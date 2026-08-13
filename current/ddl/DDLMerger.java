package ddl;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeSet;

import cardinality.DynamicDemiLog;
import map.IntObjectMap;
import parse.Parse;
import parse.Parser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;

/**
 * Loads multiple DDL files, optionally merges records sharing a TID
 * (e.g., merging mitochondrial/plastid/plasmid DDLs into their host
 * genome), sorts, renumbers, and writes a single combined DDL file.
 *
 * Source categories are detected from filenames:
 *   *mito* = mitochondrial, *plastid* = plastid, *plasmid* = plasmid,
 *   everything else = main genome.
 *
 * Usage: DDLMerger in=a.ddl.gz,b.ddl.gz out=combined.ddl.gz
 *
 * @author Ady
 * @date April 18, 2026
 */
public class DDLMerger {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		DDLMerger dm=new DDLMerger(args);
		dm.process(t);
	}

	public DDLMerger(String[] args){
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		//Genome sketching defaults to exponent 5, which trades unreachable high cardinalities for an
		//extra mantissa bit.  Set before parsing so an 'exponent=' flag still wins, and before loading
		//so a sketch file's #exponent header still wins over both.
		DynamicDemiLog.setExponent(5);

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("buckets")){
				targetBuckets=Integer.parseInt(b);
			}else if(a.equals("exponent") || a.equals("ebits")){
				DynamicDemiLog.setExponent(Integer.parseInt(b));
			}else if(a.equals("mergemito") || a.equals("mergemitochondrion")){
				mergeMito=Parse.parseBoolean(b);
			}else if(a.equals("mergeplastid")){
				mergePlastid=Parse.parseBoolean(b);
			}else if(a.equals("mergeplasmid")){
				mergePlasmid=Parse.parseBoolean(b);
			}else if(a.equals("merge") || a.equals("mergeall")){
				mergeMain=mergeMito=mergePlastid=mergePlasmid=Parse.parseBoolean(b);
			}else if(a.equals("minsize") || a.equals("originfileminsizefilter")){
				parseMinSize(b);
			}else if(a.equals("tidsout") || a.equals("outtids")){
				tidsOut=b;
			}else if(a.equals("blacklist")){
				blacklistFile=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null && Tools.looksLikeInputStream(arg)){
				//Accumulate bare input files
				inList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{
			Parser.processQuality();
			overwrite=parser.overwrite;
			out=parser.out1;
			if(parser.in1!=null){
				for(String s : parser.in1.split(",")){inList.add(s);}
			}
		}

		inFiles=inList.toArray(new String[0]);
		assert(inFiles.length>0) : "No input files specified.";
		assert(out!=null) : "No output file specified.";
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		//Load all records from all files, tagged with source category
		final ArrayList<TaggedRecord> all=new ArrayList<TaggedRecord>();
		int condensed=0;
		//Establish the authoritative k from the inputs before anything uses it (blacklist load,
		//record parsing, output header): k PROPAGATES from the input #k headers, and setFileK
		//(called per file during load) fails loud if the inputs disagree.  An explicit k= must
		//match the inputs.
		{
			int inputK=DDLLoader.peekK(inFiles[0]);
			assert(inputK>0) : KillSwitch.assertDie("No #k header in "+inFiles[0]+"; specify k= explicitly.");
			assert(k<0 || k==inputK) : KillSwitch.assertDie("Requested k="+k+" but input "+inFiles[0]+" has k="+inputK+".");
			k=inputK;
		}
		//Load the blacklist once so condense can prefer non-blacklisted kmers (bake the blacklist
		//into the condensed DB, equivalent to re-sketching with blacklist= but without reading genomes).
		if(blacklistFile!=null){DynamicDemiLog.loadBlacklist(blacklistFile, k);}
		for(String path : inFiles){
			int cat=categorize(path);
			String origin=originFromFilename(path);
			ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(path, k);
			for(DDLRecord rec : records){
				if(rec.origin==null){rec.origin=origin;}
				//Filter AT LOAD, before anything accumulates: a dropped record never enters the
				//list, so filtering makes this tool cheaper rather than more expensive (it holds
				//every surviving record in memory, and the 26% being cut is 26% never held).
				if(!keep(rec)){filtered++; continue;}
				if(targetBuckets>0 && rec.ddl.buckets>targetBuckets){
					DynamicDemiLog cddl=DynamicDemiLog.condenseBlacklisted(rec.ddl, targetBuckets);
					DDLRecord crec=new DDLRecord(cddl, rec.id, rec.taxID, rec.name);
					crec.filename=rec.filename;
					crec.bases=rec.bases;
					crec.contigs=rec.contigs;
					crec.gc=rec.gc;
					crec.origin=rec.origin;
					crec.lineage=rec.lineage;
					crec.cardinality=cddl.cardinality();
					rec=crec;
					condensed++;
				}
				all.add(new TaggedRecord(rec, cat));
			}
			if(verbose){outstream.println("Loaded "+records.size()+" records from "+path+" ("+catName(cat)+")");}
		}
		outstream.println("Loaded "+all.size()+" total records from "+inFiles.length+" files."
			+(condensed>0 ? " Condensed "+condensed+" records to "+targetBuckets+" buckets." : "")
			+(filtered>0 ? " Dropped "+filtered+" below the size filter." : ""));
		reportFilter();

		//Merge by TID where appropriate
		final IntObjectMap<TaggedRecord> tidMap=new IntObjectMap<TaggedRecord>(TaggedRecord.class);
		final IntObjectMap<long[]> gcMap=new IntObjectMap<long[]>(long[].class); //{gcBases, atBases}
		final ArrayList<DDLRecord> unmerged=new ArrayList<DDLRecord>();

		for(TaggedRecord tr : all){
			DDLRecord rec=tr.rec;
			int tid=rec.taxID;

			//Decide whether this record should merge with existing
			boolean shouldMerge=(tid>=0) && shouldMerge(tr.category);

			if(shouldMerge && tidMap.get(tid)!=null){
				//Merge into existing
				TaggedRecord existing=tidMap.get(tid);
				long[] existGC=gcMap.get(tid);
				//Weighted GC merge
				long gcA=(existGC!=null ? existGC[0] : 0);
				long atA=(existGC!=null ? existGC[1] : 0);
				long gcB=gcBases(rec);
				long atB=atBases(rec);
				existing.rec.ddl.add(rec.ddl);
				existing.rec.bases+=rec.bases;
				existing.rec.contigs+=rec.contigs;
				if(existGC==null){existGC=new long[2]; gcMap.put(tid, existGC);}
				existGC[0]=gcA+gcB;
				existGC[1]=atA+atB;
			}else if(shouldMerge){
				//First record for this TID
				tidMap.put(tid, tr);
				long gc=gcBases(rec);
				long at=atBases(rec);
				gcMap.put(tid, new long[]{gc, at});
			}else{
				//Don't merge — keep separate
				unmerged.add(rec);
			}
		}

		//Collect merged records and recompute GC
		ArrayList<DDLRecord> results=new ArrayList<DDLRecord>();
		for(int tid : tidMap.toArray()){
			DDLRecord rec=tidMap.get(tid).rec;
			long[] gc=gcMap.get(tid);
			if(gc!=null && gc[0]+gc[1]>0){rec.gc=gc[0]*1f/(gc[0]+gc[1]);}
			rec.cardinality=rec.ddl.cardinality();
			results.add(rec);
		}

		//Add unmerged records
		for(DDLRecord rec : unmerged){
			rec.cardinality=rec.ddl.cardinality();
			results.add(rec);
		}

		//Sort: taxID, then name, then cardinality
		Collections.sort(results, (a, b) -> {
			int x=a.taxID-b.taxID;
			if(x!=0){return x;}
			int y=(a.name==null ? "" : a.name).compareTo(b.name==null ? "" : b.name);
			if(y!=0){return y;}
			return Long.compare(a.cardinality, b.cardinality);
		});

		//Renumber sequentially — create new DDLRecords with final id
		ArrayList<DDLRecord> numbered=new ArrayList<DDLRecord>(results.size());
		for(int i=0; i<results.size(); i++){
			DDLRecord old=results.get(i);
			DDLRecord fresh=new DDLRecord(old.ddl, i, old.taxID, old.name);
			fresh.filename=old.filename;
			fresh.bases=old.bases;
			fresh.contigs=old.contigs;
			fresh.cardinality=old.cardinality;
			fresh.gc=old.gc;
			fresh.origin=old.origin;
			fresh.lineage=old.lineage;//[ddl/DDLMerger#001] was omitted here -> #lineage silently dropped from every merged record
			numbered.add(fresh);
		}

		final String bl=(blacklistFile!=null ? blacklistFile : DDLLoader.lastBlacklistHeader);
		DDLLoader.writeFile(numbered, out, overwrite, k, 12345L, bl);
		writeTids(numbered);

		t.stop();
		int merged=all.size()-numbered.size();
		outstream.println("Wrote "+numbered.size()+" DDL records to "+out
			+(merged>0 ? " (merged "+merged+" duplicates)" : ""));
		outstream.println("Time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Size Filter         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses minsize=key,bases+key,bases -- a per-origin minimum genome size.
	 * Keyed on the SOURCE FILE rather than on a clade, because clade membership is exactly the
	 * thing that is ambiguous here: 'plastid' vs 'chloroplast' naming differs between releases,
	 * and organelle taxonomy is confusing enough that hardcoding it would be its own bug.  The
	 * file split already encodes the domain, so the file is the honest key -- and filenames can
	 * change, so they belong in a flag rather than in the source.
	 *
	 * ⚠ THE SEPARATOR IS '+', AND THAT IS NOT A STYLE CHOICE.  Every BBTools launcher (294 of
	 * them) ends in `eval $CMD`, so bash re-parses the assembled command line AFTER the caller's
	 * quoting is gone.  Candidates were tested through a real launcher rather than assumed:
	 *   ';' and '|'  SEVER the command.  minsize='a,2000;b,100000' tidsout=x silently applied
	 *                only rule 'a', never set tidsout, printed "b,100000: command not found",
	 *                and exited 0 with plausible output.  Semicolons are unusable in ANY
	 *                BBTools argument for this reason.
	 *   '~'          EXPANDS, silently.  Bash tilde-expands after '=' in an assignment-shaped
	 *                word: minsize=~/x became minsize=/home/<user>/x.  Silent substitution is
	 *                worse than severance because nothing is printed.
	 *   ':'          survives alone, but ENABLES tilde expansion after itself (x=a,1:~root/y
	 *                became x=a,1:/root/y), and already means something else in BBTools --
	 *                vectorutils out=f1:0.9,f2:0.1 uses it for partition fractions.
	 *   '+'          inert under every test, through the real launcher.  Chosen.
	 * ';' is still accepted for direct `java ddl.DDLMerger` invocation, where no eval intervenes.
	 * The fully eval-proof form needs no separator at all -- repeat the flag, entries accumulate:
	 *   minsize=refseq_bacteria.fa.gz,100000 minsize=refseq_invertebrate.fa.gz,1000000
	 * @param s '+'- or ';'-delimited key,bases pairs
	 */
	private void parseMinSize(String s){
		if(s==null || s.length()<1){return;}
		for(String pair : s.split("[+;]")){
			if(pair.trim().length()<1){continue;}
			int comma=pair.lastIndexOf(',');
			assert(comma>0) : "minsize entries must be key,bases -- got '"+pair+"'";
			String key=normalizeOrigin(pair.substring(0, comma).trim());
			long min=Parse.parseKMG(pair.substring(comma+1).trim());
			assert(key.length()>0) : "Empty key in minsize entry '"+pair+"'";
			minSize.put(key, min);
			ruleHits.put(key, 0L);
			ruleDrops.put(key, 0L);
		}
	}

	/** Strips directory, compression and format suffixes so a rule key and a record origin can
	 * be compared regardless of which layer the name came from (source fasta vs sketch file). */
	static String normalizeOrigin(String s){
		if(s==null){return "";}
		String name=new File(s).getName().toLowerCase();
		for(String ext : new String[]{".gz", ".bz2", ".zip"}){
			if(name.endsWith(ext)){name=name.substring(0, name.length()-ext.length());}
		}
		for(String ext : new String[]{".tsv", ".ddl", "_ddl", ".fna", ".fasta", ".fa", ".txt"}){
			if(name.endsWith(ext)){name=name.substring(0, name.length()-ext.length());}
		}
		//Collapse every run of . _ - to a single _ so a rule key and a record value match regardless
		//of which separator each side used: the #file header is dot-separated (refseq.bacteria.fna.gz
		//-> refseq_bacteria) while a hand-typed key or a 32k #origin (bacteria_b64k) is underscored.
		//Both sides pass through this one function, so keys and values normalize identically.
		name=name.replaceAll("[._-]+", "_").replaceAll("^_+|_+$", "");
		return name;
	}

	/** Rule key matching this record, or null if none.  Matches the record's SOURCE FILE (#file
	 * header) first and its origin second -- either hitting is enough.  Filename is the load-bearing
	 * field: in a merged sketch file the per-record #file always survives, while #origin is either
	 * absent (the 4k file, where rec.origin then gets overwritten with the input filename) or a
	 * clade+build tag (the 32k file, e.g. bacteria_b64k).  So origin alone would match nothing. */
	private String ruleFor(DDLRecord rec){
		String key=ruleForValue(rec.filename);
		if(key==null){key=ruleForValue(rec.origin);}
		return key;
	}

	/** Rule key matching a single origin/filename string, or null.  Exact match wins; otherwise a
	 * containment match either way, so 'refseq_bacteria' matches 'refseq.bacteria.genomic' variants. */
	private String ruleForValue(String value){
		final String o=normalizeOrigin(value);
		if(o.length()<1){return null;}
		if(minSize.containsKey(o)){return o;}
		for(String key : minSize.keySet()){
			if(o.contains(key) || key.contains(o)){return key;}
		}
		return null;
	}

	/**
	 * Whether this record survives the size filter.
	 * An origin named by NO rule is kept whole.  That default is what makes organelles, plasmids
	 * and viruses exempt without any clade logic: you simply do not list those files.  It also
	 * makes the failure mode safe -- a new or renamed input file passes through visibly (and is
	 * reported below) rather than silently vanishing from the database.
	 */
	private boolean keep(DDLRecord rec){
		//Report against the SOURCE FILE when present (the clade-specific field); fall back to origin.
		//On the 4k file rec.origin is the input filename for every record, so the filename is what
		//makes the "unfiltered (kept whole)" line name the actual exempt clades (plasmid, viral, ...).
		final String label=(rec.filename!=null && rec.filename.length()>0 ? rec.filename
			: (rec.origin!=null ? rec.origin : ""));
		originSeen.merge(label, 1L, Long::sum);
		String key=ruleFor(rec);
		if(key==null){
			unruledOrigins.add(label);
			return true;
		}
		ruleHits.merge(key, 1L, Long::sum);
		if(rec.bases>=minSize.get(key)){return true;}
		ruleDrops.merge(key, 1L, Long::sum);
		originDropped.merge(label, 1L, Long::sum);
		return false;
	}

	/** Reports what the filter did, including the cases that are silent failures otherwise:
	 * a rule that matched nothing (renamed or missing input) and an origin no rule covered. */
	private void reportFilter(){
		if(minSize.isEmpty()){return;}
		outstream.println("\nSize filter (per source file):");
		for(Map.Entry<String, Long> e : minSize.entrySet()){
			final String key=e.getKey();
			final long seen=ruleHits.get(key), dropped=ruleDrops.get(key);
			outstream.println("  "+key+"\tmin="+e.getValue()+"\tmatched="+seen
				+"\tdropped="+dropped+"\tkept="+(seen-dropped)
				+(seen==0 ? "\t<-- WARNING: matched NO records; renamed or missing input?" : ""));
		}
		if(!unruledOrigins.isEmpty()){
			ArrayList<String> list=new ArrayList<String>(unruledOrigins);
			Collections.sort(list);
			outstream.println("  Unfiltered (no rule, kept whole): "+list.size()+" origin(s): "
				+String.join(", ", list));
		}
	}

	/** Writes the surviving taxIDs, sorted and unique, one per line -- the whitelist for the
	 * downstream clade filter, so the two databases can be made congruent without a second scan. */
	private void writeTids(ArrayList<DDLRecord> records){
		if(tidsOut==null){return;}
		TreeSet<Integer> tids=new TreeSet<Integer>();
		for(DDLRecord rec : records){
			if(rec.taxID>0){tids.add(rec.taxID);}
		}
		FileFormat ff=FileFormat.testOutput(tidsOut, FileFormat.TXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder(1<<16);
		for(int tid : tids){
			bb.append(tid).nl();
			if(bb.length()>=(1<<16)){bsw.print(bb); bb.clear();}
		}
		if(bb.length()>0){bsw.print(bb);}
		bsw.poisonAndWait();
		outstream.println("Wrote "+tids.size()+" distinct taxIDs to "+tidsOut);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Categorize a file by its name. */
	static int categorize(String path){
		String name=new File(path).getName().toLowerCase();
		if(name.contains("mito")){return CAT_MITO;}
		if(name.contains("plastid")){return CAT_PLASTID;}
		if(name.contains("plasmid")){return CAT_PLASMID;}
		return CAT_MAIN;
	}

	static String originFromFilename(String path){
		String name=new File(path).getName();
		if(name.endsWith(".gz")){name=name.substring(0, name.length()-3);}
		if(name.endsWith(".tsv")){name=name.substring(0, name.length()-4);}
		if(name.endsWith("_ddl")){name=name.substring(0, name.length()-4);}
		if(name.endsWith(".ddl")){name=name.substring(0, name.length()-4);}
		return name;
	}

	/** Whether records of this category should be merge-eligible. */
	boolean shouldMerge(int category){
		switch(category){
			case CAT_MITO: return mergeMito;
			case CAT_PLASTID: return mergePlastid;
			case CAT_PLASMID: return mergePlasmid;
			default: return mergeMain;
		}
	}

	static String catName(int cat){
		switch(cat){
			case CAT_MITO: return "mitochondrion";
			case CAT_PLASTID: return "plastid";
			case CAT_PLASMID: return "plasmid";
			default: return "main";
		}
	}

	/** Estimate GC base count from gc fraction and total bases. */
	static long gcBases(DDLRecord rec){
		if(rec.gc<0 || rec.bases<1){return 0;}
		return (long)(rec.gc*rec.bases);
	}

	/** Estimate AT base count from gc fraction and total bases. */
	static long atBases(DDLRecord rec){
		if(rec.gc<0 || rec.bases<1){return 0;}
		return (long)((1f-rec.gc)*rec.bases);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	static class TaggedRecord {
		TaggedRecord(DDLRecord rec_, int category_){
			rec=rec_;
			category=category_;
		}
		final DDLRecord rec;
		final int category;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<String> inList=new ArrayList<String>();
	private String[] inFiles;
	private String out;
	private int k=-1;//-1 = adopt from input #k headers (all inputs must agree); k= forces/validates it
	private int targetBuckets=-1;
	private boolean overwrite=false;
	private boolean verbose=false;
	private String blacklistFile=null;
	private String tidsOut=null;

	/** Per-source-file minimum genome size, from minsize=.  Insertion-ordered so the report
	 * lists rules in the order they were given. */
	private final LinkedHashMap<String, Long> minSize=new LinkedHashMap<String, Long>();
	/** Records each rule matched, and of those, how many it dropped. */
	private final LinkedHashMap<String, Long> ruleHits=new LinkedHashMap<String, Long>();
	private final LinkedHashMap<String, Long> ruleDrops=new LinkedHashMap<String, Long>();
	/** Origins seen, and origins covered by no rule (kept whole, and reported). */
	private final LinkedHashMap<String, Long> originSeen=new LinkedHashMap<String, Long>();
	private final LinkedHashMap<String, Long> originDropped=new LinkedHashMap<String, Long>();
	private final TreeSet<String> unruledOrigins=new TreeSet<String>();
	/** Records removed by the size filter. */
	private long filtered=0;

	private boolean mergeMain=true;
	private boolean mergeMito=true;
	private boolean mergePlastid=true;
	private boolean mergePlasmid=true;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int CAT_MAIN=0, CAT_MITO=1, CAT_PLASTID=2, CAT_PLASMID=3;
	private static final PrintStream outstream=System.err;
}
