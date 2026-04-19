package ddl;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import cardinality.DynamicDemiLog;
import map.IntObjectMap;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import fileIO.ReadWrite;

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

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("mergemito") || a.equals("mergemitochondrion")){
				mergeMito=Parse.parseBoolean(b);
			}else if(a.equals("mergeplastid")){
				mergePlastid=Parse.parseBoolean(b);
			}else if(a.equals("mergeplasmid")){
				mergePlasmid=Parse.parseBoolean(b);
			}else if(a.equals("merge")){
				mergeMito=mergePlastid=mergePlasmid=Parse.parseBoolean(b);
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
		for(String path : inFiles){
			int cat=categorize(path);
			ArrayList<DDLRecord> records=DDLLoader.loadFile(path, k);
			for(DDLRecord rec : records){
				all.add(new TaggedRecord(rec, cat));
			}
			if(verbose){outstream.println("Loaded "+records.size()+" records from "+path+" ("+catName(cat)+")");}
		}
		outstream.println("Loaded "+all.size()+" total records from "+inFiles.length+" files.");

		//Merge by TID where appropriate
		final IntObjectMap<TaggedRecord> tidMap=new IntObjectMap<TaggedRecord>();
		final IntObjectMap<long[]> gcMap=new IntObjectMap<long[]>(); //{gcBases, atBases}
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
			numbered.add(fresh);
		}

		DDLLoader.writeFile(numbered, out, overwrite, k, 12345L);

		t.stop();
		int merged=all.size()-numbered.size();
		outstream.println("Wrote "+numbered.size()+" DDL records to "+out
			+(merged>0 ? " (merged "+merged+" duplicates)" : ""));
		outstream.println("Time: \t"+t);
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

	/** Whether records of this category should be merge-eligible. */
	boolean shouldMerge(int category){
		switch(category){
			case CAT_MITO: return mergeMito;
			case CAT_PLASTID: return mergePlastid;
			case CAT_PLASMID: return mergePlasmid;
			default: return true; //Main genome always merge-eligible
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
	private int k=31;
	private boolean overwrite=false;
	private boolean verbose=false;

	private boolean mergeMito=true;
	private boolean mergePlastid=true;
	private boolean mergePlasmid=true;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int CAT_MAIN=0, CAT_MITO=1, CAT_PLASTID=2, CAT_PLASMID=3;
	private static final PrintStream outstream=System.err;
}
