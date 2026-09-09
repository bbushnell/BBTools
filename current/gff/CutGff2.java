package gff;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import aligner.Alignment;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.PGMTools;
import prok.ProkObject;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;
import tax.GiToTaxid;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Streaming clone of CutGff (Brian, 2026-09-09: "if it loads both into memory and thus is
 * inefficient, it would probably be best to clone it as CutGff2 and make that one stream the
 * sequences while holding the Gff in memory"). Holds the GFF in memory grouped by seqid and
 * STREAMS the FASTA one batch at a time through the Streamer/Writer interface
 * (template.A_SampleStreamerMT pattern, per Brian) — memory is O(GFF + largest batch), not
 * O(genome), so whole-clade extraction needs a few GB instead of hundreds
 * (measured motivator: plant tRNA extraction held a ~178Gbp scaffold payload for 12Mbp of
 * output features, MaxRSS ~358GiB).
 *
 * <p>Per-feature emission semantics (filters, flank, GC/filename/tid header decoration,
 * strand handling, maxNs, alignRibo validation) are CLONED from CutGff.processLines and run
 * per-contig — for any feature whose contig appears in the FASTA the output record is
 * identical to CutGff's. Deliberate deltas, all loud or documented:
 * <ul>
 * <li>Output order is FASTA-major (features of one contig together, GFF order within a
 *   contig), not GFF-major.</li>
 * <li>Features whose seqid never streams by are counted and cause a RuntimeException at end
 *   of input (CutGff asserts per-line instead); allowmissingseqids=t downgrades to a loud
 *   stderr report — used when the annotation legitimately covers records absent from the
 *   sequence file, and the count is part of the caller's receipt.</li>
 * <li>Duplicate contig ids: FIRST occurrence consumes the features (CutGff's map.put was
 *   last-wins); later occurrences are counted and reported.</li>
 * <li>Sequence headers of the form tid|taxid|ACC fall back to the bare ACC (text after the
 *   second pipe) for GFF-seqid lookup, so clade bulk FASTAs work without a rewrite pass.</li>
 * <li>UNSUPPORTED (throws, with the reason — use cutgff.sh instead): pickbest, oneperfile
 *   (file-scoped selection over small curation inputs, where whole-file load is fine),
 *   multiple fna/gff pairs, and renamebytaxid with network tax modes (accession/gi/header
 *   would make one server query per contig when streaming); renamebytaxid with tid|-style
 *   headers or taxmode=taxid is supported (local parsing).</li>
 * </ul>
 *
 * @author Brian Bushnell
 * @author G11
 * @date 2026-09-09
 */
public class CutGff2 implements Accumulator<CutGff2.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		CutGff2 x=new CutGff2(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CutGff2(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Shared.TRIM_READ_DESCRIPTION=Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;
		Read.VALIDATE_IN_CONSTRUCTOR=true;
		GffLine.parseAttributes=true;

		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			workers=parser.workers();
			threadsIn=parser.threadsIn;
			threadsOut=parser.threadsOut;
			maxReads=parser.maxReads;
			out=parser.out1;
		}

		if(alignRibo){
			ProkObject.loadConsensusSequenceFromFile(false, false);
		}

		if(pickBest || onePerFile){
			throw new RuntimeException("pickbest/oneperfile are file-scoped selection modes over small "
				+"curation inputs where CutGff's whole-file load is appropriate; CutGff2 does not "
				+"support them. Use cutgff.sh.");
		}
		if(fnaList.size()!=1 || gffList.size()!=1){
			throw new RuntimeException("CutGff2 processes exactly one fna/gff pair per invocation "
				+"(got "+fnaList.size()+"/"+gffList.size()+"); loop invocations or use cutgff.sh.");
		}
		if(renameByTaxID && taxMode!=TAXID_MODE){
			throw new RuntimeException("renamebytaxid with taxmode accession/gi/header would query the "
				+"tax server once per contig when streaming; CutGff2 supports it only for tid|-style "
				+"headers or taxmode=taxid (local parsing). Use cutgff.sh for network modes.");
		}

		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffin=FileFormat.testInput(fnaList.get(0), FileFormat.FA, null, true, true);
		ffout=FileFormat.testOutput(out, FileFormat.FA, null, true, overwrite, append, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse arguments from the command line */
	private Parser parse(String[] args){

		Parser parser=new Parser();
		parser.overwrite=overwrite;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null) : "Bad parameter: "+arg;
				Tools.addFiles(b, fnaList);
			}else if(a.equals("gff") || a.equals("ingff") || a.equals("gffin")){
				assert(b!=null) : "Bad parameter: "+arg;
				Tools.addFiles(b, gffList);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("alignribo") || a.equals("align")){
				alignRibo=Parse.parseBoolean(b);
			}else if(a.equals("adjustendpoints")){
				adjustEndpoints=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("slop16s") || a.equalsIgnoreCase("16sslop") || a.equalsIgnoreCase("ssuslop")){
				ssuSlop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("slop23s") || a.equalsIgnoreCase("23sslop") || a.equalsIgnoreCase("lsuslop")){
				lsuSlop=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxns") || a.equalsIgnoreCase("maxundefined")){
				maxNs=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("maxnrate") || a.equalsIgnoreCase("maxnfraction")){
				maxNFraction=Integer.parseInt(b);
			}else if(a.equals("invert")){
				invert=Parse.parseBoolean(b);
			}else if(a.equals("type") || a.equals("types")){
				types=b;
			}else if(a.equals("attributes") || a.equals("requiredattributes")){
				requiredAttributes=b.split(",");
			}else if(a.equals("banattributes") || a.equals("bannedattributes")){
				bannedAttributes=b.split(",");
			}else if(a.equals("banpartial")){
				banPartial=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("allowmissingseqids") || a.equalsIgnoreCase("allowmissing")){
				allowMissingSeqids=Parse.parseBoolean(b);
			}

			else if(a.equalsIgnoreCase("renameByTaxID")){
				renameByTaxID=Parse.parseBoolean(b);
			}else if(a.equals("taxmode")){
				if("accession".equalsIgnoreCase(b)){
					taxMode=ACCESSION_MODE;
				}else if("header".equalsIgnoreCase(b)){
					taxMode=HEADER_MODE;
				}else if("gi".equalsIgnoreCase(b)){
					taxMode=GI_MODE;
				}else if("taxid".equalsIgnoreCase(b)){
					taxMode=TAXID_MODE;
				}else{
					assert(false) : "Bad tax mode: "+b;
				}
			}else if(a.equals("requirepresent")){
				requirePresent=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("onePerFile")){
				onePerFile=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("pickBest") || a.equalsIgnoreCase("findBest") || a.equalsIgnoreCase("keepBest")){
				pickBest=Parse.parseBoolean(b);
			}

			else if(a.equals("minlen")){
				minLen=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxLen=Integer.parseInt(b);
			}else if(a.equals("flank") || a.equals("pad")){
				flank=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("gccontig")){
				appendGC=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("filename")){
				appendFilename=Parse.parseBoolean(b);
			}

			else if(ProkObject.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		ArrayList<String> banned=new ArrayList<String>();
		if(banPartial){banned.add("partial=true");}
		if(bannedAttributes!=null){
			for(String s : bannedAttributes){banned.add(s);}
		}
		bannedAttributes=banned.isEmpty() ? null : banned.toArray(new String[0]);

		if(gffList.isEmpty()){
			for(String s : fnaList){
				String prefix=ReadWrite.stripExtension(s);
				String gff=prefix+".gff";
				File f=new File(gff);
				if(!f.exists()){
					String gz=gff+".gz";
					f=new File(gz);
					assert(f.exists() && f.canRead()) : "Can't read file "+gff;
					gff=gz;
				}
				gffList.add(gff);
			}
		}
		assert(gffList.size()==fnaList.size()) : "Number of fna and gff files do not match: "+fnaList.size()+", "+gffList.size();
		return parser;
	}

	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		gffList=Tools.fixExtension(gffList);
	}

	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(gffList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");
		}
		if(out!=null && !Tools.testOutputFiles(overwrite, append, false, out)){
			throw new RuntimeException("\nCan't write to output file "+out+"\n");
		}
	}

	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public void process(Timer t){
		//Load and group the annotation first; this is the part that stays in memory.
		{
			ArrayList<GffLine> lines=GffLine.loadGffFile(gffList.get(0), types, false);
			featuresLoaded=lines.size();
			for(GffLine gline : lines){
				ArrayList<GffLine> sub=bySeqid.get(gline.seqid);
				if(sub==null){bySeqid.put(gline.seqid, sub=new ArrayList<GffLine>(4));}
				sub.add(gline);
			}
		}
		outstream.println("Loaded "+featuresLoaded+" features of type "+types+" on "+bySeqid.size()+" seqids.");

		//Default 1 worker (Brian, 2026-09-09: template-MT shape is fine but "default to running
		//it with 1 worker" — keeps resources low and output fully ordered/deterministic; the
		//pipeline is I/O-bound so extra workers buy little). workers= overrides.
		final int threads=Tools.max(1, workers>0 ? workers : 1);
		Read.VALIDATE_IN_CONSTRUCTOR=(threads<2 && threadsIn<2);

		Streamer st=StreamerFactory.makeStreamer(ffin, null, true, maxReads, false, true, threadsIn);
		st.start();

		Writer fw=(ffout==null ? null : WriterFactory.makeWriter(ffout, null, threadsOut, null, false));
		if(fw!=null){fw.start();}

		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(st, fw, i));
		}

		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;

		if(fw!=null){
			if(errorState){
				//A dead processing worker may have left an ordered batch-id gap, and the
				//graceful poison/drain path could wait forever on the missing id — exactly
				//the case Writer.finishError() exists for (non-blocking abandon; see its
				//javadoc). Qiqi review, 2026-09-09. The throw below still fires.
				fw.finishError();
			}else{
				errorState|=fw.poisonAndWait();//@return errorState (true=error), per Writer's javadoc
			}
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}
		errorState|=ReadWrite.closeStream(st);

		//End-of-input accounting: everything left in bySeqid never met its contig.
		if(!bySeqid.isEmpty()){
			long leftover=0;
			for(ArrayList<GffLine> sub : bySeqid.values()){leftover+=sub.size();}
			StringBuilder sb=new StringBuilder();
			int shown=0;
			for(String s : bySeqid.keySet()){
				if(shown>=5){sb.append(", ..."); break;}
				if(shown>0){sb.append(", ");}
				sb.append(s); shown++;
			}
			String msg="Features whose seqid was absent from the sequence input: "+leftover
				+" (on "+bySeqid.size()+" seqids, e.g. "+sb+")";
			if(allowMissingSeqids){
				outstream.println("WARNING: "+msg);
				outstream.println("missing_seqid_features="+leftover+" missing_seqids="+bySeqid.size());
			}else{
				throw new RuntimeException(msg+"; the annotation does not match the sequence file. "
					+"Pass allowmissingseqids=t only if this is expected and the count belongs in "
					+"your receipt.");
			}
		}
		if(duplicateContigs.get()>0){
			outstream.println("duplicate_feature_bearing_contig_ids_skipped="+duplicateContigs.get()
				+" (first occurrence consumed the features)");
		}

		t.stop();
		if(ffout!=null){outstream.println("Wrote "+out);}
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		if(alignRibo){
			outstream.println(Tools.number("Flipped:           ", flipped.get(), 8));
			outstream.println(Tools.number("Failed Alignment:  ", failed.get(), 8));
		}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the id up to (not including) the first space or tab; the id itself if none. */
	private static String firstToken(String id){
		if(id==null){return null;}
		for(int i=0; i<id.length(); i++){
			final char c=id.charAt(i);
			if(c==' ' || c=='\t'){return id.substring(0, i);}
		}
		return id;
	}

	/** Claims the feature sublist for a streamed contig, or null if it has none (or a
	 * same-id contig already claimed them — counted). Thread-safe; each sublist is handed
	 * out exactly once. Headers of the form tid|taxid|ACC also try bare ACC. */
	private ArrayList<GffLine> claimFeatures(String id){
		String key=id;
		synchronized(bySeqid){
			ArrayList<GffLine> sub=bySeqid.remove(key);
			if(sub==null && id.startsWith("tid|")){
				int pipe2=id.indexOf('|', 4);
				if(pipe2>0 && pipe2+1<id.length()){
					key=id.substring(pipe2+1);
					sub=bySeqid.remove(key);
				}
			}
			if(sub!=null){consumedIds.add(key);}
			else if(consumedIds.contains(key)){duplicateContigs.incrementAndGet();}
			return sub;
		}
	}

	/** Sets r.obj (taxid) and normalizes r.id, mirroring CutGff.renameByTaxID's local
	 * (non-network) paths; the constructor rejects network tax modes for streaming. */
	private void renameByTaxID(Read r){
		if(r.id!=null && r.id.startsWith("tid|")){
			r.obj=TaxTree.parseHeaderStatic(r.id);
		}else{
			int id=GiToTaxid.parseTaxidNumber(r.id, '|');
			assert(id>=0 || !requirePresent) : "Can't find taxID for header: "+id+", "+r.name();
			r.obj=id;
			r.id="tid|"+id+"|"+id;
		}
	}

	/**
	 * Checks if a GFF line passes length and attribute filter criteria.
	 * (Cloned from CutGff.)
	 */
	private boolean hasAttributes(GffLine gline){
		if(gline.attributes==null){return false;}
		int len=gline.length();
		if(len<minLen || len>maxLen){return false;}
		if(hasAttributes(gline, bannedAttributes)){return false;}
		return requiredAttributes==null || hasAttributes(gline, requiredAttributes);
	}

	/** Tests if a GFF line contains any of the specified attributes. (Cloned from CutGff.) */
	private static boolean hasAttributes(GffLine gline, String[] attributes){
		if(attributes==null){return false;}
		for(String s : attributes){
			if(gline.attributes.contains(s)){
				return true;
			}
		}
		return false;
	}

	/**
	 * Processes one contig's GFF lines and extracts or masks corresponding sequences.
	 * CLONED from CutGff.processLines with the map lookup replaced by the single already-
	 * claimed scaffold — per-feature emission is otherwise identical. gcCache is per-thread
	 * (a contig is processed by exactly one thread).
	 */
	private ArrayList<Read> processLines(ArrayList<GffLine> lines, Read scaf, boolean invertSelection,
			HashMap<String, Float> gcCache, String sourceBasename){
		ArrayList<Read> list=null;
		for(GffLine gline : lines){
			if(hasAttributes(gline)){
				assert(scaf!=null) : "Can't find "+gline.seqid;

				boolean pass=true;
				Float identity=null;
				if(alignRibo && gline.inbounds(scaf.length())){
					int type=gline.prokType();
					identity=align(gline, scaf.bases, type);
					if(identity==null){pass=false;}
				}

				if(pass){
					final int start=gline.start-1;
					final int stop=gline.stop-1;

					if(invertSelection){
						byte[] bases=scaf.bases;
						for(int i=start; i<=stop; i++){
							if(i>=0 && i<bases.length){
								bases[i]='N';
							}
						}
					}else{
						if(start>=0 && stop<scaf.length()){
							final int extStart=Tools.max(0, start-flank);
							final int extStop=Tools.min(scaf.length()-1, stop+flank);
							final int leftFlank=start-extStart;//genomic-left flank actually applied (may be < flank at a contig edge)
							final int rightFlank=extStop-stop;//genomic-right flank actually applied
							assert(leftFlank>=0 && leftFlank<=flank) : leftFlank+", "+flank;
							assert(rightFlank>=0 && rightFlank<=flank) : rightFlank+", "+flank;

							String id=gline.attributes;
							if(renameByTaxID){
								id="tid|"+scaf.obj+"|"+id;
							}
							if(flank>0){
								//lflank/rflank describe the OUTPUT record's 5'->3' orientation: on a minus-strand
								//feature the genomic-right flank becomes the record's 5' leader after the
								//reverseComplement below, so left/right swap. Plus strand: no swap needed.
								final int lflank=(gline.strand==GffLine.MINUS ? rightFlank : leftFlank);
								final int rflank=(gline.strand==GffLine.MINUS ? leftFlank : rightFlank);
								id=id+" lflank="+lflank+" rflank="+rflank;
							}
							if(appendGC){
								id=id+" contig_gc="+String.format("%.4f", contigGC(scaf, gcCache));
							}
							if(appendFilename){
								id=id+" source="+sourceBasename;
							}
							Read r=new Read(Arrays.copyOfRange(scaf.bases, extStart, extStop+1), null, id, 1);
							r.obj=identity;

							assert(!r.containsLowercase()) : r.toFasta()+"\n"
							+ "validated="+r.validated()+", scaf.validated="+scaf.validated()+", tuc="+Read.TO_UPPER_CASE+", vic="+Read.VALIDATE_IN_CONSTRUCTOR;
							if(maxNs>=0 || maxNFraction>=0){
								long allowed=Tools.min(maxNs>=0 ? maxNs : r.length(), (long)(r.length()*(maxNFraction>=0 ? maxNFraction : 1)));
								if(r.countUndefined()>allowed){r=null;}
							}

							if(r!=null){
								if(gline.strand==1){r.reverseComplement();}
								if(list==null){list=new ArrayList<Read>(8);}
								list.add(r);
							}
						}
					}
				}
			}
		}
		return list;
	}

	/** GC fraction of a scaffold's FULL bases, cached per contig id. (Cloned from CutGff.) */
	private static float contigGC(Read scaf, HashMap<String, Float> cache){
		Float f=cache.get(scaf.id);
		if(f==null){
			f=Tools.calcGC(scaf.bases);
			cache.put(scaf.id, f);
		}
		return f;
	}

	/** Consensus-alignment validation of an rRNA feature. (Cloned from CutGff.) */
	private Float align(GffLine gline, byte[] scaf, int type){
		Read[] consensusReads=ProkObject.consensusReads(type);
		if(consensusReads==null || consensusReads.length==0){
			assert(false) : type+"\n"+gline.toString();
			return null;
		}
		byte[] universal=consensusReads[0].bases;
		float minIdentity=ProkObject.minID(type)*ID_MULT;
		if(universal==null){assert(false); return 1F;}

		int start=gline.start-1;
		int stop=gline.stop-1;
		assert(start<=stop) : start+", "+stop+", "+scaf.length;
		assert(start>=0 && start<scaf.length) : start+", "+stop+", "+scaf.length;
		assert(stop>=0 && stop<scaf.length) : start+", "+stop+", "+scaf.length;
		final int a=Tools.max(0, start);
		final int b=Tools.min(scaf.length-1, stop);

		byte[] ref=Arrays.copyOfRange(scaf, a, b+1);
		Read r=new Read(ref, null, 0);
		if(gline.strand==GffLine.MINUS){r.reverseComplement();}

		Alignment plus=new Alignment(r);
		plus.align(universal);

		r.reverseComplement();
		Alignment minus=new Alignment(r);
		minus.align(universal);

		Alignment best=null;
		if(plus.id>=minus.id){
			best=plus;
		}else{
			best=minus;
			if(minus.id>=minIdentity){
				if(verbose) {System.err.println("Flipped: "+plus.id+" \t"+minus.id+"");}
				flipped.incrementAndGet();
				gline.strand=Shared.MINUS;
			}
		}
		if(best.id>=minIdentity){
			return best.id;
		}else{
			if(verbose) {System.err.println("Failed alignment: "+plus.id+" \t"+minus.id);}
			failed.incrementAndGet();
			return null;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			errorState|=(!pt.success);
		}
	}

	@Override
	public final boolean success(){return !errorState;}

	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	class ProcessThread extends Thread {

		ProcessThread(Streamer st_, Writer fw_, int tid_){
			st=st_;
			fw=fw_;
			tid=tid_;
			gcCacheT=appendGC ? new HashMap<String, Float>() : null;
			sourceBasenameT=appendFilename ? new File(fnaList.get(0)).getName() : null;
		}

		@Override
		public void run(){
			ListNum<Read> ln=st.nextList();
			while(ln!=null && ln.size()>0){
				//One output list per input batch, re-added under the batch's own id so an
				//ordered Writer sees every id exactly once (the A_SampleStreamerMT contract).
				ArrayList<Read> outList=new ArrayList<Read>();
				for(Read r : ln){
					readsProcessedT++;
					basesProcessedT+=r.length();
					//FastaStreamer ignores Shared.TRIM_READ_DESCRIPTION (see the TODO there),
					//so r.id may carry the description. Trim to the first token HERE so the
					//GFF lookup, tid| parsing, GC cache key and invert emission all see the
					//same trimmed id CutGff's legacy input path produced.
					r.id=firstToken(r.id);
					ArrayList<GffLine> sub=claimFeatures(r.id);
					if(sub!=null){
						if(renameByTaxID){renameByTaxID(r);}//After claim: lookup uses the original id
						ArrayList<Read> cut=processLines(sub, r, invert, gcCacheT, sourceBasenameT);
						if(cut!=null){outList.addAll(cut);}
					}
					if(invert){outList.add(r);}//invert emits every (possibly masked) scaffold
				}
				if(fw!=null){fw.add(outList, ln.id);}
				ln=st.nextList();
			}
			success=true;
		}

		long readsProcessedT=0;
		long basesProcessedT=0;
		boolean success=false;

		private final HashMap<String, Float> gcCacheT;
		private final String sourceBasenameT;
		private final Streamer st;
		private final Writer fw;
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** List of input FASTA file paths (exactly one supported) */
	private ArrayList<String> fnaList=new ArrayList<String>();
	/** List of input GFF annotation file paths (exactly one supported) */
	private ArrayList<String> gffList=new ArrayList<String>();
	/** Output file path for extracted sequences */
	private String out=null;
	/** Comma-separated list of feature types to extract (default: "CDS") */
	private String types="CDS";
	/** If true, mask regions instead of extracting them (all scaffolds are emitted) */
	private boolean invert=false;
	/** If true, exclude features marked as partial */
	private boolean banPartial=true;
	/** Minimum feature length to retain */
	private int minLen=1;
	/** Maximum feature length to retain */
	private int maxLen=Integer.MAX_VALUE;
	/** Bases of genomic-space flank to add to each side of an extracted feature (0=off) */
	private int flank=0;
	/** If true, append contig_gc=X.XXXX (GC fraction of the FULL source contig) to headers */
	private boolean appendGC=false;
	/** If true, append source=<basename of input fna> to headers */
	private boolean appendFilename=false;
	/** If true, features referencing seqids absent from the FASTA are reported, not fatal */
	private boolean allowMissingSeqids=false;

	/** Array of attributes that must be present for feature retention */
	private String[] requiredAttributes;
	/** Array of attributes that exclude features when present */
	private String[] bannedAttributes;

	/*--------------------------------------------------------------*/

	/** GFF features grouped by seqid; sublists are claimed (removed) as contigs stream by,
	 * so leftovers at EOF are exactly the never-matched features. Guarded by its own lock. */
	private final HashMap<String, ArrayList<GffLine>> bySeqid=new HashMap<String, ArrayList<GffLine>>();
	/** Lookup keys already claimed — detects duplicate feature-bearing contig ids. */
	private final HashSet<String> consumedIds=new HashSet<String>();
	/** Count of later same-id contigs whose features were already claimed (first-wins). */
	private final AtomicLong duplicateContigs=new AtomicLong(0);
	/** Features loaded from the GFF (post type-filter). */
	private long featuresLoaded=0;

	/** If true, rename sequences using taxonomic identifiers (local modes only; see javadoc) */
	private boolean renameByTaxID=false;
	/** Taxonomic parsing mode; only TAXID_MODE (or tid|-style headers) supported here */
	private int taxMode=ACCESSION_MODE;
	/** If true, require taxonomic ID to be found */
	private boolean requirePresent=false;
	/** If true, validate rRNA features using consensus alignment */
	private boolean alignRibo=false;
	/** Parsed for CutGff flag compatibility; endpoint adjustment is inactive there too */
	private boolean adjustEndpoints=false;
	/** Unsupported; constructor throws if set (see javadoc) */
	private boolean onePerFile=false;
	/** Unsupported; constructor throws if set (see javadoc) */
	private boolean pickBest=false;
	/** Tolerance for small subunit rRNA endpoint adjustment (flag compatibility) */
	private int ssuSlop=999;
	/** Tolerance for large subunit rRNA endpoint adjustment (flag compatibility) */
	private int lsuSlop=999;

	/** Identity multiplier for alignment scoring */
	private float ID_MULT=0.96f;

	/** Maximum number of N bases allowed in extracted sequences */
	private int maxNs=-1;
	/** Maximum fraction of N bases allowed in extracted sequences */
	private double maxNFraction=-1;

	private static int ACCESSION_MODE=0, GI_MODE=1, HEADER_MODE=2, TAXID_MODE=3;

	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of reads retained (from the Writer) */
	protected long readsOut=0;
	/** Number of bases retained (from the Writer) */
	protected long basesOut=0;

	/** Thread-safe counter for features with corrected strand orientation */
	protected AtomicLong flipped=new AtomicLong(0);
	/** Thread-safe counter for features that failed alignment validation */
	protected AtomicLong failed=new AtomicLong(0);

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	/** Worker thread count (workers= flag; default min(4, machine threads)) */
	private int workers=-1;
	/** Input decompression/parsing threads (threadsin= flag) */
	private int threadsIn=-1;
	/** Output compression threads (threadsout= flag) */
	private int threadsOut=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input file format specification */
	private FileFormat ffin;
	/** Output file format specification */
	private FileFormat ffout;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Output stream for logging and status messages */
	private PrintStream outstream=System.err;
	/** Global verbose mode flag for detailed logging */
	public static boolean verbose=false;
	/** Flag indicating whether processing encountered errors */
	public boolean errorState=false;
	/** Flag allowing overwrite of existing output files */
	private boolean overwrite=true;
	/** Flag controlling whether to append to existing output files */
	private boolean append=false;

}
