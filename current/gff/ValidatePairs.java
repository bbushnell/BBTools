package gff;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Validates directories of paired genome (.fna) + annotation (.gff) files.
 * For each basename that has both an .fna(.gz) and a .gff(.gz), checks:
 *   (1) every GFF feature line has exactly 9 tab-separated fields -- a short field
 *       count flags a truncated/partial record (a common download-truncation signature);
 *   (2) the set of sequence IDs the GFF references (feature seqids + ##sequence-region
 *       pragmas) matches the set of sequence headers in the FNA, in BOTH directions.
 * Also flags unpaired files (fna without gff or vice versa) and empty files.
 * Multithreaded across pairs, so a directory of tens of thousands of pairs validates fast.
 *
 * This complements a gzip-integrity check: gzip -t catches truncated/corrupt compression,
 * but a file can decompress cleanly and still carry partial GFF records or an FNA/GFF
 * version mismatch (RefSeq drift) -- which only a content-level check like this detects.
 *
 * @author Brian Bushnell, Amber
 * @date August 16, 2026
 */
public class ValidatePairs implements Accumulator<ValidatePairs.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Program entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		ValidatePairs x=new ValidatePairs(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor; parses arguments and collects input directories.
	 * @param args Command line arguments
	 */
	public ValidatePairs(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("dir") || a.equals("path")){
				if(b!=null){for(String s : b.split(",")){dirs.add(s);}}
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("printpass") || a.equals("showpass")){
				printPass=Parse.parseBoolean(b);
			}else if(a.equals("maxbad") || a.equals("maxreport")){
				maxReport=Integer.parseInt(b);
			}else if(a.equals("checkfnaingff") || a.equals("fnaingff") || a.equals("bidirectional")){
				checkFnaInGff=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//Standard flags in the parser (t=, ow=, etc.)
			}else if(b==null && new File(arg).isDirectory()){
				dirs.add(arg);//Bare directory as a positional argument
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(out==null && parser.out1!=null){out=parser.out1;}
		overwrite=parser.overwrite;
		append=parser.append;

		if(dirs.isEmpty()){
			throw new RuntimeException("At least one input directory is required (in=<dir>).");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Scans the directories, validates every pair, and writes the report.
	 * @param t Timer */
	void process(Timer t){

		ArrayList<Pair> pairList=scanPairs();
		pairs=pairList.toArray(new Pair[pairList.size()]);
		outstream.println("Found "+pairs.length+" basename(s) across "+dirs.size()+
				" director"+(dirs.size()==1 ? "y" : "ies")+".");

		spawnThreads();

		writeReport();

		t.stop();
		outstream.println();
		outstream.println("Pairs checked:          \t"+pairsChecked);
		outstream.println("Passed:                 \t"+passCount);
		outstream.println("Failed:                 \t"+failCount);
		outstream.println("  unpaired:             \t"+unpairedCount);
		outstream.println("  empty file:           \t"+emptyCount);
		outstream.println("  partial GFF records:  \t"+partialCount);
		outstream.println("  gff-seqid-not-in-fna: \t"+gffNotInFnaCount);
		outstream.println("  fna-seqid-not-in-gff: \t"+fnaNotInGffCount);
		outstream.println();
		outstream.println("Time:                   \t"+t);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	/** Scans input directories and pairs .fna files to .gff files by basename.
	 * @return List of pairs (either side may be null if unpaired) */
	private ArrayList<Pair> scanPairs(){
		HashMap<String, String> fnaMap=new HashMap<String, String>();
		HashMap<String, String> gffMap=new HashMap<String, String>();
		for(String dir : dirs){
			File d=new File(dir);
			File[] files=d.listFiles();
			if(files==null){
				outstream.println("Warning: not a readable directory: "+dir);
				continue;
			}
			for(File f : files){
				if(!f.isFile()){continue;}
				String name=f.getName();
				String base=strip(name, FNA_EXT);
				if(base!=null){fnaMap.put(base, f.getAbsolutePath()); continue;}
				base=strip(name, GFF_EXT);
				if(base!=null){gffMap.put(base, f.getAbsolutePath());}
			}
		}
		HashSet<String> allBases=new HashSet<String>(fnaMap.keySet());
		allBases.addAll(gffMap.keySet());
		ArrayList<Pair> list=new ArrayList<Pair>(allBases.size());
		for(String base : allBases){
			list.add(new Pair(base, fnaMap.get(base), gffMap.get(base)));
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawns worker threads that pull pairs from a shared atomic index. */
	private void spawnThreads(){
		final int threads=Tools.max(1, Tools.min(Shared.threads(), pairs.length));
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){alpt.add(new ProcessThread(i));}
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
	}

	/** Merges one finished thread's results into the shared list.
	 * @param pt Completed worker thread */
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			results.addAll(pt.resultsT);
			errorState|=!pt.success;
			errorState|=pt.errorStateT;
		}
	}

	/** @return true if no error occurred */
	@Override
	public final boolean success(){return !errorState;}

	/** @return Lock for the accumulator */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------          Reporting           ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes the per-pair failure report (to out= or the output stream) and tallies counts. */
	private void writeReport(){
		Result[] arr=results.toArray(new Result[results.size()]);
		Arrays.sort(arr, (x, y) -> x.name.compareTo(y.name));//Deterministic order

		ByteStreamWriter bsw=(out==null ? null : makeBSW());
		ByteBuilder bb=new ByteBuilder();
		for(Result r : arr){
			pairsChecked++;
			if(r.pass){
				passCount++;
				if(printPass){
					bb.append("PASS\t").append(r.name).append('\t').append(r.fnaSeqs).append(" seqs, ").append(r.features).append(" features").nl();
				}
			}else{
				failCount++;
				if(r.unpaired){unpairedCount++;}
				if(r.empty){emptyCount++;}
				if(r.partial){partialCount++;}
				if(r.gffNotInFna){gffNotInFnaCount++;}
				if(r.fnaNotInGff){fnaNotInGffCount++;}
				bb.append("FAIL\t").append(r.name).nl();
				for(String reason : r.reasons){
					bb.append('\t').append(reason).nl();
				}
			}
			if(bb.length()>16384){
				flush(bb, bsw);
			}
		}
		flush(bb, bsw);
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}

	/** Flushes the buffer to the writer or the output stream and clears it. */
	private void flush(ByteBuilder bb, ByteStreamWriter bsw){
		if(bb.length()<=0){return;}
		if(bsw!=null){bsw.print(bb.toBytes());}
		else{outstream.print(bb.toString());}
		bb.clear();
	}

	private ByteStreamWriter makeBSW(){
		FileFormat ff=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Validation (static)     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Validates a single fna+gff pair.
	 * @param p Pair to validate
	 * @param maxReport Max example seqids/lines to list per reason
	 * @param checkFnaInGff Whether to also report fna seqids absent from the gff
	 * @return Result describing pass/fail and reasons
	 */
	static Result validatePair(Pair p, int maxReport, boolean checkFnaInGff){
		Result r=new Result(p.name);
		if(p.fna==null){r.unpaired=true; r.fail("UNPAIRED: gff present but no matching fna"); return r;}
		if(p.gff==null){r.unpaired=true; r.fail("UNPAIRED: fna present but no matching gff"); return r;}

		//--- Read FNA sequence headers ---
		final HashSet<String> fnaSeqids=new HashSet<String>();
		long fnaSeqs=0;
		{
			FileFormat ff=FileFormat.testInput(p.fna, FileFormat.FASTA, null, false, false);
			ByteFile bf=ByteFile.makeByteFile(ff);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0 && line[0]=='>'){
					fnaSeqs++;
					fnaSeqids.add(fnaSeqid(line));
				}
			}
			if(bf.close()){r.fail("FNA_READ_ERROR");}
		}
		if(fnaSeqs==0){r.empty=true; r.fail("EMPTY_FNA: no sequences");}

		//--- Read GFF ---
		final HashSet<String> gffSeqids=new HashSet<String>();//union of feature seqids + sequence-region pragmas
		long features=0, malformed=0, lineNum=0;
		boolean inFasta=false;
		final ArrayList<String> malformedEx=new ArrayList<String>();
		{
			FileFormat ff=FileFormat.testInput(p.gff, FileFormat.GFF, null, false, false);
			ByteFile bf=ByteFile.makeByteFile(ff);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNum++;
				if(line.length==0 || inFasta){continue;}
				if(line[0]=='#'){
					if(matchesPrefix(line, FASTA_DIRECTIVE)){inFasta=true;}
					else if(matchesPrefix(line, SEQREGION_DIRECTIVE)){
						String sr=sequenceRegionSeqid(line);
						if(sr!=null){gffSeqids.add(sr);}
					}
					continue;
				}
				//Feature line
				features++;
				final int fields=countTabs(line)+1;
				if(fields!=9){
					malformed++;
					if(malformedEx.size()<maxReport){malformedEx.add("line "+lineNum+"="+fields+"f");}
				}
				gffSeqids.add(gffSeqid(line));
			}
			if(bf.close()){r.fail("GFF_READ_ERROR");}
		}
		if(features==0){r.empty=true; r.fail("EMPTY_GFF: no feature lines");}
		if(malformed>0){r.partial=true; r.fail("PARTIAL_RECORDS ("+malformed+"): "+joinExamples(malformedEx, maxReport));}

		//--- Cross-check headers (both directions) ---
		final ArrayList<String> gffNotFna=notIn(gffSeqids, fnaSeqids);
		if(!gffNotFna.isEmpty()){
			r.gffNotInFna=true;
			r.fail("GFF_SEQID_NOT_IN_FNA ("+gffNotFna.size()+"): "+joinExamples(gffNotFna, maxReport));
		}
		if(checkFnaInGff){
			final ArrayList<String> fnaNotGff=notIn(fnaSeqids, gffSeqids);
			if(!fnaNotGff.isEmpty()){
				r.fnaNotInGff=true;
				r.fail("FNA_SEQID_NOT_IN_GFF ("+fnaNotGff.size()+"): "+joinExamples(fnaNotGff, maxReport));
			}
		}

		r.fnaSeqs=fnaSeqs;
		r.features=features;
		return r;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Static Byte Helpers     ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts tab bytes in a line. */
	static int countTabs(byte[] line){
		int c=0;
		for(int i=0; i<line.length; i++){if(line[i]=='\t'){c++;}}
		return c;
	}

	/** Extracts the FNA sequence id: the token after '>' up to the first whitespace. */
	static String fnaSeqid(byte[] line){
		int a=1, b=1;
		while(b<line.length && line[b]!=' ' && line[b]!='\t'){b++;}
		return new String(line, a, b-a, StandardCharsets.US_ASCII);
	}

	/** Extracts the GFF seqid: field 0, up to the first tab. */
	static String gffSeqid(byte[] line){
		int b=0;
		while(b<line.length && line[b]!='\t'){b++;}
		return new String(line, 0, b, StandardCharsets.US_ASCII);
	}

	/** Returns true if the line starts with the given ASCII prefix. */
	static boolean matchesPrefix(byte[] line, String pre){
		if(line.length<pre.length()){return false;}
		for(int i=0; i<pre.length(); i++){
			if(line[i]!=pre.charAt(i)){return false;}
		}
		return true;
	}

	/** Parses the seqid from a "##sequence-region <seqid> <start> <end>" pragma line. */
	static String sequenceRegionSeqid(byte[] line){
		int i=SEQREGION_DIRECTIVE.length();
		while(i<line.length && (line[i]==' ' || line[i]=='\t')){i++;}
		int a=i;
		while(i<line.length && line[i]!=' ' && line[i]!='\t'){i++;}
		if(i<=a){return null;}
		return new String(line, a, i-a, StandardCharsets.US_ASCII);
	}

	/** Returns members of a that are not in b. */
	static ArrayList<String> notIn(HashSet<String> a, HashSet<String> b){
		ArrayList<String> list=new ArrayList<String>();
		for(String s : a){if(!b.contains(s)){list.add(s);}}
		return list;
	}

	/** Joins up to max examples, appending an overflow count. */
	static String joinExamples(ArrayList<String> list, int max){
		StringBuilder sb=new StringBuilder();
		final int n=Math.min(max, list.size());
		for(int i=0; i<n; i++){
			if(i>0){sb.append(", ");}
			sb.append(list.get(i));
		}
		if(list.size()>n){sb.append(", ...(+").append(list.size()-n).append(" more)");}
		return sb.toString();
	}

	/** Strips the first matching extension (case-insensitive); returns the basename or null. */
	static String strip(String name, String[] exts){
		for(String e : exts){
			if(name.length()>e.length() && name.regionMatches(true, name.length()-e.length(), e, 0, e.length())){
				return name.substring(0, name.length()-e.length());
			}
		}
		return null;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** A basename with its (possibly-null) fna and gff paths. */
	static class Pair{
		Pair(String name_, String fna_, String gff_){name=name_; fna=fna_; gff=gff_;}
		final String name, fna, gff;
	}

	/** Per-pair validation result. */
	static class Result{
		Result(String name_){name=name_;}
		void fail(String reason){pass=false; reasons.add(reason);}
		final String name;
		boolean pass=true;
		boolean unpaired=false, empty=false, partial=false, gffNotInFna=false, fnaNotInGff=false;
		long fnaSeqs=0, features=0;
		final ArrayList<String> reasons=new ArrayList<String>();
	}

	/** Worker thread; pulls pairs from the shared atomic index and validates each. */
	class ProcessThread extends Thread {

		ProcessThread(final int tid_){tid=tid_;}

		@Override
		public void run(){
			final int n=pairs.length;
			for(int i=nextIndex.getAndIncrement(); i<n; i=nextIndex.getAndIncrement()){
				Result r;
				try{
					r=validatePair(pairs[i], maxReport, checkFnaInGff);
				}catch(Throwable e){
					r=new Result(pairs[i].name);
					r.fail("EXCEPTION: "+e);
					errorStateT=true;
				}
				resultsT.add(r);
			}
			success=true;
		}

		final int tid;
		boolean success=false;
		boolean errorStateT=false;
		final ArrayList<Result> resultsT=new ArrayList<Result>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<String> dirs=new ArrayList<String>();
	private String out=null;
	private boolean printPass=false;
	private int maxReport=10;
	private boolean checkFnaInGff=true;

	private Pair[] pairs;
	private final AtomicInteger nextIndex=new AtomicInteger(0);
	private final ArrayList<Result> results=new ArrayList<Result>();

	private long pairsChecked=0, passCount=0, failCount=0;
	private long unpairedCount=0, emptyCount=0, partialCount=0, gffNotInFnaCount=0, fnaNotInGffCount=0;

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

	private static final String[] FNA_EXT={".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa"};
	private static final String[] GFF_EXT={".gff3.gz", ".gff.gz", ".gtf.gz", ".gff3", ".gff", ".gtf"};
	private static final String FASTA_DIRECTIVE="##FASTA";
	private static final String SEQREGION_DIRECTIVE="##sequence-region";

}
