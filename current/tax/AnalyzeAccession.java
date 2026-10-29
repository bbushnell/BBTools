package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import map.LongHashMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringNum;
import template.Accumulator;

/**
 * Counts patterns in Accessions and handles hashing for Accession to TaxID lookups.
 * Analyzes accession strings from input files to determine pattern distributions
 * and calculates the number of possible combinations for each pattern.
 * Also owns the static pattern tables used by disk-table builders, TaxServer,
 * and AccessionToTaxid to turn accession text into compact numeric keys.
 * Supports both single-threaded and multi-threaded processing with per-file options.
 *
 * @author Brian Bushnell
 * @date May 9, 2018
 */
public class AnalyzeAccession implements Accumulator<AnalyzeAccession.ProcessThread>{

	/**
	 * Runs accession-pattern analysis from the command line.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		AnalyzeAccession x=new AnalyzeAccession(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses command-line arguments and prepares input and output formats.
	 * @param args Command-line arguments
	 */
	public AnalyzeAccession(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				if(b==null){in.clear();}
				else{
					String[] split2=b.split(",");
					for(String s2 : split2){
						in.add(s2);
					}
				}
			}else if(a.equals("perfile")){
				perFile=Parse.parseBoolean(b);
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}

		assert(FastaReadInputStream.settingsOK());

		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}

		//		if(!ByteFile.FORCE_MODE_BF2){
		//			ByteFile.FORCE_MODE_BF2=false;
		//			ByteFile.FORCE_MODE_BF1=true;
		//		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}

		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ffina=new FileFormat[in.size()];
		for(int i=0; i<in.size(); i++){
			ffina[i]=FileFormat.testInput(in.get(i), FileFormat.TXT, null, true, false);
		}
	}

	/**
	 * Counts accession patterns from all inputs, optionally writes a sorted summary,
	 * and throws after worker, reader or output-writer failures.
	 * @param t Runtime timer started by the caller
	 */
	void process(Timer t){

		if(perFile){
			process_perFile();
		}else{
			for(FileFormat ffin : ffina){
				process_inner(ffin);
			}
		}

		if(ffout!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffout);
			/* FIXED [tax/AnalyzeAccession#001]: start() must precede any
			 * print; addJob requires a started writer once enough pattern
			 * lines fill its internal buffer. */
			bsw.start();
			bsw.println("#Pattern\tCount\tCombos\tBits");
			ArrayList<StringNum> list=new ArrayList<StringNum>();
			list.addAll(countMap.values());
			Collections.sort(list);
			Collections.reverse(list);
			for(StringNum sn : list){
				double combos=1;
				for(int i=0; i<sn.s.length(); i++){
					char c=sn.s.charAt(i);
					if(c=='D'){combos*=10;}
					else if(c=='L'){combos*=26;}
				}
				bsw.print(sn.toString().getBytes());
				bsw.println("\t"+(long)combos+"\t"+Tools.format("%.2f", Tools.log2(combos)));
			}
			errorState|=bsw.poisonAndWait();
		}

		t.stop();

		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));

		outstream.println();
		outstream.println("Valid Lines:       \t"+linesOut);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Counts one input file with a shared reader and up to eight worker threads.
	 * Worker threads accumulate private pattern maps that merge after all workers
	 * terminate and the shared reader closes.
	 * @param ffin Input file format descriptor
	 */
	void process_inner(FileFormat ffin){

		final int threads=Tools.min(8, Shared.threads());
		final ArrayList<ProcessThreadGroup> groups=
				new ArrayList<ProcessThreadGroup>(1);
		Throwable failure=null;
		try{
			failure=startGroup(makeGroup(ffin, threads, groups));
		}catch(Throwable t){
			failure=recordFailure(failure, t);
		}
		failure=finishGroups(groups, failure);
		throwFailure(failure);
	}

	/**
	 * Counts all inputs concurrently, assigning one reader and a worker group to
	 * each file.  Each reader is shared only by the workers created for that file.
	 */
	void process_perFile(){
		final ArrayList<ProcessThreadGroup> groups=
				new ArrayList<ProcessThreadGroup>(ffina.length);
		Throwable failure=null;
		for(FileFormat ffin : ffina){
			try{
				final int threads=Tools.min(16, Shared.threads());
				failure=startGroup(makeGroup(ffin, threads, groups));
			}catch(Throwable t){
				failure=recordFailure(failure, t);
			}
			if(failure!=null){break;}
		}
		failure=finishGroups(groups, failure);
		throwFailure(failure);
	}

	/**
	 * Registers a new ownership group, opens its reader and attaches workers.
	 * @param ffin Input file format descriptor
	 * @param threads Intended workers sharing this reader
	 * @param groups Ownership groups to close after worker termination
	 * @return Registered group with workers ready to start
	 */
	private static ProcessThreadGroup makeGroup(final FileFormat ffin,
			final int threads, final ArrayList<ProcessThreadGroup> groups){
		assert(threads>0) : "AnalyzeAccession requires at least one worker per "+
			"opened ByteFile; threads="+threads+" file="+ffin;
		final ProcessThreadGroup group=new ProcessThreadGroup(threads);
		groups.add(group);
		group.bf=ByteFile.makeByteFile(ffin);
		for(int i=0; i<threads; i++){group.workers.add(new ProcessThread(group.bf));}
		return group;
	}

	/**
	 * Starts workers one at a time and records only successfully-started threads.
	 * @param group Group whose workers should start
	 * @return First startup failure, or null if all starts succeeded
	 */
	private static Throwable startGroup(final ProcessThreadGroup group){
		for(ProcessThread pt : group.workers){
			try{
				pt.start();
				group.started.add(pt);
			}catch(Throwable t){
				return t;
			}
			assert(group.started.size()<=group.workers.size()) :
				"Started-worker tracking exceeded intended capacity; started="+
				group.started.size()+" workers="+group.workers.size();
		}
		return null;
	}

	/**
	 * Joins workers, closes owned readers, then merges successful worker output.
	 * @param groups Reader/worker groups to finish
	 * @param failure First failure already observed during open or start
	 * @return First failure across open, start, join, worker, close and merge
	 */
	private Throwable finishGroups(final ArrayList<ProcessThreadGroup> groups,
			Throwable failure){
		failure=joinGroups(groups, failure);
		failure=checkWorkers(groups, failure);
		failure=closeGroups(groups, failure);
		failure=accumulateGroups(groups, failure);
		errorState|=(failure!=null);
		return failure;
	}

	/**
	 * Joins every started worker without touching never-started threads.
	 * @param groups Groups whose started workers should terminate
	 * @param failure First prior failure, if any
	 * @return First prior or interrupted-join failure
	 */
	private static Throwable joinGroups(final ArrayList<ProcessThreadGroup> groups,
			Throwable failure){
		for(ProcessThreadGroup group : groups){
			for(ProcessThread pt : group.started){
				while(pt.getState()!=Thread.State.TERMINATED){
					try{
						pt.join();
					}catch(InterruptedException e){
						failure=recordFailure(failure, e);
					}
				}
			}
		}
		return failure;
	}

	/**
	 * Records abnormal worker exits before partial results can be merged.
	 * @param groups Joined worker groups to inspect
	 * @param failure First prior failure, if any
	 * @return First prior or worker-completion failure
	 */
	private static Throwable checkWorkers(final ArrayList<ProcessThreadGroup> groups,
			Throwable failure){
		for(ProcessThreadGroup group : groups){
			for(ProcessThread pt : group.started){
				if(pt.error!=null){
					failure=recordFailure(failure, pt.error);
				}else if(!pt.complete){
					failure=recordFailure(failure, new RuntimeException(
							"AnalyzeAccession worker failed to consume "+
							group.bf.name()+" to EOF without recording an error."));
				}
			}
		}
		return failure;
	}

	/**
	 * Closes every owned reader after all workers that share it have terminated.
	 * @param groups Groups whose readers should close
	 * @param failure First prior failure, if any
	 * @return First prior or close failure
	 */
	private static Throwable closeGroups(final ArrayList<ProcessThreadGroup> groups,
			Throwable failure){
		for(ProcessThreadGroup group : groups){
			if(group.bf==null){continue;}
			try{
				if(group.bf.close()){
					failure=recordFailure(failure, new RuntimeException(
							"Error reading or closing "+group.bf.name()));
				}
			}catch(Throwable t){
				failure=recordFailure(failure, t);
			}
		}
		return failure;
	}

	/**
	 * Merges only complete worker counts after all owned readers have closed.
	 * @param groups Joined and closed worker groups
	 * @param failure First prior failure, if any
	 * @return First prior or accumulation failure
	 */
	private Throwable accumulateGroups(final ArrayList<ProcessThreadGroup> groups,
			Throwable failure){
		for(ProcessThreadGroup group : groups){
			for(ProcessThread pt : group.started){
				try{
					accumulate(pt);
				}catch(Throwable t){
					failure=recordFailure(failure, t);
				}
			}
		}
		return failure;
	}

	/**
	 * Preserves the first failure and attaches later cleanup failures.
	 * @param primary First failure seen so far, or null
	 * @param secondary Later failure, or null
	 * @return Primary when present, otherwise secondary
	 */
	private static Throwable recordFailure(final Throwable primary,
			final Throwable secondary){
		if(secondary==null){return primary;}
		if(primary==null){return secondary;}
		if(primary!=secondary){primary.addSuppressed(secondary);}
		return primary;
	}

	/**
	 * Rethrows a recorded failure without wrapping unchecked failures.
	 * @param failure Failure to throw, or null
	 */
	private static void throwFailure(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		throw new RuntimeException(failure);
	}

	/*--------------------------------------------------------------*/

	/** Reader ownership group for one ByteFile and its workers. */
	private static class ProcessThreadGroup{

		/**
		 * Creates preallocated worker lists for one owned reader.
		 * @param threads Intended worker count
		 */
		ProcessThreadGroup(final int threads){
			workers=new ArrayList<ProcessThread>(threads);
			started=new ArrayList<ProcessThread>(threads);
		}

		/** Owned reader to close after all started workers terminate. */
		ByteFile bf;
		/** Workers intended for this reader. */
		final ArrayList<ProcessThread> workers;
		/** Successfully-started workers safe to join. */
		final ArrayList<ProcessThread> started;

	}

	/*--------------------------------------------------------------*/

	/** Worker that reads accession lines and counts remapped patterns privately. */
	static class ProcessThread extends Thread{

		/**
		 * Creates a worker over a shared ByteFile.
		 * @param bf_ Reader to consume
		 */
		ProcessThread(ByteFile bf_){
			bf=bf_;
		}

		/**
		 * Consumes accession lines from the shared reader and counts patterns.
		 * Records uncaught failures so the owner can skip partial counts.
		 */
		@Override
		public void run(){
			try{
				final StringBuilder buffer=new StringBuilder(128);
				for(ListNum<byte[]> lines=bf.nextList(); lines!=null; lines=bf.nextList()){
					assert(lines.size()>0);
					if(lines.id==0){
						//Header may start with "accession" (accession2taxid) or "#" (assembly_summary).
						assert(Tools.startsWith(lines.get(0), "accession") ||
							Tools.startsWith(lines.get(0), "#")) : bf.name()+"[0]: "+new String(lines.get(0));
					}else{
						assert(!Tools.startsWith(lines.get(0), "accession")) : bf.name()+"["+lines.id+"]: "+new String(lines.get(0));
					}
					for(byte[] line : lines){
						if(line.length>0){
							linesProcessedT++;
							bytesProcessedT+=(line.length+1);

							if(Tools.startsWith(line, "#")){continue;}//Skip comment/header lines (assembly_summary format)

							boolean valid=lines.id>0 || !(Tools.startsWith(line, "accession")); //Skips test for most lines

							if(valid){
								linesOutT++;
								increment(line, buffer);
							}
						}
					}
				}
				complete=true;
			}catch(Throwable t){
				error=t;
				throw t;
			}
		}

		/**
		 * Remaps one accession token into a pattern and increments this thread's
		 * private count for that pattern.
		 * @param line Input line beginning with the accession token
		 * @param buffer Reused scratch buffer
		 */
		void increment(byte[] line, StringBuilder buffer){
			buffer.setLength(0);
			for(int i=0; i<line.length; i++){
				final byte b=line[i];
				if(b==' ' || b=='\t' || b=='.' || b==':' || b==','){break;}
				final char b2=(char)remap[b];
				assert(b2!='?' || b=='+') : "unprocessed symbol in "+new String(line)+"\n"+"'"+(char)b+"'";
				buffer.append(b2);
			}
			if(buffer.length()>=3 && line.length>=3
					&& (line[0]=='N' || line[0]=='n') && (line[1]=='Z' || line[1]=='z')
					&& (line[2]=='_' || line[2]=='-')
					&& buffer.charAt(0)=='L' && buffer.charAt(1)=='L' && buffer.charAt(2)=='-'){
				buffer.setCharAt(0, 'N');
				buffer.setCharAt(1, 'Z');
			}
			String key=buffer.toString();
			StringNum value=countMapT.get(key);
			if(value!=null){value.increment();}
			else{countMapT.put(key, new StringNum(key, 1));}
		}

		/** Pattern counts produced by this worker. */
		private HashMap<String, StringNum> countMapT=new HashMap<String, StringNum>();
		/** Shared reader consumed by this worker. */
		private final ByteFile bf;
		/** Nonempty input lines seen by this worker. */
		long linesProcessedT=0;
		/** Non-header accession rows counted by this worker. */
		long linesOutT=0;
		/** Approximate input bytes seen by this worker, including line terminators. */
		long bytesProcessedT=0;
		/** True only after this worker consumes the shared reader to EOF. */
		boolean complete=false;
		/** Uncaught failure that killed this worker before normal EOF. */
		Throwable error=null;

	}

	/*--------------------------------------------------------------*/

	/**
	 * Merges one worker's private pattern counts and counters.
	 * @param t Finished worker to merge
	 */
	@Override
	public void accumulate(ProcessThread t){
		if(!t.complete || t.error!=null){
			errorState=true;
			return;
		}
		linesProcessed+=t.linesProcessedT;
		linesOut+=t.linesOutT;
		bytesProcessed+=t.bytesProcessedT;
		for(Entry<String, StringNum> e : t.countMapT.entrySet()){
			StringNum value=e.getValue();
			final String key=e.getKey();
			StringNum old=countMap.get(key);
			if(old==null){countMap.put(key, value);}
			else{old.add(value);}
		}
	}

	/**
	 * Reports whether a worker, reader or writer has raised the shared error state.
	 * @return {@code true} if processing completed without an error state
	 */
	@Override
	public boolean success(){
		return !errorState;
	}

	/*--------------------------------------------------------------*/

	/**
	 * Counts possible concrete accessions represented by a pattern.
	 * @param s Pattern using {@code D} for digits and {@code L} for letters
	 * @return Product of all variable positions, clamped to Long.MAX_VALUE
	 */
	public static long combos(String s){
		double combos=1;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='D'){combos*=10;}
			else if(c=='L'){combos*=26;}
		}
		return (combos>=Long.MAX_VALUE ? Long.MAX_VALUE : (long)Math.ceil(combos));
	}

	/**
	 * Counts possible concrete accessions represented by a byte pattern.
	 * @param s Pattern using {@code D} for digits and {@code L} for letters
	 * @return Product of all variable positions, or -1 if it overflows a long
	 */
	public static long combos(byte[] s){
		double combos=1;
		for(int i=0; i<s.length; i++){
			byte c=s[i];
			if(c=='D'){combos*=10;}
			else if(c=='L'){combos*=26;}
		}
		return (combos>=Long.MAX_VALUE ? -1 : (long)Math.ceil(combos));
	}

	/*--------------------------------------------------------------*/

	/**
	 * Loads tab-delimited accession patterns and prepares compact digitization
	 * tables.  The first column of each non-comment row is interpreted as a
	 * pattern; encodable patterns receive nonnegative dense indexes, while
	 * patterns whose combination space is too large are recorded as overflow.
	 * This method initializes the process-wide {@link #codeMap} and packed
	 * pattern table once; later disk lookup calls read those static tables.
	 * The pattern file is closed before global table construction or publication;
	 * close failures identify the file and are suppressed behind a primary read
	 * failure.
	 * @param fname Pattern table path
	 * @return Pattern to dense-index map, with -1 for overflow patterns
	 */
	public static HashMap<String, Integer> loadCodeMap(String fname){
		assert(codeMap==null);
		assert(patternCodeMap==null);
		TextFile tf=new TextFile(fname);
		ArrayList<String> list=new ArrayList<String>();
		Throwable failure=null;
		try{
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(!line.startsWith("#")){
					String[] split=line.split("\t");
					list.add(split[0]);
				}
			}
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			try{
				if(tf.close()){throw new RuntimeException("Error closing "+fname);}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}

		final int numPatterns=list.size();
		long[] combosArr=new long[numPatterns];
		for(int i=0; i<numPatterns; i++){
			String s=list.get(i);
			longestPattern=Tools.max(longestPattern, s.length());
			combosArr[i]=combos(s);
		}

		//Iteratively assign patterns to table 0; overflow goes to table 1
		boolean[] isTable1=new boolean[numPatterns];
		boolean changed=true;
		while(changed){
			changed=false;
			int count0=0;
			for(int i=0; i<numPatterns; i++){
				if(!isTable1[i]){count0++;}
			}
			if(count0==0){break;}
			int cb=Math.max(1, (int)Math.ceil(Tools.log2(count0)));
			long maxCombos0=(1L<<(62-cb));
			for(int i=0; i<numPatterns; i++){
				if(!isTable1[i] && (combosArr[i]<0 || combosArr[i]>=maxCombos0)){
					isTable1[i]=true;
					changed=true;
				}
			}
		}

		int count0=0, count1=0;
		for(int i=0; i<numPatterns; i++){
			if(isTable1[i]){count1++;}
			else{count0++;}
		}

		codeBits0=(count0>0) ? Math.max(1, (int)Math.ceil(Tools.log2(count0))) : 1;
		codeBits1=(count1>0) ? Math.max(1, (int)Math.ceil(Tools.log2(count1))) : 1;
		long maxCombos1=(1L<<(62-codeBits1));

		HashMap<String, Integer> map=new HashMap<String, Integer>(numPatterns*3);
		LongHashMap map2=new LongHashMap(numPatterns*3);
		patternShift=new int[numPatterns];
		patternLowBits=new int[numPatterns];
		patternList=list.toArray(new String[numPatterns]);

		int code0=0, code1=0;
		for(int i=0; i<numPatterns; i++){
			String s=list.get(i);
			if(isTable1[i]){
				if(combosArr[i]<0 || combosArr[i]>=maxCombos1){
					map.put(s, -1);
					patternShift[i]=-1;
					long key=patternKey(s);
					if(key>=0){map2.put(key, 0);}
				}else{
					int pc=code1++;
					patternShift[i]=codeBits1+1;
					patternLowBits[i]=(pc<<1)|1;
					map.put(s, i);
					long key=patternKey(s);
					if(key>=0){map2.put(key, i+1);}
				}
			}else{
				int pc=code0++;
				patternShift[i]=codeBits0+1;
				patternLowBits[i]=(pc<<1)|0;
				map.put(s, i);
				long key=patternKey(s);
				if(key>=0){map2.put(key, i+1);}
			}
		}

		codeMap=map;
		patternCodeMap=map2;
		return map;
	}

	/**
	 * Converts an ASCII accession string to the compact positive key used by
	 * accession lookup tables after {@link #loadCodeMap(String)} has initialized
	 * the static pattern tables.
	 * @param s Accession string; parsing stops at space, tab, '.', ':', or ','
	 * @return Positive key, -1 for overflow patterns, or -2 for no matching pattern
	 */
	public static long digitize(String s){
		final long key=remapKey(s);
		if(key<0){return digitize(remap(s), s);}
		final int code=patternCodeMap.get(key);
		if(code<0){return -2;}
		if(code==0){return -1;}
		final int idx=code-1;
		final String pattern=patternList[idx];

		long number=0;
		for(int i=0; i<pattern.length(); i++){
			char c=s.charAt(i);
			char p=pattern.charAt(i);
			if(p=='-' || p=='?' || p=='N' || p=='Z'){
				//do nothing - constant symbols, zero bits
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : s;
			}
		}
		number=((number+1)<<patternShift[idx])|patternLowBits[idx];
		return number;
	}

	/**
	 * Converts an ASCII accession byte array to the compact positive key used by
	 * accession lookup tables after {@link #loadCodeMap(String)} has initialized
	 * the static pattern tables.
	 * @param s Accession bytes; parsing stops at space, tab, '.', ':', or ','
	 * @return Positive key, -1 for overflow patterns, or -2 for no matching pattern
	 */
	public static long digitize(byte[] s){
		final long key=remapKey(s);
		if(key<0){return digitize(remap(s), s);}
		final int code=patternCodeMap.get(key);
		if(code<0){return -2;}
		if(code==0){return -1;}
		final int idx=code-1;
		final String pattern=patternList[idx];

		long number=0;
		for(int i=0; i<pattern.length(); i++){
			byte c=s[i];
			char p=pattern.charAt(i);
			if(p=='-' || p=='?' || p=='N' || p=='Z'){
				//do nothing - constant symbols, zero bits
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : new String(s);
			}
		}
		number=((number+1)<<patternShift[idx])|patternLowBits[idx];
		return number;
	}

	/**
	 * Slow-path digitization for accessions whose remapped pattern is too long for
	 * the packed pattern-code table.
	 * @param pattern Remapped accession pattern
	 * @param s Original accession string
	 * @return Positive key, -1 for overflow patterns, or -2 for no matching pattern
	 */
	private static long digitize(String pattern, String s){
		Integer idx=codeMap.get(pattern);
		if(idx==null){return -2;}
		if(idx.intValue()<0){return -1;}

		long number=0;
		for(int i=0; i<pattern.length(); i++){
			char c=s.charAt(i);
			char p=pattern.charAt(i);
			if(p=='-' || p=='?' || p=='N' || p=='Z'){
				//do nothing - constant symbols, zero bits
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : s;
			}
		}
		number=((number+1)<<patternShift[idx])|patternLowBits[idx];
		return number;
	}

	/**
	 * Slow-path digitization for byte accessions whose remapped pattern is too
	 * long for the packed pattern-code table.
	 * @param pattern Remapped accession pattern
	 * @param s Original accession bytes
	 * @return Positive key, -1 for overflow patterns, or -2 for no matching pattern
	 */
	private static long digitize(String pattern, byte[] s){
		Integer idx=codeMap.get(pattern);
		if(idx==null){return -2;}
		if(idx.intValue()<0){return -1;}

		long number=0;
		for(int i=0; i<pattern.length(); i++){
			byte c=s[i];
			char p=pattern.charAt(i);
			if(p=='-' || p=='?' || p=='N' || p=='Z'){
				//do nothing - constant symbols, zero bits
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : new String(s);
			}
		}
		number=((number+1)<<patternShift[idx])|patternLowBits[idx];
		return number;
	}

	/**
	 * Converts an accession prefix into its D/L/-/?/NZ pattern.
	 * @param s ASCII accession string; parsing stops at space, tab, '.', ':', or ','
	 * @return Remapped pattern, or an empty pattern for null and empty inputs
	 */
	public static String remap(String s){
		if(s==null || s.length()<1){return "";}
		ByteBuilder buffer=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char b=s.charAt(i);
			if(b==' ' || b=='\t' || b=='.' || b==':' || b==','){break;}
			buffer.append((char)remap[b]);
		}
		String result=buffer.toString();
		if(result.length()>=3 && s.length()>=3
				&& (s.charAt(0)=='N' || s.charAt(0)=='n')
				&& (s.charAt(1)=='Z' || s.charAt(1)=='z')
				&& (s.charAt(2)=='_' || s.charAt(2)=='-')
				&& result.startsWith("LL-")){
			result="NZ"+result.substring(2);
		}
		return result;
	}

	/**
	 * Converts accession bytes into their D/L/-/?/NZ pattern.
	 * @param s ASCII accession bytes; parsing stops at space, tab, '.', ':', or ','
	 * @return Remapped pattern
	 */
	public static String remap(byte[] s){
		ByteBuilder buffer=new ByteBuilder(s.length);
		for(int i=0; i<s.length; i++){
			final byte b=s[i];
			if(b==' ' || b=='\t' || b=='.' || b==':' || b==','){break;}
			buffer.append((char)remap[b]);
		}
		String result=buffer.toString();
		if(result.length()>=3 && s.length>=3
				&& (s[0]=='N' || s[0]=='n')
				&& (s[1]=='Z' || s[1]=='z')
				&& (s[2]=='_' || s[2]=='-')
				&& result.startsWith("LL-")){
			result="NZ"+result.substring(2);
		}
		return result;
	}

	/**
	 * Packs a remapped pattern into three-bit symbols for direct map lookup.
	 * @param s Remapped pattern
	 * @return Packed pattern key, or -1 if the pattern is too long to pack
	 */
	private static long patternKey(String s){
		if(s.length()>MAX_CODED_PATTERN_LENGTH){return -1;}
		long key=0;
		for(int i=0; i<s.length(); i++){
			key|=((long)patternCode(s.charAt(i)))<<(PATTERN_BITS*i);
		}
		return key;
	}

	/**
	 * Remaps and packs an accession string without allocating the intermediate
	 * pattern when it fits in a packed key.
	 * @param s ASCII accession string; parsing stops at space, tab, '.', ':', or ','
	 * @return Packed pattern key, or -1 if the pattern is too long to pack
	 */
	private static long remapKey(String s){
		if(s==null || s.length()<1){return 0;}
		long key=0;
		final boolean nz=(s.length()>=3 && (s.charAt(0)=='N' || s.charAt(0)=='n')
				&& (s.charAt(1)=='Z' || s.charAt(1)=='z') && (s.charAt(2)=='_' || s.charAt(2)=='-'));
		for(int i=0; i<s.length(); i++){
			final char b=s.charAt(i);
			if(b==' ' || b=='\t' || b=='.' || b==':' || b==','){break;}
			if(i>=MAX_CODED_PATTERN_LENGTH){return -1;}
			final char p=(nz && i<2) ? b : (char)remap[b];
			key|=((long)patternCode(p))<<(PATTERN_BITS*i);
		}
		return key;
	}

	/**
	 * Remaps and packs accession bytes without allocating the intermediate pattern
	 * when it fits in a packed key.
	 * @param s ASCII accession bytes; parsing stops at space, tab, '.', ':', or ','
	 * @return Packed pattern key, or -1 if the pattern is too long to pack
	 */
	private static long remapKey(byte[] s){
		long key=0;
		final boolean nz=(s.length>=3 && (s[0]=='N' || s[0]=='n')
				&& (s[1]=='Z' || s[1]=='z') && (s[2]=='_' || s[2]=='-'));
		for(int i=0; i<s.length; i++){
			final byte b=s[i];
			if(b==' ' || b=='\t' || b=='.' || b==':' || b==','){break;}
			if(i>=MAX_CODED_PATTERN_LENGTH){return -1;}
			final char p=(nz && i<2) ? (char)b : (char)remap[b];
			key|=((long)patternCode(p))<<(PATTERN_BITS*i);
		}
		return key;
	}

	/**
	 * Converts one remapped pattern character to its three-bit packed code.
	 * @param c Pattern character
	 * @return Packed code in the range 1..6, or 7 for an unrecognized symbol when
	 * assertions are disabled
	 * @throws AssertionError If assertions are enabled and the symbol is unrecognized
	 */
	private static int patternCode(char c){
		if(c=='D'){return 1;}
		if(c=='L' || c=='l'){return 2;}
		if(c=='-'){return 3;}
		if(c=='?'){return 4;}
		if(c=='N' || c=='n'){return 5;}
		if(c=='Z' || c=='z'){return 6;}
		assert(false) : c;
		return 7;
	}

	/*--------------------------------------------------------------*/

	/** Input filenames from {@code in=} or bare existing-file arguments. */
	private ArrayList<String> in=new ArrayList<String>();
	/** Optional summary output filename. */
	private String out=null;
	/** True to read all input files concurrently instead of serially. */
	private boolean perFile=true;

	/*--------------------------------------------------------------*/

	/** Aggregated pattern counts across all workers. */
	private HashMap<String, StringNum> countMap=new HashMap<String, StringNum>();
	/**
	 * Pattern-to-index table populated by {@link #loadCodeMap(String)} for
	 * slow-path patterns and for callers that test whether disk tables are loaded.
	 */
	public static HashMap<String, Integer> codeMap;
	/** Packed-pattern lookup table with stored indexes offset by one. */
	private static LongHashMap patternCodeMap;
	/** Dense-index to pattern lookup table. */
	private static String[] patternList;
	/** Bits assigned to table-0 pattern codes. */
	static int codeBits0=-1;
	/** Bits assigned to table-1 pattern codes. */
	static int codeBits1=-1;
	/** Per-pattern bit shift used to leave room for the pattern code. */
	static int[] patternShift;
	/** Per-pattern low bits containing the pattern code and table flag. */
	static int[] patternLowBits;
	/** Longest loaded pattern length. */
	private static int longestPattern=-1;
	/** Bits per pattern character in packed pattern keys. */
	private static final int PATTERN_BITS=3;
	/** Longest accession pattern that fits in a 63-bit packed pattern key. */
	private static final int MAX_CODED_PATTERN_LENGTH=21;

	/** Total nonempty input lines seen. */
	private long linesProcessed=0;
	/** Total non-header lines counted as accession rows. */
	private long linesOut=0;
	/** Approximate input bytes including line terminators. */
	private long bytesProcessed=0;
	/** Output bytes; currently unused. */
	private long bytesOut=0;

	/*--------------------------------------------------------------*/

	/** Input file descriptors. */
	private final FileFormat[] ffina;
	/** Optional summary output descriptor. */
	private final FileFormat ffout;

	/**
	 * Provides the lock used by the public Accumulator contract.
	 * @return Accumulator merge lock
	 */
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}
	/** Accumulator lock retained for the public Accumulator contract. */
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/** ASCII byte-to-pattern remap table. */
	private static final byte[] remap=makeRemap();

	/**
	 * Builds the ASCII byte-to-pattern lookup table used by the remap helpers.
	 * @return Remap table for byte values below 128
	 */
	private static byte[] makeRemap(){
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)'?');
		for(int i='A'; i<='Z'; i++){array[i]='L';}
		for(int i='a'; i<='z'; i++){array[i]='L';}
		for(int i='0'; i<='9'; i++){array[i]='D';}
		array['_']=array['-']='-';
		return array;
	}

	/*--------------------------------------------------------------*/

	/** Status and progress output stream. */
	private PrintStream outstream=System.err;
	/** Enables verbose logging in this class and the delegated I/O classes. */
	public static boolean verbose=false;
	/** Shared worker, reader and writer failure latch. */
	public boolean errorState=false;
	/** Allow output files to be overwritten. */
	private boolean overwrite=true;
	/** Append to existing output files instead of replacing them. */
	private boolean append=false;

}
