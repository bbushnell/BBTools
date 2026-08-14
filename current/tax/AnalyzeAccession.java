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
import template.ThreadWaiter;

/**
 * Counts patterns in Accessions and handles hashing for Accession to TaxID lookups.
 * Analyzes accession strings from input files to determine pattern distributions
 * and calculates the number of possible combinations for each pattern.
 * Supports both single-threaded and multi-threaded processing with per-file options.
 *
 * @author Brian Bushnell
 * @date May 9, 2018
 */
public class AnalyzeAccession implements Accumulator<AnalyzeAccession.ProcessThread> {

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

	void process(Timer t){

		if(perFile) {
			process_perFile();
		}else{
			for(FileFormat ffin : ffina){
				process_inner(ffin);
			}
		}

		if(ffout!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffout);
			bsw.start();//FIXED [tax/AnalyzeAccession#001]: start() MUST precede any print. print/println buffer into an ArrayBlockingQueue(5) whose addJob asserts(started) (ByteStreamWriter:292) and put()-blocks; writing before start() crashes under -ea (BBTools default) once buffered output reaches maxLen=32768 (i.e. >~1300 distinct L/D/- patterns). Latent today (real accession data has few distinct shapes -> <32KB) but a real contract violation. Matches the FindAncestor start()-first idiom. Twin: AnalyzeAccession_ST#001.
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

	void process_inner(FileFormat ffin){

		ByteFile bf=ByteFile.makeByteFile(ffin);

		final int threads=Tools.min(8, Shared.threads());
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){alpt.add(new ProcessThread(bf));}
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
	}


	void process_perFile(){
		ArrayList<ArrayList<ProcessThread>> perFileList=new ArrayList<ArrayList<ProcessThread>>(ffina.length);
		for(FileFormat ffin : ffina) {
			ByteFile bf=ByteFile.makeByteFile(ffin);

			final int threads=Tools.min(16, Shared.threads());
			ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
			for(int i=0; i<threads; i++){alpt.add(new ProcessThread(bf));}
			perFileList.add(alpt);
			ThreadWaiter.startThreads(alpt);
		}
		for(ArrayList<ProcessThread> alpt : perFileList){
			boolean success=ThreadWaiter.waitForThreadsToFinish(alpt, this);
			errorState|=!success;
		}
	}

	/*--------------------------------------------------------------*/

	static class ProcessThread extends Thread {

		ProcessThread(ByteFile bf_){
			bf=bf_;
		}

		@Override
		public void run() {
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
		}

		void increment(byte[] line, StringBuilder buffer){
			buffer.setLength(0);
			for(int i=0; i<line.length; i++){
				final byte b=line[i];
				if(b==' ' || b=='\t' || b=='.' || b==':'){break;}
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

		private HashMap<String, StringNum> countMapT=new HashMap<String, StringNum>();
		private final ByteFile bf;
		long linesProcessedT=0;
		long linesOutT=0;
		long bytesProcessedT=0;

	}

	/*--------------------------------------------------------------*/

	@Override
	public void accumulate(ProcessThread t) {
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

	@Override
	public boolean success() {
		return !errorState;
	}

	/*--------------------------------------------------------------*/

	public static long combos(String s){
		double combos=1;
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='D'){combos*=10;}
			else if(c=='L'){combos*=26;}
		}
		return (combos>=Long.MAX_VALUE ? Long.MAX_VALUE : (long)Math.ceil(combos));
	}

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

	public static HashMap<String, Integer> loadCodeMap(String fname){
		assert(codeMap==null);
		TextFile tf=new TextFile(fname);
		ArrayList<String> list=new ArrayList<String>();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(!line.startsWith("#")){
				String[] split=line.split("\t");
				list.add(split[0]);
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
		patternShift=new int[numPatterns];
		patternLowBits=new int[numPatterns];

		int code0=0, code1=0;
		for(int i=0; i<numPatterns; i++){
			String s=list.get(i);
			if(isTable1[i]){
				if(combosArr[i]<0 || combosArr[i]>=maxCombos1){
					map.put(s, -1);
					patternShift[i]=-1;
				}else{
					int pc=code1++;
					patternShift[i]=codeBits1+1;
					patternLowBits[i]=(pc<<1)|1;
					map.put(s, i);
				}
			}else{
				int pc=code0++;
				patternShift[i]=codeBits0+1;
				patternLowBits[i]=(pc<<1)|0;
				map.put(s, i);
			}
		}

		codeMap=map;
		return map;
	}

	public static long digitize(String s){
		String pattern=remap(s);
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

	public static long digitize(byte[] s){
		String pattern=remap(s);
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

	public static String remap(String s){
		if(s==null || s.length()<1){return "";}
		ByteBuilder buffer=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char b=s.charAt(i);
			if(b==' ' || b=='\t' || b=='.' || b==':'){break;}
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

	public static String remap(byte[] s){
		ByteBuilder buffer=new ByteBuilder(s.length);
		for(int i=0; i<s.length; i++){
			final byte b=s[i];
			if(b==' ' || b=='\t' || b=='.' || b==':'){break;}
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

	/*--------------------------------------------------------------*/

	private ArrayList<String> in=new ArrayList<String>();
	private String out=null;
	private boolean perFile=true;

	/*--------------------------------------------------------------*/

	private HashMap<String, StringNum> countMap=new HashMap<String, StringNum>();
	public static HashMap<String, Integer> codeMap;
	static int codeBits0=-1;
	static int codeBits1=-1;
	static int[] patternShift;
	static int[] patternLowBits;
	private static int longestPattern=-1;

	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;

	/*--------------------------------------------------------------*/

	private final FileFormat[] ffina;
	private final FileFormat ffout;

	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	private static final byte[] remap=makeRemap();

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

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

}
