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
import shared.Parse;
import shared.Parser;
import shared.PreParser;
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
 * Counts patterns in Accessions.
 * Handles hashing for Accession to TaxID lookups.
 * @author Brian Bushnell
 * @date May 9, 2018
 *
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
			bsw.start();
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
					//This one is not really important; the header could be missing.
					assert(Tools.startsWith(lines.get(0), "accession")) : bf.name()+"[0]: "+new String(lines.get(0));
				}else{
					assert(!Tools.startsWith(lines.get(0), "accession")) : bf.name()+"["+lines.id+"]: "+new String(lines.get(0));
				}
				for(byte[] line : lines){
					if(line.length>0){
						linesProcessedT++;
						bytesProcessedT+=(line.length+1);
						
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
		HashMap<String, Integer> map=new HashMap<String, Integer>(list.size()*3);
		codeBits=(int)Math.ceil(Tools.log2(list.size()));
		final int patternBits=63-codeBits;
		final long maxCombos=((1L<<(patternBits-1))-1);
		for(int i=0; i<list.size(); i++){
			String s=list.get(i);
			longestPattern=Tools.max(longestPattern, s.length());
			long combos=combos(s);
			if(combos<0 || combos>=maxCombos){map.put(s, -1);}
			else{map.put(s, i);}
		}
		codeMap=map;
		return map;
	}
	
	public static long digitize(String s){
		String pattern=remap(s);
		Integer code=codeMap.get(pattern);
		if(code==null){return -2;}
		if(code.intValue()<0){return -1;}
		
		long number=0;
		for(int i=0; i<pattern.length(); i++){
			char c=s.charAt(i);
			char p=pattern.charAt(i);
			if(p=='-' || p=='?'){
				//do nothing
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : s;
			}
		}
		number=(number<<codeBits)+code;
		return number;
	}
	
	public static long digitize(byte[] s){
		String pattern=remap(s);
		Integer code=codeMap.get(pattern);
		if(code==null){return -2;}
		if(code.intValue()<0){return -1;}
		
		long number=0;
		for(int i=0; i<pattern.length(); i++){
			byte c=s[i];
			char p=pattern.charAt(i);
			if(p=='-' || p=='?'){
				//do nothing
			}else if(p=='D'){
				number=(number*10)+(c-'0');
			}else if(p=='L'){
				number=(number*26)+(Tools.toUpperCase(c)-'A');
			}else{
				assert(false) : new String(s);
			}
		}
		number=(number<<codeBits)+code;
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
		return buffer.toString();
	}
	
	public static String remap(byte[] s){
		ByteBuilder buffer=new ByteBuilder(s.length);
		for(int i=0; i<s.length; i++){
			final byte b=s[i];
			if(b==' ' || b=='\t' || b=='.' || b==':'){break;}
			buffer.append((char)remap[b]);
		}
		return buffer.toString();
	}
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	private String out=null;
	private boolean perFile=true;
	
	/*--------------------------------------------------------------*/

	private HashMap<String, StringNum> countMap=new HashMap<String, StringNum>();
	public static HashMap<String, Integer> codeMap;
	private static int codeBits=-1;
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
