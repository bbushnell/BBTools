package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

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
import structures.StringNum;

/**
 * Analyzes patterns in biological sequence accession numbers for efficient hashing.
 * Counts accession pattern frequencies, calculates combinatorial possibilities, and supports digitization for TaxID lookup systems.
 * @author Brian Bushnell
 * @date May 9, 2018
 */
public class AnalyzeAccession_ST {
	
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AnalyzeAccession_ST x=new AnalyzeAccession_ST(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public AnalyzeAccession_ST(String[] args){
		
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

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
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
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

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

		for(FileFormat ffin : ffina){
			process_inner(ffin);
		}
		
		if(ffout!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffout);
			bsw.start();//FIXED [tax/AnalyzeAccession_ST#001] (twin of AnalyzeAccession#001): start() MUST precede any print. addJob asserts(started) (ByteStreamWriter:292) + ArrayBlockingQueue(5) put()-blocks; writing before start() crashes under -ea once buffered output reaches maxLen=32768. Latent (few distinct L/D/- shapes) but a real contract violation; matches the FindAncestor start()-first idiom.
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
		
		byte[] line=bf.nextLine();
		StringBuilder buffer=new StringBuilder(32);
		
		for(int lineNum=0; line!=null; lineNum++){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				assert((lineNum==0)==(Tools.startsWith(line, "accession"))) : "Line "+lineNum+": "+new String(line);
//				final boolean valid=(line[0]!='#');
				
				if(true){
					linesOut++;
					bytesOut+=(line.length+1);
					increment(line, buffer);
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
	}
	
	void increment(byte[] line, StringBuilder buffer){
		buffer.setLength(0);
		for(int i=0; i<line.length; i++){
			final byte b=line[i];
			if(b==' ' || b=='\t' || b=='.'){break;}
			buffer.append((char)remap[b]);
		}
		if(buffer.length()>=3 && line.length>=3
				&& (line[0]=='N' || line[0]=='n') && (line[1]=='Z' || line[1]=='z')
				&& (line[2]=='_' || line[2]=='-')
				&& buffer.charAt(0)=='L' && buffer.charAt(1)=='L' && buffer.charAt(2)=='-'){
			buffer.setCharAt(0, 'N');
			buffer.setCharAt(1, 'Z');
		}
		String key=buffer.toString();
		StringNum value=countMap.get(key);
		if(value!=null){value.increment();}
		else{countMap.put(key, new StringNum(key, 1));}
	}
	
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
		ByteBuilder buffer=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			final char b=s.charAt(i);
			if(b==' ' || b=='\t' || b=='.'){break;}
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
			if(b==' ' || b=='\t' || b=='.'){break;}
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
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat[] ffina;
	private final FileFormat ffout;
	
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
