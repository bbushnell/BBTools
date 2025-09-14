package ml;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.EntropyTracker;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2023
 *
 */
public class SequenceToVector {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SequenceToVector x=new SequenceToVector(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public SequenceToVector(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		int k_=0, d_=Integer.MAX_VALUE;
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("result")){
				result0=Float.parseFloat(b);
			}else if(a.equals("width")){
				width=Integer.parseInt(b);
			}else if(a.equals("dimensions") || a.equals("dims")){
				d_=Integer.parseInt(b);
				if(d_<1) {d_=Integer.MAX_VALUE;}
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("parse") || a.equals("parseheader")){
				parseHeader=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
				mode=SPECTRUM;//(k_>0 ? SPECTRUM : RAW);
			}else if(a.equals("spectrum") || a.equals("spectra") || a.equals("kfreq") || a.equals("kmerfrequency")){
				boolean c=Parse.parseBoolean(b);
				mode=(c ? SPECTRUM : RAW);
			}else if(a.equals("raw")){
				boolean c=Parse.parseBoolean(b);
				mode=(!c ? SPECTRUM : RAW);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		

		setDimensions(d_);
		if(mode==SPECTRUM) {
			assert(k_>0) : "k must be in range 1-8 in spectrum mode: "+k_;
			k=k_;
			fullSpace=(1<<(2*k));
			kSpace=calcKSpace(k);
			kmap=kmap(k, maxDimensions);
		}else {
			k=fullSpace=kSpace=0;
			kmap=null;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		assert(result0!=-1 || parseHeader);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		
		long readsProcessed=0, basesProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			final int bufWidth=(k<1 ? 0 : 4+Tools.min(kspaceArray[k], maxDimensions));
			int inputs=(k<1 ? width*4+4 : bufWidth);
			bsw.print("#dims\t"+inputs+"\t1\n");
			final ByteBuilder bb=new ByteBuilder();
			final float[] buffer=(k<1 ? null : new float[bufWidth]);
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

					float result=(parseHeader ? Parse.parseFloat(r1.id, "result=", '\t') : result0);
					toVector(r1, bb, bsw, width, k, result, buffer, rcomp);
					toVector(r2, bb, bsw, width, k, result, buffer, rcomp);
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}

		errorState=bsw.poisonAndWait() | errorState;
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		assert(!errorState) : "An error was encountered.";
	}
	
	/*--------------------------------------------------------------*/
	
	private static void toVector(Read r, ByteBuilder bb, ByteStreamWriter bsw, int width, int k,
			float result, float[] buffer, boolean rcomp) {
		if(r==null) {return;}
		toVector(r.bases, bb, width, k, result, buffer);
		if(bsw!=null) {bsw.println(bb);}
		bb.clear();
		if(!rcomp) {return;}
		toVector(r.reverseComplement(), bb, bsw, width, k, result, buffer, false);
	}
	
	private static ByteBuilder toVector(byte[] bases, ByteBuilder bb, int width, int k, float result, float[] buffer) {
		float len=bases.length/(float)(width+5);
		float gc=Tools.calcGC(bases);
		EntropyTracker[] eTrackers=localETrackers.get();
		final int tnum=Tools.min(bases.length, eTrackers.length-1);
		EntropyTracker eTracker=eTrackers[tnum];
		assert(eTracker!=null) : tnum+", "+eTrackers.length;
		float entropy=eTracker.averageEntropy(bases, true);
		float poly=Read.longestHomopolymer(bases);
		poly=poly/(poly+5);
		
		if(k<1) {
			bb.append(len, 4).tab().append(gc, 4).tab().append(entropy, 4).tab().append(poly, 4);
			appendRaw(bases, bb, width);
		}else {
			Arrays.fill(buffer, 0);
//			int lim=Tools.min(width, bases.length);
			int count=fillSpectrum(bases, buffer, 4, 0, bases.length, k);
			len=(count*0.25f)/kspaceArray[k];
			buffer[0]=len; buffer[1]=gc; buffer[2]=entropy; buffer[3]=poly;
			appendSpectrum(bb, buffer);
		}
		
		if(result==(int)result) {
			bb.tab().append((int)result);
		}else {
			bb.tab().append(result, 4);
		}
		return bb;
	}
	
	private static void appendRaw(byte[] bases, ByteBuilder bb, int width) {
		for(int i=0; i<bases.length && i<width; i++) {
			byte b=bases[i];
			int n=AminoAcid.baseToNumber4[b];
			bb.tab().append(hotCodes[n]);
		}
		for(int i=bases.length; i<width; i++) {//Pad with zeros
			bb.tab().append("0\t0\t0\t0");
		}
	}
	
	private static void appendSpectrum(ByteBuilder bb, float[] vec) {
		for(float f : vec) {
			bb.append(f, 5, true).tab();
		}
		bb.length--;
	}
	
	/*--------------------------------------------------------------*/

	public static float[] fillVector(byte[] bases, float[] vec, int k) {
		assert(vec!=null);
		return fillVector(bases, vec, k, 0, bases.length-1);
	}
	
	public static float[] fillVector(byte[] bases, float[] vec, int k, int from, int to) {//To is inclusive
		assert(vec!=null);
		final int len=(to-from+1);
		final int width=(vec.length-4)/4;
		float flen=len/(float)(width+5);
		float gc=Tools.calcGC(bases, from, to);
		EntropyTracker[] eTrackers=localETrackers.get();
		EntropyTracker eTracker=eTrackers[Tools.min(len, eTrackers.length-1)];
		float entropy=eTracker.averageEntropy(bases, true, from, to);
		float poly=Read.longestHomopolymer(bases, from, to);
		poly=poly/(poly+5);
		Arrays.fill(vec, 0);
		vec[0]=flen;
		vec[1]=gc;
		vec[2]=entropy;
		vec[3]=poly;
//		System.err.println(Arrays.toString(vec));
		
		int lim=Tools.min(to+1, from+width, bases.length);
		if(k<1) {
			fillRaw(bases, vec, 4, from, lim);//To is exclusive
			assert(vec.length==width*4+4);
		}else {
			fillSpectrum(bases, vec, 4, from, lim, k);
		}
		
//		vec[width*4+4]=result;
		return vec;
	}
	
	public static void fillRaw(byte[] bases, float[] vec, int offset, int from, int to) {
		for(int i=from, j=offset; i<to; i++, j+=4) {
			byte b=bases[i];
			int n=AminoAcid.baseToNumber0[b];
			vec[j+n]=1;
		}
	}
	
	public static int fillSpectrum(byte[] bases, float[] vec, int offset, int from, int to, int k) {
		final int[] map=kmapArray[k];

//		assert(false) : maxDimensions+", "+Arrays.toString(map);
		final int kspace=Tools.min(kspaceArray[k], maxDimensions);
		int count=0, len=0;
		int kmer=0;
		int mask=~((-1)<<(2*k));
		for(int i=from; i<to; i++) {
			byte b=bases[i];
			int x=AminoAcid.baseToNumber0[b];
			kmer=((kmer<<2)|x)&mask;
			len=(AminoAcid.isFullyDefined(b) ? len+1 : 0);
			if(len>=k){
				vec[map[kmer]+offset]++;
				count++;
			}
		}
		float mult=(kspace*0.25f)/count;//Relative fractions with an average of 0.25.
		for(int i=offset; i<vec.length; i++) {vec[i]*=mult;}
		return count;
	}
	
	public static float score(byte[] bases, float[] vec, CellNet net, boolean rcomp) {
		net.applyInput(vec);
		final float r, f=net.feedForward();
		if(rcomp) {
			AminoAcid.reverseComplementBasesInPlace(bases);
			net.applyInput(vec);
			r=net.feedForward();
			AminoAcid.reverseComplementBasesInPlace(bases);
		}else {r=f;}
		return Tools.max(r, f);
	}
	
	/*--------------------------------------------------------------*/
	
	public static int calcKSpace(int k) {
		assert(k<16 && k>=0) : k;
		final int fullSpace=1<<(2*k);
//		if((k&1)==1) {return fullSpace/2;}
		int count=0;
		for(int kmer=0; kmer<fullSpace; kmer++) {//This is slow; a closed-form solution could be used.
			int rcomp=AminoAcid.reverseComplementBinaryFast(kmer, k);
			count+=(rcomp>=kmer ? 1 : 0);
		}
		int palindromes=(((k&1)==1) ? 0 : 1<<k);
		int closedForm=(fullSpace+palindromes)/2;
		assert(count==closedForm) : count+", "+closedForm;
		return closedForm;
	}
	
	/** Produces a condensed kmer space */
	public static int[] kmap(int k, int maxDims) {
		assert(k<16 && k>=0) : k;
		final int fullSpace=1<<(2*k);
		final int[] map=new int[fullSpace];
		
		int count=0;
		final Random randy=(maxDims<fullSpace ? new Random(k) : null);
		for(int kmer=0; kmer<fullSpace; kmer++) {
			int rcomp=AminoAcid.reverseComplementBinaryFast(kmer, k);
			if(kmer<=rcomp) {
				if(count<maxDims) {
					map[kmer]=count;
				}else {
					map[kmer]=randy.nextInt(maxDims);
				}
				count++;
			}else {
				map[kmer]=map[rcomp];
			}
		}
		return map;
	}
	
	public static void setDimensions(int maxDims) {
		if(maxDimensions==maxDims) {return;}
		maxDimensions=maxDims;
		fillArrays(maxDimensions);
	}
	
	private static void fillArrays(int maxDims) {
		for(int k=1; k<=kMax; k++) {
			fullspaceArray[k]=(1<<(2*k));
			kspaceArray[k]=calcKSpace(k);
			kmapArray[k]=kmap(k, maxDims);
		}
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;

	private static final ThreadLocal<EntropyTracker[]> localETrackers=new ThreadLocal<EntropyTracker[]>(){
        @Override protected EntropyTracker[] initialValue() {
        	EntropyTracker[] array=new EntropyTracker[maxWindow+1];
        	for(int i=minWindow; i<array.length; i++) {
        		array[i]=new EntropyTracker(ke, i, false);
        	}
        	return array;
        }
    };
//	public EntropyTracker eTracker=new EntropyTracker(3, 20, false);
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private boolean errorState=false;
	
	/*--------------------------------------------------------------*/

	private boolean rcomp=false;
	private boolean parseHeader=false;
	private int width=55;
	float result0=-1;
	
	/*--------------------------------------------------------------*/
	
	final int k;//=4;
	final int fullSpace;//=(1<<(2*k));
	final int kSpace;//=calcKSpace(k);
	final int[] kmap;

	private static final int[] fullspaceArray;
	private static final int[] kspaceArray;
	private static final int[][] kmapArray;
	
	private static int maxDimensions=Integer.MAX_VALUE;
	private static final int kMax=8;
	
	static {
		fullspaceArray=new int[kMax+1];
		kspaceArray=new int[kMax+1];
		kmapArray=new int[kMax+1][];
		fillArrays(maxDimensions);
	}
	
	/*--------------------------------------------------------------*/
	
	private int mode=RAW;
	static final int RAW=0, SPECTRUM=1;
	
	/*--------------------------------------------------------------*/

	public static final int minWindow=16;
	public static final int maxWindow=40;
	public static final int ke=3;
	
	/*--------------------------------------------------------------*/
	
	public static final String[] hotCodes=new String[] {"1\t0\t0\t0", "0\t1\t0\t0", "0\t0\t1\t0", "0\t0\t0\t1", "1\t0\t0\t0"};
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
