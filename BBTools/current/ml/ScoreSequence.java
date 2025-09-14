package ml;

import java.util.ArrayList;

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
import structures.LongList;
import tracker.PalindromeTracker;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class ScoreSequence {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		ScoreSequence x=new ScoreSequence(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public ScoreSequence(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("net") || a.equals("nn")){
				netFile=b;
			}else if(a.equals("hist")){
				histFile=b;
			}else if(a.equals("width")){
				width=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("parse")){
				parseHeader=Parse.parseBoolean(b);
			}else if(a.equals("cutoff")){
				cutoff=Float.parseFloat(b);
				filter=true;
			}else if(a.equals("highpass")){
				highpass=Parse.parseBoolean(b);
				filter=true;
			}else if(a.equals("lowpass")){
				highpass=!Parse.parseBoolean(b);
				filter=true;
			}else if(a.equals("filter")){
				filter=Parse.parseBoolean(b);
			}else if(a.equals("annotate") || a.equals("rename")){
				annotate=Parse.parseBoolean(b);
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffnet=FileFormat.testOutput(netFile, FileFormat.BBNET, null, true, true, false, false);
		net=CellNetParser.load(netFile);
		assert(net!=null) : netFile;
		if(width<0) {width=(net.numInputs()-4)/4;}
		else {assert(width==(net.numInputs()-4)/4) : width+", "+net.numInputs();}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}

		int inputs=width*4+4;
		final ByteStreamWriter bsw=(ffout1==null ? null : new ByteStreamWriter(ffout1));
		if(bsw!=null) {
			bsw.start();
//			bsw.print("#dims\t"+inputs+"\t1\n");
		}
		
		final float[] vec=new float[width*4+4];
		long readsProcessed=0, basesProcessed=0;
		final ByteBuilder bb=new ByteBuilder();
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

//					float result=(parseHeader ? Parse.parseFloat(r1.id, "result=", '\t') : 0);
					processRead(r1, bb, bsw, vec);
					processRead(r2, bb, bsw, vec);
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

		if(bsw!=null) {errorState=bsw.poisonAndWait() | errorState;}
		
		if(histFile!=null) {
			bb.clear();
			PalindromeTracker.append(bb, "#count\tneg\tpos", new LongList[] {mhist, phist}, 101);
			ReadWrite.writeStringInThread(bb, histFile, false);
		}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Reads Out:          "+readsOut);
		assert(!errorState) : "An error was encountered.";
	}
	
	/*--------------------------------------------------------------*/
	
	private boolean processRead(Read r, ByteBuilder bb, ByteStreamWriter bsw, float[] vec) {
		if(r==null) {return false;}
		float result=(parseHeader ? Parse.parseFloat(r.id, "result=", '\t') : 0);
		float score=score(r.bases, vec, k, net, rcomp);
		int iscore=Tools.mid(0, Math.round(score*100), 100);
		boolean pass=(!filter ? true : (score>=cutoff)==highpass);
		readsOut+=(pass ? 1 : 0);
		if(result<0.5f) {mhist.increment(iscore);}
		else {phist.increment(iscore);}
		if(bsw!=null && pass) {
			if(annotate) {r.id+=String.format("\tscore=%.4f", score);}
			bsw.print(r.toFasta(bb.clear()).nl());
		}
		return pass;
	}
	
	public static float score(byte[] bases, float[] vec, int k, CellNet net, boolean rcomp) {
		SequenceToVector.fillVector(bases, vec, k);
		net.applyInput(vec);
		final float r;
		final float f=net.feedForward();
		if(rcomp) {
			AminoAcid.reverseComplementBasesInPlace(bases);
			SequenceToVector.fillVector(bases, vec, k);
			net.applyInput(vec);
			r=net.feedForward();
			AminoAcid.reverseComplementBasesInPlace(bases);
		}else {r=f;}
		return Tools.max(r, f);
	}
	
	public static float score(byte[] bases, float[] vec, int k, CellNet net) {
		assert(vec!=null);
		SequenceToVector.fillVector(bases, vec, k);
		net.applyInput(vec);
		final float f=net.feedForward();
		return f;
	}
	
	public static float score(byte[] bases, float[] vec, int k, CellNet net, int from, int to) {
		SequenceToVector.fillVector(bases, vec, k, from, to);
		net.applyInput(vec);
		final float f=net.feedForward();
		return f;
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String netFile=null;
	private String histFile=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffnet;

	private LongList phist=new LongList(101);
	private LongList mhist=new LongList(101);
	
	/*--------------------------------------------------------------*/

	private long readsOut=0;
	private long maxReads=-1;
	private boolean errorState=false;

	private int k=0;
	private boolean rcomp=false;
	private boolean parseHeader=false;
	private int width=-1;
	private final CellNet net;

	private boolean filter=false;
	private boolean highpass=true;
	private boolean annotate=true;
	private float cutoff=0.5f;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
