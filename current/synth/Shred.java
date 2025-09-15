package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class Shred {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		Shred x=new Shred(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Shred(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		Shared.capBufferLen(100);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		long seed=-1;
		Parser parser=new Parser();
		parser.minReadLength=parser.maxReadLength=-1;
		int increment_=0;
		boolean even=false;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("length") || a.equals("len") || a.equals("shredlen") || a.equals("shredlength")){
				shredLength=Parse.parseIntKMG(b);
			}else if(a.equals("overlap")){
				overlap=Parse.parseIntKMG(b);
			}else if(a.equals("increment")){
				increment_=Parse.parseIntKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("even") || a.equals("equal")){
				even=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("median")){
				median=Parse.parseIntKMG(b);
			}else if(a.equals("variance")){
				variance=Parse.parseIntKMG(b);
			}else if(a.equals("mode")){
				mode=Tools.find(b.toUpperCase(), modes);
			}else if(Tools.find(a.toUpperCase(), modes)>=0){
				mode=Tools.find(a.toUpperCase(), modes);
				assert(Parse.parseBoolean(b)) : "Modes should only be true: '"+arg+"'";
			}else if(a.equals("filetid")){
				parseFileTID=Parse.parseBoolean(b);
			}else if(a.equals("sequencetid") || a.equals("headertid")){
				parseSequenceTID=Parse.parseBoolean(b);
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		evenLengths=even;
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;

			minLength=parser.minReadLength;
			maxLength=parser.maxReadLength;
			maxNs=parser.maxNs;
		}
		
		minLength=Tools.mid(1, minLength, shredLength);
		if(increment_>0) {
			assert(shredLength>0);
			increment=increment_;
			overlap=(maxLength<0 ? shredLength : minLength)-increment;
		} else {
			assert(shredLength>0);
			assert(shredLength>overlap);
			increment=shredLength-overlap;
		}
		assert(increment>0);
		incMult=1.0/increment;
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		
		randy=Shared.threadLocalRandom(seed);
		if(median>0 && variance>0 && (minLength<0 || maxLength<0)) {
			minLength=Tools.max(1, median-variance);
			maxLength=median+variance;
		}
		System.err.println("minlen="+minLength+", maxlen="+maxLength);
		if(minLength>-1 && maxLength>-1) {
			assert(maxLength>=minLength);
			range=maxLength-minLength+1;
		}else {
			range=0;
		}
	}
	
	public boolean parseArgument(String arg, String a, String b){
		if(a.equals("reads") || a.equals("maxreads")){
			maxReads=Parse.parseKMG(b);
			return true;
		}else if(a.equals("some_argument")){
			maxReads=Parse.parseKMG(b);
			return true;
		}
		return false;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=2;
			
			if(cris.paired()){KillSwitch.kill("This program does not support paired reads.");}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;
		
		readsOut=0;
		basesOut=0;
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}
		
		int tid=-1;
		if(parseFileTID) {
			tid=bin.BinObject.parseTaxID(in1);
		}

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
			
			ArrayList<Read> listOut=new ArrayList<Read>();
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);

				final int initialLength1=r1.length();

				readsProcessed++;
				basesProcessed+=initialLength1;
				
				if(parseSequenceTID) {
					tid=bin.BinObject.parseTaxID(r1.id);
				}
				
				if(range>1){
					processRandomly(r1, listOut, tid);
				}else if(evenLengths){
					processEvenly(r1, listOut, tid);
				}else{
					processUnevenly(r1, listOut, tid);
				}
			}

			if(ros!=null){ros.add(listOut, ln.id);}

			cris.returnList(ln);
			if(verbose){outstream.println("Returned a list.");}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	void processRead(final Read r1, final ArrayList<Read> list, int tid){
//		if(r1.length()<minLength){return;}
//		if(r1.length()<=shredLength){
//			r1.numericID=readsOut;
//			list.add(r1);
//			readsOut++;
//			basesOut+=r1.length();
//			return;
//		}
//		if(evenLengths){
//			processEvenly(r1, list, tid);
//		}else{
//			processUnevenly(r1, list, tid);
//		}
//	}
	
	void processEvenly(final Read r1, final ArrayList<Read> list, int tid){
		final byte[] bases=r1.bases;
		final byte[] quals=r1.quality;
		final String name=(prefix==null ? r1.id : prefix+r1.numericID);
		
		final int chunks=(int)Math.ceil((bases.length-overlap)*incMult);
		assert(chunks>0);
		double inc2=bases.length/(double)chunks;
		
		for(int chunk=0; chunk<chunks; chunk++){
			int a=(int)Math.floor(inc2*chunk);
			int b=(chunk==chunks-1 ? bases.length : overlap+(int)Math.floor(inc2*(chunk+1)));
			b=Tools.min(b, a+shredLength);
			final int length=b-a;
			if(length<minLength){return;}
			final byte[] bases2=KillSwitch.copyOfRange(bases, a, b);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, a, b));
			Read shred=new Read(bases2, quals2, toName(name, a, (b-1), tid), readsOut);
			boolean pass=true;
			if(pass && maxNs>=0 && shred.countUndefined()>maxNs) {pass=false;}
			if(pass) {
				readsOut++;
				basesOut+=shred.length();
				list.add(shred);
			}
		}
	}
	
	void processUnevenly(final Read r1, final ArrayList<Read> list, int tid){
		final byte[] bases=r1.bases;
		final byte[] quals=r1.quality;
		final String name=(prefix==null ? r1.id : prefix+r1.numericID);
		for(int i=0; i<bases.length; i+=increment){
			final int limit=Tools.min(i+shredLength, bases.length);
			final int length=limit-i;
			if(length<minLength){return;}
			final byte[] bases2=KillSwitch.copyOfRange(bases, i, limit);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, i, limit));
			Read shred=new Read(bases2, quals2, toName(name, i, (limit-1), tid), readsOut, r1.flags);
			boolean pass=true;
			if(pass && maxNs>=0 && shred.countUndefined()>maxNs) {pass=false;}
			if(pass) {
				readsOut++;
				basesOut+=shred.length();
				list.add(shred);
			}
			if(limit==bases.length){return;}
			assert(limit<bases.length);
		}
	}
	
	void processRandomly(final Read r1, final ArrayList<Read> list, int tid){
		final byte[] bases=r1.bases;
		final byte[] quals=r1.quality;
		final String name=(prefix==null ? r1.id : prefix+r1.numericID);
		for(int i=0; i<bases.length;){
//			int rand=Tools.min(randy.nextInt(2*variance), randy.nextInt(3*variance), 2*variance);
//			final int limit=Tools.max(i+minLength, Tools.min(i+rand+median-variance, bases.length));
//			final int rand=randy.nextInt(range);
//			final int limit=Tools.min(i+rand+minLength, bases.length);
//			final int length=limit-i;
			final int length=randomLength(bases.length-i);
			if(length<minLength){return;}
			final int limit=i+length;
			assert(limit<=bases.length) : "i="+i+", length="+length+", limit="+limit+", blen="+bases.length;
			final byte[] bases2=KillSwitch.copyOfRange(bases, i, limit);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, i, limit));
			Read shred=new Read(bases2, quals2, toName(name, i, (limit-1), tid), readsOut);
			readsOut++;
			basesOut+=shred.length();
			list.add(shred);
			if(limit==bases.length){return;}
			assert(limit<bases.length);
			i=limit;
		}
	}
	
	final String toName(String name, int start, int stop, int tid) {
		return (name==null ? "" : name+"_")+start+"-"+stop+(tid>0 ? "_tid_"+tid : "");
	}
	
	final int randomLength(final int remainder) {
		if(mode==LINEAR) {return randomLengthLinear(remainder);}
		else if(mode==LOG) {return logUniformLength(remainder);}
		else if(mode==EXP) {return randomLengthExp(remainder);}
		else {
			throw new RuntimeException("Unknown mode "+mode);
		}
	}
	
	final int randomLengthLinear(final int remainder) {
		return Tools.min(minLength+randy.nextInt(range), remainder);
	}
	
	final int randomLengthExp(final int remainder) {
		if(remainder<=minLength) {return remainder;}
		long max=Tools.min(maxLength, remainder);
		float invlambda=(float)(1/Math.sqrt(minLength*max));
//		System.err.println("1/Lambda="+invlambda);
		int rand=(int)Tools.exponential(randy, invlambda);
		for(int i=0; i<20 && (rand>max || rand<minLength); i++) {
			rand=(int)Tools.exponential(randy, invlambda);
		}
//		System.err.println("Rand.exp="+rand);
		return (rand<minLength || rand>max ? remainder : rand);
	}
	
    /**
     * Samples a length from a continuous log-uniform distribution.
     * @author Isla
     */
    private int logUniformLength(int remainder) {
        double u = randy.nextDouble(); // Uniform between 0 and 1
        
        // Calculate log-uniform value between minLength and maxLength
        double logMin = Math.log(minLength);
        double logMax = Math.log(maxLength);
        double logValue = logMin + u * (logMax - logMin);
        
        // Convert back from log space and round to integer
        int len=(int) Math.round(Math.exp(logValue));
        return Math.min(remainder, len);
    }
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	boolean parseFileTID=false;
	boolean parseSequenceTID=false;
	String prefix=null;
	
	/*--------------------------------------------------------------*/

	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long readsOut=0;
	protected long basesOut=0;
	
	private long maxReads=-1;
	
	private int median=-1;
	private int variance=-1;
	
	private int shredLength=500;
	private int minLength=-1;
	private int maxLength=-1;
	private int maxNs=-1;
	private final int range;
	private int overlap=0;
	private final int increment;
	private final double incMult;
	
	private final boolean evenLengths;
	
	private final Random randy;
	
	int mode=LINEAR;
	static final int LINEAR=0, LOG=1, EXP=2;
	static final String[] modes={"LINEAR", "LOG", "EXP"};
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
