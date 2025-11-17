package synth;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class SynthMDA {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point for MDA simulation.
	 * @param args Command-line arguments for MDA configuration */
	public static void main(String[] args){
		Timer t=new Timer();
		SynthMDA x=new SynthMDA(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs SynthMDA instance and parses command-line arguments.
	 * Configures input/output files, amplification parameters, and validation settings.
	 * @param args Command-line arguments containing MDA simulation parameters
	 */
	public SynthMDA(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		parser.build=7;
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
			}else if(a.equals("minlen") || a.equals("ml")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("maxlen") || a.equals("mxl")){
				maxlen=Integer.parseInt(b);
			}else if(a.equals("cycles")){
				cycles=Integer.parseInt(b);
			}else if(a.equals("initialratio")){
				initialRatio=Float.parseFloat(b);
			}else if(a.equals("ratio")){
				ratio=Float.parseFloat(b);
			}else if(a.equals("refout")){
				out1=b;
			}else if(a.equals("perfect")){
				perfectrate=Float.parseFloat(b);
			}else if(a.equals("length")){
				readlength=Integer.parseInt(b);
			}else if(a.equals("paired")){
				paired=Parse.parseBoolean(b);
			}else if(a.equals("amp")){
				amp=Integer.parseInt(b);
			}
//			else if(a.equals("build")){
//				assert(false) : "Build should have been parsed by parser.";
//				build=Integer.parseInt(b);
//			}
			else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		minlen2=Tools.min(minlen2, minlen);
		
		{//Process parser fields
			Parser.processQuality();
			
			if(parser.maxReads>0){reads=parser.maxReads;}
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			if(ref==null){ref=parser.in1;}

			readsOut=parser.out1;
			
			extref=parser.extin;
			extout=parser.extout;
			build=parser.build;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(ref==null){throw new RuntimeException("Error - input reference must be specified.");}
		
		if(out1==null){
			out1=ReadWrite.stripToCore(ref)+"_"+Long.toHexString(new Random().nextLong()&Long.MAX_VALUE)+".fa";
		}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffref=FileFormat.testInput(ref, FileFormat.FASTQ, extref, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that executes the MDA simulation pipeline.
	 * Reads input reference sequences, performs iterative amplification cycles,
	 * writes amplified reference output, and optionally generates synthetic reads.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		ByteBuilder bb=new ByteBuilder();
		bb.append('$');
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ffref, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		assert(!cris.paired());
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffref==null || ffref.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					
					bb.append(r1.bases);
					bb.append('$');
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}
				
				final ArrayList<Read> listOut=reads;

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStream(cris);
		
		ByteBuilder dest=amplify(bb, false, minlen, maxlen, initialRatio);
		bb=null;
		for(int i=0; i<cycles; i++){
			dest=amplify(dest, i<1, minlen, maxlen, ratio);
//			if(dest.length()*ratio>1500000000){break;}
		}
//		assert(false) : cycles+", "+dest.length();
		
		TextStreamWriter tsw=(ffout1==null ? null : new TextStreamWriter(ffout1));
		if(tsw!=null){tsw.start();}
		
		bb=new ByteBuilder();
		for(int i=0, id=1; i<dest.length(); i++){
			byte b=dest.get(i);
			if(b=='$'){
				if(bb.length()>0){
					tsw.print(">"+id+"\n");
					tsw.println(bb.toString());
					id++;
				}
				bb.setLength(0);
			}else{
				bb.append(b);
			}
		}
		dest=null;
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(readsOut!=null){
			FileFormat ff=FileFormat.testOutput(readsOut, FileFormat.FASTQ, null, true, overwrite, false, false);
			assert(ff!=null);
			ArrayList<String> list=new ArrayList<String>();
			list.add("reads="+reads);
			list.add("length="+readlength);
			list.add("amp="+amp);
			if(paired){
				list.add("paired="+paired);
				list.add("interleaved="+paired);
			}
			list.add("build="+build);
			list.add("out="+readsOut);
			list.add("ow="+overwrite);
			list.add("minq="+16);
			list.add("midq="+25);
			list.add("maxq="+38);
			list.add("adderrors");
			list.add("snprate="+0.02);
			list.add("delrate="+0.005);
			list.add("insrate="+0.005);
			list.add("nrate="+0.005);
			list.add("maxinslen="+3);
			list.add("maxdellen="+3);
			list.add("maxnlen="+3);
			list.add("maxinss="+2);
			list.add("maxdels="+2);
			list.add("maxns="+2);
			list.add("maxsnps="+2);
			list.add("seed=-1");
			list.add("ref="+out1);
			if(prefix!=null){list.add("prefix="+prefix);}
			if(perfectrate>0){
				list.add("perfect="+perfectrate);
			}
			RandomReads3.main(list.toArray(new String[list.size()]));
		}
		
		boolean deleteRef=(readsOut!=null);
		if(deleteRef){
			if(verbose){System.err.println("Trying to delete "+out1);}
			try {
				File f=new File(out1);
				if(f.exists()){f.delete();}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Performs single amplification cycle by random fragment sampling.
	 * Randomly selects fragments from source sequence in forward or reverse orientation,
	 * amplifies content by the specified ratio until target size is reached.
	 *
	 * @param source Source sequence data to amplify from
	 * @param retain Whether to retain original source data in output
	 * @param minlen Minimum fragment length for amplification
	 * @param maxlen Maximum fragment length for amplification
	 * @param ratio Target amplification factor for this cycle
	 * @return ByteBuilder containing amplified sequence data
	 */
	private ByteBuilder amplify(ByteBuilder source, boolean retain, int minlen, int maxlen, float ratio){
		assert(minlen<=maxlen && minlen>0 && maxlen>0);
		final int range=maxlen-minlen+1;
		final int slen=source.length();
		if(slen<minlen2*1.1f){
			KillSwitch.kill("Input ("+slen+") must be at least 10% longer than minlen ("+minlen2+").");
		}
		if(source.length()>=500000000){retain=false;}
		ByteBuilder dest=(retain ? source : new ByteBuilder());
		int goal=(int)Tools.min((long)(slen*ratio), 600000000);
		while(dest.length()<goal){
			final long initialLength=dest.length();
			final int start=randy.nextInt(slen);
			final int len0=minlen+randy.nextInt(range);
			final boolean forward=Tools.nextBoolean(randy);
			if(initialLength+(long)len0>1500000000){break;}
//			System.err.println(forward+", "+start+", "+len0);
			if(forward){
				final int stop=Tools.min(source.length(), start+len0);
//				System.err.println("stop="+stop);
				for(int i=start; i<stop; i++){
					byte b=source.get(i);
					if(b=='$'){
//						System.err.println("b="+(char)b);
						break;
					}
					dest.append(b);
				}
			}else{
				final int stop=Tools.max(0, start-len0);
//				System.err.println("stop="+stop);
				for(int i=start; i>=stop; i--){
					byte b=source.get(i);
					if(b=='$'){
//						System.err.println("b="+(char)b);
						break;
					}
					dest.append(AminoAcid.baseToComplementExtended[b]);
				}
			}
			dest.append('$');
			long added=dest.length()-initialLength;
//			System.err.println("added "+added+"/"+len0+" ("+initialLength+" -> "+dest.length()+")");
//			if(added<Tools.min(200, minlen) || (added<Tools.min(1000, minlen) && Tools.nextBoolean(randy))){dest.setLength(initialLength);}
			if(added<Tools.min(minlen2, minlen)){dest.setLength((int)initialLength);}
		}
		return dest;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input reference file path */
	private String ref=null;
	/** Output amplified reference file path */
	private String out1=null;
	
	/** File extension override for reference input */
	private String extref=null;
	/** File extension override for amplified reference output */
	private String extout=null;
	
	/** File format handler for reference input */
	private final FileFormat ffref;
	/** File format handler for amplified reference output */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Minimum fragment length for amplification sampling */
	private int minlen=10000;
	/** Secondary minimum length threshold for fragment retention */
	private int minlen2=4000;
	/** Maximum fragment length for amplification sampling */
	private int maxlen=150000;
	/** Number of amplification cycles to perform */
	private int cycles=9;
	/** Amplification ratio for the first cycle */
	private float initialRatio=1.3f;
	/** Amplification ratio for subsequent cycles */
	private float ratio=2;
	
	/** Prefix for generated sequence identifiers */
	private String prefix=null;
	
	/*--------------------------------------------------------------*/
	
	/** Number of synthetic reads to generate if specified */
	private long reads=12000000;
	/** Length of synthetic reads to generate */
	private int readlength=150;
	/** Amplification factor for synthetic read generation */
	private int amp=200;
	/** Whether to generate paired-end synthetic reads */
	private boolean paired=true;
	/** Genome build version for synthetic read generation */
	private int build=7;
	/** Output file path for synthetic reads */
	private String readsOut=null;
	/** Fraction of synthetic reads to generate without errors */
	private float perfectrate=0;
	
	/** Thread-local random number generator for fragment sampling */
	private final Random randy=Shared.threadLocalRandom();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Global flag controlling verbose output across components */
	public static boolean verbose=false;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	
}
