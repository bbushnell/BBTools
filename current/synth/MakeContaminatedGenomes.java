package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import dna.AminoAcid;
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
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;


/**
 * Fuses files together randomly to make chimeric genomes.
 * @author Brian Bushnell
 * @date Oct 7, 2014
 *
 */
public class MakeContaminatedGenomes {

	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		MakeContaminatedGenomes x=new MakeContaminatedGenomes(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs MakeContaminatedGenomes with command-line arguments.
	 * Parses input parameters, initializes file formats, and validates configuration.
	 * @param args Command-line arguments containing input files and options
	 */
	public MakeContaminatedGenomes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("chimeras") || a.equals("count")){
				chimeras=Integer.parseInt(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("exp")){
				exponent1=exponent2=Double.parseDouble(b);
			}else if(a.equals("exp1")){
				exponent1=Double.parseDouble(b);
			}else if(a.equals("exp2")){
				exponent2=Double.parseDouble(b);
			}else if(a.equals("delimiter")){
				delimiter=b;
			}else if(a.equals("regex")){
				regex=b;
			}else if(a.equals("subrate")){
				subRate=Double.parseDouble(b);
			}else if(a.equals("indelrate")){
				indelRate=Double.parseDouble(b);
			}else if(a.equals("id") || a.equals("ani") || a.equals("identity")){
				errorRate=Double.parseDouble(b);
				subRate=0.99*errorRate;
				indelRate=0.01*errorRate;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		errorRate=subRate+indelRate;
		
		{//Process parser fields
			Parser.processQuality();
			
			fofn=parser.in1;

			outPattern=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(fofn==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(outPattern!=null && outPattern.equalsIgnoreCase("null")){outPattern=null;}
		
		fffofn=FileFormat.testInput(fofn, FileFormat.TXT, null, true, true);
	}
	
	/**
	 * Main processing method that generates contaminated genomes.
	 * Creates specified number of chimeric sequences by randomly sampling and fusing input files.
	 * @param t Timer for tracking execution time and performance statistics
	 */
	void process(Timer t){
		final String[] in=TextFile.toStringLines(fffofn);
//		final long sizes[]=calcSizes(in);
		final Random randy=Shared.threadLocalRandom(seed);
		
		final StringBuilder sb=new StringBuilder();
		for(int cid=0; cid<chimeras; cid++){
			String s=makeOne(in, randy, cid);
			sb.append(s).append('\n');
		}
		if(outNames!=null){
			ReadWrite.writeString(sb, outNames);
		}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/* Calculates the size of each fasta file in an array */
//	private long[] calcSizes(String[] in){
//		long[] sizes=new long[in.length];
//		for(int i=0; i<in.length; i++){
//			ByteFile bf=ByteFile.makeByteFile(in[i], true);
//			long[] symbols=new long[255];
//			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
//				if(line.length>0 && line[0]!='>'){
//					assert(Tools.isLetter(line[0]));
//					for(byte b : line){
//						symbols[b]++;
//					}
//				}
//			}
//			final long defined=symbols['A']+symbols['C']+symbols['G']+symbols['T'];
//			final long total=shared.Vector.sum(symbols);
//			sizes[i]=total;
//		}
//		return sizes;
//	}
	
	/**
	 * Creates a single contaminated genome by fusing two randomly selected input files.
	 * Selects two different files, determines sampling fractions, and generates chimeric output.
	 *
	 * @param in Array of input file paths
	 * @param randy Random number generator for stochastic sampling
	 * @param cid Chimera identifier number
	 * @return Output file path of the generated contaminated genome
	 */
	String makeOne(String[] in, Random randy, int cid){
//		System.err.println("A");
		int a=randy.nextInt(in.length);
		int b=a;
		while(b==a){
			b=randy.nextInt(in.length);
		}
		double fracA=Math.pow(randy.nextDouble(), exponent1);
		double fracB=Math.pow(randy.nextDouble(), exponent2);
//		System.err.println("B: "+fracA+", "+fracB);

		FileFormat ffa=FileFormat.testInput(in[a], FileFormat.FASTA, null, true, true);
		FileFormat ffb=FileFormat.testInput(in[b], FileFormat.FASTA, null, true, true);
//		System.err.println("B1");

		ArrayList<Read> readsA=ConcurrentReadInputStream.getReads(-1, false, ffa, null, null, null);
//		System.err.println("B2: "+readsA.size());
		ArrayList<Read> readsB=ConcurrentReadInputStream.getReads(-1, false, ffb, null, null, null);
//		System.err.println("B3: "+readsB.size());
		
		return writeChimera(in[a], in[b], readsA, readsB, fracA, fracB, randy, cid);
	}
	
	/**
	 * Writes a contaminated genome by processing and combining reads from two sources.
	 * Applies sampling fractions to reads and writes output in order of dominance.
	 *
	 * @param inA First input file path
	 * @param inB Second input file path
	 * @param readsA Reads from first source
	 * @param readsB Reads from second source
	 * @param fracA Sampling fraction for first source
	 * @param fracB Sampling fraction for second source
	 * @param randy Random number generator for processing
	 * @param cid Chimera identifier number
	 * @return Output file path containing the contaminated genome
	 */
	String writeChimera(String inA, String inB, ArrayList<Read> readsA, ArrayList<Read> readsB, double fracA, double fracB, Random randy, int cid){
		ByteBuilder bb=new ByteBuilder();
		long sizeA=0, sizeB=0;
		for(Read r : readsA){
			readsProcessed++;
			basesProcessed+=r.length();
			processRead(r, bb, fracA, randy);
			sizeA+=r.length();
		}
//		System.err.println("D: "+sizeA);
		for(Read r : readsB){
			readsProcessed++;
			basesProcessed+=r.length();
			processRead(r, bb, fracB, randy);
			sizeB+=r.length();
		}
//		System.err.println("E: "+sizeB);
		
		final String out;
		if(fracA>=fracB){
//			System.err.println("F");
			out=outPattern.replaceFirst(regex, delimiter+sizeA+delimiter+Tools.format("%.3f", fracA)+delimiter+ReadWrite.stripToCore(inA)+delimiter+sizeB+delimiter+Tools.format("%.3f", fracB)+delimiter+ReadWrite.stripToCore(inB)+delimiter+cid+delimiter);
			ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, true);
			bsw.start();
//			System.err.println("G");
			for(Read r : readsA){bsw.println(r);}
			for(Read r : readsB){bsw.println(r);}
//			System.err.println("H");
			bsw.poisonAndWait();
		}else{
//			System.err.println("I");
			out=outPattern.replaceFirst(regex, delimiter+sizeB+delimiter+Tools.format("%.3f", fracB)+delimiter+ReadWrite.stripToCore(inB)+delimiter+sizeA+delimiter+Tools.format("%.3f", fracA)+delimiter+ReadWrite.stripToCore(inA)+delimiter+cid+delimiter);
			ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, true);
			bsw.start();
//			System.err.println("J");
			for(Read r : readsB){bsw.println(r);}
			for(Read r : readsA){bsw.println(r);}
//			System.err.println("K");
			bsw.poisonAndWait();
		}
//		System.err.println("L");
		return out;
	}
	
	/**
	 * Processes a single read by applying genome sampling and introducing mutations.
	 * Randomly samples a fraction of the genome and optionally adds substitutions and indels.
	 *
	 * @param r The read to process (modified in place)
	 * @param bb Byte builder for sequence manipulation
	 * @param genomeFraction Fraction of the genome to retain (0.0-1.0)
	 * @param randy Random number generator for stochastic operations
	 */
	public void processRead(Read r, ByteBuilder bb, double genomeFraction, Random randy){
		
		//Setup
		bb.clear();
		r.quality=null;
		
		long mutationsAdded=0;
		
		//Handle genomeFraction
		if(genomeFraction<1){
			final byte[] bases0=r.bases;
			int retain=(int)(bases0.length*(genomeFraction));
//			System.err.println("retain: "+retain);
			if(retain<bases0.length){
				final int start=randy.nextInt(bases0.length);
				int i=0, j=start;
				for(; i<retain && j<bases0.length; i++, j++){
					bb.append(bases0[j]);
				}
				j=0;
				
				if(i<retain){mutationsAdded++;} //Chimeric junction
				
				for(; i<retain; i++, j++){
					bb.append(bases0[j]);
				}
				r.bases=bb.toBytes();
				bb.clear();
			}
		}
		
		//Handle mutations
		//Not really the point of this tool but easy to add
		//Here, subs+indels=errors
		if(errorRate>0){
			final byte[] bases=r.bases;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				float x=randy.nextFloat();
				if(x<errorRate && AminoAcid.isFullyDefined(b)){
					mutationsAdded++;
					if(x<subRate){
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt(3)+1)&3)];
						bb.append(b);
					}else if(randy.nextBoolean()){//del
						//do nothing
					}else{//ins
						i--;
						b=AminoAcid.numberToBase[randy.nextInt(4)];
						bb.append(b);
					}
				}else{
					bb.append(b);
				}
			}
			//Modify read
			r.bases=bb.toBytes();
		}
//		if(prefix!=null){
//			r.id=prefix+r.numericID;
//		}
		basesRetained+=r.bases.length;
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** File of filenames containing input genome paths */
	private String fofn=null;

	/** Output filename pattern with placeholders for contamination metadata */
	private String outPattern=null;
	/** Output file for writing generated filename list */
	private String outNames=null;
	
	/** Number of contaminated genomes to generate */
	private int chimeras=1;
	/** Random seed for reproducible contamination patterns */
	private long seed=-1;
	/** Power law exponent for sampling first genome fraction */
	double exponent1=1;
	/** Power law exponent for sampling second genome fraction */
	double exponent2=1;
	/** Delimiter character for output filename metadata fields */
	String delimiter="_";
	/** Regular expression pattern for filename substitution */
	String regex="#";
	
	/** Rate of substitution mutations to introduce */
	double subRate=0;
	/** Rate of insertion/deletion mutations to introduce */
	double indelRate=0;
	/** Total error rate combining substitutions and indels */
	double errorRate=0;
	/** Total number of bases retained after processing */
	long basesRetained=0;

	/** Total number of reads processed */
	long readsProcessed=0;
	/** Total number of bases processed */
	long basesProcessed=0;
	
	/** Output stream for logging and status messages */
	private PrintStream outstream=System.err;
	
	/*--------------------------------------------------------------*/
	
	/** File format object for the input file of filenames */
	private final FileFormat fffofn;
	
	/*--------------------------------------------------------------*/
	
	/** Enables verbose output for debugging and detailed logging */
	public static boolean verbose=false;
	/** Indicates whether the program encountered errors during execution */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
