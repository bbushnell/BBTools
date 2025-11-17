package jgi;

import java.io.PrintStream;

import dna.AminoAcid;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class MakePolymers {
	
	/** Program entry point for generating k-mer polymer sequences.
	 * @param args Command-line arguments including k-mer range and output settings */
	public static void main(String[] args){
		Timer t=new Timer();
		MakePolymers x=new MakePolymers(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes the polymer generator.
	 * Processes parameters for k-mer range (mink, maxk), minimum sequence length,
	 * and output file settings.
	 * @param args Command-line arguments array
	 */
	public MakePolymers(String[] args){
		
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
			}else if(a.equals("k")){
				mink=maxk=Integer.parseInt(b);
			}else if(a.equals("mink")){
				mink=Integer.parseInt(b);
			}else if(a.equals("maxk")){
				maxk=Integer.parseInt(b);
			}else if(a.equals("len") || a.equals("minlen")){
				minLen=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);
	}
	
	/**
	 * Main processing method that generates polymer sequences for each k-mer length.
	 * Creates output writer, iterates through k-mer sizes from mink to maxk,
	 * and writes sequences containing all possible k-mers of each length.
	 * @param t Timer for tracking execution time and reporting performance
	 */
	void process(Timer t){

		final ByteStreamWriter bsw;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		for(int i=mink; i<=maxk; i++){
			writeSequence(i, bsw);
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Writes polymer sequences containing all possible k-mers of length k.
	 * Generates sequences by enumerating all 4^k possible k-mers and repeating
	 * each k-mer enough times to meet the minimum sequence length requirement.
	 *
	 * @param k Length of k-mers to generate sequences for
	 * @param bsw Output writer for FASTA formatted sequences
	 */
	private void writeSequence(int k, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		
		final int minLen2=((minLen+k-1)/k)*k;
		final int minCount;
		if(minLen2-minLen>=k-1){
			minCount=minLen2/k;
		}else{
			minCount=minLen2/k+1;
		}
		
		final long max=(1<<(2*k))-1;
		for(long kmer=0; kmer<=max; kmer++){
			bb.append('>').append(k).append('_').append(kmer).append('\n');
			for(int i=0; i<minCount; i++){
				basesProcessed+=k;
				toBytes(kmer, k, bb);
			}
			readsProcessed++;
			bb.nl();
			if(bb.length>=16384){
				bsw.print(bb);
				bb.clear();
			}
		}
		if(bb.length>0){
			bsw.print(bb);
			bb.clear();
		}
	}
	
	/**
	 * Converts a k-mer encoded as a long integer to DNA sequence bytes.
	 * Decodes the 2-bit per base encoding back to ACGT nucleotide characters.
	 *
	 * @param kmer K-mer encoded as long with 2 bits per base
	 * @param k Length of the k-mer to decode
	 * @param bb ByteBuilder to append the decoded sequence to
	 */
	public static final void toBytes(long kmer, int k, ByteBuilder bb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			bb.append(AminoAcid.numberToBase[x]);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Number of sequences (reads) generated */
	private long readsProcessed=0;
	/** Total number of bases written to output */
	private long basesProcessed=0;
	
	/** Maximum k-mer length to generate sequences for */
	/** Minimum k-mer length to generate sequences for */
	private int mink=1, maxk=1;
	
	/** Minimum length of generated polymer sequences */
	private int minLen=31;
	
	/** Output file path for generated FASTA sequences */
	private String out1=null;

	/** File format specification for output file */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and detailed logging */
	public static boolean verbose=false;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
