package jgi;

import java.io.PrintStream;

import dna.AminoAcid;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tracker.ReadStats;

/**
 * Generates synthetic repeated-motif FASTA records.
 * For each k in the requested range, this emits one record for each nucleotide
 * k-mer and repeats that k-base motif to the smallest multiple of k at least
 * minlen+k-1, exposing every starting phase of a minlen-base window in the motif.
 *
 * @author Brian Bushnell
 * @date Oct 17, 2014
 */
public class MakePolymers{
	
	/** Construct, validate parameters, generate records, then close redirected diagnostics on success. */
	public static void main(String[] args){
		Timer t=new Timer();
		MakePolymers x=new MakePolymers(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parse settings, validate them, and initialize the primary output descriptor. */
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
			String[] split=Parse.splitOnFirst(arg, '=');
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
		if(out1==null || out1.trim().isEmpty()){
			throw new IllegalArgumentException("An output file is required; set out=<file> or out=stdout.fa.");
		}
		
		validateParams();
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);
	}
	
	/** Generate all requested records, wait for output completion, and report final counts. */
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
	
	/** Write every repeated k-base motif for this k. */
	private void writeSequence(int k, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		
		final int minCount=repeatCount(k, minLen);
		final long max=maxKey(k);
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
	
	/** Append one 2-bit encoded nucleotide motif of length k. */
	public static final void toBytes(long kmer, int k, ByteBuilder bb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			bb.append(AminoAcid.numberToBase[x]);
		}
	}
	
	/** Return the inclusive maximum 2-bit key for a requested polymer length. */
	public static final long maxKey(final int k){
		assert(k>0 && k<=MAX_K) : "MakePolymers stores 2-bit kmers in a positive signed long; k="+k+", max="+MAX_K;
		return (1L<<(2*k))-1;
	}
	
	/** Return the repeat count needed to expose every minLen-base window phase. */
	public static final int repeatCount(final int k, final int minLen){
		assert(k>0) : "MakePolymers validates positive k before repeatCount; k="+k;
		assert(minLen>0) : "MakePolymers validates positive minLen before repeatCount; minLen="+minLen;
		final long target=(long)minLen+k-1;
		return (int)((target+k-1)/k);
	}
	
	/**
	 * Return the largest FASTA record ByteBuilder will need to hold for one key.
	 * Near-MAX records start from zero because same-k sequence lengths are constant
	 * and each large completed record triggers the existing 16384-byte flush.
	 */
	public static final long maxRecordBytes(final int k, final int minLen){
		final long seqBases=repeatCount(k, minLen)*(long)k;
		return seqBases+digits(k)+digits(maxKey(k))+FASTA_RECORD_OVERHEAD;
	}
	
	private static final int digits(final long x){
		assert(x>=0) : "MakePolymers FASTA record size formats only nonnegative k and kmer ids; x="+x;
		int digits=1;
		for(long y=x; y>=10; y/=10){digits++;}
		return digits;
	}
	
	/** Validate numeric bounds before opening the primary output file. */
	private void validateParams(){
		if(minLen<1){
			throw new IllegalArgumentException("minlen must be at least 1: "+minLen);
		}
		if(mink<1 || maxk<1 || mink>maxk || maxk>MAX_K){
			throw new IllegalArgumentException("Require 1<=mink<=maxk<="+MAX_K+"; mink="+mink+", maxk="+maxk);
		}
		for(int k=mink; k<=maxk; k++){
			final long recordBytes=maxRecordBytes(k, minLen);
			if(recordBytes>ByteBuilder.MAXLEN){
				throw new IllegalArgumentException("FASTA record for k="+k+", minlen="+minLen+
					" needs "+recordBytes+" bytes; max is "+ByteBuilder.MAXLEN);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	public static final int MAX_K=31;
	/** Fixed bytes for '>', '_', header newline, and sequence newline. */
	private static final int FASTA_RECORD_OVERHEAD=4;
	
	/** Number of FASTA records generated. */
	private long readsProcessed=0;
	/** Number of sequence bases generated. */
	private long basesProcessed=0;
	
	/** Minimum and maximum repeated motif length; valid values are 1 through MAX_K. */
	private int mink=1, maxk=1;
	
	/** Minimum requested sequence length before adding k-1 bases and rounding up to a multiple of k. */
	private int minLen=31;
	
	/** Required output filename. */
	private String out1=null;
	/** Primary output descriptor. */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	/** Diagnostic output stream, possibly redirected by PreParser. */
	private PrintStream outstream=System.err;
	/** Print verbose parser and I/O diagnostics. */
	public static boolean verbose=false;
	/** True if any writer reports a terminal error. */
	public boolean errorState=false;
	/** Allow overwriting an existing output file. */
	private boolean overwrite=true;
	/** Append to an existing output file instead of replacing it. */
	private boolean append=false;
	
}
