package jgi;

import java.io.PrintStream;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Reads a fasta/fastq file, encodes every kmer as a 2-bit-packed long via the
 * standard shift-mask approach, anonymizes it with Tools.hash64shift, and
 * dumps the resulting hashcodes one per line (A48-encoded) so the sequence
 * content cannot be recovered from the output - only its statistical
 * structure (repeats, duplicates, complexity) survives.
 * Built for Amber's cardinality-estimator work, which needs realistic hash
 * streams but may not touch or reason about biological sequence directly.
 *
 * @author UMP45
 * @date July 5, 2026
 */
public class KmerHashDump {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) {
		Timer t=new Timer();
		KmerHashDump x=new KmerHashDump(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public KmerHashDump(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
		}

		assert(k>0 && k<=31) : "k must be 1-31 to fit in a 2-bit-packed long; k="+k;

		if(in1==null){throw new RuntimeException("Error - an input file is required.");}
		if(out1==null){throw new RuntimeException("Error - an output file (out=) is required.");}
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);

		shift=2*k;
		mask=(shift>=64 ? -1L : ~((-1L)<<shift));
	}

	/** Parse arguments from the command line. */
	private Parser parse(String[] args){

		Parser parser=new Parser();

		for(int i=0; i<args.length; i++){
			String arg=args[i];

			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null && parser.out1==null){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private void process(Timer t) {
		Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, true, maxReads,
			false, false, 1);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		st.start();
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
			for(Read r : ln) {
				dumpKmerHashes(r, bb);
				if(r.mate!=null){dumpKmerHashes(r.mate, bb);}
				if(bb.length()>=BUFFER_LIMIT){
					bsw.print(bb);
					bb=new ByteBuilder();
				}
			}
		}
		st.close();
		if(bb.length()>0){bsw.print(bb);}
		errorState|=bsw.poisonAndWait();

		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		System.err.println("Kmers written:  \t"+kmersOut);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Slide a k-mer window across a single read using the standard shift-mask
	 * encoding, hash each kmer with Tools.hash64shift, and append the
	 * A48-encoded hash (one per line) to the buffer.
	 * @param r Read to scan
	 * @param bb Output buffer
	 */
	private void dumpKmerHashes(Read r, ByteBuilder bb) {
		readsIn++;
		final byte[] bases=r.bases;
		if(bases==null){return;}
		basesIn+=bases.length;

		long kmer=0;
		int len=0;
		for(int i=0, blen=bases.length; i<blen; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				continue;
			}
			kmer=((kmer<<2)|x)&mask;
			len++;
			if(len>=k){
				final long hash=Tools.hash64shift(kmer);
				bb.appendA48(hash).nl();
				kmersOut++;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String in2=null;
	private String out1=null;

	private FileFormat ffin1;
	private FileFormat ffin2;
	private FileFormat ffout1;

	/** Kmer length; must be 1-31 so the 2-bit-packed kmer fits in a long. */
	private int k=31;
	/** 2*k; number of bits used in the packed kmer. */
	private final int shift;
	/** Mask selecting the low 'shift' bits, isolating the current kmer. */
	private final long mask;

	/** Quit after processing this many input reads; -1 means no limit. */
	private long maxReads=-1;
	/** Flush the output buffer once it reaches this many bytes. */
	private static final int BUFFER_LIMIT=1<<20;

	private long readsIn=0;
	private long basesIn=0;
	private long kmersOut=0;

	private boolean overwrite=false;
	private boolean append=false;
	private boolean errorState=false;

	private PrintStream outstream=System.err;
	public static boolean verbose=false;

}
