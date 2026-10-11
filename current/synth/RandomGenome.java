package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import shared.Random;

import bin.AdjustEntropy;
import clade.Clade;
import clade.CladeLoader;
import clade.CladeObject;
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
import structures.ByteBuilder;

/**
 * Writes random DNA or amino acid records in wrapped FASTA form.
 * Each record receives totalLength/chroms generated symbols using integer division;
 * any remainder is omitted. Padding is added separately at both ends as N for DNA
 * or X for amino acids. Record IDs are one-based contig or gene numbers.
 *
 * DNA sampling uses the configured GC fraction or, when an input is supplied,
 * conditional base probabilities derived from its kmer counts. Input-derived sampling
 * cannot be combined with noPoly; GC-only and amino sampling honor that option for
 * generated symbols. Padding can contain repeated symbols regardless of noPoly.
 *
 * Generation runs on the calling thread and each invocation creates a writer thread.
 * The RNG and accumulated writer error state belong to the instance and persist across
 * invocations. Construction also changes shared I/O settings; concurrent use of an
 * instance is unsupported. Writer completion errors are reported by process.
 *
 * @author Brian Bushnell
 * @date Jan 3, 2013
 */
public class RandomGenome {
	
	/** Runs generation and closes redirected diagnostics on normal completion.
	 * @param args Tool and shared BBTools command-line options
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		RandomGenome x=new RandomGenome(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parses output dimensions and sampling options, optionally loading an input model.
	 * Shared parser and compression settings may be changed. Negative padding is clamped
	 * to zero; assertions check wrapping, chromosome count and per-record length.
	 * For clade-file input only the first loaded clade is used. Numerical sampling ranges
	 * and the compatibility of every option combination are not comprehensively checked.
	 * @param args Tool options and shared BBTools arguments, including config-file options
	 */
	public RandomGenome(String[] args){
		
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
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("chroms")){
				chroms=Parse.parseIntKMG(b);
			}else if(a.equals("len") || a.equals("length") || a.equals("size")){
				totalLength=Parse.parseKMG(b);
			}else if(a.equals("pad")){
				pad=Tools.max(0, Parse.parseIntKMG(b));
			}else if(a.equals("gc")){
				gc=Float.parseFloat(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("nohomopolymers") || a.equals("banhomopolymers") || a.equals("nopoly")){
				noPoly=Parse.parseBoolean(b);
			}else if(a.equals("includestop") || a.equals("stop")){
				includeStop=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
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
			in=parser.in1;
		}

		wrap=Shared.FASTA_WRAP;
		assert(wrap>0) : "Wrap is too small.";
		assert(chroms>0) : "Chroms must be greater than 0.";
		assert(totalLength>=chroms) : "Length must be at least chroms.";
		assert(totalLength/chroms<Shared.MAX_ARRAY_LEN-2L*pad) : "Length per chrom must be at most "+Shared.MAX_ARRAY_LEN;
		chromLength=(int)(totalLength/chroms);

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		if(in!=null && noPoly){
			throw new IllegalArgumentException("nopoly is incompatible with input-derived generation; omit in= or set nopoly=f");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffin=FileFormat.testInput(in, FileFormat.CLADE, null, true, false, false);
		ffout=FileFormat.testOutput(out, FileFormat.FA, null, true, overwrite, append, false);
		
		if(ffin!=null){
			Clade clade;
			AdjustEntropy.load();
			if(ffin.clade()){
				ArrayList<Clade> clades=CladeLoader.loadCladesFromClade(ffin);
				clade=clades.get(0);
			}else{
				clade=CladeLoader.loadCladeFromSequence(ffin);
			}
			long[] counts=clade.counts[k];
			if(k>2){counts=unfold(counts, k);}
			prefixMatrix=countsToPrefixProb(counts, k);
		}else{
			prefixMatrix=null;
		}
		
		randy=Shared.threadLocalRandom(seed);
	}
	
	/** Selects the current shared amino/DNA mode, waits for output and reports writer errors.
	 * RNG position and prior errorState are not reset between invocations.
	 * @param t Caller timer; currently passed through without being stopped or reported
	 */
	void process(Timer t){
		if(Shared.AMINO_IN){
			processAmino(t);
		}else{
			processNucleotide(t);
		}
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Writes contig records using GC-only or input-derived conditional base sampling.
	 * Each input-model record starts with the all-A encoded prefix; N padding resets
	 * that prefix. GC is not applied in the input-model branch.
	 * Waits for the writer and accumulates its returned error state.
	 * @param t Caller timer, currently unused
	 */
	void processNucleotide(Timer t){
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		
		for(long chrom=1; chrom<=chroms; chrom++){
			bsw.print('>').print("contig").println(chrom);
			ByteBuilder bb=new ByteBuilder(wrap+1);
			byte prev='N';
			final int max=chromLength+2*pad;
			final int pad2=chromLength+pad;
			if(prefixMatrix!=null){
				if(prefixMatrix!=null){
					final int mask=(1<<(2*(k-1)))-1; // Precalculate mask
					int prefix=0; // Start with empty (k-1)-mer
					for(int i=0; i<max; ){
						for(int j=0; j<wrap && i<max; i++, j++){
							byte b;
							if(i<pad || i>=pad2){
								b='N';
								prefix=0; // Reset on N
							}else{
								b=nextBase(prefix, prefixMatrix, randy);
								// Update prefix: shift left, add new base, mask to k-1 length
								prefix=((prefix<<2)|AminoAcid.baseToNumber[b])&mask;
							}
							bb.append(b);
							prev=b;
						}
						bb.nl();
						bsw.print(bb);
						bb.clear();
					}
				}
			}else if(gc==0.5f){
				for(int i=0; i<max; ){
					for(int j=0; j<wrap && i<max; i++, j++){
						byte b;
						if(i<pad || i>=pad2){b='N';}
						else{
							b=AminoAcid.numberToBase[randy.nextInt()&3];
							while(noPoly && b==prev){b=AminoAcid.numberToBase[randy.nextInt()&3];}
						}
						bb.append(b);
						prev=b;
					}
					bb.nl();
					bsw.print(bb);
					bb.clear();
				}
			}else{
				for(int i=0; i<max; ){
					for(int j=0; j<wrap && i<max; i++, j++){
						boolean at=randy.nextFloat()>=gc;
						char b;
						if(i<pad || i>=pad2){b='N';}
						else{
							boolean low=randy.nextBoolean();
							if(at){
								b=low ? 'A' : 'T';
							}else{
								b=low ? 'C' : 'G';
							}
							while(noPoly && b==prev){
								low=randy.nextBoolean();
								if(at){
									b=low ? 'A' : 'T';
								}else{
									b=low ? 'C' : 'G';
								}
							}
						}
						bb.append(b);
						prev=(byte)b;
					}
					bb.nl();
					bsw.print(bb);
					bb.clear();
				}
			}
		}
		errorState|=bsw.poisonAndWait();
	}
	
	/** Writes gene records sampled uniformly from the configured amino acid alphabet.
	 * includeStop adds the stop symbol; noPoly rejects repeated generated symbols.
	 * GC and the optional input-derived DNA model do not affect this branch.
	 * Waits for the writer and accumulates its returned error state.
	 * @param t Caller timer, currently unused
	 */
	void processAmino(Timer t){
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		final byte[] acids=AminoAcid.numberToAcid;
		final int limit=(includeStop ? acids.length : acids.length-1);
		for(long chrom=1; chrom<=chroms; chrom++){
			bsw.print('>').print("gene").println(chrom);
			ByteBuilder bb=new ByteBuilder(wrap+1);
			byte prev='X';
			final int max=chromLength+2*pad;
			final int pad2=chromLength+pad;
			for(int i=0; i<max; ){
				for(int j=0; j<wrap && i<max; i++, j++){
					byte b;
					if(i<pad || i>=pad2){b='X';}
					else{
						b=acids[randy.nextInt(limit)];
						while(noPoly && b==prev){b=acids[randy.nextInt(limit)];}
					}
					bb.append(b);
					prev=b;
				}
				bb.nl();
				bsw.print(bb);
				bb.clear();
			}
		}
		errorState|=bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	
	/** Creates cumulative A/C/G/T probabilities for every encoded (k-1)-base prefix.
	 * Rows with no positive total use a uniform distribution. Counts are read only.
	 * @param counts Full oriented-kmer counts indexed by two-bit A/C/G/T encoding
	 * @param k Kmer length used to address counts
	 * @return New prefix-by-four cumulative probability matrix, each row ending at one
	 */
	static float[][] countsToPrefixProb(long[] counts, int k){
		final int prefixes=1<<(2*(k-1)); // 4^(k-1) possible (k-1)-mer prefixes
		float[][] matrix=new float[prefixes][4];
		
		// For each prefix (k-1)-mer
		for(int prefix=0; prefix<prefixes; prefix++){
			long[] baseCounts=new long[4];
			
			// Count occurrences of each base following this prefix
			for(int base=0; base<4; base++){
				int kmer=(prefix<<2)|base; // Append base to prefix
				baseCounts[base]=counts[kmer];
			}
			
			// Convert to cumulative probabilities
			long total=baseCounts[0]+baseCounts[1]+baseCounts[2]+baseCounts[3];
			if(total>0){
				matrix[prefix][0]=(float)baseCounts[0]/total;
				matrix[prefix][1]=matrix[prefix][0]+(float)baseCounts[1]/total;
				matrix[prefix][2]=matrix[prefix][1]+(float)baseCounts[2]/total;
				matrix[prefix][3]=1.0f; // Always 1.0 for last
			}else{
				// No data for this prefix, use uniform
				matrix[prefix][0]=0.25f;
				matrix[prefix][1]=0.50f;
				matrix[prefix][2]=0.75f;
				matrix[prefix][3]=1.00f;
			}
		}
		
		return matrix;
	}

	/** Samples one DNA base from the selected cumulative-probability row.
	 * @param prefix Encoded (k-1)-mer selecting a matrix row
	 * @param prefixMatrix Cumulative A/C/G/T probabilities; read only
	 * @param randy Generator supplying a draw in [0,1)
	 * @return A, C, G or T; T is the final cumulative bucket
	 */
	static byte nextBase(int prefix, float[][] prefixMatrix, Random randy){
		float[] probs=prefixMatrix[prefix];
		float r=randy.nextFloat();
		
		if(r<probs[0]){return (byte)'A';}
		if(r<probs[1]){return (byte)'C';}
		if(r<probs[2]){return (byte)'G';}
		return (byte)'T';
	}

	/** Expands dense canonical counts into oriented-kmer counts using the clade remap.
	 * Nonpalindromic entries copy the canonical count into each orientation; a
	 * reverse-complement palindrome receives twice its count. Input counts are unchanged.
	 * @param counts Counts indexed by CladeObject.remapMatrix[k]
	 * @param k Kmer length with an available canonical remap; constructor uses k greater than two
	 * @return New array indexed by every full two-bit kmer encoding
	 */
	static long[] unfold(long[] counts, int k){
		
		final int[] remap=CladeObject.remapMatrix[k];
		
		final int max=(1<<(2*k))-1;
		long[] unfolded=new long[max+1];
		
		for(int kmer=0; kmer<=max; kmer++){
			int rc=AminoAcid.reverseComplementBinaryFast(kmer, k);
			long count=counts[remap[kmer]];
			
			// The remap above supplies the dense canonical index shared by both orientations.
			if(kmer==rc){
				// Palindrome - double it
				unfolded[kmer]=count*2;
			}else{
				// Non-palindrome - use canonical count
				unfolded[kmer]=count;
			}
		}
		
		return unfolded;
	}
	
	/*--------------------------------------------------------------*/

	/** Optional clade or sequence input used to build the DNA prefix model. */
	private String in=null;
	/** Output destination selected by the shared parser. */
	private String out=null;
	
	/** Requested number of output records. */
	int chroms=1;
	/** Requested generated-symbol total before integer division and padding. */
	long totalLength=1000000;
	/** GC fraction for DNA generation without an input-derived model. */
	float gc=0.5f;
	/** Generated symbols per record, truncated from totalLength/chroms. */
	final int chromLength;
	/** FASTA line width captured from Shared.FASTA_WRAP during construction. */
	final int wrap;
	/** N or X symbols added at each end of every record. */
	int pad=0;
	/** Rejects repeated generated symbols in GC-only DNA and amino modes. */
	boolean noPoly=false;
	/** Includes the stop symbol in the amino sampling alphabet. */
	boolean includeStop=false;
	/** Generator seed; a negative value requests a nonfixed seed. */
	long seed=-1;
	
	/** Kmer length used for input-derived DNA probabilities; launcher documents two through five. */
	int k=5;
	/** Conditional cumulative A/C/G/T rows, or null when no input was supplied. */
	final float[][] prefixMatrix;
	
	/*--------------------------------------------------------------*/

	/** Instance sampling generator; state advances across process calls. */
	final Random randy;
	
	/** Reserved output counter; currently neither updated nor reported. */
	private long linesOut=0;
	/** Reserved output counter; currently neither updated nor reported. */
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/

	/** Optional model-input format; null when no input was supplied. */
	private final FileFormat ffin;
	/** Output format and overwrite/append policy passed to the writer. */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	
	/** Diagnostics selected by PreParser; main closes it on normal completion. */
	private PrintStream outstream=System.err;
	/** Shared verbose flag propagated to ReadWrite while parsing. */
	public static boolean verbose=false;
	/** Accumulated writer completion errors, reported by process after generation. */
	public boolean errorState=false;
	/** Permission to replace an existing output, copied from Parser. */
	private boolean overwrite=true;
	/** Whether to append to an existing output, copied from Parser. */
	private boolean append=false;
	
}
