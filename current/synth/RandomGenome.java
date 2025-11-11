package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import bin.AdjustEntropy;
import clade.Clade;
import clade.CladeLoader;
import clade.CladeObject;
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
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @contributor Isla
 * @date Jan 3, 2013
 *
 */
public class RandomGenome {
	
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
		assert(2*pad+totalLength/chroms<Shared.MAX_ARRAY_LEN) : "Length per chrom must be at most "+Shared.MAX_ARRAY_LEN;
		chromLength=(int)(totalLength/chroms);

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffin=FileFormat.testInput(in, FileFormat.CLADE, null, true, false, false);
		ffout=FileFormat.testOutput(out, FileFormat.FA, null, true, overwrite, append, false);
		
		if(ffin!=null) {
			Clade clade;
			AdjustEntropy.load();
			if(ffin.clade()) {
				ArrayList<Clade> clades=CladeLoader.loadCladesFromClade(ffin);
				clade=clades.get(0);
			}else {
				clade=CladeLoader.loadCladeFromSequence(ffin);
			}
			long[] counts=clade.counts[k];
			if(k>2) {counts=unfold(counts, k);}
			prefixMatrix=countsToPrefixProb(counts, k);
		}else {
			prefixMatrix=null;
		}
		
		randy=Shared.threadLocalRandom(seed);
	}
	
	void process(Timer t){
		if(Shared.AMINO_IN){
			processAmino(t);
		}else{
			processNucleotide(t);
		}
	}
	
	void processNucleotide(Timer t){
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		
		for(int chrom=1; chrom<=chroms; chrom++){
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
							b=AminoAcid.numberToBase[randy.nextInt(4)];
							while(noPoly && b==prev){b=AminoAcid.numberToBase[randy.nextInt(4)];}
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
		bsw.poisonAndWait();
	}
	
	void processAmino(Timer t){
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		final byte[] acids=AminoAcid.numberToAcid;
		final int limit=(includeStop ? acids.length : acids.length-1);
		for(int chrom=1; chrom<=chroms; chrom++){
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
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Convert kmer counts to prefix probability matrix.
	 * @param counts Array of kmer counts (may be folded/canonical form)
	 * @param k Kmer length
	 * @return float[4^(k-1)][4] where [prefix][base] = cumulative probability
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

	/**
	 * Choose next base based on prefix probabilities.
	 * @param prefix The (k-1)-mer prefix as binary encoding
	 * @param prefixMatrix Cumulative probability matrix
	 * @param randy Random number generator
	 * @return Next base (A/C/G/T)
	 */
	static byte nextBase(int prefix, float[][] prefixMatrix, Random randy){
		float[] probs=prefixMatrix[prefix];
		float r=randy.nextFloat();
		
		if(r<probs[0]){return (byte)'A';}
		if(r<probs[1]){return (byte)'C';}
		if(r<probs[2]){return (byte)'G';}
		return (byte)'T';
	}

	/**
	 * Unfold canonical kmer counts to forward orientation.
	 * Palindromes get doubled, non-palindromes get their count plus their RC's count.
	 * @param counts Canonical (folded) counts
	 * @param k Kmer length
	 * @return Unfolded counts array
	 */
	static long[] unfold(long[] counts, int k){
		
		final int[] remap=CladeObject.remapMatrix[k];
		
		final int max=(1<<(2*k))-1;
		long[] unfolded=new long[max+1];
		
		for(int kmer=0; kmer<=max; kmer++){
			int rc=AminoAcid.reverseComplementBinaryFast(kmer, k);
			long count=counts[remap[kmer]];
			
			// Find canonical index (this is simplified - you'd use the actual remap)
			// For now assuming counts[canon] exists
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

	private String in=null;
	private String out=null;
	
	int chroms=1;
	long totalLength=1000000;
	float gc=0.5f;
	final int chromLength;
	final int wrap;
	int pad=0;
	boolean noPoly=false;
	boolean includeStop=false;
	long seed=-1;
	
	int k=5;
	final float[][] prefixMatrix;
	
	/*--------------------------------------------------------------*/

	final Random randy;
	
	private long linesOut=0;
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/

	private final FileFormat ffin;
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
