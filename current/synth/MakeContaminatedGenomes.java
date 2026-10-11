package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.regex.Matcher;
import shared.Random;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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
 * Builds synthetic contaminated genomes from pairs of entries in an input file list.
 * Each output combines the retained records from two loaded genomes, with the larger
 * sampled fraction first. Retention is applied separately to each read as a circular
 * fragment; optional substitutions and indels follow retention. Read IDs remain intact
 * and qualities are discarded. The generated filename records retained base counts,
 * sampled fractions, source filename cores and the output index.
 *
 * Selection and mutation run on the calling thread; each output uses a writer thread.
 * Instances retain mutable counters and error state, and construction changes shared
 * I/O settings. Concurrent use of one instance is unsupported. Writer completion errors
 * are accumulated and reported after the requested generation loop finishes.
 *
 * @author Brian Bushnell
 * @date Oct 7, 2014
 */
public class MakeContaminatedGenomes {

	/** Runs generation and closes a redirected diagnostic stream on normal completion.
	 * @param args Command-line arguments accepted by the constructor
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		MakeContaminatedGenomes x=new MakeContaminatedGenomes(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parses the input list, output pattern, sampling settings and mutation rates.
	 * Also configures shared input, compression and diagnostic settings. The legacy
	 * id/ani/identity aliases are currently parsed directly as combined error rates.
	 * Output patterns and numerical ranges are not comprehensively validated here.
	 * @param args Tool options and shared BBTools arguments, including config-file options
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
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(fofn==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(outPattern!=null && outPattern.equalsIgnoreCase("null")){outPattern=null;}
		
		fffofn=FileFormat.testInput(fofn, FileFormat.TXT, null, true, true);
	}
	
	/** Loads input entries and generates the requested outputs with one seeded generator.
	 * Positive output counts require at least two entries before selection or output creation;
	 * nonpositive counts skip generation. Duplicate path strings remain distinct entries.
	 * Reports counts of input records and pre-mutation bases, then throws for writer errors.
	 * Counters and prior error state are retained if this instance is invoked again.
	 * @param t Running timer to stop after generation
	 */
	void process(Timer t){
		final String[] in=TextFile.toStringLines(fffofn);
		if(chimeras>0 && in.length<2){
			throw new IllegalArgumentException("At least two input entries are required for chimeras="+
					chimeras+"; found "+in.length);
		}
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
	
	/** Selects two distinct entry indices, loads both genomes and writes their modified reads.
	 * Samples each fraction separately and raises it to the respective configured exponent.
	 * @param in At least two genome-path entries; duplicate path values are permitted
	 * @param randy Generator shared by selection, retention and mutation
	 * @param cid Zero-based output index used in the generated filename
	 * @return Generated output path
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
	
	/** Modifies both read lists in place, derives an output name and writes every retained record.
	 * Larger sampled fraction goes first; ties place A first. Filename base counts are measured
	 * after retention and mutation. Waits for the writer and accumulates its completion error.
	 * Empty records are written too.
	 * @param inA First source path, used for its filename core
	 * @param inB Second source path, used for its filename core
	 * @param readsA First genome's reads; bases and qualities may be replaced
	 * @param readsB Second genome's reads; bases and qualities may be replaced
	 * @param fracA Fraction applied to every read from A before mutation
	 * @param fracB Fraction applied to every read from B before mutation
	 * @param randy Generator used for fragment starts and mutations
	 * @param cid Output index embedded in the filename
	 * @return Output path, including when writer completion sets errorState
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
			out=outPattern.replaceFirst(regex, Matcher.quoteReplacement(delimiter+sizeA+delimiter+Tools.format("%.3f", fracA)+delimiter+ReadWrite.stripToCore(inA)+delimiter+sizeB+delimiter+Tools.format("%.3f", fracB)+delimiter+ReadWrite.stripToCore(inB)+delimiter+cid+delimiter));
			ByteStreamWriter bsw=new ByteStreamWriter(out, overwrite, append, true);
			bsw.start();
//			System.err.println("G");
			for(Read r : readsA){bsw.println(r);}
			for(Read r : readsB){bsw.println(r);}
//			System.err.println("H");
			errorState|=bsw.poisonAndWait();
		}else{
//			System.err.println("I");
			out=outPattern.replaceFirst(regex, Matcher.quoteReplacement(delimiter+sizeB+delimiter+Tools.format("%.3f", fracB)+delimiter+ReadWrite.stripToCore(inB)+delimiter+sizeA+delimiter+Tools.format("%.3f", fracA)+delimiter+ReadWrite.stripToCore(inA)+delimiter+cid+delimiter));
			ByteStreamWriter bsw=new ByteStreamWriter(out, overwrite, append, true);
			bsw.start();
//			System.err.println("J");
			for(Read r : readsB){bsw.println(r);}
			for(Read r : readsA){bsw.println(r);}
//			System.err.println("K");
			errorState|=bsw.poisonAndWait();
		}
//		System.err.println("L");
		return out;
	}
	
	/** Retains a circular fragment, optionally mutates it, and counts the resulting bases.
	 * Clears qualities even when bases are unchanged. Fractions below one truncate the retained
	 * length to an integer and choose a start within the original read; copying can wrap to its
	 * beginning. Mutations affect fully defined bases: substitutions choose another base,
	 * deletions consume the position, and insertions emit a base and retry that position.
	 * @param r Read with nonnull bases; modified in place while preserving its ID
	 * @param bb Reusable scratch builder; cleared on entry and its final contents are unspecified
	 * @param genomeFraction Fraction retained before mutation; values at least one keep all bases
	 * @param randy Generator for fragment placement and mutation choices
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
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt3()+1)&3)];
						bb.append(b);
					}else if(randy.nextBoolean()){//del
						//do nothing
					}else{//ins
						i--;
						b=AminoAcid.numberToBase[randy.nextInt()&3];
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
	
	/** File containing one genome input path per line. */
	private String fofn=null;

	/** Output path whose first regex match is replaced with generated source/count information. */
	private String outPattern=null;
	/** Optional output-path list; no constructor option currently sets this field. */
	private String outNames=null;
	
	/** Requested number of mixed-genome output files. */
	private int chimeras=1;
	/** Generator seed; a negative value requests a nonfixed seed. */
	private long seed=-1;
	/** Exponent applied to the first genome's uniform fraction draw. */
	double exponent1=1;
	/** Exponent applied to the second genome's uniform fraction draw. */
	double exponent2=1;
	/** Separator inserted between generated filename components. */
	String delimiter="_";
	/** Regular expression whose first match in outPattern is replaced. */
	String regex="#";
	
	/** Substitution threshold within the combined per-position mutation draw. */
	double subRate=0;
	/** Indel contribution to errorRate; an indel draw chooses insertion or deletion. */
	double indelRate=0;
	/** Combined mutation threshold, initialized from subRate plus indelRate after parsing. */
	double errorRate=0;
	/** Bases retained after both fragment selection and mutation; not the displayed input total. */
	long basesRetained=0;

	/** Source records processed across all outputs, counting repeated input selections again. */
	long readsProcessed=0;
	/** Source bases counted before each read is modified. */
	long basesProcessed=0;
	
	/** Diagnostic stream selected by PreParser; closed by main only on normal completion. */
	private PrintStream outstream=System.err;
	
	/*--------------------------------------------------------------*/
	
	/** Input format for the text file listing genome paths. */
	private final FileFormat fffofn;
	
	/*--------------------------------------------------------------*/
	
	/** Shared verbose setting propagated to the selected I/O helpers during parsing. */
	public static boolean verbose=false;
	/** Accumulated writer completion failures; process reports these after generation. */
	public boolean errorState=false;
	/** Permission to replace an existing output, copied from Parser. */
	private boolean overwrite=true;
	/** Whether to append output records, copied from Parser. */
	private boolean append=false;
	
}
