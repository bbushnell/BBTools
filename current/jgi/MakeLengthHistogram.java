package jgi;

import java.io.PrintStream;
import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.MetadataWriter;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Generates a length histogram, counting each mate with a nonnull base array separately.
 * Minimum, maximum, average and base totals use actual lengths. Median, mode and
 * population standard deviation use the read-count histogram's binned, capped lengths.
 * Bin labels and all length statistics are in bases; bin base totals retain actual
 * bases rather than the product of count and label. Cumulative columns include the
 * current bin and every larger bin, not exact raw-length thresholds.
 * <p>With positive width B and nonnegative cap M, main uses highest bin floor(M/B).
 * Lengths are floored by default, or rounded to the nearest bin with halfway ties
 * upward, then capped. The highest label is floor(M/B)*B; longer reads remain counted.
 * <p>Settings and counters are mutable static state and are not reset here. These
 * contracts assume one invocation in a fresh JVM, not concurrent or repeated use.
 * Numeric validation, mode ties and I/O status handling remain unresolved where
 * annotated below.
 *
 * @author Brian Bushnell
 * @date Jul 16, 2012
 */
public class MakeLengthHistogram{

	/**
	 * Parses settings, converts the requested base cap to a bin index and runs calc.
	 * On normal completion reports elapsed time, requests optional execution metadata
	 * and asks Shared to close a redirected status stream. Shared settings are not restored.
	 * @param args Input/output paths and flags; bin width and cap are integer bases
	 */
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object(){}.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Timer t=new Timer();

		String in1=null, in2=null;
		String out=null;

		Data.GENOME_BUILD=-1;
		ReadWrite.USE_UNPIGZ=true;
		Shared.capBufferLen(20);

		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out") || a.equals("hist") || a.equals("lhist")){
				out=b;
			}else if(a.equals("max") || a.equals("maxlength")){
				MAX_LENGTH=Integer.parseInt(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				NON_ZERO_ONLY=Parse.parseBoolean(b);
			}else if(a.startsWith("mult") || a.startsWith("div") || a.startsWith("bin")){
				MULT=Integer.parseInt(b);
			}else if(a.equals("round")){
				ROUND_BINS=Parse.parseBoolean(b);
			}else if(i==0 && !arg.contains("=")){
				in1=arg;
			}else if(i==1 && !arg.contains("=")){
				in2=arg;
			}else if(i==3 && !arg.contains("=")){
				//TODO: Unresolved - Bare output uses argument index 3, but not 2; positional syntax needs a separate decision.
				out=arg;
			}else{
				throw new RuntimeException("Unknown argument: "+arg);
			}
		}

		{//Process parser fields
			Parser.processQuality();
		}

		//TODO: Unresolved - Repeated main calls rescale MAX_LENGTH again, and read/base counters are not reset.
		//TODO: Unresolved - Width/cap/rounding arithmetic is not fully validated before reader startup and allocation.
		MAX_LENGTH/=MULT;

		calc(in1, in2, out);
		t.stop();
		outstream.println("Time: \t"+t);
		MetadataWriter.write(null, readsProcessed, basesProcessed, 0, 0, false);

		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}

	/**
	 * Consumes legacy concurrent input batches and emits one complete histogram report.
	 * Each nonnull base array from either mate contributes a read and its actual length.
	 * For ordinary FASTA/FASTQ, maxReads limits fragments (single reads or pairs), while
	 * report counts include individual mates; the generic reader rejects a zero quota.
	 * <p>Direct callers must prepare MAX_LENGTH as the highest bin INDEX, with positive
	 * MULT and nonnegative, representable allocation/rounding arithmetic. This method
	 * does not perform main's bases-to-bin conversion and does not reset static counters.
	 * <p>Returns spent batches, requests input closure and joins the output writer. Throws
	 * when writer finalization reports an error. The retained input-close result is still
	 * ignored, and no cleanup-on-exception or complete-output guarantee is made.
	 * <p>Rows contain label, read count/percent/upper-tail count/percent, then actual
	 * base count/percent/upper-tail bases/percent. Upper tails include the current bin.
	 * Percentages have three decimal places; average and population deviation have one.
	 * Empty bins may be suppressed; scanning stops at the first label at least the raw
	 * maximum, or at the cap. With flooring, the final displayed bin can therefore be empty.
	 * @param in1 Required primary input path or standard-input alias
	 * @param in2 Optional second mate path, or null
	 * @param out Histogram path, or null for standard output
	 */
	public static void calc(String in1, String in2, String out){

		FastaReadInputStream.MIN_READ_LEN=1;

		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
//			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		if(verbose){outstream.println("Paired: "+paired);}

		final int max=MAX_LENGTH;
		long[] readHist=new long[max+1];
		long[] baseHist=new long[max+1];

		int maxFound=0;
		int minFound=Integer.MAX_VALUE;

		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(paired==(r.mate!=null));
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//outstream.println("reads.size()="+reads.size());
				for(Read r1 : reads){
					final Read r2=r1.mate;

//					System.out.println("Processing read "+r.numericID);

					if(r1!=null && r1.bases!=null){
						readsProcessed++;
						final int x=r1.length();
						basesProcessed+=x;
						final int y=Tools.min(max, ((ROUND_BINS ? x+MULT/2 : x))/MULT);
						readHist[y]++;
						baseHist[y]+=x;
						maxFound=Tools.max(maxFound, x);
						minFound=Tools.min(minFound, x);
					}

					if(r2!=null && r2.bases!=null){
						readsProcessed++;
						final int x=r2.length();
						basesProcessed+=x;
						final int y=Tools.min(max, ((ROUND_BINS ? x+MULT/2 : x))/MULT);
						readHist[y]++;
						baseHist[y]+=x;
						maxFound=Tools.max(maxFound, x);
						minFound=Tools.min(minFound, x);
					}

				}
				//outstream.println("returning list");
				cris.returnList(ln);
				//outstream.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(verbose){outstream.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){outstream.println("Returned list");}
			//TODO: Probable bug - Ignoring the native close error can publish a histogram after a failed input read.
			ReadWrite.closeStream(cris);
			if(verbose){outstream.println("Closed stream");}
			outstream.println("Processed "+readsProcessed+" reads.");
		}

		//TODO: Unresolved - Empty-population average/percentage divisors can be zero; no finite-value policy is selected.
		if(readsProcessed<1){minFound=0;}
		double stdev=Tools.standardDeviationHistogram(readHist)*MULT;
		int median=medianBin(readHist)*MULT;
		//TODO: Unresolved - Tied modes still use the shared percentile rank, which can differ from medianBin.
		int mode=Tools.calcModeHistogram(readHist)*MULT;

		double[] readHistF=new double[max+1];
		long[] readHistC=new long[max+1];
		double[] readHistCF=new double[max+1];

		double[] baseHistF=new double[max+1];
		long[] baseHistC=new long[max+1];
		double[] baseHistCF=new double[max+1];

		readHistC[max]=readHist[max];
		baseHistC[max]=baseHist[max];
		for(int i=max; i>0; i--){
			readHistC[i-1]=readHistC[i]+readHist[i-1];
			baseHistC[i-1]=baseHistC[i]+baseHist[i-1];
		}
		for(int i=0; i<=max; i++){
			readHistCF[i]=readHistC[i]*100d/readHistC[0];
			readHistF[i]=readHist[i]*100d/readHistC[0];
			baseHistCF[i]=baseHistC[i]*100d/baseHistC[0];
			baseHistF[i]=baseHist[i]*100d/baseHistC[0];
		}

		TextStreamWriter tsw=new TextStreamWriter(out==null ? "stdout" : out, overwrite, append, false);
		tsw.start();
		tsw.println("#Reads:\t"+readsProcessed);
		tsw.println("#Bases:\t"+baseHistC[0]);
		tsw.println("#Max:\t"+maxFound);
		tsw.println("#Min:\t"+minFound);
		tsw.println("#Avg:\t"+Tools.format("%.1f", (baseHistC[0]*1d/readsProcessed)));
		tsw.println("#Median:\t"+median);
		tsw.println("#Mode:\t"+mode);
		tsw.println("#Std_Dev:\t"+Tools.format("%.1f", stdev));
		tsw.println("#Read Length Histogram:");
		tsw.println("#Length\treads\tpct_reads\tcum_reads\tcum_pct_reads\tbases\tpct_bases\tcum_bases\tcum_pct_bases");
		for(int i=0; i<=max; i++){
			if(readHist[i]>0 || !NON_ZERO_ONLY){
				tsw.println((i*MULT)+"\t"+readHist[i]+Tools.format("\t%.3f%%", readHistF[i])+"\t"+readHistC[i]+Tools.format("\t%.3f%%", readHistCF[i])+
						"\t"+baseHist[i]+Tools.format("\t%.3f%%", baseHistF[i])+"\t"+baseHistC[i]+Tools.format("\t%.3f%%", baseHistCF[i]));
			}
			if(i*MULT>=maxFound){break;}
		}
		if(tsw.poisonAndWait()){
			throw new RuntimeException("MakeLengthHistogram terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Select the exact middle bin for odd totals and the lower middle bin for
	 * even totals. Empty populations return zero; totals beyond the signed long
	 * range retain the shared helper's behavior.
	 * @param histogram Nonnegative read counts per length bin
	 * @return Selected bin index, before multiplication by the bin width
	 */
	private static int medianBin(final long[] histogram){
		if(histogram.length<1){return 0;}
		long total=0;
		for(final long count : histogram){
			assert(count>=0) : "Median cumulative ranks require nonnegative read-bin counts: "+count;
			if(count>Long.MAX_VALUE-total){return Tools.percentileHistogram(histogram, 0.5);}
			total+=count;
		}
		final long target=total/2+(total&1L);
		long cumulative=0;
		for(int i=0; i<histogram.length; i++){
			cumulative+=histogram[i];
			if(cumulative>=target){return i;}
		}
		return histogram.length-1;
	}

	/** Input fragment quota; positive limits reads/pairs, negative is unlimited; zero is unsupported. */
	public static long maxReads=-1;
	/** Individual nonnull-base reads accumulated by calc; not reset between calls. */
	public static long readsProcessed=0;
	/** Actual-base accumulator used for metadata; report bases instead come from the current histogram. */
	public static long basesProcessed=0;
	/** Requested CLI cap in bases before main scales it; highest bin index during calc. */
	public static int MAX_LENGTH=80000;
	/** Positive bin width in bases, also used to scale bin-valued statistics. */
	public static int MULT=10;
	/** False floors lengths; true adds integer half-width before division to select the nearest bin. */
	public static boolean ROUND_BINS=false;
	/** Suppresses zero-count output rows without changing accumulation or statistics. */
	public static boolean NON_ZERO_ONLY=true;

	/** Status and elapsed-time destination; PreParser may redirect it from standard error. */
	static PrintStream outstream=System.err;

	/** Appends a complete report in supported formats; requires overwrite=false. */
	public static boolean append=false;
	/** Allows replacement of an existing histogram file; mutually exclusive with append. */
	public static boolean overwrite=true;
	/** Programmatic extra reader-status messages; there is no dedicated verbose CLI branch here. */
	public static boolean verbose=false;

}
