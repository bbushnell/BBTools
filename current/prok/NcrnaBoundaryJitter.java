package prok;

import java.util.Random;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Observed member-level endpoint errors conditioned on a small radius.
 * Consumes NcrnaEndpointSummary's signed histograms without calling or matching.
 * Out-of-radius mass is reported, never moved to edge bins. The retained native
 * distributions are immutable after loading; each sampling worker owns its RNG.
 * An empty window remains explicitly unavailable rather than inventing jitter.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryJitter{

	/** Reads both complete signed error distributions for one emitted member.
	 * Absolute-error rows are validated as the named alternate measure and ignored.
	 * @param path Java-produced family/end/measure/offset_bp/loci histogram
	 * @param member Exact emitted-member label in every input row
	 * @param radius Positive conditioning radius in bases
	 */
	public static NcrnaBoundaryJitter load(final String path, final String member, final int radius){
		if(path==null || member==null || member.isEmpty() || radius<1){
			throw new IllegalArgumentException("Conditioned jitter requires a histogram, member and positive radius");
		}
		final NcrnaBoundaryJitter result=new NcrnaBoundaryJitter(member, radius);
		final LineParser1 parser=new LineParser1('\t');
		final boolean[] seen=new boolean[2];
		final int[] last=new int[2];
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, HEADER, path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==5 && parser.termEquals(member, 0), "Jitter histogram must bind exactly the selected emitted member");
				final int end=parser.termEquals("5prime", 1) ? 0 : parser.termEquals("3prime", 1) ? 1 : -1;
				require(end>=0, "Jitter histogram must name a measured gene-sense end");
				final boolean signed=parser.termEquals("signed", 2);
				require(signed || parser.termEquals("absolute", 2), "Unknown histogram measure cannot define signed jitter");
				final int error=parser.parseInt(3);
				final long count=parser.parseLong(4);
				require(count>0 && (signed || error>=0), "Sparse histogram bins require positive counts and valid signed/absolute units");
				if(!signed){continue;}
				require(!seen[end] || error>last[end], "Signed bins must be unique and strictly increasing within each end");
				seen[end]=true;
				last[end]=error;
				result.observed[end]=Math.addExact(result.observed[end], count);
				if(error<-radius){result.low[end]=Math.addExact(result.low[end], count);}
				else if(error>radius){result.high[end]=Math.addExact(result.high[end], count);}
				else{result.distributions[end].add(error, count);}
			}
		}
		require(result.observed[0]==result.observed[1], "Both measured ends must retain the same matched-locus denominator");
		for(int end=0; end<2; end++){
			require(result.observed[end]==Math.addExact(result.retained(end), Math.addExact(result.low[end], result.high[end])),
				"Every observed error must remain inside the conditioning window or in exactly one excluded tail");
		}
		return result;
	}

	private NcrnaBoundaryJitter(final String member_, final int radius_){
		member=member_;
		radius=radius_;
		distributions=new RrnaEndpointVector.Distribution[]{
			new RrnaEndpointVector.Distribution(radius), new RrnaEndpointVector.Distribution(radius)};
	}

	/** No fallback is defined for an end whose conditioning window has zero mass. */
	public boolean available(final int end){return retained(end)>0;}
	public long retained(final int end){return distributions[end].total;}
	public long observed(final int end){return observed[end];}
	public long excludedLow(final int end){return low[end];}
	public long excludedHigh(final int end){return high[end];}

	/** Draws an actual observed signed call-minus-reference error, in bases. */
	public int draw(final int end, final Random random){
		require(end>=0 && end<2 && random!=null && available(end), "Cannot draw from an absent conditioned endpoint distribution");
		return distributions[end].draw(random);
	}

	/** Fresh small accounting table; this is written even when a window is empty. */
	public void write(final String path){
		try(Output output=new Output(path)){
			output.writer.println("member\tend\tradius_bp\tobserved\tretained\texcluded_low\texcluded_high\tstate");
			final ByteBuilder buffer=new ByteBuilder();
			for(int end=0; end<2; end++){
				output.writer.print(buffer.clear().append(member).tab().append(end==0 ? "5prime" : "3prime")
					.tab().append(radius).tab().append(observed[end]).tab().append(retained(end))
					.tab().append(low[end]).tab().append(high[end])
					.tab().append(available(end) ? "CONDITIONED_OBSERVED_BINS" : "EMPTY_CONDITIONING_WINDOW").nl());
			}
		}
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	public final String member;
	public final int radius;
	private final RrnaEndpointVector.Distribution[] distributions;
	private final long[] observed=new long[2], low=new long[2], high=new long[2];
	private static final String HEADER="family\tend\tmeasure\toffset_bp\tloci";
}
