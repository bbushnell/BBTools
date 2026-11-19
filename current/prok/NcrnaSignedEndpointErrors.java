package prok;

import java.util.Random;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import shared.Shared;

/** Complete observed signed call-minus-truth errors for one member and both ends.
 * Unlike conditioned jitter, this immutable distribution preserves every bin,
 * including zero and distant errors. Sampling uses integer mass without clamping
 * or an invented empty-bin fallback. Each sampling worker owns its random source.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSignedEndpointErrors{

	/** Native small-report CLI; no contexts, tables, training or caller is opened. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSignedEndpointErrors.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String member=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("family")){member=b;}
				else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown signed-error report parameter: "+arg);}
			}
			require(parser.in1!=null && parser.out1!=null && member!=null && parser.in2==null && parser.out2==null
				&& parser.maxReads<0 && parser.samplerate==1 && !parser.overwrite && !parser.append,
				"Use in=member.histogram.tsv out=fresh.bands.tsv family=exact_member, without limits, sampling or overwrite");
			load(parser.in1, member).writeBands(parser.out1);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Reads the existing endpoint-summary histogram for exactly one member.
	 * Absolute rows are validated as such and skipped; signed bins define draws.
	 */
	public static NcrnaSignedEndpointErrors load(final String path, final String member){
		if(path==null || member==null || member.isEmpty()){
			throw new IllegalArgumentException("Full endpoint errors require an explicit histogram and member");
		}
		final IntList[] offsets={new IntList(), new IntList()};
		final LongList[] counts={new LongList(), new LongList()};
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, HEADER, path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==5 && parser.termEquals(member, 0),
					"Full error histogram must bind exactly the selected member");
				final int end=parser.termEquals("5prime", 1) ? 0 : parser.termEquals("3prime", 1) ? 1 : -1;
				require(end>=0, "Full error histogram must name a measured gene-sense end");
				final boolean signed=parser.termEquals("signed", 2);
				require(signed || parser.termEquals("absolute", 2), "Unknown measure cannot define signed endpoint draws");
				final int offset=parser.parseInt(3);
				final long count=parser.parseLong(4);
				require(count>0 && (signed || offset>=0), "Sparse measured bins need positive counts and valid signed/absolute units");
				if(!signed){continue;}
				offsets[end].add(offset);
				counts[end].add(count);
			}
		}
		return new NcrnaSignedEndpointErrors(member,
			new Distribution(offsets[0].toArray(), counts[0].toArray()),
			new Distribution(offsets[1].toArray(), counts[1].toArray()));
	}

	/** Owns immutable distributions; the two measured ends share a locus denominator. */
	public NcrnaSignedEndpointErrors(final String member_, final Distribution five, final Distribution three){
		require(member_!=null && !member_.isEmpty() && five!=null && three!=null && five.total==three.total,
			"Both full endpoint distributions must bind one member and the same matched-locus population");
		member=member_;
		distributions=new Distribution[]{five, three};
	}

	/** Returns an immutable distribution; its private primitive arrays are never exposed. */
	public Distribution distribution(final int end){
		require(end>=0 && end<2, "Endpoint axis is zero for start and one for stop");
		return distributions[end];
	}

	/** Reports directions separately; beyond-radius counts never disappear from draws. */
	public void writeBands(final String path){
		try(Output output=new Output(path)){
			output.writer.println("member\tend\tradius_bp\tnegative_beyond\tnegative_inside\tzero\tpositive_inside\tpositive_beyond\tobserved");
			final ByteBuilder buffer=new ByteBuilder();
			for(int end=0; end<2; end++){
				final Distribution distribution=distributions[end];
				for(final int radius:RADII){
					long negativeBeyond=0, negativeInside=0, positiveInside=0, positiveBeyond=0;
					for(int bin=0; bin<distribution.size(); bin++){
						final int offset=distribution.offset(bin);
						final long count=distribution.count(bin);
						if(offset<-radius){negativeBeyond+=count;}
						else if(offset<0){negativeInside+=count;}
						else if(offset>radius){positiveBeyond+=count;}
						else if(offset>0){positiveInside+=count;}
					}
					assert(negativeBeyond+negativeInside+distribution.zero+positiveInside+positiveBeyond==distribution.total)
						: "Signed bands and zero must conserve every measured matched locus";
					output.writer.print(buffer.clear().append(member).tab().append(end==0 ? "5prime" : "3prime")
						.tab().append(radius).tab().append(negativeBeyond).tab().append(negativeInside)
						.tab().append(distribution.zero).tab().append(positiveInside).tab().append(positiveBeyond)
						.tab().append(distribution.total).nl());
				}
			}
		}
	}

	/** Exact integer sampling support, with no per-draw allocations. */
	public static final class Distribution{
		/** Copies strictly increasing offsets and positive counts, including zero if observed. */
		public Distribution(final int[] offsets_, final long[] counts_){
			require(offsets_!=null && counts_!=null && offsets_.length>0 && offsets_.length==counts_.length,
				"Measured signed support must contain paired, nonempty offset and count arrays");
			offsets=offsets_.clone();
			counts=counts_.clone();
			cumulative=new long[counts.length];
			long sum=0, atZero=0;
			for(int bin=0; bin<counts.length; bin++){
				require(counts[bin]>0 && (bin==0 || offsets[bin]>offsets[bin-1]),
					"Measured signed bins must have positive counts and unique increasing offsets");
				sum=Math.addExact(sum, counts[bin]);
				cumulative[bin]=sum;
				if(offsets[bin]==0){atZero=counts[bin];}
			}
			total=sum;
			zero=atZero;
		}

		public int size(){return offsets.length;}
		public int offset(final int bin){return offsets[bin];}
		public long count(final int bin){return counts[bin];}
		public long incorrect(){return total-zero;}

		/** Nonzero-conditioned background weight; zero is solely the positive anchor. */
		public double backgroundWeight(final int bin){
			require(incorrect()>0, "An all-correct population has no measured wrong-end background");
			return offsets[bin]==0 ? 0 : counts[bin]/(double)incorrect();
		}

		/** Unbiased bounded-long draw, then binary search on the original integer mass. */
		public int draw(final Random random){
			require(random!=null, "Full endpoint draws require a worker-owned random source");
			long bits, value;
			do{
				bits=random.nextLong()>>>1;
				value=bits%total;
			}while(bits-value+(total-1)<0);
			int low=0, high=offsets.length-1;
			while(low<high){
				final int middle=(low+high)>>>1;
				if(value<cumulative[middle]){high=middle;}
				else{low=middle+1;}
			}
			assert(value>=0 && value<cumulative[low] && (low==0 || value>=cumulative[low-1]))
				: "Each sampled integer mass must land in exactly one unchanged measured bin";
			return offsets[low];
		}

		/** Counts refer to matched loci, not original training loci or vector rows. */
		public final long total, zero;
		private final int[] offsets;
		private final long[] counts, cumulative;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	public final String member;
	private final Distribution[] distributions;
	private static final int[] RADII={0, 6, 10, 16, 20};
	private static final String HEADER="family\tend\tmeasure\toffset_bp\tloci";
}
