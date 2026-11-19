package test.prok;

import java.nio.charset.StandardCharsets;

import prok.NcrnaSignedEndpointErrors;
import prok.NcrnaSignedEndpointErrors.Distribution;
import prok.NcrnaPositionalLeaveOneOut;
import prok.NcrnaWeightedEndpointTable;
import prok.RrnaPositionalKmerTable;
import prok.RrnaPositionalKmerTable.Table;

/** Independent per-offset count oracle for signed, nonuniform backgrounds.
 * A repeated word occurs at offsets with different valid denominators; one
 * locus also has unavailable offsets. No fitted network or biological input.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWeightedEndpointTableTest{

	public static void runFixtures(){
		final byte[][] sequences={
			bases("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
			bases("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
			bases("AAAAAAAANAAAAAACGTACGTACGTACGTACGTAAAAAAA"),
			bases("AAAAAAAAAAAAAA")};
		for(final String end:new String[]{"5prime", "3prime"}){
			check(sequences, end, 5);
			check(sequences, end, 9);
		}
	}

	private static void check(final byte[][] sequences, final String end, final int k){
		final Distribution distribution=new Distribution(ERROR_OFFSETS, ERROR_COUNTS);
		final NcrnaSignedEndpointErrors errors=new NcrnaSignedEndpointErrors("fixture", distribution, distribution);
		final Table geometry=new Table("fixture", end, k, 0, .5, 3);
		final NcrnaWeightedEndpointTable table=new NcrnaWeightedEndpointTable(geometry, errors);
		for(final byte[] sequence:sequences){
			geometry.add(sequence, ENDPOINT);
			table.add(sequence, ENDPOINT);
		}
		geometry.finish();
		table.finish();
		final NcrnaWeightedEndpointTable.View view=table.leaveOneOut();
		final NcrnaPositionalLeaveOneOut windows=new NcrnaPositionalLeaveOneOut(geometry, table);
		for(int omitted=-1; omitted<sequences.length; omitted++){
			if(omitted>=0){
				view.set(sequences[omitted], ENDPOINT);
				windows.set(sequences[omitted], ENDPOINT);
			}
			final int limit=k==5 ? geometry.space : 5;
			for(int index=0; index<limit; index++){
				final int word=k==5 ? index : index==0 ? 0 : index==1 ? geometry.space/3
					: index==2 ? geometry.space-1 : index==3 ? 12345 : 54321;
				final double expected=literal(sequences, omitted, word, k);
				final double actual=omitted<0 ? table.scoreKey(word) : view.scoreKey(word);
				require(Double.doubleToLongBits(actual)==Double.doubleToLongBits(expected),
					"Weighted full/LOO scores must exactly match independently counted per-offset frequencies: "+omitted+" key="+word);
			}
			checkWindows(sequences, omitted, k, table, omitted<0 ? null : windows);
		}
		// A failed binding must not leave the preceding locus usable.
		try{view.set(null, ENDPOINT);}
		catch(IllegalArgumentException expected){/* The next lookup must fail. */}
		boolean rejected=false;
		try{view.scoreKey(0);}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Weighted LOO must invalidate a preceding locus before a failed bind");
	}

	/** Both cached and distant paths use the independently counted weighted view. */
	private static void checkWindows(final byte[][] sequences, final int omitted, final int k,
			final NcrnaWeightedEndpointTable table, final NcrnaPositionalLeaveOneOut view){
		final byte[] bases=sequences[Math.max(0, omitted)];
		final float[] output=new float[7];
		final float[] expected=new float[7];
		for(final int delta:new int[]{-17, -4, 0, 4, 17}){
			boolean complete=true;
			for(int site=0; site<output.length; site++){
				final int position=ENDPOINT+delta+site-3;
				final int word=position<0 || position+k>bases.length ? -1 : RrnaPositionalKmerTable.encode(bases, position, k);
				if(word<0){complete=false;}
				else{expected[site]=(float)literal(sequences, omitted, word, k);}
			}
			final boolean available=view==null ? table.scoreWindow(bases, ENDPOINT+delta, output, 0)
				: view.scoreDistantWindow(ENDPOINT+delta, output, 0);
			require(available==complete, "Explicit weighted window availability must match literal clipped/ambiguous word positions");
			if(available){
				for(int site=0; site<output.length; site++){
					require(Float.floatToIntBits(output[site])==Float.floatToIntBits(expected[site]),
						"Every available weighted word must match the independent full/omitted-locus count oracle");
				}
			}
		}
	}

	/** Literal scan recomputes each queried word's counts without using table rows. */
	private static double literal(final byte[][] sequences, final int omitted, final int word, final int k){
		final long[] valid=new long[4], counts=new long[4];
		for(int locus=0; locus<sequences.length; locus++){
			if(locus==omitted){continue;}
			final byte[] bases=sequences[locus];
			for(int bin=0; bin<ERROR_OFFSETS.length; bin++){
				final int position=ENDPOINT+ERROR_OFFSETS[bin];
				if(position<0 || position+k>bases.length){continue;}
				final int observed=RrnaPositionalKmerTable.encode(bases, position, k);
				if(observed<0){continue;}
				valid[bin]++;
				if(observed==word){counts[bin]++;}
			}
		}
		final double prior=.5*(1<<(2*k));
		final double anchor=(counts[1]+.5)/(valid[1]+prior);
		double background=0;
		for(int bin=0; bin<ERROR_OFFSETS.length; bin++){
			if(bin==1){continue;}
			background+=(ERROR_COUNTS[bin]/6.0)*((counts[bin]+.5)/(valid[bin]+prior));
		}
		return Math.log(anchor/background);
	}

	private static byte[] bases(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
	private static final int ENDPOINT=10;
	private static final int[] ERROR_OFFSETS={-7, 0, 2, 8};
	private static final long[] ERROR_COUNTS={2, 4, 1, 3};
}
