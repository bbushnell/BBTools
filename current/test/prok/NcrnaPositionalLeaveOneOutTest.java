package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.PreParser;
import prok.NcrnaPositionalLeaveOneOut;
import prok.RrnaPositionalKmerTable.Table;
import shared.Shared;

/** Literal omitted-locus rebuild oracle for immutable positional training views.
 * Includes repeated words, unequal valid denominators, ambiguous and clipped
 * windows, both endpoint labels and all three probe placements.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalLeaveOneOutTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPositionalLeaveOneOutTest.class, false);
		try{
			if(pp.args.length!=0){throw new IllegalArgumentException("This finite oracle accepts no parameters: "+pp.args[0]);}
			runFixtures();
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void runFixtures(){
		final byte[][] sequences={
			bases("AAAAAAAAACGTACGTACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
			bases("AAAAAAAAACGTACGTACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
			bases("AAAAAAAANCGTACGTACGTNCGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
			bases("AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG")};
		for(final String end:new String[]{"5prime", "3prime"}){
			for(final int offset:new int[]{-9, -3, 1}){
				for(final int endpoint:new int[]{2, 21, 54}){checkRebuild(sequences, end, offset, endpoint, 5, 3);}
			}
		}
		checkRebuild(sequences, "3prime", -9, 21, 9, 10);
		checkEmptyRemainder();
		checkBindingFailure();
		NcrnaWeightedEndpointTableTest.runFixtures();
		System.out.println("PASS NcrnaPositionalLeaveOneOutTest literal_rebuild all_keys both_ends three_placements repeated_words clipped_ambiguous denominators cached_windows distant_windows immutable_parent reuse failed_rebinding empty_remainder");
	}

	/** Compare every packed key, not only keys occurring in the held-out sequence. */
	private static void checkRebuild(final byte[][] sequences, final String end, final int offset,
			final int endpoint, final int k, final int radius){
		final double alpha=k==9 ? 1.0/(1<<18) : .5;
		final Table full=new Table("fixture", end, k, offset, alpha, radius);
		// A clean interior background keeps every remaining offset observed even
		// when the target locus is clipped or ambiguous.
		final byte[] background=sequences[3];
		full.add(background, 25);
		full.add(background, 25);
		for(final byte[] sequence:sequences){full.add(sequence, endpoint);}
		full.finish();
		final long[][] before=new long[full.sites][];
		for(int site=0; site<before.length; site++){before[site]=full.counts[site].clone();}
		final long[] totals=full.valid.clone();
		final NcrnaPositionalLeaveOneOut view=new NcrnaPositionalLeaveOneOut(full);
		for(int omitted=0; omitted<sequences.length; omitted++){
			final Table rebuilt=new Table("fixture", end, k, offset, alpha, radius);
			rebuilt.add(background, 25);
			rebuilt.add(background, 25);
			for(int locus=0; locus<sequences.length; locus++){
				if(locus!=omitted){rebuilt.add(sequences[locus], endpoint);}
			}
			rebuilt.finish();
			view.set(sequences[omitted], endpoint);
			for(int word=0; word<full.space; word++){
				check(Double.doubleToLongBits(view.scoreKey(word))==Double.doubleToLongBits(rebuilt.scoreKey(word)),
					"Leave-one-out must exactly match literal omitted-locus normalization for every key");
			}
			final float[] expected=new float[full.sites], actual=new float[full.sites+2];
			for(int delta=-radius; delta<=radius; delta++){
				rebuilt.scoreWindow(sequences[omitted], endpoint+delta+offset, expected);
				Arrays.fill(actual, Float.POSITIVE_INFINITY);
				final boolean available=view.scoreWindow(endpoint+delta, actual, 1);
				checkWindow(expected, actual, available);
				if(!available){
					for(final float value:actual){check(value==Float.POSITIVE_INFINITY, "An unavailable cached window must not copy unused scores");}
				}
				checkWindow(expected, actual, view.scoreDistantWindow(endpoint+delta, actual, 1));
			}
			for(final int delta:new int[]{-sequences[omitted].length, -2*radius, -radius-1,
					radius+1, 2*radius, sequences[omitted].length}){
				rebuilt.scoreWindow(sequences[omitted], endpoint+delta+offset, expected);
				Arrays.fill(actual, Float.POSITIVE_INFINITY);
				checkWindow(expected, actual, view.scoreDistantWindow(endpoint+delta, actual, 1));
			}
			boolean rejected=false;
			try{view.scoreWindow(endpoint+radius+1, actual, 1);}
			catch(IllegalArgumentException expectedFailure){rejected=true;}
			check(rejected, "The legacy cached-window API must still reject distant candidates");
		}
		check(Arrays.equals(totals, full.valid), "Worker exclusions must never mutate shared valid denominators");
		for(int site=0; site<before.length; site++){
			check(Arrays.equals(before[site], full.counts[site]), "Worker exclusions must never mutate shared positional counts");
		}
	}

	/** The legacy rebuilt table supplies independent missing-word and exact-score oracles. */
	private static void checkWindow(final float[] expected, final float[] actual, final boolean available){
		boolean complete=true;
		for(final float value:expected){complete&=Float.isFinite(value);}
		check(available==complete, "Explicit availability must match the literal table's clipped/ambiguous word positions");
		if(available){
			for(int site=0; site<expected.length; site++){
				check(Float.floatToIntBits(expected[site])==Float.floatToIntBits(actual[site+1]),
					"Available windows must exactly equal the literal omitted-locus table");
			}
		}
		check(actual[0]==Float.POSITIVE_INFINITY && actual[actual.length-1]==Float.POSITIVE_INFINITY,
			"A distant probe must not overwrite adjacent feature slots");
	}

	/** An invalid next binding cannot leave the preceding locus usable by accident. */
	private static void checkBindingFailure(){
		final Table table=new Table("fixture", "5prime", 5, 0, .5, 3);
		final byte[] clean=bases("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		table.add(clean, 15);
		table.add(clean, 15);
		table.finish();
		final NcrnaPositionalLeaveOneOut view=new NcrnaPositionalLeaveOneOut(table);
		final float[] output=new float[table.sites];
		boolean rejected=false;
		try{view.scoreDistantWindow(25, output, 0);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A distant lookup requires a successfully bound counted locus");
		view.set(clean, 15);
		view.scoreDistantWindow(25, output, 0);
		try{view.set(null, 15);}
		catch(IllegalArgumentException expected){/* The old binding must already be invalid. */}
		rejected=false;
		try{view.scoreDistantWindow(25, output, 0);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Failed rebinding must release the old bases and invalidate distant lookup");
		view.set(clean, 15);
		view.scoreDistantWindow(25, output, 0);
	}

	/** Literal table finalization rejects a prior-only offset; the view must agree. */
	private static void checkEmptyRemainder(){
		final Table table=new Table("fixture", "5prime", 5, 0, .5, 3);
		final byte[] clean=bases("AAAAAAAAAAAAAAAAAAAA"), missing=clean.clone();
		Arrays.fill(missing, (byte)'N');
		table.add(clean, 8);
		table.add(missing, 8);
		table.finish();
		boolean rejected=false;
		try{new NcrnaPositionalLeaveOneOut(table).set(clean, 8);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An offset with no remaining observation must fail before emitting pseudocount-only features");
	}

	private static byte[] bases(final String sequence){return sequence.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean condition, final String why){
		if(!condition){throw new AssertionError(why);}
	}
}
