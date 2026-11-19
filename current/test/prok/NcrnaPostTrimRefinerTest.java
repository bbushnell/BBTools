package test.prok;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import idaligner.EndClippedAligner;
import prok.NcrnaBoundaryRefiner;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;

/** Finite oracle for refiner timing, eligibility and retained caller state.
 * The trim spy supplies literal independent coordinates; no alignment outcome is
 * inferred from the new branch. Existing real trim behavior is tested elsewhere.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPostTrimRefinerTest{

	public static void main(final String[] args) throws Exception{
		check(args.length==0, "This finite oracle takes no data or configuration overrides");
		final byte[] bases=new byte[160];
		Arrays.fill(bases, (byte)'A');
		final byte[] original=bases.clone();
		final NcrnaBoundaryRefiner defaultMode=new NcrnaBoundaryRefiner(){
			@Override
			public boolean applies(final int model){return false;}
			@Override
			public boolean refine(final Orf orf, final byte[] sequence, final int model, final float gc,
					final float identity, final boolean quantum, final int fixed, final int minimum, final int maximum){
				throw new AssertionError("The default-mode probe does not own any model");
			}
		};
		check(!defaultMode.afterAlignmentTrim(), "Existing refiner implementations must inherit pre-trim mode without source changes");
		for(int strand=0; strand<2; strand++){
			final TrimSpy off=new TrimSpy(true);
			final Orf raw=candidate(bases, strand, 40, 89);
			check(finish(off, raw, bases, 0) && off.trims==1 && raw.start==44 && raw.stop==85,
				"Null refiner must retain the ordinary successful trim and its true return value");
			final TrimSpy legacy=new TrimSpy(true);
			final RecordingRefiner before=new RecordingRefiner(false, true, 40, 89, 0);
			legacy.setBoundaryRefiner(before);
			final Orf old=candidate(bases, strand, 40, 89);
			check(!finish(legacy, old, bases, 0) && legacy.trims==0 && before.calls==1 && old.start==41 && old.stop==88,
				"Legacy pre-trim replacement must still bypass ordinary trim and return false");
			final TrimSpy after=new TrimSpy(true);
			final RecordingRefiner post=new RecordingRefiner(true, true, 44, 85, 0);
			after.setBoundaryRefiner(post);
			final int[] observed={-1, -1, 0};
			final boolean[] flags={false, false};
			after.setInstrumentSink((contig, orientation, model, ws, we, left, right, sequence, origin,
					trimmed, invoked, source, identity, length)->{
				observed[0]=left;
				observed[1]=right;
				observed[2]++;
				flags[0]=trimmed;
				flags[1]=invoked;
			});
			final Orf shifted=candidate(bases, strand, 40, 89);
			check(finish(after, shifted, bases, 0) && after.trims==1 && post.calls==1 && shifted.start==45 && shifted.stop==84,
				"Post-trim replacement must see the actual trim anchor once and preserve trimSucceeded");
			final ArrayList<Orf> candidates=new ArrayList<Orf>();
			candidates.add(shifted);
			after.commitSnapshotCandidates(candidates, new ArrayList<int[]>(), new ArrayList<Orf>());
			check(observed[0]==44 && observed[1]==85 && observed[2]==1 && flags[0] && flags[1],
				"Instrumentation must retain one pre-refinement trim anchor and truthful trim/invocation flags");
			final TrimSpy skipped=new TrimSpy(true);
			final RecordingRefiner absent=new RecordingRefiner(true, false, -1, -1, 0);
			skipped.setBoundaryRefiner(absent);
			final Orf unchanged=candidate(bases, strand, 40, 89);
			check(finish(skipped, unchanged, bases, 0) && skipped.trims==1 && absent.calls==0
				&& unchanged.start==raw.start && unchanged.stop==raw.stop,
				"A refiner that does not own the selected model cannot alter the ordinary OFF path");
			final TrimSpy failed=new TrimSpy(false);
			final RecordingRefiner inactive=new RecordingRefiner(true, true, -1, -1, 0);
			failed.setBoundaryRefiner(inactive);
			final Orf noTrim=candidate(bases, strand, 40, 89);
			check(!finish(failed, noTrim, bases, 0) && failed.trims==1 && inactive.calls==0 && noTrim.start==40 && noTrim.stop==89,
				"A failed/no-op trim without raw or fixed-end eligibility cannot invoke post-trim inference");
			final TrimSpy rawEligible=new TrimSpy(false);
			final java.lang.reflect.Field rawMode=NcrnaScavenger.class.getDeclaredField("boundaryOnRawEndpoints");
			rawMode.setAccessible(true);
			rawMode.setBoolean(rawEligible, true);
			final RecordingRefiner rawRefiner=new RecordingRefiner(true, true, 40, 89, 0);
			rawEligible.setBoundaryRefiner(rawRefiner);
			final Orf retainedAnchor=candidate(bases, strand, 40, 89);
			check(!finish(rawEligible, retainedAnchor, bases, 0) && rawRefiner.calls==1
				&& retainedAnchor.start==41 && retainedAnchor.stop==88,
				"Explicit raw-anchor eligibility allows inference without misreporting trim success");
			for(int mask=1; mask<=3; mask++){
				final int start=(mask&1)==0 ? 40 : 0, stop=(mask&2)==0 ? 89 : bases.length-1;
				final TrimSpy fixed=new TrimSpy(false);
				final RecordingRefiner interior=new RecordingRefiner(true, true, start, stop, mask);
				fixed.setBoundaryRefiner(interior);
				final Orf partial=candidate(bases, strand, start, stop);
				check(!finish(fixed, partial, bases, mask), "Post-trim fixed-end inference cannot pretend a trim succeeded");
				check(interior.calls==(mask==3 ? 0 : 1) && fixed.trims==(mask==3 ? 0 : 1),
					"One fixed endpoint enables its interior partner; two fixed endpoints bypass every refinement");
				check(partial.start==start+((mask&1)==0 ? 1 : 0) && partial.stop==stop-((mask&2)==0 ? 1 : 0),
					"Interior refinement must leave the physical truncation coordinate unchanged on both strands");
			}
		}
		final TrimSpy competing=new TrimSpy(true);
		competing.voteEnds=true;
		boolean rejected=false;
		try{competing.setBoundaryRefiner(new RecordingRefiner(true, true, 44, 85, 0));}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Competing vote-centred and post-trim alignment-centred refinement must fail during installation");
		check(Arrays.equals(bases, original), "No refiner or trim spy may mutate the caller's source sequence");
		System.out.println("NCRNA_POST_TRIM_REFINER_TEST_PASS");
	}

	private static Orf candidate(final byte[] bases, final int strand, final int start, final int stop){
		final Orf result=new Orf("fixture", start, stop, strand, 0, bases, false, ProkObject.RNA);
		result.ncrnaFamily="cobalamin_rf00174";
		result.trnaModel="model_description n=1";
		result.orfScore=83.25f;
		return result;
	}

	/** Reflection reaches only the existing private physical-edge overload; the
	 * public no-partial caller seam is used for every ordinary case.
	 */
	private static boolean finish(final NcrnaScavenger caller, final Orf orf, final byte[] bases, final int mask) throws Exception{
		final boolean result;
		if(mask==0){
			result=caller.finishAcceptedBoundary(orf, bases, 0, 0, bases.length-1, 0, .8f, orf.stop-orf.start+1, .8f, true, null, null);
		}else{
			final Method method=NcrnaScavenger.class.getDeclaredMethod("finishAcceptedBoundary", Orf.class, byte[].class,
				int.class, int.class, int.class, int.class, float.class, int.class, float.class,
				boolean.class, int[].class, long[].class, int.class, EndClippedAligner.Result.class);
			method.setAccessible(true);
			result=(Boolean)method.invoke(caller, orf, bases, 0, 0, bases.length-1, 0, .8f,
				orf.stop-orf.start+1, .8f, true, null, null, mask, null);
		}
		check(orf.orfScore==83.25f && "model_description n=1".equals(orf.trnaModel)
			&& "cobalamin_rf00174".equals(orf.ncrnaFamily), "Endpoint refinement must preserve detection score and original family/model identity");
		return result;
	}

	/** Supplies a known trim outcome independently of the dispatch being checked. */
	private static final class TrimSpy extends NcrnaScavenger{
		TrimSpy(final boolean success_){
			super(new byte[][]{"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".getBytes(StandardCharsets.US_ASCII)},
				null, new String[]{"model_description n=1"}, null, 17, 16, 10);
			success=success_;
		}
		@Override
		public boolean trimToAlignmentExtent(final Orf orf, final byte[] bases, final int model,
				final int start, final int stop, final float identity){
			trims++;
			if(success){
				orf.start=44;
				orf.stop=85;
			}
			return success;
		}
		final boolean success;
		int trims;
	}

	private static final class RecordingRefiner implements NcrnaBoundaryRefiner{
		RecordingRefiner(final boolean post_, final boolean applies_, final int start_, final int stop_, final int fixed_){
			post=post_;
			active=applies_;
			start=start_;
			stop=stop_;
			fixed=fixed_;
		}
		@Override
		public boolean afterAlignmentTrim(){return post;}
		@Override
		public boolean applies(final int model){return active;}
		@Override
		public boolean refine(final Orf orf, final byte[] bases, final int model, final float gc,
				final float identity, final boolean quantum, final int fixedMask, final int minimum, final int maximum){
			check(orf.start==start && orf.stop==stop && fixedMask==fixed && model==0 && identity==.8f
				&& gc==0 && quantum && minimum==16 && maximum==Integer.MAX_VALUE,
				"Refiner received the wrong anchor, fixed end, selected model, contig GC, identity or length interval");
			calls++;
			if((fixedMask&1)==0){orf.start++;}
			if((fixedMask&2)==0){orf.stop--;}
			return true;
		}
		final boolean post, active;
		final int start, stop, fixed;
		int calls;
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
}
