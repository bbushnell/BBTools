package test.prok;

import java.util.ArrayList;
import java.util.Arrays;

import map.LongHashSet;
import prok.NcrnaScavenger;

/** Checks the default-off observation seam before an intentionally rejecting index.
 * No model alignment is needed to establish the capture event's position.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWindowCaptureTest{

	public static void main(final String[] args){
		final byte[] bases=new byte[32];
		Arrays.fill(bases, (byte)'a');
		final byte[] before=bases.clone();
		final NcrnaScavenger caller=fixture();
		final int[] events={0};
		assert(caller.scavenge("contig", bases, 0, new ArrayList<int[]>(), new int[]{8}).isEmpty())
			: "The deliberately high index floor must reject this seed-supported window";
		caller.setWindowCaptureSink((name, strand, pass, start, stop, source, contig, window)->{
			events[0]++;
			assert(name.equals("contig") && strand==0 && pass==1 && start==0 && stop==31)
				: "The observer must retain the actual scheduled event and oriented inclusive coordinates";
			assert(contig==bases && window!=bases && window.length==32)
				: "Capture borrows the original contig and the caller's separate uppercase window";
			for(byte base:window){assert(base=='A') : "Capture must see the same uppercase bytes used by the index";}
		});
		assert(caller.scavenge("contig", bases, 0, new ArrayList<int[]>(), new int[]{8}).isEmpty())
			: "Observation must not bypass the rejecting index";
		assert(events[0]==1 && caller.alignmentCount()==0)
			: "Capture must fire once before index rejection, without creating an alignment";
		assert(Arrays.equals(before, bases)) : "Observing an uppercase private window must not alter soft-masked source bases";
		caller.setWindowCaptureSink(null);
		caller.scavenge("contig", bases, 0, new ArrayList<int[]>(), new int[]{8});
		assert(events[0]==1) : "Disarming capture must remove observation from subsequent windows";
		System.err.println("NCRNA_WINDOW_CAPTURE_TEST_PASS");
	}

	/** A32 has one distinct9-mer, far below the explicit99-key index floor. */
	private static NcrnaScavenger fixture(){
		final byte[] model=new byte[32];
		Arrays.fill(model, (byte)'A');
		final LongHashSet seeds=new LongHashSet();
		seeds.add(0L);
		final NcrnaScavenger caller=new NcrnaScavenger(new byte[][]{model}, null,
			new String[]{"fixture"}, seeds, 17, 20, 16,
			9, 1, false, 0f, 0f, 0f, 99, 0f, 1f, .64f, .64f);
		caller.family="cobalamin_rf00174";
		caller.strictIndexCutoff=true;
		return caller;
	}
}
