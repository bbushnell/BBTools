package test.prok;

import java.util.ArrayList;
import java.util.Arrays;

import ml.CellNet;
import prok.NcrnaBoundaryScorer;
import prok.NcrnaBoundaryScorer.RefinementResult;
import prok.NcrnaCdsStopConstraint;
import prok.Orf;
import prok.ProkObject;
import prok.TrnaBoundaryFeatures;

/** Finite coordinate and sequential-search oracles for the CDS stop constraint.
 * Uses the actual V3 feature builder with deterministic length-scoring networks.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCdsStopConstraintTest {

	public static void main(String[] args){
		boolean assertions=false;
		assert(assertions=true);
		check(assertions, "Run finite invariants with assertions enabled");
		coordinates();
		sequential();
		geometry();
		NcrnaBoundarySweepGeometryTest.main(new String[0]);
		System.out.println("PASS NcrnaCdsStopConstraintTest coordinates_reset_workers_both_orders_edges_legal_context");
	}

	/** Mirrors native candidate creation: sense coordinates, then one minus flip. */
	private static Orf cds(int start, int stop, int strand, byte[] bases){
		final Orf result=new Orf("fixture", start, stop, strand, 0, bases, false, ProkObject.CDS);
		if(strand==1){result.flip();}
		return result;
	}

	private static void coordinates(){
		final byte[] bases=new byte[1000];
		Arrays.fill(bases, (byte)'A');
		@SuppressWarnings("unchecked")
		final ArrayList<Orf>[] candidates=new ArrayList[6];
		candidates[0]=new ArrayList<Orf>();
		candidates[3]=new ArrayList<Orf>();
		for(int strand=0; strand<2; strand++){
			candidates[strand*3].add(cds(300, 500, strand, bases));
			candidates[strand*3].add(cds(205, 400, strand, bases));
			candidates[strand*3].add(cds(50, 150, strand, bases));
		}
		final NcrnaCdsStopConstraint.Factory factory=new NcrnaCdsStopConstraint.Factory(null);
		final NcrnaCdsStopConstraint first=factory.create(), second=factory.create();
		first.prepare(candidates, bases.length);
		second.prepare(null, bases.length);
		for(int strand=0; strand<2; strand++){
			check(first.ceiling(strand, 100, 199, 1000)==204, "Next CDS begins at205 in either sense frame");
			check(first.ceiling(strand, 100, 210, 1000)==210, "Existing overlap must keep offset zero");
			check(first.ceiling(strand, 205, 220, 1000)==299, "CDS starts at the RNA start are not downstream");
			check(first.ceiling(strand, 501, 600, 1000)==Integer.MAX_VALUE, "Upstream CDS imposes no stop ceiling");
			check(second.ceiling(strand, 100, 199, 1000)==Integer.MAX_VALUE, "Workers must not share mutable contig state");
		}
		candidates[0].clear();
		first.prepare(candidates, 1000);
		check(first.ceiling(0, 100, 199, 1000)==Integer.MAX_VALUE, "Opposite-strand CDS must not constrain a plus RNA");
		first.prepare(null, 800);
		check(first.ceiling(1, 100, 199, 800)==Integer.MAX_VALUE, "Every new contig clears both strand indices");
		check(factory.totals()[0]==0, "Preparing CDS indices does not report NN sweep events");
	}

	private static void sequential(){
		for(boolean stopFirst : new boolean[]{false, true}){
			final RefinementResult raw=score(new LengthNet(stopFirst), 100, 199, 500, 32,
				Integer.MAX_VALUE, false, false, 0);
			check(stopFirst ? raw.stopOffset==20 && raw.startOffset==20 : raw.startOffset==10 && raw.stopOffset==30,
				"Unbounded length targets have independently calculated sequential offsets");
			final RefinementResult limited=score(new LengthNet(stopFirst), 100, 199, 500, 32,
				204, false, false, 0);
			check(limited.stopOffset==5 && limited.startOffset==(stopFirst ? 5 : 10),
				"Both orders respect the ceiling; a stop-first start sees the legal105base span");
			check(limited.stopSitesExcluded==27 && limited.stopSitesScored==38 && limited.startSitesScored==65,
				"Only stop offsets6through32 are excluded; all65start sites remain");
			final RefinementResult overlap=score(new LengthNet(stopFirst), 100, 199, 500, 32,
				199, false, false, 0);
			check(overlap.stopOffset==0 && overlap.stopSitesExcluded==32, "Existing overlap forbids only outward extension");
		}
	}

	private static void geometry(){
		final RefinementResult shortSpan=score(new FlatNet(), 100, 115, 400, 10, 115, false, false, 0);
		check(shortSpan.startSitesScored==11 && shortSpan.stopSitesScored==1 && shortSpan.stopSitesExcluded==10,
			"A16base candidate preserves minimum length, zero stop and all legal start sites");
		final RefinementResult fixedStart=score(new FlatNet(), 100, 199, 400, 10, 199, false, false, 1);
		check(fixedStart.startSitesScored==0 && fixedStart.stopSitesScored==11 && fixedStart.stopSitesExcluded==10,
			"A physically fixed start leaves the constrained stop search active");
		final RefinementResult fixedStop=score(new FlatNet(), 100, 199, 400, 10, 199, false, false, 2);
		check(fixedStop.stopSitesScored==0 && fixedStop.stopSitesExcluded==0 && fixedStop.startSitesScored==21,
			"A physically fixed stop has no search exclusions to count");
		final RefinementResult fixedBoth=score(new FlatNet(), 100, 199, 400, 10, 199, false, false, 3);
		check(fixedBoth.startSitesScored==0 && fixedBoth.stopSitesScored==0, "Both physical ends stay untouched");
		final RefinementResult edge=score(new FlatNet(), 100, 195, 200, 10, 197, false, true, 0);
		check(edge.stopSitesScored==13 && edge.stopSitesExcluded==2,
			"Physical edge clips six sites; only two available sites are deliberately excluded");
		final RefinementResult allowed=score(new FlatNet(), 100, 195, 200, 10, 197, false, false, 0);
		check(allowed.stopSitesScored==13, "Unavailable forbidden offsets need no sequence context");
		boolean failed=false;
		try{score(new FlatNet(), 100, 195, 200, 10, 201, false, false, 0);}
		catch(AssertionError e){failed=e.getMessage().contains("missing context for a valid candidate");}
		check(failed, "Missing context at a permitted offset must still fail loudly");
		failed=false;
		try{score(new FlatNet(), 100, 195, 400, 10, 194, false, false, 0);}
		catch(IllegalArgumentException e){failed=e.getMessage().contains("offset zero");}
		check(failed, "The ceiling must never remove the unchanged RNA stop");
	}

	private static RefinementResult score(CellNet net, int start, int stop, int length, int radius,
			int ceiling, boolean left, boolean right, int fixedMask){
		final byte[] bases=new byte[length];
		Arrays.fill(bases, (byte)'A');
		final int[] offsets=new int[2*radius+1];
		for(int i=0; i<offsets.length; i++){offsets[i]=i-radius;}
		final TrnaBoundaryFeatures.NinemerTable table=(seq, pos)->0;
		return NcrnaBoundaryScorer.refineBoundariesWithCounts(net, net, bases, start, stop,
			new byte[100], null, table, table, 9, 0, 9, 0, .5f, 100f, offsets, offsets,
			0, 0, left, right, NcrnaBoundaryScorer.FEATURES_V3, .8f, null, null,
			null, null, fixedMask, ceiling);
	}

	private static class FlatNet extends CellNet {
		@Override
		public CellNet applyInput(float[] values){
			NcrnaCdsStopConstraintTest.check(values.length==10, "Fixture consumes the real V3 feature layout");
			input=values;
			return this;
		}
		@Override
		public float feedForward(){return 0;}
		float[] input;
	}

	private static final class LengthNet extends FlatNet {
		LengthNet(boolean stopFirst_){stopFirst=stopFirst_;}
		@Override
		public float feedForward(){
			final boolean weak=(input[3]==1)==stopFirst;
			final float target=stopFirst ? (weak ? 120 : 100) : (weak ? 90 : 120);
			return (weak ? 1 : 10)-Math.abs(input[4]*100-target)/100;
		}
		final boolean stopFirst;
	}

	private static void check(boolean ok, String reason){
		if(!ok){throw new AssertionError(reason);}
	}
}
