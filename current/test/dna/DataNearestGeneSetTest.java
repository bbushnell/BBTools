package test.dna;

import java.util.ArrayList;

import dna.Data;
import dna.Exon;
import dna.Gene;
import dna.GeneSet;
import shared.Shared;
import structures.Range;

/** Literal inclusive-interval checks against a private replacement of one public range-cache slot.
 * No annotation resources are loaded. Restores the exact prior cache reference on
 * success or failure. Expected payloads are explicit references, not computed by
 * the production range builder or nearest-distance algorithm.
 * @author Sucrose
 */
public final class DataNearestGeneSetTest{

	/** Runs case=controls|gaps in a fresh JVM under either assertion mode. */
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("case="), "Require one case= argument");
		final String mode=args[0].substring(5);
		require("controls".equals(mode) || "gaps".equals(mode), "Unknown interval case: "+mode);
		final Range[] prior=Data.geneSetRangeMatrix[CHROM];
		try{
			final DataNearestGeneSetTest test=new DataNearestGeneSetTest();
			Data.geneSetRangeMatrix[CHROM]=test.ranges;
			if("controls".equals(mode)){test.controls();}
			else{test.gaps();}
		}finally{
			Data.geneSetRangeMatrix[CHROM]=prior;
		}
		require(Data.geneSetRangeMatrix[CHROM]==prior, "Fixture must restore the prior range-cache reference");
		System.out.println("PASS "+mode);
	}

	/** Creates valid tiny records and independently specified padded-range bounds. */
	private DataNearestGeneSetTest(){
		left=new GeneSet[]{geneSet("L1", 1, 200, 209), geneSet("L2", 2, 205, 209)};
		right=new GeneSet[]{geneSet("R1", 3, 800, 809), geneSet("R2", 4, 805, 809)};
		tail=new GeneSet[]{geneSet("T", 5, 1400, 1409)};
		ranges=new Range[]{new Range(0, 409, left), new Range(600, 1009, right), new Range(1200, 1609, tail)};
	}

	/** Preserves containment, inclusive overlap, end cases, payload order, emptiness and the point overload. */
	private void controls(){
		boolean ok=true;
		ok&=check("containing", 201, 205, left, true);
		ok&=check("left-endpoint", 409, 409, left, true);
		ok&=check("right-endpoint", 600, 600, right, true);
		ok&=check("right-overlap-only", 580, 620, right, false);
		ok&=check("two-overlaps", 400, 610, new GeneSet[]{left[0], left[1], right[0], right[1]}, false);
		ok&=check("all-overlaps", 0, 1609, new GeneSet[]{left[0], left[1], right[0], right[1], tail[0]}, false);
		ok&=check("before-first", -10, -1, left, true);
		ok&=check("after-last", 1610, 1620, tail, true);
		final GeneSet[] point=Data.getNearestGeneSets(CHROM, 585);
		final boolean pointOk=point==left && sameSets(point, left);
		System.out.println("point585="+labels(point)+" pass="+pointOk);
		ok&=pointOk;
		Data.geneSetRangeMatrix[CHROM]=new Range[0];
		ok&=check("empty", 10, 20, null, false);
		Data.geneSetRangeMatrix[CHROM]=ranges;
		require(ok, "Interval controls must retain complete range payloads and point-overload behavior");
	}

	/** Explicit distances: left 1 vs 180, right 171 vs 10, and tie 91 vs 91. */
	private void gaps(){
		boolean ok=true;
		ok&=check("gap-nearer-left", 410, 420, left, true);
		ok&=check("gap-nearer-right", 580, 590, right, true);
		ok&=check("gap-equal-left-tie", 500, 509, left, true);
		require(ok, "Gap queries must choose the nearer merged range, preserving left ties and whole shared payloads");
	}

	/** Compares literal payload references and order; selected single ranges retain their shared array. */
	private static boolean check(String label, int start, int stop, GeneSet[] expected, boolean shared){
		final GeneSet[] actual=Data.getNearestGeneSets(CHROM, start, stop);
		final boolean ok=sameSets(actual, expected) && (!shared || actual==expected);
		System.out.println(label+" ["+start+","+stop+"] actual="+labels(actual)+
				" expected="+labels(expected)+" pass="+ok);
		return ok;
	}

	/** Reference equality prevents silently cloning, dropping, reordering or selecting individual members. */
	private static boolean sameSets(GeneSet[] actual, GeneSet[] expected){
		if(actual==null || expected==null){return actual==expected;}
		if(actual.length!=expected.length){return false;}
		for(int i=0; i<actual.length; i++){
			if(actual[i]!=expected[i]){return false;}
		}
		return true;
	}

	/** Prints compact payload labels without using annotation formatting as an oracle. */
	private static String labels(GeneSet[] sets){
		if(sets==null){return "null";}
		final StringBuilder text=new StringBuilder();
		for(GeneSet set : sets){
			if(text.length()>0){text.append(',');}
			text.append(set==null ? "null" : set.name);
		}
		return text.toString();
	}

	/** Constructs a single fully coding transcript; coordinates and expected padded bounds are separate literals. */
	private static GeneSet geneSet(String name, int id, int start, int stop){
		final Exon[] exons={new Exon(start, stop, (byte)CHROM, Shared.PLUS, false, true)};
		final Gene gene=new Gene((byte)CHROM, Shared.PLUS, start, stop, start, stop, id,
				name, name+"_tx", null, null, null, exons, false, false, true,
				"seqGene", "fixture", "fixture");
		final ArrayList<Gene> genes=new ArrayList<Gene>(1);
		genes.add(gene);
		return new GeneSet(name, genes);
	}

	/** Keeps all finite-oracle checks active under -da. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Cache slot replaced only for this fresh-JVM test and restored in finally. */
	private static final int CHROM=1;
	/** Two-member left range payload, including a member not individually overlapping every query. */
	private final GeneSet[] left;
	/** Two-member right range payload used to catch truncation to one nearest gene. */
	private final GeneSet[] right;
	/** Last-range payload for boundary controls. */
	private final GeneSet[] tail;
	/** Sorted, disjoint inclusive ranges with known independent gap distances. */
	private final Range[] ranges;
}
