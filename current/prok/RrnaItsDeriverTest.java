package prok;

import java.util.ArrayList;

/** Focused, Java-8-safe checks for RrnaItsDeriver's pure interval math and the ported
 * CallITS mutual-nearest-neighbor pairing algorithm. No file I/O, no HBM/genome loading,
 * no production wiring -- exercises the standalone class directly. Covers the design
 * doc's F01/F02 (plus/minus ordered triple), F03-F06 (mechanical rejections), F09/F10
 * (partial anchors -- ITS1-only, ITS2-only), F12 (missing anchors, no fabricated gap),
 * and an ambiguity case (competing anchors -> no pairing, per the ported tie-detection).
 * @author G11
 */
public class RrnaItsDeriverTest{

	public static void main(String[] args){
		testPlusOrderedTriple();
		testMinusOrderedTriple();
		testCrossContig();
		testMixedStrand();
		testOverlap();
		testZeroLength();
		testIts1Only();
		testIts2Only();
		testMissingAnchors();
		testAmbiguousCompetingAnchors();
		testMaxSpan();
		testLongDescriptionSeqidTruncation();
		System.out.println("PASS RrnaItsDeriverTest");
	}

	/*--------------------------------------------------------------*/

	private static RrnaItsDeriver.Anchor anchor(String seqid, int start, int stop, byte strand, String fam, String id){
		return new RrnaItsDeriver.Anchor(seqid, start, stop, strand, fam, id, 0, 0f);
	}

	/** Builds a plus-strand synthetic contig: 18S[0,19], gap1(ITS1) sentinel 'A'x10,
	 * 5.8S[30,49], gap2(ITS2) sentinel 'C'x10, LSU[60,79]. Distinct sentinel bases per
	 * gap so a byte-level check proves the extracted interval, not just its coordinates. */
	private static byte[] plusGenome(){
		byte[] b=new byte[100];
		java.util.Arrays.fill(b, (byte)'T');
		for(int i=20; i<30; i++){b[i]='A';}//ITS1
		for(int i=50; i<60; i++){b[i]='C';}//ITS2
		return b;
	}

	private static void testPlusOrderedTriple(){
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor s18=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "s18_a");
		RrnaItsDeriver.Anchor r58=anchor("ctg1", 30, 49, (byte)0, RrnaItsDeriver.R58, "r58_a");
		RrnaItsDeriver.Anchor lsu=anchor("ctg1", 60, 79, (byte)0, RrnaItsDeriver.LSU, "lsu_a");
		ArrayList<RrnaItsDeriver.Anchor> group=new ArrayList<RrnaItsDeriver.Anchor>();
		group.add(s18); group.add(r58); group.add(lsu);

		ArrayList<RrnaItsDeriver.Derived> derived=RrnaItsDeriver.resolveGroup(group, bases, true, true, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		RrnaItsDeriver.Derived its1=find(derived, "ITS1"), its2=find(derived, "ITS2"), its=find(derived, "ITS");
		if(its1.start!=20 || its1.stop!=29) throw new AssertionError("ITS1 coords "+its1.start+"-"+its1.stop);
		if(its2.start!=50 || its2.stop!=59) throw new AssertionError("ITS2 coords "+its2.start+"-"+its2.stop);
		if(its.start!=20 || its.stop!=59) throw new AssertionError("combined ITS coords "+its.start+"-"+its.stop);
		if(!allBases(its1.sequence, (byte)'A')) throw new AssertionError("ITS1 sequence not all sentinel A");
		if(!allBases(its2.sequence, (byte)'C')) throw new AssertionError("ITS2 sequence not all sentinel C");
		if(its.sequence.length!=40) throw new AssertionError("combined ITS length "+its.sequence.length);
	}

	private static void testMinusOrderedTriple(){
		//Same biological order (18S->5.8S->LSU) represented at DESCENDING genomic
		//coordinates: LSU genomic-left, 18S genomic-right, both on strand=1.
		byte[] bases=plusGenome();//reuse: LSU at [0,19], gap 'A' at [20,29], 5.8S [30,49], gap 'C' [50,59], 18S [60,79]
		RrnaItsDeriver.Anchor lsu=anchor("ctg1", 0, 19, (byte)1, RrnaItsDeriver.LSU, "lsu_m");
		RrnaItsDeriver.Anchor r58=anchor("ctg1", 30, 49, (byte)1, RrnaItsDeriver.R58, "r58_m");
		RrnaItsDeriver.Anchor s18=anchor("ctg1", 60, 79, (byte)1, RrnaItsDeriver.S18, "s18_m");
		ArrayList<RrnaItsDeriver.Anchor> group=new ArrayList<RrnaItsDeriver.Anchor>();
		group.add(lsu); group.add(r58); group.add(s18);

		ArrayList<RrnaItsDeriver.Derived> derived=RrnaItsDeriver.resolveGroup(group, bases, true, true, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		RrnaItsDeriver.Derived its1=find(derived, "ITS1"), its2=find(derived, "ITS2");
		//ITS1 = 18S(genomic-right, upstream in transcript order on minus)-to-5.8S gap = genomic [50,59], sentinel 'C'
		if(its1.start!=50 || its1.stop!=59) throw new AssertionError("minus ITS1 coords "+its1.start+"-"+its1.stop);
		//ITS2 = 5.8S-to-LSU gap = genomic [20,29], sentinel 'A'
		if(its2.start!=20 || its2.stop!=29) throw new AssertionError("minus ITS2 coords "+its2.start+"-"+its2.stop);
		//Coordinates stay ascending-genomic (no inversion); orientation is applied only to sequence.
		if(!allBases(revcomp(its1.sequence), (byte)'C')) throw new AssertionError("minus ITS1 sequence not reverse-complemented sentinel C");
		if(!allBases(revcomp(its2.sequence), (byte)'A')) throw new AssertionError("minus ITS2 sequence not reverse-complemented sentinel A");
	}

	private static void testCrossContig(){
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor a=anchor("ctgA", 0, 19, (byte)0, RrnaItsDeriver.S18, "a");
		RrnaItsDeriver.Anchor b=anchor("ctgB", 30, 49, (byte)0, RrnaItsDeriver.R58, "b");
		RrnaItsDeriver.Result r=RrnaItsDeriver.derive("ITS1", a, b, bases, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(r.reason!=RrnaItsDeriver.Reason.CROSS_CONTIG) throw new AssertionError("expected CROSS_CONTIG, got "+r.reason);
	}

	private static void testMixedStrand(){
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor a=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "a");
		RrnaItsDeriver.Anchor b=anchor("ctg1", 30, 49, (byte)1, RrnaItsDeriver.R58, "b");
		RrnaItsDeriver.Result r=RrnaItsDeriver.derive("ITS1", a, b, bases, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(r.reason!=RrnaItsDeriver.Reason.MIXED_STRAND) throw new AssertionError("expected MIXED_STRAND, got "+r.reason);
	}

	private static void testOverlap(){
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor a=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "a");
		RrnaItsDeriver.Anchor b=anchor("ctg1", 10, 29, (byte)0, RrnaItsDeriver.R58, "b");//overlaps a
		RrnaItsDeriver.Result r=RrnaItsDeriver.derive("ITS1", a, b, bases, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(r.reason!=RrnaItsDeriver.Reason.OVERLAP) throw new AssertionError("expected OVERLAP, got "+r.reason);
	}

	private static void testZeroLength(){
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor a=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "a");
		RrnaItsDeriver.Anchor b=anchor("ctg1", 20, 39, (byte)0, RrnaItsDeriver.R58, "b");//adjacent, no gap
		RrnaItsDeriver.Result r=RrnaItsDeriver.derive("ITS1", a, b, bases, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(r.reason!=RrnaItsDeriver.Reason.ZERO_LENGTH) throw new AssertionError("expected ZERO_LENGTH, got "+r.reason);
	}

	private static void testIts1Only(){
		//18S+5.8S only: ITS1 permitted, no ITS2/combined.
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor s18=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "s18");
		RrnaItsDeriver.Anchor r58=anchor("ctg1", 30, 49, (byte)0, RrnaItsDeriver.R58, "r58");
		ArrayList<RrnaItsDeriver.Anchor> group=new ArrayList<RrnaItsDeriver.Anchor>();
		group.add(s18); group.add(r58);
		ArrayList<RrnaItsDeriver.Derived> derived=RrnaItsDeriver.resolveGroup(group, bases, true, true, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(derived.size()!=1 || !derived.get(0).product.equals("ITS1")) throw new AssertionError("expected exactly 1 ITS1, got "+derived.size());
	}

	private static void testIts2Only(){
		//5.8S+LSU only: ITS2 permitted, no ITS1; no combined (needs 18S).
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor r58=anchor("ctg1", 30, 49, (byte)0, RrnaItsDeriver.R58, "r58");
		RrnaItsDeriver.Anchor lsu=anchor("ctg1", 60, 79, (byte)0, RrnaItsDeriver.LSU, "lsu");
		ArrayList<RrnaItsDeriver.Anchor> group=new ArrayList<RrnaItsDeriver.Anchor>();
		group.add(r58); group.add(lsu);
		ArrayList<RrnaItsDeriver.Derived> derived=RrnaItsDeriver.resolveGroup(group, bases, true, true, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(derived.size()!=1 || !derived.get(0).product.equals("ITS2")) throw new AssertionError("expected exactly 1 ITS2, got "+derived.size());
	}

	private static void testMissingAnchors(){
		byte[] bases=plusGenome();
		ArrayList<RrnaItsDeriver.Anchor> empty=new ArrayList<RrnaItsDeriver.Anchor>();
		if(!RrnaItsDeriver.resolveGroup(empty, bases, true, true, RrnaItsDeriver.DEFAULT_MAX_SPAN).isEmpty())
			throw new AssertionError("empty group fabricated a derived record");
		ArrayList<RrnaItsDeriver.Anchor> oneOnly=new ArrayList<RrnaItsDeriver.Anchor>();
		oneOnly.add(anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "s18"));
		if(!RrnaItsDeriver.resolveGroup(oneOnly, bases, true, true, RrnaItsDeriver.DEFAULT_MAX_SPAN).isEmpty())
			throw new AssertionError("single anchor fabricated a derived record");
	}

	private static void testAmbiguousCompetingAnchors(){
		//Two overlapping 5.8S candidates -> same-family overlap -> both selfAmbiguous ->
		//excluded from ANY pairing, per the ported CallITS semantics (unresolved copy
		//identity, not resolved by picking the nearer one).
		byte[] bases=new byte[200]; java.util.Arrays.fill(bases, (byte)'T');
		RrnaItsDeriver.Anchor s18=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "s18");
		RrnaItsDeriver.Anchor r58a=anchor("ctg1", 30, 39, (byte)0, RrnaItsDeriver.R58, "r58a");
		RrnaItsDeriver.Anchor r58b=anchor("ctg1", 35, 44, (byte)0, RrnaItsDeriver.R58, "r58b");//overlaps r58a
		ArrayList<RrnaItsDeriver.Anchor> group=new ArrayList<RrnaItsDeriver.Anchor>();
		group.add(s18); group.add(r58a); group.add(r58b);
		ArrayList<RrnaItsDeriver.Derived> derived=RrnaItsDeriver.resolveGroup(group, bases, false, true, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		for(RrnaItsDeriver.Derived d : derived){
			if(d.product.equals("ITS1")) throw new AssertionError("same-family-overlap ambiguity was not excluded from pairing");
		}
	}

	private static void testMaxSpan(){
		byte[] bases=new byte[20000]; java.util.Arrays.fill(bases, (byte)'T');
		RrnaItsDeriver.Anchor a=anchor("ctg1", 0, 19, (byte)0, RrnaItsDeriver.S18, "a");
		RrnaItsDeriver.Anchor b=anchor("ctg1", 15000, 15019, (byte)0, RrnaItsDeriver.R58, "b");//gap=14980 > DEFAULT_MAX_SPAN
		RrnaItsDeriver.Result r=RrnaItsDeriver.derive("ITS1", a, b, bases, RrnaItsDeriver.DEFAULT_MAX_SPAN);
		if(r.reason!=RrnaItsDeriver.Reason.MAX_SPAN) throw new AssertionError("expected MAX_SPAN, got "+r.reason);
		//maxSpan=-1 means unlimited -- the same pair must succeed.
		RrnaItsDeriver.Result r2=RrnaItsDeriver.derive("ITS1", a, b, bases, -1);
		if(r2.reason!=RrnaItsDeriver.Reason.OK) throw new AssertionError("maxSpan=-1 should be unlimited, got "+r2.reason);
	}

	/** Real-genome-shaped seqid, e.g. "NC_026746.1 Cryptococcus neoformans var. grubii H99
	 * chromosome 2, complete sequence" (an actual header this bug was caught against on a
	 * real Dori canary run, 2026-09-17, G11+Ganyu): descriptive text after the accession,
	 * INCLUDING a comma, which would corrupt the GFF's own comma-delimited attribute column
	 * if not truncated at the accession. GFF column 1 and the outits FASTA header must agree
	 * on the SAME truncated seqid so a consumer can join them by contig name. */
	private static void testLongDescriptionSeqidTruncation(){
		final String longSeqid="ctg1 Some Organism, chromosome 2, complete sequence";
		byte[] bases=plusGenome();
		RrnaItsDeriver.Anchor s18=anchor(longSeqid, 0, 19, (byte)0, RrnaItsDeriver.S18, "s18");
		RrnaItsDeriver.Anchor r58=anchor(longSeqid, 30, 49, (byte)0, RrnaItsDeriver.R58, "r58");
		ArrayList<RrnaItsDeriver.Anchor> group=new ArrayList<RrnaItsDeriver.Anchor>();
		group.add(s18); group.add(r58);
		RrnaItsDeriver.Derived its1=find(RrnaItsDeriver.resolveGroup(group, bases, false, true, RrnaItsDeriver.DEFAULT_MAX_SPAN), "ITS1");

		structures.ByteBuilder bb=new structures.ByteBuilder();
		its1.appendGff(bb);
		String gffCol1=bb.toString().split("\t")[0];
		if(!gffCol1.equals("ctg1")) throw new AssertionError("GFF seqid not truncated at first space/comma: '"+gffCol1+"'");

		String fastaHeader=its1.toRead().id;
		String fastaCol1=fastaHeader.split("\t")[0];
		if(!fastaCol1.equals("ctg1")) throw new AssertionError("FASTA header contig field not truncated: '"+fastaCol1+"'");
		if(!fastaCol1.equals(gffCol1)) throw new AssertionError("GFF seqid and FASTA header contig field disagree: '"+gffCol1+"' vs '"+fastaCol1+"'");
	}

	/*--------------------------------------------------------------*/

	private static RrnaItsDeriver.Derived find(ArrayList<RrnaItsDeriver.Derived> list, String product){
		for(RrnaItsDeriver.Derived d : list){if(d.product.equals(product)) return d;}
		throw new AssertionError("no "+product+" in result ("+list.size()+" total)");
	}

	private static boolean allBases(byte[] seq, byte b){
		for(byte x : seq){if(x!=b) return false;}
		return seq.length>0;
	}

	private static byte[] revcomp(byte[] seq){
		return dna.AminoAcid.reverseComplementBases(seq);
	}
}
