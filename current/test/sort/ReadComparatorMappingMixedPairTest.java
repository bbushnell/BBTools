package test.sort;

import sort.ReadComparatorMapping;
import shared.Shared;
import stream.Read;

/** Baseline and fixed oracle for coordinate sorting mixed paired and unpaired roots.
 * This mirrors the root shape produced by odd single-file interleaved FASTA input:
 * one linked pair followed by an unmatched final first mate.
 * @author Sucrose
 */
public final class ReadComparatorMappingMixedPairTest{

	/** Runs expect=baseline or expect=fixed in a fresh test JVM. */
	public static void main(String[] args){
		if(args.length!=1 || (!args[0].equals("expect=baseline") && !args[0].equals("expect=fixed"))){
			throw new IllegalArgumentException("Expected expect=baseline or expect=fixed");
		}
		checkOddInterleavedFastaShape(args[0].equals("expect=fixed"));
		checkMappedPairVsMateOnlyCycle(args[0].equals("expect=fixed"));
		checkMatchStringResidual(args[0].equals("expect=fixed"));
		if(args[0].equals("expect=fixed")){
			checkFixedLegalMappingGrid();
			checkFixedTiePreservation();
		}
	}

	/** Current baseline dereferences the unpaired read's null mate when the paired
	 * root is the first compare argument. Fixed behavior must compare both orders. */
	private static void checkOddInterleavedFastaShape(boolean fixed){
		final Read paired=read("pair/1", 0);
		final Read mate=read("pair/2", 0);
		final Read unpaired=read("odd-tail/1", 1);
		paired.mate=mate;
		mate.mate=paired;
		mate.setPairnum(1);

		final ReadComparatorMapping comparator=new ReadComparatorMapping();
		Integer forward=null;
		boolean nullPointer=false;
		try{
			forward=comparator.compare(paired, unpaired);
		}catch(NullPointerException e){
			nullPointer=true;
		}

		if(!fixed){
			if(!nullPointer){
				throw new AssertionError("Expected baseline NullPointerException, observed "+forward);
			}
			System.out.println("PASS\tbaseline\todd-interleaved-fasta-shape\tNullPointerException");
			return;
		}

		if(nullPointer){throw new AssertionError("Fixed comparator still dereferenced a null mate");}
		final int reverse=comparator.compare(unpaired, paired);
		if(Integer.signum(forward)!=-Integer.signum(reverse)){
			throw new AssertionError("Mixed paired/unpaired comparison is not antisymmetric: "+
					forward+", "+reverse);
		}
		System.out.println("PASS\tfixed\todd-interleaved-fasta-shape\t"+forward+","+reverse);
	}

	/** Baseline switches anchors: fully mapped pairs compare by primary against each
	 * other but by mate against an unmapped-root/mapped-mate fragment. */
	private static void checkMappedPairVsMateOnlyCycle(boolean fixed){
		final Read a=mappedPairRoot("a", 100, 300, 0);
		final Read b=mappedPairRoot("b", 150, 200, 1);
		final Read c=mateMappedRoot("c", 250, 2);

		final ReadComparatorMapping comparator=new ReadComparatorMapping();
		if(!fixed){
			checkCompare("baseline", "mapped-pair-mate-only-a-b", -1, comparator, a, b);
			checkCompare("baseline", "mapped-pair-mate-only-b-c", -1, comparator, b, c);
			checkCompare("baseline", "mapped-pair-mate-only-c-a", -1, comparator, c, a);
			System.out.println("PASS\tbaseline\tmapped-pair-mate-only-cycle");
		}else{
			checkOrder("fixed", "mapped-pair-mate-only", comparator, a, b, c);
		}
	}

	/** Fixed ordering first compares stable anchor coordinate and then ranks every
	 * legal mapped/unmapped and linked/unlinked root shape at equal anchors. */
	private static void checkFixedLegalMappingGrid(){
		final Read[] reads={
				mappedPairRoot("plus100-both-paired", Shared.PLUS, 100, 900, 70),
				bothMappedUnpairedRoot("plus100-both-unpaired", Shared.PLUS, 100, 900, 60),
				primaryMappedLinkedRoot("plus100-primary-linked", Shared.PLUS, 100, 50),
				mappedRead("plus100-primary-unlinked", 40, Shared.PLUS, 100),
				mateMappedRoot("plus100-mate-only", Shared.PLUS, 100, 30),

				mappedPairRoot("plus101-both-paired", Shared.PLUS, 101, 901, 29),
				bothMappedUnpairedRoot("plus101-both-unpaired", Shared.PLUS, 101, 901, 28),
				primaryMappedLinkedRoot("plus101-primary-linked", Shared.PLUS, 101, 27),
				mappedRead("plus101-primary-unlinked", 26, Shared.PLUS, 101),
				mateMappedRoot("plus101-mate-only", Shared.PLUS, 101, 25),

				mappedPairRoot("minus100-both-paired", Shared.MINUS, 100, 900, 24),
				bothMappedUnpairedRoot("minus100-both-unpaired", Shared.MINUS, 100, 900, 23),
				primaryMappedLinkedRoot("minus100-primary-linked", Shared.MINUS, 100, 22),
				mappedRead("minus100-primary-unlinked", 21, Shared.MINUS, 100),
				mateMappedRoot("minus100-mate-only", Shared.MINUS, 100, 20),

				unmappedLinkedRoot("unmapped-linked", 19),
				read("unmapped-unlinked", 18)
		};

		final ReadComparatorMapping comparator=new ReadComparatorMapping();
		checkAllOrdered("fixed", "legal-shape-grid", comparator, reads);
	}

	/** Shape ties keep the historical detail read instead of always consulting the mate. */
	private static void checkFixedTiePreservation(){
		final ReadComparatorMapping comparator=new ReadComparatorMapping();

		checkCompare("fixed", "both-mapped-unpaired-primary-detail", -1, comparator,
				bothMappedUnpairedRoot("both-a", Shared.PLUS, 100, 900, 1),
				bothMappedUnpairedRoot("both-b", Shared.PLUS, 100, 100, 2));

		final Read primaryLinkedA=primaryMappedLinkedRoot("primary-a", Shared.PLUS, 100, 1);
		final Read primaryLinkedB=primaryMappedLinkedRoot("primary-b", Shared.PLUS, 100, 2);
		primaryLinkedB.mate.setStrand(Shared.PLUS);
		checkCompare("fixed", "primary-linked-primary-detail", -1, comparator,
				primaryLinkedA, primaryLinkedB);

		checkCompare("fixed", "unmapped-linked-mate-detail", -1, comparator,
				unmappedLinkedRoot("unmapped-a", 10, 1),
				unmappedLinkedRoot("unmapped-b", 1, 10));

		final Read paired=mappedRead("paired", 2, Shared.PLUS, 100);
		final Read unpaired=mappedRead("unpaired", 1, Shared.PLUS, 100);
		paired.setPaired(true);
		checkCompare("fixed", "primary-unlinked-paired-bit", -1, comparator,
				paired, unpaired);

		System.out.println("PASS\tfixed\ttie-preservation");
	}

	/** SRT-205 baseline for match-string null mixing and prefix equality.
	 * BBTools expanded and short match text can contain m/S/I/D; the comparator's
	 * secondary match path is reached by mapped reads tied at the same coordinate. */
	private static void checkMatchStringResidual(boolean fixed){
		final Read insertion=mappedRead("insertion", 1, Shared.PLUS, 100);
		final Read noMatch=mappedRead("no-match", 2, Shared.PLUS, 100);
		final Read deletion=mappedRead("deletion", 3, Shared.PLUS, 100);
		insertion.match=new byte[] {'I'};
		deletion.match=new byte[] {'D'};

		final ReadComparatorMapping comparator=new ReadComparatorMapping();
		if(fixed){
			checkOrder("fixed", "match-null", comparator, deletion, insertion, noMatch);
		}else{
			checkCompare("baseline", "match-null-insertion-nomatch", -1,
					comparator, insertion, noMatch);
			checkCompare("baseline", "match-null-nomatch-deletion", -1,
					comparator, noMatch, deletion);
			checkCompare("baseline", "match-null-deletion-insertion", -1,
					comparator, deletion, insertion);
		}

		final Read prefix=mappedRead("prefix", 4, Shared.PLUS, 100);
		final Read longer=mappedRead("prefix-longer", 5, Shared.PLUS, 100);
		prefix.match=new byte[] {'m'};
		longer.match=new byte[] {'m', 'S'};
		if(fixed){
			checkCompare("fixed", "match-prefix-length", -1, comparator, prefix, longer);
		}else{
			if(comparator.compareMatchStrings(prefix.match, longer.match)!=0){
				throw new AssertionError("Expected prefix-equal match strings to compare equal");
			}
			checkCompare("baseline", "match-prefix-later-id", -1, comparator, prefix, longer);
		}

		final Read insertionPrefix=mappedRead("prefix-insertion", 1, Shared.PLUS, 100);
		final Read noEventPrefix=mappedRead("prefix", 2, Shared.PLUS, 100);
		final Read deletionPrefix=mappedRead("prefix-deletion", 3, Shared.PLUS, 100);
		insertionPrefix.match=new byte[] {'m', 'I'};
		noEventPrefix.match=new byte[] {'m'};
		deletionPrefix.match=new byte[] {'m', 'D'};
		if(fixed){
			checkOrder("fixed", "match-prefix", comparator,
					noEventPrefix, deletionPrefix, insertionPrefix);
		}else{
			checkCompare("baseline", "match-prefix-insertion-prefix", -1,
					comparator, insertionPrefix, noEventPrefix);
			checkCompare("baseline", "match-prefix-prefix-deletion", -1,
					comparator, noEventPrefix, deletionPrefix);
			checkCompare("baseline", "match-prefix-deletion-insertion", -1,
					comparator, deletionPrefix, insertionPrefix);
		}

		System.out.println("PASS\t"+(fixed ? "fixed" : "baseline")+"\tmatch-string-residual");
	}

	/** Mapping fields stay at defaults; only mate links distinguish the roots. */
	private static Read read(String id, long numericID){
		return new Read(new byte[]{'A'}, null, id, numericID);
	}

	/** Builds a primary-plus/mate-minus mapped pair on the same scaffold. */
	private static Read mappedPairRoot(String id, int primaryStart, int mateStop,
			long numericID){
		return mappedPairRoot(id, Shared.PLUS, primaryStart, mateStop, numericID);
	}

	/** Builds an opposite-strand mapped pair on the same scaffold. */
	private static Read mappedPairRoot(String id, int primaryStrand, int primaryCoord,
			int mateCoord, long numericID){
		final int mateStrand=(primaryStrand==Shared.PLUS ? Shared.MINUS : Shared.PLUS);
		final Read primary=mappedRead(id+"/1", numericID, primaryStrand, primaryCoord);
		final Read mate=mappedRead(id+"/2", numericID, mateStrand, mateCoord);
		primary.mate=mate;
		mate.mate=primary;
		mate.setPairnum(1);
		primary.setPaired(true);
		mate.setPaired(true);
		return primary;
	}

	/** Builds an unmapped primary with a mapped second mate, preserving the root shape
	 * that triggers mate-vs-primary anchor switching. */
	private static Read mateMappedRoot(String id, int mateStop, long numericID){
		return mateMappedRoot(id, Shared.MINUS, mateStop, numericID);
	}

	/** Builds an unmapped primary with a mapped second mate at the requested anchor. */
	private static Read mateMappedRoot(String id, int mateStrand, int mateCoord,
			long numericID){
		final int primaryStrand=(mateStrand==Shared.PLUS ? Shared.MINUS : Shared.PLUS);
		final Read primary=read(id+"/1", numericID);
		primary.setStrand(primaryStrand);
		final Read mate=mappedRead(id+"/2", numericID, mateStrand, mateCoord);
		primary.mate=mate;
		mate.mate=primary;
		mate.setPairnum(1);
		return primary;
	}

	/** Builds a two-end mapped root whose mate is present but whose paired-alignment
	 * bit is false, matching the preexisting no-mate-tie-break path. */
	private static Read bothMappedUnpairedRoot(String id, int primaryStart, int mateStop,
			long numericID){
		return bothMappedUnpairedRoot(id, Shared.PLUS, primaryStart, mateStop, numericID);
	}

	/** Builds a mapped root with a mapped mate but no paired-alignment flag. */
	private static Read bothMappedUnpairedRoot(String id, int primaryStrand, int primaryCoord,
			int mateCoord, long numericID){
		final Read primary=mappedPairRoot(id, primaryStrand, primaryCoord, mateCoord, numericID);
		primary.setPaired(false);
		primary.mate.setPaired(false);
		return primary;
	}

	/** Builds a mapped primary with a present, unmapped second mate. */
	private static Read primaryMappedLinkedRoot(String id, int primaryStart, long numericID){
		return primaryMappedLinkedRoot(id, Shared.PLUS, primaryStart, numericID);
	}

	/** Builds a mapped primary at the requested anchor with a present, unmapped second mate. */
	private static Read primaryMappedLinkedRoot(String id, int primaryStrand, int primaryCoord,
			long numericID){
		final int mateStrand=(primaryStrand==Shared.PLUS ? Shared.MINUS : Shared.PLUS);
		final Read primary=mappedRead(id+"/1", numericID, primaryStrand, primaryCoord);
		final Read mate=read(id+"/2", numericID);
		mate.setStrand(mateStrand);
		primary.mate=mate;
		mate.mate=primary;
		mate.setPairnum(1);
		return primary;
	}

	/** Builds two linked unmapped reads on opposite strands. */
	private static Read unmappedLinkedRoot(String id, long numericID){
		return unmappedLinkedRoot(id, numericID, numericID);
	}

	/** Builds two linked unmapped reads on opposite strands with independent IDs. */
	private static Read unmappedLinkedRoot(String id, long primaryID, long mateID){
		final Read primary=read(id+"/1", primaryID);
		final Read mate=read(id+"/2", mateID);
		primary.setStrand(Shared.PLUS);
		mate.setStrand(Shared.MINUS);
		primary.mate=mate;
		mate.mate=primary;
		mate.setPairnum(1);
		return primary;
	}

	/** Builds one mapped read of length 1 whose plus start or minus stop is coord. */
	private static Read mappedRead(String id, long numericID, int strand, int coord){
		final Read r=read(id, numericID);
		r.setMapped(true);
		r.chrom=0;
		r.setStrand(strand);
		r.start=coord;
		r.stop=coord;
		return r;
	}

	/** Requires {@code first < second < third} and all reverse comparisons. */
	private static void checkOrder(String expect, String label, ReadComparatorMapping comparator,
			Read first, Read second, Read third){
		checkCompare(expect, label+"-first-second", -1, comparator, first, second);
		checkCompare(expect, label+"-second-third", -1, comparator, second, third);
		checkCompare(expect, label+"-first-third", -1, comparator, first, third);
		checkCompare(expect, label+"-second-first", 1, comparator, second, first);
		checkCompare(expect, label+"-third-second", 1, comparator, third, second);
		checkCompare(expect, label+"-third-first", 1, comparator, third, first);
		System.out.println("PASS\t"+expect+"\t"+label+"-order");
	}

	/** Requires a full ordered list, covering antisymmetry and transitivity edges. */
	private static void checkAllOrdered(String expect, String label,
			ReadComparatorMapping comparator, Read[] reads){
		for(int i=0; i<reads.length; i++){
			for(int j=0; j<reads.length; j++){
				checkCompare(expect, label+"-"+i+"-"+j, Integer.compare(i, j),
						comparator, reads[i], reads[j]);
			}
		}
		System.out.println("PASS\t"+expect+"\t"+label);
	}

	/** Checks one comparison sign. */
	private static void checkCompare(String expect, String label, int expected,
			ReadComparatorMapping comparator, Read a, Read b){
		final int observed=sign(comparator.compare(a, b));
		if(observed!=expected){
			throw new AssertionError(expect+" "+label+" expected "+expected+" observed "+observed);
		}
	}

	/** Returns the sign of a comparator result. */
	private static int sign(int value){return value<0 ? -1 : value>0 ? 1 : 0;}
}
