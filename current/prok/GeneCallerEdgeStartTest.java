package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import dna.AminoAcid;
import shared.Shared;

/** Regression tests for synthetic truncated starts at a contig's 5-prime edge. */
public class GeneCallerEdgeStartTest {

	public static void main(String[] args){
		int tests=0;
		tests+=checkCodonPolicy();
		tests+=checkEdgeStopSuppressed();
		tests+=checkAmbiguousEdgeSuppressed();
		tests+=checkSenseEdgeRetained();
		tests+=checkRealStartRetained();
		System.out.println("GeneCallerEdgeStartTest PASS: "+tests+" checks");
	}

	private static int checkCodonPolicy(){
		check(!GeneCaller.validSyntheticEdgeStart(codon("TAA")), "TAA accepted as an edge start");
		check(!GeneCaller.validSyntheticEdgeStart(codon("TAG")), "TAG accepted as an edge start");
		check(!GeneCaller.validSyntheticEdgeStart(codon("TGA")), "TGA accepted as an edge start");
		check(!GeneCaller.validSyntheticEdgeStart(codon("NNN")), "ambiguous codon accepted as an edge start");
		check(GeneCaller.validSyntheticEdgeStart(codon("AAA")), "sense codon rejected as a truncated edge start");
		check(GeneCaller.validSyntheticEdgeStart(codon("ATG")), "ATG rejected as an edge start");
		return 6;
	}

	private static int checkEdgeStopSuppressed(){
		int checks=0;
		for(String stop : new String[] {"TAA", "TAG", "TGA"}){
			for(int frame=0; frame<3; frame++){
				String prefix=(frame==0 ? "" : frame==1 ? "C" : "CC");
				ArrayList<Orf> orfs=orfs(prefix+stop+"ATGAAATAA", frame);
				check(orfs.size()==1, "expected one ORF after leading "+stop+
						" in frame "+frame+", found "+orfs.size());
				check(orfs.get(0).start==frame+3 && orfs.get(0).stop==frame+11,
						"leading "+stop+" was retained in frame "+frame+" coordinates: "+orfs.get(0));
				checks+=2;
			}
		}
		return checks;
	}

	private static int checkAmbiguousEdgeSuppressed(){
		ArrayList<Orf> orfs=orfs("NNNATGAAATAA");
		check(orfs.size()==1, "expected one ORF after leading NNN, found "+orfs.size());
		check(orfs.get(0).start==3 && orfs.get(0).stop==11,
				"leading ambiguous codon was retained in ORF coordinates: "+orfs.get(0));
		return 2;
	}

	private static int checkSenseEdgeRetained(){
		ArrayList<Orf> orfs=orfs("AAAAAATAA");
		check(orfs.size()==1, "expected one truncated edge ORF, found "+orfs.size());
		check(orfs.get(0).start==0 && orfs.get(0).stop==8,
				"valid truncated edge start was lost: "+orfs.get(0));
		return 2;
	}

	private static int checkRealStartRetained(){
		ArrayList<Orf> orfs=orfs("ATGAAATAA");
		check(orfs.size()==1, "expected one ATG-started ORF, found "+orfs.size());
		check(orfs.get(0).start==0 && orfs.get(0).stop==8,
				"ATG-started edge ORF was lost: "+orfs.get(0));
		return 2;
	}

	private static ArrayList<Orf> orfs(String bases){
		return orfs(bases, 0);
	}

	private static ArrayList<Orf> orfs(String bases, int frame){
		return GeneCaller.makeOrfsForFrame("fixture", bases.getBytes(StandardCharsets.US_ASCII),
				frame, Shared.PLUS, 9);
	}

	private static int codon(String bases){
		int code=0;
		for(byte b : bases.getBytes(StandardCharsets.US_ASCII)){
			code=(code<<2)|AminoAcid.baseToNumber[b];
		}
		return code;
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
