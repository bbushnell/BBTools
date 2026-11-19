package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import dna.GeneticCode;
import prok.StopTailCandidateScan;
import prok.StopTailCandidateScan.Candidate;

/** Literal synonymous coding-DNA and strand/view fixtures; no biological panel.
 * @author Qiqi
 */
public final class StopTailGcTest{
	public static void main(String[] args){
		final byte[] low=dna("GCT",100,"TAA"),high=dna("GCC",100,"TGA");
		final Candidate a=find(low,300,0),b=find(high,300,0);
		check(Arrays.equals(a.copyTail(),b.copyTail()),"Synonymous DNA has identical peptides");
		check(a.tailGc()==2.0/3 && b.tailGc()==1,"GC must come from actual DNA, excluding different stops");
		check(a.gcAvailable && b.gcAvailable,"Measured availability is explicit");
		final byte[] mixed=dna("AAA",100,"TGA");mixed[0]='G';mixed[1]='G';mixed[2]='G';
		final Candidate c=find(mixed,300,0);
		check(c.tailGc()==.01 && c.shortTailGc()==0,"LONG100 and SHORT99 have different actual intervals");
		final Candidate minus=find(reverseComplement(mixed),0,1);
		check(Arrays.equals(c.copyTail(),minus.copyTail()) && minus.tailGc()==.01 && minus.shortTailGc()==0,
			"Minus-strand coordinates and complement retain the same coding interval");
		final byte[] shorter=dna("AAA",40,"TGA");final Candidate zero=find(shorter,120,0);
		check(zero.gcAvailable && zero.tailGc()==0 && zero.shortTailGc()==0,"Measured all-AT is not missing");
		final Candidate ordinary=findOrdinary(mixed,300);
		check(!ordinary.gcAvailable && Arrays.equals(ordinary.copyTail(),c.copyTail()),
			"Default scan preserves peptide and explicitly lacks GC");
		shorter[15]='N';
		boolean failed=false;
		try{
			StopTailCandidateScan.senseGcCount(shorter,120,0,40);
		}catch(IllegalArgumentException expected){failed=true;}
		check(failed,"An ambiguous measured interval fails rather than changing its denominator");
		System.out.println("STOP_TAIL_GC_TEST_PASS checks="+checks);
	}
	private static Candidate find(byte[] dna,int stop,int strand){
		for(Candidate c:StopTailCandidateScan.scanCappedWithGc(
				"literal",dna,GeneticCode.forTable(11),40,false).candidates){
			if(c.stopLeft==stop && c.strand==strand){return c;}
		}
		throw new AssertionError("Missing expected stop");
	}
	private static Candidate findOrdinary(byte[] dna,int stop){
		for(Candidate c:StopTailCandidateScan.scanCapped(
				"literal",dna,GeneticCode.forTable(11),40,100,false).candidates){
			if(c.stopLeft==stop && c.strand==0){return c;}
		}
		throw new AssertionError("Missing ordinary stop");
	}
	private static byte[] dna(String codon,int n,String stop){
		final StringBuilder b=new StringBuilder();
		for(int i=0;i<n;i++){b.append(codon);}
		return b.append(stop).toString().getBytes(StandardCharsets.US_ASCII);
	}
	private static byte[] reverseComplement(byte[] dna){
		final byte[] out=new byte[dna.length];
		for(int i=0;i<dna.length;i++){
			final byte b=dna[i];
			out[dna.length-1-i]=(byte)(b=='A' ? 'T' : b=='T' ? 'A' : b=='G' ? 'C' : 'G');
		}
		return out;
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}checks++;}
	private static int checks;
}
