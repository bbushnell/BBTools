package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.GeneticCode;
import prok.StopTailCandidateScan;
import prok.StopTailCandidateScan.Candidate;
import prok.StopTailCandidateScan.Scan;

/** Literal first-start query oracles; never invokes a classifier or infers expected residues from it.
 * @author Keqing
 */
public final class StopTailFirstStartTest{

	public static void main(String[] args){
		check(args.length==0, "This scanner fixture needs no files");
		for(int table:new int[]{4,11,25}){for(int frame=0; frame<3; frame++){
			for(int prefix:new int[]{0,17,217}){for(int length:new int[]{39,40,99,100,124,125,149,150,199,200,201,601}){
				final String text=repeat("N", frame)+repeat("CAA", prefix)+"ATG"+repeat("AAA", length-1)+"TAA";
				final String peptide="M"+repeat("K", length-1);
				both(text, table, frame, prefix+length, length,
					length<40 ? null : peptide.substring(Math.max(0, length-200)),
					StopTailCandidateScan.BOUNDARY_EDGE);
			}}
		}}
		// The later ATG must not replace the first; the full80 residues remain in the query.
		both(repeat("CAA",17)+"ATG"+repeat("CAA",39)+"ATG"+repeat("AAA",39)+"TAA",11,0,97,80,
			"M"+repeat("Q",39)+"M"+repeat("K",39),StopTailCandidateScan.BOUNDARY_EDGE);
		for(String reset:new String[]{"TAA","NNN"}){
			final String text="ATG"+repeat("CAA",216)+reset+repeat("CAA",17)+"ATG"+repeat("AAA",39)+"TAA";
			both(text,11,0,57,40,"M"+repeat("K",39),
				reset.equals("TAA") ? StopTailCandidateScan.BOUNDARY_STOP : StopTailCandidateScan.BOUNDARY_AMBIGUITY);
			both("ATG"+repeat("CAA",216)+reset+repeat("CAA",80)+"TAA",11,0,80,0,null,
				reset.equals("TAA") ? StopTailCandidateScan.BOUNDARY_STOP : StopTailCandidateScan.BOUNDARY_AMBIGUITY);
		}
		// CTG can initiate in table11 but not25; elongation translation remains L.
		both("CTG"+repeat("CAA",79)+"TAA",11,0,80,80,"L"+repeat("Q",79),StopTailCandidateScan.BOUNDARY_EDGE);
		both("CTG"+repeat("CAA",79)+"TAA",25,0,80,0,null,StopTailCandidateScan.BOUNDARY_EDGE);
		for(int n:new int[]{39,40,200,601}){
			final Scan absent=scan(repeat("CAA",n)+"TAA",11,true);
			check(absent.candidates.isEmpty() && absent.noStartStops==(n>=40 ? 1 : 0),
				"No-start losses count otherwise eligible real stops, excluding the sub40 interval");
		}
		final Scan gap=scan("NNNATG"+repeat("CAA",39)+"TAA",11,true);
		check(gap.ambiguityBoundedQueries==1 && gap.edgeBoundedQueries==0,
			"An explicit ambiguous codon bounds the candidate separately from a contig edge");
		final String disabled=repeat("CAA",200)+"ATG"+repeat("CAA",124)+"TAA";
		final Scan bins=StopTailCandidateScan.scanBins(
			"disabled",bytes(disabled),GeneticCode.forTable(11),true,true,false,true,40,true);
		check(target(bins,disabled.length()-3,0)==null,
			"Actual125-aa query must respect a disabled125–149 bin despite raw325 interval");
		check(scan("ATG"+repeat("CAA",300),11,true).candidates.isEmpty(),"A contig end is never a fabricated stop");
		System.out.println("STOP_TAIL_FIRST_START_TEST_PASS checks="+checks);
	}

	private static void both(String text,int table,int frame,int upstream,int querySpan,String expected,int boundary){
		for(int strand=0; strand<2; strand++){
			final String oriented=strand==0 ? text : reverseComplement(text);
			final byte[] dna=bytes(oriented), saved=dna.clone();
			final GeneticCode code=GeneticCode.forTable(table);
			final Scan actual=StopTailCandidateScan.scanBins("fixture",dna,code,true,true,true,true,40,true);
			final int left=strand==0 ? text.length()-3 : 0;
			final Candidate found=target(actual,left,strand);
			check((found!=null)==(expected!=null),"First-start eligibility must match the literal selected span");
			if(found!=null){
				check(found.frame==frame && found.upstreamCodons==upstream && found.queryCodons==querySpan &&
					found.intervalBoundary==boundary && Arrays.equals(found.copyTail(),bytes(expected)),
					"Earliest start, selected span, genomic frame, boundary and exact capped residues must agree");
			}
			check(Arrays.equals(saved,dna),"Scanning must never alter the nucleotide input");
			final Scan old=StopTailCandidateScan.scanBins("fixture",dna,code,true,true,true,true,40);
			final Scan explicitOff=StopTailCandidateScan.scanBins("fixture",dna,code,true,true,true,true,40,false);
			check(old.candidates.size()==explicitOff.candidates.size() && old.realStops()==actual.realStops(),
				"Optional start filtering preserves physical-stop census and default candidate population");
			for(int i=0; i<old.candidates.size(); i++){
				final Candidate a=old.candidates.get(i),b=explicitOff.candidates.get(i);
				check(a.key()==b.key() && a.upstreamCodons==b.upstreamCodons && a.queryCodons==a.upstreamCodons &&
					Arrays.equals(a.copyTail(),b.copyTail()),"Explicit OFF must preserve default query bytes and keys");
			}
		}
	}
	private static Candidate target(Scan scan,int left,int strand){
		Candidate found=null;
		for(Candidate candidate:scan.candidates){if(candidate.stopLeft==left && candidate.strand==strand){
			check(found==null,"Every physical stop/strand has at most one query"); found=candidate;
		}}
		return found;
	}
	private static Scan scan(String dna,int table,boolean first){
		return StopTailCandidateScan.scanBins("fixture",bytes(dna),GeneticCode.forTable(table),true,true,true,true,40,first);
	}
	private static byte[] bytes(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String repeat(String text,int count){
		final StringBuilder out=new StringBuilder(); for(int i=0; i<count; i++){out.append(text);} return out.toString();
	}
	private static String reverseComplement(String text){
		final StringBuilder out=new StringBuilder();
		for(int i=text.length()-1; i>=0; i--){
			final char c=text.charAt(i);
			out.append(c=='A' ? 'T' : c=='T' ? 'A' : c=='C' ? 'G' : c=='G' ? 'C' : 'N');
		}
		return out.toString();
	}
	private static void check(boolean value,String message){checks++;if(!value){throw new AssertionError(message);}}
	private static int checks;
}
