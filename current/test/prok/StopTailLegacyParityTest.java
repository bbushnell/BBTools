package test.prok;

import java.nio.charset.StandardCharsets;
import dna.GeneticCode;
import prok.StopTailCandidateScan;
import prok.StopTailCandidateScan.Candidate;
import prok.StopTailCandidateScan.Scan;

/** Emit deterministic existing-API observations under old and new classpaths.
 * The build recipe compares complete bytes, not just candidate counts.
 * @author Keqing
 */
public final class StopTailLegacyParityTest{
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		for(int table:new int[]{4, 11, 25}){
			for(int frame=0; frame<3; frame++){
				for(int length:new int[]{0, 39, 40, 59, 60, 79, 80, 99, 100, 124, 125, 149, 150, 199, 200, 201, 601}){
					final String text=repeat("N", frame)+repeat("CAA", length)+"TGA"
						+repeat("GCT", 100)+"NNN"+repeat("AAA", length)+"TAA";
					for(int strand=0; strand<2; strand++){
						final byte[] bases=(strand==0 ? text : reverse(text)).getBytes(StandardCharsets.US_ASCII);
						final GeneticCode code=GeneticCode.forTable(table);
						final String id=table+"_"+frame+"_"+length+"_"+strand;
						emit("default", StopTailCandidateScan.scan(id, bases, code));
						emit("bins", StopTailCandidateScan.scanBins(id, bases, code, true, false, true));
						for(int minimum:new int[]{100, 125, 150, 177, 200}){
							emit("minimum"+minimum, StopTailCandidateScan.scan(id, bases, code, minimum));
						}
					}
				}
			}
		}
	}
	private static void emit(String mode, Scan scan){
		final StringBuilder s=new StringBuilder(mode).append('\t').append(scan.contig).append('\t').append(scan.length)
			.append('\t').append(scan.under100).append('\t').append(scan.from100To199).append('\t').append(scan.atLeast200)
			.append('\t').append(scan.ambiguousCodons).append('\t').append(scan.unterminatedIntervals)
			.append('\t').append(scan.minimumTailLength).append('\t').append(scan.lengthBins)
			.append('\t').append(scan.shortQueries).append('\n');
		for(Candidate c:scan.candidates){
			s.append(c.key()).append('\t').append(c.stopLeft).append('\t').append(c.strand).append('\t').append(c.frame)
				.append('\t').append(c.upstreamCodons).append('\t')
				.append(new String(c.copyTail(), StandardCharsets.US_ASCII)).append('\n');
		}
		System.out.print(s);
	}
	private static String repeat(String unit, int n){
		final StringBuilder s=new StringBuilder();
		for(int i=0; i<n; i++){s.append(unit);}
		return s.toString();
	}
	private static String reverse(String dna){
		final StringBuilder s=new StringBuilder(dna.length());
		for(int i=dna.length()-1; i>=0; i--){
			final char c=dna.charAt(i);
			s.append(c=='A' ? 'T' : c=='T' ? 'A' : c=='C' ? 'G' : c=='G' ? 'C' : 'N');
		}
		return s.toString();
	}
}
