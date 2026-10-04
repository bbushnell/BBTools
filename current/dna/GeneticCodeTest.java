package dna;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

/** Synthetic genetic-code checks: biology, native indexing, isolation, malformed input.
 * @author Keqing
 */
public final class GeneticCodeTest {

	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("GeneticCodeTest takes no arguments");}
		final byte[] legacy=AminoAcid.codeToByte.clone();
		testCodons();
		testInitiation();
		testCustom();
		check(Arrays.equals(legacy, AminoAcid.codeToByte), "GeneticCode must not mutate global AminoAcid translation");
		System.out.println("PASS GeneticCodeTest: 64-codon parity, tables 4/11/25, initiation, custom TSV and isolation");
	}

	private static void testCodons(){
		final GeneticCode bacterial=GeneticCode.forTable(11), four=GeneticCode.forTable(4), twentyFive=GeneticCode.forTable(25);
		int stops11=0, stops4=0, stops25=0;
		for(int i=0; i<64; i++){
			final String triplet=AminoAcid.codonToString(i);
			check(GeneticCode.codon(triplet)==i, "Native codon order mismatch for "+triplet);
			check(GeneticCode.codon(triplet.toLowerCase(java.util.Locale.ROOT))==i, "Lowercase codon mismatch");
			final byte expected=(byte)AminoAcid.toChar(triplet);
			check(bacterial.aminoAcid(i)==expected, "Table11 elongation mismatch for "+triplet);
			check(four.aminoAcid(i)==(triplet.equals("TGA") ? 'W' : expected), "Table4 reassignment mismatch for "+triplet);
			check(twentyFive.aminoAcid(i)==(triplet.equals("TGA") ? 'G' : expected), "Table25 reassignment mismatch for "+triplet);
			stops11+=bacterial.isStop(i) ? 1 : 0;
			stops4+=four.isStop(i) ? 1 : 0;
			stops25+=twentyFive.isStop(i) ? 1 : 0;
		}
		check(stops11==3 && stops4==2 && stops25==2, "Stop counts must reflect TGA reassignment");
		check(GeneticCode.codon("AAN")==-1 && GeneticCode.codon("AUG")==-1, "Ambiguous/RNA keys must not be DNA codons");
		check(GeneticCode.codon("A\u0080G")==-1, "Non-ASCII codons must not index outside base lookup");
		check(GeneticCode.codon(new byte[]{'A', -1, 'G'}, 0)==-1, "Signed input bytes must not index outside base lookup");
		check(bacterial.aminoAcid(-1)=='X' && !bacterial.isStart(-1) && !bacterial.isStop(-1), "Unknown codon must remain unknown");
		reject("Unsupported translation table", ()->GeneticCode.forTable(999));
		reject("Codon index", ()->bacterial.aminoAcid(64));
		reject("Codon index", ()->bacterial.isStop(-2));
		reject("exactly three", ()->GeneticCode.codon("AT"));
		reject("three bases", ()->GeneticCode.codon(new byte[]{'A', 'T'}, 0));
		reject("three bases", ()->GeneticCode.codon(new byte[]{'A', 'T', 'G'}, Integer.MAX_VALUE));
	}

	private static void testInitiation(){
		final int gtg=GeneticCode.codon("GTG"), ctg=GeneticCode.codon("CTG");
		final GeneticCode eleven=GeneticCode.forTable(11);
		check(eleven.translateCodon(gtg, true)=='M', "Complete GTG initiation must produce M");
		check(eleven.translateCodon(gtg, false)=='V', "Internal/edge GTG must retain V");
		check(eleven.translateCodon(ctg, true)=='M' && eleven.translateCodon(ctg, false)=='L', "CTG initiation and elongation differ");
		final String[][] starts={{"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"},
			{"TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}, {"TTG", "ATG", "GTG"}};
		final int[] tables={11, 4, 25};
		for(int t=0; t<tables.length; t++){
			final GeneticCode code=GeneticCode.forTable(tables[t]);
			for(int i=0; i<64; i++){
				final boolean expected=Arrays.asList(starts[t]).contains(AminoAcid.codonToString(i));
				check(code.isStart(i)==expected, "Unexpected initiation set for table "+tables[t]+" codon "+i);
			}
		}
		reject("Cannot initiate", ()->eleven.translateCodon(GeneticCode.codon("AGT"), true));
		reject("Cannot initiate", ()->eleven.translateCodon(-1, true));
	}

	private static void testCustom() throws Exception{
		final Path file=Files.createTempFile("genetic-code-", ".tsv");
		try{
			final String source=customTable();
			Files.write(file, source.getBytes(StandardCharsets.US_ASCII));
			final GeneticCode custom=GeneticCode.load(file.toString());
			check(custom.id()==0, "Custom code must not claim an NCBI table ID");
			check(custom.isStart(GeneticCode.codon("AGT")), "Custom serine initiation permission lost");
			check(custom.aminoAcid(GeneticCode.codon("AGT"))=='S', "Custom start must not change internal serine");
			check(custom.translateCodon(GeneticCode.codon("AGT"), true)=='M', "Custom initiation must produce M");
			check(!GeneticCode.forTable(11).isStart(GeneticCode.codon("AGT")), "Custom table changed shared built-in starts");
			check(custom.aminoAcid(GeneticCode.codon("TGA"))=='W', "Custom table lost TGA=W");
			final boolean oldSimd=shared.Shared.SIMD;
			try{
				for(boolean vector : new boolean[]{false, simd.Vector.simd256}){
					shared.Shared.SIMD=vector;
					for(String text : new String[]{source.replace("\n", "\r\n"), source.substring(0, source.length()-1)}){
						Files.write(file, text.getBytes(StandardCharsets.US_ASCII));
						final GeneticCode reread=GeneticCode.load(file.toString());
						for(int i=0; i<64; i++){
							check(reread.aminoAcid(i)==custom.aminoAcid(i) && reread.isStart(i)==custom.isStart(i),
								"CRLF/unterminated input changed codon "+i+", SIMD="+vector);
						}
					}
				}
			}finally{shared.Shared.SIMD=oldSimd;}
			Files.write(file, "destroyed after load".getBytes(StandardCharsets.US_ASCII));
			check(custom.aminoAcid(GeneticCode.codon("TGA"))=='W', "Loaded code must not depend on later file contents");
			badFile(file, source+"AAA\tK\t0\n", "duplicate codon");
			badFile(file, source.replace("AAA\tK\t0\n", ""), "64 unique codons");
			badFile(file, source.replace("AAA\tK\t0", "AAN\tK\t0"), "only A/C/G/T");
			badFile(file, source.replace("TTT\tF\t0", "UUU\tF\t0"), "only A/C/G/T");
			badFile(file, source.replace("AAA\tK\t0", "AAA\tX\t0"), "canonical uppercase");
			badFile(file, source.replace("AAA\tK\t0", "AAA\tK\t2"), "start must be");
			badFile(file, source.replace("TAA\t*\t0", "TAA\t*\t1"), "cannot also initiate");
			badFile(file, source.replace("AAA\tK\t0", "AAA\tK\t0\textra"), "separated by tabs");
			badFile(file, source.replace("codon\tamino_acid\tstart", "wrong\theader"), "expected header");
			badFile(file, "", "64 unique codons");
			badFile(file, source.replace("\t1\n", "\t0\n"), "needs starts and stops");
			badFile(file, source.replace("\t*\t", "\tW\t"), "needs starts and stops");
		}finally{Files.delete(file);}
	}

	/** Hand-specified difference from table11: TGA=W, and AGT can initiate. */
	private static String customTable(){
		final StringBuilder out=new StringBuilder("# Synthetic test code\ncodon\tamino_acid\tstart\n");
		for(int i=0; i<64; i++){
			final String codon=AminoAcid.codonToString(i);
			final char aa=(codon.equals("TGA") ? 'W' : AminoAcid.toChar(codon));
			out.append(codon).append('\t').append(aa).append('\t');
			out.append(codon.equals("ATG") || codon.equals("AGT") ? '1' : '0').append('\n');
		}
		return out.toString();
	}

	private static void badFile(final Path path, final String content, final String diagnostic) throws Exception{
		Files.write(path, content.getBytes(StandardCharsets.US_ASCII));
		reject(diagnostic, ()->GeneticCode.load(path.toString()));
	}

	private static void reject(final String diagnostic, final Runnable action){
		try{action.run();}
		catch(IllegalArgumentException e){
			check(e.getMessage()!=null && e.getMessage().contains(diagnostic), "Wrong failure reason: "+e);
			return;
		}
		throw new AssertionError("Expected diagnostic: "+diagnostic);
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
