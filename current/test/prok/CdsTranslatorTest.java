package test.prok;

import prok.CdsTranslator;
import prok.Orf;
import prok.ProkObject;

import java.nio.charset.StandardCharsets;

import dna.AminoAcid;
import gff.GffLine;

/** Independent literal translation oracles for the AnnotateGenes input contract.
 * @author Keqing
 */
public final class CdsTranslatorTest {

	public static void main(String[] args){
		peptide("ATGGTTTAA", "+", 0, ".", 11, "MV");
		peptide("TTAAACCAT", "-", 0, ".", 11, "MV");
		peptide("GTGAAATAA", "+", 0, ".", 11, "MK");
		peptide("TTGAAATAA", "+", 0, ".", 11, "MK");
		peptide("GTGAAATAA", "+", 0, "partial=10", 11, "VK");
		peptide("TTATTTCAC", "-", 0, "partial=01", 11, "VK");
		peptide("GTGAAAA", "+", 0, "partial=01", 11, "MK");
		peptide("TTTTCAC", "-", 0, "partial=10", 11, "MK");
		peptide("GTGAAAA", "+", 0, "partial=11", 11, "VK");
		peptide("GTGAAAA", "+", 0, "partial=true", 11, "VK");
		peptide("AGTGAAATAA", "+", 1, "partial=10", 11, "VK");
		peptide("TTATTTCACT", "-", 1, "partial=01", 11, "VK");
		peptide("CCGTGAAATAA", "+", 2, "partial=10", 11, "VK");
		peptide("TTATTTCACGG", "-", 2, "partial=01", 11, "VK");
		peptide("ATGTGATAA", "+", 0, ".", 4, "MW");
		peptide("ATGTGATAA", "+", 0, ".", 25, "MG");
		peptide("ATGTGATAA", "+", 0, "transl_table=4", 11, "MW");
		unsupported("ATGTGATAA", "+", 0, ".", 11, "internal_stop");
		unsupported("ATGAAATGA", "+", 0, ".", 4, "missing_terminal_stop");
		unsupported("ATGAAATGA", "+", 0, ".", 25, "missing_terminal_stop");
		unsupported("ATGNNNTAA", "+", 0, ".", 11, "ambiguous_codon");
		unsupported("AAAAAATAA", "+", 0, ".", 11, "invalid_start_codon");
		unsupported("ATGAAATAA", "+", 0, "transl_except=(pos:4..6%2Caa:Sec)", 11, "translation_exception");
		unsupported("ATGAAATAA", "+", -1, ".", 11, "unknown_phase");
		unsupported("ATGAAATAA", "+", 0, "transl_table=2", 11, "unsupported_translation_table");
		final GffLine offset=line(3, 11, "+", 0, "ID=x;Note=partial=true;partial=00;anticodon=TCC");
		final String attributes=offset.attributes;
		final CdsTranslator.Translation result=CdsTranslator.translate(offset, bytes("CCATGGTTTAAGG"), 11);
		check("MV".equals(new String(result.peptide, StandardCharsets.US_ASCII)), "Offset interval oracle");
		check(attributes.equals(offset.attributes) && offset.start==3 && offset.stop==11, "Translation must not mutate GFF metadata");
		boolean failed=false;
		try{CdsTranslator.translate(line(1, 10, "+", 0, "."), bytes("ATGAAATAA"), 11);}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Out-of-bounds intervals must fail before extraction");
		callerPartial("ATGAAATAA", GffLine.PLUS, false, false, null, "MK");
		callerPartial("TTATTTCAT", GffLine.MINUS, false, false, null, "MK");
		callerPartial("GTGAAATAA", GffLine.PLUS, true, false, "10", "VK");
		callerPartial("TTATTTCAC", GffLine.MINUS, true, false, "01", "VK");
		callerPartial("ATGAAAC", GffLine.PLUS, false, true, "01", "MK");
		callerPartial("GTTTCAT", GffLine.MINUS, false, true, "10", "MK");
		callerPartial("GTGAAAC", GffLine.PLUS, true, true, "11", "VK");
		callerPartial("GTTTCAC", GffLine.MINUS, true, true, "11", "VK");
		final Orf rna=new Orf("ctg", 0, 8, GffLine.PLUS, 0, bytes("ATGAAATAA"), false, ProkObject.r5S, true, true);
		check(!rna.toGff().contains(";partial="), "CDS truncation flags alone must not invent RNA clipping provenance");
		System.err.println("CDS_TRANSLATION_PASS literal_peptides=18 unsupported=8 preservation=true bounds=true caller_gff_cases=9");
	}

	/** Literal peptides/edge bits exercise the actual producer-to-consumer handoff. */
	private static void callerPartial(String dna, int strand, boolean five, boolean three, String expectedBits, String expectedPeptide){
		final byte[] oriented=bytes(dna);
		if(strand==GffLine.MINUS){AminoAcid.reverseComplementBasesInPlace(oriented);}
		// A terminal frame-zero candidate ends at the last complete codon, even
		// when the physical scaffold has one or two remaining bases.
		final int stop=three ? oriented.length-1-oriented.length%3 : oriented.length-1;
		final Orf orf=new Orf("ctg", 0, stop, strand, 0, oriented, true, ProkObject.CDS, five, three);
		final String annotated=orf.toGff();
		final int partial=annotated.indexOf(";partial=");
		check(expectedBits==null ? partial<0 : partial>=0 && annotated.endsWith(";partial="+expectedBits)
			&& partial==annotated.lastIndexOf(";partial="), "Caller must emit exactly the expected genomic edge bits");
		final CdsTranslator.Translation translated=CdsTranslator.translate(new GffLine(bytes(annotated)), bytes(dna), 11);
		check(translated.reason==null && expectedPeptide.equals(new String(translated.peptide, StandardCharsets.US_ASCII)),
			"Caller partial metadata must yield the independently specified peptide on both strands; reason="+translated.reason);
	}

	private static void peptide(String dna, String strand, int phase, String attrs, int table, String expected){
		final CdsTranslator.Translation result=CdsTranslator.translate(line(1, dna.length(), strand, phase, attrs), bytes(dna), table);
		check(result.reason==null && result.peptide!=null, "Expected peptide for "+dna+" "+strand+" "+attrs+", got "+result.reason);
		check(expected.equals(new String(result.peptide, StandardCharsets.US_ASCII)), "Literal peptide differs for "+dna+" "+strand+" "+attrs);
	}
	private static void unsupported(String dna, String strand, int phase, String attrs, int table, String reason){
		final CdsTranslator.Translation result=CdsTranslator.translate(line(1, dna.length(), strand, phase, attrs), bytes(dna), table);
		check(result.peptide==null && reason.equals(result.reason), "Expected "+reason+" for "+dna+", got "+result.reason);
	}
	private static GffLine line(int start, int end, String strand, int phase, String attrs){
		return new GffLine(bytes("ctg\tfixture\tCDS\t"+start+"\t"+end+"\t.\t"+strand+"\t"+(phase<0 ? "." : phase)+"\t"+attrs));
	}
	private static byte[] bytes(String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}
}
