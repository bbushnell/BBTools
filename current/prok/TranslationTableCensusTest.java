package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPOutputStream;

/** Exact finite metadata oracles; no NCBI data or expected output derived from the census.
 * @author Keqing
 */
public final class TranslationTableCensusTest {

	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected a new fixture directory");}
		final Path dir=Paths.get(args[0]);
		Files.createDirectory(dir);
		final Path gff=dir.resolve("mixed.gff");
		final String content="##gff-version 3\n"
			+row("CDS", "ID=a;transl_table=4")+row("CDS", "transl_table=25;ID=b")
			+row("CDS", "Note=transl_table=11")+row("tRNA", "transl_table=999")
			+row("CDS", "transl_table=4;transl_except=(pos:3..5%2Caa:Sec)")
			+"##FASTA\n>sequence\nATGTAA\n";
		Files.write(gff, content.getBytes(StandardCharsets.US_ASCII));
		final TranslationTableCensus.Counts counts=TranslationTableCensus.scan(gff.toFile());
		check(counts.total==4 && counts.missing==1 && counts.exceptions==1, "CDS/missing/exception counts must be 4/1/1");
		check(counts.codes.size==3 && counts.codes.get(0)==4 && counts.rows.get(0)==2
			&& counts.codes.get(1)==25 && counts.rows.get(1)==1 && counts.codes.get(2)==0 && counts.rows.get(2)==1,
			"Expected exact table histogram {4:2,25:1,missing:1}; non-CDS and Note substrings must not count");
		check(TranslationTableCensus.pairedFasta(gff.toFile())==null, "Absent pair must remain absent");
		final Path fasta=dir.resolve("mixed.fna.gz");
		try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(fasta))){out.write(">sequence\nATGTAA\n".getBytes(StandardCharsets.US_ASCII));}
		check(TranslationTableCensus.pairedFasta(gff.toFile()).equals(fasta.toFile()), "Compressed exact pair must resolve");
		Files.write(dir.resolve("mixed.fa"), ">s\nATG\n".getBytes(StandardCharsets.US_ASCII));
		reject(()->TranslationTableCensus.pairedFasta(gff.toFile()), "Ambiguous FASTA pair");
		Files.move(dir.resolve("mixed.fa"), dir.resolve("not_a_pair.fa"));
		for(String invalid : new String[]{"transl_table=0", "transl_table=-1", "transl_table=", "transl_table=x",
			"transl_table=2147483648", "transl_table=4;transl_table=4", "transl_table=4;transl_table=25"}){
			final byte[] bytes=invalid.getBytes(StandardCharsets.US_ASCII);
			reject(()->new TranslationTableCensus.Counts().add(bytes, 0, bytes.length), "transl_table");
		}
		final Path compressed=dir.resolve("uniform.gff.gz");
		try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(compressed))){
			out.write((row("CDS", "transl_table=11")+row("CDS", "transl_table=11")).getBytes(StandardCharsets.US_ASCII));
		}
		final TranslationTableCensus.Counts compressedCounts=TranslationTableCensus.scan(compressed.toFile());
		check(compressedCounts.total==2 && compressedCounts.missing==0 && compressedCounts.codes.get(0)==11,
			"Gzip input must preserve declared table11 counts");
		Files.write(dir.resolve("empty.gff"), row("tRNA", "ID=rna").getBytes(StandardCharsets.US_ASCII));
		System.out.println("PASS TranslationTableCensusTest: exact CDS metadata, gzip, missing/mixed codes, pairs and rejection diagnostics");
	}

	private static String row(String type, String attributes){return "sequence\tRefSeq\t"+type+"\t1\t6\t.\t+\t0\t"+attributes+"\n";}
	private static void reject(Runnable action, String message){
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(message), "Wrong failure: "+e); return;}
		throw new AssertionError("Expected rejection: "+message);
	}
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}
}
