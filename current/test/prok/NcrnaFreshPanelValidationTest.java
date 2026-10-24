package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaFreshPanelValidation;
import prok.RrnaResourceIO;

/** Native FASTA observation, assembly-scoped lookup, complete GFF evidence and bad-binding fixtures.
 * @author Brian Bushnell, Raiden */
public final class NcrnaFreshPanelValidationTest {
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("outdir="), "Use outdir=fresh_directory");
		final String d=args[0].substring(7);check(new File(d).mkdir(), "Fixture directory must be fresh");
		final String fa1=d+"/a.fa", fa2=d+"/b.fa", gff1=d+"/a.gff", gff2=d+"/b.gff";
		final String fasta=">shared first assembly\nACGTACGT\n>tid|11|tagged native identity\nACGT\n";
		write(fa1, fasta);write(fa2, ">shared second assembly\nACGT\n");
		final String valid="shared\tRefSeq\tgene\t1\t8\t.\t+\t.\tID=x;Note=spaces and escaped%3Bsemicolon";
		final String[] rows={valid, valid,
			"shared\tRefSeq\tRNA\t8\t9\t.\t+\t.\tID=range",
			"shared\tRefSeq\tRNA\tbad\t4\t.\t+\t.\tID=bad",
			"absent\tRefSeq\tRNA\t1\t4\t.\t-\t.\tID=missing",
			"tagged\tRefSeq\ttRNA\t1\t4\t.\t-\t.\tID=alias"};
		final StringBuilder original=new StringBuilder("##gff-version 3\n\n");
		for(String row:rows){original.append(row).append('\n');}original.append("##FASTA\n>shared\nACGT\n");
		write(gff1, original.toString());
		write(gff2, "shared\tRefSeq\tRNA\t1\t5\t.\t+\t.\tID=assembly_specific_range\n");
		final String selection=NcrnaFreshPanelValidation.SELECTION_HEADER+'\n'+selected("GCF_1.1", 11, "bacteria", 12)+selected("GCF_2.1", 22, "archaea", 4);
		final String first=download(1, "GCF_1.1", fa1, gff1), second=download(2, "GCF_2.1", fa2, gff2);
		write(d+"/selection.tsv", selection);write(d+"/downloads.tsv", NcrnaFreshPanelValidation.DOWNLOAD_HEADER+'\n'+first+second);
		final String bindings=HASH+'\t'+d+"/selection.tsv\n"+HASH+'\t'+d+"/downloads.tsv\n"+HASH+'\t'+fa1+'\n'+HASH+'\t'+fa2+'\n'+HASH+'\t'+gff1+'\n'+HASH+'\t'+gff2+'\n';
		write(d+"/sources.tsv", bindings);run(d, "good");
		final ArrayList<byte[]> features=read(d+"/good.features.tsv.gz");check(features.size()==8, "Header plus every source feature, including duplicates, must survive");
		final String[] statuses={"VALID_INTERVAL", "VALID_INTERVAL", "OUT_OF_RANGE", "INVALID_COORDINATES", "UNRESOLVED_SEQID", "VALID_INTERVAL"};
		final LineParser1 p=new LineParser1('\t');
		for(int i=0; i<rows.length; i++){
			p.set(features.get(i+1));check(p.terms()==23 && p.termEquals("bacteria", 0) && p.termEquals("GCF_1.1", 1)
				&& p.termEquals(statuses[i], 5) && p.parseLong(4)==i+3, "Every original row keeps its assembly, status and source line");
			check(new String(features.get(i+1), StandardCharsets.UTF_8).endsWith('\t'+rows[i]), "The nine original GFF fields must remain byte-exact");
		}
		p.set(features.get(6));check(p.termEquals("tid|11|tagged", 8) && p.termEquals("PROKARYOTIC", 10), "Native-tag alias resolves without calling a prokaryote nuclear");
		p.set(features.get(7));check(p.termEquals("archaea", 0) && p.termEquals("GCF_2.1", 1)
			&& p.termEquals("OUT_OF_RANGE", 5) && p.parseLong(9)==4, "Shared seqid must use the second assembly's own four-base contig");
		final ArrayList<byte[]> contigs=read(d+"/good.contigs.tsv");check(contigs.size()==4, "Observe every actual contig");
		p.set(contigs.get(1));check(p.parseLong(3)==8 && p.parseInt(5)==0 && p.parseInt(6)==11 && p.termEquals("NOT_TAGGED", 7), "Catalog TaxID must not be invented as a native FASTA tag");
		final ArrayList<byte[]> summary=read(d+"/good.summary.tsv");check(summary.size()==3, "One summary per selected assembly");
		p.set(summary.get(1));check(p.parseLong(5)==11 && p.parseLong(6)==2 && p.parseLong(7)==1 && p.parseLong(8)==2
			&& p.parseLong(9)==6 && p.parseLong(10)==3 && p.parseLong(11)==1 && p.parseLong(12)==1 && p.parseLong(13)==1,
			"Source comments, blank lines, embedded FASTA, duplicate features and every status must conserve");
		write(fa1, ">shared\nACGTACGT\n>shared\nACGT\n");bad(d, "duplicate", "duplicate FASTA contig");
		write(fa1, ">tagged\nACGTACGT\n>tid|11|tagged\nACGT\n");bad(d, "ambiguous", "Ambiguous normalized FASTA accession");
		write(fa1, fasta.replace("ACGTACGT", "ACGTACG"));bad(d, "length", "Observed FASTA bases differ");write(fa1, fasta);
		write(d+"/downloads.tsv", NcrnaFreshPanelValidation.DOWNLOAD_HEADER+'\n'+first.replace("GCF_1.1", "GCF_2.1")+second);bad(d, "swapped", "assembly/index disagreement");
		write(d+"/downloads.tsv", NcrnaFreshPanelValidation.DOWNLOAD_HEADER+'\n'+first+first);bad(d, "duplicate_index", "Duplicate index");
		write(d+"/downloads.tsv", NcrnaFreshPanelValidation.DOWNLOAD_HEADER+'\n'+first);bad(d, "missing_pair", "exactly one downloaded pair");
		write(d+"/downloads.tsv", NcrnaFreshPanelValidation.DOWNLOAD_HEADER+'\n'+first+second);
		write(d+"/sources.tsv", bindings.replace(HASH+'\t'+fa1, "11111111111111111111\t"+fa1));bad(d, "hash", "conflicting verified source binding");
		write(d+"/sources.tsv", bindings);
		System.out.println("FRESH_PANEL_VALIDATION_TEST_PASS native_streaming assembly_scoped_seqids original_rows_preserved duplicate_and_length_rejection frozen_index_hash_bindings no_nuclear_relabeling");
	}
	private static String selected(String id, int tid, String domain, int bases){
		return id+"\tna\t"+tid+'\t'+tid+'\t'+tid+"\t2\tfixture_phylum\t"+domain+"\t2026-10-01\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/fixture/\tfixture\t"+bases+"\tKNOWN_SOURCES_EXCLUDED_NEEDS_PROVENANCE\n";
	}
	private static String download(int i, String id, String fasta, String gff){return i+"\t"+id+'\t'+fasta+'\t'+HASH+'\t'+gff+'\t'+HASH+'\t'+HASH+'\n';}
	private static void run(String d, String label){
		NcrnaFreshPanelValidation.main(new String[]{"selection="+d+"/selection.tsv", "downloads="+d+"/downloads.tsv", "sources="+d+"/sources.tsv", "expected=2", "outprefix="+d+'/'+label});
	}
	private static void bad(String d, String label, String reason){
		boolean rejected=false;
		try{run(d, label);}catch(IllegalArgumentException e){check(e.getMessage().contains(reason), "Wrong failure for "+label+": "+e);rejected=true;}
		check(rejected && !new File(d+'/'+label+".summary.tsv").exists(), "Binding/FASTA failure must not publish a successful GFF summary: "+label);
	}
	private static ArrayList<byte[]> read(String path){
		final ArrayList<byte[]> rows=new ArrayList<byte[]>();final ByteFile in=RrnaResourceIO.open(path);
		try{for(byte[] line; (line=in.nextLine())!=null; ){rows.add(line);}}finally{check(!in.close(), "Fixture output read failed");}return rows;
	}
	private static void write(String path, String text){
		// Fixture inputs are deliberately changed between negative cases; product outputs stay fresh.
		final ByteStreamWriter out=new ByteStreamWriter(path, true, false, false);out.start();out.print(text);RrnaResourceIO.finish(out);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String HASH="00000000000000000000";
}
