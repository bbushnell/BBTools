package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CobalaminNativePanels;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Independent tiny metadata oracle: a negative input leads the permutation,
 * repeated contig/coordinate names belong to distinct inputs, and a truncated
 * main-family locus must enter protection without becoming a construction input.
 * FASTA paths deliberately do not exist, proving this adapter does not read DNA.
 * @author Raiden
 */
public final class CobalaminNativePanelsTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CobalaminNativePanelsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture option: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null, "Supply out=fresh_fixture_directory");
			check(new File(parser.out1).mkdir(), "Fixture directory must be fresh");
			new CobalaminNativePanelsTest(parser.out1).run();
			pp.outstream.println("COBALAMIN_NATIVE_PANELS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private CobalaminNativePanelsTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		write("order.tsv", "assembly_id\tinput_id\tchunk\tadmitted_loci\tlinked_fasta\tsource_fasta\n"
			+"g7\t7\t0\t0\t/absent/g7.fa\t/absent/g7.fa\n"
			+"g2\t2\t0\t1\t/absent/g2.fa\t/absent/g2.fa\n"
			+"g9\t9\t1\t2\t/absent/g9.fa\t/absent/g9.fa\n");
		write("manifest.tsv", "assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role\n"
			+manifest(2)+manifest(9)+manifest(7));
		write("counts.tsv", "assembly_id\tinput_id\tcohort\tmember\tnative_events\tnative_loci\tconstruction_loci\ttruncation_only_loci\n"
			+counts(2, "RF00174", 1, 0)+counts(2, "RF01482", 1, 1)+counts(2, "RF01689", 0, 0)
			+counts(9, "RF00174", 1, 1)+counts(9, "RF01482", 0, 0)+counts(9, "RF01689", 1, 1)
			+counts(7, "RF00174", 0, 0)+counts(7, "RF01482", 0, 0)+counts(7, "RF01689", 0, 0));
		final String raw="assembly_id\tinput_id\tcohort\tsource_fasta\tmember\tcontig\tstart\tstop\tstrand\tnative_events\tuntruncated_events\tconstruction_eligible\texact_scored_links\n"
			+locus(2, "RF00174", "+", false)+locus(2, "RF01482", "-", true)
			+locus(9, "RF00174", "+", true)+locus(9, "RF01689", "-", true);
		write("loci.tsv", raw);
		final String sources="source_id\tmember\tassembly_id\tinput_id\tcontig\tstart1\tstop1\tstrand\tlength\ttaxid\tchunk\n"
			+source(2, "RF01482", "-", 1)+source(9, "RF00174", "+", 2)+source(9, "RF01689", "-", 3);
		write("sources.tsv", sources);
		invoke("valid", "loci.tsv", "sources.tsv");
		check(read("valid/panel1.genomes.txt").equals("/absent/g7.fa\n"), "The negative first genome must not be filtered out");
		check(read("valid/panel3.genomes.txt").equals("/absent/g7.fa\n/absent/g2.fa\n/absent/g9.fa\n"), "Both panels must retain the exact saved prefix order");
		final String header="seqid\tfamily\tstart\tstop\tstrand\tinput_ordinal0\n";
		check(read("valid/panel1.raw_native.truth.tsv").equals(header), "A negative panel retains an explicit empty truth view");
		check(read("valid/panel3.raw_native.truth.tsv").equals(header
			+"same\tRF00174\t5\t8\t+\t1\n"+"same\tRF01482\t5\t8\t-\t1\n"
			+"same\tRF00174\t5\t8\t+\t2\n"+"same\tRF01689\t5\t8\t-\t2\n"),
			"Truncated raw loci and repeated coordinates must retain their own panel input ordinal");
		final String summary=read("valid/panels.tsv");
		check(summary.contains("1\tRF00174\t0\t0\t0\t0\tNA_UNION_REFERENCE_ONLY\n")
			&& summary.contains("3\tRF00174\t2\t1\t2\t1\tNA_UNION_REFERENCE_ONLY\n"), "Raw and construction positive genomes have different denominators; scored subtype counts must not be fabricated");
		final String coverage=read("valid/seed_coverage.tsv");
		check(coverage.contains("RF01482_input_2_locus_1\tg2\tg2\tRF01482\tRF01482_input_2_locus_1\tTRAIN\tfalse\n")
			&& coverage.split("\n").length==4 && !coverage.contains("RF00174_input_2"), "Only admitted original IDs enter all-resource seed coverage");
		final String binding="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff\n";
		write("completed.tsv", binding+"0\tg7\t7\t/absent/g7.fa\t/absent/g7.gff\n");
		invoke("suffix", "loci.tsv", "sources.tsv", "chunkstart=1", "chunksize=1", "completedprefix="+path("completed.tsv"));
		check(read("suffix/chunk1.genomes.txt").equals("/absent/g2.fa\n")
			&& read("suffix/chunk2.genomes.txt").equals("/absent/g9.fa\n"), "Suffix chunks retain exactly the uncalled saved order");
		check(read("suffix/chunk1.sources.tsv").equals(binding+"0\tg2\t2\t/absent/g2.fa\t/absent/g2.gff\n")
			&& read("suffix/chunk2.sources.tsv").equals(binding+"0\tg9\t9\t/absent/g9.fa\t/absent/g9.gff\n"),
			"Each chunk restarts caller ordinals while retaining original input IDs and source paths");
		check(read("suffix/panel3.sources.tsv").equals(read("valid/panel3.sources.tsv")),
			"Optional suffix output cannot alter the full ordered source bindings");
		invoke("tail", "loci.tsv", "sources.tsv", "chunkstart=1", "chunksize=3", "completedprefix="+path("completed.tsv"));
		check(read("tail/chunk1.genomes.txt").equals("/absent/g2.fa\n/absent/g9.fa\n")
			&& !new File(path("tail/chunk2.genomes.txt")).exists(), "The last chunk stops at the roster, without padding or duplication");
		write("wrong_prefix.tsv", binding+"0\tg7\t99\t/absent/g7.fa\t/absent/g7.gff\n");
		rejectPrefix("wrong_prefix", "wrong_prefix.tsv", 1, "Completed prefix must match");
		rejectPrefix("short_prefix", "completed.tsv", 2, "Completed prefix ended");
		write("long_prefix.tsv", read("valid/panel3.sources.tsv"));
		rejectPrefix("long_prefix", "long_prefix.tsv", 1, "exactly chunkstart");
		write("duplicate_sources.tsv", sources+source(2, "RF01482", "-", 1));
		reject("duplicate", "loci.tsv", "duplicate_sources.tsv", "Original ID/geometry");
		write("changed_loci.tsv", raw.replace("\t0\tfalse\t0\n", "\t1\ttrue\t0\n"));
		reject("changed", "changed_loci.tsv", "sources.tsv", "Raw truth must retain truncations");
		write("foreign_loci.tsv", raw.replace("g2\t2\tfixture", "g2\t99\tfixture"));
		reject("foreign", "foreign_loci.tsv", "sources.tsv", "Metadata outside");
	}

	private void reject(final String output, final String loci, final String sources, final String reason){
		boolean rejected=false;
		try{invoke(output, loci, sources);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(reason);}
		check(rejected && !new File(path(output+"/complete.tsv")).exists(), "Invalid metadata must fail without completion: "+output);
	}
	private void rejectPrefix(final String output, final String prefix, final int skip, final String reason){
		boolean rejected=false;
		try{invoke(output, "loci.tsv", "sources.tsv", "chunkstart="+skip, "chunksize=1", "completedprefix="+path(prefix));}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(reason);}
		check(rejected && !new File(path(output)).exists(), "A changed completed population must fail before creating outputs: "+output);
	}
	private void invoke(final String output, final String loci, final String sources, final String... options){
		final String[] base={"in="+path("order.tsv"), "manifest="+path("manifest.tsv"), "counts="+path("counts.tsv"),
			"loci="+path(loci), "sourceids="+path(sources), "out="+path(output), "expectedgenomes=3", "expected=1,1,1", "rawexpected=2,1,1", "panels=1,3", "t=1"};
		final String[] args=java.util.Arrays.copyOf(base, base.length+options.length);
		System.arraycopy(options, 0, args, base.length, options.length);
		CobalaminNativePanels.main(args);
	}
	private String manifest(final int id){
		return "g"+id+'\t'+id+"\tfixture\t/absent/g"+id+".fa\t/absent/g"+id+".gff\t.\t.\t/absent/g"+id+".fa\t/absent/g"+id+".gff\tfasta\tncrna\tfixture\n";
	}
	private String counts(final int id, final String member, final int raw, final int eligible){
		return "g"+id+'\t'+id+"\tfixture\t"+member+'\t'+raw+'\t'+raw+'\t'+eligible+'\t'+(raw-eligible)+'\n';
	}
	private String locus(final int id, final String member, final String strand, final boolean eligible){
		return "g"+id+'\t'+id+"\tfixture\t/absent/g"+id+".fa\t"+member+"\tsame\t5\t8\t"+strand+"\t1\t"+(eligible ? 1 : 0)+'\t'+eligible+"\t0\n";
	}
	private String source(final int id, final String member, final String strand, final int ordinal){
		return member+"_input_"+id+"_locus_"+ordinal+'\t'+member+"\tg"+id+'\t'+id+"\tsame\t5\t8\t"+strand+"\t4\t.\t0\n";
	}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean condition, final String reason){if(!condition){throw new AssertionError(reason);}}
	private final String root;
}
