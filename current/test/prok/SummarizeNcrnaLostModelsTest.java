package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaLostModels;
import shared.Shared;

/** Literal two-panel oracle: deduplicated removed winners, retained winners,
 * a gained call and an overlap-present loss that must never request restoration.
 * @author Raiden
 */
public final class SummarizeNcrnaLostModelsTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SummarizeNcrnaLostModelsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown loss-tally fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			new SummarizeNcrnaLostModelsTest(parser.out1).run();
			pp.outstream.println("NCRNA_LOST_MODELS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private SummarizeNcrnaLostModelsTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		final String absent="LOST_GRADED_REFERENCE\tNO_OVERLAPPING_CLAN_CALL\t";
		final String common="RF00174\tmain%25key n=42\tfalse\t.\t.\n";
		final String first=HEADER
			+"a\tr\tcontig\t1\t99\t+\t"+absent+common
			+"b\tr\tcontig\t1\t99\t+\t"+absent+common
			+"c\tr\tcontig\t1\t99\t+\t"+absent+"RF01689\tkept n=901\ttrue\t.\t.\n"
			+"d\tr\tcontig\t1\t99\t+\tLOST_GRADED_REFERENCE\tOVERLAPPING_CLAN_CALL_PRESENT\tRF01482\tignored\tfalse\tRF00174\tother\n"
			+"e\tr\tcontig\t1\t99\t+\tGAINED_GRADED_REFERENCE\tNO_OVERLAPPING_CLAN_CALL\t.\t.\tNA\tRF00174\tnew\n";
		final String second=HEADER
			+"x\tr\tcontig\t1\t99\t+\t"+absent+common
			+"y\tr\tcontig\t1\t99\t+\t"+absent+"RF01482\tminor n=8\tfalse\t.\t.\n";
		write("one.tsv", first);
		write("two.tsv", second);
		invoke("one.tsv,two.tsv", "one,two", "valid");
		final String summary=read("valid.summary.tsv");
		check(summary.contains("one\tALL\t4\t3\t2\t1\t1\t1\t"), "First panel has four graded losses, three absences and only one distinct removed winner");
		check(summary.contains("two\tALL\t2\t2\t2\t0\t2\t0\t"), "Second panel must retain its own locus and model counts");
		check(summary.contains("UNION\tALL\t6\t5\t4\t1\t2\t1\t"), "Union model count deduplicates across panels while evidence rows remain additive");
		check(summary.contains("UNION\tRF00174\t3\t3\t3\t0\t1\t0\t"), "Three main-member losses across panels identify one unique removed model");
		check(read("valid.removed.ids.txt").equals("main%key\nminor\n"), "Only distinct removed native IDs belong in the name list; metadata and GFF percent escaping are removed");
		check(read("valid.RF00174.ids.txt").equals("main%key\n") && read("valid.RF01482.ids.txt").equals("minor\n")
			&& read("valid.RF01689.ids.txt").isEmpty(), "Member-specific lists must partition the unique removed winners, including empty members");
		check(read("valid.already_present.loci.tsv").contains("one\tc\tr\t") && !read("valid.removed.ids.txt").contains("kept"), "Retained-winner absences must remain individually visible and never be added twice");
		write("missing.tsv", first.replace("\ttrue\t", "\tNA\t"));
		reject("missing.tsv", "missing", "Require measured");
		reject("one.tsv,one.tsv", "duplicate", "Duplicate panel");
		write("conflict.tsv", second.replace("main%25key n=42\tfalse", "main%25key n=42\ttrue"));
		reject("one.tsv,conflict.tsv", "conflict", "conflicting");
	}

	private void reject(final String input, final String output, final String reason){
		boolean rejected=false;
		try{invoke(input, input.indexOf(',')>=0 ? "one,two" : "one", output);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(reason);}
		check(rejected && !new File(path(output+".removed.ids.txt")).exists(), "Invalid evidence must fail before an add-back list exists: "+reason);
	}
	private void invoke(final String input, final String panels, final String output){
		final String[] names=input.split(",");
		for(int i=0; i<names.length; i++){names[i]=path(names[i]);}
		SummarizeNcrnaLostModels.main(new String[]{"in="+String.join(",", names), "panels="+panels, "out="+path(output)});
	}
	private void write(final String name, final String value){try(Output output=new Output(path(name))){output.writer.print(value);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final String root;
	private static final String HEADER="input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_member\tbaseline_model\tbaseline_model_in_subset\tcandidate_member\tcandidate_model\n";
}
