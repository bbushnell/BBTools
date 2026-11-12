package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.GradeNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.WriteNcrnaReference;

/** Ordered aliases, zero-call populations and invalid selection metadata. These
 * tests do not claim registry authentication; the project Bash fixture tests it.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaEvaluationPanelTest{

	NcrnaEvaluationPanelTest(final File root_){root=root_;}
	void run(){
		check(root.mkdir(), "Evaluation metadata fixture root must be new");
		write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/never/a.fa\t10\t100\t1\n2\t/never/b.fa\t20\t200\t1\n");
		write("members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\npanel\t1\t1\t/alias/a.fa\tA\tASM_A\t1\tprok\tNA\tNA\n"
			+"panel\t2\t2\t/never/b.fa\tB\tASM_B\t2\tprok\tNA\tNA\n");
		write("caller.config", "# B first; no query calls are needed to select A\nin=/never/b.fa\nin=/alias/a.fa\ncallcds=t\n");
		write("reference.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\n1\t/never/a.fa\t/never/a.gff\n2\t/never/b.fa\t/never/b.gff\n");
		write("observed.tsv", WriteNcrnaReference.CONTIG_HEADER+"\n1\t%2Fnever%2Fa.fa\tcontig%25a\t100\n2\t/never/b.fa\tcontig_b\t200\n");
		invoke("good");
		check(Arrays.equals(lines("good.reference.tsv"), new String[]{GradeNcrnaReference.MANIFEST_HEADER,
			"2\t/never/b.fa\t/never/b.gff", "1\t/never/a.fa\t/never/a.gff"}),
			"Selection must preserve caller order and both canonical inputs without reading missing payloads");
		check(lines("good.genomes.tsv")[2].equals("2\t1\tASM_A\t1\tNA\tNA\tNA\t/alias/a.fa\tNA\t100\t1"),
			"Caller alias must survive; unknown locus count must not become zero truth");
		check(lines("good.contigs.tsv")[1].equals("1\tASM_A\tcontig%a\t100"), "Saved contig escaping must decode exactly once");
		write("single.config", "in=/alias/a.fa\n");
		invoke("single", "in="+path("single.config"));
		check(lines("single.reference.tsv").length==2, "An explicit one-input selection must retain exactly that input");
		write("duplicate.config", "in=/alias/a.fa\nin=/never/a.fa\n");
		bad("duplicate", "Repeated catalog input", "in="+path("duplicate.config"));
		write("unknown.config", "in=/never/foreign.fa\n");
		bad("unknown", "not a member", "in="+path("unknown.config"));
		write("nested.config", "config=elsewhere\nin=/alias/a.fa\n");
		bad("nested", "flat unpaired", "in="+path("nested.config"));
		write("missing_reference.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\n1\t/never/a.fa\t/never/a.gff\n");
		bad("missing", "complete registered catalog", "reference="+path("missing_reference.tsv"));
		write("raw.tsv", "seqid\tfamily\tstart\tstop\tstrand\nx\tF\t1\t6\t+\n");
		bad("raw", null, "reference="+path("raw.tsv"));
		write("bad_contigs.tsv", WriteNcrnaReference.CONTIG_HEADER+"\n1\t/never/a.fa\tx\t99\n2\t/never/b.fa\ty\t200\n");
		bad("dimensions", "entire input", "observed="+path("bad_contigs.tsv"));
		write("collision.tsv", WriteNcrnaReference.CONTIG_HEADER+"\n1\t/never/a.fa\tx\t100\n2\t/never/b.fa\tx\t200\n");
		bad("collision", "ambiguous", "observed="+path("collision.tsv"));
		bad("source", "Unknown or empty", "sourcepanel=other");
		System.err.println("NCRNA_EVALUATION_PANEL_TEST_PASS");
	}
	private void invoke(final String name, final String... extra){
		final String[] base={"mode=panel", "in="+path("caller.config"), "census="+path("census.tsv"),
			"memberships="+path("members.tsv"), "reference="+path("reference.tsv"), "observed="+path("observed.tsv"),
			"sourcepanel=panel", "out="+path(name), "t=1"};
		final String[] args=Arrays.copyOf(base, base.length+extra.length);
		System.arraycopy(extra, 0, args, base.length, extra.length);
		BuildNcrnaReference.main(args);
	}
	private void bad(final String name, final String fragment, final String... extra){
		try{invoke(name, extra);}
		catch(IllegalArgumentException failure){
			check(fragment==null || failure.getMessage().contains(fragment), "Unexpected failure: "+failure);
			check(!new File(path(name+".reference.tsv")).exists(), "Failed selection cannot publish a completed reference manifest");
			return;
		}
		throw new AssertionError("Expected invalid panel rejection: "+name);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private void write(final String name, final String text){
		try(Output output=new Output(path(name))){output.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private String[] lines(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input input=new Input(path(name))){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				result.add(new String(row, StandardCharsets.UTF_8));
			}
		}
		return result.toArray(new String[result.size()]);
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private final File root;
}
