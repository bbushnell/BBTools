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
		sourceUnion();
		System.err.println("NCRNA_EVALUATION_PANEL_TEST_PASS");
	}
	/** A caller panel may cross original catalog sources; exact input ownership,
	 * aliases and roster metadata remain unambiguous across their explicit union. */
	private void sourceUnion(){
		final String first="left\t1\t1\t/alias/a.fa\tA\tASM_A\t1\tprok\tNA\tNA\n";
		final String second="right\t1\t2\t/never/b.fa\tB\tASM_B\t2\tprok\tNA\tNA\n";
		final String overlap="right\t2\t1\t/second/a.fa\tA\tASM_A\t1\tprok\tNA\tNA\n";
		write("union.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\n"+first+second+overlap);
		final String membership="memberships="+path("union.tsv");
		invoke("union", membership, "sourcepanel=left,right");
		check(Arrays.equals(lines("good.reference.tsv"), lines("union.reference.tsv")),
			"Disjoint and overlapping sources must preserve the ordered original input denominator");
		check(lines("union.selection.tsv")[1].startsWith("2\tleft,right\t"),
			"Selection provenance must name the requested source union, not invent one original source");
		write("union_alias.config", "in=/second/a.fa\nin=/never/b.fa\n");
		invoke("union_alias", membership, "sourcepanel=right,left", "in="+path("union_alias.config"));
		check(lines("union_alias.genomes.tsv")[1].contains("\t/second/a.fa\t"),
			"A second source's literal alias must remain usable without rereading sequence payloads");
		bad("union_missing", "Unknown or empty", membership, "sourcepanel=left,missing");
		bad("union_empty", "nonempty, distinct", membership, "sourcepanel=left,");
		bad("union_duplicate_source", "nonempty, distinct", membership, "sourcepanel=left,left");
		bad("union_restricted", "not a member", membership, "sourcepanel=left");
		bad("union_duplicate_input", "Repeated catalog input", membership, "sourcepanel=left,right",
			"in="+path("duplicate.config"));
		write("union_conflict.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\n"+first+second
			+overlap.replace("\tASM_A\t1\t", "\tASM_A\t9\t"));
		bad("union_conflict", "Conflicting known TaxID", "memberships="+path("union_conflict.tsv"),
			"sourcepanel=left,right");
		metadataResolution(first, second, overlap);
	}

	private void metadataResolution(final String first, final String second, final String overlap){
		final String header=NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\n";
		final String unknown=first.replace("\tASM_A\t1\t", "\tNA\tNA\t");
		write("mixed.tsv", header+unknown+second+overlap);
		invoke("mixed", "memberships="+path("mixed.tsv"), "sourcepanel=left,right");
		check(Arrays.equals(lines("good.genomes.tsv"), lines("mixed.genomes.tsv")),
			"A known assembly and TaxID must resolve NA across every selected-input membership");
		write("mixed_reverse.tsv", header+second+overlap+unknown);
		write("reversed_census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER
			+"\n2\t/never/b.fa\t20\t200\t1\n1\t/never/a.fa\t10\t100\t1\n");
		invoke("mixed_reverse", "memberships="+path("mixed_reverse.tsv"), "sourcepanel=right,left",
			"census="+path("reversed_census.tsv"));
		check(Arrays.equals(lines("mixed.genomes.tsv"), lines("mixed_reverse.genomes.tsv")),
			"Known metadata resolution must be independent of membership, catalog and source-list order");
		invoke("mixed_alias", "memberships="+path("mixed.tsv"), "sourcepanel=left,right",
			"in="+path("union_alias.config"));
		check(Arrays.equals(lines("union_alias.genomes.tsv"), lines("mixed_alias.genomes.tsv")),
			"A selected alias must not choose which source's metadata is considered");
		check(lines("mixed.selection.tsv")[1].endsWith("\tSELECTED_UNION_KNOWN_VALUES_NA_UNKNOWN"),
			"The selected-only known-value rule must be explicit in provenance");
		write("all_unknown.tsv", header+unknown+second+overlap.replace("\tASM_A\t1\t", "\tNA\tNA\t"));
		invoke("all_unknown", "memberships="+path("all_unknown.tsv"), "sourcepanel=left,right");
		check(lines("all_unknown.genomes.tsv")[2].equals("2\t1\tNA\tNA\tNA\tNA\tNA\t/alias/a.fa\tNA\t100\t1"),
			"All unknown claims stay NA; no assembly or TaxID may be inferred");
		write("assembly_conflict.tsv", header+first+second+overlap.replace("ASM_A", "ASM_OTHER"));
		bad("assembly_conflict", "Conflicting known assembly", "memberships="+path("assembly_conflict.tsv"),
			"sourcepanel=left,right");
		write("local_names.tsv", header+first+second+overlap.replace("\tA\tASM_A", "\tsource_local_name\tASM_A"));
		invoke("local_names", "memberships="+path("local_names.tsv"), "sourcepanel=left,right");
		check(Arrays.equals(lines("good.genomes.tsv"), lines("local_names.genomes.tsv")),
			"Source-local names are retained claims, not another canonical-input identity");
		write("alias_owner.tsv", header+first+second.replace("/never/b.fa", "/alias/a.fa")+overlap);
		bad("alias_owner", "Ambiguous FASTA ownership", "memberships="+path("alias_owner.tsv"), "sourcepanel=left,right");
		write("offpanel_census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER
			+"\n1\t/never/a.fa\t10\t100\t1\n2\t/never/b.fa\t20\t200\t1\n3\t/never/c.fa\t30\t300\t1\n");
		write("offpanel_reference.tsv", GradeNcrnaReference.MANIFEST_HEADER
			+"\n1\t/never/a.fa\t/never/a.gff\n2\t/never/b.fa\t/never/b.gff\n3\t/never/c.fa\t/never/c.gff\n");
		final String third="left\t2\t3\t/never/c.fa\tC\tASM_C\t3\tprok\tNA\tNA\n";
		final String thirdUnknown="right\t3\t3\t/alias/c.fa\tC\tNA\tNA\tprok\tNA\tNA\n";
		write("offpanel.tsv", header+first+second+overlap+third+thirdUnknown);
		invoke("offpanel", "memberships="+path("offpanel.tsv"), "sourcepanel=left,right",
			"census="+path("offpanel_census.tsv"), "reference="+path("offpanel_reference.tsv"));
		write("offpanel_conflict.tsv", header+first+second+overlap+third
			+thirdUnknown.replace("\tNA\tNA\tprok", "\tASM_OTHER\t999\tprok"));
		invoke("offpanel_conflict", "memberships="+path("offpanel_conflict.tsv"), "sourcepanel=left,right",
			"census="+path("offpanel_census.tsv"), "reference="+path("offpanel_reference.tsv"));
		check(Arrays.equals(lines("good.genomes.tsv"), lines("offpanel.genomes.tsv"))
			&& Arrays.equals(lines("good.genomes.tsv"), lines("offpanel_conflict.genomes.tsv")),
			"Off-panel NA/known and incompatible known metadata must not add, reorder or choose selected roster values");
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
