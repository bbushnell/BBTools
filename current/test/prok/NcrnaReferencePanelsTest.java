package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;

import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Exact panel input counts despite repeated and interleaved source memberships.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferencePanelsTest{
	NcrnaReferencePanelsTest(final File root_){root=root_;}
	void run(){
		check(root.mkdir(), "Panel-count fixture directory must be new");
		write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER
			+"\n10\t/missing/one.fa\t101\t100\t1\n7\t/missing/two.fa\t201\t200\t1\n");
		write("members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
			+member("A", 1, 10, "one")+member("B", 1, 10, "one")
			+member("A", 2, 10, "one")+member("A", 3, 7, "two")+member("B", 2, 10, "one"));
		BuildNcrnaReference.main(new String[]{"mode=panels", "census="+path("census.tsv"),
			"memberships="+path("members.tsv"), "out="+path("counts"), "t=1"});
		try(Input in=new Input(path("counts.panel_population.tsv"))){
			check(text(in).equals("source\tmembership_rows\tdistinct_input_ids\tnucleotide_bases\tstatus"), "Panel count schema changed");
			check(text(in).equals("A\t3\t2\t300\tCATALOG_INPUTS_NOT_INDEPENDENT_ORGANISMS"),
				"Three A memberships cover two exact inputs and 300 bases");
			check(text(in).equals("B\t2\t1\t100\tCATALOG_INPUTS_NOT_INDEPENDENT_ORGANISMS"),
				"Repeated B memberships must not double its input or base denominator");
			check(in.reader.nextLine()==null, "Each declared source is counted once");
		}
		System.err.println("NCRNA_REFERENCE_PANELS_TEST_PASS");
	}
	private static String member(final String source, final int index, final int input, final String name){
		return source+'\t'+index+'\t'+input+"\t/missing/"+name+".fa\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n";
	}
	private static String text(final Input in){
		final byte[] row=in.reader.nextLine();
		check(row!=null, "Missing expected panel count row");
		return new String(row, StandardCharsets.UTF_8);
	}
	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
}
