package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BindNcrnaCalls;
import prok.BuildNcrnaReference;
import prok.GradeNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.WriteNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Literal saved-call ownership, conservation and budget cases; no genome payloads.
 * @author Brian Bushnell, Raiden
 */
final class BindNcrnaCallsTest{

	BindNcrnaCallsTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Call-binding fixture root must be new");
		prepare();
		invoke("result", "maxlines=5", "maxfiles=10");
		check(Arrays.equals(lines("result/input_1.gff.gz"), new String[]{"##gff-version 3", a.trim(), a.trim(), unknown.trim()}),
			"Interleaving must preserve per-input bytes/order, duplicate calls and unknown models");
		check(Arrays.equals(lines("result/input_2.gff.gz"), new String[]{"##gff-version 3", b.trim(), excluded.trim()}),
			"Excluded-family RNA must remain in the projected population");
		check(Arrays.equals(lines("result/input_3.gff.gz"), new String[]{"##gff-version 3"}), "Input with only tRNA must have an explicit empty RNA query");
		check(Arrays.equals(lines("result/line_map.tsv.gz"), new String[]{"source_physical_line\tinput_id\tquery_physical_line",
			"4\t1\t2", "6\t1\t3", "9\t1\t4", "5\t2\t2", "10\t2\t3"}), "Physical-line witnesses must include header/blank offsets and preserve interleaving");
		final String[] counts=lines("result/original_counts.tsv");
		final String qualification="\tPROJECTED_ZERO_OTHER_TYPES_IS_NOT_SOURCE_ABSENCE";
		check(counts.length==4 && counts[1].equals("1\t4\t1\t0\t0\t3\t3"+qualification)
			&& counts[2].equals("2\t3\t0\t0\t1\t2\t2"+qualification) && counts[3].equals("3\t1\t0\t1\t0\t0\t0"+qualification),
			"All eight original features must be counted, including three omitted non-RNA features");
		final String[] summary=lines("result/summary.tsv");
		check(summary.length==2 && summary[1].startsWith("fixture\t3\t10\t1\t1\t1\t5\t1\t1\t1\t"), "Original line/type inventory and one empty input must be explicit");
		check(Arrays.equals(lines("result/calls.tsv"), new String[]{GradeNcrnaReference.MANIFEST_HEADER,
			"1\t/canonical/a%name.fa\t"+path("result/input_1.gff.gz"),
			"2\t/canonical/b.fa\t"+path("result/input_2.gff.gz"),
			"3\t/canonical/c.fa\t"+path("result/input_3.gff.gz")}), "Query manifest must bind canonical paths including the empty input");
		check(Arrays.equals(lines("result/contigs.tsv.gz"), lines("observed.tsv")), "Escaped literal percent names must round-trip into selected observations");
		check(lines("result/bindings.tsv")[1].contains("\t/original/a.fa\t/canonical/a%25name.fa\t2\t"), "Original alias and canonical FASTA must remain separately visible");
		check(lines("result/sources.tsv")[1].endsWith("\t00000000000000000001\tfixture\tDECLARED_METADATA_IDENTITY_NOT_SEQUENCE_VERIFICATION"),
			"Supplied source fingerprint is a declaration, not an invented independently verified digest");
		check(new File(path("result")).list().length==10, "Three payloads and seven metadata files exhaust the declared file budget");
		bad("row_budget", "line budget", "maxlines=4");
		bad("file_budget", "file budget", "maxfiles=9");
		bad("buffer_budget", "capacity budget", "maxbuffer=255");
		bad("output_budget", "output-byte budget", "maxoutput=1");
		badInput("foreign", original.replace("b\tBBTools\tRNA", "foreign\tBBTools\tRNA"), "unknown ownership");
		badInput("source", original.replace("a%1\tBBTools\tCDS", "a%1\tother\tCDS"), "nine-field CallGenes");
		badInput("type", original.replace("\tCDS\t", "\tgene\t"), "Unknown feature type");
		badInput("bounds", original.replace("\t1\t6\t", "\t1\t101\t"), "invalid contig bounds");
		final String duplicateContigs=write("duplicate_contigs.tsv", panelContigText.replace("B\tASM_B\tb\t", "B\tASM_B\ta%1\t"));
		bad("duplicate_contig", "repeated literal contig", "panelcontigs="+duplicateContigs);
		final String wrongObserved=write("wrong_observed.tsv", observedText.replace("\ta%251\t100", "\ta%251\t99"));
		bad("observed_mismatch", "differs from its original panel", "observed="+wrongObserved);
		final String wrongOrder=write("wrong_order.config", "in=/canonical/b.fa\nin=/original/a.fa\nin=/canonical/c.fa\n");
		bad("config_order", "exact original panel FASTAs", "callerconfig="+wrongOrder);
		final String missingEmpty=write("missing_empty.config", "in=/original/a.fa\nin=/canonical/b.fa\n");
		bad("config_incomplete", "complete panel", "callerconfig="+missingEmpty);
		final StringBuilder many=new StringBuilder("##gff-version 3\n");
		for(int i=0; i<33; i++){many.append(a);}
		final String growth=write("growth.gff", many.toString());
		invoke("growth", "in="+growth, "maxlines=33");
		check(lines("growth/input_1.gff.gz").length==34 && lines("growth/line_map.tsv.gz").length==34,
			"Byte and primitive-line arrays must grow without losing duplicate rows or witnesses");
		System.err.println("BIND_NCRNA_CALLS_TEST_PASS");
	}

	private void prepare(){
		census=write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/canonical/a%name.fa\t123\t100\t1\n2\t/canonical/b.fa\t123\t100\t1\n3\t/canonical/c.fa\t123\t100\t1\n");
		memberships=write("memberships.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
			+"panel\t1\t1\t/original/a.fa\tA\tASM_A\t1\tprok\tNA\tNA\n"
			+"panel\t2\t2\t/canonical/b.fa\tB\tASM_B\t2\tprok\tNA\tNA\n"
			+"panel\t3\t3\t/canonical/c.fa\tC\tASM_C\t3\tprok\tNA\tNA\n");
		genomes=write("genomes.tsv", "index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n"
			+"1\tA\tASM_A\t1\tbacteria\tP\tNA\t/original/a.fa\t0\t100\t1\n"
			+"2\tB\tASM_B\t2\tbacteria\tP\tNA\t/canonical/b.fa\t0\t100\t1\n"
			+"3\tC\tASM_C\t3\tbacteria\tP\tNA\t/canonical/c.fa\t0\t100\t1\n");
		panelContigText=BindNcrnaCalls.PANEL_CONTIGS+"\nA\tASM_A\ta%1\t100\nB\tASM_B\tb\t100\nC\tASM_C\tc\t100\n";
		panelContigs=write("panel_contigs.tsv", panelContigText);
		observedText=WriteNcrnaReference.CONTIG_HEADER+'\n'+fields("1", "/canonical/a%name.fa", "a%1", "100")
			+fields("2", "/canonical/b.fa", "b", "100")+fields("3", "/canonical/c.fa", "c", "100");
		observed=write("observed.tsv", observedText);
		callerConfig=write("caller.config", "# Original ordered input list\nin=/original/a.fa\nin=/canonical/b.fa\nin=/canonical/c.fa\ncallcds=t\n");
		a=call("a%1", "RNA", "model:F_1,family:F");
		b=call("b", "RNA", "model:G_1,family:G");
		unknown=call("a%1", "RNA", "model:unknown_prefix");
		excluded=call("b", "RNA", "model:X_1,family:X");
		original="##gff-version 3\n\n"+call("a%1", "CDS", ".")+a+b+a+call("c", "tRNA", ".")+call("b", "rRNA", ".")+unknown+excluded;
		input=write("original.gff", original);
	}
	private void badInput(final String name, final String text, final String expected){bad(name, expected, "in="+write(name+".gff", text));}
	private void bad(final String name, final String expected, final String... args){
		try{invoke(name, args);}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(expected), "Unexpected binding failure: "+failure);
			check(!new File(path(name+"/calls.tsv")).exists(), "Failed binding must not publish the final query manifest");
			return;
		}
		throw new AssertionError("Expected binding failure: "+expected);
	}
	private void invoke(final String name, final String... extra){
		final String[] base={"in="+input, "genomes="+genomes, "panelcontigs="+panelContigs, "observed="+observed,
			"census="+census, "memberships="+memberships, "callerconfig="+callerConfig, "baseline=fixture", "inputsha80=00000000000000000001",
			"out="+path(name), "t=1", "bgzip=f", "zipthreads=1"};
		final String[] args=Arrays.copyOf(base, base.length+extra.length);
		System.arraycopy(extra, 0, args, base.length, extra.length);
		BindNcrnaCalls.main(args);
	}
	private static String call(final String contig, final String type, final String attributes){return contig+"\tBBTools\t"+type+"\t1\t6\t.\t+\t.\t"+attributes+'\n';}
	private static String fields(final String... fields){
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<fields.length; i++){
			if(i>0){b.tab();}
			BuildNcrnaReference.field(b, fields[i]);
		}
		b.nl();
		return new String(b.array, 0, b.length(), StandardCharsets.UTF_8);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private String write(final String name, final String text){
		final String path=path(name);
		try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return path;
	}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private String input, genomes, panelContigs, observed, census, memberships, callerConfig;
	private String a, b, unknown, excluded, original, panelContigText, observedText;
}
