package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BindNcrnaCalls;
import prok.BuildNcrnaReference;
import prok.NcrnaCalibrationPrep;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PartitionNcrnaBenchmarkInputs;
import prok.WriteNcrnaReference;
import shared.Shared;
import structures.ByteBuilder;

/** Independent finite partitions and actual saved-call binder compatibility.
 * No sequence payload exists at the fixture's declared FASTA paths.
 * @author Brian Bushnell, Raiden
 */
public final class PartitionNcrnaBenchmarkInputsTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, PartitionNcrnaBenchmarkInputsTest.class, false);
		final Parser parser=new Parser();
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		check(parser.out1!=null, "Supply out=fresh_fixture_directory");
		final File root=new File(parser.out1).getAbsoluteFile();
		check(!root.exists(), "Fixture output must be fresh");
		Files.createDirectory(root.toPath());
		try{test(new Fixture(root));}
		finally{Shared.closeStream(pp.outstream);}
		System.err.println("NCRNA_BENCHMARK_PARTITION_TEST_PASS");
	}

	private static void test(final Fixture f) throws Exception{
		f.prepare();
		f.run("partition", f.observed, f.selection);
		final String prefix=f.path("partition")+File.separator;
		final String[] assignment={"input_id\tfasta\tbatch\tbatch_input_index",
			"1\t"+fasta(1)+"\t0\t1", "2\t"+fasta(2)+"\t0\t2",
			"3\t"+fasta(3)+"\t1\t1", "4\t"+fasta(4)+"\t2\t1"};
		check(Arrays.equals(lines(prefix+"assignments.tsv"), assignment),
			"The ordered path graph {a},{c},{a,b},{b,c} must yield first-fit 0,0,1,2 without renaming inputs");
		check(lines(prefix+"summary.tsv")[1].equals(
			"4\t600\t6\t3\t3\t4\t6\t2\t3\tSaved metadata; greedy batch count is not a proven minimum"),
			"The literal six-record population needs achieved3/lowerbound2, not a false minimum claim");
		check(Arrays.equals(lines(prefix+"batch_0.inputs.config"), new String[]{"in="+fasta(1), "in="+fasta(2)}),
			"Batch inputs must retain original census order and exact UTF-8/percent/equals paths");
		check(Arrays.equals(lines(prefix+"batches.tsv"), new String[]{"batch\tinputs\tnucleotide_bases\tfasta_records\tbinder_maxfiles",
			"0\t2\t200\t2\t9", "1\t1\t200\t2\t8", "2\t1\t200\t2\t8"}),
			"Each batch must retain its literal input/base/record dimensions and binder budget for grading-cost planning");
		check(Arrays.equals(lines(prefix+"batch_0.contigs.tsv"), new String[]{BindNcrnaCalls.PANEL_CONTIGS,
			"1\tNA\t"+NAME+"\t100", "2\tNA\tc\t100"}),
			"Escaped names must become literal binder names, with canonical IDs and consistent NA assembly");
		check(lines(prefix+"batch_0.genomes.tsv")[0].equals(NcrnaCalibrationPrep.ROSTER_HEADER)
			&& lines(prefix+"batch_0.binder.config")[0].equals("maxfiles=9"),
			"Generated rosters must use the existing schema and include all seven binder metadata files");
		bind(f, 0, true);
		bind(f, 1, true);
		bind(f, 2, false);
		check(Arrays.equals(lines(f.path("bound_0/input_1.gff.gz")), new String[]{GFF_HEADER, CALL}),
			"A literal percent/Unicode call must remain attributed to input1");
		check(Arrays.equals(lines(f.path("bound_0/input_2.gff.gz")), new String[]{GFF_HEADER}),
			"An input with no calls must retain an explicit empty query in a nonempty batch");
		check(Arrays.equals(lines(f.path("bound_1/input_1.gff.gz")), new String[]{GFF_HEADER, CALL})
			&& lines(f.path("bound_1/calls.tsv"))[1].startsWith("3\t"+fasta(3)+"\t"),
			"The same name in another batch must bind independently to input3");
		check(lines(f.path("bound_2/calls.tsv"))[1].startsWith("4\t"+fasta(4)+"\t")
			&& Arrays.equals(lines(f.path("bound_2/input_1.gff.gz")), new String[]{GFF_HEADER}),
			"A wholly empty-call batch must still bind every selected input");

		final String reordered=f.write("reordered.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'
			+f.blocks[3]+f.blocks[2]+f.blocks[1]+f.blocks[0]);
		f.run("reordered", reordered, f.selection);
		for(final String file:new File(prefix).list()){
			check(Arrays.equals(Files.readAllBytes(new File(prefix+file).toPath()),
				Files.readAllBytes(new File(f.path("reordered/"+file)).toPath())),
				"Saved input-block order must not change census-ordered partitions or output bytes");
		}
		final String interleaved=f.write("interleaved.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'
			+fieldRow("3", fasta(3), NAME, "100")+f.blocks[0]
			+fieldRow("4", fasta(4), "b", "100")+f.blocks[1]
			+fieldRow("3", fasta(3), "b", "100")+fieldRow("4", fasta(4), "c", "100"));
		f.run("interleaved", interleaved, f.selection);
		for(final String file:new File(prefix).list()){
			check(Arrays.equals(Files.readAllBytes(new File(prefix+file).toPath()),
				Files.readAllBytes(new File(f.path("interleaved/"+file)).toPath())),
				"Interleaving worker-buffered inputs must preserve the same census-ordered partition and binding bytes");
		}
		final String subset=f.write("subset.tsv", SELECTION_HEADER+'\n'+f.selected(1, 1)+f.selected(3, 2));
		f.run("subset", f.observed, subset);
		check(lines(f.path("subset/summary.tsv"))[1].startsWith("2\t300\t3\t2\t1\t2\t2\t2\t2\t"),
			"Subset conflicts and totals must exclude unselected inputs while preserving full observation validation");

		final String duplicate=f.write("duplicate.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'
			+f.blocks[0]+f.blocks[0]+f.blocks[1]+f.blocks[2]+f.blocks[3]);
		f.reject("duplicate", duplicate, f.selection, "repeats within input");
		final String interleavedDuplicate=f.write("interleaved_duplicate.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'
			+fieldRow("3", fasta(3), NAME, "100")+f.blocks[0]
			+fieldRow("3", fasta(3), NAME, "100")+f.blocks[1]+f.blocks[3]);
		f.reject("interleaved_duplicate", interleavedDuplicate, f.selection, "repeats within input");
		final String badDimension=f.write("bad_dimensions.tsv", f.observedText.replace("\t100\n", "\t99\n"));
		f.reject("dimensions", badDimension, f.selection, "dimensions disagree");
		final String missing=f.write("missing.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'+f.blocks[0]+f.blocks[1]+f.blocks[2]);
		f.reject("missing", missing, f.selection, "missing catalog input");
		final String delimiter=f.write("delimiter.tsv", f.observedText.replace("α%251", "x%09y"));
		f.reject("delimiter", delimiter, f.selection, "Decoded contig names");
		final String malformed=f.write("malformed.tsv", f.observedText.replace("α%251", "x%Q0"));
		f.reject("malformed", malformed, f.selection, "Invalid GFF percent escape");
		final String wrongPath=f.write("wrong_path.tsv", f.observedText.replace("/absent/", "/other/"));
		f.reject("path", wrongPath, f.selection, "exact catalog binding");
		f.reject("partition", f.observed, f.selection, "fresh directory");
		clique(f);
	}

	/** The existing binder, rather than another partition implementation, tests interoperability. */
	private static void bind(final Fixture f, final int batch, final boolean hasCall){
		final String gff=f.write("batch_"+batch+".gff", GFF_HEADER+'\n'+(hasCall ? CALL+'\n' : ""));
		final String prefix=f.path("partition/batch_"+batch);
		BindNcrnaCalls.main(new String[]{"in="+gff, "genomes="+prefix+".genomes.tsv", "panelcontigs="+prefix+".contigs.tsv",
			"observed="+f.observed, "census="+f.census, "memberships="+f.memberships, "callerconfig="+prefix+".inputs.config",
			"baseline=fixture", "inputsha80="+DECLARED_FIXTURE_IDENTITY, "config="+prefix+".binder.config", "out="+f.path("bound_"+batch)});
	}

	private static void clique(final Fixture f){
		final String census=f.write("clique_census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'
			+"1\t"+fasta(1)+"\t1\t100\t1\n2\t"+fasta(2)+"\t1\t100\t1\n3\t"+fasta(3)+"\t1\t100\t1\n");
		final String memberships=f.write("clique_memberships.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+f.members(1)+f.members(2)+f.members(3));
		final String selection=f.write("clique_selection.tsv", SELECTION_HEADER+'\n'+f.selected(1, 1)+f.selected(2, 1)+f.selected(3, 1));
		final String contigs=f.write("clique_contigs.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'
			+fieldRow("1", fasta(1), "shared", "100")+fieldRow("2", fasta(2), "shared", "100")+fieldRow("3", fasta(3), "shared", "100"));
		PartitionNcrnaBenchmarkInputs.main(new String[]{"in="+selection, "census="+census, "memberships="+memberships,
			"contigs="+contigs, "out="+f.path("clique")});
		check(lines(f.path("clique/summary.tsv"))[1].startsWith("3\t300\t3\t1\t1\t3\t3\t3\t3\t"),
			"One shared name in three complete inputs requires exactly three distinct batches");
	}

	private static final class Fixture{
		Fixture(final File root_){root=root_;}
		String path(final String name){return new File(root, name).toString();}
		String write(final String name, final String text){
			final String path=path(name);
			try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
			return path;
		}
		void prepare(){
			final ByteBuilder genomes=new ByteBuilder().append(NcrnaEvidenceCatalog.CENSUS_HEADER).nl();
			final ByteBuilder members=new ByteBuilder().append(NcrnaEvidenceCatalog.MEMBERSHIP_HEADER).nl();
			final ByteBuilder chosen=new ByteBuilder().append(SELECTION_HEADER).nl();
			for(int id=1; id<=4; id++){
				final int records=id<=2 ? 1 : 2;
				genomes.append(id).tab().append(fasta(id).getBytes(StandardCharsets.UTF_8)).append("\t1\t").append(100*records).tab().append(records).nl();
				members.append(members(id).getBytes(StandardCharsets.UTF_8));
				chosen.append(selected(id, records).getBytes(StandardCharsets.UTF_8));
			}
			census=write("census.tsv", utf8(genomes));
			memberships=write("memberships.tsv", utf8(members));
			selection=write("selection.tsv", utf8(chosen));
			blocks=new String[]{fieldRow("1", fasta(1), NAME, "100"), fieldRow("2", fasta(2), "c", "100"),
				fieldRow("3", fasta(3), NAME, "100").replace("α%251", "%CE%B1%251")+fieldRow("3", fasta(3), "b", "100"),
				fieldRow("4", fasta(4), "b", "100")+fieldRow("4", fasta(4), "c", "100")};
			observedText=WriteNcrnaReference.CONTIG_HEADER+'\n'+String.join("", blocks);
			observed=write("observed.tsv", observedText);
		}
		String members(final int id){return "source\t"+id+'\t'+id+'\t'+fasta(id)+"\tNA\tNA\tNA\tNA\tNA\tNA\n";}
		String selected(final int id, final int records){return id+"\t"+fasta(id)+"\t"+(100*records)+'\t'+records+"\tUNKNOWN\tUNKNOWN\t0\t0\n";}
		void run(final String output, final String contigs, final String selected){
			PartitionNcrnaBenchmarkInputs.main(new String[]{"in="+selected, "census="+census,
				"memberships="+memberships, "contigs="+contigs, "out="+path(output)});
		}
		void reject(final String output, final String contigs, final String selected, final String diagnostic){
			try{run(output, contigs, selected);}
			catch(IllegalArgumentException failure){
				check(failure.getMessage().contains(diagnostic), "Wrong rejection: "+failure.getMessage());
				return;
			}
			throw new AssertionError("Invalid partition inputs accepted: "+diagnostic);
		}
		final File root;
		String census, memberships, selection, observed, observedText;
		String[] blocks;
	}

	private static String fieldRow(final String... fields){
		final ByteBuilder out=new ByteBuilder();
		for(int i=0; i<fields.length; i++){
			if(i>0){out.tab();}
			BuildNcrnaReference.field(out, fields[i]);
		}
		return utf8(out.nl());
	}
	/** ByteBuilder's String append/toString paths are ASCII; fixtures must preserve literal UTF-8. */
	private static String utf8(final ByteBuilder bytes){return new String(bytes.array, 0, bytes.length(), StandardCharsets.UTF_8);}
	private static String[] lines(final String path){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private static String fasta(final int id){return "/absent/α="+id+"%seq.fa";}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private static final String NAME="α%1", GFF_HEADER="##gff-version 3";
	private static final String CALL=NAME+"\tBBTools\tRNA\t1\t10\t.\t+\t.\tmodel=fixture";
	/** The binder records supplied identities as declarations; this fixture authenticates no file. */
	private static final String DECLARED_FIXTURE_IDENTITY="00000000000000000001";
	private static final String SELECTION_HEADER="input_id\tfasta\tnucleotide_bases\tfasta_records\tdomain\tphylum\tbase_log2\trecord_log2";
}
