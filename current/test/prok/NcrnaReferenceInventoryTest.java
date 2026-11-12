package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.GradeNcrnaReference;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Hand-specified inventory outcomes and parity with the existing subtype grader.
 * FASTA paths intentionally do not exist: a metadata-only inventory must not open
 * them. Three original inputs include a masked-only and an empty reference.
 * @author Ganyu
 */
public final class NcrnaReferenceInventoryTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaReferenceInventoryTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){
					throw new IllegalArgumentException("Unknown parameter: "+arg);
				}
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append,
				"Supply out=fresh_fixture_directory under an existing parent; no input or overwrite");
			new NcrnaReferenceInventoryTest(new File(parser.out1).getAbsoluteFile()).run();
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}
	private NcrnaReferenceInventoryTest(final File root_){
		root=root_;
		check(root.mkdir(), "Fixture directory must be fresh");
	}

	private void run(){
		write("members.tsv", "family\tmember\treference_model\ncobalamin\tRF00174\tRF00174\ncobalamin\tRF01482\tRF01482\ncobalamin\tRF01689\tRF01689\n6S\tRF00013\tRF00013\n6S\tRF01685\tRF01685\ntmRNA\tRF00023\tRF00023\ntmRNA\tRF01851\tRF01851\nglycine\tRF00504\tRF00504\n");
		final ByteBuilder census=new ByteBuilder("input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records\n");
		final ByteBuilder manifest=new ByteBuilder("input_id\tfasta\tgff\n");
		for(int i=1; i<=3; i++){
			census.append(i).tab().append(path("unopened"+i+".fa")).append("\t100\t10000\t1\n");
			manifest.append(i).tab().append(path("unopened"+i+".fa")).tab().append(path("ref"+i+".gff")).nl();
		}
		write("census.tsv", census.toString());
		write("reference.tsv", manifest.toString());
		final ByteBuilder memberships=new ByteBuilder("source\tmember_index\tinput_id\tsupplied_fasta\tsource_genome\tsource_assembly\ttaxid\tdeclared_compartment\treported_gff\treported_sha80\n");
		membership(memberships, "challenge", 1, 1);
		membership(memberships, "challenge", 2, 1); // Same input in a source twice is one positive input.
		membership(memberships, "challenge", 3, 2);
		membership(memberships, "control", 1, 1);
		membership(memberships, "control", 2, 3);
		write("memberships.tsv", memberships.toString());

		final ArrayList<Row> first=new ArrayList<Row>(), second=new ArrayList<Row>();
		final Row multi=parent(first, "1", "cobalamin", "RF00174", 100, false, '+');
		hit(multi, "RF00174", "CONFIDENT", GROUP, -1);
		hit(multi, "RF01482", "CONFIDENT", GROUP, -1);
		final Row marginal=parent(first, "1", "cobalamin", "RF00174", 300, false, '+');
		hit(marginal, "RF00174", "CONFIDENT", GROUP, -1);
		hit(marginal, "RF01689", "MARGINAL", GROUP, -1);
		final Row prior=parent(first, "1", "cobalamin", "RF01689", 500, false, '-');
		prior(prior);
		hit(prior, "RF00174", "COVERED_BY_PROTECTED_PRIOR", PRIOR, 9);
		final Row six=parent(first, "1", "6S", "RF00013", 700, false, '+');
		hit(six, "RF00013", "CONFIDENT", GROUP, -1);
		hit(six, "RF01685", "MARGINAL", GROUP, -1);
		final Row tm=parent(first, "1", "tmRNA", "RF00023", 900, false, '+');
		hit(tm, "RF00023", "CONFIDENT", GROUP, -1);
		hit(tm, "RF01851", "MARGINAL", BRIDGE, -1);
		final Row left=parent(first, "1", "glycine", "RF00504", 1100, false, '-');
		final Row right=parent(first, "1", "glycine", "RF00504", 1200, false, '-');
		hit(left, "RF00504", "CONFIDENT", GROUP, -1);
		hit(right, "RF00504", "CONFIDENT", GROUP, -1);
		parent(first, "1", "glycine_GGAnGA", "RF03166", 1400, true, '-');
		final Row masked=parent(second, "2", "cobalamin", "RF00174", 100, true, '+');
		hit(masked, "RF00174", "MARGINAL", GROUP, -1);
		final Row tmMask=parent(second, "2", "tmRNA", "RF01851", 300, true, '+');
		hit(tmMask, "RF01851", "MARGINAL", GROUP, -1);
		codec.write(path("ref1.gff"), "1", first);
		codec.write(path("ref2.gff"), "2", second);
		codec.write(path("ref3.gff"), "3", new ArrayList<Row>());
		write("events.tsv", events.toString());
		write("witnesses.tsv", witnesses.toString());
		inventory("result", "reference.tsv", "members.tsv", "witnesses.tsv");

		check(lines("result.inputs.tsv").size()==13, "Three original inputs times four unions, including zeros");
		check(lines("result.members.tsv").size()==25, "Three original inputs times eight nonadditive members");
		final ArrayList<String> loci=lines("result.loci.tsv.gz");
		check(loci.size()==10, "Nine distinct selected parents; separate related-class row excluded");
		check(find(loci, "1\t"+left.id+"\t")!=null && find(loci, "1\t"+right.id+"\t")!=null,
			"Adjacent minus-strand aptamers retain separate parent IDs and intervals");
		check(find(loci, "1\t"+tm.id+"\t").contains("\tfalse\tfalse\t1\t"), "Bridge evidence is retained without becoming member support");
		check(row("result.summary.tsv", "UNION\tcobalamin\t.\t").startsWith("UNION\tcobalamin\t.\t3\t1\t1\t1\t3\t1\t1\t1\t"),
			"Three scored parents in one positive genome; masked-only and empty genomes remain in denominator");
		check(row("result.strata.tsv", "challenge\tUNION\tcobalamin\t.\t").startsWith("challenge\tUNION\tcobalamin\t.\t2\t1\t1\t0\t"),
			"Repeated source memberships do not inflate the two distinct source inputs");
		check(row("result.members.tsv", "1\t6S\tRF01685\t").startsWith("1\t6S\tRF01685\tRF01685\t0\t1\tfalse\tNA"),
			"Marginal-only member evidence is not a positive member genome");
		check(row("result.members.tsv", "1\ttmRNA\tRF01851\t").startsWith("1\ttmRNA\tRF01851\tRF01851\t0\t0\tfalse\tNA"),
			"Bridge witness does not promote a second tmRNA member");
		check(row("result.summary.tsv", "EXCLUSIVE\tcobalamin\tRF01482\t").startsWith("EXCLUSIVE\tcobalamin\tRF01482\t3\t0\t2\t1\t0\t2\t"),
			"Multi-support masks exclusive accuracy but not recovery; original masks remain in all views");
		gradeParity();
		final ArrayList<String> mapping=lines("members.tsv");
		final ByteBuilder reordered=new ByteBuilder(mapping.get(0)).nl();
		for(int i=4; i<mapping.size(); i++){reordered.append(mapping.get(i)).nl();}
		for(int i=1; i<4; i++){reordered.append(mapping.get(i)).nl();}
		write("reordered_members.tsv", reordered.toString());
		inventory("reordered", "reference.tsv", "reordered_members.tsv", "witnesses.tsv");
		for(final String line:lines("result.summary.tsv")){
			if(line.startsWith("view\t")){continue;}
			final String[] fields=line.split("\t");
			check(line.equals(row("reordered.summary.tsv", fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t")),
				"Member bit positions and family order must not change coverage or cobalamin reference views");
		}

		write("duplicate.tsv", manifest.toString()+manifest.toString().split("\n")[1]+"\n");
		reject("duplicate", "duplicate.tsv", "members.tsv", "witnesses.tsv", "Duplicate input_id");
		write("incomplete_members.tsv", String.join("\n", lines("members.tsv")).replace("tmRNA\tRF01851\tRF01851\n", "")+"\n");
		reject("incomplete", "reference.tsv", "incomplete_members.tsv", "witnesses.tsv", "unmapped");
		write("substituted.tsv", witnesses.toString());
		reject("substituted", "reference.tsv", "members.tsv", "substituted.tsv", "evidence path differs");
		check(!new File(path("unopened1.fa")).exists(), "Test must not fabricate the missing FASTA inputs");
		System.err.println("NCRNA_REFERENCE_INVENTORY_TEST_PASS");
	}

	private void gradeParity(){
		write("grade_map.tsv", "family\tmodel_prefix\tdisplay_label\tinclude\ncobalamin\tcbA\tA\ttrue\ncobalamin\tcbB\tB\ttrue\ncobalamin\tcbC\tC\ttrue\n");
		write("grade_subtypes.tsv", "subtype\tmodel_prefix\treference_model\nRF00174\tcbA\tRF00174\nRF01482\tcbB\tRF01482\nRF01689\tcbC\tRF01689\n");
		final ByteBuilder queries=new ByteBuilder("input_id\tfasta\tgff\n");
		for(int i=1; i<=3; i++){
			write("query"+i+".gff", "##gff-version 3\n");
			queries.append(i).tab().append(path("unopened"+i+".fa")).tab().append(path("query"+i+".gff")).nl();
		}
		write("queries.tsv", queries.toString());
		GradeNcrnaReference.main(new String[]{"in="+path("queries.tsv"), "reference="+path("reference.tsv"),
			"familymap="+path("grade_map.tsv"), "restrictreference=t", "subtypes="+path("grade_subtypes.tsv"),
			"witnesses="+path("witnesses.tsv"), "events="+path("events.tsv"), "out="+path("graded")});
		for(final String view:new String[]{"RECOVERY", "EXCLUSIVE"}){
			for(final String model:new String[]{"RF00174", "RF01482", "RF01689"}){
				final String[] inventory=row("result.summary.tsv", view+"\tcobalamin\t"+model+"\t").split("\t");
				final String nativeView=view.equals("RECOVERY") ? "MODEL_SUPPORT_RECOVERY" : "EXCLUSIVE_SUBTYPE_ACCURACY";
				final String[] grade=row("graded.subtype_summary.tsv", nativeView+"\t"+model+"\t").split("\t");
				check(inventory[7].equals(grade[6]) && inventory[8].equals(grade[7]),
					"Inventory reference denominators must equal the existing grader in "+view+"/"+model);
			}
		}
	}

	private void inventory(final String output, final String reference, final String mapping, final String witness){
		BuildNcrnaReference.main(new String[]{"mode=inventory", "census="+path("census.tsv"),
			"memberships="+path("memberships.tsv"), "reference="+path(reference), "inventorymembers="+path(mapping),
			"witnesses="+path(witness), "events="+path("events.tsv"), "out="+path(output), "t=1"});
	}
	private void reject(final String output, final String reference, final String mapping, final String witness, final String expected){
		try{
			inventory(output, reference, mapping, witness);
			throw new AssertionError("Invalid inventory succeeded: "+output);
		}catch(IllegalArgumentException error){
			check(error.getMessage().contains(expected), "Wrong rejection: "+error);
		}
		check(!new File(path(output+".provenance.tsv")).exists(), "Failed inventory cannot publish its final completion record");
	}
	private Row parent(final ArrayList<Row> rows, final String input, final String family, final String model,
		final int start, final boolean mask, final char strand){
		final Row row=codec.create(input, new Locus("seq", family, start, start+99, strand), mask,
			mask ? "retained_partial" : null, model, path("witnesses.tsv"));
		rows.add(row);
		return row;
	}
	private void hit(final Row row, final String model, final String tier, final String association, final int priorLine){
		final int id=++nextEvent;
		events.append(id).tab().append(row.input).append("\tseq\t").append(row.locus.start).tab().append(row.locus.stop)
			.tab().append(row.locus.strand).tab().append(model).tab().append(row.locus.family)
			.append("\t1\t100\t100\tno\tfalse\tfalse\t50.0\t40.0\t125.000000\ttrue\tfalse\t90\ttrue\tfalse\t")
			.append(tier).append("\tsynthetic\t.\t").append(association.equals(BRIDGE) ? 16 : 0)
			.append("\t/synthetic/native.tbl\t").append(id+1).append("\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n");
		BuildNcrnaReference.field(witnesses, row.id).tab().append(row.input).append("\tCM_HIT\t").append(id)
			.append("\t/synthetic/native.tbl\t").append(id+1).tab().append(model).tab().append(row.locus.start)
			.tab().append(row.locus.stop).tab().append(priorLine).tab().append(association).nl();
	}
	private void prior(final Row row){
		BuildNcrnaReference.field(witnesses, row.id).tab().append(row.input).append("\tPROTECTED_PRIOR\t.\t.\t.\t")
			.append(row.model).tab().append(row.locus.start).tab().append(row.locus.stop)
			.append("\t9\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n");
	}
	private void membership(final ByteBuilder out, final String source, final int index, final int input){
		out.append(source).tab().append(index).tab().append(input).tab().append(path("unopened"+input+".fa"))
			.append("\tgenome").append(input).append("\tNA\tNA\tprokaryotic\tNA\tNA\n");
	}
	private void write(final String name, final String value){
		try(Output out=new Output(path(name))){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}
	}
	private ArrayList<String> lines(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				result.add(new String(line, StandardCharsets.UTF_8));
			}
		}
		return result;
	}
	private String row(final String name, final String prefix){
		final String found=find(lines(name), prefix);
		check(found!=null, "Missing row "+prefix+" in "+name);
		return found;
	}
	private static String find(final ArrayList<String> rows, final String prefix){
		String found=null;
		for(final String row:rows){
			if(row.startsWith(prefix)){
				check(found==null, "Repeated fixture row "+prefix);
				found=row;
			}
		}
		return found;
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean ok, final String why){
		if(!ok){throw new AssertionError(why);}
	}
	private final File root;
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private int nextEvent;
	private final ByteBuilder events=new ByteBuilder("event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision\n");
	private final ByteBuilder witnesses=new ByteBuilder("annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation\n");
	private static final String GROUP="QUALIFYING_SAME_FAMILY_GROUP_MEMBER", PRIOR="SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS", BRIDGE="TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE";
}
