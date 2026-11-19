package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import fileIO.FileFormat;
import dna.Data;
import prok.CobalaminSelectedPopulation;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaProcess;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Source collisions, duplicate native events, truncations, masked/prior rows and
 * coordinate disagreement must remain distinct. Uses three synthetic inputs;
 * their deliberately invalid DNA payloads ensure this metadata tool reads no DNA.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSelectedPopulationTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CobalaminSelectedPopulationTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String mode="full";
			int expected=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("mode")){mode=b;}
				else if(a.equals("expectedgenomes")){expected=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown fixture argument: "+arg);}
			}
			require(parser.out1!=null && !parser.overwrite && !parser.append, "Supply out=fresh_fixture_directory");
			if(mode.equals("missingstate")){
				require(parser.in1!=null && expected>0, "Missing-state control needs in=dropped_ledger expectedgenomes=N");
				verifyDroppedLedger(parser.in1, expected);
				throw new AssertionError("A dropped input must fail the production conservation guard");
			}
			require(mode.equals("full"), "Unknown fixture mode: "+mode);
			fixture(new File(parser.out1));
			parallel(new File(parser.out1, "parallel"));
			orderedGap(new File(parser.out1, "missing_job.tsv"));
			producerFailure(new File(parser.out1, "failed_producer"));
			pp.outstream.println("COBALAMIN_SELECTED_POPULATION_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** The native writer must finish and throw for a deliberately missing job1. */
	private static void orderedGap(final File path){
		final Output output=new Output(FileFormat.testOutput(path.getAbsolutePath(), FileFormat.TEXT, null, false, false, false, true));
		boolean failed=false;
		try(Output owner=output){
			owner.writer.add(new ByteBuilder().append("header\n"), 0);
			owner.writer.add(new ByteBuilder().append("after_gap\n"), 2);
		}catch(RuntimeException expected){failed=expected.getMessage().contains("jobID 1 was never supplied");}
		require(failed && output.writer.getState()==Thread.State.TERMINATED,
			"A missing ordered ID must fail loudly after the writer has terminated");
	}

	/** One bad second input precedes more valid empty jobs than the writer's128-job
	 * flow-control threshold. The failed worker must fill its ID and drain, so
	 * other producers can finish and the error reaches the caller without a hang. */
	private static void producerFailure(final File root) throws Exception{
		require(root.mkdir(), "Fresh producer-failure fixture required");
		final String path=root.getAbsolutePath()+"/", fasta=path+"source.fa";
		write(fasta, "Not DNA.\n");
		final String header="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope";
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		try(Output selected=new Output(path+"selected.tsv"); Output provenance=new Output(path+"provenance.tsv"); Output reference=new Output(path+"reference.tsv")){
			selected.writer.println(header+"\tselection_role");
			provenance.writer.println(header);
			reference.writer.println("input_id\tfasta\tgff");
			for(int input=1; input<=512; input++){
				final String gff=path+input+".gff";
				if(input==2){write(gff, "Not a GFF header.\n");}
				else{codec.write(gff, Integer.toString(input), new ArrayList<Row>());}
				final String row="assembly_"+input+"\t"+input+"\tfixture\t"+fasta+"\t"+gff+"\t.\t.\t"+fasta+"\t"+gff+"\tfasta\tfixture";
				selected.writer.println(row+"\tretained");
				provenance.writer.println(row);
				reference.writer.println(input+"\t"+fasta+"\t"+gff);
			}
		}
		write(path+"events.tsv", EVENT_HEADER+"\n");
		write(path+"witnesses.tsv", "annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation\n");
		write(path+"cohorts.tsv", "cohort\treference\tevents\twitnesses\nfixture\t"+path+"reference.tsv\t"+path+"events.tsv\t"+path+"witnesses.tsv\n");
		boolean failed=false;
		try{
			CobalaminSelectedPopulation.main(new String[]{"in="+path+"selected.tsv", "provenance="+path+"provenance.tsv", "cohorts="+path+"cohorts.tsv", "out="+path+"result", "expectedgenomes=512", "t=8"});
		}catch(IllegalStateException expected){failed=expected.getMessage().startsWith("Reference producers failed:")
			&& expected.getMessage().contains("input ordinal0=1 assembly=assembly_2");}
		require(failed, "A parallel producer error must throw with its exact input ordinal and reason");
		verifyFailureLedger(new File(path+"result"), false);
		final String launcher=Data.findPath("?dev/cobalaminselectedpopulation.sh", false);
		require(launcher!=null, "Standalone negative control requires the actual installed launcher");
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(launcher);
		command.add("in="+path+"selected.tsv");
		command.add("provenance="+path+"provenance.tsv");
		command.add("cohorts="+path+"cohorts.tsv");
		command.add("out="+path+"result_serial");
		command.add("expectedgenomes=512");
		command.add("t=1");
		NcrnaProcess.addJvmOptions(command, 1024);
		boolean nonzero=false;
		try{NcrnaProcess.run(command, null, "Expected population failure");}
		catch(IllegalStateException expected){nonzero=expected.getMessage().equals("Expected population failure exited 1");}
		require(nonzero, "Standalone invocation must exit1 for a recorded producer failure");
		write(path+"standalone.exit", "1\n");
		verifyFailureLedger(new File(path+"result_serial"), true);
		try(Input in=new Input(path+"result_serial/producer_status.tsv"); Output out=new Output(path+"dropped_status.tsv")){
			out.writer.println(in.reader.nextLine());
			for(int index=0; index<511; index++){out.writer.println(in.reader.nextLine());}
			require(in.reader.nextLine()!=null && in.reader.nextLine()==null, "Drop exactly one of the512producer-state rows");
		}
		final String testLauncher=Data.findPath("?dev/cobalaminselectedpopulation_test.sh", false);
		require(testLauncher!=null, "Dropped-row control requires the installed test launcher");
		command.clear();
		command.add("bash");
		command.add(testLauncher);
		command.add("mode=missingstate");
		command.add("in="+path+"dropped_status.tsv");
		command.add("out="+path+"unused");
		command.add("expectedgenomes=512");
		NcrnaProcess.addJvmOptions(command, 1024);
		boolean conservationFailed=false;
		try{NcrnaProcess.run(command, null, "Expected state conservation failure");}
		catch(IllegalStateException expected){conservationFailed=expected.getMessage().equals("Expected state conservation failure exited 1");}
		require(conservationFailed, "A dropped ledger row must exit1 with the production conservation mismatch");
		write(path+"missing_state.exit", "1\n");
	}

	/** Reads an actually shortened ledger, then invokes the same always-on guard
	 * used before the producer writes any count. Independent oracle:511versus512. */
	private static void verifyDroppedLedger(final String path, final int expected){
		final long[] counts=new long[4];
		final String[] states={"COMPLETE_EMPTY_OUTPUT", "COMPLETE_NONEMPTY_OUTPUT", "FAILED", "SKIPPED_AFTER_FAILURE"};
		final LineParser1 parser=new LineParser1('\t');
		try(Input in=new Input(path)){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				parser.set(line);
				require(parser.terms()==7, "Dropped-state fixture must retain complete remaining records");
				boolean found=false;
				for(int state=0; state<states.length; state++){
					if(parser.termEquals(states[state], 4)){
						counts[state]++;
						found=true;
					}
				}
				require(found, "Every retained fixture row has one declared terminal state");
			}
		}
		CobalaminSelectedPopulation.validateProducerCounts(counts[0], counts[1], counts[2], counts[3], expected);
	}

	/** The deterministic serial control has one genuine empty, one failed and510
	 * skipped inputs. The parallel control may distribute successful/skipped work
	 * differently, but must retain the same failed input and all512status rows. */
	private static void verifyFailureLedger(final File root, final boolean serial){
		for(final String name:new String[]{"summary.tsv", "inputs.tsv", "native_loci.tsv.gz", "native_events.tsv.gz", "result.tsv"}){
			require(!new File(root, name).exists(), "A producer failure must precede every population-count file: "+name);
		}
		final LineParser1 parser=new LineParser1('\t');
		int rows=0, empty=0, failed=0, skipped=0;
		try(Input in=new Input(new File(root, "producer_status.tsv").getPath())){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				parser.set(line);
				require(parser.terms()==7 && parser.parseInt(0)==rows && parser.parseInt(1)==rows+1,
					"Every input needs its original ordinal and exact ordered job identity");
				if(parser.termEquals("FAILED", 4)){
					require(rows==1 && parser.termEquals("assembly_2", 2) && parser.parseInt(5)==-1
						&& parser.parseString(6).contains("version directive"), "A failed producer must retain its exact input/reason, never measured zero bytes");
					failed++;
				}else if(parser.termEquals("SKIPPED_AFTER_FAILURE", 4)){
					require(rows>1 && parser.parseInt(5)==-1 && parser.termEquals("after_failed_input_ordinal0=1", 6),
						"Unprocessed inputs must retain skipped state and causing failure, never biological zero");
					skipped++;
				}else{
					require(parser.termEquals("COMPLETE_EMPTY_OUTPUT", 4) && parser.parseInt(5)==0 && parser.termEquals(".", 6),
						"A legitimate empty output requires successful complete processing");
					empty++;
				}
				rows++;
			}
		}
		require(rows==512 && failed==1 && empty>=1 && empty+failed+skipped==rows,
			"Successful empty, failed and skipped inputs must partition all512 original inputs");
		require(!serial || empty==1 && skipped==510, "Serial control must expose all three distinct categories deterministically");
	}

	/** Every source reuses the same contig/coordinates. Parallel evidence must
	 * retain all2048distinct source/event pairs across multiple writer buffers. */
	private static void parallel(final File root) throws Exception{
		require(root.mkdir(), "Fresh parallel fixture directory required");
		final String path=root.getAbsolutePath()+"/";
		final String header="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope";
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		try(Output selected=new Output(path+"selected.tsv"); Output provenance=new Output(path+"provenance.tsv");
			Output reference=new Output(path+"reference.tsv"); Output events=new Output(path+"events.tsv"); Output witnesses=new Output(path+"witnesses.tsv")){
			selected.writer.println(header+"\tselection_role");
			provenance.writer.println(header);
			reference.writer.println("input_id\tfasta\tgff");
			events.writer.println(EVENT_HEADER);
			witnesses.writer.println("annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation");
			for(int source=0; source<16; source++){
				final int input=100+source;
				final String fasta=path+input+".fa", gff=path+input+".gff";
				write(fasta, "Not DNA.\n");
				final String row="assembly_"+input+"\t"+input+"\tfixture\t"+fasta+"\t"+gff+"\t.\t.\t"+fasta+"\t"+gff+"\tfasta\tfixture";
				selected.writer.println(row+"\tretained");
				provenance.writer.println(row);
				reference.writer.println(input+"\t"+fasta+"\t"+gff);
				final ArrayList<Row> rows=new ArrayList<Row>();
				for(int locus=0; locus<128; locus++){
					final int start=10+locus*20, event=1+source*128+locus;
					final Row parent=codec.create(Integer.toString(input), new Locus("shared", "cobalamin", start, start+10, '+'), false, null, "RF00174", "fixture");
					rows.add(parent);
					events.writer.print(event(event, input, "RF00174", start, start+10, 0, true));
					witnesses.writer.print(witness(parent, event, "RF00174", start, start+10));
				}
				codec.write(gff, Integer.toString(input), rows);
			}
		}
		write(path+"cohorts.tsv", "cohort\treference\tevents\twitnesses\nfixture\t"+path+"reference.tsv\t"+path+"events.tsv\t"+path+"witnesses.tsv\n");
		CobalaminSelectedPopulation.main(new String[]{"in="+path+"selected.tsv", "provenance="+path+"provenance.tsv", "cohorts="+path+"cohorts.tsv", "out="+path+"result", "expectedgenomes=16", "t=8"});
		final boolean[] seen=new boolean[2049];
		final LineParser1 parser=new LineParser1('\t');
		int count=0;
		try(Input in=new Input(path+"result/raw_reference_links.tsv.gz")){
			in.reader.nextLine();
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				parser.set(line);
				require(parser.terms()==16, "Parallel correspondence rows must remain complete");
				final int id=parser.parseInt(3);
				require(id>0 && id<seen.length && !seen[id], "No native correspondence may disappear or be duplicated");
				require(id==count+1, "Native ordered jobs must preserve complete selected-input and within-input order");
				seen[id]=true;
				final int input=100+(id-1)/128, start=10+((id-1)%128)*20;
				require(parser.parseInt(1)==input && parser.termEquals("assembly_"+input, 0)
					&& parser.parseInt(8)==start && parser.parseInt(9)==start+10
					&& parser.parseInt(12)==start && parser.parseInt(13)==start+10
					&& parser.termEquals("true", 15), "A shared contig name cannot mix coordinates or sources across worker buffers");
				count++;
			}
		}
		require(count==2048, "Parallel output must conserve the complete source/event multiset");
	}

	private static void fixture(final File root) throws Exception{
		require(root.mkdir(), "Fresh fixture directory required");
		final String path=root.getAbsolutePath()+"/";
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> first=new ArrayList<Row>();
		first.add(codec.create("1", new Locus("shared", "cobalamin", 10, 20, '+'), false, null, "RF00174", "fixture"));
		first.add(codec.create("1", new Locus("shared", "cobalamin", 30, 45, '+'), false, null, "RF00174", "fixture"));
		first.add(codec.create("1", new Locus("shared", "cobalamin", 100, 130, '+'), false, null, "RF00174", "fixture"));
		first.add(codec.create("1", new Locus("shared", "cobalamin", 150, 160, '+'), true, "fixture_mask", "RF00174", "fixture"));
		codec.write(path+"a.gff", "1", first);
		final ArrayList<Row> second=new ArrayList<Row>();
		second.add(codec.create("2", new Locus("shared", "cobalamin", 10, 20, '+'), false, null, "RF01482", "fixture"));
		codec.write(path+"b.gff", "2", second);
		codec.write(path+"c.gff", "3", new ArrayList<Row>());
		for(final String name:new String[]{"a", "b", "c"}){write(path+name+".fa", "This is not a FASTA file.\n");}
		final String header="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role\n";
		try(Output out=new Output(path+"manifest.tsv"); Output provenance=new Output(path+"provenance.tsv")){
			out.writer.print(header);
			provenance.writer.print(header.substring(0, header.lastIndexOf('\t'))+"\n");
			for(final String name:new String[]{"a", "c", "b"}){
				final int id=name.equals("a") ? 1 : name.equals("b") ? 2 : 3;
				final String fields=path+name+".fa\t"+path+name+".gff\t.\t.\t"+path+name+".fa\t"+path+name+".gff\tfasta\tfixture";
				out.writer.print(name+"\t"+id+"\t"+(id==2 ? "older_ANI95_addition" : "fixture")+"\t"+fields+"\tretained\n");
				provenance.writer.print(name+"\t"+id+"\tfixture\t"+fields+"\n");
			}
		}
		write(path+"reference.tsv", "input_id\tfasta\tgff\n1\t"+path+"a.fa\t"+path+"a.gff\n2\t"+path+"b.fa\t"+path+"b.gff\n3\t"+path+"c.fa\t"+path+"c.gff\n");
		write(path+"cohorts.tsv", "cohort\treference\tevents\twitnesses\nfixture\t"+path+"reference.tsv\t"+path+"events.tsv\t"+path+"witnesses.tsv\n");
		write(path+"events.tsv", EVENT_HEADER+"\n"
			+event(1, 1, "RF00174", 10, 20, 0, true)
			+event(2, 1, "RF00174", 10, 20, 1, true)
			+event(3, 1, "RF00174", 30, 40, 1, true)
			+event(4, 2, "RF01482", 10, 20, 0, true)
			+event(5, 2, "RF01689", 50, 60, 0, false));
		write(path+"witnesses.tsv", "annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation\n"
			+witness(first.get(0), 1, "RF00174", 10, 20)
			+witness(first.get(1), 3, "RF00174", 30, 40)
			+prior(first.get(2))+prior(first.get(3))
			+witness(second.get(0), 4, "RF01482", 10, 20));
		CobalaminSelectedPopulation.main(new String[]{"in="+path+"manifest.tsv", "provenance="+path+"provenance.tsv", "cohorts="+path+"cohorts.tsv", "out="+path+"result", "expectedgenomes=3", "t=2"});
		try(Input in=new Input(path+"result/summary.tsv")){
			in.reader.nextLine();
			final String[] expected={"RF00174\t3\t3\t2\t1\t1\t1", "RF01482\t3\t1\t1\t1\t0\t1", "RF01689\t3\t0\t0\t0\t0\t0"};
			for(final String row:expected){require(row.equals(new String(in.reader.nextLine(), "UTF-8")), "Native events/loci/truncations and distinct sources must conserve fixture counts");}
			require(in.reader.nextLine()==null, "Three members only");
		}
		try(Input in=new Input(path+"result/reference_categories.tsv")){
			in.reader.nextLine();
			require("cobalamin_union\t3\t1\t1\tREFERENCE_PARENTS_NOT_RAW_NATIVE_CALLS".equals(new String(in.reader.nextLine(), "UTF-8")), "Prior and mask labels must not select native population members");
		}
		try(Input in=new Input(path+"result/result.tsv")){
			in.reader.nextLine();
			require("raw_scored_witness_links\t3".equals(new String(in.reader.nextLine(), "UTF-8")), "Every native/scored association retained");
			require("nonidentical_raw_scored_links\t1".equals(new String(in.reader.nextLine(), "UTF-8")), "Changed GFF coordinates must stop construction");
		}
		int rows=0;
		try(Input in=new Input(path+"result/inputs.tsv")){
			in.reader.nextLine();
			while(in.reader.nextLine()!=null){rows++;}
		}
		require(rows==9, "Every input, including the zero-hit genome, must have every member row");
		boolean failed=false;
		try{
			CobalaminSelectedPopulation.main(new String[]{"in="+path+"manifest.tsv", "provenance="+path+"provenance.tsv", "cohorts="+path+"cohorts.tsv", "out="+path+"bad", "expectedgenomes=4"});
		}catch(IllegalArgumentException expected){failed=true;}
		require(failed && !new File(path+"bad").exists(), "A missing selected genome must fail before publication");
	}

	private static String event(final int id, final int input, final String model, final int start, final int stop, final int truncation, final boolean nativeGa){
		return id+"\t"+input+"\tshared\t"+start+"\t"+stop+"\t+\t"+model+"\tcobalamin\t1\t11\t11\t"+truncation
			+"\tfalse\tfalse\t100.0\t50\t200\t"+nativeGa+"\tfalse\t100\ttrue\tfalse\tCONFIDENT\tfixture\t.\t0\t/saved/native.tbl\t"+id+"\tfixture\tfixture\n";
	}
	private static String witness(final Row row, final int event, final String model, final int start, final int stop){
		return row.id+"\t"+row.input+"\tCM_HIT\t"+event+"\t/saved/native.tbl\t"+event+"\t"+model+"\t"+start+"\t"+stop+"\t-1\tQUALIFYING_SAME_FAMILY_GROUP_MEMBER\n";
	}
	private static String prior(final Row row){
		return row.id+"\t"+row.input+"\tPROTECTED_PRIOR\t.\t.\t.\tRF00174\t"+row.locus.start+"\t"+row.locus.stop+"\t1\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n";
	}
	private static void write(final String path, final String value){
		try(Output out=new Output(path)){out.writer.print(value);}
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
}
