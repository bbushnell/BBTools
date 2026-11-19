package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareNcrnaSeedInputs;
import shared.Shared;
import structures.ByteBuilder;

/** Protects exact-source requests and the all-original denominator without DNA I/O.
 * Uses deliberately nonexistent FASTA paths so the metadata-only contract is tested.
 * Adapted from template.A_SampleBasic2 and PrepareNcrnaRosterBaselineTest.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaSeedInputsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PrepareNcrnaSeedInputsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown seed-input fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory without input or overwrite");
			new PrepareNcrnaSeedInputsTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_SEED_INPUT_PREPARATION_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixture            ----------------*/
	/*--------------------------------------------------------------*/

	private PrepareNcrnaSeedInputsTest(final File root_){
		root=root_;
		check(root.mkdir(), "Preserve fixture evidence by requiring a fresh directory");
	}

	private void run(){
		final String[][] originals=genes();
		writePopulation("accepted", originals, false, 2);
		prepare("accepted", "accepted-plan");
		try(Input in=new Input(path("accepted-plan.requests.tsv"))){
			check(text(in.reader.nextLine()).equals(PrepareNcrnaSeedInputs.REQUEST_HEADER), "Request schema changed");
			final String expected="2\t2\tframe_1\t8\t"+path("unread-missing.fa")
				+"\tsame\t5\t19\t-\t1\ttrue\tAVAILABLE_RECORDED_FRAME";
			check(text(in.reader.nextLine()).equals(expected),
				"Retain original ID, exact source, minus strand and inclusive bounds; do not retrieve the eligible duplicate");
			check(in.reader.nextLine()==null, "Only the excluded original needs new DNA");
		}
		try(Input in=new Input(path("accepted-plan.sources.tsv"))){
			in.reader.nextLine();
			check(text(in.reader.nextLine()).equals("frame_1\t8\t"+path("unread-missing.fa")
				+"\t500\t80\t1\t1\tAVAILABLE_RECORDED_FRAME"),
				"Cost only the exact missing-gene source, not all frames or an assembly alias");
			check(in.reader.nextLine()==null, "One missing source must be read once");
		}
		check(metric("accepted-plan", "original_genes")==3
			&& metric("accepted-plan", "recorded_source_frames")==3
			&& metric("accepted-plan", "originals_with_saved_sequences")==2
			&& metric("accepted-plan", "distinct_saved_representatives")==1
			&& metric("accepted-plan", "missing_original_sequences")==1
			&& metric("accepted-plan", "requested_evaluation_bases")==15
			&& metric("accepted-plan", "planned_source_file_bytes")==500
			&& metric("accepted-plan", "planned_source_nucleotide_bases")==80,
			"Keep originals, construction representatives, empty frames and requested bases as distinct units");

		final String[][] unavailable=genes();
		unavailable[2][30]="UNRESOLVED_SOURCE_UNAVAILABLE";
		writePopulation("unavailable", unavailable, true, 2);
		prepare("unavailable", "unavailable-plan");
		check(metric("unavailable-plan", "recorded_unavailable_missing_sources")==1
			&& metric("unavailable-plan", "missing_original_sequences")==1,
			"An unavailable source remains a request, never a removed recall gene");

		final String[][] duplicate=genes();
		duplicate[1][0]="0";
		reject("duplicate", duplicate, 2, "component merges or duplicate IDs");
		final String[][] merged=genes();
		merged[1][13]="0";
		reject("collapsed", merged, 2, "component merges or duplicate IDs");
		final String[][] crossGroup=genes();
		crossGroup[1][1]="8";
		crossGroup[1][2]=path("unread-missing.fa");
		crossGroup[1][29]="frame_1";
		reject("cross-group", crossGroup, 1, "SAME accepted group");
		final String[][] indirect=genes();
		indirect[0][14]="1";
		indirect[0][15]="EXACT_SEQUENCE_DUPLICATE_FOR_CONSTRUCTION_ONLY";
		indirect[1][14]="0";
		reject("cycle", indirect, 2, "saved-member counts");
		final String[][] notTruncated=genes();
		notTruncated[2][17]="0";
		reject("not-truncated", notTruncated, 2, "Only resolved excluded truncations");
		final String[][] wrongFrame=genes();
		wrongFrame[2][29]="frame_0";
		reject("wrong-frame", wrongFrame, 2, "exact canonical frame");
		reject("wrong-count", genes(), 1, "counts must conserve");
		boolean fresh=false;
		try{prepare("accepted", "accepted-plan");}
		catch(IllegalArgumentException expected){fresh=expected.getMessage().contains("must be fresh");}
		check(fresh, "A second invocation must not overwrite the accepted requests");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input Builders       ----------------*/
	/*--------------------------------------------------------------*/

	private String[][] genes(){
		return new String[][]{
			gene(0, 5, 0, 0, "CONSTRUCTION_REPRESENTATIVE", 0),
			gene(1, 5, 0, 0, "EXACT_SEQUENCE_DUPLICATE_FOR_CONSTRUCTION_ONLY", 0),
			gene(2, 8, 1, -1, "EXCLUDED_NATIVE_PARTIAL_OR_TRUNCATED", 1)
		};
	}

	private String[] gene(final int id, final int group, final int frame, final int representative,
		final String state, final int truncation){
		return new String[]{""+id, ""+group, path(frame==0 ? "unread-eligible.fa" : "unread-missing.fa"),
			"same", "5", "19", id==2 ? "-" : "+", "true", "false", "1", "0", "0",
			"QUALIFIED", ""+id, ""+representative, state, "1", ""+truncation, "1", "15", "15",
			"false", truncation>0 ? "NATIVE_TRUNCATION" : "ELIGIBLE_NATIVE_HIT", "0", "0", "0", "0",
			"CONTIG_END_REACH", "50.0", "frame_"+frame, "AVAILABLE_RECORDED_FRAME"};
	}

	private void writePopulation(final String prefix, final String[][] genes, final boolean unavailable,
		final int firstFrameGenes){
		final ByteBuilder frames=new ByteBuilder(PrepareNcrnaSeedInputs.FRAME_HEADER).nl();
		frames.append("frame_0\t5\t").append(path("unread-eligible.fa")).tab().append(path("unread-eligible.fa"))
			.append("\t1000\t100\t1\t").append(firstFrameGenes).append("\tAVAILABLE_RECORDED_FRAME\n");
		frames.append("frame_1\t8\t").append(path("unread-missing.fa")).tab().append(path("unread-missing.fa"))
			.append("\t500\t80\t1\t").append(3-firstFrameGenes).tab()
			.append(unavailable ? "UNRESOLVED_SOURCE_UNAVAILABLE" : "AVAILABLE_RECORDED_FRAME").nl();
		frames.append("frame_2\t10\t").append(path("unread-zero.fa")).tab().append(path("unread-zero.fa"))
			.append("\t200\t60\t1\t0\tUNRESOLVED_SOURCE_UNAVAILABLE\n");
		write(prefix+".frames.tsv", frames);
		final ByteBuilder rows=new ByteBuilder(PrepareNcrnaSeedInputs.ORIGINAL_HEADER).nl();
		for(final String[] gene:genes){
			for(int i=0; i<gene.length; i++){
				if(i>0){rows.tab();}
				rows.append(gene[i]);
			}
			rows.nl();
		}
		write(prefix+".original_calls.tsv.gz", rows);
	}

	private void prepare(final String input, final String output){
		PrepareNcrnaSeedInputs.main(new String[]{"in="+path(input), "out="+path(output),
			"expectedgenes=3", "expectedframes=3", "expectedmembers=1", "expectedmissing=1"});
	}

	private void reject(final String name, final String[][] genes, final int firstFrameGenes, final String reason){
		writePopulation(name, genes, false, firstFrameGenes);
		boolean rejected=false;
		try{prepare(name, name+"-plan");}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(reason);}
		check(rejected, "Reject malformed source/representative metadata before requesting DNA: "+name+" "+reason);
		check(!new File(path(name+"-plan.requests.tsv")).exists(),
			"Invalid metadata must not publish partial retrieval requests: "+name);
	}

	private long metric(final String prefix, final String key){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path(prefix+".summary.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals(key, 0)){return p.parseLong(1);}
			}
		}
		throw new AssertionError("Missing required denominator/cost metric: "+key);
	}

	private String path(final String name){return new File(root, name).getPath();}
	private static String text(final byte[] bytes){return bytes==null ? "" : new String(bytes, StandardCharsets.UTF_8);}
	private void write(final String name, final ByteBuilder bytes){try(Output out=new Output(path(name))){out.writer.print(bytes);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
}
