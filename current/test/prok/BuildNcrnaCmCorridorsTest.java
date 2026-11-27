package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import aligner.CovarianceModel;
import aligner.CovarianceModelParser;
import aligner.CovarianceModelStateCorridors;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaCmCorridors;
import prok.NcrnaCmCorridorResource;
import prok.PrepareNcrnaCmTraining;
import prok.TraceNcrnaCmTraining;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal saved paths, not another alignment run.102finite searches contain
 *101usable anchors; one greatest-departure trajectory is trimmed, leaving100
 *from10resolved groups. A separate raw member remains on the original axis.
 *All cross-validation fits are sparse and use their explicit uniform fallback.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaCmCorridorsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BuildNcrnaCmCorridorsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Use out=fresh_fixture_directory");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(),
				"Literal corridor fixtures require fresh owned data files");
			test.aligner.CovarianceModelCykTest.main(new String[0]);
			test.aligner.CovarianceModelStateCorridorsTest.main(new String[0]);
			NcrnaCmCorridorSelectionTest.main(new String[0]);
			run(new File(parser.out1).getAbsolutePath());
			System.err.println("NCRNA_CM_CORRIDOR_FIT_FIXTURE_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixture            ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final String root){
		write(root+"/literal.cm", "INFERNAL1/a\nNAME LiteralAC\nACC RF00174\nSTATES 3\nNODES 3\nCLEN 2\nW 6\nALPH RNA\nNULL 0 0 0 0\nELSELF -1\nEFP7GF 0 1\nGA 100\nQDBBETA1 1e-7\nQDBBETA2 1e-15\nCM\n"
			+"[ ROOT 0 ] - - - - - -\nS 0 -1 0 1 1 0 0 6 6 0\n"
			+"[ MATP 1 ] - - - - - -\nMP 1 0 1 2 1 0 0 6 6 0 * 0 * * * * * * * * * * * * * *\n"
			+"[ END 2 ] - - - - - -\nE 2 1 1 -1 0 0 0 0 0\n//\n");
		final ByteBuilder audit=new ByteBuilder(AUDIT_HEADER).nl();
		for(int i=0; i<21; i++){
			final int group=i==20 ? 0 : i/2;
			audit.append('a').append(i).tab().append('i').append(i).append("\t/unused").append(i).append(".fa\t")
				.append(group).tab().append(group).append("\tfalse\tfalse\tfalse\n");
		}
		write(root+"/audit.tsv", audit.toString());
		write(root+"/models.tsv", MODEL_HEADER+"\n0\tRF00174\tliteral\t102\t10\ttrue\n1\tRAW\traw\t1\t1\tfalse\n");
		final ByteBuilder original=new ByteBuilder(PrepareNcrnaCmTraining.CALL_HEADER).nl();
		final ByteBuilder summaries=new ByteBuilder(TraceNcrnaCmTraining.CALL_HEADER).nl();
		final ByteBuilder paths=new ByteBuilder(TraceNcrnaCmTraining.PATH_HEADER).nl();
		for(int call=0; call<103; call++){
			final int source=call%20;
			original.append(call).append("\ta").append(source).tab().append(call==102 ? 1 : 0).append("\ti").append(source)
				.append("\tC%2525\t2\t5\t+\t").append(call==102 ? "RAW\traw" : "RF00174\tliteral n=4")
				.append("\tTP\ttrue\tr").append(call).nl();
			if(call==102){continue;}
			final boolean anchored=call!=101;
			final int shift=call==0 ? 1 : 0;
			summaries.append(call).append("\ta").append(source).append("\t0\ti").append(source)
				.append("\tC%2525\t2\t5\t+\t1\t6\t6\tOK\t4\t").append(3+shift).tab().append(4+shift)
				.tab().append(anchored).tab().append(anchored ? shift : -1).append("\t1\t1\t1\t147\tTP\tr").append(call).nl();
			for(int state=0; state<3; state++){
				paths.append(call).append("\t0\t").append(state).tab().append(state).tab().append((state==2 ? 3 : 2)+shift)
					.tab().append((state==2 ? 3 : 4)+shift).tab().append(state==2 ? -1 : 0)
					.tab().append(state==2 ? -1 : state+1).append("\t-1");
				if(anchored){
					paths.tab().append(state==2 ? 3 : 2).tab().append(state==2 ? 3 : 2)
						.tab().append(state==2 ? 3 : 4).tab().append(state==2 ? 3 : 4);
				}else{paths.append("\t\t\t\t");}
				paths.nl();
			}
		}
		write(root+"/original.tsv", original.toString());
		write(root+"/calls.tsv.gz", summaries.toString());
		write(root+"/paths.tsv.gz", paths.toString());
		write(root+"/parts.tsv", "calls\tpaths\n"+root+"/calls.tsv.gz\t"+root+"/paths.tsv.gz\n");
		invoke(root, root+"/audit.tsv", root+"/fit");
		checkSummary(root+"/fit/summary.tsv");
		checkContainment(root+"/fit/containment.tsv");
		final CovarianceModel cm=CovarianceModelParser.read(root+"/literal.cm");
		final String[] names={"literal", "raw"}, members={"RF00174", "RAW"};
		final CovarianceModel[] models={cm, null};
		final NcrnaCmCorridorResource.Data loaded=NcrnaCmCorridorResource.read(root+"/fit/RF00174.corridors.tsv.gz",
			new String[]{"literal"}, new String[]{"RF00174"}, new CovarianceModel[]{cm});
		check(loaded.retainedCalls[0]==100 && loaded.resolvedGroups[0]==10 && loaded.profiles.length==1,
			"The public reader must retain100calls/10groups on exactly the selected member's axis");
		check(!new File(root+"/fit/RAW.corridors.tsv.gz").exists(), "An unrefined member remains in reports without an invented fitted resource");
		for(int state=0; state<3; state++){
			check(loaded.profiles[0].observed(state), "The hundred retained literal paths visit all three grammar states");
			for(int bound=0; bound<4; bound++){check(loaded.profiles[0].offset(state, bound)==0, "The removed shifted trajectory must not widen retained zero-residual bounds");}
		}
		invoke(root, root+"/audit.tsv", root+"/explicit_default", "trimpercent=1");
		check(text(root+"/fit/RF00174.corridors.tsv.gz").equals(text(root+"/explicit_default/RF00174.corridors.tsv.gz")),
			"Omitted and explicit default must produce identical serving resources, including fitted support");
		invoke(root, root+"/audit.tsv", root+"/five_percent", "trimpercent=5", "mincalls=95");
		final NcrnaCmCorridorResource.Data trimmed=NcrnaCmCorridorResource.read(root+"/five_percent/RF00174.corridors.tsv.gz",
			new String[]{"literal"}, new String[]{"RF00174"}, new CovarianceModel[]{cm});
		check(trimmed.retainedCalls[0]==96 && trimmed.resolvedGroups[0]==10,
			"Explicit five-percent CLI must reach the selector and serving resource: floor101*5/100 removes five complete paths");
		boolean invalidTrimRejected=false;
		try{invoke(root, root+"/audit.tsv", root+"/invalid_trim", "trimpercent=100");}
		catch(IllegalArgumentException expected){invalidTrimRejected=true;}
		check(invalidTrimRejected && !new File(root+"/invalid_trim").exists(),
			"Invalid CLI trimming must fail before creating output or reading training streams");
		// The unused alias20belongs to group0. Excluding it must also exclude
		// used source0; a check of only each used row's flags would miss this.
		final String badAudit=audit.toString().replace("a20\ti20\t/unused20.fa\t0\t0\tfalse", "a20\ti20\t/unused20.fa\t0\t0\ttrue");
		write(root+"/blocked_alias.tsv", badAudit);
		boolean rejected=false;
		try{invoke(root, root+"/blocked_alias.tsv", root+"/bad_alias_fit");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An excluded alias must block its whole resolved source group before learning");
		resourceFailures(root, cm, names, members, models);
		contextFit(root, cm);
	}

	/** Shuffled references bind names, including the unrefined member. */
	private static void contextFit(final String root, final CovarianceModel cm){
		write(root+"/context_refs.fa", ">raw\nGG\n>literal n=102\nAC\n");
		invoke(root, root+"/audit.tsv", root+"/context_fit", "refs="+root+"/context_refs.fa", "anchorpadding=100", "margin=2");
		final NcrnaCmCorridorResource.Data data=NcrnaCmCorridorResource.read(root+"/context_fit/RF00174.corridors.tsv.gz",
			new String[]{"literal"}, new String[]{"RF00174"}, new CovarianceModel[]{cm},
			new byte[][]{{'A', 'C'}});
		check(data.anchorPadding==100 && data.unclippedAnchors && data.retainedCalls[0]==100 && data.resolvedGroups[0]==10,
			"Context export preserves the selected hundred whole paths/ten groups and explicit100-base anchor contract");
		for(int state=0; state<3; state++){
			for(int bound=0; bound<4; bound++){
				check(data.profiles[0].offset(state, bound)==((bound&1)==0 ? -2 : 2),
					"Selected margin2 expands retained zero-residual bounds once, without changing trajectory support");
			}
		}
		check(!new File(root+"/context_fit/RAW.corridors.tsv.gz").exists(),
			"A bound unrefined sequence does not acquire a CM or state profile");
		contextOrder(root, cm);
		write(root+"/unknown_context_refs.fa", ">raw\nGG\n>unknown\nAC\n");
		write(root+"/missing_context_refs.fa", ">literal\nAC\n");
		rejectContextFit(root, "unknown_context", "refs="+root+"/unknown_context_refs.fa", "anchorpadding=100");
		rejectContextFit(root, "missing_context", "refs="+root+"/missing_context_refs.fa", "anchorpadding=100");
		rejectContextFit(root, "duplicate_context", "refs="+root+"/context_refs.fa,"+root+"/context_refs.fa", "anchorpadding=100");
		rejectContextFit(root, "undeclared_context", "refs="+root+"/context_refs.fa");
		rejectContextFit(root, "unbound_context", "anchorpadding=100");
		rejectContextFit(root, "wrong_context_geometry", "refs="+root+"/context_refs.fa", "anchorpadding=0");
		rejectContextFit(root, "negative_margin", "margin=-1");
		System.err.println("NCRNA_CM_CORRIDOR_CONTEXT_FIT_FIXTURE_PASS");
	}

	/** Same-member descriptors reverse the census order, with distinct sequences
	 * and supported/sparse states so order cannot be inferred from dimensions. */
	private static void contextOrder(final String root, final CovarianceModel cm){
		write(root+"/context_order_models.tsv", MODEL_HEADER
			+"\n0\tRF00174\tliteral\t102\t10\ttrue\n1\tRAW\traw\t1\t1\tfalse\n2\tRF00174\tsparse\t0\t0\tfalse\n");
		write(root+"/context_order_refs.fa", ">raw\nGG\n>sparse\nCC\n>literal n=102\nAC\n");
		invoke(root, root+"/audit.tsv", root+"/context_order_fit", "models="+root+"/context_order_models.tsv", "expectedmodels=3",
			"refs="+root+"/context_order_refs.fa", "anchorpadding=100", "margin=2");
		// The actual reader checks each descriptor's ordered name AND bases.
		// A writer falling back to census order fails before these support checks.
		final NcrnaCmCorridorResource.Data data=NcrnaCmCorridorResource.read(root+"/context_order_fit/RF00174.corridors.tsv.gz",
			new String[]{"sparse", "literal"}, new String[]{"RF00174", "RF00174"}, new CovarianceModel[]{cm, cm},
			new byte[][]{{'C', 'C'}, {'A', 'C'}});
		check(data.profiles.length==2 && data.anchorPadding==100 && data.unclippedAnchors
			&& data.profiles[0]==null && data.retainedCalls[0]==0 && data.resolvedGroups[0]==0
			&& data.profiles[1]!=null && data.retainedCalls[1]==100 && data.resolvedGroups[1]==10,
			"Within-member FASTA order must bind sparse CC first and supported AC second, preserving their separate support without transferring a profile");
		System.err.println("NCRNA_CM_CORRIDOR_SAME_MEMBER_ORDER_FIXTURE_PASS");
	}

	private static void rejectContextFit(final String root, final String name, final String... extra){
		final String out=root+'/'+name;
		boolean rejected=false;
		try{invoke(root, root+"/audit.tsv", out, extra);}catch(final IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(out).exists(), "A malformed context contract must fail before exporting resources: "+name);
	}

	private static void checkSummary(final String path){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			input.reader.nextLine();
			p.set(input.reader.nextLine());
			check(p.terms()==11 && p.parseInt(0)==103 && p.parseInt(1)==102 && p.parseInt(2)==102
				&& p.parseInt(3)==101 && p.parseInt(4)==1 && p.parseInt(5)==306 && p.parseInt(6)==1
				&& p.parseInt(7)==102 && p.parseInt(8)==1 && p.parseLong(9)>0 && input.reader.nextLine()==null,
				"Literal all-call/search/anchor/visit/census/fitted denominators must stay separate and conserve the declared population");
		}
	}

	private static void checkContainment(final String path){
		final LineParser1 p=new LineParser1('\t');
		int checks=0;
		try(Input input=new Input(path)){
			input.reader.nextLine();
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				if(p.parseInt(0)!=0){continue;}
				if(p.termEquals("FULL_LEARNED_RETAINED", 3)){
					check(p.parseInt(4)==100 && p.parseInt(5)==100, "All hundred retained full-fit trajectories must be contained");
					checks++;
				}else if(p.termEquals("FULL_LEARNED_TRIMMED", 3)){
					check(p.parseInt(4)==1 && p.parseInt(5)==0, "The deliberately removed shifted whole path lies outside zero-width learned residuals");
					checks++;
				}else if(p.termEquals("CV_UNIFORM_SPARSE", 3)){
					check(p.parseInt(4)==101 && p.parseInt(5)==100, "Every held-out path uses the sparse uniform0control; only the shifted path misses it");
					checks++;
				}
			}
		}
		check(checks==3, "The literal retained, trimmed and sparse-CV strata must all be reported");
	}

	private static void resourceFailures(final String root, final CovarianceModel cm, final String[] names,
			final String[] members, final CovarianceModel[] models){
		final CovarianceModelStateCorridors signed=new CovarianceModelStateCorridors(cm,
			new short[]{1, -2, -1, 2, 0, 0, 0, 0, 0, 0, 0, 0}, new boolean[]{true, false, false});
		NcrnaCmCorridorResource.write(root+"/signed.tsv.gz", names, members, models,
			new CovarianceModelStateCorridors[]{signed, null}, new int[]{100, 0}, new int[]{10, 0});
		final NcrnaCmCorridorResource.Data loaded=NcrnaCmCorridorResource.read(root+"/signed.tsv.gz", names, members, models);
		check(loaded.profiles[0].offset(0, 0)==1 && loaded.profiles[0].offset(0, 1)==-2
			&& !loaded.profiles[0].observed(1), "The versioned reader must preserve signed unordered residuals and unseen-state fallback");
		final String text=text(root+"/signed.tsv.gz");
		reject(root+"/overflow.tsv", text.replace("\nS\t0\t1", "\nS\t0\t40000"), names, members, models);
		reject(root+"/duplicate.tsv", text.replace("\nM\t1", "\nS\t0\t1\t-2\t-1\t2\nM\t1"), names, members, models);
		reject(root+"/truncated.tsv", text.substring(0, text.indexOf("\nM\t1")+1), names, members, models);
		reject(root+"/wrong_model.tsv", text.replace("\tliteral\tRF00174", "\tother\tRF00174"), names, members, models);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private static void invoke(final String root, final String audit, final String out, final String... extra){
		final String[] baseline={"in="+root+"/parts.tsv", "calls="+root+"/original.tsv",
			"audit="+audit, "models="+root+"/models.tsv", "cms="+root+"/literal.cm", "out="+out,
			"expectedcalls=103", "expectedsearches=102", "expectedsteps=306", "expectedmodels=2", "expectedassemblies=21",
			"mincalls=100", "mingroups=10", "folds=5", "uniform=0", "t=1", "ow=f", "unpigz=f", "pigz=f"};
		final String[] args=java.util.Arrays.copyOf(baseline, baseline.length+extra.length);
		System.arraycopy(extra, 0, args, baseline.length, extra.length);
		BuildNcrnaCmCorridors.main(args);
	}
	private static void reject(final String path, final String value, final String[] names, final String[] members, final CovarianceModel[] models){
		write(path, value);
		boolean rejected=false;
		try{NcrnaCmCorridorResource.read(path, names, members, models);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Malformed resource fixture must fail: "+path);
	}
	private static String text(final String path){
		final ByteBuilder rows=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] line; (line=input.reader.nextLine())!=null;){rows.append(line).nl();}
		}
		return new String(rows.toBytes(), StandardCharsets.UTF_8);
	}
	private static void write(final String path, final String text){
		assert(path!=null && text!=null) : "Every literal fixture file has explicit owned content and destination";
		try(Output output=new Output(path)){output.writer.print(text);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
	private static final String AUDIT_HEADER="assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled";
	private static final String MODEL_HEADER="model_ordinal0\tmember\tmodel\taccepted_calls\tresolved_assembly_groups\tlearned_corridor_supported";
}
