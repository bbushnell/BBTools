package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.GradeNcrnaReference;
import prok.NcrnaDdlTruth;
import prok.NcrnaDdlTruth.Group;
import prok.NcrnaDdlTruth.Parent;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal reference qualification for the DDL driver, without caller/resources.
 * One selected input has six physical parents; an empty input and an unselected
 * nonexistent reference disprove denominator loss and accidental broad reads.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlTruthTest{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlTruthTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown truth fixture parameter: "+arg);
			}
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use -ea and out=fresh_directory; this fixture owns complete synthetic metadata");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Truth fixture directory must be fresh");
			new NcrnaDdlTruthTest(root).run();
			pp.outstream.println("NCRNA_DDL_TRUTH_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaDdlTruthTest(final File root_){root=root_;}

	/*--------------------------------------------------------------*/
	/*----------------       Literal Coverage       ----------------*/
	/*--------------------------------------------------------------*/

	private void run(){
		write("families.tsv", GradeNcrnaReference.RULE_HEADER+"\n"
			+"cobalamin\tcbA\tA\ttrue\ncobalamin\tcbB\tB\ttrue\n");
		write("subtypes.tsv", "subtype\tmodel_prefix\treference_model\nA\tcbA\tRF00174\nB\tcbB\tRF01482\n");
		events.append(EVENT_HEADER).nl();
		witnesses.append(WITNESS_HEADER).nl();
		final Row main=parent(100, 199, '+', false, 0);
		hit(main, 0, "CONFIDENT", GROUP, -1);
		hit(main, 0, "CONFIDENT", GROUP, -1);
		final Row multi=parent(210, 309, '+', false, 0);
		hit(multi, 0, "CONFIDENT", GROUP, -1);
		hit(multi, 1, "CONFIDENT", GROUP, -1);
		final Row prior=parent(320, 419, '+', false, 0);
		BuildNcrnaReference.field(witnesses, prior.id).append("\tI\tPROTECTED_PRIOR\t.\t.\t.\tRF00174\t320\t419\t2\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n");
		hit(prior, 1, "COVERED_BY_PROTECTED_PRIOR", PRIOR, 2);
		final Row mask=parent(430, 529, '+', true, 0);
		hit(mask, 0, "MARGINAL", GROUP, -1);
		final Row other=parent(540, 639, '+', false, 1);
		hit(other, 1, "CONFIDENT", GROUP, -1);
		final Row minus=parent(650, 749, '-', false, 0);
		hit(minus, 0, "CONFIDENT", GROUP, -1);
		codec.write(path("reference.gff"), "I", parents);
		codec.write(path("empty.gff"), "ZERO", new ArrayList<Row>());
		write("events.tsv", events.toString());
		write("witnesses.tsv", witnesses.toString());
		write("reference.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\n"
			+"I\t"+path("input.fa")+'\t'+path("reference.gff")+'\n'
			+"ZERO\t"+path("zero.fa")+'\t'+path("empty.gff")+'\n'
			+"UNSELECTED\t"+path("absent.fa")+'\t'+path("must_not_be_opened.gff")+'\n');
		final NcrnaDdlTruth.Input[] inputs=load(new String[]{path("zero.fa"), path("input.fa")});
		check(inputs.length==2 && inputs[0].binding.id.equals("ZERO") && inputs[0].parents.length==0
			&& inputs[0].targets==0 && inputs[0].group("chr")==null,
			"Caller order and explicitly empty input must survive without opening an unselected GFF");
		final NcrnaDdlTruth.Input input=inputs[1];
		check(input.binding.id.equals("I") && input.parents.length==6 && input.targets==4 && input.masked==1,
			"Six distinct physical parents include exactly four RF00174 scored targets and one mask");
		final Group group=input.group("chr");
		check(group!=null && group.parents.length==6 && group.targets.length==4 && group.scored.length==5,
			"Other-subtype scored truth joins attribution; masks and repeated subtype witnesses do not add parents");
		check(input.group("chr description")==null, "The truth index uses the literal native GFF seqid, never an invented full title");
		final Parent[] bound=input.parents;
		check(bound[0].target && !bound[0].multi && bound[0].scoredModels.equals("RF00174"),
			"Two witnesses for the same subtype still describe one single-supported physical target");
		check(bound[1].target && bound[1].multi && bound[1].scoredModels.equals("RF00174,RF01482"),
			"Multiple scored subtype support remains target-eligible and explicitly ambiguous");
		check(bound[2].target && bound[2].prior && !bound[2].multi && bound[2].scoredModels.equals("RF00174"),
			"An overlapping current-model observation cannot replace the protected prior's original subtype authority");
		check(!bound[3].target && bound[3].scoredIndex==-1 && bound[3].scoredModels.equals(".")
			&& bound[3].marginalModels.equals("RF00174"), "A marginal mask remains unscored and cannot consume unique attribution");
		check(!bound[4].target && bound[4].scoredIndex==3 && bound[4].scoredModels.equals("RF01482")
			&& bound[5].target && bound[5].row.locus.strand=='-' && bound[5].scoredIndex==4,
			"Other-subtype and minus-strand parents retain their original physical attribution indices");
		for(final Parent p:bound){
			check(p.scoredIndex<0 || group.scored[p.scoredIndex]==p.row.locus,
				"Every scored attribution index must bind its original physical locus object");
		}
		reject(new String[]{path("input.fa"), path("input.fa")});
		reject(new String[]{path("not_in_manifest.fa")});
	}

	private NcrnaDdlTruth.Input[] load(final String[] selected){
		return NcrnaDdlTruth.load(path("reference.tsv"), selected, path("subtypes.tsv"),
			path("witnesses.tsv"), path("events.tsv"), GradeNcrnaReference.loadRules(path("families.tsv")));
	}
	private void reject(final String[] selected){
		boolean rejected=false;
		try{load(selected);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Duplicate or unbound selected inputs must not alter the original truth denominator");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixture Encoding       ----------------*/
	/*--------------------------------------------------------------*/

	private Row parent(final int start, final int stop, final char strand, final boolean mask, final int type){
		final Row row=codec.create("I", new Locus("chr", "cobalamin", start, stop, strand), mask,
			mask ? "partial_model_support" : null, MODELS[type], path("witnesses.tsv"));
		parents.add(row);
		return row;
	}
	private void hit(final Row row, final int type, final String tier, final String association, final int priorLine){
		final int event=++nextEvent;
		events.append(event).append("\tI\tchr\t").append(row.locus.start).tab().append(row.locus.stop)
			.tab().append(row.locus.strand).tab().append(MODELS[type])
			.append("\tcobalamin\t1\t100\t100\t0\tfalse\tfalse\t50.0\t40.0\t125.000000\ttrue\tfalse\t90\ttrue\tfalse\t")
			.append(tier).append("\tsynthetic\t.\t0\t/synthetic/native.tbl\t").append(event+1)
			.append("\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n");
		BuildNcrnaReference.field(witnesses, row.id).append("\tI\tCM_HIT\t").append(event)
			.append("\t/synthetic/native.tbl\t").append(event+1).tab().append(MODELS[type]).tab()
			.append(row.locus.start).tab().append(row.locus.stop).tab().append(priorLine).tab().append(association).nl();
	}
	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final File root;
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private final ArrayList<Row> parents=new ArrayList<Row>();
	private final ByteBuilder events=new ByteBuilder(), witnesses=new ByteBuilder();
	private int nextEvent;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] MODELS={"RF00174", "RF01482"};
	private static final String GROUP="QUALIFYING_SAME_FAMILY_GROUP_MEMBER";
	private static final String PRIOR="SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS";
	private static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation";
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
}
