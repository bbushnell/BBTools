package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import parse.LineParser1;
import prok.BuildNcrnaReference;
import prok.GradeNcrnaReference;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Literal subtype outcomes and union parity; no caller or CM searches.
 * Each case has one input plus an explicitly empty second input. Evidence is
 * synthetic, with fixed roles and tiers rather than derived scientific scores.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaSubtypeGradeTest{

	NcrnaSubtypeGradeTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Subtype fixture root must be fresh");
		final Fixture single=new Fixture("single");
		final Row a=single.parent(100, 199, false, 0);
		single.hit(a, 0, "CONFIDENT", GROUP, -1);
		single.call(0, 100, 199);
		single.run();
		single.counts(0, 0, 1, 0, 0, 0, 0, 0);
		single.counts(1, 0, 1, 0, 0, 0, 0, 0);

		final Fixture multi=new Fixture("multi");
		final Row ab=multi.parent(100, 199, false, 0);
		multi.hit(ab, 0, "CONFIDENT", GROUP, -1);
		multi.hit(ab, 1, "CONFIDENT", GROUP, -1);
		multi.call(0, 100, 199);
		multi.call(2, 100, 199);
		multi.run();
		multi.counts(0, 0, 1, 0, 0, 0, 0, 0);
		multi.counts(0, 1, 0, 0, 1, 0, 0, 0);
		multi.counts(1, 0, 0, 0, 0, 1, 1, 0);
		multi.counts(1, 1, 0, 0, 0, 1, 0, 0);
		for(int view=0; view<2; view++){multi.counts(view, 2, 0, 1, 0, 0, 0, 0);}
		check(lines(multi.path("dual.subtype_counts.tsv"))[1].startsWith("2\t1\t0\t1\t0\t0\t"),
			"Two model supports count as one ambiguous parent, including the empty input in population");

		final Fixture both=new Fixture("both");
		final Row bothParent=both.parent(100, 199, false, 0);
		both.hit(bothParent, 0, "CONFIDENT", GROUP, -1);
		both.hit(bothParent, 1, "CONFIDENT", GROUP, -1);
		both.call(0, 100, 199);
		both.call(1, 100, 199);
		both.run();
		for(int type=0; type<2; type++){
			both.counts(0, type, 1, 0, 0, 0, 0, 0);
			both.counts(1, type, 0, 0, 0, 1, 1, 0);
		}
		check(lines(both.path("dual.summary.tsv"))[1].startsWith("cobalamin\t2\t1\t1\t0\t"),
			"Two subtype TPs cannot create a second union reference; duplicate union call is FP");

		final Fixture duplicate=new Fixture("duplicate");
		final Row da=duplicate.parent(100, 199, false, 0);
		duplicate.hit(da, 0, "CONFIDENT", GROUP, -1);
		duplicate.hit(da, 0, "CONFIDENT", GROUP, -1);
		duplicate.call(0, 100, 199);
		duplicate.call(0, 100, 199);
		duplicate.run();
		for(int view=0; view<2; view++){duplicate.counts(view, 0, 1, 1, 0, 0, 0, 1);}
		check(lines(duplicate.path("dual.subtype_counts.tsv"))[1].startsWith("2\t1\t1\t0\t"),
			"Repeated hits from one subtype cannot make a parent multi-supported");

		final Fixture mask=new Fixture("mask");
		final Row ma=mask.parent(100, 199, true, 0);
		mask.hit(ma, 0, "MARGINAL", GROUP, -1);
		mask.call(0, 100, 199);
		mask.call(0, 100, 199);
		mask.call(2, 100, 199);
		mask.run();
		for(int view=0; view<2; view++){
			mask.counts(view, 0, 0, 0, 0, 1, 2, 1);
			mask.counts(view, 2, 0, 0, 0, 1, 1, 0);
		}

		final Fixture marginal=new Fixture("marginal");
		final Row mb=marginal.parent(100, 199, false, 0);
		marginal.hit(mb, 0, "CONFIDENT", GROUP, -1);
		marginal.hit(mb, 1, "MARGINAL", GROUP, -1);
		marginal.call(0, 100, 199);
		marginal.call(1, 100, 199);
		marginal.run();
		for(int view=0; view<2; view++){
			marginal.counts(view, 0, 1, 0, 0, 0, 0, 0);
			marginal.counts(view, 1, 0, 0, 0, 1, 1, 0);
		}
		final Fixture thirdMarginal=new Fixture("third_marginal");
		final Row abc=thirdMarginal.parent(100, 199, false, 0);
		thirdMarginal.hit(abc, 0, "CONFIDENT", GROUP, -1);
		thirdMarginal.hit(abc, 1, "CONFIDENT", GROUP, -1);
		thirdMarginal.hit(abc, 2, "MARGINAL", GROUP, -1);
		thirdMarginal.call(2, 100, 199);
		thirdMarginal.run();
		for(int view=0; view<2; view++){thirdMarginal.counts(view, 2, 0, 0, 0, 1, 1, 0);}

		for(final String role:new String[]{WEAK, BRIDGE}){
			final Fixture weak=new Fixture(role.equals(WEAK) ? "weak" : "bridge");
			final Row wa=weak.parent(100, 199, false, 0);
			weak.hit(wa, 0, "CONFIDENT", GROUP, -1);
			weak.hit(wa, 1, role.equals(WEAK) ? "BELOW_MARGINAL" : "CONFIDENT", role, -1);
			weak.call(1, 100, 199);
			weak.run();
			for(int view=0; view<2; view++){
				weak.counts(view, 0, 0, 0, 1, 0, 0, 0);
				weak.counts(view, 1, 0, 1, 0, 0, 0, 0);
			}
		}
		for(int hit=0; hit<2; hit++){
			final Fixture prior=new Fixture("prior"+hit);
			final Row pa=prior.parent(100, 199, false, 0);
			prior.prior(pa, 0);
			if(hit==1){prior.hit(pa, 1, "COVERED_BY_PROTECTED_PRIOR", PRIOR, 2);}
			prior.call(1, 100, 199);
			prior.run();
			for(int view=0; view<2; view++){
				prior.counts(view, 0, 0, 0, 1, 0, 0, 0);
				prior.counts(view, 1, 0, 1, 0, 0, 0, 0);
			}
		}

		final Fixture competition=new Fixture("competition");
		final Row ambiguous=competition.parent(100, 199, false, 0);
		final Row exclusive=competition.parent(150, 249, false, 0);
		competition.hit(ambiguous, 0, "CONFIDENT", GROUP, -1);
		competition.hit(ambiguous, 1, "CONFIDENT", GROUP, -1);
		competition.hit(exclusive, 0, "CONFIDENT", GROUP, -1);
		competition.call(0, 100, 199);
		competition.run();
		competition.counts(0, 0, 1, 0, 1, 0, 0, 0);
		competition.counts(1, 0, 1, 0, 0, 1, 0, 0);
		competition.paired(0, ambiguous.id);
		competition.paired(1, exclusive.id);

		final Fixture precedence=new Fixture("precedence");
		final Row maskFirst=precedence.parent(100, 199, true, 0);
		final Row scored=precedence.parent(150, 249, false, 0);
		precedence.hit(maskFirst, 0, "MARGINAL", GROUP, -1);
		precedence.hit(scored, 0, "CONFIDENT", GROUP, -1);
		precedence.call(0, 100, 199);
		precedence.run();
		for(int view=0; view<2; view++){precedence.counts(view, 0, 1, 0, 0, 1, 0, 0);}

		final Fixture zero=new Fixture("zero");
		final Row za=zero.parent(100, 199, false, 0);
		zero.hit(za, 0, "CONFIDENT", GROUP, -1);
		zero.run();
		for(int view=0; view<2; view++){zero.counts(view, 0, 0, 0, 1, 0, 0, 0);}

		final Fixture reverse=new Fixture("reverse");
		final Row ra=reverse.parent(100, 199, false, 0, '-');
		reverse.hit(ra, 0, "CONFIDENT", GROUP, -1);
		reverse.call(0, 101, 199, '-');
		reverse.run();
		for(int view=0; view<2; view++){
			reverse.counts(view, 0, 1, 0, 0, 0, 0, 0);
			reverse.offsets(view);
		}

		final Fixture missing=new Fixture("missing_event");
		final Row missingParent=missing.parent(100, 199, false, 0);
		missing.hit(missingParent, 0, "CONFIDENT", GROUP, -1);
		missing.events.clear().append(EVENT_HEADER).nl();
		missing.reject("binding mismatch");
		final Fixture unsupported=new Fixture("no_scored_support");
		final Row ua=unsupported.parent(100, 199, false, 0);
		unsupported.hit(ua, 0, "MARGINAL", GROUP, -1);
		unsupported.reject("lacks justified scored subtype support");
		final Fixture badBinding=new Fixture("bad_binding");
		final Row ba=badBinding.parent(100, 199, false, 0);
		badBinding.hit(ba, 0, "CONFIDENT", GROUP, -1);
		badBinding.wrongInput=true;
		badBinding.reject("binding mismatch");
		System.err.println("NCRNA_SUBTYPE_GRADE_TEST_PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Synthetic Inputs       ----------------*/
	/*--------------------------------------------------------------*/

	private final class Fixture{
		Fixture(final String name){
			directory=new File(root, name);
			check(directory.mkdir(), "Each subtype fixture must be fresh");
			events.append(EVENT_HEADER).nl();
			witnesses.append(WITNESS_HEADER).nl();
		}
		Row parent(final int start, final int stop, final boolean mask, final int representative){
			return parent(start, stop, mask, representative, '+');
		}
		Row parent(final int start, final int stop, final boolean mask, final int representative, final char strand){
			final Row row=codec.create("I", new Locus("seq", "cobalamin", start, stop, strand), mask,
				mask ? "partial_model_support" : null, MODELS[representative], path("witnesses.tsv"));
			parents.add(row);
			return row;
		}
		void hit(final Row row, final int type, final String tier, final String association, final int priorLine){
			final int id=++nextEvent;
			events.append(id).append("\tI\tseq\t").append(row.locus.start).tab().append(row.locus.stop)
				.tab().append(row.locus.strand).tab().append(MODELS[type]).append("\tcobalamin\t1\t100\t100\tno\tfalse\tfalse\t50.0\t40.0\t125.000000\ttrue\tfalse\t80\ttrue\tfalse\t")
				.append(tier).append("\tsynthetic\t.\t").append(association.equals(BRIDGE) ? 16 : 0)
				.append("\t/synthetic/native.tbl\t").append(id+1).append("\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n");
			BuildNcrnaReference.field(witnesses, row.id).append("\tI\tCM_HIT\t").append(id)
				.append("\t/synthetic/native.tbl\t").append(id+1).tab().append(MODELS[type]).tab()
				.append(row.locus.start).tab().append(row.locus.stop).tab().append(priorLine).tab().append(association).nl();
		}
		void prior(final Row row, final int type){
			BuildNcrnaReference.field(witnesses, row.id).append("\tI\tPROTECTED_PRIOR\t.\t.\t.\t")
				.append(MODELS[type]).tab().append(row.locus.start).tab().append(row.locus.stop)
				.append("\t2\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n");
		}
		void call(final int type, final int start, final int stop){
			call(type, start, stop, '+');
		}
		void call(final int type, final int start, final int stop, final char strand){
			calls.append("seq\tBBTools\tRNA\t").append(start).tab().append(stop)
				.append("\t.\t").append(strand).append("\t.\tmodel:cb").append(NAMES[type]).append("_c1,family:").append(NAMES[type]).nl();
		}
		void prepare(){
			codec.write(path("reference.gff"), "I", parents);
			codec.write(path("zero.ref.gff"), "ZERO", new ArrayList<Row>());
			write("events.tsv", events.toString());
			write("witnesses.tsv", wrongInput ? witnesses.toString().replace("\tI\tCM_HIT", "\tZERO\tCM_HIT") : witnesses.toString());
			write("calls.gff", calls.toString());
			write("zero.calls.gff", "");
			write("reference.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\n"+binding("I", "reference.gff")+binding("ZERO", "zero.ref.gff"));
			write("queries.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\n"+binding("I", "calls.gff")+binding("ZERO", "zero.calls.gff"));
			write("families.tsv", GradeNcrnaReference.RULE_HEADER+"\ncobalamin\tcbA\tA\ttrue\ncobalamin\tcbB\tB\ttrue\ncobalamin\tcbC\tC\ttrue\n");
			write("subtypes.tsv", "subtype\tmodel_prefix\treference_model\nA\tcbA\tRF00174\nB\tcbB\tRF01482\nC\tcbC\tRF01689\n");
		}
		void run(){
			prepare();
			invoke(false);
			invoke(true);
			for(final String suffix:new String[]{".raw.tsv", ".inputs.tsv", ".summary.tsv"}){
				check(Arrays.equals(lines(path("union"+suffix)), lines(path("dual"+suffix))), "Subtype option must leave every ordinary union row unchanged: "+directory+suffix);
			}
			check(lines(path("dual.subtype_summary.tsv")).length==7, "Both views must contain all three configured subtypes, including zero counts");
			final LineParser1 p=new LineParser1('\t');
			for(final String line:lines(path("dual.subtype_raw.tsv"))){
				p.set(line.getBytes(StandardCharsets.UTF_8));
				check(p.terms()==22, "Every subtype raw row must have twenty-two columns");
			}
			for(final String line:lines(path("dual.subtype_support.tsv"))){
				p.set(line.getBytes(StandardCharsets.UTF_8));
				check(p.terms()==12, "Each parent support row retains twelve fields");
			}
		}
		void counts(final int view, final int type, final int tp, final int fp, final int fn,
			final int maskRefs, final int masked, final int duplicates){
			final LineParser1 p=new LineParser1('\t');
			int found=0;
			for(final String line:lines(path("dual.subtype_summary.tsv"))){
				p.set(line.getBytes(StandardCharsets.UTF_8));
				if(!p.termEquals(VIEWS[view], 0) || !p.termEquals(NAMES[type], 1)){continue;}
				found++;
				check(p.terms()==17 && p.parseInt(2)==2 && p.parseLong(3)==tp && p.parseLong(4)==fp
					&& p.parseLong(5)==fn && p.parseLong(6)==tp+fn && p.parseLong(7)==maskRefs
					&& p.parseLong(8)==masked && p.parseLong(9)==duplicates, "Literal subtype counts disagree: "+directory+" / "+line);
				check(line.contains(view==0 ? "model-support recovery, not a missed gene, not wrong-subtype proof"
					: "not unique-biological-subtype proof"), "Scientific interpretation must be attached to each summary row");
			}
			check(found==1, "Subtype summary row must exist exactly once");
		}
		void paired(final int view, final String parent){
			final LineParser1 p=new LineParser1('\t');
			int found=0;
			for(final String line:lines(path("dual.subtype_raw.tsv"))){
				p.set(line.getBytes(StandardCharsets.UTF_8));
				if(p.termEquals(VIEWS[view], 0) && p.termEquals("A", 1) && p.termEquals("QUERY", 3)){
					check(p.termEquals("TP", 11) && p.termEquals(parent, 12), "Each view must match independently, allowing reassignment after ambiguity masking");
					found++;
				}
			}
			check(found==1, "Competition fixture has exactly one A query");
		}
		void offsets(final int view){
			final LineParser1 p=new LineParser1('\t');
			int found=0;
			for(final String line:lines(path("dual.subtype_raw.tsv"))){
				p.set(line.getBytes(StandardCharsets.UTF_8));
				if(p.termEquals(VIEWS[view], 0) && p.termEquals("A", 1) && p.termEquals("QUERY", 3)){
					check(p.parseInt(20)==0 && p.parseInt(21)==-1, "Minus-strand subtype offsets must use transcript direction");
					found++;
				}
			}
			check(found==1, "Reverse fixture has exactly one A query per view");
		}
		void reject(final String expected){
			prepare();
			try{invoke(true);}
			catch(IllegalArgumentException failure){
				check(failure.getMessage().contains(expected), "Unexpected subtype rejection: "+failure);
				check(!new File(path("dual.summary.tsv")).exists(), "Failed subtype validation must not publish overall completion");
				return;
			}
			throw new AssertionError("Expected subtype rejection: "+expected);
		}
		void invoke(final boolean dual){
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+path("queries.tsv"), "reference="+path("reference.tsv"),
				"familymap="+path("families.tsv"), "out="+path(dual ? "dual" : "union"), "t=1", "bgzip=f", "zipthreads=1"));
			if(dual){
				args.add("subtypes="+path("subtypes.tsv"));
				args.add("witnesses="+path("witnesses.tsv"));
				args.add("events="+path("events.tsv"));
			}
			GradeNcrnaReference.main(args.toArray(new String[args.size()]));
		}
		String binding(final String id, final String file){return id+"\t/synthetic/"+id+".fa\t"+path(file)+"\n";}
		String path(final String file){return new File(directory, file).getAbsolutePath();}
		void write(final String file, final String text){try(Output out=new Output(path(file))){out.writer.print(text);}}
		final File directory;
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> parents=new ArrayList<Row>();
		final ByteBuilder events=new ByteBuilder(), witnesses=new ByteBuilder(), calls=new ByteBuilder();
		int nextEvent;
		boolean wrongInput;
	}

	private static String[] lines(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){result.add(new String(line, StandardCharsets.UTF_8));}
		}
		return result.toArray(new String[result.size()]);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private static final String[] MODELS={"RF00174", "RF01482", "RF01689"}, NAMES={"A", "B", "C"};
	private static final String[] VIEWS={"MODEL_SUPPORT_RECOVERY", "EXCLUSIVE_SUBTYPE_ACCURACY"};
	private static final String GROUP="QUALIFYING_SAME_FAMILY_GROUP_MEMBER", PRIOR="SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS";
	private static final String WEAK="WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT", BRIDGE="TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE";
	private static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation";
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
}
