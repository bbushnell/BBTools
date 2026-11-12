package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;

/** Describes an existing tiered reference without changing its decisions.
 * Raw events, physical annotations and protected priors have distinct counters.
 * Primary marginal reasons partition their own unit; cause incidences overlap.
 * Physical masks use the serialized reason, never a witness-derived replacement.
 * One thread streams events and retains only small counters, bounded examples,
 * protected-prior identities and one input's decoded GFF rows at a time.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceTierReport{

	static void process(final String input, final String output, final String sampleModel, final int examples){
		final NcrnaReferenceTierReport report=new NcrnaReferenceTierReport(input, output, sampleModel, examples);
		report.loadExpected();
		report.readEvents();
		report.readPriors();
		report.readAnnotations();
		report.writeOverlaps();
		report.writeCounts();
	}
	private NcrnaReferenceTierReport(final String input_, final String output_, final String sampleModel_, final int examples_){
		input=input_;
		output=output_;
		sampleModel=sampleModel_;
		examples=examples_;
		require(examples>=0 && examples<=10, "examplespergroup must be within0..10 to bound retained examples");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved populations      ----------------*/
	/*--------------------------------------------------------------*/

	private void loadExpected(){
		boolean schema=false;
		try(Input in=new Input(path("parameters.tsv"))){
			RrnaResourceIO.header(in.reader, "parameter\tvalue", input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==2, "Tier parameter rows require two fields");
				if(parser.termEquals("schema", 0)){
					require(!schema && parser.termEquals("NCRNA_REFERENCE_TIERS_V1", 1), "Unknown or repeated reference tier schema");
					schema=true;
				}
			}
		}
		require(schema, "A completed tier parameter schema is required");
		try(Input in=new Input(path("family_counts.tsv"))){
			RrnaResourceIO.header(in.reader, FAMILY_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==12, "Family accounting requires all twelve fields");
				probe.set(parser, 0);
				require(!familyKeys.contains(probe), "Repeated reference family accounting");
				final Family family=new Family(text(0));
				familyKeys.put(probe.copy(), family);
				familyNames.put(family.name, family);
				families.add(family);
				readExpected(family.events, 2);
				family.expectedScored=nonnegative(9);
				family.expectedMasked=nonnegative(10);
				family.expectedPriors=nonnegative(11);
			}
		}
		boolean total=false;
		try(Input in=new Input(path("event_counts.tsv"))){
			RrnaResourceIO.header(in.reader, MODEL_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==9, "Model accounting requires all nine fields");
				if(parser.termEquals("ALL_QUERIED", 0)){
					require(!total, "Repeated global event accounting");
					total=true;
					readExpected(all, 2);
					continue;
				}
				probe.set(parser, 0);
				require(!modelKeys.contains(probe), "Repeated model accounting");
				final Model model=new Model(text(0));
				modelKeys.put(probe.copy(), model);
				models.add(model);
				readExpected(model.events, 2);
			}
		}
		require(total && !families.isEmpty() && !models.isEmpty(), "Family/model/global population records are required");
		if(sampleModel!=null){
			final byte[] key=sampleModel.getBytes(StandardCharsets.UTF_8);
			probe.set(key, 0, key.length);
			require(modelKeys.contains(probe), "Requested samplemodel is not in the saved model population");
		}
	}
	private void readExpected(final Counts counts, final int first){
		long sum=0;
		for(int i=0; i<STATES.length; i++){
			counts.expected[i]=nonnegative(first+i);
			sum=Math.addExact(sum, counts.expected[i]);
		}
		require(sum==nonnegative(1), "Saved tier partition must conserve its raw-event total");
	}

	/** Reads each event exactly once for counters; examples retain original rows. */
	private void readEvents(){
		try(Input in=new Input(path("event_provenance.tsv.gz")); Output out=new Output(output+".examples.tsv")){
			RrnaResourceIO.header(in.reader, EVENT_HEADER, input);
			out.writer.println(EVENT_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==30, "Tier event provenance requires thirty fields");
				probe.set(parser, 6);
				final Model model=modelKeys.get(probe);
				probe.set(parser, 7);
				final Family family=familyKeys.get(probe);
				require(model!=null && family!=null && (model.family==null || model.family==family), "Event model/family differs from saved accounting");
				model.family=family;
				final int tier=state();
				final int decision=parser.parseInt(25), source=NcrnaReferenceStrata.readFlags(parser, 24);
				require(decision>=0 && (decision&~29)==0, "Unknown tier decision flag cannot silently disappear from reporting");
				final boolean excluded=!bool(17);
				final boolean partial=parser.parseInt(8)!=1 || parser.parseInt(9)!=parser.parseInt(10) || parser.parseInt(11)!=0;
				final boolean cross=(decision&1)!=0, warning=(source&SOURCE_BLOCKERS)!=0, bridge=(decision&16)!=0;
				final int reason=tier==1 ? reason(23) : -1;
				count(all, tier, reason, excluded, partial, cross, warning, bridge);
				count(family.events, tier, reason, excluded, partial, cross, warning, bridge);
				count(model.events, tier, reason, excluded, partial, cross, warning, bridge);
				if(sampleModel!=null && parser.termEquals(sampleModel, 6)){
					final int bucket=(reason<0 ? 0 : reason)+(cross ? 4 : 0);
					if(sampleCounts[tier][bucket]++<examples){
						final Sample sample=new Sample(parser);
						probe.set(parser, 1);
						IntList list=samplesByInput.get(probe);
						if(list==null){
							list=new IntList();
							samplesByInput.put(probe.copy(), list);
						}
						list.add(sampleRows.size());
						sampleRows.add(sample);
						samples++;
						out.writer.println(row);
					}
				}
			}
		}
		checkCounts(all);
		for(final Family family:families){checkCounts(family.events);}
		for(final Model model:models){checkCounts(model.events);}
	}
	private static void count(final Counts counts, final int tier, final int reason, final boolean excluded,
		final boolean partial, final boolean cross, final boolean warning, final boolean bridge){
		assert(tier>=0 && tier<STATES.length && (tier!=1 || reason>=0 && reason<4)) : "Only known saved tiers and marginal reasons can enter the report partition";
		counts.actual[tier]++;
		if(reason>=0){counts.reasons[reason]++;}
		final long[] incidence=counts.incidence[tier];
		if(excluded){incidence[0]++;}
		if(partial){incidence[1]++;}
		if(cross){incidence[2]++;}
		if(warning){incidence[3]++;}
		if(bridge){incidence[4]++;}
	}
	private static void checkCounts(final Counts counts){
		for(int i=0; i<STATES.length; i++){
			require(counts.expected[i]==counts.actual[i], "Observed event tier differs from its saved family/model accounting");
		}
		long marginal=0;
		for(final long n:counts.reasons){marginal+=n;}
		require(marginal==counts.actual[1], "Exclusive primary reasons must partition raw marginal events");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Physical annotations    ----------------*/
	/*--------------------------------------------------------------*/

	private void readPriors(){
		try(Input in=new Input(path("annotation_witnesses.tsv.gz"))){
			RrnaResourceIO.header(in.reader, WITNESS_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==11, "Annotation witnesses require eleven fields");
				if(parser.termEquals("CM_HIT", 2)){continue;}
				require(parser.termEquals("PROTECTED_PRIOR", 2), "Unknown annotation witness kind");
				final String id=decoded(0);
				require(!priorIds.contains(id), "A protected physical annotation needs exactly one prior authority witness");
				priorIds.put(id, new Prior(decoded(1), decoded(6), parser.parseInt(7), parser.parseInt(8)));
			}
		}
	}
	/** Reuses the canonical codec, retaining decoded rows for only one input. */
	private void readAnnotations(){
		final ObjectSet<String> inputs=new ObjectSet<String>(String.class), files=new ObjectSet<String>(String.class);
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		long priorsSeen=0;
		try(Input in=new Input(path("reference.tsv"))){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.MANIFEST_HEADER, input);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==3, "Reference manifest requires three literal fields");
				final String id=text(0), gff=text(2);
				require(inputs.add(id) && files.add(gff) && new File(gff).isAbsolute(), "Reference input IDs/GFFs must be unique with absolute paths");
				for(final Row row:codec.read(gff, id)){
					final Family family=familyNames.get(row.locus.family);
					require(family!=null && row.scores!=null, "Every tier annotation needs a declared family and complete score metadata");
					if(row.mask){
						final int reason=reason(row.reason);
						require(reason<3, "A bridge witness cannot be an independent physical mask");
						family.masks[reason]++;
						family.masked++;
					}else{family.scored++;}
					final Prior prior=priorIds.get(row.id);
					if(prior!=null){
						require(!prior.seen && !row.mask && id.equals(prior.input) && prior.model.equals(row.model)
							&& row.locus.start==prior.start && row.locus.stop==prior.stop, "Protected prior must retain one scored physical identity and its original interval/model");
						prior.seen=true;
						priorsSeen++;
						if(row.scores.missing){family.missingPrior++;}
						else{family.measuredPrior++;}
					}else{require(!row.scores.missing, "Missing native score is permitted only for a named protected prior");}
				}
			}
		}
		inputCount=inputs.size();
		require(inputCount>0 && priorsSeen==priorIds.size(), "Every prior authority witness must resolve exactly once in the primary annotation files");
		for(final Family family:families){
			require(family.scored==family.expectedScored && family.masked==family.expectedMasked
				&& family.missingPrior+family.measuredPrior==family.expectedPriors, "Physical annotation/prior counts disagree with saved family totals");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Bounded witnesses      ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent coordinate lookup for sampled other-family evidence. This does
	 * not reinterpret dispositions: only the constructor's reported qualification
	 * fields select witnesses, and both original rows retain their exact values. */
	private void writeOverlaps(){
		try(Output out=new Output(output+".example_overlaps.tsv")){
			out.writer.println("sample_event_id\t"+EVENT_HEADER);
			if(samples==0){return;}
			try(Input in=new Input(path("event_provenance.tsv.gz"))){
				RrnaResourceIO.header(in.reader, EVENT_HEADER, input);
				for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
					parser.set(row);
					probe.set(parser, 1);
					final IntList selected=samplesByInput.get(probe);
					if(selected==null || parser.termEquals("INVALID_INTERVAL", 22) || (parser.parseInt(25)&8)!=0){continue;}
					if(!bool(17) && !(bool(20) && !bool(21))){continue;}
					probe.set(parser, 2);
					familyProbe.set(parser, 7);
					for(int i=0; i<selected.size; i++){
						final Sample sample=sampleRows.get(selected.get(i));
						if(!sample.valid || sample.witnesses>=examples || !sample.target.equals(probe) || sample.family.equals(familyProbe)
							|| !parser.termEquals(sample.minus ? "-" : "+", 5)
							|| parser.parseInt(3)>sample.stop || parser.parseInt(4)<sample.start){continue;}
						sample.witnesses++;
						witnesses++;
						out.writer.print(buffer.clear().append(sample.id).tab().append(row).nl());
					}
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Reports          ----------------*/
	/*--------------------------------------------------------------*/

	private void writeCounts(){
		models.sort((a, b)->a.name.compareTo(b.name));
		try(Output events=new Output(output+".events.tsv"); Output reasons=new Output(output+".reasons.tsv")){
			events.writer.println("scope\tfamily\tmodel\ttier\traw_events\tnative_excluded_incidence\tpartial_incidence\tcross_family_incidence\tsource_warning_incidence\tbridge_incidence\tqualification");
			reasons.writer.println("scope\tfamily\tmodel\traw_marginal_events\tband_primary\tdiscordant_primary\tpartial_primary\tbridge_witness_primary\tqualification");
			writeCounts(events, reasons, "ALL", ".", ".", all);
			for(final Family family:families){writeCounts(events, reasons, "FAMILY", family.name, ".", family.events);}
			for(final Model model:models){writeCounts(events, reasons, "MODEL", model.family==null ? "." : model.family.name, model.name, model.events);}
		}
		try(Output out=new Output(output+".physical.tsv")){
			out.writer.println("family\tscored_annotations\tmasked_annotations\tband_primary\tdiscordant_primary\tpartial_primary\tprotected_priors\tmissing_score_priors\tmeasured_score_priors\tqualification");
			for(final Family family:families){
				buffer.clear().append(family.name).tab().append(family.scored).tab().append(family.masked);
				for(final long n:family.masks){buffer.tab().append(n);}
				out.writer.print(buffer.tab().append(family.missingPrior+family.measuredPrior).tab().append(family.missingPrior)
					.tab().append(family.measuredPrior).append("\tPHYSICAL_ANNOTATIONS_EXCLUSIVE_PRIMARY_MASK_REASONS\n"));
			}
		}
		long total=0, scored=0, masked=0, missing=0, measured=0;
		for(final long count:all.actual){total+=count;}
		for(final Family family:families){
			scored+=family.scored;
			masked+=family.masked;
			missing+=family.missingPrior;
			measured+=family.measuredPrior;
		}
		try(Output out=new Output(output+".summary.tsv")){
			out.writer.println("inputs\traw_events\tscored_annotations\tmasked_annotations\tprotected_priors\tmissing_score_priors\tmeasured_score_priors\tsample_events\tsample_cross_family_witnesses\tstatus");
			out.writer.print(buffer.clear().append(inputCount).tab().append(total).tab().append(scored).tab().append(masked)
				.tab().append(missing+measured).tab().append(missing).tab().append(measured).tab().append(samples).tab().append(witnesses)
				.append("\tSAVED_REFERENCE_REPORT_COMPLETE_NO_DECISIONS_CHANGED\n"));
		}
	}
	private void writeCounts(final Output events, final Output reasons, final String scope, final String family, final String model, final Counts counts){
		for(int i=0; i<STATES.length; i++){
			buffer.clear().append(scope).tab().append(family).tab().append(model).tab().append(STATES[i]).tab().append(counts.actual[i]);
			for(final long n:counts.incidence[i]){buffer.tab().append(n);}
			events.writer.print(buffer.append("\tCAUSE_INCIDENCES_OVERLAP_AND_ARE_NOT_ADDITIVE\n"));
		}
		buffer.clear().append(scope).tab().append(family).tab().append(model).tab().append(counts.actual[1]);
		for(final long n:counts.reasons){buffer.tab().append(n);}
		reasons.writer.print(buffer.append("\tRAW_EVENTS_EXCLUSIVE_PRIMARY_REASONS\n"));
	}
	private int state(){
		for(int i=0; i<STATES.length; i++){if(parser.termEquals(STATES[i], 22)){return i;}}
		throw new IllegalArgumentException("Unknown saved reference event tier");
	}
	private int reason(final int column){
		for(int i=0; i<REASONS.length; i++){if(parser.termEquals(REASONS[i], column)){return REASON_CLASS[i];}}
		throw new IllegalArgumentException("Unknown serialized marginal reason");
	}
	private static int reason(final String value){
		for(int i=0; i<REASONS.length; i++){if(REASONS[i].equals(value)){return REASON_CLASS[i];}}
		throw new IllegalArgumentException("Unknown serialized physical mask reason: "+value);
	}
	private long nonnegative(final int field){
		final long value=parser.parseLong(field);
		require(value>=0, "Saved accounting counts must be nonnegative");
		return value;
	}
	private boolean bool(final int field){
		require(parser.termEquals("true", field) || parser.termEquals("false", field), "Saved evidence booleans must be explicit");
		return parser.termEquals("true", field);
	}
	private String text(final int field){
		parser.setBounds(field);
		return new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
	}
	private String decoded(final int field){
		parser.setBounds(field);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private String path(final String name){return new File(input, name).getPath();}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}

	private static final class Counts{
		final long[] expected=new long[7], actual=new long[7], reasons=new long[4];
		final long[][] incidence=new long[7][5];
	}
	private static final class Family{
		Family(final String name_){name=name_;}
		final String name;
		final Counts events=new Counts();
		final long[] masks=new long[3];
		long expectedScored, expectedMasked, expectedPriors, scored, masked, missingPrior, measuredPrior;
	}
	private static final class Model{
		Model(final String name_){name=name_;}
		final String name;
		final Counts events=new Counts();
		Family family;
	}
	private static final class Prior{
		Prior(final String input_, final String model_, final int start_, final int stop_){
			input=input_;
			model=model_;
			start=start_;
			stop=stop_;
		}
		final String input, model;
		final int start, stop;
		boolean seen;
	}
	private static final class Sample{
		Sample(final LineParser1 parser){
			id=parser.parseLong(0);
			start=parser.parseInt(3);
			stop=parser.parseInt(4);
			minus=parser.termEquals("-", 5);
			valid=!parser.termEquals("INVALID_INTERVAL", 22);
			final Target key=new Target();
			key.set(parser, 2);
			target=key.copy();
			key.set(parser, 7);
			family=key.copy();
		}
		final long id;
		final int start, stop;
		final boolean minus, valid;
		final Target target, family;
		int witnesses;
	}
	private final String input, output, sampleModel;
	private final int examples;
	private int inputCount, samples, witnesses;
	private final Counts all=new Counts();
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ObjectMap<Target,Family> familyKeys=new ObjectMap<Target,Family>(Target.class, Family.class);
	private final ObjectMap<String,Family> familyNames=new ObjectMap<String,Family>(String.class, Family.class);
	private final ObjectMap<Target,Model> modelKeys=new ObjectMap<Target,Model>(Target.class, Model.class);
	private final ObjectMap<String,Prior> priorIds=new ObjectMap<String,Prior>(String.class, Prior.class);
	private final ArrayList<Sample> sampleRows=new ArrayList<Sample>();
	private final ObjectMap<Target,IntList> samplesByInput=new ObjectMap<Target,IntList>(Target.class, IntList.class);
	private final int[][] sampleCounts=new int[7][8];
	private final Target probe=new Target(), familyProbe=new Target();
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private static final String[] STATES={"CONFIDENT", "MARGINAL", "BELOW_MARGINAL", "UNCERTAIN_FLOOR", "REJECTED_WRONG_FAMILY", "COVERED_BY_PROTECTED_PRIOR", "INVALID_INTERVAL"};
	private static final String[] REASONS={"below_GA_above_marginal_floor", "qualifying_competing_CM", "discordant_source_evidence", "partial_model_support", "tandem_bridge_witness_not_representative"};
	private static final int[] REASON_CLASS={0, 1, 1, 2, 3};
	private static final int SOURCE_BLOCKERS=NcrnaReferenceStrata.OTHER_FAMILY_SOURCE|NcrnaReferenceStrata.UNKNOWN_RFAM
		|NcrnaReferenceStrata.PSEUDO_TRUE|NcrnaReferenceStrata.PSEUDOGENE_TAG|NcrnaReferenceStrata.PSEUDO_BIOTYPE
		|NcrnaReferenceStrata.PARTIAL_TRUE|NcrnaReferenceStrata.EXCEPTION|NcrnaReferenceStrata.UNKNOWN_BOOLEAN|NcrnaReferenceStrata.INCOMPATIBLE_ID;
	private static final String FAMILY_HEADER="family\traw_events\tCONFIDENT\tMARGINAL\tBELOW_MARGINAL\tUNCERTAIN_FLOOR\tREJECTED_WRONG_FAMILY\tCOVERED_BY_PROTECTED_PRIOR\tINVALID_INTERVAL\tscored_annotations\tmasked_annotations\tprotected_priors";
	private static final String MODEL_HEADER="model\traw_events\tCONFIDENT\tMARGINAL\tBELOW_MARGINAL\tUNCERTAIN_FLOOR\tREJECTED_WRONG_FAMILY\tCOVERED_BY_PROTECTED_PRIOR\tINVALID_INTERVAL";
	private static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation";
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
	static final String[] INPUTS={"parameters.tsv", "family_counts.tsv", "event_counts.tsv", "event_provenance.tsv.gz", "annotation_witnesses.tsv.gz", "reference.tsv"};
	static final String[] SUFFIXES={".events.tsv", ".reasons.tsv", ".physical.tsv", ".examples.tsv", ".example_overlaps.tsv", ".summary.tsv"};
}
