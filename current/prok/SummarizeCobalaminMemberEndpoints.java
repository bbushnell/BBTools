package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Joins closed member/support assignments to the same capture's saved endpoint
 * errors. Neither support nor endpoints are recalculated. Comparators retain all
 * matching calls in that capture, without applying the selected member's filter.
 * @author Brian Bushnell, Diona
 */
public final class SummarizeCobalaminMemberEndpoints{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeCobalaminMemberEndpoints.class, false);
		try{
			new SummarizeCobalaminMemberEndpoints(pp.args).process();
			pp.outstream.println("COBALAMIN_MEMBER_ENDPOINTS_SUMMARIZED");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeCobalaminMemberEndpoints(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("witnesses")){witnesses=b;}
			else if(a.equals("metrics")){metrics=b;}
			else if(a.equals("capture")){capture=b;}
			else if(a.equals("expectedkept")){expectedKept=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedselected")){expectedSelected=NcrnaNumbers.intValue(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown endpoint summary argument: "+arg);}
		}
		details=parser.in1;
		annotations=parser.in2;
		out=parser.out1;
		require(details!=null && annotations!=null && witnesses!=null && metrics!=null && out!=null
			&& capture!=null && capture.matches("[A-Za-z0-9_]+") && expectedKept>0
			&& expectedSelected>0 && expectedSelected<=expectedKept && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Use complete accepted details/annotations, same-capture witnesses/metrics, explicit cohort counts and fresh out=");
		require(Tools.testInputFiles(false, true, details, annotations, witnesses, metrics)
			&& Tools.testOutputFiles(false, false, false, out), "Read closed inputs and write one fresh table");
		Shared.setThreads(1);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Processing         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		readMetrics();
		readAnnotations();
		readSelected();
		readWitnesses();
		write();
	}

	private void readMetrics(){
		try(Input input=new Input(metrics)){
			final byte[] header=input.reader.nextLine();
			require(header!=null, "Capture metrics need a header");
			p.set(header);
			final int width=p.terms();
			final int scored=column("scored_references");
			final int masks=column("masked_references");
			final int matched=column("matched");
			final int exact=column("both_exact");
			final int inputs=column("inputs");
			final byte[] line=input.reader.nextLine();
			require(line!=null, "Capture metrics need one complete row");
			p.set(line);
			require(p.terms()==width && input.reader.nextLine()==null, "Use one complete capture metric row");
			scoredReferences=count(scored);
			maskedReferences=count(masks);
			captureMatches=count(matched);
			captureBoth=count(exact);
			captureInputs=count(inputs);
			require(captureInputs>0 && captureMatches<=scoredReferences && captureBoth<=captureMatches,
				"Capture denominators must be consistent");
		}
	}

	private void readAnnotations(){
		long scored=0;
		long masked=0;
		try(Input input=new Input(annotations)){
			RrnaResourceIO.header(input.reader, ANNOTATIONS, annotations);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==11 && (p.termEquals("true", 2) || p.termEquals("false", 2)), "Use complete accepted annotation rows");
				final int support=support(6);
				final String label=support==0 ? "UNKNOWN" : Integer.bitCount(support)>1 ? "AMBIGUOUS" : MEMBERS[Integer.numberOfTrailingZeros(support)];
				require(p.termEquals(label, 5), "Accepted support set and label disagree");
				final Reference ref=new Reference(p.parseString(0), p.parseString(1), p.termEquals("true", 2), support);
				key(ref.input, ref.id);
				require(!references.contains(probe), "Duplicate input-qualified annotation");
				references.put(probe.copy(), ref);
				order.add(ref);
				if(ref.scored){scored++;}else{masked++;}
			}
		}
		require(scored==scoredReferences && masked==maskedReferences, "Annotation population differs from the capture");
	}

	private void readSelected(){
		int kept=0;
		int selected=0;
		try(Input input=new Input(details)){
			RrnaResourceIO.header(input.reader, DETAILS, details);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==14, "Use complete accepted member-assignment details");
				if(!p.termEquals("MEMBER_PREDICTED_ALL_MATCHES", 0)){continue;}
				final Reference ref=reference(1, 4);
				require(ref!=null && ref.scored && ref.keptModel==null, "Selected member references must be unique and scored");
				require(p.termEquals("RF01482", 2) && support(8)==ref.support, "Selected member/support differs from accepted annotations");
				ref.keptModel=p.parseString(3);
				require(ref.keptModel.startsWith("cobalamin_RF01482_"), "Selected exact model must retain the emitted member");
				kept++;
				if(ref.support==1){selected++;}
			}
		}
		require(kept==expectedKept && selected==expectedSelected, "The closed selected-member cohort changed");
	}

	private void readWitnesses(){
		long matched=0;
		long both=0;
		try(Input input=new Input(witnesses)){
			RrnaResourceIO.header(input.reader, WITNESSES, witnesses);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==12, "Use complete normalized one-to-one endpoint witnesses");
				final Reference ref=reference(0, 1);
				require(ref!=null && ref.scored && !ref.seen, "Each endpoint witness needs one unique scored reference");
				ref.seen=true;
				final String member=p.parseString(2);
				require(member.equals(MEMBERS[0]) || member.equals(MEMBERS[1]) || member.equals(MEMBERS[2]), "Unknown emitted witness member");
				final String model=RrnaResourceIO.first(p.parseString(3));
				require(model.startsWith("cobalamin_"+member+"_"), "Saved witness model/member mismatch");
				final long d5=p.parseLong(10);
				final long d3=p.parseLong(11);
				require(d5>=-Integer.MAX_VALUE && d5<=Integer.MAX_VALUE && d3>=-Integer.MAX_VALUE && d3<=Integer.MAX_VALUE,
					"Saved endpoint differences must be finite coordinate differences");
				final boolean exact=d5==0 && d3==0;
				matched++;
				if(exact){both++;}
				if(ref.keptModel!=null){
					require(member.equals("RF01482") && ref.keptModel.equals(model), "Selected model differs from its original saved capture assignment");
					if(ref.support==1){add(0, exact);}
				}
				if(member.equals("RF00174") && (ref.support&1)!=0){
					add(2, exact);
					if(ref.support==1){add(1, exact);}
				}
			}
		}
		for(Reference ref:order){require(ref.keptModel==null || ref.seen, "A selected member reference is missing from the saved capture");}
		require(matched==captureMatches && both==captureBoth, "Complete capture matched/both-exact counts disagree with saved metrics");
		require(totals[0]==expectedSelected, "Every requested selected call needs its original endpoint witness");
	}

	/*--------------------------------------------------------------*/
	/*----------------             Output           ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		try(Output output=new Output(out)){
			output.writer.println("cohort\tcapture\temitted_member\treference_support\tselection\tmatched\tboth_exact\tboth_exact_fraction\tcapture_inputs\tcapture_matched\tcapture_both_exact\twitnesses\tqualification");
			for(int i=0; i<3; i++){
				row.clear().append(COHORTS[i]).tab().append(capture).tab().append(i==0 ? "RF01482" : "RF00174")
					.tab().append(i==2 ? "INCLUDES_RF00174" : "RF00174_ONLY")
					.tab().append(i==0 ? "FROZEN_MEMBER_PREDICTED_KEPT" : "ALL_MATCHED_CALLS_IN_SAME_CAPTURE")
					.tab().append(totals[i]).tab().append(exacts[i]).tab();
				if(totals[i]>0){row.append((double)exacts[i]/totals[i], 10);}else{row.append("NA");}
				output.writer.print(row.tab().append(captureInputs).tab().append(captureMatches).tab().append(captureBoth)
					.tab().append(witnesses).append("\tSAVED_ENDPOINT_ERRORS;SUPPORT_NOT_SUBTYPE_TRUTH;COMPARATORS_UNFILTERED;MAIN_COHORTS_NESTED_NOT_ADDITIVE;NO_REMATCH\n"));
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private void add(int cohort, boolean exact){
		totals[cohort]++;
		if(exact){exacts[cohort]++;}
	}
	private int column(String name){
		int found=-1;
		for(int i=0; i<p.terms(); i++){
			if(p.termEquals(name, i)){
				require(found<0, "Duplicate metric column");
				found=i;
			}
		}
		require(found>=0, "Missing metric column: "+name);
		return found;
	}
	private long count(int column){
		final long value=p.parseLong(column);
		require(value>=0, "Negative metric count");
		return value;
	}
	private int support(int column){
		for(int bits=0; bits<SUPPORT.length; bits++){if(p.termEquals(SUPPORT[bits], column)){return bits;}}
		throw new IllegalArgumentException("Unknown or noncanonical retained support set");
	}
	private Reference reference(int inputColumn, int idColumn){
		key(p.parseString(inputColumn), p.parseString(idColumn));
		return references.get(probe);
	}
	private void key(String input, String id){
		require(!input.isEmpty() && !id.isEmpty() && !id.equals("NA"), "Reference keys must be explicit");
		key.clear().append(input).tab().append(id);
		probe.set(key.array, 0, key.length());
	}
	private static void require(boolean value, String why){if(!value){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Reference{
		Reference(String input_, String id_, boolean scored_, int support_){
			input=input_;
			id=id_;
			scored=scored_;
			support=support_;
		}
		final String input;
		final String id;
		final boolean scored;
		final int support;
		String keptModel;
		boolean seen;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private String details;
	private String annotations;
	private String witnesses;
	private String metrics;
	private String capture;
	private String out;
	private int expectedKept=-1;
	private int expectedSelected=-1;
	private long scoredReferences;
	private long maskedReferences;
	private long captureMatches;
	private long captureBoth;
	private long captureInputs;
	private final long[] totals=new long[3];
	private final long[] exacts=new long[3];
	private final ObjectMap<Target, Reference> references=new ObjectMap<Target, Reference>(Target.class, Reference.class);
	private final ArrayList<Reference> order=new ArrayList<Reference>();
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder();
	private final ByteBuilder row=new ByteBuilder();
	private final LineParser1 p=new LineParser1('\t');

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] SUPPORT={"NA", "RF00174", "RF01482", "RF00174,RF01482", "RF01689", "RF00174,RF01689", "RF01482,RF01689", "RF00174,RF01482,RF01689"};
	private static final String[] COHORTS={"SECOND_PREDICTED_KEPT_MAIN_ONLY", "MAIN_ALL_CAPTURE_MAIN_ONLY", "MAIN_ALL_CAPTURE_MAIN_SUPPORTED"};
	public static final String ANNOTATIONS="input_id\treference_id\tscored\trepresentative_raw\trepresentative_label\tannotation_label\tconfident_members\tmarginal_members\tprotected_prior\tmeasurement_models\twitnesses";
	public static final String DETAILS="cohort\tinput_id\temitted_member\tmodel\treference_id\trepresentative_raw\trepresentative_label\tannotation_label\tconfident_members\tmarginal_members\tprotected_prior\tbaseline_reference_matched\tbaseline_same_member\tqualification";
	public static final String WITNESSES="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
}
