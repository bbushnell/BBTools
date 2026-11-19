package prok;

import java.util.ArrayList;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Conditional endpoint replay on the existing comparison's validated paired references.
 * Reads only saved ungated scores, never genomes, networks or a new matching rule.
 * One thread owns retained selected candidates and source-qualified keys. Rejected
 * candidates remain in conservation counts; unknown eligibility remains unavailable.
 * Scalar replay chooses the captured second-end context after gating the first.
 * Empty numeric report cells retain their presence evidence: unknown eligibility
 * has no replayed error, zero replayable/subset counts have no RMS, unchanged
 * errors have no better fraction, and an unchosen margin has no numeric value.
 * Existing population counts and status columns distinguish each case.
 * This diagnostic cannot establish detection or a later path-selection outcome.
 * @author Brian Bushnell, Raiden
 */
final class CobalaminBoundaryReplay{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	CobalaminBoundaryReplay(final String trace, final String bindings, final int radius_, final float[] margins_){
		require(radius_>0, "Replay needs the same declared positive radius as its captured caller");
		radius=radius_;
		margins=margins_==null ? DEFAULT_MARGINS.clone() : margins_.clone();
		marginText=new String[margins.length];
		require(margins.length>1 && margins[0]==0, "Diagnostic margins must begin with zero and include a positive comparison");
		for(int i=0; i<margins.length; i++){
			require(Float.isFinite(margins[i]) && margins[i]>=0 && (i==0 || margins[i]>margins[i-1]),
				"Diagnostic margins must be finite, nonnegative and strictly increasing");
			marginText[i]=i==0 ? "0" : Float.toString(margins[i]);
		}
		readBindings(bindings);
		readTrace(trace);
		for(final Source source:sources){
			require(source.inputSeen, "Trace omitted a declared input, including an input with no calls: "+source.id);
			for(int i=0; i<5; i++){
				require(source.candidates[i]==source.totals[i], "Candidate rows disagree with contig conservation totals for "+source.id);
			}
		}
	}

	/** The normalized grade must retain exactly the declared source population. */
	void input(final String id){
		final Source source=byInput.get(id);
		require(source!=null && !source.gradeSeen, "Normalized input is absent from replay bindings or repeated: "+id);
		source.gradeSeen=true;
	}

	/** Binds every ON query, including masked/unmatched calls, by exact final identity. */
	Call query(final String input, final String contig, final char strand, final String member,
			final String model, final int start, final int stop){
		final Source source=byInput.get(input);
		require(source!=null, "A replay query needs its declared input identity");
		callKey(source.ordinal, contig, strand, member, CompareCobalaminBenchmarkGrades.modelId(model), start, stop);
		final Call call=byCall.get(probe);
		require(call!=null && !call.linked, "ON query lacks a unique exact selected trace: "+input+" "+contig+" "+start+".."+stop);
		call.linked=true;
		return call;
	}

	/** Completes exact trace-to-grade conservation before any diagnostic selection. */
	void linked(){
		for(final Source source:sources){require(source.gradeSeen, "Normalized grade omitted replay input "+source.id);}
		for(final Call call:calls){require(call.linked, "Selected trace lacks its exact normalized ON query");}
	}

	/** Retains one reference already matched in BOTH grades by the unchanged matcher. */
	void pair(final Call call, final String reference, final int cohort, final int refStart, final int refStop,
			final long off5, final long off3){
		require(call!=null && call.linked && !call.paired && cohort>=0 && cohort<3 && refStart>0 && refStop>=refStart,
			"Replay pairs require distinct validated winners and their fixed OFF-member cohort");
		call.paired=true;
		pairs.add(new Pair(call, reference, cohort, refStart, refStop, off5, off3));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved Inputs         ----------------*/
	/*--------------------------------------------------------------*/

	private void readBindings(final String path){
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff", path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==5 && p.parseInt(0)==sources.size(), "Replay bindings require complete contiguous source ordinals");
				final Source source=new Source(sources.size(), p.parseString(2), p.parseString(3));
				require(!byInput.contains(source.id), "Replay input IDs must be unique");
				byInput.put(source.id, source);
				sources.add(source);
			}
		}
		require(!sources.isEmpty(), "Replay cannot guess an absent source population");
	}

	/** C and its S rows are written atomically by the producer. Only selected calls
	 * retain score arrays. No transient candidate or per-field String is allocated
	 * for the rejected population; its fates still contribute to conservation.
	 */
	private void readTrace(final String path){
		Call active=null;
		try(Input input=new Input(path)){
			final byte[] version=input.reader.nextLine();
			require(version!=null, "Boundary trace requires a version header");
			p.set(version);
			require(p.terms()==2 && p.termEquals("#ncrna_boundary_trace", 0)
				&& (p.termEquals("1", 1) || p.termEquals("2", 1) || p.termEquals("3", 1)), "Unsupported boundary-trace version");
			final boolean legacy=p.termEquals("1", 1);
			final boolean legacyMaximum=!p.termEquals("3", 1);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){continue;}
				p.set(line);
				if(!p.termEquals("S", 0)){
					if(active!=null){active.validate();}
					active=null;
				}
				final Source source=source(p.parseInt(1));
				if(p.termEquals("I", 0)){
					require(p.terms()==3 && !source.inputSeen && p.termEquals(source.path, 2), "Trace assembly path differs from its actual ordered source binding");
					source.inputSeen=true;
				}else if(p.termEquals("C", 0)){
					require(p.terms()==37 && p.parseInt(29)==radius, "Unsupported trace candidate schema or radius");
					final int maximum=NcrnaTextIO.maximumLength(p, 22, legacyMaximum);
					final boolean positional=p.termEquals("1-HOT", 28);
					require(positional || p.termEquals("SLIDING", 28), "Unknown captured network kind");
					if(kind==null){kind=positional ? "1-HOT" : "SLIDING";}
					require(p.termEquals(kind, 28), "One replay compares one complete network kind/radius");
					final boolean returned=flag(31), selected=flag(32), known=flag(14), eligible=flag(16);
					require(!selected || returned, "Selected trace must be a scavenger return");
					source.candidates[0]++;
					if(returned){source.candidates[1]++;}
					if(selected){source.candidates[2]++;}
					if(!known){source.candidates[3]++;}
					if(eligible && p.parseInt(30)==0){source.candidates[4]++;}
					if(selected){
						active=new Call(source, p, positional, radius, legacy, maximum);
						callKey(source.ordinal, active.contig, active.strand, active.member, active.model, active.finalStart, active.finalStop);
						require(!byCall.contains(probe), "Exact selected trace identity is ambiguous; do not guess a final-call association");
						byCall.put(probe.copy(), active);
						calls.add(active);
					}
				}else if(p.termEquals("S", 0)){
					require(p.terms()==9, "Unsupported trace score schema");
					if(active!=null){active.scores(p);}
				}else if(p.termEquals("T", 0)){
					require(p.terms()==11, "Unsupported trace contig-total schema");
					for(int i=0; i<5; i++){
						final long value=p.parseLong(i+6);
						require(value>=0, "Trace conservation totals cannot be negative");
						source.totals[i]=Math.addExact(source.totals[i], value);
					}
				}else{throw new IllegalArgumentException("Unknown boundary trace record");}
			}
			if(active!=null){active.validate();}
		}
	}

	private Source source(final int ordinal){
		require(ordinal>=0 && ordinal<sources.size(), "Trace source ordinal is outside its bound population");
		return sources.get(ordinal);
	}
	private boolean flag(final int column){
		require(p.termEquals("0", column) || p.termEquals("1", column), "Trace flags must be explicit0/1");
		return p.termEquals("1", column);
	}
	private void callKey(final int ordinal, final String contig, final char strand, final String member,
			final String model, final int start, final int stop){
		key.clear().append(ordinal).tab().append(contig).tab().append(strand).tab().append(member).tab().append(model)
			.tab().append(start).tab().append(stop);
		probe.set(key.array, 0, key.length());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reports            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reports the same reference cohort at all declared diagnostic margins. Unknown
	 * eligibility makes selection UNAVAILABLE, never an implicitly unchanged end.
	 * Selection uses pooled full squared errors across both ends, not mean RMS.
	 */
	void write(final String prefix){
		final Stats[][][] stats=new Stats[margins.length][4][2];
		for(int margin=0; margin<margins.length; margin++){
			for(int cohort=0; cohort<4; cohort++){
				for(int end=0; end<2; end++){stats[margin][cohort][end]=new Stats();}
			}
		}
		long unknown=0;
		final ByteBuilder b=new ByteBuilder();
		try(Output output=new Output(prefix+".replay_loci.tsv")){
			output.writer.println("input_id\treference_id\toff_member\tcaptured_member\tmargin\tstatus\toff_d5\treplay_d5\toff_d3\treplay_d3\tqualification");
			for(final Pair pair:pairs){
				final Call call=pair.call;
				if(!call.known){unknown++;}
				for(int margin=0; margin<margins.length; margin++){
					b.clear().append(call.source.id).tab().append(pair.reference).tab().append(MEMBERS[pair.cohort]).tab()
						.append(call.member).tab().append(marginText[margin]).tab().append(call.known ? "REPLAYED" : "UNKNOWN_ELIGIBILITY").tab()
						.append(pair.five).tab();
					if(call.known){
						final long ends=call.endpoints(margins[margin]);
						final int five=(int)(ends>>>32), three=(int)ends;
						final long z5=call.strand=='+' ? (long)five+1-pair.start : (long)pair.stop-call.length+five;
						final long z3=call.strand=='+' ? (long)three+1-pair.stop : (long)pair.start-call.length+three;
						stats[margin][pair.cohort][0].add(pair.five, z5, radius, five!=call.start);
						stats[margin][3][0].add(pair.five, z5, radius, five!=call.start);
						stats[margin][pair.cohort][1].add(pair.three, z3, radius, three!=call.stop);
						stats[margin][3][1].add(pair.three, z3, radius, three!=call.stop);
						b.append(z5).tab().append(pair.three).tab().append(z3);
					}else{b.tab().append(pair.three).tab();}
					output.writer.print(b.append("\tFIXED_PAIRED_REFERENCE_COHORT_CONDITIONAL_ON_CAPTURED_CALLS\n"));
				}
			}
		}
		try(Output output=new Output(prefix+".replay_endpoints.tsv")){
			output.writer.println("kind\tradius_bp\tmargin\toff_member_cohort\tend\tpaired_total\treplayable_n\tunknown_n\toff_rms_nt\treplay_rms_nt\tbetter\tworse\tunchanged\tbetter_fraction\texact\twithin3\twithin10\tmoved_from_anchor\tunmoved_from_anchor\tworse_off_within_radius\tworse_off_beyond_radius\toff_within_radius_n\toff_within_radius_rms_nt\treplay_same_subset_rms_nt\tqualification");
			for(int margin=0; margin<margins.length; margin++){
				for(int cohort=0; cohort<4; cohort++){
					long total=0;
					for(final Pair pair:pairs){if(cohort==3 || pair.cohort==cohort){total++;}}
					for(int end=0; end<2; end++){
						final Stats s=stats[margin][cohort][end];
						require(s.n==s.better+s.worse+s.unchanged && s.worse==s.worseWithin+s.worseBeyond,
							"Every replayable end requires exactly one error-change and worsening-radius disposition");
						b.clear().append(kind==null ? "NO_CAPTURED_CANDIDATES" : kind).tab().append(radius).tab().append(marginText[margin]).tab()
							.append(cohort==3 ? "cobalamin" : MEMBERS[cohort]).tab().append(end==0 ? "5prime" : "3prime")
							.tab().append(total).tab().append(s.n).tab().append(total-s.n).tab();
						rms(b, s.offSquares, s.n);
						rms(b.tab(), s.onSquares, s.n);
						b.tab().append(s.better).tab().append(s.worse).tab().append(s.unchanged).tab();
						ratio(b, s.better, s.worse);
						b.tab().append(s.exact).tab().append(s.near3).tab().append(s.near10).tab().append(s.moved).tab().append(s.n-s.moved)
							.tab().append(s.worseWithin).tab().append(s.worseBeyond).tab().append(s.withinN).tab();
						rms(b, s.offWithinSquares, s.withinN);
						rms(b.tab(), s.onWithinSquares, s.withinN);
						output.writer.print(b.append("\tFULL_PAIRED_RMS_PRIMARY_UNKNOWN_EXPLICIT_NOT_DETECTION_EVIDENCE\n"));
					}
				}
			}
		}
		int chosen=-1;
		final long ungatedBetter=better(stats[0]), ungatedWorse=worse(stats[0]);
		if(unknown==0 && ungatedBetter+ungatedWorse>0){
			for(int margin=1; margin<margins.length; margin++){
				final long better=better(stats[margin]), worse=worse(stats[margin]);
				if(better+worse==0 || (double)better*(ungatedBetter+ungatedWorse)<(double)ungatedBetter*(better+worse)){continue;}
				if(chosen<0 || squares(stats[margin])<squares(stats[chosen])
					|| squares(stats[margin])==squares(stats[chosen]) && margins[margin]>margins[chosen]){chosen=margin;}
			}
		}
		try(Output output=new Output(prefix+".replay_selection.tsv")){
			output.writer.println("kind\tradius_bp\tpaired_references\tunknown_references\tselected_margin\tstatus\tqualification");
			b.clear().append(kind==null ? "NO_CAPTURED_CANDIDATES" : kind).tab().append(radius).tab().append(pairs.size()).tab().append(unknown).tab();
			if(chosen>=0){b.append(marginText[chosen]);}
			output.writer.print(b.tab().append(unknown>0 ? "UNAVAILABLE_INCOMPLETE_CAPTURE" : chosen>=0 ? "SELECTED_DIAGNOSTIC_MARGIN" : "NONE_QUALIFIES")
				.append("\tPOOLED_FULL_SQUARED_ERROR_DEFINED_RATIO_AT_LEAST_UNGATED_TIES_LARGER_MARGIN_CONFIRM_IN_CALLER\n"));
		}
		try(Output output=new Output(prefix+".replay_accounting.tsv")){
			output.writer.println("input_id\tcandidates\tscavenger_returned\tselected_and_query_linked\tunknown_candidates\teligible_without_scores");
			for(final Source source:sources){
				b.clear().append(source.id);
				for(final long count:source.candidates){b.tab().append(count);}
				output.writer.print(b.nl());
			}
		}
	}

	private static long better(final Stats[][] stats){return stats[3][0].better+stats[3][1].better;}
	private static long worse(final Stats[][] stats){return stats[3][0].worse+stats[3][1].worse;}
	private static double squares(final Stats[][] stats){return stats[3][0].onSquares+stats[3][1].onSquares;}
	private static void rms(final ByteBuilder b, final double squares, final long n){
		require(Double.isFinite(squares) && squares>=0 && n>=0, "Measured coordinate-square sums and counts must be finite/nonnegative");
		if(n>0){b.append(Math.sqrt(squares/n), 8);}
	}
	private static void ratio(final ByteBuilder b, final long better, final long worse){
		if(better+worse>0){b.append(better/(double)(better+worse), 8);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Score Replay        ----------------*/
	/*--------------------------------------------------------------*/

	/** One retained final selected candidate. Sense coordinates and score buffers
	 * belong to this result; finalStart/finalStop are one-based genomic positions.
	 */
	static final class Call{
		Call(final Source source_, final LineParser1 p, final boolean positional_, final int radius_, final boolean legacy_, final int maximum){
			source=source_;
			legacy=legacy_;
			positional=positional_;
			radius=radius_;
			worker=p.parseInt(2);
			id=p.parseLong(3);
			contig=p.parseString(5);
			final int orientation=p.parseInt(6);
			require(orientation==0 || orientation==1, "Trace strand must preserve the original sense-to-genomic transform");
			strand=orientation==0 ? '+' : '-';
			member=NcrnaFamilyRegistry.sharedBoundaryMember(p.parseString(7)).member;
			model=p.parseString(8);
			length=p.parseInt(9);
			start=p.parseInt(10);
			stop=p.parseInt(11);
			known=p.parseInt(14)==1;
			fixed=known ? p.parseInt(12) : 0;
			eligible=p.parseInt(16)==1;
			min=p.parseInt(21);
			max=maximum;
			observedStart=p.parseInt(23);
			observedStop=p.parseInt(24);
			finalStart=Math.addExact(p.parseInt(25), 1);
			finalStop=Math.addExact(p.parseInt(26), 1);
			mask=p.parseInt(30);
			require(length>0 && start>=0 && stop>=start && stop<length && finalStart>0 && finalStop>=finalStart && finalStop<=length,
				"Replay needs valid original sense anchors and final genomic bounds");
			require(min>0 && max>=min && (known ? fixed>=0 && fixed<=3 : legacy ? p.parseInt(12)==-1 : p.length(12)==0), "Known physical-end masks and accepted length bounds must remain explicit");
		}
		void scores(final LineParser1 p){
			require(known && eligible && count<3 && p.parseInt(1)==source.ordinal && p.parseInt(2)==worker && p.parseLong(3)==id,
				"Retained score window belongs to a different or ineligible selected candidate");
			final int end=p.parseInt(4), context=p.parseInt(5);
			require(end>=0 && end<2 && (fixed&(1<<end))==0 && context>=0 && context<=3,
				"Score windows cannot move fixed physical ends or invent a context");
			p.setBounds(8);
			final byte[] line=p.line();
			final float[] values=new float[2*radius+1];
			final boolean[] available=new boolean[values.length];
			int from=p.a(), index=0;
			for(int to=from; to<=p.b(); to++){
				if(to==p.b() || line[to]==','){
					require(index<values.length && (!legacy || to>from), "Score window width must equal2R+1 with versioned optional cells");
					if(to>from){
						final float value=Float.intBitsToFloat(Parse.parseInt(line, from, to));
						if(!(legacy && Float.isNaN(value))){
							require(Float.isFinite(value), "A measured trace score must be finite");
							values[index]=value;
							available[index]=true;
						}
					}
					index++;
					from=to+1;
				}
			}
			require(index==values.length && available[radius], "Captured score window must retain its finite unchanged position");
			windows[count++]=new Window(end, context, p.parseInt(6), p.parseInt(7), values, available);
		}
		void validate(){
			int actualMask=0;
			for(int i=0; i<count; i++){actualMask|=1<<windows[i].end;}
			require(actualMask==mask, "Selected candidate's saved window mask differs from its score rows");
			if(!known){require(!eligible && count==0, "Unknown eligibility cannot synthesize replay scores");return;}
			final long zero=endpoints(0);
			require((int)(zero>>>32)==observedStart && (int)zero==observedStop,
				"Decoded zero-margin replay must exactly reconstruct observed ungated inference");
			final int a=strand=='+' ? observedStart+1 : length-observedStop;
			final int z=strand=='+' ? observedStop+1 : length-observedStart;
			require(a==finalStart && z==finalStop, "Final trace must be the recorded original sense-to-genomic transformation");
		}
		long endpoints(final float margin){
			require(known, "Unknown eligibility has no replayed endpoint");
			int five=start, three=stop;
			if(eligible && count>0){
				if(positional){
					for(int i=0; i<count; i++){
						final Window w=windows[i];
						require(w.context==NcrnaBoundaryTrace.INDEPENDENT && w.start==start && w.stop==stop,
							"Positional ends use the same original post-trim span");
						final int offset=select(w.values, w.available, radius, margin, true);
						if(w.end==0){five+=offset;}else{three+=offset;}
					}
					final long span=(long)three-five+1;
					if(five<0 || three>=length || span<min || span>max){
						five=start;
						three=stop;
					}
				}else{
					final Window first=windows[0];
					require(first.context==NcrnaBoundaryTrace.FIRST && first.start==start && first.stop==stop,
						"Scalar replay uses the captured first-end order and original span");
					final int offset=select(first.values, first.available, radius, margin, false);
					if(first.end==0){five+=offset;}else{three+=offset;}
					boolean found=false;
					for(int i=1; i<count; i++){
						final Window w=windows[i];
						if(w.context==(offset==0 ? NcrnaBoundaryTrace.HELD : NcrnaBoundaryTrace.MOVED)){
							require(w.end!=first.end && w.start==five && w.stop==three, "Second-end features must use the gated first-end context");
							final int second=select(w.values, w.available, radius, margin, false);
							if(w.end==0){five+=second;}else{three+=second;}
							found=true;
							break;
						}
					}
					require(fixed!=0 || found, "Two-end scalar replay needs its captured held/moved context");
				}
			}
			return ((long)five<<32)|(three&0xffffffffL);
		}
		final Source source;
		final String contig, member, model;
		final char strand;
		final boolean known, eligible, positional, legacy;
		final int worker, radius, length, start, stop, fixed, min, max, observedStart, observedStop, finalStart, finalStop, mask;
		final long id;
		final Window[] windows=new Window[3];
		int count;
		boolean linked, paired;
	}

	/** Uses the same shared confidence helper as inference; legacy scalar ties hold
	 * current while legacy positional ties use the first maximum.
	 */
	static int select(final float[] values, final boolean[] available, final int radius, final float margin, final boolean positional){
		int best=-1;
		for(int i=0; i<values.length; i++){
			if(!available[i]){continue;}
			require(Float.isFinite(values[i]), "Infinite score cannot establish confidence");
			if(best<0 || values[i]>values[best]){best=i;}
		}
		require(best>=0, "A captured window must have an available score");
		final int legacy=positional || values[best]>values[radius] ? best : radius;
		return NcrnaEndpointConfidence.select(values, radius, legacy, margin, available)-radius;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Retained State      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final int ordinal_, final String id_, final String path_){
			ordinal=ordinal_;
			id=id_;
			path=path_;
		}
		final int ordinal;
		final String id, path;
		final long[] candidates=new long[5], totals=new long[5];
		boolean inputSeen, gradeSeen;
	}
	private static final class Window{
		Window(final int end_, final int context_, final int start_, final int stop_, final float[] values_, final boolean[] available_){
			end=end_;
			context=context_;
			start=start_;
			stop=stop_;
			values=values_;
			available=available_;
		}
		final int end, context, start, stop;
		final float[] values;
		final boolean[] available;
	}
	private static final class Pair{
		Pair(final Call call_, final String reference_, final int cohort_, final int start_, final int stop_, final long five_, final long three_){
			call=call_;
			reference=reference_;
			cohort=cohort_;
			start=start_;
			stop=stop_;
			five=five_;
			three=three_;
		}
		final Call call;
		final String reference;
		final int cohort, start, stop;
		final long five, three;
	}
	private static final class Stats{
		void add(final long before, final long after, final int radius, final boolean moved_){
			final double a=Math.abs((double)before), z=Math.abs((double)after);
			n++;
			offSquares+=a*a;
			onSquares+=z*z;
			if(z<a){better++;}
			else if(z>a){
				worse++;
				if(a<=radius){worseWithin++;}else{worseBeyond++;}
			}else{unchanged++;}
			if(z==0){exact++;}
			if(z<=3){near3++;}
			if(z<=10){near10++;}
			if(moved_){moved++;}
			if(a<=radius){
				withinN++;
				offWithinSquares+=a*a;
				onWithinSquares+=z*z;
			}
			assert(Double.isFinite(offSquares) && Double.isFinite(onSquares)) : "Coordinate differences must be widened before squaring for finite full-cohort RMS";
		}
		long n, better, worse, unchanged, exact, near3, near10, moved, worseWithin, worseBeyond, withinN;
		double offSquares, onSquares, offWithinSquares, onWithinSquares;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final int radius;
	private final float[] margins;
	/** Preformatted once, preserving requested float values without per-row formatting. */
	private final String[] marginText;
	private String kind;
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Call> calls=new ArrayList<Call>();
	private final ArrayList<Pair> pairs=new ArrayList<Pair>();
	private final ObjectMap<String, Source> byInput=new ObjectMap<String, Source>(String.class, Source.class);
	private final ObjectMap<Target, Call> byCall=new ObjectMap<Target, Call>(Target.class, Call.class);
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder();
	private final LineParser1 p=new LineParser1('\t');
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final float[] DEFAULT_MARGINS={0, .005f, .02f, .1f};
}
