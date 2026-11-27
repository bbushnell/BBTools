package prok;

import java.util.ArrayList;
import java.util.Arrays;
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
import structures.DoubleList;

/** Fits a finite OR grid on frozen emitted assignments, selection-only by
* default or pooling both original partitions with explicit refit=t. A 3-D
* rejection prefix counts calls failing every CM without scanning calls per
* triple. Optional mode=fixed predicts declared policies without fitting.
* Neither mode predicts a newly optimized caller path.
* @author Brian Bushnell, Ganyu
*/
public final class SelectCobalaminClanAcceptance{

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, SelectCobalaminClanAcceptance.class, false);
		try{new SelectCobalaminClanAcceptance(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SelectCobalaminClanAcceptance(String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){
				require(b!=null && (b.equalsIgnoreCase("fit") || b.equalsIgnoreCase("fixed")), "Use mode=fit or mode=fixed");
				fixed=b.equalsIgnoreCase("fixed");
			}
			else if(a.equals("fixedpolicy")){fixedPolicyPath=b;}
			else if(a.equals("gains")){gains=b;}
			else if(a.equals("expectedgains")){expectedGains=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("raw")){raw=b;}
			else if(a.equals("baseline")){baseline=b;}
			else if(a.equals("capture")){capture=b;}
			else if(a.equals("sources")){sources=b;}
			else if(a.equals("audit")){audit=b;}
			else if(a.equals("annotations")){annotations=b;}
			else if(a.equals("commoncutoff")){common=NcrnaNumbers.doubleValue(a, b);}
			else if(a.equals("expectedinputs")){expectedInputs=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedreferences")){referenceN=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("expectedmasked")){maskedN=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("baselineunmatched")){unmatchedN=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("quantiles")){quantiles=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("refit")){refit=Parse.parseBoolean(b);}
			else if(a.equals("upperidentity")){upperIdentity=NcrnaNumbers.floatValue(a, b);}
			else{require(parser.parse(arg, a, b), "Unknown OR fitter argument: "+arg);}
		}
		scores=parser.in1;
		assignments=parser.in2;
		out=parser.out1;
		require(scores!=null && assignments!=null && out!=null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1 && expectedInputs>0 && referenceN>0 && maskedN>=0 && unmatchedN>0
			&& (fixed || quantiles>=2 && quantiles<=128 && Double.isFinite(common)),
			"Bind complete saved inputs and positive denominators; fitting additionally requires a finite common cutoff and a bounded 2..128 quantile grid");
		require(fixed ? fixedPolicyPath!=null && gains!=null && expectedGains>=0
			: fixedPolicyPath==null && gains==null && expectedGains<0,
			"Fixed mode requires the frozen policy, actual original-gain roster and explicit gain denominator; fit mode cannot consume them");
		require(upperIdentity==HIGH || upperIdentity==.75f,
			"Use upperidentity=.778 for the historical fit or .75 for the fixed product band; saved score coverage remains .60 to below .778");
		require(!fixed || !refit && upperIdentity==HIGH,
			"refit=t and a changed upperidentity apply only to mode=fit; fixed comparisons retain their original bands and cutoffs");
		require(Tools.testInputFiles(false, true, scores, assignments, raw, sources, audit, annotations), "All event and annotation bindings must exist");
		for(String prefix:new String[]{baseline, capture}){
			require(prefix!=null && Tools.testInputFiles(false, true, prefix+".queries.tsv", prefix+".references.tsv", prefix+".genomes.tsv"),
				"Both arms require complete normalized one-to-one grades");
		}
		if(fixed){
			require(Tools.testInputFiles(false, true, fixedPolicyPath, gains), "Fixed policy and original gain roster must exist");
			RrnaResourceIO.fresh(out, ".fixed_policies.tsv", ".fixed_comparison.tsv", ".fixed_references.tsv.gz",
				".fixed_genomes.tsv.gz", ".fixed_events.tsv.gz", ".fixed_identity.tsv");
		}else{
			RrnaResourceIO.fresh(out, ".axes.tsv", ".surface.tsv.gz", ".choice.tsv", ".policy.tsv", ".comparison.tsv", ".references.tsv", ".genomes.tsv", ".breakdown.tsv");
		}
		Shared.setThreads(1);
	}

	private void process(){
		readSources();
		readAudit();
		readReferences(0, baseline);
		readReferences(1, capture);
		readAnnotations();
		readQueries(0, baseline);
		readQueries(1, capture);
		readRaw();
		readAssignments();
		readScores();
		if(fixed){
			readFixedPolicy();
			readOriginalGains();
			reportFixed();
			System.err.println("COBALAMIN_FIXED_ACCEPTANCE_PREDICTED");
			return;
		}
		fit();
		report();
		System.err.println("COBALAMIN_CLAN_ACCEPTANCE_SELECTED");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input Binding        ----------------*/
	/*--------------------------------------------------------------*/

	private void readSources(){
		try(Input input=new Input(sources)){
			RrnaResourceIO.header(input.reader, "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff", sources);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==5 && p.parseInt(0)==genomes.size(), "Source ordinals must remain consecutive");
				final Genome g=new Genome(p.parseString(2), p.parseString(1));
				require(!byInput.contains(g.id) && !g.id.isEmpty(), "Source IDs must be unique");
				byInput.put(g.id, g);
				genomes.add(g);
			}
		}
		require(genomes.size()==expectedInputs, "The complete input roster, including negatives, is required");
	}
	private void readAudit(){
		try(Input input=new Input(audit)){
			RrnaResourceIO.header(input.reader, "assembly_id\tinput_id\tcanonical_fasta\thistorical_alias_group\tresolved_group\theld_out\tduplicate_original_group\tsampled", audit);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==8, "Resolved-group audit schema changed");
				final Genome g=byInput.get(p.parseString(1));
				if(g==null){continue;}
				require(g.group<0 && p.termEquals(g.assembly, 0) && p.parseInt(4)>=0, "Group audit must bind each selected assembly exactly once");
				g.group=p.parseInt(4);
			}
		}
		for(Genome g:genomes){require(g.group>=0, "Every positive and negative input needs its frozen group");}
	}
	private void readReferences(int side, String prefix){
		int count=0;
		int masks=0;
		try(Input input=new Input(prefix+".references.tsv")){
			RrnaResourceIO.header(input.reader, "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored", prefix);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==7, "Reference schema changed");
				final Genome g=genome(0);
				final String id=p.parseString(1);
				refKey(g.id, id);
				final boolean scored=bool(6);
				if(side==0){
					require(!byRef.contains(probe), "Duplicate source-qualified reference");
					final Ref ref=new Ref(g, id, p.parseString(2), p.parseInt(3), p.parseInt(4), strand(5), scored);
					require(ref.start>0 && ref.stop>=ref.start, "Ordered reference bounds required");
					byRef.put(probe.copy(), ref);
					refs.add(ref);
				}else{
					final Ref ref=byRef.get(probe);
					require(ref!=null && !ref.otherSeen && p.termEquals(ref.seqid, 2) && p.parseInt(3)==ref.start
						&& p.parseInt(4)==ref.stop && strand(5)==ref.strand && scored==ref.scored, "Both arms must retain identical reference identities, geometry and masks");
					ref.otherSeen=true;
				}
				if(scored){
					count++;
				}else{
					masks++;
				}
			}
		}
		require(count==referenceN && masks==maskedN, "Reference denominators changed");
	}
	private void readAnnotations(){
		try(Input input=new Input(annotations)){
			RrnaResourceIO.header(input.reader, SummarizeCobalaminMemberEndpoints.ANNOTATIONS, annotations);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==11, "Use closed annotation-support rows");
				refKey(p.parseString(0), p.parseString(1));
				final Ref ref=byRef.get(probe);
				require(ref!=null && ref.support<0 && ref.scored==bool(2), "Every annotation must bind one exact reference and status");
				ref.support=-1;
				for(int bits=0; bits<SETS.length; bits++){
					if(p.termEquals(SETS[bits], 6)){
						ref.support=bits;
					}
				}
				require(ref.support>=0 && p.termEquals(supportLabel(ref.support), 5), "Annotation label and complete support set must agree");
			}
		}
		for(Ref ref:refs){require(ref.support>=0 && ref.otherSeen, "Complete annotation and paired reference populations required");}
	}
	private void readQueries(int side, String prefix){
		try(Input input=new Input(prefix+".queries.tsv")){
			RrnaResourceIO.header(input.reader, SummarizeCobalaminMemberAssignments.QUERY_HEADER, prefix);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==10, "Normalized query schema changed");
				final Genome g=genome(0);
				final int kind=kind(7);
				Ref ref=null;
				if(kind==0){
					refKey(g.id, p.parseString(9));
					ref=byRef.get(probe);
					require(ref!=null && ref.scored && ref.calls[side]==null && bool(8), "Recorded TP references are eligible and consumed once, never rematched");
				}else{require(p.termEquals(".", 9) || p.termEquals("NA", 9), "Unmatched or masked calls cannot acquire a reference");}
				final Query q=new Query(g, p.parseString(1), p.parseInt(2), p.parseInt(3), strand(4), member(5), first(6), kind, ref);
				require(q.start>0 && q.stop>=q.start && q.model.startsWith("cobalamin_"+MEMBERS[q.member]+"_"), "Query geometry/model must be valid");
				if(ref!=null){
					require(q.seqid.equals(ref.seqid) && q.strand==ref.strand && q.start<=ref.stop && ref.start<=q.stop, "Recorded TP must overlap its bound reference");
					ref.calls[side]=q;
				}
				g.counts[side][kind]++;
				queries[side].add(q);
				if(side==1){
					queryKey(q.g.id, q.seqid, q.start, q.stop, q.strand, q.model);
					require(!byQuery.contains(probe), "Duplicate exact emitted query is ambiguous; use its native assignment rather than silently choosing");
					byQuery.put(probe.copy(), q);
				}
			}
		}
		int seen=0;
		try(Input input=new Input(prefix+".genomes.tsv")){
			RrnaResourceIO.header(input.reader, SummarizeCobalaminMemberAssignments.GENOME_HEADER, prefix);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==5, "Genome accounting schema changed");
				final Genome g=genome(0);
				require(!g.seen[side], "Duplicate genome accounting");
				g.seen[side]=true;
				long total=0;
				for(int k=0; k<3; k++){
					require(g.counts[side][k]==p.parseLong(k+2), "Per-input query outcome counts must reconcile");
					total+=g.counts[side][k];
				}
				require(total==p.parseLong(1), "Genome calls must partition into recorded outcomes");
				seen++;
			}
		}
		require(seen==expectedInputs, "Negative genomes cannot disappear from accounting");
		long unmatched=0;
		for(Genome g:genomes){unmatched+=g.counts[side][2];}
		if(side==0){require(unmatched==unmatchedN, "Objective unmatched denominator must equal the selected ordinary reference");}
	}
	private void readRaw(){
		try(Input input=new Input(raw)){
			RrnaResourceIO.header(input.reader, CobalaminCmAcceptance.HEADER, raw);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==17 && bool(12), "Use the unfiltered native capture with complete event identity");
				final int ordinal=p.parseInt(0);
				require(ordinal>=0 && ordinal<genomes.size(), "Native input ordinal must bind the exact roster");
				final Event e=new Event(line, genomes.get(ordinal), ordinal, p.parseInt(15), p.parseLong(16), p.parseString(1), strand(2), member(3), first(4), p.parseInt(6)+1, p.parseInt(7)+1);
				e.identity=bits(8, false);
				e.own=bits(9, true);
				e.ga=bits(10, false);
				require(e.worker>=0 && e.trace>=0 && Float.isFinite(e.identity) && e.identity>=LOW && e.identity<=1
					&& Float.isFinite(e.ga) && e.ga>0 && e.start>0 && e.stop>=e.start && e.stop<=p.parseInt(5), "Capture identity, geometry and score scale must be valid");
				events.add(e);
			}
		}
	}
	private void readAssignments(){
		try(Input input=new Input(assignments)){
			RrnaResourceIO.header(input.reader, AnalyzeCobalaminCmProposals.EVENT_ASSIGNMENT_HEADER, assignments);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==17, "Use the exact native-event companion export");
				final Event e=event(0);
				require(e.query==null && p.termEquals(e.g.id, 1) && p.parseInt(2)==e.g.group && p.parseInt(3)==e.ordinal
					&& p.parseInt(4)==e.worker && p.parseLong(5)==e.trace && p.termEquals(e.seqid, 6) && strand(7)==e.strand
					&& p.parseInt(8)==e.start && p.parseInt(9)==e.stop && member(10)==e.member && first(11).equals(e.model)
					&& sameBits(bits(14, false), e.identity) && sameBits(bits(15, true), e.own) && sameBits(bits(16, false), e.ga),
					"Assignment must bind the raw event index, native worker/trace ID, model, coordinates and score bits");
				queryKey(e.g.id, e.seqid, e.start, e.stop, e.strand, e.model);
				final Query q=byQuery.get(probe);
				require(q!=null && q.event==null && q.kind==kind(12) && (q.ref==null ? p.termEquals(".", 13) || p.termEquals("NA", 13) : p.termEquals(q.ref.id, 13)),
					"An already-linked native selected event must bind exactly one saved query/reference assignment");
				q.event=e;
				e.query=q;
			}
		}
		for(Query q:queries[1]){require(q.event!=null, "Every emitted query needs its exact selected native event");}
	}
	private void readScores(){
		// HIGH binds the original saved-score capture, independent of the fitted
		// upperIdentity. Lowering bypass must not waive missing/duplicate tail rows.
		try(Input input=new Input(scores)){
			RrnaResourceIO.header(input.reader, ScoreCobalaminClanCandidates.SCORE_HEADER, scores);
			for(byte[] line; (line=input.reader.nextLine())!=null; ){
				p.set(line);
				require(p.terms()==28, "Cross-CM output must contain index, raw17 fields and ten target fields");
				final Event e=event(0);
				require(e.identity<HIGH, "Identity-bypass events must not be substituted into the scored band");
				p.setBounds(1);
				final int from=p.a();
				p.setBounds(17);
				final int to=p.b();
				require(to-from==e.raw.length, "Cross-score raw event width changed");
				for(int i=0; i<e.raw.length; i++){require(e.raw[i]==line[from+i], "Cross-score native event bytes changed");}
				final int target=member(18);
				require((e.seen&(1<<target))==0 && first(19).startsWith("cobalamin_"+MEMBERS[target]+"_"), "Each target CM occurs once with its own anchor library");
				final float score=bits(20, true);
				final float ga=bits(21, false);
				require(Float.isFinite(ga) && ga>0, "Target CM needs its own finite positive GA");
				if(target==e.member){require(sameBits(score, e.own) && sameBits(ga, e.ga) && first(19).equals(e.model), "Own-CM scores must be reused bit for bit");}
				e.margins[target]=Float.isNaN(score) ? Double.NaN : (double)score-ga;
				e.seen|=1<<target;
			}
		}
		for(Event e:events){require(e.seen==(e.identity<HIGH ? 7 : 0), "Exactly three target scores are required for every middle-band event, even unselected candidates");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Fixed Policy Prediction ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable policy axes shared with finite fixtures and saved cost readers.
	 * Floors apply before own support. Cross floors gate only fallback after an
	 * own failure. Identity at the ceiling accepts without a CM. Each CM retains
	 * its frozen cutoff; none of these axes is selected from the current data.
	 */
	public static final class FixedPolicy{
		private FixedPolicy(final String name_, final float floor_, final boolean nativeFloor_,
				final float cross_, final float minorCross_, final float ceiling_, final boolean nativeCeiling_,
				final int directions_, final boolean nativeOrder_){
			name=name_;
			lower=floor_;
			nativeFloor=nativeFloor_;
			cross=cross_;
			minorCross=minorCross_;
			upper=ceiling_;
			nativeCeiling=nativeCeiling_;
			directions=directions_;
			nativeOrder=nativeOrder_;
		}
		public float floor(final int member){
			validateMember(member);
			return nativeFloor ? NATIVE[member] : lower;
		}
		public float crossFloor(final int member){
			validateMember(member);
			return member==0 ? cross : minorCross;
		}
		public float ceiling(final int member){
			validateMember(member);
			return nativeCeiling ? NATIVE[member] : upper;
		}
		public int targetAt(final int position){
			require(position>=0 && position<3, "A fixed cross priority contains exactly three targets");
			return nativeOrder ? position : position==0 ? 1 : position==1 ? 2 : 0;
		}
		public boolean keeps(final int emitted, final float identity, final int passTargets){
			validateMember(emitted);
			require(Float.isFinite(identity) && identity>=0 && identity<=1 && passTargets>=0 && passTargets<=7,
				"Fixed retention requires actual finite identity and three finite-CM pass bits");
			return identity>=floor(emitted) && CobalaminCrossDirectionMasks.retained(passTargets, emitted,
				directions, identity, crossFloor(emitted), ceiling(emitted));
		}
		public final String name;
		public final int directions;
		private final float lower, cross, minorCross, upper;
		private final boolean nativeFloor, nativeCeiling, nativeOrder;
	}
	public static FixedPolicy fixedPolicy(final String name){
		for(FixedPolicy policy:FIXED_POLICIES){if(policy.name.equals(name)){return policy;}}
		throw new IllegalArgumentException("Unknown fixed acceptance policy: "+name);
	}
	public static int finitePassTargets(final double[] margins, final double[] cutoffs){
		require(margins!=null && cutoffs!=null && margins.length==3 && cutoffs.length==3,
			"A fixed OR decision needs exactly three bound score-minus-GA/cutoff axes");
		int bits=0;
		for(int m=0; m<3; m++){
			require(Double.isFinite(cutoffs[m]) && (Double.isFinite(margins[m]) || Double.isNaN(margins[m])
				|| margins[m]==Double.NEGATIVE_INFINITY), "Only unavailable or impossible scores may be nonfinite; cutoffs must remain finite");
			if(Double.isFinite(margins[m]) && margins[m]>=cutoffs[m]){bits|=1<<m;}
		}
		return bits;
	}
	private static void validateMember(final int member){require(member>=0 && member<3, "Member must use the original ordered three-family axis");}

	private void readFixedPolicy(){
		try(Input input=new Input(fixedPolicyPath)){
			RrnaResourceIO.header(input.reader, CobalaminCmAcceptance.OR_POLICY_HEADER, fixedPolicyPath);
			for(int m=0; m<3; m++){
				final byte[] line=input.reader.nextLine();
				require(line!=null, "Frozen policy must retain all three ordered CM cutoffs");
				p.set(line);
				require(p.terms()==6 && p.termEquals("CLAN_OR_ANY_FINITE_PASS", 0) && p.termEquals(MEMBERS[m], 1)
					&& p.parseFloat(2)==LOW && p.parseFloat(3)==HIGH && p.termEquals("FAIL_THIS_CM", 5),
					"Use the original .60-to-.778 finite OR policy, never refitted or altered score coverage");
				fixedCutoffs[m]=p.parseDouble(4);
				require(Double.isFinite(fixedCutoffs[m]), "Frozen target cutoffs must be finite");
			}
			require(input.reader.nextLine()==null, "Additional rows cannot redefine the fixed CM policy");
		}
	}
	private void readOriginalGains(){
		int count=0;
		try(Input input=new Input(gains)){
			RrnaResourceIO.header(input.reader, GAIN_HEADER, gains);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==13, "Use the actual original caller's complete reference-turnover schema");
				if(!p.termEquals("GAINED_GRADED_REFERENCE", 6)){continue;}
				final Genome g=genome(0);
				refKey(g.id, p.parseString(1));
				final Ref ref=byRef.get(probe);
				require(ref!=null && ref.scored && !ref.originalGain && ref.calls[0]==null
					&& p.termEquals(ref.seqid, 2) && p.parseInt(3)==ref.start && p.parseInt(4)==ref.stop && strand(5)==ref.strand,
					"Actual gain IDs must uniquely bind scored references absent from the same ordinary baseline");
				ref.gainMember=member(11);
				ref.gainModel=first(12);
				require(ref.gainModel.startsWith("cobalamin_"+MEMBERS[ref.gainMember]+"_"), "Original gained model must belong to its declared emitted member");
				ref.originalGain=true;
				count++;
			}
		}
		require(count==expectedGains, "Retain the complete declared original gain roster; do not infer it from fitted survivors");
	}

	private void reportFixed(){
		final int n=FIXED_POLICIES.length;
		final long[][][] counts=new long[n][2][3], turnover=new long[n][2][3], presence=new long[n][2][4];
		final long[][][] eventCounts=new long[n][2][4], changes=new long[n][2][6];
		final long[][][] identityCounts=new long[n+1][2][7];
		final long[][][] gainCounts=new long[n][2][5];
		final long[][] ordinary=new long[2][3];
		for(Genome g:genomes){
			g.fixedCounts=new long[n][3];
			for(int k=0; k<3; k++){ordinary[g.group&1][k]+=g.counts[0][k];}
		}
		for(Event e:events){
			e.fixedPass=finitePassTargets(e.margins, fixedCutoffs);
			for(int i=0; i<n; i++){
				final boolean keep=FIXED_POLICIES[i].keeps(e.member, e.identity, e.fixedPass);
				if(keep){e.fixedKeep|=1<<i;}
				final int split=e.g.group&1;
				eventCounts[i][split][0]++;
				if(keep){eventCounts[i][split][1]++;}
				if(e.query==null){
					eventCounts[i][split][2]++;
					if(keep){eventCounts[i][split][3]++;}
				}
			}
		}
		for(Query q:queries[1]){
			final int split=q.g.group&1;
			addIdentity(identityCounts[n][split], q.event.identity, q.kind);
			final boolean pruned=kept(q, 1);
			for(int i=0; i<n; i++){
				final boolean keep=kept(q, i);
				if(keep){
					counts[i][split][q.kind]++;
					q.g.fixedCounts[i][q.kind]++;
					addIdentity(identityCounts[i][split], q.event.identity, q.kind);
				}
				if(keep && !pruned){changes[i][split][q.kind]++;}
				if(!keep && pruned){changes[i][split][3+q.kind]++;}
			}
		}
		try(Output loci=new Output(out+".fixed_references.tsv.gz"); Output gs=new Output(out+".fixed_genomes.tsv.gz")){
			loci.writer.println("policy\tinput_id\treference_id\tsplit\ttransition\toriginal_or_gain\toriginal_gain_member\toriginal_gain_model\tcapture_model\tcapture_matches_original_gain_model\tqualification");
			gs.writer.println("policy\tinput_id\tsplit\tordinary_matched\tpredicted_matched\tordinary_calls\tpredicted_calls\tfirst_detection\tlast_detection\tqualification");
			for(int i=0; i<n; i++){
				for(Ref ref:refs){
					if(!ref.scored){continue;}
					final boolean before=ref.calls[0]!=null;
					final Query q=ref.calls[1];
					final boolean after=q!=null && kept(q, i);
					final int split=ref.g.group&1;
					if(ref.originalGain){
						gainCounts[i][split][0]++;
						if(q!=null){gainCounts[i][split][1]++;}
						if(after){gainCounts[i][split][2]++;}
						if(q!=null && q.member==ref.gainMember && q.model.equals(ref.gainModel)){
							gainCounts[i][split][3]++;
							if(after){gainCounts[i][split][4]++;}
						}
					}
					if(!before && !after && !ref.originalGain){continue;}
					final int state=before && after ? 0 : before ? 1 : after ? 2 : -1;
					if(state>=0){turnover[i][split][state]++;}
					loci.writer.print(row.clear().append(FIXED_POLICIES[i].name).tab().append(ref.g.id).tab().append(ref.id)
						.tab().append(SPLITS[split]).tab().append(state<0 ? "ABSENT_BOTH" : STATES[state]).tab().append(ref.originalGain)
						.tab().append(ref.originalGain ? MEMBERS[ref.gainMember] : "NA").tab().append(ref.originalGain ? ref.gainModel : "NA")
						.tab().append(q==null ? "NA" : q.model).tab().append(!ref.originalGain ? "NA" : q!=null
							&& q.member==ref.gainMember && q.model.equals(ref.gainModel) ? "true" : "false").tab().append(FIXED_QUALIFICATION).nl());
				}
				for(Genome g:genomes){
					final long a=g.counts[0][0], b=g.fixedCounts[i][0];
					final long ac=all(g.counts[0]), bc=all(g.fixedCounts[i]);
					final int split=g.group&1;
					if(a==0 && b>0){presence[i][split][0]++;}
					if(a>0 && b==0){presence[i][split][1]++;}
					if(ac==0 && bc>0){presence[i][split][2]++;}
					if(ac>0 && bc==0){presence[i][split][3]++;}
					gs.writer.print(row.clear().append(FIXED_POLICIES[i].name).tab().append(g.id).tab().append(SPLITS[split])
						.tab().append(a).tab().append(b).tab().append(ac).tab().append(bc).tab().append(a==0 && b>0)
						.tab().append(a>0 && b==0).tab().append(FIXED_QUALIFICATION).nl());
				}
			}
		}
		try(Output output=new Output(out+".fixed_comparison.tsv")){
			output.writer.println("policy\tsplit\tmatched\tunmatched\tmasked\tordinary_matched\tordinary_unmatched\tordinary_masked\tretained\tlost\tgained\toriginal_or_gains\toriginal_gains_in_capture\toriginal_gains_retained\tfirst_detection_genomes\tlast_detection_genomes\tfirst_call_genomes\tlast_call_genomes\tnew_matched_vs_pruned\tnew_masked_vs_pruned\tnew_unmatched_vs_pruned\tdropped_matched_vs_pruned\tdropped_masked_vs_pruned\tdropped_unmatched_vs_pruned\traw_events\traw_accepted\tunselected_raw_events\tunselected_raw_accepted\toriginal_gain_models_equal_capture\toriginal_gain_models_equal_and_retained\tqualification");
			for(int i=0; i<n; i++){
				for(int split=0; split<3; split++){
					require(sum(counts[i], split, 0)==sum(turnover[i], split, 0)+sum(turnover[i], split, 2)
						&& sum(ordinary, split, 0)==sum(turnover[i], split, 0)+sum(turnover[i], split, 1),
						"Frozen one-to-one references must reconcile with independently counted accepted queries");
					row.clear().append(FIXED_POLICIES[i].name).tab().append(SPLITS[split])
						.tab().append(sum(counts[i], split, 0)).tab().append(sum(counts[i], split, 2)).tab().append(sum(counts[i], split, 1))
						.tab().append(sum(ordinary, split, 0)).tab().append(sum(ordinary, split, 2)).tab().append(sum(ordinary, split, 1));
					for(int k=0; k<3; k++){row.tab().append(sum(turnover[i], split, k));}
					for(int k=0; k<3; k++){row.tab().append(sum(gainCounts[i], split, k));}
					for(int k=0; k<4; k++){row.tab().append(sum(presence[i], split, k));}
					for(int k=0; k<6; k++){row.tab().append(sum(changes[i], split, k));}
					for(int k=0; k<4; k++){row.tab().append(sum(eventCounts[i], split, k));}
					row.tab().append(sum(gainCounts[i], split, 3)).tab().append(sum(gainCounts[i], split, 4));
					output.writer.print(row.tab().append(FIXED_QUALIFICATION).nl());
				}
			}
		}
		writeFixedIdentity(identityCounts);
		writeFixedEvents();
		writeFixedPolicies();
	}
	private static boolean kept(final Query q, final int policy){return (q.event.fixedKeep&(1<<policy))!=0;}
	private static void addIdentity(final long[] counts, final float identity, final int kind){
		assert(counts.length==7 && kind>=0 && kind<3) : "Identity tails retain matched/unmatched/masked separately at both thresholds";
		if(kind==0){counts[0]++;}
		final int offset=kind==0 ? 0 : kind==1 ? 2 : 1;
		if(identity>=.75f){counts[1+offset]++;}
		if(identity>=HIGH){counts[4+offset]++;}
	}
	private void writeFixedIdentity(final long[][][] counts){
		try(Output output=new Output(out+".fixed_identity.tsv")){
			output.writer.println("population\tsplit\tall_matched\tmatched_at_least_075\tunmatched_at_least_075\tmasked_at_least_075\tunmasked_at_least_075\tmatched_share_of_unmasked_tail_075\tmatched_coverage_at_least_075\tmatched_at_least_0778\tunmatched_at_least_0778\tmasked_at_least_0778\tunmasked_at_least_0778\tmatched_share_of_unmasked_tail_0778\tmatched_coverage_at_least_0778\tqualification");
			for(int i=0; i<counts.length; i++){
				for(int split=0; split<3; split++){
					final long total=sum(counts[i], split, 0);
					require(sum(counts[i], split, 4)<=sum(counts[i], split, 1) && sum(counts[i], split, 1)<=total,
						"Nested matched tails must conserve the full matched denominator");
					row.clear().append(i==FIXED_POLICIES.length ? "CAPTURE_SELECTED" : FIXED_POLICIES[i].name)
						.tab().append(SPLITS[split]).tab().append(total);
					for(int offset=1; offset<=4; offset+=3){
						final long matched=sum(counts[i], split, offset);
						final long unmatched=sum(counts[i], split, offset+1);
						final long masked=sum(counts[i], split, offset+2);
						final long unmasked=Math.addExact(matched, unmatched);
						row.tab().append(matched).tab().append(unmatched).tab().append(masked).tab().append(unmasked).tab();
						if(unmasked>0){row.appendSlow(matched/(double)unmasked);}
						row.tab();
						if(total>0){row.appendSlow(matched/(double)total);}
					}
					output.writer.print(row.tab().append(FIXED_QUALIFICATION).nl());
				}
			}
		}
	}
	private void writeFixedEvents(){
		try(Output output=new Output(out+".fixed_events.tsv.gz")){
			row.clear().append("event_index\tinput_id\tinput_ordinal0\tworker\ttrace_id\tsplit\tcontig\tstrand\tstart1\tstop1\temitted_member\tmodel\tidentity_bits\tcapture_selected\tfrozen_outcome\tfrozen_reference\tfinite_cm_pass_mask");
			for(FixedPolicy policy:FIXED_POLICIES){row.tab().append(policy.name);}
			output.writer.print(row.append("\tqualification\n"));
			for(int index=0; index<events.size(); index++){
				final Event e=events.get(index);
				final Query q=e.query;
				row.clear().append(index).tab().append(e.g.id).tab().append(e.ordinal).tab().append(e.worker).tab().append(e.trace)
					.tab().append(SPLITS[e.g.group&1]).tab().append(e.seqid).tab().append(e.strand).tab().append(e.start).tab().append(e.stop)
					.tab().append(MEMBERS[e.member]).tab().append(e.model).tab().append(Float.floatToIntBits(e.identity))
					.tab().append(q!=null).tab().append(q==null ? "UNSELECTED_RAW_EVENT" : q.kind==0 ? "TP" : q.kind==1 ? "MASKED" : "UNMATCHED")
					.tab().append(q==null || q.ref==null ? "NA" : q.ref.id).tab().append(e.fixedPass);
				for(int i=0; i<FIXED_POLICIES.length; i++){row.tab().append((e.fixedKeep&(1<<i))!=0);}
				output.writer.print(row.tab().append(FIXED_QUALIFICATION).nl());
			}
		}
	}
	private void writeFixedPolicies(){
		try(Output output=new Output(out+".fixed_policies.tsv")){
			output.writer.println("policy\tmember\tlower_identity\tcross_fallback_floor\tidentity_bypass_ceiling\tallowed_direction_mask\tcross_target_priority\tcm_minus_ga\tqualification");
			for(FixedPolicy policy:FIXED_POLICIES){
				for(int m=0; m<3; m++){
					row.clear().append(policy.name).tab().append(MEMBERS[m]).tab().appendSlow(policy.floor(m))
						.tab().appendSlow(policy.crossFloor(m)).tab().appendSlow(policy.ceiling(m)).tab().append(policy.directions).tab();
					for(int pos=0; pos<3; pos++){
						if(pos>0){row.append(',');}
						row.append(MEMBERS[policy.targetAt(pos)]);
					}
					output.writer.print(row.tab().appendSlow(fixedCutoffs[m]).tab().append(FIXED_QUALIFICATION).nl());
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Grid Fitting        ----------------*/
	/*--------------------------------------------------------------*/

	private void fit(){
		final int n=queries[1].size();
		final double[][] margins=new double[n][];
		final int[] splits=new int[n];
		final int[] kinds=new int[n];
		final boolean[] band=new boolean[n];
		for(int i=0; i<n; i++){
			final Query q=queries[1].get(i);
			margins[i]=q.event.margins;
			splits[i]=q.g.group&1;
			kinds[i]=q.kind;
			band[i]=q.event.identity<upperIdentity;
		}
		final long[][] base=new long[2][3];
		for(Genome g:genomes){
			for(int k=0; k<3; k++){
				base[g.group&1][k]+=g.counts[0][k];
			}
		}
		grid=new Grid(margins, splits, kinds, band, base, referenceN, unmatchedN, quantiles, common, refit);
		final String qualification=fitQualification();
		try(Output axes=new Output(out+".axes.tsv"); Output surface=new Output(out+".surface.tsv.gz"); Output choice=new Output(out+".choice.tsv"); Output policy=new Output(out+".policy.tsv")){
			axes.writer.println("member\tindex\tcutoff\tqualification");
			for(int m=0; m<3; m++){
				for(int i=0; i<grid.axes[m].length; i++){
					axes.writer.print(row.clear().append(MEMBERS[m]).tab().append(i).tab().appendSlow(grid.axes[m][i])
						.tab().append(refit ? "ALL_ELIGIBLE_REFIT_UNMASKED_FINITE_MARGIN_QUANTILES_PLUS_COMMON_AND_REJECT_FINITE_CONTROL"
							: "SELECTION_UNMASKED_FINITE_MARGIN_QUANTILES_PLUS_COMMON_AND_REJECT_FINITE_CONTROL").nl());
				}
			}
			final String header="cutoff_main\tcutoff_second\tcutoff_third\tsplit\tmatched\tunmatched\tmasked\tnormalized_matched_gain\tnormalized_unmatched_drop\tselected\tqualification";
			surface.writer.println(header);
			choice.writer.println(header);
			for(int a=0; a<grid.axes[0].length; a++){
				for(int b=0; b<grid.axes[1].length; b++){
					for(int c=0; c<grid.axes[2].length; c++){
						final boolean best=a==grid.best[0] && b==grid.best[1] && c==grid.best[2];
						for(int split=0; split<3; split++){
							final long matched=grid.kept(a, b, c, split, 0);
							final long unmatched=grid.kept(a, b, c, split, 2);
							row.clear().appendSlow(grid.axes[0][a]).tab().appendSlow(grid.axes[1][b]).tab().appendSlow(grid.axes[2][c])
							.tab().append(SPLITS[split]).tab().append(matched).tab().append(unmatched).tab().append(grid.kept(a, b, c, split, 1))
							.tab().appendSlow((matched-sum(base, split, 0))/(double)referenceN)
							.tab().appendSlow((sum(base, split, 2)-unmatched)/(double)unmatchedN).tab().append(best).tab().append(qualification).nl();
							surface.writer.print(row);
							if(best){choice.writer.print(row);}
						}
					}
				}
			}
			policy.writer.println("acceptance_policy\tmember\tlower_identity\tupper_identity\tcm_minus_ga\tunavailable_score");
			// A nonpositive fitted-population tradeoff does not produce a loadable recommendation.
			if(grid.bestMinimum>0){
				for(int m=0; m<3; m++){
					policy.writer.print(row.clear().append("CLAN_OR_ANY_FINITE_PASS").tab().append(MEMBERS[m])
						.tab().appendSlow(LOW).tab().appendSlow(upperIdentity).tab().appendSlow(grid.axes[m][grid.best[m]])
						.append("\tFAIL_THIS_CM\n"));
				}
			}
		}
	}

	/** Public finite-grid kernel for independent exhaustive fixture comparison. */
	public static final class Grid{
		public Grid(double[][] margins, int[] split, int[] kind, boolean[] band, long[][] baseline, int references, int unmatched, int quantiles, double common){
			this(margins, split, kind, band, baseline, references, unmatched, quantiles, common, false);
		}
		public Grid(double[][] margins, int[] split, int[] kind, boolean[] band, long[][] baseline, int references, int unmatched, int quantiles, double common, boolean refit){
			require(margins.length==split.length && split.length==kind.length && kind.length==band.length && references>0 && unmatched>0
				&& quantiles>=2 && quantiles<=128 && Double.isFinite(common), "Grid inputs and objective scales must be explicit and bounded");
			for(int m=0; m<3; m++){
				final DoubleList values=new DoubleList();
				for(int i=0; i<margins.length; i++){
					require(margins[i].length==3 && split[i]>=0 && split[i]<2 && kind[i]>=0 && kind[i]<3, "Each frozen call has three scores, one group split and one outcome");
					require(Double.isFinite(margins[i][m]) || Double.isNaN(margins[i][m]) || margins[i][m]==Double.NEGATIVE_INFINITY,
						"Positive infinity is malformed, not a passing CM score");
					if(band[i] && (refit || split[i]==0) && kind[i]!=1 && Double.isFinite(margins[i][m])){values.add(margins[i][m]);}
				}
				Arrays.sort(values.array, 0, values.size);
				final DoubleList points=new DoubleList();
				points.add(common);
				if(values.size>0){
					for(int q=0; q<quantiles; q++){points.add(values.get((int)((long)q*(values.size-1)/(quantiles-1))));}
					points.add(Math.nextUp(values.get(values.size-1)));
				}
				Arrays.sort(points.array, 0, points.size);
				int size=0;
				for(int i=0; i<points.size; i++){
					if(size==0 || Double.compare(points.array[i], points.array[size-1])!=0){
						points.array[size++]=points.array[i];
					}
				}
				axes[m]=Arrays.copyOf(points.array, size);
			}
			final int cells=Math.multiplyExact(Math.multiplyExact(axes[0].length, axes[1].length), axes[2].length);
			rejected=new long[6][cells];
			for(int i=0; i<margins.length; i++){
				totals[split[i]][kind[i]]++;
				if(!band[i]){continue;}
				final int a=bin(axes[0], margins[i][0]);
				final int b=bin(axes[1], margins[i][1]);
				final int c=bin(axes[2], margins[i][2]);
				if(a<axes[0].length && b<axes[1].length && c<axes[2].length){rejected[split[i]*3+kind[i]][index(a, b, c)]++;}
			}
			for(long[] cube:rejected){
				for(int a=0; a<axes[0].length; a++){
					for(int b=0; b<axes[1].length; b++){
						for(int c=0; c<axes[2].length; c++){
							final int at=index(a, b, c);
							if(a>0){cube[at]+=cube[index(a-1, b, c)];}
						}
					}
				}
				for(int a=0; a<axes[0].length; a++){
					for(int b=1; b<axes[1].length; b++){
						for(int c=0; c<axes[2].length; c++){
							cube[index(a, b, c)]+=cube[index(a, b-1, c)];
						}
					}
				}
				for(int a=0; a<axes[0].length; a++){
					for(int b=0; b<axes[1].length; b++){
						for(int c=1; c<axes[2].length; c++){
							cube[index(a, b, c)]+=cube[index(a, b, c-1)];
						}
					}
				}
			}
			final int fitSplit=refit ? 2 : 0;
			for(int a=0; a<axes[0].length; a++){
				for(int b=0; b<axes[1].length; b++){
					for(int c=0; c<axes[2].length; c++){
						final double gain=(kept(a, b, c, fitSplit, 0)-sum(baseline, fitSplit, 0))/(double)references;
						final double drop=(sum(baseline, fitSplit, 2)-kept(a, b, c, fitSplit, 2))/(double)unmatched;
						final double minimum=Math.min(gain, drop);
						final double difference=Math.abs(gain-drop);
						if(!chosen || minimum>bestMinimum || minimum==bestMinimum && difference<bestDifference){
							chosen=true;
							bestMinimum=minimum;
							bestDifference=difference;
							best[0]=a;
							best[1]=b;
							best[2]=c;
						}
					}
				}
			}
			require(chosen, "At least one finite grid triple must be evaluated");
		}
		public long kept(int a, int b, int c, int split, int kind){
			require(split>=0 && split<=2 && kind>=0 && kind<3, "Grid queries use selection, judging or their sum");
			if(split==2){return kept(a, b, c, 0, kind)+kept(a, b, c, 1, kind);}
			final long n=totals[split][kind]-rejected[split*3+kind][index(a, b, c)];
			require(n>=0 && n<=totals[split][kind], "Rejection prefixes conserve every frozen outcome");
			return n;
		}
		private int index(int a, int b, int c){return (a*axes[1].length+b)*axes[2].length+c;}
		private static int bin(double[] axis, double value){
			if(!Double.isFinite(value)){return 0;}// Missing and impossible scores fail every finite cutoff.
			int lo=0;
			int hi=axis.length;
			while(lo<hi){
				int mid=(lo+hi)>>>1;
				if(axis[mid]<=value){
					lo=mid+1;
				}else{
					hi=mid;
				}
			}
			return lo;
		}
		public final double[][] axes=new double[3][];
		public final int[] best=new int[3];
		public double bestMinimum;
		private double bestDifference;
		private boolean chosen;
		private final long[][] totals=new long[2][3];
		private final long[][] rejected;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Selected Reports       ----------------*/
	/*--------------------------------------------------------------*/

	private void report(){
		final String qualification=fitQualification();
		final long[][][] counts=new long[3][2][3];
		final long[][][][] breakdown=new long[3][3][3][6];
		for(Query q:queries[0]){q.keep[0]=true;}
		for(Query q:queries[1]){
			final Event e=q.event;
			q.keep[1]=e.identity>=upperIdentity || (Float.isNaN(e.own) ? e.identity>=NATIVE[q.member] : (double)e.own-e.ga>=common);
			q.keep[2]=e.identity>=upperIdentity;
			if(!q.keep[2]){
				for(int m=0; m<3; m++){
					q.keep[2]|=Double.isFinite(e.margins[m]) && e.margins[m]>=grid.axes[m][grid.best[m]];
				}
			}
		}
		for(int policy=0; policy<3; policy++){
			for(Query q:queries[policy==0 ? 0 : 1]){
				if(!q.keep[policy]){continue;}
				counts[policy][q.g.group&1][q.kind]++;
				q.g.kept[policy][q.kind]++;
				if(q.kind!=1){breakdown[policy][0][q.member][q.ref==null ? 5 : supportIndex(q.ref.support)]++;}
				if(q.ref!=null){q.g.supportCounts[policy][q.member][supportIndex(q.ref.support)]++;}
			}
		}
		final long[][][] turnoverCounts=new long[2][2][3];
		final long[][][] genomeCounts=new long[2][2][4];
		try(Output loci=new Output(out+".references.tsv"); Output gs=new Output(out+".genomes.tsv")){
			loci.writer.println("policy\tinput_id\treference_id\ttransition\tbaseline_emitted\tcandidate_emitted\treference_support\tqualification");
			gs.writer.println("policy\tinput_id\tsplit\tbaseline_matched\tcandidate_matched\tbaseline_calls\tcandidate_calls\tfirst_detection\tlast_detection\tqualification");
			for(int policy=1; policy<3; policy++){
				for(Ref ref:refs){
					if(!ref.scored){continue;}
					final Query old=ref.calls[0];
					final Query next=ref.calls[1];
					final boolean a=old!=null;
					final boolean b=next!=null && next.keep[policy];
					if(!a && !b){continue;}
					final int state=a && b ? 0 : a ? 1 : 2;
					turnoverCounts[policy-1][ref.g.group&1][state]++;
					if(state>0){breakdown[policy][state][state==1 ? old.member : next.member][supportIndex(ref.support)]++;}
					loci.writer.print(row.clear().append(POLICIES[policy]).tab().append(ref.g.id).tab().append(ref.id).tab().append(STATES[state])
						.tab().append(a ? MEMBERS[old.member] : "NA").tab().append(b ? MEMBERS[next.member] : "NA")
						.tab().append(SETS[ref.support]).tab().append(qualification).nl());
				}
				for(Genome g:genomes){
					final long a=g.kept[0][0];
					final long b=g.kept[policy][0];
					final long ac=all(g.kept[0]);
					final long bc=all(g.kept[policy]);
					if(a==0 && b>0){genomeCounts[policy-1][g.group&1][0]++;}
					if(a>0 && b==0){genomeCounts[policy-1][g.group&1][1]++;}
					if(ac==0 && bc>0){genomeCounts[policy-1][g.group&1][2]++;}
					if(ac>0 && bc==0){genomeCounts[policy-1][g.group&1][3]++;}
					gs.writer.print(row.clear().append(POLICIES[policy]).tab().append(g.id).tab().append(SPLITS[g.group&1])
						.tab().append(a).tab().append(b).tab().append(ac).tab().append(bc).tab().append(a==0 && b>0).tab().append(a>0 && b==0).tab().append(qualification).nl());
				}
			}
		}
		try(Output summary=new Output(out+".comparison.tsv"); Output bd=new Output(out+".breakdown.tsv")){
			summary.writer.println("policy\tsplit\tmatched\tunmatched\tmasked\tretained\tlost\tgained\tfirst_detection_genomes\tlast_detection_genomes\tfirst_call_genomes\tlast_call_genomes\tqualification");
			bd.writer.println("policy\tpopulation\temitted_member\tannotation_support_label\tcount\tfirst_category_detection_genomes\tlast_category_detection_genomes\tqualification");
			for(int policy=0; policy<3; policy++){
				for(int split=0; split<3; split++){
					row.clear().append(POLICIES[policy]).tab().append(SPLITS[split]).tab().append(sum(counts[policy], split, 0))
					.tab().append(sum(counts[policy], split, 2)).tab().append(sum(counts[policy], split, 1));
					for(int i=0; i<3; i++){row.tab().append(policy==0 ? i==0 ? sum(counts[0], split, 0) : 0 : sum(turnoverCounts[policy-1], split, i));}
					for(int i=0; i<4; i++){row.tab().append(policy==0 ? 0 : sum(genomeCounts[policy-1], split, i));}
					summary.writer.print(row.tab().append(qualification).nl());
				}
				for(int population=0; population<3; population++){
					for(int m=0; m<3; m++){
						for(int label=0; label<6; label++){
							long first=0;
							long last=0;
							if(population==0 && policy>0 && label<5){
								for(Genome g:genomes){
									final long before=g.supportCounts[0][m][label];
									final long after=g.supportCounts[policy][m][label];
									if(before==0 && after>0){first++;}
									if(before>0 && after==0){last++;}
								}
							}
							bd.writer.print(row.clear().append(POLICIES[policy]).tab().append(population==0 ? "ALL_KEPT_UNMASKED" : STATES[population])
								.tab().append(MEMBERS[m]).tab().append(LABELS[label]).tab().append(breakdown[policy][population][m][label])
								.tab().append(first).tab().append(last).tab().append(qualification)
								.append(";CATEGORY_PRESENCE_IS_NOT_WHOLE_GENOME_DETECTION;LABEL_SWITCHES_CAN_CHANGE_CATEGORY\n"));
						}
					}
				}
			}
		}
		for(int split=0; split<2; split++){
			for(int kind=0; kind<3; kind++){
				require(counts[2][split][kind]==grid.kept(grid.best[0], grid.best[1], grid.best[2], split, kind), "Selected-grid counts must match independently enumerated retained calls");
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private String fitQualification(){
		assert(!fixed) : "Only fitted reports describe selection-only or pooled refit scope; fixed mode has its own qualification";
		if(!refit && upperIdentity==HIGH){return QUALIFICATION;}
		return (refit ? REFIT_QUALIFICATION : QUALIFICATION)+";UPPER_IDENTITY="+Float.toString(upperIdentity);
	}

	private Event event(int col){
		final long i=p.parseLong(col);
		require(i>=0 && i<events.size(), "Capture event index outside complete raw file");
		return events.get((int)i);
	}
	private Genome genome(int col){
		final Genome g=byInput.get(p.parseString(col));
		require(g!=null, "Undeclared input ID");
		return g;
	}
	private String first(int col){return RrnaResourceIO.first(p.parseString(col));}
	private boolean bool(int col){
		require(p.termEquals("true", col) || p.termEquals("false", col), "Explicit boolean required");
		return p.termEquals("true", col);
	}
	private char strand(int col){
		require(p.termEquals("+", col) || p.termEquals("-", col), "Explicit strand required");
		return p.termEquals("+", col) ? '+' : '-';
	}
	private int member(int col){
		for(int i=0; i<3; i++){
			if(p.termEquals(MEMBERS[i], col)){
				return i;
			}
		}
		throw new IllegalArgumentException("Unknown emitted/scoring member");
	}
	private int kind(int col){
		if(p.termEquals("TP", col)){
			return 0;
		}
		if(p.termEquals("MASKED", col)){
			return 1;
		}
		require(p.termEquals("UNMATCHED_OTHER", col) || p.termEquals("UNMATCHED_DUPLICATE", col), "Unknown grade outcome");
		return 2;
	}
	private float bits(int col, boolean optional){
		if(optional && p.termEquals("", col)){
			return Float.NaN;
		}
		final float f=Float.intBitsToFloat(p.parseInt(col));
		require(Float.isFinite(f) || optional && f==Float.NEGATIVE_INFINITY, "Invalid encoded score; only blank is unavailable");
		return f;
	}
	private static boolean sameBits(float a, float b){return Float.floatToIntBits(a)==Float.floatToIntBits(b);}
	private void refKey(String input, String id){
		require(!input.isEmpty() && !id.isEmpty() && !id.equals("."), "Explicit source-qualified identity required");
		key.clear().append(input).tab().append(id);
		probe.set(key.array, 0, key.length());
	}
	private void queryKey(String input, String seqid, int start, int stop, char strand, String model){
		key.clear().append(input).tab().append(seqid).tab().append(start).tab().append(stop).tab().append(strand).tab().append(model);
		probe.set(key.array, 0, key.length());
	}
	private static long sum(long[][] values, int split, int kind){return split==2 ? values[0][kind]+values[1][kind] : values[split][kind];}
	private static long all(long[] values){return values[0]+values[1]+values[2];}
	private static int supportIndex(int bits){return bits==0 ? 4 : Integer.bitCount(bits)>1 ? 3 : Integer.numberOfTrailingZeros(bits);}
	private static String supportLabel(int bits){return LABELS[supportIndex(bits)];}
	private static void require(boolean ok, String why){
		if(!ok){
			throw new IllegalArgumentException(why);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Records            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Genome{
		Genome(String id_, String assembly_){
			id=id_;
			assembly=assembly_;
		}
		final String id;
		final String assembly;
		int group=-1;
		final long[][] counts=new long[2][3];
		final long[][] kept=new long[3][3];
		final long[][][] supportCounts=new long[3][3][5];
		final boolean[] seen=new boolean[2];
		long[][] fixedCounts;
	}
	private static final class Ref{
		Ref(Genome g_, String id_, String seqid_, int start_, int stop_, char strand_, boolean scored_){
			g=g_;
			id=id_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			scored=scored_;
		}
		final Genome g;
		final String id;
		final String seqid;
		final int start;
		final int stop;
		final char strand;
		final boolean scored;
		final Query[] calls=new Query[2];
		int support=-1;
		boolean otherSeen;
		boolean originalGain;
		int gainMember=-1;
		String gainModel;
	}
	private static final class Query{
		Query(Genome g_, String seqid_, int start_, int stop_, char strand_, int member_, String model_, int kind_, Ref ref_){
			g=g_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			member=member_;
			model=model_;
			kind=kind_;
			ref=ref_;
		}
		final Genome g;
		final String seqid;
		final String model;
		final int start;
		final int stop;
		final int member;
		final int kind;
		final char strand;
		final Ref ref;
		final boolean[] keep=new boolean[3];
		Event event;
	}
	private static final class Event{
		Event(byte[] raw_, Genome g_, int ordinal_, int worker_, long trace_, String seqid_, char strand_, int member_, String model_, int start_, int stop_){
			raw=raw_.clone();
			g=g_;
			ordinal=ordinal_;
			worker=worker_;
			trace=trace_;
			seqid=seqid_;
			strand=strand_;
			member=member_;
			model=model_;
			start=start_;
			stop=stop_;
		}
		final byte[] raw;
		final Genome g;
		final int ordinal;
		final int worker;
		final int member;
		final int start;
		final int stop;
		final long trace;
		final String seqid;
		final String model;
		final char strand;
		float identity;
		float own;
		float ga;
		int seen;
		int fixedPass, fixedKeep;
		final double[] margins={Double.NaN, Double.NaN, Double.NaN};
		Query query;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	private String scores;
	private String assignments;
	private String raw;
	private String baseline;
	private String capture;
	private String sources;
	private String audit;
	private String annotations;
	private String out;
	private boolean fixed;
	private boolean refit;
	private float upperIdentity=HIGH;
	private String fixedPolicyPath, gains;
	private int expectedGains=-1;
	private final double[] fixedCutoffs=new double[3];
	private double common=Double.NaN;
	private int expectedInputs=-1;
	private int referenceN=-1;
	private int maskedN=-1;
	private int unmatchedN=-1;
	private int quantiles=64;
	private Grid grid;
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Ref> refs=new ArrayList<Ref>();
	private final ArrayList<Event> events=new ArrayList<Event>();
	@SuppressWarnings("unchecked")
	private final ArrayList<Query>[] queries=new ArrayList[]{new ArrayList<Query>(), new ArrayList<Query>()};
	private final ObjectMap<String, Genome> byInput=new ObjectMap<String, Genome>(String.class, Genome.class);
	private final ObjectMap<Target, Ref> byRef=new ObjectMap<Target, Ref>(Target.class, Ref.class);
	private final ObjectMap<Target, Query> byQuery=new ObjectMap<Target, Query>(Target.class, Query.class);
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder();
	private final ByteBuilder row=new ByteBuilder();
	private final LineParser1 p=new LineParser1('\t');

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final float LOW=.60f;
	private static final float HIGH=.778f;
	private static final float[] NATIVE={.63f, .66f, .66f};
	private static final FixedPolicy[] FIXED_POLICIES={
		new FixedPolicy("FULL_OR", LOW, false, LOW, LOW, HIGH, false, 63, true),
		new FixedPolicy("PRUNED", LOW, false, LOW, LOW, HIGH, false, 30, false),
		new FixedPolicy("A", LOW, false, .63f, .66f, HIGH, false, 30, false),
		new FixedPolicy("B", LOW, false, LOW, LOW, .75f, false, 30, false),
		new FixedPolicy("C", LOW, false, .63f, .66f, .75f, false, 30, false),
		new FixedPolicy("D1", LOW, false, LOW, LOW, HIGH, true, 30, false),
		new FixedPolicy("D2", LOW, true, LOW, LOW, HIGH, false, 30, false),
		new FixedPolicy("E", LOW, false, .63f, .63f, .75f, false, 63, false),
		new FixedPolicy("E_MINOR_66", LOW, false, .63f, .66f, .75f, false, 63, false)
	};
	private static final String GAIN_HEADER="input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_member\tbaseline_model\tbaseline_model_in_subset\tcandidate_member\tcandidate_model";
	private static final String FIXED_QUALIFICATION="FROZEN_CUTOFFS_NO_FIT;FULL_CAPTURE_SCORES_060_0778;CAPTURE_SELECTED_ASSIGNMENTS_ONLY_FOR_GRADE;RAW_UNSELECTED_EVENTS_SEPARATE;ORIGINAL_GAIN_REFERENCE_ROSTER;NO_REMATCH_OR_PATH_RESELECTION;NOT_ACTUAL_CALLER_RECALL";
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] SETS={"NA", "RF00174", "RF01482", "RF00174,RF01482", "RF01689", "RF00174,RF01689", "RF01482,RF01689", "RF00174,RF01482,RF01689"};
	private static final String[] LABELS={"RF00174", "RF01482", "RF01689", "AMBIGUOUS", "UNKNOWN", "NO_MATCHED_REFERENCE"};
	private static final String[] SPLITS={"SELECTION", "JUDGING", "ALL"};
	private static final String[] POLICIES={"ORDINARY", "COMMON_OWN_CM_WITH_NATIVE_MISSING_FALLBACK", "BEST_OF_CLAN_OR_UNAVAILABLE_REJECTS"};
	private static final String[] STATES={"RETAINED", "LOST", "GAINED"};
	private static final String QUALIFICATION="SELECTION_ONLY_QUANTILE_GRID_APPROXIMATION;JUDGING_UNCHANGED;FINITE_CUTOFFS;MISSING_AND_IMPOSSIBLE_FAIL_OR;FIXED_EMITTED_ASSIGNMENTS_NOT_NEW_CALLER_PATH;NO_REMATCH";
	private static final String REFIT_QUALIFICATION="ALL_ELIGIBLE_SAVED_REFIT_QUANTILE_GRID_APPROXIMATION;SELECTION_AND_JUDGING_POOLED_FOR_FIT;ORIGINAL_GROUPS_AND_PARTITIONS_RETAINED_NOT_HELDOUT;FINITE_CUTOFFS;MISSING_AND_IMPOSSIBLE_FAIL_OR;FIXED_EMITTED_ASSIGNMENTS_NOT_NEW_CALLER_PATH;NO_REMATCH";
}
