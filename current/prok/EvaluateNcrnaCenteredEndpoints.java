package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Compares saved uncentered and centered predictions on every original trial.
 * Adapted from template.A_SampleByteFile. No inference, fitting, replacement
 * draws or model selection occurs. Final invalid proposals are counted, NOT
 * held: this preserves the historical independent-end assay. Caller joint holds
 * are a separate integration test. Saved context bounds are not contig bounds.
 * @author Brian, Ganyu
 */
public final class EvaluateNcrnaCenteredEndpoints{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, EvaluateNcrnaCenteredEndpoints.class, false);
		try{
			final EvaluateNcrnaCenteredEndpoints tool=new EvaluateNcrnaCenteredEndpoints(pp.args);
			tool.process();
			pp.outstream.println(tool.geometryOnly ? "NCRNA_CENTERED_GEOMETRY_COMPLETE_NO_SCORES_NO_INFERENCE"
				: "NCRNA_CENTERED_ENDPOINT_COMPARISON_COMPLETE_SAVED_SCORES_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Explicit accepted counts prevent a truncated population from becoming a result. */
	public EvaluateNcrnaCenteredEndpoints(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("centered")){centered=b;}
			else if(a.equals("ledger")){ledger=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("expected")){expected=Parse.parseKMG(b);}
			else if(a.equals("usable")){usable=Parse.parseKMG(b);}
			else if(a.equals("geometryonly")){geometryOnly=Parse.parseBoolean(b);}
			else if(a.equals("end")){end="5prime".equals(b) ? 0 : "3prime".equals(b) ? 1 : -1;}
			else{require(parser.parse(arg, a, b), "Unknown centered-comparison argument: "+arg);}
		}
		baseline=parser.in1;
		out=parser.out1;
		if(geometryOnly){
			require(baseline==null && centered==null && usable==-1 && ledger!=null && loci!=null && out!=null && end>=0 && expected>0
				&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
				"geometryonly=t requires only ledger= loci= end= expected= out=prefix and complete unsampled populations");
			require(Tools.testInputFiles(false, true, ledger, loci) && Tools.testOutputFiles(false, false, false, out+".geometry.tsv")
				&& Tools.testForDuplicateFiles(true, ledger, loci, out+".geometry.tsv"), "Geometry inputs and fresh output must be distinct and readable/writable");
			return;
		}
		require(baseline!=null && ledger!=null && loci!=null && out!=null && end>=0 && expected>0
			&& usable>=0 && usable<=expected && (usable==0 ? centered==null : centered!=null),
			"Require in=saved_baseline centered=saved_centered ledger=all_trials loci=native_loci end= expected= usable= out=prefix; omit centered only when usable=0");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Compare complete original populations into fresh outputs without sampling");
		require(Tools.testInputFiles(false, true, baseline, centered, ledger, loci)
			&& Tools.testOutputFiles(false, false, false, out+".summary.tsv", out+".trials.tsv")
			&& Tools.testForDuplicateFiles(true, baseline, centered, ledger, loci, out+".summary.tsv", out+".trials.tsv"),
			"Comparison inputs and fresh outputs must be readable, writable and distinct");
	}

	/** All rows are joined in their original order; held trials consume no centered score. */
	public void process(){
		if(geometryOnly){
			processGeometry();
			return;
		}
		final Geometry geometry=new Geometry(loci);
		final Stats[] results=new Stats[POLICIES.length];
		for(int i=0; i<results.length; i++){results[i]=new Stats();}
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder text=new ByteBuilder();
		final long[] starts=new long[POLICIES.length], stops=new long[POLICIES.length];
		final int[] invalid=new int[POLICIES.length];
		long count=0, scored=0, missing=0, reverted=0, unsupported=0, rawInvalid=0;
		try(Input audit=new Input(ledger); Scores oldScores=new Scores(baseline, end);
				Scores newScores=centered==null ? null : new Scores(centered, end);
				Output trials=new Output(out+".trials.tsv")){
			RrnaResourceIO.header(audit.reader, NcrnaFixedEndpointTrials.AUDIT_HEADER, ledger);
			trials.writer.println("original_row1\tlocus_id\tsample\tend\tdisposition\toffset_supported\tpolicy\tfinal_start0\tfinal_stop0\tsigned_error_bp\tinvalid_proposal");
			for(byte[] line=audit.reader.nextLine(); line!=null; line=audit.reader.nextLine()){
				p.set(line);
				count++;
				final Bounds bounds=validateTrial(p, geometry, count, expected);
				final long rawStart=p.parseInt(7), rawStop=p.parseInt(8);
				final long truth=end==0 ? p.parseInt(9) : p.parseInt(10);
				final long centerStart=p.parseLong(12), centerStop=p.parseLong(13);
				final int minimum=p.parseInt(20);
				final boolean supported=p.parseInt(19)==1;
				final int centerInvalid=invalid(centerStart, centerStop, bounds.length, minimum);
				final boolean vector=p.termEquals("VECTORIZED", 16);
				final boolean noFeature=p.termEquals("CENTER_MISSING_FEATURE_HOLD", 16);
				final boolean held=centerInvalid!=0;
				if(!supported){unsupported++;}
				if(held){reverted++;}
				if(noFeature){missing++;}
				if(invalid(rawStart, rawStop, bounds.length, minimum)!=0){rawInvalid++;}
				oldScores.next(p, false);
				if(vector){
					require(newScores!=null && scored<usable, "A vectorized trial needs one accepted centered score row");
					newScores.next(p, true);
					scored++;
				}
				Arrays.fill(starts, rawStart);
				Arrays.fill(stops, rawStop);
				for(int i=3; i<POLICIES.length; i++){
					starts[i]=held ? rawStart : centerStart;
					stops[i]=held ? rawStop : centerStop;
				}
				move(starts, stops, 1, oldScores.selected0-10, end);
				move(starts, stops, 2, oldScores.selected01-10, end);
				if(vector){
					move(starts, stops, 4, newScores.selected0-10, end);
					move(starts, stops, 5, newScores.selected01-10, end);
				}
				final long rawError=(end==0 ? rawStart : rawStop)-truth;
				for(int i=0; i<POLICIES.length; i++){
					invalid[i]=invalid(starts[i], stops[i], bounds.length, minimum);
					final long error=(end==0 ? starts[i] : stops[i])-truth;
					final int oldArm=i==4 ? 1 : i==5 ? 2 : 0;
					final long oldError=(end==0 ? starts[oldArm] : stops[oldArm])-truth;
					results[i].add(error, rawError, oldError, invalid[i]);
					if(i>=3){
						if(held){results[i].reverted++;}
						else if(starts[i]!=starts[3] || stops[i]!=stops[3]){results[i].centerPlusNetwork++;}
						else{results[i].centerOnly++;}
					}
					text.clear().append(count).tab();
					field(text, p, 0).tab();
					field(text, p, 6).tab().append(end==0 ? "5prime" : "3prime").tab();
					field(text, p, 16).tab().append(supported ? 1 : 0).tab().append(POLICIES[i])
						.tab().append(starts[i]).tab().append(stops[i]).tab().append(error).tab().append(INVALID[invalid[i]]).nl();
					trials.writer.print(text);
				}
			}
			require(count==expected && scored==usable && scored+missing+reverted==count,
				"Scored, missing-feature and reverted original trials must partition the exact accepted denominator");
			oldScores.finish(expected);
			if(newScores!=null){newScores.finish(usable);}
			trials.writer.print(text.clear().append("#COMPLETE\t").append(count).nl());
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("policy\tend\ttrials\texact\twithin3\twithin10\texact_pct\twithin3_pct\trms_bp\traw_rms_bp\trms_minus_raw_bp\tcomparator\tcomparator_rms_bp\trms_minus_comparator_bp\tbetter_than_raw\tworse_than_raw\ttied_raw\tbetter_than_comparator\tworse_than_comparator\ttied_comparator\tinvalid_proposals\tcrossed\toutside_context\tbelow_minimum\traw_invalid_trials\tvectorized\tmissing_features\tcenter_reverted\tunsupported\tcenter_only\tcenter_plus_network\treverted_to_raw\tqualification");
			for(int i=0; i<POLICIES.length; i++){
				results[i].row(text.clear(), POLICIES[i], end, i==4 ? POLICIES[1] : i==5 ? POLICIES[2] : POLICIES[0],
					rawInvalid, scored, missing, reverted, unsupported);
				summary.writer.print(text);
			}
		}
	}

	/** Joint raw/center counts distinguish inherited short spans, newly invalid
	 * centers and repairs. This reads no scores and repeats no inference or metric.
	 */
	private void processGeometry(){
		final Geometry geometry=new Geometry(loci);
		final long[][] counts=new long[4][4];
		final LineParser1 p=new LineParser1('\t');
		long count=0;
		try(Input input=new Input(ledger)){
			RrnaResourceIO.header(input.reader, NcrnaFixedEndpointTrials.AUDIT_HEADER, ledger);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				count++;
				final Bounds bounds=validateTrial(p, geometry, count, expected);
				final int minimum=p.parseInt(20);
				final int raw=invalid(p.parseLong(7), p.parseLong(8), bounds.length, minimum);
				final int center=invalid(p.parseLong(12), p.parseLong(13), bounds.length, minimum);
				counts[raw][center]++;
			}
		}
		require(count==expected, "Geometry cross-tab must retain the complete accepted original trial population");
		long total=0;
		try(Output output=new Output(out+".geometry.tsv")){
			output.writer.println("end\traw_geometry\tcentered_geometry\ttrials\ttrial_pct");
			final ByteBuilder text=new ByteBuilder();
			for(int raw=0; raw<4; raw++){
				for(int center=0; center<4; center++){
					final long n=counts[raw][center];
					total+=n;
					output.writer.print(text.clear().append(end==0 ? "5prime" : "3prime").tab().append(INVALID[raw])
						.tab().append(INVALID[center]).tab().append(n).tab().append(100.0*n/count, 8).nl());
				}
			}
		}
		assert(total==count) : "The disjoint raw/center geometry cross-tab must conserve every original trial";
	}

	/** Shared validation keeps the diagnostic and prediction paths on one contract. */
	private static Bounds validateTrial(final LineParser1 p, final Geometry geometry, final long count, final long expected){
		require(count<=expected && p.terms()==21 && p.parseLong(11)==count,
			"The all-trial ledger must preserve every accepted original row exactly once in order");
		for(int i=0; i<21; i++){require(p.length(i)>0, "Centered ledger fields cannot be empty");}
		final Bounds bounds=geometry.get(p);
		final long rawStart=p.parseInt(7), rawStop=p.parseInt(8);
		final long centerStart=p.parseLong(12), centerStop=p.parseLong(13);
		final int m5=p.parseInt(14), m3=p.parseInt(15), minimum=p.parseInt(20);
		final long support=p.parseLong(17), floor=p.parseLong(18);
		final boolean supported=p.parseInt(19)==1;
		require((p.termEquals("0", 19) || p.termEquals("1", 19)) && support>=0 && floor>0
			&& supported==(support>=floor) && (supported || m5==0 && m3==0) && minimum>0
			&& centerStart==rawStart-m5 && centerStop==rawStop-m3,
			"Ledger centers and support must agree with its declared applied training medians");
		final int bad=invalid(centerStart, centerStop, bounds.length, minimum);
		require(bad!=0 ? p.termEquals(HOLDS[bad], 16) : p.termEquals("VECTORIZED", 16) || p.termEquals("CENTER_MISSING_FEATURE_HOLD", 16),
			"Centered disposition must agree with the retained joint geometry, or explicitly retain a valid feature-missing center");
		return bounds;
	}

	/** Exclusive geometry categories, in exactly the input producer's order. */
	private static int invalid(final long start, final long stop, final int length, final int minimum){
		assert(length>0 && minimum>0) : "Validated source geometry and family minimum bound all proposal classifications";
		return stop<start ? 1 : start<0 || stop>=length ? 2 : stop-start+1<minimum ? 3 : 0;
	}
	private static void move(final long[] starts, final long[] stops, final int arm, final int shift, final int end){
		assert(shift>=-10 && shift<=10 && (end==0 || end==1)) : "Only one reachable endpoint moves in the historical conditional assay";
		if(end==0){starts[arm]+=shift;}else{stops[arm]+=shift;}
	}
	private static ByteBuilder field(final ByteBuilder b, final LineParser1 p, final int term){
		p.setBounds(term);
		return b.append(p.line(), p.a(), p.b()-p.a());
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	/** Only small native metadata are retained; no sequence or score population is copied. */
	private static final class Geometry{
		Geometry(final String path){
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, NcrnaLearningPopulation.NATIVE_HEADER, path);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					require(p.terms()==14, "Centered replay requires native occurrence geometry with observed sense-oriented flanks");
					for(int i=0; i<14; i++){require(p.length(i)>0, "Native occurrence fields must be explicit rather than silently interpreted as zero");}
					final String id=p.parseString(0);
					probe.set(p, 0);
					require(!id.isEmpty() && p.termEquals(id, 1) && !index.contains(probe), "Native original IDs must be explicit and unique");
					final int start=p.parseInt(9), stop=p.parseInt(10), core=p.parseInt(11);
					final int left=p.parseInt(12), right=p.parseInt(13);
					require(start>0 && stop>=start && (long)stop-start+1==core && left>=0 && right>=0
						&& (long)left+core+right<=Integer.MAX_VALUE && (p.termEquals("+", 8) || p.termEquals("-", 8))
						&& (p.termEquals("-", 8) ? right : left)<start && p.length(2)>0 && p.length(3)>0,
						"Native context geometry must fit source coordinates without inventing a contig length");
					index.put(probe.copy(), rows.size());
					rows.add(new Bounds(left, core, left+core+right, p.parseString(2), p.parseString(3)));
				}
			}
			require(!rows.isEmpty(), "Native geometry cannot be empty");
		}
		Bounds get(final LineParser1 p){
			probe.set(p, 0);
			final int id=index.get(probe);
			require(id>=0, "Every original trial must join a known native context");
			final Bounds b=rows.get(id);
			require(p.termEquals(b.source, 1) && p.termEquals(b.assembly, 2) && p.termEquals("evaluation", 3)
				&& p.parseInt(4)>=0 && p.length(5)>0 && p.parseInt(6)>=0 && p.parseInt(7)>=0
				&& p.parseInt(8)>=p.parseInt(7) && p.parseInt(8)<b.length && p.parseInt(9)==b.left
				&& (long)p.parseInt(10)-b.left+1==b.core, "Original trial identity and raw/truth coordinates must match the unchanged native context");
			return b;
		}
		final NcrnaReferenceGeometry.Target probe=new NcrnaReferenceGeometry.Target();
		final ObjectIntMap<NcrnaReferenceGeometry.Target> index=new ObjectIntMap<NcrnaReferenceGeometry.Target>(NcrnaReferenceGeometry.Target.class);
		final ArrayList<Bounds> rows=new ArrayList<Bounds>();
	}
	private static final class Bounds{
		Bounds(final int left_, final int core_, final int length_, final String source_, final String assembly_){
			left=left_;
			core=core_;
			length=length_;
			source=source_;
			assembly=assembly_;
		}
		final int left, core, length;
		final String source, assembly;
	}

	/** Validates saved FP32 scores and rederives stored policies without another inference. */
	private static final class Scores implements AutoCloseable{
		Scores(final String path, final int end_){
			end=end_;
			input=new Input(path);
			try{
				RrnaResourceIO.header(input.reader, "#NCRNA_FULL_ERROR_SCORES_V1", path);
				RrnaResourceIO.header(input.reader, "#coding\tFP32_RAW_SIGNED_INT", path);
				for(final String name:new String[]{"#matrix", "#rows", "#network", "#end"}){
					final byte[] line=input.reader.nextLine();
					require(line!=null, "Saved score provenance is incomplete");
					p.set(line);
					require(p.terms()==2 && p.termEquals(name, 0) && p.length(1)>0
						&& (!name.equals("#end") || p.termEquals(end==0 ? "5prime" : "3prime", 1)), "Saved score provenance/end mismatch");
				}
				byte[] line=input.reader.nextLine();
				require(line!=null, "Missing saved score header");
				p.set(line);
				if(p.termEquals("#featurevariant", 0)){
					require(p.terms()==2 && p.termEquals("baseline", 1), "Centering changes no first-stage feature layout");
					line=input.reader.nextLine();
				}
				final ByteBuilder header=new ByteBuilder(NcrnaFixedEndpointTrials.ROW_HEADER)
					.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
				for(int i=0; i<22; i++){header.append("\tscore_bits_").append(i);}
				require(Arrays.equals(line, header.toBytes()), "Saved score columns must match the qualified 22-output evaluator");
			}catch(RuntimeException | Error failure){
				try{input.close();}catch(RuntimeException | Error close){failure.addSuppressed(close);}
				throw failure;
			}
		}
		void next(final LineParser1 ledger, final boolean centered){
			final byte[] line=input.reader.nextLine();
			require(line!=null, "Saved scores ended before their all-trial ledger");
			p.set(line);
			require(p.terms()==38, "Saved score rows require all 38 fields");
			for(int i=0; i<38; i++){require(p.length(i)>0, "Saved score fields cannot be empty");}
			for(int i=0; i<11; i++){
				if(centered && (i==7 || i==8)){
					require(p.parseLong(i)==ledger.parseLong(i+5), "Centered score coordinates disagree with the retained center");
				}else{
					p.setBounds(i);
					ledger.setBounds(i);
					require(p.length(i)==ledger.length(i), "Saved score identity width differs from the original trial");
					for(int j=0; j<p.length(i); j++){
						require(line[p.a()+j]==ledger.line()[ledger.a()+j], "Saved score and original trial identities differ");
					}
				}
			}
			final long correction=end==0 ? p.parseLong(9)-p.parseLong(7) : p.parseLong(10)-p.parseLong(8);
			final int target=correction<-10 || correction>10 ? 21 : (int)correction+10;
			require(p.parseInt(11)==target && p.parseLong(12)==correction, "Saved score targets must retain real signed distances, including outside the window");
			for(int i=0; i<22; i++){
				scores[i]=Float.intBitsToFloat(p.parseInt(16+i));
				require(Float.isFinite(scores[i]), "Saved scores must be finite");
			}
			selected0=NcrnaEndpointConfidence.selectFull(scores, 10, 0);
			selected01=NcrnaEndpointConfidence.selectFull(scores, 10, .1f);
			require(p.parseInt(13)==selected0 && p.parseInt(14)==selected01 && p.parseInt(15)==Tools.maxIndex(scores),
				"Saved policy choices must agree with the preserved score bits");
			count++;
		}
		void finish(final long expected){
			final byte[] line=input.reader.nextLine();
			require(line!=null, "Saved scores require a terminal completion marker");
			p.set(line);
			require(count==expected && p.terms()==2 && p.termEquals("#COMPLETE", 0)
				&& p.parseLong(1)==count && input.reader.nextLine()==null, "Saved score completion count and terminal position must match the full accepted stream");
		}
		@Override public void close(){input.close();}
		final Input input;
		final int end;
		final LineParser1 p=new LineParser1('\t');
		final float[] scores=new float[22];
		long count;
		int selected0, selected01;
	}

	/** Primitive sufficient statistics, with full untruncated distances in bases. */
	private static final class Stats{
		void add(final long error, final long raw, final long comparator, final int bad){
			assert(bad>=0 && bad<4) : "Every final proposal has one geometry disposition";
			n++;
			final long absolute=Math.abs(error);
			if(absolute==0){exact++;}
			if(absolute<=3){within3++;}
			if(absolute<=10){within10++;}
			if(absolute<Math.abs(raw)){better++;}else if(absolute>Math.abs(raw)){worse++;}else{tied++;}
			if(absolute<Math.abs(comparator)){compareBetter++;}else if(absolute>Math.abs(comparator)){compareWorse++;}else{compareTied++;}
			squared+=(double)error*error;
			rawSquared+=(double)raw*raw;
			compareSquared+=(double)comparator*comparator;
			invalid[bad]++;
			assert(better+worse+tied==n && compareBetter+compareWorse+compareTied==n && exact<=within3 && within3<=within10 && within10<=n)
				: "Matched comparisons and nested tolerances must retain every original trial";
		}
		void row(final ByteBuilder b, final String policy, final int end, final String comparator,
				final long rawInvalid, final long scored, final long missing, final long held, final long unsupported){
			assert(n>0 && invalid[0]+invalid[1]+invalid[2]+invalid[3]==n) : "Final geometry counts must partition the complete assay";
			final double rms=Math.sqrt(squared/n), rawRms=Math.sqrt(rawSquared/n), compareRms=Math.sqrt(compareSquared/n);
			b.append(policy).tab().append(end==0 ? "5prime" : "3prime").tab().append(n).tab().append(exact).tab().append(within3).tab().append(within10)
				.tab().append(100.0*exact/n, 8).tab().append(100.0*within3/n, 8).tab().append(rms, 8).tab().append(rawRms, 8).tab().append(rms-rawRms, 8)
				.tab().append(comparator).tab().append(compareRms, 8).tab().append(rms-compareRms, 8).tab().append(better).tab().append(worse).tab().append(tied)
				.tab().append(compareBetter).tab().append(compareWorse).tab().append(compareTied).tab().append(n-invalid[0]).tab().append(invalid[1])
				.tab().append(invalid[2]).tab().append(invalid[3]).tab().append(rawInvalid).tab().append(scored).tab().append(missing).tab().append(held)
				.tab().append(unsupported).tab().append(centerOnly).tab().append(centerPlusNetwork).tab().append(reverted)
				.append("\tINDEPENDENT_ENDS_NO_POST_NETWORK_HOLD_CONTEXT_BOUNDS_NOT_CONTIG_BOUNDS").nl();
		}
		long n, exact, within3, within10, better, worse, tied, compareBetter, compareWorse, compareTied;
		long centerOnly, centerPlusNetwork, reverted;
		double squared, rawSquared, compareSquared;
		final long[] invalid=new long[4];
	}

	private String baseline, centered, ledger, loci, out;
	private int end=-1;
	private long expected=-1, usable=-1;
	private boolean geometryOnly;
	private static final String[] POLICIES={"raw", "uncentered_margin0", "uncentered_margin0.1", "center_only", "centered_margin0", "centered_margin0.1"};
	private static final String[] HOLDS={"", "CENTER_CROSSED_HOLD", "CENTER_OUTSIDE_CONTEXT_HOLD", "CENTER_BELOW_MINIMUM_HOLD"};
	private static final String[] INVALID={"VALID", "CROSSED", "OUTSIDE_CONTEXT", "BELOW_MINIMUM"};
}
