package prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ReadWrite;
import map.LongHashSet;
import map.ObjectIntMap;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.DoubleList;
import structures.FloatList;

/** Scores a move/hold gate once per proposed move and replays thresholds without
 * reevaluating either network. All original trials, including deterministic holds,
 * remain in every denominator. Invalid proposed spans are counted, not clipped;
 * these independent-end trials do not replace joint caller validation.
 * Optional numeric cells are empty: held trials have no gate score, controls
 * have no gate threshold, unchanged errors have no better fraction, and empty
 * move/classes have no percentages, AUC or quantiles. Existing policy/count
 * columns determine presence; computed nonfinite values always fail.
 * @author Brian, Ganyu, Raiden
 */
public final class NcrnaEndpointGateEvaluation{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointGateEvaluation.class, false);
		try{
			new NcrnaEndpointGateEvaluation(pp.args).process();
			pp.outstream.println("NCRNA_ENDPOINT_GATE_EVALUATION_COMPLETE_NO_FIT_NO_FIRST_STAGE_INFERENCE");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaEndpointGateEvaluation(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("replay")){replay=Parse.parseBoolean(b);}
			else if(a.equals("controls")){controls=Parse.parseBoolean(b);}
			else if(a.equals("pool")){pool=Parse.parseBoolean(b);}
			else if(a.equals("matchedrates")){matchedRates=Parse.parseBoolean(b);}
			else if(a.equals("oormargins")){
				final String[] values=b.split(",");
				require(values.length>0, "A declared out-of-range margin list must contain at least one comparison");
				oorMargins=new float[values.length];
				oorHeadUsed=new boolean[values.length];
				for(int i=0; i<values.length; i++){
					oorHeadUsed[i]=!values[i].equalsIgnoreCase("ignore");
					oorMargins[i]=oorHeadUsed[i] ? Float.parseFloat(values[i]) : 0;
					require(Float.isFinite(oorMargins[i]) && (oorHeadUsed[i]
						? i==0 || !oorHeadUsed[i-1] || oorMargins[i]>oorMargins[i-1] : i==0),
						"Out-of-range margins must be finite and ascend; optional ignore belongs first");
				}
			}
			else if(a.equals("rows")){rows=b;}
			else if(a.equals("scores")){scoresPath=b;}
			else if(a.equals("net")){network=b;}
			else if(a.equals("firstnet")){firstNet=b;}
			else if(a.equals("firstmatrix")){firstMatrix=b;}
			else if(a.equals("firstrows")){firstRows=b;}
			else if(a.equals("expected")){expected=Long.parseLong(b);}
			else if(a.equals("proposals")){expectedProposals=Long.parseLong(b);}
			else if(a.equals("minlen")){minimum=Integer.parseInt(b);}
			else if(a.equals("end")){
				require("5prime".equals(b) || "3prime".equals(b), "Specify end=5prime or3prime");
				end=b.equals("5prime") ? 0 : 1;
			}else if(a.equals("thresholds")){
				final String[] values=b.split(",");
				thresholds=new float[values.length];
				for(int i=0; i<values.length; i++){
					thresholds[i]=Float.parseFloat(values[i]);
					require(Float.isFinite(thresholds[i]) && (i==0 || thresholds[i]>thresholds[i-1]),
						"Thresholds must be finite and strictly increasing; scores are not probabilities");
				}
			}else{require(parser.parse(arg, a, b), "Unknown gate-evaluation argument: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		if((controls || pool) && thresholds==null){thresholds=new float[0];}
		require(input!=null && output!=null && end>=0 && expected>0 && expectedProposals>=0
			&& expectedProposals<=expected && minimum>0 && thresholds!=null && (controls || pool || thresholds.length>0),
			"Bind in= out=fresh_directory end= expected= proposals= minlen= thresholds=explicit_list");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && !new File(output).exists(),
			"Evaluation requires complete populations and fresh outputs without sampling");
		require(!(replay && controls), "Choose saved gate replay or first-stage controls, not both");
		require(!matchedRates || replay, "Matched movement rates require a completed saved gate cache; no additional inference");
		if(pool){
			require(!controls && !replay && rows==null && scoresPath==null && network==null && firstNet==null
				&& firstMatrix==null && firstRows==null && thresholds.length==0, "pool=t consumes only a five-fold control manifest");
		}else if(controls){
			require(network==null && scoresPath==null && rows!=null && firstNet!=null && firstMatrix!=null && firstRows!=null
				&& thresholds.length==0, "Controls use in=first_stage_scores rows=all_trial_audit and original bindings; no gate or cutoffs");
			scoresPath=input;
		}else if(replay){
			require(rows==null && scoresPath==null && network==null && firstNet==null && firstMatrix==null && firstRows==null,
				"Replay consumes only its complete saved gate cache, never vectors or networks");
		}else{
			require(rows!=null && scoresPath!=null && network!=null && firstNet!=null && firstMatrix!=null && firstRows!=null,
				"Inference requires rows= scores= net= firstnet= firstmatrix= firstrows= exact original bindings");
		}
		require(Tools.testInputFiles(false, true, input, rows, controls ? null : scoresPath, network), "Evaluation inputs must be readable");
		require(Tools.testForDuplicateFiles(true, input, rows, controls ? null : scoresPath, network), "Distinct vector/audit/score/network inputs are required");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){process(null);}

	/** Injection is only for finite fixtures; replay cannot invoke a network. */
	public void process(final CellNet fixture){
		require(!(replay || controls || pool) || fixture==null, "Replay and controls never evaluate a gate");
		final CellNet net=replay || controls || pool ? null : fixture==null ? CellNetParser.load(network, false) : fixture;
		require(net==null || net.numInputs()==32 && net.numOutputs()==1, "The move gate must have32inputs and one binary-target output");
		require(new File(output).mkdir(), "Cannot create fresh evaluation directory");
		final boolean oldBgzf=ReadWrite.USE_BGZF;
		final boolean oldForce=ReadWrite.FORCE_BGZIP;
		ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=true;
		try{
			if(pool){poolControls();}
			else if(replay){
				if(matchedRates){
					require(expectedProposals<=Integer.MAX_VALUE, "Exact score ranking requires a Java-indexable proposal population");
					ranking=new float[(int)expectedProposals];
					positiveScores=new FloatList();
					negativeScores=new FloatList();
				}
				readCache();
				if(matchedRates){
					prepareMatchedRates();
					matchingPass=true;
					readCache();
					for(int i=0; i<matchedCount; i++){
						require(matches[i].counts.moved==matches[i].attainable,
							"Ranked cutoff and replay disagree; the frozen gate cache may have changed between passes");
					}
				}
			}else{infer(net);}
		}finally{ReadWrite.USE_BGZF=oldBgzf;ReadWrite.FORCE_BGZIP=oldForce;}
		// A report exists only after all input and cached-output owners close.
		try(Output report=new Output(path("summary.tsv"))){
			final ByteBuilder text=new ByteBuilder(HEADER).nl();
			for(int i=0; i<results.length; i++){
				final String policy=i==0 ? "do_nothing" : i==1 ? "margin0" : i==2 ? "margin0.1"
					: i==3 ? "all_proposals" : (i-4)%2==0 ? "gate" : "gate_hard_oor_veto";
				results[i].row(text, policy, i<4 ? null : Float.toString(thresholds[(i-4)/2]), end, results[0].exact);
			}
			report.writer.print(text);
		}
		if(matchedRates){writeMatchedRates();writeRanking();}
		if(gaps!=null){writeGaps();}
		if(oorMargins!=null){writeOorSweep();}
	}

	private void initialize(){
		results=new Counts[4+2*thresholds.length];
		for(int i=0; i<results.length; i++){results[i]=new Counts();}
		gaps=controls || pool || oorMargins!=null ? new DoubleList() : null;
		if(oorMargins!=null){
			oorResults=new Counts[2*oorMargins.length];
			for(int i=0; i<oorResults.length; i++){oorResults[i]=new Counts();}
		}
		seen=new LongHashSet();
		genes=new ObjectIntMap<String>(String.class);
		lastGene=null;
		lastIndex=-1;
		total=proposals=0;
	}

	/** Pools primitive counts and squared errors, never rounded percentages/RMS.
	 * Each fold is read once and also gets its individual report. The manifest
	 * names previously qualified disjoint fold streams; no new split is inferred.
	 */
	private void poolControls(){
		initialize();
		final LineParser1 p=new LineParser1('\t');
		final String[] usedScores=new String[5];
		final String[] usedRows=new String[5];
		int count=0;
		try(Input manifest=new Input(input)){
			RrnaResourceIO.header(manifest.reader, "fold\tscores\trows\tfirstnet\tfirstmatrix\tfirstrows\ttrials\tproposals", input);
			for(byte[] line=manifest.reader.nextLine(); line!=null; line=manifest.reader.nextLine()){
				p.set(line);
				require(p.terms()==8 && count<5 && p.parseInt(0)==count, "Pool manifest must name folds0..4 exactly once in order");
				for(int i=0; i<8; i++){require(p.length(i)>0, "Empty pooled input binding");}
				usedScores[count]=p.parseString(1);
				usedRows[count]=p.parseString(2);
				for(int i=0; i<count; i++){
					require(!usedScores[i].equals(usedScores[count]) && !usedRows[i].equals(usedRows[count]), "A fold stream cannot be pooled twice");
				}
				final NcrnaEndpointGateEvaluation child=new NcrnaEndpointGateEvaluation(new String[]{
					"controls=t", "in="+usedScores[count], "rows="+usedRows[count], "firstnet="+p.parseString(3),
					"firstmatrix="+p.parseString(4), "firstrows="+p.parseString(5), "expected="+p.parseLong(6),
					"proposals="+p.parseLong(7), "end="+endName(), "minlen="+minimum, "out="+path("fold"+count)});
				child.oorMargins=oorMargins;
				child.oorHeadUsed=oorHeadUsed;
				child.process();
				require(!child.outer, "Only internal fold-held streams belong in the five-fold pool");
				for(int i=0; i<results.length; i++){results[i].add(child.results[i]);}
				gaps.addAll(child.gaps);
				if(oorMargins!=null){
					for(int i=0; i<oorResults.length; i++){oorResults[i].add(child.oorResults[i]);}
				}
				total+=child.total;
				proposals+=child.proposals;
				count++;
			}
		}
		require(count==5, "A pooled comparison requires all five completed folds");
		complete();
	}

	private void infer(final CellNet net){
		initialize();
		final LineParser1 audit=new LineParser1('\t');
		final LineParser1 p=new LineParser1('\t');
		final LineParser1 v=new LineParser1('\t');
		final float[] heads=new float[22];
		final float[] features=new float[32];
		final ByteBuilder text=new ByteBuilder();
		try(Input vectors=controls ? null : new Input(input); Input identities=new Input(rows); Input scores=new Input(scoresPath);
			Output cache=controls ? null : new Output(path("scores.tsv.gz"))){
			if(vectors!=null){RrnaResourceIO.header(vectors.reader, "#dims\t32\t1\t0", input);}
			RrnaResourceIO.header(identities.reader, AUDIT_HEADER, rows);
			final int offset=readFirstHeader(scores);
			if(cache!=null){cache.writer.print(new ByteBuilder(CACHE_MAGIC).nl().append("#coding\tFP32_RAW_SIGNED_INT\n")
				.append("#vectors\t").append(input).nl().append("#audit\t").append(rows).nl()
				.append("#first_scores\t").append(scoresPath).nl().append("#gate_network\t").append(network).nl()
				.append("#end\t").append(endName()).nl().append("#minlen\t").append(minimum).nl()
				.append("#trials\t").append(expected).nl().append("#proposals\t").append(expectedProposals).nl()
				.append(cacheHeader()).nl());}
			for(byte[] line=identities.reader.nextLine(); line!=null; line=identities.reader.nextLine()){
				audit.set(line);
				final byte[] raw=scores.reader.nextLine();
				require(raw!=null, "First-stage scores ended before the all-trial audit");
				p.set(raw);
				require(p.terms()==offset+22, "Unexpected first-stage score column count");
				join(audit, p, 11);
				readHeads(p, offset, heads);
				final Trial trial=validate(audit, heads);
				if(offset==16){
					final int label=Math.abs(trial.correction)>10 ? 21 : (int)trial.correction+10;
					require(p.parseInt(11)==label && p.parseLong(12)==trial.correction
						&& p.parseInt(13)==NcrnaEndpointConfidence.selectFull(heads, 10, 0)
						&& p.parseInt(14)==NcrnaEndpointConfidence.selectFull(heads, 10, .1f)
						&& p.parseInt(15)==Tools.maxIndex(heads), "Saved baseline decisions disagree with their raw scores");
				}
				float score=0; // Unused unless this non-control trial has a proposal.
				if(trial.proposal>=0 && !controls){
					final byte[] vector=vectors.reader.nextLine();
					require(vector!=null, "Proposal vector stream is truncated");
					v.set(vector);
					require(v.terms()==33, "Gate vectors require32features and one binary target");
					for(int i=0; i<32; i++){
						features[i]=number(v, i);
						require(Float.isFinite(features[i]), "Gate features must be finite");
						// GateVectors writes18decimal places: tiny scores can round to zero.
						if(i<22){require(features[i]==heads[i] || Math.abs((double)features[i]-heads[i])<=5.1e-19,
							"Gate features must preserve their paired first-stage scores");}
					}
					require(v.termEquals(trial.better ? "1" : "0", 32), "Gate target must be strict distance improvement");
					net.applyInput(features);
					net.feedForward();
					score=net.getOutput(0);
					require(Float.isFinite(score), "Gate inference produced a nonfinite decision score");
				}
				accept(trial, heads, score);
				if(cache!=null){
					text.clear().append(line);
					for(final float head:heads){text.tab().append(Float.floatToRawIntBits(head));}
					text.tab();
					if(trial.proposal>=0){text.append(Float.floatToRawIntBits(score));}
					cache.writer.print(text.tab().append(endName()).tab().append(Math.abs(trial.correction))
						.tab().append(Math.abs(trial.correction-(trial.proposal<0 ? 0 : trial.proposal-10))).nl());
				}
			}
			require(vectors==null || vectors.reader.nextLine()==null, "Proposal vectors extend beyond the all-trial audit");
			terminal(scores, p);
			complete();
			if(cache!=null){cache.writer.print(text.clear().append("#COMPLETE\t").append(total).tab().append(proposals).nl());}
		}
	}

	private int readFirstHeader(final Input scores){
		final LineParser1 p=new LineParser1('\t');
		final byte[] magic=scores.reader.nextLine();
		require(magic!=null, "Empty first-stage score stream");
		p.set(magic);
		final boolean full=p.terms()==1 && p.termEquals("#NCRNA_FULL_ERROR_SCORES_V1", 0);
		require(full || p.terms()==1 && p.termEquals("#NCRNA_GATE_FIRST_STAGE_V1", 0), "Unrecognized first-stage score schema");
		RrnaResourceIO.header(scores.reader, "#coding\tFP32_RAW_SIGNED_INT", scoresPath);
		RrnaResourceIO.header(scores.reader, "#matrix\t"+firstMatrix, scoresPath);
		RrnaResourceIO.header(scores.reader, "#rows\t"+firstRows, scoresPath);
		RrnaResourceIO.header(scores.reader, "#network\t"+firstNet, scoresPath);
		RrnaResourceIO.header(scores.reader, "#end\t"+endName(), scoresPath);
		byte[] header=scores.reader.nextLine();
		require(header!=null, "Missing first-stage column header");
		p.set(header);
		if(full && p.termEquals("#featurevariant", 0)){
			require(p.terms()==2 && p.termEquals("baseline", 1), "Outer comparison requires the frozen baseline feature layout");
			header=scores.reader.nextLine();
		}
		final ByteBuilder expectedHeader=new ByteBuilder(NcrnaPositionalFeatureVariants.ROW_HEADER);
		if(full){expectedHeader.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");}
		appendHeads(expectedHeader);
		require(header!=null && Arrays.equals(header, expectedHeader.toBytes()), "First-stage column order differs from the qualified schema");
		outer=full;
		return full ? 16 : 11;
	}

	private void readCache(){
		initialize();
		final LineParser1 p=new LineParser1('\t');
		final float[] heads=new float[22];
		try(Input cache=new Input(input)){
			final byte[] magic=cache.reader.nextLine();
			require(magic!=null, "Gate cache requires a version header");
			p.set(magic);
			final boolean legacy=p.terms()==1 && p.termEquals("#NCRNA_GATE_EVALUATION_V1", 0);
			require(legacy || p.terms()==1 && p.termEquals(CACHE_MAGIC, 0), "Unsupported gate-evaluation cache version");
			RrnaResourceIO.header(cache.reader, "#coding\tFP32_RAW_SIGNED_INT", input);
			for(final String field:new String[]{"#vectors", "#audit", "#first_scores", "#gate_network"}){
				final byte[] line=cache.reader.nextLine();
				require(line!=null, "Missing cache provenance");
				p.set(line);
				require(p.terms()==2 && p.termEquals(field, 0) && p.length(1)>0, "Incomplete cache provenance: "+field);
			}
			RrnaResourceIO.header(cache.reader, "#end\t"+endName(), input);
			RrnaResourceIO.header(cache.reader, "#minlen\t"+minimum, input);
			RrnaResourceIO.header(cache.reader, "#trials\t"+expected, input);
			RrnaResourceIO.header(cache.reader, "#proposals\t"+expectedProposals, input);
			RrnaResourceIO.header(cache.reader, cacheHeader(legacy).toString(), input);
			boolean closed=false;
			for(byte[] line=cache.reader.nextLine(); line!=null; line=cache.reader.nextLine()){
				p.set(line);
				if(p.termEquals("#COMPLETE", 0)){
					require(p.terms()==3 && p.parseLong(1)==total && p.parseLong(2)==proposals
						&& cache.reader.nextLine()==null, "Gate cache must close at its exact trial/proposal counts");
					closed=true;
					break;
				}
				require(p.terms()==44, "Gate cache needs18audit fields,22first-stage scores, gate score, end and two distances");
				if(total==0){outer=p.termEquals("evaluation", 3);}
				readHeads(p, 18, heads);
				final Trial trial=validate(p, heads);
				require(p.termEquals(endName(), 41) && p.parseLong(42)==Math.abs(trial.correction)
					&& p.parseLong(43)==Math.abs(trial.correction-(trial.proposal<0 ? 0 : trial.proposal-10)),
					"Cached end and before/proposed distances must agree with original coordinates");
				final float score;
				if(trial.proposal<0){
					require(legacy ? p.termEquals("NA", 40) : p.length(40)==0, "A bypassed trial cannot carry a gate inference");
					score=0; // Trial.proposal carries presence, never this unused finite slot.
				}else{
					require(p.length(40)>0 && !p.termEquals("NA", 40), "A proposed move requires its saved gate score");
					score=Float.intBitsToFloat(p.parseInt(40));
					require(Float.isFinite(score), "Saved gate scores must be finite");
					if(matchedRates && !matchingPass){
						ranking[(int)proposals-1]=score;
						if(trial.better){positiveScores.add(score);}else{negativeScores.add(score);}
					}
				}
				accept(trial, heads, score);
			}
			require(closed, "Truncated gate cache has no terminal marker");
			complete();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Stream Validation     ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates truth, deterministic proposal and every stored audit decision. */
	private Trial validate(final LineParser1 p, final float[] heads){
		require(p.terms()==18 || replay && p.terms()==44, "Unexpected all-trial audit schema");
		for(int i=0; i<18; i++){require(p.length(i)>0, "Empty audit field at trial "+(total+1));}
		require(p.termEquals(outer ? "evaluation" : "FOLD_HELD", 3), "A report cannot mix outer and fold-held populations");
		if(lastGene==null || !p.termEquals(lastGene, 0)){
			lastGene=p.parseString(0);
			lastIndex=genes.get(lastGene);
			if(lastIndex<0){lastIndex=genes.size();genes.put(lastGene, lastIndex);}
		}
		final int sample=p.parseInt(6);
		require(sample>=0 && p.parseInt(4)>=0 && seen.add(((long)lastIndex<<32)|sample), "Duplicate or invalid trial identity");
		final int start=p.parseInt(7);
		final int stop=p.parseInt(8);
		final int trueStart=p.parseInt(9);
		final int trueStop=p.parseInt(10);
		require(start>=0 && stop>=start && trueStart>=0 && trueStop>=trueStart, "Audit coordinates require valid raw and true spans");
		final Trial t=trial;
		t.correction=end==0 ? (long)trueStart-start : (long)trueStop-stop;
		t.proposal=NcrnaEndpointGateVectors.proposal(heads);
		final int shift=t.proposal<0 ? 0 : t.proposal-10;
		final long length=(long)stop-start+1+(end==0 ? -shift : shift);
		t.better=t.proposal>=0 && Math.abs(t.correction-shift)<Math.abs(t.correction);
		t.crossed=t.proposal>=0 && length<=0;
		t.shortSpan=t.proposal>=0 && length<minimum;
		require(p.parseInt(12)==t.proposal && p.parseLong(13)==t.correction
			&& p.termEquals(t.better ? "1" : "0", 14) && p.termEquals(t.crossed ? "1" : "0", 15)
			&& p.termEquals(t.shortSpan ? "1" : "0", 16)
			&& p.termEquals(t.proposal>=0 && heads[21]>heads[t.proposal] ? "1" : "0", 17),
			"Audit proposal/truth/geometry must agree with raw coordinates and saved first-stage scores");
		if(t.proposal>=0){proposals++;}
		require(p.parseLong(11)==(t.proposal>=0 ? proposals : 0), "Gate row indices must cover proposals exactly once in order");
		total++;
		require(total<=expected && proposals<=expectedProposals, "Input exceeds accepted producer counts");
		return t;
	}

	private void accept(final Trial t, final float[] heads, final float score){
		assert(t!=null && heads.length==22) : "Every policy must consume the same validated trial and original scores";
		results[0].accept(t, 0);
		results[1].accept(t, NcrnaEndpointConfidence.selectFull(heads, 10, 0)-10);
		results[2].accept(t, NcrnaEndpointConfidence.selectFull(heads, 10, .1f)-10);
		results[3].accept(t, t.proposal<0 ? 0 : t.proposal-10);
		if(gaps!=null && t.proposal>=0){gaps.add((double)heads[t.proposal]-heads[21]);}
		if(oorMargins!=null){
			for(int p=0; p<2; p++){
				for(int m=0; m<oorMargins.length; m++){
					final int selected=NcrnaEndpointConfidence.selectFull(heads, 10, p==0 ? 0 : .1f, oorMargins[m], oorHeadUsed[m]);
					oorResults[p*oorMargins.length+m].accept(t, selected-10);
				}
			}
		}
		for(int i=0; i<thresholds.length; i++){
			final int shift=t.proposal>=0 && score>=thresholds[i] ? t.proposal-10 : 0;
			results[4+2*i].accept(t, shift);
			results[5+2*i].accept(t, shift!=0 && heads[t.proposal]>heads[21] ? shift : 0);
		}
		if(matchingPass){
			for(int i=0; i<matchedCount; i++){
				final Matched m=matches[i];
				final boolean take=t.proposal>=0 && (m.strict ? score>m.cutoff : score>=m.cutoff);
				m.counts.accept(t, take ? t.proposal-10 : 0);
			}
		}
	}

	/** Uses only score ranks and the two control move counts, never truth labels.
	 * A tie spanning the requested count is reported on both sides; no arbitrary
	 * tie quota is silently presented as a deployable threshold.
	 */
	private void prepareMatchedRates(){
		Arrays.sort(ranking);
		matches=new Matched[5];
		matchedCount=0;
		addMatched("margin0", results[1].moved);
		addMatched("margin0.1", results[2].moved);
		final float allCutoff=ranking.length==0 ? 0 : Math.min(0, ranking[0]);
		matches[matchedCount++]=new Matched("all_proposals", proposals, proposals, allCutoff, false, 0, "EXACT");
	}

	private void addMatched(final String control, final long target){
		require(target>=0 && target<=ranking.length, "A control cannot accept more than the available proposals");
		if(target==0){
			final float cutoff=ranking.length==0 ? 0 : ranking[ranking.length-1];
			matches[matchedCount++]=new Matched(control, 0, 0, cutoff, true, 0, "EXACT");
			return;
		}
		final int boundary=ranking.length-(int)target;
		final float cutoff=ranking[boundary];
		int first=boundary;
		int after=boundary+1;
		while(first>0 && ranking[first-1]==cutoff){first--;}
		while(after<ranking.length && ranking[after]==cutoff){after++;}
		final long atLeast=ranking.length-first;
		final long above=ranking.length-after;
		final long tied=after-first;
		if(atLeast==target){
			matches[matchedCount++]=new Matched(control, target, target, cutoff, false, tied, "EXACT");
		}else{
			assert(above<target && atLeast>target) : "The ranked boundary tie brackets the unattainable exact move count";
			matches[matchedCount++]=new Matched(control, target, above, cutoff, true, tied, "TIE_LOWER");
			matches[matchedCount++]=new Matched(control, target, atLeast, cutoff, false, tied, "TIE_UPPER");
		}
	}

	private void writeMatchedRates(){
		try(Output report=new Output(path("matched_rates.tsv"))){
			final ByteBuilder text=new ByteBuilder("matched_to\ttarget_moves\ttarget_move_pct\tmatch_status\tcutoff_operator\ttied_at_cutoff\t")
				.append(HEADER).nl();
			for(int i=0; i<matchedCount; i++){
				final Matched m=matches[i];
				text.append(m.control).tab().append(m.target).tab().append(100.0*m.target/total, 6).tab().append(m.status).tab();
				if(ranking.length>0){text.append(m.strict ? ">" : ">=");}
				text.tab().append(m.tied).tab();
				m.counts.row(text, "gate_matched_rate", ranking.length>0 ? Float.toString(m.cutoff) : null, end, results[0].exact);
			}
			report.writer.print(text);
		}
	}

	/** Labels are used only to describe ranking quality, never to choose rate-matched cutoffs. */
	private void writeRanking(){
		require((long)positiveScores.size()+negativeScores.size()==proposals, "Every scored proposal has exactly one strict-improvement label");
		final float[] positive=positiveScores.toArray();
		final float[] negative=negativeScores.toArray();
		Arrays.sort(positive);
		Arrays.sort(negative);
		final boolean bothClasses=positive.length>0 && negative.length>0;
		try(Output quality=new Output(path("ranking.tsv")); Output distribution=new Output(path("gate_score_quantiles.tsv"))){
			final ByteBuilder text=new ByteBuilder();
			quality.writer.println("end\tproposals\timproving\tnon_improving\timproving_pct\tauc_higher_score_predicts_improvement\tqualification");
			text.append(endName()).tab().append(proposals).tab().append(positive.length).tab().append(negative.length).tab();
			if(proposals>0){text.append(100.0*positive.length/proposals, 6);}
			text.tab();
			if(bothClasses){text.append(auc(positive, negative), 10).append("\tTIES_COUNT_HALF_NO_NEW_INFERENCE");}
			else{text.append("\tSINGLE_CLASS_OR_NO_PROPOSALS");}
			quality.writer.print(text.nl());
			distribution.writer.println("end\tlabel\tproposals\tquantile\tgate_score\tdefinition");
			for(int label=0; label<2; label++){
				final float[] values=label==0 ? negative : positive;
				if(values.length==0){continue;}
				for(final double fraction:new double[]{0, .01, .05, .25, .5, .75, .95, .99, 1}){
					text.clear().append(endName()).tab().append(label==0 ? "non_improving" : "improving")
						.tab().append(values.length).tab().append(fraction, 6).tab();
					text.append(Float.toString(values[Math.max(0, (int)Math.ceil(fraction*values.length)-1)]));
					distribution.writer.print(text.append("\tPOSITIONAL_NEAREST_RANK\n"));
				}
			}
		}
	}

	/** Exact pair-count AUC for finite ascending arrays; equal scores earn half.
	 * Both classes must be nonempty; a missing class has no AUC measurement.
	 * Input arrays are borrowed and unchanged. */
	public static double auc(final float[] positive, final float[] negative){
		require(positive!=null && negative!=null && positive.length>0 && negative.length>0,
			"AUC requires a nonempty positive and negative score population");
		for(final float[] values:new float[][]{positive, negative}){
			for(int i=0; i<values.length; i++){
				require(Float.isFinite(values[i]) && (i==0 || values[i]>=values[i-1]), "AUC score arrays must be finite and ascending");
			}
		}
		int below=0;
		int above=0;
		long twiceWins=0;
		for(final float score:positive){
			while(below<negative.length && negative[below]<score){below++;}
			above=Math.max(above, below);
			while(above<negative.length && negative[above]==score){above++;}
			twiceWins+=2L*below+above-below;
		}
		final double result=twiceWins/(2.0*positive.length*negative.length);
		assert(result>=0 && result<=1) : "Concordant and half-tied positive/negative pairs cannot exceed all pairs";
		return result;
	}

	/** Positional nearest-rank quantiles; DoubleList.percentile is sum-weighted
	 * and would be incorrect for this signed score-difference distribution. */
	private void writeGaps(){
		require(gaps.size()==proposals, "The gap distribution must contain each unique movement proposal once");
		gaps.sort();
		try(Output report=new Output(path("gap_distribution.tsv"))){
			report.writer.println("end\tproposals\tquantile\tproposal_minus_oor_score\tdefinition");
			final ByteBuilder text=new ByteBuilder();
			if(proposals==0){return;}
			for(final double fraction:new double[]{0, .001, .01, .05, .1, .25, .5, .75, .9, .95, .99, .999, 1}){
				text.clear().append(endName()).tab().append(proposals).tab().append(fraction, 6).tab();
				final int index=Math.max(0, (int)Math.ceil(fraction*gaps.size())-1);
				text.append(Double.toString(gaps.get(index)));
				report.writer.print(text.append("\tPOSITIONAL_NEAREST_RANK_ALL_UNIQUE_PROPOSALS\n"));
			}
		}
	}

	private void writeOorSweep(){
		long bestNet=oorResults[0].better-oorResults[0].worse;
		double bestSquared=oorResults[0].squaredAfter;
		for(final Counts c:oorResults){
			require(c.trials==total, "Every signed-head-margin row must retain the complete original trial cohort");
			bestNet=Math.max(bestNet, c.better-c.worse);
			bestSquared=Math.min(bestSquared, c.squaredAfter);
		}
		try(Output report=new Output(path("oor_margin_sweep.tsv"))){
			final ByteBuilder text=new ByteBuilder("positional_margin\toor_margin\tbest_net_on_grid\tbest_rms_on_grid\t").append(HEADER).nl();
			for(int p=0; p<2; p++){
				for(int m=0; m<oorMargins.length; m++){
					final Counts c=oorResults[p*oorMargins.length+m];
					final String margin=oorHeadUsed[m] ? Float.toString(oorMargins[m]) : "ignore";
					text.append(p==0 ? "0" : "0.1").tab().append(margin)
						.tab().append(c.better-c.worse==bestNet ? 1 : 0).tab().append(c.squaredAfter==bestSquared ? 1 : 0).tab();
					c.row(text, "signed_oor_margin", margin, end, results[0].exact);
				}
			}
			report.writer.print(text);
		}
	}

	private void terminal(final Input stream, final LineParser1 p){
		final byte[] line=stream.reader.nextLine();
		require(line!=null, "First-stage cache lacks its completion marker");
		p.set(line);
		require(p.terms()==2 && p.termEquals("#COMPLETE", 0) && p.parseLong(1)==total && stream.reader.nextLine()==null,
			"First-stage scores must close at the complete all-trial count");
	}
	private void complete(){
		require(total==expected && proposals==expectedProposals, "Expected complete trial/proposal counts: "+expected+"/"+expectedProposals+" got "+total+"/"+proposals);
		for(final Counts c:results){require(c.trials==total, "Every policy must preserve the original denominator");}
	}
	private static void join(final LineParser1 a, final LineParser1 b, final int terms){
		require(a.terms()>=terms && b.terms()>=terms, "Short trial identity");
		for(int i=0; i<terms; i++){
			a.setBounds(i);
			b.setBounds(i);
			require(a.b()-a.a()==b.b()-b.a(), "Mismatched first-stage/audit identity width");
			for(int j=0; j<a.b()-a.a(); j++){require(a.line()[a.a()+j]==b.line()[b.a()+j], "Mismatched first-stage/audit identity");}
		}
	}
	private static void readHeads(final LineParser1 p, final int offset, final float[] heads){
		assert(heads.length==22) : "The frozen first stage has21positions and one non-positional head";
		for(int i=0; i<22; i++){
			require(p.length(offset+i)>0, "Empty first-stage score bits");
			heads[i]=Float.intBitsToFloat(p.parseInt(offset+i));
			require(Float.isFinite(heads[i]), "Nonfinite first-stage score");
		}
	}
	private static float number(final LineParser1 p, final int i){
		require(p.length(i)>0, "Empty gate feature");
		p.setBounds(i);
		for(int j=p.a(); j<p.b(); j++){
			if(p.line()[j]=='e' || p.line()[j]=='E'){return (float)Parse.parseDoubleSlow(p.line(), p.a(), p.b());}
		}
		return p.parseFloat(i);
	}
	private static void appendHeads(final ByteBuilder b){
		assert(b!=null) : "A raw-score schema must begin with its trial identity fields";
		for(int i=0; i<22; i++){b.append("\tscore_bits_").append(i);}
	}
	private static ByteBuilder cacheHeader(){return cacheHeader(false);}
	private static ByteBuilder cacheHeader(final boolean legacy){
		final ByteBuilder b=new ByteBuilder(AUDIT_HEADER);
		appendHeads(b);
		return b.append(legacy ? "\tgate_score_bits" : "\tgate_score_bits_if_proposal").append("\tend\tdistance_before\tdistance_proposed");
	}
	private String endName(){return end==0 ? "5prime" : "3prime";}
	private String path(final String name){return new File(output, name).getPath();}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------          Trial Counts        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Trial{
		long correction;
		int proposal;
		boolean better, crossed, shortSpan;
	}

	private static final class Matched{
		Matched(final String control_, final long target_, final long attainable_, final float cutoff_,
				final boolean strict_, final long tied_, final String status_){
			assert(control_!=null && status_!=null && target_>=0 && attainable_>=0 && Float.isFinite(cutoff_))
				: "Matched-rate rows bind a real control count and a finite, explicit score comparison";
			control=control_;
			target=target_;
			attainable=attainable_;
			cutoff=cutoff_;
			strict=strict_;
			tied=tied_;
			status=status_;
		}
		final String control, status;
		final long target, attainable, tied;
		final float cutoff;
		final boolean strict;
		final Counts counts=new Counts();
	}

	/** Primitive accounting over the complete independent-end trial population. */
	private static final class Counts{
		void add(final Counts c){
			assert(c!=null && c.trials>0) : "Pooling requires a complete nonempty fold, before rounding its derived metrics";
			trials+=c.trials;
			proposed+=c.proposed;
			moved+=c.moved;
			exact+=c.exact;
			within1+=c.within1;
			within3+=c.within3;
			better+=c.better;
			worse+=c.worse;
			tied+=c.tied;
			correct+=c.correct;
			harmed+=c.harmed;
			reachable+=c.reachable;
			fixed+=c.fixed;
			oor+=c.oor;
			oorBetter+=c.oorBetter;
			oorWorse+=c.oorWorse;
			oorTied+=c.oorTied;
			oorHeld+=c.oorHeld;
			crossed+=c.crossed;
			acceptedCrossed+=c.acceptedCrossed;
			shortSpan+=c.shortSpan;
			acceptedShort+=c.acceptedShort;
			missedGains+=c.missedGains;
			squaredBefore+=c.squaredBefore;
			squaredAfter+=c.squaredAfter;
		}
		void accept(final Trial t, final int shift){
			assert(shift>=-10 && shift<=10 && (shift==0 || t.proposal>=0 && shift==t.proposal-10))
				: "All policies either hold or accept the same unique reachable proposal";
			final long before=Math.abs(t.correction);
			final long after=Math.abs(t.correction-shift);
			trials++;
			if(t.proposal>=0){proposed++;}
			if(t.crossed){crossed++;}
			if(t.shortSpan){shortSpan++;}
			if(shift!=0){
				moved++;
				if(t.crossed){acceptedCrossed++;}
				if(t.shortSpan){acceptedShort++;}
			}
			if(t.better && shift==0){missedGains++;}
			if(after==0){exact++;}
			if(after<=1){within1++;}
			if(after<=3){within3++;}
			if(after<before){better++;}else if(after>before){worse++;}else{tied++;}
			if(before==0){correct++;if(shift!=0){harmed++;}}
			else if(before<=10){reachable++;if(after==0){fixed++;}}
			else{
				oor++;
				if(after<before){oorBetter++;}else if(after>before){oorWorse++;}else{oorTied++;}
				if(shift==0){oorHeld++;}
			}
			squaredBefore+=(double)before*before;
			squaredAfter+=(double)after*after;
			assert(better+worse+tied==trials && correct+reachable+oor==trials && oorBetter+oorWorse+oorTied==oor
				&& exact<=within1 && within1<=within3 && within3<=trials && acceptedCrossed<=crossed && acceptedShort<=shortSpan)
				: "Coordinate cohorts, tolerances and invalid-accept counts must conserve all trials";
		}
		void row(final ByteBuilder b, final String policy, final String threshold, final int end, final long baselineExact){
			assert(trials>0) : "Complete trial denominators are required for endpoint percentages and RMS";
			final double before=Math.sqrt(squaredBefore/trials);
			final double after=Math.sqrt(squaredAfter/trials);
			require(Double.isFinite(before) && Double.isFinite(after), "Endpoint RMS must be computed from finite squared-error totals");
			b.append(policy).tab();
			if(threshold!=null){b.append(threshold);}
			b.tab().append(end==0 ? "5prime" : "3prime")
				.tab().append(trials).tab().append(proposed).tab().append(moved).tab().append(exact).tab().append(within1).tab().append(within3)
				.tab().append(100.0*exact/trials, 6).tab().append(100.0*within3/trials, 6)
				.tab().append(before, 6).tab().append(after, 6).tab().append(after-before, 6)
				.tab().append(better).tab().append(worse).tab().append(tied).tab();
			if(better+worse>0){b.append(better/(double)(better+worse), 6);}
			b.tab().append(correct).tab().append(harmed).tab().append(reachable).tab().append(fixed)
				.tab().append(oor).tab().append(oorBetter).tab().append(oorWorse).tab().append(oorTied).tab().append(oorHeld)
				.tab().append(crossed).tab().append(acceptedCrossed).tab().append(shortSpan).tab().append(acceptedShort)
				.tab().append(exact-baselineExact).tab().append(100.0*(exact-baselineExact)/trials, 6)
				.tab().append(100.0*moved/trials, 6);
			percentage(b.tab(), better, moved);
			percentage(b.tab(), worse, moved);
			b.tab().append(moved-better-worse);
			percentage(b.tab(), moved-better-worse, moved);
			b.tab().append(missedGains).tab().append(100.0*(better-worse)/trials, 6).nl();
		}
		private static void percentage(final ByteBuilder b, final long n, final long d){
			assert(n>=0 && n<=d) : "Move outcome percentages use disjoint subsets of accepted moves";
			if(d>0){b.append(100.0*n/d, 6);}
		}
		long trials, proposed, moved, exact, within1, within3, better, worse, tied;
		long correct, harmed, reachable, fixed, oor, oorBetter, oorWorse, oorTied, oorHeld;
		long crossed, acceptedCrossed, shortSpan, acceptedShort, missedGains;
		double squaredBefore, squaredAfter;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields and Schemas     ----------------*/
	/*--------------------------------------------------------------*/

	public static final String AUDIT_HEADER=NcrnaPositionalFeatureVariants.ROW_HEADER+
		"\tgate_row1\tproposal_class\ttrue_minus_raw_bp\tstrictly_better\tcrossed\tbelow_minimum\toor_head_above_proposal";
	private static final String CACHE_MAGIC="#NCRNA_GATE_EVALUATION_V2";
	private static final String HEADER="policy\tthreshold\tend\ttrials\tproposals\tmoved\texact\twithin1\twithin3\texact_pct\twithin3_pct"
		+"\trms_before\trms_after\trms_delta\tbetter\tworse\ttied\tbetter_fraction\tcorrect_before\tcorrect_harmed\treachable_wrong\treachable_fixed"
		+"\toor\toor_improved\toor_worsened\toor_tied\toor_held\tcrossed_proposals\taccepted_crossed\tshort_proposals\taccepted_short\texact_delta\texact_delta_pp"
		+"\tmoved_pct_trials\tbetter_pct_moves\tworse_pct_moves\ttied_moves\ttied_pct_moves\tmissed_gains\tnet_better_minus_worse_per100_trials";
	private String input, output, rows, scoresPath, network, firstNet, firstMatrix, firstRows;
	private int end=-1, minimum=-1;
	private long expected=-1, expectedProposals=-1, total, proposals;
	private boolean replay, outer, controls, pool, matchedRates, matchingPass;
	private float[] ranking;
	private FloatList positiveScores, negativeScores;
	private float[] oorMargins;
	/** Ignore is a comparison mode, never an infinite signed numeric margin. */
	private boolean[] oorHeadUsed;
	private Counts[] oorResults;
	private DoubleList gaps;
	private Matched[] matches;
	private int matchedCount;
	private float[] thresholds;
	private Counts[] results;
	private LongHashSet seen;
	private ObjectIntMap<String> genes;
	private String lastGene;
	private int lastIndex;
	private final Trial trial=new Trial();
}
