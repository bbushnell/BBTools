package test.prok;

import java.io.File;
import java.util.Locale;

import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEndpointConfidence;
import prok.NcrnaEndpointGateEvaluation;
import prok.NcrnaEndpointGateVectors;
import prok.NcrnaPositionalFeatureVariants;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Literal independent-end outcomes plus complete inference/replay joins.
 * These tiny fixtures qualify bookkeeping, not biological performance.
 * @author Brian, Ganyu, Raiden
 */
public final class NcrnaEndpointGateEvaluationTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointGateEvaluationTest.class, false);
		try{
			final Parser p=new Parser();
			p.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(p.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown fixture argument");
			}
			check(p.out1!=null && p.in1==null && p.in2==null && p.out2==null && !p.overwrite
				&& !p.append && p.maxReads<0 && p.samplerate==1, "Use out=fresh_directory and native runtime options only");
			final File root=new File(p.out1);
			check(root.mkdir(), "Fresh fixture root required");
			check(NcrnaEndpointGateEvaluation.auc(new float[]{2}, new float[]{1})==1, "Perfect higher-positive ranking has AUC1");
			check(NcrnaEndpointGateEvaluation.auc(new float[]{0}, new float[]{1})==0, "Reversed ranking has AUC0");
			check(NcrnaEndpointGateEvaluation.auc(new float[]{1, 1}, new float[]{1, 1})==.5, "Equal scores earn half credit");
			check(Math.abs(NcrnaEndpointGateEvaluation.auc(new float[]{1, 2, 2}, new float[]{0, 2, 3})-4.0/9)<1e-12,
				"Repeated positive and negative ties contribute literal pair counts");
			fails(()->NcrnaEndpointGateEvaluation.auc(new float[0], new float[]{1}));
			fails(()->NcrnaEndpointGateEvaluation.auc(new float[]{1}, new float[0]));
			fails(()->NcrnaEndpointGateEvaluation.auc(new float[]{2, 1}, new float[]{0}));
			for(int end=0; end<2; end++){
				for(final boolean full:new boolean[]{false, true}){run(new File(root, end+"_"+full), end, full);}
			}
			pp.outstream.println("NCRNA_ENDPOINT_GATE_EVALUATION_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Finite Fixtures       ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final File root, final int end, final boolean full){
		check(root.mkdir(), "Fixture case must be fresh");
		final String net=write(root, "network.txt", "Finite injected network; never loaded.\n");
		final ByteBuilder audit=new ByteBuilder(NcrnaEndpointGateEvaluation.AUDIT_HEADER).nl();
		final ByteBuilder vectors=new ByteBuilder("#dims\t32\t1\t0\n");
		final ByteBuilder scores=new ByteBuilder(full ? "#NCRNA_FULL_ERROR_SCORES_V1\n" : "#NCRNA_GATE_FIRST_STAGE_V1\n")
			.append("#coding\tFP32_RAW_SIGNED_INT\n#matrix\toriginal_matrix\n#rows\toriginal_rows\n#network\toriginal_net\n")
			.append("#end\t").append(end==0 ? "5prime" : "3prime").nl();
		if(full){scores.append("#featurevariant\tbaseline\n");}
		scores.append(NcrnaPositionalFeatureVariants.ROW_HEADER);
		if(full){scores.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");}
		for(int i=0; i<22; i++){scores.append("\tscore_bits_").append(i);}
		scores.nl();
		final int[] corrections={2, 0, 20, -20, 1, 0, 2, 10, 1};
		final int[] shifts={2, 2, 5, 5, 1, 0, 2, 10, 2};
		int proposed=0;
		for(int row=0; row<9; row++){
			final int sign=end==0 ? 1 : -1;
			final int correction=sign*corrections[row];
			final float[] heads=new float[22];
			heads[10+sign*shifts[row]]=1;
			if(row==0){heads[21]=2;}
			if(row==4){heads[10+2*sign]=1;}
			if(row==6){heads[21]=1;}
			final int choice=NcrnaEndpointGateVectors.proposal(heads);
			final int stop=row==7 ? 104 : 199;
			final ByteBuilder identity=new ByteBuilder("gene").append(row).append("\tsource\tassembly\t")
				.append(full ? "evaluation" : "FOLD_HELD").append("\t0\tmodel\t0\t100\t").append(stop)
				.tab().append(end==0 ? 100+correction : 0).tab().append(end==0 ? 300 : stop+correction);
			final boolean better=choice>=0 && Math.abs(correction-(choice-10))<Math.abs(correction);
			if(choice>=0){proposed++;}
			audit.append(identity).tab().append(choice>=0 ? proposed : 0).tab().append(choice).tab().append(correction)
				.tab().append(better ? 1 : 0).tab().append(row==7 ? 1 : 0).tab().append(row==7 ? 1 : 0)
				.tab().append(row==0 ? 1 : 0).nl();
			scores.append(identity);
			if(full){
				scores.tab().append(Math.abs(correction)>10 ? 21 : correction+10).tab().append(correction)
					.tab().append(NcrnaEndpointConfidence.selectFull(heads, 10, 0))
					.tab().append(NcrnaEndpointConfidence.selectFull(heads, 10, .1f)).tab().append(Tools.maxIndex(heads));
			}
			for(final float head:heads){scores.tab().append(Float.floatToRawIntBits(head));}
			scores.nl();
			if(choice>=0){
				for(final float head:heads){vectors.append(head, 18).tab();}
				vectors.append("1\t1\t0.5\t0.5\t0.5\t0\t0.7\t1e-4\t0\t0\t").append(better ? 1 : 0).nl();
			}
		}
		check(proposed==7, "Literal fixture has seven proposals and two bypasses, including one distance-tied move");
		scores.append("#COMPLETE\t9\n");
		final String vectorPath=write(root, "vectors.tsv", vectors.toString());
		final String auditPath=write(root, "rows.tsv", audit.toString());
		final String scoresPath=write(root, "first.tsv", scores.toString());
		final String output=path(root, "inference");
		final String[] args={"in="+vectorPath, "rows="+auditPath, "scores="+scoresPath, "net="+net,
			"firstnet=original_net", "firstmatrix=original_matrix", "firstrows=original_rows",
			"end="+(end==0 ? "5prime" : "3prime"), "expected=9", "proposals=7", "minlen=90",
			"thresholds=0.5,0.75,0.8", "out="+output};
		final FixedNet fixed=new FixedNet(.75f);
		new NcrnaEndpointGateEvaluation(args).process(fixed);
		check(fixed.calls==7, "One gate inference per proposal, none for the two deterministic bypasses");
		verifySummary(output);
		final String[] replay={"replay=t", "in="+output+"/scores.tsv.gz", "out="+path(root, "replay"),
			args[7], args[8], args[9], args[10], args[11], "matchedrates=t"};
		new NcrnaEndpointGateEvaluation(replay).process();
		check(read(output+"/summary.tsv").equals(read(path(root, "replay/summary.tsv"))),
			"Saved raw gate-score replay must reproduce every policy byte without any network");
		checkMatched(path(root, "replay/matched_rates.tsv"), true);
		checkRanking(path(root, "replay/ranking.tsv"), .5);
		final String cache=read(output+"/scores.tsv.gz");
		final String[] legacyReplay=replay.clone();
		legacyReplay[1]="in="+write(root, "legacy_cache.tsv", legacy(cache));
		legacyReplay[2]="out="+path(root, "legacy_replay");
		new NcrnaEndpointGateEvaluation(legacyReplay).process();
		for(final String suffix:new String[]{"summary", "matched_rates", "ranking", "gate_score_quantiles"}){
			check(read(path(root, "replay/"+suffix+".tsv")).equals(read(path(root, "legacy_replay/"+suffix+".tsv"))),
				"Complete finite results must survive old-cache normalization: "+suffix);
		}
		if(end==0 && !full){emptyClasses(root, cache);}
		final String[] varying=args.clone();
		varying[12]="out="+path(root, "varying");
		new NcrnaEndpointGateEvaluation(varying).process(new FixedNet(.1f, true));
		final String[] varyingReplay=replay.clone();
		varyingReplay[1]="in="+path(root, "varying/scores.tsv.gz");
		varyingReplay[2]="out="+path(root, "varying_replay");
		new NcrnaEndpointGateEvaluation(varyingReplay).process();
		checkMatched(path(root, "varying_replay/matched_rates.tsv"), false);
		checkRanking(path(root, "varying_replay/ranking.tsv"), 5.0/12);
		final String[] control={"controls=t", "in="+scoresPath, "rows="+auditPath, args[4], args[5], args[6],
			args[7], args[8], args[9], args[10], "out="+path(root, "controls"), "oormargins=ignore,-0.1,0,0.1"};
		new NcrnaEndpointGateEvaluation(control).process();
		final String[] emptyMargins=java.util.Arrays.copyOf(control, control.length+1);
		emptyMargins[10]="out="+path(root, "empty_margins");
		emptyMargins[control.length]="oormargins=,";
		fails(()->new NcrnaEndpointGateEvaluation(emptyMargins));
		check(!new File(path(root, "empty_margins")).exists(), "An empty declared comparison grid must fail before output creation");
		try(Input c=new Input(path(root, "controls/summary.tsv")); Input all=new Input(output+"/summary.tsv")){
			for(int i=0; i<5; i++){
				check(java.util.Arrays.equals(c.reader.nextLine(), all.reader.nextLine()), "Saved first-stage controls equal the same all-trial comparison after gate inference");
			}
			check(c.reader.nextLine()==null, "Controls emit exactly unchanged, two margins and all proposals");
		}
		checkOor(path(root, "controls"));
		if(!full){
			final ByteBuilder manifest=new ByteBuilder("fold\tscores\trows\tfirstnet\tfirstmatrix\tfirstrows\ttrials\tproposals\n");
			for(int fold=0; fold<5; fold++){
				final String foldScores=write(root, "fold"+fold+".scores.tsv", scores.toString().replace("gene", "fold"+fold+"_gene"));
				final String foldRows=write(root, "fold"+fold+".rows.tsv", audit.toString().replace("gene", "fold"+fold+"_gene"));
				manifest.append(fold).tab().append(foldScores).tab().append(foldRows)
					.append("\toriginal_net\toriginal_matrix\toriginal_rows\t9\t7\n");
			}
			final String manifestPath=write(root, "folds.tsv", manifest.toString());
			new NcrnaEndpointGateEvaluation(new String[]{"pool=t", "in="+manifestPath, "out="+path(root, "pool"),
				args[7], "expected=45", "proposals=35", args[10]}).process();
			try(Input pooled=new Input(path(root, "pool/summary.tsv"))){
				pooled.reader.nextLine();
				final LineParser1 p=new LineParser1('\t').set(pooled.reader.nextLine());
				check(p.parseLong(3)==45 && p.parseLong(4)==35 && p.parseLong(6)==10 && p.parseLong(38)==20,
					"Five-fold aggregation sums exact counts including missed gains without dropping deterministic holds");
			}
		}

		final String truncated=write(root, "truncated.tsv", scores.toString().replace("#COMPLETE\t9\n", ""));
		fails(args, 2, "scores="+truncated, root, "truncated_first", new FixedNet(.75f));
		final String wrong=write(root, "wrong.tsv", scores.toString().replace("gene0\t", "other0\t"));
		fails(args, 2, "scores="+wrong, root, "wrong_identity", new FixedNet(.75f));
		final String badTarget=write(root, "bad_target.tsv", vectors.toString().replaceFirst("1e-4\t0\t0\t1", "1e-4\t0\t0\t0"));
		fails(args, 0, "in="+badTarget, root, "wrong_target", new FixedNet(.75f));
		fails(args, 5, "firstmatrix=wrong_matrix", root, "wrong_provenance", new FixedNet(.75f));
		fails(args, 9, "proposals=8", root, "wrong_count", new FixedNet(.75f));
		fails(args, 11, "thresholds=0.8,0.5", root, "unsorted_thresholds", new FixedNet(.75f));
		fails(args, 11, "thresholds=NaN", root, "nonfinite_threshold", new FixedNet(.75f));
		fails(args, 8, "expected=9", root, "nonfinite_net", new FixedNet(Float.NaN));
		final String missingTerminal=write(root, "cache_truncated.tsv", cache.replace("#COMPLETE\t9\t7\n", ""));
		final String[] badReplay=replay.clone();
		badReplay[1]="in="+missingTerminal;
		badReplay[2]="out="+path(root, "bad_replay");
		fails(()->new NcrnaEndpointGateEvaluation(badReplay).process());
		check(!new File(path(root, "bad_replay/summary.tsv")).exists(), "Truncated caches must not publish a summary");
	}

	private static void verifySummary(final String output){
		try(Input in=new Input(output+"/summary.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			int rows=0;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==40 && p.parseLong(3)==9 && p.parseLong(4)==7, "Every policy retains nine trials and seven proposals");
				if(p.termEquals("do_nothing", 0)){
					check(p.parseLong(6)==2 && p.parseLong(16)==9 && p.length(17)==0 && p.length(1)==0 && p.parseLong(38)==4,
						"Unchanged baseline, zero-movement ratio and four missed improving proposals");
				}else if(p.termEquals("gate", 0) && p.parseFloat(1)<=.75f){
					check(p.parseLong(5)==7 && p.parseLong(6)==4 && p.parseLong(14)==4 && p.parseLong(15)==2
						&& p.parseLong(16)==3 && p.parseLong(19)==1 && p.parseLong(21)==3
						&& p.parseLong(23)==1 && p.parseLong(24)==1 && p.parseLong(28)==1 && p.parseLong(30)==1
						&& Math.abs(p.parseFloat(33)-700f/9)<.0001 && p.parseLong(36)==1 && p.parseLong(38)==0,
						"Gate fixes three reachable ends, partly improves one distant end, harms two and retains invalid acceptance counts");
					check(Math.abs(p.parseFloat(39)-200f/9)<.0001, "Net movement benefit is better minus worse per100original trials");
				}else if(p.termEquals("gate_hard_oor_veto", 0) && p.parseFloat(1)<=.75f){
					check(p.parseLong(5)==5 && p.parseLong(6)==2, "Hard head veto holds on both a winning and tied non-positional head");
				}else if(p.termStartsWith("gate", 0)){
					check(p.parseLong(5)==0 && p.parseLong(6)==2, "Threshold above gate score reproduces do-nothing outcomes");
				}
				rows++;
			}
			check(rows==10, "Four controls plus two policies for each of three thresholds");
		}
	}

	/** One literal hold, improving move or harmful move isolates missing-class state. */
	private static void emptyClasses(final File root, final String cache){
		for(int kind=0; kind<3; kind++){
			final boolean proposal=kind>0, improving=kind==1;
			final int correction=improving ? 2 : 0;
			final String prefix=cache.substring(0, cache.indexOf("\ngene0\t")+1)
				.replace("#trials\t9\n", "#trials\t1\n").replace("#proposals\t7\n", "#proposals\t"+(proposal ? 1 : 0)+"\n");
			final ByteBuilder b=new ByteBuilder(prefix).append("one\tsource\tassembly\tFOLD_HELD\t0\tmodel\t0\t100\t199\t")
				.append(100+correction).append("\t199\t").append(proposal ? 1 : 0).tab().append(proposal ? 12 : -1)
				.tab().append(correction).tab().append(improving ? 1 : 0).append("\t0\t0\t0");
			for(int i=0; i<22; i++){b.tab().append(Float.floatToIntBits(i==(proposal ? 12 : 10) ? 1 : 0));}
			b.tab();
			if(proposal){b.append(Float.floatToIntBits(0));}
			b.append("\t5prime\t").append(correction).tab().append(Math.abs(correction-(proposal ? 2 : 0)))
				.append("\n#COMPLETE\t1\t").append(proposal ? 1 : 0).nl();
			final String current=b.toString();
			final String output=path(root, "one_"+kind);
			final String[] args={"replay=t", "in="+write(root, "one_"+kind+".tsv", current), "out="+output,
				"end=5prime", "expected=1", "proposals="+(proposal ? 1 : 0), "minlen=90", "thresholds=0", "matchedrates=t", "oormargins=ignore,0"};
			new NcrnaEndpointGateEvaluation(args).process();
			try(Input input=new Input(output+"/ranking.tsv")){
				input.reader.nextLine();
				final LineParser1 p=new LineParser1('\t').set(input.reader.nextLine());
				check(p.parseLong(1)==(proposal ? 1 : 0) && p.parseLong(2)==(improving ? 1 : 0)
					&& p.parseLong(3)==(kind==2 ? 1 : 0) && p.length(5)==0 && (proposal || p.length(4)==0),
					"An absent class has no AUC; absence of all proposals also removes its improving percentage");
			}
			check(read(output+"/gate_score_quantiles.tsv").split("\n").length==(proposal ? 10 : 1),
				"Only a populated score class has nine measured quantiles");
			check(read(output+"/gap_distribution.tsv").split("\n").length==(proposal ? 14 : 1),
				"No proposal population has no measured gap quantiles");
			for(final String name:new String[]{"summary", "matched_rates", "ranking", "gate_score_quantiles", "gap_distribution", "oor_margin_sweep"}){
				for(final String cell:read(output+"/"+name+".tsv").split("[\t\n]")){
					check(!cell.equals("NA") && !cell.equals("NaN") && !cell.contains("Infinity"), "Missing statistics are omitted, computed values stay finite");
				}
			}
			if(!proposal){
				try(Input input=new Input(output+"/matched_rates.tsv")){
					input.reader.nextLine();
					final LineParser1 p=new LineParser1('\t').set(input.reader.nextLine());
					check(p.length(4)==0 && p.length(7)==0 && p.length(23)==0 && p.length(40)==0,
						"No observed scores supply no matched cutoff/operator, unchanged fraction or accepted-move percentage");
				}
				final String[] bad=args.clone();
				bad[1]="in="+write(root, "bad_current_hold.tsv", current.replace("\t\t5prime\t", "\tNA\t5prime\t"));
				bad[2]="out="+path(root, "bad_current_hold");
				fails(()->new NcrnaEndpointGateEvaluation(bad).process());
				final String[] old=args.clone();
				old[1]="in="+write(root, "bad_legacy_hold.tsv", legacy(current).replace("\tNA\t5prime\t", "\t\t5prime\t"));
				old[2]="out="+path(root, "bad_legacy_hold");
				fails(()->new NcrnaEndpointGateEvaluation(old).process());
			}else{
				try(Input input=new Input(output+"/summary.tsv")){
					for(int i=0; i<5; i++){input.reader.nextLine();}
					final LineParser1 p=new LineParser1('\t').set(input.reader.nextLine());
					check(p.termEquals("gate", 0) && p.parseFloat(1)==0 && p.parseLong(5)==1,
						"A measured zero gate score is accepted at its real inclusive zero cutoff");
				}
				final String[] missing=args.clone();
				missing[1]="in="+write(root, "bad_missing_"+kind+".tsv", current.replace("\t0\t5prime\t", "\t\t5prime\t"));
				missing[2]="out="+path(root, "bad_missing_"+kind);
				fails(()->new NcrnaEndpointGateEvaluation(missing).process());
				for(final float invalid:new float[]{Float.NaN, Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY}){
					final String suffix="bad_"+kind+"_"+Float.floatToIntBits(invalid);
					final String[] bad=args.clone();
					bad[1]="in="+write(root, suffix+".tsv", current.replace("\t0\t5prime\t", "\t"+Float.floatToIntBits(invalid)+"\t5prime\t"));
					bad[2]="out="+path(root, suffix);
					fails(()->new NcrnaEndpointGateEvaluation(bad).process());
				}
			}
		}
	}

	/** The old version alone declares NA for a bypassed gate, never for inference. */
	private static String legacy(final String current){
		final ByteBuilder b=new ByteBuilder();
		for(final String line:current.split("\n")){
			if(line.equals("#NCRNA_GATE_EVALUATION_V2")){b.append("#NCRNA_GATE_EVALUATION_V1\n");}
			else if(line.contains("\tgate_score_bits_if_proposal\t")){b.append(line.replace("\tgate_score_bits_if_proposal\t", "\tgate_score_bits\t")).nl();}
			else if(!line.startsWith("#")){
				final String[] fields=line.split("\t", -1);
				check(fields.length==44, "Legacy conversion must retain every cache field");
				if(fields[40].isEmpty()){fields[40]="NA";}
				for(int i=0; i<fields.length; i++){
					if(i>0){b.tab();}
					b.append(fields[i]);
				}
				b.nl();
			}else{b.append(line).nl();}
		}
		return b.toString();
	}

	private static void checkRanking(final String path, final double expected){
		try(Input input=new Input(path)){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t').set(input.reader.nextLine());
			check(p.terms()==7 && p.parseLong(1)==7 && p.parseLong(2)==4 && p.parseLong(3)==3
				&& Math.abs(p.parseDouble(5)-expected)<1e-9 && input.reader.nextLine()==null,
				"Saved-score ranking uses every proposed row and its literal strict-improvement label");
		}
	}

	private static void checkMatched(final String path, final boolean tied){
		try(Input in=new Input(path)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			int rows=0;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==46 && p.parseLong(9)==9, "Matched rows preserve the complete cohort and all metric columns");
				if(p.termEquals("all_proposals", 0)){
					check(p.parseLong(1)==7 && p.parseLong(11)==7 && p.termEquals("EXACT", 3) && p.parseFloat(7)==0,
						"Nonnegative gate outputs accept every proposal at cutoff0");
				}else{
					check(p.parseLong(1)==5, "Each matched rule targets its own cohort's exact control move count");
					if(tied){
						check(p.parseLong(5)==7 && p.parseFloat(7)==.75f, "All seven equal scores form one explicit cutoff tie");
						if(p.termEquals("TIE_LOWER", 3)){
							check(p.termEquals(">", 4) && p.parseLong(11)==0, "Strict cutoff reports the lower attainable count, without hidden tie quota");
						}else{
							check(p.termEquals("TIE_UPPER", 3) && p.termEquals(">=", 4) && p.parseLong(11)==7,
								"Inclusive cutoff reports the upper attainable count");
						}
					}else{
						check(p.termEquals("EXACT", 3) && p.parseLong(11)==5 && p.parseLong(5)==1,
							"Unique ranked scores match the control count exactly without consulting labels");
					}
				}
				rows++;
			}
			check(rows==(tied ? 5 : 3), "Each matched rule gets one exact row or both attainable tie bounds, plus all proposals");
		}
	}

	/** Original policies occupy exact cells of the decomposed signed-margin grid. */
	private static void checkOor(final String directory){
		final String[] controlRows=new String[4];
		try(Input in=new Input(directory+"/summary.tsv")){
			in.reader.nextLine();
			for(int i=0; i<4; i++){controlRows[i]=new String(in.reader.nextLine(), java.nio.charset.StandardCharsets.UTF_8);}
		}
		int checked=0;
		try(Input in=new Input(directory+"/oor_margin_sweep.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			final LineParser1 control=new LineParser1('\t');
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==44, "Signed-margin report carries both independent margins and all metrics");
				int index=-1;
				if(p.termEquals("0", 0) && p.termEquals("0.0", 1)){index=1;}
				if(p.termEquals("0.1", 0) && p.termEquals("0.1", 1)){index=2;}
				if(p.termEquals("0", 0) && p.termEquals("ignore", 1)){index=3;}
				if(index>=0){
					control.set(controlRows[index].getBytes(java.nio.charset.StandardCharsets.UTF_8));
					for(int i=3; i<40; i++){
						check(p.termEquals(control.parseString(i), i+4), "Old policy counts/metrics must reproduce exactly at their signed-margin grid cell");
					}
					checked++;
				}
				if(p.termEquals("0", 0) && p.termEquals("-0.1", 1)){
					check(p.parseLong(9)==6, "A negative head margin admits the head-tied proposal but not the worse-than-head proposal");
				}
			}
		}
		check(checked==3, "Both old strict margins and all proposals have exact grid equivalents");
		try(Input in=new Input(directory+"/gap_distribution.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t').set(in.reader.nextLine());
			check(p.parseLong(1)==7 && p.parseDouble(2)==0 && p.parseDouble(3)==-1,
				"Observed signed-gap distribution counts seven unique proposals and retains its negative minimum");
		}
	}
	/*--------------------------------------------------------------*/
	/*----------------        Fixture Helpers       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class FixedNet extends CellNet{
		FixedNet(final float value_){this(value_, false);}
		FixedNet(final float value_, final boolean varying_){
			value=value_;
			varying=varying_;
		}
		@Override public int numInputs(){return 32;}
		@Override public int numOutputs(){return 1;}
		@Override public CellNet applyInput(final float[] values){
			NcrnaEndpointGateEvaluationTest.check(values.length==32 && values[29]==.0001f, "Complete gate feature layout accepts exponent values");
			return this;
		}
		@Override public float feedForward(){calls++;return value;}
		@Override public float getOutput(final int index){return value+(varying ? calls*.1f : 0);}
		final float value;
		final boolean varying;
		int calls;
	}
	private static void fails(final String[] original, final int index, final String replacement,
			final File root, final String label, final CellNet net){
		final String[] args=original.clone();
		args[index]=replacement;
		args[12]="out="+path(root, label);
		fails(()->new NcrnaEndpointGateEvaluation(args).process(net));
		check(!new File(path(root, label+"/summary.tsv")).exists(), "Failed input cannot publish a complete report");
	}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed gate evaluation input must fail loudly");
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static String write(final File root, final String name, final String value){
		final String path=path(root, name);
		try(Output out=new Output(path)){out.writer.print(value);}
		return path;
	}
	private static String read(final String path){
		final ByteBuilder text=new ByteBuilder();
		try(Input in=new Input(path)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){text.append(line).nl();}
		}
		return text.toString();
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
