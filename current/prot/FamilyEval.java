package prot;

import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import ddl.DDLFamilySearch.StructuralInputs;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prot.FamilyPanelRunner.Panel;
import prot.FamilyPanelRunner.PanelCounts;
import shared.Resources;
import shared.Shared;
import structures.ByteBuilder;

/** Evaluates labeled protein panels through the existing deployed classifier.
 * One loaded index and worker pool serve every panel; counts do not need replay TSVs.
 * Optional benchmarkrepeats reuses the classifier and reports repeated work separately from scientific counts.
 * @author Brian, Keqing, Fischl
 */
public final class FamilyEval{
	private String requestedMatrix;

	public static void main(String[] args) throws Exception{
		try{
			final FamilyEval tool=new FamilyEval(args);
			try{tool.process();}finally{Shared.closeStream(tool.outstream);}
		}catch(Resources.MissingResourceException e){
			System.err.print(e.getMessage());
			System.exit(1);
		}
	}

	public FamilyEval(String[] args){
		final PreParser pp=new PreParser(args,getClass(),false);
		outstream=pp.outstream;
		ReadWrite.ZIPLEVEL=9; ReadWrite.FORCE_BGZIP=true; ReadWrite.USE_BGZIP=true;
		final Parser parser=new Parser();
		String positive=null,negatives=null,names=null,positivePin=null,negativePins=null;
		String verification=null,positiveCount=null,negativeCounts=null,evaluationGc=null;boolean metadataSpecified=false;
		for(String arg:pp.args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT),value=pair.length>1 ? pair[1] : null;
			if(Arrays.asList("lengthmodel","lengthmodelsha80","lengthmedians","lengthmedianssha80",
				"familyhydro","familyhydrosha80","pairedscales").contains(key)){
				metadataSpecified=true;
			}
			if(key.equals("verification")){
				require(verification==null && ("sha80".equals(value) || "structural".equals(value)),
					"verification requires one explicit sha80 or structural value");
				verification=value;
			}else if(key.equals("evaluationgc")){
				require(evaluationGc==null && ("package".equals(value) || "off".equals(value)),
					"evaluationgc requires one explicit package or off value");
				evaluationGc=value;
			}else if(key.equals("querygc")){
				require(queryGcManifest==null && value!=null && !value.isEmpty(),"querygc requires one explicit manifest path");
				queryGcManifest=Paths.get(value);
				require(queryGcManifest.isAbsolute() && queryGcManifest.normalize().equals(queryGcManifest),
					"querygc requires a normalized absolute path");
			}else if(key.equals("matrix")){
				require(requestedMatrix==null && ProteinSubstitutionMatrix.supportedSelector(value),
					"matrix requires one supported explicit selector: "+value);requestedMatrix=value;
			}else if(key.equals("poscount")){
				require(positiveCount==null && value!=null && !value.isEmpty(),
					"poscount requires one explicit population");
				positiveCount=value;
			}else if(key.equals("negcounts")){
				require(negativeCounts==null && value!=null && !value.isEmpty(),
					"negcounts requires one explicit population list");
				negativeCounts=value;
			}else if(key.equals("pos")){positive=value;}
			else if(key.equals("neg")){negatives=value;}
			else if(key.equals("negnames")){names=value;}
			else if(key.equals("possha80")){positivePin=value;}
			else if(key.equals("negsha80")){negativePins=value;}
			else if(key.equals("insha80")){floorParentPin=floorArgument(value,key);}
			else if(key.equals("floors")){floorPath=floorArgument(value,key);}
			else if(key.equals("floorssha80")){floorPin=floorArgument(value,key);}
			else if(key.equals("fitting")){fittingPath=floorArgument(value,key);}
			else if(key.equals("fittingsha80")){fittingPin=floorArgument(value,key);}
			else if(key.equals("marginscale")){marginScale=Double.parseDouble(value);marginSpecified=true;}
			// Default-off November 2026 probe: M scales stored member lowers instead
			// of using negative uppers. The accepted 42-row grid in
			// results/member_multiplier_report_20261116 still had 2527.8..12807.75
			// negative ppm at M1.75 and was not adopted.
			else if(key.equals("membermultiplier")){memberMultiplier=Double.parseDouble(value);memberMultiplierSpecified=true;}
			else if(key.equals("paired")){paired=Parse.parseBoolean(value);}
			else if(key.equals("lengthmodel")){model=value;}
			else if(key.equals("lengthmodelsha80")){modelPin=value;}
			else if(key.equals("lengthmedians")){medians=value;}
			else if(key.equals("lengthmedianssha80")){medianPin=value;}
			else if(key.equals("top")){topOverride=Parse.parseIntKMG(value);require(topOverride>0,"top must be positive");}
			else if(key.equals("ddlreject")){
				rejectOverride=Parse.parseIntKMG(value);
				require(rejectOverride>=0,"ddlreject must be nonnegative");
			}
			else if(key.equals("heuristic")){heuristic=value;heuristicSpecified=true;}
			else if(key.equals("admission")){admission=value;}
			else if(key.equals("admissionstats")){admissionStats=Parse.parseBoolean(value);}
			else if(key.equals("acceptingrankstats")){acceptingRankStats=Parse.parseBoolean(value);}
			else if(key.equals("overhangdiagnostics")){overhangDiagnostics=Parse.parseBoolean(value);}
			else if(key.equals("stagecpu")){stageCpu=Parse.parseBoolean(value);}
			else if(key.equals("benchmarkrepeats")){
				require(value!=null && !value.isEmpty(),"benchmarkrepeats requires an explicit positive count");
				final long repeats;
				try{repeats=Parse.parseKMG(value);}
				catch(RuntimeException invalid){throw new IllegalArgumentException("Invalid benchmarkrepeats: "+value,invalid);}
				require(repeats>=1 && repeats<=Integer.MAX_VALUE,
					"benchmarkrepeats must be within1..2147483647 and include the ordinary first pass");
				benchmarkRepeats=(int)repeats;
			}
			else if(key.equals("familyhydro")){familyHydro=value;}
			else if(key.equals("familyhydrosha80")){hydroPin=value;}
			else if(key.equals("minqueryaa")){
				minQueryAa=Parse.parseIntKMG(value);
				require(minQueryAa>=1,"minqueryaa must be positive");
			}
			else if(key.equals("maxqueryaa")){
				maxQueryAa=Parse.parseIntKMG(value);
				require(maxQueryAa>=1,"maxqueryaa must be positive");
			}
			else if(key.equals("pairedscales")){scaleText=value;}
			else if(key.equals("lengthbins")){lengthEdges=parseLengthBins(value);}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		verificationMode="sha80".equals(verification) ? FamilyCalibration.VerificationMode.SHA80 : FamilyCalibration.VerificationMode.STRUCTURAL;
		require(queryGcManifest==null || structural() && !paired,"Query GC currently requires ordinary structural evaluation with exact panel counts");
		require(parser.in1!=null && parser.out1!=null && (positive!=null || negatives!=null),
			"Require in=classifier_package out=result_directory and at least one of pos= or neg=");
		require(Double.isFinite(marginScale) && marginScale>=0,"marginscale must be finite and nonnegative");
		require(!memberMultiplierSpecified || Double.isFinite(memberMultiplier) && memberMultiplier>=0,
			"membermultiplier must be finite and nonnegative");
		require(!memberMultiplierSpecified || !marginSpecified,
			"membermultiplier and explicit marginscale are mutually exclusive replay modes");
		require("off".equals(heuristic) || "hitfourth_h020_c150".equals(heuristic),
			"heuristic must be off or the released hitfourth_h020_c150 pivot");
		require("fixed".equals(admission) || "short_scaled".equals(admission) || "family_mean".equals(admission),
			"admission must be fixed, short_scaled or family_mean");
		fitted=floorPath!=null || floorPin!=null || floorParentPin!=null || fittingPath!=null || fittingPin!=null;
		evaluationGcOff="off".equals(evaluationGc);
		require(!evaluationGcOff || queryGcManifest==null && !paired && !fitted,
			"evaluationgc=off requires ordinary evaluation without querygc, paired ranking or saved floors");
		if(!structural()){require(positiveCount==null && negativeCounts==null,"poscount/negcounts require verification=structural");}
		if(fitted){
			require(floorPath!=null && floorPin!=null && floorParentPin!=null && fittingPath!=null && fittingPin!=null && marginSpecified,
				"Saved floors require floors/floorssha80, fitting/fittingsha80, insha80 and explicit absolute marginscale together");
			require(!paired && benchmarkRepeats==1 && topOverride==0 && rejectOverride<0 && admission.equals("fixed") &&
				!heuristicSpecified && model==null && modelPin==null && medians==null && medianPin==null &&
				familyHydro==null && hydroPin==null && scaleText==null,
				"A saved floor vector permits one frozen-policy evaluation, without scoring overrides or benchmark repeats");
			require(minQueryAa<0 && maxQueryAa<0 || minQueryAa==40 && maxQueryAa==99,
				"Saved short-family floors require the unchanged40..99 query bounds");
			minQueryAa=40;maxQueryAa=99;admissionStats=true;admission="family_fitted";
			verifyFloorFile(floorPath,floorPin);verifyFloorFile(fittingPath,fittingPin);
			require(floorParentPin.matches("[0-9a-f]{20}") || structural() && floorParentPin.equals("unmeasured"),
				"The fitted parent requires a lowercase sha80 or structural unmeasured label");
		}
		if(admission.equals("family_mean")){
			require(!paired && (!heuristicSpecified || heuristic.equals("hitfourth_h020_c150")) && (topOverride==0 || topOverride==50) &&
				(rejectOverride<0 || rejectOverride==6),
				"Family-mean admission requires the released new pivot, top50 and fixed6 parent anchor");
			admissionStats=true;
		}
		if(admission.equals("short_scaled")){
			require(!paired && (!heuristicSpecified || heuristic.equals("hitfourth_h020_c150")) && (topOverride==0 || topOverride==50),
				"Length-scaled short evaluation requires the released new pivot and top50, not learned paired ranking");
			require(minQueryAa<0 && maxQueryAa<0 || minQueryAa==40 && maxQueryAa==99,
				"short_scaled validates every normalized query in40..99; other query ranges are unsupported");
			minQueryAa=40;maxQueryAa=99;admissionStats=true;
			require(rejectOverride<0 || rejectOverride==6,"The scaled floor is anchored at six hits; a different fixed override is a separate comparison");
		}
		require((minQueryAa<0)==(maxQueryAa<0),"minqueryaa and maxqueryaa must be supplied together");
		queryBounds=minQueryAa>0;
		require(!queryBounds || maxQueryAa>=minQueryAa && maxQueryAa<=200,"Query bounds must be ordered within1..200; no input filtering is performed");
		require(!paired || benchmarkRepeats==1,"Repeated-work timing supports ordinary evaluation, not paired ranking");
		if(paired){
			require(!acceptingRankStats,"Accepting-rank telemetry supports ordinary single-policy evaluation only");
			require(!overhangDiagnostics,"Overhang diagnostics require ordinary single-policy frozen evaluation");
			require(!admissionStats,"Admission telemetry is supported by ordinary single-policy evaluation only");
			require(topOverride==0 && rejectOverride<0 && heuristic.equals("off") && familyHydro==null && hydroPin==null && !queryBounds,
				"Paired learned-ranking mode cannot receive ordinary cap, heuristic, KD or query-bound overrides");
			require(lengthEdges==null, "lengthbins currently reports ordinary single-policy evaluation, not paired ranking arms");
			require(!marginSpecified && !memberMultiplierSpecified && (positive==null)!=(negatives==null),
				"Paired evaluation requires positive-only or negative-only panels and pairedscales instead of marginscale or membermultiplier");
			require(model!=null && medians!=null && modelPin!=null && medianPin!=null && scaleText!=null,
				"Paired evaluation requires pinned lengthmodel/lengthmedians and three explicit pairedscales");
			final String[] fields=split(scaleText,"pairedscales");require(fields.length==3,"pairedscales must name OFF,multiply,add1 in that order");
			pairedScales=new double[3];
			for(int a=0;a<3;a++){
				pairedScales[a]=Double.parseDouble(fields[a]);
				require(Double.isFinite(pairedScales[a]) && pairedScales[a]>=0,"Paired margins must be finite and nonnegative");
			}
			verifyLengthInputs();
		}else{
			require(model==null && modelPin==null && scaleText==null,"Learned length-model inputs require explicit paired=t");
			if(heuristic.equals("off")){
				require(medians==null && medianPin==null && familyHydro==null && hydroPin==null,
					"Member medians/KD require an explicit heuristic or the separate paired model");
			}else{
				final boolean supplied=medians!=null || medianPin!=null || familyHydro!=null || hydroPin!=null;
				if(supplied){requireHeuristicInputs();}
				// An otherwise complete package may supply all metadata and query bounds.
			}
		}
		require(positive!=null || positivePin==null,"possha80 requires a positive panel");
		require(negatives!=null || names==null && negativePins==null,"Negative labels and pins require negative panels");
		require(parser.in2==null && parser.out2==null && !parser.append && parser.maxReads<0,
			"Evaluation consumes complete unpaired panels; paired input, append and read limits are unsupported");
		packagePath=Paths.get(parser.in1).toAbsolutePath().normalize();
		output=Paths.get(parser.out1).toAbsolutePath().normalize();
		threads=parser.workers()>0 ? parser.workers() : Shared.threads();
		inputThreads=parser.threadsIn>0 ? parser.threadsIn : 1;
		require(threads>0,"Evaluation needs at least one scoring worker");
		final String[] negativePaths=negatives==null ? new String[0] : split(negatives,"neg"),labels=names==null ? null : split(names,"negnames");
		final String[] pins=negativePins==null ? null : split(negativePins,"negsha80");
		require(labels==null || labels.length==negativePaths.length,"negnames must match the negative panel count");
		require(pins==null || pins.length==negativePaths.length,"negsha80 must match the negative panel count");
		require(positive!=null || positiveCount==null,"poscount requires a positive panel");
		require(negatives!=null || negativeCounts==null,"negcounts requires negative panels");
		final String[] populations=negativeCounts==null ? null : split(negativeCounts,"negcounts");
		require(populations==null || populations.length==negativePaths.length,"negcounts must match the negative panel count");
		final int offset=positive==null ? 0 : 1;
		panels=new Panel[negativePaths.length+offset];
		if(positive!=null){
			panels[0]=new Panel("positive",positive,positivePin,true,verificationMode,structural()?optionalPopulation(positiveCount):-1);
		}
		final HashSet<String> uniqueNames=new HashSet<String>(); uniqueNames.add("positive");
		final HashSet<Path> uniquePaths=new HashSet<Path>();
		if(positive!=null){uniquePaths.add(panels[0].path);}
		for(int i=0;i<negativePaths.length;i++){
			final String label=labels==null ? "negative_"+(i+1) : labels[i];
			require(uniqueNames.add(label),"Duplicate panel label: "+label);
			panels[i+offset]=new Panel(label,negativePaths[i],pins==null ? null : pins[i],false,verificationMode,
				structural() && populations!=null ? population(populations[i]) : -1);
			require(uniquePaths.add(panels[i+offset].path),"A sequence file cannot be used as two evaluation panels");
		}
		if(fitted){for(Panel panel:panels){require(panel.pin!=null,"Saved-floor evaluation requires an explicit pin for every panel");}}
		if(structural()){
			try{
				for(Panel panel:panels){panelInputs.track(panel.path);}
				if(queryGcManifest!=null){panelInputs.track(queryGcManifest);}
			}
			catch(java.io.IOException e){throw new java.io.UncheckedIOException("Cannot track structural panels",e);}
		}
		// A failed run never overwrites an accepted result directory, regardless of writer defaults.
		require(!Files.exists(output),"Evaluation output directory already exists: "+output);
	}

	public void process() throws Exception{
		final long start=System.nanoTime();
		Shared.AMINO_IN=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		ReadWrite.setZipThreads(threads);
		final FamilyScoringPackage source=new FamilyScoringPackage(packagePath,verificationMode);
		source.protectCompactOutput(output);
		require(!fitted || structural() || source.fingerprint.equals(floorParentPin),
			"Saved floors require the exact original parent recipe");
		verifyPanels();
		final FamilyHitFloorBinding floorBinding=fitted ? new FamilyHitFloorBinding(packagePath,floorParentPin,marginScale,
			Paths.get(fittingPath),fittingPin,threads,verificationMode) : null;
		final FamilyTailCascade.LoadedClassifier base=evaluationGcOff ?
			source.loadForEvaluationWithoutGc(threads,topOverride,rejectOverride) : fitted ? floorBinding.parent() : topOverride>0 || rejectOverride>=0 ?
			source.loadWithOverrides(threads,topOverride,rejectOverride) : source.load(threads);
		require(!base.usesFamilyGcPenalty() || queryGcManifest!=null,
			"GC-ranked evaluation requires querygc for every panel; use explicit evaluationgc=off only for an evaluation without GC ranking or filtering");
		require(requestedMatrix==null || requestedMatrix.equals(base.config().alignmentMetric.matrix.name),
			"Requested evaluation matrix differs from calibrated package; matrix selection cannot reuse another matrix's gates");
		final boolean inheritedRanking=inheritRanking(source);
		final boolean compact=source.compactHeader()!=null;
		require(!compact || medians==null && medianPin==null && familyHydro==null && hydroPin==null,
			"Compact packages embed member medians/KD; external evaluator member metadata cannot match the packaged ranking");
		if(paired){processPaired(source,base,start);return;}
		require(admission.equals("fixed") || heuristic.equals("hitfourth_h020_c150"),
			"Length-based admission requires an explicit or inherited common pivot");
		if(!heuristic.equals("off")){requireHeuristicInputs();}
		final String compactMargin=compact ? source.compactHeader().get("marginscale") : null;
		effectiveMarginScale=compact && !marginSpecified && compactMargin!=null ?
			Double.parseDouble(compactMargin) : marginScale;
		FamilyTailCascade.LoadedClassifier classifier=fitted ? floorBinding.load(Paths.get(floorPath),floorPin) :
			memberMultiplierSpecified ? base.withMemberCutoffMultiplier(memberMultiplier) :
			base.withEvaluationMarginScale(effectiveMarginScale);
		if(admission.equals("short_scaled")){classifier=classifier.withLengthScaledShortAdmission();}
		if(admission.equals("family_mean")){classifier=classifier.withFamilyMeanAdmission();}
		if(!heuristic.equals("off") && !inheritedRanking){
			classifier=classifier.withHeuristicRanking(new FamilyHeuristicScore(classifier,Paths.get(medians),medianPin,
				Paths.get(familyHydro),hydroPin,minQueryAa,maxQueryAa,heuristic,verificationMode));
		}
		final int families=classifier.familyCount();
		require(families>0,"The classifier package contains no live families");
		queryGc=queryGcManifest==null ? null : FamilyQueryGc.loadPanels(queryGcManifest, panels,
			new FamilyQueryGc.InputContext(panelInputs));
		final FamilyPanelRunner runner=newRunner(classifier,"Family evaluation");
		final long scoringStart=System.nanoTime(),cpuStart=processCpuNanos();
		runner.process();
		final long scoringCpu=processCpuNanos()-cpuStart,scoringWall=System.nanoTime()-scoringStart;
		require(scoringCpu>=0 && scoringWall>=0,"Evaluation scoring clocks must not move backward");
		counts=runner.counts;
		long total=0;
		for(PanelCounts panel:counts){total=Math.addExact(total,panel.queries);}
		final long queryEvaluations=Math.multiplyExact(total,(long)benchmarkRepeats);
		for(int repeat=1;repeat<benchmarkRepeats;repeat++){
			final int pass=repeat+1;
			final FamilyPanelRunner repeated=newRunner(classifier,"Family benchmark pass "+pass);
			repeated.process();
			verifyRepeatedCounts(counts,repeated.counts,pass);
		}
		final long benchmarkCpu=benchmarkRepeats==1 ? scoringCpu : processCpuNanos()-cpuStart;
		final long benchmarkWall=benchmarkRepeats==1 ? scoringWall : System.nanoTime()-scoringStart;
		require(benchmarkCpu>=scoringCpu && benchmarkWall>=scoringWall,
			"Repeated-work clocks must include the complete first scoring pass");
		source.verify();
		if(fitted){
			floorBinding.verify();verifyFloorFile(floorPath,floorPin);verifyFloorFile(fittingPath,fittingPin);
		}
		classifier.verifyHeuristicRanking();
		verifyPanels();
		if(queryGc!=null){FamilyQueryGc.verifyAll(queryGc);}
		Files.createDirectories(output);
		// A parent-only recipe must never masquerade as the complete evaluation policy.
		if(source.compactHeader()!=null){
			write(output.resolve("classifier_reference.tsv"),source.compactReferenceRecipe());
		}else{
			write(output.resolve(fitted || evaluationGcOff ? "parent_classifier.config" : "classifier.config"),
				topOverride>0 || rejectOverride>=0 ? source.recipeWithOverrides(topOverride,rejectOverride) : source.recipe());
		}
		// Persist the altered gate formula so M-mode results cannot be mistaken
		// for legacy marginscale replays.
		final ByteBuilder policy=memberMultiplierSpecified ?
			new ByteBuilder("setting\tvalue\nmember_multiplier_mode\ton\nmembermultiplier\t").append(Double.toString(memberMultiplier))
				.append("\nformula\tmax(runtime_floor,membermultiplier*stored_member_lower)\n")
				.append("identity_runtime_floor\t").append(Float.toString(classifier.config().identityFloor))
				.append("\nblosum_runtime_floor\t").append(Float.toString(classifier.config().blosumFloor)).nl()
				.append("identity_lower_rejection\tunchanged_stored_member_lower\n")
				.append("negative_derived_upper_contribution\tnone\n")
				.append("tie_definition\tfloor_binding_when_runtime_floor_ge_scaled_member_lower\n")
				.append("identity_floor_binding_families\t").append(classifier.memberMultiplierIdentityFloorBindings())
				.append("\nblosum_floor_binding_families\t").append(classifier.memberMultiplierBlosumFloorBindings())
				.append("\nzero_member_families\t").append(classifier.memberMultiplierZeroMemberFamilies())
				.append("\nalignment_unavailable_families\t").append(classifier.memberMultiplierAlignmentUnavailableFamilies())
				.append("\ncomparison\tstrict_greater;double_threshold\n")
				.append(topOverride==0 && rejectOverride<0 && heuristic.equals("off") && admission.equals("fixed") ?
					"unchanged\tranking,shortlists,scalar_gates,DDL_and_cosine_shortcuts,reference_bytes\n" :
					"unchanged\tscalar_gates,DDL_and_cosine_shortcuts,reference_bytes\n") :
			new ByteBuilder("setting\tvalue\nmarginscale\t").append(Double.toString(effectiveMarginScale))
				.append("\nformula\tbase+scale*(current-base)\nbase\tmax(member_lower,global_floor)\ncurrent\tmax(frozen_upper,global_floor)\n")
				.append("scale_one\texact_original_threshold\ncomparison\tstrict_greater;double_threshold\n")
				.append(topOverride==0 && rejectOverride<0 && heuristic.equals("off") && admission.equals("fixed") ?
					"unchanged\tmember_rejection,ranking,shortlists,scalar_gates,DDL_and_cosine_shortcuts,reference_bytes\n" :
					"unchanged\tmember_rejection,scalar_gates,DDL_and_cosine_shortcuts,reference_bytes\n");
		if(topOverride>0 || rejectOverride>=0 || !heuristic.equals("off") || queryBounds || admissionStats){
			policy.append("top\t").append(classifier.config().top).append("\nheuristic\t").append(heuristic)
				.append("\ncosine_weight\t").append(heuristic.equals("off") ? "NA" : "1.5")
				.append("\nhydrophobicity_weight\t").append(heuristic.equals("off") ? "NA" : "0.2")
				.append("\nlengthmedians\t").append(compactRanking ? "embedded" : medians==null ? "NA" : medians)
				.append("\nlengthmedians_sha80\t").append(compactRanking ? "embedded" : medianPin==null ? "NA" : medianPin)
				.append("\nfamilyhydro\t").append(compactRanking ? "embedded" : familyHydro==null ? "NA" : familyHydro)
				.append("\nfamilyhydro_sha80\t").append(compactRanking ? "embedded" : hydroPin==null ? "NA" : hydroPin)
				.append("\nmin_query_aa\t").append(queryBounds ? Integer.toString(minQueryAa) : "NA")
				.append("\nmax_query_aa\t").append(queryBounds ? Integer.toString(maxQueryAa) : "NA")
				.append("\nquery_length_policy\tvalidate_only_no_filter_no_truncate\nraw_hit_admission\t")
				.append(fitted ? "ddl_ge_pinned_family_floor" : admission.equals("fixed") ? "ddl_ge_reject" : admission.equals("family_mean") ?
					"ddl_ge_family_mean_band_4_5_6" : "ddl_ge_ceil_6_times_length_minus_4_over_95")
				.append("\nddl_reject\t")
				.append(classifier.config().ddlReject).append("\nscalar_filter_order\tafter_ranked_cap\nno_backfill\ttrue\n");
		}
		if(admissionStats){
			policy.append("admission\t").append(admission).append("\nadmission_length\t")
				.append(fitted ? "frozen_live_roster_not_length_bands" : admission.equals("family_mean") ? "retained_family_residues_per_member_exact_integer_comparison" : "normalized_residues_not_distinct_or_unmasked_kmers").nl()
				.append("pool_scope\tall_raw_hit_admitted_families_before_ranking_cap_and_scalar_gates\n");
			if(admission.equals("family_mean")){
				policy.append("family_mean_bands\tmean<65:4;65<=mean<100:5;mean>=100:6\n")
					.append("family_mean_source\tpinned_reference_members_and_residues\n")
					.append("floor_range_scope\tconfigured_live_family_roster;not_one_query_floor\n");
			}
		}
		if(fitted){
			policy.append("floor_parent\t").append(packagePath.toString()).append("\nfloor_parent_sha80\t").append(floorParentPin)
				.append("\nfloors\t").append(floorPath).append("\nfloors_sha80\t").append(floorPin)
				.append("\nfitting_manifest\t").append(fittingPath).append("\nfitting_manifest_sha80\t").append(fittingPin)
				.append("\nclassifier_recipe_role\tparent_only;floors_and_margin_are_required\n")
				.append("initial_support\tfrozen_from_fitting;never_reestimated_from_evaluation\n");
		}
		if(structural()){policy.append("verification\tstructural\ninput_identities\tinherited_or_unmeasured_not_recomputed\n");}
		if(evaluationGcOff){policy.append("evaluation_gc\toff\n")
			.append("family_gc_policy\texplicitly_disabled_for_this_evaluation;package_unchanged\n")
			.append("classifier_recipe_role\tparent_only;evaluationgc=off_is_required\n");}
		if(queryGc!=null){policy.append("query_gc_manifest\t").append(queryGcManifest.toString())
			.append("\nquery_gc_loads\t1\nquery_gc_identity\texact_declared_panel;count;original_ordinal_and_first_token\nmissing_gc\texplicit_neutral_no_inference\n");}
		write(output.resolve("evaluation_policy.tsv"),policy);
		if(acceptingRankStats){
			write(output.resolve("accepting_rank_policy.tsv"),new ByteBuilder("#family_accepting_rank_policy\t1\nsetting\tvalue\n"+
				"rank\tone_based_final_sorted_shortlist_before_scalar_and_GC_filtering\n"+
				"denominator\tranked_accepts_only;unranked_shortcut_accepts_reported_separately\n"+
				"rejected_queries\tno_accepting_rank\nempty_ranked_population\tNA\n"+
				"median\taverage_of_middle_two_order_statistics\np90\tnearest_rank_ceil_0.9N\n"+
				"scoring\tordinary_single_pass_no_extra_alignments\n"));
		}
		if(overhangDiagnostics){
			write(output.resolve("overhang_diagnostic_policy.tsv"),new ByteBuilder(
				"metric\tdefinition\ncompared_dp_pairs\tactual candidate DPs, including failed spans; excludes unavailable calibration and short-input skips\n"+
				"trim_positive_overhang_nonpositive\tsame query/reference/matrix, maximal-trim raw>0 and overhang raw<=0; denominator compared_dp_pairs (also report maximal_trim_positive)\n"+
				"identity_path_removed\tspan-qualified candidate passing the identity route but raw<=0; denominator identity_path_eligible (also report span_qualified_pairs)\n"+
				"scope\tone frozen evaluation pass per view/panel, not unique queries or calibration-round totals\n"+
				"timing\tdiagnostic scoring includes one extra maximal-trim DP per compared pair\n"));
			runner.writeOverhangDiagnostics(output.resolve("overhang_diagnostics.tsv"));
		}
		runner.writeSummary(output.resolve("summary.tsv"));
		if(queryGc!=null){runner.writeGcCoverage(output.resolve("gc_coverage.tsv"));}
		if(lengthEdges!=null){runner.writeLengths(output.resolve("length_summary.tsv"));}
		if(admissionStats){runner.writeAdmission(output.resolve("admission_summary.tsv"));runner.writeAlignments(output.resolve("alignment_counts.tsv"));}
		if(acceptingRankStats){runner.writeAcceptingRanks(output.resolve("accepting_rank_summary.tsv"));}
		if(stageCpu){runner.writeStageCpu(output.resolve("stage_cpu.tsv"));}
		runner.writeFamilies(output.resolve("family_false_accepts.tsv.gz"));writeInputs(source);
		final double seconds=(System.nanoTime()-start)*1e-9;
		write(output.resolve("timing.tsv"),new ByteBuilder("metric\tvalue\nwall_seconds\t").appendSlow(seconds)
			.nl().append("threads\t").append(threads).nl());
		write(output.resolve("scoring_time.tsv"),new ByteBuilder("metric\tvalue\nwall_seconds\t").appendSlow(scoringWall*1e-9)
			.append("\nprocess_cpu_seconds\t").appendSlow(scoringCpu*1e-9).append("\nqueries\t").append(total)
			.append("\nworker_threads\t").append(threads).nl());
		if(benchmarkRepeats>1){
			write(output.resolve("benchmark_time.tsv"),new ByteBuilder("metric\tvalue\nwall_seconds\t").appendSlow(benchmarkWall*1e-9)
				.append("\nprocess_cpu_seconds\t").appendSlow(benchmarkCpu*1e-9)
				.append("\npasses\t").append(benchmarkRepeats).append("\nunique_queries\t").append(total)
				.append("\nquery_evaluations\t").append(queryEvaluations).append("\nworker_threads\t").append(threads)
				.append("\nclassifier_loads\t1\nscientific_counts\tfirst_pass_only\nper_pass_count_parity\texact\n")
				.append("timing_scope\tordinary_scoring_passes_including_repeat_setup_and_parity\n")
				.append("repeated_queries_are_independent_samples\tfalse\n"));
		}
		if(queryGc!=null){FamilyQueryGc.verifyAll(queryGc);}
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_EVAL_PASS panels=").append(panels.length)
			.append(" queries=").append(total).append(" live_families=").append(families).nl());
		outstream.println("Family evaluation completed: "+total+" queries, "+panels.length+" panels, "+families+" live families.");
	}

	private boolean structural(){return verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL;}
	private static long population(String value){
		require(value!=null && value.matches("[0-9]+"),"Panel population must be a positive decimal integer");
		final long count=Long.parseLong(value);require(count>0,"Panel population must be positive");return count;
	}
	private static long optionalPopulation(String value){return value==null ? -1 : population(value);}
	/** Structural panel identities and counts are explicit; metadata stability is checked
	 * from declaration through the final scoring pass, without changing denominators. */
	private void verifyPanels()throws Exception{
		for(Panel panel:panels){panel.verify();}
		if(structural()){panelInputs.verify();}
	}

	/** Reuse the immutable loaded classifier; each pass retains ordinary input and worker lifecycle. */
	private FamilyPanelRunner newRunner(FamilyTailCascade.LoadedClassifier classifier,String label){
		assert(classifier!=null && classifier.familyCount()>0) : "Benchmark passes require the validated live classifier used by the first pass";
		final FamilyPanelRunner runner=new FamilyPanelRunner(classifier,panels,threads,inputThreads,null,outstream,label,lengthEdges,
			queryBounds ? minQueryAa : 0,queryBounds ? maxQueryAa : Integer.MAX_VALUE,admissionStats,stageCpu,false,1,queryGc);
		if(overhangDiagnostics){runner.enableOverhangDiagnostics();}
		if(acceptingRankStats){runner.enableAcceptingRankStats();}
		return runner;
	}

	/** Repeated work must reproduce each panel, terminal stage and winning-family count exactly. */
	private void verifyRepeatedCounts(PanelCounts[] expected,PanelCounts[] observed,int pass){
		assert(expected!=null && observed!=null && pass>1) : "Parity compares a completed extra pass against the retained first pass";
		require(expected.length==observed.length,"Benchmark pass "+pass+" changed the panel roster");
		for(int panel=0;panel<expected.length;panel++){
			final PanelCounts a=expected[panel],b=observed[panel];
			require(a.queries==b.queries && a.residues==b.residues && a.accepted==b.accepted &&
				a.measuredGcQueries==b.measuredGcQueries && a.neutralGcQueries==b.neutralGcQueries &&
				a.measuredGcAccepts==b.measuredGcAccepts && a.neutralGcAccepts==b.neutralGcAccepts &&
				Arrays.equals(a.stages,b.stages) && Arrays.equals(a.winners,b.winners),
				"Benchmark pass "+pass+" changed panel "+panels[panel].name+" counts: queries "+a.queries+"/"+b.queries+
				", accepted "+a.accepted+"/"+b.accepted+"; terminal stages and family winners must also agree");
		}
	}

	/** A calibrated package's ranking is part of its recipe, not the evaluator's default CLI flag.
	 * Explicit duplicate bindings must agree; silently replacing fitted ranking is unsupported.
	 */
	private boolean inheritRanking(FamilyScoringPackage source){
		final java.util.HashMap<String,String> binding=source.rankingOptions();
		if(binding.isEmpty()){return false;}
		require(!paired,"A packaged heuristic cannot be replaced by the unrelated learned paired-ranking study");
		require(!heuristicSpecified || heuristic.equals(binding.get("heuristic")),
			"Explicit evaluator heuristic must agree with the calibrated package");
		final int min=Integer.parseInt(binding.get("minqueryaa")),max=Integer.parseInt(binding.get("maxqueryaa"));
		require(!queryBounds || minQueryAa==min && maxQueryAa==max,"Evaluator query bounds must agree with the packaged ranking contract");
		if(source.compactHeader()!=null){
			require(medians==null && medianPin==null && familyHydro==null && hydroPin==null,
				"Compact packages embed member medians/KD; external evaluator member metadata cannot match the packaged ranking");
			heuristic=binding.get("heuristic");
			minQueryAa=min; maxQueryAa=max; queryBounds=true; compactRanking=true;
			return true;
		}
		if(medians!=null){
			require(Paths.get(medians).toAbsolutePath().normalize().toString().equals(binding.get("lengthmedians")) &&
				medianPin.equals(binding.get("lengthmedianssha80")) &&
				Paths.get(familyHydro).toAbsolutePath().normalize().toString().equals(binding.get("familyhydro")) &&
				hydroPin.equals(binding.get("familyhydrosha80")),"Evaluator member metadata must match the calibrated package binding");
		}
		heuristic=binding.get("heuristic"); medians=binding.get("lengthmedians"); medianPin=binding.get("lengthmedianssha80");
		familyHydro=binding.get("familyhydro"); hydroPin=binding.get("familyhydrosha80");
		minQueryAa=min; maxQueryAa=max; queryBounds=true;
		return true;
	}
	private void requireHeuristicInputs(){
		if(compactRanking){return;}
		require(queryBounds && medians!=null && medianPin!=null && familyHydro!=null && hydroPin!=null,
			"Heuristic evaluation needs complete CLI or packaged member medians/KD and query bounds");
		require(Paths.get(medians).isAbsolute() && Paths.get(familyHydro).isAbsolute(),
			"Heuristic member tables must use absolute paths for reproducible policy records");
	}

	/** Explicitly frozen policies; this invocation never selects a scale from a panel. */
	private void processPaired(FamilyScoringPackage source,FamilyTailCascade.LoadedClassifier base,long start) throws Exception{
		verifyLengthInputs();
		final FamilyLengthScore binding=source.bindLengthRanking(base,Paths.get(model),modelPin,Paths.get(medians),medianPin,
			verificationMode,structural() ? lengthInputs : null);
		final FamilyTailCascade.LengthPairedClassifier factory=base.lengthPaired(binding,pairedScales);
		final FamilyPairedPanelRunner runner=new FamilyPairedPanelRunner(base,panels,threads,inputThreads,
			()->new PairedScorer(factory.newWorker()),outstream,"Paired family evaluation");
		runner.process();
		long total=0;for(FamilyPairedCounts panel:runner.counts){total=Math.addExact(total,panel.queries());}
		source.verify();binding.verify();verifyLengthInputs();verifyPanels();
		Files.createDirectories(output);
		write(output.resolve("classifier.config"),source.recipe());
		final ByteBuilder policy=new ByteBuilder("setting\tvalue\npaired\ttrue\nmodel\t").append(model)
			.append("\nmodel_sha80\t").append(modelPin).append("\nmedians\t").append(medians).append("\nmedians_sha80\t").append(medianPin)
			.append("\nverification\t").append(structural() ? "structural\nidentity_semantics\tinherited_unverified" :
				"sha80\nidentity_semantics\tmeasured_sha80")
			.append("\narms\toff,multiply,add1\nscore_forms\tcurrent;raw_hits*f;raw_hits+f\nscale_selection\texternal_negative_only_frozen_before_positives\n")
			.append("oracle_scope\tunion_of_three_strict_top50_lists_before_scalar_gates;not_exhaustive_truth\n")
			.append("oracle_tie\tlower_physical_family_rank\noracle_metric\tglocal_BLOSUM_raw_score/full_query_self_score\n")
			.append("shortcut_policy\tstudy_requires_ddlaccept_1025;no_DDL_shortcut\n");
		for(int a=0;a<3;a++){policy.append("margin_").append(FamilyPairedCounts.ARM_NAMES[a]).tab().append(Double.toString(pairedScales[a])).nl();}
		policy.append(binding.recipe());
		write(output.resolve("evaluation_policy.tsv"),policy);
		runner.writeSummary(output.resolve("summary.tsv"),pairedScales);
		runner.writeFamilies(output.resolve("family_false_accepts.tsv.gz"));
		runner.writePairs(output.resolve("paired_transitions.tsv"),output.resolve("oracle_membership.tsv"));
		writeInputs(source);
		write(output.resolve("timing.tsv"),new ByteBuilder("metric\tvalue\nwall_seconds\t").appendSlow((System.nanoTime()-start)*1e-9)
			.nl().append("threads\t").append(threads).nl());
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_PAIRED_EVAL_PASS panels=").append(panels.length)
			.append(" queries=").append(total).append(" arms=3 live_families=").append(base.familyCount()).nl());
		outstream.println("Paired family evaluation completed: "+total+" queries across three frozen ranking arms.");
	}

	/** One reusable synchronous observer per scoring worker; no callback allocation per query. */
	private static final class PairedScorer implements FamilyPairedPanelRunner.Scorer,FamilyTailCascade.LengthPairedObserver{
		PairedScorer(FamilyTailCascade.LengthPairedWorker worker_){worker=worker_;}
		@Override public void score(String id,byte[] raw,FamilyPairedCounts target){
			require(current==null && target!=null,"Paired scorer must not be re-entered or lack its panel counts");
			current=target;try{worker.assignObserve(id,raw,this);}finally{current=null;}
		}
		@Override public void accept(int length,String[] terminals,int[] winners,int[] ddl,float[] identities,float[] blosum,
				int oracle,int mask,int size){current.observe(length,terminals,winners,oracle,mask,size);}
		private final FamilyTailCascade.LengthPairedWorker worker;
		private FamilyPairedCounts current;
	}

	private void verifyLengthInputs(){
		for(String[] input:new String[][]{{model,modelPin},{medians,medianPin}}){
			require(input[0]!=null && input[1]!=null &&
				(input[1].matches("[0-9a-f]{20}") || structural() && input[1].equals("unmeasured")),
				"Ranking inputs need paths and sha80 fingerprints or structural unmeasured labels");
			final Path path=Paths.get(input[0]);
			require(path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path),
				"Frozen ranking input must be an absolute regular file: "+path);
			if(structural()){
				try{lengthInputs.track(path);}
				catch(java.io.IOException e){throw new IllegalArgumentException("Frozen ranking input changed: "+path,e);}
			}
			else{require(input[1].equals(DigestSuffix.file(path.toString())),"Frozen ranking input changed: "+path);}
		}
	}

	private void writeInputs(FamilyScoringPackage source){
		final ByteBuilder out=new ByteBuilder("panel\tpath\texpected_sha80\npackage\t").append(source.path.toString()).tab().append(source.fingerprint).nl();
		if(source.compactHeader()!=null){
			if(source.compactReference()!=null){
				out.append("compact_reference\t").append(source.compactReference().toString()).append("\tunmeasured\n");
			}
			final Path[] inputs=source.compactInputs();
			final String[] names={"compact_family","compact_alphabet","compact_blacklist","compact_matrix"};
			for(int i=0;i<inputs.length;i++){
				out.append(names[i]).tab().append(inputs[i].toString()).append("\tunmeasured\n");
			}
		}
		if(fitted){out.append("floors\t").append(floorPath).tab().append(floorPin).nl()
			.append("fitting_manifest\t").append(fittingPath).tab().append(fittingPin).nl();}
		for(Panel panel:panels){out.append(panel.name).tab().append(panel.path.toString()).tab().append(panel.pin==null ? "NA" : panel.pin).nl();}
		if(queryGc!=null){out.append("query_gc_manifest\t").append(queryGcManifest.toString()).append("\tunmeasured\n");}
		write(output.resolve("inputs.tsv"),out);
		if(structural()){
			final ByteBuilder populations=new ByteBuilder("panel\texpected_queries\n");
			for(Panel panel:panels){
				populations.append(panel.name).tab();
				if(panel.expectedQueries<0){populations.append("NA");}
				else{populations.append(panel.expectedQueries);}
				populations.nl();
			}
			write(output.resolve("input_populations.tsv"),populations);
		}
		write(output.resolve("interpretation.txt"),new ByteBuilder("Positive TPR measures acceptance by any family; no family-truth labels are inferred.\n")
			.append("Each negative panel has its own observed denominator. A reused calibration panel is not independent cold validation.\n")
			.append("Per-family false accepts count terminal winners, including zero-count live families; empty reference models are excluded.\n"));
	}

	/** A named binding option cannot silently disappear through a missing value. */
	private static String floorArgument(final String value,final String key){
		require(value!=null && !value.isEmpty(),"Missing saved-floor argument: "+key);return value;
	}
	/** Verify the immutable final table/provenance without rereading fitting sequences. */
	private void verifyFloorFile(final String name,final String pin){
		require(name!=null && pin!=null && (pin.matches("[0-9a-f]{20}") || structural() && pin.equals("unmeasured")),
			"Saved-floor inputs require paths and lowercase sha80 fingerprints or structural unmeasured labels");
		final Path path=Paths.get(name);
		require(path.isAbsolute() && Files.isRegularFile(path) && !Files.isSymbolicLink(path),
			"Saved-floor input must be an absolute regular non-symlink file: "+path);
		if(structural()){
			try{floorInputs.track(path);}
			catch(java.io.IOException e){throw new IllegalArgumentException("Saved-floor input changed: "+path,e);}
		}
		else{require(pin.equals(DigestSuffix.file(name)),"Saved-floor input changed: "+path);}
	}

	private static String[] split(String text,String name){
		require(text!=null && !text.isEmpty(),"Missing "+name);
		final String[] values=text.split(",",-1);
		for(String value:values){require(!value.isEmpty(),"Empty entry in "+name);}
		return values;
	}

	/** Boundaries affect reporting only; positive increasing edges include explicit outside bins. */
	private static int[] parseLengthBins(String text){
		final String[] values=split(text, "lengthbins");
		require(values.length<=1000, "Length reports support at most1000 boundaries to bound per-worker counters");
		final int[] edges=new int[values.length];
		for(int i=0; i<edges.length; i++){
			final long value;
			try{value=Parse.parseKMG(values[i]);}
			catch(RuntimeException invalid){
				throw new IllegalArgumentException("Invalid lengthbins boundary: "+values[i], invalid);
			}
			require(value>0 && value<=Integer.MAX_VALUE && (i==0 || value>edges[i-1]),
				"lengthbins requires strictly increasing positive integer boundaries");
			edges[i]=(int)value;
		}
		return edges;
	}
	private static void write(Path path,ByteBuilder content){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);writer.start();
		writer.print(content); require(!writer.poisonAndWait(),"Output failed: "+path);
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	/** Whole-JVM CPU includes workers, GC and input threads over the same scoring interval as wall time. */
	private static long processCpuNanos(){
		final java.lang.management.OperatingSystemMXBean bean=java.lang.management.ManagementFactory.getOperatingSystemMXBean();
		require(bean instanceof com.sun.management.OperatingSystemMXBean,"A process CPU clock is required for evaluation timing");
		final long value=((com.sun.management.OperatingSystemMXBean)bean).getProcessCpuTime();
		require(value>=0,"Process CPU clock is unavailable");return value;
	}
	private final FamilyCalibration.VerificationMode verificationMode;
	private final boolean evaluationGcOff;
	private final StructuralInputs panelInputs=new StructuralInputs();
	private final StructuralInputs floorInputs=new StructuralInputs();
	private final StructuralInputs lengthInputs=new StructuralInputs();
	private final Path packagePath,output;
	private final Panel[] panels;
	private Path queryGcManifest;
	private FamilyQueryGc[] queryGc;
	private final int threads,inputThreads;
	private final PrintStream outstream;
	private PanelCounts[] counts;
	/** Only memberMultiplierSpecified activates this value; NaN poisons accidental inactive arithmetic. */
	private double marginScale=1,memberMultiplier=Double.NaN;
	private double effectiveMarginScale=1;
	private boolean paired,marginSpecified,memberMultiplierSpecified,stageCpu;
	private String model,modelPin,medians,medianPin,scaleText;
	private double[] pairedScales;
	private int[] lengthEdges;
	private String heuristic="off",familyHydro,hydroPin;
	private boolean heuristicSpecified;
	private boolean compactRanking;
	private String admission="fixed";
	private boolean admissionStats;
	private boolean acceptingRankStats;
	private boolean overhangDiagnostics;
	private int topOverride,minQueryAa=-1,maxQueryAa=-1;
	private int rejectOverride=-1;
	private int benchmarkRepeats=1;
	private boolean queryBounds;
	private boolean fitted;
	private String floorPath,floorPin,floorParentPin,fittingPath,fittingPin;
}
