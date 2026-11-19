package prot;

import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Locale;

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
import ddl.DDLFamilySearch.StructuralInputs;

/** Iterative negative-only calibration using one loaded classifier and shared panel execution.
 * Each round scores one immutable parent; all observations merge before any gate changes.
 * @author Brian, Keqing
 */
public final class FamilyCalibrate{
	private String requestedMatrix;

	public static void main(String[] args) throws Exception{
		try{
			final FamilyCalibrate tool=new FamilyCalibrate(args);
			try{tool.process();}finally{Shared.closeStream(tool.outstream);}
		}catch(Resources.MissingResourceException e){
			System.err.print(e.getMessage());
			System.exit(1);
		}
	}

	public FamilyCalibrate(String[] args){
		final PreParser pp=new PreParser(args,getClass(),false);outstream=pp.outstream;
		final Parser parser=new Parser();
		String negatives=null,names=null,pins=null,prepare=null;boolean targetRateSpecified=false;
		for(String arg:pp.args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT),value=pair.length>1 ? pair[1] : null;
			if(key.equals("prepare")){prepare=value;}
			else if(key.equals("querygc")){
				require(queryGcManifest==null && value!=null && !value.isEmpty(),"querygc requires one explicit manifest path");
				queryGcManifest=Paths.get(value);
				require(queryGcManifest.isAbsolute() && queryGcManifest.normalize().equals(queryGcManifest),"querygc requires a normalized absolute path");
			}
			else if(key.equals("matrix")){
				require(requestedMatrix==null && ProteinSubstitutionMatrix.supportedSelector(value),
					"matrix requires one supported explicit selector: "+value);requestedMatrix=value;
			}
			else if(key.equals("verification")){
				require("sha80".equals(value) || "structural".equals(value),"verification must be sha80 or structural");
				verificationMode="structural".equals(value) ? FamilyCalibration.VerificationMode.STRUCTURAL : FamilyCalibration.VerificationMode.SHA80;
			}
			else if(key.equals("heldout")){heldout=value;}
			else if(key.equals("heldoutsha80")){heldoutPin=value;}
			else if(key.equals("heldoutqueries")){heldoutQueries=Parse.parseKMG(value);}
			else if(key.equals("neg")){negatives=value;}
			else if(key.equals("negnames")){names=value;}
			else if(key.equals("negsha80")){pins=value;}
			else if(key.equals("expectedinitialaccepts")){
				final String[] counts=split(value,"expectedinitialaccepts");expectedInitialAccepts=new long[counts.length];
				for(int i=0;i<counts.length;i++){
					expectedInitialAccepts[i]=Parse.parseKMG(counts[i]);
					require(expectedInitialAccepts[i]>=0,"Expected initial accepts must be nonnegative: "+counts[i]);
				}
			}else if(key.equals("expectedpanelqueries")){
				final String[] values=split(value,"expectedpanelqueries");
				panelQueryCounts=new long[values.length];
				for(int i=0;i<values.length;i++){
					panelQueryCounts[i]=Parse.parseKMG(values[i]);
					require(panelQueryCounts[i]>0,"expectedpanelqueries must contain positive complete-panel populations");
				}
				if(values.length==1){expectedPanelQueries=panelQueryCounts[0];panelQueryCounts=null;}
			}
			else if(key.equals("targetfpr")){target=Double.parseDouble(value);targetRateSpecified=true;}
			else if(key.equals("targetaccepts")){targetAccepts=Long.parseLong(value);require(targetAccepts>=0,"targetaccepts must be nonnegative");}
			else if(key.equals("targetacceptspooled")){
				require(targetAcceptsPooled<0,"targetacceptspooled may be specified only once");
				targetAcceptsPooled=Long.parseLong(value);
				require(targetAcceptsPooled>=0,"targetacceptspooled must be a nonnegative integer");
			}
			else if(key.equals("targetacceptsbyclass")){
				require(targetAcceptsByClass==null,"targetacceptsbyclass may be specified only once");
				final String[] values=split(value,"targetacceptsbyclass");targetAcceptsByClass=new long[values.length];
				for(int i=0;i<values.length;i++){
					targetAcceptsByClass[i]=Long.parseLong(values[i]);
					require(targetAcceptsByClass[i]>=0,"Per-class acceptance ceilings must be nonnegative integers");
				}
			}
			else if(key.equals("marginscale")){marginScale=Double.parseDouble(value);}
			else if(key.equals("perclass")){perClass=Parse.parseBoolean(value);}
			else if(key.equals("rounds")){
				final long n=Parse.parseKMG(value);require(n>=1 && n<=32,"rounds must be 1..32 evaluations, including the initial one");
				rounds=(int)n;
			}else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown calibration argument: "+arg);}
		}
		require(parser.out1!=null && (prepare==null ? parser.in1!=null && negatives!=null :
			parser.in1==null && negatives==null && names==null && pins==null),
			"Require out=new_directory and either in=package neg=panels or prepare=manifest, without mixed inputs");
		require(parser.in2==null && parser.out2==null && !parser.append && parser.maxReads<0,
			"Calibration consumes complete unpaired negative panels; paired input, append and read limits are unsupported");
		require(Double.isFinite(target) && target>=0 && target<=1,"targetfpr must be finite and in [0,1]");
		require(Double.isFinite(marginScale) && marginScale>0 && (prepare==null || marginScale==1),
			"A positive finite calibration margin is supported for an existing package only");
		require(targetAccepts<0 || !targetRateSpecified && prepare==null && perClass,
			"An exact acceptance count replaces targetfpr and requires existing-package per-class calibration");
		require(targetAcceptsByClass==null || !targetRateSpecified && targetAccepts<0 && prepare==null && perClass,
			"Per-class acceptance counts replace targetfpr/targetaccepts and require existing-package per-class calibration");
		require(targetAcceptsPooled<0 || targetAcceptsByClass!=null && !targetRateSpecified && targetAccepts<0 && prepare==null && perClass,
			"targetacceptspooled requires explicit existing-package per-class acceptance counts");
		final boolean anyHeldout=heldout!=null || heldoutPin!=null || heldoutQueries!=-1;
		require(queryGcManifest==null || verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL && prepare==null && !anyHeldout,
			"Query GC currently requires structural existing-package calibration without combined preparation or held-out scoring");
		if(queryGcManifest!=null){
			try{queryGcInputs.track(queryGcManifest);}
			catch(java.io.IOException problem){throw new java.io.UncheckedIOException("Cannot track the query GC manifest", problem);}
		}
		require(!anyHeldout || verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL && prepare!=null &&
			heldout!=null && heldoutPin!=null && heldoutPin.matches("[0-9a-f]{20}") && heldoutQueries>0 &&
			Paths.get(heldout).isAbsolute() && target==.0001 && rounds<=4,
			"Resident held-out scoring requires structural preparation, an absolute inherited panel, exact query count, at most four rounds and targetfpr=.0001");
		packagePath=prepare==null ? Paths.get(parser.in1).toAbsolutePath().normalize() : null;
		output=Paths.get(parser.out1).toAbsolutePath().normalize();
		require(!Files.exists(output),"Calibration output directory already exists: "+output);
		threads=parser.workers()>0 ? parser.workers() : Shared.threads();
		inputThreads=parser.threadsIn>0 ? parser.threadsIn : 1;
		require(threads>0,"Calibration requires at least one scoring worker");
		preparation=prepare==null ? null : new FamilyPreparation(Paths.get(prepare),verificationMode);
		require(preparation==null || requestedMatrix==null || requestedMatrix.equals(preparation.matrixName()),
			"Requested matrix differs from preparation recipe; old member bounds cannot select custom scores");
		if(preparation!=null){panels=preparation.panels;validateExpectedCounts();return;}
		final String[] paths=split(negatives,"neg"),labels=names==null ? null : split(names,"negnames");
		final String[] expected=pins==null ? null : split(pins,"negsha80");
		require(labels==null || labels.length==paths.length,"negnames must match the negative panel count");
		require(expected==null || expected.length==paths.length,"negsha80 must match the negative panel count");
		panels=new Panel[paths.length];
		require(panelQueryCounts==null || panelQueryCounts.length==paths.length,
			"Expected panel counts must be one shared count or one count per panel");
		require(targetAccepts<0 || paths.length==1 && expectedPanelQueries>0 && targetAccepts<=expectedPanelQueries,
			"Exact-count repair requires one explicitly sized negative class");
		require(targetAcceptsByClass==null || labels!=null && targetAcceptsByClass.length==paths.length &&
			(panelQueryCounts!=null || paths.length==1 && expectedPanelQueries>0),
			"Per-class counts require ordered explicit names, one target and one population per negative panel");
		final HashSet<String> uniqueNames=new HashSet<String>();uniqueNames.add("positive");
		final HashSet<Path> uniquePaths=new HashSet<Path>();
		for(int i=0;i<paths.length;i++){
			final String label=labels==null ? "negative_"+(i+1) : labels[i];
			require(uniqueNames.add(label),"Duplicate or reserved negative panel label: "+label);
			panels[i]=new Panel(label,paths[i],expected==null ? defaultPanelPin(paths[i]) : expected[i],false,
				verificationMode,expectedQueries(i));
			require(uniquePaths.add(panels[i].path),"A negative file cannot be counted as two panels");
		}
		validateExpectedCounts();
	}

	/** Positive panels never guide calibration; optional held-out scoring happens after the saved negative gate. */
	public void process() throws Exception{
		final long start=System.nanoTime();
		Shared.AMINO_IN=true;ReadWrite.ZIPLEVEL=9;ReadWrite.FORCE_BGZIP=true;ReadWrite.USE_BGZIP=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;ReadWrite.setZipThreads(threads);
		final boolean structural=verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL;
		final FamilyScoringPackage source=preparation==null ? new FamilyScoringPackage(packagePath,verificationMode) : null;
		//TODO: Compact recalibration needs a compatible writer; recipe() serializes legacy packages only.
		require(source==null || source.compactHeader()==null,
			"Compact family packages are not supported for recalibration. Use in=classifier.config with its legacy components, "+
			"or prepare= with a complete member-calibration recipe.");
		final StructuralInputs panelState=structural ? queryGcManifest==null ? new StructuralInputs() : queryGcInputs : null;
		if(queryGcManifest!=null){panelState.verify();}
		for(Panel panel:panels){if(panelState!=null){panelState.track(panel.path);}panel.verify();}
		FamilyTailCascade.LoadedClassifier current=preparation==null ? source.load(threads) : preparation.build(output,threads,outstream);
		final FamilyReference[] initialFamilies=targetAcceptsPooled<0 ? null : current.config().families;
		require(!current.usesFamilyGcPenalty() || queryGcManifest!=null,"GC-ranked calibration requires querygc for every negative panel");
		final FamilyQueryGc[] queryGc=queryGcManifest==null ? null : FamilyQueryGc.loadPanels(queryGcManifest, panels,
			new FamilyQueryGc.InputContext(panelState));
		require(requestedMatrix==null || requestedMatrix.equals(current.config().alignmentMetric.matrix.name),
			"Requested calibration matrix differs from the calibrated package");
		final java.util.HashMap<String,String> ranking=source==null ? null : source.rankingOptions();
		final int minimumQuery=ranking!=null && ranking.containsKey("minqueryaa") ? Integer.parseInt(ranking.get("minqueryaa")) : 0;
		final int maximumQuery=ranking!=null && ranking.containsKey("maxqueryaa") ? Integer.parseInt(ranking.get("maxqueryaa")) : Integer.MAX_VALUE;
		FamilyPanelRunner last=null;
		final ByteBuilder history=new ByteBuilder("round\tqueries\tfalse_accepts\tfpr\tddl_accepts\tchanged_families\tstatus\twall_seconds\n");
		final ByteBuilder cohorts=new ByteBuilder("round\tpanel\tqueries\tfalse_accepts\tfpr\tddl_accepts\tidentity_accepts\tblosum_accepts\n");
		final ByteBuilder gcHistory=queryGc==null ? null : new ByteBuilder("round\tpanel\tgc_coverage\tqueries\taccepted\n");
		final long[] observedPanelQueries=new long[panels.length];
		String status=null;int evaluations=0,updateRounds=0;long familyUpdates=0;
		for(int round=1;round<=rounds;round++){
			final long roundStart=System.nanoTime();
			final FamilyCalibrationUpdates updates=current.newCalibrationUpdates();
			last=new FamilyPanelRunner(current,panels,threads,inputThreads,updates,outstream,"Calibration round "+round,
				null,minimumQuery,maximumQuery,false,false,false,marginScale,queryGc);
			last.process();evaluations=round;
			if(panelState!=null){panelState.verify();}
			long queries=0,accepted=0,shortcuts=0;boolean classesPass=true,shortcutClassFails=false;
			for(int p=0;p<panels.length;p++){
				final PanelCounts c=last.counts[p];
				final long expectedQueries=expectedQueries(p);
				require(expectedQueries<0 || c.queries==expectedQueries,
					"Calibration panel population changed before gate update: round="+round+" panel="+panels[p].name+
					" expected="+expectedQueries+" observed="+c.queries);
				if(round==1){observedPanelQueries[p]=c.queries;}
				else{require(c.queries==observedPanelQueries[p],
					"Calibration panel population changed after first observation: round="+round+" panel="+panels[p].name+
					" expected="+observedPanelQueries[p]+" observed="+c.queries);}
				require(round!=1 || expectedInitialAccepts==null || c.accepted==expectedInitialAccepts[p],
					"Initial calibration accepts changed before gate update: panel="+panels[p].name+
					" expected="+(expectedInitialAccepts==null ? -1 : expectedInitialAccepts[p])+" observed="+c.accepted);
				queries=Math.addExact(queries,c.queries);accepted=Math.addExact(accepted,c.accepted);shortcuts=Math.addExact(shortcuts,c.stages[1]);
				classesPass&=withinClass(c.accepted,c.queries,p);shortcutClassFails|=!withinClass(c.stages[1],c.queries,p);
				cohorts.append(round).tab().append(panels[p].name).tab().append(c.queries).tab().append(c.accepted)
					.tab().appendSlow(c.accepted/(double)c.queries).tab().append(c.stages[1]).tab().append(c.stages[2]).tab().append(c.stages[3]).nl();
				if(gcHistory!=null){
					for(int measured=1; measured>=0; measured--){
						gcHistory.append(round).tab().append(panels[p].name).tab().append(measured==1 ? "measured" : "neutral")
							.tab().append(last.gcQueries(p, measured==1)).tab().append(last.gcAccepts(p, measured==1)).nl();
					}
				}
			}
			int proposed=0;
			for(int f=0;f<current.familyCount();f++){if(updates.changed(f)){proposed++;}}
			status=((perClass ? classesPass : within(accepted,queries)) &&
				(targetAcceptsPooled<0 || accepted<=targetAcceptsPooled)) ? "converged" :
				((perClass ? shortcutClassFails : !within(shortcuts,queries)) ||
				(targetAcceptsPooled>=0 && shortcuts>targetAcceptsPooled)) ? "ddl_floor_exceeds_target" :
				round==rounds ? "round_limit" : proposed==0 ? "no_progress" : "raise";
			final boolean raise=status.equals("raise");
			history.append(round).tab().append(queries).tab().append(accepted).tab().appendSlow(accepted/(double)queries)
				.tab().append(shortcuts).tab().append(raise ? proposed : 0).tab().append(status)
				.tab().appendSlow((System.nanoTime()-roundStart)*1e-9).nl();
			outstream.println("Calibration round "+round+": "+accepted+" accepts / "+queries+" negatives; "+status+".");
			if(!raise){break;}
			current=current.withCalibrationUpdates(updates);updateRounds++;familyUpdates=Math.addExact(familyUpdates,proposed);
		}
		require(last!=null && status!=null && !status.equals("raise"),"Every bounded calibration ends with an evaluated terminal parent");
		if(preparation==null){Files.createDirectories(output.getParent());Files.createDirectory(output);}
		final Path config=output.resolve("classifier.config");
		if(preparation!=null){preparation.writeFinal(output,current);}
		else if(updateRounds==0){write(config,source.recipe());}
		else{
			final Path references=output.resolve("references.tsv.gz");current.writeReferences(references);
			write(config,source.recipe(references,structural ? "unmeasured" : DigestSuffix.file(references.toString())));
		}
		last.writeSummary(output.resolve("summary.tsv"));last.writeFamilies(output.resolve("family_false_accepts.tsv.gz"));
		if(queryGc!=null){last.writeGcCoverage(output.resolve("gc_coverage.tsv"));write(output.resolve("gc_rounds.tsv"), gcHistory);}
		write(output.resolve("rounds.tsv"),history);write(output.resolve("cohorts.tsv"),cohorts);
		write(output.resolve("calibration.tsv"),new ByteBuilder("metric\tvalue\nstatus\t").append(status)
			.nl().append("evaluations\t").append(evaluations).nl().append("max_evaluations\t").append(rounds)
			.nl().append(targetAcceptsByClass!=null ? "target_accepts_by_class\t"+countVector() :
				targetAccepts<0 ? "target_fpr\t"+Double.toString(target) : "target_accepts\t"+targetAccepts)
			.nl().append("target_scope\t").append(targetAcceptsPooled>=0 ? "pooled_and_per_class_counts" : targetAcceptsByClass!=null ? "per_class_counts" :
				targetAccepts>=0 ? "per_class_count" : perClass ? "per_class" : "pooled")
			.nl().append("marginscale\t").appendSlow(marginScale)
			.nl().append("update_rounds\t").append(updateRounds).nl().append("family_updates\t").append(familyUpdates)
			.nl().append("threads\t").append(threads).nl().append("wall_seconds\t").appendSlow((System.nanoTime()-start)*1e-9).nl());
		if(targetAcceptsPooled>=0){
			long queries=0,accepted=0;
			for(PanelCounts counts:last.counts){queries=Math.addExact(queries,counts.queries);accepted=Math.addExact(accepted,counts.accepted);}
			write(output.resolve("pooled_target.tsv"),new ByteBuilder("queries\tmaximum_accepts\tobserved_accepts\tmaximum_fpr\tobserved_fpr\n")
				.append(queries).tab().append(targetAcceptsPooled).tab().append(accepted)
				.tab().appendSlow(targetAcceptsPooled/(double)queries).tab().appendSlow(accepted/(double)queries).nl());
			writeGateChanges(initialFamilies,current);
		}
		if(targetAcceptsByClass!=null){
			final ByteBuilder ceilings=new ByteBuilder("panel\tqueries\tmaximum_accepts\tobserved_accepts");
			if(targetAcceptsPooled>=0){ceilings.append("\tmaximum_fpr\tobserved_fpr");}ceilings.nl();
			for(int p=0;p<panels.length;p++){
				ceilings.append(panels[p].name).tab().append(expectedQueries(p)).tab().append(targetAcceptsByClass[p])
					.tab().append(last.counts[p].accepted);
				if(targetAcceptsPooled>=0){ceilings.tab().appendSlow(targetAcceptsByClass[p]/(double)last.counts[p].queries)
					.tab().appendSlow(last.counts[p].accepted/(double)last.counts[p].queries);}ceilings.nl();
			}
			write(output.resolve("targets.tsv"),ceilings);
		}
		// This reads only the final recipe and verifies file pins; prepared arrays are never reloaded for scoring.
		final FamilyScoringPackage provenance=structural ? null : source==null ? new FamilyScoringPackage(config) : source;
		final ByteBuilder inputs=new ByteBuilder("panel\tpath\tsha80\npackage\t")
			.append(structural ? config.toString() : provenance.path.toString()).tab().append(structural ? "unmeasured" : provenance.fingerprint).nl();
		for(Panel p:panels){
			inputs.append(p.name).tab().append(p.path.toString()).tab().append(p.pin==null ? "unmeasured" : p.pin).nl();
		}
		if(queryGc!=null){inputs.append("query_gc_manifest\t").append(queryGcManifest.toString()).append("\tunmeasured\n");}
		write(output.resolve("inputs.tsv"),inputs);
		write(output.resolve("interpretation.txt"),new ByteBuilder("These are errors on reused calibration panels, not independent cold false-positive rates.\n")
			.append("Each round scores one immutable parent; worker/panel maxima merge before gates change. Positive data is never used.\n")
			.append(targetAcceptsByClass==null ? "" : "Per-class counts are ceilings, not equality promises. Updates still merge accepted winners from every panel; realized counts can undershoot their ceilings.\n")
			.append(targetAcceptsPooled<0 ? "" : "The pooled count ceiling and every per-class count ceiling must both pass; either DDL shortcut-floor violation prevents convergence.\n")
			.append("Identity increments are 1 percentage point; BLOSUM increments are 0.005, using terminal float scores before upward rounding.\n")
			.append("Increments apply to effective gates at marginscale=").appendSlow(marginScale)
			.append("; stored uppers are increased by inverse margin conversion. Reload evaluation must use this same margin.\n")
			.append("PASS means a completed calculation. Only status=converged creates FROZEN; other statuses retain the last evaluated package without claiming convergence.\n")
			.append(preparation==null ? "The package reuses pinned ranking dependencies; no member-based preparation was requested.\n" :
				"Preparation passes references and modal profiles directly in memory; only final package/audit files are retained. Only the explicitly named panel supplies raw cosine maxima.\n"));
		if(provenance!=null){provenance.verify();}for(Panel panel:panels){panel.verify();}
		if(structural && source!=null){source.verify();}
		if(panelState!=null){panelState.verify();}
		if(queryGc!=null){FamilyQueryGc.verifyAll(queryGc);}
		if(preparation!=null){preparation.verify();}
		if(structural){
			final String blacklistIdentity=preparation==null ? null : preparation.blacklistIdentity();
			final ByteBuilder verification=new ByteBuilder("key\tvalue\nverification\tstructural\ninput_identity\t")
				.append(preparation==null ? "inherited_or_unmeasured_not_recomputed" : preparation.inputIdentity())
				.append("\noutput_byte_identity\tunmeasured\nserialized_reload_validated\tfalse\n")
				.append(preparation==null ? "parent_recipe\t" : "preparation_recipe\t")
				.append(preparation==null ? source.path.toString() : preparation.recipePath().toString())
				.append(queryGc==null ? "\nquery_identity_checks\tnonempty;population_count;reader_order_preserved;no_independent_ID_sidecar_join\n" :
					"\nquery_identity_checks\texact_declared_panel;population_count;original_ordinal_and_first_token_GC_join\nquery_gc_loads\t1\nmissing_gc\texplicit_neutral_no_inference\n");
			if(blacklistIdentity!=null){verification.append("blacklist_input_identity\t").append(blacklistIdentity).nl();}
			write(output.resolve("verification.tsv"),verification);
		}
		if(status.equals("converged")){write(output.resolve("FROZEN"),new ByteBuilder(config.toString()).tab().append(structural ? "unmeasured" : DigestSuffix.file(config.toString())).nl());}
		if(heldout!=null){evaluateHeldout(current,last,status);}
		if(queryGc!=null){FamilyQueryGc.verifyAll(queryGc);}
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_CALIBRATE_PASS status=").append(status).append(" evaluations=").append(evaluations).nl());
		outstream.println("Calibration completed: "+status+" after "+evaluations+" evaluations; package "+config+".");
	}

	/** One final roster comparison counts each changed family once, separately from round updates. */
	private void writeGateChanges(FamilyReference[] initial,FamilyTailCascade.LoadedClassifier current){
		final FamilyReference[] result=current.config().families;
		require(initial!=null && initial.length==result.length,"Gate summary requires the same initial/final live roster");
		int changed=0,identityChanged=0,blosumChanged=0;
		final double[] sums=new double[4],maxima=new double[4],deltas=new double[4];
		for(int i=0;i<initial.length;i++){
			final FamilyReference before=initial[i],after=result[i];
			require(before.id.equals(after.id) && before.rank==after.rank && before.members==after.members &&
				before.identityLower==after.identityLower && before.blosumLower==after.blosumLower &&
				after.identityUpper>=before.identityUpper && after.blosumFinal>=before.blosumFinal,
				"Gate summary found changed roster/member lowers or decreasing gates");
			final double id=(double)after.identityUpper-before.identityUpper,bl=(double)after.blosumFinal-before.blosumFinal;
			if(id>0 || bl>0){changed++;}if(id>0){identityChanged++;}if(bl>0){blosumChanged++;}
			final double idEffective=FamilyTailCascade.evaluationGate(after.identityLower,after.identityUpper,current.config().identityFloor,marginScale)-
				FamilyTailCascade.evaluationGate(before.identityLower,before.identityUpper,current.config().identityFloor,marginScale);
			final double blEffective=FamilyTailCascade.evaluationGate(after.blosumLower,after.blosumFinal,current.config().blosumFloor,marginScale)-
				FamilyTailCascade.evaluationGate(before.blosumLower,before.blosumFinal,current.config().blosumFloor,marginScale);
			deltas[0]=id;deltas[1]=bl;deltas[2]=idEffective;deltas[3]=blEffective;
			for(int d=0;d<deltas.length;d++){sums[d]+=deltas[d];maxima[d]=Math.max(maxima[d],deltas[d]);}
		}
		final ByteBuilder text=new ByteBuilder("metric\tvalue\nlive_families\t").append(initial.length)
			.append("\nunique_changed_families\t").append(changed).append("\nidentity_changed_families\t").append(identityChanged)
			.append("\nblosum_changed_families\t").append(blosumChanged).append("\nmarginscale\t").appendSlow(marginScale)
			.append("\nmember_lowers_unchanged\ttrue\n");
		final String[] names={"identity_stored","blosum_stored","identity_effective","blosum_effective"};
		for(int d=0;d<names.length;d++){
			text.append(names[d]).append("_delta_sum\t").appendSlow(sums[d]).nl();
			text.append(names[d]).append("_delta_max\t").appendSlow(maxima[d]).nl();
		}
		write(output.resolve("gate_change_summary.tsv"),text);
	}

	/** The saved negative gate precedes all positive-file access; this path never changes the classifier. */
	private void evaluateHeldout(FamilyTailCascade.LoadedClassifier frozen,FamilyPanelRunner negative,String calibrationStatus)throws Exception{
		require(verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL && preparation!=null && heldout!=null,
			"Resident held-out evaluation must belong to its completed structural preparation");
		boolean allowed=calibrationStatus.equals("converged");
		for(PanelCounts panel:negative.counts){allowed&=within(panel.accepted,panel.queries);}
		final String reason=!calibrationStatus.equals("converged") ? "calibration_not_converged" : allowed ? "every_class_within_ceiling" : "negative_class_exceeds_ceiling";
		write(output.resolve("heldout_gate.tsv"),new ByteBuilder("key\tvalue\nallowed\t").append(allowed)
			.append("\nreason\t").append(reason).append("\ncalibration_status\t").append(calibrationStatus)
			.append("\nceiling_per_class_per_million\t100\nmarginscale\t1\npositive_guided_selection\tfalse\n"));
		if(!allowed){
			write(output.resolve("heldout_status.tsv"),new ByteBuilder("key\tvalue\nstatus\tskipped\nreason\t").append(reason).nl());
			return;
		}
		// Construction itself checks the path; even this is deliberately after the persisted gate.
		final Panel positive=new Panel("positive",heldout,heldoutPin,true,verificationMode,heldoutQueries);
		positive.verify();
		final FamilyPanelRunner runner=new FamilyPanelRunner(frozen,new Panel[]{positive},threads,inputThreads,
			outstream,"Resident held-out evaluation",preparation.minimumHeldoutQueryAa(),preparation.maximumHeldoutQueryAa(),false);
		final long wallStart=System.nanoTime(),cpuStart=processCpuNanos();
		runner.process();
		final long wall=System.nanoTime()-wallStart,cpu=processCpuNanos()-cpuStart;
		require(wall>=0 && cpu>=0,"Resident scoring clocks must not move backward");
		positive.verify();preparation.verify();frozen.verifyHeuristicRanking();
		final Path result=output.resolve("heldout");Files.createDirectory(result);
		runner.writeSummary(result.resolve("summary.tsv"));
		write(result.resolve("scoring_time.tsv"),new ByteBuilder("metric\tvalue\nwall_seconds\t").appendSlow(wall*1e-9)
			.append("\nprocess_cpu_seconds\t").appendSlow(cpu*1e-9).append("\nqueries\t").append(runner.completedQueries(0))
			.append("\nworker_threads\t").append(threads).nl());
		write(result.resolve("inputs.tsv"),new ByteBuilder("panel\tpath\tinherited_sha80\npositive\t").append(positive.path.toString()).tab().append(positive.pin).nl());
		write(result.resolve("evaluation_policy.tsv"),new ByteBuilder("setting\tvalue\nmarginscale\t1\nverification\tstructural\nclassifier_source\tfrozen_resident_calibration\npositive_guided_selection\tfalse\nserialized_reload_validated\tfalse\n"));
		write(result.resolve("PASS"),new ByteBuilder("FAMILY_EVAL_PASS panels=1 queries=").append(runner.completedQueries(0))
			.append(" live_families=").append(frozen.familyCount()).nl());
		write(output.resolve("heldout_status.tsv"),new ByteBuilder("key\tvalue\nstatus\tcompleted\nreason\t").append(reason).nl());
	}
	private static long processCpuNanos(){
		final java.lang.management.OperatingSystemMXBean bean=java.lang.management.ManagementFactory.getOperatingSystemMXBean();
		require(bean instanceof com.sun.management.OperatingSystemMXBean,"A process CPU clock is required for resident scoring");
		final long value=((com.sun.management.OperatingSystemMXBean)bean).getProcessCpuTime();
		require(value>=0,"Process CPU clock is unavailable");return value;
	}

	/** An expected first-round vector binds the complete panel order before expensive scoring. */
	private void validateExpectedCounts(){
		require(panelQueryCounts==null || panelQueryCounts.length==panels.length,
			"Expected panel counts must match the complete ordered panel roster");
		require(expectedInitialAccepts==null || expectedInitialAccepts.length==panels.length,
			"expectedinitialaccepts must have one count per negative panel, in panel order");
		if(expectedInitialAccepts!=null){
			for(int i=0;i<expectedInitialAccepts.length;i++){require(expectedQueries(i)<0 || expectedInitialAccepts[i]<=expectedQueries(i),
				"Expected initial accepts cannot exceed the declared complete-panel population");}
		}
		if(targetAcceptsByClass!=null){
			long pooledQueries=0;
			for(int i=0;i<targetAcceptsByClass.length;i++){
				require(expectedQueries(i)>0 && targetAcceptsByClass[i]<=expectedQueries(i),
					"Per-class acceptance ceilings cannot exceed their declared complete-panel populations");
				if(targetAcceptsPooled>=0){pooledQueries=Math.addExact(pooledQueries,expectedQueries(i));}
			}
			require(targetAcceptsPooled<0 || targetAcceptsPooled<=pooledQueries,
				"Pooled acceptance ceiling cannot exceed the declared complete-panel population sum");
		}
	}
	private long expectedQueries(int panel){return panelQueryCounts==null ? expectedPanelQueries : panelQueryCounts[panel];}
	private String defaultPanelPin(String path){return verificationMode==FamilyCalibration.VerificationMode.SHA80 ? DigestSuffix.file(path) : null;}
	/** Compare one complete named class with its integer ceiling, retaining scalar behavior when absent. */
	private boolean withinClass(long accepted,long queries,int panel){
		assert(queries>0 && accepted>=0 && accepted<=queries) : "Calibration targets compare observed finite panel counts";
		return targetAcceptsByClass==null ? within(accepted,queries) : accepted<=targetAcceptsByClass[panel];
	}
	/** Render the validated count vector in panel order for the final provenance table only. */
	private String countVector(){
		assert(targetAcceptsByClass!=null && targetAcceptsByClass.length>0) : "Count provenance requires a nonempty validated vector";
		final StringBuilder text=new StringBuilder();
		for(long count:targetAcceptsByClass){if(text.length()>0){text.append(',');}text.append(count);}
		return text.toString();
	}
	private boolean within(long accepted,long queries){
		assert(queries>0 && accepted>=0 && accepted<=queries) : "Calibration targets compare observed finite panel counts";
		return targetAccepts>=0 ? accepted<=targetAccepts : accepted/(double)queries<=target;
	}
	private static String[] split(String text,String name){
		require(text!=null && !text.isEmpty(),"Missing "+name);final String[] values=text.split(",",-1);
		for(String value:values){require(!value.isEmpty(),"Empty entry in "+name);}return values;
	}
	private static void write(Path path,ByteBuilder text){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);writer.start();
		writer.print(text);require(!writer.poisonAndWait(),"Calibration output failed: "+path);
	}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final Path packagePath,output;
	private final FamilyPreparation preparation;
	private final Panel[] panels;
	private Path queryGcManifest;
	private final StructuralInputs queryGcInputs=new StructuralInputs();
	private final int threads,inputThreads;
	private final PrintStream outstream;
	private double target=.0001;
	private double marginScale=1;
	private long targetAccepts=-1;
	private long targetAcceptsPooled=-1;
	private long[] targetAcceptsByClass;
	private int rounds=4;
	private boolean perClass;
	private FamilyCalibration.VerificationMode verificationMode=FamilyCalibration.VerificationMode.STRUCTURAL;
	private String heldout,heldoutPin;
	private long heldoutQueries=-1;
	private long[] expectedInitialAccepts;
	private long[] panelQueryCounts;
	private long expectedPanelQueries=-1;
}
