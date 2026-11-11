package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import parse.PreParser;
import parse.Parse;
import shared.Shared;
import shared.Tools;

/**
 * Runs the development-only multi-K comparison with optional learned rescue on
 * the first assembly arm.
 * <p>
 * The initial live branch gate runs while the seed k-mer table is still owned by
 * the initial {@link Tadpole} instance.  The optional postpass rescue transforms
 * the initial raw contigs before ordinary graph cleanup and later-K table
 * construction.  All arms retain formed short contigs and use the same cleanup
 * order, so neither network observes truth labels, later-K tables or filtered
 * FASTA intermediates.
 * @author Fischl
 */
public final class MultiKBranchRescue implements TadpoleMulti.InitialAssembly{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Expands standard configuration files before selecting the two optional bindings.
	 * @param args command line arguments
	 */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MultiKBranchRescue.class, false);
		try{new MultiKBranchRescue(pp.args).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses experiment-specific model bindings, then forwards the remaining
	 * arguments into the ordinary multi-K configuration.
	 * <p>
	 * The constructor owns output-mode and cleanup-order switches that would make
	 * the raw neural-model comparison ambiguous, but leaves ordinary Tadpole and
	 * TadpoleMulti assembly knobs visible for paired fixture recipes.
	 * @param args command line arguments after {@link PreParser} expansion
	 */
	public MultiKBranchRescue(final String[] args){
		if(args==null){throw new IllegalArgumentException("Expected reads, out=FASTA and a multi-K schedule.");}
		final ArrayList<String> nativeArgs=new ArrayList<String>(Arrays.asList(
				"retainshortcontigs=t", "trimends=0", "trimcircular=f", "simdfma=f"));
		boolean overlapGuard=false;
		boolean branchGuard=false;
		boolean initialBranch=false;
		boolean shortEvidence=false;
		float initialRatio=0;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			String key=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			while(key.startsWith("-")){key=key.substring(1);}
			final String raw=eq<0 ? null : arg.substring(eq+1);
			final String value=raw!=null && raw.equalsIgnoreCase("null") ? null : raw;
			if(key.equals("livenet")){liveModel=value;}
			else if(key.equals("livecutoff")){liveCutoff=Float.parseFloat(value);}
			else if(key.equals("postnet")){postModel=value;}
			else if(key.equals("postcutoff")){postCutoff=Float.parseFloat(value);}
			else if(key.equals("graphoverlapguard")){overlapGuard=Parse.parseBoolean(value);}
			else if(key.equals("graphbranchguard")){branchGuard=Parse.parseBoolean(value);}
			else if(key.equals("initialbranchguard")){initialBranch=Parse.parseBoolean(value);}
			else if(key.equals("fuseshortevidence")){shortEvidence=Parse.parseBoolean(value);}
			else if(key.equals("initialmergecoverageratio")){initialRatio=Float.parseFloat(value);}
			else if((key.startsWith("out") && !key.equals("out")) || key.equals("dot") || key.equals("gfa")
					|| key.equals("hist") || key.equals("khist") || key.equals("dump") || key.equals("ihist")
					|| key.equals("inserthistogram") || key.equals("fusionvectors") || key.equals("processcontigs")
					|| key.startsWith("contigrescue")){
				throw new IllegalArgumentException("This experiment supports out=FASTA only, and owns initial cleanup ordering: "+key);
			}else{nativeArgs.add(arg);}
		}
		validateBinding(liveModel, liveCutoff, "live");
		validateBinding(postModel, postCutoff, "post");
		if(initialRatio!=0 && (!Float.isFinite(initialRatio) || initialRatio<1)){
			throw new IllegalArgumentException("initialmergecoverageratio must be zero (disabled) or finite and at least 1.");
		}
		initialCoverageRatio=initialRatio;
		initialBranchGuard=initialBranch;
		config=new TadpoleMulti.Config(nativeArgs.toArray(new String[0]));
		if(overlapGuard && !config.finalGraphNeeded()){
			throw new IllegalArgumentException("graphoverlapguard requires an existing final graph phase; it never adds a table load.");
		}
		config.graphMergeOverlapGuard=overlapGuard;
		if(branchGuard && !config.finalGraphNeeded()){
			throw new IllegalArgumentException("graphbranchguard requires an existing final graph phase; it never adds a table load.");
		}
		config.graphMergeBranchGuard=branchGuard;
		config.fuseShortEvidence=shortEvidence;
		if(!config.ordered || config.assembleK<24 || config.assembleK>127 || config.out==null
				|| !FileFormat.isFastaExt(ReadWrite.rawExtension(config.out))){
			throw new IllegalArgumentException("Use korder=input, assemblek24-127 and out=FASTA for the frozen branch models.");
		}
	}

	/** Requires a model and a finite cutoff to appear as a pair.
	 * @param model model path, or {@code null} when this rescue stage is disabled
	 * @param cutoff inclusive probability cutoff for the model
	 * @param name stage name used in error messages
	 */
	private static void validateBinding(final String model, final float cutoff, final String name){
		if(model==null ? !Float.isNaN(cutoff) : model.isEmpty() || !Float.isFinite(cutoff) || cutoff<0 || cutoff>1){
			throw new IllegalArgumentException(name+"net and "+name+"cutoff in[0,1] must be specified together.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Executes the multi-K run while preserving the process-wide direct graph
	 * gates that the initial assembly temporarily overrides.
	 * <p>
	 * Concurrent invocations in the same JVM are unsupported because those gates
	 * are static {@link BubblePopper} hooks.
	 */
	private void run(){
		previousCoverageGate=BubblePopper.directCoverageGate;
		previousBranchGate=BubblePopper.directBranchGate;
		try{new TadpoleMulti(config, this).process();}
		finally{
			BubblePopper.directCoverageGate=previousCoverageGate;
			BubblePopper.directBranchGate=previousBranchGate;
		}
	}

	/**
	 * Validates raw-unitig invariants, installs the initial hooks and binds all
	 * rescue models before assembly read streams or output writers open.
	 * <p>
	 * Only the initial assembly table may own these hooks.  The live branch gate
	 * must see the initial table's k-mers and ownership bits directly, while the
	 * postpass transform must receive the same raw contig list immediately before
	 * ordinary graph cleanup prepares contigs for the later-K schedule.
	 * @param tad initial Tadpole instance that owns the seed table
	 * @param schedule multi-K configuration that created this callback
	 */
	@Override
	public void prepare(final Tadpole tad, final TadpoleMulti.Config schedule){
		assert(schedule==config && tad!=null) : "This binding belongs only to this schedule's initial assembler.";
		if(tad.hashKmerMode!=Tadpole.HASH_EXPLICIT || tad.processingMode!=Tadpole.contigMode ||
				!tad.tables().rcomp() || AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN || Shared.SIMD_FMA ||
				!tad.retainShortContigs || tad.trimEnds!=0 || tad.trimCircular || tad.minExtension<2 ||
				!tad.useOwnership || tad.minCountExtend<1 || tad.minContigLen<1 || !tad.processContigs){
			throw new IllegalArgumentException("Require explicit canonical DNA, retained untrimmed unitigs, minextension>=2, mce>=1, ownership=t, simdfma=f and ordinary initial cleanup.");
		}
		if(tad.ecc || tad.ecco || tad.merge || tad.markErrors || tad.hpIndel || tad.localEdit ||
				tad.extendLeft>0 || tad.extendRight>0 || tad.tossJunk || tad.discardLowDepthReads>=0 || tad.discardUncorrectable ||
				tad.evictLowDepthContigs || Tadpole.append){
			throw new IllegalArgumentException("Preprocess reads separately; no early contig eviction or append in this comparison.");
		}
		final ArrayList<String> inputs=new ArrayList<String>(tad.tables().in1);
		inputs.addAll(tad.tables().in2); inputs.addAll(tad.tables().extra);
		if(inputs.isEmpty()){throw new IllegalArgumentException("At least one read input is required.");}
		if(liveModel!=null){inputs.add(liveModel);}
		if(postModel!=null){inputs.add(postModel);}
		if(config.fuseNet!=null){inputs.add(config.fuseNet);}
		TadpoleGraph.checkPaths(inputs.get(0), inputs, config.out);
		if(!Tools.testInputFiles(true, true, inputs.toArray(new String[0])) ||
				!Tools.testOutputFiles(false, false, false, config.out)){
			throw new IllegalArgumentException("All inputs must be readable; output must be new and distinct from every input/model.");
		}
		assert(tad.liveBranchGate==null && tad.rawContigTransform==null) : "Initial rescue hooks must have one owner.";
		live=liveModel==null ? null : new LiveBranchGate(tad, liveModel, liveCutoff);
		final SameKJoinRescue post=postModel==null ? null : new SameKJoinRescue(tad, postModel, postCutoff);
		tad.liveBranchGate=live;
		if(post!=null){
			tad.rawContigTransform=contigs -> post.mergeLoaded(contigs,
					live==null ? SameKJoinFeatures.moments(tad) : live.populationMoments());
		}
		if(initialCoverageRatio>0){
			initialCoverageGate=new BubblePopper.CoverageGate(initialCoverageRatio);
			BubblePopper.directCoverageGate=initialCoverageGate;
		}
		if(initialBranchGuard){
			initialBranchGate=new GraphBranchGate();
			BubblePopper.directBranchGate=initialBranchGate;
		}
		Tadpole.outstream.println("MULTIK_BRANCH_BINDING live="+(live!=null)+" postpass="+(post!=null)+
				" assemblek="+tad.k()+" cleanup="+(initialCoverageGate==null ? "unchanged" : "coverage_guard")+
				" retained_short=true");
		if(initialBranchGate!=null){Tadpole.outstream.println("MULTIK_INITIAL_BRANCH_BINDING enabled=true");}
	}

	/**
	 * Reports learned-gate counters after the initial process2 cleanup completes
	 * and before the initial table is released for the later-K arms.
	 * <p>
	 * Initial direct graph gates are restored here as soon as their protected
	 * cleanup has finished.  {@link #run()} repeats the same restoration in a
	 * {@code finally} block so aborted runs also leave the process-wide gates as
	 * they were before this experiment started.
	 * @param tad initial Tadpole instance that has just finished cleanup
	 */
	@Override
	public void finish(final Tadpole tad){
		assert(tad!=null && tad.k()==config.assembleK) : "Only the initial assembly may report live rescue counters.";
		if(initialBranchGate!=null){
			assert(BubblePopper.directBranchGate==initialBranchGate) :
					"Initial branch policy must remain installed until finish(), before any later-K phase.";
			BubblePopper.directBranchGate=previousBranchGate;
			Tadpole.outstream.println("MULTIK_INITIAL_BRANCH evaluations="+initialBranchGate.evaluations+
					" rejected="+initialBranchGate.rejected);
		}
		if(initialCoverageGate!=null){
			assert(BubblePopper.directCoverageGate==initialCoverageGate) :
					"Initial cleanup must retain its own policy until TadpoleMulti invokes finish(), before later-K phases.";
			BubblePopper.directCoverageGate=previousCoverageGate;
			Tadpole.outstream.println("MULTIK_INITIAL_COVERAGE ratio="+initialCoverageGate.ratio+
					" evaluations="+initialCoverageGate.evaluations+" rejected="+initialCoverageGate.rejected);
		}
		final long[] totals=live==null ? new long[2] : live.totals();
		Tadpole.outstream.println("MULTIK_LIVE_BRANCH evaluations="+totals[0]+" rescued_encounters="+totals[1]);
		if(live!=null){Tadpole.outstream.println("LIVE_BRANCH_MOMENTS seconds="+(live.momentNanos()*1e-9));}
		live=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Frozen multi-K schedule after private experiment arguments are removed. */
	private final TadpoleMulti.Config config;
	/** Maximum high/low coverage ratio for optional initial bubble popping, or zero. */
	private final float initialCoverageRatio;
	/** True when initial graph cleanup rejects branch-conflicting merges. */
	private final boolean initialBranchGuard;
	/** Active and saved initial-coverage gates for process-wide BubblePopper state. */
	private BubblePopper.CoverageGate initialCoverageGate, previousCoverageGate;
	/** Active and saved initial-branch gates for process-wide BubblePopper state. */
	private GraphBranchGate initialBranchGate, previousBranchGate;
	/** Optional neural model paths for initial live rescue and raw-contig postpass rescue. */
	private String liveModel, postModel;
	/** Inclusive cutoffs paired with {@link #liveModel} and {@link #postModel}. */
	private float liveCutoff=Float.NaN, postCutoff=Float.NaN;
	/** Installed live gate, retained only long enough to report encounter counts. */
	private LiveBranchGate live;
}
