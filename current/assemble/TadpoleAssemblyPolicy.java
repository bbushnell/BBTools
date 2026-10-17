package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.Data;
import kmer.AbstractKmerTableSet;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * User-facing assembly defaults, separate from constructor and experimental APIs.
 * The postpass borrows the initial table; later phases retain their own policies.
 * Explicit options override defaults regardless of argument order.
 * @author Fischl
 */
final class TadpoleAssemblyPolicy implements TadpoleMulti.InitialAssembly, AutoCloseable {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves policy options without opening streams or loading a count table. */
	TadpoleAssemblyPolicy(final String[] input){this(input, false);}

	/** Also handles the explicit tadpolemulti.sh entry point with a single listed K. */
	private TadpoleAssemblyPolicy(final String[] input, final boolean forceMulti){
		assert(input!=null) : "Policy selection follows config and PacBio expansion.";
		assembly=isAssembly(input);
		final ArrayList<String> explicit=new ArrayList<String>();
		for(String arg : input){
			final String name=key(arg), value=value(arg);
			if(name.equals("assemblydefaults")){defaults=Parse.parseBoolean(value);}
			else if(name.equals("branchrescue")){
				rescue=("auto".equalsIgnoreCase(value) ? -1 : Parse.parseBoolean(value) ? 1 : 0);
			}else if(name.equals("branchrescuenet")){
				if(value==null || value.isEmpty()){
					throw new IllegalArgumentException("branchrescuenet requires a model path or null.");
				}
				postModel="null".equalsIgnoreCase(value) ? null : value;
				modelSet=true;
			}else if(name.equals("branchrescuecutoff")){
				postCutoff=Float.parseFloat(value);
			}else if(name.equals("initialmergecoverageratio")){
				initialRatio=Float.parseFloat(value);
				if(!Float.isFinite(initialRatio)){
					throw new IllegalArgumentException("initialmergecoverageratio must be finite.");
				}
			}else{explicit.add(arg);}
		}
		if(!Float.isFinite(postCutoff) || postCutoff<0 || postCutoff>1){
			throw new IllegalArgumentException("branchrescuecutoff must be finite and in [0,1].");
		}
		if(Float.isNaN(initialRatio)){initialRatio=assembly && defaults ? 1.75f : 0;}
		if(initialRatio!=0){BubblePopper.CoverageGate.validateRatio(initialRatio);}
		if(!assembly && (rescue==1 || modelSet || initialRatio!=0)){
			throw new IllegalArgumentException("Branch rescue and initial merge coverage apply only to contig assembly.");
		}
		rescueRequested=assembly && (rescue==1 || (rescue<0 && (defaults || modelSet))) && postModel!=null;
		multi=forceMulti || TadpoleMulti.hasMultipleK(explicit.toArray(new String[0]));
		final ArrayList<String> resolved=new ArrayList<String>();
		if(assembly && (defaults || rescueRequested)){
			resolved.addAll(Arrays.asList("retainshortcontigs=t", "trimends=0", "trimcircular=f",
					"minprob=0", "minprobmain=f", "simdfma=f"));
		}
		if(assembly && defaults && multi){
			resolved.addAll(Arrays.asList("fusecoverageratio=1.75", "graphmergecoverageratio=1.75",
					"fuseshortevidence=t", "graphbranchguard=t"));
		}
		if(assembly && multi){
			String fusion=defaults ? "auto" : null;
			boolean cutoff=false;
			for(String arg : explicit){
				if(key(arg).equals("fusenet")){fusion=value(arg);}
				if(key(arg).equals("fusencutoff")){cutoff=true;}
			}
			if("auto".equalsIgnoreCase(fusion)){
				bundledFusion=true;
				resolved.add("fusenet="+modelPath("?tadpole_fusion.bbnet"));
				if(!cutoff){resolved.add("fusencutoff=0.667098");}
				// Auto is a public shorthand, not a literal filename for Config.
				for(int i=explicit.size()-1; i>=0; i--){
					if(key(explicit.get(i)).equals("fusenet")){explicit.remove(i);}
				}
			}
		}
		resolved.addAll(explicit);
		args=resolved.toArray(new String[0]);
		// Protect bundled/custom models even from outstream's early PreParser open.
		if(rescueRequested){postModel=modelPath(postModel); checkModelOutputs(postModel, input);}
		for(String arg : args){
			if(key(arg).equals("fusenet") && !"null".equalsIgnoreCase(value(arg))){
				checkModelOutputs(value(arg), input);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Execution          ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs public single/multi-K entry points without changing lower-level callers. */
	static void run(final String[] input, final boolean forceMulti){
		try(TadpoleAssemblyPolicy policy=new TadpoleAssemblyPolicy(input, forceMulti)){
			if(forceMulti || policy.multi){
				if(TadpoleGraph.requested(policy.args)){
					throw new IllegalArgumentException("Graph-only mode requires one K, not a multi-K schedule.");
				}
				new TadpoleMulti(policy.multiConfig(), policy).process();
			}else{
				final Timer total=new Timer(), initialization=new Timer();
				final Tadpole tad=Tadpole.makeTadpole(policy.args, true);
				policy.prepare(tad, null);
				tad.printExecutionPlan();
				initialization.stop();
				Tadpole.outstream.println("Initialization Time:      \t"+initialization);
				try{tad.process(total);}
				finally{policy.finish(tad); tad.rawContigTransform=null;}
			}
		}finally{Shared.closeStream(Tadpole.outstream);}
	}

	/** Installs initial-only policy before workers start, never on bridge tables. */
	@Override
	public void prepare(final Tadpole tad, final TadpoleMulti.Config schedule){
		assert(tad!=null && !prepared) : "Each policy owns at most one initial assembly.";
		if(!assembly){return;}
		if(tad.processingMode!=Tadpole.contigMode){
			throw new IllegalStateException("Assembly policy disagrees with Tadpole's resolved mode.");
		}
		prepared=true;
		previousCoverage=BubblePopper.directCoverageGate;
		if(schedule!=null && schedule.automaticFusion
				&& (!tad.tables().rcomp() || AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN)){
			schedule.fuseNet=null;
			Tadpole.outstream.println("Automatic fusion NN disabled: requires canonical DNA kmers.");
		}
		if(initialRatio>0){
			coverage=new BubblePopper.CoverageGate(initialRatio);
			BubblePopper.directCoverageGate=coverage;
		}
		final String unsupported=unsupported(tad);
		if(rescueRequested && unsupported!=null){
			if(rescue==1 || modelSet){throw new IllegalArgumentException("branchrescue: "+unsupported);}
			Tadpole.outstream.println("Branch postpass disabled: "+unsupported);
		}else if(rescueRequested){
			assert(tad.rawContigTransform==null) : "The public postpass must not replace another raw-contig transform.";
			final SameKJoinRescue post=new SameKJoinRescue(tad, postModel, postCutoff);
			tad.rawContigTransform=contigs -> post.mergeLoaded(contigs, SameKJoinFeatures.moments(tad));
			Tadpole.outstream.println("Branch postpass: model="+postModel+" cutoff="+postCutoff+
					"; original assembly table, before ordinary cleanup.");
		}
	}

	/** Reports and releases initial cleanup policy before any later-K operation. */
	@Override
	public void finish(final Tadpole tad){
		assert(tad!=null) : "Finish belongs to the original table owner.";
		if(coverage!=null){
			Tadpole.outstream.println("Initial merge coverage: ratio="+coverage.ratio+
					", evaluations="+coverage.evaluations+", rejected="+coverage.rejected+".");
		}
		close();
	}

	/** Restores caller state on success, validation failure, or worker failure. */
	@Override
	public void close(){
		if(prepared){BubblePopper.directCoverageGate=previousCoverage; prepared=false;}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	/** Carries automatic-model eligibility separately from a custom explicit binding. */
	TadpoleMulti.Config multiConfig(){
		assert(multi) : "Only the multi-K entry point resolves a phase schedule.";
		final TadpoleMulti.Config config=new TadpoleMulti.Config(args);
		config.automaticFusion=bundledFusion;
		return config;
	}

	/** Auto abstains outside the frozen model's domain rather than extrapolating. */
	private static String unsupported(final Tadpole tad){
		assert(tad!=null) : "Model eligibility uses the fully resolved native assembler.";
		if(tad.k()<24 || tad.k()>127){return "the bundled postpass supports assembly K24-127.";}
		if(tad.hashKmerMode!=Tadpole.HASH_EXPLICIT || !tad.tables().rcomp()
				|| AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN){return "requires explicit canonical DNA kmers.";}
		if(!tad.retainShortContigs || tad.trimEnds!=0 || tad.trimCircular || tad.minExtension<2
				|| !tad.useOwnership || tad.minCountExtend<1 || tad.minContigLen<1){
			return "requires retained untrimmed unitigs, minextension>=2, mce>=1, mincontig>=1 and ownership.";
		}
		if(Shared.SIMD_FMA){return "requires simdfma=f for the validated scoring arithmetic.";}
		return null;
	}

	/** Checks output aliases before native constructors can open redirected streams. */
	private static void checkModelOutputs(final String model, final String[] input){
		assert(model!=null) : "Only enabled models need path protection.";
		final ArrayList<String> outputs=new ArrayList<String>();
		for(String arg : input){
			final String name=key(arg), val=value(arg);
			if(val==null || "null".equalsIgnoreCase(val)){continue;}
			if(name.startsWith("out") || name.equals("dot") || name.equals("gfa")
					|| name.equals("dump") || name.equals("hist") || name.equals("khist")
					|| name.equals("ihist") || name.equals("inserthistogram")){
				outputs.addAll(Arrays.asList(val.split(",")));
			}
		}
		// Other writers validate their own collisions. Here only model preservation
		// matters; repeated output options are not themselves a model-path conflict.
		for(String output : outputs){TadpoleGraph.checkPaths(model, new ArrayList<String>(), output);}
	}

	/** Resolves shipped resources by installation path, not the caller's directory. */
	private static String modelPath(final String name){
		assert(name!=null) : "Disabled neural policies never resolve a resource.";
		final String path=Data.findPath(name);
		if(path==null || !Tools.testInputFiles(false, true, path)){
			throw new IllegalArgumentException("Unreadable assembly network: "+name);
		}
		return path;
	}

	/** Mirrors mode-selection triggers; preprocessing flags alone do not imply assembly. */
	static boolean isAssembly(final String[] args){
		assert(args!=null) : "Mode selection needs expanded arguments.";
		String mode=null;
		boolean contigs=false, ecc=false, mark=false, hp=false, local=false, toss=false, uncorrectable=false;
		int left=-1, right=-1, depth=-1;
		for(String arg : args){
			final String name=key(arg), val=value(arg);
			if(name.equals("mode")){mode=val;}
			else if(name.equals("contigs")){contigs=val!=null && !"null".equalsIgnoreCase(val);}
			else if(name.equals("ecc") || name.equals("ecct")){ecc=Parse.parseBoolean(val);}
			else if(name.equals("aecc") || name.equals("aec") || name.equals("aggressive")
					|| name.equals("cecc") || name.equals("cec") || name.equals("conservative")){
				if(Parse.parseBoolean(val)){ecc=true;}
			}else if(name.equals("markerrors")){mark=Parse.parseBoolean(val);}
			else if(name.equals("hpindel")){hp=Parse.parseBoolean(val);}
			else if(name.equals("fixindels") || name.equals("localedit")){local=Parse.parseBoolean(val);}
			else if(name.equals("tossjunk") || name.equals("tj")){toss=Parse.parseBoolean(val);}
			else if(name.equals("tossuncorrectable") || name.equals("tu")){uncorrectable=Parse.parseBoolean(val);}
			else if(name.equals("tossdepth") || name.equals("discarddepth")){depth=Integer.parseInt(val);}
			else if(name.equals("extendleft") || name.equals("el")){left=Parse.parseIntKMG(val);}
			else if(name.equals("extendright") || name.equals("er")){right=Parse.parseIntKMG(val);}
			else if(name.equals("maxextension") || name.equals("maxe")){left=right=Parse.parseIntKMG(val);}
		}
		if(mode!=null){
			// Match Tadpole's numeric parser, including leading zeros and KMG suffixes.
			return mode.equalsIgnoreCase("contig") || (mode.length()>0 && Tools.isDigit(mode.charAt(0))
					&& Parse.parseIntKMG(mode)==Tadpole.contigMode);
		}
		return !contigs && !ecc && !mark && !hp && !local && !toss && !uncorrectable
				&& depth<=0 && left<=0 && right<=0;
	}

	/** BBTools accepts case-insensitive names and strips optional leading hyphens. */
	private static String key(final String arg){
		final int eq=arg.indexOf('=');
		String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
		while(name.startsWith("-")){name=name.substring(1);}
		return name;
	}

	/** Keeps path case and everything after the first equals sign intact. */
	private static String value(final String arg){
		final int eq=arg.indexOf('=');
		return eq<0 ? null : arg.substring(eq+1);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final String[] args;
	final boolean assembly, multi, rescueRequested;
	private boolean defaults=true, modelSet=false, prepared=false, bundledFusion=false;
	private int rescue=-1;
	private String postModel="?tadpole_contig_rescue.bbnet";
	private float postCutoff=0.7f, initialRatio=Float.NaN;
	private BubblePopper.CoverageGate coverage, previousCoverage;
}
