package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import ml.CellNet;
import ml.CellNetParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Opt-in single-K assembly with learned joining of retained raw unitigs.
 * Replaces ordinary post-contig graph cleanup, not the conservative contig walk.
 * Features and inference borrow the original table; no read reload or truth input.
 * @author Fischl
 */
final class SameKJoinRescue {

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Routes partial bindings here too, so missing model/cutoff fails explicitly. */
	static boolean requested(final String[] args){
		assert(args!=null) : "Rescue dispatch follows Tadpole config expansion.";
		for(String arg : args){if(key(arg).startsWith("contigrescue")){return true;}}
		return false;
	}

	/** Enters through tadpole.sh; the ordinary constructor still owns native options. */
	static void run(final String[] args){new SameKJoinRescue(args).process();}

	/** Development-only composition; the ordinary tadpole.sh binding remains unchanged. */
	static void runCombined(final String[] args, final String liveModel, final float liveCutoff){
		if(liveModel==null || liveModel.isEmpty() || !Float.isFinite(liveCutoff) || liveCutoff<0 || liveCutoff>1){
			throw new IllegalArgumentException("Combined rescue requires a live model and explicit cutoff in[0,1].");
		}
		final SameKJoinRescue rescue=new SameKJoinRescue(args);
		// The existing constructor checks every other input. This extra binding must
		// also be protected before the postpass can open output or score-trace files.
		TadpoleGraph.checkPaths(liveModel, new ArrayList<String>(Arrays.asList(liveModel)), rescue.output, rescue.scores);
		if(!Tools.testInputFiles(true, true, liveModel)){
			throw new IllegalArgumentException("The live branch model must be readable.");
		}
		rescue.process(new LiveBranchGate(rescue.tad, liveModel, liveCutoff));
	}

	/** Handles standard optional leading hyphens without changing filename case. */
	private static String key(final String arg){
		final int eq=arg.indexOf('=');
		String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
		while(name.startsWith("-")){name=name.substring(1);}
		return name;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prepended defaults permit explicit overrides, which validation then checks. */
	private SameKJoinRescue(final String[] args){
		final ArrayList<String> nativeArgs=new ArrayList<String>(Arrays.asList(
				"mode=contig", "processcontigs=f", "pop=f", "sweeplen=0",
				"retainshortcontigs=t", "trimends=0", "trimcircular=f", "hashkmers=f"));
		for(String arg : args){
			final String name=key(arg);
			final int eq=arg.indexOf('=');
			final String raw=eq<0 ? null : arg.substring(eq+1);
			final String value=raw!=null && raw.equalsIgnoreCase("null") ? null : raw;
			if(name.equals("contigrescuenet")){model=value;}
			else if(name.equals("contigrescuecutoff")){cutoff=Float.parseFloat(value);}
			else if(name.equals("contigrescuescores")){scores=value;}
			else if(name.equals("out") || name.equals("out1") || name.equals("oute") || name.equals("oute1")){
				output=value;
			}else if(name.startsWith("out") || name.equals("dot") || name.equals("gfa") || name.equals("hist")
					|| name.equals("khist") || name.equals("dump") || name.equals("ihist") || name.equals("inserthistogram")){
				// In particular, don't let another PreParser open/truncate an unchecked log path.
				throw new IllegalArgumentException("Contig rescue supports out= and contigrescuescores= only; not "+name+". Redirect stderr in the shell.");
			}else{nativeArgs.add(arg);}
		}
		if(model==null || model.isEmpty() || output==null || output.isEmpty() || output.indexOf(',')>=0
				|| !Float.isFinite(cutoff) || cutoff<0 || cutoff>1){
			throw new IllegalArgumentException("Contig rescue requires contigrescuenet=, explicit contigrescuecutoff in[0,1], and one out=FASTA.");
		}
		if(!FileFormat.isFastaExt(ReadWrite.rawExtension(output))){
			throw new IllegalArgumentException("Contig rescue output must use a FASTA extension (including stdout.fa).");
		}
		if(TadpoleMulti.hasMultipleK(args) || TadpoleGraph.requested(args)){
			throw new IllegalArgumentException("Contig rescue requires single-K contig assembly, not a multi-K schedule or supplied graph.");
		}
		nativeArgs.add("out=null");
		tad=Tadpole.makeTadpole(nativeArgs.toArray(new String[0]), true);
		k=tad.k();
		if(k<24 || k>127 || tad.processingMode!=Tadpole.contigMode || tad.hashKmerMode!=Tadpole.HASH_EXPLICIT
				|| !tad.tables().rcomp() || AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN){
			throw new IllegalArgumentException("Contig rescue requires single-K24–127 explicit canonical DNA contig assembly.");
		}
		if(tad.processContigs || tad.evictLowDepthContigs || !tad.retainShortContigs || tad.trimEnds!=0
				|| tad.trimCircular || tad.minExtension<2 || !tad.useOwnership || tad.minContigLen<1){
			throw new IllegalArgumentException("Contig rescue requires raw unitigs: processcontigs=f pop=f sweeplen=0 retainshortcontigs=t trimends=0 trimcircular=f minextension>=2 ownership=t.");
		}
		if(tad.ecc || tad.ecco || tad.merge || tad.markErrors || tad.hpIndel || tad.localEdit
				|| tad.extendLeft>0 || tad.extendRight>0 || tad.tossJunk || tad.discardLowDepthReads>=0
				|| tad.discardUncorrectable){
			throw new IllegalArgumentException("Contig rescue does not run read correction, extension, merging or discard modes; preprocess separately.");
		}
		if(tad.contigIDOffset!=0 || tad.taxID>=0){
			throw new IllegalArgumentException("Contig rescue writes samek_ IDs and depth headers; contig ID offsets and taxID headers are unsupported.");
		}
		if(Tadpole.append){throw new IllegalArgumentException("Contig rescue does not support append=t.");}
		final ArrayList<String> inputs=new ArrayList<String>(tad.tables().in1);
		inputs.addAll(tad.tables().in2); inputs.addAll(tad.tables().extra); inputs.add(model);
		TadpoleGraph.checkPaths(model, inputs, output, scores);
		if(!Tools.testInputFiles(true, true, inputs.toArray(new String[0]))
				|| !Tools.testOutputFiles(Tadpole.overwrite, false, false, output, scores)){
			throw new IllegalArgumentException("Rescue inputs must be readable and outputs writable without input/output aliases.");
		}
		net=CellNetParser.load(model, false);
		if(net==null || net.numInputs()!=SameKJoinFeatures.INPUTS || net.numOutputs()!=1){
			throw new IllegalArgumentException("Contig rescue requires a27-input, one-output samek27-v1 model.");
		}
		net.simdFF=false;// Matches validation's layer setting; dot products still follow simdfma.
		// Vector.feedForwardDense falls back to Vector.fma even here. For bit-stable
		// scalar dot products use simdfma=f; default SIMD can vary by an ULP with JIT state.
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows an already constructed assembler; never loads reads or owns its table. */
	SameKJoinRescue(final Tadpole tad_, final String model_, final float cutoff_){
		if(tad_==null || model_==null || !Float.isFinite(cutoff_) || cutoff_<0 || cutoff_>1){
			throw new IllegalArgumentException("Borrowed postpass requires a table owner, model and cutoff in[0,1].");
		}
		tad=tad_; k=tad.k(); model=model_; cutoff=cutoff_;
		net=CellNetParser.load(model, false);
		if(net==null || net.numInputs()!=SameKJoinFeatures.INPUTS || net.numOutputs()!=1){
			throw new IllegalArgumentException("Borrowed postpass requires a27-input, one-output samek27-v1 model.");
		}
		net.simdFF=false;
	}

	/** Scores before ordinary cleanup and returns every product, including short contigs. */
	ArrayList<Contig> mergeLoaded(final ArrayList<Contig> contigs, final SameKJoinFeatures.Moments moments){
		assert(output==null && scores==null && contigs!=null && moments!=null) :
				"Borrowed postpass has no writer; its caller owns final filtering and table lifetime.";
		final SameKJoinMerge.Choices choices=score(contigs, moments);
		final SameKJoinMerge.Result result=SameKJoinMerge.merge(contigs, choices.mates(), k);
		final Kmer word=new Kmer(k);
		for(Contig c : result.contigs){
			if(tad instanceof Tadpole1){((Tadpole1)tad).tables.calcCoverage(c);}
			else{((Tadpole2)tad).tables.calcCoverage(c, word);}
			if(c.minCov<1 || !Float.isFinite(c.coverage)){
				throw new IllegalStateException("Exact K-1 joining must preserve present kmers before multi-K traversal.");
			}
			c.name=null;// Restore native depth/end annotations for later graph emission.
			if(!c.canonical()){
				// Contig.rcomp swaps native endpoint state but not graph exploration limits.
				final boolean limited=c.graphLeftLimited;
				c.graphLeftLimited=c.graphRightLimited; c.graphRightLimited=limited;
				c.rcomp();
			}
		}
		Tadpole.outstream.println("MULTIK_POSTPASS input="+contigs.size()+" candidates="+candidates+
				" selected="+choices.supported+" joined="+result.joins+" cycle_unitigs="+result.cycleUnitigs+
				" output="+result.contigs.size());
		return result.contigs;
	}

	/** Keeps the table live through scoring and metadata, clearing it even after failure. */
	private void process(){process(null);}

	/** Runs both decisions against one table, sharing immutable moments when a live gate exists. */
	private void process(final LiveBranchGate liveGate){
		assert(tad!=null && net!=null) : "Validate the complete binding before loading reads or opening outputs.";
		final long start=System.nanoTime();
		Tadpole.outstream.println("Contig rescue: raw-unitig postpass replaces ordinary graph cleanup; model="+model+" cutoff="+cutoff);
		tad.printExecutionPlan();
		try{
			assert(tad.liveBranchGate==null) : "This driver owns the live gate for exactly one raw assembly.";
			tad.liveBranchGate=liveGate;
			tad.process2(Tadpole.contigMode);
			if(tad.errorState || tad.tables().errorState){throw new IllegalStateException("Contig assembly failed before rescue.");}
			final ArrayList<Contig> contigs=tad.detachContigs();
			final long momentStart=System.nanoTime();
			final SameKJoinFeatures.Moments moments=liveGate==null ? SameKJoinFeatures.moments(tad) : liveGate.populationMoments();
			final long scoreStart=System.nanoTime();
			final SameKJoinMerge.Choices choices=score(contigs, moments);
			final long mergeStart=System.nanoTime();
			final SameKJoinMerge.Result result=SameKJoinMerge.merge(contigs, choices.mates(), k);
			final int emitted=writeContigs(result.contigs);
			final long end=System.nanoTime();
			Tadpole.outstream.println("SAME_K_RESCUE_PASS input="+contigs.size()+" candidates="+candidates+
					" selected="+choices.supported+" joined="+result.joins+" cycle_unitigs="+result.cycleUnitigs+
					" output="+result.contigs.size()+" emitted="+emitted+" mincontig="+tad.minContigLen);
			Tadpole.outstream.println("SAME_K_RESCUE_SECONDS build="+((momentStart-start)*1e-9)+
					" moments="+((scoreStart-momentStart)*1e-9)+" score="+((mergeStart-scoreStart)*1e-9)+
					" merge_write="+((end-mergeStart)*1e-9)+" total="+((end-start)*1e-9));
			if(liveGate!=null){
				final long[] totals=liveGate.totals();
				Tadpole.outstream.println("LIVE_BRANCH_RESCUE_PASS k="+k+" evaluations="+totals[0]+" rescued_encounters="+totals[1]);
				Tadpole.outstream.println("LIVE_BRANCH_MOMENTS seconds="+(liveGate.momentNanos()*1e-9));
			}
		}finally{
			tad.liveBranchGate=null;
			tad.clearData();
		}
		if(tad.showStats && !FileFormat.isStdio(output)){
			jgi.AssemblyStats2.main(new String[]{"in="+output, "printextended"});
		}
	}

	/** Scores all physical candidates before selecting reciprocal unique partners. */
	private SameKJoinMerge.Choices score(final ArrayList<Contig> contigs, final SameKJoinFeatures.Moments moments){
		assert(contigs!=null) : "Native contig mode must finish before endpoint enumeration.";
		final SameKJoinMerge.Choices choices=new SameKJoinMerge.Choices(contigs.size());
		final ByteStreamWriter trace=scores==null ? null : new ByteStreamWriter(scores, Tadpole.overwrite, false, false);
		final float[] vector=new float[SameKJoinFeatures.INPUTS];
		final ByteBuilder row=new ByteBuilder();
		final TadpoleGraph.Counts counts=new TadpoleGraph.Counts(){
			@Override public int count(final Kmer word){return tad.bridgeCount(word);}
		};
		if(trace!=null){trace.start(); trace.print("event\tk\tscore_bits\tscore\taccepted\n");}
		try{
			candidates=new SameKJoinCollector(contigs, counts, k).collect(new SameKJoinCollector.Sink(){
				@Override
				public void accept(final long event, final int se, final int de, final Contig a, final boolean ar,
						final Contig b, final boolean br, final int[] ac, final int[] bc){
					if(a.length()<k+2 || b.length()<k+2){
						throw new IllegalStateException("Rescue features require formed K+2 unitigs; preserve minextension=2.");
					}
					final int ai=AminoAcid.baseToNumber[SameKJoinCollector.baseAt(a, ar, a.length()-k-1)];
					final int bi=AminoAcid.baseToNumber[SameKJoinCollector.baseAt(b, br, k)];
					SameKJoinFeatures.fill(vector, k, ac, bc, ai, bi, a.coverage, b.coverage,
							a.length()-k+1, b.length()-k+1, moments.uniqueMean(), moments.observationMean());
					final float probability=net.applyInput(vector).feedForward();
					if(!Float.isFinite(probability) || probability<0 || probability>1){
						throw new IllegalStateException("Model output is not a finite probability: "+probability);
					}
					final boolean accepted=probability>=cutoff;
					if(accepted){choices.add(se, de);}
					if(trace!=null){
						trace.print(row.clear().append(event).tab().append(k).tab().append(Float.floatToRawIntBits(probability))
								.tab().appendSlow(probability).tab().append(accepted ? 1 : 0).nl());
					}
				}
			});
		}finally{if(trace!=null && trace.poisonAndWait()){throw new IllegalStateException("Rescue score trace failed.");}}
		return choices;
	}

	/** Applies final mincontig only now; depths are recomputed from the borrowed table. */
	private int writeContigs(final ArrayList<Contig> contigs){
		assert(tad.minContigLen>0) : "Only the final emission stage applies mincontig to rescued paths.";
		final Kmer word=new Kmer(k);
		final ByteBuilder row=new ByteBuilder();
		final ByteStreamWriter writer=new ByteStreamWriter(output, Tadpole.overwrite, false, true);
		writer.start();
		int emitted=0;
		try{
			for(Contig c : contigs){
				if(c.length()<tad.minContigLen){continue;}
				if(tad instanceof Tadpole1){((Tadpole1)tad).tables.calcCoverage(c);}
				else{((Tadpole2)tad).tables.calcCoverage(c, word);}
				if(c.minCov<1 || !Float.isFinite(c.coverage)){
					throw new IllegalStateException("Exact K-1 joining must preserve present kmers; missing count in "+c.name);
				}
				// These are depth annotations, not invented native branch/end classifications.
				writer.print(row.clear().append('>').append(c.name).append(" length=").append(c.length())
						.append(" cov=").append(c.coverage, 6).append(" min=").append(c.minCov)
						.append(" max=").append(c.maxCov).nl().append(c.bases).nl());
				emitted++;
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Rescued assembly output failed.");}}
		return emitted;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Tadpole tad;
	private final CellNet net;
	private final int k;
	private String model, output, scores;
	private float cutoff=Float.NaN;
	private long candidates;
}
