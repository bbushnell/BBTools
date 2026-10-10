package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Runs assembly followed by ordered, shared-table fusion and bridging phases.
 * The last phase lends its count table to final graph processing when K matches.
 *
 * @author Brian Bushnell, Noelle
 */
public class TadpoleMulti {

	public static void main(String[] args){
		args=Tadpole.expandConfigArgs(args);
		final Config config=new Config(args);
		new TadpoleMulti(config).process();
	}

	TadpoleMulti(Config config_){config=config_;}

	void process(){
		final Tadpole longest=Tadpole.makeTadpole(makeArgs(config.assembleK, true), true);
		fusionNeural=makeNeuralGate(longest);
		try(FusionJoinCollector collector=makeCollector(longest)){
			fusionCollector=collector;
			processInner(longest);
		}finally{
			fusionCollector=null;
			clearFusionSupport();
		}
	}

	/** Validates model/output aliases before opening any collector or assembly writer. */
	private FusionNeuralGate makeNeuralGate(final Tadpole longest){
		if(config.fuseNet==null){return null;}
		final ArrayList<String> inputs=new ArrayList<String>(longest.tables().in1);
		inputs.addAll(longest.tables().in2);
		inputs.addAll(longest.tables().extra);
		inputs.add(config.fuseNet);
		TadpoleGraph.checkPaths(config.fuseNet, inputs, config.out, config.outGfa);
		if(!Tools.testInputFiles(false, true, config.fuseNet)){
			throw new IllegalArgumentException("Unreadable fusion network: "+config.fuseNet);
		}
		return new FusionNeuralGate(config.fuseNet, config.fuseCutoff);
	}

	/** Creates optional evidence sidecars, protecting all source files and final graph outputs. */
	private FusionJoinCollector makeCollector(final Tadpole longest){
		if(config.fusionVectorPrefix==null){return null;}
		final ArrayList<String> inputs=new ArrayList<String>(longest.tables().in1);
		inputs.addAll(longest.tables().in2);
		inputs.addAll(longest.tables().extra);
		if(config.fuseNet!=null){inputs.add(config.fuseNet);}
		return new FusionJoinCollector(config.fusionVectorPrefix, inputs, config.out, config.outGfa);
	}

	/** Constructs initial contigs, then applies the selected multi-K schedule. */
	private void processInner(final Tadpole longest){
		assert(longest!=null) : "The initial assembler owns source arguments and original contigs.";
		config.printExecutionPlan(longest);
		if(!Tools.testOutputFiles(Tadpole.overwrite, Tadpole.append, false, config.out)){
			throw new RuntimeException("Can't write output file "+config.out+"; overwrite="+Tadpole.overwrite);
		}
		longest.process2(Tadpole.contigMode);
		checkErrorState(longest);
		final EarlyLowDepthTracker earlyLowDepth;
		if(config.earlyLowDepthDiag()){
			config.applyLowDepthThresholds(longest);
			earlyLowDepth=new EarlyLowDepthTracker(longest.diagnoseLowDepthContigs("assemble-k"));
		}else{earlyLowDepth=null;}
		longest.markBridgeEndpoints();
		longest.clearContigEdges();
		ArrayList<Contig> contigs=longest.detachContigs();
		final int minContig=longest.minContigLen;
		final int idOffset=longest.contigIDOffset;
		if(config.ordered){
			contigs=processOrdered(contigs, longest, minContig);
		}else{
			longest.tables().clear();
			System.gc();
			contigs=processLegacy(contigs, longest, minContig);
		}
		if(earlyLowDepth!=null && finalLowDepthDiagnostic!=null){
			earlyLowDepth.reportFates(finalLowDepthDiagnostic);
		}
		writeContigs(contigs, config.out, minContig, idOffset);
		if(config.showStats && FileFormat.isFastaExt(ReadWrite.rawExtension(config.out)) && !FileFormat.isStdio(config.out)){
			System.err.println();
			jgi.AssemblyStats2.main(new String[] {"in="+config.out, "printextended"});
		}
	}

	/** Keeps the historical ordering available for controlled comparisons. */
	private ArrayList<Contig> processLegacy(ArrayList<Contig> contigs, final Tadpole longest, final int minContig){

		/* Preserve fusion-before-bridging order. Each fusion K owns at most one table. */
		Tadpole reusableGraphTadpole=null;
		for(int i=0; i<config.fuseKs.length; i++){
			final int k=config.fuseKs[i], before=contigs.size();
			try{
				loadFusionSupport(k);
				final Timer timer=new Timer();
				longest.setContigs(contigs);
				longest.clearContigEdges();
				final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, k,
						config.assembleK-1, false, minContig, config.fuseMaxMismatches,
						config.fuseDeadEndsOnly, config.fuseConflicts);
				overlapper.support=fusionSupport;
				overlapper.collector=fusionCollector;
				overlapper.neural=fusionNeural;
				overlapper.allowTrim=config.fuseTrim;
				overlapper.maxCoverageRatio=config.fuseCoverageRatio;
				if(overlapper.addEdges()>0){mergeCrossK(longest);}
				contigs=longest.detachContigs();
				checkErrorState(longest);
				timer.stop();
				System.err.println("Cross-k overlaps "+k+": "+before+" -> "+contigs.size()+" contigs; "+timer);
			}finally{clearFusionSupport();}
		}

		/* Bridge tables are needed only for unbranched paths across actual sequence gaps. */
		for(int i=0; i<config.bridgeKs.length; i++){
			final int k=config.bridgeKs[i], before=contigs.size();
			final int endpoints=countBridgeEndpoints(contigs, k);
			if(endpoints<2){
				System.err.println("Cross-k bridges "+k+": skipped; only "+endpoints+" eligible endpoint"+
						(endpoints==1 ? "." : "s."));
				continue;
			}
			final Tadpole tad=Tadpole.makeTadpole(makeArgs(k, false), true);
			setCrossKGraph(tad, true);
			tad.loadKmers(new Timer());
			tad.pruneLoadedKmers();
			tad.setContigs(contigs);
			tad.cleanLoadedKmers(contigs);
			checkErrorState(tad);
			tad.clearContigEdges();
			final boolean resolveRepeats=tad.resolveRepeats;
			tad.resolveRepeats=false;
			try{
				tad.processContigs();
			}finally{
				tad.resolveRepeats=resolveRepeats;
			}
			mergeCrossK(tad);

			contigs=tad.detachContigs();
			checkErrorState(tad);
			if(config.finalGraphNeeded() && k==config.graphK && i==config.bridgeKs.length-1){
				reusableGraphTadpole=tad;
			}else{tad.tables().clear();}
			System.gc();
			System.err.println("Cross-k bridges "+k+": "+before+" -> "+contigs.size()+" contigs.");
		}

		if(config.finalGraphNeeded()){
			contigs=extractFinalGraph(contigs, reusableGraphTadpole, minContig);
		}
		return contigs;
	}

	/** Uses one table per occurrence; fusion sees raw counts before bridge cleaning. */
	private ArrayList<Contig> processOrdered(ArrayList<Contig> contigs, final Tadpole initial, final int minContig){
		assert(config.phaseKs.length>0 && config.phaseKs[0]==config.assembleK) :
				"The initial assembly table must be the first ordered phase owner.";
		boolean finalDone=false;
		for(int i=0; i<config.phaseKs.length; i++){
			final int k=config.phaseKs[i];
			final boolean fuse=Config.contains(config.fuseKs, k);
			final boolean bridge=Config.contains(config.bridgeKs, k) && (i>0 || config.bridgeInitial);
			final boolean finish=config.finalGraphNeeded() && i==config.phaseKs.length-1 && k==config.graphK;
			final boolean evidence=fuse && (config.fusePath || fusionCollector!=null || fusionNeural!=null);
			Tadpole owner=(i==0 ? initial : null);
			try{
				if(owner==null && (evidence || finish || (bridge && countBridgeEndpoints(contigs, k)>=2))){
					owner=Tadpole.makeTadpole(evidence ? makeFusionSupportArgs(k) : makeArgs(k, false), true);
					System.err.println("Loading ordered phase "+(i+1)+" at k="+k+
							(evidence ? " (explicit fusion evidence shared with bridging)." : "."));
					owner.loadKmers(new Timer());
					checkErrorState(owner);
				}
				if(fuse){
					if(evidence){
						if(config.fusePath){attachFusionSupport(owner);}
						if(fusionCollector!=null){fusionCollector.beginPhase(k, config.assembleK-1, false, countProvider(owner));}
						if(fusionNeural!=null){fusionNeural.beginPhase(k, countProvider(owner));}
					}
					try{contigs=fuseContigs(contigs, initial, k, minContig);}
					finally{clearFusionSupport();}// Borrower only; owner survives until this phase ends.
				}
				if(owner!=null && i>0){
					setCrossKGraph(owner, bridge);
					owner.pruneLoadedKmers();
					owner.setContigs(contigs);
					owner.cleanLoadedKmers(contigs);
					checkErrorState(owner);
				}
				if(bridge){
					final int endpoints=countBridgeEndpoints(contigs, k);
					if(endpoints>=2){
						assert(owner!=null) : "Bridge traversal requires the current phase's live count table.";
						contigs=bridgeContigs(contigs, owner);
					}else{System.err.println("Cross-k bridges "+k+": skipped; only "+endpoints+" eligible endpoints.");}
				}
				if(finish){
					contigs=extractFinalGraph(contigs, owner, minContig);
					owner=null;// extractFinalGraph released the table successfully.
					finalDone=true;
				}
			}finally{
				clearFusionSupport();
				if(owner!=null){owner.tables().clear();}
			}
			System.gc();
		}
		if(config.finalGraphNeeded() && !finalDone){contigs=extractFinalGraph(contigs, null, minContig);}
		return contigs;
	}

	/** Adds shorter-K overlap proposals using the current optional borrowed evidence. */
	private ArrayList<Contig> fuseContigs(final ArrayList<Contig> contigs, final Tadpole initial,
			final int k, final int minContig){
		assert(k<config.assembleK) : "CrossKTipOverlapper's overlap ceiling remains assembleK-1.";
		final int before=contigs.size();
		final Timer timer=new Timer();
		initial.setContigs(contigs);
		initial.clearContigEdges();
		final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, k, config.assembleK-1,
				false, minContig, config.fuseMaxMismatches, config.fuseDeadEndsOnly, config.fuseConflicts);
		overlapper.support=fusionSupport;
		overlapper.collector=fusionCollector;
		overlapper.neural=fusionNeural;
		overlapper.allowTrim=config.fuseTrim;
		overlapper.maxCoverageRatio=config.fuseCoverageRatio;
		if(overlapper.addEdges()>0){mergeCrossK(initial);}
		final ArrayList<Contig> result=initial.detachContigs();
		checkErrorState(initial);
		timer.stop();
		System.err.println("Cross-k overlaps "+k+": "+before+" -> "+result.size()+" contigs; "+timer);
		return result;
	}

	/** Traverses gaps on the cleaned phase table without enabling final graph operations. */
	private ArrayList<Contig> bridgeContigs(final ArrayList<Contig> contigs, final Tadpole tad){
		assert(tad!=null) : "Bridge graph discovery must borrow a live phase table.";
		final int before=contigs.size();
		setCrossKGraph(tad, true);
		tad.setContigs(contigs);
		tad.clearContigEdges();
		final boolean oldResolve=tad.resolveRepeats, oldPop=tad.popBubbles;
		tad.resolveRepeats=false;
		tad.popBubbles=false;
		try{tad.processContigs();}
		finally{tad.resolveRepeats=oldResolve; tad.popBubbles=oldPop;}
		mergeCrossK(tad);
		final ArrayList<Contig> result=tad.detachContigs();
		checkErrorState(tad);
		System.err.println("Cross-k bridges "+tad.k()+": "+before+" -> "+result.size()+" contigs.");
		return result;
	}

	/** Builds one complete graph after all cross-k merges, then emits the requested terminal representation. */
	private ArrayList<Contig> extractFinalGraph(final ArrayList<Contig> contigs,
			Tadpole tad, final int minContig){
		if(tad==null){
			System.err.println("Loading graph-k table at k="+config.graphK+".");
			tad=Tadpole.makeTadpole(makeArgs(config.graphK, false), true);
			tad.loadKmers(new Timer());
			tad.pruneLoadedKmers();
			tad.setContigs(contigs);
			tad.cleanLoadedKmers(contigs);
			checkErrorState(tad);
		}else{
			System.err.println("Reusing bridge-k table for final graph at k="+config.graphK+".");
		}
		setCrossKGraph(tad, false);
		tad.minContigLen=minContig;
		tad.refreshGraphEndpoints=true;
		tad.popBubbles=false;
		// The borrowed initial table must not repeat its initial-only sweep.
		tad.sweepContigLen=0;
		tad.classifyGraphContigs=false;
		/* Refresh graph-k endpoint topology before resolving exact overlaps that were
		 * ineligible under the original longest-k endpoint classifications. */
		final boolean resolveRepeats=tad.resolveRepeats;
		tad.resolveRepeats=false;
		tad.simpleOmnitigs=false;
		tad.graphCover=false;
		tad.lowDepthContigDiag=false;
		tad.setContigs(contigs);
		tad.clearContigEdges();
		tad.processContigs();
		tad.clearContigEdges();
		int graphOverlapBefore=-1;
		final Timer graphOverlapTimer=new Timer();
		if(config.graphK<config.assembleK){
			//Borrow the current graph counts, which may already be pruned/washed.
			//Never load another table while this graph table remains live.
			if(config.fusePath){attachFusionSupport(tad);}
			if(fusionCollector!=null){fusionCollector.beginPhase(tad.k(), config.assembleK-1, true, countProvider(tad));}
			if(fusionNeural!=null){fusionNeural.beginPhase(tad.k(), countProvider(tad));}
			graphOverlapBefore=contigs.size();
			try{
				final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, config.graphK,
						config.assembleK-1, true, minContig, config.fuseMaxMismatches,
						config.fuseDeadEndsOnly, config.fuseConflicts);
				overlapper.support=fusionSupport;
				overlapper.collector=fusionCollector;
				overlapper.neural=fusionNeural;
				overlapper.allowTrim=config.fuseTrim;
				overlapper.maxCoverageRatio=config.fuseCoverageRatio;
				if(overlapper.addEdges()>0){mergeCrossK(tad);}
			}finally{clearFusionSupport();}
		}
		final ArrayList<Contig> merged=tad.detachContigs();
		graphOverlapTimer.stop();
		if(graphOverlapBefore>=0){
			System.err.println("Graph-k overlaps "+config.graphK+": "+graphOverlapBefore+
					" -> "+merged.size()+" contigs; "+graphOverlapTimer);
		}
		tad.resolveRepeats=resolveRepeats;
		tad.simpleOmnitigs=config.simpleOmnitigs;
		tad.graphCover=config.graphCover;
		config.applyLowDepthDiagnostic(tad);
		config.applyFinalGraphClassification(tad);
		config.applyFinalGraphOutput(tad);
		tad.setContigs(merged);
		tad.clearContigEdges();
		final BubblePopper.CoverageGate previousGate=BubblePopper.directCoverageGate;
		final BubblePopper.CoverageGate gate=config.graphMergeCoverageRatio==0 ? null
				: new BubblePopper.CoverageGate(config.graphMergeCoverageRatio);
		BubblePopper.directCoverageGate=gate;
		try{
			tad.processContigs();
		}finally{
			BubblePopper.directCoverageGate=previousGate;
			if(gate!=null){
				System.err.println("Graph merge coverage: ratio="+gate.ratio+", evaluations="+
						gate.evaluations+", rejected="+gate.rejected+".");
			}
		}
		finalLowDepthDiagnostic=tad.lastLowDepthDiagnostic;
		final ArrayList<Contig> extracted=tad.detachContigs();
		checkErrorState(tad);
		tad.tables().clear();
		System.gc();
		return extracted;
	}

	private static void setCrossKGraph(final Tadpole tad, final boolean value){
		if(tad instanceof Tadpole1){((Tadpole1)tad).crossKGraph=value;}
		else{((Tadpole2)tad).crossKGraph=value;}
	}

	private static int countBridgeEndpoints(ArrayList<Contig> contigs, final int k){
		int count=0;
		for(Contig c : contigs){
			if(c.length()<k){continue;}
			if(c.leftBridgeEndpoint){count++;}
			if(c.rightBridgeEndpoint){count++;}
		}
		return count;
	}

	private void mergeCrossK(Tadpole tad){
		final boolean oldDirect=BubblePopper.popDirect;
		final boolean oldIndirect=BubblePopper.popIndirect;
		final boolean oldCrossK=BubblePopper.crossKMerge;
		final float oldDepthRatio=BubblePopper.crossKMaxDepthRatio;
		final int oldMismatches=BubblePopper.crossKMaxMismatches;
		final FusionKmerSupport oldSupport=BubblePopper.crossKSupport;
		BubblePopper.popDirect=true;
		BubblePopper.popIndirect=false;
		BubblePopper.crossKMerge=true;
		BubblePopper.crossKMaxDepthRatio=config.maxDepthRatio;
		BubblePopper.crossKMaxMismatches=config.fuseMaxMismatches;
		BubblePopper.crossKSupport=fusionSupport;
		try{
			for(int pass=0, merged=1; pass<config.passes && merged>0; pass++){
				merged=tad.popBubbles(false);
			}
		}finally{
			BubblePopper.popDirect=oldDirect;
			BubblePopper.popIndirect=oldIndirect;
			BubblePopper.crossKMerge=oldCrossK;
			BubblePopper.crossKMaxDepthRatio=oldDepthRatio;
			BubblePopper.crossKMaxMismatches=oldMismatches;
			BubblePopper.crossKSupport=oldSupport;
		}
	}

	/** Loads the sole live table for an initial same-K fusion pass. */
	private void loadFusionSupport(final int k){
		if(!config.fusePath && fusionCollector==null && fusionNeural==null){return;}
		assert(fusionSupport==null && fusionSupportTadpole==null) :
				"Each serial fusion phase must release its support table before another load.";
		final Tadpole evidence=Tadpole.makeTadpole(makeFusionSupportArgs(k), true);
		fusionSupportTadpole=evidence;
		if(evidence.k()!=k){
			throw new IllegalArgumentException("Fusion support K changed while constructing its count table.");
		}
		System.err.println("Loading fusion path at k="+k+", depth="+config.fusePathDepth+".");
		evidence.loadKmers(new Timer());
		checkErrorState(evidence);
		if(config.fusePath){attachFusionSupport(evidence);}
		if(fusionCollector!=null){fusionCollector.beginPhase(k, config.assembleK-1, false, countProvider(evidence));}
		if(fusionNeural!=null){fusionNeural.beginPhase(k, countProvider(evidence));}
	}

	/** Borrows counts only; feature collection must not implicitly enable a path veto. */
	private static TadpoleGraph.Counts countProvider(final Tadpole evidence){
		assert(evidence!=null) : "Fusion features require a live-phase count owner.";
		if(!evidence.tables().rcomp() || kmer.AbstractKmerTableSet.MASK_MIDDLE || shared.Shared.AMINO_IN){
			throw new IllegalArgumentException("Fusion features require rcomp=t, maskmiddle=f, amino=f.");
		}
		return new TadpoleGraph.Counts(){
			@Override
			public int count(final Kmer word){return evidence.bridgeCount(word);}
		};
	}

	/** Installs a serial checker; table ownership is separate so final graphs can lend theirs. */
	private void attachFusionSupport(final Tadpole evidence){
		assert(fusionSupport==null) : "A fusion checker must be detached before changing its count table.";
		if(!evidence.tables().rcomp() || kmer.AbstractKmerTableSet.MASK_MIDDLE || shared.Shared.AMINO_IN){
			throw new IllegalArgumentException("Fusion paths require rcomp=t, maskmiddle=f, amino=f.");
		}
		fusionSupport=new FusionKmerSupport(new FusionKmerSupport.Evidence(){
			@Override
			public int count(final Kmer key){return evidence.bridgeCount(key);}
			@Override
			public boolean isJunction(final int max, final int second){return evidence.isJunction(max, second);}
		}, evidence.k(), config.fusePathDepth, config.fusePathFlank);
	}

	/** Evidence counts must not lose singleton words to a probabilistic prefilter. */
	String[] makeFusionSupportArgs(final int k){
		assert((config.fusePath || config.fusionVectorPrefix!=null || config.fuseNet!=null) && k<config.assembleK) :
				"Only a shorter-K support or collection pass may load a separate, sole evidence table.";
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(makeArgs(k, false)));
		args.add("hashkmers=explicit");
		args.add("prefilter=f");
		args.add("prepasses=1");
		return args.toArray(new String[args.size()]);
	}

	/** Releases an owned table, never a borrowed final-graph table, even after a failed load. */
	private void clearFusionSupport(){
		if(fusionCollector!=null){fusionCollector.endPhase();}
		if(fusionNeural!=null){fusionNeural.endPhase();}
		if(fusionSupport!=null){
			System.err.println("Fusion path: k="+fusionSupport.k+", evaluations="+fusionSupport.evaluations+
					", rejected="+fusionSupport.rejected+", words="+fusionSupport.words+
					", flank="+fusionSupport.flank+".");
		}
		fusionSupport=null;
		if(fusionSupportTadpole!=null){
			fusionSupportTadpole.tables().clear();
			fusionSupportTadpole=null;
			System.gc();
		}
	}

	private static void checkErrorState(Tadpole tad){
		if(tad.errorState || tad.tables().errorState){
			throw new RuntimeException(tad.getClass().getSimpleName()+" terminated in an error state; no assembly will be written.");
		}
	}

	String[] makeArgs(final int k, final boolean longest){
		final ArrayList<String> list=new ArrayList<String>(config.common.size()+6);
		list.addAll(config.common);
		final int hashMode=(longest ? Tadpole.HASH_EXPLICIT : config.hashMode(k));
		list.add("hashkmers="+Tadpole.hashModeName(hashMode));
		list.add("seedfromtable="+longest);
		if(!longest){list.add("tipseededwash=t");}
		list.add("k="+k);
		list.add("out=null");
		list.add("mode=contig");
		list.add("showstats=f");
		list.add("processcontigs="+(longest ? "t" : "f"));
		/* Sweep once after initial contig construction; bridge tables never discard contigs. */
		list.add("sweeplen="+(longest ? config.sweepContigLen : 0));
		if(longest && config.evictLowDepthContigs){
			list.add("evictlowdepthcontigs=t");
			config.addLowDepthThresholdArgs(list);
		}
		if(longest && config.retainShortContigsSet){list.add("retainshortcontigs="+config.retainShortContigs);}
		else if(longest && config.graphClassificationRequested()){list.add("retainshortcontigs=t");}
		if(!longest){list.add("pop=f");}
		return list.toArray(new String[list.size()]);
	}

	private void writeContigs(ArrayList<Contig> contigs, String out, int minContig, int idOffset){
		if(out==null){return;}
		if(!Tools.testOutputFiles(Tadpole.overwrite, Tadpole.append, false, out)){
			throw new RuntimeException("Can't write output file "+out+"; overwrite="+Tadpole.overwrite);
		}
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.FA, 0, 0, true,
				Tadpole.overwrite, Tadpole.append, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		int id=idOffset;
		for(Contig c : contigs){
			c.id=id++;
			if(config.emitContig(c, minContig)){bsw.println(c);}
		}
		if(bsw.poisonAndWait()){throw new RuntimeException("Error writing "+out);}
	}

	static boolean hasMultipleK(final String[] args){
		for(String arg : args){
			final int equals=arg.indexOf('=');
			String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			while(a.startsWith("-")){a=a.substring(1);}
			final String b=(equals<0 ? null : arg.substring(equals+1));
			if(a.equals("k") && equals>=0 && arg.indexOf(',', equals+1)>=0){return true;}
			if(a.equals("assemblek") || a.equals("fusek") || a.equals("joink")
					|| a.equals("bridgek") || a.equals("graphk") || a.equals("fusemaxmismatches")
					|| a.equals("fusedeadends") || a.equals("fuseconflicts") || a.equals("fusecoverageratio")
					|| a.equals("graphmergecoverageratio")
					|| a.equals("fusepath") || a.equals("fusepathdepth") || a.equals("fusepathflank")
					|| a.equals("fusetrim") || a.equals("fusesupportk") || a.equals("fusionvectors")
					|| a.equals("fusenet") || a.equals("fusencutoff") || a.equals("korder")){return true;}
			if(a.equals("lowdepthcontigdiagstage") || a.equals("ldcdstage")){return true;}
			if((a.equals("lowdepthcontigdiag") || a.equals("diagnoselowdepthcontigs") || a.equals("ldcd"))
					&& isLowDepthDiagStage(b)){return true;}
		}
		return false;
	}

	static class Config {

		Config(String[] args){
			String kList=null, assembleText=null, fuseText=null, bridgeText=null;
			boolean assembleExplicit=false, fuseExplicit=false, bridgeExplicit=false, graphExplicit=false;
			for(String arg : args){
				final int equals=arg.indexOf('=');
				String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
				while(a.startsWith("-")){a=a.substring(1);}
				if(a.equals("packed")){Kmer.PACKED=Parse.parseBoolean(equals<0 ? null : arg.substring(equals+1));}
			}
			for(String arg : args){
				final int equals=arg.indexOf('=');
				String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
				final String b=(equals<0 ? null : arg.substring(equals+1));
				while(a.startsWith("-")){a=a.substring(1);}
				if(a.equals("k")){kList=b;}
				else if(a.equals("korder")){
					if("input".equalsIgnoreCase(b)){ordered=true;}
					else if("legacy".equalsIgnoreCase(b)){ordered=false;}
					else{throw new IllegalArgumentException("korder must be input or legacy: "+b);}
				}
				else if(a.equals("assemblek")){
					assembleText=b;
					assembleExplicit=(b==null || !b.equalsIgnoreCase("auto"));
				}else if(a.equals("fusek") || a.equals("joink")){
					fuseText=b;
					fuseExplicit=(b==null || !b.equalsIgnoreCase("auto"));
				}else if(a.equals("bridgek")){
					bridgeText=b;
					bridgeExplicit=(b==null || !b.equalsIgnoreCase("auto"));
				}else if(a.equals("fusemaxmismatches")){
					fuseMaxMismatches=Integer.parseInt(b);
				}else if(a.equals("fusedeadends")){
					fuseDeadEndsOnly=Parse.parseBoolean(b);
				}else if(a.equals("fuseconflicts")){
					fuseConflicts=Parse.parseBoolean(b);
				}else if(a.equals("fusecoverageratio")){
					fuseCoverageRatio=Float.parseFloat(b);
				}else if(a.equals("graphmergecoverageratio")){
					graphMergeCoverageRatio=Float.parseFloat(b);
				}else if(a.equals("fusesupportk")){
					throw new IllegalArgumentException("Experimental fusesupportk was replaced by fusepath=t (same-K, one table).");
				}else if(a.equals("fusepath")){
					fusePath=Parse.parseBoolean(b);
				}else if(a.equals("fusepathdepth")){
					fusePathDepth=Integer.parseInt(b);
				}else if(a.equals("fusepathflank")){
					fusePathFlank=("auto".equalsIgnoreCase(b) ? 0 : Parse.parseIntKMG(b));
				}else if(a.equals("fusetrim")){
					fuseTrim=Parse.parseBoolean(b);
				}else if(a.equals("fusionvectors")){
					if(b==null || b.length()==0){throw new IllegalArgumentException("fusionvectors requires a file prefix.");}
					fusionVectorPrefix=b.equalsIgnoreCase("null") ? null : b;
				}else if(a.equals("fusenet")){
					if(b==null || b.length()==0){throw new IllegalArgumentException("fusenet requires a model path or null.");}
					fuseNet=b.equalsIgnoreCase("null") ? null : b;
				}else if(a.equals("fusencutoff")){
					fuseCutoff=Float.parseFloat(b);
					fuseCutoffSet=true;
				}else if(a.equals("hashkmers") || a.equals("kmerhash") || a.equals("bridgehash") || a.equals("hashbridges")
						|| a.equals("hashonly") || a.equals("hashedkmers")){
					hashModeRequested=Tadpole.parseHashMode(b);
				}else if(a.equals("out") || a.equals("out1") || a.equals("oute") || a.equals("oute1")){out=b;}
				else if(a.equals("crosskmaxdepthratio") || a.equals("ckmdr")){maxDepthRatio=Float.parseFloat(b);}
				else if(a.equals("crosskpasses") || a.equals("ckpasses")){passes=Integer.parseInt(b);}
				else if(a.equals("graphk")){
					graphK=(b==null || b.equalsIgnoreCase("auto") ? -1 : parseK(b, "graphk"));
					graphExplicit=true;
				}
				else if(a.equals("simpleomnitigs") || a.equals("omnitigs")){
					simpleOmnitigs=Parse.parseBoolean(b);
				}else if(a.equals("graphcover") || a.equals("pathcover") || a.equals("nonredundantpaths")){
					graphCover=Parse.parseBoolean(b);
				}else if(a.equals("pop") || a.equals("popbubbles")){
					popBubbles=Parse.parseBoolean(b);
					common.add(arg);
				}else if(a.equals("gfa") || a.equals("outgfa")){
					outGfa=b;
				}else if(a.equals("lowdepthcontigdiag") || a.equals("diagnoselowdepthcontigs") || a.equals("ldcd")){
					if(isLowDepthDiagStage(b)){lowDepthContigDiag=true; lowDepthContigDiagStage=parseLowDepthDiagStage(b);}
					else{lowDepthContigDiag=Parse.parseBoolean(b);}
				}else if(a.equals("lowdepthcontigdiagstage") || a.equals("ldcdstage")){
					lowDepthContigDiag=true;
					lowDepthContigDiagStage=parseLowDepthDiagStage(b);
				}else if(a.equals("lowdepthcontigmaxlen") || a.equals("ldcmaxlen")){
					lowDepthContigMaxLen=(b==null || b.equalsIgnoreCase("auto") ? -1 : Parse.parseIntKMG(b));
				}else if(a.equals("lowdepthcontigmaxcov") || a.equals("ldcmaxcov")){
					lowDepthContigMaxCov=Float.parseFloat(b);
				}else if(a.equals("lowdepthcontigfraction") || a.equals("ldcfrac")){
					lowDepthContigFraction=Float.parseFloat(b);
				}else if(a.equals("lowdepthcontigtopology") || a.equals("ldctopology")){
					lowDepthContigTopology=Tadpole.parseLowDepthTopology(b);
				}else if(a.equals("evictlowdepthcontigs") || a.equals("removelowdepthcontigs") || a.equals("ldce")){
					evictLowDepthContigs=Parse.parseBoolean(b);
				}else if(a.equals("retainshortcontigs") || a.equals("retainshortgraph")){
					if(b==null || b.equalsIgnoreCase("auto")){retainShortContigsSet=false;}
					else{retainShortContigs=Parse.parseBoolean(b); retainShortContigsSet=true;}
				}else if(a.equals("sweeplen") || a.equals("graphsweeplen")){
					sweepContigLen=Parse.parseIntKMG(b);
				}else if(a.equals("classifygraphcontigs") || a.equals("classifycontigs") || a.equals("graphclassify")){
					classifyGraphContigs=Parse.parseBoolean(b);
				}else if(a.equals("graphclasslowmaxcov") || a.equals("gclmc")){
					graphClassLowMaxCov=Float.parseFloat(b);
				}else if(a.equals("graphclasslowfraction") || a.equals("gclf")){
					graphClassLowFraction=Float.parseFloat(b);
				}else if(a.equals("graphclassmediumfraction") || a.equals("gcmf")){
					graphClassMediumFraction=Float.parseFloat(b);
				}else if(a.equals("graphclasshighfraction") || a.equals("gchf")){
					graphClassHighFraction=Float.parseFloat(b);
				}else if(a.equals("emitsuspect") || a.equals("suspect")){
					final boolean x=Parse.parseBoolean(b);
					emitTerminal=x; emitBranchedTerminal=x; emitUnanchored=x; emitLoopback=x;
				}else if(a.equals("emitterminal")){
					emitTerminal=Parse.parseBoolean(b);
				}else if(a.equals("emitbranchedterminal")){
					emitBranchedTerminal=Parse.parseBoolean(b);
				}else if(a.equals("emitunanchored")){
					emitUnanchored=Parse.parseBoolean(b);
				}else if(a.equals("emitloopback")){
					emitLoopback=Parse.parseBoolean(b);
				}else if(a.equals("emitbranchedconnected")){
					emitBranchedConnected=Parse.parseBoolean(b);
				}else if(a.equals("emitmulticonnected")){
					emitMultiConnected=Parse.parseBoolean(b);
				}else if(a.equals("emitselfloop")){
					emitSelfLoop=Parse.parseBoolean(b);
				}else if(a.equals("emitconnectedmax") || a.equals("ecm")){
					emitConnectedMax=(b==null || b.equalsIgnoreCase("all") || b.equalsIgnoreCase("auto") ? -1 : Integer.parseInt(b));
					classifyGraphContigs=true;
				}else if(a.equals("evictsuspect") || a.equals("es")){
					final boolean x=Parse.parseBoolean(b);
					evictTerminal=x; evictBranchedTerminal=x; evictUnanchored=x; evictLoopback=x;
				}else if(a.equals("evictterminal")){
					evictTerminal=Parse.parseBoolean(b);
				}else if(a.equals("evictbranchedterminal")){
					evictBranchedTerminal=Parse.parseBoolean(b);
				}else if(a.equals("evictunanchored")){
					evictUnanchored=Parse.parseBoolean(b);
				}else if(a.equals("evictloopback")){
					evictLoopback=Parse.parseBoolean(b);
				}else if(a.equals("evictbranchedconnected")){
					evictBranchedConnected=Parse.parseBoolean(b);
				}else if(a.equals("evictmulticonnected")){
					evictMultiConnected=Parse.parseBoolean(b);
				}else if(a.equals("evictgraphclass") || a.equals("evictgraphtopology")){
					evictGraphTopologyMask=Tadpole.parseGraphTopologyMask(b);
				}else if(a.equals("evictgraphdepth") || a.equals("evictdepth")){
					evictGraphDepthMask=Tadpole.parseGraphDepthMask(b);
				}else if(a.equals("evictconnectedabove") || a.equals("eca")){
					evictConnectedAbove=(b==null || b.equalsIgnoreCase("none") || b.equalsIgnoreCase("false") ? -1 : Integer.parseInt(b));
				}else if(a.equals("showstats")){showStats=Parse.parseBoolean(b); common.add(arg);}
				else if(a.equals("dot") || a.equals("outdot")){
					throw new RuntimeException("DOT output is not yet supported by TadpoleMulti.");
				}else if(a.equals("mode") && !"contig".equalsIgnoreCase(b)){
					throw new RuntimeException("TadpoleMulti requires mode=contig.");
				}else{
					if(requiresExplicitTable(a, b)){explicitTableRequired=true;}
					if(a.equals("prealloc") || a.equals("preallocate")){preallocRequested=parsePrealloc(b);}
					common.add(arg);
				}
			}
			final int[] shorthand=parseKList(kList, "k", true);
			if(assembleExplicit){assembleK=parseK(assembleText, "assemblek");}
			else if(shorthand.length>0){assembleK=shorthand[0];}
			else{throw new RuntimeException("TadpoleMulti requires k=<list> or assemblek=<value>.");}
			fuseKs=(fuseExplicit ? parseKList(fuseText, "fusek", true)
					: selectBelow(shorthand, assembleK));
			bridgeKs=(bridgeExplicit ? parseKList(bridgeText, "bridgek", true)
					: (ordered || assembleExplicit ? shorthand : selectBelow(shorthand, assembleK)));
			bridgeInitial=bridgeExplicit && contains(bridgeKs, assembleK);
			phaseKs=makePhaseKs(shorthand, kList!=null);
			for(int k : fuseKs){
				if(k>=assembleK){
					throw new RuntimeException("fusek values must be shorter than assemblek="+assembleK+": "+k);
				}
			}
			if(out==null){throw new RuntimeException("TadpoleMulti requires an output file.");}
			if(maxDepthRatio<0){throw new RuntimeException("crosskmaxdepthratio must be nonnegative.");}
			if(passes<1){throw new RuntimeException("crosskpasses must be positive.");}
			if(fuseMaxMismatches<-1){throw new IllegalArgumentException("fusemaxmismatches must be -1 or nonnegative.");}
			CrossKTipOverlapper.validateCoverageRatio(fuseCoverageRatio);
			if(graphMergeCoverageRatio!=0){BubblePopper.CoverageGate.validateRatio(graphMergeCoverageRatio);}
			if(fuseConflicts && fuseMaxMismatches<0){
				throw new IllegalArgumentException("fuseconflicts requires fusemaxmismatches>=0.");
			}
			if(fusePathDepth<1){throw new IllegalArgumentException("fusepathdepth must be positive.");}
			if(fusePathFlank<0){throw new IllegalArgumentException("fusepathflank must be auto, zero or positive.");}
			if(fuseCutoffSet && (!Float.isFinite(fuseCutoff) || fuseCutoff<0 || fuseCutoff>1)){
				throw new IllegalArgumentException("fusencutoff must be finite and in [0,1].");
			}
			if((fuseNet!=null)!=fuseCutoffSet){
				throw new IllegalArgumentException("Specify fusenet and an explicit fusencutoff together.");
			}
			if(simpleOmnitigs && graphCover){throw new RuntimeException("simpleOmnitigs and graphCover are mutually exclusive output modes.");}
			if(lowDepthContigMaxLen==0 || lowDepthContigMaxLen< -1){throw new RuntimeException("lowDepthContigMaxLen must be positive or auto.");}
			if(lowDepthContigMaxCov<0){throw new RuntimeException("lowDepthContigMaxCov must be nonnegative.");}
			if(lowDepthContigFraction<0 || lowDepthContigFraction>1){throw new RuntimeException("lowDepthContigFraction must be from 0 to 1.");}
			if(graphClassificationRequested()){
				if(graphClassLowMaxCov<0){throw new RuntimeException("graphClassLowMaxCov must be nonnegative.");}
				if(graphClassLowFraction<0 || graphClassLowFraction>1){throw new RuntimeException("graphClassLowFraction must be from 0 to 1.");}
				if(graphClassMediumFraction<graphClassLowFraction){
					throw new RuntimeException("graphClassMediumFraction must be at least graphClassLowFraction.");
				}
				if(graphClassHighFraction<graphClassMediumFraction){
					throw new RuntimeException("graphClassHighFraction must be at least graphClassMediumFraction.");
				}
			}
			if(emitConnectedMax==0 || emitConnectedMax< -1){throw new RuntimeException("emitConnectedMax must be positive or all.");}
			if(evictConnectedAbove==0 || evictConnectedAbove< -1){throw new RuntimeException("evictConnectedAbove must be positive or none.");}
			if(evictGraphTopologyMask<0 || evictGraphTopologyMask>255){throw new RuntimeException("Invalid graph eviction topology mask.");}
			if(sweepContigLen<0){throw new RuntimeException("sweeplen must be nonnegative.");}
			if(graphClassificationRequested()){classifyGraphContigs=true;}
			if(finalGraphNeeded()){
				if(graphK<0){
					graphK=(!ordered ? assembleK : shorthand.length>0 ? shorthand[shorthand.length-1]
							: phaseKs[phaseKs.length-1]);
				}
			}else if(graphExplicit){throw new RuntimeException("graphk requires graph operations, graph classification, or lowDepthContigDiag=t.");}
			else{graphK=assembleK;}
		}

		private static int parseK(final String text, final String name){
			if(text==null || text.length()<1){throw new RuntimeException(name+" requires a positive kmer length.");}
			final int requested=Integer.parseInt(text);
			if(requested<1){throw new RuntimeException(name+" requires a positive kmer length: "+text);}
			return Kmer.getKbig(requested);
		}

		private int[] parseKList(final String text, final String name, final boolean nullIsEmpty){
			if(text==null){
				if(nullIsEmpty){return new int[0];}
				throw new RuntimeException(name+" requires a comma-delimited kmer list.");
			}
			if(text.length()<1 || text.equalsIgnoreCase("none") || text.equalsIgnoreCase("false")
					|| text.equalsIgnoreCase("f")){return new int[0];}
			final String[] split=text.split(",", -1);
			final int[] array=new int[split.length];
			for(int i=0; i<split.length; i++){array[i]=parseK(split[i], name);}
			if(ordered){return array;}
			Arrays.sort(array);
			final int[] descending=new int[array.length];
			int unique=0, last=-1;
			for(int i=array.length-1; i>=0; i--){
				final int value=array[i];
				if(unique==0 || value!=last){descending[unique++]=value; last=value;}
			}
			return Arrays.copyOf(descending, unique);
		}

		/** Reserves the first assemble-K occurrence for assembly; later repeats remain real phases. */
		private int[] makePhaseKs(final int[] shorthand, final boolean hasList){
			assert(assembleK>0) : "Phase construction requires the resolved initial assembly K.";
			final structures.IntList phases=new structures.IntList();
			phases.add(assembleK);
			boolean consumed=false;
			for(int k : shorthand){
				if(!consumed && k==assembleK){consumed=true;}
				else{phases.add(k);}
			}
			for(int[] requested : new int[][] {fuseKs, bridgeKs}){
				for(int k : requested){
					if(contains(phases.array, phases.size, k)){continue;}
					if(ordered && hasList){
						throw new IllegalArgumentException("Phase K="+k+" is absent from k; include it in the ordered k list.");
					}
					phases.add(k);
				}
			}
			return Arrays.copyOf(phases.array, phases.size);
		}

		/** Tests membership without sorting away the caller's phase order. */
		static boolean contains(final int[] values, final int k){return contains(values, values.length, k);}

		/** Searches only the populated prefix of a phase list. */
		private static boolean contains(final int[] values, final int length, final int k){
			assert(length<=values.length) : "Phase membership cannot inspect unused list capacity.";
			for(int i=0; i<length; i++){if(values[i]==k){return true;}}
			return false;
		}

		private static int[] selectBelow(final int[] source, final int ceiling){
			int count=0;
			for(int k : source){if(k<ceiling){count++;}}
			final int[] selected=new int[count];
			int next=0;
			for(int k : source){if(k<ceiling){selected[next++]=k;}}
			return selected;
		}

		boolean graphOperations(){return simpleOmnitigs || graphCover || outGfa!=null;}
		boolean earlyLowDepthDiag(){return lowDepthContigDiag && (lowDepthContigDiagStage&DIAG_EARLY)!=0;}
		boolean finalLowDepthDiag(){return lowDepthContigDiag && (lowDepthContigDiagStage&DIAG_FINAL)!=0;}
		boolean finalGraphNeeded(){
			return graphOperations() || finalLowDepthDiag() || graphClassificationRequested() || graphMergeCoverageRatio>0;
		}
		boolean explicitGraphClassificationRequested(){
			return classifyGraphContigs || emitTerminal || emitBranchedTerminal || emitUnanchored || emitLoopback
					|| emitBranchedConnected || emitMultiConnected || emitSelfLoop || emitConnectedMax>0
					|| graphTopologyEvictionRequested() || evictConnectedAbove>0;
		}
		boolean graphClassificationRequested(){
			return explicitGraphClassificationRequested() || sweepContigLen>0;
		}
		boolean graphTopologyEvictionRequested(){
			return evictGraphTopologyMask!=0 || evictUnanchored || evictTerminal || evictBranchedTerminal || evictLoopback ||
					evictBranchedConnected || evictMultiConnected;
		}
		boolean useHashBridgeTables(){return hashModeForAny(bridgeKs)>0;}
		int bridgeHashMode(){return hashModeForAny(bridgeKs);}
		int hashMode(final int k){
			if(explicitTableRequired || k<64){return Tadpole.HASH_EXPLICIT;}
			if(hashModeRequested>=0){return hashModeRequested;}
			if(k==64 && !preallocRequested){return Tadpole.HASH_EXPLICIT;}
			return preallocRequested ? Tadpole.HASH_FIXED : Tadpole.HASH_PAIR;
		}
		private int hashModeForAny(final int[] array){
			for(int k : array){
				final int mode=hashMode(k);
				if(mode>0){return mode;}
			}
			return Tadpole.HASH_EXPLICIT;
		}
		private static boolean parsePrealloc(String b){
			if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){return Parse.parseBoolean(b);}
			return Double.parseDouble(b)>0;
		}
		private static boolean requiresExplicitTable(String a, String b){
			if(a.equals("maxcountretain") || a.equals("maxcr") || a.equals("maxdepthretain") || a.equals("maxdr")){
				return b!=null && !b.equalsIgnoreCase("inf") && !b.equalsIgnoreCase("infinity");
			}
			if(a.equals("outkmers") || a.equals("outk") || a.equals("dump")){return b!=null && !b.equalsIgnoreCase("null");}
			if(a.equals("gchist")){return !(b!=null && (b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false") || b.equals("0")));}
			return false;
		}
		void applyFinalGraphOutput(final Tadpole tad){
			if(outGfa!=null){tad.setGfaOutput(outGfa);}
		}
		void applyLowDepthThresholds(final Tadpole tad){
			tad.lowDepthContigMaxLen=lowDepthContigMaxLen;
			tad.lowDepthContigMaxCov=lowDepthContigMaxCov;
			tad.lowDepthContigFraction=lowDepthContigFraction;
			tad.lowDepthContigTopology=lowDepthContigTopology;
		}
		void addLowDepthThresholdArgs(final ArrayList<String> list){
			list.add("lowdepthcontigmaxlen="+(lowDepthContigMaxLen<0 ? "auto" : lowDepthContigMaxLen));
			list.add("lowdepthcontigmaxcov="+lowDepthContigMaxCov);
			list.add("lowdepthcontigfraction="+lowDepthContigFraction);
			list.add("lowdepthcontigtopology="+Tadpole.lowDepthTopologyName(lowDepthContigTopology));
		}
		void applyLowDepthDiagnostic(final Tadpole tad){
			tad.lowDepthContigDiag=finalLowDepthDiag();
			applyLowDepthThresholds(tad);
		}
		void applyFinalGraphClassification(final Tadpole tad){
			/* The automatic sweep already ran before cross-k fusion and bridging. */
			tad.classifyGraphContigs=explicitGraphClassificationRequested();
			tad.emitTerminal=emitTerminal;
			tad.emitBranchedTerminal=emitBranchedTerminal;
			tad.emitUnanchored=emitUnanchored;
			tad.emitLoopback=emitLoopback;
			tad.emitBranchedConnected=emitBranchedConnected;
			tad.emitMultiConnected=emitMultiConnected;
			tad.emitSelfLoop=emitSelfLoop;
			tad.emitConnectedMax=emitConnectedMax;
			tad.evictTerminal=evictTerminal;
			tad.evictBranchedTerminal=evictBranchedTerminal;
			tad.evictUnanchored=evictUnanchored;
			tad.evictLoopback=evictLoopback;
			tad.evictBranchedConnected=evictBranchedConnected;
			tad.evictMultiConnected=evictMultiConnected;
			tad.evictGraphTopologyMask=evictGraphTopologyMask;
			tad.evictGraphDepthMask=evictGraphDepthMask;
			tad.sweepContigLen=0;
			tad.popBubbles=popBubbles;
			tad.evictConnectedAbove=evictConnectedAbove;
			tad.graphClassLowMaxCov=graphClassLowMaxCov;
			tad.graphClassLowFraction=graphClassLowFraction;
			tad.graphClassMediumFraction=graphClassMediumFraction;
			tad.graphClassHighFraction=graphClassHighFraction;
			applyLowDepthThresholds(tad);
		}
		boolean emitContig(final Contig c, final int minContig){
			if(c.length()>=minContig){return true;}
			if(!classifyGraphContigs || c.graphClass<0){return false;}
			if(c.graphClass==Contig.GRAPH_CONNECTED){return emitConnectedMax<0 || c.graphClassHop<=emitConnectedMax;}
			if(c.graphClass==Contig.GRAPH_TERMINAL){return emitTerminal;}
			if(c.graphClass==Contig.GRAPH_BRANCHED_TERMINAL){return emitBranchedTerminal;}
			if(c.graphClass==Contig.GRAPH_UNANCHORED){return emitUnanchored;}
			if(c.graphClass==Contig.GRAPH_LOOPBACK){return emitLoopback;}
			if(c.graphClass==Contig.GRAPH_BRANCHED_CONNECTED){return emitBranchedConnected;}
			if(c.graphClass==Contig.GRAPH_MULTI_CONNECTED){return emitMultiConnected;}
			if(c.graphClass==Contig.GRAPH_SELF_LOOP){return emitSelfLoop;}
			throw new IllegalStateException("Unknown graph class "+c.graphClass+" on contig "+c.id+".");
		}

		void printExecutionPlan(final Tadpole tad){
			final ByteBuilder extras=new ByteBuilder();
			tad.appendExecutionPlanExtras(extras);
			if(simpleOmnitigs){Tadpole.appendPlanWord(extras, "simpleomnitigs");}
			if(graphCover){Tadpole.appendPlanWord(extras, "graphcover");}
			if(lowDepthContigDiag){Tadpole.appendPlanWord(extras, "lowdepthcontigdiag");}
			if(classifyGraphContigs && !tad.classifyGraphContigs){Tadpole.appendPlanWord(extras, "classifygraphcontigs");}
			if(evictTerminal){Tadpole.appendPlanWord(extras, "evictterminal");}
			if(evictBranchedTerminal){Tadpole.appendPlanWord(extras, "evictbranchedterminal");}
			if(evictUnanchored){Tadpole.appendPlanWord(extras, "evictunanchored");}
			if(evictLoopback){Tadpole.appendPlanWord(extras, "evictloopback");}
			if(evictBranchedConnected){Tadpole.appendPlanWord(extras, "evictbranchedconnected");}
			if(evictMultiConnected){Tadpole.appendPlanWord(extras, "evictmulticonnected");}
			if(evictGraphTopologyMask!=0){
				Tadpole.appendPlanWord(extras, "evictgraphclass="+Tadpole.graphTopologyMaskName(evictGraphTopologyMask));
			}
			if(graphTopologyEvictionRequested() || evictConnectedAbove>0){
				Tadpole.appendPlanWord(extras, "evictgraphdepth="+Tadpole.graphDepthMaskName(evictGraphDepthMask));
			}
			if(evictConnectedAbove>0){Tadpole.appendPlanWord(extras, "evictconnectedabove="+evictConnectedAbove);}
			Tadpole.printPlanLine("mode", "assemble");
			if(extras.length()>0){Tadpole.printPlanLine("extra", extras.toString());}
			Tadpole.printPlanLine("assemblek", assembleK);
			Tadpole.printPlanLine("korder", ordered ? "input" : "legacy");
			if(ordered){Tadpole.printPlanLine("phasek", toKList(phaseKs));}
			if(fuseKs.length>0){Tadpole.printPlanLine("fusek", toKList(fuseKs));}
			if(fuseMaxMismatches>=0){Tadpole.printPlanLine("fusemaxmismatches", fuseMaxMismatches);}
			if(fuseDeadEndsOnly){Tadpole.printPlanLine("fusedeadends", "true");}
			if(fuseConflicts){Tadpole.printPlanLine("fuseconflicts", "true");}
			if(fuseCoverageRatio>0){Tadpole.printPlanLine("fusecoverageratio", ""+fuseCoverageRatio);}
			if(graphMergeCoverageRatio>0){Tadpole.printPlanLine("graphmergecoverageratio", ""+graphMergeCoverageRatio);}
			if(!fuseTrim){Tadpole.printPlanLine("fusetrim", "false");}
			if(fusionVectorPrefix!=null){Tadpole.printPlanLine("fusionvectors", fusionVectorPrefix);}
			if(fuseNet!=null){
				Tadpole.printPlanLine("fusenet", fuseNet);
				Tadpole.printPlanLine("fusencutoff", ""+fuseCutoff);
			}
			if(fusePath){
				Tadpole.printPlanLine("fusepath", "true");
				Tadpole.printPlanLine("fusepathdepth", fusePathDepth);
				Tadpole.printPlanLine("fusepathflank", fusePathFlank==0 ? "auto" : ""+fusePathFlank);
			}
			if(bridgeKs.length>0){Tadpole.printPlanLine("bridgek", toKList(bridgeKs));}
			final int displayedHashMode=displayedHashMode();
			if(displayedHashMode>0){
				Tadpole.printPlanLine("hashkmers", Tadpole.hashModeName(displayedHashMode));
			}
			if(finalGraphNeeded()){Tadpole.printPlanLine("graphk", graphK);}
			if(lowDepthContigDiag){Tadpole.printPlanLine("diagstage", lowDepthDiagStageName(lowDepthContigDiagStage));}
			Tadpole.outstream.println();
		}

		private static String toKList(final int[] array){
			final ByteBuilder bb=new ByteBuilder(array.length*4);
			for(int i=0; i<array.length; i++){
				if(i>0){bb.comma();}
				bb.append(array[i]);
			}
			return bb.toString();
		}

		private int displayedHashMode(){
			final int bridgeMode=hashModeForAny(bridgeKs);
			if(bridgeMode>0){return bridgeMode;}
			return finalGraphNeeded() ? hashMode(graphK) : Tadpole.HASH_EXPLICIT;
		}

		final ArrayList<String> common=new ArrayList<String>();
		final int assembleK;
		final int[] fuseKs, bridgeKs, phaseKs;
		final boolean bridgeInitial;
		boolean ordered=true;
		String out, outGfa;
		String fusionVectorPrefix;
		String fuseNet;
		float fuseCutoff=Float.NaN;
		boolean fuseCutoffSet=false;
		float maxDepthRatio=3;
		int passes=10;
		int fuseMaxMismatches=-1;
		boolean fuseDeadEndsOnly=false;
		boolean fuseConflicts=false;
		/** Optional experimental whole-contig depth ratio; zero disables the veto. */
		float fuseCoverageRatio=0;
		/** Optional direct-merge depth policy, scoped only to normal final-graph simplification. */
		float graphMergeCoverageRatio=0;
		boolean fusePath=false;
		int fusePathDepth=1;
		int fusePathFlank=0;
		boolean fuseTrim=true;
		int graphK=-1;
		boolean simpleOmnitigs=false, graphCover=false, lowDepthContigDiag=false, evictLowDepthContigs=false, popBubbles=true;
		boolean classifyGraphContigs=false;
		boolean emitTerminal=false, emitBranchedTerminal=false, emitUnanchored=false, emitLoopback=false;
		boolean emitBranchedConnected=false, emitMultiConnected=false, emitSelfLoop=false;
		boolean evictTerminal=false, evictBranchedTerminal=false, evictUnanchored=false, evictLoopback=false;
		boolean evictBranchedConnected=false, evictMultiConnected=false;
		int evictGraphTopologyMask=0;
		int lowDepthContigDiagStage=DIAG_FINAL;
		int lowDepthContigMaxLen=-1;
		float lowDepthContigMaxCov=3, lowDepthContigFraction=0.2f;
		int lowDepthContigTopology=0, emitConnectedMax=-1, evictConnectedAbove=-1;
		int evictGraphDepthMask=1<<Contig.DEPTH_LOW;
		int sweepContigLen=500;
		float graphClassLowMaxCov=4, graphClassLowFraction=0.2f;
		float graphClassMediumFraction=0.4f, graphClassHighFraction=2.5f;
		boolean explicitTableRequired=false, preallocRequested=false;
		boolean retainShortContigs=false, retainShortContigsSet=false;
		int hashModeRequested=Tadpole.HASH_AUTO;
		boolean showStats=true;
	}

	private final Config config;
	private Tadpole fusionSupportTadpole;
	private FusionKmerSupport fusionSupport;
	private FusionJoinCollector fusionCollector;
	private FusionNeuralGate fusionNeural;
	private Tadpole.LowDepthDiagnostic finalLowDepthDiagnostic;

	private static boolean isLowDepthDiagStage(final String s){
		return s!=null && (s.equalsIgnoreCase("early") || s.equalsIgnoreCase("final") || s.equalsIgnoreCase("both"));
	}

	private static int parseLowDepthDiagStage(final String s){
		if(s==null || s.equalsIgnoreCase("final")){return DIAG_FINAL;}
		if(s.equalsIgnoreCase("early")){return DIAG_EARLY;}
		if(s.equalsIgnoreCase("both")){return DIAG_EARLY|DIAG_FINAL;}
		throw new RuntimeException("lowDepthContigDiagStage must be early, final, or both: "+s);
	}

	private static String lowDepthDiagStageName(final int stage){
		return stage==(DIAG_EARLY|DIAG_FINAL) ? "both" : stage==DIAG_EARLY ? "early" : "final";
	}

	private static final class EarlyLowDepthTracker {
		EarlyLowDepthTracker(final Tadpole.LowDepthDiagnostic diagnostic){
			candidates=new ArrayList<Contig>(diagnostic.candidates);
			lengths=new int[candidates.size()];
			for(int i=0; i<lengths.length; i++){lengths[i]=candidates.get(i).length(); bases+=lengths[i];}
		}
		void reportFates(final Tadpole.LowDepthDiagnostic diagnostic){
			final IdentityHashMap<Contig, Boolean> live=new IdentityHashMap<Contig, Boolean>(diagnostic.live.size()*2+1);
			final IdentityHashMap<Contig, Boolean> finalCandidates=new IdentityHashMap<Contig, Boolean>(diagnostic.candidates.size()*2+1);
			for(Contig c : diagnostic.live){live.put(c, Boolean.TRUE);}
			for(Contig c : diagnostic.candidates){finalCandidates.put(c, Boolean.TRUE);}
			int present=0, absorbed=0, grown=0, stillCandidate=0;
			long grownBases=0;
			for(int i=0; i<candidates.size(); i++){
				final Contig c=candidates.get(i);
				if(!live.containsKey(c)){absorbed++; continue;}
				present++;
				if(c.length()>lengths[i]){grown++; grownBases+=c.length()-lengths[i];}
				if(finalCandidates.containsKey(c)){stillCandidate++;}
			}
			Tadpole.outstream.println("Early low-depth candidate fates: initial="+candidates.size()+"/"+bases+
					" contigs/bases, present="+present+", retiredOrAbsorbed="+absorbed+", grown="+grown+
					"/+"+grownBases+" bases, stillCandidates="+stillCandidate+
					", reclassified="+(present-stillCandidate)+".");
		}
		final ArrayList<Contig> candidates;
		final int[] lengths;
		long bases;
	}

	private static final int DIAG_EARLY=1, DIAG_FINAL=2;
}
