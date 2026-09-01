package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.Kmer;

/**
 * Runs a Tadpole assembly, exhausts exact tip overlaps at requested shorter kmer
 * lengths, then bridges graph-disconnected ends with independent kmer tables.
 *
 * @author Brian Bushnell, Noelle
 */
public class TadpoleMulti {

	public static void main(String[] args){
		final Config config=new Config(args);
		new TadpoleMulti(config).process();
	}

	TadpoleMulti(Config config_){config=config_;}

	void process(){
		final Tadpole longest=Tadpole.makeTadpole(makeArgs(config.assembleK, true), true);
		config.printExecutionPlan(longest);
		if(!Tools.testOutputFiles(Tadpole.overwrite, Tadpole.append, false, config.out)){
			throw new RuntimeException("Can't write output file "+config.out+"; overwrite="+Tadpole.overwrite);
		}
		longest.process2(Tadpole.contigMode);
		checkErrorState(longest);
		longest.markBridgeEndpoints();
		longest.clearContigEdges();
		ArrayList<Contig> contigs=longest.detachContigs();
		final int minContig=longest.minContigLen;
		final int idOffset=longest.contigIDOffset;
		longest.tables().clear();
		System.gc();

		/* Exact terminal overlaps need only contig sequence; exhaust them before rereading the reads. */
		Tadpole reusableGraphTadpole=null;
		for(int i=0; i<config.fuseKs.length; i++){
			final int k=config.fuseKs[i], before=contigs.size();
			final Timer timer=new Timer();
			longest.setContigs(contigs);
			longest.clearContigEdges();
			final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, k,
					config.assembleK-1, false, minContig);
			if(overlapper.addEdges()>0){mergeCrossK(longest);}
			contigs=longest.detachContigs();
			checkErrorState(longest);
			timer.stop();
			System.err.println("Cross-k overlaps "+k+": "+before+" -> "+contigs.size()+" contigs; "+timer);
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

		writeContigs(contigs, config.out, minContig, idOffset);
		if(config.showStats && FileFormat.isFastaExt(ReadWrite.rawExtension(config.out)) && !FileFormat.isStdio(config.out)){
			System.err.println();
			jgi.AssemblyStats2.main(new String[] {"in="+config.out, "printextended"});
		}
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
			graphOverlapBefore=contigs.size();
			final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, config.graphK,
					config.assembleK-1, true, minContig);
			if(overlapper.addEdges()>0){mergeCrossK(tad);}
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
		config.applyFinalGraphOutput(tad);
		tad.setContigs(merged);
		tad.clearContigEdges();
		tad.processContigs();
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
		BubblePopper.popDirect=true;
		BubblePopper.popIndirect=false;
		BubblePopper.crossKMerge=true;
		BubblePopper.crossKMaxDepthRatio=config.maxDepthRatio;
		try{
			for(int pass=0, merged=1; pass<config.passes && merged>0; pass++){
				merged=tad.popBubbles(false);
			}
		}finally{
			BubblePopper.popDirect=oldDirect;
			BubblePopper.popIndirect=oldIndirect;
			BubblePopper.crossKMerge=oldCrossK;
			BubblePopper.crossKMaxDepthRatio=oldDepthRatio;
		}
	}

	private static void checkErrorState(Tadpole tad){
		if(tad.errorState || tad.tables().errorState){
			throw new RuntimeException(tad.getClass().getSimpleName()+" terminated in an error state; no assembly will be written.");
		}
	}

	private String[] makeArgs(final int k, final boolean longest){
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
		if(!longest){list.add("pop=f");}
		return list.toArray(new String[list.size()]);
	}

	private static void writeContigs(ArrayList<Contig> contigs, String out, int minContig, int idOffset){
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
			if(c.length()>=minContig){bsw.println(c);}
		}
		if(bsw.poisonAndWait()){throw new RuntimeException("Error writing "+out);}
	}

	static boolean hasMultipleK(final String[] args){
		for(String arg : args){
			final int equals=arg.indexOf('=');
			String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			while(a.startsWith("-")){a=a.substring(1);}
			if(a.equals("k") && equals>=0 && arg.indexOf(',', equals+1)>=0){return true;}
			if(a.equals("assemblek") || a.equals("fusek") || a.equals("joink")
					|| a.equals("bridgek") || a.equals("graphk")){return true;}
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
				else if(a.equals("assemblek")){
					assembleText=b;
					assembleExplicit=(b==null || !b.equalsIgnoreCase("auto"));
				}else if(a.equals("fusek") || a.equals("joink")){
					fuseText=b;
					fuseExplicit=(b==null || !b.equalsIgnoreCase("auto"));
				}else if(a.equals("bridgek")){
					bridgeText=b;
					bridgeExplicit=(b==null || !b.equalsIgnoreCase("auto"));
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
				}else if(a.equals("gfa") || a.equals("outgfa")){
					outGfa=b;
				}else if(a.equals("lowdepthcontigdiag") || a.equals("diagnoselowdepthcontigs") || a.equals("ldcd")){
					lowDepthContigDiag=Parse.parseBoolean(b);
				}else if(a.equals("lowdepthcontigmaxlen") || a.equals("ldcmaxlen")){
					lowDepthContigMaxLen=(b==null || b.equalsIgnoreCase("auto") ? -1 : Parse.parseIntKMG(b));
				}else if(a.equals("lowdepthcontigmaxcov") || a.equals("ldcmaxcov")){
					lowDepthContigMaxCov=Float.parseFloat(b);
				}else if(a.equals("lowdepthcontigfraction") || a.equals("ldcfrac")){
					lowDepthContigFraction=Float.parseFloat(b);
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
					: (assembleExplicit ? shorthand : selectBelow(shorthand, assembleK)));
			for(int k : fuseKs){
				if(k>=assembleK){
					throw new RuntimeException("fusek values must be shorter than assemblek="+assembleK+": "+k);
				}
			}
			if(out==null){throw new RuntimeException("TadpoleMulti requires an output file.");}
			if(maxDepthRatio<0){throw new RuntimeException("crosskmaxdepthratio must be nonnegative.");}
			if(passes<1){throw new RuntimeException("crosskpasses must be positive.");}
			if(simpleOmnitigs && graphCover){throw new RuntimeException("simpleOmnitigs and graphCover are mutually exclusive output modes.");}
			if(lowDepthContigMaxLen==0 || lowDepthContigMaxLen< -1){throw new RuntimeException("lowDepthContigMaxLen must be positive or auto.");}
			if(lowDepthContigMaxCov<0){throw new RuntimeException("lowDepthContigMaxCov must be nonnegative.");}
			if(lowDepthContigFraction<0 || lowDepthContigFraction>1){throw new RuntimeException("lowDepthContigFraction must be from 0 to 1.");}
			if(finalGraphNeeded()){
				if(graphK<0){graphK=assembleK;}
			}else if(graphExplicit){throw new RuntimeException("graphk requires graph operations or lowDepthContigDiag=t.");}
			else{graphK=assembleK;}
		}

		private static int parseK(final String text, final String name){
			if(text==null || text.length()<1){throw new RuntimeException(name+" requires a positive kmer length.");}
			final int requested=Integer.parseInt(text);
			if(requested<1){throw new RuntimeException(name+" requires a positive kmer length: "+text);}
			return Kmer.getKbig(requested);
		}

		private static int[] parseKList(final String text, final String name, final boolean nullIsEmpty){
			if(text==null){
				if(nullIsEmpty){return new int[0];}
				throw new RuntimeException(name+" requires a comma-delimited kmer list.");
			}
			if(text.length()<1 || text.equalsIgnoreCase("none") || text.equalsIgnoreCase("false")
					|| text.equalsIgnoreCase("f")){return new int[0];}
			final String[] split=text.split(",");
			final int[] array=new int[split.length];
			for(int i=0; i<split.length; i++){array[i]=parseK(split[i], name);}
			Arrays.sort(array);
			final int[] descending=new int[array.length];
			int unique=0, last=-1;
			for(int i=array.length-1; i>=0; i--){
				final int value=array[i];
				if(unique==0 || value!=last){descending[unique++]=value; last=value;}
			}
			return Arrays.copyOf(descending, unique);
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
		boolean finalGraphNeeded(){return graphOperations() || lowDepthContigDiag;}
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
		void applyLowDepthDiagnostic(final Tadpole tad){
			tad.lowDepthContigDiag=lowDepthContigDiag;
			tad.lowDepthContigMaxLen=lowDepthContigMaxLen;
			tad.lowDepthContigMaxCov=lowDepthContigMaxCov;
			tad.lowDepthContigFraction=lowDepthContigFraction;
		}

		void printExecutionPlan(final Tadpole tad){
			final ByteBuilder extras=new ByteBuilder();
			tad.appendExecutionPlanExtras(extras);
			if(simpleOmnitigs){Tadpole.appendPlanWord(extras, "simpleomnitigs");}
			if(graphCover){Tadpole.appendPlanWord(extras, "graphcover");}
			if(lowDepthContigDiag){Tadpole.appendPlanWord(extras, "lowdepthcontigdiag");}
			Tadpole.printPlanLine("mode", "assemble");
			if(extras.length()>0){Tadpole.printPlanLine("extra", extras.toString());}
			Tadpole.printPlanLine("assemblek", assembleK);
			if(fuseKs.length>0){Tadpole.printPlanLine("fusek", toKList(fuseKs));}
			if(bridgeKs.length>0){Tadpole.printPlanLine("bridgek", toKList(bridgeKs));}
			final int displayedHashMode=displayedHashMode();
			if(displayedHashMode>0){
				Tadpole.printPlanLine("hashkmers", Tadpole.hashModeName(displayedHashMode));
			}
			if(finalGraphNeeded()){Tadpole.printPlanLine("graphk", graphK);}
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
		final int[] fuseKs, bridgeKs;
		String out, outGfa;
		float maxDepthRatio=3;
		int passes=10;
		int graphK=-1;
		boolean simpleOmnitigs=false, graphCover=false, lowDepthContigDiag=false;
		int lowDepthContigMaxLen=-1;
		float lowDepthContigMaxCov=3, lowDepthContigFraction=0.2f;
		boolean explicitTableRequired=false, preallocRequested=false;
		int hashModeRequested=Tadpole.HASH_AUTO;
		boolean showStats=true;
	}

	private final Config config;
}
