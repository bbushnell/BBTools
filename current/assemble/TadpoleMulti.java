package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Timer;
import shared.Tools;
import ukmer.Kmer;

/**
 * Runs a long-k Tadpole assembly, exhausts exact tip overlaps at every requested
 * shorter k, then bridges graph-disconnected ends with progressively shorter tables.
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
		final Tadpole longest=Tadpole.makeTadpole(makeArgs(config.kmers[0], true), true);
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
		for(int i=1; i<config.kmers.length; i++){
			final int k=config.kmers[i], before=contigs.size();
			final Timer timer=new Timer();
			longest.setContigs(contigs);
			longest.clearContigEdges();
			final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, k, config.kmers[0]-1);
			if(overlapper.addEdges()>0){mergeCrossK(longest);}
			contigs=longest.detachContigs();
			checkErrorState(longest);
			timer.stop();
			System.err.println("Cross-k overlaps "+k+": "+before+" -> "+contigs.size()+" contigs; "+timer);
		}

		/* Short-k tables are needed only for unbranched paths across actual sequence gaps. */
		for(int i=1; i<config.kmers.length; i++){
			final int k=config.kmers[i], before=contigs.size();
			final int endpoints=countBridgeEndpoints(contigs);
			if(endpoints<2){
				System.err.println("Cross-k bridges "+k+": skipped; only "+endpoints+" eligible endpoint"+
						(endpoints==1 ? "." : "s."));
				break;
			}
			final Tadpole tad=Tadpole.makeTadpole(makeArgs(k, false), true);
			setCrossKGraph(tad, true);
			tad.loadKmers(new Timer());
			tad.cleanLoadedKmers();
			checkErrorState(tad);
			tad.setContigs(contigs);
			tad.clearContigEdges();
			tad.processContigs();
			mergeCrossK(tad);

			contigs=tad.detachContigs();
			checkErrorState(tad);
			if(config.graphOperations() && k==config.graphK && i==config.kmers.length-1){
				reusableGraphTadpole=tad;
			}else{tad.tables().clear();}
			System.gc();
			System.err.println("Cross-k bridges "+k+": "+before+" -> "+contigs.size()+" contigs.");
		}

		if(config.graphOperations()){
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
			tad.cleanLoadedKmers();
			checkErrorState(tad);
		}else{System.err.println("Reusing shortest-k table for final graph at k="+config.graphK+".");}
		setCrossKGraph(tad, false);
		tad.minContigLen=minContig;
		tad.refreshGraphEndpoints=true;
		/* Refresh graph-k endpoint topology before resolving exact overlaps that were
		 * ineligible under the original longest-k endpoint classifications. */
		final boolean resolveRepeats=tad.resolveRepeats;
		tad.resolveRepeats=false;
		tad.simpleOmnitigs=false;
		tad.graphCover=false;
		tad.setContigs(contigs);
		tad.clearContigEdges();
		tad.processContigs();
		tad.clearContigEdges();
		int graphOverlapBefore=-1;
		final Timer graphOverlapTimer=new Timer();
		if(config.graphK<config.kmers[0]){
			graphOverlapBefore=contigs.size();
			final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, config.graphK,
					config.kmers[0]-1, true);
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

	private static int countBridgeEndpoints(ArrayList<Contig> contigs){
		int count=0;
		for(Contig c : contigs){
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
		}
		return false;
	}

	static class Config {

		Config(String[] args){
			String kList=null;
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
				else if(a.equals("out") || a.equals("out1") || a.equals("oute") || a.equals("oute1")){out=b;}
				else if(a.equals("crosskmaxdepthratio") || a.equals("ckmdr")){maxDepthRatio=Float.parseFloat(b);}
				else if(a.equals("crosskpasses") || a.equals("ckpasses")){passes=Integer.parseInt(b);}
				else if(a.equals("graphk")){graphK=(b==null || b.equalsIgnoreCase("auto") ? -1 : Integer.parseInt(b));}
				else if(a.equals("simpleomnitigs") || a.equals("omnitigs")){
					simpleOmnitigs=Parse.parseBoolean(b);
				}else if(a.equals("graphcover") || a.equals("pathcover") || a.equals("nonredundantpaths")){
					graphCover=Parse.parseBoolean(b);
				}else if(a.equals("showstats")){showStats=Parse.parseBoolean(b); common.add(arg);}
				else if(a.equals("dot") || a.equals("outdot")){
					throw new RuntimeException("DOT output is not yet supported by TadpoleMulti.");
				}else if(a.equals("mode") && !"contig".equalsIgnoreCase(b)){
					throw new RuntimeException("TadpoleMulti requires mode=contig.");
				}else{common.add(arg);}
			}
			if(kList==null || kList.indexOf(',')<0){throw new RuntimeException("TadpoleMulti requires at least two comma-delimited k values.");}
			final String[] split=kList.split(",");
			kmers=new int[split.length];
			for(int i=0; i<split.length; i++){
				final int requested=Integer.parseInt(split[i]);
				if(requested<1){throw new RuntimeException("Kmer lengths must be positive: "+kList);}
				kmers[i]=Kmer.getKbig(requested);
			}
			if(graphK>0){graphK=Kmer.getKbig(graphK);}
			Arrays.sort(kmers);
			for(int i=0, j=kmers.length-1; i<j; i++, j--){
				final int temp=kmers[i];
				kmers[i]=kmers[j];
				kmers[j]=temp;
			}
			for(int i=1; i<kmers.length; i++){
				if(kmers[i]==kmers[i-1]){throw new RuntimeException("Duplicate kmer length: "+kmers[i]);}
			}
			if(out==null){throw new RuntimeException("TadpoleMulti requires an output file.");}
			if(maxDepthRatio<0){throw new RuntimeException("crosskmaxdepthratio must be nonnegative.");}
			if(passes<1){throw new RuntimeException("crosskpasses must be positive.");}
			if(simpleOmnitigs && graphCover){throw new RuntimeException("simpleOmnitigs and graphCover are mutually exclusive output modes.");}
			if(graphOperations()){
				if(graphK<0){graphK=kmers[kmers.length-1];}
				if(graphK<kmers[kmers.length-1] || graphK>kmers[0]){
					throw new RuntimeException("graphk must be between the shortest and longest assembly k: "+graphK);
				}
			}else if(graphK>=0){throw new RuntimeException("graphk requires simpleomnitigs=t or graphcover=t.");}
		}

		boolean graphOperations(){return simpleOmnitigs || graphCover;}

		final ArrayList<String> common=new ArrayList<String>();
		final int[] kmers;
		String out;
		float maxDepthRatio=3;
		int passes=10;
		int graphK=-1;
		boolean simpleOmnitigs=false, graphCover=false;
		boolean showStats=true;
	}

	private final Config config;
}
