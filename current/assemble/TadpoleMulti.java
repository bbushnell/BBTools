package assemble;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import shared.Timer;
import shared.Tools;

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
			if(tad instanceof Tadpole1){((Tadpole1)tad).crossKGraph=true;}
			else{((Tadpole2)tad).crossKGraph=true;}
			tad.loadKmers(new Timer());
			checkErrorState(tad);
			tad.setContigs(contigs);
			tad.clearContigEdges();
			tad.processContigs();
			mergeCrossK(tad);

			contigs=tad.detachContigs();
			checkErrorState(tad);
			tad.tables().clear();
			System.gc();
			System.err.println("Cross-k bridges "+k+": "+before+" -> "+contigs.size()+" contigs.");
		}

		writeContigs(contigs, config.out, minContig, idOffset);
		if(config.showStats && FileFormat.isFastaExt(ReadWrite.rawExtension(config.out)) && !FileFormat.isStdio(config.out)){
			System.err.println();
			jgi.AssemblyStats2.main(new String[] {"in="+config.out, "printextended"});
		}
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

	private static class Config {

		Config(String[] args){
			String kList=null;
			for(String arg : args){
				final int equals=arg.indexOf('=');
				String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
				final String b=(equals<0 ? null : arg.substring(equals+1));
				while(a.startsWith("-")){a=a.substring(1);}
				if(a.equals("k")){kList=b;}
				else if(a.equals("out") || a.equals("out1") || a.equals("oute") || a.equals("oute1")){out=b;}
				else if(a.equals("crosskmaxdepthratio") || a.equals("ckmdr")){maxDepthRatio=Float.parseFloat(b);}
				else if(a.equals("crosskpasses") || a.equals("ckpasses")){passes=Integer.parseInt(b);}
				else if(a.equals("showstats")){showStats=Parse.parseBoolean(b); common.add(arg);}
				else if(a.equals("dot") || a.equals("outdot")){
					throw new RuntimeException("DOT output is not yet supported by TadpoleMulti.");
				}else if(a.equals("mode") && !"contig".equalsIgnoreCase(b)){
					throw new RuntimeException("TadpoleMulti requires mode=contig.");
				}else{common.add(arg);}
			}
			if(kList==null || kList.indexOf(',')<0){throw new RuntimeException("TadpoleMulti requires descending comma-delimited k values, such as k=140,75,31.");}
			final String[] split=kList.split(",");
			kmers=new int[split.length];
			for(int i=0; i<split.length; i++){
				kmers[i]=Integer.parseInt(split[i]);
				if(kmers[i]<1 || (i>0 && kmers[i]>=kmers[i-1])){
					throw new RuntimeException("Kmer lengths must be positive and strictly descending: "+kList);
				}
			}
			if(out==null){throw new RuntimeException("TadpoleMulti requires an output file.");}
			if(maxDepthRatio<0){throw new RuntimeException("crosskmaxdepthratio must be nonnegative.");}
			if(passes<1){throw new RuntimeException("crosskpasses must be positive.");}
		}

		final ArrayList<String> common=new ArrayList<String>();
		final int[] kmers;
		String out;
		float maxDepthRatio=3;
		int passes=10;
		boolean showStats=true;
	}

	private final Config config;
}
