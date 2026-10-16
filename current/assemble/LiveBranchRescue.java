package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import kmer.AbstractKmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;

/** Development-only raw contig comparison with optional live branch rescue. @author Fischl */
public final class LiveBranchRescue {

	/** Expands standard configs; no reference or labels are available to the live decision. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchRescue.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Forces identical raw-unitig/scalar settings in the no-model and model arms. */
	private static void run(final String[] args){
		if(args==null){throw new IllegalArgumentException("Expected read, output and optional network arguments.");}
		String output=null, model=null;
		float cutoff=Float.NaN;
		boolean verifyMoments=false;
		final ArrayList<String> nativeArgs=new ArrayList<String>();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(name.equals("out")){output=value;}
			else if(name.equals("net")){model=value;}
			else if(name.equals("cutoff")){cutoff=Float.parseFloat(value);}
			else if(name.equals("verifymoments")){verifyMoments=Parse.parseBoolean(value);}
			else if(name.startsWith("out") || name.equals("dot") || name.equals("gfa") || name.equals("khist") ||
					name.equals("hist") || name.equals("dump") || name.equals("ihist") || name.equals("inserthistogram")){
				throw new IllegalArgumentException("Only out=raw_contigs.fa is supported by this comparison driver.");
			}else{nativeArgs.add(arg);}
		}
		if(output==null || output.equalsIgnoreCase("null") || TadpoleMulti.hasMultipleK(args) ||
				(model==null ? !Float.isNaN(cutoff) : !Float.isFinite(cutoff) || cutoff<0 || cutoff>1)){
			throw new IllegalArgumentException("Required out=NEW.fa and single K; optional net=model.bbnet requires cutoff in[0,1].");
		}
		nativeArgs.addAll(Arrays.asList("mode=contig", "showstats=f", "processcontigs=f", "pop=f", "sweeplen=0",
				"retainshortcontigs=t", "trimends=0", "trimcircular=f", "hashkmers=f", "simdfma=f", "out="+output));
		final Tadpole tad=Tadpole.makeTadpole(nativeArgs.toArray(new String[0]), true);
		if(tad.k()<24 || tad.k()>127 || tad.processContigs || tad.evictLowDepthContigs || !tad.useOwnership ||
				!tad.tables().rcomp() || AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN || tad.minCountExtend<1){
			throw new IllegalArgumentException("Live comparison requires raw single-K24-127 canonical DNA unitigs with ownership and mce>=1.");
		}
		final ArrayList<String> inputs=new ArrayList<String>(tad.tables().in1);
		inputs.addAll(tad.tables().in2); inputs.addAll(tad.tables().extra);
		if(inputs.isEmpty()){throw new IllegalArgumentException("At least one read input is required.");}
		if(model!=null){inputs.add(model);}
		TadpoleGraph.checkPaths(inputs.get(0), inputs, output);
		if(!Tools.testInputFiles(true, true, inputs.toArray(new String[0])) || !Tools.testOutputFiles(false, false, false, output)){
			throw new IllegalArgumentException("Inputs must be readable and output must be new and distinct.");
		}
		final LiveBranchGate gate=model==null ? null : new LiveBranchGate(tad, model, cutoff);
		try{
			tad.liveBranchGate=gate;
			tad.process2(Tadpole.contigMode);
			if(tad.errorState || tad.tables().errorState){throw new IllegalStateException("Live branch assembly failed.");}
			if(verifyMoments){verifyMoments(tad);}
			final long[] totals=gate==null ? new long[2] : gate.totals();
			if(gate!=null){Tadpole.outstream.println("LIVE_BRANCH_MOMENTS seconds="+(gate.momentNanos()*1e-9));}
			Tadpole.outstream.println("LIVE_BRANCH_RESCUE_PASS k="+tad.k()+" evaluations="+totals[0]+" rescued_encounters="+totals[1]);
		}finally{
			tad.liveBranchGate=null;
			tad.clearData();
		}
	}

	/** Optional post-run oracle adds extra scans; it must not be used for end-to-end timing claims. */
	private static void verifyMoments(final Tadpole tad){
		final SameKJoinFeatures.Moments direct=SameKJoinFeatures.moments(tad);
		final SameKJoinFeatures.Moments walked=SameKJoinFeatures.momentsByWalker(tad);
		if(direct.unique!=walked.unique ||
				Double.doubleToRawLongBits(direct.observations)!=Double.doubleToRawLongBits(walked.observations) ||
				Double.doubleToRawLongBits(direct.squares)!=Double.doubleToRawLongBits(walked.squares)){
			throw new IllegalStateException("Count-only traversal changed exact depth moments; frozen model inputs must remain identical.");
		}
		Tadpole.outstream.println("LIVE_BRANCH_MOMENTS_PARITY_PASS unique="+direct.unique);
	}
}
