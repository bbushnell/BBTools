package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import kmer.AbstractKmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;

/**
 * Development-only collection at actual native contig branch checks. The same
 * conservative walk runs with or without the observer. Reference labeling and
 * training are separate processes, after the read count table is released.
 * @author Fischl
 */
public final class LiveBranchDiagnostic {

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Expands config files through the standard preparser, then forwards native arguments. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchDiagnostic.class, false);
		try{collect(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Installs the observer only on this single-K raw-unitig diagnostic instance. */
	private static void collect(final String[] args){
		if(args==null){throw new IllegalArgumentException("Expected diagnostic arguments.");}
		String output=null, raw=null;
		int controlEvery=64;
		final ArrayList<String> nativeArgs=new ArrayList<String>();
		for(String arg : args){
			final int eq=arg.indexOf('=');
			String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=eq<0 ? null : arg.substring(eq+1);
			if(name.equals("out")){output=value;}
			else if(name.equals("rawcontigs")){raw=value;}
			else if(name.equals("controlevery")){controlEvery=Parse.parseIntKMG(value);}
			else if(name.startsWith("out") || name.equals("dot") || name.equals("gfa") || name.equals("khist")
					|| name.equals("hist") || name.equals("dump") || name.equals("ihist") || name.equals("inserthistogram")){
				throw new IllegalArgumentException("Use only out=trace and rawcontigs=FASTA; unsupported output: "+name);
			}else{nativeArgs.add(arg);}
		}
		if(output==null || output.equalsIgnoreCase("null") || controlEvery<1 || TadpoleMulti.hasMultipleK(args)){
			throw new IllegalArgumentException("Required out=trace.tsv, controlevery>=1, and one K.");
		}
		// Override only the documented raw-unitig collection mode, not branch thresholds.
		nativeArgs.addAll(Arrays.asList("mode=contig", "showstats=f", "processcontigs=f", "pop=f", "sweeplen=0",
				"retainshortcontigs=t", "trimends=0", "trimcircular=f", "hashkmers=f", "out="+(raw==null ? "null" : raw)));
		final Tadpole tad=Tadpole.makeTadpole(nativeArgs.toArray(new String[0]), true);
		if(tad.k()<24 || tad.k()>127 || tad.processContigs || tad.evictLowDepthContigs || !tad.useOwnership ||
				!tad.tables().rcomp() || AbstractKmerTableSet.MASK_MIDDLE || Shared.AMINO_IN || tad.minCountExtend<1){
			throw new IllegalArgumentException("Trace requires raw single-K24–127 canonical DNA unitigs with ownership and mce>=1.");
		}
		final ArrayList<String> inputs=new ArrayList<String>(tad.tables().in1);
		inputs.addAll(tad.tables().in2); inputs.addAll(tad.tables().extra);
		if(inputs.isEmpty()){throw new IllegalArgumentException("At least one read input is required.");}
		TadpoleGraph.checkPaths(inputs.get(0), inputs, output, raw);
		if(!Tools.testInputFiles(true, true, inputs.toArray(new String[0])) ||
				!Tools.testOutputFiles(false, false, false, output, raw)){
			throw new IllegalArgumentException("Diagnostic inputs must be readable; trace/FASTA outputs must be new and distinct.");
		}
		try(final LiveBranchObserver observer=new LiveBranchObserver(output, tad.k(), tad.minCountExtend, controlEvery)){
			tad.liveBranchObserver=observer;
			tad.process2(Tadpole.contigMode);
			if(tad.errorState || tad.tables().errorState){throw new IllegalStateException("Assembly failed during branch observation.");}
			observer.complete(SameKJoinFeatures.moments(tad));
		}finally{
			tad.liveBranchObserver=null;
			tad.clearData();
		}
		Tadpole.outstream.println("LIVE_BRANCH_COLLECT_PASS k="+tad.k()+" controlevery="+controlEvery);
	}
}
