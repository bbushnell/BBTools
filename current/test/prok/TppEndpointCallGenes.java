package test.prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaCdsStopConstraint;
import prok.NcrnaFamily;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/** Opt-in TPP stop-search experiment using the normal CallGenes pipeline.
 * Argument parsing follows the native template; CallGenes owns input, worker
 * scheduling, gene calling and ordinary outputs. This driver owns only the
 * optional constraint log. Each invocation requires a fresh JVM.
 * @author Brian Bushnell, Raiden
 */
public final class TppEndpointCallGenes {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, TppEndpointCallGenes.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Experiment options are removed; every other argument reaches CallGenes. */
	private static void run(final String[] args){
		final ArrayList<String> forwarded=new ArrayList<String>();
		final ArrayList<String> paths=new ArrayList<String>();
		boolean enabled=false;
		String log=null;
		for(String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("tppcdsstop")){enabled=Parse.parseBoolean(b);}
			else if(a.equals("constraintlog")){log=b;}
			else{
				forwarded.add(arg);
				if(b!=null && (a.startsWith("out") || a.equals("in") || a.equals("in1") || a.equals("in2"))){paths.add(b);}
			}
		}
		if(!enabled && log!=null){throw new IllegalArgumentException("constraintlog requires tppcdsstop=t");}
		if(log!=null){
			// Fresh-only diagnostic output cannot overwrite a source or prior result.
			if(!Tools.testOutputFiles(false, false, false, log)){
				throw new IllegalArgumentException("Constraint log must be a new writable file: "+log);
			}
			paths.add(log);
			Tools.testForDuplicateFiles(true, paths.toArray(new String[0]));
		}
		final Timer timer=new Timer();
		final ByteStreamWriter writer=log==null ? null : new ByteStreamWriter(
			FileFormat.testOutput(log, FileFormat.TXT, null, true, false, false, false));
		if(writer!=null){writer.start();}
		final NcrnaCdsStopConstraint.Factory factory=enabled ? new NcrnaCdsStopConstraint.Factory(writer) : null;
		try{
			final CallGenes caller=new CallGenes(forwarded.toArray(new String[0]));
			caller.run(timer, factory==null ? null : () -> configure(factory));
		}finally{
			if(writer!=null && writer.poisonAndWait()){
				throw new RuntimeException("Failed writing constraint log: "+log);
			}
		}
		final long[] counts=factory==null ? new long[3] : factory.totals();
		System.err.println("TPP stop constraint: enabled="+enabled+", candidate sweeps="+counts[0]
			+", limited sweeps="+counts[1]+", excluded available stop positions="+counts[2]);
	}

	/** Runs after family loading and before any caller worker starts. */
	private static void configure(final NcrnaCdsStopConstraint.Factory factory){
		NcrnaFamily selected=null;
		for(NcrnaFamily family : GeneCaller.ncrnaFamilies){
			if("tpp".equals(family.name)){
				if(selected!=null){throw new IllegalStateException("Only one TPP family may own the experiment");}
				selected=family;
			}
		}
		if(selected==null){throw new IllegalArgumentException("tppcdsstop=t requires tpp=t and its explicit family resources");}
		if(selected.boundaryNetsByModel==null && (selected.boundary5NetTemplate==null || selected.boundary3NetTemplate==null)){
			throw new IllegalArgumentException("TPP stop constraint requires loaded endpoint networks");
		}
		selected.cdsStopConstraintFactory=factory;
	}
}
