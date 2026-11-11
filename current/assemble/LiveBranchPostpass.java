package assemble;

import java.util.ArrayList;
import java.util.Locale;

import parse.PreParser;
import shared.Shared;

/**
 * Development-only comparison of the shipped postpass with and without live rescue.
 * Both modes still require the ordinary contigrescuenet and contigrescuecutoff
 * binding; livenet/livecutoff optionally add a separate gate during the native walk.
 * Each invocation runs one selected mode, not both arms of a comparison. The ordinary
 * postpass takes a 27-input model; the optional live gate takes a separate 12-input model.
 * The combined run shares the original count table and population moments with the
 * subsequent postpass, without reloading reads. This wrapper appends simdfma=f and
 * does not restore shared settings. Output validation, writing and cleanup are delegated
 * to SameKJoinRescue; partial output and score files are not rolled back on failure.
 * @author Fischl
 */
public final class LiveBranchPostpass{

	/**
	 * Expands standard configs while keeping this experimental binding out of tadpole.sh.
	 * Closes the preparser's owned stream, preserving standard output/error.
	 * @param args Ordinary postpass options plus optional livenet/livecutoff
	 */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchPostpass.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Consumes livenet/livecutoff and delegates all other options to the existing postpass.
	 * Omitted, empty and literal "null" livenet values select the control with the
	 * NaN cutoff sentinel. The ordinary postpass model/cutoff remain required.
	 * A non-NaN livecutoff without a model is rejected; NaN with a model is rejected
	 * by the combined driver. The two cutoff options apply to their respective models.
	 * With a model, runCombined validates a finite cutoff in [0,1], model
	 * readability and output aliases in addition to the ordinary postpass checks.
	 * The delegated postpass requires one FASTA output, honors native overwrite and
	 * rejects append. Values preserve text after the first equals sign; underlying
	 * drivers load models, write outputs and own cleanup.
	 * @param args Expanded comparison and ordinary postpass arguments
	 */
	private static void run(final String[] args){
		if(args==null){throw new IllegalArgumentException("Expected postpass arguments and optional livenet/livecutoff.");}
		String model=null;
		float cutoff=Float.NaN;
		final ArrayList<String> postpassArgs=new ArrayList<String>();
		for(String arg : args){
			final int eq=arg.indexOf('=');
			String name=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=eq<0 ? null : arg.substring(eq+1);
			// PreParser rewrites literal null to an empty value before this parser.
			if(name.equals("livenet")){model=value==null || value.isEmpty() || value.equalsIgnoreCase("null") ? null : value;}
			else if(name.equals("livecutoff")){cutoff=Float.parseFloat(value);}
			else{postpassArgs.add(arg);}
		}
		if(model==null && !Float.isNaN(cutoff)){
			throw new IllegalArgumentException("livecutoff requires livenet; omit both for the postpass-only control.");
		}
		postpassArgs.add("simdfma=f");
		final String[] nativeArgs=postpassArgs.toArray(new String[0]);
		if(model==null){SameKJoinRescue.run(nativeArgs);}
		else{SameKJoinRescue.runCombined(nativeArgs, model, cutoff);}
	}
}
