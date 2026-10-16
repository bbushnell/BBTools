package assemble;

import java.util.ArrayList;
import java.util.Locale;

import parse.PreParser;
import shared.Shared;

/** Development-only comparison of the shipped postpass with and without live rescue. @author Fischl */
public final class LiveBranchPostpass {

	/** Uses the standard preparser while keeping the experimental binding out of tadpole.sh. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, LiveBranchPostpass.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Both arms use scalar inference and the existing postpass's path/mode validation. */
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
			if(name.equals("livenet")){model=value==null || value.equalsIgnoreCase("null") ? null : value;}
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
