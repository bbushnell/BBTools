package prok;

import java.util.Map;

import prot.MagQCCompositeInference;

/**
 * Assembly FASTA completeness/contamination client, including concurrent batches.
 * Implementation shares the MAG-QC caller, assignment and inference components.
 * @author Yoimiya, Nilou
 */
public final class ProkCC {

	/** Passes the public arguments to the shared native assembly engine. */
	public static void main(String[] args) throws Exception{
		//Embedded callers bypass prokcc.sh, so report installation help before config parsing fails.
		for(String arg:args){
			if(arg.regionMatches(true, 0, "config=", 0, 7)){
				final String config=arg.substring(7);
				if(!new java.io.File(config).isFile()){
					throw new IllegalArgumentException("Missing ProkCC release config: "+config+"\n"+
						shared.Resources.prokccDownloadInstructions());
				}
			}
		}
		prot.MagQCAssemblyBatch.main(args);
	}

	private ProkCC(){}

	/** Resolves human-summary routing after input expansion; literal file names retain their case. */
	public static String summaryDestination(String selection, int inputFiles){
		if(inputFiles<1){throw new IllegalArgumentException("Summary routing requires at least one input file");}
		if(selection==null || selection.isEmpty()){
			throw new IllegalArgumentException("summary requires auto, stderr, stdout, null, or a file path");
		}
		if("auto".equalsIgnoreCase(selection)){return inputFiles==1 ? "stderr" : null;}
		if("null".equalsIgnoreCase(selection)){return null;}
		if("stderr".equalsIgnoreCase(selection)){return "stderr";}
		if("stdout".equalsIgnoreCase(selection)){return "stdout";}
		return selection;
	}

	/**
	 * Scores the shipped default composite, loading its release and model on first
	 * call inside this monitor. Input is the raw composite vector; output order and
	 * width follow the selected model's declared contract. Input is copied and
	 * outputs are private.
	 * No gene calling, assignment, taxonomy query, or subnet inference is performed.
	 * A JVM already bound to another configured release is rejected, never reused.
	 */
	public static synchronized float[] calcCCSynced(float[] input){
		if(input==null || input.length==0){throw new IllegalArgumentException("Composite input must be nonempty");}
		final MagQCCompositeInference.Binding binding=defaultBinding==null ?
			MagQCCompositeInference.shippedDefault(input.length) : defaultBinding;
		final float[] result=calcCCSynced(input, binding);
		defaultBinding=binding;
		return result;
	}

	/**
	 * Batch overload with an immutable release binding. The first successful load
	 * fixes one configured model for this JVM; mismatched paths, pins or widths
	 * fail before scoring. Loading failure leaves the cache unbound and retryable.
	 */
	public static synchronized float[] calcCCSynced(float[] input, MagQCCompositeInference.Binding binding){
		if(binding==null || input==null || input.length!=binding.inputs()){
			throw new IllegalArgumentException("Composite input and binding width must agree");
		}
		if(compositeBinding!=null && !compositeBinding.equals(binding)){
			throw new IllegalArgumentException("Synchronized composite is already bound to a different release");
		}
		if(composite==null){
			final MagQCCompositeInference loaded=new MagQCCompositeInference(binding);
			composite=loaded; compositeBinding=binding;
		}
		return composite.score(input);
	}

	/**
	 * Optionally preloads a known-needed batch model inside the same cache monitor.
	 * The ordinary calcCCSynced API remains lazy. Width comes from the parsed model;
	 * the batch validates it against its independently loaded subnet layout before
	 * starting bin workers. A failed load leaves the cache unbound.
	 */
	public static synchronized MagQCCompositeInference.Binding preloadCCSynced(Map<String,String> options){
		if(compositeBinding!=null){
			final MagQCCompositeInference.Binding requested=
				new MagQCCompositeInference.Binding(options, compositeBinding.inputs());
			if(!compositeBinding.equals(requested)){
				throw new IllegalArgumentException("Synchronized composite is already bound to a different release");
			}
			return compositeBinding;
		}
		final MagQCCompositeInference loaded=MagQCCompositeInference.loadConfigured(options);
		composite=loaded; compositeBinding=loaded.binding();
		return compositeBinding;
	}

	/** All cache access, including first publication, is guarded by ProkCC.class. */
	private static MagQCCompositeInference composite;
	private static MagQCCompositeInference.Binding compositeBinding, defaultBinding;
}
