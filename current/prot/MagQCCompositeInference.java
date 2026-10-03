package prot;

import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

import dna.Data;

/**
 * Resource-bound six-head composite inference for the synchronized ProkCC API.
 * A binding contains no loaded networks. Loading validates the same contract as
 * the native harness, and discards the construction master after making one worker.
 * @author Nilou
 */
public final class MagQCCompositeInference {

	/** Loads once for an owning monitor; a failed load publishes no cached state. */
	public MagQCCompositeInference(Binding binding){
		this.binding=binding;
		scorer=load(binding.options, binding.width);
	}

	/** Loads a known-needed model before its independent subnet resources finish loading. */
	private MagQCCompositeInference(HashMap<String,String> options){
		scorer=load(options, 0);
		binding=new Binding(options, scorer.input.length);
	}

	/** Takes an owned option snapshot and validates before publishing an inference instance. */
	public static MagQCCompositeInference loadConfigured(Map<String,String> options){
		if(options==null){throw new IllegalArgumentException("Configured composite requires release options");}
		return new MagQCCompositeInference(new HashMap<String,String>(options));
	}

	/** Validated identity, including model-declared width, for the synchronized cache. */
	public Binding binding(){return binding;}

	/** Width zero is internal loading without a vector; other callers supply a validated width. */
	private static MagQCNetworkHarness.Scorer load(HashMap<String,String> options, int width){
		try{
			for(String key:RESOURCES){MagQCNetworkHarness.pinned(options, key, key+"sha80");}
			final MagQCNetworkHarness.Scorer loaded=width==0 ?
				MagQCNetworkHarness.Scorer.loadConfigured(options, false) :
				new MagQCNetworkHarness.Scorer(options, width, false);
			if(loaded.dummy || loaded.worker.getTag("magqc_fixture")!=null){
				throw new IllegalArgumentException("The synchronized API requires a real six-output composite");
			}
			return loaded;
		}catch(RuntimeException e){throw e;}
		catch(Exception e){throw new IllegalArgumentException("Could not load the bound composite", e);}
	}

	/**
	 * Copies caller inputs and returns six private outputs. The instance lock also
	 * protects direct callers; ProkCC additionally serializes binding and first load.
	 * Callers must not modify their input concurrently with this call.
	 */
	public synchronized float[] score(float[] input){return scorer.evaluateCopy(input);}

	/** Resolves the release next to this installation, without loading a network. */
	public static Binding shippedDefault(int width){
		final String config=Paths.get(Data.NETWORKS(), "magqc", "release.config").toString();
		return new Binding(MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}), width);
	}

	/** Immutable configured model/resource identity; no caller-owned map is retained. */
	public static final class Binding {
		/** Takes the native release options and the composite input width. */
		public Binding(Map<String,String> source, int width_){
			if(source==null || width_<1){throw new IllegalArgumentException("Composite binding requires options and positive width");}
			width=width_; options=new HashMap<String,String>();
			for(String key:RESOURCES){
				final String path=source.get(key), pin=source.get(key+"sha80");
				if(path==null || path.isEmpty() || pin==null || !pin.matches("[0-9a-f]{20}")){
					throw new IllegalArgumentException("Composite binding requires path and sha80 for "+key);
				}
				options.put(key, Paths.get(path).toAbsolutePath().normalize().toString());
				options.put(key+"sha80", pin);
			}
		}

		/** The expected caller vector width, checked again against the parsed model. */
		public int inputs(){return width;}

		/** Model plus all frozen resource paths/pins must agree before cache reuse. */
		@Override public boolean equals(Object other){
			if(this==other){return true;}
			if(!(other instanceof Binding)){return false;}
			final Binding binding=(Binding)other;
			return width==binding.width && options.equals(binding.options);
		}
		@Override public int hashCode(){return 31*width+options.hashCode();}

		private final HashMap<String,String> options;
		private final int width;
	}

	private final MagQCNetworkHarness.Scorer scorer;
	private final Binding binding;
	private static final String[] RESOURCES={"net", "bundle", "familylist", "subnetmanifest",
		"expectedcopytable", "subnetpopulations"};
}
