package prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectSet;
import parse.Parse;
import prok.NcrnaTextIO.ConfigInput;

/** Defaults declared by one family's boundary_config manifest role.
 * Reads plain key=value records without changing global command-line state.
 * Explicit loader choices win. Payload existence and library normalization matter
 * only when refinement is selected; the ordinary resource loaders validate payloads.
 * @author Brian Bushnell, Raiden
 */
public final class PackagedBoundaryConfig{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads one complete packaged configuration; unknown or duplicate fields fail loudly. */
	public static PackagedBoundaryConfig read(final String path){return new PackagedBoundaryConfig(path);}

	private PackagedBoundaryConfig(final String path){
		if(path==null){throw new IllegalArgumentException("Packaged boundary configuration needs a path");}
		file=new File(path).getAbsoluteFile();
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(ConfigInput input=new ConfigInput(path)){
			while(input.next()){
				final String a=input.key.toLowerCase(Locale.ROOT), b=input.value;
				require(!a.isEmpty() && !b.isEmpty(), "Empty configuration key or value: "+a);
				require(seen.add(a), "Duplicate key: "+a);
				parse(a, b);
			}
		}
		for(String key:REQUIRED){require(seen.contains(key), "Missing key: "+key);}
		validate();
	}

	/** Configuration is small startup state; strings survive only for family/path bindings. */
	private void parse(final String a, final String b){
		assert(a!=null && b!=null) : "Record framing must finish before interpreting a boundary setting";
		switch(a){
			case "family":
				family=CallGenes.parseNcrnaFamily(b);
				break;
			case "enabled":
				require(Parse.isBoolean(b) || b.equals("0") || b.equals("1"), "Enabled must be a boolean value");
				enabled=Parse.parseBoolean(b);
				break;
			case "features":
				require(b.equalsIgnoreCase("v1"), "Only feature version v1 is supported");
				break;
			case "featureorder":
				require(b.equals(FEATURE_ORDER), "Feature order does not match the shared scorer");
				break;
			case "meanlen":
				normalization=NcrnaNumbers.intValue(a, b);
				break;
			case "startradius":
				startRadius=NcrnaNumbers.intValue(a, b);
				break;
			case "stopradius":
				stopRadius=NcrnaNumbers.intValue(a, b);
				break;
			case "marginstart":
				startMargin=NcrnaNumbers.floatValue(a, b);
				break;
			case "marginstop":
				stopMargin=NcrnaNumbers.floatValue(a, b);
				break;
			case "cutoff":
				cutoff=b.equalsIgnoreCase("none") ? Float.NaN : NcrnaNumbers.floatValue(a, b);
				require(b.equalsIgnoreCase("none") || Float.isFinite(cutoff), "Numeric cutoff must be finite; none disables it");
				break;
			case "net":
				net=resolve(b);
				break;
			case "tablestart":
				startTable=resolve(b);
				break;
			case "tablestop":
				stopTable=resolve(b);
				break;
			default: throw new IllegalArgumentException("Unknown packaged boundary key "+a+" in "+file);
		}
	}

	private void validate(){
		NcrnaFamilyRegistry.sharedBoundaryMember(family);
		require(normalization>0 && startRadius>=0 && stopRadius>=0 && startRadius<=normalization && stopRadius<=normalization
			&& startRadius<(Integer.MAX_VALUE-1)/2 && stopRadius<(Integer.MAX_VALUE-1)/2,
			"Correction radii must fit the declared upper-median model length");
		require(Float.isFinite(startMargin) && startMargin>=0 && Float.isFinite(stopMargin) && stopMargin>=0,
			"Confidence margins must be finite and nonnegative");
		require(enabled ? net!=null && startTable!=null && stopTable!=null : net==null && startTable==null && stopTable==null,
			"ON requires a network and both tables; OFF must omit their bindings");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Defaults           ----------------*/
	/*--------------------------------------------------------------*/

	/** Applies only missing choices; boundaryOverride already resolves family/global CLI precedence. */
	public void applyDefaults(final String name, final NcrnaFamilyLoader.Options options, final byte[][] library){
		assert(options!=null) : "Defaults belong to an existing caller-owned options object";
		require(family.equals(name), "Configuration cannot alter a different family: "+name);
		options.boundary=options.boundaryOverride==null ? enabled : options.boundaryOverride;
		if(options.boundary){
			require(upperMedianLength(library)==normalization, "Packaged normalization differs from the loaded model library");
		}
		if(options.resources.boundary==null){options.resources.boundary=net;}
		if(options.resources.boundaryStart==null){options.resources.boundaryStart=startTable;}
		if(options.resources.boundaryStop==null){options.resources.boundaryStop=stopTable;}
		if(options.boundaryStartOffsets==null){options.boundaryStartOffsets=offsets(startRadius);}
		if(options.boundaryStopOffsets==null){options.boundaryStopOffsets=offsets(stopRadius);}
		if(Float.isNaN(options.boundaryMarginStart)){options.boundaryMarginStart=startMargin;}
		if(Float.isNaN(options.boundaryMarginStop)){options.boundaryMarginStop=stopMargin;}
	}

	/** Binds the score cutoff after the native loader establishes the selected scorer's geometry. */
	void loaded(final NcrnaFamily loaded){
		require(family.equals(loaded.name), "Loaded family differs from its packaged configuration");
		if(loaded.boundary5NetTemplate==null){
			require(loaded.boundary3NetTemplate==null, "A shared adapter cannot have only one enabled end");
			return;
		}
		require(loaded.boundaryFeatureVersion==NcrnaBoundaryScorer.FEATURES_V1 && loaded.boundaryNetsByModel==null,
			"Packaged v1 settings cannot reinterpret another feature layout or per-model networks");
		require(loaded.boundaryMeanLen==normalization, "Loaded normalization differs from the package");
		require(loaded.boundary3NetTemplate!=null && loaded.boundary5NetTemplate.numOutputs()==1
			&& loaded.boundary3NetTemplate.numOutputs()==1, "Boundary ranking requires one scalar output per end");
		loaded.nnCutoff=cutoff;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	/** Upper median, matching the shared scorer's length normalization contract. */
	public static int upperMedianLength(final byte[][] library){
		if(library==null || library.length==0){throw new IllegalArgumentException("Cannot normalize an empty model library");}
		final int[] lengths=new int[library.length];
		for(int i=0; i<library.length; i++){
			if(library[i]==null || library[i].length==0){throw new IllegalArgumentException("Empty model at index "+i);}
			lengths[i]=library[i].length;
		}
		Arrays.sort(lengths);
		return lengths[lengths.length/2];
	}

	private static int[] offsets(final int radius){
		assert(radius>=0 && radius<(Integer.MAX_VALUE-1)/2) : "Validated radius must fit a centered offset array";
		final int[] offsets=new int[2*radius+1];
		for(int i=0; i<offsets.length; i++){offsets[i]=i-radius;}
		return offsets;
	}

	private String resolve(final String path){
		assert(path!=null && !path.isEmpty()) : "Resource records require a nonempty value";
		final File resource=new File(path);
		return (resource.isAbsolute() ? resource : new File(file.getParentFile(), path)).getAbsolutePath();
	}

	private void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message+" in "+file);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final File file;
	private String family, net, startTable, stopTable;
	private boolean enabled;
	private int normalization, startRadius, stopRadius;
	private float startMargin, stopMargin, cutoff;
	/** Serialized v1 feature order; must match the shared scorer's ten network inputs. */
	public static final String FEATURE_ORDER="ani,prof0,prof1,prof2,isStop,fuzz0,fuzz1,fuzz2,lengthRatio,contigGC";
	private static final String[] REQUIRED={"family", "enabled", "features", "featureorder", "meanlen",
		"startradius", "stopradius", "marginstart", "marginstop", "cutoff"};
}
