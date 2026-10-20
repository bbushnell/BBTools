package prok;

import parse.Parse;

/** Per-subtype 6S boundary overrides; resource loading remains in the family loader.
 * @author Ganyu */
public final class SixsBoundaryConfig {

	public SixsBoundaryConfig(String familyId){
		if(!"sixs_rf00013".equals(familyId) && !"sixs_rf01685".equals(familyId)){
			throw new IllegalArgumentException("6S boundary configuration requires sixs_rf00013 or sixs_rf01685: "+familyId);
		}
		family=familyId;prefix=familyId.replace("_","");
	}

	/** Returns false for another subtype's keys, without changing this instance. */
	public boolean parse(String key, String value){
		require(key!=null, "A boundary option requires a key");
		if(key.equalsIgnoreCase(prefix+"boundarynn")){enabled=Parse.parseBoolean(value);return true;}
		if(key.equalsIgnoreCase(prefix+"boundarynet")){net=path(value);return true;}
		if(key.equalsIgnoreCase(prefix+"boundarystarttable")){startTable=path(value);return true;}
		if(key.equalsIgnoreCase(prefix+"boundarystoptable")){stopTable=path(value);return true;}
		if(key.equalsIgnoreCase(prefix+"boundarystartrange")){start=R58BoundaryConfig.range(value);return true;}
		if(key.equalsIgnoreCase(prefix+"boundarystoprange")){stop=R58BoundaryConfig.range(value);return true;}
		return false;
	}

	public boolean enabled(boolean global){return enabled==null?global:enabled.booleanValue();}

	/** Validate after all arguments, so resource and enable order is immaterial. */
	public void validate(boolean ncrna, boolean sixs, boolean global){
		require(!(Boolean.TRUE.equals(enabled) || hasResources()) || ncrna && sixs,
			"Explicit "+family+" boundary enable/paths/ranges require ncrna=t and sixs=t");
		validateResources(global);
	}

	/** Changes only this subtype's boundary selections, retaining other loader options. */
	public void apply(String familyId, NcrnaFamilyLoader.Options options, boolean global){
		require(family.equals(familyId), "6S boundary controls for "+family+" cannot alter "+familyId);
		require(options!=null, "Boundary overrides require the target family's loader options");
		validateResources(global);
		options.boundary=enabled(global);
		if(net!=null){options.resources.boundary=net;}
		if(startTable!=null){options.resources.boundaryStart=startTable;}
		if(stopTable!=null){options.resources.boundaryStop=stopTable;}
		if(start!=null){options.boundaryStartOffsets=start.clone();}
		if(stop!=null){options.boundaryStopOffsets=stop.clone();}
	}

	public void clear(){enabled=null;net=startTable=stopTable=null;start=stop=null;}

	private void validateResources(boolean global){
		require(enabled(global) || !hasResources(),
			"Explicit "+family+" boundary paths/ranges require effective boundary refinement enabled");
	}
	private boolean hasResources(){return net!=null || startTable!=null || stopTable!=null || start!=null || stop!=null;}
	private static String path(String value){
		require(value!=null && !value.trim().isEmpty(), "Explicit boundary resource path must be nonempty");
		return value;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	private final String family,prefix;
	private Boolean enabled;
	private String net,startTable,stopTable;
	private int[] start,stop;
}
