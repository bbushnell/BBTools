package prok;

import java.util.ArrayList;
import parse.Parse;

/** Run-local23S endpoint inference; disabled families and explicit off do no resource I/O.
 * @author Brian, Raiden */
final class Prok23sEndpointConfig {
	boolean parse(String key,String value){
		if(key.equalsIgnoreCase("rrna23sendpointnn")){infer=Parse.parseBoolean(value);explicit=true;}
		else if(key.equalsIgnoreCase("rrna23sendpointtables")){tables=path(key,value);}
		else if(key.equalsIgnoreCase("rrna23sendpointnets")){nets=path(key,value);}
		else{return false;}return true;
	}
	void apply(ArrayList<NcrnaFamily> families){
		require(families!=null,"Endpoint controls require initialized family registration");NcrnaFamily selected=null;
		for(NcrnaFamily f:families){if(f.name.equals("23S")){require(selected==null,"Only one23S family may own endpoint resources");selected=f;}}
		if(selected==null){require((!explicit || !infer) && tables==null && nets==null,"23S endpoint controls require rrna17=23s or rrna17=prok");return;}
		require(infer || (tables==null && nets==null),"Explicit resources require rrna23sendpointnn=t");if(!infer){return;}
		require(CallGenes.RRNA23S_ENDPOINT_MODE==CallGenes.RRNA23S_ENDPOINT_OFF && selected.boundary5NetTemplate==null && selected.boundary3NetTemplate==null && selected.boundaryNetsByModel==null && !selected.voteEnds && !selected.trimAlignmentExtent,
			"Current-recipe per-end nets must not combine with historical23S voting/scoring/trim modes");
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.tables=tables;overrides.endpoints=nets;
		final RrnaEndpointCallerFeatures.Resources inference=NcrnaFamilyLoader.loadEndpoints(selected, overrides, resolver);
		selected.setRrnaEndpointInference(inference,null);
		System.err.println("23S endpoints: enabled; tables="+resolver.tables(selected.name, overrides).path
			+" nets="+resolver.resolve(selected.name, "endpoints", overrides).path+" radius="+inference.radius);
	}
	private static String path(String key,String value){require(value!=null && !value.isEmpty(),key+" requires an explicit resource path");return value;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private boolean infer=true,explicit=false;private String tables,nets;
}
