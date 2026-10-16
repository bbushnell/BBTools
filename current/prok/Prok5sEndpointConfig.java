package prok;

import java.util.ArrayList;
import parse.Parse;

/** Loads the measured prokaryotic5S endpoint resources for the seeded caller.
 * Disabled or PGM families perform no resource I/O; eukaryotic5S is separate.
 * @author Brian, Raiden */
final class Prok5sEndpointConfig {
	boolean parse(String key,String value){
		if(key.equalsIgnoreCase("rrna5sendpointnn")){infer=Parse.parseBoolean(value);explicit=true;}
		else if(key.equalsIgnoreCase("rrna5sendpointtables")){tables=path(key,value);}
		else if(key.equalsIgnoreCase("rrna5sendpointnets")){nets=path(key,value);}
		else{return false;}return true;
	}
	void apply(ArrayList<NcrnaFamily> families){
		require(families!=null,"Endpoint controls require initialized family registration");
		NcrnaFamily selected=null;
		for(NcrnaFamily family:families){if(family.name.equals("prok5S")){
			require(selected==null,"Only one seeded prok5S family may own endpoint resources");selected=family;
		}}
		if(selected==null){
			require((!explicit || !infer) && tables==null && nets==null,"5S endpoint options require the seeded prok5S family, not the PGM or euk5S family");
			return;
		}
		require(infer || (tables==null && nets==null),"Explicit5S endpoint resources require rrna5sendpointnn=t");
		if(!infer){System.err.println("prok5S endpoints: disabled");return;}
		require(selected.reuseConsensusAlignment && !selected.trimAlignmentExtent && !selected.voteEnds && !selected.voteWindows
			&& selected.boundary5NetTemplate==null && selected.boundary3NetTemplate==null && selected.boundaryNetsByModel==null,
			"5S endpoint features require raw accepted alignment endpoints and cannot combine with other boundary policies");
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.tables=tables;overrides.endpoints=nets;
		final RrnaEndpointCallerFeatures.Resources inference=NcrnaFamilyLoader.loadEndpoints(selected, overrides, resolver);
		selected.setRrnaEndpointInference(inference,null);
		System.err.println("prok5S endpoints: enabled; tables="+resolver.tables(selected.name, overrides).path
			+" nets="+resolver.resolve(selected.name, "endpoints", overrides).path+" radius="+inference.radius);
	}
	private static String path(String key,String value){require(value!=null && !value.isEmpty(),key+" requires an explicit resource path");return value;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private boolean infer=true,explicit=false;
	private String tables,nets;
}
