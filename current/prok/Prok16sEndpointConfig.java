package prok;

import java.util.ArrayList;
import parse.Parse;

/** Run-local loading of the measured16S endpoint resources; disabled families do no I/O.
 * @author Brian, Raiden */
final class Prok16sEndpointConfig {
	boolean parse(String key,String value){
		if(key.equalsIgnoreCase("rrna16sendpoint")){infer=Parse.parseBoolean(value);explicit=true;}
		else if(key.equalsIgnoreCase("rrna16sendpointtables")){tables=path(key,value);}
		else if(key.equalsIgnoreCase("rrna16sendpointnets")){nets=path(key,value);}
		else{return false;}
		return true;
	}
	void apply(ArrayList<NcrnaFamily> families){
		require(families!=null,"Endpoint configuration requires the registered caller families");
		NcrnaFamily selected=null;
		for(NcrnaFamily family:families){if(family.name.equals("16S")){
			require(selected==null,"Only one16S family may own endpoint resources");selected=family;
		}}
		if(selected==null){
			require((!explicit || !infer) && tables==null && nets==null,"16S endpoint options require rrna17=16s or rrna17=prok");
			return;
		}
		require(infer || (tables==null && nets==null),"Endpoint resource overrides require rrna16sendpoint=t");
		if(!infer){System.err.println("16S endpoints: disabled");return;}
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.tables=tables;overrides.endpoints=nets;
		final RrnaEndpointCallerFeatures.Resources inference=NcrnaFamilyLoader.loadEndpoints(selected, overrides, resolver);
		selected.setRrnaEndpointInference(inference,null);
		System.err.println("16S endpoints: enabled; mitochondrial model retains raw ends; tables="+resolver.tables(selected.name, overrides).path
			+" nets="+resolver.resolve(selected.name, "endpoints", overrides).path);
	}
	private static String path(String key,String value){
		require(value!=null && !value.isEmpty(),key+" requires a resource path");return value;
	}
	private static void require(boolean value,String why){if(!value){throw new IllegalArgumentException(why);}}
	private boolean infer=true,explicit=false;
	private String tables=null,nets=null;
}
