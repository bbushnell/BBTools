package prok;

import java.io.File;
import java.util.ArrayList;
import dna.Data;
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
		final String tableDir=tables==null ? new File(resource("16S_endpoint_tables/TABLES.tsv")).getParent() : tables;
		final String netPath=nets==null ? resource("16S_endpoints.bbnets.gz") : nets;
		require(new File(tableDir,"TABLES.tsv").isFile(),"Missing16S endpoint table manifest: "+tableDir);
		require(new File(netPath).isFile() || new File(netPath,"MANIFEST.tsv").isFile(),"Missing16S endpoint network file/manifest: "+netPath);
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load("16S",tableDir,
			selected.modelNames,selected.library,new String[]{"mitochondria"});
		final RrnaEndpointCallerFeatures.Resources inference=RrnaEndpointNetworkLoader.load(netPath,loaded);
		require(inference.acceptedVerifierIdentity,"16S endpoint networks must declare the calibrated accepted-verifier identity feature");
		selected.setRrnaEndpointInference(inference,null);
		System.err.println("16S endpoints: enabled; mitochondrial model retains raw ends; tables="+tableDir+" nets="+netPath);
	}
	private static String resource(String name){
		final String path=Data.findPath("?"+name,false);
		require(path!=null && new File(path).isFile(),"Missing16S endpoint resource: "+name);return path;
	}
	private static String path(String key,String value){
		require(value!=null && !value.isEmpty(),key+" requires a resource path");return value;
	}
	private static void require(boolean value,String why){if(!value){throw new IllegalArgumentException(why);}}
	private boolean infer=true,explicit=false;
	private String tables=null,nets=null;
}
