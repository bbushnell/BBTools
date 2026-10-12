package prok;

import java.io.File;
import java.util.ArrayList;
import dna.Data;
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
		final String tableDir=tables==null ? new File(resource("23S_endpoint_tables/TABLES.tsv")).getParent() : tables;
		final String netPath=nets==null ? resource("23S_endpoints.bbnets.gz") : nets;
		require(new File(tableDir,"TABLES.tsv").isFile() && (new File(netPath).isFile() || new File(netPath,"MANIFEST.tsv").isFile()),"23S endpoint resources require a table manifest and a selected-network directory or framed bundle file");
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load("23S",tableDir,selected.modelNames,selected.library,new String[0]);
		final RrnaEndpointCallerFeatures.Resources inference=RrnaEndpointNetworkLoader.load(netPath,loaded);
		require(inference.acceptedVerifierIdentity,"23S nets must explicitly carry accepted-verifier identity metadata");selected.setRrnaEndpointInference(inference,null);
		System.err.println("23S endpoints: enabled; tables="+tableDir+" nets="+netPath+" radius="+inference.five[0].radius);
	}
	private static String resource(String name){
		final String path=Data.findPath("?"+name,false);
		require(path!=null && new File(path).isFile(),"Missing23S endpoint resource: "+name);return path;
	}
	private static String path(String key,String value){require(value!=null && !value.isEmpty(),key+" requires an explicit resource path");return value;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private boolean infer=true,explicit=false;private String tables,nets;
}
