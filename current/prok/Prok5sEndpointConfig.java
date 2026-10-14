package prok;

import java.io.File;
import java.util.ArrayList;
import dna.Data;
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
		final String tableDir=tables==null ? new File(resource("5S_endpoint_tables/TABLES.tsv")).getParent() : tables;
		final String netPath=nets==null ? resource("5S_endpoints.bbnets.gz") : nets;
		require(selected.reuseConsensusAlignment && !selected.trimAlignmentExtent && !selected.voteEnds && !selected.voteWindows
			&& selected.boundary5NetTemplate==null && selected.boundary3NetTemplate==null && selected.boundaryNetsByModel==null,
			"5S endpoint features require raw accepted alignment endpoints and cannot combine with other boundary policies");
		require(new File(tableDir,"TABLES.tsv").isFile() && (new File(netPath).isFile() || new File(netPath,"MANIFEST.tsv").isFile()),
			"5S endpoint inference needs a complete table manifest and selected-network directory or bundle");
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load("prok5S",tableDir,selected.modelNames,selected.library,new String[0]);
		final RrnaEndpointCallerFeatures.Resources inference=RrnaEndpointNetworkLoader.load(netPath,loaded);
		require(inference.acceptedVerifierIdentity,"5S nets must preserve the trained accepted-verifier-window identity feature");
		selected.setRrnaEndpointInference(inference,null);
		System.err.println("prok5S endpoints: enabled; tables="+tableDir+" nets="+netPath+" radius="+inference.radius);
	}
	private static String resource(String name){
		final String path=Data.findPath("?"+name,false);
		require(path!=null && new File(path).isFile(),"Missing5S endpoint resource: "+name);return path;
	}
	private static String path(String key,String value){require(value!=null && !value.isEmpty(),key+" requires an explicit resource path");return value;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private boolean infer=true,explicit=false;
	private String tables,nets;
}
