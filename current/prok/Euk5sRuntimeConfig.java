package prok;

import java.io.File;
import java.util.ArrayList;
import parse.Parse;

/** Ordinary CallGenes loading of selected euk5S cutoffs and endpoint resources.
 * The disabled family performs no resource I/O. Explicit scalar sweeps retain
 * their historical meaning unless an explicit model table is supplied.
 * @author Raiden
 */
final class Euk5sRuntimeConfig {
	boolean parse(String key,String value){
		if(key.equalsIgnoreCase("euk5sendpoint")){
			infer=Parse.parseBoolean(value);inferExplicit=true;
		}else if(key.equalsIgnoreCase("euk5sendpointtables")){
			tables=path(key,value);
		}else if(key.equalsIgnoreCase("euk5sendpointnets")){
			nets=path(key,value);
		}else if(key.equalsIgnoreCase("euk5smodelcutoffs")){
			require(value!=null && !value.isEmpty(),key+" requires a path or f");
			cutoffsDisabled=value.equalsIgnoreCase("f") || value.equalsIgnoreCase("false");
			cutoffs=cutoffsDisabled?null:value;
		}else{return false;}
		return true;
	}

	void apply(ArrayList<NcrnaFamily> families){
		require(families!=null,"Runtime configuration requires the registered caller families");
		NcrnaFamily selected=null;
		for(NcrnaFamily f:families){if(f.name.equals("euk5S")){
			require(selected==null,"Only one euk5S family may own runtime resources");selected=f;
		}}
		if(selected==null){
			require((!inferExplicit || !infer) && tables==null && nets==null && cutoffs==null,
				"euk5S runtime options require euk5s=t or an euk5S rrna17 profile");
			return;
		}
		require(infer || (tables==null && nets==null),"Endpoint table/network overrides require euk5sendpoint=t");
		final boolean scalar=CallGenes.ncrnaSweepTarget("euk5S")
			&& (!Float.isNaN(CallGenes.NCRNA_ID_PASS_OVERRIDE) || !Float.isNaN(CallGenes.NCRNA_ID_BORDERLINE_OVERRIDE));
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.cutoffs=cutoffs;overrides.tables=tables;overrides.endpoints=nets;
		if(!cutoffsDisabled && (cutoffs!=null || !scalar)){
			NcrnaFamilyLoader.loadCutoffs(selected,overrides,resolver);
		}
		if(infer){
			selected.setRrnaEndpointInference(NcrnaFamilyLoader.loadEndpoints(selected,overrides,resolver),null);
			final String tableDir=new File(resolver.tables(selected.name,overrides).path).getParent();
			final String netDir=resolver.resolve(selected.name,"endpoints",overrides).path;
			System.err.println("euk5S endpoints: enabled models="+selected.library.length+" tables="+tableDir+" nets="+netDir);
		}else{System.err.println("euk5S endpoints: disabled");}
	}
	private static String path(String key,String value){
		require(value!=null && !value.isEmpty(),key+" requires a resource path");return value;
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private boolean infer=true,cutoffsDisabled=false,inferExplicit=false;
	private String tables=null,nets=null,cutoffs=null;
}
