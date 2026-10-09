package prok;

import java.io.File;
import fileIO.ByteFile;
import map.ObjectIntMap;
import parse.LineParser1;

/** Explicit per-model endpoint table loader. Launchers verify the manifest and
 * payload sha80s; this loader validates the biological model/end/geometry binding.
 * Ordinary CallGenes installs these only after the corresponding family is enabled.
 * @author Raiden
 */
public final class RrnaEndpointResourceLoader {
	/** Existing euk5S callers retain the published three-argument entry point. */
	public static RrnaEndpointCallerFeatures.Resources load(String directory,String[] names,byte[][] refs){
		return load(RrnaEndpointCallerFeatures.Resources.legacyFamily(names),directory,names,refs);
	}
	public static RrnaEndpointCallerFeatures.Resources load(String family,String directory,String[] names,byte[][] refs){
		require(directory!=null&&names!=null&&refs!=null&&names.length==refs.length&&names.length>0,"Explicit table directory and caller library required");
		final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);final int n=names.length;
		for(int i=0;i<n;i++){require(names[i]!=null&&names[i].matches("[A-Za-z0-9_]+")&&!index.contains(names[i]),"Unique safe caller model names must match the table bundle");index.put(names[i],i);}
		final RrnaPositionalKmerTable.Table[] five=new RrnaPositionalKmerTable.Table[n],three=new RrnaPositionalKmerTable.Table[n];
		final String manifest=new File(directory,"TABLES.tsv").toString();final ByteFile in=RrnaResourceIO.open(manifest);final LineParser1 p=new LineParser1('\t');int rows=0;
		try{RrnaResourceIO.header(in,HEADER,manifest);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==8,"Table manifest must retain all geometry/provenance columns");final String model=p.parseString(0),end=p.parseString(1);require(index.contains(model),"Table model is absent from the caller library: "+model);
				final int i=index.get(model);final RrnaPositionalKmerTable.Table[] dest=end.equals("5prime")?five:end.equals("3prime")?three:null;
				require(dest!=null&&dest[i]==null&&p.parseInt(3)>0&&p.parseString(7).matches("[0-9a-f]{20}"),"Exactly one sha80-bound table per model and end");
				final RrnaPositionalKmerTable.Table t=RrnaPositionalKmerTable.load(new File(directory,model+'.'+end+".tsv.gz").toString(),model);
				require(t.radius==p.parseInt(3)&&t.end.equals(end)&&t.k==p.parseInt(2)&&t.offset==p.parseInt(4)&&t.alpha==RrnaPositionalKmerTable.savedDouble(p,5)&&t.records==p.parseLong(6),"Loaded table differs from its declared model/end/k/radius/offset/smoothing/population");dest[i]=t;rows++;
			}require(rows==2*n,"Every caller model requires both endpoint tables");
		}finally{require(!in.close(),"Endpoint table manifest read failed");}
		return new RrnaEndpointCallerFeatures.Resources(family,names,refs,five,three);
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final String HEADER="model\tend\tk\tradius\toffset\tpseudocount\tmembers\ttable_sha80";
}
