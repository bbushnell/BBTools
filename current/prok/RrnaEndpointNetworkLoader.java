package prok;

import java.io.File;
import fileIO.ByteFile;
import map.ObjectIntMap;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;

/** Bind explicit selected networks to actual caller models. Producers verify all
 * sha80 payloads; this loader rejects incomplete or ambiguous biological mappings.
 * @author Raiden
 */
public final class RrnaEndpointNetworkLoader {
	public static RrnaEndpointCallerFeatures.Resources load(String directory,RrnaEndpointCallerFeatures.Resources tables){
		require(directory!=null && tables!=null,"Explicit network directory and verified table binding required");
		final int n=tables.names.length;final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);
		for(int i=0;i<n;i++){index.put(tables.names[i],i);}
		final CellNet[] five=new CellNet[n],three=new CellNet[n];final String manifest=new File(directory,"MANIFEST.tsv").toString();
		final ByteFile in=RrnaResourceIO.open(manifest);final LineParser1 p=new LineParser1('\t');int rows=0;
		try{RrnaResourceIO.header(in,"model\tend\tseed\tfile\tsha80",manifest);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==5,"Network manifest requires model/end/seed/file/sha80");
				final String model=p.parseString(0),end=p.parseString(1),path=p.parseString(3);require(index.contains(model),"Network model absent from caller library: "+model);
				final int i=index.get(model),seed=p.parseInt(2);final CellNet[] dest=end.equals("5prime")?five:end.equals("3prime")?three:null;
				require(dest!=null && dest[i]==null && seed>=0 && p.parseString(4).matches("[0-9a-f]{20}"),"Exactly one selected network per model/end is required");
				require(path.equals("networks/"+model+'.'+end+".s"+seed+".bbnet"),"Network filename must preserve its declared model/end/seed within the bundle");
				final CellNet net=CellNetParser.load(new File(directory,path).toString());RrnaEndpointCallerFeatures.Resources.validateNet(net,tables.inputs,tables.sites);dest[i]=net;rows++;
			}require(rows==2*n,"Every model must have both selected endpoint networks");
		}finally{require(!in.close(),"Network manifest read failed");}
		return new RrnaEndpointCallerFeatures.Resources(tables.family,tables.names,tables.refs,tables.five,tables.three,five,three);
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
