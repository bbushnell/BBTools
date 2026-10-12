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
		if(new File(directory).isFile()){return RrnaEndpointNetBundle.load(directory,tables);}
		final int n=tables.names.length;final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);
		int active=0;for(boolean raw:tables.raw){if(!raw){active++;}}
		for(int i=0;i<n;i++){index.put(tables.names[i],i);}
		final CellNet[] five=new CellNet[n],three=new CellNet[n];final String manifest=new File(directory,"MANIFEST.tsv").toString();
		final ByteFile in=RrnaResourceIO.open(manifest);final LineParser1 p=new LineParser1('\t');int rows=0;
		try{RrnaResourceIO.header(in,"model\tend\tseed\tfile\tsha80",manifest);
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==5,"Network manifest requires model/end/seed/file/sha80");
				final String model=p.parseString(0),end=p.parseString(1),path=p.parseString(3);require(index.contains(model),"Network model absent from caller library: "+model);
				final int i=index.get(model),seed=p.parseInt(2);final CellNet[] dest=end.equals("5prime")?five:end.equals("3prime")?three:null;
				require(!tables.raw[i],"An explicitly raw model must have no selected networks: "+model);
				require(dest!=null && dest[i]==null && seed>=0 && p.parseString(4).matches("[0-9a-f]{20}"),"Exactly one selected network per model/end is required");
				require(path.equals("networks/"+model+'.'+end+".s"+seed+".bbnet"),"Network filename must preserve its declared model/end/seed within the bundle");
				final CellNet net=CellNetParser.load(new File(directory,path).toString());RrnaEndpointCallerFeatures.Resources.validateNet(net,tables.inputs,tables.sites);dest[i]=net;rows++;
			}require(rows==2*active,"Every refined model must have both selected endpoint networks");
		}finally{require(!in.close(),"Network manifest read failed");}
		return new RrnaEndpointCallerFeatures.Resources(tables.family,tables.names,tables.refs,tables.five,tables.three,five,three,tables.raw,acceptedIdentity(directory,tables.family));
	}
	/** A new feature definition must travel with its trained network bundle.
	 * Missing metadata means the unchanged legacy candidate-span contract. */
	static boolean acceptedIdentity(String directory,String family){
		final File file=new File(directory,"IDENTITY_SOURCE.txt");if(!file.exists()){return false;}
		final ByteFile in=RrnaResourceIO.open(file.toString());try{final byte[] line=in.nextLine();require(line!=null && in.nextLine()==null,"Identity metadata must contain exactly one explicit feature definition");final String value=new String(line,java.nio.charset.StandardCharsets.US_ASCII);
			if(value.equals("candidate_span")){return false;}require(value.equals("accepted_verifier_window") && (family.equals("16S") || family.equals("23S")),"Only explicitly trained16S/23S bundles may select accepted-verifier identity;18S/euk5S stay unchanged");return true;
		}finally{require(!in.close(),"Network identity metadata read failed");}
	}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
