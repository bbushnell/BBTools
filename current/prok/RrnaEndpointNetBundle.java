package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import map.ObjectSet;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** One line-count-framed .bbnets file retaining model/end/seed and identity semantics.
 * Numerical CellNet lines are copied, not reserialized or rounded. Only a recognized
 * legacy checkpoint digest is normalized to sha80. Member sha80 fields identify the
 * original source files, not the normalized framed payloads. Source sha80 fields
 * are provenance; the producer verifies input and output payload hashes separately.
 * @author Brian, Raiden */
public final class RrnaEndpointNetBundle {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaEndpointNetBundle.class,false);String in=null,out=null;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Use flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("out")){out=v;}else{throw new IllegalArgumentException(k);}}
		require(in!=null && out!=null,"Explicit selected-network directory and fresh bundle output required");
		pack(in,out);Shared.closeStream(pp.outstream);
	}
	static void pack(String directory,String output){
		require(RrnaEndpointNetworkLoader.acceptedIdentity(directory,"16S"),"This release bundle requires accepted-verifier identity metadata");
		final ByteFile in=RrnaResourceIO.open(new File(directory,"MANIFEST.tsv").toString());
		final LineParser1 p=new LineParser1('\t');final ArrayList<Entry> entries=new ArrayList<Entry>();final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,"model\tend\tseed\tfile\tsha80",directory);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==5,"Selected network manifest must retain model/end/seed/file/sha80");
				final String model=p.parseString(0),end=p.parseString(1),file=p.parseString(3),digest=p.parseString(4);final int seed=p.parseInt(2);
				validateEntry(model,end,seed,digest);require(seen.add(model+'\t'+end),"Duplicate selected network member");
				require(file.equals("networks/"+model+'.'+end+".s"+seed+".bbnet"),"Selected network must resolve within its declared bundle directory");
				final ByteFile net=RrnaResourceIO.open(new File(directory,file).toString());final ArrayList<byte[]> lines=new ArrayList<byte[]>();
				try{for(byte[] line;(line=net.nextLine())!=null;){lines.add(normalizeLegacyProvenance(line));}}finally{require(!net.close(),"Network payload read failed");}
				require(!lines.isEmpty(),"Cannot package an empty network");entries.add(new Entry(model,end,seed,digest,lines));
			}
		}finally{require(!in.close(),"Selected network manifest read failed");}
		require(!entries.isEmpty(),"A bundle needs selected network members");RrnaResourceIO.fresh(output,"");
		final ByteStreamWriter out=RrnaResourceIO.writer(output);final ByteBuilder b=new ByteBuilder();
		try{out.println(SCHEMA);out.println("family\t16S");out.println("identity_source\taccepted_verifier_window");out.print(b.clear().append("networks\t").append(entries.size()).nl());
			for(Entry entry:entries){
				out.print(b.clear().append("network\t").append(entry.model).tab().append(entry.end).tab().append(entry.seed).tab().append(entry.sha80).tab().append(entry.lines.size()).nl());
				for(byte[] line:entry.lines){out.print(b.clear().append(line).nl());}
			}out.println(END);
		}finally{RrnaResourceIO.finish(out);}
		System.err.println("RRNA_ENDPOINT_NET_BUNDLE_PACK_PASS members="+entries.size()+" identity=accepted_verifier_window");
	}
	static RrnaEndpointCallerFeatures.Resources load(String path,RrnaEndpointCallerFeatures.Resources tables){
		require(tables!=null && "16S".equals(tables.family),"This endpoint bundle must bind to the declared 16S tables");
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);int active=0;
		for(int i=0; i<tables.names.length; i++){indices.put(tables.names[i],i);if(!tables.raw[i]){active++;}}
		final CellNet[] five=new CellNet[tables.names.length],three=new CellNet[tables.names.length];
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,SCHEMA,path);RrnaResourceIO.header(in,"family\t16S",path);RrnaResourceIO.header(in,"identity_source\taccepted_verifier_window",path);
			p.set(in.nextLine());require(p.terms()==2 && p.termEquals("networks",0),"Bundle must declare its member count");
			final int count=p.parseInt(1);require(count==active*2,"Every refined model needs exactly two bundled networks");
			for(int member=0; member<count; member++){
				p.set(in.nextLine());require(p.terms()==6 && p.termEquals("network",0),"Each network requires framed model/end/seed/sha80/line-count metadata");
				final String model=p.parseString(1),end=p.parseString(2),digest=p.parseString(4);final int seed=p.parseInt(3),size=p.parseInt(5);
				validateEntry(model,end,seed,digest);require(indices.contains(model) && size>0 && size<=1000000,"Network must name a caller model and a bounded nonempty payload");
				final int index=indices.get(model);final CellNet[] dest=end.equals("5prime") ? five : three;
				require(!tables.raw[index] && dest[index]==null,"Raw or duplicate model/end members are not valid endpoint networks");
				final ArrayList<byte[]> lines=new ArrayList<byte[]>(size);
				for(int i=0; i<size; i++){final byte[] line=in.nextLine();require(line!=null,"Truncated framed network payload");lines.add(line);}
				final CellNet net=CellNetParser.loadFromLines(lines);RrnaEndpointCallerFeatures.Resources.validateNet(net,tables.inputs,tables.sites);dest[index]=net;
			}
			RrnaResourceIO.header(in,END,path);require(in.nextLine()==null,"Unexpected trailing bundle content");
		}finally{require(!in.close(),"Endpoint bundle read failed");}
		return new RrnaEndpointCallerFeatures.Resources("16S",tables.names,tables.refs,tables.five,tables.three,five,three,tables.raw,true);
	}
	static void validateEntry(String model,String end,int seed,String digest){
		require(model.matches("[A-Za-z0-9_.:-]+") && (end.equals("5prime") || end.equals("3prime")) && seed>=0 && digest.matches("[0-9a-f]{20}"),"Network metadata must preserve an unambiguous model/end/seed and sha80");
	}
	/** Version1 migration of the known GPU export tag; all numerical bytes survive. */
	static byte[] normalizeLegacyProvenance(byte[] line){
		assert(line!=null):"A framed member must contain complete source lines";
		final String prefix="##checkpoint_sha256 ";
		if(!shared.Tools.startsWith(line,prefix)){return line;}
		final String digest=new String(line,prefix.length(),line.length-prefix.length(),java.nio.charset.StandardCharsets.US_ASCII);
		require(digest.matches("[0-9a-f]{64}"),"Legacy checkpoint metadata must be exactly64 lowercase hex digits before normalization");
		return new ByteBuilder("##checkpoint_sha80 ").append(digest.substring(44)).toBytes();
	}
	static void require(boolean value,String why){if(!value){throw new IllegalArgumentException(why);}}
	static final class Entry{
		Entry(String m,String e,int s,String h,ArrayList<byte[]> l){model=m;end=e;seed=s;sha80=h;lines=l;}
		final String model,end,sha80;final int seed;final ArrayList<byte[]> lines;
	}
	static final String SCHEMA="#RRNA_ENDPOINT_BBNETS_V1",END="#END_RRNA_ENDPOINT_BBNETS";
}
