package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Portable caller envelope around existing position annotations. Exact consensus
 * bytes and full key conservation bind the resource without a runtime hash tool
 * or any CM calculation. The inner version1 table remains unchanged.
 * @author Brian Bushnell, Raiden */
public final class SeedModelPositionResource {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,SeedModelPositionResource.class,false);
		final ObjectMap<String,String> flags=new ObjectMap<String,String>(String.class,String.class);
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT);
			require(KNOWN.contains("|"+key+"|") && flags.put(key,arg.substring(eq+1))==null,"Unknown or repeated binding: "+key);}
		for(String key:KNOWN.substring(1,KNOWN.length()-1).split("\\|")){require(flags.get(key)!=null,"Missing resource binding: "+key);}
		Shared.setThreads(1);final ArrayList<String> names=new ArrayList<String>();final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		RrnaResourceIO.stream(flags.get("models"),r->{names.add(RrnaResourceIO.first(r.id));sequences.add(r.bases.clone());});
		final String[] ids=names.toArray(new String[0]);final byte[][] bases=sequences.toArray(new byte[0][]);final int[] lengths=lengths(bases);
		final LongHashSet seeds=ProkObject.loadLongKmers(flags.get("seeds"),17);final int clen=Integer.parseInt(flags.get("clen"));final int[] sites=sites(flags.get("sites"));
		final String[] hashes=new String[4];for(int i=0;i<hashes.length;i++){hashes[i]=flags.get(SeedModelPositionTable.BINDINGS[i].toLowerCase(Locale.ROOT));}
		final SeedModelPositionTable table=SeedModelPositionTable.load(flags.get("positions"),seeds,ids,lengths,clen,sites,hashes);
		write(flags.get("out"),table,bases);
		final SeedModelPositionTable reread=load(flags.get("out"),seeds,ids,bases);
		require(reread.seedCount()==seeds.size(),"Round-trip must retain every mapped and unmapped caller key");
		System.err.println("SEED_POSITION_RESOURCE_PASS models="+ids.length+" seeds="+seeds.size()+" unchanged_annotations=true exact_consensus_binding=true");Shared.closeStream(pp.outstream);
	}
	static SeedModelPositionTable load(String path,LongHashSet seeds,String[] names,byte[][] sequences){
		require(names!=null && sequences!=null && names.length==sequences.length,"Runtime model names and complete sequences must stay index-aligned");
		final int[] lengths=lengths(sequences);final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			final byte[] header=in.nextLine();require(header!=null,"Missing bound resource header");p.set(header);
			require(p.terms()==8 && p.termEquals("#SeedModelPositionResource",0) && p.parseInt(1)==1,"Expected the version1 exact-consensus envelope");
			final int clen=p.parseInt(2);final int[] sites=sites(p.parseString(3));final String[] hashes=new String[4];
			for(int i=0;i<hashes.length;i++){hashes[i]=p.parseString(i+4);require(SeedOffsetTable.validSha80(hashes[i]),"Resource provenance must use sha80");}
			for(int m=0;m<names.length;m++){
				final byte[] row=in.nextLine();require(row!=null,"Every actual consensus needs its complete bound sequence");p.set(row);
				require(p.terms()==3 && p.termEquals("#Sequence",0) && p.termEquals(names[m],1),"Bound consensus names must retain runtime index order");
				final byte[] saved=p.parseByteArray(2);require(java.util.Arrays.equals(saved,sequences[m]),"Seed-position map belongs to a different consensus sequence: "+names[m]);
			}
			final SeedModelPositionTable table=SeedModelPositionTable.read(in,path,seeds,names,lengths,clen,sites,hashes);
			validatePositions(table,sequences);return table;
		}finally{require(!in.close(),"Bound position resource read failed");}
	}
	static void write(String path,SeedModelPositionTable table,byte[][] sequences){
		require(table!=null && sequences!=null && table.modelCount()==sequences.length,"Envelope preserves the existing annotation model order");
		RrnaResourceIO.fresh(path,"");final ByteStreamWriter out=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{
			b.append("#SeedModelPositionResource\t1\t").append(table.clen).tab();for(int s=0;s<table.sites.length;s++){if(s>0){b.append(',');}b.append(table.sites[s]);}
			for(String hash:table.hashes){b.tab().append(hash);}out.print(b.nl());
			for(int m=0;m<table.modelCount();m++){
				require(sequences[m].length==table.modelLength(m),"Bound sequence length must agree with annotations");
				out.print(b.clear().append("#Sequence\t").append(table.modelName(m)).tab().append(sequences[m]).nl());
			}
			table.write(out);
		}finally{RrnaResourceIO.finish(out);}
	}
	static int[] lengths(byte[][] sequences){require(sequences!=null && sequences.length>0,"A position resource must bind a nonempty library");final int[] lengths=new int[sequences.length];for(int i=0;i<lengths.length;i++){require(sequences[i]!=null && sequences[i].length>=17,"Consensus cannot be shorter than its17-mer anchors");lengths[i]=sequences[i].length;}return lengths;}
	static void validatePositions(SeedModelPositionTable table,byte[][] sequences){
		for(long key:table.keys){final SeedModelPositionTable.Entry entry=table.get(key);for(int i=0;i<entry.size();i++){
			final SeedModelPositionTable.Position p=entry.get(i);long word=0;for(int j=0;j<17;j++){final int base=RrnaResourceIO.base(sequences[p.model][p.start+j]);require(base>=0,"An exact mapped seed cannot contain an ambiguous consensus base");word=(word<<2)|base;}
			require(word==key,"Every annotated seed must equal the actual caller consensus bases at its declared position");
		}}
	}
	static int[] sites(String value){final String[] words=value.split(",");final int[] sites=new int[words.length];for(int i=0;i<sites.length;i++){sites[i]=Integer.parseInt(words[i]);}return sites;}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final String KNOWN="|positions|models|seeds|clen|sites|seedfastasha80|consensussha80|anchorssha80|cmsha80|out|";
}
