package prok;

import java.util.ArrayList;
import java.util.Arrays;
import aligner.CovarianceModelConsensusMap;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.LongObjectMap;
import parse.LineParser1;
import structures.ByteBuilder;

/** Versioned exact-consensus seed positions, with explicit unmapped keys.
 * All coordinates are boundaries or starts on the supplied sense-oriented
 * consensus, never pooled endpoint estimates. No alignment or CM scoring.
 * @author Brian Bushnell, Raiden */
public final class SeedModelPositionTable {
	/** Build from ALL existing caller keys and every bound consensus. A shared
	 * key can have positions in several models; repeats within one model remain
	 * separate postings and are marked ambiguous by the coverage report. */
	static SeedModelPositionTable build(LongHashSet seeds, String[] names, byte[][] sequences,
			CovarianceModelConsensusMap[] maps, int[] sites, String[] hashes){
		require(seeds!=null && seeds.size()>0 && names!=null && names.length>0 && sequences.length==names.length
			&& maps.length==names.length && maps[0]!=null, "Seed positions require the full nonempty seed set and index-aligned consensus maps");
		final int[] lengths=new int[names.length];final int clen=maps[0].modelLength();
		for(int m=0; m<names.length; m++){
			require(sequences[m]!=null && maps[m]!=null && sequences[m].length==maps[m].consensusLength && maps[m].modelLength()==clen,
				"Every model must use its actual complete consensus and the same bound CM");lengths[m]=sequences[m].length;
		}
		final SeedModelPositionTable table=new SeedModelPositionTable(seeds, names, lengths, clen, sites, hashes);
		for(int m=0; m<names.length; m++){
			for(int s=0; s<sites.length; s++){table.siteLow[m][s]=maps[m].boundaryLow(sites[s]);table.siteHigh[m][s]=maps[m].boundaryHigh(sites[s]);}
			final int[] inverseLow=new int[lengths[m]+1], inverseHigh=new int[lengths[m]+1];Arrays.fill(inverseLow, -1);
			for(int b=0; b<=clen; b++){
				for(int q=maps[m].boundaryLow(b); q<=maps[m].boundaryHigh(b); q++){
					if(inverseLow[q]<0){inverseLow[q]=b;}inverseHigh[q]=b;
				}
			}
			for(int q=0; q<inverseLow.length; q++){require(inverseLow[q]>=0, "The complete global map must represent every consensus boundary");}
			long word=0;int valid=0;
			for(int q=0; q<sequences[m].length; q++){
				final int base=RrnaResourceIO.base(sequences[m][q]);if(base<0){word=0;valid=0;continue;}
				word=((word<<2)|base)&MASK;valid++;if(valid<K){continue;}final Entry e=table.index.get(word);if(e==null){continue;}
				final int start=q-K+1, end=q+1;
				e.add(new Position(m, start, inverseLow[start], inverseHigh[start], inverseLow[end], inverseHigh[end]), lengths, clen);
			}
		}
		return table;
	}
	private SeedModelPositionTable(LongHashSet seeds, String[] names_, int[] lengths_, int clen_, int[] sites_, String[] hashes_){
		require(seeds!=null && seeds.size()>0 && names_!=null && lengths_!=null && names_.length>0
			&& names_.length==lengths_.length && clen_>0 && sites_!=null && sites_.length>0 && hashes_!=null && hashes_.length==4,
			"Position resource binds seed FASTA, consensus FASTA, Stockholm anchors and CM plus declared sites");
		for(String hash:hashes_){require(SeedOffsetTable.validSha80(hash), "Position resource provenance uses canonical sha80");}
		names=names_.clone();lengths=lengths_.clone();clen=clen_;sites=sites_.clone();hashes=hashes_.clone();
		for(int m=0; m<names.length; m++){
			require(names[m]!=null && names[m].matches("[A-Za-z0-9_.:-]+") && lengths[m]>0, "Resource model IDs and lengths must be explicit and safe TSV fields");
			for(int j=0; j<m; j++){require(!names[m].equals(names[j]), "Duplicate consensus names would merge distinct model coordinates");}
		}
		for(int s=0; s<sites.length; s++){require(sites[s]>0 && sites[s]<clen && (s==0 || sites[s]>sites[s-1]), "Insertion sites are unique ordered internal CM boundaries");}
		keys=seeds.toArray();Arrays.sort(keys);index=new LongObjectMap<Entry>(Math.max(16, keys.length*2), Entry.class);
		for(long key:keys){require(key>=0 && key<=MASK, "Packed forward seeds must contain exactly the K17 key space");index.put(key, new Entry());}
		siteLow=new int[names.length][sites.length];siteHigh=new int[names.length][sites.length];
	}
	/** Reads with actual caller resources supplied by the orchestrator. Exact
	 * metadata equality and full key conservation prevent loading another library. */
	public static SeedModelPositionTable load(String path, LongHashSet seeds, String[] names, int[] lengths, int clen, int[] sites, String[] hashes){
		final ByteFile in=RrnaResourceIO.open(path);
		try{return read(in,path,seeds,names,lengths,clen,sites,hashes);}
		finally{require(!in.close(), "Seed-position resource read failed");}
	}
	/** Shared parser also reads the unchanged table following an exact-consensus resource envelope. */
	static SeedModelPositionTable read(ByteFile in,String path,LongHashSet seeds,String[] names,int[] lengths,int clen,int[] sites,String[] hashes){
		final SeedModelPositionTable t=new SeedModelPositionTable(seeds,names,lengths,clen,sites,hashes);
		final LineParser1 p=new LineParser1('\t');
		{
			RrnaResourceIO.header(in, "#SeedModelPositionTable\t1\tK17\tEXACT_CONSENSUS", path);
			RrnaResourceIO.header(in, "#Counts\t"+t.keys.length+'\t'+names.length+'\t'+clen+'\t'+sites.length, path);
			for(int i=0; i<4; i++){RrnaResourceIO.header(in, "#Binding\t"+BINDINGS[i]+'\t'+hashes[i], path);}
			for(int m=0; m<names.length; m++){
				RrnaResourceIO.header(in, "#Model\t"+m+'\t'+names[m]+'\t'+lengths[m], path);
				for(int s=0; s<sites.length; s++){
					final byte[] line=in.nextLine();require(line!=null, "Every site needs a model-specific projection");p.set(line);
					require(p.terms()==5 && p.termEquals("#Site", 0) && p.parseInt(1)==m && p.parseInt(2)==sites[s], "Projected site identity must match the requested model and RF boundary");
					t.siteLow[m][s]=p.parseInt(3);t.siteHigh[m][s]=p.parseInt(4);
					require(t.siteLow[m][s]>=0 && t.siteHigh[m][s]>=t.siteLow[m][s] && t.siteHigh[m][s]<=lengths[m], "Site projection must fit its declared consensus");
				}
			}
			RrnaResourceIO.header(in, HEADER, path);int keyIndex=-1;long prior=-1;boolean unmapped=false;
			for(byte[] line; (line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==8, "Seed position rows preserve all eight fields");final long key=p.parseLongA48(0);
				if(key!=prior){keyIndex++;require(keyIndex<t.keys.length && key==t.keys[keyIndex], "Every original seed must appear once in sorted groups, including unmapped keys");prior=key;unmapped=false;}
				final Entry e=t.index.get(key);final int model=p.parseInt(1);
				if(model==-1){
					require(e.size()==0 && !unmapped, "An unmapped key has one sentinel row and no position rows");
					for(int c=2; c<8; c++){require(p.parseInt(c)==-1, "Unmapped keys cannot fabricate coordinate values");}unmapped=true;
				}else{
					require(!unmapped && p.parseInt(3)==p.parseInt(2)+K-1, "Mapped K17 seeds have one actual inclusive17-base consensus span");
					e.add(new Position(model, p.parseInt(2), p.parseInt(4), p.parseInt(5), p.parseInt(6), p.parseInt(7)), lengths, clen);
				}
			}
			require(keyIndex+1==t.keys.length, "Trailing unmapped keys cannot disappear from a position resource");
		}return t;
	}
	void write(String path){
		RrnaResourceIO.fresh(path, "");final ByteStreamWriter out=RrnaResourceIO.writer(path);
		try{write(out);}finally{RrnaResourceIO.finish(out);}
	}
	void write(ByteStreamWriter out){
		require(out!=null,"Position serialization requires a live writer owned by its caller");final ByteBuilder b=new ByteBuilder();
		{
			out.println("#SeedModelPositionTable\t1\tK17\tEXACT_CONSENSUS");out.println("#Counts\t"+keys.length+'\t'+names.length+'\t'+clen+'\t'+sites.length);
			for(int i=0; i<4; i++){out.println("#Binding\t"+BINDINGS[i]+'\t'+hashes[i]);}
			for(int m=0; m<names.length; m++){
				out.println("#Model\t"+m+'\t'+names[m]+'\t'+lengths[m]);
				for(int s=0; s<sites.length; s++){out.println("#Site\t"+m+'\t'+sites[s]+'\t'+siteLow[m][s]+'\t'+siteHigh[m][s]);}
			}
			out.println(HEADER);
			for(long key:keys){final Entry e=index.get(key);
				if(e.size()==0){out.print(b.clear().appendA48(key).append("\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n"));}
				for(int i=0; i<e.size(); i++){final Position p=e.get(i);
					out.print(b.clear().appendA48(key).tab().append(p.model).tab().append(p.start).tab().append(p.start+K-1)
						.tab().append(p.rfStartLow).tab().append(p.rfStartHigh).tab().append(p.rfEndLow).tab().append(p.rfEndHigh).nl());
				}
			}
		}
	}
	public Entry get(long key){return index.get(key);}
	public int seedCount(){return keys.length;}
	public int modelCount(){return names.length;}
	public String modelName(int model){return names[model];}
	public int modelLength(int model){return lengths[model];}
	public int siteCount(){return sites.length;}
	public int siteBoundary(int site){return sites[site];}
	/** -1=left,1=right,0=crosses or lies in the uncertain boundary interval. */
	public int side(Position p, int site){
		require(p!=null && site>=0 && site<sites.length, "Site attribution requires an actual stored model placement");
		return p.start+K<=siteLow[p.model][site] ? -1 : p.start>=siteHigh[p.model][site] ? 1 : 0;
	}
	public static final class Position {
		Position(int m, int s, int a, int b, int c, int d){model=m;start=s;rfStartLow=a;rfStartHigh=b;rfEndLow=c;rfEndHigh=d;}
		public final int model, start, rfStartLow, rfStartHigh, rfEndLow, rfEndHigh;
	}
	public static final class Entry {
		public int size(){return positions.size();}
		public Position get(int i){return positions.get(i);}
		/** Multiple models are conditional mappings, not repeated positions within one model. */
		public boolean ambiguous(){for(int i=1; i<size(); i++){if(get(i).model==get(i-1).model){return true;}}return false;}
		void add(Position p, int[] lengths, int clen){
			require(p.model>=0 && p.model<lengths.length && p.start>=0 && (long)p.start+K<=lengths[p.model]
				&& p.rfStartLow>=0 && p.rfStartHigh>=p.rfStartLow && p.rfEndLow>=p.rfStartLow
				&& p.rfEndHigh>=p.rfStartHigh && p.rfEndHigh>=p.rfEndLow && p.rfEndHigh<=clen,
				"Seed positions must retain ordered consensus and CM-boundary intervals");
			if(size()>0){final Position old=get(size()-1);require(p.model>old.model || p.model==old.model && p.start>old.start, "Postings are model/start sorted and may never duplicate one seed occurrence");}
			positions.add(p);
		}
		private final ArrayList<Position> positions=new ArrayList<Position>();
	}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final int K=17;
	static final long MASK=(1L<<(2*K))-1;
	static final String[] BINDINGS={"seedFastaSha80","consensusSha80","anchorsSha80","cmSha80"};
	static final String HEADER="seedA48\tmodel\tconsensusStart0\tconsensusStop0\tcmStartBoundaryLow\tcmStartBoundaryHigh\tcmEndBoundaryLow\tcmEndBoundaryHigh";
	final long[] keys;
	final String[] names, hashes;
	final int[] lengths, sites;
	final int clen;
	private final int[][] siteLow, siteHigh;
	private final LongObjectMap<Entry> index;
}
