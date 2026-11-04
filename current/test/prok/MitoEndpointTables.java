package test.prok;

import prok.NcrnaBoundaryScorer;
import prok.NcrnaBoundaryVectorGen;
import prok.NcrnaTextIO;
import prok.RrnaResourceIO;
import prok.SeedOffsetTable;
import prok.TrnaBoundaryFeatures;

import java.util.Arrays;
import java.util.function.Function;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;

/** One sparse resource for model-specific and pooled k9 endpoint counts.
 * Uses the existing V3 profile definition; no new scoring formula or aligner.
 * Construction membership/split decisions belong to the caller of add().
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointTables {
	public MitoEndpointTables(int models, int startInside, int stopInside, float meanLength){
		require(models>0 && models<1024 && startInside>=0 && startInside<=K && stopInside>=0 && stopInside<=K,
			"Packed count keys require bounded model slots and explicit k9 placement");
		require(Float.isFinite(meanLength) && meanLength>0,"A single family-wide length normalization must survive table serialization and scoring");
		this.models=models;this.meanLength=meanLength;inside=new int[]{startInside,stopInside};views=new TrnaBoundaryFeatures.NinemerTable[models+1][2];
		for(int model=-1; model<models; model++){for(int end=0; end<2; end++){
			final int slot=slot(model,end);views[model+1][end]=(seq,pos)->{
				final int key=pack(seq,pos);return key<0?0:Math.max(0,counts.get(slot|key));
			};
		}}
	}
	/** Count one original member once, in its assigned model and the pooled table. */
	public void add(byte[] bases, int start, int stop, int model){
		require(model>=0 && model<models && start>=0 && stop>=start && stop<bases.length,
			"Training observations require a known model and true genomic endpoints");
		for(int end=0; end<2; end++){
			final int pos=TrnaBoundaryFeatures.kmerWindowStart(end==0?start:stop,type(end),inside[end],K-inside[end]);
			final int key=pack(bases,pos);if(key<0){continue;}
			counts.increment(slot(model,end)|key);counts.increment(slot(-1,end)|key);
		}
	}
	public TrnaBoundaryFeatures.NinemerTable table(int model, int end){
		require(model>=-1 && model<models && end>=0 && end<2,"Endpoint lookup must bind a valid model or explicit family pool");return views[model+1][end];
	}
	public int insideCount(int end){require(end>=0 && end<2,"Endpoint placement requires start or stop");return inside[end];}
	/** Match the scorer's independently normalized50/50 family/model profiles. */
	public float[] features(byte[] bases, int start, int stop, int model, int end, float gc, float identity, int trueEnd){
		require(model>=0 && model<models && start>=0 && stop>=start && stop<bases.length && end>=0 && end<2,
			"V3 vectors must retain the caller's candidate span, selected model and length definition");
		final int pos=end==0?start:stop;
		final float[] local=profile(bases,pos,model,end,trueEnd),family=profile(bases,pos,-1,end,trueEnd);
		NcrnaBoundaryVectorGen.blendProfiles(local,family);
		return NcrnaBoundaryScorer.buildV3Input(local,end==1,(stop-start+1)/meanLength,gc,identity);
	}
	/** Development harness, not installed mitochondrial caller dispatch. Both ends
	 * use the real confidence-ordered sequential scorer and the stored family mean. */
	public NcrnaBoundaryScorer.RefinementResult refine(CellNet net,byte[] bases,int start,int stop,int model,
			byte[] consensus,float gc,float identity,int[] startOffsets,int[] stopOffsets,int windowEdgeMask){
		require(model>=0 && model<models && consensus!=null && consensus.length>0,
			"Development scoring must bind a real model and both pooled/model tables");
		require((windowEdgeMask&~3)==0,"Physical-window edge provenance uses only the two sense-end bits");
		return NcrnaBoundaryScorer.refineBoundariesWithCounts(net,net,bases,start,stop,consensus,null,
			table(model,0),table(model,1),inside[0],K-inside[0],inside[1],K-inside[1],gc,meanLength,
			startOffsets,stopOffsets,0,0,(windowEdgeMask&1)!=0,(windowEdgeMask&2)!=0,NcrnaBoundaryScorer.FEATURES_V3,
			identity,null,null,table(-1,0),table(-1,1));
	}
	/** Recover physical-edge provenance from sealed source coordinates, never from
	 * a predicted anchor. Explicit origin-crossing windows have no physical ends. */
	public static int windowEdgeMask(MitoEndpointGeometry source,int genomeLength,int truthStart,int truthStop,int windowLength){
		require(source!=null && !source.mixedStrand && genomeLength>0 && truthStart>=0 && truthStop>=truthStart && truthStop<windowLength,
			"Physical edges require a supported source window with valid true coordinates");
		final int last=source.lo.length-1;
		final int first=source.reverse[0]?genomeLength-source.hi[0]:source.lo[0]-1;
		int end=source.reverse[last]?genomeLength-source.lo[last]:source.hi[last]-1;
		if(source.crossesOrigin){end+=genomeLength;}
		require(end-first==truthStop-truthStart,"True genomic span must reproduce the sealed source-segment geometry");
		if(source.crossesOrigin){require(windowLength<=genomeLength,"Circular windows cannot repeat a genome");return 0;}
		final int left=first-truthStart,right=end+windowLength-1-truthStop;
		require(left>=0 && right<genomeLength,"Sealed window must fit the original physical contig");
		return (left==0?1:0)|(right==genomeLength-1?2:0);
	}
	/** One fixed median for the added library; never a selected-model denominator. */
	public static float medianLength(byte[][] library,int first){
		require(library!=null && first>=0 && first<library.length,"Length normalization requires the explicit mitochondrial model subset");
		final int[] lengths=new int[library.length-first];
		for(int i=first;i<library.length;i++){require(library[i]!=null && library[i].length>0,"Empty consensus cannot define endpoint length normalization");lengths[i-first]=library[i].length;}
		Arrays.sort(lengths);final int mid=lengths.length/2;
		return (lengths.length&1)==1?lengths[mid]:((long)lengths[mid-1]+lengths[mid])/2f;
	}
	private float[] profile(byte[] bases,int pos,int model,int end,int trueEnd){
		assert(end>=0 && end<2):"Inside/outside placement is endpoint-specific";
		return trueEnd<0 ? TrnaBoundaryFeatures.enrichmentProfile(bases,pos,type(end),inside[end],K-inside[end],table(model,end))
			:TrnaBoundaryFeatures.enrichmentProfile(bases,pos,trueEnd,type(end),inside[end],K-inside[end],table(model,end));
	}
	public void write(String path,String consensusSha80){
		write(path, consensusSha80, pathName -> new ByteStreamWriter(fileIO.FileFormat.testOutput(pathName,
			fileIO.FileFormat.TEXT, null, true, false, false, false)));
	}
	/** Own startup, scratch and serialization; cleanup must not replace a write failure. */
	private void write(String path, String consensusSha80, Function<String, ByteStreamWriter> factory){
		require(SeedOffsetTable.validSha80(consensusSha80),"The packed tables must bind the exact consensus library");validatePool();RrnaResourceIO.fresh(path,"");
		final IntList keys=new IntList();final int[] all=counts.keys();for(int key:all){if(key!=counts.invalid()){keys.add(key);}}
		Arrays.sort(keys.array,0,keys.size);
		try(NcrnaTextIO.Output output=new NcrnaTextIO.Output(factory.apply(path))){
			final ByteStreamWriter out=output.writer;final ByteBuilder b=new ByteBuilder();
			out.println(HEADER);out.print(b.clear().append(models).tab().append(K).tab().append(inside[0]).tab().append(inside[1]).tab().append(consensusSha80).tab().appendSlow(meanLength).nl());
			out.println("model_index\tend\tpacked_kmer\tcount");
			for(int i=0;i<keys.size;i++){final int key=keys.array[i],s=key>>>BITS;
				out.print(b.clear().append((s>>>1)-1).tab().append(s&1).tab().append(key&MASK).tab().append(counts.get(key)).nl());}
		}
	}
	public static MitoEndpointTables load(String path,int models,String consensusSha80){
		require(SeedOffsetTable.validSha80(consensusSha80),"Inference must supply an exact library pin");final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,HEADER,path);final byte[] row=in.nextLine();require(row!=null,"Missing packed table metadata");p.set(row);
			require(p.terms()==6 && p.parseInt(0)==models && p.parseInt(1)==K && p.termEquals(consensusSha80,4),"Table count/k/library/mean must match the explicit development feature contract");
			final MitoEndpointTables t=new MitoEndpointTables(models,p.parseInt(2),p.parseInt(3),p.parseFloat(5));RrnaResourceIO.header(in,"model_index\tend\tpacked_kmer\tcount",path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==4,"Packed table schema changed");final int m=p.parseInt(0),end=p.parseInt(1),key=p.parseInt(2),count=p.parseInt(3);
				require(m>=-1 && m<models && end>=0 && end<2 && key>=0 && key<=MASK && count>0,"Invalid packed table observation");
				final int index=t.slot(m,end)|key;require(t.counts.get(index)<0,"Duplicate packed observation");t.counts.put(index,count);}
			t.validatePool();return t;
		}finally{require(!in.close(),"Packed table read failed");}
	}
	private int slot(int model,int end){
		assert(model>=-1 && model<models && end>=0 && end<2):"Family pool and per-model endpoints require disjoint key slots";return ((model+1)*2+end)<<BITS;
	}
	/** The family pool must be exactly the sum of model observations, per end/key. */
	private void validatePool(){
		final IntHashMap sum=new IntHashMap();final int[] keys=counts.keys(),values=counts.values();
		for(int i=0;i<keys.length;i++){
			if(keys[i]==counts.invalid()){continue;}
			final int s=keys[i]>>>BITS;if(s<2){continue;}
			final int key=((s&1)<<BITS)|(keys[i]&MASK),before=Math.max(0,sum.get(key));
			require(values[i]>0 && before<=Integer.MAX_VALUE-values[i],"Endpoint observation totals must not overflow");sum.put(key,before+values[i]);
		}
		for(int i=0;i<keys.length;i++){
			if(keys[i]!=counts.invalid() && (keys[i]>>>BITS)<2){require(sum.get(keys[i])==values[i],"Stored family count differs from its model-member sum");}
		}
		for(int key:sum.keys()){if(key!=sum.invalid()){require(counts.get(key)==sum.get(key),"Every model observation requires an equal pooled contribution");}}
	}
	private static int pack(byte[] bases,int pos){
		if(pos<0 || pos>bases.length-K){return -1;}int key=0;
		for(int i=0;i<K;i++){final int v=RrnaResourceIO.base(bases[pos+i]);if(v<0){return -1;}key=(key<<2)|v;}
		assert(key>=0 && key<=MASK):"Forward k9 codes must fit the reserved18bits";return key;
	}
	private static TrnaBoundaryFeatures.BoundaryType type(int end){return end==0?TrnaBoundaryFeatures.BoundaryType.START:TrnaBoundaryFeatures.BoundaryType.STOP;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final int models;public final float meanLength;private final int[] inside;
	private final IntHashMap counts=new IntHashMap();
	private final TrnaBoundaryFeatures.NinemerTable[][] views;
	public static final int K=9;
	private static final int BITS=18,MASK=(1<<BITS)-1;
	private static final String HEADER="models\tk\tstart_inside\tstop_inside\tconsensus_sha80\tmean_length";
}
