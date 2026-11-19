package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import prot.FamilyTailCascade.Comparison;
import prot.FamilyTailCascade.Counters;
import prot.FamilyTailCascade.Result;
import prot.FamilyTailCascade.SearchConfig;
import prot.FamilyTailCascade.Searcher;
import prot.FamilyTailCascade.Truth;

/** Experimental finite-grid measurement, never a replacement for production timing.
 * One DDL score per query; production candidate preparation per reject; one alignment
 * per distinct pair. First-winner B intervals update A/B difference arrays.
 * @author Sayu
 */
final class FamilyThresholdGrid{

	private FamilyThresholdGrid(){}
	static final int[] R={4,5,6,7,8,9,10};
	static final int[] A={12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,1025};
	static final float[] I={30,31,32,33,34,35,36,37,38,39,40,50,100};
	static final float[] B=blosumFloors(new double[]{0,.118,.125,.130,.135,.140,.145,.150,.160,.180,.200,0x1.fffffep127});
	static final GridDefinition DEFAULT_GRID=new GridDefinition("default",R,A,I,B);
	static final GridDefinition BOUNDED_GRID=new GridDefinition("bounded",new int[]{6},new int[]{1025},boundedIdentities(),
		blosumFloors(new double[]{.225,.250,.300,.350,.400,.500,0x1.fffffep127}));
	static final int TUPLES=DEFAULT_GRID.tuples;
	static final int WITNESS_LIMIT=40;
	static final int ASSIGNED=0,CORRECT=1,DDL_ACCEPT=2,IDENTITY=3,BLOSUM=4,
		SCALAR_REJECT=5,IDENTITY_REJECT=6,BLOSUM_REJECT=7,DP=8,ID_REJECTS=9,B_REJECTS=10,
		RANK_FAMILIES=11,RANK_COSINES=12,RANK_ZERO=13,DDL_ACTIVE=14,WINDOW=15,
		SHORTLIST=16,SCALAR_CHECKED=17,SCALAR_REMOVED=18;
	static final String[] METRICS={"assigned","correct","ddl_accept","identity_accept","blosum_accept",
		"scalar_reject","identity_reject_terminal","blosum_reject_terminal","dp_attempts","identity_rejects","blosum_rejects",
		"rank_families","rank_cosines","rank_zero_norm","ddl_active","cosine_window","dp_shortlist","scalar_checked","scalar_removed"};
	private static float[] boundedIdentities(){
		final float[] result=new float[50];
		for(int i=0;i<49;i++){result[i]=51+i;}
		result[49]=100;
		return result;
	}
	private static float[] blosumFloors(final double[] values){
		final float[] result=new float[values.length];
		for(int j=0;j<result.length;j++){result[j]=FamilyGateCodec.roundUpper(values[j]);}
		assert(result[result.length-1]==Float.MAX_VALUE) : "B-off must remain finite and exact through the CLI's directed rounding";
		return result;
	}
	static GridDefinition gridDefinition(String name){
		if(name==null || name.equals("default")){return DEFAULT_GRID;}
		if(name.equals("bounded")){return BOUNDED_GRID;}
		throw new IllegalArgumentException("gridaxes must be default or bounded");
	}

	/** For each I, unaccepted B settings are a suffix. A candidate consumes a prefix
	 * of that suffix (or all of it via identity). Emit only first-winner intervals. */
	static void evaluateFloors(GridDefinition grid,int count,int[] families,float[] identity,float[] blosum,
			float[] lower,float[] upper,float[] finals,FloorSink sink){
		assert(count>0 && count<=families.length) : "Floor evaluation requires the already capped, scalar-filtered DP list";
		for(int i=0;i<grid.I.length;i++){
			int start=0,ir=0,br=0,last=BLOSUM_REJECT;
			for(int c=0;c<count && start<grid.B.length;c++){
				if(identity[c]>Math.max(upper[c],grid.I[i])){
					sink.accept(i,start,grid.B.length,families[c],IDENTITY,c+1,ir,br);
					start=grid.B.length;
				}else if(identity[c]<lower[c]){ir++; last=IDENTITY_REJECT;}
				else{
					int end=start;
					if(blosum[c]>finals[c]){
						while(end<grid.B.length && blosum[c]>grid.B[end]){end++;}
					}
					if(end>start){sink.accept(i,start,end,families[c],BLOSUM,c+1,ir,br); start=end;}
					// Remaining B suffix rejected this candidate; accepted prefixes have left.
					br++; last=BLOSUM_REJECT;
				}
			}
			if(start<grid.B.length){sink.accept(i,start,grid.B.length,-1,last,count,ir,br);}
		}
	}
	interface FloorSink{
		void accept(int identityIndex,int bStart,int bEnd,int family,int stage,int dpAttempts,int identityRejects,int blosumRejects);
	}

	static Summary run(SearchConfig config,Path inputPath,Path output,Path truthPath,
			HashMap<String,String> options,int threads,long expected) throws Exception{
		assert(config.truth==null && !config.union) : "Grid correctness uses a separate truth map, never ranking diagnostics";
		final GridDefinition grid=gridDefinition(options.get("gridaxes"));
		final String cohort=options.get("gridcohort");
		final boolean negative="negative".equals(cohort);
		if(!negative && !"positive".equals(cohort)){throw new IllegalArgumentException("gridcohort must be negative or positive");}
		final int[] selected=selectedWitnessTuples(options.get("gridwitnessselected"),negative,grid);
		if(negative== (truthPath!=null) || !negative && !options.containsKey("truthsha80")){
			throw new IllegalArgumentException("Positive grid requires pinned truth; negative grid forbids truth");
		}
		final long offset=Long.parseLong(options.getOrDefault("gridoffset","0"));
		if(offset<0 || expected<1 || expected>1000000 || offset>10000000-expected){
			throw new IllegalArgumentException("Grid shards require 1..1M reads and ordinals within the10M cohort");
		}
		final Path witnesses=options.containsKey("gridwitnessout") ? Paths.get(options.get("gridwitnessout")).toAbsolutePath().normalize() : null;
		if(negative && witnesses==null || !negative && witnesses!=null){throw new IllegalArgumentException("Only negative grids require gridwitnessout");}
		if(witnesses!=null && (Files.exists(witnesses) || Files.isSymbolicLink(witnesses) || witnesses.toFile().getCanonicalFile().equals(output.toFile().getCanonicalFile()))){
			throw new IllegalArgumentException("Witness output exists or aliases grid output");
		}
		final HashMap<String,Truth> truth=negative ? null : FamilyTailCascade.loadTruth(truthPath,config.families);
		final FileFormat format=FileFormat.testInput(inputPath.toString(),FileFormat.FASTA,null,false,true);
		final Streamer stream=StreamerFactory.makeStreamer(format,null,true,expected,false,true,1);
		final Source source=new Source(stream,offset);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Worker[] workers=new Worker[threads];
		final long start=System.nanoTime();
		try{
			for(int w=0;w<threads;w++){workers[w]=new Worker(source,config,truth,negative,selected,grid,failure,w);}
			stream.start();
			for(Worker w:workers){w.start();}
			for(Worker w:workers){w.join();}
		}catch(Throwable t){
			t.printStackTrace();KillSwitch.kill("Grid startup/join failure cannot leave query producers blocked");
			throw new IllegalStateException("Grid startup/join failed",t);
		}finally{if(ReadWrite.closeStream(stream)){throw new IllegalStateException("Grid query read failed");}}
		if(failure.get()!=null){throw new IllegalStateException("Grid worker failed",failure.get());}
		final Counts total=new Counts(grid,negative,selected);
		long alignments=0,reads=0;
		for(Worker w:workers){total.merge(w.engine.counts); alignments+=w.engine.alignments; reads+=w.engine.reads;}
		if(reads!=expected || source.next!=offset+expected){throw new IllegalStateException("Grid query count/ordinal coverage mismatch");}
		total.finish();
		total.write(output,witnesses,cohort,offset,reads);
		System.out.println("FAMILY_THRESHOLD_GRID_TIMING query_aggregate_write_seconds="+((System.nanoTime()-start)*1e-9));
		return new Summary(reads,alignments,offset,grid.tuples);
	}
	/** Complete evidence is retained only for an explicitly bounded set of frozen tuples. */
	private static int[] selectedWitnessTuples(String text,boolean negative,GridDefinition grid){
		if(text==null){return new int[0];}
		if(!negative){throw new IllegalArgumentException("gridwitnessselected requires a negative grid");}
		final String[] fields=text.split(",",-1);
		if(fields.length<1 || fields.length>5){throw new IllegalArgumentException("Select one to five witness tuples, never the whole grid");}
		final int[] selected=new int[fields.length];
		for(int i=0;i<fields.length;i++){
			selected[i]=Integer.parseInt(fields[i]);
			if(selected[i]<0 || selected[i]>=grid.tuples){throw new IllegalArgumentException("Witness tuple outside measured grid: "+selected[i]);}
		}
		Arrays.sort(selected);
		for(int i=1;i<selected.length;i++){
			if(selected[i]==selected[i-1]){throw new IllegalArgumentException("Duplicate selected witness tuple: "+selected[i]);}
		}
		return selected;
	}
	static final class Summary{
		Summary(long reads_,long alignments_,long offset_,int tuples_){reads=reads_;alignments=alignments_;offset=offset_;tuples=tuples_;}
		final long reads,alignments,offset;
		final int tuples;
	}
	private static final class Source{
		Source(Streamer stream_,long offset){stream=stream_;next=offset;}
		synchronized Batch next(){
			final ListNum<Read> list=stream.nextList();
			if(list==null){return null;}
			final Batch result=new Batch(list,next); next+=list.list.size(); return result;
		}
		final Streamer stream;
		long next;
	}
	private static final class Batch{
		Batch(ListNum<Read> list_,long start_){list=list_;start=start_;}
		final ListNum<Read> list;
		final long start;
	}
	private static final class Worker extends Thread{
		Worker(Source source_,SearchConfig config,HashMap<String,Truth> truth_,boolean negative,int[] selected,GridDefinition grid,
				AtomicReference<Throwable> failure_,int id){
			super("FamilyThresholdGrid-"+id);source=source_;truth=truth_;failure=failure_;engine=new Engine(config,negative,selected,grid);
		}
		@Override public void run(){
			try{
				for(Batch batch=source.next();batch!=null;batch=source.next()){
					if(failure.get()!=null){continue;}
					try{
						long ordinal=batch.start;
						for(Read read:batch.list.list){
							if(read==null || read.bases==null || read.mate!=null){throw new IllegalArgumentException("Grid requires unpaired proteins");}
							final ProteinSequence query=new ProteinSequence(FamilyTailCascade.firstToken(read.id),read.bases);
							final Truth t=truth==null ? null : truth.get(query.id);
							if(truth!=null && t==null){throw new IllegalArgumentException("Missing positive truth: "+query.id);}
							engine.observe(Arrays.copyOf(read.bases,query.length()),query,ordinal++,t==null ? -1 : t.ordinal);
						}
					}catch(Throwable t){failure.compareAndSet(null,t);}
				}
			}catch(Throwable t){t.printStackTrace();KillSwitch.kill("Unrecoverable grid input dispatch failure");}
		}
		final Source source;
		final HashMap<String,Truth> truth;
		final AtomicReference<Throwable> failure;
		final Engine engine;
	}

	/** Worker-local state; score arrays survive only until the next query generation. */
	static final class Engine{
		Engine(SearchConfig config,boolean negative){
			this(config,negative,new int[0],DEFAULT_GRID);
		}
		Engine(SearchConfig config,boolean negative,int[] selected,GridDefinition grid_){
			if(config.ddlAccept!=1025 || config.truth!=null || config.top!=50 || !config.scalarGates ||
					config.cosineMode!=FamilyTailCascade.CosineMode.RANKING || config.features!=FamilyTailCascade.RankingFeatures.COMBINED25 || config.union){
				throw new IllegalArgumentException("Grid engine requires fixed top50, combined25/scalars and DDL-off preparation");
			}
			grid=grid_; searcher=new Searcher(config); counts=new Counts(grid,negative,selected);
			seen=new long[config.families.length]; cachedIdentity=new float[seen.length];cachedBlosum=new float[seen.length];
		}
		void observe(byte[] raw,ProteinSequence query,long ordinal,int truth){
			assert(query!=null && raw!=null && ordinal>=0) : "Input ordinals identify exact selected-negative witnesses";
			reads++; generation++;
			if(generation==0){Arrays.fill(seen,0);generation=1;}
			final Witness witness=counts.negative ? new Witness(ordinal,query.id) : null;
			int bestScore=-1,bestFamily=-1,self=0; boolean selfReady=false;
			for(int r=0;r<grid.R.length;r++){
				if(r>0 && bestScore<grid.R[r]){counts.ddlRejected[r]++;continue;}
				searcher.ddlReject=grid.R[r];
				final long oldRank=prep.rankFamilies,oldCos=prep.rankCosines,oldZero=prep.rankZeroNorm,
					oldChecked=prep.scalarPairsChecked,oldRemoved=prep.scalarPairsRemoved;
				final boolean ready=searcher.prepare(raw,query,result,prep,null,r>0);
				if(r==0){
					for(int j=0;j<searcher.work.size;j++){
						final int f=searcher.work.touched[j],score=searcher.work.scores[f];
						if(score>bestScore || score==bestScore && (bestFamily<0 || searcher.families[f].rank<searcher.families[bestFamily].rank)){
							bestScore=score;bestFamily=f;
						}
					}
				}
				if(bestScore<grid.R[r]){counts.ddlRejected[r]++;continue;}
				int firstA=0;while(firstA<grid.A.length && grid.A[firstA]<=bestScore){firstA++;}
				for(int i=0;i<grid.I.length;i++){
					counts.range(DDL_ACTIVE,r,i,0,grid.A.length,0,grid.B.length,result.ddlActive);
					if(firstA>0){counts.outcome(r,i,0,firstA,0,grid.B.length,bestFamily,DDL_ACCEPT,0,0,0,truth,witness);}
					counts.range(RANK_FAMILIES,r,i,firstA,grid.A.length,0,grid.B.length,prep.rankFamilies-oldRank);
					counts.range(RANK_COSINES,r,i,firstA,grid.A.length,0,grid.B.length,prep.rankCosines-oldCos);
					counts.range(RANK_ZERO,r,i,firstA,grid.A.length,0,grid.B.length,prep.rankZeroNorm-oldZero);
					counts.range(WINDOW,r,i,firstA,grid.A.length,0,grid.B.length,result.cosineWindow);
					counts.range(SHORTLIST,r,i,firstA,grid.A.length,0,grid.B.length,result.dpShortlist);
					counts.range(SCALAR_CHECKED,r,i,firstA,grid.A.length,0,grid.B.length,prep.scalarPairsChecked-oldChecked);
					counts.range(SCALAR_REMOVED,r,i,firstA,grid.A.length,0,grid.B.length,prep.scalarPairsRemoved-oldRemoved);
				}
				if(!ready){
					if(!"scalar_reject".equals(result.terminal)){throw new IllegalStateException("Unexpected grid preparation terminal: "+result.terminal);}
					for(int i=0;i<grid.I.length;i++){counts.outcome(r,i,firstA,grid.A.length,0,grid.B.length,-1,SCALAR_REJECT,0,0,0,truth,witness);}
					continue;
				}
				if(!selfReady){self=ProteinSearcher.selfScore(query);selfReady=true;}
				for(int c=0;c<result.dpShortlist;c++){
					final Comparison candidate=searcher.shortlist.get(c);final int f=candidate.ordinal;
					if(seen[f]!=generation){
						FamilyCalibration.align(query,candidate.family.consensus,self,searcher.metrics);
						cachedIdentity[f]=(float)searcher.metrics.identityPct;cachedBlosum[f]=(float)searcher.metrics.blosum;
						seen[f]=generation;alignments++;
					}
					families[c]=f;identity[c]=cachedIdentity[f];blosum[c]=cachedBlosum[f];
					lower[c]=candidate.family.identityLower;upper[c]=candidate.family.identityUpper;finals[c]=candidate.family.blosumFinal;
				}
				final int rr=r,aa=firstA;
				evaluateFloors(grid,result.dpShortlist,families,identity,blosum,lower,upper,finals,
					(i,b0,b1,f,stage,dp,ir,br)->counts.outcome(rr,i,aa,grid.A.length,b0,b1,f,stage,dp,ir,br,truth,witness));
			}
		}
		final GridDefinition grid;
		final Searcher searcher;
		final Counts counts;
		final Result result=new Result();
		final Counters prep=new Counters();
		final long[] seen;
		final float[] cachedIdentity,cachedBlosum;
		final int[] families=new int[50];
		final float[] identity=new float[50],blosum=new float[50],lower=new float[50],upper=new float[50],finals=new float[50];
		long generation,reads,alignments;
	}

	static final class Witness implements Comparable<Witness>{
		Witness(long ordinal_,String id_){ordinal=ordinal_;id=id_;}
		@Override public int compareTo(Witness other){return Long.compare(ordinal,other.ordinal);}
		final long ordinal;
		final String id;
	}

	/** Difference arrays make A/B rectangle aggregation independent of rectangle area. */
	static final class Counts{
		Counts(boolean negative_){this(DEFAULT_GRID,negative_,new int[0]);}
		Counts(GridDefinition grid_,boolean negative_,int[] selected){
			assert(selected!=null && selected.length<=5 && (negative_ || selected.length==0)) : "Complete witnesses are bounded to five selected negative tuples";
			grid=grid_; negative=negative_;witnesses=negative ? new Witness[grid.tuples][] : null; sizes=negative ? new byte[grid.tuples] : null;
			complete=new ArrayList<ArrayList<Witness>>(selected.length);
			selectedTuples=selected.clone();
			selectedSlots=selected.length==0 ? null : new int[grid.tuples];
			if(selectedSlots!=null){Arrays.fill(selectedSlots,-1);}
			delta=new long[METRICS.length][grid.cells];
			ddlRejected=new long[grid.R.length];
			for(int t:selectedTuples){
				assert(t>=0 && t<grid.tuples) : "Selected witness tuple must belong to the finite measured grid";
				assert(selectedSlots[t]<0) : "A selected tuple has exactly one complete witness list";
				selectedSlots[t]=complete.size();
				complete.add(new ArrayList<Witness>());
			}
		}
		void range(int metric,int r,int i,int a0,int a1,int b0,int b1,long value){
			assert(!finished && a0>=0 && a1<=grid.A.length && b0>=0 && b1<=grid.B.length) : "Difference updates precede prefix reconstruction and use half-open grid intervals";
			if(value==0 || a0==a1 || b0==b1){return;}
			final int base=(r*grid.I.length+i)*grid.plane,stride=grid.B.length+1;
			final long[] d=delta[metric];
			d[base+a0*stride+b0]+=value; d[base+a1*stride+b0]-=value;
			d[base+a0*stride+b1]-=value; d[base+a1*stride+b1]+=value;
		}
		void outcome(int r,int i,int a0,int a1,int b0,int b1,int family,int stage,int dp,int ir,int br,int truth,Witness witness){
			assert(stage>=DDL_ACCEPT && stage<=BLOSUM_REJECT) : "Each grid interval has one disjoint production terminal stage";
			range(stage,r,i,a0,a1,b0,b1,1);
			range(DP,r,i,a0,a1,b0,b1,dp);range(ID_REJECTS,r,i,a0,a1,b0,b1,ir);range(B_REJECTS,r,i,a0,a1,b0,b1,br);
			if(family>=0){
				range(ASSIGNED,r,i,a0,a1,b0,b1,1);
				if(family==truth){range(CORRECT,r,i,a0,a1,b0,b1,1);}
				if(negative){
					for(int a=a0;a<a1;a++){for(int b=b0;b<b1;b++){remember(grid.tuple(r,a,i,b),witness);}}
				}
			}
		}
		void remember(int t,Witness witness){
			assert(negative && witness!=null) : "Only observed negative input ordinals may become FP witnesses";
			final int selected=selectedSlots==null ? -1 : selectedSlots[t];
			if(selected>=0){complete.get(selected).add(witness);return;}
			final int size=sizes[t];
			if(size>WITNESS_LIMIT){return;}
			if(size==WITNESS_LIMIT){witnesses[t]=null;sizes[t]=WITNESS_LIMIT+1;return;}
			if(size==0){witnesses[t]=new Witness[WITNESS_LIMIT];}
			witnesses[t][size]=witness;sizes[t]++;
		}
		void merge(Counts other){
			assert(!finished && !other.finished && grid==other.grid && negative==other.negative) : "Merge raw grid differences before prefix summation";
			if(!Arrays.equals(selectedTuples,other.selectedTuples)){throw new IllegalArgumentException("Grid workers disagree about selected witness tuples");}
			for(int m=0;m<delta.length;m++){for(int c=0;c<grid.cells;c++){delta[m][c]+=other.delta[m][c];}}
			for(int r=0;r<grid.R.length;r++){ddlRejected[r]+=other.ddlRejected[r];}
			if(negative){
				for(int i=0;i<complete.size();i++){complete.get(i).addAll(other.complete.get(i));}
				for(int t=0;t<grid.tuples;t++){
					if(sizes[t]+other.sizes[t]>WITNESS_LIMIT){sizes[t]=WITNESS_LIMIT+1;witnesses[t]=null;}
					else{for(int j=0;j<other.sizes[t];j++){remember(t,other.witnesses[t][j]);}}
				}
			}
		}
		void finish(){
			assert(!finished) : "A grid difference array must be reconstructed exactly once";
			final int stride=grid.B.length+1;
			for(long[] d:delta){for(int r=0;r<grid.R.length;r++){for(int i=0;i<grid.I.length;i++){
				final int base=(r*grid.I.length+i)*grid.plane;
				for(int a=0;a<grid.A.length;a++){for(int b=0;b<grid.B.length;b++){
					final int c=base+a*stride+b;
					if(a>0){d[c]+=d[c-stride];}if(b>0){d[c]+=d[c-1];}if(a>0 && b>0){d[c]-=d[c-stride-1];}
					if(d[c]<0){throw new IllegalStateException("Negative reconstructed grid count");}
				}}
			}}}
			finished=true;
		}
		long value(int m,int r,int a,int i,int b){
			assert(finished) : "Grid counts are readable only after prefix reconstruction";
			return delta[m][(r*grid.I.length+i)*grid.plane+a*(grid.B.length+1)+b];
		}
		void write(Path output,Path witnessOut,String cohort,long offset,long reads){
			assert(finished) : "Only validated reconstructed aggregate counts may be published";
			final ByteStreamWriter writer=aggregateWriter(output);
			final ByteStreamWriter ww=witnessOut==null ? null : aggregateWriter(witnessOut);
			writer.start();if(ww!=null){ww.start();}
			try{
				final ByteBuilder header=new ByteBuilder("tuple\treject\taccept\tidentity_floor\tblosum_floor\tblosum_bits\tcohort\tstart\tcount\tddl_reject\tfalse_accept\twrong_family");
				for(String m:METRICS){header.tab().append(m);}writer.print(header.nl());
				if(ww!=null){ww.print("tuple\tordinal\tquery_id\n");}
				for(int r=0;r<grid.R.length;r++){for(int a=0;a<grid.A.length;a++){for(int i=0;i<grid.I.length;i++){for(int b=0;b<grid.B.length;b++){
					final int t=grid.tuple(r,a,i,b);
					final long assigned=value(ASSIGNED,r,a,i,b),correct=value(CORRECT,r,a,i,b);
					long terminals=ddlRejected[r];for(int m=DDL_ACCEPT;m<=BLOSUM_REJECT;m++){terminals+=value(m,r,a,i,b);}
					if(terminals!=reads || correct>assigned || assigned!=value(DDL_ACCEPT,r,a,i,b)+value(IDENTITY,r,a,i,b)+value(BLOSUM,r,a,i,b)){
						throw new IllegalStateException("Grid terminal/assignment counts fail conservation: tuple"+t);
					}
					final ByteBuilder row=new ByteBuilder().append(t).tab().append(grid.R[r]).tab().append(grid.A[a]).tab().append(Float.toString(grid.I[i])).tab()
						.append(Float.toHexString(grid.B[b])).tab().append(Integer.toHexString(Float.floatToRawIntBits(grid.B[b]))).tab().append(cohort).tab()
						.append(offset).tab().append(reads).tab().append(ddlRejected[r]).tab().append(negative ? assigned : 0).tab().append(negative ? 0 : assigned-correct);
					for(int m=0;m<METRICS.length;m++){row.tab().append(value(m,r,a,i,b));}writer.print(row.nl());
					final int selected=selectedSlots==null ? -1 : selectedSlots[t];
					if(negative && selected>=0){
						final ArrayList<Witness> list=complete.get(selected);
						if(list.size()!=assigned){throw new IllegalStateException("Missing complete selected witnesses: tuple"+t);}
						Collections.sort(list);
						long previous=-1;
						for(Witness w:list){
							if(w.ordinal<offset || w.ordinal>=offset+reads || w.ordinal<=previous){throw new IllegalStateException("Invalid or duplicate selected witness ordinal");}
							previous=w.ordinal;
							ww.print(new ByteBuilder().append(t).tab().append(w.ordinal).tab().append(w.id).nl());
						}
					}else if(negative && assigned<=WITNESS_LIMIT){
						if(sizes[t]!=assigned){throw new IllegalStateException("Missing bounded negative witnesses: tuple"+t);}
						if(assigned>0){
							Arrays.sort(witnesses[t],0,sizes[t]);
							for(int j=0;j<sizes[t];j++){
								final Witness w=witnesses[t][j];
								if(w.ordinal<offset || w.ordinal>=offset+reads || j>0 && witnesses[t][j-1].ordinal==w.ordinal){throw new IllegalStateException("Invalid negative witness ordinal");}
								ww.print(new ByteBuilder().append(t).tab().append(w.ordinal).tab().append(w.id).nl());
							}
						}
					}
				}}}}
			}finally{
				boolean bad=writer.poisonAndWait();if(ww!=null){bad=ww.poisonAndWait() || bad;}
				if(bad){throw new IllegalStateException("Grid aggregate/witness write failed");}
			}
		}
		/** A single producer publishes these tables serially through print(), not batch IDs. */
		private static ByteStreamWriter aggregateWriter(Path path){
			return new ByteStreamWriter(FileFormat.testOutput(path.toString(),FileFormat.TEXT,null,false,false,false,false));
		}
		final GridDefinition grid;
		final long[][] delta;
		final long[] ddlRejected;
		final Witness[][] witnesses;
		final byte[] sizes;
		final int[] selectedTuples,selectedSlots;
		final ArrayList<ArrayList<Witness>> complete;
		final boolean negative;
		boolean finished;
	}
	static int tuple(int r,int a,int i,int b){return DEFAULT_GRID.tuple(r,a,i,b);}
	static final class GridDefinition{
		GridDefinition(String name_,int[] R_,int[] A_,float[] I_,float[] B_){
			name=name_; R=R_; A=A_; I=I_; B=B_;
			plane=(A.length+1)*(B.length+1);
			cells=R.length*I.length*plane;
			tuples=R.length*A.length*I.length*B.length;
			assert(tuples>0 && B[B.length-1]==Float.MAX_VALUE) : "Each threshold grid must include a finite B-off control";
		}
		int tuple(int r,int a,int i,int b){return ((r*A.length+a)*I.length+i)*B.length+b;}
		final String name;
		final int[] R,A;
		final float[] I,B;
		final int plane,cells,tuples;
	}
}
