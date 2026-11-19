package prot;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.IntHashMap;
import structures.ListNum;

/** Workers pull aligned ten-line compact family records and retain private parsed batches.
 * The master combines them in physical file order only after all workers join.
 * Header census keys are {@code families}, {@code live}, and {@code empty}.
 * This loader checks record syntax, unique identities, dense rank coverage and
 * count conservation. It does not open alphabet/blacklist/matrix resources or
 * enforce their header bindings; the callable package performs those checks.
 * Rank order need not equal physical order. Returned all/live arrays and their
 * parallel metadata use physical order, with zero-member rows omitted from live.
 * Loading temporarily changes process-wide ByteFile batch limits. Callers must
 * serialize loads and other readers that depend on those limits; this class
 * does not provide a process-wide lock. Returned arrays are exposed mutable
 * storage and require caller ownership despite their final field references.
 * @author Yelan
 */
public final class FamilyCompactLoader{

	/** Static loading entry point only. */
	private FamilyCompactLoader(){}

	/** Uses 200-line ByteFile1 batches, including the ten-line header record.
	 * The caller must configure process-wide reader/decompression settings before
	 * loading, without concurrent readers changing those settings. Batch limits
	 * are restored after every load, including failures.
	 * The input is closed after worker joins. No partial Loaded result is published
	 * on failure; a close failure is suppressed on an earlier failure when present.
	 * @param path Nonempty input path interpreted by the ordinary ByteFile/FileFormat APIs.
	 * @param threads Positive requested parser count, capped using maxFamilies; one runs inline.
	 * @param maxFamilies Nonnegative physical-family cap, including empty records.
	 * @return Parsed families, parallel metadata, computed totals and timing in nanoseconds.
	 */
	public static Loaded load(String path,int threads,int maxFamilies){
		if(path==null || path.isEmpty() || threads<1 || maxFamilies<0){
			throw new IllegalArgumentException("Require a compact reference path, positive threads and nonnegative maxFamilies");
		}
		final long start=System.nanoTime();
		final Metrics metrics=new Metrics();
		final int oldLines=ByteFile.TARGET_LIST_SIZE,oldBytes=ByteFile.TARGET_LIST_BYTES;
		ByteFile input=null;
		Loaded loaded;
		Throwable failure=null;
		try{
			ByteFile.TARGET_LIST_SIZE=FamilyCompactCodec.LINES_PER_BATCH;
			ByteFile.TARGET_LIST_BYTES=100_000_000;
			// BF1 honors the fixed line limit in both scalar and SIMD paths.
			// Its inherited nextList serializes claims and supplies firstRecordNum.
			input=new ByteFile1(FileFormat.testInput(path,FileFormat.TEXT,null,false,true));
			final InputStream backend=input.is();
			metrics.inputClass=backend==null ? "unavailable" : backend.getClass().getName();
			metrics.open=System.nanoTime()-start;
			final int count=Math.min(threads,1+maxFamilies/FamilyCompactCodec.FAMILIES_PER_BATCH);
			final ProcessThread[] workers=readWorkers(input,count,maxFamilies,metrics);
			final long mergeStart=System.nanoTime();
			FamilyCompactCodec.Header header=null;
			int blocks=0;
			for(ProcessThread worker : workers){
				blocks=Math.addExact(blocks,worker.blocks.size());
				if(worker.header!=null){
					if(header!=null){throw new IllegalArgumentException("Duplicate compact header batch");}
					header=worker.header;
				}
			}
			if(header==null){throw new IllegalArgumentException("Missing ten-line compact header");}
			final int total=headerInt(header,"families");
			if(total>maxFamilies){throw new IllegalArgumentException("Compact census exceeds maxFamilies="+maxFamilies);}
			final LoadedFamily[] all=new LoadedFamily[total];
			for(ProcessThread worker : workers){
				for(ParsedBlock block : worker.blocks){
					if(block.first>all.length-block.records.length){
						throw new IllegalArgumentException("Compact family block exceeds header census at ordinal "+block.first);
					}
					for(int i=0; i<block.records.length; i++){
						final int ordinal=block.first+i;
						if(all[ordinal]!=null){throw new IllegalArgumentException("Duplicate compact family ordinal "+ordinal);}
						all[ordinal]=block.records[i];
					}
				}
			}
			metrics.merge=System.nanoTime()-mergeStart;
			final long censusStart=System.nanoTime();
			loaded=validate(header,all,blocks);
			metrics.census=System.nanoTime()-censusStart;
		}catch(Throwable problem){failure=problem; throw propagate(problem);}
		finally{
			final long closeStart=System.nanoTime();
			try{
				if(input!=null && input.close()){
					throw new IllegalStateException("Error reading or closing compact reference: "+path);
				}
			}catch(Throwable closeFailure){
				if(failure!=null){failure.addSuppressed(closeFailure);}else{throw propagate(closeFailure);}
			}finally{
				metrics.close=System.nanoTime()-closeStart;
				ByteFile.TARGET_LIST_SIZE=oldLines; ByteFile.TARGET_LIST_BYTES=oldBytes;
			}
		}
		return new Loaded(loaded,new Timing(metrics,System.nanoTime()-start));
	}

	/** All worker state is private until join; a first failure stops further claims.
	 * Joins every started worker after startup/parse failure. Interrupted joins are
	 * retried and the calling thread's interrupt flag is restored before propagating.
	 * Already-claimed reads/parses may finish after another worker reports failure.
	 */
	private static ProcessThread[] readWorkers(ByteFile input,int count,int maxFamilies,Metrics metrics){
		assert(input!=null && count>0) : "Compact worker claims require an open reader and at least one parser";
		final long start=System.nanoTime();
		final AtomicReference<Throwable> failure=new AtomicReference<>();
		final ProcessThread[] workers=new ProcessThread[count];
		for(int i=0; i<count; i++){workers[i]=new ProcessThread(input,maxFamilies,failure,i);}
		int started=0;
		boolean interrupted=false;
		try{
			if(count==1){workers[0].run();}
			else{for(ProcessThread worker : workers){worker.start(); started++;}}
		}catch(Throwable problem){failure.compareAndSet(null,problem);}
		finally{
			for(int i=0; i<started; i++){
				boolean joined=false;
				while(!joined){
					try{workers[i].join(); joined=true;}
					catch(InterruptedException problem){interrupted=true; failure.compareAndSet(null,problem);}
				}
			}
			for(ProcessThread worker : workers){
				metrics.readWorker+=worker.readNanos;
				metrics.parseWorker+=worker.parseNanos;
			}
			metrics.workers=System.nanoTime()-start;
			if(interrupted){Thread.currentThread().interrupt();}
		}
		if(failure.get()!=null){throw propagate(failure.get());}
		return workers;
	}

	/** Checks the whole-file identity/census contract and builds the live-only view.
	 * Ranks uniquely cover 0..families-1; family IDs and nonnegative stable IDs are
	 * unique. Totals are computed from records, with GC-measured equal to members;
	 * this does not compare arbitrary package-total header fields or external maps.
	 */
	private static Loaded validate(FamilyCompactCodec.Header header,LoadedFamily[] records,int blocks){
		final int total=headerInt(header,"families"),liveCount=headerInt(header,"live"),emptyCount=headerInt(header,"empty");
		if(total!=records.length || liveCount<0 || emptyCount<0 || liveCount+emptyCount!=total){
			throw new IllegalArgumentException("Compact header census mismatch");
		}
		final boolean[] ranks=new boolean[records.length];
		final HashSet<String> ids=new HashSet<String>();
		final IntHashMap stableIds=new IntHashMap(records.length*2+1);
		final FamilyReference[] all=new FamilyReference[records.length],live=new FamilyReference[liveCount];
		final int[] allStableIds=new int[records.length],liveStableIds=new int[liveCount];
		final double[] lengthMedians=new double[liveCount],kdMeans=new double[liveCount];
		final double[] gcMeans=new double[liveCount],gcSampleSds=new double[liveCount];
		final long[] gcMeasured=new long[liveCount];
		int next=0;
		long members=0,residues=0,words=0,gc=0;
		for(int i=0; i<records.length; i++){
			final LoadedFamily record=records[i];
			if(record==null){
				throw new IllegalArgumentException("Missing compact family at file ordinal "+i+
					" of "+records.length);
			}
			final FamilyReference family=record.family;
			if(family.rank>=ranks.length || ranks[family.rank]){
				throw new IllegalArgumentException("Duplicate or out-of-range compact family rank: "+family.rank);
			}
			ranks[family.rank]=true;
			if(!ids.add(family.id)){throw new IllegalArgumentException("Duplicate compact family ID: "+family.id);}
			if(stableIds.get(record.stableId)>=0){
				throw new IllegalArgumentException("Duplicate compact stable ID: "+record.stableId);
			}
			stableIds.put(record.stableId,i);
			allStableIds[i]=record.stableId;
			all[i]=family;
			if(family.members>0){
				if(next==live.length){throw new IllegalArgumentException("Live compact family count exceeds header census");}
				live[next]=all[i];
				liveStableIds[next]=record.stableId;
				lengthMedians[next]=record.lengthMedian;
				kdMeans[next]=record.kdMean;
				gcMeasured[next]=record.measured;
				gcMeans[next]=record.meanGc;
				gcSampleSds[next]=record.sampleSd;
				next++;
			}
			try{
				members=Math.addExact(members,family.members);
				residues=Math.addExact(residues,family.residues);
				words=Math.addExact(words,family.validWords);
				gc=Math.addExact(gc,record.measured);
			}catch(ArithmeticException problem){
				throw new IllegalArgumentException("Compact family count overflow at "+family.id,problem);
			}
		}
		if(next!=liveCount || all.length-next!=emptyCount){
			throw new IllegalArgumentException("Compact live/empty census does not match header");
		}
		if(gc!=members){throw new IllegalArgumentException("Compact GC measured total does not equal member total");}
		return new Loaded(header,all,live,allStableIds,liveStableIds,lengthMedians,kdMeans,
			gcMeasured,gcMeans,gcSampleSds,members,residues,words,blocks);
	}

	/** Maps explicit wire calibration availability to the in-memory family enum. */
	private static FamilyReference.AlignmentCalibration alignment(FamilyCompactCodec.AlignmentCalibration value){
		if(value==FamilyCompactCodec.AlignmentCalibration.ENABLED){return FamilyReference.AlignmentCalibration.ENABLED;}
		if(value==FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE){
			return FamilyReference.AlignmentCalibration.UNAVAILABLE;
		}
		throw new IllegalArgumentException("Missing compact alignment state");
	}

	/** Requires a nonempty header value containing decimal digits within the int range. */
	private static int headerInt(FamilyCompactCodec.Header header,String key){
		final String value=header.get(key);
		if(value==null){throw new IllegalArgumentException("Compact header missing "+key);}
		long x=0;
		for(int i=0; i<value.length(); i++){
			final int digit=value.charAt(i)-'0';
			if(digit<0 || digit>9 || x>(Integer.MAX_VALUE-digit)/10){
				throw new IllegalArgumentException("Compact header "+key+" must be a nonnegative int");
			}
			x=x*10+digit;
		}
		return (int)x;
	}

	/** Preserves Error/runtime exception identity and wraps checked worker/I/O failures. */
	private static RuntimeException propagate(Throwable problem){
		if(problem instanceof Error){throw (Error)problem;}
		if(problem instanceof RuntimeException){return (RuntimeException)problem;}
		return new IllegalStateException("Compact family reference loading failed",problem);
	}

	/** Published arrays: all follows physical order; live excludes only zero-member rows.
	 * The live entries reference the same FamilyReference objects as all. Stable IDs
	 * are parallel to their corresponding all/live array; lengths, KD and GC arrays
	 * are parallel to live only. None is indexed implicitly by family rank.
	 * Arrays are not cloned on access. Callers keep shared results stable while in use.
	 */
	public static final class Loaded{
		/** Adds completed timing without copying any validated family or metadata array. */
		private Loaded(Loaded data,Timing timing){
			header=data.header; all=data.all; live=data.live; allStableIds=data.allStableIds;
			liveStableIds=data.liveStableIds; lengthMedians=data.lengthMedians; kdMeans=data.kdMeans;
			gcMeasured=data.gcMeasured; gcMeans=data.gcMeans; gcSampleSds=data.gcSampleSds;
			members=data.members; residues=data.residues; validWords=data.validWords; blocks=data.blocks;
			this.timing=timing;
		}
		/** Takes ownership of merged arrays; timing is attached after input closure. */
		private Loaded(FamilyCompactCodec.Header header,FamilyReference[] all,FamilyReference[] live,
				int[] allStableIds,int[] liveStableIds,double[] lengthMedians,double[] kdMeans,
				long[] gcMeasured,double[] gcMeans,double[] gcSampleSds,
				long members,long residues,long validWords,int blocks){
			this.header=header; this.all=all; this.live=live; this.allStableIds=allStableIds;
			this.liveStableIds=liveStableIds; this.lengthMedians=lengthMedians; this.kdMeans=kdMeans;
			this.gcMeasured=gcMeasured; this.gcMeans=gcMeans; this.gcSampleSds=gcSampleSds;
			this.members=members; this.residues=residues; this.validWords=validWords;
			this.blocks=blocks; timing=null;
		}
		public final FamilyCompactCodec.Header header;
		public final FamilyReference[] all,live;
		public final int[] allStableIds,liveStableIds;
		public final double[] lengthMedians,kdMeans,gcMeans,gcSampleSds;
		public final long[] gcMeasured;
		public final long members,residues,validWords;
		/** Physical 200-line batches, including the batch containing the header. */
		public final int blocks;
		public final Timing timing;
	}

	/** Master phases are disjoint. Worker read/parse sums overlap each other and
	 * workerPhaseWallNanos; reading includes the ByteFile monitor and decompression waits.
	 * All numeric timings are nanoseconds. Worker sums are elapsed worker intervals,
	 * not CPU time, and may exceed total wall time when workers overlap.
	 */
	public static final class Timing{
		/** Captures completed loader counters after closure. */
		private Timing(Metrics m,long total){
			totalWallNanos=total; openWallNanos=m.open; workerPhaseWallNanos=m.workers;
			readBatchWorkerNanos=m.readWorker; parseWorkerNanos=m.parseWorker;
			mergeWallNanos=m.merge; censusWallNanos=m.census; closeWallNanos=m.close;
			inputClass=m.inputClass;
			assert(total>=m.open+m.workers+m.merge+m.census+m.close) :
				"Disjoint compact loader master phases must fit within wall time";
		}
		public final long totalWallNanos,openWallNanos,workerPhaseWallNanos,mergeWallNanos;
		public final long censusWallNanos,closeWallNanos;
		public final long readBatchWorkerNanos,parseWorkerNanos;
		public final String inputClass;
	}

	/** Per-load counters; worker contributions are read only after joins. */
	private static final class Metrics{
		long open,workers,readWorker,parseWorker,merge,census,close;
		String inputClass;
	}

	/** One parser owns its decoded blocks; the shared reader serializes batch claims. */
	private static final class ProcessThread extends Thread{
		/** Binds shared input/failure signaling and the physical-family cap. */
		ProcessThread(ByteFile input,int maxFamilies,AtomicReference<Throwable> failure,int id){
			super("FamilyCompactLoader-"+id);
			this.input=input; this.maxFamilies=maxFamilies; this.failure=failure;
		}
		/** Claims until EOF or first failure and retains only decoded block payloads. */
		@Override
		public void run(){
			try{
				final LineParser1 parser=new LineParser1('\t');
				while(failure.get()==null){
					if(Thread.currentThread().isInterrupted()){throw new InterruptedException("Compact parser interrupted");}
					final long readStart=System.nanoTime();
					final ListNum<byte[]> lines;
					try{lines=input.nextList();}finally{readNanos+=System.nanoTime()-readStart;}
					if(lines==null){break;}
					final long parseStart=System.nanoTime();
					try{parse(lines,parser);}finally{parseNanos+=System.nanoTime()-parseStart;}
				}
			}catch(Throwable problem){failure.compareAndSet(null,problem);}
		}
		/** Validates physical alignment and derives family ordinals from the first line.
		 * The first batch contains the header and at most nineteen families; later
		 * full batches contain twenty. Raw line arrays are not kept in parsed blocks.
		 */
		private void parse(ListNum<byte[]> batch,LineParser1 parser){
			final int size=batch.size();
			if(size<1 || size>FamilyCompactCodec.LINES_PER_BATCH ||
					size%FamilyCompactCodec.LINES_PER_RECORD!=0 || batch.firstRecordNum<0 ||
					batch.firstRecordNum%FamilyCompactCodec.LINES_PER_RECORD!=0){
				throw new IllegalArgumentException("Incomplete or misaligned ten-line batch at line "+
					batch.firstRecordNum+", lines="+size);
			}
			final byte[][] lines=batch.list.toArray(new byte[size][]);
			final boolean first=batch.firstRecordNum==0;
			final int offset=first ? FamilyCompactCodec.LINES_PER_RECORD : 0;
			if(first){header=FamilyCompactCodec.parseHeader(lines,0,parser);}
			final long ordinal=first ? 0 : batch.firstRecordNum/FamilyCompactCodec.LINES_PER_RECORD-1;
			final int count=(size-offset)/FamilyCompactCodec.LINES_PER_RECORD;
			if(ordinal>maxFamilies-count){
				throw new IllegalArgumentException("Compact data exceeds maxFamilies at ordinal "+ordinal);
			}
			final LoadedFamily[] records=new LoadedFamily[count];
			for(int i=0; i<count; i++){
				try{
					records[i]=FamilyCompactCodec.parseRecord(lines,
						offset+i*FamilyCompactCodec.LINES_PER_RECORD,parser,LoadedFamily.SINK);
				}catch(IllegalArgumentException problem){
					throw new IllegalArgumentException("Invalid compact family at file ordinal "+(ordinal+i)+": "+
						problem.getMessage(),problem);
				}
			}
			blocks.add(new ParsedBlock((int)ordinal,records));
		}
		private final ByteFile input;
		private final int maxFamilies;
		private final AtomicReference<Throwable> failure;
		private final ArrayList<ParsedBlock> blocks=new ArrayList<ParsedBlock>();
		private FamilyCompactCodec.Header header;
		private long readNanos,parseNanos;
	}

	/** Worker-owned decoded family segment positioned by physical file ordinal. */
	private static final class ParsedBlock{
		/** Retains the decoded array for the master's ordered merge. */
		ParsedBlock(int first,LoadedFamily[] records){
			assert(first>=0 && records!=null) : "A compact batch retains its physical family origin and owned records";
			this.first=first; this.records=records;
		}
		final int first;
		final LoadedFamily[] records;
	}

	/** Temporary owner of one family and the metadata not stored in FamilyReference. */
	private static final class LoadedFamily{
		/** Retains the decoded family and its stable ID/ranking/GC values. */
		LoadedFamily(FamilyReference family,int stableId,double lengthMedian,double kdMean,
				long measured,double meanGc,double sampleSd){
			assert(family!=null && stableId>=0) : "Compact family loading needs a family and stable ID";
			this.family=family; this.stableId=stableId; this.lengthMedian=lengthMedian;
			this.kdMean=kdMean; this.measured=measured; this.meanGc=meanGc; this.sampleSd=sampleSd;
		}
		/** Constructs final family storage directly from codec callbacks, avoiding a Record wrapper. */
		static final FamilyCompactCodec.RecordSink<LoadedFamily> SINK=
			new FamilyCompactCodec.RecordSink<LoadedFamily>(){
				@Override
				public LoadedFamily empty(String familyId,int rank,int stableId){
					return new LoadedFamily(FamilyReference.empty(familyId,rank),stableId,0,0,0,0,0);
				}
				@Override
				public LoadedFamily live(String familyId,int rank,int stableId,long members,
						long residues,long validWords,float cosineLower,float identityLower,
						float blosumLower,float cosineAccept,float identityUpper,float blosumFinal,
						int scalarMask,FamilyCompactCodec.AlignmentCalibration alignmentCalibration,
						double[] scalarBounds,double lengthMedian,double kdMean,long measured,
						double meanGc,double sampleSd,byte[] consensusRaw,char[] ddl,float[] monomers){
					final FamilyReference family=new FamilyReference(familyId,rank,members,residues,
						validWords,new ProteinSequence(familyId,consensusRaw),ddl,monomers,cosineLower,
						identityLower,blosumLower,cosineAccept,identityUpper,blosumFinal,
						consensusRaw,scalarMask,scalarBounds,alignment(alignmentCalibration));
					return new LoadedFamily(family,stableId,lengthMedian,kdMean,measured,meanGc,sampleSd);
				}
			};
		final FamilyReference family;
		final int stableId;
		final double lengthMedian,kdMean,meanGc,sampleSd;
		final long measured;
	}
}
