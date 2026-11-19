package prot;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.ListNum;

/** Workers pull aligned five-line records and retain private parsed batches.
 * The master combines them in physical file order only after all workers join.
 * No additional producer, work queue, reference hashing or second input pass.
 * @author Keqing
 */
public final class FamilyReferenceLoader{

	private FamilyReferenceLoader(){}

	/** Validates an already constructed roster without serialization or parsing.
	 * Records own immutable payloads; copy only the outer array so caller reuse cannot reorder it.
	 */
	public static Loaded fromMemory(FamilyReferenceCodec.Header header, FamilyReference[] records){
		return fromMemory(header, records, false);
	}

	/** Caller explicitly promises to honor the returned header's geometry and span at every scoring site. */
	public static Loaded fromMemoryMetricBound(FamilyReferenceCodec.Header header, FamilyReference[] records){
		return fromMemory(header, records, true);
	}

	private static Loaded fromMemory(FamilyReferenceCodec.Header header, FamilyReference[] records, boolean metricAware){
		if(header==null || records==null || records.length!=header.totalFamilies || records.length>1000000){
			throw new IllegalArgumentException("Resident references require a matching header and at most1000000 families");
		}
		return validate(header, records.clone(), 0, metricAware);
	}

	/** Uses 200-line ByteFile1 batches, including the five-line header record.
	 * The caller must configure process-wide reader/decompression settings before
	 * loading, without concurrent readers changing those settings. Batch limits
	 * are restored after every load, including failures.
	 */
	public static Loaded load(String path, int threads, int maxFamilies){
		return load(path, threads, maxFamilies, false);
	}

	/** Accepts explicit package metrics as well as legacy files. The consumer must use header.alignmentMetric;
	 * this entry point is deliberately separate so old scoring clients fail closed on new packages.
	 */
	public static Loaded loadMetricBound(String path, int threads, int maxFamilies){
		return load(path, threads, maxFamilies, true);
	}

	private static Loaded load(String path, int threads, int maxFamilies, boolean metricAware){
		if(path==null || path.isEmpty() || threads<1 || maxFamilies<0){
			throw new IllegalArgumentException("Require a reference path, positive threads and nonnegative maxFamilies");
		}
		final long start=System.nanoTime();
		final Metrics metrics=new Metrics();
		final int oldLines=ByteFile.TARGET_LIST_SIZE, oldBytes=ByteFile.TARGET_LIST_BYTES;
		ByteFile input=null;
		Loaded loaded;
		Throwable failure=null;
		try{
			ByteFile.TARGET_LIST_SIZE=LINES_PER_BLOCK;
			ByteFile.TARGET_LIST_BYTES=100_000_000;
			// BF1 honors the fixed line limit in both scalar and SIMD paths.
			// Its inherited nextList serializes claims and supplies firstRecordNum.
			input=new ByteFile1(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));
			final InputStream backend=input.is();
			metrics.inputClass=backend==null ? "unavailable" : backend.getClass().getName();
			metrics.open=System.nanoTime()-start;
			final int count=Math.min(threads, 1+maxFamilies/FAMILIES_PER_BLOCK);
			final ProcessThread[] workers=readWorkers(input, count, maxFamilies, metrics, metricAware);
			final long mergeStart=System.nanoTime();
			FamilyReferenceCodec.Header header=null;
			int blocks=0;
			for(ProcessThread worker : workers){
				blocks=Math.addExact(blocks, worker.blocks.size());
				if(worker.header!=null){
					if(header!=null){throw new IllegalArgumentException("Duplicate reference header batch");}
					header=worker.header;
				}
			}
			if(header==null){throw new IllegalArgumentException("Missing five-line reference header");}
			final FamilyReference[] all=new FamilyReference[header.totalFamilies];
			for(ProcessThread worker : workers){
				for(ParsedBlock block : worker.blocks){
					if(block.first>all.length-block.records.length){
						throw new IllegalArgumentException("Family block exceeds header census at ordinal "+block.first);
					}
					for(int i=0; i<block.records.length; i++){
						final int ordinal=block.first+i;
						if(all[ordinal]!=null){throw new IllegalArgumentException("Duplicate family ordinal "+ordinal);}
						all[ordinal]=block.records[i];
					}
				}
			}
			metrics.merge=System.nanoTime()-mergeStart;
			final long censusStart=System.nanoTime();
			loaded=validate(header, all, blocks, metricAware);
			metrics.census=System.nanoTime()-censusStart;
		}catch(Throwable problem){failure=problem; throw propagate(problem);}
		finally{
			final long closeStart=System.nanoTime();
			try{
				if(input!=null && input.close()){throw new IllegalStateException("Error reading or closing family reference: "+path);}
			}catch(Throwable closeFailure){
				if(failure!=null){failure.addSuppressed(closeFailure);}else{throw propagate(closeFailure);}
			}finally{
				metrics.close=System.nanoTime()-closeStart;
				ByteFile.TARGET_LIST_SIZE=oldLines; ByteFile.TARGET_LIST_BYTES=oldBytes;
			}
		}
		return new Loaded(loaded, new Timing(metrics, System.nanoTime()-start));
	}

	/** All worker state is private until join; a first failure stops further claims. */
	private static ProcessThread[] readWorkers(ByteFile input, int count, int maxFamilies, Metrics metrics, boolean metricAware){
		assert(input!=null && count>0) : "Worker claims require an open reader and at least one parser";
		final long start=System.nanoTime();
		final AtomicReference<Throwable> failure=new AtomicReference<>();
		final ProcessThread[] workers=new ProcessThread[count];
		for(int i=0; i<count; i++){workers[i]=new ProcessThread(input, maxFamilies, failure, i, metricAware);}
		int started=0;
		boolean interrupted=false;
		try{
			if(count==1){workers[0].run();}
			else{for(ProcessThread worker : workers){worker.start(); started++;}}
		}catch(Throwable problem){failure.compareAndSet(null, problem);}
		finally{
			for(int i=0; i<started; i++){
				boolean joined=false;
				while(!joined){
					try{workers[i].join(); joined=true;}
					catch(InterruptedException problem){interrupted=true; failure.compareAndSet(null, problem);}
				}
			}
			for(ProcessThread worker : workers){
				metrics.readWorker+=worker.readNanos;
				metrics.parseWorker+=worker.parseNanos;
				metrics.codec.add(worker.codecTiming);
			}
			metrics.workers=System.nanoTime()-start;
			if(interrupted){Thread.currentThread().interrupt();}
		}
		if(failure.get()!=null){throw propagate(failure.get());}
		return workers;
	}

	/** Checks the whole-file identity contract and builds the live-only search view. */
	private static Loaded validate(FamilyReferenceCodec.Header header, FamilyReference[] all, int blocks, boolean metricAware){
		assert(all.length==header.totalFamilies) : "Result storage is allocated from the validated header census";
		if(!metricAware && !header.alignmentMetric.same(FamilyAlignmentMetric.LEGACY)){
			throw new IllegalArgumentException("Metric-bound references require an explicitly metric-aware scoring consumer");
		}
		final boolean[] ranks=new boolean[all.length];
		final boolean explicitAlignment=FamilyReferenceCodec.alignmentStateSchema(header.schema);
		final HashSet<String> ids=new HashSet<>();
		final FamilyReference[] live=new FamilyReference[header.liveFamilies];
		int next=0;
		long members=0, residues=0, words=0;
		for(int i=0; i<all.length; i++){
			final FamilyReference family=all[i];
			if(family==null){throw new IllegalArgumentException("Missing family at file ordinal "+i+" of "+all.length);}
			if(header.scalarZ==0 && family.scalarMask!=0){throw new IllegalArgumentException("Scalar-enabled family in a z=0 package: "+family.id);}
			if(family.rank>=ranks.length || ranks[family.rank]){throw new IllegalArgumentException("Duplicate or out-of-range family rank: "+family.rank);}
			ranks[family.rank]=true;
			if(!ids.add(family.id)){throw new IllegalArgumentException("Duplicate family ID: "+family.id);}
			if(family.members>0){
				if(family.ddl.length!=header.width){
					throw new IllegalArgumentException("Reference DDL width disagrees with package header: "+family.id+
						" observed="+family.ddl.length+" expected="+header.width);
				}
				if(family.explicitAlignmentCalibration!=explicitAlignment){
					throw new IllegalArgumentException("Reference alignment-state row disagrees with package header: "+family.id);
				}
				if(family.spectrum.length!=header.spectrumWidth){
					throw new IllegalArgumentException("Reference spectrum width disagrees with package header: "+family.id+
						" observed="+family.spectrum.length+" expected="+header.spectrumWidth);
				}
				if(next==live.length){throw new IllegalArgumentException("Live family count exceeds header census");}
				live[next++]=family;
			}
			try{
				members=Math.addExact(members, family.members);
				residues=Math.addExact(residues, family.residues);
				words=Math.addExact(words, family.validWords);
			}catch(ArithmeticException problem){throw new IllegalArgumentException("Family count total overflow at "+family.id, problem);}
		}
		if(next!=header.liveFamilies || all.length-next!=header.emptyFamilies){
			throw new IllegalArgumentException("Reference live/empty census does not match header");
		}
		return new Loaded(header, all, live, members, residues, words, blocks);
	}

	private static RuntimeException propagate(Throwable problem){
		if(problem instanceof Error){throw (Error)problem;}
		if(problem instanceof RuntimeException){return (RuntimeException)problem;}
		return new IllegalStateException("Family reference loading failed", problem);
	}

	/** Published arrays: all follows physical order; live excludes only zero-member rows. */
	public static final class Loaded{
		private Loaded(Loaded data, Timing timing){
			header=data.header; all=data.all; live=data.live; members=data.members;
			residues=data.residues; validWords=data.validWords; blocks=data.blocks; this.timing=timing;
		}
		private Loaded(FamilyReferenceCodec.Header header, FamilyReference[] all, FamilyReference[] live,
				long members, long residues, long validWords, int blocks){
			this.header=header; this.all=all; this.live=live; this.members=members;
			this.residues=residues; this.validWords=validWords; this.blocks=blocks; timing=null;
		}
		public final FamilyReferenceCodec.Header header;
		public final FamilyReference[] all, live;
		public final long members, residues, validWords;
		/** Physical 200-line batches, including the batch containing the header. */
		public final int blocks;
		public final Timing timing;
	}

	/** Master phases are disjoint. Worker read/parse sums overlap each other and
	 * workerPhaseWallNanos; reading includes the ByteFile monitor and decompression waits. */
	public static final class Timing{
		private Timing(Metrics m, long total){
			totalWallNanos=total; openWallNanos=m.open; workerPhaseWallNanos=m.workers;
			readBatchWorkerNanos=m.readWorker; parseWorkerNanos=m.parseWorker;
			mergeWallNanos=m.merge; censusWallNanos=m.census; closeWallNanos=m.close;
			parsedRecords=m.codec.records; metadataGatesWorkerNanos=m.codec.metadataGatesNanos;
			consensusDdlWorkerNanos=m.codec.consensusDdlNanos; spectrumWorkerNanos=m.codec.spectrumNanos;
			dtoNormWorkerNanos=m.codec.dtoNormNanos;
			inputClass=m.inputClass;
			assert(total>=m.open+m.workers+m.merge+m.census+m.close) : "Disjoint master phases must fit within loader wall time";
		}
		public final long totalWallNanos,openWallNanos,workerPhaseWallNanos,mergeWallNanos,censusWallNanos,closeWallNanos;
		public final long readBatchWorkerNanos,parseWorkerNanos,parsedRecords;
		public final long metadataGatesWorkerNanos,consensusDdlWorkerNanos,spectrumWorkerNanos,dtoNormWorkerNanos;
		public final String inputClass;
	}

	private static final class Metrics{
		long open,workers,readWorker,parseWorker,merge,census,close;
		final FamilyReferenceCodec.ParseTiming codec=new FamilyReferenceCodec.ParseTiming();
		String inputClass;
	}

	private static final class ProcessThread extends Thread{
		ProcessThread(ByteFile input, int maxFamilies, AtomicReference<Throwable> failure, int id, boolean metricAware){
			super("FamilyReferenceLoader-"+id);
			this.input=input; this.maxFamilies=maxFamilies; this.failure=failure;
			this.metricAware=metricAware;
		}
		@Override
		public void run(){
			try{
				final LineParser1 parser=new LineParser1('\t');
				while(failure.get()==null){
					if(Thread.currentThread().isInterrupted()){throw new InterruptedException("Reference parser interrupted");}
					final long readStart=System.nanoTime();
					final ListNum<byte[]> lines;
					try{lines=input.nextList();}finally{readNanos+=System.nanoTime()-readStart;}
					if(lines==null){break;}
					final long parseStart=System.nanoTime();
					try{parse(lines, parser);}finally{parseNanos+=System.nanoTime()-parseStart;}
				}
			}catch(Throwable problem){failure.compareAndSet(null, problem);}
		}
		private void parse(ListNum<byte[]> batch, LineParser1 parser){
			final int size=batch.size();
			if(size<1 || size>LINES_PER_BLOCK || size%FamilyReferenceCodec.LINES_PER_FAMILY!=0 ||
				batch.firstRecordNum<0 || batch.firstRecordNum%FamilyReferenceCodec.LINES_PER_FAMILY!=0){
				throw new IllegalArgumentException("Incomplete or misaligned five-line batch at line "+batch.firstRecordNum+", lines="+size);
			}
			final byte[][] lines=batch.list.toArray(new byte[size][]);
			final boolean first=batch.firstRecordNum==0;
			final int offset=first ? FamilyReferenceCodec.LINES_PER_FAMILY : 0;
			if(first){
				header=metricAware ? FamilyReferenceCodec.parseMetricBoundHeader(lines, 0, parser) :
					FamilyReferenceCodec.parseHeader(lines, 0, parser);
				if(header.totalFamilies>maxFamilies){throw new IllegalArgumentException("Reference census exceeds maxFamilies="+maxFamilies);}
			}
			final long ordinal=first ? 0 : batch.firstRecordNum/FamilyReferenceCodec.LINES_PER_FAMILY-1;
			final int count=(size-offset)/FamilyReferenceCodec.LINES_PER_FAMILY;
			if(ordinal>maxFamilies-count){throw new IllegalArgumentException("Reference data exceeds maxFamilies at ordinal "+ordinal);}
			final FamilyReference[] records=new FamilyReference[count];
			for(int i=0; i<count; i++){
				try{records[i]=FamilyReferenceCodec.parseFiveLinesAnySpectrum(lines, offset+i*FamilyReferenceCodec.LINES_PER_FAMILY, parser, codecTiming);}
				catch(IllegalArgumentException problem){throw new IllegalArgumentException("Invalid family at file ordinal "+(ordinal+i)+": "+problem.getMessage(), problem);}
			}
			blocks.add(new ParsedBlock((int)ordinal, records));
		}
		private final ByteFile input;
		private final int maxFamilies;
		private final AtomicReference<Throwable> failure;
		private final boolean metricAware;
		private final ArrayList<ParsedBlock> blocks=new ArrayList<>();
		private FamilyReferenceCodec.Header header;
		private long readNanos,parseNanos;
		private final FamilyReferenceCodec.ParseTiming codecTiming=new FamilyReferenceCodec.ParseTiming();
	}

	private static final class ParsedBlock{
		ParsedBlock(int first, FamilyReference[] records){
			assert(first>=0 && records!=null) : "A parsed batch retains its physical family origin and owned records";
			this.first=first; this.records=records;
		}
		final int first;
		final FamilyReference[] records;
	}

	public static final int FAMILIES_PER_BLOCK=40;
	public static final int LINES_PER_BLOCK=200;
}
