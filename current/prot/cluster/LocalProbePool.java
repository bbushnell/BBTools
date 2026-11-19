package prot.cluster;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.FamilySketchBuilder;
import ddl.TailPoolIndex;
import prot.ProteinSequence;
import template.Accumulator;
import template.ThreadWaiter;

/** Parallel construction of the frozen residual-tail index.
 * Sequence storage and baseline must be the same immutable inputs used by the
 * selected-baseline loader. Ordinary builds release temporary sketches after
 * index creation. Recursive fission waves may retain the current residual sketch
 * rows through a bounded cache and rebuild only the sorted ordinal delta.
 * @author Brian Bushnell, Keqing
 */
public final class LocalProbePool{

	public static TailPoolIndex build(final LocalProbeBaseline baseline,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final int threads){
		if(baseline==null || sequences==null || builder==null || threads<1){
			throw new IllegalArgumentException("Residual index requires the frozen baseline, input store, sketch recipe and positive worker count");
		}
		final ClusterSketchCodec.Header h=baseline.sourceHeader;
		if(!builder.matchesGeometry(h.width,h.geometry,h.sketchSeed,h.positionMode,h.partitions,h.overlap,h.blacklistSha80)){
			throw new IllegalArgumentException("Residual-tail sketch geometry/blacklist differs from the unchanged source models");
		}
		return build(baseline.residualOrdinals,sequences,builder,threads);
	}

	/** Builds the next pass's pool from exact sorted original ordinals. The caller
	 * binds the sketch recipe to its committed model payloads before this call.
	 */
	public static TailPoolIndex build(final int[] residualOrdinals,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final int threads){
		final char[][] rows=sketch(residualOrdinals,sequences,builder,threads);
		return new TailPoolIndex(builder.width(),residualOrdinals,rows,threads);
	}

	/** Sketches exact sorted original ordinals into rows parallel to the ordinal vector. */
	private static char[][] sketch(final int[] ordinals,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final int threads){
		validate(ordinals,sequences,builder,threads);
		final Construction construction=new Construction(ordinals,sequences,builder,threads);
		construction.run();
		return construction.rows;
	}

	/** Rejects null context, nonpositive workers and unsorted, duplicate or out-of-range ordinals. */
	private static void validate(final int[] ordinals,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final int threads){
		if(ordinals==null || sequences==null || builder==null || threads<1){
			throw new IllegalArgumentException("Residual construction needs ordered ordinals, sequences, sketch recipe and workers");
		}
		for(int i=0;i<ordinals.length;i++){
			if(ordinals[i]<0 || ordinals[i]>=sequences.length || i>0 && ordinals[i]<=ordinals[i-1]){
				throw new IllegalArgumentException("Residual ordinals must be unique, increasing and inside the original input at slot "+i);
			}
		}
	}

	/** Reuses sketches across one exclusive recursive-fission loop.
	 * The caller owns immutable sequence storage, an immutable sketch recipe, and
	 * sorted unique residual ordinals. The retained cursor is one ordinal clone,
	 * one row-vector, and one char row per current residual; a delta transition can
	 * transiently hold old rows, new rows, rows for additions, and the TailPoolIndex
	 * constructor's defensive clones until publication finishes.
	 */
	public static final class ResidualCache implements AutoCloseable{
		/** Binds a cache to one frozen input array, one sketch recipe object, and one worker count. */
		public ResidualCache(final ProteinSequence[] sequences_,final FamilySketchBuilder builder_,final int threads_){
			if(sequences_==null || builder_==null || threads_<1){
				throw new IllegalArgumentException("Residual cache needs frozen sequences, a sketch recipe and positive workers");
			}
			sequences=sequences_;builder=builder_;threads=threads_;
		}
		/** Builds from a live state that belongs to this cache's frozen input population. */
		TailPoolIndex build(final LiveProbeState state){
			if(state==null || state.sequenceCount()!=sequences.length){
				throw new IllegalArgumentException("Residual cache state must match the frozen input store");
			}
			return build(state.sourceHeader(),state.residualOrdinals());
		}
		/** True only for the exact objects and worker count bound at cache construction. */
		public boolean matchesContext(final ProteinSequence[] sequences_,final FamilySketchBuilder builder_,final int threads_){
			return sequences==sequences_ && builder==builder_ && threads==threads_;
		}
		/** Publishes an immutable pool for these exact sorted original ordinals.
		 * The first call cold-builds every row. Later calls drop removed ordinals,
		 * keep intersection rows by ordinal, and sketch only ordinals newly present
		 * in the residual set.
		 */
		public TailPoolIndex build(final ClusterSketchCodec.Header header,final int[] nextOrdinals){
			if(header==null){throw new IllegalArgumentException("Residual cache needs the exact source sketch header");}
			if(!builder.matchesGeometry(header.width,header.geometry,header.sketchSeed,header.positionMode,
					header.partitions,header.overlap,header.blacklistSha80)){
				throw new IllegalArgumentException("Residual cache sketch geometry/blacklist differs from the live state");
			}
			if(closed){throw new IllegalStateException("Residual cache is closed");}
			validate(nextOrdinals,sequences,builder,threads);
			if(ordinals==null){
				final int[] next=nextOrdinals.clone();
				final char[][] nextRows=sketch(next,sequences,builder,threads);
				ordinals=next;rows=nextRows;
				sketched+=ordinals.length;coldBuilds++;
			}else{
				merge(nextOrdinals);
			}
			return new TailPoolIndex(builder.width(),ordinals,rows,threads);
		}
		/** Retains old rows for the sorted intersection and sketches only sorted additions. */
		private void merge(final int[] nextOrdinals){
			final int[] next=nextOrdinals.clone();
			int old=0, additions=0;
			for(final int nextOrdinal:next){
				while(old<ordinals.length && ordinals[old]<nextOrdinal){old++;}
				if(old>=ordinals.length || ordinals[old]!=nextOrdinal){additions++;}
			}
			final int[] added=new int[additions];
			old=0;additions=0;
			for(final int nextOrdinal:next){
				while(old<ordinals.length && ordinals[old]<nextOrdinal){old++;}
				if(old>=ordinals.length || ordinals[old]!=nextOrdinal){added[additions++]=nextOrdinal;}
			}
			final char[][] addedRows=additions==0 ? null : sketch(added,sequences,builder,threads);
			final char[][] nextRows=new char[next.length][];
			old=0;int fresh=0;
			for(int i=0;i<next.length;i++){
				final int nextOrdinal=next[i];
				while(old<ordinals.length && ordinals[old]<nextOrdinal){old++;}
				nextRows[i]=(old<ordinals.length && ordinals[old]==nextOrdinal) ? rows[old] : addedRows[fresh++];
			}
			ordinals=next;rows=nextRows;sketched+=additions;deltaBuilds++;
		}
		/** Number of residual rows sketched by this cache, excluding reused rows. */
		public long sketched(){return sketched;}
		/** Number of full cold rebuilds retained by this cache. */
		public long coldBuilds(){return coldBuilds;}
		/** Number of delta transitions, including transitions with zero additions. */
		public long deltaBuilds(){return deltaBuilds;}
		/** Releases retained row arrays and permanently rejects future builds. */
		@Override
		public void close(){
			ordinals=null;rows=null;closed=true;
		}
		private final ProteinSequence[] sequences;
		private final FamilySketchBuilder builder;
		private final int threads;
		private int[] ordinals;
		private char[][] rows;
		private long sketched,coldBuilds,deltaBuilds;
		private boolean closed;
	}

	private static final class Construction implements Accumulator<Worker>{
		Construction(int[] ordinals_,ProteinSequence[] sequences_,FamilySketchBuilder builder_,int threads){
			ordinals=ordinals_; sequences=sequences_; rows=new char[ordinals.length][];
			for(int i=0;i<Math.min(threads,ordinals.length);i++){
				workers.add(new Worker(this,builder_.fork()));
			}
		}
		void run(){
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Residual sketch workers did not complete"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null || completed!=ordinals.length){
				throw new IllegalStateException("Residual sketches failed before index publication",failure);
			}
		}
		int claim(){return (int)Math.min(next.getAndAdd(BATCH),ordinals.length);}
		synchronized void fail(Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(Worker worker){
			if(!worker.success){fail(new IllegalStateException("A residual sketch worker stopped without completing its assigned work"));}
			completed=Math.addExact(completed,worker.completed);
		}
		@Override public boolean success(){return failure==null && completed==ordinals.length;}
		@Override public ReadWriteLock rwlock(){return null;}
		final int[] ordinals;
		final ProteinSequence[] sequences;
		final char[][] rows;
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		int completed;
	}

	private static final class Worker extends Thread{
		Worker(Construction owner_,FamilySketchBuilder builder_){owner=owner_;builder=builder_;}
		@Override public void run(){
			try{
				for(int start=owner.claim();start<owner.ordinals.length && owner.failure==null;start=owner.claim()){
					final int end=(int)Math.min((long)start+BATCH,owner.ordinals.length);
					for(int slot=start;slot<end && owner.failure==null;slot++){
						final int ordinal=owner.ordinals[slot];
						if(ordinal<0 || ordinal>=owner.sequences.length || owner.sequences[ordinal]==null){
							throw new IllegalArgumentException("Residual ordinal is absent from frozen input: "+ordinal);
						}
						final char[] sketch=new char[builder.width()];
						builder.sketchEncoded(owner.sequences[ordinal].enc,sketch);
						owner.rows[slot]=sketch; completed++;
					}
				}
				success=owner.failure==null;
			}catch(Throwable cause){owner.fail(cause);}
		}
		final Construction owner;
		final FamilySketchBuilder builder;
		int completed;
		boolean success;
	}

	private LocalProbePool(){}
	private static final int BATCH=200;
}
