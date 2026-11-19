package prot.cluster;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.atomic.AtomicBoolean;

/** Append-only primitive member list shared by regenerated family objects. */
final class FamilyMembers {

	FamilyMembers(){this(8);}

	FamilyMembers(int initial){
		if(initial<1){initial=1;}
		chunks=new AtomicReference<AtomicIntegerArray[]>(new AtomicIntegerArray[Math.max(1, tableLength(initial))]);
	}

	int size(){checkFailure(); return published.get();}

	int get(int index){
		final int size=published.get();
		checkFailure();
		if(index<0 || index>=size){throw new ArrayIndexOutOfBoundsException(index);}
		return raw(index);
	}

	int[] toArray(){
		final int size=published.get();
		checkFailure();
		final int[] copy=new int[size];
		for(int i=0; i<size; i++){copy[i]=raw(i);}
		return copy;
	}

	int add(int value){
		checkFailure();
		if(value<0 || value==Integer.MAX_VALUE){throw new IllegalArgumentException("Invalid member index");}
		final int slot=claimed.getAndIncrement();
		try{
			if(slot<0 || slot>=MAX_MEMBERS){throw new IllegalStateException("Family membership count exceeded supported capacity");}
			ensure(slot);
			chunks.get()[slot>>>CHUNK_SHIFT].set(slot&CHUNK_MASK, value+1);
			advancePublished();
		}catch(Throwable problem){
			failure.compareAndSet(null, problem);
			throw propagate(problem);
		}
		return slot;
	}

	/** Serial reassignment clears ownership after every worker from the prior pass has joined. */
	void clear(){
		checkFailure();
		final RankingState state=ranking;
		if(claimed.get()!=published.get() || state!=null && state.publishing.get()){
			throw new IllegalStateException("Cannot clear while a member append is unpublished");
		}
		if(state!=null){state.statistics.clear();}// Published descriptors retain last-nonempty provenance.
		chunks.set(new AtomicIntegerArray[1]);
		claimed.set(0);
		published.set(0);
	}

	private int raw(int index){
		return chunks.get()[index>>>CHUNK_SHIFT].get(index&CHUNK_MASK)-1;
	}

	private void ensure(int slot){
		final int chunk=slot>>>CHUNK_SHIFT;
		for(;;){
			final AtomicIntegerArray[] table=chunks.get();
			if(chunk<table.length && table[chunk]!=null){return;}
			final int newLength=chunk<table.length ? table.length : tableLength(slot+1);
			final AtomicIntegerArray[] copy=new AtomicIntegerArray[newLength];
			System.arraycopy(table, 0, copy, 0, table.length);
			if(copy[chunk]==null){copy[chunk]=new AtomicIntegerArray(CHUNK_SIZE);}
			if(chunks.compareAndSet(table, copy)){return;}
		}
	}

	private void advancePublished(){
		final RankingState state=ranking;
		if(state!=null){advanceRanked(state); return;}
		for(;;){
			checkFailure();
			final int size=published.get();
			if(!ready(size)){return;}
			published.compareAndSet(size, size+1);
		}
	}

	/** Attach once before this membership is shared with workers; exported descriptors do not add ordinals. */
	void enableRanking(final MemberRankInputs inputs, final int layer, final MemberRankStatistics.Snapshot seed){
		checkFailure();
		if(inputs==null || layer<0 || ranking!=null || claimed.get()!=published.get()){
			throw new IllegalStateException("Attach ranking once to quiescent membership with its original input cache and nonnegative layer");
		}
		final RankingState state=new RankingState(inputs,layer);
		final int count=published.get();
		for(int i=0; i<count; i++){inputs.addTo(state.statistics,raw(i));}
		final MemberRankStatistics.Snapshot initial=state.statistics.freeze(layer,seed);
		state.members=initial.members; state.median=initial.medianLength;
		state.mean=initial.meanKd; state.layer=initial.layer;
		ranking=state;
	}

	boolean hasRanking(){return ranking!=null;}
	/** Identity, not matching array length, binds ordinal-derived statistics to their population. */
	boolean rankingUses(final MemberRankInputs inputs){
		checkFailure();
		final RankingState state=ranking;
		return inputs!=null && state!=null && state.inputs==inputs;
	}

	/** The coordinator selects the next measurement layer only after previous workers join. */
	void rankingLayer(final int layer){
		checkFailure();
		final RankingState state=requireRanking();
		if(layer<0 || claimed.get()!=published.get() || state.publishing.get()){
			throw new IllegalStateException("Measurement layer changes require quiescent membership and a nonnegative layer");
		}
		state.nextLayer=layer;
	}

	/** All observed fields are volatile: their ordering makes the version check safe under the Java memory model. */
	void readRanking(final MemberRankStatistics.View out){
		if(out==null){throw new IllegalArgumentException("Ranking reads need query-private reusable scratch");}
		final RankingState state=requireRanking();
		for(;;){
			checkFailure();
			final long before=state.version;
			if((before&1)!=0){continue;}
			final int members=state.members, layer=state.layer;
			final double median=state.median, mean=state.mean;
			if(before!=state.version){continue;}
			out.members=members; out.layer=layer; out.medianLength=median; out.meanKd=mean;
			return;
		}
	}

	/** Allocated only at a frozen phase/seed handoff; candidate ranking uses readRanking's reused view. */
	MemberRankStatistics.Snapshot rankingSnapshot(){
		final MemberRankStatistics.View view=new MemberRankStatistics.View();
		readRanking(view);
		return new MemberRankStatistics.Snapshot(view.members,view.medianLength,view.meanKd,view.layer);
	}

	/** Publish at most one native chunk per batch so a continuously growing family remains visible. */
	private void advanceRanked(final RankingState state){
		for(;;){
			if(!state.publishing.compareAndSet(false,true)){return;}
			try{
				checkFailure();
				final int before=published.get(), limit=Math.min(MAX_MEMBERS,before+CHUNK_SIZE);
				int next=before;
				while(next<limit && ready(next)){
					state.inputs.addTo(state.statistics,raw(next)); next++;
				}
				if(next!=before){
					final double median=state.statistics.medianLength(), mean=state.statistics.meanKd();
					assert(state.statistics.count()==next) : "Published descriptors cover exactly the contiguous ready member prefix";
					state.version++;
					state.members=next; state.median=median; state.mean=mean; state.layer=state.nextLayer;
					published.set(next);
					state.version++;
				}
			}catch(Throwable problem){
				failure.compareAndSet(null,problem);
				throw propagate(problem);
			}finally{state.publishing.set(false);}
			// Release before checking: an append that lost the CAS either appears here or becomes the next publisher.
			if(!ready(published.get())){return;}
		}
	}

	private RankingState requireRanking(){
		checkFailure();
		final RankingState state=ranking;
		if(state==null){throw new IllegalStateException("No retained-member ranking descriptors are attached");}
		return state;
	}

	private static final class RankingState{
		RankingState(final MemberRankInputs inputs_, final int layer_){inputs=inputs_; nextLayer=layer_;}
		final MemberRankInputs inputs;
		final MemberRankStatistics statistics=new MemberRankStatistics();
		final AtomicBoolean publishing=new AtomicBoolean();
		volatile long version;
		volatile int members,layer,nextLayer;
		volatile double median,mean;
	}

	private boolean ready(int slot){
		final AtomicIntegerArray[] table=chunks.get();
		final int chunk=slot>>>CHUNK_SHIFT;
		return chunk<table.length && table[chunk]!=null && table[chunk].get(slot&CHUNK_MASK)>0;
	}

	private void checkFailure(){
		final Throwable problem=failure.get();
		if(problem==null){return;}
		throw propagate(problem);
	}

	private static RuntimeException propagate(Throwable problem){
		if(problem instanceof Error){throw (Error)problem;}
		if(problem instanceof RuntimeException){throw (RuntimeException)problem;}
		return new IllegalStateException("Family member append failed", problem);
	}

	private static int tableLength(int capacity){
		final int needed=(capacity+CHUNK_SIZE-1)>>>CHUNK_SHIFT;
		int length=1;
		while(length<needed){
			if(length>=(1<<30)){throw new IllegalArgumentException("Requested family membership capacity is too large");}
			length<<=1;
		}
		return length;
	}

	private final AtomicReference<AtomicIntegerArray[]> chunks;
	private final AtomicInteger claimed=new AtomicInteger(), published=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private volatile RankingState ranking;
	private static final int CHUNK_SHIFT=8;
	private static final int CHUNK_SIZE=1<<CHUNK_SHIFT;
	private static final int CHUNK_MASK=CHUNK_SIZE-1;
	private static final int MAX_MEMBERS=1<<30;
}
