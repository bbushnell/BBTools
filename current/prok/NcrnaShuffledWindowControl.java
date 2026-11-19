package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import prok.ReplayNcrnaWindows.OriginalInput;
import prok.ReplayNcrnaWindows.Window;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Seed-bypassed null control on a declared reservoir of real saved candidates.
 * Candidate bases alone are permuted; original outside-window refinement flanks
 * and physical contig ends remain unchanged. No whole caller or path selection.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaShuffledWindowControl{

	static void run(final CallGenes caller, final ArrayList<OriginalInput> inputs,
			final String prefix, final int count, final long seed, final boolean cached, final int threads){
		final Window[] sample=new Window[count];
		final Random random=new Random(seed);
		int seen=0;
		for(final OriginalInput input:inputs){
			for(final Window window:input.windows){
				final int pick=seen<count ? seen : random.nextInt(seen+1);
				if(pick<count){sample[pick]=window;}
				seen++;
			}
		}
		if(count<1 || seen<count){throw new IllegalArgumentException("A shuffled sample needs the complete declared window universe");}
		final Item[] items=new Item[count];
		for(int i=0; i<count; i++){items[i]=new Item(sample[i]);}
		Arrays.sort(items, (a, b)->Long.compare(a.window.key, b.window.key));
		final int[] starts=new int[inputs.size()+1];
		for(final Item item:items){starts[item.window.input.slot+1]++;}
		for(int i=1; i<starts.length; i++){starts[i]+=starts[i-1];}
		// Source contexts belong to an input-file worker, so claim input files,
		// never individual windows that could race their original-source fallback.
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		for(int i=0; i<Math.min(threads, inputs.size()); i++){
			workers.add(new Worker(caller, inputs, items, starts, seed, cached, next, failure));
		}
		ThreadWaiter.startAndWait(workers);
		if(failure.get()!=null){throw new IllegalStateException("Shuffled verifier control failed", failure.get());}
		for(final Worker worker:workers){if(!worker.success){throw new IllegalStateException("Every shuffled-control worker must finish");}}
		write(prefix, items, inputs, seed, cached, seen);
	}

	private static final class Worker extends Thread implements NcrnaModelAttemptInstrumentSink{
		Worker(final CallGenes caller, final ArrayList<OriginalInput> inputs_, final Item[] items_, final int[] starts_,
				final long seed_, final boolean cached, final AtomicInteger next_, final AtomicReference<Throwable> failure_){
			inputs=inputs_;
			items=items_;
			starts=starts_;
			seed=seed_;
			next=next_;
			failure=failure_;
			byte[][] library=null;
			for(final NcrnaFamily family:GeneCaller.ncrnaFamilies){if("cobalamin_rf00174".equals(family.name)){library=family.library;}}
			if(library==null){throw new IllegalArgumentException("Shuffled control needs the configured RF00174 library");}
			cache=cached ? new NcrnaAlignmentCache(library, 1, Integer.MAX_VALUE) : null;
			for(int span=0; span<2; span++){
				for(int cut=0; cut<CUTS.length; cut++){verifiers[span][cut]=caller.shuffledCobalaminWindowWorker(CUTS[cut], span==1);}
			}
		}

		@Override
		public void run(){
			try{
				for(int slot=next.getAndIncrement(); slot<inputs.size(); slot=next.getAndIncrement()){
					if(failure.get()!=null){return;}
					for(int i=starts[slot]; i<starts[slot+1]; i++){measure(items[i]);}
				}
				success=true;
			}catch(RuntimeException | Error invalid){failure.compareAndSet(null, invalid);}
		}

		private void measure(final Item item){
			final Window window=item.window;
			active=item;
			final byte[] shuffled=window.context.bases.clone();
			final int from=window.start-window.context.origin;
			final int length=window.stop-window.start+1;
			random.setSeed(seed^window.key);
			for(int i=length-1; i>0; i--){
				final int j=random.nextInt(i+1);
				final byte base=shuffled[from+i];
				shuffled[from+i]=shuffled[from+j];
				shuffled[from+j]=base;
			}
			assert(sameCounts(window.context.bases, shuffled, from, length)) : "A shuffled null must preserve every candidate byte count and its length";
			final NcrnaScavenger.SavedContextSource original=window.input.contexts.sourceFor(window.context);
			final NcrnaScavenger.SavedContextSource flank=(context, origin, left, right)->{
				original.fill(window.context.bases, origin, left, right);
				for(int pos=left-origin; pos<=right-origin; pos++){
					if(context[pos]==0){context[pos]=window.context.bases[pos];}
				}
			};
			if(cache!=null){cache.clear();}
			final long before=cache==null ? 0 : cache.fills();
			final long hits=cache==null ? 0 : cache.hits();
			for(int span=0; span<2; span++){
				for(int cut=0; cut<CUTS.length; cut++){
					final NcrnaScavenger worker=verifiers[span][cut];
					final long work=worker.alignmentCount();
					final NcrnaDdlOrf accepted=worker.verifyShuffledWindow(window.name, shuffled, window.context.origin,
						window.length, window.strand, window.pass, window.start, window.stop, window.source, flank, this, cache);
					item.logical+=worker.alignmentCount()-work;
					if(accepted!=null){item.accepted|=1<<(span*CUTS.length+cut);}
				}
			}
			item.primaryFills=cache==null ? primaryAttempts : cache.fills()-before;
			item.physical=item.logical-(cache==null ? 0 : cache.hits()-hits);
			assert(item.physical>=item.primaryFills && item.physical<=item.logical)
				: "Only reused primary identities avoid kernels; actual refinement and edge work remain counted";
			primaryAttempts=0;
			item.done=true;
			active=null;
		}

		@Override
		public void alignmentSpan(final String name, final int strand, final int pass, final int start, final int stop,
				final int model, final float identity, final int left, final int right, final int stage,
				final boolean clip, final boolean edgeFloor){
			if(stage!=PRIMARY){return;}
			if(!Float.isFinite(identity)){throw new IllegalStateException("Actual primary alignment identity must be finite");}
			primaryAttempts++;
			active.best[0]=Math.max(active.best[0], identity);
			if(left>=0 && right-left+1>=109){active.best[1]=Math.max(active.best[1], identity);}
		}

		@Override
		public void modelAttempt(final String name, final int strand, final int pass, final int start, final int stop,
				final int rank, final int model, final int shared, final float identity, final int alignedLength, final boolean accepted){}

		private final ArrayList<OriginalInput> inputs;
		private final Item[] items;
		private final int[] starts;
		private final long seed;
		private final AtomicInteger next;
		private final AtomicReference<Throwable> failure;
		private final NcrnaScavenger[][] verifiers=new NcrnaScavenger[2][CUTS.length];
		private final NcrnaAlignmentCache cache;
		private final Random random=new Random();
		private Item active;
		private long primaryAttempts;
		boolean success;
	}

	private static boolean sameCounts(final byte[] original, final byte[] shuffled, final int from, final int length){
		final int[] counts=new int[256];
		for(int i=from; i<from+length; i++){
			counts[original[i]&255]++;
			counts[shuffled[i]&255]--;
		}
		for(final int count:counts){if(count!=0){return false;}}
		return true;
	}

	private static void write(final String prefix, final Item[] items,
			final ArrayList<OriginalInput> inputs, final long seed, final boolean cached, final int universe){
		final ByteBuilder row=new ByteBuilder();
		final long[] accepted=new long[12];
		long logical=0, physical=0, fills=0, sourceFiles=0, sourceBases=0;
		try(NcrnaTextIO.Output out=new NcrnaTextIO.Output(prefix+".shuffle.windows.tsv")){
			out.writer.println("input_id\tcaptured_event\tbest_primary_no_span\tbest_primary_span109\taccepted_setting_mask\tlogical_alignments\tphysical_primary_fills\tphysical_total_alignments");
			for(final Item item:items){
				if(!item.done){throw new IllegalStateException("Each selected shuffled window must finish before publication");}
				for(int setting=0; setting<accepted.length; setting++){if((item.accepted&(1<<setting))!=0){accepted[setting]++;}}
				logical+=item.logical;
				physical+=item.physical;
				fills+=item.primaryFills;
				out.writer.print(row.clear().append(item.window.input.name).tab().append(item.window.event).tab()
					.append(Float.toString(item.best[0])).tab().append(Float.toString(item.best[1])).tab()
					.append(item.accepted).tab().append(item.logical).tab().append(item.primaryFills).tab().append(item.physical).nl());
			}
		}
		for(final OriginalInput input:inputs){
			sourceFiles+=input.contexts.sourceLoads();
			sourceBases+=input.contexts.sourceBases();
		}
		try(NcrnaTextIO.Output out=new NcrnaTextIO.Output(prefix+".shuffle.summary.tsv")){
			out.writer.println("span_min\tidentity_cut\tsampled_windows\tprimary_distribution_windows\tmedian\tp95\tp99\tp999\tmaximum\taccepted_windows\taccepted_fraction");
			for(int span=0; span<2; span++){
				final float[] distribution=new float[items.length];
				int n=0;
				for(final Item item:items){if(item.best[span]>=0){distribution[n++]=item.best[span];}}
				Arrays.sort(distribution, 0, n);
				for(int cut=0; cut<CUTS.length; cut++){
					row.clear().append(span==0 ? 1 : 109).tab().append(Float.toString(CUTS[cut])).tab().append(items.length).tab().append(n);
					for(final int quantile:QUANTILES){
						final int rank=(int)(((long)n*quantile+999)/1000);
						row.tab().append(n==0 ? "NA" : Float.toString(distribution[rank-1]));
					}
					final long count=accepted[span*CUTS.length+cut];
					out.writer.print(row.tab().append(count).tab().append((double)count/items.length, 8).nl());
				}
			}
			out.writer.print(row.clear().append("# source=G11_ADDITION; shuffle=FISHER_YATES_CANDIDATE_ONLY; outside_flanks=ORIGINAL_REAL_BASES; seeds=BYPASSED; dinucleotide=UNAVAILABLE; seed=")
				.append(seed).append("; captured_universe=").append(universe).append("; exact_primary_cache=").append(cached)
				.append("; logical_alignments=").append(logical).append("; physical_primary_fills=").append(fills)
				.append("; physical_total_alignments=").append(physical)
				.append("; original_source_files=").append(sourceFiles).append("; original_source_bases=").append(sourceBases).nl());
		}
	}

	private static final class Item{
		Item(final Window window_){window=window_;}
		final Window window;
		final float[] best={-1,-1};
		int accepted;
		long logical, primaryFills, physical;
		boolean done;
	}
	private static final float[] CUTS={.60f,.62f,.64f,.66f,.68f,.70f};
	private static final int[] QUANTILES={500,950,990,999,1000};
}
