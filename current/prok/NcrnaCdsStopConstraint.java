package prok;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import structures.ByteBuilder;
import structures.IntList;

/** Worker-owned CDS start positions for an opt-in ncRNA endpoint experiment.
 * Uses scored/pruned CDS candidates before RNA generation or gene-path selection.
 * Coordinates in each list are zero-based and oriented to that CDS strand.
 * Existing RNA/CDS overlap is preserved; only additional stop extension is limited.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCdsStopConstraint {

	/** Shared experiment owner; creates one private coordinate index per scavenger.
	 * The optional writer is borrowed and closed by the experiment driver. */
	public static final class Factory {
		public Factory(final ByteStreamWriter writer_){
			writer=writer_;
			if(writer!=null){writer.println(HEADER);}
		}

		/** Called during worker initialization, never per base or per seed. */
		public synchronized NcrnaCdsStopConstraint create(){
			final NcrnaCdsStopConstraint worker=new NcrnaCdsStopConstraint(this);
			workers.add(worker);
			return worker;
		}

		/** Aggregate only after CallGenes has joined every calling worker. */
		public synchronized long[] totals(){
			final long[] totals=new long[3];
			for(NcrnaCdsStopConstraint worker : workers){
				totals[0]+=worker.sweeps;
				totals[1]+=worker.limitedSweeps;
				totals[2]+=worker.excludedSites;
			}
			assert(totals[1]<=totals[0]) : "Limited sweeps are a subset of observed NN sweeps";
			return totals;
		}

		/** Rare accepted-locus output; serialize native byte-writer handoff only. */
		private synchronized void write(final ByteBuilder row){
			assert(writer!=null) : "Only an explicitly requested constraint log receives rows";
			writer.print(row);
		}

		private final ByteStreamWriter writer;
		private final ArrayList<NcrnaCdsStopConstraint> workers=new ArrayList<NcrnaCdsStopConstraint>();
	}

	public NcrnaCdsStopConstraint(){this(null);}

	private NcrnaCdsStopConstraint(final Factory owner_){owner=owner_;}

	/** Replaces this contig's coordinates without modifying borrowed CDS objects.
	 * Null candidates mean no CDS context, as in an isolated makeRnas invocation. */
	public void prepare(final ArrayList<Orf>[] candidates, final int length){
		assert(length>=0) : "GeneCaller supplies the current complete contig length";
		contigLength=length;
		starts[0].clear();
		starts[1].clear();
		if(candidates==null){return;}
		for(ArrayList<Orf> list : candidates){
			if(list==null){continue;}
			for(Orf cds : list){
				assert(cds!=null && cds.type==ProkObject.CDS && (cds.strand==0 || cds.strand==1)
					&& cds.flipped()==cds.strand && cds.scaflen==length
					&& cds.start>=0 && cds.stop>=cds.start && cds.stop<length)
					: KillSwitch.assertDie("GeneCaller.breakOrfs flips minus candidates once into genomic coordinates before RNA merging; length="+length);
				final int start=cds.strand==0 ? cds.start : length-1-cds.stop;
				starts[cds.strand].add(start);
			}
		}
		starts[0].sort();
		starts[1].sort();
	}

	/** Earliest downstream CDS start in the RNA's sense frame, with raw overlap kept.
	 * The result is frozen before either endpoint moves. MAX_VALUE means no ceiling. */
	public int ceiling(final int strand, final int start, final int stop, final int length){
		assert((strand==0 || strand==1) && length==contigLength && start>=0
			&& stop>=start && stop<length)
			: KillSwitch.assertDie("Constraint lookup must use the prepared contig and unchanged sense-oriented NN anchors");
		final IntList list=starts[strand];
		int lo=0, hi=list.size;
		while(lo<hi){
			final int mid=(lo+hi)>>>1;
			if(list.array[mid]<=start){lo=mid+1;}
			else{hi=mid;}
		}
		return lo==list.size ? Integer.MAX_VALUE : Math.max(stop, list.array[lo]-1);
	}

	/** Records candidate-sweep events, not an inferred number of unique final loci. */
	public void record(final Orf orf, final int model, final int start, final int stop,
			final int ceiling, final NcrnaBoundaryScorer.RefinementResult result,
			final int appliedStartOffset, final int appliedStopOffset){
		assert(result!=null && result.stopSitesExcluded>=0 && ceiling>=stop)
			: KillSwitch.assertDie("Observed exclusions must belong to a completed legal endpoint sweep");
		sweeps++;
		if(result.stopSitesExcluded>0){limitedSweeps++;}
		excludedSites+=result.stopSitesExcluded;
		if(owner==null || owner.writer==null){return;}
		bb.clear();
		for(int i=0; i<orf.scafName.length(); i++){
			final char c=orf.scafName.charAt(i);
			bb.append(c=='\t' || c=='\n' || c=='\r' ? ' ' : c);
		}
		bb.tab().append(orf.strand).tab().append(contigLength).tab().append(model)
			.tab().append(start).tab().append(stop).tab().append(ceiling)
			.tab().append(result.startOffset).tab().append(result.stopOffset)
			.tab().append(appliedStartOffset).tab().append(appliedStopOffset)
			.tab().append(result.startSitesScored).tab().append(result.stopSitesScored)
			.tab().append(result.stopSitesExcluded).nl();
		owner.write(bb);
	}

	private final Factory owner;
	private final IntList[] starts={new IntList(), new IntList()};
	private final ByteBuilder bb=new ByteBuilder();
	private int contigLength=-1;
	private long sweeps, limitedSweeps, excludedSites;
	private static final String HEADER="contig\tstrand\tcontig_length\tmodel\tanchor_start0\tanchor_stop0"
		+"\tstop_ceiling0\tsearch_start_offset\tsearch_stop_offset\tapplied_start_offset\tapplied_stop_offset"
		+"\tstart_sites\tstop_sites\tstop_sites_excluded";
}
