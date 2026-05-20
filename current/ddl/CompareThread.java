package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import idaligner.QuantumAligner;

/**
 * Worker thread for parallel DDL comparison.  Distributes queries via
 * AtomicLong work-stealing, fills per-query result lists, and optionally
 * performs inline SSU alignment with a per-thread QuantumAligner.
 *
 * Supports buffered alignment: heap collects top `buffer` candidates by
 * composite score, aligns SSU on those, re-sorts by ANI, then trims to
 * `maxRecords`.  This bounds alignment cost while ensuring the best-by-ANI
 * hit is likely captured.
 *
 * @author Noire, Brian Bushnell
 * @date May 19, 2026
 */
class CompareThread extends Thread {
	final ArrayList<DDLRecord> queries, refs;
	final AtomicLong nextQuery;
	final ArrayList<DDLComparison>[] results;
	final int nQueries, nRefs, k, maxRecords, buffer, minHits;
	final boolean useIndex, alignSSU, banSelf;
	final DDLIndex idx;
	final AtomicLong totalComparisonsPerformed;
	long alignCountT;

	CompareThread(ArrayList<DDLRecord> queries, ArrayList<DDLRecord> refs,
			AtomicLong nextQuery, ArrayList<DDLComparison>[] results,
			int nQueries, int nRefs, int k, int maxRecords, int minHits,
			boolean useIndex, DDLIndex idx, AtomicLong totalComparisonsPerformed,
			boolean alignSSU){
		this(queries, refs, nextQuery, results, nQueries, nRefs, k, maxRecords, 0,
			minHits, useIndex, idx, totalComparisonsPerformed, alignSSU, false);
	}

	CompareThread(ArrayList<DDLRecord> queries, ArrayList<DDLRecord> refs,
			AtomicLong nextQuery, ArrayList<DDLComparison>[] results,
			int nQueries, int nRefs, int k, int maxRecords, int buffer, int minHits,
			boolean useIndex, DDLIndex idx, AtomicLong totalComparisonsPerformed,
			boolean alignSSU, boolean banSelf){
		this.queries=queries; this.refs=refs;
		this.nextQuery=nextQuery; this.results=results;
		this.nQueries=nQueries; this.nRefs=nRefs;
		this.k=k; this.maxRecords=maxRecords;
		this.buffer=Math.max(buffer, 20+2*maxRecords);
		this.minHits=minHits;
		this.useIndex=useIndex; this.idx=idx;
		this.totalComparisonsPerformed=totalComparisonsPerformed;
		this.alignSSU=alignSSU;
		this.banSelf=banSelf;
	}

	@Override
	public void run(){
		long localComparisons=0;
		long localAligns=0;
		DDLComparison working=new DDLComparison();
		QuantumAligner aligner=(alignSSU ? new QuantumAligner() : null);
		while(true){
			int qi=(int)nextQuery.getAndIncrement();
			if(qi>=nQueries){break;}
			DDLComparisonHeap heap=new DDLComparisonHeap(buffer);
			DDLRecord qRec=queries.get(qi);
			final int qTid=qRec.taxID;
			if(useIndex){
				int[] counts=idx.query(qRec.ddl);
				for(int ri=0; ri<nRefs; ri++){
					if(counts[ri]<minHits){continue;}
					if(banSelf && qTid>0 && refs.get(ri).taxID==qTid){continue;}
					working.compare(qRec, refs.get(ri), k);
					heap.offer(working);
					localComparisons++;
				}
				if(heap.size()<buffer && minHits>1){
					for(int ri=0; ri<nRefs; ri++){
						if(counts[ri]>=minHits || counts[ri]<1){continue;}
						if(banSelf && qTid>0 && refs.get(ri).taxID==qTid){continue;}
						working.compare(qRec, refs.get(ri), k);
						heap.offer(working);
						localComparisons++;
					}
				}
			}else{
				for(int ri=0; ri<nRefs; ri++){
					if(banSelf && qTid>0 && refs.get(ri).taxID==qTid){continue;}
					working.compare(qRec, refs.get(ri), k);
					heap.offer(working);
					localComparisons++;
				}
			}
			ArrayList<DDLComparison> list=heap.toList();
			if(aligner!=null){
				for(DDLComparison c : list){
					if(c.queryRecord==null || c.refRecord==null){continue;}
					int qt=c.queryRecord.riboType(), rt=c.refRecord.riboType();
					if(qt!=DDLRecord.RIBO_NONE && qt==rt){
						c.ssuIdentity=aligner.align(c.queryRecord.riboBytes(), c.refRecord.riboBytes());
						localAligns++;
					}
				}
				list.sort((a, b) -> {
					float sa=1000f*(a.ssuIdentity>=0 ? a.ssuIdentity : 0)+a.cachedScore;
					float sb=1000f*(b.ssuIdentity>=0 ? b.ssuIdentity : 0)+b.cachedScore;
					return sa>sb ? -1 : sa<sb ? 1 : 0;
				});
			}
			if(list.size()>maxRecords){
				list=new ArrayList<>(list.subList(0, maxRecords));
			}
			results[qi]=list;
		}
		totalComparisonsPerformed.addAndGet(localComparisons);
		alignCountT=localAligns;
	}
}
