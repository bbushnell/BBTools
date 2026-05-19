package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import idaligner.QuantumAligner;

/**
 * Worker thread for parallel DDL comparison.  Distributes queries via
 * AtomicLong work-stealing, fills per-query result lists, and optionally
 * performs inline SSU alignment with a per-thread QuantumAligner.
 *
 * @author Noire, Brian Bushnell
 * @date May 19, 2026
 */
class CompareThread extends Thread {
	final ArrayList<DDLRecord> queries, refs;
	final AtomicLong nextQuery;
	final ArrayList<DDLComparison>[] results;
	final int nQueries, nRefs, k, maxRecords, minHits;
	final boolean useIndex, alignSSU;
	final DDLIndex idx;
	final AtomicLong totalComparisonsPerformed;
	long alignCountT;

	CompareThread(ArrayList<DDLRecord> queries, ArrayList<DDLRecord> refs,
			AtomicLong nextQuery, ArrayList<DDLComparison>[] results,
			int nQueries, int nRefs, int k, int maxRecords, int minHits,
			boolean useIndex, DDLIndex idx, AtomicLong totalComparisonsPerformed,
			boolean alignSSU){
		this.queries=queries; this.refs=refs;
		this.nextQuery=nextQuery; this.results=results;
		this.nQueries=nQueries; this.nRefs=nRefs;
		this.k=k; this.maxRecords=maxRecords; this.minHits=minHits;
		this.useIndex=useIndex; this.idx=idx;
		this.totalComparisonsPerformed=totalComparisonsPerformed;
		this.alignSSU=alignSSU;
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
			DDLComparisonHeap heap=new DDLComparisonHeap(maxRecords);
			DDLRecord qRec=queries.get(qi);
			if(useIndex){
				int[] counts=idx.query(qRec.ddl);
				for(int ri=0; ri<nRefs; ri++){
					if(counts[ri]<minHits){continue;}
					working.compare(qRec, refs.get(ri), k);
					heap.offer(working);
					localComparisons++;
				}
			}else{
				for(int ri=0; ri<nRefs; ri++){
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
			}
			results[qi]=list;
		}
		totalComparisonsPerformed.addAndGet(localComparisons);
		alignCountT=localAligns;
	}
}
