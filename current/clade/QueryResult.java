package clade;

import java.util.ArrayList;
import java.util.Collections;

import idaligner.IDAligner;

/**
 * Centralizes the per-query search pipeline: clade search, sketch lookup,
 * SSU alignment, sorting, and trimming. Produces the display list and
 * records top hits for future composite scoring.
 *
 * @author Chloe
 */
public class QueryResult {

	ArrayList<Comparison> displayList;
	int topCladeTID=-1;
	int topSketchTID=-1;
	int topSSUTID=-1;

	/**
	 * Build a QueryResult using the same pipeline as ProcessThread.run().
	 * Output is identical to the inline code it replaces.
	 */
	static QueryResult build(Clade query, CladeIndex index, int maxHits,
			int maxHitsToPrint, boolean showRecords, IDAligner ssa){

		QueryResult qr=new QueryResult();
		ArrayList<Comparison> list=index.findBest(query, maxHits);

		if(list!=null){
			if(Clade.callSSU && ssa!=null){
				for(Comparison comp : list){comp.align(ssa);}
			}
			Collections.sort(list);
			if(list.size()>maxHitsToPrint){
				ArrayList<Comparison> refined=new ArrayList<Comparison>();
				int cladeCount=0, sketchCount=0;
				for(Comparison c : list){
					if(c.isSketchHit){
						if(sketchCount<CladeIndex.maxSketchHits){refined.add(c); sketchCount++;}
					}else if(cladeCount<maxHitsToPrint){
						refined.add(c); cladeCount++;
					}
				}
				list.clear();
				list.addAll(refined);
			}
			if(showRecords){
				for(Comparison c : list){c.cacheConfidence();}
			}else if(!list.isEmpty()){
				list.get(0).cacheConfidence();
			}

			// Record top TIDs for future composite scoring
			for(Comparison c : list){
				if(!c.isSketchHit && c.ref!=null && qr.topCladeTID<0){
					qr.topCladeTID=c.ref.taxID;
				}
			}
			if(!list.isEmpty() && list.get(0).sketchTaxID>0){
				qr.topSketchTID=list.get(0).sketchTaxID;
			}
			float bestSSUDif=1f;
			for(Comparison c : list){
				if(c.ref!=null && c.ssudif<bestSSUDif){
					bestSSUDif=c.ssudif;
					qr.topSSUTID=c.ref.taxID;
				}
			}
		}

		qr.displayList=list;
		return qr;
	}
}
