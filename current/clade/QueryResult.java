package clade;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import idaligner.IDAligner;
import tax.TaxTree;

/**
 * Centralizes the per-query search pipeline: clade search, sketch lookup,
 * SSU alignment, composite scoring, and display list assembly.
 *
 * The display list is the union of:
 *   - Top X records by kmer distance (k5dif or 2*k4dif)
 *   - Top X records by DDL sketch similarity (matches * wkid)
 *   - Top 1 record by SSU ANI
 * Deduplicated by TID, sorted by composite score.
 *
 * @author Chloe
 */
public class QueryResult {

	ArrayList<Comparison> displayList;
	int topCladeTID=-1;
	int topSketchTID=-1;
	int topSSUTID=-1;

	static QueryResult build(Clade query, CladeIndex index, int maxHits,
			int maxHitsToPrint, boolean showRecords, IDAligner ssa){

		QueryResult qr=new QueryResult();
		ArrayList<Comparison> all=index.findBest(query, maxHits);

		if(all==null || all.isEmpty()){
			qr.displayList=all;
			return qr;
		}

		// SSU align all entries
		if(Clade.callSSU && ssa!=null){
			for(Comparison comp : all){comp.align(ssa);}
		}

		// Separate clade hits from sketch-only hits; skip null refs
		ArrayList<Comparison> cladeHits=new ArrayList<Comparison>();
		ArrayList<Comparison> sketchHits=new ArrayList<Comparison>();
		for(Comparison c : all){
			if(c.ref==null) continue;
			if(c.isSketchHit) sketchHits.add(c);
			else cladeHits.add(c);
		}

		// Sort clade hits by effective k5dif (lower = better match)
		cladeHits.sort((a, b) -> Float.compare(effectiveKdif(a), effectiveKdif(b)));

		// Sort sketch hits by matches*wkid (higher = better match)
		sketchHits.sort((a, b) -> Float.compare(sketchScore(b), sketchScore(a)));

		// Collect kept entries: union of top X kmer, top X sketch, top 1 SSU
		LinkedHashMap<Integer, Comparison> kept=new LinkedHashMap<Integer, Comparison>();

		// Top clade hits by kmer
		int cladeCount=0;
		for(Comparison c : cladeHits){
			if(cladeCount>=maxHitsToPrint) break;
			if(!kept.containsKey(c.ref.taxID)){
				kept.put(c.ref.taxID, c);
			}
			cladeCount++;
		}
		Clade topCladeRef=cladeHits.isEmpty() ? null : cladeHits.get(0).ref;
		qr.topCladeTID=(topCladeRef!=null ? topCladeRef.taxID : -1);

		// Top sketch hits by matches*wkid
		int sketchCount=0;
		for(Comparison c : sketchHits){
			if(sketchCount>=CladeIndex.maxSketchHits) break;
			if(!kept.containsKey(c.ref.taxID)){
				kept.put(c.ref.taxID, c);
			}
			sketchCount++;
		}
		Clade topSketchRef=null;
		if(!sketchHits.isEmpty()){
			topSketchRef=sketchHits.get(0).ref;
			qr.topSketchTID=topSketchRef.taxID;
		}else if(topCladeRef!=null && cladeHits.get(0).sketchTaxID>0){
			qr.topSketchTID=cladeHits.get(0).sketchTaxID;
		}

		// Top 1 SSU hit
		Comparison bestSSU=null;
		for(Comparison c : all){
			if(c.ref!=null && c.ssudif<1f){
				if(bestSSU==null || c.ssudif<bestSSU.ssudif) bestSSU=c;
			}
		}
		if(bestSSU!=null){
			qr.topSSUTID=bestSSU.ref.taxID;
			if(!kept.containsKey(bestSSU.ref.taxID)){
				kept.put(bestSSU.ref.taxID, bestSSU);
			}
		}

		// Build display list sorted by composite score (higher = better)
		ArrayList<Comparison> display=new ArrayList<Comparison>(kept.values());
		final Clade fTopClade=topCladeRef;
		final Clade fTopSketch=topSketchRef;
		//CLEVER [verified in-file]: display list is a UNION of three rankings (top-X by kmer dist, top-X by DDL sketch matches*wkid, top-1 by SSU ANI), deduped by taxID via LinkedHashMap, then sorted by compositeScore. ALL sorts here use Float.compare (NaN-safe -- NaN sorts as max, consistently) unlike Comparison.compareTo's raw `<`, so a NaN dif won't break the sort contract.
		display.sort((a, b) -> {
			float sa=a.compositeScore(lcaFor(a, fTopClade, fTopSketch));
			float sb=b.compositeScore(lcaFor(b, fTopClade, fTopSketch));
			return Float.compare(sb, sa);
		});

		// Cache confidence
		if(showRecords){
			for(Comparison c : display){c.cacheConfidence();}
		}else if(!display.isEmpty()){
			display.get(0).cacheConfidence();
		}

		qr.displayList=display;
		return qr;
	}

	private static float effectiveKdif(Comparison c){
		return (c.k5dif<1f) ? c.k5dif : 2f*c.k4dif;
	}

	private static float sketchScore(Comparison c){
		return Math.max(c.wkid, 0f)*Math.max(c.kmerMatches, 0);
	}

	/** Compute LCA level using lineage strings. Always returns canonical levels. */
	private static int lcaFor(Comparison c, Clade topCladeRef, Clade topSketchRef){
		if(c.sketchLCA>0) return c.sketchLCA;

		Clade other=(c.isSketchHit ? topCladeRef : topSketchRef);
		if(c.ref!=null && other!=null){
			CharSequence refLin=c.ref.lineage();
			CharSequence otherLin=other.lineage();
			if(refLin!=null && otherLin!=null){
				int level=CladeIndex.lineageLCA(refLin, otherLin);
				if(level>0) return level;
			}
		}
		return TaxTree.LIFE;
	}
}
