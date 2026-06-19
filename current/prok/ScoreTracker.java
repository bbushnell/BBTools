package prok;

import java.util.ArrayList;

import json.JsonObject;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Accumulates scoring statistics for Open Reading Frames (ORFs) of a specific type.
 * Tracks gene start scores, stop scores, inner k-mer scores, and sequence lengths
 * to calculate average scores and genome coverage estimates for prokaryotic gene prediction.
 * @author Brian Bushnell
 */
public class ScoreTracker {
	
	/** Creates a ScoreTracker for ORFs of the specified type.
	 * @param type_ The ORF type to track (e.g., gene, start codon, stop codon) */
	public ScoreTracker(int type_){
		type=type_;
	}
	
	/**
	 * Merges statistics from another ScoreTracker into this one.
	 * Adds all accumulated sums and counts from the source tracker.
	 * @param st The ScoreTracker to merge into this one
	 */
	public void add(ScoreTracker st){
		geneStartScoreSum+=st.geneStartScoreSum;
		geneStopScoreSum+=st.geneStopScoreSum;
		geneInnerScoreSum+=st.geneInnerScoreSum;
		lengthSum+=st.lengthSum;
		
		geneStartScoreCount+=st.geneStartScoreCount;
		geneStopScoreCount+=st.geneStopScoreCount;
		geneInnerScoreCount+=st.geneInnerScoreCount;
		lengthCount+=st.lengthCount;
	}
	
	/**
	 * Adds statistics from an array of ORF lists.
	 * Processes each list in the array by calling add(ArrayList<Orf>).
	 * @param array Array of ORF lists to process
	 */
	public void add(ArrayList<Orf>[] array){
		for(ArrayList<Orf> list : array){add(list);}
	}
	
	/**
	 * Adds statistics from all ORFs in the list that match this tracker's type.
	 * Only ORFs with the same type as this tracker contribute to the statistics.
	 * @param list List of ORFs to process
	 */
	public void add(ArrayList<Orf> list){
		if(list==null){return;}
		for(Orf orf : list){
			//only matching-type ORFs contribute; add(Orf) re-checks type+null. A null list *element* would NPE at orf.type (out of scope: malformed list).
			if(orf.type==type){add(orf);}
		}
	}
	
	/**
	 * Adds statistics from a single ORF if it matches this tracker's type.
	 * Accumulates start score, stop score, average k-mer score, and length.
	 * @param orf The ORF to add (null or wrong type ORFs are ignored)
	 */
	public void add(Orf orf){
		if(orf==null || orf.type!=type){return;}
		geneStartScoreSum+=orf.startScore;
		geneStopScoreSum+=orf.stopScore;
		geneInnerScoreSum+=orf.averageKmerScore();
		lengthSum+=orf.length();
		
		//lockstep with the sums above: each count rises with its sum, so count = #contributing ORFs (count==0 only when sum==0)
		geneStartScoreCount++;
		geneStopScoreCount++;
		geneInnerScoreCount++;
		lengthCount++;
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		bb.append("Start Score:          \t ").append(geneStartScoreSum/geneStartScoreCount, 4).nl();
		bb.append("Stop Score:           \t ").append(geneStopScoreSum/geneStopScoreCount, 4).nl();
		bb.append("Inner Score:          \t ").append(geneInnerScoreSum/geneInnerScoreCount, 4).nl();
		bb.append("Length:               \t ").append(lengthSum/(double)lengthCount, 2).nl();
		if(genomeSize>0){
			bb.append("Approx Genic Fraction:\t ").append(Tools.min(1.0, lengthSum/(double)genomeSize), 4).nl();
		}
		return bb.toString();
	}
	
	/**
	 * Converts the accumulated statistics to a JSON object.
	 * Includes average scores and genome coverage fraction if genome size is set.
	 * @return JSON object containing all computed statistics
	 */
	public JsonObject toJson(){
		JsonObject jo=new JsonObject();
		//TODO: Possible bug [prok/ScoreTracker#001] (MEDIUM; reachable) - sum/count is non-finite when count==0 (empty tracker: 0.0/0=NaN) or when geneInnerScoreSum absorbed NaN/Inf from Orf.averageKmerScore() on a short ORF (Orf.java:85 divides by length()-kInnerCDS-2). addLiteral formats NaN via String.format -> a bare "NaN" token = INVALID JSON (RFC 8259). (toString instead prints "0.0000": ByteBuilder casts (long)NaN=0 -- silently wrong but valid text.) Serialized at CallGenes.printStatsJson (L540-573, gated by callX) -> reachable when an enabled type finds zero features + stats=json. Fix is a design call (Brian): guard count==0 here, or guard non-finite in json/JsonLiteral.
		jo.addLiteral("Start Score", geneStartScoreSum/geneStartScoreCount, 4);
		jo.addLiteral("Stop Score", geneStopScoreSum/geneStopScoreCount, 4);
		jo.addLiteral("Inner Score", geneInnerScoreSum/geneInnerScoreCount, 4);
		jo.addLiteral("Length", lengthSum/(double)lengthCount, 2);
		if(genomeSize>0){
			jo.addLiteral("Approx Genic Fraction", Tools.min(1.0, lengthSum/(double)genomeSize), 4);
		}
		return jo;
	}
	
	/** Number of gene start scores accumulated */
	long geneStartScoreCount=0;
	/** Number of gene stop scores accumulated */
	long geneStopScoreCount=0;
	/** Number of gene inner k-mer scores accumulated */
	long geneInnerScoreCount=0;
	/** Number of ORF lengths accumulated */
	long lengthCount=0;
	
	/** Sum of all gene start scores */
	double geneStartScoreSum=0;
	/** Sum of all gene stop scores */
	double geneStopScoreSum=0;
	/** Sum of all gene inner k-mer scores */
	double geneInnerScoreSum=0;
	/** Sum of all ORF lengths in bases */
	long lengthSum=0;
	
	/** Total genome size for calculating genic fraction coverage */
	long genomeSize=0;
	
	/** ORF type that this tracker monitors */
	final int type;
	
}
