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
		if(geneStartScoreCount>0){bb.append("Start Score:          \t ").append(average(geneStartScoreSum, geneStartScoreCount, "start"), 4).nl();}
		if(geneStopScoreCount>0){bb.append("Stop Score:           \t ").append(average(geneStopScoreSum, geneStopScoreCount, "stop"), 4).nl();}
		if(geneInnerScoreCount>0){bb.append("Inner Score:          \t ").append(average(geneInnerScoreSum, geneInnerScoreCount, "inner"), 4).nl();}
		if(lengthCount>0){bb.append("Length:               \t ").append(average(lengthSum, lengthCount, "length"), 2).nl();}
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
		if(geneStartScoreCount>0){jo.addLiteral("Start Score", average(geneStartScoreSum, geneStartScoreCount, "start"), 4);}
		if(geneStopScoreCount>0){jo.addLiteral("Stop Score", average(geneStopScoreSum, geneStopScoreCount, "stop"), 4);}
		if(geneInnerScoreCount>0){jo.addLiteral("Inner Score", average(geneInnerScoreSum, geneInnerScoreCount, "inner"), 4);}
		if(lengthCount>0){jo.addLiteral("Length", average(lengthSum, lengthCount, "length"), 2);}
		if(genomeSize>0){
			jo.addLiteral("Approx Genic Fraction", Tools.min(1.0, lengthSum/(double)genomeSize), 4);
		}
		return jo;
	}

	private static double average(double sum, long count, String label){
		if(count<=0 || !Double.isFinite(sum)){
			throw new IllegalStateException("ScoreTracker "+label+" average needs a positive count and finite sum: count="+count+" sum="+sum);
		}
		return sum/count;
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
