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
