package hmm;

import java.util.HashMap;

/**
 * Tracks protein sequence hit lengths across HMM search results.
 * Maintains a mapping of reference model names to their maximum hit lengths,
 * allowing tracking of best protein sequence matches for a query.
 * @author Brian Bushnell
 */
public class ProteinSummary {
	
	/** Creates a protein summary for tracking hits against a named query sequence.
	 * @param name Query sequence identifier */
	public ProteinSummary(String name){
		this.name=name;
	}
	
	/** Returns true if anything changed */
	public boolean add(HMMSearchLine line){
		Integer old=map.get(line.name);
		if(old==null || old<line.length){
			map.put(line.name, line.length);
			return true;
		}
		return false;
	}
	
	/** Name of query sequence */
	public final String name;
	
	//This could alternatively be a map of names to lines, to include all data
	/** Map of reference model name to hit length */
	public HashMap<String, Integer> map=new HashMap<String, Integer>();
	
}
