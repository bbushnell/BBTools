package sketch;

import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import shared.Tools;

/**
 * Collection of taxonomic classification records for a single query sequence.
 * Provides sorting, deduplication, SSU processing, and accuracy testing.
 * @author Brian Bushnell
 * @date January 1, 2020
 */
class RecordSet {

	/** Creates a RecordSet for the specified query ID.
	 * @param qID_ Query sequence identifier */
	RecordSet(int qID_){
		qID=qID_;
	}
	
	/**
	 * Sorts records by ANI and removes duplicate taxonomic levels.
	 * Uses bit masking to track seen taxonomic levels and removes duplicates.
	 * Records are sorted in descending order of ANI (best matches first).
	 */
	public void sortAndSweep() {
		if(verbose){System.err.println("RecordSet.sortAndSweep(): qID="+qID+", sorted="+sorted+", swept="+swept+", ssuProcessed="+ssuProcessed);}
		if(sorted && swept){return;}
		
		Collections.sort(records);
		sorted=true;
		
		long seen=0;
		int removed=0;
		for(int i=0; i<records.size(); i++){
			Record r=records.get(i);
			long mask=1L<<r.taxLevelExtended;
			if((seen&mask)!=0){
				records.set(i, null);
				removed++;
			}
			seen|=mask;
		}
		if(removed>0){
			Tools.condenseStrict(records);
		}
		swept=true;
	}
	
	/**
	 * Processes 16S SSU data for all records if available.
	 * Calls processSSU() on each record to compute SSU alignments from the SSUMap.
	 * Only processes once and skips if SSUMap is not initialized.
	 */
	public void processSSU(){
		if(verbose){System.err.println("RecordSet.processSSU(): qID="+qID+", sorted="+sorted+", swept="+swept+", ssuProcessed="+ssuProcessed);}
		if(ssuProcessed || SSUMap.r16SMap==null){return;}//TODO
		for(Record r : records){
			r.processSSU();
		}
		ssuProcessed=true;
	}

	/**
	 * Tests taxonomic classification accuracy at different taxonomic levels.
	 * Compares predicted taxonomy against expected taxonomy using both
	 * NCBI taxonomy and SSU data when available.
	 *
	 * @param bswBad Writer for failed classification details
	 * @return Status array indexed by taxonomic level with classification results
	 */
	int[] test(ByteStreamWriter bswBad){
		if(verbose){System.err.println("RecordSet.test(): qID="+qID+", sorted="+sorted+", swept="+swept+", ssuProcessed="+ssuProcessed);}
		boolean failed=false;
		int[] status=new int[AnalyzeSketchResults.taxLevels];
		for(int level=1; level<AnalyzeSketchResults.taxLevels; level++){
			Record first=null;
			boolean correctTax=true;
			boolean correctSSU=true;
			boolean missingSSU=false;
			for(Record r : records){
				if(r.taxLevelExtended<level){
					//Ignore
				}else{
					if(first==null){
						first=r;
						correctTax=true;
						correctSSU=true;
						missingSSU=first.ssu()<=0;
					}else{
						if(r.taxLevelExtended>=first.taxLevelExtended){
							//OK
						}else{
							//Incorrect NCBI
							correctTax=false;
							if(r.ssu()<=0 || first.ssu()<=0){
								//Missing SSU
								missingSSU=true;
							}else if(first.ssu()>=r.ssu()){
								//Correct SSU
							}else{
								//Incorrect SSU
								correctSSU=false;
								failed=true;
							}
						}
					}
				}
			}
			
			int x;
			if(first==null){
				//No hit
				x=AnalyzeSketchResults.NOHIT;
			}else if(correctTax){
				x=AnalyzeSketchResults.CORRECT;
			}else{
				x=AnalyzeSketchResults.INCORRECT_TAX;
				if(!correctSSU){x|=AnalyzeSketchResults.INCORRECT_SSU;}
				else if(missingSSU){x|=AnalyzeSketchResults.MISSING_SSU;}
			}
			status[level]=x;
			if(first==null){break;}//array is initialized to zero anyway
		}
		
		if(failed && bswBad!=null){
			bswBad.println();
			for(Record r : records){
				bswBad.println(r.text);
			}
		}
		
		return status;
	}
	
	/** Collection of taxonomic classification records for this query */
	ArrayList<Record> records=new ArrayList<Record>(8);
	
	/** Marks a taxonomic level as present in this record set.
	 * @param level Taxonomic level to add */
	void addLevel(int level){
		long mask=1L<<level;
		assert((levels&mask)==0);
		levels|=mask;
	}
	
	/**
	 * Checks if a taxonomic level is present in this record set.
	 * @param level Taxonomic level to check
	 * @return True if level is present
	 */
	boolean hasLevel(int level){
		long mask=1L<<level;
		return (levels&mask)==mask;
	}
	
	/** Bit mask tracking which taxonomic levels are present */
	long levels;
	/** Query sequence identifier */
	final int qID;
	
	/** Whether records have been sorted by ANI */
	boolean sorted=false;
	/** Whether duplicate taxonomic levels have been removed */
	boolean swept=false;
	/** Whether SSU processing has been completed */
	boolean ssuProcessed=false;
	
	/** Enable verbose debugging output */
	static boolean verbose;
}