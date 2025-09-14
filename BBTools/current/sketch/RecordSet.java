package sketch;

import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import shared.Tools;

class RecordSet {

	RecordSet(int qID_){
		qID=qID_;
	}
	
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
	
	public void processSSU(){
		if(verbose){System.err.println("RecordSet.processSSU(): qID="+qID+", sorted="+sorted+", swept="+swept+", ssuProcessed="+ssuProcessed);}
		if(ssuProcessed || SSUMap.r16SMap==null){return;}//TODO
		for(Record r : records){
			r.processSSU();
		}
		ssuProcessed=true;
	}

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
	
	ArrayList<Record> records=new ArrayList<Record>(8);
	
	void addLevel(int level){
		long mask=1L<<level;
		assert((levels&mask)==0);
		levels|=mask;
	}
	
	boolean hasLevel(int level){
		long mask=1L<<level;
		return (levels&mask)==mask;
	}
	
	long levels;
	final int qID;
	
	boolean sorted=false;
	boolean swept=false;
	boolean ssuProcessed=false;
	
	static boolean verbose;
}