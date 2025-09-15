package align2;

import java.util.ArrayList;

import fileIO.TextFile;

public class ReformatBatchOutput2 {
	
//	Elapsed:	31.7
//
//	Mapping Statistics for 0s_default.sam:
//	mapped:                	100.00%
//	retained:              	96.06%
//	discarded:             	0.00%
//	ambiguous:             	3.94%
//
//	Strict correctness (both ends exactly correct):
//	true positive:         	96.06%
//	false positive:        	0.00%
//
//	Loose correctness (one end approximately correct):
//	true positive:         	96.06%
//	false positive:        	0.00%
//
//	false negative:        	0.00%
//	Elapsed:	2.34
//	Elapsed:	20.51
	
	
	public static void main(String[] args){
		TextFile tf=new TextFile(args[0], false);
		String[] lines=tf.toStringLines();
		ArrayList<String> list=new ArrayList<String>();
		
		int mode=0;
		
		System.out.println(header());
		
		for(String s : lines){
			if(s.startsWith("Elapsed:")){mode++;}
			if(mode>1){
				mode=0;
			}else{
//				list.add(s);
				if(s.startsWith("Mapping Statistics for ")){
					System.out.println(s.replace("Mapping Statistics for ", "").replace(".sam:", "")+"\t");
				}else if(s.startsWith("Mapping:")){
					s=s.replace("Mapping:", "").replace("seconds.", "").trim();
					System.out.print(s+"\t");
				}
			}
		}
	}
	
	
	public static String header() {
		return("name\tcount\ttime\tmapTime\tmapped\tretained\tdiscarded\tambiguous\ttruePositive\t" +
				"falsePositive\ttruePositiveL\tfalsePositiveL\tfalseNegative");
	}
	
}
