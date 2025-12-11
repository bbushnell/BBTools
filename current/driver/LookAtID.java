package driver;

import java.util.Arrays;

import fileIO.TextFile;
import shared.Tools;
import stream.SiteScoreR;

/**
 * Utility for analyzing numeric IDs in SiteScoreR text files.
 * Reads tab-separated files containing SiteScoreR objects and checks for
 * numeric ID overflow conditions (IDs >= Integer.MAX_VALUE).
 * Reports the maximum ID found and any overflow instances.
 *
 * @author Brian Bushnell
 * @date Dec 3, 2012
 */
public class LookAtID {
	
	/**
	 * Program entry point for analyzing SiteScoreR ID values in a tab-separated file.
	 * Reads each line, converts entries to SiteScoreR objects, tracks the maximum numeric ID, and reports any IDs that overflow Integer.MAX_VALUE.
	 * @param args args[0] is the input file path
	 */
	public static void main(String[] args){
		
		TextFile tf=new TextFile(args[0], true);
		
		long max=0;
		
		long line=0;
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			SiteScoreR[] array=SiteScoreR.fromTextArray(s);
			String[] split=s.split("\t");
			for(int i=0; i<array.length; i++){
				SiteScoreR ssr=array[i];
				String s2=split[i];
				max=Tools.max(ssr.numericID, max);
				if(ssr.numericID>=Integer.MAX_VALUE){
					System.out.println("Found overflow ID "+ssr.numericID+" at line "+line);
					System.out.println("ssr="+ssr.toText());
					System.out.println("raw="+s2);
					System.out.println("All:\n"+Arrays.toString(split));
					System.out.println();
					break;
				}
			}
			line++;
		}
		tf.close();
		System.out.println("Max ID was "+max);
		
	}
	
}
