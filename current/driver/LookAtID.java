package driver;

import java.util.Arrays;

import fileIO.TextFile;
import shared.Tools;
import stream.SiteScoreR;

/**
 * @author Brian Bushnell
 * @date Dec 3, 2012
 *
 */
public class LookAtID {
	
	/**
	 * Program entry point for analyzing SiteScoreR ID values.
	 * Reads a tab-separated text file and examines each SiteScoreR object
	 * for numeric ID overflow conditions. Reports the maximum ID found
	 * and details of any IDs that exceed Integer.MAX_VALUE.
	 *
	 * @param args Command-line arguments where args[0] is the input file path
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
