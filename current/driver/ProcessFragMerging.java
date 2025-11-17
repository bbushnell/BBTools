package driver;

import fileIO.TextFile;
import shared.Tools;

/**
 * For BBMerge comparison data collation
 * @author Brian Bushnell
 * @date Mar 15, 2016
 *
 */
public class ProcessFragMerging {
	
	/**
	 * Program entry point that processes a BBMerge log file.
	 * Parses the input file line by line, extracting specific statistics
	 * and reformatting them as tab-separated output for data analysis.
	 * @param args Command-line arguments where args[0] is the input filename
	 */
	public static void main(String[] args){
		
		String sym="\t";
		
		String fname=args[0];
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			String[] split=line.split("\\p{javaWhitespace}+");
			if(line.startsWith("***")){
				System.out.print("\n"+split[1]+sym);
//				System.out.println("\n"+line);
			}else if(line.startsWith("real")){
				String time=line.split("\t")[1];
				double seconds=toSeconds(time);
				System.out.print(Tools.format("%.3f", seconds)+sym);
			}else if(line.startsWith("Reads Used:")){
				System.out.print(split[2]+sym+split[3].substring(1)+sym);
			}else if(line.startsWith("mapped:")){
				System.out.print(split[2]+sym+split[4]+sym);
			}else if(line.startsWith("Error Rate:")){
				System.out.print(split[3]+sym+split[5]+sym);
			}else if(line.startsWith("Sub Rate:")){
				System.out.print(split[3]+sym+split[5]+sym);
			}else if(line.startsWith("Del Rate:")){
				System.out.print(split[3]+sym+split[5]+sym);
			}else if(line.startsWith("Ins Rate:")){
				System.out.print(split[3]+sym+split[5]+sym);
			}
//				Del Rate:        	  0.0161% 	      168 	  0.0276% 	       64385
//				Ins Rate:        	  0.0017% 	       18 	  0.0002% 	         366
			
		}
		
	}
	
	/**
	 * Converts time string in "Xm Ys" format to total seconds.
	 * Parses minutes and seconds from time output and returns combined value.
	 * @param s Time string in format like "1m 23.456s"
	 * @return Total time in seconds as a double
	 */
	public static double toSeconds(String s){
		s=s.replaceAll("s", "");
		String[] split=s.split("m");
		String seconds=split[1], minutes=split[0];
		return 60*Double.parseDouble(minutes)+Double.parseDouble(seconds);
	}
	
}
