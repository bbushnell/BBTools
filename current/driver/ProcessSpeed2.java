package driver;

import fileIO.TextFile;
import shared.Tools;

/**
 * For generic data collation
 * @author Brian Bushnell
 * @date December 6, 2016
 *
 */
public class ProcessSpeed2 {
	
	/**
	 * Parses timing data file and outputs converted times in seconds.
	 * Reads a file specified in args[0] (with optional "in=" prefix) and processes
	 * lines starting with "real", "user", or "sys" followed by time measurements.
	 * Outputs three columns: real, user, sys times converted to decimal seconds.
	 *
	 * @param args Command-line arguments where args[0] is the input file path
	 */
	public static void main(String[] args){
		
		System.out.println("#real\tuser\tsys");
		
		String fname=args[0].replace("in=", "");
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("real\t")){
				String time=line.split("\t")[1];
				double seconds=toSeconds(time);
				System.out.print(Tools.format("%.3f\t", seconds));
			}else if(line.startsWith("user\t")){
				String time=line.split("\t")[1];
				double seconds=toSeconds(time);
				System.out.print(Tools.format("%.3f\t", seconds));
			}else if(line.startsWith("sys\t")){
				String time=line.split("\t")[1];
				double seconds=toSeconds(time);
				System.out.print(Tools.format("%.3f\n", seconds));
			}
			
		}
		
	}
	
	/**
	 * Converts time string in minutes:seconds format to total seconds.
	 * Parses strings like "2m15.5s" and converts to decimal seconds (135.5).
	 * Removes 's' suffix and splits on 'm' to extract minutes and seconds.
	 *
	 * @param s Time string in format "Nm##.##s" (e.g., "2m15.5s")
	 * @return Total time in seconds as a double
	 */
	public static double toSeconds(String s){
		s=s.replaceAll("s", "");
		String[] split=s.split("m");
		String seconds=split[1], minutes=split[0];
		return 60*Double.parseDouble(minutes)+Double.parseDouble(seconds);
	}
	
}
