package driver;

import fileIO.TextFile;
import shared.Tools;

/**
 * Utility for processing BBMerge comparison data and timing measurements.
 * Parses timing and accuracy statistics from BBMerge output files and formats
 * them into tab-separated values for analysis and comparison.
 *
 * @author Brian Bushnell
 * @date Feb 28, 2016
 */
public class ProcessSpeed {
	
	public static void main(String[] args){
		
		System.out.println("#real\tuser\tsys\tcorrect\tincorrect\tSNR");
		
		//NOTE [driver/ProcessSpeed#002] LOW/dev: args[0] used with no length guard → AIOOBE if invoked with no args.
		//.replace("in=","") strips the literal "in=" anywhere (fine for a filename; not a prefix-only strip). Crash-loud.
		String fname=args[0].replace("in=", "");
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("***")){
				//FIXED [driver/ProcessSpeed#001] was `line.replace("\\*\\*\\*","")`: String.replace is LITERAL, so it searched
				//for the 6-char sequence \*\*\* (backslash-star ×3), which never occurs in a line that startsWith("***").
				//=> the *** decoration was NEVER stripped (author confused replace with replaceAll). Now strips literal "***".
				System.out.println(line.replace("***", "").trim());
			}else if(line.startsWith("real\t")){
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
				System.out.print(Tools.format("%.3f\t", seconds));
			}else if(line.startsWith("Correct:")){
				System.out.print(line.split("\\p{javaWhitespace}+")[2]+"\t");
			}else if(line.startsWith("Incorrect:")){
				System.out.print(line.split("\\p{javaWhitespace}+")[2]+"\t");
			}else if(line.startsWith("SNR:")){
				System.out.print(line.split("\\p{javaWhitespace}+")[1]+"\n");
			}
//				Correct:                	99.72071%	15941011 reads
//				Incorrect:              	0.27929%	44646 reads
//				Too Short:              	0.02666%	4262 reads
//				Too Long:               	0.25263%	40384 reads
//				SNR:                    	25.539
			
			
			
		}
		
	}
	
	public static double toSeconds(String s){
		//Parses bash `time` format "MmSS.sss s" (e.g. "1m30.500s") → 60*min+sec. Duplicated in ProcessSpeed2.toSeconds.
		//Format-contract: input WITHOUT an 'm' (e.g. plain "30.5s") → split("m") length 1 → split[1] AIOOBE. Crash-loud, dev tool.
		s=s.replaceAll("s", "");
		String[] split=s.split("m");
		String seconds=split[1], minutes=split[0];
		return 60*Double.parseDouble(minutes)+Double.parseDouble(seconds);
	}
	
}
