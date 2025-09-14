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
	
	public static double toSeconds(String s){
		s=s.replaceAll("s", "");
		String[] split=s.split("m");
		String seconds=split[1], minutes=split[0];
		return 60*Double.parseDouble(minutes)+Double.parseDouble(seconds);
	}
	
}
