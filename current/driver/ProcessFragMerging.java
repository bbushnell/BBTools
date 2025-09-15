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
	
	public static double toSeconds(String s){
		s=s.replaceAll("s", "");
		String[] split=s.split("m");
		String seconds=split[1], minutes=split[0];
		return 60*Double.parseDouble(minutes)+Double.parseDouble(seconds);
	}
	
}
