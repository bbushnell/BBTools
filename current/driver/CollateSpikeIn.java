package driver;

import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date Oct 10, 2013
 *
 */
public class CollateSpikeIn {
	
	/**
	 * Parses BBMap log files and outputs tabulated alignment statistics.
	 * Reads log files line by line looking for mapping percentages, accuracy rates,
	 * and job identifiers, then outputs formatted results to stdout.
	 * Expected output format: "jobID\t%Control (BBMap)\t%Accuracy (BBMap)"
	 *
	 * @param args Command-line arguments where args[0] is the log file path to parse
	 */
	public static void main(String[] args){
//		Executing align2.BBMapPacBio [minratio=0.40, fastareadlen=500, out=null, in=/projectb/shared/pacbio/jobs/026/026437/data/filtered_subreads.fasta]

		System.out.println("jobID\t%Control (BBMap)\t%Accuracy (BBMap)");
		
		TextFile tf=new TextFile(args[0], false);
		String file=null, mapped=null, acc=null;
		String line=tf.nextLine();
		while(line!=null){
			if(line.startsWith("mapped:")){
				String[] split=line.split("\\p{javaWhitespace}+");
				mapped=split[1].replace("%", "");
			}else if(line.startsWith("Match Rate:")){
				String[] split=line.split("\\p{javaWhitespace}+");
				acc=split[2].replace("%", "");
				System.out.println(file+"\t"+mapped+"\t"+acc);
				file=acc=mapped=null;
			}else if(line.startsWith("Executing align2.BBMap")){
				String[] split=line.split("\\p{javaWhitespace}+");
				for(String s : split){
					if(s.startsWith("in=")){
						file=s.replace("in=", "").replace("]", "").replace(",", "");
						file=file.replace("/projectb/shared/pacbio/jobs/", "").replace("/data/filtered_subreads.fasta", "");
						file=file.substring(file.indexOf('/')+1);
						mapped=acc=null;
						break;
					}
				}
			}
			line=tf.nextLine();
		}
		
	}
	
}
