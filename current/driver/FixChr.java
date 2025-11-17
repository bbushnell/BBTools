package driver;

import fileIO.TextFile;
import fileIO.TextStreamWriter;

/** One-off program for converting grch38 sam files to hg19 */
public class FixChr {
	
	/**
	 * Program entry point for chromosome name conversion.
	 * Reads input SAM file line by line and adds "chr" prefix to chromosome names.
	 * Processes both data lines and contig header lines appropriately.
	 * @param args Command-line arguments: [0] input file path, [1] output file path
	 */
	public static void main(String[] args){
		
		String in=args[0];
		String out=args[1];
		
		TextFile tf=new TextFile(in);
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, true);
		tsw.start();
		
		String s=null;
		while((s=tf.nextLine())!=null){
			if(!s.startsWith("#")){s="chr"+s;}
			else if(s.startsWith("##contig=<ID=")){
				s="##contig=<ID=chr"+s.substring("##contig=<ID=".length());
			}
			tsw.println(s);
		}
		tf.close();
		tsw.poisonAndWait();
	}
	
}
