package driver;

import fileIO.TextFile;
import fileIO.TextStreamWriter;

/**
 * One-off utility program for converting GRCH38 SAM files to hg19 format.
 * Adds "chr" prefix to chromosome names in SAM headers and data lines.
 * Specifically handles contig headers and non-comment lines in SAM files.
 * @author Brian Bushnell
 */
public class FixChr {
	
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
