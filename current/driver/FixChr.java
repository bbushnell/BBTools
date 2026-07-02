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
			//NOTE [driver/FixChr#001] LOW/dev DOC-vs-CODE MISMATCH (no .sh, no callers, "one-off" per javadoc): the javadoc
			//claims SAM files ("adds chr prefix to chromosome names in SAM headers and data lines"), but this logic is
			//VCF-shaped, not SAM: (a) it keys header lines on '#'/'##contig=<ID=' — SAM header lines start with '@' (@HD/@SQ/@PG),
			//not '#'; (b) for a non-'#' line it prepends "chr" to the WHOLE line, which for a real SAM data line lands on QNAME
			//(col 1), NOT the RNAME chromosome (col 3). So on actual SAM it corrupts @-headers ("chr@HD...") and mis-prefixes
			//read names. It only does the right thing on a chrom-FIRST format (e.g. a stripped TSV/BED-like file). Either the
			//javadoc is mislabeled or the tool is broken-for-SAM. args[0]/args[1] also unguarded (AIOOBE). Dead dev toy → LOW.
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
