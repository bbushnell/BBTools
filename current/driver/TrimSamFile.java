package driver;

import java.util.HashSet;

import fileIO.ByteFile;
import shared.LineParser1;
import shared.Tools;
import stream.SamLine;

/**
 * Utility for filtering SAM files by removing reads that overlap with specified genomic regions.
 * Identifies reads that map to or span a target scaffold region and excludes them from output.
 * Useful for removing contaminated or problematic genomic regions from alignment data.
 * @author Brian Bushnell
 */
public class TrimSamFile {
	
	/**
	 * Program entry point for SAM file trimming utility.
	 * Expects 4 arguments: filename, scaffold name, start position, and end position.
	 * @param args Command-line arguments: [filename] [scaffold] [start] [end]
	 */
	public static void main(String[] args){
		String fname=args[0];
		String scaf=args[1];
		int from=Integer.parseInt(args[2]);
		int to=Integer.parseInt(args[3]);
		ByteFile tf=ByteFile.makeByteFile(fname, false);
		HashSet<String> set=findBadLines(tf, scaf, from, to);
		tf.reset();
		printExcludingSet(tf, set);
	}
	
	
	/**
	 * Identifies read names that should be excluded from the SAM file.
	 * Finds reads that overlap the specified genomic region or have mapping problems.
	 * Includes reads where either mate maps within the target region, reads that span
	 * the region boundary, and unpaired or unmapped reads.
	 *
	 * @param tf Input SAM file to scan
	 * @param scafS Target scaffold/chromosome name
	 * @param from Start position of region to exclude (inclusive)
	 * @param to End position of region to exclude (inclusive)
	 * @return Set of read names to exclude from output
	 */
	public static HashSet<String> findBadLines(ByteFile tf, String scafS, int from, int to){
		byte[] scaf=scafS.getBytes();
		HashSet<String> set=new HashSet<String>(16000);
		LineParser1 lp=new LineParser1('\t');
		
		for(byte[] s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s[0]!='@'){//header
				SamLine sl=new SamLine(lp.set(s));
				
				if(sl.pos>=from && sl.pos<=to && Tools.equals(sl.rname(), scaf)){
					set.add(sl.qname);
				}else if(sl.pnext>=from && sl.pnext<=to && Tools.equals(sl.rnext(), scaf)){
					set.add(sl.qname);
				}else if(Tools.equals(sl.rname(), scaf) && Tools.equals(sl.rnext(), scaf) && (sl.pos<from != sl.pnext<from)){
					set.add(sl.qname);
				}else if(!sl.mapped() || !sl.nextMapped() || !sl.pairedOnSameChrom()){
					set.add(sl.qname);
				}
			}
		}
		return set;
	}
	
	
	/**
	 * Prints SAM file content while excluding reads with names in the provided set.
	 * Preserves all header lines and only filters alignment records.
	 * @param tf Input SAM file to process
	 * @param set Set of read names to exclude from output
	 */
	public static void printExcludingSet(ByteFile tf, HashSet<String> set){
		LineParser1 lp=new LineParser1('\t');
		for(byte[] s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s[0]=='@'){//header
				System.out.println(s);
			}else{
				SamLine sl=new SamLine(lp.set(s));
				
				if(!set.contains(sl.qname)){
					System.out.println(s);
				}
			}
		}
	}
	
	
}

