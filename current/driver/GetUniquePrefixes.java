package driver;

import java.util.HashSet;

import fileIO.TextFile;

/**
 * Extracts unique sequence prefixes from FASTA files.
 * Reads a FASTA file and outputs only sequences with unique prefixes of specified length.
 * Useful for generating non-redundant sequence datasets based on initial subsequences.
 * @author Brian Bushnell
 */
public class GetUniquePrefixes {
	
	/**
	 * Program entry point for extracting unique sequence prefixes.
	 * Processes a FASTA file and outputs sequences with unique prefixes.
	 * @param args Command-line arguments: [0] input FASTA filename, [1] prefix length
	 */
	public static void main(String[] args){
		
		
		String fname=args[0];
		int prefix=Integer.parseInt(args[1]);
		TextFile tf=new TextFile(fname);
		
		HashSet<String> set=new HashSet<String>();
		
		String header=null;
		StringBuilder sequence=new StringBuilder();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith(">")){
				if(sequence.length()>0){
					if(sequence.length()>prefix){sequence.setLength(prefix);}
					String s=sequence.toString();
					if(set.contains(s)){}
					else{
						set.add(s);
						System.out.println(header+"\n"+s);
					}
				}
				sequence.setLength(0);
				header=line;
			}else{
				sequence.append(line);
			}
		}
		//TODO: Possible bug [driver/GetUniquePrefixes#001] LOW/dev (no .sh, no callers): the LAST FASTA record is never
		//emitted — dedup+print only fires when the NEXT ">" header is seen (L32-41), so after EOF the final accumulated
		//`sequence`/`header` is dropped. Missing final flush (same shape as GenerateNoCallsFromCoverage#001). Also args[0]/
		//parseInt(args[1]) unguarded (AIOOBE/NFE) and tf is never closed (leak). Dead one-off → LOW; the dropped-last-record
		//is the real one. Fix: replicate the L32-41 flush block once more after the loop.
	}
	
}
