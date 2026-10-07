package stream.bam;

import java.io.FileInputStream;
import java.io.IOException;

/**
 * Reports returned byte counts while consuming input through BgzfInputStreamMT.
 * Output uses the historical label "Block" for each nonnegative read result;
 * these counts describe read calls, not necessarily compressed BGZF blocks.
 * Does not compare byte contents against expected data.
 */
public class DebugBgzfRead{
	/**
	 * Reads with an 8192-byte buffer and prints each count, EOF and the final total.
	 * Passes one as the reader's thread argument; counts zero-byte results as calls.
	 * The file stream and reader are owned by the try-with-resources statement.
	 * @param args Optional input path, defaulting to test_mt_output.bam; later arguments ignored
	 * @throws IOException If opening, reading or closing the input fails
	 */
	public static void main(String[] args) throws IOException{
		String file=args.length>0 ? args[0] : "test_mt_output.bam";

		try(FileInputStream fis=new FileInputStream(file);
			BgzfInputStreamMT in=new BgzfInputStreamMT(fis, 1)){

			byte[] buf=new byte[8192];
			long total=0;
			int blockNum=0;

			while(true){
				int n=in.read(buf);
				if(n<0){
					System.out.println("Block "+blockNum+": EOF");
					break;
				}
				System.out.println("Block "+blockNum+": "+n+" bytes (total: "+(total+n)+")");
				total+=n;
				blockNum++;

				// Historical inspection conditions; neither branch changes reading.
				if(blockNum>=374 && blockNum<=380){
					// No action for this inspection range.
				}else if(blockNum>374){
					// No action outside the inspection range.
				}
			}

			System.out.println("\nTotal: "+total+" bytes in "+blockNum+" blocks");
		}
	}
}
