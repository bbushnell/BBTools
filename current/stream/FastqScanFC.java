package stream;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

import shared.Timer;
import shared.Tools;
import simd.Vector;
import structures.IntList;

/**
 * Retained FileChannel experiment counting four-line FASTQ-like raw files.
 * Groups every four LF positions into a record and counts the bytes between
 * the first two LF positions as sequence bytes. Does not validate headers,
 * separators or quality lengths, decode compression, or support wrapped records.
 *
 * Historical note attributed to Brian: this approach was slower than the
 * InputStream scanner on HPC and was deliberately kept as a negative result.
 * The current repository search finds only this class's diagnostic main as a
 * caller; do not remove the experiment or infer a new performance result here.
 * Existing limitations #001–#003 remain below. The earlier note referenced
 * bug_reports/stream/FastqScanFC.md, which is not included in this clone.
 * @author Collei
 */
public class FastqScanFC{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Scans the first argument as a raw filename and prints timing/counts to stderr.
	 * @param args At least one element; args[0] is the filename
	 * @throws RuntimeException Wraps an IOException from the scan
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		String fname=args[0];
		// We don't need full FileFormat logic for this raw scan, just the name
		FastqScanFC fqs=new FastqScanFC(fname);
		try{fqs.scan();}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop();
		String s=Tools.timeReadsBasesProcessed(t, fqs.totalRecords, fqs.totalBases, 12);
		System.err.println(s);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Stores the filename without opening it.
	 * @param fname_ Raw file to open when scan is called */
	public FastqScanFC(String fname_){
		fname=fname_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds this file's four-LF groups and sequence-byte counts to the existing totals.
	 * Uses a fixed 262144-byte buffer, shifting residual bytes after complete groups.
	 * Totals are not reset. The known limitations below remain unfixed in this
	 * retained experiment; this method does not validate record structure.
	 * @throws IOException Propagates file opening, reading or closing errors
	 */
	void scan() throws IOException{
		//[stream/FastqScanFC#003 LOW latent] Resources close only on the normal path,
		//not in finally; a scan that exits exceptionally can leave them open. Retained dormant limitation.
		@SuppressWarnings("resource")
		RandomAccessFile raf=new RandomAccessFile(fname, "r");
		FileChannel channel=raf.getChannel();
		
		// Fixed 256 KiB experimental buffer; retained without a new tuning claim.
		final int bufSize=262144; 
		byte[] buffer=new byte[bufSize];
		ByteBuffer bb=ByteBuffer.wrap(buffer);
		
		IntList newlines=new IntList(4096);
		int bstop=0, residue=0, bstart=0;
		
		// FileChannel read loop
		while(true){
			// Read into the buffer, respecting the residue (bytes moved to front)
			bb.position(residue);
			bb.limit(buffer.length);
			int r=channel.read(bb);
			
			if(r<=0 && residue==0){break;} // EOF and no residue
			
			if(r<0){r=0;}// EOF
			bstop=residue+r;
			
			// Scan for newlines
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines);
			
			int records=newlines.size/4;
			totalRecords+=records;
			
			// Process records
			for(int i=0, j=0; i<records; i++, j+=4){
				// int headerEnd=newlines.get(j);
				// int basesEnd=newlines.get(j+1); // Not needed for simple counting if just skipping
				// int plusEnd=newlines.get(j+2);
				int recordEnd=newlines.get(j+3);
				
				// Calculate bases length if needed for stats
				// bases = basesEnd - headerEnd - 1
				//[stream/FastqScanFC#002 LOW latent] No trailing-CR strip: the LF distance
				//includes CR, counting one extra sequence byte per CRLF record. Retained limitation.
				int bases=newlines.get(j+1)-newlines.get(j)-1;
				totalBases+=bases;
				
				bstart=recordEnd+1;
			}
			
			residue=bstop-bstart;
			if(residue>0){
				// Shift residue to beginning of array
				System.arraycopy(buffer, bstart, buffer, 0, residue);
			}
			bstart=0;
			newlines.clear();
			
			if(r==0 && residue>0){
				// No new bytes and residue remains: stop without counting the partial group.
				// A full buffer can also cause a zero-byte read, as described by #001.
				//[stream/FastqScanFC#001 LOW latent] No expansion: a record exceeding the buffer
				//can fill the residue without completing a group. A zero-byte read then takes
				//this break, dropping that record and the remaining file. Retained limitation.
				break;
			}
		}
		
		channel.close();
		raf.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Raw filename reopened by each scan call. */
	private final String fname;
	/** Cumulative complete four-LF groups counted across scan calls. */
	long totalRecords;
	/** Cumulative sequence-byte distances; includes trailing CR when present. */
	long totalBases;
}
