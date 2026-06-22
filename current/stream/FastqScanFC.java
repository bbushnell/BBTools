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
 * FastqScan using FileChannel for IO.
 * Uses a large buffer and FileChannel.read() to bypass InputStream synchronization overhead.
 * NOTE (Brian): KEPT as a deliberate negative result — this FileChannel approach proved SLOWER than the
 * InputStream path in HPC, so it is intentionally not wired in (no caller). Its existence is the record
 * that FileChannels are not worthwhile here; do not delete it as "dead code." Its latent divergences from
 * FastqScan (no expand(), no trailing-\r strip, no try/finally close) are documented in
 * bug_reports/stream/FastqScanFC.md and left unfixed since it is not a live path.
 * @author Collei
 */
public class FastqScanFC {

	public static void main(String[] args) {
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
	
	public FastqScanFC(String fname_) {
		fname=fname_;
	}
	
	void scan() throws IOException {
		//[stream/FastqScanFC#003 LOW latent] raf/channel are closed on the happy path (L99-100) but NOT in a try/finally → leaked if scan() throws mid-loop (the @SuppressWarnings acknowledges it). Orphan class (no caller) → latent.
		@SuppressWarnings("resource")
		RandomAccessFile raf=new RandomAccessFile(fname, "r");
		FileChannel channel=raf.getChannel();
		
		// 64k buffer is small for FileChannel. Let's go bigger. 
		// 256k or 1MB is often the sweet spot for modern SSD/OS page caching.
		final int bufSize=262144; 
		byte[] buffer=new byte[bufSize];
		ByteBuffer bb=ByteBuffer.wrap(buffer);
		
		IntList newlines=new IntList(4096);
		int bstop=0, residue=0, bstart=0;
		
		// FileChannel read loop
		while(true) {
			// Read into the buffer, respecting the residue (bytes moved to front)
			bb.position(residue);
			bb.limit(buffer.length);
			int r=channel.read(bb);
			
			if(r<=0 && residue==0) {break;} // EOF and no residue
			
			if(r<0) r=0; // EOF
			bstop=residue+r;
			
			// Scan for newlines
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines);
			
			int records=newlines.size/4;
			totalRecords+=records;
			
			// Process records
			for(int i=0, j=0; i<records; i++, j+=4) {
				// int headerEnd=newlines.get(j);
				// int basesEnd=newlines.get(j+1); // Not needed for simple counting if just skipping
				// int plusEnd=newlines.get(j+2);
				int recordEnd=newlines.get(j+3);
				
				// Calculate bases length if needed for stats
				// bases = basesEnd - headerEnd - 1
				//[stream/FastqScanFC#002 LOW latent] omits the trailing-\r strip that sibling FastqScan does (slashr1) → +1 base per record on \r\n (Windows) FASTQ. Orphan → latent.
				int bases=newlines.get(j+1) - newlines.get(j) - 1;
				totalBases+=bases;
				
				bstart=recordEnd+1;
			}
			
			residue=bstop-bstart;
			if(residue>0) {
				// Shift residue to beginning of array
				System.arraycopy(buffer, bstart, buffer, 0, residue);
			}
			bstart=0;
			newlines.clear();
			
			if(r==0 && residue>0) {
				// We hit EOF but have a partial record left.
				// This usually means a truncated file or a file not ending in newline.
				// In a scanner, we might just drop it or count it as partial.
				//[stream/FastqScanFC#001 LOW latent] Doubles as the no-expand failure mode: a record > bufSize (256KB) never completes → records=0 → residue fills the buffer → next read returns 0 (no room) → this break DROPS the giant record AND the rest of the file → silent undercount on long-read (>256KB) FASTQ. Sibling FastqScan grows via expand(); this orphan variant has none.
				break;
			}
		}
		
		channel.close();
		raf.close();
	}
	
	private final String fname;
	long totalRecords;
	long totalBases;
}