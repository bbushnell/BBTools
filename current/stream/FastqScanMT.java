package stream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.bam.BgzfSettings;
import structures.IntList;

/**
 * Counts reads and bases in a sequence file, with low overhead and multithreading.
 * Uses a "stateless" synchronization approach where threads grab chunks,
 * find the FASTQ frame alignment (@...+) locally, and count based on that.
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 7, 2025
 */
public final class FastqScanMT {

	public static void main(String[] args){
		Timer t=new Timer();
		if(args.length<1) {throw new RuntimeException("Usage: fastqscan.sh filename");}
		String fname=args[0];
		while(fname.startsWith("-")) {fname=fname.substring(1);}
		if(fname.startsWith("in=")) {fname=fname.substring(3);}
		int threads=Math.min(2, Shared.threads());
		BgzfSettings.READ_THREADS=Tools.mid(1, 18, Shared.threads());
		for(int i=1; i<args.length; i++) {
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("t") || a.equals("threads")) {threads=Integer.parseInt(b);}
			else if(a.equalsIgnoreCase("simd")) {Shared.SIMD&=Parse.parseBoolean(b);}
			else if(Tools.isNumeric(arg)) {threads=Integer.parseInt(arg);}
			else if(Parser.parseCommonStatic(arg, a, b)) {}
			else if(Parser.parseZip(arg, a, b)) {}
			else {assert(false) : "Unknown parameter "+arg;}
		}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		if(ff.stdin()){
			//Do nothing
		}else{
			File f=new File(fname);
			if(!f.isFile() || !f.canRead()){
				throw new RuntimeException("Can't read "+fname);
			}
		}
		FastqScanMT fqs=new FastqScanMT(ff);
		try{fqs.read(threads);}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:   \t");
		System.err.println("Records:\t"+fqs.totalRecords);
		System.err.println("Bases:  \t"+fqs.totalBases);
	}

	public static long[] countReadsAndBases(String fname, boolean halveInterleaved, int readThreads, int zipThreads) {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		return countReadsAndBases(ff, halveInterleaved, readThreads, zipThreads);
	}

	/** Returns molecules, reads, bases, file headers */
	public static long[] countReadsAndBases(FileFormat ff, boolean halveInterleaved, int readThreads, int zipThreads) {
		final int oldZT=BgzfSettings.READ_THREADS;
		if(ff.compressed()) {
			BgzfSettings.READ_THREADS=(zipThreads>1 ? zipThreads : Tools.mid(1, Shared.threads(), 18));
		}
		int recordsPerRead=1;
		if(ff.fastq() && halveInterleaved) {
			int[] iq=FileFormat.testInterleavedAndQuality(ff.name(), false);
			recordsPerRead=(iq[1]==FileFormat.INTERLEAVED ? 2 : 1);
		}
		FastqScanMT fqs=new FastqScanMT(ff);
		try{fqs.read(readThreads);}
		catch(IOException e){
			e.printStackTrace();
			//throw new RuntimeException(e);
			return null;
		}finally {BgzfSettings.READ_THREADS=oldZT;}
		long[] ret=new long[] {fqs.totalRecords/recordsPerRead, fqs.totalRecords, 
			fqs.totalBases, fqs.totalRecords};
		return ret;
	}

	FastqScanMT(FileFormat ff_){ff=ff_;}

	void read(int threads) throws IOException {
		if(!ff.fastq()){throw new RuntimeException("FastqScanMT only supports FASTQ.");}
		
		final InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		threads=(threads<1 ? Math.min(2, Shared.threads()) : threads);//Peaks at 2
		final ArrayList<ScanThread> alst=new ArrayList<ScanThread>(threads);
		
		for(int i=0; i<threads; i++){
			ScanThread st=new ScanThread(is);
			alst.add(st);
			st.start();
		}
		
		boolean success=true;
		for(ScanThread st : alst){
			while(st.getState()!=Thread.State.TERMINATED){
				try {st.join();}
				catch(InterruptedException e){e.printStackTrace();}
			}
			synchronized(st){
				success&=st.success;
				totalRecords+=st.recordsT;
				totalBases+=st.basesT;
			}
		}
		
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
		if(!success){throw new RuntimeException("Scanning failed.");}
	}
	
	private class ScanThread extends Thread {
		
		ScanThread(InputStream is_){
			is=is_;
		}
		
		@Override
		public void run(){
			try {
				process();
				success=true;
			}catch (IOException e){
				e.printStackTrace();
			}
		}
		
		private void process() throws IOException {
			while(true){
				final int len=fillBuffer();
				if(len<1){break;}
				scanBuffer(len);
			}
		}
		
		/** * Fills buffer synchronized from IS. 
		 * Does not return until buffer is full or stream is empty. 
		 */
		private int fillBuffer() throws IOException {
			synchronized(is){
				int len=0;
				while(len<buffer.length){
					int r=is.read(buffer, len, buffer.length-len);
					if(r<0){break;}
					len+=r;
				}
				return len;
			}
		}
		
		private void scanBuffer(final int len){
			if(len<1){return;}
			
			// 1. Find newlines
			newlines.clear();
			Vector.findSymbols(buffer, 0, len, (byte)'\n', newlines);
			
			// Handle missing terminal newline on the very last block
			if(len<buffer.length && (len==0 || buffer[len-1]!='\n')){
				newlines.add(len); // Pretend there is a newline at the very end
			}
			
			final int lines=newlines.size;
			if(lines==0){
				// Special case: Huge block with no newlines (single sequence line?)
				// We can't identify it, so we assume it's a Sequence line if we can't prove otherwise.
				// But 1MB without newlines is weird. Assuming sequence base count.
				assert(false) : "Record exceeded buffer length";
				basesT+=len;
				return;
			}
			
			// 2. Determine Frame (Self-Stabilization)
			int frameStart=-1; // Index in newlines of the first confirmed Header line
			
			// Scan for @...+.  i is the index of the PREVIOUS newline.
			// Line i starts at newlines.get(i-1)+1.
			// We iterate through lines to find @Header (Frame 0) and +Plus (Frame 2)
			for(int i=0; i<lines-2; i++){
				final int startHeader=(i==0 ? 0 : newlines.get(i-1)+1);
				if(buffer[startHeader]=='@'){
					final int startPlus=newlines.get(i+1)+1;
					if(buffer[startPlus]=='+'){
						frameStart=i;
						break;
					}
				}
			}
			
			// 4a. Edge case: No markers found (End of file, or weird small buffer)
			// Assume standard FASTQ structure relative to end: Last line is Qual (Frame 3)
			if(frameStart<0){
				// If we are at EOF (len < buffer.length), the last line is a Qual line.
				// If we are NOT at EOF, this is a very weird buffer (all sequence?), but
				// with 1MB buffers this shouldn't happen for valid FASTQ.
				// We will assume the last line is Frame 3.
				frameStart=(lines-1)-3; 
				// frameStart might be negative, which is fine for the modulo math below.
			}
			
			// Helper: Frame 0=Head, 1=Seq, 2=Plus, 3=Qual
			// We calculate frame relative to frameStart (which is Frame 0)
			
			int lineStart=0;
			for(int i=0; i<lines; i++){
				final int lineEnd=newlines.get(i);
				// Distance from known header line
				final int dist=i-frameStart;
				// Modulo 4, handling negatives
				final int frame=(dist & 3); 
				
				if(frame==1){ // Sequence Line
					int lineLen=lineEnd-lineStart; // Exclude newline
					if(lineLen>0 && buffer[lineEnd-1]=='\r'){lineLen--;}
					basesT+=lineLen;
				}else if(frame==0){ // Header Line
					recordsT++;
				}
				
				lineStart=lineEnd+1;
			}
			
			// Handle residue (bytes after last newline) if any
			if(lineStart < len){
				// The residue is the start of the NEXT line.
				// Current line i=lines.
				final int dist=lines-frameStart;
				final int frame=(dist & 3);
				
				if(frame==1){ // Partial Sequence Line
					int lineLen=len-lineStart;
					basesT+=lineLen;
				}else if(frame==0){ 
					// Partial Header Line. 
					// Do NOT count record yet; it will be counted by the next thread 
					// which sees the newline terminating this header.
				}
			}
		}
		
		private final InputStream is;
		private final byte[] buffer=new byte[1024*1024]; // 1MB buffer
		private final IntList newlines=new IntList(1024*16);
		
		long recordsT=0;
		long basesT=0;
		boolean success=false;
	}

	private final FileFormat ff;
	long totalRecords=0;
	long totalBases=0;

}