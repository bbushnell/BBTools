package stream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.bam.BgzfSettings;

/**
 * Counts literal LF bytes and bytes delivered by a ReadWrite input stream.
 * An unterminated final line contributes bytes but no LF; CR alone adds no line.
 * Decoding depends on the selected input stream, so the byte count need not equal
 * the compressed file size. This is a byte scanner, not a record-format validator.
 * Scan workers share input reads under a lock and count their private chunks outside it.
 * The static count API temporarily changes global BGZF settings; callers must coordinate
 * such configuration changes rather than assuming independent per-call settings.
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 15, 2025
 */
public final class FileScanMT{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Scans one filename (optionally prefixed by in=) and prints time, LF count and byte count.
	 * Optional t/threads controls scan workers; it does not update Shared's thread setting.
	 * Initializes global BGZF read threads before parsing options and does not restore them.
	 * @param args Filename followed by optional scan-thread, SIMD and common/zip settings
	 */
	public static void main(String[] args){
		Timer t=new Timer(System.out);
		if(args.length<1){throw new RuntimeException("Usage: filescan.sh filename");}
		String fname=args[0];
		while(fname.startsWith("-")){fname=fname.substring(1);}
		if(fname.startsWith("in=")){fname=fname.substring(3);}
		int threads=1;
		BgzfSettings.READ_THREADS=Tools.mid(1, 18, Shared.threads());
		for(int i=1; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equalsIgnoreCase("simd")){Shared.SIMD&=Parse.parseBoolean(b);}
			else if(Tools.isNumeric(arg)){threads=Integer.parseInt(arg);}
			else if(Parser.parseCommonStatic(arg, a, b)){}
			else if(Parser.parseZip(arg, a, b)){}
			else{assert(false) : "Unknown parameter "+arg;}
		}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		if(ff.stdin()){
			//Do nothing
		}else{
			File f=new File(fname);
			if(!f.isFile() || !f.canRead()){throw new RuntimeException("Can't read "+fname);}
		}
		FileScanMT fqs=new FileScanMT(ff);
		try{fqs.read(threads);}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:   \t");
		System.out.println("Lines:  \t"+fqs.totalLines);
		System.out.println("Bytes:  \t"+fqs.totalBytes);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an input descriptor and delegates to the descriptor-based count method.
	 * @param fname Input name accepted by FileFormat and ReadWrite
	 * @param readThreads Positive scan-worker count, or below one for automatic selection
	 * @param zipThreads Requested BGZF worker setting; see the descriptor overload
	 * @return New [LF count, byte count] array, or null for a caught IOException
	 */
	public static long[] countLinesAndBytes(String fname, int readThreads, int zipThreads){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		return countLinesAndBytes(ff, readThreads, zipThreads);
	}
	
	/** Counts decoded stream bytes and LF occurrences, restoring the prior global BGZF setting.
	 * A compressed descriptor sets the BGZF worker count to zipThreads when greater
	 * than one, or Shared's thread count bounded to 1 through 18 otherwise. An uncompressed
	 * descriptor does not override the setting before scanning.
	 * The finally block restores the setting but does not isolate overlapping calls.
	 * IOException from read() is printed and converted to null; unchecked failures propagate.
	 * @param ff Nonnull input descriptor; read() opens by its filename
	 * @param readThreads Positive worker count, or below one for min(2, Shared.threads())
	 * @param zipThreads BGZF setting for compressed descriptors; ignored otherwise
	 * @return New two-element array ordered [LF count, byte count], or null as described
	 */
	public static long[] countLinesAndBytes(FileFormat ff, int readThreads, int zipThreads){
		final int oldZT=BgzfSettings.READ_THREADS;
		if(ff.compressed()){
			BgzfSettings.READ_THREADS=(zipThreads>1 ? zipThreads : Tools.mid(1, Shared.threads(), 18));
		}
		FileScanMT fqs=new FileScanMT(ff);
		try{fqs.read(readThreads);}
		catch(IOException e){
			e.printStackTrace();
			return null;
		}finally{BgzfSettings.READ_THREADS=oldZT;}
		long[] ret=new long[]{fqs.totalLines, fqs.totalBytes};
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains the input descriptor; no stream is opened until read() is called. */
	FileScanMT(FileFormat ff_){ff=ff_;}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens by filename, waits for scan workers and adds their totals to this instance.
	 * Calls ReadWrite.getInputStream(name, false, false), rather than forwarding all
	 * descriptor settings. Totals are not reset, so normal repeated calls accumulate.
	 * finishReading's returned status is currently ignored; worker failure instead throws.
	 * @param threads Positive worker count, or below one for min(2, Shared.threads())
	 * @throws IOException Retained by the API; worker IOExceptions instead cause an unchecked scan failure
	 */
	void read(int threads) throws IOException{
		final InputStream is=ReadWrite.getInputStream(ff.name(), false, false);
		threads=(threads<1 ? Math.min(2, Shared.threads()) : threads);//Automatic scan-worker cap is two.
		final ArrayList<ScanThread> alst=new ArrayList<ScanThread>(threads);
		
		for(int i=0; i<threads; i++){
			ScanThread st=new ScanThread(is);
			alst.add(st);
			st.start();
		}
		
		boolean success=true;
		for(ScanThread st : alst){
			while(st.getState()!=Thread.State.TERMINATED){
				try{st.join();}
				catch(InterruptedException e){e.printStackTrace();}
			}
			synchronized(st){
				success&=st.success;
				totalLines+=st.linesT;
				totalBytes+=st.bytesT;
			}
		}
		
		//TODO: Possible bug [stream/FileScanMT#001]: the reported finish/close status is
		//discarded, so totals may be returned despite a reported cleanup error. Not runtime-verified.
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
		if(!success){throw new RuntimeException("Scanning failed.");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns one 1MiB buffer and counters while sharing the serialized input reads. */
	private class ScanThread extends Thread{
		
		/** Retains the shared stream; the outer scanner owns final stream closure. */
		ScanThread(InputStream is_){is=is_;}
		
		/** Marks success only after normal completion; caught IOExceptions leave it false. */
		@Override
		public void run(){
			try{
				process();
				success=true;
			}catch(IOException e){
				e.printStackTrace();
			}
		}
		
		/** Repeatedly fills and counts private chunks until no more bytes are obtained. */
		private void process() throws IOException{
			while(true){
				final int len=fillBuffer();
				if(len<1){break;}
				scanBuffer(len);
			}
		}
		
		/** Fills the private buffer under the shared input lock, stopping at EOF or capacity.
		 * Adds only bytes obtained by this worker to bytesT.
		 * @return Number of bytes read, including zero at EOF
		 */
		private int fillBuffer() throws IOException{
			//Parallel-scan invariant: the synchronized(is) read is the ONLY shared access; each thread reads a
			//DISTINCT chunk into its own 1MiB buffer, then counts '\n' using Vector's selected path outside the lock. The chunks
			//partition the stream exactly -> sum of per-thread linesT/bytesT == total newlines/bytes. The lock
			//serializes only the reads, leaving counting of the private chunks outside the lock.
			synchronized(is){
				int len=0;
				while(len<buffer.length){
					int r=is.read(buffer, len, buffer.length-len);
					if(r<0){break;}
					len+=r;
				}
				bytesT+=len;
				return len;
			}
		}
		
		/** Adds the LF count in buffer[0..len); uses Vector's scalar/SIMD selection.
		 * @param len Valid prefix length; values below one contribute nothing
		 */
		private void scanBuffer(final int len){
			if(len<1){return;}
			int count=Vector.countSymbols(buffer, 0, len, (byte)'\n');
			linesT+=count;
		}
		
		/** Shared stream; all worker reads synchronize on this object. */
		private final InputStream is;
		/** Private chunk storage, counted after releasing the stream lock. */
		private final byte[] buffer=new byte[1024*1024];
		
		/** LF bytes found by this worker. */
		long linesT=0;
		/** InputStream bytes obtained by this worker. */
		long bytesT=0;
		/** True only after process() completes normally. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained input descriptor, used for the name and finishReading policy. */
	private final FileFormat ff;
	/** LF counts already folded from workers into this instance, including partial work. */
	long totalLines=0;
	/** InputStream byte counts already folded from workers, including partial work. */
	long totalBytes=0;

}
