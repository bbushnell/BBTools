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
 * Counts lines and bytes in a file.
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 15, 2025
 */
public final class FileScanMT {

	public static void main(String[] args){
		Timer t=new Timer(System.out);
		if(args.length<1) {throw new RuntimeException("Usage: fastqscan.sh filename");}
		String fname=args[0];
		while(fname.startsWith("-")) {fname=fname.substring(1);}
		if(fname.startsWith("in=")) {fname=fname.substring(3);}
		int threads=1;//Math.min(2, Shared.threads());
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
		FileScanMT fqs=new FileScanMT(ff);
		try{fqs.read(threads);}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:   \t");
		System.out.println("Lines:  \t"+fqs.totalLines);
		System.out.println("Bytes:  \t"+fqs.totalBytes);
	}

	public static long[] countLinesAndBytes(String fname, int readThreads, int zipThreads) {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		return countLinesAndBytes(ff, readThreads, zipThreads);
	}
	
	public static long[] countLinesAndBytes(FileFormat ff, int readThreads, int zipThreads) {
		final int oldZT=BgzfSettings.READ_THREADS;
		if(ff.compressed()) {
			BgzfSettings.READ_THREADS=(zipThreads>1 ? zipThreads : Tools.mid(1, Shared.threads(), 18));
		}
		FileScanMT fqs=new FileScanMT(ff);
		try{fqs.read(readThreads);}
		catch(IOException e){
			e.printStackTrace();
			//throw new RuntimeException(e);
			return null;
		}finally {BgzfSettings.READ_THREADS=oldZT;}
		long[] ret=new long[] {fqs.totalLines, fqs.totalBytes};
		return ret;
	}

	FileScanMT(FileFormat ff_){ff=ff_;}

	void read(int threads) throws IOException {
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
				totalLines+=st.linesT;
				totalBytes+=st.bytesT;
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
				bytesT+=len;
				return len;
			}
		}
		
		private void scanBuffer(final int len){
			if(len<1){return;}
			int count=Vector.countSymbols(buffer, 0, len, (byte)'\n');
			linesT+=count;
		}
		
		private final InputStream is;
		private final byte[] buffer=new byte[1024*1024]; // 1MB buffer
		
		long linesT=0;
		long bytesT=0;
		boolean success=false;
	}

	private final FileFormat ff;
	long totalLines=0;
	long totalBytes=0;

}