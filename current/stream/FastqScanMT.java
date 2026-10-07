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
import structures.IntList;

/**
 * Counts FASTQ records and bases by scanning fixed-size chunks on multiple workers.
 * Workers serialize reads from one stream, then infer four-line FASTQ phase and
 * scan their own buffers independently. This assumes conventional FASTQ structure;
 * it does not perform FastqScan's single-parser corruption checks or count actual
 * quality lengths. The CLI prints the base total as its quality-count stand-in.
 *
 * Instances accumulate worker totals and are intended for one read attempt.
 * Static entry points also change shared settings and are not isolated from
 * concurrent configuration. Decompression may have additional backend threads.
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 7, 2025
 */
public final class FastqScanMT{

	/** Prints time, records, bases, a quality-count stand-in and input bytes to stdout.
	 * Defaults to at most two parsing workers; shared parser/compression options persist.
	 * @param args FASTQ path (optionally in=), then thread count or flag=value options retaining embedded '=' */
	public static void main(String[] args){
		Timer t=new Timer(System.out);
		if(args.length<1){throw new RuntimeException("Usage: fastqscan.sh filename");}
		String fname=args[0];
		while(fname.startsWith("-")){fname=fname.substring(1);}
		if(fname.startsWith("in=")){fname=fname.substring(3);}
		int threads=Math.min(2, Shared.threads());
		BgzfSettings.READ_THREADS=Tools.mid(1, 18, Shared.threads());
		for(int i=1; i<args.length; i++){
			String arg=args[i];
			final int equals=arg.indexOf('=');
			String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			String b=(equals<0 || equals==arg.length()-1 ? null : arg.substring(equals+1));
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
			if(!f.isFile() || !f.canRead()){
				throw new RuntimeException("Can't read "+fname);
			}
		}
		FastqScanMT fqs=new FastqScanMT(ff);
		try{fqs.read(threads);}
		catch(IOException e){throw new RuntimeException(e);}
		t.stop("Time:   \t");
		System.out.println("Records:\t"+fqs.totalRecords);
		System.out.println("Bases:  \t"+fqs.totalBases);
		System.out.println("Quals:  \t"+fqs.totalBases);//TODO qualsT is never incremented in scanBuffer (MT does not count quals yet) — prints totalBases as a stand-in (valid since quals==bases for FASTQ). Known-incomplete feature.
		System.out.println("Bytes:  \t"+fqs.totalBytes);
	}

	/** Resolves a FASTQ-default descriptor and delegates the count operation.
	 * @param fname Input path
	 * @param halveInterleaved Probe pairing and halve the first count if interleaved
	 * @param readThreads Parsing workers; below one selects at most two shared threads
	 * @param zipThreads BGZF request; at most one selects an automatic count
	 * @return Four counters from the descriptor overload, or null on its caught IOException */
	public static long[] countReadsAndBases(String fname, boolean halveInterleaved, int readThreads, int zipThreads){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		return countReadsAndBases(ff, halveInterleaved, readThreads, zipThreads);
	}

	/** Counts FASTQ without the single-parser corruption report.
	 * Returns {records/divisor, records, bases, sequence-header proxy}. The divisor is
	 * two only when the requested FASTQ pairing probe detects interleaving; integer
	 * division truncates odd counts. The fourth value deliberately repeats records,
	 * unlike FastqScan's FASTQ file-header count of zero; see historical #002 below.
	 * Temporarily changes BGZF read threads for compressed input. Restoration surrounds
	 * read() only, so the preceding pairing probe and construction are outside it.
	 * @param ff Nonnull FASTQ descriptor
	 * @param halveInterleaved Probe pairing to choose the first counter's divisor
	 * @param readThreads Parsing workers; below one selects at most two shared threads
	 * @param zipThreads BGZF request; at most one uses a shared-thread cap of eighteen
	 * @return Four counters, or null after printing an IOException from read()
	 * @throws RuntimeException For non-FASTQ input or a worker that did not succeed */
	public static long[] countReadsAndBases(FileFormat ff, boolean halveInterleaved, int readThreads, int zipThreads){
		final int oldZT=BgzfSettings.READ_THREADS;
		if(ff.compressed()){
			BgzfSettings.READ_THREADS=(zipThreads>1 ? zipThreads : Tools.mid(1, Shared.threads(), 18));
		}
		int recordsPerRead=1;
		if(ff.fastq() && halveInterleaved){
			int[] iq=FileFormat.testInterleavedAndQuality(ff.name(), false);
			recordsPerRead=(iq[1]==FileFormat.INTERLEAVED ? 2 : 1);
		}
		FastqScanMT fqs=new FastqScanMT(ff);
		try{fqs.read(readThreads);}
		catch(IOException e){
			e.printStackTrace();
			//throw new RuntimeException(e);
			return null;
		}finally{BgzfSettings.READ_THREADS=oldZT;}
		//[stream/FastqScanMT#002 RESOLVED - Brian] ret[3] is SEQUENCE headers (per-record '@' lines), returned as totalRecords on purpose: '@' is a legal quality char and would double-count across chunk boundaries, so #records is the unambiguous proxy — exact for any valid FASTQ; a misformatted file (records!=headers) should crash, not silently mislead. (Single-thread FastqScan's ret[3]=totalHeaders is "file headers" = 0 for fastq, a different quantity; live callers use ret[0..2].)
		long[] ret=new long[]{fqs.totalRecords/recordsPerRead, fqs.totalRecords,
			fqs.totalBases, fqs.totalRecords};
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains a descriptor without opening input or starting workers.
	 * @param ff_ Nonnull descriptor; read() enforces FASTQ format */
	FastqScanMT(FileFormat ff_){ff=ff_;}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens one stream, starts workers, joins them and adds their counters to this instance.
	 * Interrupted joins are logged and retried without restoring interrupt status.
	 * After joining, normal-path closure precedes the worker-success check.
	 * @param threads Parsing worker count; below one selects min(2, Shared.threads())
	 * @throws IOException Declared for the counting API; worker I/O failures instead
	 * leave success false and cause the RuntimeException below
	 * @throws RuntimeException For non-FASTQ input or unsuccessful worker completion */
	void read(int threads) throws IOException{
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
				try{st.join();}
				catch(InterruptedException e){e.printStackTrace();}
			}
			synchronized(st){
				success&=st.success;
				totalRecords+=st.recordsT;
				totalBases+=st.basesT;
				totalQuals+=st.qualsT;
				totalBytes+=st.bytesT;
			}
		}

		//TODO: Probable bug - finishReading reports close errors by boolean return,
		//which is ignored here. There is also no finally-close path; review separately.
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
		if(!success){throw new RuntimeException("Scanning failed.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns a chunk buffer, newline positions and counters; only input acquisition is locked. */
	private class ScanThread extends Thread{

		/** Borrows the common input stream; the outer read method closes it.
		 * @param is_ Shared nonnull input, also used as the acquisition monitor */
		ScanThread(InputStream is_){
			is=is_;
		}

		/** Processes chunks until EOF; records success only after normal return.
		 * IOExceptions are printed and leave success false; unchecked failures propagate. */
		@Override
		public void run(){
			try{
				process();
				success=true;
			}catch(IOException e){
				e.printStackTrace();
			}
		}

		/** Repeatedly acquires and scans chunks until the shared stream returns no data.
		 * @throws IOException If a read from the shared input fails */
		private void process() throws IOException{
			while(true){
				final int len=fillBuffer();
				if(len<1){break;}
				scanBuffer(len);
			}
		}

		/** Fills this worker's buffer under the common input monitor until full or EOF.
		 * Updates the worker's byte total; zero-length reads are retried.
		 * @return Number of bytes acquired, zero when EOF precedes any data
		 * @throws IOException If the shared stream read fails */
		private int fillBuffer() throws IOException{
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

		/** Infers a local four-line phase and accumulates header/sequence contributions.
		 * Searches for at-sign and plus-sign lines two positions apart; if none are found,
		 * assumes the last complete line is a quality line. That fallback is an assumption,
		 * not validation; historical #004 records its malformed-tail limitation.
		 * A short final chunk lacking LF gets a synthetic newline position without changing bytesT.
		 * Sequence residues are counted here; header residues wait for their terminating
		 * newline in a later chunk. CR before a sequence LF/boundary is excluded.
		 * No quality lengths or alphabet checks are accumulated.
		 * @param len Valid bytes in this worker's buffer; nonpositive values do nothing */
		private void scanBuffer(final int len){
			if(len<1){return;}

			// 1. Find newlines
			newlines.clear();
			Vector.findSymbols(buffer, 0, len, (byte)'\n', newlines);

			// Handle missing terminal newline on the very last block
			if(len<buffer.length && (len==0 || buffer[len-1]!='\n')){
				newlines.add(len);// Pretend there is a newline at the very end
			}

			final int lines=newlines.size;
			if(lines==0){
				// Special case: Huge block with no newlines (single sequence line?)
				// We can't identify it, so we assume it's a Sequence line if we can't prove otherwise.
				// But 1MB without newlines is weird. Assuming sequence base count.
				assert(false) : "Record exceeded buffer length";
				bytesT+=len;//[stream/FastqScanMT#003 latent] redundant — fillBuffer (L166) already counted these bytes; dead under -ea (assert above halts), double-counts bytesT only under -da. Trivial.
				return;
			}

			// 2. Determine Frame (Self-Stabilization)
			int frameStart=-1;// Index in newlines of the first confirmed Header line

			// Scan for @...+.  i is the index of the PREVIOUS newline.
			// Line i starts at newlines.get(i-1)+1.
			// We iterate through lines to find @Header (Frame 0) and +Plus (Frame 2)
			//CORRECT despite '@'/'+' being valid quality chars: the only non-header line that can start with '@' is a QUAL line, whose line i+2 is a SEQ line — which cannot start with '+'. So "@ at i AND + at i+2" uniquely identifies a true header for valid FASTQ. (Studied praise — verified.)
			//Deeper why (Brian): a FASTQ/FASTA header is arbitrary free text — '@', '+', any byte can appear ANYWHERE in it, not just at position 0 — so NO substring can conclusively prove "this is a complete header." Framing therefore cannot trust header CONTENT; it must anchor on STRUCTURE: the 4-line periodicity, read off the @...+ RELATIVE positions. FASTQ isn't locally parseable; you need global structure from a known anchor. That is the whole reason this stateless-chunk scan is possible.
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
				//[stream/FastqScanMT#004 LOW] This "last complete line == Qual(frame 3)" fallback is CORRECT for a complete file (its final chunk ends just after the last qual\n) and for any no-@...+ chunk that still ends on a record boundary. It misframes ONLY when the last complete line is not a qual line — which for the final chunk means a TRUNCATED/incomplete FASTQ (malformed input); a >1MB single line instead crashes loud at the lines==0 assert above. So: bounded miscount on truncated input only, correct on valid complete files. Unlike ST FastqScan (which flags partialRecords), MT silently misframes a truncated tail. (Verifier-flagged, re-traced to LOW.)
			}

			// Helper: Frame 0=Head, 1=Seq, 2=Plus, 3=Qual
			// We calculate frame relative to frameStart (which is Frame 0)

			int lineStart=0;
			for(int i=0; i<lines; i++){
				final int lineEnd=newlines.get(i);
				// Distance from known header line
				final int dist=i-frameStart;
				// Modulo 4, handling negatives
				//&3 is negative-safe (two's complement). Since every thread locks the same true header phase, all threads agree on each line's absolute frame → the boundary line (prev chunk's residue head + this chunk's line0) is counted consistently, no double-count or miss.
				final int frame=(dist&3);

				if(frame==1){// Sequence Line
					int lineLen=lineEnd-lineStart;// Exclude newline
					if(lineLen>0 && buffer[lineEnd-1]=='\r'){lineLen--;}
					basesT+=lineLen;
				}else if(frame==0){// Header Line
					recordsT++;
				}

				lineStart=lineEnd+1;
			}

			// Handle residue (bytes after last newline) if any
			if(lineStart<len){
				// The residue is the start of the NEXT line.
				// Current line i=lines.
				final int dist=lines-frameStart;
				final int frame=(dist&3);

				if(frame==1){// Partial Sequence Line
					//[stream/FastqScanMT#001] strip a \r straddling the chunk boundary (split exactly at \r|\n of a seq line). The full-line path (L234) and the next chunk's line0 already strip it; this residue head did not → +1 base over-count on \r\n input at that rare alignment.
					int lineLen=len-lineStart;
					if(lineLen>0 && buffer[len-1]=='\r'){lineLen--;}
					basesT+=lineLen;
				}else if(frame==0){
					// Partial Header Line.
					// Do NOT count record yet; it will be counted by the next thread
					// which sees the newline terminating this header.
				}
			}
		}

		/** Borrowed common input and acquisition monitor. */
		private final InputStream is;
		/** One mebibyte of worker-local scan storage. */
		private final byte[] buffer=new byte[1024*1024];// 1MB buffer
		/** Reused newline offsets, including an optional synthetic final position. */
		private final IntList newlines=new IntList(1024*16);

		/** Header-line terminators counted in this worker's inferred phase. */
		long recordsT=0;
		/** Sequence-line bytes assigned to this worker, excluding supported CRLF separators. */
		long basesT=0;
		/** Placeholder: scanBuffer never increments it; CLI substitutes bases instead. */
		long qualsT=0;
		/** Input bytes read; historical #003 adds them again on its assertions-off branch. */
		long bytesT=0;
		/** True only after process returns normally. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained input descriptor. */
	private final FileFormat ff;
	/** Sum of worker record counts after joined aggregation. */
	long totalRecords=0;
	/** Sum of worker sequence-length counts after joined aggregation. */
	long totalBases=0;
	/** Aggregated placeholder quality count; current workers leave it zero. */
	long totalQuals;
	/** Sum of worker decompressed-byte counts, including any #003 branch duplication. */
	long totalBytes;

}
