package stream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import shared.Random;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/**
 * Produces reads and bases from a sequence file, with low overhead.
 * @author Brian Bushnell
 * @contributor Collei
 * @date November 29, 2025
 */
public final class FastqScanStreamer implements Streamer{

	public static void main(String[] args) {
		Timer t=new Timer();
		if(args.length<1 || args.length>2) {
			System.err.println("Usage: FastqScanStreamer filename");
		}
		String fname=args[0];
		if(args.length>1) {Read.SKIP_SLOW_VALIDATION=!Parse.parseBoolean(args[1]);}
		while(fname.startsWith("-")) {fname=fname.substring(1);}
		if(fname.startsWith("in=")) {fname=fname.substring(3);}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false);
		if(ff.stdin()) {
			//Do nothing
		}else{
			File f=new File(fname);
			if(!f.isFile() || !f.canRead()) {
				throw new RuntimeException("Can't read "+fname);
			}
		}
		FastqScanStreamer fqs=new FastqScanStreamer(ff, 0, -1);
		fqs.start();
		for(ListNum<Read> ln=fqs.nextList(); ln!=null; ln=fqs.nextList()) {
			//Do nothing
		}
		fqs.close();
		t.stop("Time:   \t");
		System.err.println("Records:\t"+fqs.totalRecords);
		System.err.println("Bases:  \t"+fqs.totalBases);
		if(ff.samOrBam()) {System.err.println("Headers:\t"+fqs.totalHeaders);}
		ByteBuilder bb=fqs.corruption();
		if(fqs.slashrLines>0) {
			System.err.println("Contained Windows-style \r\n");
		}
		if(bb!=null) {
			System.err.print(bb);
			System.exit(1);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	FastqScanStreamer(String fname, int pairnum_, long maxReads_) {
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, true, false), pairnum_, maxReads_);
	}

	FastqScanStreamer(FileFormat ff_, int pairnum_, long maxReads_) {
		ff=ff_;
		interleaved=ff.interleaved();
		pairnum=pairnum_;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	void readFastq() throws IOException {
		is=ReadWrite.getInputStream(ff.name(), false, false);
		IntList newlines=new IntList(8192);
		for(int r=is.read(buffer); r>0 || bstop>0; r=is.read(buffer, bstop, buffer.length-bstop)) {
			assert(bstart==0);
			r=Math.max(r, 0);
			bstop+=r;
			if(r==0 && buffer[bstop-1]!='\n') {
				if(bstop>=buffer.length) {expand();}
				buffer[bstop++]='\n';
				missingTerminalNewline=true;
			}
			Vector.findSymbols(buffer, 0, bstop, (byte)'\n', newlines.clear());
			final int records=(!interleaved ? newlines.size/4 : (newlines.size/4)&(~1));
			totalRecords+=records;
			final ArrayList<Read> reads;
			if(interleaved) {
				reads=makeReadsInterleaved(records, newlines);
			}else {
				reads=makeReadsSingle(records, newlines);
			}
			if(reads!=null && !reads.isEmpty()) {
				if(samplerate<1) {sample(reads);}
				ListNum<Read> ln=new ListNum<Read>(reads, nextLID++);
				boolean b=add(ln);
				if(!b) {break;}
			}
			
			final int residue=bstop-bstart;
			if(residue>0) {
				if(bstart>0) {
					System.arraycopy(buffer, bstart, buffer, 0, residue);
				}else if(r>0){
					expand();
				}
			}
			bstart=0;
			bstop=residue;
			if(r<1) {
				if(residue>0) {partialRecords++;}
				break;
			}
			if(nextRID>=maxReads) {break;}
		}
	}
	
	private ArrayList<Read> makeReadsSingle(int records, IntList newlines){
		int recordStart=0;
		long maxRecords=Math.min(records, maxReads-nextRID);
		records=(int)maxRecords;
		ArrayList<Read> reads=new ArrayList<Read>(records);
		final int offset=-FASTQ.ASCII_OFFSET;
		for(int i=0, j=0; i<records; i++, j+=4) {
			final int headerEnd=newlines.get(j);
			final int basesEnd=newlines.get(j+1);
			final int plusEnd=newlines.get(j+2);
			final int recordEnd=newlines.get(j+3);
			final int slashr0=(buffer[headerEnd-1]=='\r') ? 1 : 0;
			final int slashr1=(buffer[basesEnd-1]=='\r') ? 1 : 0;
			final int slashr2=(buffer[recordEnd-1]=='\r') ? 1 : 0;
			final int headerLen=headerEnd-recordStart-1-slashr0;
			final int basesLen=basesEnd-headerEnd-1-slashr1;
			final int qualsLen=recordEnd-plusEnd-1-slashr2;
			
			final String header=new String(buffer, recordStart+1, headerLen);
			final byte[] bases=Arrays.copyOfRange(buffer, headerEnd+1, basesEnd-slashr1);
			final byte[] quals=Arrays.copyOfRange(buffer, plusEnd+1, recordEnd-slashr2);
			Vector.applyQualOffset(quals, bases, offset);
			final Read r=new Read(bases, quals, header, nextRID);
			r.setPairnum(pairnum);
			reads.add(r);
			nextRID++;
			
			slashrLines+=slashr1+slashr2;
			totalBases+=basesLen;
			bstart=recordEnd+1;
			qualMismatch|=(qualsLen!=basesLen);
			missingAt|=(buffer[recordStart]!='@');
			missingPlus|=(buffer[basesEnd+1]!='+');
			recordStart=recordEnd+1;
		}
		return reads;
	}
	
	private ArrayList<Read> makeReadsInterleaved(int records, IntList newlines){
		int recordStart=0;
		final long remainingRecords=(maxReads-nextRID);
		long maxRecords=(records<remainingRecords ? records : Math.min(records, remainingRecords*2));
		records=(int)maxRecords;
		ArrayList<Read> reads=new ArrayList<Read>(records);
		final int offset=-FASTQ.ASCII_OFFSET;
		for(int i=0, j=0; i<records; i++, j+=4) {
			final int headerEnd=newlines.get(j);
			final int basesEnd=newlines.get(j+1);
			final int plusEnd=newlines.get(j+2);
			final int recordEnd=newlines.get(j+3);
			final int slashr0=(buffer[headerEnd-1]=='\r') ? 1 : 0;
			final int slashr1=(buffer[basesEnd-1]=='\r') ? 1 : 0;
			final int slashr2=(buffer[recordEnd-1]=='\r') ? 1 : 0;
			final int headerLen=headerEnd-recordStart-1-slashr0;
			final int basesLen=basesEnd-headerEnd-1-slashr1;
			final int qualsLen=recordEnd-plusEnd-1-slashr2;
			
			final String header=new String(buffer, recordStart+1, headerLen);
			final byte[] bases=Arrays.copyOfRange(buffer, headerEnd+1, basesEnd-slashr1);
			final byte[] quals=Arrays.copyOfRange(buffer, plusEnd+1, recordEnd-slashr2);
			Vector.applyQualOffset(quals, bases, offset);
			final Read r=new Read(bases, quals, header, nextRID);
			r.setPairnum(pairnum);
			reads.add(r);
			nextRID+=(i&1);
			
			slashrLines+=slashr1+slashr2;
			totalBases+=basesLen;
			bstart=recordEnd+1;
			qualMismatch|=(qualsLen!=basesLen);
			missingAt|=(buffer[recordStart]!='@');
			missingPlus|=(buffer[basesEnd+1]!='+');
			recordStart=recordEnd+1;
		}
		pair(reads);
		return reads;
	}
	
	private void pair(ArrayList<Read> reads){
		final int lim=reads.size()&~1;
		if(reads.size()!=lim) {
			System.err.println("Warning: File "+ff.name()+" was processed"
				+ " as interleaved but had an odd number of reads.");
		}
		if(lim<1) {return;}
		for(int i=0; i<lim; i+=2) {
			Read r1=reads.get(i), r2=reads.set(i+1, null);
			r1.mate=r2;
			r2.mate=r1;
			r2.numericID=r1.numericID;
			r2.setPairnum(1);
		}
		Tools.condenseStrict(reads);
	}
	
	private void sample(ArrayList<Read> reads) {
		assert(samplerate<1);
		for(int i=0; i<reads.size(); i++) {
			if(randy.nextFloat()>samplerate) {
				reads.set(i, null);
			}
		}
		Tools.condenseStrict(reads);
	}
	
	private boolean add(ListNum<Read> ln) {
		while(ln!=null && !closed) {
			try{
				queue.put(ln);
				ln=null;
			}catch(InterruptedException e){
				if(closed) {return false;}
				throw new RuntimeException(e);
			}
		}
		return true;
	}

	private void expand() {
		long newlen=Math.min(buffer.length*2L, Shared.MAX_ARRAY_LEN);
		assert(newlen>buffer.length) : "Record "+totalRecords+" is too long.";
		buffer=Arrays.copyOf(buffer, (int)newlen);
	}
	
	public ByteBuilder corruption() {
		if(partialRecords<1 && !qualMismatch && !missingTerminalNewline && !missingPlus && !missingAt) {
			return null;
		}
		ByteBuilder bb=new ByteBuilder();
		if(partialRecords>0 || missingAt || missingPlus || qualMismatch) {
			bb.appendln("At least "+Math.max(partialRecords, 1)+" corrupt records.");
		}
		if(partialRecords>0) {bb.appendln("At least "+partialRecords+" incomplete records.");}
		if(qualMismatch) {bb.appendln("At least "+1+" base/quality mismatches.");}
		if(missingAt) {bb.appendln("At least "+1+" missing @ symbols.");}
		if(missingPlus) {bb.appendln("At least "+1+" missing + symbols.");}
		if(missingTerminalNewline) {bb.appendln("Missing terminal newline.");}
		assert(bb.length()>0);
		return bb;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Reader Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ReaderRunnable implements Runnable{
		@Override
		public void run(){
			try{
				readFastq();
			}catch(Throwable e){
				e.printStackTrace();
				synchronized(FastqScanStreamer.this) {errorState=true;}
			}
			try{
				poisonAndClose();
			}catch(Throwable e){
				e.printStackTrace();
				synchronized(FastqScanStreamer.this) {errorState=true;}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Thread Control        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public synchronized void start(){
		if(started) {return;}
		started=true;
		readerThread=new Thread(new ReaderRunnable());
		readerThread.start();
		if(verbose){System.err.println("Started "+getClass().getName());}
	}

	@Override
	public synchronized void close(){
		if(closed) {return;}
		readerThread.interrupt();
		ReadWrite.finishReading(is, ff.name(), ff.allowSubprocess());
		closed=true;
	}
	
	public synchronized void poisonAndClose(){
		poison();
		close();
	}
	
	private synchronized void poison() {
		if(poisoned) {return;}
		add(new ListNum<Read>(null, nextLID++, ListNum.POISON));
		poisoned=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public String fname(){return ff.name();}

	@Override
	public boolean paired(){return interleaved;}

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public long readsProcessed(){return totalRecords/(interleaved ? 2 : 1);}

	@Override
	public long basesProcessed(){return totalBases;}

	@Override
	public synchronized void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.random(seed));
	}

	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> ln=null;
		while(ln==null) {
			try{
				ln=queue.take();
				if(ln.poison()) {queue.put(ln);}
			}catch(InterruptedException e){
				if(drained) {return null;}
			}
		}
		if(ln.poison()) {
			synchronized(this) {
				drained=true;
				close();
			}
		}
		return ln.poison() ? null : ln;
	}

	@Override
	public ListNum<SamLine> nextLines(){
		assert(!interleaved);
		throw new RuntimeException("Not supported");
	}

	@Override
	public boolean hasMore(){return !drained;}

	@Override
	public synchronized boolean errorState(){
		return errorState || partialRecords>0 || qualMismatch || missingPlus || missingAt;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Final Fields          ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ff;
	private final boolean interleaved;
	private final long maxReads;
	private final int pairnum;
	private final ArrayBlockingQueue<ListNum<Read>> queue=new ArrayBlockingQueue<ListNum<Read>>(4);

	/*--------------------------------------------------------------*/
	/*----------------        Mutable Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** For emergency shutdown */
	private volatile boolean closed=false;
	private volatile boolean started=false;
	private volatile boolean poisoned=false;
	private volatile boolean drained=false;
	private boolean errorState=false;
	private int bufferLen=262144;
	private byte[] buffer=new byte[bufferLen];
	private InputStream is;
	private Thread readerThread;
	
	private int bstart=0, bstop=0;
	private long nextRID=0;
	private long nextLID=0;
	
	private float samplerate=1f;
	private Random randy;
	
	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/
	
	public long totalHeaders;
	public long totalRecords;
	public long totalBases;
	
	public long partialRecords;
	public long slashrLines;
	public boolean qualMismatch;
	public boolean missingTerminalNewline;
	public boolean missingPlus;
	public boolean missingAt;

	/*--------------------------------------------------------------*/
	/*----------------            Stats             ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final boolean verbose=false;
	
}
