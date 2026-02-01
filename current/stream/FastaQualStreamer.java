package stream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Single-threaded loader for Fasta + Quality files.
 * Reads two files in lockstep: .fa (bases) and .qual (quality scores).
 * @author Collei, Brian Bushnell
 * @date November 21, 2025
 */
public class FastaQualStreamer implements Streamer {

	public static void main(String[] args) {
		// Usage: java stream.FastaQualStreamer reads.fa reads.qual
		Timer t=new Timer();
		String faName=args[0];
		String qualName=args[1];

		if(args.length>2) {Shared.SIMD=Parse.parseBoolean(args[2]);}

		FileFormat ffFa=FileFormat.testInput(faName, FileFormat.FASTA, null, true, true);
		// FileFormat doesn't explicitly detect .qual usually, but we treat it as generic input here
		FileFormat ffQual=FileFormat.testInput(qualName, FileFormat.UNKNOWN, null, true, true);

		FastaQualStreamer st=new FastaQualStreamer(ffFa, ffQual, 0, -1);
		st.start();
		long reads=0, bases=0;
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
			for(Read r : ln) {
				reads+=r.pairCount();
				bases+=r.pairLength();
			}
		}
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public FastaQualStreamer(FileFormat ff, String qf, int pairnum_, long maxReads_){
		this(ff, FileFormat.testInput(qf, FileFormat.QUAL, null, ff.allowSubprocess(), false), 
			pairnum_, maxReads_);
	}

	/** Constructor. */
	public FastaQualStreamer(FileFormat ffinFa_, FileFormat ffinQual_, int pairnum_, long maxReads_){
		ffinFa=ffinFa_;
		ffinQual=ffinQual_;
		fname=ffinFa_.name(); // Primary name is the fasta file
		pairnum=pairnum_;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		// Simple output queue
		outputQueue=new ArrayBlockingQueue<ListNum<Read>>(QUEUE_SIZE);

		if(verbose){outstream.println("Made FastaQualStreamer for "+fname+" + "+ffinQual_.name());}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void start(){
		if(verbose){outstream.println("FastaQualStreamer.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Start processing thread
		thread=new ProcessThread();
		thread.start();

		if(verbose){outstream.println("FastaQualStreamer started.");}
	}

	@Override
	public void close(){
		if(bfFa!=null) {bfFa.close(); bfFa=null;}
		if(bfQual!=null) {bfQual.close(); bfQual=null;}
	}

	@Override
	public String fname() {return fname;}

	@Override
	public boolean hasMore(){
		return !finished;
	}

	@Override
	public boolean errorState() {return errorState;}

	@Override
	public boolean paired(){return false;} 

	@Override
	public int pairnum(){return pairnum;}

	@Override
	public long readsProcessed() {return readsProcessed;}

	@Override
	public long basesProcessed() {return basesProcessed;}

	@Override
	public void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	@Override
	public ListNum<Read> nextList(){
		try{
			ListNum<Read> list=outputQueue.take();
			if(list==null || list.last()){
				finished=true;
				readsProcessed=thread.readsProcessedT;
				basesProcessed=thread.basesProcessedT;
				errorState=!thread.success;
				if(list!=null) {outputQueue.add(list);}//Re-inject terminal
				return null;
			}
			return list;
		}catch(InterruptedException e){
			errorState=true;
			return null;
		}
	}

	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FastaQualStreamer does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private class ProcessThread extends Thread {

		ProcessThread(){
			setName("FastaQualStreamer-Worker");
		}

		@Override
		public void run(){
			try{
				process();
				success=true;
			}catch(Exception e){
				e.printStackTrace();
				errorState=true;
			}finally{
				// Send terminal list
				try{
					ListNum<Read> terminal=new ListNum<Read>(null, -1, ListNum.LAST);
					outputQueue.put(terminal);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}

		void process() throws InterruptedException{
			bfFa=ByteFile.makeByteFile(ffinFa);
			bfQual=ByteFile.makeByteFile(ffinQual);

			long listNumber=0;
			long readID=0;
			int bytes=0;

			final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
			ListNum<Read> ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);

			// Using LineParser to split space-delimited qualities
			LineParser1 lp=new LineParser1(' ');
			ByteBuilder bbBases=new ByteBuilder();
			ByteBuilder bbQuals=new ByteBuilder();

			// Lookahead buffers
			byte[] nextFaLine=bfFa.nextLine();
			byte[] nextQualLine=bfQual.nextLine();

			while(readID<maxReads && nextFaLine!=null){

				// 1. Process Header
				byte[] faHeader=nextFaLine;
				if(faHeader==null || faHeader.length==0 || faHeader[0]!='>'){
					// Skip garbage/blank lines between records if any, or sync issues
					if(nextFaLine==null) break; //EOF
					nextFaLine=bfFa.nextLine();
					continue;
				}

				// 2. Process Bases
				bbBases.clear();
				while(true){
					nextFaLine=bfFa.nextLine();
					if(nextFaLine==null) break;
					if(nextFaLine.length>0 && nextFaLine[0]=='>'){
						break; // Found next header, break loop, hold in nextFaLine
					}
					bbBases.append(nextFaLine);
				}

				// 3. Process Qual Header
				// Skip garbage until we find a header
				while(nextQualLine!=null && (nextQualLine.length==0 || nextQualLine[0]!='>')){
					nextQualLine=bfQual.nextLine();
				}

				if(nextQualLine==null){
					throw new RuntimeException("FASTA file has more records than QUAL file (EOF reached in QUAL).");
				}

				// Validate Headers match (loose check on ID length or content could go here)
				// For speed, we assume sync if ordered correctly.

				// 4. Process Qual Scores
				bbQuals.clear();
				while(true){
					nextQualLine=bfQual.nextLine();
					if(nextQualLine==null) break;
					if(nextQualLine.length>0 && nextQualLine[0]=='>'){
						break; // Found next header
					}

					// Parse line of space-delimited integers
					lp.set(nextQualLine);
					int terms=lp.terms();
					for(int i=0; i<terms; i++){
						int q=lp.parseInt(i);
						// Store as raw byte value (numeric Phred score)
						bbQuals.append((byte)q);
					}
				}

				// 5. Integrity Check
				if(bbBases.length != bbQuals.length) {
					throw new RuntimeException("Read "+readID+": Sequence length ("+bbBases.length+
						") != Quality length ("+bbQuals.length+").\nHeader: "+new String(faHeader));
				}

				if(samplerate>=1f || randy.nextFloat()<samplerate){
					// Fasta headers start with '>', strip it for the Read object
					String id=new String(faHeader, 1, faHeader.length-1);
					Read r=new Read(bbBases.toBytes(), bbQuals.toBytes(), id, readID);
					r.setPairnum(pairnum);
					if(!r.validated()){r.validate(true);} //Validate to be safe

					ln.add(r);
					bytes+=r.length();
				}
				readID++;

				readsProcessedT++;
				basesProcessedT+=bbBases.length;

				if(ln.size()>=slimit || bytes>=blimit){
					outputQueue.put(ln);
					ln=new ListNum<Read>(new ArrayList<Read>(slimit), listNumber++);
					bytes=0;
				}
			}

			if(ln.size()>0){
				outputQueue.put(ln);
			}
			bfFa.close();
			bfQual.close();
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** True only if this thread has completed successfully */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final String fname;
	final FileFormat ffinFa;
	final FileFormat ffinQual;

	final ArrayBlockingQueue<ListNum<Read>> outputQueue;
	private ProcessThread thread;

	private ByteFile bfFa;
	private ByteFile bfQual;

	final int pairnum;

	protected long readsProcessed=0;
	protected long basesProcessed=0;

	final long maxReads;

	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	public static int TARGET_LIST_BYTES=262144;
	private static final int QUEUE_SIZE=4;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	protected PrintStream outstream=System.err;
	public static final boolean verbose=false;
	public boolean errorState=false;
	private float samplerate=1f;
	private shared.Random randy=null;

}