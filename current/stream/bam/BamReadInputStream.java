package stream.bam;

import java.util.ArrayList;

import fileIO.FileFormat;
import stream.FASTQ;
import stream.Read;
import stream.ReadInputStream;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Reads BAM files and outputs Read objects.
 * Converts delegate SamLine batches and optionally links adjacent mate records.
 * Construction does not start the delegate; callers must invoke {@link #start()}.
 * Mutable buffers and counters require a single caller or external synchronization.
 * Close and empty-input limitations are documented at their implementation sites.
 *
 * @author Brian Bushnell
 * @date October 2025
 */
public class BamReadInputStream extends ReadInputStream{

	/**
	 * Prints the first read and its SAM representation from a BAM-default input.
	 * This demonstration assumes a filename argument and a nonempty first batch.
	 * @param args Input filename as the first argument
	 */
	public static void main(String[] args){
		BamReadInputStream bris=new BamReadInputStream(args[0], false, false, false, true);
		bris.start();

		//TODO: Probable bug - stream/bam/BamReadInputStream#004: nextList returns
		//null for an empty first batch, so this demonstration dereferences null.
		//Source-confirmed only; demo repair remains separate.
		Read r=bris.nextList().get(0);
		System.out.println(r.toText(false));
		System.out.println();
		if(r.samline!=null){
			System.out.println(r.samline.toText());
			System.out.println();
		}

		bris.close();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves a BAM-default input with an unlimited delegate read limit.
	 * @param fname Input filename
	 * @param loadHeader_ Whether to request shared header publication
	 * @param ordered_ Whether to request ordered delegate batches
	 * @param interleaved_ Whether to link adjacent first/second mate flags within each batch
	 * @param allowSubprocess_ Whether input format resolution permits subprocesses
	 */
	public BamReadInputStream(String fname, boolean loadHeader_, boolean ordered_, boolean interleaved_, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, allowSubprocess_, false), loadHeader_, ordered_, interleaved_, -1);
	}

	/**
	 * Resolves a BAM-default input and forwards the requested limit unchanged.
	 * @param fname Input filename
	 * @param loadHeader_ Whether to request shared header publication
	 * @param ordered_ Whether to request ordered delegate batches
	 * @param interleaved_ Whether to link adjacent first/second mate flags within each batch
	 * @param allowSubprocess_ Whether input format resolution permits subprocesses
	 * @param maxReads_ Delegate limit; negative requests unlimited input
	 */
	public BamReadInputStream(String fname, boolean loadHeader_, boolean ordered_, boolean interleaved_, boolean allowSubprocess_, long maxReads_){
		this(FileFormat.testInput(fname, FileFormat.BAM, null, allowSubprocess_, false), loadHeader_, ordered_, interleaved_, maxReads_);
	}

	/**
	 * Creates an unstarted reader with a default worker hint and Read conversion enabled.
	 * Warns for a non-BAM descriptor but still asks the factory to select a reader.
	 * This adapter itself consumes SamLine batches and converts them to Reads.
	 * @param ff Nonnull input descriptor compatible with SamLine batch access
	 * @param loadHeader_ Whether to request shared header publication
	 * @param ordered_ Whether to request ordered delegate batches
	 * @param interleaved_ Whether to link adjacent first/second mate flags within each batch
	 * @param maxReads_ Delegate limit forwarded unchanged
	 */
	public BamReadInputStream(FileFormat ff, boolean loadHeader_, boolean ordered_, boolean interleaved_, long maxReads_){
		loadHeader=loadHeader_;
		ordered=ordered_;
		interleaved=interleaved_;
		maxReads=maxReads_;

		stdin=ff.stdio();
		if(!ff.bam()){
			System.err.println("Warning: Did not find expected bam file extension for filename "+ff.name());
		}

		fname=ff.name();
		header=new ArrayList<byte[]>();

		bls=StreamerFactory.makeSamOrBamStreamer(ff, -1, loadHeader, ordered_, maxReads, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the nonnull delegate; this adapter does not guard repeated starts. */
	public void start(){
		if(bls!=null){bls.start();}
	}

	/**
	 * Refills an exhausted buffer unless finished has been set.
	 * On the finished/exhausted path, asserts that at least one entry was generated.
	 * @return Whether the current buffer has an unread entry
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			if(!finished){
				fillBuffer();
			}else{
				//TODO: Probable bug - stream/bam/BamReadInputStream#003: after empty
				//input sets finished with generated==0, this later query asserts instead
				//of returning false. Source-confirmed only; EOF repair is separate.
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	/**
	 * Transfers the current batch, refilling when absent or exhausted.
	 * Empty batches become null. Counts returned list entries, with linked mates
	 * represented by their first read; does not use finished to suppress refills.
	 * @return Next batch of first reads/singletons, or null for a null/empty delegate batch
	 * @throws RuntimeException If the local per-read cursor is nonzero
	 */
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}

	/**
	 * Obtains one delegate batch and converts it, resetting the local cursor first.
	 * A null or empty batch sets finished and installs an empty buffer. IDs and
	 * generated count advance by the number of returned first reads/singletons.
	 */
	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());

		buffer=null;
		next=0;

		// Get next batch of SamLines from Streamer
		ListNum<SamLine> samLines=bls.nextLines();
		if(samLines==null || samLines.size()==0){
			finished=true;
			buffer=new ArrayList<Read>(0);
			return;
		}

		// Convert SamLine objects to Read objects
		buffer=toReadList(samLines.list, nextReadID);
		nextReadID+=buffer.size();
		generated+=buffer.size();
	}

	/**
	 * Converts nonnull SAM entries and attaches their source SamLine objects.
	 * Interleaved mode links adjacent first/second mate flags within this batch only;
	 * it does not compare names or carry a pending mate across batch boundaries.
	 * Each output entry receives one numeric ID; its linked mate shares that ID.
	 * Linked second reads are marked with pair number one for downstream Read consumers.
	 * @param samLines Nonnull batch of nonnull SAM alignments
	 * @param nextReadID2 First numeric ID to assign
	 * @return New list containing singletons and first reads of linked pairs
	 */
	private final ArrayList<Read> toReadList(ArrayList<SamLine> samLines, long nextReadID2){
		ArrayList<Read> list=new ArrayList<Read>(samLines.size());
		for(int idx=0; idx<samLines.size(); idx++){
			SamLine sl1=samLines.get(idx);

			Read r1=sl1.toRead(FASTQ.PARSE_CUSTOM);
			r1.samline=sl1;
			r1.numericID=nextReadID2;
			list.add(r1);

			// Handle paired reads if interleaved
			//TODO: Possible bug [stream/bam/BamReadInputStream#001] - interleaved pairing trusts ADJACENCY +
			//FLAGS only: sl1 is first-of-pair (0x1|0x40) → it pairs sl1 with sl1+1 if sl1+1 is ANY second-of-
			//pair (0x1|0x80), WITHOUT verifying sl1.qname==sl2.qname. Correct for genuinely interleaved input
			//(mates adjacent, the mode's contract), but on a COORDINATE-sorted BAM fed with interleaved=true
			//an unrelated second-of-pair record can be linked. Interleaved mode assumes interleaved data.
			//A 2026-09-30 Java-source search found construction only in this class's main;
			//that does not rule out external callers. A qname-equality guard could check that assumption. NOTE: the pair shares ONE numericID (both r1,r2 = nextReadID2, ++ once) -
			//correct paired semantics; only r1 enters the list (r1.mate=r2), so list.size()=pairs+singletons.
			if(interleaved && (sl1.flag&0x1)!=0 && (sl1.flag&0x40)!=0){
				// This is first of pair, next should be second
				if(idx+1<samLines.size()){
					SamLine sl2=samLines.get(idx+1);
					if((sl2.flag&0x1)!=0 && (sl2.flag&0x80)!=0){
						Read r2=sl2.toRead(FASTQ.PARSE_CUSTOM);
						r2.samline=sl2;
						r2.numericID=nextReadID2;
						r2.setPairnum(1);//STR379: SamLine.toRead leaves mate numbering to its caller.
						r1.mate=r2;
						r2.mate=r1;
						idx++; // Skip the mate since we processed it
					}
				}
			}

			nextReadID2++;
		}
		return list;
	}

	/**
	 * Sets finished without closing the delegate or clearing buffered reads.
	 * @return Always true; this does not reflect the delegate's error state
	 */
	@Override
	public boolean close(){
		//TODO: Probable bug - stream/bam/BamReadInputStream#002: the parent close
		//contract releases input resources and returns an error status. This method
		//does not close bls and always reports true without consulting its state.
		//Source-confirmed contract mismatch; runtime/resource repair is separate.
		finished=true;
		return true;
	}

	/** @throws RuntimeException Always; this reader does not implement restart. */
	@Override
	public synchronized void restart(){
		throw new RuntimeException("BamReadInputStream.restart() not supported - BAM streams cannot be reset");
	}

	/** Returns the filename captured from the input descriptor. */
	@Override
	public String fname(){return fname;}

	/** Returns the requested interleaved mode, not a count of successfully linked pairs. */
	@Override
	public boolean paired(){return interleaved;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Current converted batch; transferred to the caller by nextList. */
	private ArrayList<Read> buffer=null;
	/** Unused local header list; header publication is requested through the delegate. */
	private ArrayList<byte[]> header=null;
	/** Local read cursor, reset to zero; this implementation does not advance it. */
	private int next=0;

	/** Factory-selected reader, started only by start. */
	private final Streamer bls;
	/** Whether within-batch adjacency/flag pairing is requested. */
	private final boolean interleaved;
	/** Requested shared-header retention policy. */
	private final boolean loadHeader;
	/** Captured ordering request, also forwarded to the delegate. */
	private final boolean ordered;
	/** Captured input filename. */
	private final String fname;
	/** Captured delegate limit, forwarded unchanged at construction. */
	private final long maxReads;
	/** Null/empty batch observed, or close called; nextList does not consult this flag. */
	private boolean finished=false;

	/** Converted output entries, counting a linked pair as one. */
	public long generated=0;
	/** Entries transferred by nextList, counting a linked pair as one. */
	public long consumed=0;
	/** Next numeric ID assigned to a converted output entry. */
	private long nextReadID=0;

	/** Whether the input descriptor denotes a standard stream. */
	public final boolean stdin;

}
