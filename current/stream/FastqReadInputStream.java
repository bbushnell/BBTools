package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import hiseq.IlluminaHeaderParser2;
import shared.Shared;
import structures.ByteBuilder;

/** Legacy batch FASTQ reader backed by ByteFile and the shared FASTQ parser.
 * Captures descriptor pairing, amino mode, header shrinking and batch capacity at
 * construction. Other FASTQ parsing options and its error flag remain shared globals.
 * Use one consumer; synchronization of batch/reset methods does not make every
 * public method safe for concurrent use or configuration changes.
 * Counters measure list entries, normally one per pair in interleaved mode; parser
 * options can alter mate retention/orientation. Call close for backend error folding,
 * then errorState when the shared parser status is also required.
 * @author Brian Bushnell, Shinobu
 */
public class FastqReadInputStream extends ReadInputStream{
	
	/** Legacy demo that prints the first entry of one batch; assumes nonempty input. */
	public static void main(final String[] args){
		
		FastqReadInputStream fris=new FastqReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a descriptor with FASTQ as the fallback format and delegates.
	 * @param fname Input path or standard-input name
	 * @param allowSubprocess_ Permit subprocess-assisted input through the descriptor
	 */
	public FastqReadInputStream(final String fname, final boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false));
	}
	
	/** Creates the byte backend and captures this wrapper's configuration.
	 * Retains a legacy unused filename split when PARSE_CUSTOM is enabled; actual
	 * record parsing remains delegated to FASTQ, including its live global settings.
	 * @param ff Nonnull descriptor supplying name, pairing and backend options
	 */
	public FastqReadInputStream(final FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		shrinkHeaders=FASTQ.SHRINK_HEADERS;
		if(!ff.fastq()){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+ff.name());
		}
		
		if(FASTQ.PARSE_CUSTOM){
			try{
				final String[] s=ff.name().split("_");
			}catch(final Exception e){
				//Preserve historical tolerance of errors from this unused filename hook.
			}
		}
		interleaved=ff.interleaved();
		tf=ByteFile.makeByteFile(ff);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Prefetches through the backend when necessary and checks for a buffered entry.
	 * If no buffered entry remains and the backend is closed, retains the historical
	 * generated>0 assertion. That counter resets on restart; a repeated hasMore on an
	 * exhausted empty input can therefore assert.
	 * @return true if an unread entry is buffered
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}
	
	/** Transfers the next batch without reusing the returned list and counts its entries.
	 * Reads and mates are retained references, not copies. EOF's empty batch becomes null.
	 * @return Next nonempty list, or null when the delegated parser returns no entries
	 * @throws RuntimeException If the legacy single-entry cursor indicates mixed access
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
	
	/** Fills at most BUF_LEN entries, advances generated/ID counters and optionally shrinks IDs.
	 * Closes the byte backend after a short batch; an exact full final batch requires
	 * another fill to observe EOF. The backend close result is not folded here; callers
	 * must still call public close. A defensive null result latches a local error,
	 * whereas the current parser's ordinary EOF result is an empty list.
	 */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		buffer=FASTQ.toReadList(tf, BUF_LEN, nextReadID, interleaved, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		//A short batch means exhaustion. An exact multiple needs one extra empty fill before closure.
		if(bsize<BUF_LEN){tf.close();}

		generated+=bsize;
		//Preserve the null-versus-empty distinction: only a null result trips the defensive error branch.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}else if(shrinkHeaders){//Regenerates header strings only when the captured option is enabled.
			for(Read r : buffer){shrinkHeader(r, r.mate);}
		}
	}
	
	/** Rewrites shrinkable Illumina IDs as coordinates with first/second read markers.
	 * Tests only r1's header with the delegated parser. Both new IDs use r1's lane,
	 * tile and x/y coordinates; an unshrinkable first header leaves both reads unchanged.
	 * Reuses instance parsing/building scratch while fillBuffer owns the reader monitor.
	 * @param r1 Nonnull first or unpaired read
	 * @param r2 Optional mate whose ID is rewritten from the same coordinates
	 */
	private void shrinkHeader(final Read r1, final Read r2){
		ihp.parse(r1.id);
		if(!ihp.canShrink()){return;}
		bbh.clear().colon().colon().colon();
		ihp.appendCoordinates(bbh).space().append(1).colon();
		r1.id=bbh.toString();
		if(r2!=null){
			bbh.set(bbh.length-2, (byte)'2');
			r2.id=bbh.toString();
		}
	}
	
	/** Closes the byte backend and folds its reported status into the sticky local flag.
	 * Does not clear a prefetched batch or incorporate FASTQ's shared error flag;
	 * errorState queries that flag separately. This method is not synchronized here.
	 * @return Local error state after including the backend's close result
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/** Clears counters/cursor/batch state and delegates reopening to the byte backend.
	 * Retains local/shared error flags, captured configuration and header scratch.
	 * Backend reset semantics determine support for standard streams or other sources.
	 */
	@Override
	public synchronized void restart(){
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	/** Returns captured input pairing mode; shared parser options still govern retained mates. */
	@Override
	public boolean paired(){return interleaved;}
	
	/** Returns the byte backend's input name. */
	@Override
	public String fname(){return tf.name();}
	
	/** Reports the local flag or FASTQ's shared flag, which can reflect other readers.
	 * Does not directly query the byte backend; public close folds its reported errors.
	 */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Prefetched entries; detached on nextList so caller-owned lists remain intact. */
	private ArrayList<Read> buffer=null;
	/** Legacy single-entry cursor; current batch methods reset and require zero. */
	private int next=0;
	
	/** Backend selected by ByteFile.makeByteFile from the descriptor and global settings. */
	private final ByteFile tf;
	/** Pairing mode captured from the descriptor. */
	private final boolean interleaved;
	/** Additional Read flags captured from Shared.AMINO_IN. */
	private final int flag;
	/** Header-shrinking option captured once at construction. */
	private boolean shrinkHeaders;
	
	/** Reused builder for replacement first/second mate IDs. */
	private final ByteBuilder bbh=new ByteBuilder(128);
	/** Reused first-header parser; accessed by shrinkHeader under fillBuffer's monitor. */
	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();

	/** Maximum list entries per batch, captured at construction. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Unused captured data-cap placeholder; fillBuffer limits only entry count.
	 * TODO: The original proposal for very long FASTQ reads required disabling the
	 * cap for paired ends. No such base/data limit is implemented here.
	 */
	private final long MAX_DATA=Shared.bufferData();

	/** Parsed list entries since construction/restart, including prefetched entries. */
	public long generated=0;
	/** Entries transferred through nextList since construction/restart. */
	public long consumed=0;
	/** Starting numeric ID passed to the next parser call; advances by list entries. */
	private long nextReadID=0;
	
	/** Descriptor's captured standard-stream indicator; does not imply rewind support. */
	public final boolean stdin;

	/*--------------------------------------------------------------*/
	/*----------------       Shared Settings        ----------------*/
	/*--------------------------------------------------------------*/

	/** Enables wrapper diagnostics; configure before concurrent use. */
	public static boolean verbose=false;

}
