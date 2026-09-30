package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * Reads unpaired nucleotide sequences from FASTA with a separate QUAL file.
 * Construction opens both ByteFiles and captures the shared batch length. Parsing
 * runs on the calling thread; the selected ByteFile implementation may use workers.
 * Configure shared input options before use and coordinate iteration, close and restart.
 *
 * Only the initial headers are asserted equal. Later records follow file order, with
 * asserted sequence/quality length equality and subsequent Read validation.
 * Numeric QUAL parsing assumes its historical decimal/space grammar; this reader
 * does not validate all header identities, token syntax or trailing QUAL records.
 * Names use the default charset and sequence bytes are uppercased before Read creation.
 *
 * Iteration is blockwise through {@link #nextList()}; {@link #hasMore()} may prefetch.
 * Generated and consumed counts refer to read entries, not batches. Short fills
 * close the inputs before returning the buffered reads; callers should still close
 * explicitly on early termination or failure. Parsing exceptions are propagated
 * separately from the cached {@link #errorState()} flag.
 *
 * @contributor Shinobu (EOF handling and documentation)
 */
public class FastaQualReadInputStream extends ReadInputStream{
	
	/**
	 * Prints the first read from up to five batches, then closes both inputs.
	 * @param args Sequence filename followed by quality filename
	 * @throws RuntimeException If normal processing completes and input closure reports an error
	 */
	public static void main(final String[] args){
		final FastaQualReadInputStream fris=new FastaQualReadInputStream(args[0], args[1], true);
		boolean error=false;
		try{
			//STR-056: Test EOF before indexing and stop before fetching a sixth sample batch.
			for(int i=0; i<5; i++){
				final ArrayList<Read> list=fris.nextList();
				if(list==null){break;}
				System.out.println(list.get(0).toText(false));
			}
		}finally{
			error=fris.close();
		}
		if(error){throw new RuntimeException("Error closing FASTA/QUAL inputs.");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves the sequence format and immediately opens both inputs.
	 * @param fname Sequence input name; FASTA is the default format
	 * @param qfname Required quality input name; QUAL is the default format
	 * @param allowSubprocess_ Whether both inputs may use subprocess decompression
	 */
	public FastaQualReadInputStream(String fname, String qfname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false), qfname);
	}
	
	/**
	 * Opens automatically selected ByteFiles for sequences and qualities.
	 * Shared buffer length and backend preferences must be set before construction.
	 * A non-FASTA, non-stdio sequence descriptor produces a warning, not conversion.
	 * @param ff Sequence input descriptor, also supplying the subprocess policy
	 * @param qfname Required quality input name
	 */
	public FastaQualReadInputStream(FileFormat ff, String qfname){
		if(!ff.fasta() && !ff.stdio()){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+ff.name());
		}
		
		btf=ByteFile.makeByteFile(ff);
		qtf=ByteFile.makeByteFile(FileFormat.testInput(qfname, FileFormat.QUAL, null, ff.allowSubprocess(), false));
		//QUAL records supply scores, not mates. The private pairing branch is unreachable.
		interleaved=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Iteration           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests for buffered reads, filling a batch if needed while the sequence input is open.
	 * Prefetch may increase generated without increasing consumed. Coordinate this
	 * unsynchronized method with lifecycle changes; it is not a multicaller iterator.
	 * @return Whether a buffered read remains, including a batch retained after close
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			//STR-055: Zero generated reads is normal for empty input; closed readers stay at EOF.
			if(btf.isOpen()){fillBuffer();}
		}
		return (buffer!=null && next<buffer.size());
	}
	
	/**
	 * Transfers the next nonempty batch to the caller, or returns null at EOF.
	 * Returned lists and their copied sequence/quality arrays survive subsequent fills.
	 * If parsing fails while building a batch, its partial list is neither returned nor counted.
	 * @return Next batch of unpaired reads, or null when no records remain
	 */
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
		return r;
	}

	/** Fills one captured-size batch; short fills close inputs and fold their error flags. */
	private synchronized void fillBuffer(){
		if(builder==null){builder=new ByteBuilder(2000);}
		if(verbose){System.err.println("Filling buffer.  buffer="+(buffer==null ? null : buffer.size()));}
		assert(buffer==null || next>=buffer.size());
		buffer=null;
		next=0;
		buffer=toReads(BUF_LEN, nextReadID, interleaved);
		final int count=(buffer==null ? 0 : buffer.size());

		if(verbose){System.err.println("Filled buffer.  size="+count);}
		
		nextReadID+=count;
		if(count<BUF_LEN){
			if(verbose){System.err.println("Closing tf");}
			errorState|=close();
		}
		generated+=count;
		if(verbose){System.err.println("generated="+generated);}
	}
	
	/**
	 * Builds a list and checks the parser's terminal and size invariants.
	 * @param maxReadsToReturn Maximum list entries
	 * @param numericID ID assigned to the first entry
	 * @param interleaved Private pairing switch; always false for this reader
	 * @return Parsed list, possibly empty, or null after termination
	 */
	private ArrayList<Read> toReads(int maxReadsToReturn, long numericID, boolean interleaved){
		ArrayList<Read> list=toReadList(maxReadsToReturn, numericID, interleaved);
		if(list==null){assert(finished);}else{assert(list.size()<=maxReadsToReturn);}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Primes initial headers and constructs a bounded list in FASTA record order.
	 * Only the first full headers are asserted equal; later QUAL headers are discarded.
	 * FASTA termination drives iteration, without checking for extra QUAL records.
	 * @param maxReadsToReturn Maximum list entries
	 * @param numericID ID assigned to the first entry
	 * @param interleaved Private pairing switch; unreachable as true through this class
	 * @return List of entries, possibly empty, or null if finished or no initial FASTA header
	 */
	private ArrayList<Read> toReadList(int maxReadsToReturn, long numericID, boolean interleaved){
		if(finished){return null;}
		if(verbose){System.err.println("FastaQualRIS fetching a list.");}
		
		if(currentHeader==null && numericID==0){
			//Prime both headers. Copies preceding them are normally empty and are discarded;
			//leading non-header content is not retained. The side effects save the headers.
			nextBases(btf, builder);
			nextQualities(qtf, builder);
			if(nextHeaderB==null){
				finish();
				return null;
			}
			assert(Tools.equals(nextHeaderB, nextHeaderQ)) : "Quality and Base headers differ for read "+numericID;
			currentHeader=nextHeaderB;
			nextHeaderB=nextHeaderQ=null;
			if(currentHeader==null){
				finish();
				return null;
			}
		}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(1000, maxReadsToReturn));
		
		int added=0;
		Read prev=null;
		
		while(added<maxReadsToReturn){
			Read r=makeRead(numericID);
			if(verbose){System.err.println("Made "+r);}
			if(r==null){
				finish();
				if(verbose){System.err.println("makeRead returned null.");}
				break;
			}
			if(interleaved){
				if(prev==null){prev=r;}else{
					prev.mate=r;
					r.mate=prev;
					list.add(prev);
					added++;
					numericID++;
					prev=null;
				}
			}else{
				list.add(r);
				added++;
				numericID++;
			}
		}
		
		assert(list.size()<=maxReadsToReturn);
		if(verbose){System.err.println("FastaQualRIS returning a list.  Size="+list.size());}
		return list;
	}
	
	/**
	 * Concatenates sequence lines up to the next header or EOF and saves header lookahead.
	 * Case normalization occurs later in makeRead. EOF leaves lookahead unchanged.
	 * @param btf Sequence ByteFile, also used to find the initial header during priming
	 * @param bb Empty reusable builder, cleared before returning
	 * @return Newly copied bytes, possibly empty; never null
	 */
	private final byte[] nextBases(ByteFile btf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=btf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			bb.append(line);
			line=btf.nextLine();
		}
		if(line==null){}else{
			assert(line.length>0);
			assert(line[0]==carrot);
			nextHeaderB=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/**
	 * Accumulates scores through the next QUAL header or EOF and saves header lookahead.
	 * Numeric mode expects unsigned decimal values separated by single literal spaces,
	 * without leading/trailing spaces; token syntax and integer/byte ranges are not
	 * generally checked here. ASCII mode subtracts FASTQ.ASCII_OFFSET from each byte.
	 * @param qtf Quality ByteFile, also used to find the initial header during priming
	 * @param bb Empty reusable builder, cleared before returning
	 * @return Newly copied score bytes, possibly empty; never null
	 */
	private final byte[] nextQualities(ByteFile qtf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=qtf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			if(NUMERIC_QUAL && line.length>0){
				int x=0;
				for(int i=0; i<line.length; i++){
					byte b=line[i];
					if(b==space){
						assert(i>0);
						bb.append((byte)x);
						x=0;
					}else{
						x=10*x+(b-zero);
					}
				}
				bb.append((byte)x);
			}else{
				for(byte b : line){bb.append((byte)(b-FASTQ.ASCII_OFFSET));}
			}
			line=qtf.nextLine();
		}
		if(line==null){}else{
			assert(line.length>0);
			assert(line[0]==carrot);
			nextHeaderQ=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/**
	 * Reads one sequence and its scores, advances header state and constructs a Read.
	 * Length equality asserts before Read validation. With assertions enabled, empty
	 * records reach the first-base assertion operand and fail before Read construction.
	 * @param numericID ID to assign without incrementing the reader counters
	 * @return New unpaired read, or null if finished or no current FASTA header
	 */
	private Read makeRead(long numericID){
		if(finished){
			if(verbose){System.err.println("Returning null because finished.");}
			return null;
		}
		if(currentHeader==null){return null;}
		assert(nextHeaderB==null);
		assert(nextHeaderQ==null);
		
		final byte[] bases=nextBases(btf, builder);
		final byte[] quals=nextQualities(qtf, builder);
		final byte[] header=currentHeader;
		
		currentHeader=nextHeaderB;
		nextHeaderB=nextHeaderQ=null;
		
		//Defensive only: nextBases returns a copy, never null. Normal termination uses
		//finished/currentHeader above after the last header has been consumed.
		if(bases==null){
			if(verbose){System.err.println("Returning null because tf.nextLine()==null: A");}
			return null;
		}
		
		assert(bases.length==quals.length) :
			"\nFor sequence "+numericID+", name "+new String(header)+":\n"+
					"The bases and quality scores are different lengths, "+bases.length+" and "+quals.length;
		
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Tools.toUpperCase(bases[i]);
		}
		assert(bases[0]!=carrot) : new String(bases)+"\n"+numericID+"\n"+header[0];
		String hd=new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1);
		Read r=new Read(bases, quals, hd, numericID);
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Lifecycle           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Closes both inputs on the first call, finishes parsing and releases the builder.
	 * Already buffered reads remain available. The first return is not itself cached
	 * in errorState; fillBuffer folds it on short fills. Later calls return the cached
	 * reader flag, which may differ from an earlier direct close result.
	 * @return Backend error flags on first close, or cached reader errorState if already closed
	 */
	@Override
	public synchronized boolean close(){
		if(closed){return errorState;}
		if(verbose){System.err.println("FastaQualRIS closing.");}
		builder=null;
		finish();
		boolean a=btf.close();
		boolean b=qtf.close();
		closed=true;
		//Report backend flags rather than a literal false (historical FastaReadInputStream#001).
		//The old equivalence to errorState|a|b assumes a fresh lifetime in which fillBuffer
		//first sets errorState from this return. Restart retains errorState, so that is conditional.
		return a|b;
	}

	/**
	 * Resets reader counters, buffers and headers, then delegates reset to both ByteFiles.
	 * Replay depends on the backends and inputs. Captured batch size and ByteFile
	 * instances are retained; inherited errorState is not cleared.
	 */
	@Override
	public synchronized void restart(){
		if(verbose){System.err.println("FastaQualRIS restarting.");}
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		
		finished=false;
		closed=false;

		buffer=null;
		nextHeaderB=null;
		nextHeaderQ=null;
		currentHeader=null;
		builder=null;

		btf.reset();
		qtf.reset();
	}

	/**
	 * Reports this reader's fixed unpaired mode.
	 * @return false
	 */
	@Override
	public boolean paired(){return interleaved;}
	
	/** Stops subsequent parsing without closing backends or discarding buffered reads. */
	private synchronized void finish(){
		if(verbose){System.err.println("FastaQualRIS setting finished "+finished+" -> "+true);}
		finished=true;
	}
	
	/**
	 * Returns a diagnostic label for both inputs.
	 * @return Sequence name, comma, quality name
	 */
	@Override
	public String fname(){return btf.name()+","+qtf.name();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Prefetched batch, relinquished when nextList transfers it to the caller. */
	private ArrayList<Read> buffer=null;
	/** Vestigial element index; blockwise access never increments it. */
	private int next=0;

	/** Sequence input, opened during construction. */
	private final ByteFile btf;
	/** Quality input, opened during construction using the sequence subprocess policy. */
	private final ByteFile qtf;
	/** Always false; QUAL supplies scores rather than a second mate. */
	private final boolean interleaved;

	/** Batch capacity captured from Shared.bufferLen before constructor input opening. */
	private final int BUF_LEN=Shared.bufferLen();

	/** Entries successfully buffered since construction or restart, including prefetched entries. */
	public long generated=0;
	/** Entries handed to callers since construction or restart; this is not a batch count. */
	public long consumed=0;
	/** Numeric ID for the first entry of the next successfully filled batch. */
	private long nextReadID=0;
	
	/** Shared parser mode: true for decimal scores, false for ASCII scores; keep stable during parsing. */
	public static boolean NUMERIC_QUAL=true;
	
	/** Enables reader diagnostics on standard error. */
	public static boolean verbose=false;

	/** Next FASTA header found while copying the current sequence. */
	private byte[] nextHeaderB=null;
	/** Next QUAL header; compared during initial priming, discarded for later records. */
	private byte[] nextHeaderQ=null;
	
	/** FASTA header for the next Read, including the leading marker. */
	private byte[] currentHeader=null;
	
	/** Lazily allocated scratch shared by sequence and quality parsing, each returning a copy. */
	private ByteBuilder builder=null;
	
	/** Gates further parsing after finish or close, until restart. */
	private boolean finished=false;
	/** Whether this reader's backend-closing branch has completed, until restart. */
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** FASTA/QUAL header marker. */
	private final byte carrot='>';
	/** Numeric QUAL token separator. */
	private final byte space=' ';
	/** Decimal digit offset used by numeric QUAL parsing. */
	private final byte zero='0';
	
}
