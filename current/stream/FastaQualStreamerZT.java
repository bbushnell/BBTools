package stream;

import java.io.PrintStream;
import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/** FASTA plus QUAL reader that parses on the calling thread.
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * Adapted from FastaQualReadInputStream. Construction immediately opens both inputs
 * with explicit ByteFile type2 selection, which may start input workers; start is
 * a no-op. Configure backend settings before construction and instance settings
 * before the first nextList call. Always close explicitly after consumption or a
 * parsing failure: a final data batch can mark finished before it is returned,
 * causing the following terminal call to bypass nextList's close path.
 *
 * Limits and counters count constructed records before sampling; returned batches
 * contain retained records only. The first full headers are compared, while later
 * records follow FASTA order and check sequence/quality lengths. Numeric mode assumes
 * decimal integers separated by single literal spaces; it is not a general QUAL
 * validator. Headers use the default charset. Parsing failures propagate separately
 * from the backend-close status returned by errorState.
 * @author Brian Bushnell
 * @contributor Collei
 * @contributor Shinobu (documentation and formatting)
 * @date November 22, 2025
 */
public class FastaQualStreamerZT implements Streamer{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Captures configuration and immediately opens the FASTA and QUAL backends.
	 * The QUAL descriptor inherits the FASTA subprocess permission. Allowed explicit
	 * ByteFile type2 selection takes precedence over global force flags.
	 * @param ffFa FASTA input descriptor
	 * @param qf QUAL input path
	 * @param pairnum_ Side marker assigned to retained reads, normally 0 or 1
	 * @param maxReads_ Maximum records constructed before sampling; negative for unlimited
	 */
	public FastaQualStreamerZT(FileFormat ffFa, String qf, int pairnum_, long maxReads_){
		fname=ffFa.name();
		qfname=qf;
		FileFormat ffQual=FileFormat.testInput(qfname, FileFormat.QUAL, null, ffFa.allowSubprocess(), false);
		pairnum=pairnum_;
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		
		btf=ByteFile.makeByteFile(ffFa, 2);
		qtf=ByteFile.makeByteFile(ffQual, 2);
		
		// Legacy logic used a static flag; we make it instance-level here
		numericQual=true;
		
		if(verbose){outstream.println("Made FastaQualStreamerZT for "+fname);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Performs no additional startup; construction already opened both inputs. */
	@Override
	public void start(){}

	/** Parses retained reads inline and assigns a batch ID, or returns terminal null.
	 * A nonempty final batch may already have finished=true. An empty/null parse result
	 * closes both inputs, but the initial finished fast path does not; callers must
	 * close explicitly. firstRecordNum remains unspecified. Parsing failures escape
	 * without automatically setting errorState or finished; close rather than resume.
	 * @return Next nonempty retained batch, or null when terminal
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		if(finished){return null;}
		
		// Generate the list on the calling thread (Host-Driven)
		ArrayList<Read> list=toReadList(TARGET_LIST_SIZE);
		
		if(list==null || list.isEmpty()){
			finished=true;
			close();
			return null;
		}
		
		return new ListNum<Read>(list, listID++);
	}
	
	/** Folds both backend-close results and marks closed after both calls return.
	 * Repeated normal close is a no-op. Does not mark finished, reset parser state,
	 * or convert a parsing exception into an error flag. Restart is not implemented.
	 */
	@Override
	public synchronized void close(){
		if(closed){return;}
		errorState|=btf.close();//Fold the reader's error state (truncated/corrupt fasta) — the 2b reader-fold, previously MISSING here so truncation was silently accepted on the live fasta+qual path
		errorState|=qtf.close();//Fold the qual reader's error state too (two-file streamer)
		closed=true;
	}

	/** Returns whether parser completion is still pending; a returned final batch may already clear it. */
	@Override
	public synchronized boolean hasMore(){return !finished;}

	/** Returns accumulated backend-close status, independently of thrown parsing failures. */
	@Override
	public synchronized boolean errorState(){
		return errorState; // Now propagates the fasta/qual reader truncation error folded in close() (was previously never set — a 2b reader-fold gap, fixed 2026-06-22)
	}

	/** Reports false: this reader returns one unpaired record stream. */
	@Override
	public boolean paired(){return false;}

	/** Returns the configured side marker for retained reads. */
	@Override
	public int pairnum(){return pairnum;}

	/** Returns successfully constructed records before sampling, including unreturned batch data on failure. */
	@Override
	public synchronized long readsProcessed(){return generated;}

	/** Returns bases in successfully constructed reads before sampling. */
	@Override
	public synchronized long basesProcessed(){return consumedBases;}

	/** Configures random retention after construction, before the first parsing call.
	 * Rejected reads still consume the read limit and contribute to both counters.
	 * Retained numeric IDs can therefore have gaps. This is not positional sampling.
	 * @param rate Retention threshold, normally in [0,1]
	 * @param seed Seed passed to Shared.threadLocalRandom when sampling is active
	 */
	@Override
	public synchronized void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : Shared.threadLocalRandom(seed));
	}

	/** Rejects the unsupported SAM-line view.
	 * @return Never returns normally
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException();}
	
	/** Returns both configured paths as a parenthesized comma-separated name. */
	@Override
	public String fname(){return "("+fname+","+qfname+")";}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Logic          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primes initial headers and constructs reads until the retained target, limit or EOF.
	 * Initial priming discards preceding parsed content. Only the first full headers
	 * are compared; byte equality is the fast path, followed by default-charset text
	 * comparison. Construction updates base totals, then sampling and generated-count
	 * advancement occur here. EOF can mark finished while the returned list is nonempty.
	 * @param maxReadsToReturn Maximum retained records in this batch; normally positive
	 * @return Retained records, possibly empty, or null when already finished/no FASTA header
	 */
	private synchronized ArrayList<Read> toReadList(int maxReadsToReturn){
		if(finished){return null;}
		if(builder==null){builder=new ByteBuilder(2000);}
		
		if(currentHeader==null && generated==0){
			nextBases(btf, builder);
			nextQualities(qtf, builder);
			if(nextHeaderB==null){
				finish();
				return null;
			}
			//Legacy first-header validation: compare decoded full text if the byte arrays differ.
			if(!Tools.equals(nextHeaderB, nextHeaderQ)){
				String hb=new String(nextHeaderB);
				String hq=new String(nextHeaderQ);
				if(!hb.equals(hq)){
					//Descriptions are part of this initial comparison; there is no relaxed-ID fallback.
					throw new RuntimeException("Quality and Base headers differ:\n"+hb+"\n"+hq);
				}
			}
			currentHeader=nextHeaderB;
			nextHeaderB=nextHeaderQ=null;
			if(currentHeader==null){
				finish();
				return null;
			}
		}
		
		ArrayList<Read> list=new ArrayList<Read>(Data.min(1000, maxReadsToReturn));
		int added=0;
		
		while(added<maxReadsToReturn && generated<maxReads){
			Read r=makeRead(generated);
			if(r==null){
				finish();
				break;
			}
			
			if(samplerate>=1f || randy.nextFloat()<samplerate){
				r.setPairnum(pairnum);
				list.add(r);
				added++;
			}
			generated++;
		}
		
		return list;
	}
	
	/** Consumes the current FASTA record and corresponding qualities, then advances headers.
	 * Later QUAL header names are deliberately ignored. Length mismatch throws before
	 * construction/counter updates. Bases are uppercased; the default charset decodes
	 * the FASTA name. Successful construction contributes its length before sampling.
	 * @param numericID Generated-record index before sampling
	 * @return Constructed read, possibly empty; null for no header or a null base array
	 */
	private Read makeRead(long numericID){
		if(currentHeader==null){return null;}
		
		final byte[] bases=nextBases(btf, builder);
		final byte[] quals=nextQualities(qtf, builder);
		final byte[] header=currentHeader;

		//currentHeader follows the FASTA header chain; nextHeaderQ (the qual header) is intentionally
		//discarded here. fa/qual headers are compared only for record 1 (toReadList priming); records
		//2..N trust file ordering + the per-read length check below. Legacy-faithful (FastaQualReadInputStream).
		currentHeader=nextHeaderB;
		nextHeaderB=nextHeaderQ=null;
		
		if(bases==null){return null;}
		
		if(bases.length!=quals.length){
			throw new RuntimeException("\nFor sequence "+numericID+", name "+new String(header)+":\n"+
					"The bases and quality scores are different lengths, "+bases.length+" and "+quals.length);
		}
		
		for(int i=0; i<bases.length; i++){
			bases[i]=(byte)Tools.toUpperCase(bases[i]);
		}
		
		String hd=new String(header, 1, ReadHeader.end(header, 1, Shared.TRIM_READ_DESCRIPTION)-1); //Strip '>'
		Read r=new Read(bases, quals, hd, numericID);
		consumedBases+=r.length();
		return r;
	}
	
	/** Concatenates nonheader lines and saves the next FASTA header, following legacy logic.
	 * Returns an empty array for no collected content; the shared builder is cleared.
	 * @param btf Sequence input backend
	 * @param bb Empty shared accumulation buffer
	 * @return Copied sequence bytes, possibly empty
	 */
	private final byte[] nextBases(ByteFile btf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=btf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			bb.append(line);
			line=btf.nextLine();
		}
		
		if(line!=null){
			nextHeaderB=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/** Collects qualities through the next header and clears the shared builder afterward.
	 * Numeric mode assumes decimal digits with single literal-space separators and
	 * narrows accumulated ints to bytes; it does not generally validate digits/ranges.
	 * Legacy ASCII mode subtracts FASTQ.ASCII_OFFSET from each nonheader byte.
	 * @param qtf Quality input backend
	 * @param bb Empty shared accumulation buffer
	 * @return Copied quality bytes, possibly empty
	 */
	private final byte[] nextQualities(ByteFile qtf, ByteBuilder bb){
		assert(bb.length()==0);
		byte[] line=qtf.nextLine();
		while(line!=null && (line.length==0 || line[0]!=carrot)){
			if(numericQual && line.length>0){
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
				// ASCII encoded support
				for(byte b : line){bb.append((byte)(b-FASTQ.ASCII_OFFSET));}
			}
			line=qtf.nextLine();
		}
		
		if(line!=null){
			nextHeaderQ=line;
		}
		final byte[] r=bb.toBytes();
		bb.setLength(0);
		
		return r;
	}
	
	/** Marks parsing finished without closing either backend. */
	private void finish(){finished=true;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Configured FASTA path. */
	public final String fname;
	/** Configured QUAL path. */
	public final String qfname;
	/** Side marker assigned to retained reads. */
	final int pairnum;
	/** Construction limit before sampling; negative constructor values become Long.MAX_VALUE. */
	final long maxReads;
	
	/** FASTA backend opened during construction. */
	private final ByteFile btf;
	/** QUAL backend opened during construction. */
	private final ByteFile qtf;
	
	/** Successfully constructed records counted after the sampling decision. */
	private long generated=0;
	/** Bases contributed by successful Read construction, before sampling. */
	private long consumedBases=0;
	/** ID assigned to the next nonempty returned batch. */
	private long listID=0;
	
	/** Lazily allocated scratch buffer, emptied between sequence and quality parsing. */
	private ByteBuilder builder;
	/** Current FASTA header; drives record order after initial priming. */
	private byte[] currentHeader=null;
	/** Lookahead FASTA header promoted to currentHeader by makeRead. */
	private byte[] nextHeaderB=null;
	/** Lookahead QUAL header, compared initially and discarded for later records. */
	private byte[] nextHeaderQ=null;
	
	/** Parser terminal flag, which may precede publication of the last data batch. */
	private boolean finished=false;
	/** Set after both backend-close calls return; prevents repeated normal closure. */
	private boolean closed=false;
	/** Accumulated backend-close status; parsing exceptions propagate independently. */
	public boolean errorState=false;
	
	/** Retention threshold applied after each read is constructed. */
	private float samplerate=1f;
	/** PRNG assigned by setSampleRate for rates below one. */
	private shared.Random randy=null;
	
	/** True for numeric QUAL scores; false selects legacy per-byte ASCII-offset subtraction. */
	public boolean numericQual=true;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Retained-read target per batch; configure before reading. */
	public static int TARGET_LIST_SIZE=shared.Shared.bufferLen();
	/** Header marker, numeric separator and decimal-zero bytes. */
	private final byte carrot='>', space=' ', zero='0';
	
	/** Destination for optional diagnostics. */
	protected PrintStream outstream=System.err;
	/** Compile-time optional diagnostics switch. */
	public static final boolean verbose=false;

}
