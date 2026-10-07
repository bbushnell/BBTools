package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Buffered, unpaired GenBank reader that extracts sequence letters from ORIGIN sections.
 * Records receive decimal numeric names; locus and accession metadata are not retained.
 * Sequence letters are uppercased before Read construction, whose usual validation
 * and normalization settings still apply. The amino input flag is captured at construction.
 * Callers must coordinate access: synchronized batch methods do not synchronize
 * hasMore, close, external counter reads or shared configuration changes.
 * @author Brian Bushnell
 */
public class GbkReadInputStream extends ReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints the first record from the first batch using Read's text representation.
	 * Requires an input containing a record; this method does not explicitly close the reader.
	 * @param args Command-line arguments; args[0] is the input filename
	 */
	public static void main(String[] args){
		
		GbkReadInputStream fris=new GbkReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a GenBank reader from a filename, with optional subprocess support.
	 * @param fname Input GenBank filename
	 * @param allowSubprocess_ Allow subprocess decompression if needed
	 */
	public GbkReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.GBK, null, allowSubprocess_, false));
	}
	
	/**
	 * Opens the ByteFile selected for the descriptor and captures the amino input flag.
	 * Warns if the descriptor is not marked GenBank, but still uses this parser.
	 * Buffer settings are captured during instance initialization.
	 * @param ff Nonnull descriptor of the input source
	 */
	public GbkReadInputStream(FileFormat ff){
		if(verbose){System.err.println("GbkReadInputStream("+ff+")");}//#001-fix [stream/GbkReadInputStream#001]: was "FastqReadInputStream(" (copy-paste wrong class name).
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		if(!ff.gbk()){
			System.err.println("Warning: Did not find expected gbk file extension for filename "+ff.name());//#001-fix: was "fastq file extension".
		}
		bf=ByteFile.makeByteFile(ff);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks the buffered records and fills a new batch if the source is still open.
	 * Can advance input and the generated count without handing records to the caller.
	 * A closed source with no buffered records remains at EOF, including empty input.
	 * @return Whether a buffered record is available
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			if(bf.isOpen()){
				fillBuffer();
			}
			//Zero generated records is normal for empty input; repeated EOF queries stay false.
		}
		return (buffer!=null && next<buffer.size());
	}
	
	/**
	 * Hands off the next batch, filling it as needed and incrementing consumed.
	 * Returned lists and sequence arrays are not reused by the reader. Batch capacity
	 * limits records rather than sequence bytes; this reader provides blockwise access.
	 * @return Nonempty list of unpaired reads, or null when no records remain
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
	
	/** Fills one record-limited batch, advances numeric IDs and closes input on a short batch. */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		buffer=toReadList(bf, BUF_LEN, nextReadID, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		//TODO: Review [stream/GbkReadInputStream#003] - the short-batch close result is ignored; final status propagation needs separate review before changing lifecycle behavior.
		if(bsize<BUF_LEN){bf.close();}
		
		generated+=bsize;
		//defensive/unreachable: toReadList always returns a (possibly empty) list, never null. An empty buffer -> nextList() maps empty->null as the EOF signal.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in GbkReadInputStream.");//#001-fix: was "...FastqReadInputStream." (copy-paste wrong class name).
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses records from the current position without closing or resetting the source.
	 * Each line beginning with ORIGIN starts a sequence. Sequence lines must be nonempty;
	 * their ASCII letters contribute until a line begins with slash or EOF is reached.
	 * A present terminating line is expected to begin with //. Letters are uppercased
	 * and passed to a Read with no qualities and a decimal numeric name; other metadata
	 * is ignored. Sequence scratch state is local to this call, and Read constructor
	 * validation/normalization still applies under its usual settings.
	 * @param bf Nonnull source positioned at or before the next record
	 * @param maxReadsToReturn Positive maximum number of records to return
	 * @param numericID Starting numeric ID, incremented once per constructed record
	 * @param flag Read flags passed to the constructor, such as the amino input mask
	 * @return Newly allocated, possibly empty list; each record receives a separate base array
	 */
	public static ArrayList<Read> toReadList(final ByteFile bf, final int maxReadsToReturn, long numericID, final int flag){
		ArrayList<Read> list=new ArrayList<Read>(Data.min(8192, maxReadsToReturn));
		
		int added=0;
		
		String idLine=null;
		ByteBuilder bb=new ByteBuilder();
		for(byte[] s=bf.nextLine(); s!=null; s=bf.nextLine()){
			if(Tools.startsWith(s, "ORIGIN")){
				byte[] line=null;
				for(line=bf.nextLine(); line!=null && line[0]!='/'; line=bf.nextLine()){
					for(byte b : line){
						if(Tools.isLetter(b)){
							bb.append(Tools.toUpperCase(b));
						}
					}
				}
				assert(line==null || Tools.startsWith(line, "//")) : new String(line);

				//#002 [stream/GbkReadInputStream#002] Retained naming limitation: idLine stays null, so reads use numeric names and omit locus/accession metadata. Adding names needs separate format review. The ORIGIN loop assumes nonempty lines.
				Read r=new Read(bb.toBytes(), null, idLine==null ? ""+numericID : idLine, numericID, flag);
				list.add(r);
				added++;
				numericID++;
				
				bb.clear();
				idLine=null;
				
				if(added>=maxReadsToReturn){break;}
			}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Lifecycle          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Closes the ByteFile and retains its returned status in the reader's local flag.
	 * Does not discard an existing buffered batch or merge the global FASTQ flag.
	 * @return Accumulated local error status after this close call
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+bf.name()+"; errorState="+errorState);}
		errorState|=bf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+bf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/**
	 * Resets counters, numeric IDs and buffer, then delegates reopening/rewinding to ByteFile.
	 * Retains the local error flag and captured format/flag/buffer settings. Whether the
	 * underlying source supports rereading belongs to the selected ByteFile implementation.
	 */
	@Override
	public synchronized void restart(){
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		bf.reset();
	}

	/** Indicates whether reads are paired; always false for GenBank input.
	 * @return false */
	@Override
	public boolean paired(){return false;}
	
	/** Returns the name of the underlying input file.
	 * @return Input filename */
	@Override
	public String fname(){return bf.name();}
	
	/**
	 * Combines the local flag with FASTQ's process-wide flag, even though this reader
	 * parses GenBank directly. Does not query the ByteFile here and may differ from close().
	 * @return Local error status or the current global FASTQ status
	 */
	@Override
	//GenBank parses directly in toReadList, but this method also includes FASTQ's process-global status; behavior retained.
	public boolean errorState(){return errorState || FASTQ.errorState();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Generated batch retained until the next handoff; cleared on restart. */
	private ArrayList<Read> buffer=null;
	/** Legacy block cursor; current batch-only access keeps it at zero. */
	private int next=0;//vestigial: never incremented (blockwise-only reader via nextList) -> always 0.
	
	/** Selected byte reader, retained across restart. */
	private final ByteFile bf;
	/** Amino input flag captured at construction and passed to every Read. */
	private final int flag;

	/** Maximum records per generated batch, captured at instance initialization. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Captured legacy byte-budget setting; currently unused by the GenBank parser. */
	private final long MAX_DATA=Shared.bufferData();//TODO: A per-batch byte limit is not implemented; batches are capped only by record count.

	/** Individual records placed in generated batches during the current pass. */
	public long generated=0;
	/** Individual records handed to the caller during the current pass. */
	public long consumed=0;
	/** Numeric ID assigned to the first record of the next generated batch. */
	private long nextReadID=0;
	
	/** Whether the input descriptor denotes a standard stream. */
	public final boolean stdin;
	/** Enables construction and close diagnostics on standard error. */
	public static boolean verbose=false;

}
