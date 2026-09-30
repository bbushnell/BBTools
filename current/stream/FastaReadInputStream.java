package stream;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import dna.Gene;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/** Legacy buffered FASTA reader with optional section splitting and adjacent-read pairing.
 * Batches contain unpaired reads or first mates linked to second mates. Length limits,
 * pairing and batch limits are captured at construction; quality, header and warning
 * options remain live globals. Configure these before use, not during concurrent reads.
 * A single consumer owns parsing and restart. Serialized closure does not make the
 * parser generally thread-safe. Read errors handled while open are reported through
 * {@link #errorState()} and {@link #close()}; callers must check the final status even
 * when reads were returned. Errors handled after intentional closure remain tolerated.
 * This parser retains historical tolerances and is not a strict FASTA validator.
 *
 * @author Brian Bushnell, Shinobu
 * @date Feb 13, 2013
 */
public class FastaReadInputStream extends ReadInputStream{
	
	/** Legacy timing/demo driver; consumes the first read of each fetched batch.
	 * Assumes enough nonempty batches for its requested counts; not a validation CLI.
	 * @param args Filename, optional display/count limits, minimum length and section length
	 */
	public static void main(final String[] args){
		
		int a=20, b=Integer.MAX_VALUE;
		if(args.length>1){a=Integer.parseInt(args[1]);}
		if(args.length>2){b=Integer.parseInt(args[2]);}
		if(args.length>3){MIN_READ_LEN=Integer.parseInt(args[3]);}
		if(args.length>4){TARGET_READ_LEN=Integer.parseInt(args[4]);}
		if(TARGET_READ_LEN<1){
			TARGET_READ_LEN=Integer.MAX_VALUE;
			SPLIT_READS=false;
		}else{
			SPLIT_READS=true;
		}
		
		Timer t=new Timer();
		
		FastaReadInputStream fris=new FastaReadInputStream(args[0], false, false, false, Shared.bufferData());
		Read r=fris.nextList().get(0);
		int i=0;
		
		while(r!=null){
			if(i<a){System.out.println("'"+r.toText(false)+"'");}
			r=fris.nextList().get(0);
			if(++i>=a){break;}
		}
		while(r!=null && i++<b){r=fris.nextList().get(0);}
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens a filename through a FASTA descriptor without content detection.
	 * @param fname Input path or stdin name
	 * @param interleaved_ Link adjacent generated reads as mates
	 * @param amino_ Mark generated reads as amino-acid sequences
	 * @param allowSubprocess_ Permit external decompression when opening the stream
	 * @param maxdata Positive soft base limit per batch; otherwise use Shared.bufferData()
	 */
	public FastaReadInputStream(final String fname, final boolean interleaved_, final boolean amino_, final boolean allowSubprocess_, final long maxdata){
		this(FileFormat.testInput(fname, FileFormat.FASTA, FileFormat.FASTA, 0, allowSubprocess_, false, false), interleaved_, amino_, maxdata);
	}
	
	/** Opens the descriptor's name and captures length, pairing and batch limits.
	 * Format content is not validated here; settingsOK is called through an assertion
	 * after opening. The descriptor supplies subprocess permission, not pairing.
	 * @param ff Nonnull input descriptor
	 * @param interleaved_ Link adjacent generated reads as mates
	 * @param amino_ Mark generated reads as amino-acid sequences
	 * @param maxdata Positive soft base limit per batch; otherwise use Shared.bufferData()
	 */
	public FastaReadInputStream(final FileFormat ff, final boolean interleaved_, final boolean amino_, final long maxdata){
		name=ff.name();
		amino=amino_;
		flag=(amino ? Read.AAMASK : 0);
		
		if(!fileIO.FileFormat.hasFastaExtension(name) && !name.startsWith("stdin")){
			System.err.println("Warning: Did not find expected fasta file extension for filename "+name);
		}
		
		interleaved=interleaved_;
		allowSubprocess=ff.allowSubprocess();
		minLen=MIN_READ_LEN;
		maxLen=(SPLIT_READS ? (TARGET_READ_LEN>0 ? TARGET_READ_LEN : Integer.MAX_VALUE) : Integer.MAX_VALUE);
		MAX_DATA=maxdata>0 ? maxdata : Shared.bufferData();
		
		ins=open();
		
		assert(settingsOK());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Transfers the prefetched or next batch to the caller; returned lists are not reused.
	 * Each entry is an unpaired read or first mate; a final unmatched first mate is retained.
	 * @return Next nonempty batch, or null when no reads remain
	 * @throws RuntimeException If the internal single-read index indicates mixed access
	 */
	@Override
	public ArrayList<Read> nextList(){
		if(nextReadIndex!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(currentList==null || nextReadIndex>=currentList.size()){
			boolean b=fillList();
		}
		ArrayList<Read> list=currentList;
		currentList=null;
		if(list==null || list.isEmpty()){
			list=null;
		}else{
			consumed+=list.size();
		}
		return list;
	}
	
	/** Prefetches a batch if needed and reports whether an entry remains.
	 * A false result does not guarantee closure, notably for an empty input.
	 * @return true if a prefetched batch contains an unread entry
	 */
	@Override
	public boolean hasMore(){
		if(currentList==null || nextReadIndex>=currentList.size()){
			if(open){
				fillList();
			}else{
//				assert(generated>0) : "Was the file empty?";
			}
		}
		return (currentList!=null && nextReadIndex<currentList.size());
	}
	
	/** Closes any existing input, resets parsing/counters and reopens the same name.
	 * Keeps captured settings, allocated buffer, reported-header reference and sticky
	 * error state. Reopening stdin does not rewind it. Requires exclusive parser ownership.
	 */
	@Override
	public void restart(){
		if(ins!=null){close();}
		assert(ins==null);
//		generated=0;
		consumed=0;
		nextReadIndex=0;
		nextReadID=0;
		currentList=null;
		header=null;
		bstart=0;
		bstop=0;
		
		currentSection=0;
		if(ins==null){
			ins=open();
		}else{
			assert(false) : "is should be null";
		}
	}
	
	/** Serializes closure on this reader and folds reported read/close errors.
	 * Marks the reader closed before releasing the underlying input; repeated calls
	 * return the sticky error state. A bare System.in is left open; wrappers follow
	 * their own close behavior and may close stdin. This monitor also
	 * coordinates the read-exception handler, not concurrent parsing or restart.
	 * @return true if a read or close error has been reported
	 */
	@Override
	public final boolean close(){
		synchronized(this){
			if(!open){return errorState;}
			open=false;
			assert(ins!=null);

			try{
				if(ins!=System.in){
					errorState|=ReadWrite.finishReading(ins, name, allowSubprocess);
				}
			}catch(final Exception e){
				System.err.println("Some error occured: "+e);
				errorState=true;
			}

			ins=null;
		}
		//#001-fix [stream/FastaReadInputStream#001]: this formerly always returned false.
		//ConcurrentGenericReadInputStream folds producer.close() into its error field, but its
		//errorState() override also checks retained producers. Standard ReadWrite.closeStreams
		//checked that override after closure and therefore still detected the error; a caller
		//using only the close return could miss it. Returning the sticky flag fixes that latent
		//inconsistency. Prior review corrected an initial overstatement of standard-path loss.
		return errorState;
	}
	
	/** Returns the captured pairing mode, not a guarantee that every read has a mate. */
	@Override
	public boolean paired(){return interleaved;}
	
	/** Allocates a batch bounded by entry count and a soft total-base threshold.
	 * Both mates count toward bases but only the first occupies a list slot. A whole
	 * read/pair may cross MAX_DATA. An unmatched final first mate remains in the list;
	 * its nextReadID increment is skipped when the missing-mate branch breaks.
	 * @return true if at least one list entry was loaded
	 */
	private final boolean fillList(){
//		assert(open);
		if(!open){
			currentList=null;
			return false;
		}
		assert(currentList==null || nextReadIndex>=currentList.size());
		nextReadIndex=0;
		currentList=new ArrayList<Read>(BUF_LEN);
		
		if(header==null){
			header=nextHeader();
			if(header==null){
				currentList=null;
				return false;
			}
		}
		long len=0;
		for(int i=0; i<BUF_LEN && len<MAX_DATA; i++){
			Read r=generateRead(0);
			if(r==null){break;}
			currentList.add(r);
			len+=r.length();
			if(interleaved){
				//r is already in the list: retain an unmatched final first mate. This is historical
				//FASTA tolerance, not proof of well-formed interleaving; SCARF uses a separate policy.
				Read r2=generateRead(1);
				if(r2==null){break;}
				len+=r2.length();
				r.mate=r2;
				r2.mate=r;
			}
			nextReadID++;
			if(verbose){System.err.println("Generated a read; i="+i+", BUF_LEN="+BUF_LEN);}
//			if(i==1){assert(false) : r.numericID+", "+r.mate.numericID;}
		}
		
		return currentList.size()>0;
	}
	
	/** Generates one accepted sequence or section, skipping sections below minLen.
	 * Reads current quality/header globals and retains the obsolete custom-header branch.
	 * @param pairnum Pair flag for this read (0 for first/unpaired, 1 for second)
	 * @return Generated read with the current entry ID, or null after end-of-input closure
	 */
	private final Read generateRead(final int pairnum){
		if(verbose){System.err.println("Called generateRead(); bstart="+bstart+", bstop="+bstop+", currentSection="+currentSection+", header="+header);}
		assert(header!=null) : "Null header for fasta read - input file may be corrupt: "+name;
		if(bstart<bstop && buffer[bstart]==carrot){
			header=nextHeader();
			currentSection=0;
		}
		if(header==null){
			close();
			return null;
		}
		byte[] bases=nextBases();
		
		currentSection++;
		while(bases==null){
			header=nextHeader();
			if(header==null){
				close();
				return null;
			}
			bases=nextBases();
			currentSection=1;
		}
		assert(bases!=null);
		assert(bases.length>0);
		
		byte[] quals=null;
		if(FAKE_QUALITY){
			quals=new byte[bases.length];
			Arrays.fill(quals, (byte)(Shared.FAKE_QUAL));
		}
		String hd=((!FORCE_SECTION_NAME && currentSection==1 && bases.length<=maxLen) ? header : header+"_part_"+currentSection);
		assert(currentSection==1 || bases.length>0) : "id="+hd+", section="+currentSection+", len="+bases.length+"\n"+new String(bases);
		Read r=null;
		if(FASTQ.PARSE_CUSTOM){//TODO: This appears to be legacy code from old custom header format, which won't work anymore
			assert(false) : "TODO - this code for parsing custom headers of fasta files is obsolete.";
			if(header!=null && header.indexOf('_')>0){
				String temp=header;
				if(temp.endsWith(" /1") || temp.endsWith(" /2")){temp=temp.substring(0, temp.length()-3);}
				String[] answer=temp.split("_");

				if(answer.length>=5){
					try{
						int trueChrom=Gene.toChromosome(answer[1]);
						byte trueStrand=Byte.parseByte(answer[2]);
						int trueLoc=Integer.parseInt(answer[3]);
						int trueStop=Integer.parseInt(answer[4]);
						r=new Read(bases, quals, hd, nextReadID, (flag|trueStrand), trueChrom, trueLoc, trueStop);
						r.setSynthetic(true);
					}catch(final NumberFormatException e){
						FASTQ.PARSE_CUSTOM=false;
						System.err.println("Turned off PARSE_CUSTOM because could not parse "+new String(header));
					}
				}else{
					FASTQ.PARSE_CUSTOM=false;
					System.err.println("Turned off PARSE_CUSTOM because answer="+Arrays.toString(answer));
				}
			}else{
				FASTQ.PARSE_CUSTOM=false;
				System.err.println("Turned off PARSE_CUSTOM because header="+header+", index="+header.indexOf('_'));
			}
		}
		if(r==null){
			r=new Read(bases, quals, hd, nextReadID, flag);
		}
		r.setPairnum(pairnum);
		if(verbose){System.err.println("Made read:\t"+(r.length()>1000 ? r.id : r.toString()));}
		return r;
	}
	
	/** Advances past a header and preserves its raw bytes through a Latin1 String.
	 * Honors live description trimming and legacy SOH/tab header handling. Refill
	 * expects the first non-comment record to begin with a header marker; leading
	 * blank lines are not generally accepted. See the byte-roundtrip rationale below.
	 * @return Header without the leading marker, possibly empty, or null when unavailable
	 */
	private String nextHeader(){
		if(verbose){System.err.println("Called nextHeader(); bstart="+bstart+"; bstop="+bstop);}
		assert(bstart>=bstop || buffer[bstart]=='>' || buffer[bstart]<=slashr) : bstart+", "+bstop+", '"+(char)buffer[bstart]+"'"+"\t"+name;
		while(bstart<bstop && buffer[bstart]!='>'){bstart++;}
		int x=bstart;
		assert(bstart>=bstop || buffer[x]=='>') : bstart+", "+bstop+", '"+(char)buffer[x]+"'";
		while(x<bstop && buffer[x]>slashr){x++;}
		if(verbose){System.err.println("A: x="+x);}
		if(x<bstop && (buffer[x]<=STX || buffer[x]==tab)){ //Handle deprecated 'SOH' symbol and tab
			while(x<bstop && (buffer[x]>slashr || buffer[x]<=STX || buffer[x]==tab)){
				if(buffer[x]==0x1){buffer[x]=carrot;}
				x++;
			}
		}
		if(verbose){System.err.println("B: x="+x);}
		if(x>=bstop){
			int fb=fillBuffer();
			if(verbose){System.err.println("B: fb="+fb);}
			if(fb<1){
				if(verbose){System.err.println("Returning null from nextHeader()");}
				return null;
			}
			//STR-016: fillBuffer may skip leading comments without shifting their bytes away.
			//Resume at its unread offset so comments cannot become a header or sequence.
			x=bstart;
			assert(bstart>=0 && bstart<bstop && buffer[x]=='>') : "Improperly formatted fasta file; expecting '>' symbol.\n"+
				(buffer[x]=='@' ? "If this is a fastq file, please rename it with a '.fastq' extension.\n" : "")+
					bstart+", "+bstop+", "+(int)buffer[x]+", "+(char)buffer[x]; //Note: This assertion will fire if a fasta file starts with a newline.
			while(x<bstop && buffer[x]>slashr){x++;}
			if(x<bstop && (buffer[x]<=STX || buffer[x]==tab)){ //Handle deprecated 'SOH' symbol and tab
				while(x<bstop && (buffer[x]>slashr || buffer[x]<=STX || buffer[x]==tab)){
					if(buffer[x]==0x1){buffer[x]=carrot;}
					x++;
				}
			}
		}
		if(verbose){System.err.println("C: x="+x);}
		assert(x>=bstop || buffer[x]<=slashr);
		
		int start=bstart+1, stop=x;
		if(Shared.TRIM_READ_DESCRIPTION){
			for(int i=start; i<stop; i++){
				if(Character.isWhitespace(buffer[i])){
					stop=i;
					break;
				}
			}
		}
		
		//FIXED [stream/FastaReadInputStream UTF8]: was US_ASCII, which maps every byte>127 to U+FFFD and destroys
		//UTF-8 (or any 8-bit) header bytes at READ time, irrecoverably. ISO-8859-1 is a bijective byte<->char map
		//(0-255), so raw header bytes survive read->store->write losslessly through the existing (byte)charAt writer
		//and the file reproduces byte-for-byte. Do NOT switch to UTF-8 decode here: that yields multibyte chars the
		//(byte)charAt writer would re-truncate. (Faithful passthrough; display-correct UTF-8 would be a design change.)
		String s=stop>start ? new String(buffer, start, stop-start, StandardCharsets.ISO_8859_1) : "";
		if(verbose){System.err.println("Fetched header: '"+s+"'");}
		bstart=x+1;
		
		return s;
	}
	
	/** Extracts up to maxLen sequence bytes, omitting signed byte values at or below carriage return.
	 * Refills/grows the shared byte buffer as needed and advances bstart past consumed
	 * input. Further sequence normalization belongs to Read construction. Detected empty
	 * first sections may invoke handleNoSequence; an initial refill with no new bytes
	 * returns null directly, so not every empty section reaches that policy.
	 * @return Fresh sequence array, or null for exhaustion or a section shorter than minLen
	 */
	private byte[] nextBases(){
		if(verbose){System.err.println("Called nextBases()");}
		if(bstart>=bstop){
			int bytes=fillBuffer();
			if(bytes<1 || !open){return null;}
		}
		int x=bstart;
		int bases=0;
		
		if(!(x>=bstop || buffer[x]!='>')){
			handleNoSequence(x);
		}
		
		while(x<bstop && bases<maxLen && buffer[x]!='>'){
			while(x<bstop && bases<maxLen && buffer[x]!='>'){
				if(buffer[x]>slashr){bases++;}
				x++;
			}
			assert(x==bstop || buffer[x]=='>' || bases==maxLen);
			if(x==bstop && bases<maxLen){
				int fb=fillBuffer();
				if(fb<1){
					x=bstop;
					if(verbose){System.err.println("Broke loop when fb="+fb+"; bstart="+bstart+", bstop="+bstop);}
					break;
				}
				//Recount after fillBuffer shifted unconsumed bytes to the front and set bstart=0.
				//Previously counted bases now start at buffer[0]; rescanning is intentional, not a leak.
				//Doubling the buffer keeps this amortized O(n) as it grows to fit the record.
				x=bstart;
				bases=0;
			}
		}
		
		if(bases<minLen){
			
			if(bases==0){handleNoSequence(x);}
			
			bstart=x;
			if(verbose){System.err.println("Fetched "+bases+" bases; returning null.  bstart="+bstart+", bstop="+bstop/*+"\n"+new String(buffer)*/);}
			return null;
		}
		
		byte[] r=new byte[bases];
		
		for(int i=bstart, j=0; j<bases; i++){
			assert(i<x);
			byte b=buffer[i];
			if(b>slashr){
				r[j]=b;
				j++;
			}
		}
		
		if(verbose){System.err.println("Fetched "+bases+" bases, open="+open+":\n'"+(r.length>1000 ? "*LONG*" : new String(r))+"'");}
		
		bstart=x;
		return r;
	}
	
	/** Applies the live missing-sequence warning/assertion policy to a first section.
	 * Later sections are ignored because an exact split-length multiple can enter here.
	 * Warning suppression compares header references and may disable warnings globally.
	 * When warnings are enabled, a repeated header reference returns before the abort
	 * assertion too. Otherwise ABORT_IF_NO_SEQUENCE is enforced only with assertions enabled.
	 * @param x Buffer position used to bound the diagnostic excerpt
	 */
	private void handleNoSequence(final int x){
		if(currentSection>0){return;}//This section is spuriously entered for reads that are a multiple of the target read length when splitting.
		if(WARN_IF_NO_SEQUENCE){
			synchronized(getClass()){
				if(reportedHeader==header){return;}
				reportedHeader=header;
				System.err.println("Warning: A fasta header with no sequence was encountered ("+currentSection+"):\n"+header+"\n");
				if(WARN_FIRST_TIME_ONLY){WARN_IF_NO_SEQUENCE=false;}
			}
		}
		assert(!ABORT_IF_NO_SEQUENCE) : "\n<START>"+new String(buffer, 0, Tools.min(x+1, buffer.length))+"<STOP>\n";
	}
	
	/** Compacts unread bytes and refills through a header boundary or end of input.
	 * Grows the buffer when necessary and skips initial semicolon-comment lines.
	 * Read exceptions handled while open latch errorState and are treated as EOF for
	 * buffering; callers must check the reported error. Already-closed cancellation
	 * remains tolerated. A zero result does not mean that the buffer contains no bytes.
	 * @return Newly read byte count, from the recursive call when skipping a full comment block
	 */
	private final int fillBuffer(){
//		assert(open);
		if(!open){return 0;}
		if(verbose){System.err.println("fillBuffer() : bstart="+bstart+", bstop="+bstop);}
		if(bstart<bstop){ //Shift end bytes to beginning
			if(bstart>0){
				int extra=bstop-bstart;
				for(int i=0; i<extra; i++, bstart++){
					buffer[i]=buffer[bstart];
				}
				bstop=extra;
			}
		}else{
			bstop=0;
		}
		if(verbose){System.err.println("After shift : bstart="+bstart+", bstop="+bstop);}

		
		bstart=0;
		
		int len=bstop;
		int r=-1;
		int sum=0;
		boolean seenNewline=false;
		while(len==bstop){//hit end of input without encountering a caret
			if(bstop==buffer.length){
				buffer=KillSwitch.copyOf(buffer, buffer.length*2L);
				if(verbose){System.err.println("Resized to "+buffer.length);}
			}
			if(verbose){System.err.println("A: bstop="+bstop+", len="+len);}
			try{
				r=-1;
				r=ins.read(buffer, bstop, buffer.length-bstop);
			}catch(final Exception e){
				//Use close()'s monitor to distinguish active-input errors from intentional early closure.
				//Errors handled after close remain tolerated; successful reads acquire no extra lock.
				synchronized(this){
					if(open){
						if(!errorState){System.err.println("Error reading FASTA input "+name+": "+e);}
						errorState=true;
					}
				}
			}
			if(verbose){System.err.println("B: r="+r);}
			if(r>0){
				sum+=r;
				bstop=bstop+r;
				if(bstop>0 && len==0){len=1;}//Probably to skip the first >
				
				while(len<bstop && (buffer[len]!=carrot || !seenNewline)){//I need to see a caret after newline
					seenNewline|=(buffer[len]=='\n');
					len++;
				}
				if(verbose){System.err.println("C: bstop="+bstop+", len="+len+", seenNewline="+seenNewline/*+", seenCarrot="+seenCarrot*/);}
			}else{
				len=bstop;
				if(verbose){System.err.println("D: bstop="+bstop+", len="+len);}
				break;
			}
		}
		if(verbose){System.err.println("E: bstop="+bstop+", len="+len+", r="+r);}
		
		//Skip ';'-delimited comments
		if(header==null && bstop>bstart && buffer[bstart]==';'){
			if(sum==0){return sum;}
			int lastsemi=bstart;
			assert(nextReadID==0);
			assert(bstart==0);
			while(bstop>bstart && buffer[bstart]==';'){
				while(bstop>bstart && (buffer[bstart]>slashr || buffer[bstart]<=STX || buffer[bstart]==tab)){bstart++;}
				while(bstop>bstart && buffer[bstart]<=slashr){bstart++;}
			}
			if(bstart>=bstop){ //Overflowed buffer with comments; recur
				bstart=lastsemi;
				return fillBuffer();
			}
		}
		
		assert(r==-1 || buffer[len]=='>');
		if(verbose){System.err.println("After filling: bstart="+bstart+", bstop="+bstop+", len="+len+", r="+r+", sum="+sum);}
		return sum;
	}
	
	/** Opens the configured name and resets byte-buffer bounds.
	 * Sets the open flag before delegating to ReadWrite; an opening exception propagates.
	 * Does not reset parser counters, warning history or error state.
	 * @return Underlying input stream
	 * @throws RuntimeException If already marked open, or opening fails in ReadWrite
	 */
	private final InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened fasta file "+name);
		}
		open=true;
		ins=ReadWrite.getInputStream(name, true, allowSubprocess);
		bstart=0;
		bstop=0;
		return ins;
	}
	
	/** Returns the unsynchronized open flag; not a concurrency or input-health check. */
	public boolean isOpen(){return open;}
	
	/** Validates the current global length settings, not an existing reader's captured values.
	 * Construction invokes this through an assertion after opening the input. Direct calls
	 * enforce the checks regardless of assertion mode; configure globals before use.
	 * @return true for valid settings
	 * @throws RuntimeException If the minimum or enabled split-length settings are invalid
	 */
	public static final boolean settingsOK(){
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1){
			throw new RuntimeException("Minimum FASTA read length is too long: "+MIN_READ_LEN);
		}
		if(MIN_READ_LEN<1){
			throw new RuntimeException("Minimum FASTA read length is too short: "+MIN_READ_LEN);
		}
		if(SPLIT_READS){
			if(TARGET_READ_LEN<1){
				throw new RuntimeException("Target FASTA read length is too short: "+TARGET_READ_LEN);
			}
			if(MIN_READ_LEN>TARGET_READ_LEN){
				throw new RuntimeException("Minimum FASTA read length is longer than maximum read length: "+MIN_READ_LEN+">"+TARGET_READ_LEN);
			}
		}
		if(MIN_READ_LEN>=Integer.MAX_VALUE-1 || MIN_READ_LEN<1){return false;}
		if(SPLIT_READS && (TARGET_READ_LEN<1 || MIN_READ_LEN>TARGET_READ_LEN)){return false;}
		return true;
	}
	
	/** Returns the captured input name. */
	@Override
	public String fname(){return name;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input name captured from the descriptor; also used by restart and diagnostics. */
	public final String name;
	
	/** Prefetched batch; detached on nextList so returned lists are not reused. */
	private ArrayList<Read> currentList=null;
	/** Current record header, shared by its split sections. */
	private String header=null;
	
	/** Last header reference warned about; retained across restart. */
	private String reportedHeader=null;

	/** Open state; close and the exception handler coordinate writes/checks on this monitor. */
	private boolean open=false;
	/** Reusable byte buffer, doubled when an incomplete record fills it. */
	private byte[] buffer=new byte[16384];
	/** Half-open range of unconsumed bytes in buffer. */
	private int bstart=0, bstop=0;
	/** Underlying input owned by this reader; cleared by close. Do not replace during use. */
	public InputStream ins;
	
	/** Number of list entries transferred, counting a linked pair as one. */
	private long consumed=0;
	/** Numeric ID for the next entry; both mates receive the same value. */
	private long nextReadID=0;
	/** Legacy single-entry cursor; batch access requires zero and resets it to zero. */
	private int nextReadIndex=0;
	/** One-based emitted/attempted section number within the current header, initially zero. */
	private int currentSection=0;

	/** Captured permission for subprocess-assisted input. */
	public final boolean allowSubprocess;
	/** Captured pairing mode for adjacent generated reads. */
	public final boolean interleaved;
	/** Captured amino-acid mode. */
	public final boolean amino;
	/** Read flags derived from amino mode. */
	public final int flag;
	/** Maximum list entries per batch, captured at construction. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Soft base limit checked before each entry; both mates contribute. */
	private final long MAX_DATA;
	/** Captured maximum section length and minimum accepted section length. */
	private final int maxLen, minLen;
	
	/*--------------------------------------------------------------*/
	/*----------------       Shared Settings        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Enables diagnostic output; configure before concurrent use. */
	public static boolean verbose=false;
	/** Header markers and legacy control-byte delimiters, including retained unused aliases. */
	private static final byte slashr='\r', slashn='\n', carrot='>', space=' ', tab='\t', SOH=0x1, STX=0x2;
	
	/** Whether newly constructed readers capture TARGET_READ_LEN as their section limit. */
	public static boolean SPLIT_READS=false;
	/** Global target section length; a nonpositive value is invalid when splitting is checked. */
	public static int TARGET_READ_LEN=500;
	/** Global minimum accepted section length, captured at construction. */
	public static int MIN_READ_LEN=1;
	/** Live option to allocate qualities filled with Shared.FAKE_QUAL for generated reads. */
	public static boolean FAKE_QUALITY=false;
	/** Live option to append a part suffix even to the first section. */
	public static boolean FORCE_SECTION_NAME=false;
	/** Global missing-sequence warning enablement; handleNoSequence may turn this off. */
	public static boolean WARN_IF_NO_SEQUENCE=true;
	/** Disable WARN_IF_NO_SEQUENCE globally after the first emitted warning when true. */
	public static boolean WARN_FIRST_TIME_ONLY=true;
	/** Enable handleNoSequence's guarded abort assertion; no abort effect under -da. */
	public static boolean ABORT_IF_NO_SEQUENCE=false;
	
}
