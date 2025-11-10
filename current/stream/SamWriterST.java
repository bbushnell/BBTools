package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Single-threaded SAM writer with simple buffering.
 * Simpler alternative to SamLineWriter for cases where threading overhead isn't worth it.
 * 
 * @author Isla
 * @date November 10, 2025
 */
public class SamWriterST implements Writer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public SamWriterST(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_){
		ffout=ffout_;
		fname=ffout.name();
		header=header_;
		useSharedHeader=useSharedHeader_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);
		
		if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else {
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose) {outstream2.println("Made SamWriterST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void start(){
		writeHeader();
		started=true;
	}
	
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}

	/** Add reads for writing (will be converted to SamLines). */
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		ArrayList<SamLine> lines=toSamLines(reads.list);
		writeLines(lines);
	}

	/** Add already-formatted SamLines for writing. */
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		writeLines(lines.list);
	}

	/** Write lines immediately (single-threaded). */
	private void writeLines(ArrayList<SamLine> lines){
		if(!started){start();}
		
		ByteBuilder bb=new ByteBuilder();
		try{
			for(SamLine sl : lines){
				sl.toBytes(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=sl.length();
				
				if(bb.length()>=BUFFER_SIZE){
					outstream.write(bb.toBytes());
					bb.clear();
				}
			}
			if(bb.length()>0){
				outstream.write(bb.toBytes());
			}
		}catch(IOException e){
			errorState=true;
			throw new RuntimeException("Error writing SAM lines", e);
		}
	}

	/** Signal end of input (no-op for single-threaded). */
	public final void poison(){
		finished=true;
	}

	/** Wait for all writes to complete (no-op for single-threaded). */
	public final boolean waitForFinish(){
		return errorState();
	}

	/** Convenience method - poison and wait. */
	public final boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	@Override
	public void close(){
		if(!finished){poison();}
		ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
	}

	@Override
	public long readsWritten(){return readsWritten;}

	@Override
	public long basesWritten(){return basesWritten;}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads) {
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();

		for(final Read r1 : reads){
			Read r2=(r1==null ? null : r1.mate);

			SamLine sl1=(r1==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE 
				&& r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE 
				&& r2.samline!=null ? r2.samline : new SamLine(r2, 1)));

			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || 
				!sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}
			assert(sl1!=null) : r1;
			addSamLine(r1, sl1, samLines);
			addSamLine(r2, sl2, samLines);
		}
		return samLines;
	}

	private static void addSamLine(Read r, SamLine primary, ArrayList<SamLine> samLines) {
		if(r==null || primary==null) {return;}
		
		assert(!ReadStreamWriter.ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		samLines.add(primary);

		// Handle secondary alignments
		ArrayList<SiteScore> list=r.sites;
		if(ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
			final Read clone=r.clone();
			for(int i=1; i<list.size(); i++){
				SiteScore ss=list.get(i);
				clone.match=null;
				clone.setFromSite(ss);
				clone.setSecondary(true);
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.primary());
				assert(!ReadStreamWriter.USE_ATTACHED_SAMLINE || secondary.cigar!=null) : r;
				samLines.add(secondary);
			}
		}
	}
	
	ArrayList<byte[]> getHeader(){
		if(verbose) {outstream2.println("Fetching header: "+useSharedHeader+","+(header!=null));}
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			headerLines=SamReadInputStream.getSharedHeader(true);
		}else if(header!=null){
			headerLines=header;
		}else {
			headerLines=SamHeader.makeHeaderList(supressHeaderSequences, 
				ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(headerLines==null) {
			outstream2.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		if(verbose) {outstream2.println("Fetched header: "+(headerLines==null ? "null" : headerLines.size()));}
		return headerLines;
	}

	protected void writeHeader(){
		if(headerWritten || supressHeader){return;}
		ArrayList<byte[]> headerLines=getHeader();
		
		ByteBuilder bb=new ByteBuilder();
		try{
			for(byte[] line : headerLines) {
				bb.append(line).nl();
				if(bb.length()>=16384) {
					outstream.write(bb.toBytes());
					bb.clear();
				}
			}
			if(bb.length()>=1) {
				outstream.write(bb.toBytes());
				bb.clear();
			}
		}catch(IOException e){
			throw new RuntimeException(e);
		}
		headerWritten=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Getters and Setters      ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean errorState() {return errorState;}
	
	@Override
	public boolean finishedSuccessfully() {return !errorState && finished;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output file name. */
	final String fname;
	/** Output file format. */
	final FileFormat ffout;
	/** True if header should be pulled from shared input. */
	final boolean useSharedHeader;
	/** True if header should be skipped. */
	final boolean supressHeader;
	/** True if header sequences should be written. */
	final boolean supressHeaderSequences;
	/** True after header is written */
	boolean headerWritten=false;
	/** Header lines to write. */
	final ArrayList<byte[]> header;
	/** Output stream. */
	final OutputStream outstream;
	
	/** Total reads written. */
	private long readsWritten=0;
	/** Total bases written. */
	private long basesWritten=0;
	/** Were any errors encountered */
	private boolean errorState=false;
	/** True after start() called */
	private boolean started=false;
	/** True after poison() called */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int BUFFER_SIZE=65536;
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream2=System.err;

}