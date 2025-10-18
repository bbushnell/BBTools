package stream;

import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.bam.BamWriter;

/**
 * Writes BAM files from Read objects using BamWriter.
 * Parallel to ReadStreamByteWriter but dedicated to BAM format.
 *
 * @author Brian Bushnell
 * @date October 2025
 */
public class ReadStreamBamWriter extends ReadStreamWriter {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public ReadStreamBamWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, qfname_, read1_, bufferSize, header, false, true, useSharedHeader);
		assert(OUTPUT_BAM) : "ReadStreamBamWriter requires BAM output format";
		assert(read1) : "BAM output requires read1=true (cannot write paired reads to separate files)";

		try {
			bamWriter=new BamWriter(fname);
			headerWritten=false;
			this.useSharedHeader=useSharedHeader;
		} catch(IOException e){
			throw new RuntimeException("Error creating BamWriter for "+fname, e);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Execution           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void run() {
		try {
			run2();
		} catch (IOException e) {
			finishedSuccessfully=false;
			throw new RuntimeException(e);
		}
	}

	private void run2() throws IOException{
		writeHeader();
		processJobs();
		finishWriting();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	private void writeHeader() throws IOException {
		if(!headerWritten){
			ArrayList<byte[]> headerLines;
			if(useSharedHeader){
				headerLines=SamReadInputStream.getSharedHeader(true);
				if(headerLines==null){
					System.err.println("Warning: Header was null, creating empty header");
					headerLines=new ArrayList<byte[]>();
				}
			}else{
				// Generate header from Data.scaffoldNames
				headerLines=new ArrayList<byte[]>();
				// Add @HD line
				headerLines.add("@HD\tVN:1.6\tSO:unsorted".getBytes());
				// Could add @SQ lines here if needed
			}
			bamWriter.writeHeader(headerLines);
			headerWritten=true;
		}
	}

	private void processJobs() throws IOException{
		Job job=null;
		while(job==null){
			try {
				job=queue.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		while(job!=null && !job.poison){
			if(!job.isEmpty()){
				writeBam(job);
			}

			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}

	private void finishWriting() throws IOException {
		if(bamWriter!=null){
			bamWriter.close();
		}
		finishedSuccessfully=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	private void writeBam(Job job) throws IOException {
		assert(read1);
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();

		for(final Read r1 : job.list){
			Read r2=(r1==null ? null : r1.mate);

			SamLine sl1=(r1==null ? null : (USE_ATTACHED_SAMLINE && r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (USE_ATTACHED_SAMLINE && r2.samline!=null ? r2.samline : new SamLine(r2, 1)));
			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}

			addSamLine(r1, sl1, samLines);
			addSamLine(r2, sl2, samLines);
		}

		if(samLines.size()>0){
			bamWriter.writeLines(samLines);
		}
	}

	private void addSamLine(Read r, SamLine primary, ArrayList<SamLine> samLines) {
		if(r==null || primary==null) {return;}

		assert(!ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		samLines.add(primary);

		readsWritten++;
		basesWritten+=r.length();

		// Handle secondary alignments
		ArrayList<SiteScore> list=r.sites;
		if(OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
			final Read clone=r.clone();
			for(int i=1; i<list.size(); i++){
				SiteScore ss=list.get(i);
				clone.match=null;
				clone.setFromSite(ss);
				clone.setSecondary(true);
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.primary());
				assert(!ASSERT_CIGAR || secondary.cigar!=null) : r;
				samLines.add(secondary);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Instance Fields       ----------------*/
	/*--------------------------------------------------------------*/

	private final BamWriter bamWriter;
	private boolean headerWritten;
	private final boolean useSharedHeader;

}
