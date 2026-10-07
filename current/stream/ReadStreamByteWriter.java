package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;

/**
 * Legacy byte-buffered writer for queued read batches and optional separate qualities.
 * The parent initializes output streams; this class consumes FIFO jobs through a poison
 * marker and dispatches each nonempty batch by format. Primary read mode may interleave
 * mates; mate-only mode selects mates from the submitted first-read list.
 * SAM and site output have their own pair-selection rules. Counters are maintained while
 * rendering records, before all buffered bytes necessarily reach the output stream.
 */
public class ReadStreamByteWriter extends ReadStreamWriter{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Initializes legacy outputs and the job queue without starting this thread.
	 * Header setup follows the parent constructor's format and suppression rules.
	 * @param ff Output descriptor passed to the parent
	 * @param qfname_ Optional separate quality filename
	 * @param read1_ True to render submitted reads, false to select their mates
	 * @param bufferSize Positive queue capacity
	 * @param header Optional explicit header passed to the parent
	 * @param useSharedHeader Whether parent SAM-header setup requests shared records
	 */
	public ReadStreamByteWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, qfname_, read1_, bufferSize, header, buffered, useSharedHeader);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Execution           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs format-specific output and finalization on normal completion.
	 * Caught IOExceptions clear finishedSuccessfully and are wrapped in RuntimeException.
	 * This handler does not set errorState or catch unchecked exceptions.
	 * @throws RuntimeException If an IOException escapes writing or finalization
	 */
	@Override
	public void run(){
		try{
			run2();
		}catch(IOException e){
			finishedSuccessfully=false;
//			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	/** Emits the format preamble, consumes jobs using shared buffers, then finalizes outputs. */
	private void run2() throws IOException{
		writeHeader();

		final ByteBuilder bb=new ByteBuilder(65000);
		final ByteBuilder bbq=(myQOutstream==null ? null : new ByteBuilder(65000));

		processJobs(bb, bbq);
		finishWriting(bb, bbq);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Emits the additional FASTR or native-text/site preamble selected by this writer.
	 * SAM, FASTQ, FASTA, attachment, header-only and oneline formats emit nothing here.
	 * This preamble path does not consult the parent's header-suppression flags.
	 * @throws IOException If a preamble write fails
	 */
	private void writeHeader() throws IOException{
		if(!OUTPUT_SAM && !OUTPUT_FASTQ && !OUTPUT_FASTA && !OUTPUT_ATTACHMENT && !OUTPUT_HEADER && !OUTPUT_ONELINE){
			if(OUTPUT_FASTR){
				myOutstream.write("#FASTR".getBytes());
				if(OUTPUT_INTERLEAVED){myOutstream.write("\tINT".getBytes());}
				myOutstream.write('\n');
			}else{
				if(OUTPUT_INTERLEAVED){
					//				assert(false) : OUTPUT_SAM+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+", "+OUTPUT_INTERLEAVED+", "+SITES_ONLY;
					myOutstream.write("#INTERLEAVED\n".getBytes());
				}
				if(SITES_ONLY){
					myOutstream.write(("#"+SiteScore.header()+"\n").getBytes());
				}else if(!OUTPUT_ATTACHMENT){
					myOutstream.write(("#"+Read.header()+"\n").getBytes());
				}
			}
		}
	}

	/**
	 * Consumes FIFO jobs through the first poison marker, retrying interrupted queue takes.
	 * Writes separate qualities before primary records for nonempty jobs. The close flag
	 * flushes and finalizes only the primary output; it does not terminate this loop.
	 * Job IDs are not inspected here. Both buffers remain shared across batch boundaries.
	 * @param bb Primary record buffer
	 * @param bbq Separate-quality buffer, or null when no quality output exists
	 * @throws IOException If a write fails
	 */
	private void processJobs(final ByteBuilder bb, final ByteBuilder bbq) throws IOException{

		Job job=null;
		while(job==null){
			try{
				job=queue.take();
//				job.list=queue.take();
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		while(job!=null && !job.poison){

			final OutputStream os=myOutstream;

			if(!job.isEmpty()){
				if(myQOutstream!=null){
					writeQuality(job, bbq);
				}

				if(OUTPUT_SAM){
					writeSam(job, bb, os);
				}else if(SITES_ONLY){
					writeSites(job, bb, os);
				}else if(OUTPUT_FASTQ){
					writeFastq(job, bb, os);
				}else if(OUTPUT_FASTA){
					writeFasta(job, bb, os);
				}else if(OUTPUT_ONELINE){
					writeOneline(job, bb, os);
				}else if(OUTPUT_ATTACHMENT){
					writeAttachment(job, bb, os);
				}else if(OUTPUT_HEADER){
					writeHeader(job, bb, os);
				}else if(OUTPUT_FASTR){
					writeFastr(job, bb, os);
				}else{
					writeBread(job, bb, os);
				}
			}
			if(job.close){
				if(bb.length>0){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
				boolean b=ReadWrite.finishWriting(null, myOutstream, fname, allowSubprocess);
				errorState|=b;
			}

			job=null;
			while(job==null){
				try{
					job=queue.take();
				}catch(InterruptedException e){
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Writes remaining buffered bytes and finalizes each configured output.
	 * The primary finish result is accumulated into errorState; the quality finish result
	 * is currently discarded. Normal return sets finishedSuccessfully independently of
	 * errorState; the separate source TODOs retain these known status limitations.
	 *
	 * @param bb ByteBuilder containing remaining sequence data
	 * @param bbq ByteBuilder containing remaining quality data
	 * @throws IOException if final write operations fail
	 */
	private synchronized void finishWriting(final ByteBuilder bb, final ByteBuilder bbq) throws IOException{
		if(myOutstream!=null){
			if(bb.length>0){
				myOutstream.write(bb.array, 0, bb.length);
				bb.setLength(0);
			}
			boolean b=ReadWrite.finishWriting(null, myOutstream, fname, allowSubprocess);
			errorState|=b;
		}
		if(myQOutstream!=null){
			if(bbq.length>0){
				myQOutstream.write(bbq.array, 0, bbq.length);
				bbq.setLength(0);
			}
			//TODO: Probable bug (STR229) - the quality close error result is discarded, unlike the primary result.
			ReadWrite.finishWriting(null, myQOutstream, qfname, allowSubprocess);
		}
		//TODO: Probable bug (STR230) - completion is marked successful even when errorState is true.
		finishedSuccessfully=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Appends separate quality records with greater-than-prefixed IDs and a final newline.
	 * Uses parent quality rendering only when bases are present. Primary mode selects each
	 * nonnull submitted read and optional interleaved mate, writing the buffer after each
	 * list entry. Mate-only mode requires reciprocal nonnull mates and writes at the
	 * 32768-byte threshold. Any remaining tail survives for the next job or finalization.
	 * Quality output does not advance the primary read/base counters.
	 * @param job Nonempty job containing submitted first reads
	 * @param bbq Shared separate-quality buffer, retaining any preceding tail
	 * @throws IOException If a quality-output write fails
	 */
	private void writeQuality(final Job job, final ByteBuilder bbq) throws IOException{
		//Retain an unflushed mate-quality tail from the preceding job until the next flush.
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					{
						bbq.append('>');
						bbq.append(r.id);
						bbq.append('\n');
						if(r.bases!=null){toQualityB(r.quality, r.length(), FASTA_WRAP, bbq);}
						bbq.append('\n');
					}
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						bbq.append('>');
						bbq.append(r2.id);
						bbq.append('\n');
						if(r2.bases!=null){toQualityB(r2.quality, r2.length(), FASTA_WRAP, bbq);}
						bbq.append('\n');
					}
				}
				if(bbq.length>=32768 || true){
					myQOutstream.write(bbq.array, 0, bbq.length);
					bbq.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
					assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
					bbq.append('>');
					bbq.append(r2.id);
					bbq.append('\n');
					if(r2.bases!=null){toQualityB(r2.quality, r2.length(), FASTA_WRAP, bbq);}
					bbq.append('\n');
				}
				if(bbq.length>=32768){
					myQOutstream.write(bbq.array, 0, bbq.length);
					bbq.setLength(0);
				}
			}
		}

//		if(bbq.length>0){
//			myQOutstream.write(bbq.array, 0, bbq.length);
//			bbq.setLength(0);
//		}
	}

	/**
	 * Writes reads in BREAD format (BBTools native text format).
	 * Outputs complete read information including metadata using Read.toText().
	 * Handles interleaved mode and read1/read2 selection with 32KB buffering.
	 *
	 * @param job Job containing reads to write
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if BREAD writing fails
	 */
	private void writeBread(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					r.toText(true, bb).append('\n');
					readsWritten++;
					basesWritten+=r.length();
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						r2.toText(true, bb).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}

				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
//					assert(r2!=null && r2.mate==r1 && r2!=r1) : r1.toText(false);
					if(r2!=null){
						r2.toText(true, bb).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}else{
						//TODO os.print(".\n");
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}
	}

	/**
	 * Writes attachment objects or SAM lines associated with reads.
	 * Outputs either the read's attached object (toString()) or samline data.
	 * An attached object takes precedence over a SamLine. Selected reads increment the read
	 * counter even when neither attachment exists; this path does not count bases.
	 *
	 * @param job Job containing reads with attachments
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if attachment writing fails
	 */
	private void writeAttachment(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					if(r.obj!=null){bb.append(r.obj.toString()).nl();}
					else if(r.samline!=null){r.samline.toBytes(bb).nl();}
					readsWritten++;
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						if(r2.obj!=null){bb.append(r2.obj.toString()).nl();}
						else if(r2.samline!=null){r2.samline.toBytes(bb).nl();}
						readsWritten++;
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
					if(r2!=null){
						if(r2.obj!=null){bb.append(r2.obj.toString()).nl();}
						else if(r2.samline!=null){r2.samline.toBytes(bb).nl();}
						readsWritten++;
					}else{
//						bb.append('.').append('\n');
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}
	}

	/**
	 * Writes only read headers/identifiers without sequence data.
	 * Outputs read IDs one per line, handling interleaved mode and
	 * read1/read2 selection. Counts selected reads but does not advance the base counter.
	 *
	 * @param job Job containing reads to extract headers from
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if header writing fails
	 */
	private void writeHeader(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					bb.append(r.id).append('\n');
					readsWritten++;
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						bb.append(r2.id).append('\n');
						readsWritten++;
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
					if(r2!=null){
						bb.append(r2.id).append('\n');
						readsWritten++;
					}else{
//						bb.append('.').append('\n');
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}
	}

	/**
	 * Writes reads in FASTA format with configurable line wrapping.
	 * Uses Read.toFasta() method with FASTA_WRAP setting for proper formatting.
	 * Handles interleaved output and maintains read/base counts.
	 *
	 * @param job Job containing reads to write in FASTA format
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if FASTA writing fails
	 */
	private void writeFasta(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					r.toFasta(FASTA_WRAP, bb).append('\n');
					readsWritten++;
					basesWritten+=r.length();
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						r2.toFasta(FASTA_WRAP, bb).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
					assert(ignorePairAssertions || (r2!=null && r2.mate==r1 && r2!=r1)) : "\n"+r1.toText(false)+"\n\n"+(r2==null ? "null" : r2.toText(false)+"\n");
					if(r2!=null){
						r2.toFasta(FASTA_WRAP, bb).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}
	}

	/**
	 * Writes reads in tab-delimited one-line format (ID \t SEQUENCE).
	 * Each read becomes a single line with tab-separated identifier and bases.
	 * Compact format useful for downstream processing tools.
	 *
	 * @param job Job containing reads to write in one-line format
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if one-line writing fails
	 */
	private void writeOneline(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					bb.append(r.id).append('\t').append(r.bases).append('\n');
					readsWritten++;
					basesWritten+=r.length();
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						bb.append(r2.id).append('\t').append(r2.bases).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
					assert(ignorePairAssertions || (r2!=null && r2.mate==r1 && r2!=r1)) : "\n"+r1.toText(false)+"\n\n"+(r2==null ? "null" : r2.toText(false)+"\n");
					if(r2!=null){
						bb.append(r2.id).append('\t').append(r2.bases).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}
	}

	/**
	 * Writes reads in standard FASTQ format with quality scores.
	 * Uses Read.toFastq() method for proper 4-line FASTQ formatting.
	 * Handles interleaved output and maintains read/base statistics.
	 *
	 * @param job Job containing reads to write in FASTQ format
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if FASTQ writing fails
	 */
	private void writeFastq(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		if(read1){
			for(final Read r : job.list){
				if(r!=null){
					r.toFastq(bb).append('\n');
					readsWritten++;
					basesWritten+=r.length();
					Read r2=r.mate;
					if(OUTPUT_INTERLEAVED && r2!=null){
						r2.toFastq(bb).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}else{
			for(final Read r1 : job.list){
				if(r1!=null){
					final Read r2=r1.mate;
					assert(ignorePairAssertions || (r2!=null && r2.mate==r1 && r2!=r1)) : "\n"+r1.toText(false)+"\n\n"+(r2==null ? "null" : r2.toText(false)+"\n");
					if(r2!=null){
						r2.toFastq(bb).append('\n');
						readsWritten++;
						basesWritten+=r2.length();
					}
				}
				if(bb.length>=32768){
					os.write(bb.array, 0, bb.length);
					bb.setLength(0);
				}
			}
		}
	}

	/**
	 * Writes reads in FASTR format (BBTools fast read format).
	 * Each job appends its list size, then all selected IDs, sequences and qualities.
	 * The count is the submitted list size even when interleaved mates add records.
	 * Requires nonnull entries and, in mate-only mode, nonnull mates. Read/base counters
	 * advance during sequence rendering; the buffer is written at the threshold after
	 * the complete job has been appended.
	 *
	 * @param job Job containing reads to write in FASTR format
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if FASTR writing fails
	 */
	private void writeFastr(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		bb.append(job.list.size()).append('\n');
		if(read1){
			for(final Read r : job.list){
				bb.append(r.id).append('\n');
				Read r2=r.mate;
				if(OUTPUT_INTERLEAVED && r2!=null){
					bb.append(r2.id).append('\n');
				}
			}
			for(final Read r : job.list){
				bb.append(r.bases).append('\n');
				readsWritten++;
				basesWritten+=r.length();

				Read r2=r.mate;
				if(OUTPUT_INTERLEAVED && r2!=null){
					bb.append(r2.bases).append('\n');
					readsWritten++;
					basesWritten+=r2.length();
				}
			}
			for(final Read r : job.list){
				bb.appendQuality(r.quality).append('\n');
				Read r2=r.mate;
				if(OUTPUT_INTERLEAVED && r2!=null){
					bb.appendQuality(r2.quality).append('\n');
				}
			}
		}else{
			for(final Read r1 : job.list){
				final Read r2=r1.mate;
				bb.append(r2.id).append('\n');
			}
			for(final Read r1 : job.list){
				final Read r2=r1.mate;
				bb.append(r2.bases).append('\n');
				readsWritten++;
				basesWritten+=r2.length();
			}
			for(final Read r1 : job.list){
				final Read r2=r1.mate;
				bb.appendQuality(r2.quality).append('\n');
			}
		}

		if(bb.length>=32768){
			os.write(bb.array, 0, bb.length);
			bb.setLength(0);
		}
	}

	/**
	 * Writes alignment sites information for reads.
	 * Requires primary mode. Submitted reads are rendered only when their sites list is
	 * nonnull; every nonnull mate is rendered regardless of its sites or interleaving.
	 * Read.toSites renders a dot when there are no sites. Each rendered record advances
	 * both read and base counters.
	 *
	 * @param job Job containing reads with sites data
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream to write to
	 * @throws IOException if sites writing fails
	 */
	private void writeSites(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		assert(read1);
		for(final Read r : job.list){
			Read r2=(r==null ? null : r.mate);

			if(r!=null && r.sites!=null){
				r.toSites(bb).append('\n');

				readsWritten++;
				basesWritten+=r.length();
			}
			if(r2!=null){
				r2.toSites(bb).append('\n');

				readsWritten++;
				basesWritten+=r2.length();
			}
			if(bb.length>=32768){
				os.write(bb.array, 0, bb.length);
				bb.setLength(0);
			}
		}
	}

	/**
	 * Writes each submitted read and its mate in SAM format, requiring primary mode.
	 * Reuses attached SamLines when enabled, otherwise constructs them. Unless KEEP_NAMES
	 * is set, a differing mate QNAME is replaced by the first read's QNAME, including in
	 * a reused attached line. Optional secondary records are handled by the single-read helper.
	 *
	 * @param job Job containing reads to write in SAM format
	 * @param bb ByteBuilder for output formatting
	 * @param os OutputStream used to flush buffered SAM records
	 * @throws IOException if SAM writing fails
	 */
	private void writeSam(Job job, ByteBuilder bb, OutputStream os) throws IOException{
		assert(read1);
		for(final Read r1 : job.list){
			Read r2=(r1==null ? null : r1.mate);

			SamLine sl1=(r1==null ? null : (USE_ATTACHED_SAMLINE && r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (USE_ATTACHED_SAMLINE && r2.samline!=null ? r2.samline : new SamLine(r2, 1)));
			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}

			writeSam(r1, sl1, bb);
			writeSam(r2, sl2, bb);
			if(bb.length>=32768){
				os.write(bb.array, 0, bb.length);
				bb.setLength(0);
			}
		}
	}

	/**
	 * Writes a single read and its alignments in SAM format.
	 * Outputs the primary alignment and optionally secondary alignments.
	 * Ignores null reads or null primary lines. Counts the read and bases once for the
	 * primary record; secondary lines add no counts. With secondary output enabled, a
	 * clone is reused for sites starting at index 1 and is marked secondary each time.
	 *
	 * @param r Read to write
	 * @param primary Primary alignment SamLine for this read
	 * @param bb ByteBuilder for SAM output formatting
	 */
	private void writeSam(Read r, SamLine primary, ByteBuilder bb){
		if(r==null || primary==null){return;}

		assert(!ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		primary.toBytes(bb).append('\n');

		readsWritten++;
		basesWritten+=r.length();
		ArrayList<SiteScore> list=r.sites;
		if(OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
			final Read clone=r.clone();
			for(int i=1; i<list.size(); i++){
				SiteScore ss=list.get(i);
				clone.match=null;
				clone.setFromSite(ss);
				clone.setSecondary(true);
				//TODO: Probable bug (STR372) - when Data.scaffoldLocs is null and clone.samline is nonnull,
				//SamLine copies the attachment, ignoring the new site and secondary flag. Caller reachability is unverified.
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.nonSecondary());

				assert(!ASSERT_CIGAR || secondary.cigar!=null) : r;

				secondary.toBytes(bb).append('\n');
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Requests buffered legacy streams from the parent constructor. */
	private static final boolean buffered=true;
	/** Reserved diagnostic switch; this implementation does not currently read it. */
	private static final boolean verbose=false;

}
