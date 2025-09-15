package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;

public class ReadStreamByteWriter extends ReadStreamWriter {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public ReadStreamByteWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, qfname_, read1_, bufferSize, header, false, buffered, useSharedHeader);
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
//			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
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
	
	private void writeHeader() throws IOException {
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

	private void processJobs(final ByteBuilder bb, final ByteBuilder bbq) throws IOException{
		
		Job job=null;
		while(job==null){
			try {
				job=queue.take();
//				job.list=queue.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		while(job!=null && !job.poison){

			final OutputStream os=job.outstream;
			
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
				assert(job.outstream!=null && job.outstream!=myOutstream);
				ReadWrite.finishWriting(null, job.outstream, fname, allowSubprocess); //TODO:  This should be job.fname
			}
			
			job=null;
			while(job==null){
				try {
					job=queue.take();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * @throws IOException
	 * 
	 */
	private void finishWriting(final ByteBuilder bb, final ByteBuilder bbq) throws IOException {
		if(myOutstream!=null){
			if(bb.length>0){
				myOutstream.write(bb.array, 0, bb.length);
				bb.setLength(0);
			}
			ReadWrite.finishWriting(null, myOutstream, fname, allowSubprocess);
		}
		if(myQOutstream!=null){
			if(bbq.length>0){
				myQOutstream.write(bbq.array, 0, bbq.length);
				bbq.setLength(0);
			}
			ReadWrite.finishWriting(null, myQOutstream, qfname, allowSubprocess);
		}
		finishedSuccessfully=true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeQuality(final Job job, final ByteBuilder bbq) throws IOException{
		bbq.setLength(0);
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
						if(r2.bases!=null){toQualityB(r2.quality, r2.length(), FASTA_WRAP,  bbq);}
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
					if(r2.bases!=null){toQualityB(r2.quality, r2.length(), FASTA_WRAP,  bbq);}
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeBread(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeAttachment(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeHeader(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeFasta(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeOneline(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeFastq(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeFastr(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @param os
	 * @throws IOException
	 */
	private void writeSites(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param job
	 * @param bb
	 * @throws IOException
	 */
	private void writeSam(Job job, ByteBuilder bb, OutputStream os) throws IOException {
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
	 * @param r Read to print
	 * @param primary Primary alignment of this read
	 * @param bb A buffer
	 */
	private void writeSam(Read r, SamLine primary, ByteBuilder bb) {
		if(r==null || primary==null) {return;}

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
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.primary());


				assert(!ASSERT_CIGAR || secondary.cigar!=null) : r;

				secondary.toBytes(bb).append('\n');
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final boolean buffered=true;
	private static final boolean verbose=false;
	
}
