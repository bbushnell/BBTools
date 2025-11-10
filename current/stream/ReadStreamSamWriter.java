package stream;

import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import structures.ListNum;

/**
 * Writes SAM/BAM files from Read objects using Writer.
 * Wraps the new multithreaded Writer/BamLineWriter architecture
 * to fit into the ReadStreamWriter interface.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 2025
 */
public class ReadStreamSamWriter extends ReadStreamWriter {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public ReadStreamSamWriter(FileFormat ff, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, null, true, bufferSize, header, true, useSharedHeader);
		assert(OUTPUT_SAM || OUTPUT_BAM) : "ReadStreamWriter requires SAM/BAM output format";
		assert(read1) : "SAM/BAM output requires read1=true (cannot write paired reads to separate files)";
		
		// Create header for Writer
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			headerLines=null; // Writer will pull from shared header
		}else if(header!=null){
			// Convert CharSequence header to ArrayList<byte[]>
			String headerStr=header.toString();
			String[] lines=headerStr.split("\n");
			headerLines=new ArrayList<byte[]>(lines.length);
			for(String line : lines){
				headerLines.add(line.getBytes());
			}
		}else{
			headerLines=null; // Writer will generate from Data.scaffoldNames
		}
		
		samWriter=WriterFactory.makeWriter(ff, true, true, headerLines, useSharedHeader);
		samWriter.start();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Execution           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void run() {
		try {
			run2();
		} catch (Exception e) {
			errorState=true;
			finishedSuccessfully=false;
			System.err.println("ReadStreamWriter failed: "+e.getMessage());
			throw new RuntimeException(e);
		}
	}

	private void run2() throws IOException{
		processJobs();
		finishWriting();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private void processJobs() throws IOException{
		Job job=null;
		while(job==null){
			try {
				job=queue.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		long listID=0;
		while(job!=null && !job.poison){
			if(!job.isEmpty()){
				// Convert Job to ListNum<Read>
				ListNum<Read> ln=new ListNum<Read>(job.list, listID);
				samWriter.addReads(ln);
				
				// Update statistics
				for(Read r : job.list){
					if(r!=null){
						readsWritten++;
						basesWritten+=r.length();
						if(r.mate!=null){
							readsWritten++;
							basesWritten+=r.mate.length();
						}
					}
				}
			}
			
			listID++;
			
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

	private boolean finishWriting() throws IOException {
		samWriter.poisonAndWait();
		
		// Accumulate statistics from Writer
		readsWritten=samWriter.readsWritten();
		basesWritten=samWriter.basesWritten();
		errorState|=samWriter.errorState();
		
		finishedSuccessfully=!errorState;
		return errorState;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Instance Fields       ----------------*/
	/*--------------------------------------------------------------*/

	private final Writer samWriter;

}