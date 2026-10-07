package stream;

import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import structures.ListNum;

/**
 * Adapts queued Read jobs to a SAM/BAM Writer selected by WriterFactory.
 * The delegate is started during construction; the outer writer thread consumes
 * jobs when separately started. Local batch IDs follow dequeue order.
 * The normal output selector uses this adapter when ReadWrite.USE_READ_STREAM_SAM_WRITER
 * is enabled; the parent also uses that flag to skip legacy output initialization.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 2025
 */
public class ReadStreamSamWriter extends ReadStreamWriter{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Initializes the outer job queue, prepares delegate header input and starts the delegate.
	 * Does not start this outer Thread or change the SAM-writer selection flag.
	 * Enable ReadWrite.USE_READ_STREAM_SAM_WRITER for the parent's delegating route.
	 * Shared-header mode takes precedence over explicit text for the delegate.
	 * Otherwise, text is split on newline and each segment uses the default byte encoding;
	 * trailing empty segments follow String.split semantics. Null header input leaves
	 * header selection to the concrete backend.
	 * @param ff Nonnull SAM/BAM output descriptor passed to parent and factory
	 * @param bufferSize Positive capacity of the outer job queue
	 * @param header Optional header text, also passed unchanged to the parent constructor
	 * @param useSharedHeader Whether to request the delegate's shared-header policy
	 * @throws AssertionError If enabled format or queue-capacity assertions fail
	 */
	public ReadStreamSamWriter(FileFormat ff, int bufferSize, CharSequence header, boolean useSharedHeader){
		super(ff, null, true, bufferSize, header, true, useSharedHeader);
		assert(OUTPUT_SAM || OUTPUT_BAM) : "ReadStreamWriter requires SAM/BAM output format";
		assert(read1) : "SAM/BAM output requires read1=true (cannot write paired reads to separate files)";
		
		// Create header for Writer
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			headerLines=null; // Request the selected backend's shared-header policy
		}else if(header!=null){
			// Convert CharSequence header to ArrayList<byte[]>
			String headerStr=header.toString();
			String[] lines=headerStr.split("\n");
			headerLines=new ArrayList<byte[]>(lines.length);
			for(String line : lines){
				headerLines.add(line.getBytes());
			}
		}else{
			headerLines=null; // Leave header selection to the selected backend
		}
		
		samWriter=WriterFactory.makeWriter(ff, true, true, headerLines, useSharedHeader);
		samWriter.start();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Execution           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Processes jobs and then finalizes the delegate on normal completion.
	 * Caught Exceptions set the inherited failure flags and are wrapped in RuntimeException.
	 * This handler does not catch Error instances.
	 * @throws RuntimeException If an Exception escapes job processing or finalization
	 */
	@Override
	public void run(){
		try{
			run2();
		}catch(Exception e){
			errorState=true;
			finishedSuccessfully=false;
			System.err.println("ReadStreamWriter failed: "+e.getMessage());
			throw new RuntimeException(e);
		}
	}

	/** Processes queued jobs, then performs delegate completion if processing returns normally. */
	private void run2() throws IOException{
		processJobs();
		finishWriting();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Consumes jobs through the first poison marker, retrying interrupted queue takes.
	 * Assigns local IDs starting at zero to every non-poison job, including empty jobs;
	 * the outer job ID and close flag are not used here. Nonempty lists are wrapped
	 * without copying and submitted under the delegate's payload ownership contract.
	 * Empty jobs submit an empty SamLine list to preserve the dense ID sequence.
	 * Provisional counters count each nonnull Read and its attached mate after submission;
	 * finalization replaces these counters with the delegate's reported totals.
	 */
	private void processJobs() throws IOException{
		Job job=null;
		while(job==null){
			try{
				job=queue.take();
			}catch(InterruptedException e){
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
			}else{
				//[stream/ReadStreamSamWriter#001 FIXED 2026-08-27]
				//SamWriterST2's ordered JobQueue requires dense IDs.  Omitting an empty outer batch
				//left a permanent gap: its consumer waited for that ID while later producers filled
				//the bounded queue and deadlocked.  Preserve the ID with a shared empty payload.
				samWriter.addLines(new ListNum<SamLine>(emptyLines, listID));
			}
			
			listID++;
			
			job=null;
			while(job==null){
				try{
					job=queue.take();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Calls the delegate's normal termination-and-wait operation, then reads its counters
	 * and error flag. The operation's return value is not used; the queried error flag
	 * is combined with this adapter's existing error state.
	 * Sets finishedSuccessfully to the inverse of that aggregate error flag.
	 * @return Aggregate error state; true indicates an error
	 */
	private boolean finishWriting() throws IOException{
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

	/** Empty payload reused when forwarding empty jobs with their local batch IDs. */
	private final ArrayList<SamLine> emptyLines=new ArrayList<SamLine>(0);

	/** Factory-selected delegate, started during construction and configured for both mates. */
	private final Writer samWriter;

}
