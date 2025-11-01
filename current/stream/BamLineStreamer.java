package stream;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import stream.bam.BamReader;
import stream.bam.BamToSamConverter;
import stream.bam.BgzfInputStream;
import stream.bam.BgzfInputStreamMT;
import stream.bam.BgzfSettings;
import structures.ListNum;

/**
 * Loads BAM files rapidly with multiple threads.
 * Thread 0 reads BAM binary format and converts to SAM text.
 * Worker threads convert SAM text byte[] to SamLine objects.
 *
 * @author Chloe
 * @contributor Isla
 * @date October 18, 2025
 */
public class BamLineStreamer extends SamStreamer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public BamLineStreamer(String fname_, int threads_, boolean saveHeader_, 
		boolean ordered_, long maxReads_, boolean makeReads_){
		this(FileFormat.testInput(fname_, FileFormat.BAM, null, true, false), threads_, 
			saveHeader_, ordered_, maxReads_, makeReads_);
	}

	/** Constructor. */
	public BamLineStreamer(FileFormat ffin_, int threads_, boolean saveHeader_, 
		boolean ordered_, long maxReads_, boolean makeReads_){
		super(ffin_, threads_, saveHeader_, ordered_, maxReads_, makeReads_);
		final int queueSize=2*threads+3;
		outq=new JobQueue<ListNum<SamLine>>(queueSize, ordered, true, 0);
		if(verbose) {System.err.println("Made BamLineStreamer-"+threads);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean hasMore() {return outq.hasMore();}
	
	@Override
	public boolean paired(){return false;}

	@Override
	public int pairnum(){return 0;}
	
	@Override
	public long readsProcessed() {return readsProcessed;}
	
	@Override
	public long basesProcessed() {return basesProcessed;}

	@Override
	public ListNum<SamLine> nextLines(){
		ListNum<SamLine> list=outq.take();
		if(list==null || list.last()) {
			assert(list==null || list.isEmpty());
			assert(!outq.hasMore());
			return null;
		}
		if(verbose && list!=null){outstream.println("Got list size "+list.size());}
		return list;
	}

	@Override
	public ListNum<Read> nextReads(){
		assert(makeReads);
		ListNum<SamLine> lines=nextLines();
		if(lines==null){return null;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		if(!lines.isEmpty()) {
			for(SamLine line : lines){
				assert(line.obj!=null);
				reads.add((Read)line.obj);
			}
		}
		ListNum<Read> ln=new ListNum<Read>(reads, lines.id);
		return ln;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawn process threads */
	@Override
	void spawnThreads(){
		//Determine how many threads may be used
		final int threads=this.threads+1;

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, alpt));
		}
		if(verbose){outstream.println("Spawned threads.");}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Processing thread - reads BAM records and converts to SamLine lists. */
	private class ProcessThread extends Thread {

		/** Constructor */
		ProcessThread(final int tid_, ArrayList<ProcessThread> alpt_){
			tid=tid_;
			alpt=(tid==0 ? alpt_ : null);
		}

		/** Called by start() */
		@Override
		public void run(){
			//Process the reads
			if(tid==0){
				processBamBytes();
			}else{
				makeReads();
			}

			//Indicate successful exit status
			success=true;
		}

		void processBamBytes(){
			long listNumber=0;
			try{
				FileInputStream fis=new FileInputStream(fname);
				java.io.InputStream bgzf;
				if(BgzfSettings.USE_MULTITHREADED_BGZF){
					int threads=Math.max(1, BgzfSettings.READ_THREADS);
					bgzf=new BgzfInputStreamMT(fis, threads);
					if(verbose) {System.err.println("Made bismt-"+threads);}
				}else{
					bgzf=new BgzfInputStream(fis);
				}
				BamReader reader=new BamReader(bgzf);

				//Read BAM magic
				byte[] magic=reader.readBytes(4);
				if(!Arrays.equals(magic, new byte[]{'B', 'A', 'M', 1})){
					throw new RuntimeException("Not a BAM file: "+fname);
				}

				//Read header text
				long l_text=reader.readUint32();
				byte[] text=reader.readBytes((int)l_text);

				//Parse header if requested
				if(saveHeader && header!=null){
					synchronized(header){
						//Split by newline and add to header
						int start=0;
						for(int i=0; i<text.length; i++){
							if(text[i]=='\n'){
								if(i>start){
									byte[] line=Arrays.copyOfRange(text, start, i);
									header.add(line);
								}
								start=i+1;
							}
						}
						//Add last line if not ending with newline
						if(start<text.length){
							byte[] line=Arrays.copyOfRange(text, start, text.length);
							header.add(line);
						}
						SamReadInputStream.setSharedHeader(header); //Set shared header
					}
				}

				if(verbose){outstream.println("Thread "+tid+" reading sequence lines.");}

				//Read reference sequence dictionary
				int n_ref=reader.readInt32();
				String[] refNames=new String[n_ref];
				for(int i=0; i<n_ref; i++){
					long l_name=reader.readUint32();
					refNames[i]=reader.readString((int)l_name-1); //Exclude NUL
					reader.readUint8(); //Skip NUL terminator
					long l_ref=reader.readUint32(); //Reference length (unused here)
				}

				if(verbose){outstream.println("Thread "+tid+" making converter.");}
				synchronized(BamLineStreamer.this){
					sharedConverter=new BamToSamConverter(refNames);
					BamLineStreamer.this.notifyAll();
				}
				if(verbose){outstream.println("Thread "+tid+" made converter.");}

				//Read alignment records
				ArrayList<byte[]> list=new ArrayList<byte[]>(LIST_SIZE);
				try{
					for(long reads=0; reads<maxReads; reads++){
						long block_size=reader.readUint32();
						byte[] bamRecord=reader.readBytes((int)block_size);
						list.add(bamRecord);

						if(list.size()>=LIST_SIZE){
							putBytes(new ListNum<byte[]>(list, listNumber));
							listNumber++;
							list=new ArrayList<byte[]>(LIST_SIZE);
						}
					}
				}catch(EOFException e){
					//Normal end of file
				}

				if(list.size()>0){
					putBytes(new ListNum<byte[]>(list, listNumber));
					listNumber++;
				}
				putBytes(new ListNum<byte[]>(null, listNumber, true, false));//Poison

				bgzf.close();
				fis.close();

			}catch(IOException e){
				throw new RuntimeException("Error reading BAM file: "+fname, e);
			}
			if(verbose){outstream.println("Thread "+tid+" finished.");}

			success=true;

			//Wait for completion of all threads
			boolean allSuccess=true;
			for(ProcessThread pt : alpt){

				//Wait until this thread has terminated
				if(pt!=this){
					if(verbose){outstream.println("Waiting for thread "+pt.tid);}
					while(pt.getState()!=Thread.State.TERMINATED){
						try{
							pt.join();
						}catch(InterruptedException e){
							e.printStackTrace();
						}
					}

					//Accumulate per-thread statistics
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					allSuccess&=pt.success;
				}
			}

			putReads(new ListNum<SamLine>(null, listNumber, false, true));
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning reads.");}

			//Track whether any threads failed
			if(!allSuccess){errorState=true;}
			if(verbose || verbose2){outstream.println("tid "+tid+" finished!");}
		}

		void putReads(ListNum<SamLine> list){
			if(verbose){outstream.println("tid "+tid+" putting rlist size "+list.size());}
			outq.add(list);
			if(verbose){outstream.println("tid "+tid+" done putting rlist");}
		}

		/** Iterate through the reads */
		void makeReads(){
			if(verbose){System.err.println("Thread "+tid+" waiting on converter.");}
			synchronized(BamLineStreamer.this){
				while(sharedConverter==null){
					try{
						BamLineStreamer.this.wait(100);
					}catch(InterruptedException e){
						e.printStackTrace();
					}
				}
				converter=sharedConverter;
			}

			if(verbose){outstream.println("tid "+tid+" started makeReads.");}
			ListNum<byte[]> list=takeBytes();
			while(list!=null && !list.poison()){
				ListNum<SamLine> reads=new ListNum<SamLine>(
					new ArrayList<SamLine>(list.size()), list.id);
				long readID=list.id*200;//TODO: Should be part of the listNum
				for(byte[] bamRecord : list){

//					final SamLine sl=new SamLine(converter.convertAlignment(bamRecord));//Obsolete - reparse
					final SamLine sl=converter.toSamLine(bamRecord);
					assert(sl!=null);
					if(sl!=null){
						if(makeReads){
							Read r=sl.toRead(FASTQ.PARSE_CUSTOM);
							sl.obj=r;
							r.samline=sl;
							r.numericID=readID++;
							if(!r.validated()){r.validate(true);}
						}
						reads.add(sl);

						readsProcessedT++;
						basesProcessedT+=(sl.seq==null ? 0 : sl.length());
					}
				}
				if(reads.size()>0){putReads(reads);}
				list=takeBytes();
			}
			if(verbose || verbose2){outstream.println("tid "+tid+" done making reads.");}

			putBytes(list);
			if(verbose || verbose2){outstream.println("tid "+tid+" done poisoning bytes.");}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** True only if this thread has completed successfully */
		boolean success=false;
		/** Thread ID */
		final int tid;

		ArrayList<ProcessThread> alpt;
		BamToSamConverter converter;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private volatile BamToSamConverter sharedConverter;
	final JobQueue<ListNum<SamLine>> outq;

}