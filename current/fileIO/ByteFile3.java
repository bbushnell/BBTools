package fileIO;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import shared.Parse;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.HasID;
import stream.OrderedQueueSystem;
import structures.IntList;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Multithreaded byte file reader with parallel line extraction.
 * 
 * Uses producer-consumer pattern:
 * - Producer thread reads file and finds newlines
 * - Worker threads extract individual lines from buffer chunks
 * - Output maintains original line order
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 2, 2025
 */
public final class ByteFile3 extends ByteFile {


	public static void main(String[] args){
		String fname=null;
		long first=0, last=100;
		boolean speedtest=false;
		int threads=DEFAULT_THREADS;
		int bufferSize=bufferlen;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("input")){
				fname=b;
			}else if(a.equals("speedtest")){
				speedtest=Parse.parseBoolean(b);
				if(speedtest){
					first=0;
					last=Long.MAX_VALUE;
				}
			}else if(a.equals("first") || a.equals("start")){
				first=Long.parseLong(b);
			}else if(a.equals("last") || a.equals("stop")){
				last=Long.parseLong(b);
			}else if(a.equals("lines")){
				first=0;
				last=Long.parseLong(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("buffersize") || a.equals("buffer")){
				bufferSize=(int)Parse.parseKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//Do nothing
			}else if(arg.indexOf('=')<0 && fname==null){
				// First non-flag argument is the filename
				fname=arg;
			}else{
				System.err.println("Unknown parameter "+arg);
				assert(false) : "Unknown parameter "+arg;
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		if(fname==null){fname="stdin";}
		
		ByteFile3 bf=new ByteFile3(fname, true, threads);
		ByteFile3.bufferlen=bufferSize;
		bf.start();
		speedtest(bf, first, last, !speedtest);
		
		bf.close();
		bf.reset();
		bf.close();
	}

	private static void speedtest(ByteFile3 bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length+1;
				System.out.println(new String(s));
			}

			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
//			for(long i=first; i<last; i++){
//				byte[] s=bf.nextLine();
//				if(s==null){break;}
//				lines++;
//				bytes+=s.length+1;
//			}
			for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()){
				for(byte[] line : ln.list) {
					lines++;
					bytes+=line.length+1;
				}
			}
		}
		t.stop();

		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
			System.err.println("Bytes: "+bytes);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public ByteFile3(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}

	public ByteFile3(String fname, boolean allowSubprocess_, int threads_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false), threads_);
	}

	public ByteFile3(FileFormat ff){
		this(ff, DEFAULT_THREADS);
	}
	
	public ByteFile3(FileFormat ff, int threads_){
		super(ff);
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());

		// Create OQS with prototypes
		BufferJob inputPrototype=new BufferJob(null, null, 0, BufferJob.PROTO);
		ListNum<byte[]> outputPrototype=new ListNum<byte[]>(null, 0, ListNum.PROTO);
		oqs=new OrderedQueueSystem<BufferJob, ListNum<byte[]>>(
			threads, true, inputPrototype, outputPrototype);

		if(verbose){System.err.println("Made ByteFile3 with "+threads+" threads");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public void start(){
		if(open){return;}
		open();
		spawnThreads();
	}

	@Override
	public final void reset(){
		close();
		open();
		superReset();
	}

	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open);}
		if(!open){return errorState;}
		open=false;

		if(is!=null){
			errorState|=ReadWrite.finishReading(is, name(), allowSubprocess());
			is=null;
		}

		lineNum=-1;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name());}
		return errorState;
	}

	@Override
	public final byte[] nextLine(){
		// Get next line from current list, or fetch new list
		if(currentList==null || listPos>=currentList.size()){
			currentList=oqs.getOutput();
			listPos=0;

			if(currentList==null || currentList.last()){
				if(currentList!=null && currentList.last()){
					oqs.setFinished();
				}
				return null;
			}
		}

		lineNum++;
		return currentList.get(listPos++);
	}
	
	@Override
	public final ListNum<byte[]> nextList(){
		ListNum<byte[]> ln=oqs.getOutput();
		if(ln==null || ln.last()){
			oqs.setFinished();
			return null;
		}

		lineNum++;
		return ln;
	}

	@Override
	public void pushBack(byte[] line){
		throw new UnsupportedOperationException("pushBack not supported in ByteFile3");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private void spawnThreads(){
		final int totalThreads=threads+1; // Workers plus producer

		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){
			alpt.add(new ProcessThread(i));
		}

		for(ProcessThread pt : alpt){
			pt.start();
		}
	}

	private synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened ByteFile3 "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		return is;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean isOpen(){return open;}

	@Override
	public final InputStream is(){return is;}

	@Override
	public final long lineNum(){return lineNum;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Processing thread - produces BufferJobs or processes them into lines */
	private class ProcessThread extends Thread {

		ProcessThread(final int tid_){
			tid=tid_;
			setName("ByteFile3-"+(tid==0 ? "Producer" : "Worker-"+tid));
		}

		@Override
		public void run(){
			if(tid==0){
				produceBufferJobs();
			}else{
				processBufferJobs();
			}
			success=true;
		}

		/** Producer thread - reads file and creates BufferJobs */
		void produceBufferJobs(){
			long jobId=0;
			IntList positions=new IntList();
			byte[] buffer=new byte[bufferlen];
			int bstop=0;

			try{
				while(true){
					// Fill buffer
					int r=is.read(buffer, bstop, buffer.length-bstop);
					if(r>0){
						bstop+=r;
					}

					// Find newlines
					positions.clear();
					Vector.findSymbols(buffer, 0, bstop, slashn, positions);

					// If no newlines and buffer is full, expand
					if(positions.size()==0 && bstop==buffer.length){
						buffer=Arrays.copyOf(buffer, buffer.length*2);
						continue;
					}

					// If we hit EOF with no newlines, package what we have
					if(r<0){
						if(bstop>0){
							byte[] copy=Arrays.copyOf(buffer, bstop);
							BufferJob job=new BufferJob(copy, positions.copy(), jobId++);
							oqs.addInput(job);
						}
						break;
					}

					// Package job with data up to last newline
					if(positions.size()>0){
						int lastNL=positions.get(positions.size()-1);
						byte[] copy=Arrays.copyOf(buffer, lastNL+1);
						BufferJob job=new BufferJob(copy, positions.copy(), jobId++);
						oqs.addInput(job);

						// Shift remaining bytes to start
						int remaining=bstop-lastNL-1;
						System.arraycopy(buffer, lastNL+1, buffer, 0, remaining);
						bstop=remaining;
					}
				}
			}catch(IOException e){
				e.printStackTrace();
			}

			// Signal completion
			oqs.poison();

			// Wait for workers
			ThreadWaiter.waitForThreadsToFinish(alpt);
		}

		/** Worker thread - extracts lines from BufferJobs */
		void processBufferJobs(){
			BufferJob job=oqs.getInput();
			while(job!=null && !job.poison()){
				ListNum<byte[]> output=new ListNum<byte[]>(
					new ArrayList<byte[]>(job.positions.size()+1), job.id());

				int start=0;
				for(int i=0; i<job.positions.size(); i++){
					int nlpos=job.positions.get(i);
					int limit=(nlpos>0 && job.buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;

					if(start==limit){
						output.add(blankLine);
					}else{
						byte[] line=Arrays.copyOfRange(job.buffer, start, limit);
						output.add(line);
					}
					start=nlpos+1;
				}

				// Handle any remaining bytes after last newline
				if(start<job.buffer.length){
					byte[] line=Arrays.copyOfRange(job.buffer, start, job.buffer.length);
					if(line.length>0){
						output.add(line);
					}
				}

				oqs.addOutput(output);
				job=oqs.getInput();
			}

			// Re-inject poison for other workers
			if(job!=null){oqs.addInput(job);}
		}

		boolean success=false;
		final int tid;
	}

	/** Job for passing buffer chunks between producer and workers */
	private static final class BufferJob implements HasID {

		BufferJob(byte[] buffer_, IntList positions_, long id_){
			this(buffer_, positions_, id_, NORMAL);
		}

		BufferJob(byte[] buffer_, IntList positions_, long id_, int type_){
			buffer=buffer_;
			positions=positions_;
			id=id_;
			type=type_;
		}

		@Override
		public long id(){return id;}

		@Override
		public boolean poison(){return type==POISON;}

		@Override
		public boolean last(){return type==LAST;}

		@Override
		public BufferJob makePoison(long id_){
			return new BufferJob(null, null, id_, POISON);
		}

		@Override
		public BufferJob makeLast(long id_){
			return new BufferJob(null, null, id_, LAST);
		}

		final byte[] buffer;
		final IntList positions;
		final long id;
		final int type;

		static final int PROTO=-1;
		static final int NORMAL=0;
		static final int LAST=3;
		static final int POISON=4;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final OrderedQueueSystem<BufferJob, ListNum<byte[]>> oqs;
	private final int threads;
	private ArrayList<ProcessThread> alpt;

	private boolean open=false;
	private InputStream is;
	private long lineNum=-1;

	private ListNum<byte[]> currentList=null;
	private int listPos=0;

	private static final byte[] blankLine=new byte[0];

	public static boolean verbose=false;
	public static boolean BUFFERED=false;
	public static int bufferlen=131072;
	public static int DEFAULT_THREADS=1;

	private boolean errorState=false;
}