package tax;

import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Builds disk-backed binary tables for the taxonomy server.
 * Reads NCBI accession-to-taxid files and GI tables, writes
 * memory-mappable binary files.
 *
 * Usage:
 *   BuildDiskTables accession=auto pattern=auto table=auto tree=auto \
 *     accout=accession_disk.bin giout=gi_disk.bin
 *
 * @author Brian Bushnell, Chloe
 */
public class BuildDiskTables {

	public static void main(String[] args){
		Timer t=new Timer();
		BuildDiskTables x=new BuildDiskTables(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public BuildDiskTables(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_UNPIGZ=true;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("pattern")){
				patternFile=b;
			}else if(a.equals("accout") || a.equals("accessionout")){
				accOut=b;
			}else if(a.equals("giout")){
				giOut=b;
			}else if(a.equals("bufferentries")){
				bufferEntries=Parse.parseIntKMG(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		if("auto".equalsIgnoreCase(patternFile)){patternFile=TaxTree.defaultPatternFile();}
	}

	void process(Timer t){
		if(taxTreeFile!=null){
			outstream.println("Loading tax tree.");
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, false, false);
		}

		if(patternFile!=null){
			outstream.println("Loading pattern table.");
			AnalyzeAccession.loadCodeMap(patternFile);
		}

		if(giTableFile!=null && giOut!=null){
			buildGiTable();
		}

		if(accessionFile!=null && accOut!=null){
			if(AnalyzeAccession.codeMap==null){
				throw new RuntimeException("Building a disk accession table requires pattern=...");
			}
			buildAccessionTable();
		}

		t.stop();
		outstream.println("Total time: "+t);
	}

	/*--------------------------------------------------------------*/

	private void buildGiTable(){
		Timer t=new Timer();
		outstream.println("\n--- Building GI disk table ---");

		GiToTaxid.initialize(giTableFile);

		try{
			DiskGiTable.build(giOut, GiToTaxid.array, GiToTaxid.maxGiLoaded);
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		t.stop();
		outstream.println("GI table built in "+t);
		GiToTaxid.unload();
	}

	/*--------------------------------------------------------------*/

	private void buildAccessionTable(){
		Timer t=new Timer();
		outstream.println("\n--- Building accession disk table ---");

		final String[] files=accessionFile.split(",");

		outstream.println("Pass 1: counting digitizable entries.");
		final CountStats stats=CountStats.sum(spawnCountThreads(files));
		outstream.println("Total lines: "+stats.total);
		outstream.println("Digitizable: "+stats.digitizable);
		outstream.println("Non-digitizable: "+stats.nonDigitizable);

		outstream.println("Pass 2: building disk table.");
		long collected=0, added=0;
		try(DiskAccessionTable.Builder builder=DiskAccessionTable.build(accOut, stats.digitizable)){
			final DiskAccessionTable.BufferSet buffers=builder.bufferSet(bufferEntries());
			final ProcessThread[] threads=spawnThreads(files, buffers);
			collected=waitForThreads(threads);
			buffers.flush();
			added=buffers.added();
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		outstream.println("Collected "+collected+" entries.");
		outstream.println("Added "+added+" entries.");
		t.stop();
		outstream.println("Accession table built in "+t);
	}

	private static boolean skipAccessionLine(byte[] line){
		return line.length<1 || Tools.startsWith(line, "accession") || Tools.startsWith(line, "#");
	}

	private CountThread[] spawnCountThreads(String[] files){
		final int threads=Tools.mid(1, files.length, Shared.threads());
		final CountThread[] array=new CountThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
		for(int i=0; i<array.length; i++){
			array[i]=new CountThread(files, next, error);
			array[i].start();
		}
		return array;
	}

	private ProcessThread[] spawnThreads(String[] files, DiskAccessionTable.BufferSet buffers){
		final int threads=Tools.mid(1, files.length, Shared.threads());
		final ProcessThread[] array=new ProcessThread[threads];
		final AtomicInteger next=new AtomicInteger(0);
		final AtomicReference<Throwable> error=new AtomicReference<Throwable>();
		for(int i=0; i<array.length; i++){
			array[i]=new ProcessThread(files, next, error, buffers);
			array[i].start();
		}
		return array;
	}

	private static long waitForThreads(ProcessThread[] threads){
		long collected=0;
		for(int i=0; i<threads.length; i++){
			while(threads[i].getState()!=Thread.State.TERMINATED){
				try{
					threads[i].join();
				}catch(InterruptedException e){
					throw new RuntimeException(e);
				}
			}
			if(threads[i].error!=null){throw new RuntimeException(threads[i].error);}
			collected+=threads[i].collected;
		}
		return collected;
	}

	private int bufferEntries(){
		if(bufferEntries>0){return bufferEntries;}
		final long memory=Runtime.getRuntime().maxMemory();
		final long bytes=memory/(32L*DiskAccessionTable.BUILD_PARTITIONS);
		final long entries=bytes/BUFFER_ENTRY_BYTES;
		return (int)Tools.mid(1024, MAX_BUFFER_ENTRIES, entries);
	}

	private static class CountThread extends Thread {

		CountThread(String[] files_, AtomicInteger next_, AtomicReference<Throwable> error_){
			files=files_;
			next=next_;
			errorRef=error_;
		}

		@Override
		public void run(){
			try{
				for(int i=next.getAndIncrement(); i<files.length && errorRef.get()==null; i=next.getAndIncrement()){
					processFile(files[i]);
				}
			}catch(Throwable e){
				error=e;
				errorRef.compareAndSet(null, e);
			}
		}

		private void processFile(String fname) throws IOException {
			ByteFile bf=ByteFile.makeByteFile(fname, true);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(skipAccessionLine(line)){continue;}
				total++;
				if(AnalyzeAccession.digitize(line)>0){digitizable++;}
				else{nonDigitizable++;}
			}
			if(bf.close()){throw new IOException("Error reading "+fname);}
		}

		private final String[] files;
		private final AtomicInteger next;
		private final AtomicReference<Throwable> errorRef;
		private long total=0;
		private long digitizable=0;
		private long nonDigitizable=0;
		private Throwable error=null;

	}

	private static class CountStats {

		private static CountStats sum(CountThread[] threads){
			CountStats stats=new CountStats();
			for(int i=0; i<threads.length; i++){
				while(threads[i].getState()!=Thread.State.TERMINATED){
					try{
						threads[i].join();
					}catch(InterruptedException e){
						throw new RuntimeException(e);
					}
				}
				if(threads[i].error!=null){throw new RuntimeException(threads[i].error);}
				stats.total+=threads[i].total;
				stats.digitizable+=threads[i].digitizable;
				stats.nonDigitizable+=threads[i].nonDigitizable;
			}
			return stats;
		}

		private long total=0;
		private long digitizable=0;
		private long nonDigitizable=0;

	}

	private static class ProcessThread extends Thread {

		ProcessThread(String[] files_, AtomicInteger next_, AtomicReference<Throwable> error_,
				DiskAccessionTable.BufferSet buffers_){
			files=files_;
			next=next_;
			errorRef=error_;
			buffers=buffers_;
		}

		@Override
		public void run(){
			try{
				for(int i=next.getAndIncrement(); i<files.length && errorRef.get()==null; i=next.getAndIncrement()){
					processFile(files[i]);
				}
			}catch(Throwable e){
				error=e;
				errorRef.compareAndSet(null, e);
			}
		}

		private void processFile(String fname) throws IOException {
			ByteFile bf=ByteFile.makeByteFile(fname, true);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(skipAccessionLine(line)){continue;}
				final long number=AnalyzeAccession.digitize(line);
				if(number>0){
					final int taxid=AccessionToTaxid.parseLineToTaxid(line, (byte)'\t');
					if(taxid>0){
						buffers.add(number, taxid);
						collected++;
					}
				}
			}
			if(bf.close()){throw new IOException("Error reading "+fname);}
		}

		private final String[] files;
		private final AtomicInteger next;
		private final AtomicReference<Throwable> errorRef;
		private final DiskAccessionTable.BufferSet buffers;
		private long collected=0;
		private Throwable error=null;

	}

	/*--------------------------------------------------------------*/

	private String accessionFile=null;
	private String giTableFile=null;
	private String taxTreeFile="auto";
	private String patternFile=null;
	private String accOut=null;
	private String giOut=null;
	private int bufferEntries=0;
	private TaxTree tree=null;
	private PrintStream outstream=System.err;

	private static final long BUFFER_ENTRY_BYTES=20;
	private static final int MAX_BUFFER_ENTRIES=1<<24;

}
