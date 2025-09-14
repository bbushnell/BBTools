package stream;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Set;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * This implementation allows only a single, synchronous open stream.
 * 
 * @author Brian Bushnell
 * @date May 1, 2019
 *
 */
public class MultiCros2 extends BufferedMultiCros {

	/** For testing */
	public static void main(String[] args){
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){
			names.add(args[i]);
		}
		MultiCros2 mcros=new MultiCros2(pattern, null, false, false, false, false, FileFormat.FASTQ, false);
		
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				mcros.add(r1, r1.barcode(true));
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		mcros.dumpAll();
		mcros.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** @See Details in superclass constructor */
	public MultiCros2(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, boolean threaded_){
		super(pattern1_, pattern2_, overwrite_, append_, allowSubprocess_, useSharedHeader_, defaultFormat_, threaded_, 1);
		bufferMap=new LinkedHashMap<String, Buffer>();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public boolean finishedSuccessfully(){
		return !errorState;
	}
	
	@Override
	public void add(Read r, String name){
		assert(!threaded);
		Buffer b=bufferMap.get(name);
		if(b==null){
			b=new Buffer(name);
			bufferMap.put(name, b);
			//Note: I could adjust bytesPerBuffer threshold here in response to the number of buffers.
		}
		b.add(r);
//		System.err.println("Added "+name);
	}
	
	@Override
	public long dumpResidual(ConcurrentReadOutputStream rosu){
		long dumped=0;
		
		//For each Buffer, check if it contains residual reads
		//If so, dump it into the stream
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer b=e.getValue();
			assert((b.readsIn<minReadsToDump) == (b.list!=null && !b.list.isEmpty()));
			if(b.readsIn>0 && b.readsIn<minReadsToDump){
				assert(b.list!=null && !b.list.isEmpty());
				residualReads+=b.readsIn;
				residualBases+=b.basesIn;
				if(rosu!=null){rosu.add(b.list, 0);}
			}
			b.list=null;
		}
		return dumped;
	}
	
	@Override
	public ByteBuilder report(){
		ByteBuilder bb=new ByteBuilder(1024);
		if(minReadsToDump>0){
			bb.append("Residual").tab().append(residualReads).tab().append(residualBases).nl();
		}
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			Buffer b=e.getValue();
			if(b.readsIn>=minReadsToDump){
				bb.append(b.name).tab().append(b.readsIn).tab().append(b.basesIn).nl();
			}
		}
		return bb;
	}
	
	@Override
	public Set<String> getKeys(){return bufferMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	long closeInner() {
		long x=dumpAll();
		return x;
	}
	
	@Override
	long dumpAll(){
		long dumped=0;
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump();
		}
		return dumped;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** 
	 * A Buffer holds reads destined for to a specific file.
	 * When sufficient reads are present, it opens a stream and writes them,
	 * then closes the stream.
	 */
	private class Buffer {
		
		Buffer(String name_){
			name=name_;
			String s1=pattern1.replaceFirst("%", name);
			String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);
			
			//These are created with overwrite=false append=true because 
			//the files will be appended to if the stream gets prematurely retired.
			//Therefore, files must be explicitly deleted first.
			//Alternative would be to create a new FileFormat each time.
			ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, false, true, false);
			ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, false, true, false);
			
			list=new ArrayList<Read>(readsPerBuffer);
			if(verbose){System.err.println("Made buffer for "+name);}
		}
		
		/** 
		 * Add a read to this buffer, and update all the tracking variables.
		 * This may trigger a dump.
		 */
		void add(Read r){
			list.add(r);
			long size=r.countPairBytes();
			int count=r.pairCount();
			currentBytes+=size;
			bytesInFlight+=size;
			basesIn+=size;
			readsInFlight+=count;
			readsIn+=count;
			
			if(list.size()>=readsPerBuffer || currentBytes>=bytesPerBuffer){
				if(verbose){
					System.err.println("list.size="+list.size()+"/"+readsPerBuffer+
							", bytes="+currentBytes+"/"+bytesPerBuffer+", bytesInFlight="+bytesInFlight+"/"+memLimitUpper);
				}
				dump();
			}
			
			//Too much buffered data; dump everything.
			if(bytesInFlight>=memLimitUpper){
				long dumped=dumpAll();
				if(dumped<1 && Shared.EA()){
					KillSwitch.kill("\nThis program ran out of memory."
							+ "\nTry increasing the -Xmx flag or get rid of the minreads flag,"
							+ "\nor disable assertions to skip this message and try anyway.");
				}
			}
		}
		
		/** Create a stream, dump buffered reads, and close the stream */
		long dump(){
//			System.err.println("Dumping "+name);
			if(list.isEmpty() || readsIn<minReadsToDump){return 0;}
			if(numDumps==0 && overwrite){
				delete(ff1);
				delete(ff2);
			}
//			assert(false) : counter+", "+overwrite;
			
			final long size0=list.size();
			
			ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, rswBuffers, null, useSharedHeader && numDumps==0);
			ros.start();
			ros.add(list, 0);
			errorState=ReadWrite.closeStream(ros) | errorState;
//			System.err.println("Closed stream "+name);
			
			bytesInFlight-=currentBytes;
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			numDumps++;
			list=new ArrayList<Read>(readsPerBuffer);
			return size0;
		}
		
		private void delete(FileFormat ff){
			if(ff==null){return;}
			File f=new File(ff.name());
			if(f.exists()){f.delete();}
		}
		
		/** Current list of buffered reads */
		private ArrayList<Read> list;
		
		/** Stream name, which is the variable part of the file pattern */
		private final String name;
		/** Output file 1 */
		private final FileFormat ff1;
		/** Output file 2 */
		private final FileFormat ff2;
		/** Number of reads entering the buffer */
		private long readsIn=0;
		/** Number of bases entering the buffer */
		private long basesIn=0;
		/** Number of reads written to disk */
		@SuppressWarnings("unused")
		private long readsWritten=0;//This does not count read2!
		/** Number of bytes currently in this buffer (estimated) */
		private long currentBytes=0;
		
		/** Number of dumps executed */
		private long numDumps=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Map of names to buffers */
	public final LinkedHashMap<String, Buffer> bufferMap;

}
