package test.prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import prok.CallGenesFileOutput;
import stream.Read;
import structures.ByteBuilder;

/** Ordered three-file handoff, bounded payload splitting and failure propagation.
 * @author Brian, Raiden */
public final class CallGenesFileOutputTest{
	public static void main(String[] args) throws Exception{
		ordered();bufferedCompletion();failure();
		System.out.println("CALLGENES_FILE_OUTPUT_TEST_PASS files=3 empty_batches=true failure=true");
	}
	/** Later-file producers must be able to finish a multi-MiB output before the
	 * input-order publisher reaches them, while every chunk remains bounded. */
	private static void bufferedCompletion() throws Exception{
		final CallGenesFileOutput.Failure failure=new CallGenesFileOutput.Failure();
		final int capacity=CallGenesFileOutput.queueCapacity(8L<<30, 32);
		final CallGenesFileOutput output=new CallGenesFileOutput(failure, capacity);
		final byte[] payload=new byte[5<<20];
		for(int i=0;i<payload.length;i++){payload[i]=(byte)((i>>>16)+31*i);}
		final Thread producer=new Thread(){
			@Override public void run(){
				try{output.gff(new ByteBuilder().append(payload));output.finish();}
				catch(Throwable t){output.fail(t);}
			}
		};
		producer.start();producer.join(5000);
		final boolean completed=!producer.isAlive();
		if(!completed){output.fail(new IllegalStateException("Output budget blocked the fixture producer"));producer.join();}
		assert(completed) : "An8GiB heap must let a5MiB later-file output finish without its publisher; capacity="+capacity;
		final ByteBuilder drained=new ByteBuilder(payload.length);
		for(CallGenesFileOutput.Chunk c=output.take(); c!=null; c=output.take()){
			assert(c.channel==0 && c.gff.length()<=65536) : "A larger queue must not enlarge individual GFF chunks";
			drained.append(c.gff);
		}
		assert(Arrays.equals(payload, drained.toBytes())) : "Buffered completion must preserve every output byte and chunk order";
	}
	private static void ordered() throws Exception{
		final CallGenesFileOutput.Failure failure=new CallGenesFileOutput.Failure();
		final CallGenesFileOutput[] files=new CallGenesFileOutput[3];
		final Thread[] workers=new Thread[3];
		final CountDownLatch release=new CountDownLatch(1);
		for(int i=0; i<3; i++){
			final int ordinal=i;
			files[i]=new CallGenesFileOutput(failure);
			workers[i]=new Thread(){
				@Override public void run(){
					try{
						if(ordinal==0){release.await();}
						files[ordinal].gff(new ByteBuilder().append("file").append(ordinal).nl());
						files[ordinal].reads(1, null);
						final byte[] bases=new byte[150000];Arrays.fill(bases, (byte)('A'+ordinal));
						files[ordinal].gff(new ByteBuilder().append(bases));
						final ArrayList<Read> reads=new ArrayList<Read>();
						reads.add(new Read(new byte[]{(byte)('A'+ordinal)}, null, "same-id", 0));
						files[ordinal].reads(2, reads);files[ordinal].finish();
					}catch(Throwable t){files[ordinal].fail(t);}
				}
			};
			workers[i].start();
		}
		release.countDown();
		for(int i=0; i<3; i++){
			final ByteBuilder bytes=new ByteBuilder();int empty=0, sequences=0;
			for(CallGenesFileOutput.Chunk chunk=files[i].take(); chunk!=null; chunk=files[i].take()){
				if(chunk.channel==0){
					assert(chunk.gff.length()<=65536) : "Output queue must not retain an entire GFF";
					bytes.append(chunk.gff);
				}else if(chunk.channel==1){assert(chunk.reads.isEmpty()) : "Zero-call batches retain their place";empty++;}
				else{
					assert(chunk.channel==2 && chunk.reads.size()==1) : "Sequence channel/read count changed";
					assert(chunk.reads.get(0).bases[0]=='A'+i) : "Repeated names must not mix files";sequences++;
				}
			}
			workers[i].join();
			final byte[] payload=new byte[150000];Arrays.fill(payload, (byte)('A'+i));
			final ByteBuilder expected=new ByteBuilder().append("file").append(i).nl().append(payload);
			assert(Arrays.equals(bytes.toBytes(), expected.toBytes())) : "Chunk boundaries must preserve all bytes in file order";
			assert(empty==1 && sequences==1) : "Empty and nonempty sequence batches must each occur once";
		}
	}
	private static void failure(){
		final CallGenesFileOutput.Failure failure=new CallGenesFileOutput.Failure();
		final CallGenesFileOutput output=new CallGenesFileOutput(failure);
		failure.fail(new IllegalArgumentException("fixture"));
		boolean rejected=false;
		try{output.take();}catch(IllegalStateException expected){rejected=true;}
		assert(rejected) : "Publication must fail instead of waiting for a missing terminal";
		rejected=false;
		try{output.finish();}catch(IllegalStateException expected){rejected=true;}
		assert(rejected) : "A failed publication must also release blocked producers";
	}
}
