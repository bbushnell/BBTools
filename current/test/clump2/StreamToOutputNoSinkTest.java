package test.clump2;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLongArray;

import clump2.KmerComparator;
import clump2.StreamToOutput;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.Writer;
import structures.ListNum;

/** Finite no-output consumption/count/closure oracle with no background threads.
 * @author Shinobu
 */
public class StreamToOutputNoSinkTest{

	/** Runs baseline|fixed expectations for null and empty destination arrays. */
	public static void main(String[] args) throws Exception{
		assert(args.length==1) : "Expected baseline or fixed";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown expectation: "+args[0];
		System.out.println("case\treads\tbases\treturned_lists\tclosed\terror");
		for(boolean empty : new boolean[] {false, true}){
			for(boolean inputError : new boolean[] {false, true}){
				for(boolean emptyFirst : new boolean[] {false, true}){
					final FiniteInput input=new FiniteInput(inputError, emptyFirst);
					input.start();
					final StreamToOutput tool=new StreamToOutput(input,
							empty ? new Writer[0] : null,
							empty ? new AtomicLongArray(0) : null,
							new KmerComparator(5, false, false), false, false);
					final boolean error=tool.process();
					final long reads=value(tool, "readsIn"), bases=value(tool, "basesIn");
					assert(reads==(fixed ? 2 : 0) && bases==(fixed ? 16 : 0))
							: "Count-only output must consume and count both linked mates";
					assert(input.returned==(fixed ? 2 : 0) && input.closed==fixed)
							: "Count-only processing must return the two live lists and close input";
					assert(error==(fixed && inputError)) : "Input-close errors must propagate even without output";
					System.out.println((empty ? "empty" : "null")+
							(emptyFirst ? "-emptyfirst" : "-emptylast")+"-error"+inputError+
							'\t'+reads+'\t'+bases+'\t'+input.returned+'\t'+input.closed+'\t'+error);
					if(!input.closed){input.close();}
				}
			}
		}
	}

	private static long value(StreamToOutput tool, String name) throws Exception{
		final Field field=StreamToOutput.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(tool);
	}

	/** Synchronous ready-made source with one linked pair, one live empty batch and a null EOF. */
	private static class FiniteInput implements Streamer{
		FiniteInput(boolean error_, boolean emptyFirst_){
			error=error_;
			emptyFirst=emptyFirst_;
		}
		@Override
		public void start(){ready=true;}
		@Override
		public ListNum<Read> nextList(){
			assert(ready && !closed && next<3) : "Fixture source has two live lists and a null EOF";
			if(next==2){next++; return null;}
			final ArrayList<Read> list=new ArrayList<Read>();
			if(next==payloadIndex()){
				final byte[] bases="ACGTACGT".getBytes(StandardCharsets.US_ASCII);
				final Read a=new Read(bases, null, "pair/1", 0);
				final Read b=new Read(bases.clone(), null, "pair/2", 0);
				a.mate=b; b.mate=a; b.setPairnum(1);
				list.add(a);
			}
			return new ListNum<Read>(list, next++);
		}
		@Override
		public void returnList(ListNum<Read> ln){
			assert(ln.id==returned && ln.list!=null) : "Each live list must be returned once";
			returned++;
		}
		@Override
		public void returnList(long id, boolean poison){
			throw new AssertionError("Genuine null EOF should not be returned");
		}
		@Override
		public void close(){closed=true;}
		@Override
		public boolean paired(){return true;}
		@Override
		public int pairnum(){return 0;}
		@Override
		public boolean errorState(){return error;}
		@Override
		public void setSampleRate(float rate, long seed){throw new UnsupportedOperationException("No sampling");}
		@Override
		public long readsProcessed(){return next>payloadIndex() ? 2 : 0;}
		@Override
		public long basesProcessed(){return next>payloadIndex() ? 16 : 0;}
		@Override
		public String fname(){return "finite-fixture";}
		@Override
		public boolean hasMore(){return ready && !closed && next<3;}
		@Override
		public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException("No SAM fixture");}
		private int payloadIndex(){return emptyFirst ? 1 : 0;}
		final boolean error;
		final boolean emptyFirst;
		boolean ready, closed;
		int next, returned;
	}
}
