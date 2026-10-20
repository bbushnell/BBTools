package test.clump;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import clump.KmerComparator;
import clump.StreamToOutput;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
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
				final FiniteInput input=new FiniteInput(inputError);
				input.start();
				final StreamToOutput tool=new StreamToOutput(input,
						empty ? new ConcurrentReadOutputStream[0] : null,
						new KmerComparator(5, false, false), false, false);
				final boolean error=tool.process();
				final long reads=value(tool, "readsIn"), bases=value(tool, "basesIn");
				assert(reads==(fixed ? 2 : 0) && bases==(fixed ? 16 : 0))
						: "Count-only output must consume and count both linked mates";
				assert(input.returned==(fixed ? 2 : 0) && input.closed==fixed)
						: "Count-only processing must return the payload/terminal lists and close input";
				assert(error==(fixed && inputError)) : "Input-close errors must propagate even without output";
				System.out.println((empty ? "empty" : "null")+"-error"+inputError+'\t'+reads+'\t'+bases+
						'\t'+input.returned+'\t'+input.closed+'\t'+error);
				if(!input.closed){input.close();}
			}
		}
	}

	private static long value(StreamToOutput tool, String name) throws Exception{
		final Field field=StreamToOutput.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(tool);
	}

	/** Synchronous ready-made source with one linked pair and one terminal list. */
	private static class FiniteInput extends ConcurrentReadInputStream{
		FiniteInput(boolean error_){super("finite-fixture"); error=error_;}
		@Override
		public void start(){ready=true;}
		@Override
		public ListNum<Read> nextList(){
			assert(ready && !closed && next<2) : "Fixture source has exactly one payload and one terminal";
			final ArrayList<Read> list=new ArrayList<Read>();
			if(next==0){
				final byte[] bases="ACGTACGT".getBytes(StandardCharsets.US_ASCII);
				final Read a=new Read(bases, null, "pair/1", 0);
				final Read b=new Read(bases.clone(), null, "pair/2", 0);
				a.mate=b; b.mate=a; b.setPairnum(1);
				list.add(a);
			}
			return new ListNum<Read>(list, next++);
		}
		@Override
		public void returnList(long id, boolean poison){
			assert(id==returned && poison==(id==1)) : "Each list must be returned once with its terminal flag";
			returned++;
		}
		@Override
		public void close(){closed=true;}
		@Override
		public void run(){throw new UnsupportedOperationException("No background fixture worker");}
		@Override
		public void shutdown(){closed=true;}
		@Override
		public void restart(){throw new UnsupportedOperationException("Single-use fixture");}
		@Override
		public boolean paired(){return true;}
		@Override
		public Object[] producers(){return new Object[0];}
		@Override
		public boolean errorState(){return error;}
		@Override
		public void setSampleRate(float rate, long seed){throw new UnsupportedOperationException("No sampling");}
		@Override
		public long basesIn(){return next>0 ? 16 : 0;}
		@Override
		public long readsIn(){return next>0 ? 2 : 0;}
		@Override
		public boolean verbose(){return false;}
		final boolean error;
		boolean ready, closed;
		int next, returned;
	}
}
