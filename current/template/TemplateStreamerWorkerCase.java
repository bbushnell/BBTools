package template;

import java.util.ArrayList;

import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.Writer;
import structures.ListNum;

/** Finite direct-worker check for the sample's null-only Streamer terminal contract.
 * This does not test concurrent scheduling or a concrete writer backend.
 * @author Shinobu
 */
public final class TemplateStreamerWorkerCase{

	/** Prints observed counters and submitted IDs for an ordinary or leading-empty input. */
	public static void main(final String[] args){
		if(args.length!=1 || (!args[0].equals("empty=t") && !args[0].equals("empty=f"))){
			throw new IllegalArgumentException("Expected empty=t or empty=f");
		}
		final boolean empty=args[0].equals("empty=t");
		final BatchSource source=new BatchSource(empty);
		final RecordingWriter writer=new RecordingWriter();
		final A_SampleStreamerMT.ProcessThread worker=new A_SampleStreamerMT.ProcessThread(source, writer, 0);
		worker.run();
		System.out.println(empty+"\t"+worker.readsProcessedT+"\t"+worker.basesProcessedT+"\t"+
			writer.ids+"\t"+source.calls+"\t"+worker.success);
	}

	/** Supplies an optional empty batch, one literal read and then null. */
	private static final class BatchSource implements Streamer{
		BatchSource(final boolean empty_){empty=empty_;}
		@Override
		public ListNum<Read> nextList(){
			final int index=calls++;
			if(empty && index==0){return new ListNum<Read>(new ArrayList<Read>(), 0);}
			if(index==(empty ? 1 : 0)){
				final ArrayList<Read> reads=new ArrayList<Read>();
				reads.add(new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "literal", 0));
				return new ListNum<Read>(reads, index);
			}
			return null;
		}
		@Override public String fname(){return "fixture";}
		@Override public void start(){}
		@Override public void close(){}
		@Override public boolean paired(){return false;}
		@Override public int pairnum(){return 0;}
		@Override public long readsProcessed(){return 1;}
		@Override public long basesProcessed(){return 4;}
		@Override public void setSampleRate(final float rate, final long seed){throw new UnsupportedOperationException();}
		@Override public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException();}
		@Override public boolean hasMore(){return calls<(empty ? 2 : 1);}
		@Override public boolean errorState(){return false;}
		final boolean empty;
		int calls;
	}

	/** Records submitted batch IDs without formatting or retaining mutable reads. */
	private static final class RecordingWriter implements Writer{
		@Override public void start(){}
		@Override public long readsWritten(){return 0;}
		@Override public long basesWritten(){return 0;}
		@Override public void add(final ArrayList<Read> reads, final long id){ids.add(id);}
		@Override public void addReads(final ListNum<Read> reads){add(reads.list, reads.id);}
		@Override public void addLines(final ListNum<SamLine> lines){throw new UnsupportedOperationException();}
		@Override public void poison(){}
		@Override public boolean waitForFinish(){return false;}
		@Override public boolean poisonAndWait(){return false;}
		@Override public void finishError(){throw new UnsupportedOperationException();}
		@Override public String fname(){return "fixture";}
		@Override public boolean errorState(){return false;}
		@Override public boolean finishedSuccessfully(){return true;}
		final ArrayList<Long> ids=new ArrayList<Long>();
	}
}
