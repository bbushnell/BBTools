package driver;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.ConcurrentLegacyReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Reports ordinary finite batch transfers through the legacy concurrent reader.
 * Source batches may cross an output-list boundary. Uses unlimited generation,
 * normal EOF and normal worker completion; this is not a shutdown or timing test.
 * Observations include exact Read identity/order and pre-sampling counters.
 * @author Shinobu
 * @date October 2, 2026
 */
public final class CheckLegacyCris{

	/** Runs one finite fixture and prints its observations after normal completion.
	 * @param args shape=split/aligned/repeat and paired=t/f; shared Parser flags accepted
	 * @throws InterruptedException If the caller is interrupted while joining the worker */
	public static void main(String[] args) throws InterruptedException{
		if(!CheckLegacyCris.class.desiredAssertionStatus()){
			throw new IllegalStateException("CheckLegacyCris requires -ea for its fixture checks");
		}
		final CheckLegacyCris check=new CheckLegacyCris(args);
		check.process();
		Shared.closeStream(check.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses the fixture choice before setting the buffer geometry used by the wrapper.
	 * @param args Standard flag=value arguments */
	private CheckLegacyCris(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		for(String arg : pp.args){
			final int equals=arg.indexOf('=');
			final String key=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase(Locale.ROOT);
			final String value=(equals<0 ? null : arg.substring(equals+1));
			if(key.equals("shape")){shape=value;}
			else if(key.equals("paired")){paired=Parse.parseBoolean(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		Parser.processQuality();
		if("split".equals(shape)){sizes=new int[]{3, 3};}
		else if("aligned".equals(shape)){sizes=new int[]{4, 2};}
		else if("repeat".equals(shape)){sizes=new int[]{3, 3, 3};}
		else{throw new IllegalArgumentException("Expected shape=split/aligned/repeat: "+shape);}
		Shared.setBufferLen(4);
		Shared.setBuffers(4);
		Shared.setBufferData(1000000);
		assert(Shared.bufferLen()==4 && Shared.numBuffers()==4)
			: "The fixture must force a three-entry source batch across a four-entry output boundary";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Drains ordinary EOF, joins the worker and then checks borrowed identities and totals.
	 * @throws InterruptedException If normal worker joining is interrupted */
	private void process() throws InterruptedException{
		final Fixture source=new Fixture(sizes, paired);
		final ConcurrentLegacyReadInputStream cris=new ConcurrentLegacyReadInputStream(source, -1);
		final Thread worker=new Thread(cris, "legacy-fixture");
		final ArrayList<Read> observed=new ArrayList<Read>();
		long nextID=0;
		boolean batchIDs=true;
		boolean sawTerminal=false;
		worker.start();
		for(ListNum<Read> list=cris.nextList(); list!=null; list=cris.nextList()){
			batchIDs&=(list.id==nextID++);
			final boolean terminal=list.isEmpty();
			observed.addAll(list.list);
			cris.returnList(list);
			if(terminal){sawTerminal=true; break;}
		}
		worker.join();
		cris.close();
		assert(sawTerminal && !worker.isAlive() && source.closed && !cris.errorState())
			: "Read/counter observations require an empty terminal, normal worker completion and source close";
		boolean exact=(observed.size()==source.roots.size() && batchIDs);
		final ByteBuilder ids=new ByteBuilder();
		for(int i=0; i<observed.size(); i++){
			final Read r=observed.get(i);
			if(i>0){ids.append(',');}
			ids.append(r.numericID);
			exact&=(i<source.roots.size() && r==source.roots.get(i));
			assert(r.length()==4 && (paired ? r.mate!=null && r.mate.mate==r && r.mate.length()==3 : r.mate==null))
				: "The wrapper must preserve the fixture's borrowed payload and reciprocal mate links";
		}
		final long expectedReads=source.roots.size()*(paired ? 2L : 1L);
		final long expectedBases=source.roots.size()*(paired ? 7L : 4L);
		exact&=(cris.readsIn()==expectedReads && cris.basesIn()==expectedBases);
		System.out.println("shape="+shape+" paired="+paired+" ids="+ids+
			" roots="+observed.size()+" reads="+cris.readsIn()+" bases="+cris.basesIn()+" exact="+exact);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Source Fixture       ----------------*/
	/*--------------------------------------------------------------*/

	/** Finite batch source retaining original Read references for an independent output oracle. */
	private static final class Fixture extends ReadInputStream{
		/** Creates distinct root reads and optional reciprocal mates in the specified batches.
		 * @param sizes Positive batch sizes no greater than the wrapper capacity
		 * @param paired Whether to attach one mate to each root */
		Fixture(int[] sizes, boolean paired){
			this.paired=paired;
			for(int size : sizes){
				assert(size>0 && size<=4) : "ReadInputStream fixture batches must fit the wrapper source-batch assertion";
				final ArrayList<Read> batch=new ArrayList<Read>(size);
				for(int i=0; i<size; i++){
					final int id=roots.size();
					final Read r=new Read(new byte[]{'A', 'C', 'G', 'T'}, null, "r"+id, id);
					if(paired){
						final Read mate=new Read(new byte[]{'T', 'G', 'C'}, null, "r"+id, id);
						r.mate=mate;
						mate.mate=r;
						mate.setPairnum(1);
					}
					batch.add(r);
					roots.add(r);
				}
				batches.add(batch);
			}
		}
		/** @return Next retained source batch, or null at ordinary EOF */
		@Override
		public ArrayList<Read> nextList(){return hasMore() ? batches.get(next++) : null;}
		/** @return Whether an unconsumed fixture batch remains */
		@Override
		public boolean hasMore(){return !closed && next<batches.size();}
		/** Resets only the fixture cursor and closed state, while retaining its Read objects. */
		@Override
		public void restart(){next=0; closed=false;}
		/** Marks the finite source closed. @return false, because the fixture has no I/O errors */
		@Override
		public boolean close(){closed=true; return false;}
		/** @return Whether this fixture attached mates */
		@Override
		public boolean paired(){return paired;}
		/** @return Fixed descriptive source name */
		@Override
		public String fname(){return "legacy-batch-fixture";}
		/** Original root references in intended output order. */
		final ArrayList<Read> roots=new ArrayList<Read>();
		/** Batches supplied by the source, never mutated by the fixture after construction. */
		final ArrayList<ArrayList<Read>> batches=new ArrayList<ArrayList<Read>>();
		/** Whether mates were constructed. */
		final boolean paired;
		/** Index of the next source batch. */
		int next=0;
		/** Whether close has been called. */
		boolean closed=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared parser diagnostic stream. */
	private final PrintStream outstream;
	/** Chosen finite batch layout. */
	private String shape="split";
	/** Whether roots carry mates. */
	private boolean paired=false;
	/** Source batch sizes selected by shape. */
	private final int[] sizes;
}
