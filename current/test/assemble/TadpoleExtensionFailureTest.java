package test.assemble;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.ArrayList;

import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;

/**
 * Isolates the ordinary extension worker's exception boundary using the real private worker.
 * A synthetic reader either supplies one empty terminal batch or throws a tagged exception.
 * No count table is loaded, no real read records are processed and no output stream is constructed.
 * This does not reproduce an ordered-output hang or exercise the extendReads coordinator.
 * Baseline worker failures are observed after join; a repaired fatal worker may terminate
 * the child JVM before the observation line. Run each case in its own process.
 * @author Brian Bushnell, Shinobu
 */
public final class TadpoleExtensionFailureTest{

	/** Accepts in=existing-fixture.fa, k=31|63 and fail=t|f after standard JVM-argument filtering. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TadpoleExtensionFailureTest.class, false);
		try{
			String input=null;
			int k=0;
			Boolean fail=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected in=, k= and fail=: "+arg);}
				if(pair[0].equalsIgnoreCase("in")){input=pair[1];}
				else if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(pair[1]);}
				else if(pair[0].equalsIgnoreCase("fail")){fail=Parse.parseBoolean(pair[1]);}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			if(input==null || input.length()==0 || (k!=31 && k!=63) || fail==null){
				throw new IllegalArgumentException("Required in=fixture.fa k=31|63 fail=t|f.");
			}
			run(input, k, fail.booleanValue());
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Constructs but never runs the assembler pipeline, then observes one real worker. */
	private static void run(final String input, final int k, final boolean fail) throws Exception{
		assert(k==31 || k==63) : "The two fixture lengths select the short and long table implementations.";
		final Tadpole tad=Tadpole.makeTadpole(new String[]{"in="+input, "out=null", "k="+k,
				"mode=extend", "t=1", "prefilter=0", "prealloc=f", "minprob=0",
				"ways=1", "initialsize=128", "ecc=f", "ecco=f", "merge=f", "markerrors=f",
				"hpindel=f", "fixindels=f", "hollow=f", "hollowbranchstop=f",
				"hollowbranchrollback=f", "hollowrescue=f", "el=0", "er=1", "showstats=f"}, true);
		final Input reader=new Input(fail);
		try{
			check(tad.getClass()==(k==31 ? Tadpole1.class : Tadpole2.class), "Wrong table implementation.");
			for(String name : new String[]{"hpIndel", "localEdit", "hollowBranchStop", "hollowBranchRollback", "hollowRescue"}){
				final Field field=Tadpole.class.getDeclaredField(name);
				field.setAccessible(true);
				check(!field.getBoolean(tad), "Fixture selected an experimental worker path: "+name);
			}
			final Thread worker=worker(tad, reader);
			// One writer and a completed join make this observation visible to the main thread.
			final Throwable[] uncaught=new Throwable[1];
			worker.setUncaughtExceptionHandler((thread, failure) -> {
				uncaught[0]=failure;
				failure.printStackTrace(System.err);
			});
			System.err.println("EXTENSION_FIXTURE_START k="+k+" fail="+fail+
					" implementation="+tad.getClass().getSimpleName()+" assertions="+
					TadpoleExtensionFailureTest.class.desiredAssertionStatus());
			worker.start();
			join(worker);
			check(reader.started() && reader.nextCalls==1, "Worker did not reach exactly one synthetic input request.");
			if(fail){
				check(uncaught[0]!=null && uncaught[0].getClass()==IllegalStateException.class &&
						FAILURE.equals(uncaught[0].getMessage()), "Unexpected worker failure or missing injected exception.");
				check(reader.returns==0, "A throwing nextList cannot supply a batch to return.");
			}else{
				check(uncaught[0]==null && !tad.errorState && reader.returns==1,
						"EOF control failed or did not return its terminal batch.");
			}
			System.out.println("EXTENSION_WORKER_OBSERVATION k="+k+" fail="+fail+
					" next="+reader.nextCalls+" returned="+reader.returns+" uncaught="+(uncaught[0]!=null)+
					" error_state="+tad.errorState);
		}finally{
			reader.close();
			tad.clearData();
		}
	}

	/** Uses reflection so this fixture does not widen the production worker's access. */
	private static Thread worker(final Tadpole tad, final Input reader) throws Exception{
		assert(tad!=null && reader!=null) : "The real worker needs its assembler and one synthetic reader.";
		final Class<?> type=Class.forName("assemble.Tadpole$ExtendThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(Tadpole.class, int.class,
				ConcurrentReadInputStream[].class, ConcurrentReadOutputStream[].class, ConcurrentReadOutputStream[].class);
		constructor.setAccessible(true);
		return (Thread)constructor.newInstance(tad, 0, new ConcurrentReadInputStream[]{reader}, null, null);
	}

	/** Waits for termination before cleanup; an interrupt still fails the fixture after the worker finishes. */
	private static void join(final Thread worker){
		assert(worker!=null) : "The fixture must wait for its constructed worker before releasing assembler state.";
		boolean interrupted=false;
		while(worker.isAlive()){
			try{worker.join();}
			catch(InterruptedException e){interrupted=true;}
		}
		if(interrupted){
			Thread.currentThread().interrupt();
			throw new IllegalStateException("Extension fixture was interrupted.");
		}
	}

	/** Always-on fixture checks, independent of JVM assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Synthetic read source; inherited start launches only the no-op producer below. */
	private static final class Input extends ConcurrentReadInputStream{
		/** Chooses either the tagged exception or an empty terminal batch. */
		Input(final boolean fail_){super("synthetic-extension-input"); fail=fail_;}
		@Override
		public ListNum<Read> nextList(){
			nextCalls++;
			check(nextCalls==1, "Synthetic input was requested twice.");
			if(fail){throw new IllegalStateException(FAILURE);}
			return new ListNum<Read>(new ArrayList<Read>(), 0);
		}
		@Override
		public void returnList(final long id, final boolean poison){
			check(!fail && id==0 && poison && returns==0, "Unexpected terminal return.");
			returns++;
		}
		/** Data are supplied synchronously by nextList, so the inherited producer thread has no work. */
		@Override
		public void run(){}
		/** No input resources or queue to stop. */
		@Override
		public void shutdown(){}
		/** This one-pass fixture cannot be restarted. */
		@Override
		public void restart(){throw new UnsupportedOperationException("One-pass extension fixture.");}
		/** No resources are opened by this reader. */
		@Override
		public void close(){}
		@Override
		public boolean paired(){return false;}
		@Override
		public Object[] producers(){return new Object[0];}
		@Override
		public boolean errorState(){return false;}
		@Override
		public void setSampleRate(final float rate, final long seed){
			check(rate==1, "Synthetic input does not support sampling.");
		}
		@Override
		public long basesIn(){return 0;}
		@Override
		public long readsIn(){return 0;}
		@Override
		public boolean verbose(){return false;}
		/** Exception injection policy; counters are inspected only after the worker joins. */
		final boolean fail;
		int nextCalls, returns;
	}

	/** Distinguishes the deliberate input failure from setup, reflection or unrelated crashes. */
	private static final String FAILURE="TEST_ORDINARY_EXTENSION_INPUT_FAILURE";
}
