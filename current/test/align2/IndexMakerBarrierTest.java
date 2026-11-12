package test.align2;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.Map;

import align2.Block;
import align2.IndexMaker4;
import align2.IndexMaker5;
import dna.ChromosomeArray;
import dna.Data;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
/**
 * Exercises finite IndexMaker count/allocation barrier controls.
 * Adapted from Collei's 2026-09-28 script 2892 v3 IndexMakerBarrierReviewCheck fixture.
 * @author Brian Bushnell, Collei, Jean
 */
public final class IndexMakerBarrierTest{

	/** Accepts explicit case=..., maker=4|5 for builder cases, and left/right for compare. */
	public static void main(final String[] args) throws Exception{
		final Watchdog watchdog=Watchdog.launch();
		final PreParser pp=new PreParser(args, IndexMakerBarrierTest.class, false);
		try{
			String caseName=null, makerName=null, left=null, right=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected key=value argument: "+arg);}
				if(pair[0].equalsIgnoreCase("case")){caseName=pair[1];}
				else if(pair[0].equalsIgnoreCase("maker")){makerName=pair[1];}
				else if(pair[0].equalsIgnoreCase("left")){left=pair[1];}
				else if(pair[0].equalsIgnoreCase("right")){right=pair[1];}
				else if(Parser.parseStatic(arg, pair[0].toLowerCase(), pair[1])){
					// Accept normal javasetup flags such as simd=f before fixture keys.
				}else{throw new IllegalArgumentException("Unexpected argument: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("IndexMaker barrier checks require -ea.");}
			final Case kase=Case.forName(caseName);
			if(kase==Case.COMPARE){runCompare(left, right);}
			else{kase.run(Target.forName(makerName));}
			watchdog.disarm();
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private IndexMakerBarrierTest(){}

	/** Runs a real builder over a tiny A/C chromosome and checks exact starts and sites. */
	private static void runValidTiny(final Target target) throws Exception{
		final boolean old=target.allowPolymers();
		try{
			target.setAllowPolymers(true);
			setReference(false, false);
			final Block[] blocks=target.makeIndex(-1, 1, 1, 2, 0, 1000, 0, ~0,
					Integer.MAX_VALUE, 31, false, true, null);
			if(blocks[1]==null){throw new AssertionError(target.name+" did not publish a block");}
			checkArray("starts", EXPECTED_STARTS, blocks[1].starts);
			checkArray("sites", EXPECTED_SITES, blocks[1].sites);
			System.out.println("valid\t"+target.name+"\tstarts_len="+blocks[1].starts.length+
					"\tsites_len="+blocks[1].sites.length+"\tpassed=1");
		}finally{target.setAllowPolymers(old);}
	}

	/** Runs a public count failure and verifies original failure propagation plus slot release. */
	private static void runPublicFailure(final Target target) throws Exception{
		final boolean old=target.allowPolymers();
		try{
			target.setAllowPolymers(true);
			setReference(false, false);
			boolean caught=false;
			try{
				target.makeIndex(-1, 1, 1, 2, 0, 10, 0, ~0,
						Integer.MAX_VALUE, 31, false, true, null);
			}catch(final RuntimeException expected){
				caught=expected.getMessage()!=null && expected.getMessage().startsWith("Chrom 1:");
			}
			final int active=target.activeBlocks();
			if(!caught || active!=0){
				throw new AssertionError(target.name+" count failure caught="+caught+
						", ACTIVE_BLOCKS="+active);
			}
			System.out.println("public_count_failure\t"+target.name+
					"\tpropagated\tactive_slots_0");
			runValidTiny(target);
		}finally{target.setAllowPolymers(old);}
	}

	/** Whole-fixture timeout that dumps all threads and halts hung public builders. */
	private static final class Watchdog extends Thread{
		static Watchdog launch(){
			final Watchdog watchdog=new Watchdog();
			watchdog.setDaemon(true);
			watchdog.start();
			return watchdog;
		}

		private Watchdog(){super("IndexMakerBarrierTest watchdog");}

		@Override
		public void run(){
			try{Thread.sleep(10000);}
			catch(final InterruptedException expected){return;}
			if(disarmed){return;}
			System.err.println("IndexMakerBarrierTest exceeded 10s; dumping threads.");
			for(final Map.Entry<Thread, StackTraceElement[]> entry :
					Thread.getAllStackTraces().entrySet()){
				System.err.println(entry.getKey());
				for(final StackTraceElement ste : entry.getValue()){
					System.err.println("\tat "+ste);
				}
			}
			Runtime.getRuntime().halt(70);
		}

		void disarm(){
			disarmed=true;
			interrupt();
		}

		private volatile boolean disarmed=false;
	}

	/** Injects a finite worker fault and ensures the sibling barrier waiters drain. */
	private static void runInjectedFailure(final Target target, final boolean allocation)
			throws Exception{
		final boolean old=target.allowPolymers();
		final Thread[] threads=new Thread[4];
		try{
			target.setAllowPolymers(true);
			setReference(!allocation, allocation);
			final Class<?> makerClass=Class.forName("align2."+target.className+"$BlockMaker");
			final Constructor<?> makerConstructor=makerClass.getDeclaredConstructor(int.class,
					int.class, int.class, int.class, int.class, int.class, int.class,
					int.class, int.class, boolean.class, boolean.class, Block[].class);
			makerConstructor.setAccessible(true);
			final Object maker=makerConstructor.newInstance(1, 1, 2, 0, 1000, 0, ~0,
					Integer.MAX_VALUE, 31, false, true, new Block[2]);
			final Class<?> workerClass=Class.forName("align2."+target.className+
					"$BlockMaker$CountThread");
			final Constructor<?> workerConstructor=workerClass.getDeclaredConstructor(
					makerClass, int.class, int[].class, int[].class, Block[].class);
			workerConstructor.setAccessible(true);
			final int[] sizes=new int[1], intercom=new int[4];
			final Block[] holder=new Block[1];
			for(int i=0; i<threads.length; i++){
				threads[i]=(Thread)workerConstructor.newInstance(maker, i, sizes,
						intercom, holder);
				threads[i].start();
			}
			joinUntil(threads, System.currentTimeMillis()+2000);
			final int alive=countAlive(threads);
			cleanup(threads);
			final Field failure=makerClass.getDeclaredField("failure");
			failure.setAccessible(true);
			final Throwable problem=(Throwable)failure.get(maker);
			final Class<?> expected=allocation ? AssertionError.class :
					ArrayIndexOutOfBoundsException.class;
			if(alive!=0 || !expected.isInstance(problem)){
				throw new AssertionError(target.name+" "+(allocation ? "allocation" : "count")+
						" expected alive=0 and "+expected.getSimpleName()+
						" but observed alive="+alive+", failure="+problem);
			}
			System.out.println((allocation ? "allocation" : "count")+'\t'+target.name+
					"\tall_workers_terminated\t"+problem.getClass().getSimpleName());
		}finally{
			cleanup(threads);
			target.setAllowPolymers(old);
		}
	}

	/** Compares two real cached Block files through Block.read, not container hashes. */
	private static void runCompare(final String left, final String right){
		if(left==null || right==null){
			throw new IllegalArgumentException("case=compare requires left=... and right=...");
		}
		final Block a=Block.read(left), b=Block.read(right);
		final boolean starts=Arrays.equals(a.starts, b.starts);
		final boolean sites=Arrays.equals(a.sites, b.sites);
		if(!starts || !sites){
			throw new AssertionError("Block mismatch: startsEqual="+starts+
					", sitesEqual="+sites+", leftStarts="+a.starts.length+
					", rightStarts="+b.starts.length+", leftSites="+a.sites.length+
					", rightSites="+b.sites.length);
		}
		System.out.println("compare\tstarts_len="+a.starts.length+
				"\tsites_len="+a.sites.length+
				"\tstarts_equal="+starts+"\tsites_equal="+sites);
	}

	/** Installs the tiny private reference with the padding expected by Data.getChromosome. */
	private static void setReference(final boolean countFailure,
			final boolean allocationFailure){
		final ChromosomeArray ref=new ChromosomeArray(1, Shared.PLUS, 0, 99);
		for(int i=0; i<ref.array.length; i++){
			ref.array[i]=(byte)(i%2==0 ? 'A' : 'C');
		}
		if(countFailure){
			Arrays.fill(ref.array, (byte)'A');
			ref.array[ref.maxIndex-1]='C';
		}else if(allocationFailure){
			Arrays.fill(ref.array, (byte)'N');
		}
		ref.array[0]='N';
		Data.chromosomePlusMatrix=new ChromosomeArray[] {null, ref};
		Data.chromLengths=new int[] {0, ref.array.length};
	}

	/** Requires one exact array. */
	private static void checkArray(final String name, final int[] expected,
			final int[] observed){
		if(!Arrays.equals(expected, observed)){
			throw new AssertionError(name+" expected "+Arrays.toString(expected)+
					" but observed "+Arrays.toString(observed));
		}
	}

	/** Joins all workers until a shared deadline. */
	private static void joinUntil(final Thread[] threads, final long deadline)
			throws InterruptedException{
		for(final Thread thread : threads){
			if(thread!=null){thread.join(Math.max(1, deadline-System.currentTimeMillis()));}
		}
	}

	/** Interrupts and joins any leftover workers so a failed case cannot orphan a JVM. */
	private static void cleanup(final Thread[] threads) throws InterruptedException{
		for(final Thread thread : threads){
			if(thread!=null && thread.isAlive()){thread.interrupt();}
		}
		joinUntil(threads, System.currentTimeMillis()+2000);
		final int alive=countAlive(threads);
		if(alive!=0){throw new AssertionError("Could not stop "+alive+" worker(s)");}
	}

	/** Counts live workers after a bounded join. */
	private static int countAlive(final Thread[] threads){
		int alive=0;
		for(final Thread thread : threads){
			if(thread!=null && thread.isAlive()){alive++;}
		}
		return alive;
	}

	/** One explicit fixture case. */
	private enum Case{
		COUNT("count"),
		ALLOCATION("allocation"),
		VALID("valid"),
		PUBLIC_FAILURE("publicfailure"),
		COMPARE("compare");

		Case(final String name_){name=name_;}

		static Case forName(final String name){
			for(final Case c : values()){
				if(c.name.equalsIgnoreCase(name)){return c;}
			}
			throw new IllegalArgumentException("Unknown case: "+name);
		}

		void run(final Target target) throws Exception{
			if(this==COUNT){runInjectedFailure(target, false);}
			else if(this==ALLOCATION){runInjectedFailure(target, true);}
			else if(this==VALID){runValidTiny(target);}
			else if(this==PUBLIC_FAILURE){runPublicFailure(target);}
			else{throw new IllegalStateException("Unhandled builder case: "+this);}
		}

		final String name;
	}

	/** One legacy index maker implementation. */
	private enum Target{
		FOUR("4", "IndexMaker4"),
		FIVE("5", "IndexMaker5");

		Target(final String name_, final String className_){
			name=name_;
			className=className_;
		}

		static Target forName(final String name){
			for(final Target t : values()){
				if(t.name.equalsIgnoreCase(name)){return t;}
			}
			throw new IllegalArgumentException("Unknown maker: "+name);
		}

		Block[] makeIndex(final int genome, final int minChrom, final int maxChrom,
				final int k, final int chromBits, final int maxAllowedChromIndex,
				final int chromMaskLow, final int chromMaskHigh, final int siteMask,
				final int shiftLength, final boolean write, final boolean diskInvalid,
				final Block[] index){
			if(this==FOUR){
				return IndexMaker4.makeIndex(genome, minChrom, maxChrom, k,
						chromBits, maxAllowedChromIndex, chromMaskLow, chromMaskHigh,
						siteMask, shiftLength, write, diskInvalid, index);
			}
			return IndexMaker5.makeIndex(genome, minChrom, maxChrom, k,
					chromBits, maxAllowedChromIndex, chromMaskLow, chromMaskHigh,
					siteMask, shiftLength, write, diskInvalid, index);
		}

		boolean allowPolymers(){
			return this==FOUR ? IndexMaker4.ALLOW_POLYMERS : IndexMaker5.ALLOW_POLYMERS;
		}

		void setAllowPolymers(final boolean value){
			if(this==FOUR){IndexMaker4.ALLOW_POLYMERS=value;}
			else{IndexMaker5.ALLOW_POLYMERS=value;}
		}

		int activeBlocks() throws ReflectiveOperationException{
			final Field active=Class.forName("align2."+className).getDeclaredField("ACTIVE_BLOCKS");
			active.setAccessible(true);
			return active.getInt(null);
		}

		final String name;
		final String className;
	}

	/** Exact starts from Collei's 2026-09-28 v3 valid tiny IndexMaker control. */
	private static final int[] EXPECTED_STARTS={0, 0, 48, 48, 48, 97, 97, 97,
			97, 97, 97, 97, 97, 97, 97, 97, 97};

	/** Exact sites from Collei's 2026-09-28 v3 valid tiny IndexMaker control. */
	private static final int[] EXPECTED_SITES={2, 4, 6, 8, 10, 12, 14, 16, 18,
			20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48,
			50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78,
			80, 82, 84, 86, 88, 90, 92, 94, 96, 1, 3, 5, 7, 9, 11, 13,
			15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43,
			45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73,
			75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97};
}
