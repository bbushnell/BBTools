package test.bloom;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;

import bloom.IndexCounter;
import bloom.KCountArray;
import bloom.KCountArray2;
import bloom.KCountArray3;
import bloom.KCountArray4;
import bloom.KCountArray4MT;
import bloom.KCountArray5MT;
import bloom.KCountArray6MT;
import bloom.KCountArray7MT;
import bloom.KCountArray7MTA;
import bloom.KCountArray8MT;

/** Finite contract checks for Bloom count-array boundary cases.
 * @author Brian Bushnell, Jean */
public final class BloomArrayContractTest{

	/** Exercises dual-hash prefiltering, lock modes, saturation and long-key arrays. */
	public static void main(final String[] args) throws Exception{
		String mode="all";
		String outdir=null;
		for(final String arg : args){
			if(arg.startsWith("mode=")){mode=arg.substring(5);}
			else if(arg.startsWith("outdir=")){outdir=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown argument "+arg);}
		}
		if("precise-neighbor-sums".equals(mode)){
			testPreciseNeighborSums(outdir);
			System.out.println("BLOOM_ARRAY_CONTRACT_TEST_PASS checks="+checks);
			return;
		}
		if(!"all".equals(mode)){throw new IllegalArgumentException("Unknown mode "+mode);}
		testDualHashPrefilterThreshold();
		testLockedAndUnlockedParity();
		testSaturation();
		testPreciseNeighborSums(outdir);
		testSignedPackedWordOccupancyScans();
		testLongArrayKeys();
		testIndexCounterAssertionGuard();
		System.out.println("BLOOM_ARRAY_CONTRACT_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Verifies dual-hash prefilter gating below and at the bypass threshold. */
	private static void testDualHashPrefilterThreshold(){
		final long key=17, key2=29;
		final KCountArray7MTA pre=new KCountArray7MTA(1009, 8, 3, null, 0);
		pre.increment(key, key2, 1);

		final KCountArray7MTA gated=new KCountArray7MTA(1009, 8, 3, pre, 2);
		check(gated.incrementAndReturnUnincremented(key, key2, 1)==1,
				"Prefilter counts below the limit must return without incrementing the main array.");
		check(gated.read(key, key2)==1, "Below-limit reads must report the prefilter count.");
		gated.purgeFilter();
		check(gated.read(key, key2)==0, "Below-limit prefilter increments must leave the main array untouched.");

		pre.increment(key, key2, 1);
		final KCountArray7MTA admitted=new KCountArray7MTA(1009, 8, 3, pre, 2);
		check(admitted.incrementAndReturnUnincremented(key, key2, 1)==0,
				"At-limit prefilter counts must admit the key into the main array.");
		check(admitted.read(key, key2)==1, "Admitted dual-hash keys must be readable from the main array.");
		check(admitted.read(key, key2+1)==0, "Dual-hash prefilters must not admit primary-only collisions.");
	}

	/** Ensures explicit lock modes preserve the same sequential visible contract. */
	private static void testLockedAndUnlockedParity(){
		final boolean oldLocked=KCountArray.LOCKED_INCREMENT;
		final boolean oldSet=KCountArray.SET_LOCKED_INCREMENT;
		try{
			final int[] locked=exerciseLockMode(true);
			final int[] unlocked=exerciseLockMode(false);
			final int[] expected={0, 2, 5, 255, 0, 4, 6};
			check(Arrays.equals(locked, expected), "Locked mode changed the sequential increment contract.");
			check(Arrays.equals(unlocked, expected), "Unlocked mode changed the sequential increment contract.");
			check(Arrays.equals(locked, unlocked), "Locked and unlocked modes diverged for the same sequence.");
		}finally{
			KCountArray.LOCKED_INCREMENT=oldLocked;
			KCountArray.SET_LOCKED_INCREMENT=oldSet;
		}
	}

	/** Exercises supported saturation paths without long increment loops. */
	private static void testSaturation() throws Exception{
		final long key=17;
		final long key2=0x6789A;
		final int legacyStart=0xFFFE;

		final KCountArray2 legacy=new KCountArray2(1024, 16);
		legacy.write(key, legacyStart);
		check(legacy.increment(key, 7)==0xFFFF,
				"KCountArray2 must saturate supported 16-bit cells.");
		check(legacy.read(key)==0xFFFF, "KCountArray2 must retain the saturated value.");

		final int atomicStart=Integer.MAX_VALUE-1;
		final KCountArray7MTA atomic=new KCountArray7MTA(10007, 32, 3, null, 0);
		raiseAllLanes(atomic, key, key2, atomicStart, 3);
		check(atomic.incrementAndReturnUnincremented(key, key2, 7)==atomicStart,
				"KCountArray7MTA must report the pre-increment 32-bit value.");
		check(atomic.read(key, key2)==Integer.MAX_VALUE,
				"KCountArray7MTA must saturate 32-bit cells instead of overflowing.");
	}

	/** Checks precise smoothing when a four-neighbor 32-bit side sum exceeds int range. */
	private static void testPreciseNeighborSums(final String outdir) throws Exception{
		final StringBuilder failures=new StringBuilder();
		final long[] oldSeeds=resetPreciseNeighborSeeds();
		try(final BufferedWriter writer=preciseNeighborWriter(outdir)){
			checkPreciseNeighborCase(failures, "left-overflow", 5,
					new int[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE},
					new int[]{1, 1, 1, 1}, Integer.MAX_VALUE, 4, 5, 4, writer);
			checkPreciseNeighborCase(failures, "right-overflow", 5,
					new int[]{1, 1, 1, 1},
					new int[]{Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE},
					4, Integer.MAX_VALUE, 5, 4, writer);
			checkPreciseNeighborCase(failures, "nonoverflow", 20,
					new int[]{1, 2, 3, 4}, new int[]{1, 1, 2, 4},
					10, 8, 9, 8, writer);
		}finally{
			restorePreciseNeighborSeeds(oldSeeds);
		}
		check(failures.length()==0, "Precise-neighbor failures: "+failures);
	}

	/** Verifies every signed-int scanner counts finite packed words with bit 31 clear and set. */
	private static void testSignedPackedWordOccupancyScans() throws Exception{
		final Object[] arrays={
				new KCountArray2(1024, 16),
				new KCountArray3(1024, 16),
				new KCountArray4(1024, 16, 0, 1),
				new KCountArray4MT(1024, 16, 0, 1),
				new KCountArray5MT(1024, 16, 0, 1),
				new KCountArray6MT(1024, 16, 0, 1),
				new KCountArray7MT(1024, 16, 1),
				new KCountArray8MT(1024, 16, 1, new KCountArray7MTA(1024, 16, 1, null, 0))
		};
		final Object[] arrays32={
				new KCountArray3(1024, 32),
				new KCountArray4(1024, 32, 0, 1),
				new KCountArray4MT(1024, 32, 0, 1),
				new KCountArray5MT(1024, 32, 0, 1),
				new KCountArray6MT(1024, 32, 0, 1),
				new KCountArray7MT(1024, 32, 1),
				new KCountArray8MT(1024, 32, 1, new KCountArray7MTA(1024, 32, 1, null, 0))
		};
		String comma="";
		String failures="";

		final KCountArray2 publicLegacy=new KCountArray2(1024, 16);
		publicLegacy.write(4, 0xFFFE);
		check(publicLegacy.increment(4, 7)==0xFFFF,
				"KCountArray2 must publicly saturate the high-half cell.");
		final String publicFailure=signedScanFailure(publicLegacy, 0xFFFF);
		if(publicFailure!=null){
			failures+="KCountArray2 public high="+publicFailure;
			comma=", ";
		}

		for(final Object array : arrays){
			seedPackedWord(array, 0x0000FFFF);
			final String lowFailure=signedScanFailure(array, 0xFFFF);
			if(lowFailure!=null){
				failures+=comma+array.getClass().getSimpleName()+" low="+lowFailure;
				comma=", ";
			}

			seedPackedWord(array, 0xFFFF0000);
			final String highFailure=signedScanFailure(array, 0xFFFF);
			if(highFailure!=null){
				failures+=comma+array.getClass().getSimpleName()+" high="+highFailure;
				comma=", ";
			}
		}

		for(final Object array : arrays32){
			seedPackedWord(array, Integer.MAX_VALUE);
			final String wideFailure=signedScanFailure(array, Integer.MAX_VALUE);
			if(wideFailure!=null){
				failures+=comma+array.getClass().getSimpleName()+" 32bit="+wideFailure;
				comma=", ";
			}
		}
		check(failures.length()==0,
				"Packed occupancy reflection scan controls failed: "+failures);
	}

	/** Checks that every word in a long-key array contributes to placement. */
	private static void testLongArrayKeys(){
		final KCountArray7MTA array=new KCountArray7MTA(1000003, 8, 3, null, 0);
		final long[] a={0x12345678L, 0xCAFEBABEL, 0x55AA55AAL};
		final long[] b={0x12345678L, 0xCAFEBABEL, 0x55AA55ABL};

		check(array.incrementAndReturnUnincremented(a, 2)==0, "New long-array keys must start at zero.");
		check(array.read(a)==2, "The inserted long-array key must be readable.");
		check(array.read(b)==0, "A later long-array word must affect the count-array placement.");
		check(array.incrementAndReturnUnincremented(a, 253)==2,
				"Long-array increments must return the previous minimum count.");
		check(array.read(a)==255, "Long-array increments must saturate 8-bit cells.");
		check(array.read(b)==0, "Saturating one long-array key must not leak into its neighbor.");
	}

	/** Verifies the index counter refuses non-presence counters before touching Data. */
	private static void testIndexCounterAssertionGuard(){
		boolean enabled=false;
		assert(enabled=true);
		if(!enabled){return;}
		try{
			new IndexCounter(31, true).makeKcaFromIndex(1024, 2, 2);
			check(false, "IndexCounter must assert on non-presence cell widths.");
		}catch(AssertionError expected){
			check(expected.getMessage().contains("bits=1"), "The IndexCounter assertion should explain the bits=1 limit.");
		}
	}

	/** Unit failures remain fatal without relying on the launcher's assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs a visible increment sequence in one static lock mode. */
	private static int[] exerciseLockMode(final boolean locked){
		KCountArray.SET_LOCKED_INCREMENT=true;
		KCountArray.LOCKED_INCREMENT=locked;
		final KCountArray7MTA array=new KCountArray7MTA(10007, 8, 3, null, 0);
		final long key=17;
		final long key2=KCountArray.secondaryHash(91);
		final int[] results=new int[7];
		results[0]=array.incrementAndReturnUnincremented(key, 2);
		results[1]=array.read(key);
		array.increment(key, 3);
		results[2]=array.incrementAndReturnUnincremented(key, 250);
		results[3]=array.read(key);
		results[4]=array.incrementAndReturnUnincremented(91, key2, 4);
		results[5]=array.incrementAndReturnUnincremented(91, key2, 2);
		results[6]=array.read(91, key2);
		return results;
	}

	/** Returns a cellsUsed failure for the current matrix, or null when the scan is correct. */
	private static String signedScanFailure(final Object array, final int threshold) throws Exception{
		final Method cellsUsed=array.getClass().getMethod("cellsUsed", int.class);
		final long count=((Long)cellsUsed.invoke(array, threshold)).longValue();
		return count==1 ? null : Long.toString(count);
	}

	/** Places an exact packed word directly into the first int[] row to isolate scanner behavior. */
	private static void seedPackedWord(final Object array, final int word) throws Exception{
		final Field matrixField=array.getClass().getDeclaredField("matrix");
		matrixField.setAccessible(true);
		final int[][] matrix=(int[][])matrixField.get(array);
		matrix[0]=new int[]{word};
	}

	/** Opens the optional raw lane/count evidence table for the precise-neighbor fixture. */
	private static BufferedWriter preciseNeighborWriter(final String outdir) throws Exception{
		if(outdir==null){return null;}
		final File dir=new File(outdir);
		if(!dir.exists() && !dir.mkdirs()){
			throw new RuntimeException("Could not create output directory "+outdir);
		}
		final File file=new File(dir, "precise-neighbor-sums.tsv");
		if(file.exists()){
			throw new RuntimeException("Refusing to overwrite evidence file "+file);
		}
		final BufferedWriter writer=new BufferedWriter(new FileWriter(file));
		writer.write("case\tlabel\traw_key\tlane0\texpected\tactual\n");
		return writer;
	}

	/** Appends one case failure while continuing the remaining independent cases. */
	private static void appendFailure(final StringBuilder failures, final String label,
			final String failure){
		if(failure==null){return;}
		if(failures.length()>0){failures.append("; ");}
		failures.append(label).append(": ").append(failure);
	}

	/** Finds a small key whose center, left and right neighbors occupy distinct cells. */
	private static long findPreciseNeighborKey(final KCountArray7MTA array, final int k) throws Exception{
		final long limit=1L<<(2*k);
		for(long key=0; key<limit; key++){
			if(isCollisionFreeNeighborhood(array, key, k)){return key;}
		}
		throw new AssertionError("Could not find a collision-free k="+k+" neighbor set.");
	}

	/** Returns true when all nine raw keys target distinct physical cells in this table. */
	private static boolean isCollisionFreeNeighborhood(final KCountArray7MTA array,
			final long key, final int k) throws Exception{
		final long[] rawKeys=preciseNeighborKeys(key, k);
		final long[] lanes=new long[rawKeys.length];
		for(int i=0; i<rawKeys.length; i++){
			lanes[i]=physicalLane(array, rawKeys[i]);
			for(int j=0; j<i; j++){
				if(lanes[i]==lanes[j]){return false;}
			}
		}
		return true;
	}

	/** Builds a fresh one-hash 32-bit table for independent neighbor-sum cases. */
	private static KCountArray7MTA preciseNeighborArray(){
		return new KCountArray7MTA(10007, 32, 1, null, 0);
	}

	/** Pins KCountArray7MTA's mask seed stream for deterministic fixture-local tables. */
	private static long[] resetPreciseNeighborSeeds() throws Exception{
		final long[] oldSeeds={COUNTER.getLong(null), SEEDMASK.getLong(null)};
		COUNTER.setLong(null, 0);
		SEEDMASK.setLong(null, 0);
		return oldSeeds;
	}

	/** Restores the process-wide mask seed stream after the fixture exits. */
	private static void restorePreciseNeighborSeeds(final long[] oldSeeds) throws Exception{
		COUNTER.setLong(null, oldSeeds[0]);
		SEEDMASK.setLong(null, oldSeeds[1]);
	}

	/** Returns center, the four left neighbors and the four right neighbors. */
	private static long[] preciseNeighborKeys(final long key, final int k){
		final long[] rawKeys=new long[9];
		final long leftBase=key>>>2;
		final int shift=2*(k-1);
		final long mask=(k>=32 ? -1L : ~((-1L)<<(2*k)));
		final long rightBase=(key<<2)&mask;
		rawKeys[0]=key;
		for(int i=0; i<4; i++){
			rawKeys[1+i]=leftBase|(((long)i)<<shift);
			rawKeys[5+i]=rightBase|i;
		}
		return rawKeys;
	}

	/** Seeds one precise-neighbor case and checks raw counts plus smoothed results. */
	private static void checkPreciseNeighborCase(final StringBuilder failures,
			final String label, final int centerCount, final int[] leftCounts,
			final int[] rightCounts, final int expectedLeft, final int expectedRight,
			final int expectedPrecise, final int expectedPreciseMin,
			final BufferedWriter writer) throws Exception{
		try{
			final int k=4;
			final KCountArray7MTA array=preciseNeighborArray();
			final long key=findPreciseNeighborKey(array, k);
			final long[] rawKeys=preciseNeighborKeys(key, k);
			raiseRawKey(array, rawKeys[0], centerCount);
			for(int i=0; i<4; i++){
				raiseRawKey(array, rawKeys[1+i], leftCounts[i]);
				raiseRawKey(array, rawKeys[5+i], rightCounts[i]);
			}
			String failure=checkRawCount(array, label, "center", rawKeys[0], centerCount, writer);
			for(int i=0; i<4; i++){
				failure=joinFailure(failure,
						checkRawCount(array, label, "left"+i, rawKeys[1+i], leftCounts[i], writer));
				failure=joinFailure(failure,
						checkRawCount(array, label, "right"+i, rawKeys[5+i], rightCounts[i], writer));
			}
			failure=joinFailure(failure, checkOutput(array, label, "readLeft",
					key, array.readLeft(key, k, false), expectedLeft, writer));
			failure=joinFailure(failure, checkOutput(array, label, "readRight",
					key, array.readRight(key, k, false), expectedRight, writer));
			failure=joinFailure(failure, checkOutput(array, label, "readPrecise",
					key, array.readPrecise(key, k, false), expectedPrecise, writer));
			failure=joinFailure(failure, checkOutput(array, label, "readPreciseMin",
					key, array.readPreciseMin(key, k, false), expectedPreciseMin, writer));
			appendFailure(failures, label, failure);
		}catch(Throwable e){
			appendFailure(failures, label, e.toString());
		}
	}

	/** Checks one raw key and records its lane input/count when requested. */
	private static String checkRawCount(final KCountArray7MTA array, final String caseName,
			final String label, final long rawKey, final int expected,
			final BufferedWriter writer) throws Exception{
		final long lane=array.hashForLane(rawKey, 0, 0);
		final int actual=array.read(rawKey);
		writeEvidence(writer, caseName, label, rawKey, lane, expected, actual);
		return actual==expected ? null : label+" expected "+expected+" observed "+actual;
	}

	/** Checks one visible output and records its actual and expected value. */
	private static String checkOutput(final KCountArray7MTA array, final String caseName,
			final String label, final long rawKey, final int actual, final int expected,
			final BufferedWriter writer) throws Exception{
		writeEvidence(writer, caseName, label, rawKey, physicalLane(array, rawKey), expected, actual);
		return actual==expected ? null : label+" expected "+expected+" observed "+actual;
	}

	/** Joins a new field-level failure onto an existing case failure. */
	private static String joinFailure(final String oldFailure, final String newFailure){
		if(newFailure==null){return oldFailure;}
		return oldFailure==null ? newFailure : oldFailure+", "+newFailure;
	}

	/** Writes one optional TSV evidence row. */
	private static void writeEvidence(final BufferedWriter writer, final String caseName,
			final String label, final long rawKey, final long lane,
			final int expected, final int actual) throws Exception{
		if(writer!=null){
			writer.write(caseName+'\t'+label+'\t'+rawKey+'\t'+lane+'\t'+expected+'\t'+actual+'\n');
		}
	}

	/** Returns a physical cell identifier for a one-lane raw key in this table. */
	private static long physicalLane(final KCountArray7MTA array, final long rawKey) throws Exception{
		long key=array.hashForLane(rawKey, 0, 0);
		final long arrayMask=((Number)ARRAY_MASK.get(array)).longValue();
		final int arrayBits=((Number)ARRAY_BITS.get(array)).intValue();
		final long cellMod=((Number)CELL_MOD.get(array)).longValue();
		final int indexShift=((Number)INDEX_SHIFT.get(array)).intValue();
		final int cellBits=((Number)CELL_BITS.get(array)).intValue();
		final long cellMask=((Number)CELL_MASK.get(array)).longValue();
		final long num=key&arrayMask;
		key=(key>>>arrayBits)%cellMod;
		final long index=key>>>indexShift;
		final long cellShift=cellBits*(key&cellMask);
		return (num<<48)|(index<<6)|cellShift;
	}

	/** Raises one single-hash raw key to the requested value. */
	private static void raiseRawKey(final KCountArray7MTA array, final long key, final int value) throws Exception{
		final long hashed=array.hashForLane(key, 0, 0);
		try{
			RAISE_LANE.invoke(array, hashed, value);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Seeds every lane for a dual-hash key to a high value via finite reflective setup. */
	private static void raiseAllLanes(final KCountArray7MTA array, final long key,
			final long key2, final int value, final int hashes) throws Exception{
		for(int lane=0; lane<hashes; lane++){
			final long hashed=((Long)HASH_FOR_LANE.invoke(array, key, key2, lane)).longValue();
			try{
				RAISE_LANE.invoke(array, hashed, value);
			}catch(InvocationTargetException e){
				final Throwable cause=e.getCause();
				if(cause instanceof Exception){throw (Exception)cause;}
				if(cause instanceof Error){throw (Error)cause;}
				throw new RuntimeException(cause);
			}
		}
	}

	private static final Method HASH_FOR_LANE, RAISE_LANE;
	private static final Field ARRAY_MASK, ARRAY_BITS, CELL_MOD, INDEX_SHIFT, CELL_BITS, CELL_MASK,
			COUNTER, SEEDMASK;
	static{
		try{
			HASH_FOR_LANE=KCountArray7MTA.class.getDeclaredMethod(
					"hashForLane", long.class, long.class, int.class);
			RAISE_LANE=KCountArray7MTA.class.getDeclaredMethod(
					"incrementHashedLocalAndReturnUnincremented_toAtLeast", long.class, int.class);
			ARRAY_MASK=KCountArray.class.getDeclaredField("arrayMask");
			ARRAY_BITS=KCountArray.class.getDeclaredField("arrayBits");
			CELL_MOD=KCountArray7MTA.class.getDeclaredField("cellMod");
			INDEX_SHIFT=KCountArray.class.getDeclaredField("indexShift");
			CELL_BITS=KCountArray.class.getDeclaredField("cellBits");
			CELL_MASK=KCountArray.class.getDeclaredField("cellMask");
			COUNTER=KCountArray7MTA.class.getDeclaredField("counter");
			SEEDMASK=KCountArray7MTA.class.getDeclaredField("SEEDMASK");
			HASH_FOR_LANE.setAccessible(true);
			RAISE_LANE.setAccessible(true);
			ARRAY_MASK.setAccessible(true);
			ARRAY_BITS.setAccessible(true);
			CELL_MOD.setAccessible(true);
			INDEX_SHIFT.setAccessible(true);
			CELL_BITS.setAccessible(true);
			CELL_MASK.setAccessible(true);
			COUNTER.setAccessible(true);
			SEEDMASK.setAccessible(true);
		}catch(Exception e){throw new RuntimeException(e);}
	}

	private static int checks;
}
