package test.kmer;

import java.lang.reflect.Field;

import kmer.HashArray1D;
import kmer.ScheduleMaker;

/** Finite scheduler contracts without allocating the planned table fleet.
 * Runs one case in a fresh JVM with a fixed heap, normally -Xmx256m.
 * @author Shinobu */
public final class ScheduleContractTest{

	/** Checks explicit reservations, actual probe-tail budgeting, or nearby controls. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected case=override|auto|budget|controls|invalid");}
		final String mode=args[0];
		final int extra=new HashArray1D(new int[] {5}, -1L).arrayLength()-5;
		check(extra==60, "Fixture expects the actual table's 60-cell probe tail: "+extra);
		int cases=0;
		if(mode.equals("case=override")){
			for(int passes : new int[] {1, 2}){
				for(boolean prealloc : new boolean[] {false, true}){
					final long reserved=16L*1024*1024;
					final ScheduleMaker maker=new ScheduleMaker(31, 12, prealloc, 1, 5, passes, 0.2, reserved);
					check(field(maker, "filterMemory0")==reserved && field(maker, "filterMemory1")==reserved,
						"Explicit prefilter reservation must reach both pass budgets: requested="+reserved+
						", first="+field(maker, "filterMemory0")+", second="+field(maker, "filterMemory1"));
					checkBudget(maker, 12, extra);
					cases++;
				}
			}
		}else if(mode.equals("case=auto")){
			final ScheduleMaker maker=new ScheduleMaker(31, 12, false, 1, 5, -1, 0.2, 0);
			check(field(maker, "filterMemory0")>0 && field(maker, "filterMemory1")>0,
				"Automatic prefilter passes remain enabled and must reserve both alternating filter budgets");
			checkBudget(maker, 12, extra);
			cases++;
		}else if(mode.equals("case=budget")){
			// Many small shards expose fixed probe-tail costs hidden by large default tables.
			final ScheduleMaker maker=new ScheduleMaker(3000, 41, false, 1, 5, 0, 0, 0);
			checkBudget(maker, 41, extra);
			cases++;
		}else if(mode.equals("case=controls")){
			for(int ways : new int[] {31, 127}){
				for(int bytes : new int[] {12, 41}){
					for(boolean prealloc : new boolean[] {false, true}){
						final ScheduleMaker maker=new ScheduleMaker(ways, bytes, prealloc, 1, 5, 0, 0, 16L*1024*1024);
						check(field(maker, "filterMemory0")==0 && field(maker, "filterMemory1")==0,
							"Disabled prefiltering must not reserve either filter pass");
						checkBudget(maker, bytes, extra);
						cases++;
					}
				}
			}
		}else if(mode.equals("case=invalid")){
			requireRejected(300000, 41, 1);
			requireRejected(31, 12, Double.NaN);
			requireRejected(31, 12, Double.POSITIVE_INFINITY);
			cases=3;
		}else{throw new IllegalArgumentException("Unknown case: "+mode);}
		System.out.println("PASS "+mode+" cases="+cases+" probeTail="+extra);
	}

	/** Requires a diagnostic rejection without allocating any planned arrays. */
	private static void requireRejected(final int ways, final int bytes, final double ratio){
		try{
			new ScheduleMaker(ways, bytes, false, ratio, 5, 0, 0, 0);
		}catch(IllegalArgumentException | AssertionError failure){
			final String message=failure.getMessage();
			check(message!=null && (message.contains("probe-tail reservations") || message.contains("maxSize0=")),
				"Expected a capacity-budget diagnostic, not an unrelated failure: "+failure);
			return;
		}
		throw new AssertionError("Invalid capacity request was accepted: ways="+ways+", bytes="+bytes+", ratio="+ratio);
	}

	/** Uses the actual ladder and measured tail size, rather than the scheduler's tail assumption. */
	private static void checkBudget(final ScheduleMaker maker, final int bytes, final int extra) throws Exception{
		final int[] ladder=maker.makeSchedule();
		check(ladder.length>0, "Every table needs a first allocation");
		for(int i=0; i<ladder.length; i++){
			check(ladder[i]>0 && (i==0 || ladder[i]>ladder[i-1]), "Resize capacities must increase strictly");
		}
		final long last=ladder[ladder.length-1]+(long)extra;
		final long old=ladder.length>1 ? ladder[ladder.length-2]+(long)extra : 0;
		final long peak=(last*maker.ways+old*maker.maxConcurrentResizes)*bytes;
		final long budget=field(maker, "tableMemory");
		check(peak<=budget, "Actual primary arrays plus simultaneous old copies exceed the reserved budget: peak="+
			peak+", budget="+budget+", last="+last+", old="+old+", ways="+maker.ways+", bytes="+bytes);
	}

	/** Observes package-private reservations without changing production visibility for a baseline. */
	private static long field(final ScheduleMaker maker, final String name) throws Exception{
		final Field field=ScheduleMaker.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(maker);
	}

	/** Keeps the oracle active in both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
