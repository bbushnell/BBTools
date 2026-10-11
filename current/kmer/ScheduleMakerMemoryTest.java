package kmer;

import java.util.Arrays;

/** Regression checks for bounded resize headroom without allocating large arrays.
 * @author Brian Bushnell, Fischl */
public class ScheduleMakerMemoryTest {

	/** Exercises old-constructor parity, conservative peak budgets and bounded gains. */
	public static void main(String[] args){
		boolean enabled=false;
		assert(enabled=true);
		if(!enabled){throw new IllegalStateException("Run schedule regression with -ea.");}
		int cases=0;
		for(int ways : new int[] {1, 31, 127}){
			for(int bytes : new int[] {13, 24, 41}){
				for(int initial : new int[] {5, 128000}){
					for(boolean prealloc : new boolean[] {false, true}){
						for(int passes : new int[] {0, 1, 2}){
							final double filter=(passes==0 ? 0 : 0.2);
							ScheduleMaker legacy=new ScheduleMaker(ways, bytes, prealloc, 1, initial, passes, filter, 0);
							ScheduleMaker allWays=new ScheduleMaker(ways, bytes, prealloc, 1, initial, passes, filter, 0, ways);
							assert(Arrays.equals(legacy.makeSchedule(), allWays.makeSchedule())) :
								"An explicit all-shard bound must preserve the legacy schedule.";
							for(int writers : new int[] {1, 8, 31, 256}){
								ScheduleMaker bounded=new ScheduleMaker(ways, bytes, prealloc, 1, initial, passes, filter, 0, writers);
								checkBudget(bounded);
								assert(bounded.maxSize>=legacy.maxSize) : "A lower simultaneous-copy count cannot reduce capacity.";
								if(prealloc || writers>=ways){
									assert(Arrays.equals(bounded.makeSchedule(), legacy.makeSchedule())) :
										"Preallocation and all-shard fanout must remain unchanged.";
								}
								cases++;
							}
						}
					}
				}
			}
		}
		ScheduleMaker old=new ScheduleMaker(31, 24, false, 1, 5);
		ScheduleMaker bounded=new ScheduleMaker(31, 24, false, 1, 5, 0, 0, 0, 8);
		assert(bounded.maxSize>old.maxSize) : "Eight writers across31shards should recover reserved capacity.";
		System.out.println(old);
		System.out.println(bounded);
		System.out.println("SCHEDULE_MEMORY_TEST_PASS cases="+cases);
	}

	/** Independently checks the final-resize old-array bound using the actual ladder. */
	private static void checkBudget(final ScheduleMaker maker){
		final int[] schedule=maker.makeSchedule();
		assert(schedule.length>0) : "A table must have an initial allocation.";
		final long previous=(schedule.length>1 ? schedule[schedule.length-2]+20L : 0);
		final double finalBytes=(schedule[schedule.length-1]+20.0)*maker.ways*maker.bytesPerKmer;
		final double copies=previous*(double)maker.maxConcurrentResizes*maker.bytesPerKmer;
		assert(finalBytes+copies<=maker.tableMemory*maker.memRatio) :
			"Scheduled arrays and live old copies exceed the declared table-memory budget: "+maker;
	}
}
