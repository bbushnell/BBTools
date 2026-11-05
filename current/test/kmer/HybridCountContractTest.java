package test.kmer;

import java.util.Arrays;

import kmer.HashArray;
import kmer.HashArrayHybrid;
import kmer.HashArrayHybridFast;

/**
 * Checks scalar counting in both hybrid storage tiers, including victim migration.
 * Separate set/get controls preserve reference-ID lists without mixing them with counts.
 * Each case runs in a fresh JVM under -ea or -da; literal checks remain unconditional.
 * @author Sucrose
 */
public final class HybridCountContractTest{

	/** Accepts case=hybrid-|fast- followed by one of the twelve named scenarios. */
	public static void main(String[] args){
		check(args.length==1 && args[0].startsWith("case="), "Expected exactly one case= argument");
		final String name=args[0].substring(5);
		final boolean fast=name.startsWith("fast-");
		check(fast || name.startsWith("hybrid-"), "Expected hybrid- or fast- prefix: "+name);
		final String mode=name.substring(fast ? 5 : 7);
		if(mode.equals("primary-saturate-increment") || mode.equals("victim-saturate-increment") ||
				mode.equals("primary-saturate-created") || mode.equals("victim-saturate-created")){
			verifySaturation(fast, name, mode.startsWith("victim-"), mode.endsWith("created"));
			System.out.println("PASS case="+name);
			return;
		}
		final boolean resize, association, created;
		final long target;
		switch(mode){
			case "primary-increment": resize=false; association=false; created=false; target=0; break;
			case "victim-increment": resize=false; association=false; created=false; target=VICTIM_KEY; break;
			case "primary-created": resize=false; association=false; created=true; target=0; break;
			case "victim-created": resize=false; association=false; created=true; target=VICTIM_KEY; break;
			case "resize-increment": resize=true; association=false; created=false; target=VICTIM_KEY; break;
			case "resize-created": resize=true; association=false; created=true; target=VICTIM_KEY; break;
			case "association-control": resize=false; association=true; created=false; target=VICTIM_KEY; break;
			case "association-resize": resize=true; association=true; created=false; target=VICTIM_KEY; break;
			default: throw new IllegalArgumentException("Unknown case: "+name);
		}
		final int[] schedule=resize ? new int[] {101, 211} : new int[] {101};
		final HashArray table=fast ? new HashArrayHybridFast(schedule, -1L) : new HashArrayHybrid(schedule, -1L);
		for(int i=0; i<61; i++){
			final long key=101L*i;
			final int result=association ? table.set(key, 5) : addCount(table, key, created);
			check(result==1, "Each initial key must be new and have initial count 1 where counted: "+key);
			check(table.getValue(key)==(association ? 5 : 1), "Initial scalar value differs at key "+key);
		}
		verifyTopology(table, false);
		check(table.victims().getValue(VICTIM_KEY)==(association ? 5 : 1),
			"Key 6060 must begin in the victim forest");
		if(association){
			check(table.set(0, 7)==0 && table.set(VICTIM_KEY, 7)==0,
				"Adding reference IDs must not create new keys");
			verifyList(table, 0);
			verifyList(table, VICTIM_KEY);
			if(resize){grow(table, true, false);}
			verifyTopology(table, resize);
			verifyList(table, 0);
			verifyList(table, VICTIM_KEY);
			verifyOtherKeys(table, resize, true, target);
			System.out.println("case="+name+", primaryIDs=[5, 7], formerVictimIDs=[5, 7], slots="+table.arrayLength());
		}else{
			final int updateReturn=addCount(table, target, created);
			//Do not check the bad baseline return or count before migration: both must reach the resize oracle.
			if(resize){grow(table, false, created);}
			verifyTopology(table, resize);
			final int count=table.getValue(target);
			final int[] singleton={-777}, actual=table.getValues(target, singleton);
			System.out.println("case="+name+", updateReturn="+updateReturn+", count="+count+
				", values="+Arrays.toString(actual)+", slots="+table.arrayLength());
			check(count==2, "KMR021: an existing scalar count 1 plus 1 must become 2; case="+name);
			check(updateReturn==(created ? 0 : 2), "Count API returned the wrong result; case="+name);
			check(actual!=null && actual.length>0 && actual[0]==2, "Scalar payload must contain count 2");
			for(int i=1; i<actual.length; i++){
				check(actual[i]<0, "Scalar counting must retain one logical value; case="+name);
			}
			if(resize || target==0){
				check(actual==singleton, "Primary scalar must remain direct storage rather than a list tag; case="+name);
			}
			verifyOtherKeys(table, resize, false, target);
		}
		System.out.println("PASS case="+name);
	}

	/** Starts the target at MAX-1 through its count API, then verifies both entry into and retention of saturation. */
	private static void verifySaturation(final boolean fast, final String name,
			final boolean victim, final boolean created){
		final HashArray table=fast ? new HashArrayHybridFast(new int[] {101}, -1L) :
			new HashArrayHybrid(new int[] {101}, -1L);
		final long target=victim ? VICTIM_KEY : 0;
		for(int i=0; i<61; i++){
			final long key=101L*i;
			final int initial=key==target ? Integer.MAX_VALUE-1 : 1;
			final int result=created ? table.incrementAndReturnNumCreated(key, initial) : table.increment(key, initial);
			check(result==(created ? 1 : initial), "Initial count API result differs at key "+key);
			check(table.getValue(key)==initial, "Target must be initialized directly through the selected count API");
		}
		verifyTopology(table, false);
		check(table.victims().getValue(VICTIM_KEY)==(victim ? Integer.MAX_VALUE-1 : 1),
			"The intended overflow key must exist before saturation updates");
		for(int step=1; step<=2; step++){
			final int result=addCount(table, target, created);
			final int count=table.getValue(target);
			final int[] singleton={-777}, actual=table.getValues(target, singleton);
			System.out.println("case="+name+", saturationStep="+step+", updateReturn="+result+
				", count="+count+", values="+Arrays.toString(actual));
			check(count==Integer.MAX_VALUE, "KMR021: scalar count must reach and retain MAX; case="+name+", step="+step);
			check(result==(created ? 0 : Integer.MAX_VALUE), "Saturated count API returned the wrong result; case="+name);
			check(actual!=null && actual.length>0 && actual[0]==Integer.MAX_VALUE,
				"Saturated scalar payload must contain MAX");
			for(int i=1; i<actual.length; i++){
				check(actual[i]<0, "Saturation must not append logical values; case="+name);
			}
			if(!victim){check(actual==singleton, "Saturated primary count must remain direct scalar storage");}
			verifyTopology(table, false);
			for(int i=0; i<61; i++){
				final long key=101L*i;
				check(table.getValue(key)==(key==target ? Integer.MAX_VALUE : 1),
					"Saturation changed an unrelated key: "+key);
			}
			check(!table.contains(29), "Saturation must not create key 29");
		}
	}

	/** Selects the count-returning or creation-returning API, always adding one. */
	private static int addCount(final HashArray table, final long key, final boolean created){
		return created ? table.incrementAndReturnNumCreated(key, 1) : table.increment(key, 1);
	}

	/** Adds 28 distinct keys: total size 89 crosses the initial 88-key growth threshold once. */
	private static void grow(final HashArray table, final boolean association, final boolean created){
		for(int key=1; key<=28; key++){
			final int result=association ? table.set(key, 5) : addCount(table, key, created);
			check(result==1, "Each resize-triggering key must be new: "+key);
			check(table.arrayLength()==(key<28 ? 161 : 271), "Unexpected resize point at inserted key "+key);
		}
	}

	/** Verifies actual primary/victim occupancy, rather than trusting requested capacities. */
	private static void verifyTopology(final HashArray table, final boolean resized){
		check(table.arrayLength()==(resized ? 271 : 161), "Unexpected primary storage length");
		check(table.size()==(resized ? 89 : 60) && table.victims().size()==(resized ? 0 : 1),
			"Expected 60 primary plus one victim before growth, or 89 primary and no victims after it");
		int occupied=0;
		boolean targetPrimary=false;
		for(long key : table.array()){
			if(key>=0){occupied++;}
			if(key==VICTIM_KEY){targetPrimary=true;}
		}
		check(occupied==(resized ? 89 : 60) && targetPrimary==resized,
			"Actual primary occupancy must agree with the expected victim migration");
	}

	/** Checks literal two-ID association contents, independent of backing capacity or duplicate policy. */
	private static void verifyList(final HashArray table, final long key){
		final int[] actual=table.getValues(key, new int[1]);
		check(actual!=null && actual.length>=2 && actual[0]==5 && actual[1]==7,
			"Association control must retain IDs [5, 7] at key "+key+": "+Arrays.toString(actual));
		for(int i=2; i<actual.length; i++){
			check(actual[i]<0, "Association control gained an unexpected ID at key "+key);
		}
	}

	/** Ensures every other inserted key retains its original scalar and no extra key appears. */
	private static void verifyOtherKeys(final HashArray table, final boolean resized,
			final boolean association, final long target){
		for(int i=0; i<61; i++){
			final long key=101L*i;
			final int expected=association ? 5 : key==target ? 2 : 1;
			check(table.getValue(key)==expected, "Stored key changed unexpectedly: "+key);
		}
		if(resized){
			for(int key=1; key<=28; key++){
				check(table.getValue(key)==(association ? 5 : 1), "Resize filler key changed: "+key);
			}
		}
		check(!table.contains(29), "Key 29 was never inserted");
	}

	/** Preserves literal and topology checks under -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final long VICTIM_KEY=6060;
}
