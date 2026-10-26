package test.map;

import java.util.Arrays;
import map.IntHashSetList;

/** Regression fixture for IntHashSetList sparse clear dispatch.
 * @author Sucrose
 */
public final class IntHashSetListClearTest{

	public static void main(String[] args){
		check(args.length==1, "Usage: expect=baseline|fixed");
		String[] split=args[0].split("=", 2);
		check(split.length==2 && split[0].equals("expect"), "Expected expect=baseline|fixed");
		boolean fixed;
		if(split[1].equals("baseline")){fixed=false;}
		else if(split[1].equals("fixed")){fixed=true;}
		else{throw new AssertionError("Expected expect=baseline|fixed, got "+args[0]);}

		checkSparseClear(fixed);
		checkDenseClear();
		checkDuplicateTracking();
		checkReuseAfterClear();
		checkRemoveStillRejected();
		System.err.println("PASS IntHashSetListClearTest expect="+split[1]);
	}

	private static void checkSparseClear(boolean fixed){
		IntHashSetList set=sparseSet();
		addFour(set);
		if(!fixed){
			expectUnsupported(new Runnable(){
				@Override
				public void run(){set.clear();}
			}, "baseline sparse clear should call the public remove override");
			check(set.size()==4, "throwing sparse clear should leave size unchanged");
			checkAllPresent(set);
			return;
		}

		set.clear();
		check(set.isEmpty(), "fixed sparse clear should empty the set");
		check(set.toArray().length==0, "fixed sparse clear should reset the tracked list");
		check(set.verify(), "fixed sparse clear should preserve set/list invariants");
	}

	private static void checkDenseClear(){
		IntHashSetList set=denseSet();
		addFour(set);
		check(set.size()>0.25*set.sizeLimit(), "dense clear fixture should use dense clear");
		set.clear();
		check(set.isEmpty(), "dense clear should empty the set");
		check(set.toArray().length==0, "dense clear should reset the tracked list");
		check(set.verify(), "dense clear should preserve set/list invariants");
	}

	private static void checkDuplicateTracking(){
		IntHashSetList set=denseSet();
		check(set.add(7), "first add should succeed");
		check(!set.add(7), "duplicate add should be rejected");
		check(set.size()==1, "duplicate add should not grow the set");
		int[] keys=set.toArray();
		check(keys.length==1 && keys[0]==7, "duplicate add should not duplicate the tracked key: "+
			Arrays.toString(keys));
		check(set.verify(), "duplicate add should preserve set/list invariants");
	}

	private static void checkReuseAfterClear(){
		IntHashSetList set=denseSet();
		addFour(set);
		set.clear();
		check(set.add(99), "add after dense clear should succeed");
		check(set.contains(99), "reused set should contain the new key");
		check(set.size()==1 && set.toArray().length==1, "reused set should track only the new key");
		check(set.verify(), "reuse should preserve set/list invariants");
	}

	private static void checkRemoveStillRejected(){
		IntHashSetList set=denseSet();
		check(set.add(11), "add before remove rejection should succeed");
		expectUnsupported(new Runnable(){
			@Override
			public void run(){set.remove(11);}
		}, "external remove should remain unsupported");
		check(set.size()==1 && set.contains(11), "rejected remove should leave the set unchanged");
	}

	private static IntHashSetList sparseSet(){
		return new IntHashSetList(256);
	}

	private static IntHashSetList denseSet(){
		return new IntHashSetList(8);
	}

	private static void addFour(IntHashSetList set){
		check(set.add(0), "add 0");
		check(set.add(1), "add 1");
		check(set.add(2), "add 2");
		check(set.add(3), "add 3");
		check(set.size()==4, "expected four distinct keys");
		checkAllPresent(set);
		check(set.verify(), "added keys should preserve set/list invariants");
	}

	private static void checkAllPresent(IntHashSetList set){
		for(int i=0; i<4; i++){
			check(set.contains(i), "missing key "+i);
		}
	}

	private static void expectUnsupported(Runnable runnable, String message){
		try{
			runnable.run();
		}catch(UnsupportedOperationException expected){
			return;
		}
		throw new AssertionError(message);
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private IntHashSetListClearTest(){
		throw new AssertionError();
	}
}
