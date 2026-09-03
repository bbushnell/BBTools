package align2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import stream.Read;

/**
 * Standalone driver testing BBMapS.Coordinator in isolation — the two
 * pieces of BBMapS flagged as "completely untested": (1) out-of-order
 * reassembly (t=1 never engages it; results arrive in submission order
 * with one worker) and (2) the design doc §8 "crash loud, never hang"
 * error path (forceError()/finishError() on a mid-pipeline failure,
 * Yaoyao's safety-critical finding). Does not touch BBIndex/Data/a real
 * Streamer — Coordinator's own logic has no dependency on them (verified
 * by reading BBMapS.java: the only construction-time BBMapS state a
 * Coordinator instance uses is its own constructor args).
 *
 * @author Nowi
 * @date 2026-09-03
 */
public class BBMapSCoordinatorUnitTest {

	static final class RecordingRouteWriter implements BBMapS.RouteWriter {
		final String name;
		final long throwOnId; //-1 disables
		final List<Long> submittedIds=Collections.synchronizedList(new ArrayList<Long>());
		volatile int finishCalls=0, finishErrorCalls=0;
		RecordingRouteWriter(String name){this(name, -1);}
		RecordingRouteWriter(String name, long throwOnId){this.name=name; this.throwOnId=throwOnId;}
		public void submit(long id, ArrayList<Read> reads){
			if(reads==null){throw new AssertionError(name+": submit() received null reads for id="+id+" — violates the never-null contract.");}
			if(id==throwOnId){throw new RuntimeException("Simulated writer failure on route '"+name+"' at id="+id+" (test-injected).");}
			submittedIds.add(id);
		}
		public void finish(){finishCalls++;}
		public void finishError(){finishErrorCalls++;}
	}

	static final class RecordingBBSplitterInvoker implements BBMapS.BBSplitterInvoker {
		final List<Long> invokedIds=Collections.synchronizedList(new ArrayList<Long>());
		public void invoke(long id, ArrayList<Read> pristinePrimaryList){invokedIds.add(id);}
	}

	static ArrayList<Read> emptyList(){return new ArrayList<Read>();}
	static BBMapS.MapResult newResult(long id){
		return new BBMapS.MapResult(id, emptyList(), emptyList(), emptyList(), emptyList());
	}

	public static void main(String[] args) throws Exception {
		BBMapS mapper=new BBMapS(new String[]{"build=1"});
		testOrdering(mapper);
		testErrorPath(mapper);
		System.out.println("ALL TESTS PASSED.");
	}

	/**
	 * Feeds MapResults in a deliberately scrambled arrival order and asserts
	 * every RouteWriter/BBSplitterInvoker receives them in strictly
	 * ascending id order regardless, and finish() (not finishError()) is
	 * called exactly once per route on the clean-shutdown path.
	 */
	static void testOrdering(BBMapS mapper) throws Exception {
		final BlockingQueue<BBMapS.MapResult> resultQueue=new ArrayBlockingQueue<BBMapS.MapResult>(16);
		final RecordingRouteWriter routeMapped=new RecordingRouteWriter("mapped");
		final RecordingRouteWriter routeUnmapped=new RecordingRouteWriter("unmapped");
		final RecordingRouteWriter routeBlacklisted=new RecordingRouteWriter("blacklisted");
		final RecordingRouteWriter routePrimary=new RecordingRouteWriter("primary");
		final RecordingBBSplitterInvoker splitter=new RecordingBBSplitterInvoker();

		final int numWorkers=2;
		BBMapS.Coordinator coordinator=mapper.new Coordinator(resultQueue, numWorkers,
				routeMapped, routeUnmapped, routeBlacklisted, routePrimary,
				splitter, false, false);

		Thread coordinatorThread=new Thread(coordinator, "test-coordinator-ordering");
		coordinatorThread.start();

		//Scrambled arrival: 2, 0, 4, 1, poison(A), 3, poison(B). Correct behavior: id=2 buffers
		//until id=0 arrives (drains 0, still buffers 2 since 1 is missing); id=4 buffers; id=1
		//arrives -> drains 1, then drains the already-buffered 2; 3 is still missing so 4 stays
		//buffered; id=3 arrives -> drains 3, then drains 4. Final emitted order must be 0,1,2,3,4
		//-- NOT the arrival order 2,0,4,1,3.
		resultQueue.put(newResult(2));
		resultQueue.put(newResult(0));
		resultQueue.put(newResult(4));
		resultQueue.put(newResult(1));
		resultQueue.put(BBMapS.MapResult.poison(1000));
		resultQueue.put(newResult(3));
		resultQueue.put(BBMapS.MapResult.poison(1001));

		coordinatorThread.join(10000);
		if(coordinatorThread.isAlive()){
			throw new AssertionError("FAIL testOrdering: coordinator thread did not terminate within 10s -- likely deadlocked on the reordering logic.");
		}

		List<Long> expected=java.util.Arrays.asList(0L, 1L, 2L, 3L, 4L);
		for(RecordingRouteWriter rw : new RecordingRouteWriter[]{routeMapped, routeUnmapped, routeBlacklisted, routePrimary}){
			if(!rw.submittedIds.equals(expected)){
				throw new AssertionError("FAIL testOrdering: route '"+rw.name+"' received ids in order "+rw.submittedIds+", expected strictly ascending "+expected);
			}
			if(rw.finishCalls!=1){
				throw new AssertionError("FAIL testOrdering: route '"+rw.name+"' finish() called "+rw.finishCalls+" times, expected exactly 1.");
			}
			if(rw.finishErrorCalls!=0){
				throw new AssertionError("FAIL testOrdering: route '"+rw.name+"' finishError() called "+rw.finishErrorCalls+" times, expected 0 (clean shutdown path).");
			}
		}
		if(!splitter.invokedIds.equals(expected)){
			throw new AssertionError("FAIL testOrdering: BBSplitterInvoker received ids in order "+splitter.invokedIds+", expected strictly ascending "+expected);
		}
		if(coordinator.errorState){
			throw new AssertionError("FAIL testOrdering: coordinator.errorState==true, error="+coordinator.error);
		}

		System.out.println("PASS testOrdering: Coordinator drained ids in strictly ascending order "+expected
			+" despite scrambled arrival (2,0,4,1,poison,3,poison); finish() called once per route, "
			+"zero finishError() calls, BBSplitterInvoker saw the same ascending order.");
	}

	/**
	 * Injects a writer-side failure (routeMapped.submit() throws on id=1)
	 * and asserts the design doc §8 "crash loud, never hang" contract:
	 * coordinator.errorState becomes true, forceError() reaches ALL FOUR
	 * routes' finishError() (not just the one that threw), finish() is
	 * NEVER called on the error path, and the coordinator thread still
	 * terminates promptly (does not hang).
	 */
	static void testErrorPath(BBMapS mapper) throws Exception {
		final BlockingQueue<BBMapS.MapResult> resultQueue=new ArrayBlockingQueue<BBMapS.MapResult>(16);
		final RecordingRouteWriter routeMapped=new RecordingRouteWriter("mapped", 1); //throws on id=1
		final RecordingRouteWriter routeUnmapped=new RecordingRouteWriter("unmapped");
		final RecordingRouteWriter routeBlacklisted=new RecordingRouteWriter("blacklisted");
		final RecordingRouteWriter routePrimary=new RecordingRouteWriter("primary");
		final RecordingBBSplitterInvoker splitter=new RecordingBBSplitterInvoker();

		final int numWorkers=1;
		BBMapS.Coordinator coordinator=mapper.new Coordinator(resultQueue, numWorkers,
				routeMapped, routeUnmapped, routeBlacklisted, routePrimary,
				splitter, false, false);

		Thread coordinatorThread=new Thread(coordinator, "test-coordinator-error");
		coordinatorThread.start();

		resultQueue.put(newResult(0)); //succeeds on all routes
		resultQueue.put(newResult(1)); //routeMapped.submit(1,...) throws here

		coordinatorThread.join(10000);
		if(coordinatorThread.isAlive()){
			throw new AssertionError("FAIL testErrorPath: coordinator thread did not terminate within 10s after a simulated writer failure -- this IS the hang the §8 contract exists to prevent.");
		}

		if(!coordinator.errorState){
			throw new AssertionError("FAIL testErrorPath: coordinator.errorState==false after a simulated writer failure; expected true.");
		}
		for(RecordingRouteWriter rw : new RecordingRouteWriter[]{routeMapped, routeUnmapped, routeBlacklisted, routePrimary}){
			if(rw.finishErrorCalls!=1){
				throw new AssertionError("FAIL testErrorPath: route '"+rw.name+"' finishError() called "+rw.finishErrorCalls+" times, expected exactly 1 (forceError() must reach every route, not just the one that threw).");
			}
			if(rw.finishCalls!=0){
				throw new AssertionError("FAIL testErrorPath: route '"+rw.name+"' finish() called "+rw.finishCalls+" times on the error path; expected 0 -- normal finish() must not run after a failure.");
			}
		}
		//id=0 must still have gone through cleanly on every route before the id=1 failure.
		for(RecordingRouteWriter rw : new RecordingRouteWriter[]{routeUnmapped, routeBlacklisted, routePrimary}){
			if(!rw.submittedIds.equals(java.util.Arrays.asList(0L))){
				throw new AssertionError("FAIL testErrorPath: route '"+rw.name+"' submittedIds="+rw.submittedIds+", expected [0] (id=0 should have succeeded before id=1 failed).");
			}
		}

		System.out.println("PASS testErrorPath: a simulated writer failure on id=1 set errorState=true, "
			+"reached finishError() on all 4 routes exactly once (not just the failing route), never called "
			+"finish(), and the coordinator thread terminated promptly instead of hanging.");
	}

}
