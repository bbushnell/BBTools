package test.benchmark;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import benchmark.CallGenesLeanCalls;
import benchmark.CallGenesLeanWorkers;
import benchmark.CallGenesLengthScorecard;
import benchmark.CallGenesStopAudit;
import gff.GffLine;
import prok.Orf;
import prok.ProkObject;
import structures.ByteBuilder;

/**
 * Synthetic fixtures for compact CDS reconstruction and cohort worker ownership.
 * Compares ordinary emitted GFF fields and grading counts, input-order results,
 * and worker joins after failure or interruption. Run in the Dori validation build.
 * The fixture accepts no files and does not invoke a genome caller.
 *
 * @author Sayu
 */
public final class CallGenesLeanTest{

	/** Run representation, ordering, failure and interruption cases. */
	public static void main(String[] args)throws Exception{
		check(args.length==0,"No fixture inputs required");
		compact();
		ordered();
		failureJoins(false);
		failureJoins(true);
		System.out.println("CallGenesLeanTest PASS checks="+checks);
	}

	/** Compare compact reconstruction with the real emitted GFF, including edge partials. */
	private static void compact(){
		final byte[] bases=new byte[603];
		Arrays.fill(bases,(byte)'A');
		final ArrayList<Orf> calls=new ArrayList<Orf>();
		// Literal lengths span every complete length band and both strands.
		for(int aa:new int[]{20,40,65,100,150,200}){
			for(int strand=0;strand<2;strand++){
				calls.add(new Orf("ctg description",0,3*(aa+1)-1,strand,0,bases,false,ProkObject.CDS,false,false));
			}
		}
		for(int strand=0;strand<2;strand++){
			calls.add(new Orf("ctg\tdescription",0,602,strand,0,bases,false,ProkObject.CDS,true,true));
			calls.add(new Orf("ctg",0,122,strand,0,bases,false,ProkObject.CDS,true,true));
			calls.add(new Orf("ctg",480,602,strand,0,bases,false,ProkObject.CDS,true,true));
			calls.add(new Orf("ctg",30,153,strand,0,bases,false,ProkObject.CDS,false,false));
		}
		calls.add(calls.get(2)); // Duplicate query occurrences must survive.
		calls.add(new Orf("ctg",10,90,0,0,bases,false,ProkObject.tRNA,false,false));
		final ArrayList<GffLine> ordinary=new ArrayList<GffLine>();
		for(Orf call:calls){
			if(call.type==ProkObject.CDS){
				ordinary.add(new GffLine(call.appendGff(new ByteBuilder()).toBytes()));
			}
		}
		final CallGenesLeanCalls compact=new CallGenesLeanCalls();
		compact.append("ctg",calls);
		final ArrayList<GffLine> actual=compact.lines();
		check(compact.features()==calls.size() && actual.size()==calls.size()-1,"RNA counted separately; all CDS occurrences retained");
		final boolean[] bands=new boolean[7],partials=new boolean[4];
		for(int i=0;i<ordinary.size();i++){
			final GffLine a=actual.get(i),b=ordinary.get(i);
			check(a.seqid.equals(b.seqid) && a.type.equals(b.type) && a.start==b.start && a.stop==b.stop
				&& a.strand==b.strand && a.phase==b.phase,"Exact emitted grading fields "+i);
			check(equal(CallGenesStopAudit.attribute(a,"partial"),CallGenesStopAudit.attribute(b,"partial")),"Exact partial metadata "+i);
			final int band=CallGenesLengthScorecard.calledBand(a);
			check(band==CallGenesLengthScorecard.calledBand(b),"Emitted and compact length band "+i);
			bands[band]=true;
			partials[calls.get(i).partialMask()]=true;
		}
		for(boolean present:bands){
			check(present,"Every length band covered");
		}
		for(boolean present:partials){
			check(present,"Every partial mask covered");
		}
		// A last-wins reference map, wrong starts, duplicates and missing keys exercise grading.
		final ArrayList<GffLine> refs=new ArrayList<GffLine>(ordinary.subList(1,ordinary.size()-2));
		final CallGenesStopAudit.ReferenceSet reference=CallGenesStopAudit.fromLines(refs);
		check(Arrays.equals(CallGenesStopAudit.grade(reference,actual),CallGenesStopAudit.grade(reference,ordinary)),"All eight stop/start grading counts agree");
	}

	/** Release task zero after task one starts; results must remain input ordered. */
	private static void ordered(){
		final CountDownLatch later=new CountDownLatch(1);
		final AtomicInteger active=new AtomicInteger(),peak=new AtomicInteger();
		final ArrayList<Integer> result=CallGenesLeanWorkers.run(12,3,i->{
			final int n=active.incrementAndGet();
			peak.accumulateAndGet(n,Math::max);
			try{
				if(i==0){
					await(later);
				}else if(i==1){
					later.countDown();
				}
				return i;
			}finally{
				active.decrementAndGet();
			}
		});
		check(result.size()==12 && active.get()==0 && peak.get()<=3,"Cohort and worker budget conserved");
		for(int i=0;i<result.size();i++){
			check(result.get(i)==i,"Completion order cannot reorder genome "+i);
		}
	}

	/**
	 * Keep the coordinator from returning while its sibling worker is held active.
	 * @param interrupt Select coordinator interruption instead of worker failure.
	 */
	private static void failureJoins(boolean interrupt)throws Exception{
		final CountDownLatch started=new CountDownLatch(2),release=new CountDownLatch(1),failed=new CountDownLatch(1),returned=new CountDownLatch(1);
		final AtomicBoolean siblingDone=new AtomicBoolean(),interruptPreserved=new AtomicBoolean();
		final AtomicReference<Throwable> outcome=new AtomicReference<Throwable>();
		final Thread coordinator=new Thread(()->{
			try{
				CallGenesLeanWorkers.run(2,2,i->{
					started.countDown();
					await(started);
					if(i==0 && !interrupt){
						failed.countDown();
						throw new IllegalArgumentException("injected worker failure");
					}
					try{
						await(release);
						return i;
					}finally{
						if(i==1){siblingDone.set(true);}
					}
				});
			}catch(Throwable failure){
				outcome.set(failure);
			}finally{
				interruptPreserved.set(Thread.currentThread().isInterrupted());
				returned.countDown();
			}
		},"lean-fixture-coordinator");
		coordinator.start();
		try{
			await(started);
			if(interrupt){
				coordinator.interrupt();
			}else{
				await(failed);
			}
			check(!returned.await(50,TimeUnit.MILLISECONDS),"Owner must wait for active worker after failure/interrupt");
		}finally{
			release.countDown();
			coordinator.join(10000);
		}
		check(!coordinator.isAlive() && siblingDone.get(),"Every active worker finished before coordinator exits");
		check(outcome.get() instanceof IllegalStateException,"Failure/interrupt cannot return successful cohort");
		check(interruptPreserved.get()==interrupt,"Coordinator interruption is preserved");
	}

	/** Bound fixture coordination and preserve an interrupted worker's status. */
	private static void await(CountDownLatch latch){
		try{
			if(!latch.await(10,TimeUnit.SECONDS)){
				throw new IllegalStateException("Fixture coordination timed out");
			}
		}catch(InterruptedException problem){
			Thread.currentThread().interrupt();
			throw new IllegalStateException(problem);
		}
	}

	/** Compare nullable GFF attribute values. */
	private static boolean equal(Object a,Object b){
		return a==null ? b==null : a.equals(b);
	}

	/** Count every explicit check, independently of the JVM assertion setting. */
	private static void check(boolean value,String message){
		checks++;
		if(!value){throw new AssertionError(message);}
	}

	private static int checks;
}
