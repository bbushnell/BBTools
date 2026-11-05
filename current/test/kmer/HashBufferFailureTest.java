package test.kmer;

import java.lang.reflect.Method;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReentrantLock;

import kmer.AbstractKmerTable;
import kmer.HashArray;
import kmer.HashArray1D;
import kmer.HashArray2D;
import kmer.HashBuffer;

/** Tests release of the backend lock when a buffered operation is rejected.
 * Uses an existing unsupported-operation exception, never an allocation failure.
 * Reflection observes the package-private backend lock without changing baseline APIs.
 * @author Shinobu */
public final class HashBufferFailureTest {

	/** Runs a normal count control or the rejected-operation lock-release case. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected case=controls or case=lock-release.");}
		if(args[0].equals("case=controls")){
			final HashArray1D table=new HashArray1D(new int[] {101}, -1L);
			final HashBuffer buffer=new HashBuffer(new AbstractKmerTable[] {table}, 8, 3, false, false);
			check(buffer.incrementAndReturnNumCreated(1, 1)==0, "First increment must remain buffered.");
			check(buffer.incrementAndReturnNumCreated(1, 1)==0, "Second increment must remain buffered.");
			check(buffer.flush()==1 && table.getValue(1)==2, "Flush must create one count-2 key.");
			check(buffer.flush()==0 && table.size()==1, "Empty flush must preserve counts and size.");
		}else if(args[0].equals("case=lock-release")){
			final HashArray2D table=new HashArray2D(new int[] {101}, -1L);
			final Method accessor=HashArray.class.getDeclaredMethod("getLock");
			accessor.setAccessible(true);
			final ReentrantLock lock=(ReentrantLock)accessor.invoke(table);
			final HashBuffer buffer=new HashBuffer(new AbstractKmerTable[] {table}, 1, 3, false, false);
			try{
				boolean rejected=false;
				try{buffer.incrementAndReturnNumCreated(1, 1);}
				catch(RuntimeException expected){
					rejected="Unsupported.".equals(expected.getMessage());
					if(!rejected){throw expected;}
				}
				check(rejected, "Multi-value backend must reject the count-only operation.");
				final AtomicBoolean available=new AtomicBoolean();
				final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
				final Thread probe=new Thread(new Runnable(){
					@Override
					public void run(){
						try{
							if(lock.tryLock()){
								try{available.set(true);}finally{lock.unlock();}
							}
						}catch(Throwable t){failure.set(t);}
					}
				}, "HashBufferFailureTest-lock-probe");
				probe.setDaemon(true);
				probe.start();
				probe.join(5000);
				check(!probe.isAlive() && failure.get()==null, "Nonblocking lock probe did not complete normally.");
				check(available.get(), "Buffered failure left the backend locked by the calling thread.");
			}finally{
				// Release only this fixture's own leaked hold in the failing baseline.
				while(lock.isHeldByCurrentThread()){lock.unlock();}
			}
		}else{throw new IllegalArgumentException(args[0]);}
		System.out.println("PASS "+args[0]);
	}

	/** Keeps the oracle enabled with Java assertions disabled. */
	private static void check(final boolean ok, final String message){
		if(!ok){throw new AssertionError(message);}
	}
}
