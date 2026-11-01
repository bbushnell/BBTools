package test.shared;

import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import shared.Shared;
import shared.SyncHeart;

/** Finite SyncHeart lock-release checks for invalid setters and valid controls.
 * Run one case per JVM: original unrepaired invalid -ea baselines deliberately
 * leaked a write lock, and all potentially blocked worker threads are daemon threads.
 * @author Jean
 */
public final class SyncHeartLockLeakTest{

	/** Runs exactly one isolated case: invalid-buffer-len, invalid-buffer-data, or valid-control. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected exactly one case=name argument");
		}
		final String name=args[0].substring(5);
		if("invalid-buffer-len".equals(name)){
			checkInvalidSetter(name, () -> SyncHeart.setBufferLen(0));
		}else if("invalid-buffer-data".equals(name)){
			checkInvalidSetter(name, () -> SyncHeart.setBufferData(0));
		}else if("valid-control".equals(name)){
			checkValidControl();
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		System.out.println("SyncHeartLockLeakTest case="+name+" assertions="+(Shared.EA() ? "enabled" : "disabled")+" PASS");
	}

	/**
	 * Exercises an invalid setter and verifies a later reader does not remain blocked.
	 * With assertions enabled, the setter must throw but still release SyncHeart's lock.
	 * With assertions disabled, the setter currently accepts the value; this case then
	 * only verifies that normal completion releases the lock.
	 */
	private static void checkInvalidSetter(String name, Action setter) throws Exception{
		final AtomicReference<Throwable> thrown=new AtomicReference<Throwable>();
		final Thread writer=new Thread(() -> {
			try{setter.run();}
			catch(Throwable t){thrown.set(t);}
		}, name+"-writer");
		writer.setDaemon(true);
		writer.start();
		writer.join(1000);
		require(!writer.isAlive(), name+" setter thread did not finish");
		if(Shared.EA()){
			require(thrown.get() instanceof AssertionError,
				name+" should throw AssertionError with -ea, observed "+thrown.get());
		}else{
			require(thrown.get()==null,
				name+" should not throw when assertion guards are disabled, observed "+thrown.get());
		}

		final AtomicBoolean readCompleted=new AtomicBoolean(false);
		final Thread reader=new Thread(() -> {
			SyncHeart.lowMemory();
			readCompleted.set(true);
		}, name+"-reader");
		reader.setDaemon(true);
		reader.start();
		reader.join(1000);
		require(readCompleted.get(), name+" leaked the write lock; later readers block");
	}

	/** Valid setters must release locks on the ordinary no-throw path. */
	private static void checkValidControl() throws Exception{
		final AtomicReference<Throwable> thrown=new AtomicReference<Throwable>();
		final Thread writer=new Thread(() -> {
			try{
				SyncHeart.setLowMemory(true);
				SyncHeart.setLowMemory(false);
			}catch(Throwable t){thrown.set(t);}
		}, "valid-control-writer");
		writer.setDaemon(true);
		writer.start();
		writer.join(1000);
		require(!writer.isAlive(), "valid setter thread did not finish");
		require(thrown.get()==null, "valid setter path threw "+thrown.get());

		final AtomicBoolean readCompleted=new AtomicBoolean(false);
		final Thread reader=new Thread(() -> {
			readCompleted.set(!SyncHeart.lowMemory());
		}, "valid-control-reader");
		reader.setDaemon(true);
		reader.start();
		reader.join(1000);
		require(readCompleted.get(), "valid setter path blocked a later reader");
	}

	/** Checks expectations in both assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** One setter call that may assert under -ea. */
	private interface Action{void run();}
}
