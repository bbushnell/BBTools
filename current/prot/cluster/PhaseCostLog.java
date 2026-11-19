package prot.cluster;

import java.lang.management.ManagementFactory;
import java.lang.reflect.Method;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Records nested coordinator intervals without changing clustering decisions.
 * Every interval is inclusive: children are already included in their parent.
 * CPU is the JVM process counter, including workers and GC during the interval;
 * it is neither exclusive algorithm CPU nor scheduler whole-job accounting.
 * A coordinator's context is deliberately not inherited by worker threads.
 * @author Brian Bushnell, Fischl
 */
public final class PhaseCostLog implements AutoCloseable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens one new report on the calling coordinator thread. */
	public PhaseCostLog(final String path){this(path, JVM_CLOCK);}

	/** An injectable clock permits exact interval and unavailable-counter tests. */
	public PhaseCostLog(final String path, final Clock clock_){
		if(path==null || path.isEmpty() || clock_==null || ACTIVE.get()!=null){
			throw new IllegalArgumentException("Phase costs need a new output path, a clock and no active coordinator log");
		}
		clock=clock_;
		thread=Thread.currentThread();
		writer=new ByteStreamWriter(path, false, false, true);
		writer.start();
		writer.println(HEADER);
		ACTIVE.set(this);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Scopes            ----------------*/
	/*--------------------------------------------------------------*/

	/** Begins an inclusive interval, or returns a no-op when no report is active.
	 * Call only at phase boundaries, never once per sequence or alignment.
	 */
	public static Span begin(final String name){
		final PhaseCostLog log=ACTIVE.get();
		if(log==null){return NOOP;}
		log.requireOwner();
		if(name==null || name.isEmpty() || name.indexOf('\t')>=0 || name.indexOf('\n')>=0 || name.indexOf('\r')>=0){
			throw new IllegalArgumentException("Phase names must be nonempty single TSV fields");
		}
		final Span span=new Span(log, log.top, log.nextId, name);
		log.nextId=Math.incrementExact(log.nextId);
		log.top=span;
		return span;
	}

	/** Drains the report after every interval has closed; repeated close is harmless. */
	@Override
	public void close(){
		if(closed){return;}
		requireOwner();
		if(top!=null){throw new IllegalStateException("Close the active phase before closing its cost report: "+top.name);}
		closed=true;
		ACTIVE.remove();
		if(writer.poisonAndWait()){throw new IllegalStateException("Failed writing phase costs");}
	}

	/** Prevents unrelated coordinator/worker timelines from being mixed in one tree. */
	private void requireOwner(){
		if(closed || Thread.currentThread()!=thread || ACTIVE.get()!=this){
			throw new IllegalStateException("A phase-cost report belongs to its opening coordinator until close");
		}
	}

	/** One try-with-resources interval. Closing it records elapsed time even when
	 * the operation throws; a timing row is not a success or completion marker.
	 */
	public static final class Span implements AutoCloseable{
		private Span(final PhaseCostLog owner_, final Span parent_, final int id_, final String name_){
			assert(parent_==null || parent_.id<id_) :
				"A phase parent was opened earlier; increasing IDs preserve an acyclic inclusive timing tree";
			owner=owner_;
			parent=parent_;
			id=id_;
			name=name_;
			wallStart=owner==null ? 0 : owner.clock.wallNanos();
			cpuStart=owner==null ? -1 : owner.clock.processCpuNanos();
		}

		/** Publishes an inclusive row after checking nesting and monotonic counters. */
		@Override
		public void close(){
			if(owner==null || closed){return;}
			owner.requireOwner();
			if(owner.top!=this){throw new IllegalStateException("Phase intervals must close inside-out: "+name);}
			final long cpuEnd=owner.clock.processCpuNanos();
			final long wall=owner.clock.wallNanos()-wallStart;
			final boolean cpuAvailable=cpuStart>=0 && cpuEnd>=0;
			if(wall<0 || cpuAvailable && cpuEnd<cpuStart){
				throw new IllegalStateException("Elapsed phase counters must be monotonic: "+name);
			}
			final ByteBuilder row=new ByteBuilder();
			row.append(id).tab().append(parent==null ? -1 : parent.id).tab().append(name).tab().append(wall).tab();
			if(cpuAvailable){row.append(cpuEnd-cpuStart);}
			else{row.append("NA");}
			row.tab().append(cpuAvailable ? "measured" : "unavailable");
			owner.top=parent;
			closed=true;
			owner.writer.println(row);
		}

		private final PhaseCostLog owner;
		private final Span parent;
		private final int id;
		private final String name;
		private final long wallStart, cpuStart;
		private boolean closed;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Clocks            ----------------*/
	/*--------------------------------------------------------------*/

	/** Nanosecond counters. Negative process CPU means unavailable, never zero work. */
	public interface Clock{
		long wallNanos();
		long processCpuNanos();
	}

	/** Reflective access keeps non-HotSpot JVMs and Java 8 usable without a
	 * mandatory link to a vendor management extension. Unsupported CPU is explicit.
	 */
	private static final class JvmClock implements Clock{
		JvmClock(){
			Object foundBean=null;
			Method foundMethod=null;
			try{
				foundBean=ManagementFactory.getOperatingSystemMXBean();
				final Class<?> type=Class.forName("com.sun.management.OperatingSystemMXBean");
				if(type.isInstance(foundBean)){foundMethod=type.getMethod("getProcessCpuTime");}
			}catch(ReflectiveOperationException | SecurityException unavailable){
				foundMethod=null;
			}
			bean=foundBean;
			method=foundMethod;
		}

		@Override
		public long wallNanos(){return System.nanoTime();}

		@Override
		public long processCpuNanos(){
			if(method==null){return -1;}
			try{
				final Object value=method.invoke(bean);
				if(!(value instanceof Long)){throw new IllegalStateException("Process CPU must be a nanosecond long counter");}
				return (Long)value;
			}catch(ReflectiveOperationException | SecurityException unavailable){
				return -1;
			}
		}

		private final Object bean;
		private final Method method;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Clock clock;
	private final Thread thread;
	private final ByteStreamWriter writer;
	private Span top;
	private int nextId;
	private boolean closed;
	private static final Clock JVM_CLOCK=new JvmClock();
	private static final ThreadLocal<PhaseCostLog> ACTIVE=new ThreadLocal<PhaseCostLog>();
	private static final Span NOOP=new Span(null, null, -1, null);
	public static final String HEADER="phase_id\tparent_phase_id\tphase\tinclusive_wall_ns\tinclusive_jvm_cpu_ns\tcpu_status";
}
