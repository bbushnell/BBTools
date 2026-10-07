package stream;

import parse.Parse;
import shared.Shared;

/**
 * Selects ST or ZT leaves for a resource-bounded multi-file writer.
 * ST means one background writer thread per physical output; ZT writes on the
 * submitting thread. Neither choice permits a hidden multithreaded compressor.
 * The leaf factory must enforce that contract; this class only selects a mode.
 *
 * Configure static defaults during argument parsing, then call {@link #snapshot()}.
 * Each container retains its immutable policy and selects once before opening files.
 * Supply the largest simultaneous physical-output count it may open, including
 * mates and separate quality files. For unknown fan-out, use the open-handle cap;
 * do not use the number of destinations discovered so far.
 *
 * AUTO reserves one configured thread for the caller and at most one eighth of
 * available heap for the estimated ST allowance. The allowance and thread cap are
 * provisional tuning parameters, not measured memory requirements or throughput
 * guarantees. Native stacks, descriptors and compressor buffers require separate
 * accounting by the container and factory. Rotation is independent of this policy.
 *
 * @author Shinobu
 * @date October 1, 2026
 */
public final class MultiWriterPolicy{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an immutable policy independently of the process-wide defaults.
	 * @param mode_ AUTO, ST or ZT
	 * @param maxST_ Maximum simultaneous ST outputs permitted by AUTO
	 * @param bytesPerST_ Provisional memory allowance per ST output, in bytes
	 */
	public MultiWriterPolicy(final int mode_, final int maxST_, final long bytesPerST_){
		if(mode_!=AUTO && mode_!=ST && mode_!=ZT){
			throw new IllegalArgumentException("Multi-writer mode must be AUTO, ST or ZT: "+mode_);
		}
		if(maxST_<0 || bytesPerST_<1){
			throw new IllegalArgumentException("AUTO requires maxST>=0 and bytesPerST>0: "+maxST_+", "+bytesPerST_);
		}
		mode=mode_;
		maxST=maxST_;
		bytesPerST=bytesPerST_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Selects using the configured thread count and currently available JVM heap.
	 * This is a snapshot of heap headroom, not host RAM or idle-core detection.
	 * @param outputs Upper bound on simultaneous physical outputs
	 * @return ST or ZT; never AUTO
	 */
	public int select(final int outputs){
		return select(outputs, Shared.threads(), Shared.memAvailable());
	}

	/**
	 * Selects using explicit resource estimates, allowing deterministic policy tests.
	 * Forced modes override the AUTO thresholds but do not remove handle or buffer
	 * limits elsewhere. Zero outputs selects ZT in AUTO. Negative heap headroom is
	 * treated as exhausted; configured thread counts must be positive.
	 * @param outputs Upper bound on simultaneous physical outputs, including mates
	 * @param threads Available configured threads, including the submitting thread
	 * @param availableBytes Available JVM heap headroom in bytes
	 * @return ST or ZT; never AUTO
	 */
	public int select(final int outputs, final int threads, final long availableBytes){
		if(outputs<0 || threads<1){
			throw new IllegalArgumentException("Selection requires outputs>=0 and threads>=1: "+outputs+", "+threads);
		}
		if(mode!=AUTO){return mode;}
		// Divide before comparing so even a very large output count cannot overflow.
		final long memorySlots=Math.max(0L, availableBytes)/8/bytesPerST;
		final int selected=(outputs>0 && outputs<=maxST && outputs<threads && outputs<=memorySlots ? ST : ZT);
		assert(selected!=ST || outputs<threads) : "AUTO must leave a configured thread for the submitting caller";
		return selected;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies startup defaults for a new container. Parse arguments before calling;
	 * later default changes do not alter an existing snapshot.
	 * @return Immutable policy using the current defaults
	 */
	public static MultiWriterPolicy snapshot(){
		return new MultiWriterPolicy(DEFAULT_MODE, AUTO_MAX_ST, AUTO_BYTES_PER_ST);
	}

	/**
	 * Parses multi-writer startup options. Call through Parser for ordinary CLI use.
	 * Names are case-insensitive. Recognized options require a value; an unknown
	 * name returns false without changing any setting.
	 * @param name Parameter name
	 * @param value Parameter value
	 * @return Whether this class consumed the option
	 */
	public static boolean parse(final String name, final String value){
		if("multiwritermode".equalsIgnoreCase(name)){
			DEFAULT_MODE=parseMode(value);
		}else if("multiwriterstmax".equalsIgnoreCase(name)){
			final int x=Integer.parseInt(value);
			if(x<0){throw new IllegalArgumentException("multiwriterstmax must be nonnegative: "+value);}
			AUTO_MAX_ST=x;
		}else if("multiwriterstmemory".equalsIgnoreCase(name)){
			if(value==null){throw new IllegalArgumentException("multiwriterstmemory requires a byte count");}
			final long x=Parse.parseKMG(value);
			if(x<1){throw new IllegalArgumentException("multiwriterstmemory must be positive: "+value);}
			AUTO_BYTES_PER_ST=x;
		}else{return false;}
		return true;
	}

	/**
	 * Converts an explicit mode name to its constant, rejecting unsupported modes.
	 * @param value auto, st or zt, ignoring case
	 * @return AUTO, ST or ZT
	 */
	public static int parseMode(final String value){
		if("auto".equalsIgnoreCase(value)){return AUTO;}
		if("st".equalsIgnoreCase(value)){return ST;}
		if("zt".equalsIgnoreCase(value)){return ZT;}
		throw new IllegalArgumentException("multiwritermode must be auto, st or zt: "+value);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Requested mode, including AUTO. */
	public final int mode;
	/** Maximum simultaneous ST outputs allowed by AUTO. */
	public final int maxST;
	/** Estimated bytes per ST output used only by AUTO. */
	public final long bytesPerST;

	/*--------------------------------------------------------------*/
	/*----------------      Static Configuration    ----------------*/
	/*--------------------------------------------------------------*/

	/** Startup mode, set with multiwritermode=auto|st|zt. */
	public static int DEFAULT_MODE=-1;
	/** Provisional AUTO cap, set with multiwriterstmax; zero disables AUTO ST. */
	public static int AUTO_MAX_ST=32;
	/** Provisional per-ST allowance; multiwriterstmemory accepts K/M/G suffixes. */
	public static long AUTO_BYTES_PER_ST=8_000_000L;

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Select according to output count, configured threads and heap headroom. */
	public static final int AUTO=-1;
	/** No background writer threads. */
	public static final int ZT=0;
	/** One background writer thread per physical output. */
	public static final int ST=1;
}
