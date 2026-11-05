package kmer;

/** Shared GC-by-count accumulation operation for short and long k-mer tables. */
public interface KmerTableInterface{

	/**
	 * Adds G/C-base totals to bins selected by stored k-mer count.
	 * Accumulates into existing bins without clearing them; callers coordinate
	 * access to the table and the output array. Entry filtering is implementation-specific.
	 * @param gcCounts Caller-owned accumulation array, with indices through max
	 * @param max Nonnegative highest count bin; counts above max use that bin
	 */
	public void countGC(long[] gcCounts, int max);

//	public <T extends KmerTableInterface> T[] tables();

}
