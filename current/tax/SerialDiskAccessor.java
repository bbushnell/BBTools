package tax;

/**
 * Sequential adapter around DiskAccessionTable's memory-mapped lookup path.
 * Each batch is processed in input order on the caller thread, preserving
 * DiskAccessionTable's single-key lookup semantics exactly: nonpositive or
 * absent digitized keys return -1 and present keys return their stored taxID.
 * @author Brian Bushnell, Chloe
 */
public class SerialDiskAccessor extends DiskAccessor{

	/**
	 * Creates a serial accessor over an already opened table; close ownership stays
	 * with the caller because this wrapper never opens an independent handle.
	 * @param table Source accession table
	 */
	public SerialDiskAccessor(final DiskAccessionTable table){this.table=table;}

	/**
	 * Looks up each key through {@link DiskAccessionTable#get(long)} in input order.
	 * @param keys Digitized accession keys
	 * @return One taxID or -1 per key
	 */
	@Override
	public int[] getBatch(final long[] keys){
		final int[] results=new int[keys.length];
		for(int i=0; i<keys.length; i++){
			results[i]=table.get(keys[i]);
		}
		return results;
	}

	/** Returns the label used in accessor comparisons. */
	@Override public String name(){return "Serial (MappedByteBuffer)";}
	/** Does nothing because this wrapper does not own the wrapped table. */
	@Override public void close(){}

	/** Backing table that performs each mapped lookup. */
	private final DiskAccessionTable table;
}
