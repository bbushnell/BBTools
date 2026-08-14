package tax;

/**
 * Serial disk accessor — delegates to DiskAccessionTable.get() sequentially.
 * Baseline implementation using MappedByteBuffer (page faults block per read).
 * @author Brian Bushnell, Chloe
 */
public class SerialDiskAccessor extends DiskAccessor {

	public SerialDiskAccessor(DiskAccessionTable table){this.table=table;}

	@Override
	public int[] getBatch(long[] keys){
		int[] results=new int[keys.length];
		for(int i=0; i<keys.length; i++){
			results[i]=table.get(keys[i]);
		}
		return results;
	}

	@Override public String name(){return "Serial (MappedByteBuffer)";}
	@Override public void close(){}

	private final DiskAccessionTable table;
}
