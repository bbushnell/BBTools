package driver;

import java.util.ArrayList;

import fileIO.ReadWrite;

/**
 * Memory testing utility for determining available RAM and memory allocation patterns.
 * Allocates arrays of long values to test memory limits and create temporary files
 * with serialized data for testing read/write operations.
 * @author Brian Bushnell
 */
public class ClearRam {
	
	/**
	 * Program entry point. Runs memory allocation test twice with garbage collection
	 * between attempts to determine maximum available memory.
	 * @param args Command-line arguments (unused)
	 */
	public static void main(String[] args){
		
		for(int i=0; i<2; i++){
			
			try {
				System.gc();
				attempt();
			} catch(final java.lang.OutOfMemoryError e) {
//				e.printStackTrace();
				System.err.println("Out of memory at "+((current*8)/(1<<20))+" MB");
			}
		}
	}
	
	/**
	 * Attempts to allocate maximum memory by creating long arrays.
	 * Creates arrays of 2^20 long values until OutOfMemoryError occurs.
	 * Tracks current allocation count for memory usage reporting.
	 */
	public static void attempt(){
		ArrayList<long[]> list=new ArrayList<long[]>(8000);
		current=0;
		
		while(true){
			long[] array=null;

			array=new long[1<<20];
			list.add(array);

//			for(int i=0; i<array.length; i++){
//				array[i]=current;
//				current++;
//			}
			current+=array.length;
		}
	}
	
	/**
	 * Creates and writes serialized long array data files for testing.
	 * Generates sequential long values based on current system time and writes
	 * them to disk as serialized objects. Also validates existing files.
	 * @param megs Size in megabytes for the data array to create
	 */
	public static void writeJunk(int megs){
		try {
			long[] old=(long[]) ReadWrite.readObject("JUNK"+megs+".long", false);
			for(int i=1; i<old.length; i++){
				assert(old[i]==old[i-1]+1);
			}
		} catch (Exception e) {
			
		}
		
		
		
		long[] array=new long[megs*(1<<17)];
		long current=System.nanoTime();
		for(int i=0; i<array.length; i++){
			array[i]=current+i;
		}
		ReadWrite.write(array, "JUNK"+megs+".long", false);
		System.err.println("Wrote "+((8*array.length)/(1024000))+" MB junk");
	}
	
	/** Counter tracking total number of long values allocated during memory test */
	private static long current=0;
	
}
