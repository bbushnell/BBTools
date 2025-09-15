package driver;

import java.util.ArrayList;

import fileIO.ReadWrite;

public class ClearRam {
	
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
	
	private static long current=0;
	
}
