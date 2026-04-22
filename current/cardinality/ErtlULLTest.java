package cardinality;

import shared.Tools;

/**
 * Exercises ErtlULL pack/unpack and register filling at
 * cardinalities 100, 10k, and 100k.
 *
 * @author Brian Bushnell
 */
public class ErtlULLTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final ErtlULL ull=new ErtlULL(2048, 31, -1, 0);
		// Check pack/unpack of register 0
		System.err.println("unpack(0)="+Long.toHexString(ErtlULL.unpack((byte)0)));
		System.err.println("pack(0)="+(ErtlULL.pack(0L)&0xFF));
		// Add some elements
		for(int i=1; i<=100; i++){ull.add(i);}
		// Check register fill
		int filled=0, maxR=0;
		for(int j=0; j<2048; j++){
			final int r=ull.registers[j]&0xFF;
			if(r>0){filled++;}
			if(r>maxR){maxR=r;}
		}
		System.err.println("After 100: filled="+filled+" maxR="+maxR+" card="+ull.cardinality());
		for(int i=101; i<=10000; i++){ull.add(i);}
		filled=0; maxR=0;
		for(int j=0; j<2048; j++){
			final int r=ull.registers[j]&0xFF;
			if(r>0){filled++;}
			if(r>maxR){maxR=r;}
		}
		System.err.println("After 10k: filled="+filled+" maxR="+maxR+" card="+ull.cardinality());
		for(int i=10001; i<=100000; i++){ull.add(i);}
		filled=0; maxR=0;
		for(int j=0; j<2048; j++){
			final int r=ull.registers[j]&0xFF;
			if(r>0){filled++;}
			if(r>maxR){maxR=r;}
		}
		System.err.println("After 100k: filled="+filled+" maxR="+maxR+" card="+ull.cardinality());
	}

}
