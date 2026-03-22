package cardinality;

import shared.Tools;

/**
 * Debug: compare register states between UDLL6 and ErtlULL at low cardinality.
 */
public class UDLL6Test {
	public static void main(String[] args){
		final int buckets=2048;
		final int k=31;
		final long seed=12345L;

		UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, k, seed, 0);
		ErtlULL ertl=new ErtlULL(buckets, k, seed, 0);

		// Add 10 elements and compare register states
		for(long i=1; i<=10; i++){
			udll.add(i);
			ertl.add(i);
		}

		// Compare registers for first 20 non-empty buckets
		System.out.println("Bucket\tUDLL6_reg\tErtl_reg\tUDLL6_nlzPart\tErtl_nlzPart\tUDLL6_sub\tErtl_sub");
		int shown=0;
		byte[] ertlRegs=ertl.getRegisters();
		for(int i=0; i<buckets && shown<30; i++){
			int ur=udll.getRegister(i)&0x3F;
			int er=ertlRegs[i]&0xFF;
			if(ur!=0 || er!=0){
				int unlzPart=ur>>>2;
				int usub=ur&3;
				int enlzPart=er>>>2;
				int esub=er&3;
				System.out.printf("%d\t%d\t%d\t%d\t%d\t%d\t%d%n",
					i, ur, er, unlzPart, enlzPart, usub, esub);
				shown++;
			}
		}

		// Also trace a single add to see what happens
		System.out.println("\n--- Tracing element 11 ---");
		long key11=Tools.hash64shift(11^seed);
		int nlz11=Long.numberOfLeadingZeros(key11);
		int bucket11=(int)(key11&(buckets-1));
		System.out.println("key="+Long.toHexString(key11)+" nlz="+nlz11+" bucket="+bucket11);
		System.out.println("UDLL6 bitPos=relNlz+2="+(nlz11-udll.getMinZeros())+"+2="+(nlz11-udll.getMinZeros()+2));

		// What would Ertl do?
		int p=11; // log2(2048)
		int q=64-p; // = 53
		int ertlIdx=(int)(key11>>>q);
		int ertlNlz=Long.numberOfLeadingZeros(~(~key11<<-q));
		System.out.println("Ertl: idx="+ertlIdx+" nlz="+ertlNlz+" bitPos="+(ertlNlz+p-1));
	}
}
