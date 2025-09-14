package bbmin;

import java.util.Arrays;

import structures.LongHashSet;
import structures.LongList;

/**
 * Generates an array of minimal hash codes (as positive 64-bit longs) for an input sequence.<br>
 * The resulting array is guaranteed to contain the minimal hash code<br>
 * for every window, with no duplicates.
 * On average this is expected to yield 2*(L-K)/W hash codes for sequence length L and window size W.
 * 
 * @author Brian Bushnell
 * @date October 8, 2021
 *
 */
public class Minimizer {
	
	public static void main(String[] args){
		int k=4, w=7;
		String seq="ACGTCTGAGCCTTGACACATGACT";
		try {
			k=Integer.parseInt(args[0]);
			w=Integer.parseInt(args[1]);
			seq=args[2];
		} catch (NumberFormatException e) {
			//e.printStackTrace();
			System.err.println("Usage: bbmin.Minimizer kmerlen window seq\n"
					+ "E.G.\n"
					+ "bbmin.Minimizer 4 7 ACGTCTGAGCCTTGACACATGACT");
			System.exit(1);
		}
		Minimizer minnow=new Minimizer(k, w);
		long[] array=minnow.minimize(seq.getBytes());
		System.err.println(Arrays.toString(array));
	}

	public Minimizer(int k_, int window_){this(k_, window_, 2);}
	public Minimizer(int k_, int window_, int bitsPerSymbol_){
		k=k_;
		window=window_;
		bitsPerSymbol=bitsPerSymbol_;
		shift=bitsPerSymbol*k;
		shift2=shift-bitsPerSymbol;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
	}
	
	public long[] minimize(String str){
		return minimize(str.getBytes());
	}

	public long[] minimize(byte[] bases){
		return minimize(bases, new LongList(16), new LongHashSet(16));
	}

	/** This method is typically faster since you don't need to construct a new set each time. */
	public long[] minimize(byte[] bases, LongList list, LongHashSet set){
		list.clear();
		//If the set is way too big, resize it
		if(set.capacity()*(long)window>100L+16L*bases.length){
			set.resizeDestructive(16);
		}else{
			set.clear();
		}
		
		long kmersProcessed=0;
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		long bestCode=Long.MAX_VALUE;
		int bestPosition=-1;
		long bestKmer=-1;
		long bestRkmer=-1;
		int currentWindow=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=baseToNumber[b];
			long x2=baseToComplementNumber[b];
			
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				rkmer=0;
			}else{
				len++;
			}
			
			if(len>=k){
				kmersProcessed++;
				currentWindow++;

				final long hashcode=hash(kmer, rkmer);
				System.err.println("i="+i+", code="+hashcode);

				//Track the best code in the window and its state
				if(hashcode>=minCode && hashcode<=bestCode){
					bestCode=hashcode;
					bestPosition=i;
					bestKmer=kmer;
					bestRkmer=rkmer;
				}
				
				//Once the window size is met, store the best code,
				//and backtrack to its position to start the next window
				if(currentWindow>=window && bestPosition>=0){
					if(!set.contains(bestCode)){
						set.add(bestCode);
						list.add(bestCode);
					}
					i=bestPosition;
					kmer=bestKmer;
					rkmer=bestRkmer;
					len=k;
					
					bestCode=Long.MAX_VALUE;
					bestPosition=-1;
					currentWindow=0;
				}
			}
		}
		list.sort();//optional
		return list.toArray();
	}

	public static long canon(long kmer, long rkmer){return max(kmer, rkmer);}
	public static long hash(long kmer, long rkmer){return hash(canon(kmer, rkmer));}
	public static long hash(long key) {
		key = (~key) + (key << 21); // key = (key << 21) - key - 1;
		key = key ^ (key >>> 24);
		key = (key + (key << 3)) + (key << 8); // key * 265
		key = key ^ (key >>> 14);
		key = (key + (key << 2)) + (key << 4); // key * 21
		key = key ^ (key >>> 28);
		key = key + (key << 31);
		return key;
	}
	private static final long max(long x, long y){return x>y ? x : y;}

	public final int k;
	public final int window;
	public final int bitsPerSymbol; //2 for nucleotides, 5 for amino acids.
	private final int shift;
	private final int shift2;
	private final long mask;
	private final long minCode=0;
	
	static final byte[] baseToNumber = new byte[128];
	static final byte[] baseToComplementNumber = new byte[128];
	
	static {
		Arrays.fill(baseToNumber, (byte)-1);
		Arrays.fill(baseToComplementNumber, (byte)-1);
		baseToNumber['A']=baseToNumber['a']=baseToComplementNumber['T']=baseToComplementNumber['t']=0;
		baseToNumber['C']=baseToNumber['c']=baseToComplementNumber['G']=baseToComplementNumber['c']=1;
		baseToNumber['G']=baseToNumber['g']=baseToComplementNumber['C']=baseToComplementNumber['g']=2;
		baseToNumber['T']=baseToNumber['t']=baseToComplementNumber['A']=baseToComplementNumber['a']=3;
	}
}
