package align2;

import structures.IntList;

/**
 * Legacy lossy collapse of tandem repeats in ASCII/genomic byte sequences.
 * This is not a reversible codec. Period selection and retained-copy counts
 * differ between the fixed, multiperiod, and ultra methods. Input is not changed.
 * String conversion uses the platform charset; arbitrary binary/Unicode input
 * is not the intended representation.
 *
 * @author Brian Bushnell, Collei
 * @date 2013
 */
public class CompressString{

	/** Prints the input and its period-1..3 ultra collapse. */
	public static void main(String[] args){
		if(args.length<1){throw new IllegalArgumentException("Expected a sequence to compress");}
		final String s=compressRepeatsUltra(args[0].getBytes(), 1, 3, null);
		System.out.println(args[0]+"\n"+s);
	}

	/** Applies fixed-period collapse sequentially at periods 1, 2, and 3.
	 * Later passes operate on the preceding pass's shortened output. */
	public static String compress(String s){
		final String s1=compressRepeats(s.getBytes(), 1);
		final String s2=compressRepeats(s1.getBytes(), 2);
		return compressRepeats(s2.getBytes(), 3);
	}

	/** Collapses one period. countRepeats excludes the initial occurrence.
	 * Two additional copies trigger a one-copy skip. Larger runs emit log2-sized
	 * output using the historical rule, then revisit the last original occurrence. */
	public static String compressRepeats(byte[] array, int period){
		if(period<1){throw new IllegalArgumentException("Repeat period must be positive: "+period);}
		final StringBuilder sb=new StringBuilder(array.length);
		for(int base=0; base<array.length; base++){
			final int repeats=countRepeats(array, base, period);
			if(repeats<2){
				sb.append((char)array[base]);
			}else if(repeats==2){
				base+=period-1;
			}else{
				final int log=(32-Integer.numberOfLeadingZeros(repeats+1))-1;
				assert(log>0 && log<=31) : "Positive repeat count must yield a positive int log2: "+repeats;
				for(int i=1; i<log; i++){
					for(int j=0; j<period; j++){sb.append((char)array[base+j]);}
				}
				base=base+(period*repeats)-1;
			}
		}
		return sb.toString();
	}

	/** Chooses the smallest period with at least two additional copies.
	 * Emits floor(log2(repeats+1)) copies and revisits the last original copy.
	 * If supplied, list is appended with one source index per emitted character;
	 * repeated emitted copies reuse indices, so the list need not be monotonic. */
	public static String compressRepeatsMultiperiod(byte[] array, int minPeriod, int maxPeriod, IntList list){
		validateRange(minPeriod, maxPeriod);
		final StringBuilder sb=new StringBuilder(array.length);
		for(int base=0; base<array.length; base++){
			int period=0, repeats=0;
			// Two additional complete copies require three periods in the remaining input.
			final int limit=Math.min(maxPeriod, (array.length-base)/3);
			for(int x=minPeriod; x<=limit; x++){
				final int temp=countRepeats(array, base, x);
				if(temp>1){repeats=temp; period=x; break;}
			}
			if(repeats==0){
				sb.append((char)array[base]);
				if(list!=null){list.add(base);}
			}else{
				final int log=(32-Integer.numberOfLeadingZeros(repeats+1))-1;
				assert(log>0 && log<=31) : "Selected repeat count must yield a positive int log2: "+repeats;
				for(int i=0; i<log; i++){
					for(int j=0; j<period; j++){
						sb.append((char)array[base+j]);
						if(list!=null){list.add(base+j);}
					}
				}
				base=base+(period*repeats)-1;
			}
		}
		return sb.toString();
	}

	/** Emits one copy of the first qualifying pattern, then revisits its last copy.
	 * Thus a plain isolated run generally retains two copies, not just one.
	 * Period selection and optional appended position indices follow the multiperiod method. */
	public static String compressRepeatsUltra(byte[] array, int minPeriod, int maxPeriod, IntList list){
		validateRange(minPeriod, maxPeriod);
		final StringBuilder sb=new StringBuilder(array.length);
		for(int base=0; base<array.length; base++){
			int period=0, repeats=0;
			final int limit=Math.min(maxPeriod, (array.length-base)/3);
			for(int x=minPeriod; x<=limit; x++){
				final int temp=countRepeats(array, base, x);
				if(temp>1){repeats=temp; period=x; break;}
			}
			if(repeats==0){
				sb.append((char)array[base]);
				if(list!=null){list.add(base);}
			}else{
				for(int j=0; j<period; j++){
					sb.append((char)array[base+j]);
					if(list!=null){list.add(base+j);}
				}
				base=base+(period*repeats)-1;
			}
		}
		return sb.toString();
	}

	/** Counts complete consecutive copies after the initial pattern at base.
	 * A partial matching suffix is not counted. base may equal array.length;
	 * periods longer than the remaining input return zero without index arithmetic. */
	public static int countRepeats(byte[] array, int base, int period){
		if(period<1){throw new IllegalArgumentException("Repeat period must be positive: "+period);}
		if(base<0 || base>array.length){throw new IllegalArgumentException("Repeat base outside input: "+base+" / "+array.length);}
		if(period>array.length-base){return 0;}
		final int max=array.length-period+1;
		int matches=0;
		boolean fail=false;
		for(int loc=base+period; loc<max && !fail; loc+=period){
			for(int i=0; i<period && !fail; i++){
				if(array[base+i]==array[loc+i]){matches++;}
				else{fail=true;}
			}
		}
		return matches/period;
	}

	private static void validateRange(int minPeriod, int maxPeriod){
		if(minPeriod<1 || maxPeriod<minPeriod){
			throw new IllegalArgumentException("Repeat period range must be positive and ordered: "+minPeriod+".."+maxPeriod);
		}
	}
}
