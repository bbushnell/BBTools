package test.ukmer;

import java.lang.reflect.Field;
import java.util.Arrays;

import ukmer.Kmer;

/** Checks reverse-complement cache state against independently rebuilt poly-base keys.
 * Reflection selects the existing package-private orientation mode without changing
 * production visibility merely to establish the pre-repair baseline.
 * @author Shinobu */
public final class KmerStateContractTest {

	/** Runs controls or the forward-only, already-cached reverse-complement case. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected case=controls or case=forward-warm.");}
		final boolean failing=args[0].equals("case=forward-warm");
		if(!failing && !args[0].equals("case=controls")){throw new IllegalArgumentException(args[0]);}
		final Field orientation=Kmer.class.getDeclaredField("rcomp");
		orientation.setAccessible(true);
		final boolean oldOrientation=orientation.getBoolean(null), oldPacked=Kmer.PACKED;
		int failures=0, cases=0;
		try{
			for(boolean packed : new boolean[] {false, true}){
				Kmer.PACKED=packed;
				for(int requested : new int[] {1, 31, 32, 33, 63, 64, 95}){
					if(failing){
						orientation.setBoolean(null, false);
						if(!checkCase(requested, true)){failures++;}
						cases++;
					}else{
						for(boolean canonical : new boolean[] {false, true}){
							orientation.setBoolean(null, canonical);
							if(!checkCase(requested, canonical)){failures++;}
							cases++;
						}
					}
				}
			}
		}finally{
			orientation.setBoolean(null, oldOrientation);
			Kmer.PACKED=oldPacked;
		}
		if(failures>0){throw new AssertionError("Reverse-complement state failed "+failures+" of "+cases+" cases.");}
		System.out.println("PASS "+args[0]+" cases="+cases);
	}

	/** Poly-A becomes poly-T; a rebuilt key is independent of the object's lazy caches. */
	private static boolean checkCase(final int requested, final boolean warm){
		final Kmer actual=make(requested, 0);
		assert(actual.len()==actual.kbig) : "Every case needs a complete effective window before cache evaluation.";
		if(warm){actual.key(); actual.xor(); actual.xor2(); actual.corePalindrome();}
		actual.rcomp();
		final Kmer expected=make(requested, 3);
		final boolean arrays=Arrays.equals(actual.array1(), expected.array1()) && Arrays.equals(actual.array2(), expected.array2());
		final boolean keys=Arrays.equals(actual.key(), expected.key());
		final boolean hashes=actual.xor()==expected.xor() && actual.xor2()==expected.xor2();
		final boolean palindrome=actual.corePalindrome()==expected.corePalindrome();
		final boolean valid=actual.len()==expected.len() && actual.verify(true);
		if(!(arrays && keys && hashes && palindrome && valid)){
			System.err.println("FAIL requested="+requested+", packed="+Kmer.PACKED+", warm="+warm+
					", arrays="+arrays+", key="+keys+", hashes="+hashes+", palindrome="+palindrome+", valid="+valid);
			return false;
		}
		return true;
	}

	/** Builds only through the normal rolling API; uses the actual normalized length. */
	private static Kmer make(final int requested, final int base){
		final Kmer key=new Kmer(requested);
		for(int i=0; i<key.kbig; i++){key.addRightNumeric(base);}
		return key;
	}
}
