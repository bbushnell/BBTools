package test.bloom;

import bloom.BloomFilterCorrectorWrapper;
import stream.Read;

/** Execute an enabled native invocation, then check its shared policy is restored.
 * Supply normal wrapper args with localedit=t, ecc=f, merge=f, no output required.
 * @author Fischl */
public final class BloomLocalEditScopeTest{
	/** Invokes the real wrapper with both read-policy flags initially true and checks
	 * that a successful local-edit invocation restores them. Finally restores the
	 * caller's original flags even if the wrapper or check fails; other globals remain.
	 * @param args Wrapper arguments selecting a valid local-edit invocation */
	public static void main(final String[] args){
		final boolean change=Read.CHANGE_QUALITY, header=Read.FIX_HEADER;
		try{
			Read.CHANGE_QUALITY=true; Read.FIX_HEADER=true;
			BloomFilterCorrectorWrapper.main(args);
			if(!Read.CHANGE_QUALITY || !Read.FIX_HEADER){throw new AssertionError("Experimental byte-preservation policy leaked into a later same-JVM invocation.");}
			System.out.println("BLOOM_LOCAL_EDIT_SCOPE_TEST_OK");
		}finally{Read.CHANGE_QUALITY=change; Read.FIX_HEADER=header;}
	}
}
