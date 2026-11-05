package test.kmer;

import kmer.KmerLink;

/**
 * Finite direct-link return, creation and saturation checks for KMR009/KMR019.
 * Each chain has at most three literal keys. KmerTable's separate counting path
 * is not used. Storage and return oracles are checked separately under either
 * assertion mode, with one selected case per JVM and no I/O or large allocation.
 * @author Sucrose
 */
public final class KmerLinkContractTest{

	/** Selects one failure target or a nearby control using case=name. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=name argument.");
		}
		final String mode=args[0].substring(5);
		switch(mode){
			case "new-count": insert(mode, false, 2); break;
			case "new-created": insert(mode, true, 2); break;
			case "zero-created": update(mode, true, 0, 1, 1); break;
			case "overflow-count": update(mode, false, Integer.MAX_VALUE-1, 2, Integer.MAX_VALUE); break;
			case "overflow-created": update(mode, true, Integer.MAX_VALUE-1, 2, Integer.MAX_VALUE); break;
			case "ordinary-control":
				update(mode, false, 7, 2, 9);
				update(mode, true, 7, 2, 9);
				break;
			case "new-one-control":
				insert(mode, false, 1);
				insert(mode, true, 1);
				break;
			case "saturated-control":
				update(mode, false, Integer.MAX_VALUE, 2, Integer.MAX_VALUE);
				update(mode, true, Integer.MAX_VALUE, 2, Integer.MAX_VALUE);
				break;
			default: throw new IllegalArgumentException("Unknown link case: "+mode);
		}
		System.out.println("PASS KmerLink case="+mode);
	}

	/** Creates head(10,5) and tail(20,tailCount) through public construction/set APIs. */
	private static KmerLink chain(final String mode, final int tailCount){
		final KmerLink root=new KmerLink(10, 5);
		check(root.set(20, tailCount)==1 && root.contains(20), mode, "Fixture tail must be a newly present key.");
		check(root.getValue(10)==5 && root.getValue(20)==tailCount && !root.contains(30),
			mode, "Fixture must contain only the expected head and tail before the operation.");
		return root;
	}

	/** Appends key30 after two existing links; initial count and creation flag are independent oracles. */
	private static void insert(final String mode, final boolean createdApi, final int incr){
		final KmerLink root=chain(mode, 7);
		final int result=createdApi ? root.incrementAndReturnNumCreated(30, incr) : root.increment(30, incr);
		check(root.getValue(10)==5 && root.getValue(20)==7 && root.getValue(30)==incr
			&& root.size()==3, mode, "New link must store its increment and preserve both existing counts.");
		if(createdApi){check(result==1, mode, "Expected created=1, got="+result);}
		else{check(result==incr, mode, "Expected return="+incr+", got="+result);}
	}

	/** Updates an existing tail, including zero and near-MAX counts, without reporting creation.
	 * expectedCount is supplied literally by the case rather than computed from production output. */
	private static void update(final String mode, final boolean createdApi, final int initial,
			final int incr, final int expectedCount){
		final KmerLink root=chain(mode, initial);
		final int result=createdApi ? root.incrementAndReturnNumCreated(20, incr) : root.increment(20, incr);
		final int stored=root.getValue(20);
		check(stored==expectedCount, mode, "Expected stored count="+expectedCount+", got="+stored);
		check(root.getValue(10)==5 && !root.contains(30) && root.size()==2,
			mode, "An update must preserve the head and the two-link chain.");
		if(createdApi){check(result==0, mode, "Expected created=0, got="+result);}
		else{check(result==expectedCount, mode, "Expected return="+expectedCount+", got="+result);}
	}

	/** Keeps all fixture and oracle checks active under both -ea and -da. */
	private static void check(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("KmerLink case="+mode+": "+message);}
	}
}
