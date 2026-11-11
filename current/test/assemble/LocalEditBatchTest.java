package test.assemble;

import java.util.Arrays;
import java.util.Random;
import assemble.LocalEditBatch;
import assemble.LocalSingleBaseEdit;
import stream.Read;

/**
 * Mechanical tests for original-coordinate single-base edit transactions.
 * Enumerates operation placements on reads of length zero through four with
 * fixed replacement alleles, then exercises 1500 seeded random cases. Checks
 * optional qualities, retained input buffers, budgets and selected rejection paths.
 * The differential reference applies LocalSingleBaseEdit in reverse coordinate
 * order; it shares editable-read validation with the batch implementation.
 * This does not test correction discovery, biological accuracy or concurrent use.
 * @author Fischl */
public final class LocalEditBatchTest{

	/** Runs the deterministic fixtures and prints a marker after all checks pass.
	 * Arguments are ignored; explicit checks remain active without JVM assertions. */
	public static void main(final String[] args){
		final LocalEditBatch batch=new LocalEditBatch();
		// All operation choices for prescribed alleles at every short-read position, including
		// adjacency, deletion to empty, and insertions at both terminal boundaries.
		for(int length=0; length<=4; length++){
			for(int mask=0; mask<(1<<(2*length+1)); mask++){
				final byte[] bases=new byte[length]; for(int i=0; i<length; i++){bases[i]=ALPHABET[i&3];}
				final int[] operations=new int[length+1], alleles=new int[length+1]; Arrays.fill(operations, -1);
				for(int p=0; p<length; p++){operations[p]=((mask>>>(2*p))&3)-1; alleles[p]=(p+2)&3;}
				if((mask&(1<<(2*length)))!=0){operations[length]=2; alleles[length]=3;}
				for(boolean qualities:new boolean[]{false, true}){compare(batch, bases, operations, alleles, qualities);}
			}
		}
		final Random random=new Random(902103);
		for(int trial=0; trial<1500; trial++){
			final int length=random.nextInt(128); final byte[] bases=new byte[length];
			final int[] operations=new int[length+1], alleles=new int[length+1]; Arrays.fill(operations, -1);
			for(int p=0; p<length; p++){
				bases[p]=p%17==0 ? (byte)'N' : ALPHABET[random.nextInt(4)];
				if(random.nextInt(4)==0){operations[p]=random.nextInt(3); alleles[p]=random.nextInt(4);}
			}
			if(random.nextBoolean()){operations[length]=2; alleles[length]=random.nextInt(4);}
			compare(batch, bases, operations, alleles, random.nextBoolean());
		}
		strand(); rejections(); lengths();
		System.out.println("LOCAL_EDIT_BATCH_TEST_OK checks="+checks);
	}
	/** Compares one batch against reverse-order single edits on an independent Read.
	 * Operations are -1/0/1/2 for none/substitution/deletion/insertion; alleles are
	 * A/C/G/T codes. Coordinates refer to the original bases, with a terminal gap.
	 * Also checks exact-budget acceptance and whole-list rejection one edit below it. */
	private static void compare(final LocalEditBatch batch, final byte[] bases, final int[] operations, final int[] alleles, final boolean qualities){
		final Read actual=read(bases, qualities), expected=read(bases, qualities);
		final byte[] original=actual.bases, originalQuality=actual.quality;
		final byte[] saved=original.clone(), savedQuality=originalQuality==null ? null : originalQuality.clone();
		batch.reset(actual); int count=0, consumed=0;
		for(int p=0; p<operations.length; p++){
			if(operations[p]<0){continue;}
			append(batch, operations[p], p, alleles[p]); count++; if(operations[p]!=2){consumed++;}
		}
		check(batch.size()==count, "Every accepted proposal must occupy one packed entry.");
		final LocalSingleBaseEdit oracle=new LocalSingleBaseEdit();
		for(int p=operations.length-1; p>=0; p--){
			if(operations[p]==0){oracle.apply(expected, LocalSingleBaseEdit.Operation.SUBSTITUTION, p, ALPHABET[alleles[p]]);
			}else if(operations[p]==1){oracle.apply(expected, LocalSingleBaseEdit.Operation.DELETION, p);
			}else if(operations[p]==2){oracle.apply(expected, LocalSingleBaseEdit.Operation.INSERTION, p, ALPHABET[alleles[p]]);}
		}
		check(batch.commit(count)==count, "Exactly-at-budget materialization must succeed, including an empty batch.");
		check(Arrays.equals(actual.bases, expected.bases) && Arrays.equals(actual.quality, expected.quality), "One-pass output must equal independently applied reverse-order single edits, including qualities.");
		check(Arrays.equals(original, saved) && Arrays.equals(originalQuality, savedQuality), "All original arrays must remain immutable, even on successful commit.");
		check(batch.copiedBases==(count==0 ? 0 : bases.length-consumed), "Only unchanged original bases are copied, exactly once; edited/deleted bases are not recopied.");
		check(batch.size()==0, "Commit closes and empties the transaction.");
		if(count==0){check(actual.bases==original && actual.quality==originalQuality, "Empty batch must not allocate replacement arrays.");}
		if(count>0){
			final Read rejected=read(bases, qualities); final byte[] before=rejected.bases, qbefore=rejected.quality;
			batch.reset(rejected); for(int p=0; p<operations.length; p++){if(operations[p]>=0){append(batch, operations[p], p, alleles[p]);}}
			check(batch.commit(count-1)==-1 && rejected.bases==before && rejected.quality==qbefore, "Over-budget batch rejects the complete list without replacing either original array.");
			check(batch.copiedBases==0 && batch.size()==0, "Budget rejection performs no sequence copies and leaves no pending edits.");
		}
	}
	/** Maps three prescribed edits onto the opposite strand and compares bases and qualities. */
	private static void strand(){
		final byte[] bases="ACGTACGTAACCGGTTACGTACGT".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final Read forward=read(bases, true), reverse=read(rc(bases), true);
		reverse.quality=reverse(forward.quality);
		final LocalEditBatch a=new LocalEditBatch(), b=new LocalEditBatch(); a.reset(forward); b.reset(reverse);
		a.addSubstitution(2, 0); a.addDeletion(9); a.addInsertion(18, 1);
		b.addInsertion(bases.length-18, 2); b.addDeletion(bases.length-1-9); b.addSubstitution(bases.length-1-2, 3);
		check(a.commit(3)==3 && b.commit(3)==3, "Independent strand-transformed edits must fit their identical budget.");
		check(Arrays.equals(rc(forward.bases), reverse.bases) && Arrays.equals(reverse(forward.quality), reverse.quality), "Independent edits must preserve reverse-complement sequence and reversed quality symmetry.");
	}
	/** Exercises invalid domains/order, closed batches, stale arrays and attached mates.
	 * These are selected validation paths, not allocation-failure or race injection. */
	private static void rejections(){
		final LocalEditBatch batch=new LocalEditBatch(); final Read r=read(new byte[]{'A', 'N', 'T'}, true);
		final byte[] bases=r.bases, quality=r.quality;
		batch.reset(r); batch.addDeletion(1); fails(()->batch.addInsertion(1, 0)); fails(()->batch.addSubstitution(0, 1));
		check(batch.size()==1 && r.bases==bases && r.quality==quality, "Rejected unordered/conflicting proposals must not alter the read or prior list.");
		batch.clear(); check(batch.size()==0 && r.bases==bases, "Discarding a list is rollback without a copy.");
		fails(()->batch.addDeletion(0)); fails(()->batch.commit(1));
		batch.reset(r); fails(()->batch.addDeletion(-1)); fails(()->batch.addDeletion(3)); fails(()->batch.addSubstitution(3, 0));
		fails(()->batch.addInsertion(4, 0)); fails(()->batch.addInsertion(0, -1)); fails(()->batch.addInsertion(0, 4));
		check(batch.size()==0, "Invalid positions/base codes never append a record.");
		batch.addInsertion(0, 0); fails(()->batch.commit(-1));
		check(r.bases==bases && r.quality==quality && batch.size()==0, "Commit exception must retain originals and close the transaction.");
		batch.reset(r); batch.addDeletion(1); r.bases=bases.clone(); final byte[] replacement=r.bases;
		fails(()->batch.commit(1)); check(r.bases==replacement && r.quality==quality, "A stale transaction must never overwrite another owner's replacement array.");
		r.bases=bases; batch.reset(r); batch.addDeletion(1); r.quality=quality.clone(); final byte[] qreplacement=r.quality;
		fails(()->batch.commit(1)); check(r.bases==bases && r.quality==qreplacement, "A stale quality reference must also block installation.");
		r.quality=new byte[1]; fails(()->batch.reset(r)); check(batch.size()==0, "Invalid reset cannot leave an old list active.");
		r.quality=quality; batch.reset(r); batch.addDeletion(1); r.mate=read(new byte[]{'A'}, true);
		fails(()->batch.commit(1)); check(r.bases==bases && r.quality==quality, "Attached mate metadata must reject the batch before installation.");
	}
	/** Checks length arithmetic at empty and overflow boundaries without large allocations. */
	private static void lengths(){
		check(LocalEditBatch.resultLength(0, 0)==0 && LocalEditBatch.resultLength(3, -3)==0, "Deletion to an empty sequence is representable.");
		check(LocalEditBatch.resultLength(shared.Shared.MAX_ARRAY_LEN, 0)==shared.Shared.MAX_ARRAY_LEN, "Length validation uses the shared Java array bound without allocating a giant test array.");
		fails(()->LocalEditBatch.resultLength(-1, 0)); fails(()->LocalEditBatch.resultLength(0, -1));
		fails(()->LocalEditBatch.resultLength(1, Long.MAX_VALUE)); fails(()->LocalEditBatch.resultLength(1, Long.MIN_VALUE));
		fails(()->LocalEditBatch.resultLength(shared.Shared.MAX_ARRAY_LEN, 1));
	}
	/** Dispatches one fixture operation to the batch's typed recording methods. */
	private static void append(final LocalEditBatch batch, final int operation, final int position, final int allele){
		if(operation==0){batch.addSubstitution(position, allele);}else if(operation==1){batch.addDeletion(position);}else if(operation==2){batch.addInsertion(position, allele);}else{throw new AssertionError("Invalid test operation.");}
	}
	/** Clones bases into an unpaired Read, optionally with qualities cycling from zero to forty. */
	private static Read read(final byte[] bases, final boolean withQuality){
		final byte[] q=withQuality ? new byte[bases.length] : null;
		if(q!=null){for(int i=0; i<q.length; i++){q[i]=(byte)(i%41);}}
		return new Read(bases.clone(), q, "batch-test", 0, false);
	}
	/** Returns reversed bytes without changing the supplied array. */
	private static byte[] reverse(final byte[] a){final byte[] b=a.clone(); for(int i=0; i<a.length; i++){b[i]=a[a.length-1-i];}return b;}
	/** Reverse-complements an uppercase A/C/G/T fixture; rejects other symbols. */
	private static byte[] rc(final byte[] a){final byte[] b=reverse(a); for(int i=0; i<b.length; i++){switch(b[i]){case 'A':b[i]='T'; break; case 'T':b[i]='A'; break; case 'C':b[i]='G'; break; case 'G':b[i]='C'; break; default:throw new AssertionError("Defined strand fixture expected.");}}return b;}
	/** Requires an argument/state exception; other failures propagate unchanged. */
	private static void fails(final Runnable action){boolean failed=false; try{action.run();}catch(IllegalArgumentException | IllegalStateException expected){failed=true;}check(failed, "Invalid batch operation must fail loudly.");}
	/** Counts and enforces a fixture condition regardless of assertion settings. */
	private static void check(final boolean ok, final String why){checks++; if(!ok){throw new AssertionError(why);}}
	/** Cumulative check invocations in this JVM; main does not reset this counter. */
	private static long checks;
	/** Numeric-allele decoding shared by fixture generation and the differential reference. */
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
