package assemble;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/** DEV integration controls for supplied local regions; not automatic grouping.
 * @author Fischl */
public final class LocalEditPatchCorrectorTest {

	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Deterministic patch integration test takes no inputs.");}
		final boolean packed=Kmer.PACKED,mask=Kmer.MASK_CORE;
		try{
			Kmer.PACKED=true;Kmer.MASK_CORE=false;
			for(int k:new int[]{31,62,95,127}){for(boolean reverse:new boolean[]{false,true}){for(boolean quality:new boolean[]{false,true}){pairedRegions(k,reverse,quality);}}}
			failures();
			System.out.println("LOCAL_EDIT_PATCH_CORRECTOR_TEST_OK checks="+checks+"; supplied-region integration, no automatic grouping or benchmark claim.");
		}finally{Kmer.PACKED=packed;Kmer.MASK_CORE=mask;}
	}
	private static void pairedRegions(final int k,final boolean reverse,final boolean quality){
		final byte[] truth=sequence(20*k,1603+k);final Oracle table=new Oracle(k,truth);
		final Read input=noisy(truth,k,quality);if(reverse){reverse(input);}
		final IntList regions=regions(input.length(),k,reverse);
		final Read expected=copy(input),actual=copy(input);
		final byte[] before=actual.bases,qbefore=actual.quality,saved=before.clone(),qsaved=qbefore==null ? null : qbefore.clone();
		final LocalEditPatchCorrector patch=new LocalEditPatchCorrector(k,table,1,false,8*k);
		check(new LocalEditEngine(k,table::count,1).correct(expected,8)==4,"The whole-read sequential fixture must recover all four known S/I/D errors.");
		check(patch.correct(actual,regions,4)==4 && patch.stagedPatches==2,"Two independent local regions must stage all four repairs at the exact global limit.");
		check(Arrays.equals(actual.bases,reverse ? AminoAcid.reverseComplementBases(truth) : truth),"One final patch transaction must recover the full truth sequence.");
		check(Arrays.equals(actual.bases,expected.bases) && Arrays.equals(actual.quality,expected.quality),"Local-region integration must match whole-read sequential bases and qualities, including reverse/null-quality inputs.");
		check(Arrays.equals(before,saved) && Arrays.equals(qbefore,qsaved),"Local integration must not mutate original array aliases.");
		check(patch.substitutions==2 && patch.insertions==1 && patch.deletions==1 && patch.materializedBases==actual.length(),"Only committed local operations and the one final materialization enter output counters.");
		for(int limit:new int[]{2,3}){
			final Read rejected=copy(input);final byte[] b=rejected.bases,q=rejected.quality;
			check(patch.correct(rejected,regions,limit)==0 && patch.budgetRejected && !patch.initialRejected && patch.stagedPatches==1,"Later local edits must trigger whole-read rollback, including zero remaining budget after the first patch.");
			check(rejected.bases==b && rejected.quality==q && patch.materializedBases==0 && patch.substitutions+patch.insertions+patch.deletions==0,"A read rejected after staging its first region must publish no partial sequence or edit counters.");
		}
		final Read initial=copy(input);final byte[] ib=initial.bases;
		check(patch.correct(initial,regions,1)==0 && patch.initialRejected && initial.bases==ib,"The original whole-read burden preflight must run before local patches receive budgets.");
		check(patch.correct(copy(input),regions,8)==4,"The same worker must be reusable after both preflight and later-budget rejection.");
	}
	private static void failures(){
		final int k=31,p=4*k,q=p+k+1;final byte[] truth=sequence(20*k,1603+k);final Oracle table=new Oracle(k,truth);
		final LocalEditPatchCorrector patch=new LocalEditPatchCorrector(k,table,1,false,8*k);
		final Read input=noisy(truth,k,true);final byte[] b=input.bases,quals=input.quality;
		final IntList restricted=new IntList();add(restricted,p-2*k,q+2*k,p-2,p+2);
		check(patch.correct(input,restricted,8)==0 && patch.boundaryRejected==1 && input.bases==b && input.quality==quals,"Discovering an edit in padding must discard that local candidate, not change a neighboring region.");
		final IntList huge=new IntList();add(huge,0,input.length(),k,input.length()-k);
		check(patch.correct(input,huge,8)==0 && patch.sizeRejected==1 && input.bases==b,"Scratch work is bounded by the explicit patch-size ceiling.");
		final IntList invalid=new IntList();add(invalid,10,20,Integer.MIN_VALUE,Integer.MIN_VALUE+100);
		expect(()->patch.correct(input,invalid,8),"Malformed core coordinates must fail before arithmetic can wrap into apparently valid flanks.");
		final IntList overlapping=regions(input.length(),k,false);add(overlapping,0,3*k,k,2*k);
		expect(()->patch.correct(input,overlapping,8),"Supplied copy regions must be in nonoverlapping original order.");
		table.fail=key(truth,12*k-k/2,k);
		boolean failed=false;try{patch.correct(input,regions(input.length(),k,false),8);}catch(IllegalStateException expected){failed=true;}
		check(failed && patch.stagedPatches==1,"Injected lookup failure must occur after the first patch was staged, not only in preflight.");
		check(input.bases==b && input.quality==quals && patch.materializedBases==0,"A later-region exception leaves the complete original read installed.");
		table.fail=null;check(patch.correct(input,regions(input.length(),k,false),8)==4,"A worker recovers after abandoning a failed multi-region transaction.");
	}
	private static IntList regions(final int length,final int k,final boolean reverse){
		final IntList out=new IntList();
		if(reverse){region(out,length,k,12*k,true);region(out,length,k,4*k,true);}
		else{region(out,length,k,4*k,false);region(out,length,k,12*k,false);}
		return out;
	}
	private static void region(final IntList out,final int length,final int k,final int p,final boolean reverse){
		final int a=p-4,b=p+k+6,from=a-k,to=b+k;
		if(reverse){add(out,length-to,length-from,length-b,length-a);}else{add(out,from,to,a,b);}
	}
	private static void add(final IntList out,final int from,final int to,final int coreFrom,final int coreTo){out.add(from);out.add(to);out.add(coreFrom);out.add(coreTo);}
	private static Read noisy(final byte[] truth,final int k,final boolean quality){
		final byte[] qs=quality ? new byte[truth.length] : null;if(qs!=null){for(int i=0;i<qs.length;i++){qs[i]=(byte)(15+i%25);}}
		final Read r=new Read(truth.clone(),qs,"two-local-regions",0,false);final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();
		editor.apply(r,LocalSingleBaseEdit.Operation.SUBSTITUTION,13*k+1,other(r.bases[13*k+1]));
		editor.apply(r,LocalSingleBaseEdit.Operation.DELETION,12*k);
		editor.apply(r,LocalSingleBaseEdit.Operation.INSERTION,5*k+1,r.bases[5*k+1]);
		editor.apply(r,LocalSingleBaseEdit.Operation.SUBSTITUTION,4*k,other(r.bases[4*k]));return r;
	}
	private static byte other(final byte b){return b=='A' ? (byte)'C' : (byte)'A';}
	private static Read copy(final Read r){return new Read(r.bases.clone(),r.quality==null ? null : r.quality.clone(),r.id,r.numericID,false);}
	private static void reverse(final Read r){AminoAcid.reverseComplementBasesInPlace(r.bases);if(r.quality!=null){for(int a=0,b=r.quality.length-1;a<b;a++,b--){final byte v=r.quality[a];r.quality[a]=r.quality[b];r.quality[b]=v;}}}
	private static byte[] sequence(final int length,final long seed){final Random random=new Random(seed);final byte[] b=new byte[length];for(int i=0;i<length;i++){b[i]=ALPHABET[random.nextInt(4)];}b[0]=b[length-1]='A';return b;}
	private static String key(final byte[] b,final int start,final int k){final Kmer word=new Kmer(k);for(int i=start;i<start+k;i++){word.addRight(b[i]);}return Arrays.toString(word.key());}
	private static final class Oracle implements HomopolymerIndelProposal.CountLookup {
		Oracle(final int k,final byte[] truth){for(int i=0;i<=truth.length-k;i++){counts.put(key(truth,i,k),12);}}
		@Override public int count(final Kmer word){final String s=Arrays.toString(word.key());if(s.equals(fail)){throw new IllegalStateException("Injected later-region lookup failure.");}final Integer n=counts.get(s);return n==null ? 0 : n;}
		final HashMap<String,Integer> counts=new HashMap<String,Integer>();String fail;
	}
	private static void expect(final Runnable action,final String why){boolean failed=false;try{action.run();}catch(IllegalArgumentException expected){failed=true;}check(failed,why);}
	private static void check(final boolean pass,final String why){checks++;if(!pass){throw new AssertionError(why);}}
	private static int checks;
	private static final byte[] ALPHABET={'A','C','G','T'};
}
