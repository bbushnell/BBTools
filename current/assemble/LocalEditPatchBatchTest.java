package assemble;

import java.util.Arrays;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import structures.ByteBuilder;

/** Mechanical patch-transaction tests; no kmer support or biological claims.
 * @author Fischl */
public final class LocalEditPatchBatchTest {

	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Patch transaction tests take no inputs.");}
		boundaries();reverseMaterialization();rejections();
		final Random random=new Random(441207);final LocalEditPatchBatch batch=new LocalEditPatchBatch();
		for(int trial=0;trial<1000;trial++){randomCase(random,batch,trial%2==0);}
		System.out.println("LOCAL_EDIT_PATCH_BATCH_TEST_OK checks="+checks+" random_cases=1000; mechanical splicing and whole-read transaction only.");
	}
	private static void reverseMaterialization(){
		final int n=96,subPos=12,insGap=40;final byte[] truth=sequence(n,77123);final byte inserted='G';
		for(boolean withQuality:new boolean[]{false,true}){
			final Read forward=read(truth,withQuality),reverse=reverseRead(truth,withQuality);
			final LocalEditPatchBatch fbatch=new LocalEditPatchBatch(),rbatch=new LocalEditPatchBatch();
			final byte replacement=(byte)(truth[subPos]=='A' ? 'C' : 'A');
			fbatch.reset(forward);
			fbatch.add(subPos,subPos+1,new byte[]{replacement},withQuality ? new byte[]{0} : null,1,0,0);
			fbatch.add(insGap,insGap,new byte[]{inserted},withQuality ? new byte[]{0} : null,0,1,0);
			check(fbatch.commit(2)==2,"Forward control for reverse materialization must commit both edits.");
			final int reverseSub=n-subPos-1,reverseGap=n-insGap;
			rbatch.reset(reverse);
			rbatch.add(reverseGap,reverseGap,new byte[]{complement(inserted)},withQuality ? new byte[]{0} : null,0,1,0);
			rbatch.add(reverseSub,reverseSub+1,new byte[]{complement(replacement)},withQuality ? new byte[]{0} : null,1,0,0);
			check(rbatch.commit(2)==2,"Reverse-complement materialization must commit both mapped edits.");
			check(Arrays.equals(reverse.bases,AminoAcid.reverseComplementBases(forward.bases)),"Mapped reverse patches must equal the reverse complement of forward materialization.");
			check(withQuality ? Arrays.equals(reverse.quality,reversed(forward.quality)) : reverse.quality==null,"Mapped reverse patches must preserve reverse-complement quality order and nullness.");
		}
	}
	private static void boundaries(){
		final LocalEditPatchBatch batch=new LocalEditPatchBatch();
		for(boolean withQuality:new boolean[]{false,true}){
			final Read r=read("ACGT".getBytes(),withQuality);final byte[] original=r.bases,quality=r.quality;
			batch.reset(r);check(batch.commit(0)==0 && r.bases==original && r.quality==quality,"An empty transaction must not replace original arrays.");
			batch.reset(r);
			batch.add(0,0,new byte[]{'T'},withQuality ? new byte[]{0} : null,0,1,0);
			batch.add(4,4,new byte[]{'A','C'},withQuality ? new byte[]{0,0} : null,0,2,0);
			check(batch.size()==2 && batch.editCount()==3,"Two patch intervals can contain three edits; limits count operations.");
			check(batch.commit(3)==3 && Arrays.equals(r.bases,"TACGTAC".getBytes()),"Endpoint insertions must retain every original base in order.");
			check(batch.copiedOriginalBases==4 && batch.copiedReplacementBases==3 && batch.materializedBases==7,"Copy counters distinguish original spans from staged replacements.");
			if(withQuality){check(Arrays.equals(r.quality,new byte[]{0,10,11,12,13,0,0}),"Inserted endpoint quality bytes and retained original qualities must be exact.");}
			check(Arrays.equals(original,"ACGT".getBytes()),"Commit must not modify an old base buffer retained by the caller.");
			final Read deleted=read("ACGT".getBytes(),withQuality);batch.reset(deleted);
			batch.add(0,4,new byte[0],withQuality ? new byte[0] : null,0,0,4);
			check(batch.commit(4)==4 && deleted.length()==0 && (deleted.quality==null)==!withQuality,"A complete deletion is a valid empty output, with quality nullness preserved.");
			final Read empty=read(new byte[0],withQuality);batch.reset(empty);
			batch.add(0,0,new byte[]{'G'},withQuality ? new byte[]{0} : null,0,1,0);
			check(batch.commit(1)==1 && empty.bases[0]=='G',"An insertion into an empty original interval is unambiguous.");
		}
	}
	private static void rejections(){
		final LocalEditPatchBatch batch=new LocalEditPatchBatch();final Read r=read("ACGT".getBytes(),true);
		final byte[] before=r.bases,qbefore=r.quality;
		batch.reset(r);batch.add(0,1,new byte[]{'T','G'},new byte[]{0,0},1,1,0);batch.add(3,4,new byte[]{'C'},new byte[]{0},1,0,0);
		check(batch.commit(2)==-1 && r.bases==before && r.quality==qbefore,"Exceeding the global operation budget must discard all patches, not merely the last one.");
		check(batch.materializedBases==0 && batch.substitutions+batch.insertions+batch.deletions==0 && batch.size()==0,"Rejected transactions publish no output or committed counters.");
		expect(()->batch.commit(4),"A rejected commit closes the transaction.");
		batch.reset(r);batch.add(1,3,new byte[]{'A','A'},new byte[]{0,0},2,0,0);
		expect(()->batch.add(2,4,new byte[]{'C','C'},new byte[]{0,0},2,0,0),"Overlapping original intervals must be jointly resolved, not independently staged.");
		expect(()->batch.commit(8),"An append failure must abandon the whole transaction.");
		check(r.bases==before && r.quality==qbefore,"Invalid append and closed commit must not touch the original read.");
		batch.reset(r);expect(()->batch.add(0,1,new byte[]{'A'},null,1,0,0),"Quality nullness cannot change across a replacement.");
		batch.reset(r);expect(()->batch.add(0,1,new byte[]{'A'},new byte[]{0},0,1,0),"Net replacement length must agree with insertion/deletion counts.");
		batch.reset(r);expect(()->batch.add(0,1,new byte[]{'A'},new byte[]{0},0,0,0),"Zero-operation replacements should not create phantom budget-free work.");
		batch.reset(r);batch.add(0,1,new byte[]{'A'},new byte[]{0},Integer.MAX_VALUE,0,0);
		expect(()->batch.add(2,3,new byte[]{'G'},new byte[]{0},1,0,0),"Total operation counts must not wrap at Integer.MAX_VALUE.");
		batch.reset(r);batch.add(1,2,new byte[]{'T'},new byte[]{0},1,0,0);
		final byte[] external=r.bases.clone();r.bases=external;
		expect(()->batch.commit(8),"Stale original base references must be detected before installation.");
		check(r.bases==external && r.quality==qbefore,"A stale-state error must not overwrite the external replacement that triggered it.");
		batch.reset(r);batch.add(1,2,new byte[]{'T'},new byte[]{0},1,0,0);r.chrom=1;
		expect(()->batch.commit(8),"New coordinate metadata must not be silently invalidated by commit.");
		check(r.bases==external && batch.materializedBases==0,"Metadata rejection installs no partial result.");r.chrom=-1;
		batch.reset(r);batch.add(1,2,new byte[]{'T'},new byte[]{0},1,0,0);
		check(batch.commit(1)==1,"The same worker must remain reusable after failures.");
	}
	private static void randomCase(final Random random,final LocalEditPatchBatch batch,final boolean withQuality){
		final byte[] bases=new byte[random.nextInt(41)];for(int i=0;i<bases.length;i++){bases[i]=ALPHABET[random.nextInt(4)];}
		final Read r=read(bases,withQuality);final byte[] before=r.bases,qbefore=r.quality,saved=before.clone(),qsaved=qbefore==null ? null : qbefore.clone();
		final ByteBuilder expected=new ByteBuilder(),expectedQ=new ByteBuilder();
		batch.reset(r);int cursor=0,next=0,operations=0,subs=0,ins=0,dels=0,patches=0;
		while(next<=r.length() && patches<12){
			final int from=next+random.nextInt(3);if(from>r.length()){break;}
			final int to=from+random.nextInt(Math.min(3,r.length()-from)+1);
			final Read patch=new Read(Arrays.copyOfRange(before,from,to),qbefore==null ? null : Arrays.copyOfRange(qbefore,from,to),"scratch",0,false);
			final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();int s=0,a=0,d=0;
			for(int j=1+random.nextInt(3);j>0;j--){
				final int type=patch.length()==0 ? 1 : random.nextInt(3);
				if(type==0){final int p=random.nextInt(patch.length());editor.apply(patch,LocalSingleBaseEdit.Operation.SUBSTITUTION,p,patch.bases[p]=='A' ? (byte)'C' : (byte)'A');s++;}
				else if(type==1){editor.apply(patch,LocalSingleBaseEdit.Operation.INSERTION,random.nextInt(patch.length()+1),ALPHABET[random.nextInt(4)]);a++;}
				else{editor.apply(patch,LocalSingleBaseEdit.Operation.DELETION,random.nextInt(patch.length()));d++;}
			}
			expected.append(before,cursor,from-cursor).append(patch.bases);
			if(qbefore!=null){expectedQ.append(qbefore,cursor,from-cursor).append(patch.quality);}
			batch.add(from,to,patch.bases,patch.quality,s,a,d);
			Arrays.fill(patch.bases,(byte)'N');if(patch.quality!=null){Arrays.fill(patch.quality,(byte)99);}
			operations+=s+a+d;subs+=s;ins+=a;dels+=d;patches++;cursor=to;next=Math.max(to,from+1);
		}
		expected.append(before,cursor,before.length-cursor);if(qbefore!=null){expectedQ.append(qbefore,cursor,qbefore.length-cursor);}
		check(batch.size()==patches && batch.editCount()==operations,"Staged primitive metadata and actual local operation counts must agree.");
		check(batch.commit(operations)==operations,"Exactly-at-limit patch batches must commit in full.");
		check(Arrays.equals(r.bases,expected.toBytes()) && (qbefore==null ? r.quality==null : Arrays.equals(r.quality,expectedQ.toBytes())),"Backward materialization must match the independent forward splice, despite caller scratch reuse.");
		check(Arrays.equals(before,saved) && Arrays.equals(qbefore,qsaved),"Neither successful staging nor commit may mutate original buffers.");
		check(batch.substitutions==subs && batch.insertions==ins && batch.deletions==dels,"Committed counters preserve every sequential local operation, including cancelling indels.");
		check(patches==0 ? r.bases==before : batch.copiedOriginalBases+batch.copiedReplacementBases==r.length(),"Final array coverage must be exact; an empty batch preserves identity.");
	}
	private static Read read(final byte[] bases,final boolean qualities){final byte[] q=qualities ? new byte[bases.length] : null;if(q!=null){for(int i=0;i<q.length;i++){q[i]=(byte)(10+i%30);}}return new Read(bases.clone(),q,"patch-test",0,false);}
	private static byte[] sequence(final int length,final long seed){final Random random=new Random(seed);final byte[] b=new byte[length];for(int i=0;i<length;i++){b[i]=ALPHABET[random.nextInt(4)];}b[0]='A';b[length-1]='T';return b;}
	private static Read reverseRead(final byte[] bases,final boolean qualities){final Read r=read(bases,qualities);AminoAcid.reverseComplementBasesInPlace(r.bases);if(r.quality!=null){r.quality=reversed(r.quality);}return r;}
	private static byte[] reversed(final byte[] values){final byte[] out=values.clone();for(int a=0,b=out.length-1;a<b;a++,b--){final byte v=out[a];out[a]=out[b];out[b]=v;}return out;}
	private static byte complement(final byte b){return b=='A' ? (byte)'T' : b=='C' ? (byte)'G' : b=='G' ? (byte)'C' : (byte)'A';}
	private static void expect(final Runnable action,final String why){boolean failed=false;try{action.run();}catch(IllegalArgumentException|IllegalStateException expected){failed=true;}check(failed,why);}
	private static void check(final boolean pass,final String why){checks++;if(!pass){throw new AssertionError(why);}}
	private static int checks;
	private static final byte[] ALPHABET={'A','C','G','T'};
}
