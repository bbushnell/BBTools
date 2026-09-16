package assemble;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/** DEV-only automatic region geometry and native correction integration.
 * @author Fischl */
public final class LocalEditPatchRegionsTest {

	public static void main(final String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Region tests take no inputs.");}
		final boolean packed=Kmer.PACKED,mask=Kmer.MASK_CORE;
		try{
			Kmer.PACKED=true;Kmer.MASK_CORE=false;geometry();edgeGeometry();
			for(int k:new int[]{31,62,95,127}){for(boolean reverse:new boolean[]{false,true}){for(boolean qualities:new boolean[]{false,true}){integration(k,reverse,qualities);}}}
			for(int k:new int[]{31,62,95,127}){for(boolean reverse:new boolean[]{false,true}){for(boolean right:new boolean[]{false,true}){edgeSubstitution(k,reverse,right);}}}
			System.out.println("LOCAL_EDIT_PATCH_REGIONS_TEST_OK checks="+checks+"; automatic original-coordinate grouping, no population performance claim.");
		}finally{Kmer.PACKED=packed;Kmer.MASK_CORE=mask;}
	}
	private static void edgeGeometry(){
		final int k=5,length=200;final LocalEditPatchRegions.Regions regions=new LocalEditPatchRegions.Regions(k,100);
		for(final LocalSingleBaseEdit.Operation op:LocalSingleBaseEdit.Operation.values()){
			final boolean insertion=op==LocalSingleBaseEdit.Operation.INSERTION;
			for(final int p:new int[]{0,k-1,k,k+1,length-k-1,length-k,length-k+1,length}){
				if(p==length && !insertion){continue;}
				final boolean eligible=p>=k && (insertion ? p<=length-k : p<length-k);
				int[] forward=null;
				for(final boolean reverse:new boolean[]{false,true}){
					regions.reset(length);final int contextFrom=Math.max(0,p-9),contextTo=Math.min(length,p+10);
					regions.add(op,p,contextFrom,contextTo,reverse);final IntList out=regions.finish();
					check((out.size==4)==eligible && regions.edgeDeferred==(eligible ? 0 : 1),"Optional slack may shrink, but complete S/D event or insertion gap plus K flanks must fit the original read.");
					if(!eligible){continue;}
					final int from=out.array[0],to=out.array[1],a=out.array[2],b=out.array[3],point=reverse ? length-p-(insertion ? 0 : 1) : p;
					check(from>=0 && to<=length && a-from>=k && to-b>=k && a<=point && (insertion ? point<=b : point<b),"Clipped optional padding must not remove an immutable flank or exclude its normalized event.");
					check(from<=(reverse ? length-contextTo : contextFrom) && to>=(reverse ? length-contextFrom : contextTo),"Every original native support-context base remains present after optional-slack trimming.");
					if(!reverse){forward=values(out);}else{
						check(forward[0]==length-to && forward[1]==length-from && forward[2]==length-b && forward[3]==length-a,"Edge trimming and interval reflection must commute for S/I/D, including gap endpoints.");
					}
				}
			}
		}
	}
	private static void edgeSubstitution(final int k,final boolean reverse,final boolean right){
		final Random random=new Random(2701+k);final byte[] truth=new byte[6*k];
		for(int i=0;i<truth.length;i++){truth[i]=ALPHABET[random.nextInt(4)];}truth[0]=truth[truth.length-1]='A';
		final Oracle table=new Oracle(k,truth);final byte[] q=new byte[truth.length];for(int i=0;i<q.length;i++){q[i]=(byte)(20+i%20);}
		final Read input=new Read(truth.clone(),q,"optional-edge-slack",0,false);
		final int p=right ? truth.length-k-1 : k+1;input.bases[p]=other(input.bases[p]);
		if(reverse){AminoAcid.reverseComplementBasesInPlace(input.bases);reverse(input.quality);}
		final byte[] original=input.bases,quality=input.quality,saved=original.clone(),savedQ=quality.clone();
		final Read expected=new Read(saved.clone(),savedQ.clone(),input.id,0,false);
		final LocalEditPatchRegions selector=new LocalEditPatchRegions(k,table,1,false,1,8*k);
		final LocalEditPatchCorrector executor=new LocalEditPatchCorrector(k,table,1,false,8*k);
		check(new LocalEditEngine(k,table::count,1).correct(expected,1)==1,"The native two-supported-flank kernel must accept the edge-slack fixture before the selector may recover it.");
		check(executor.correctSelected(input,selector,1)==1 && selector.regions.edgeDeferred==0,"An optional out-of-read slack base must not withhold a native-supported substitution with complete immutable K flanks.");
		check(Arrays.equals(input.bases,expected.bases) && Arrays.equals(input.quality,expected.quality),"Clipped-slack correction must match full native sequential bases AND qualities on either strand.");
		check(Arrays.equals(input.bases,reverse ? AminoAcid.reverseComplementBases(truth) : truth),"Edge-slack repairs must recover known truth, not only another implementation's output.");
		check(Arrays.equals(original,saved) && Arrays.equals(quality,savedQ),"Edge correction must leave original array aliases immutable.");
	}
	private static void geometry(){
		final LocalEditPatchRegions.Regions r=new LocalEditPatchRegions.Regions(5,100);
		r.reset(200);
		r.add(LocalSingleBaseEdit.Operation.SUBSTITUTION,100,92,109,false);
		r.add(LocalSingleBaseEdit.Operation.DELETION,20,10,29,false);
		r.add(LocalSingleBaseEdit.Operation.INSERTION,30,23,38,false);
		check(Arrays.equals(values(r.finish()),new int[]{10,38,18,32,92,109,98,103}),"Unsorted and overlapping proposals must yield sorted union copy intervals and contained mutable cores.");
		check(r.merged==1 && r.sizeDeferred==0 && r.edgeDeferred==0,"Geometry counters must distinguish merges from abstentions.");
		check(r.finish()==r.finish(),"Finishing twice must not duplicate output.");
		expect(()->r.add(LocalSingleBaseEdit.Operation.SUBSTITUTION,40,30,50,false),"A closed collection cannot silently accept more proposals.");
		r.reset(200);r.add(LocalSingleBaseEdit.Operation.INSERTION,20,10,30,true);
		check(Arrays.equals(values(r.finish()),new int[]{170,190,178,182}),"A reverse insertion maps a gap at length-p, not a base at length-p-1.");
		r.reset(200);r.add(LocalSingleBaseEdit.Operation.SUBSTITUTION,20,10,30,true);
		check(Arrays.equals(values(r.finish()),new int[]{170,190,177,182}),"Reverse base-consuming events map their whole core interval, preserving half-open endpoints.");
		r.reset(200);r.add(LocalSingleBaseEdit.Operation.DELETION,1,0,10,false);
		check(r.finish().size==0 && r.edgeDeferred==1,"Unsupported edge flanks are explicit abstentions, not clipped regions.");
		final LocalEditPatchRegions.Regions bounded=new LocalEditPatchRegions.Regions(5,25);bounded.reset(200);
		bounded.add(LocalSingleBaseEdit.Operation.SUBSTITUTION,20,10,30,false);
		bounded.add(LocalSingleBaseEdit.Operation.SUBSTITUTION,40,30,50,false);
		check(bounded.finish().size==0 && bounded.merged==1 && bounded.sizeDeferred==1,"Touching components must be merged before the size check and never split to evade the bound.");
		r.reset(Integer.MAX_VALUE);r.add(LocalSingleBaseEdit.Operation.INSERTION,Integer.MAX_VALUE,Integer.MAX_VALUE-10,Integer.MAX_VALUE,false);
		check(r.finish().size==0 && r.edgeDeferred==1,"Scratch endpoint arithmetic must not wrap near the maximum original coordinate.");
		final Random random=new Random(81273);
		for(int trial=0;trial<200;trial++){
			final LocalEditPatchRegions.Regions forward=new LocalEditPatchRegions.Regions(5,1000),reverse=new LocalEditPatchRegions.Regions(5,1000);
			forward.reset(1000);reverse.reset(1000);
			for(int i=0;i<20;i++){
				final int p=20+random.nextInt(960),from=p-10-random.nextInt(10),to=p+10+random.nextInt(10);
				final LocalSingleBaseEdit.Operation op=LocalSingleBaseEdit.Operation.values()[random.nextInt(3)];
				forward.add(op,p,from,to,false);reverse.add(op,p,from,to,true);
			}
			final IntList f=forward.finish(),v=reverse.finish();check(f.size==v.size,"Strand reflection must preserve connected component count.");
			for(int i=0;i<f.size;i+=4){final int j=v.size-4-i;
				check(f.array[i]==1000-v.array[j+1] && f.array[i+1]==1000-v.array[j] && f.array[i+2]==1000-v.array[j+3] && f.array[i+3]==1000-v.array[j+2],"Sorting and union must commute with strand reflection, including insertion cores.");
			}
		}
	}
	private static void integration(final int k,final boolean reverse,final boolean qualities){
		final Random random=new Random(1603+k);final byte[] truth=new byte[20*k];
		for(int i=0;i<truth.length;i++){truth[i]=ALPHABET[random.nextInt(4)];}truth[0]=truth[truth.length-1]='A';
		final Oracle table=new Oracle(k,truth);final byte[] qs=qualities ? new byte[truth.length] : null;
		if(qs!=null){for(int i=0;i<qs.length;i++){qs[i]=(byte)(15+i%25);}}
		final Read input=new Read(truth.clone(),qs,"automatic-regions",0,false);final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();
		editor.apply(input,LocalSingleBaseEdit.Operation.SUBSTITUTION,13*k+1,other(input.bases[13*k+1]));
		editor.apply(input,LocalSingleBaseEdit.Operation.DELETION,12*k);
		editor.apply(input,LocalSingleBaseEdit.Operation.INSERTION,5*k+1,input.bases[5*k+1]);
		editor.apply(input,LocalSingleBaseEdit.Operation.SUBSTITUTION,4*k,other(input.bases[4*k]));
		if(reverse){AminoAcid.reverseComplementBasesInPlace(input.bases);if(qs!=null){reverse(input.quality);}}
		final byte[] original=input.bases,quality=input.quality,bcopy=original.clone(),qcopy=quality==null ? null : quality.clone();
		final Read expected=new Read(bcopy.clone(),qcopy==null ? null : qcopy.clone(),input.id,0,false);
		final LocalEditPatchRegions selector=new LocalEditPatchRegions(k,table,1,false,1,8*k);
		final LocalEditPatchCorrector executor=new LocalEditPatchCorrector(k,table,1,false,8*k);
		check(selector.select(input,1).size==0 && selector.initialRejected,"Automatic selection must retain the whole-read initial burden guard.");
		final IntList selected=selector.select(input,4);
		check(selected.size==8 && selector.regions.merged==2,"Four verified original proposals should form two nearby-error regions without hand-supplied coordinates.");
		check(input.bases==original && input.quality==quality && Arrays.equals(original,bcopy) && Arrays.equals(quality,qcopy),"Selection must leave originals and their aliases intact.");
		check(new LocalEditEngine(k,table::count,1).correct(expected,8)==4,"The native sequential control must recover four fixture errors.");
		check(executor.correct(input,selected,4)==4 && executor.stagedPatches==2,"Automatically selected regions must recover all four edits with one final transaction.");
		check(Arrays.equals(input.bases,expected.bases) && Arrays.equals(input.quality,expected.quality),"Automatic local correction must preserve complete bases/qualities on both strands and with null qualities.");
		check(Arrays.equals(input.bases,reverse ? AminoAcid.reverseComplementBases(truth) : truth),"Automatically selected fixture corrections must equal known truth, not merely agree with another path.");
		check(Arrays.equals(original,bcopy) && Arrays.equals(quality,qcopy),"One final materialization cannot mutate original caller aliases.");
		final long repeatedProfile=selector.profileQueries+executor.profileQueries;
		final Read once=new Read(bcopy.clone(),qcopy==null ? null : qcopy.clone(),input.id,0,false);
		check(executor.correctSelected(once,selector,4)==4 && Arrays.equals(once.bases,input.bases) && Arrays.equals(once.quality,input.quality),"Shared preflight must preserve every output byte of the repeated-preflight control.");
		check(executor.profileQueries<repeatedProfile && executor.profileQueries>=selector.profileQueries,"Shared counters must include selection exactly once and eliminate the redundant original scan.");
		for(int limit:new int[]{1,2,3}){
			final Read rejected=new Read(bcopy.clone(),qcopy==null ? null : qcopy.clone(),input.id,0,false);final byte[] rb=rejected.bases,rq=rejected.quality;
			check(executor.correctSelected(rejected,selector,limit)==0 && (limit==1 ? executor.initialRejected : executor.budgetRejected),"Sharing preflight cannot remove initial-limit rejection or later whole-read budget rollback.");
			check(rejected.bases==rb && rejected.quality==rq && executor.materializedBases==0 && executor.substitutions+executor.insertions+executor.deletions==0,"Shared-preflight rejection must retain complete original arrays and publish no edits.");
		}
		expect(()->executor.correctSelected(once,new LocalEditPatchRegions(k,table,1,false,8,8*k),4),"A sparse selector cannot silently replace the dense original preflight contract.");
		expect(()->executor.correctSelected(once,new LocalEditPatchRegions(k,new Oracle(k,truth),1,false,1,8*k),4),"Equal-looking lookup contents do not prove identity of the original immutable evidence pool.");
		expect(()->executor.correctSelected(once,new LocalEditPatchRegions(k,table,3,false,1,8*k),4),"Different support windows cannot share preflight implicitly.");
		expect(()->executor.correctSelected(once,new LocalEditPatchRegions(k,table,1,true,1,8*k),4),"Different competition rules cannot share preflight implicitly.");
		expect(()->executor.correctSelected(once,new LocalEditPatchRegions(k,table,1,false,1,9*k),4),"Different region bounds must be explicit rather than inherited from stale selection.");
		final Read failed=new Read(bcopy.clone(),qcopy==null ? null : qcopy.clone(),input.id,0,false);final byte[] fb=failed.bases;
		table.fail=true;boolean threw=false;try{executor.correctSelected(failed,selector,4);}catch(IllegalStateException expectedFailure){threw=true;}
		check(threw && failed.bases==fb && selector.regions.finish().size==0 && executor.materializedBases==0,"A selection lookup failure must leave no usable stale regions or partially installed read.");
		table.fail=false;check(executor.correctSelected(failed,selector,4)==4,"Shared selector/executor workers remain reusable after lookup failure and guard rejection.");
		check(selector.select(input,4).size==0,"Worker reuse on the corrected read must clear all prior regions.");
	}
	private static byte other(final byte b){return b=='A' ? (byte)'C' : (byte)'A';}
	private static int[] values(final IntList x){return Arrays.copyOf(x.array,x.size);}
	private static void reverse(final byte[] q){for(int a=0,b=q.length-1;a<b;a++,b--){final byte x=q[a];q[a]=q[b];q[b]=x;}}
	private static String key(final byte[] bases,final int start,final int k){final Kmer word=new Kmer(k);for(int i=start;i<start+k;i++){word.addRight(bases[i]);}return Arrays.toString(word.key());}
	private static final class Oracle implements HomopolymerIndelProposal.CountLookup {
		Oracle(final int k,final byte[] truth){for(int i=0;i<=truth.length-k;i++){counts.put(key(truth,i,k),12);}}
		@Override public int count(final Kmer word){if(fail){throw new IllegalStateException("Injected original-selection lookup failure.");}final Integer value=counts.get(Arrays.toString(word.key()));return value==null ? 0 : value;}
		final HashMap<String,Integer> counts=new HashMap<String,Integer>();
		boolean fail;
	}
	private static void expect(final Runnable action,final String why){boolean failed=false;try{action.run();}catch(IllegalArgumentException expected){failed=true;}check(failed,why);}
	private static void check(final boolean value,final String why){checks++;if(!value){throw new AssertionError(why);}}
	private static int checks;
	private static final byte[] ALPHABET={'A','C','G','T'};
}
