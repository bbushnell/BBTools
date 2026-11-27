package prok;

import java.util.Arrays;

/** Two native choices per cobalamin member before the existing DDL order.
 * One worker owns the reusable scratch. Admission deliberately includes minor
 * members' one-hit candidates, which the shared two-hit filter can exclude.
 * @author Ganyu
 */
public final class NcrnaMemberNativePrefix{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaMemberNativePrefix(final int[] sizes){
		if(sizes==null || sizes.length!=3 || sizes[0]<2 || sizes[1]<2 || sizes[2]<2){
			throw new IllegalArgumentException("The native prefix needs three nontrivial member model axes");
		}
		ends=new int[]{0, sizes[0], Math.addExact(sizes[0], sizes[1]), Math.addExact(Math.addExact(sizes[0], sizes[1]), sizes[2])};
		seen=new boolean[ends[3]];
	}

	/*--------------------------------------------------------------*/
	/*----------------         Native Choices       ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts come from the existing query's native index, including models below
	 * its shared cutoff. Ties retain ascending original model identity. */
	public int select(final TrnaKmerIndex index){
		if(index==null){throw new IllegalArgumentException("Native prefix requires measured query counts");}
		Arrays.fill(choices, -1);
		for(int member=0; member<3; member++){
			final int floor=member==0 ? 2 : 1;
			int first=-1, second=-1, fallback=-1;
			for(int model=ends[member]; model<ends[member+1]; model++){
				final int count=index.lastSharedCount(model);
				if(count<=0){continue;}
				if(fallback<0 || count>index.lastSharedCount(fallback)){fallback=model;}
				if(count<floor){continue;}
				if(first<0 || count>index.lastSharedCount(first)){
					second=first;
					first=model;
				}else if(second<0 || count>index.lastSharedCount(second)){second=model;}
			}
			choices[member]=first>=0 ? first : fallback;
			choices[3+member]=second;
		}
		return compact();
	}

	/** Saved member-native orders provide an independent way to replay the same
	 * prefix on exact matching window geometries, without new sequence queries. */
	public int selectSaved(final int[][] memberOrders){
		if(memberOrders==null || memberOrders.length!=3){throw new IllegalArgumentException("Bind all three exact native member orders");}
		Arrays.fill(choices, -1);
		for(int member=0; member<3; member++){
			final int[] order=memberOrders[member];
			if(order==null){throw new IllegalArgumentException("An absent member window is unmeasured, not an empty shortlist");}
			for(int position=0; position<Math.min(2, order.length); position++){
				final int model=order[position];
				if(model<ends[member] || model>=ends[member+1]){throw new IllegalArgumentException("Native member order must use the complete global axis");}
				choices[position*3+member]=model;
			}
		}
		return compact();
	}

	private int compact(){
		size=0;
		for(final int model:choices){if(model>=0){prefix[size++]=model;}}
		return size;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Ordering and Admission  ----------------*/
	/*--------------------------------------------------------------*/

	/** Prefix then original order, preserving first occurrence and the total cap. */
	public int merge(final int[] original, final int count, final int[] destination, final int cap){
		if(original==null || destination==null || original==destination || count<0 || count>original.length
			|| cap<1 || destination.length<cap){throw new IllegalArgumentException("Prefix merge needs distinct bounded source/destination arrays");}
		Arrays.fill(seen, false);
		int written=0;
		for(int i=0; i<size && written<cap; i++){
			final int model=prefix[i];
			assert(model>=0 && model<seen.length && !seen[model]) : "Disjoint native members supply distinct original model IDs";
			destination[written++]=model;
			seen[model]=true;
		}
		for(int i=0; i<count && written<cap; i++){
			final int model=original[i];
			if(model<0 || model>=seen.length){throw new IllegalArgumentException("Original DDL order must use the bound model axis");}
			if(!seen[model]){
				destination[written++]=model;
				seen[model]=true;
			}
		}
		return written;
	}

	/** Extend the observed native set by the explicit member-prefix admissions.
	 * The caller supplies the original shared eligible set in this array. */
	public int extendEligibility(final int[] destination, final int count){
		if(destination==null || destination.length<seen.length || count<0 || count>seen.length){
			throw new IllegalArgumentException("Eligibility scratch must cover the full model axis");
		}
		Arrays.fill(seen, false);
		for(int i=0; i<count; i++){seen[destination[i]]=true;}
		int total=count;
		for(int i=0; i<size; i++){
			final int model=prefix[i];
			if(!seen[model]){
				destination[total++]=model;
				seen[model]=true;
			}
		}
		Arrays.sort(destination, 0, total);
		return total;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Accessors and State    ----------------*/
	/*--------------------------------------------------------------*/

	public int size(){return size;}
	public int modelAt(final int position){
		if(position<0 || position>=size){throw new IllegalArgumentException("Prefix position is outside its live entries");}
		return prefix[position];
	}
	private final int[] ends;
	private final int[] choices=new int[6], prefix=new int[6];
	private final boolean[] seen;
	private int size;
}
