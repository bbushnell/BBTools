package prot;

/**
 * Integer allocation for the D39 starting recipe in METHODS_v1 Part II.C.
 * Each group stores breakpoint, HQ, contaminated-intact and ordinary counts.
 * Global spike margins survive allocation into the 15/15/15/55 phylum groups.
 * Rounding ties favor the earlier group; ordinary rows absorb the final remainder.
 *
 * @author Yoimiya
 */
final class MagQCD39Quotas {

	private MagQCD39Quotas(){}

	/** Group names also identify the intended quota in generated bin IDs. */
	static final String[] ROLES={"main","foreign","both","remaining"};

	/**
	 * Returns four corpus counts, or sixteen counts grouped by native/foreign phylum role.
	 * Near-perfect rows use only main and remaining groups. Contaminated-intact counts
	 * are apportioned from each group's ordinary balance by largest remainder.
	 */
	static int[] allocate(int count,boolean phylum){
		if(count<0){throw new IllegalArgumentException("Negative D39 row count: "+count);}
		final int breakpoint=percent(count,15),hq=percent(count,5),isolate=percent(count,10);
		return allocate(count,phylum,new int[]{breakpoint,hq,isolate,count-breakpoint-hq-isolate});
	}

	/** Allocates externally balanced category margins while preserving this model's role margins. */
	static int[] allocate(int count,boolean phylum,int[] categories){
		if(count<0 || categories==null || categories.length!=4){throw new IllegalArgumentException("Four nonnegative D39 category counts required");}
		long totalCategories=0;
		for(int value : categories){if(value<0){throw new IllegalArgumentException("Negative D39 category quota");} totalCategories+=value;}
		if(totalCategories!=count){throw new IllegalArgumentException("D39 categories do not sum to requested rows");}
		if(!phylum){return categories.clone();}
		final int breakpoint=categories[0],hq=categories[1],isolate=categories[2];
		final int part=percent(count,15);
		final int[] margins={part,part,part,count-3*part};
		final int[] result=new int[16];
		final int eligible=margins[0]+margins[3];
		if(breakpoint+(long)hq>eligible){throw new IllegalArgumentException("Near-perfect quotas exceed eligible phylum roles");}
		final int mainBreak=eligible==0 ? 0 : Math.max(Math.max(0,breakpoint-margins[3]),
			Math.min(part,(int)((breakpoint*(long)part+eligible/2)/eligible)));
		final int remainingBreak=breakpoint-mainBreak;
		final int mainHq=eligible==0 ? 0 : Math.max(Math.max(0,hq-(margins[3]-remainingBreak)),
			Math.min(part-mainBreak,(int)((hq*(long)part+eligible/2)/eligible)));
		result[0]=mainBreak; result[1]=mainHq;
		result[12]=breakpoint-mainBreak; result[13]=hq-mainHq;
		final int[] ordinary=new int[4];
		int available=0;
		for(int role=0; role<4; role++){
			ordinary[role]=margins[role]-result[role*4]-result[role*4+1];
			if(ordinary[role]<0){throw new IllegalStateException("D39 rounded spikes exceed group size for n="+count);}
			available+=ordinary[role];
		}
		final long[] remainder=new long[4];
		int assigned=0;
		for(int role=0; role<4; role++){
			final long scaled=isolate*(long)ordinary[role];
			result[role*4+2]=available==0 ? 0 : (int)(scaled/available);
			remainder[role]=available==0 ? 0 : scaled%available;
			assigned+=result[role*4+2];
		}
		while(assigned<isolate){
			int best=-1;
			for(int role=0; role<4; role++){
				if(result[role*4+2]<ordinary[role] && (best<0 || remainder[role]>remainder[best])){best=role;}
			}
			if(best<0){throw new IllegalStateException("D39 isolate quota exceeds available mixture rows");}
			result[best*4+2]++; remainder[best]=-1; assigned++;
		}
		int total=0;
		for(int role=0; role<4; role++){
			result[role*4+3]=ordinary[role]-result[role*4+2];
			for(int category=0; category<4; category++){total+=result[role*4+category];}
		}
		assert(total==count) : "Phylum allocation must conserve the requested vector count";
		return result;
	}

	/** Half-up percentage rounding uses long arithmetic for supported int-sized pools. */
	private static int percent(int count,int pct){return (int)((count*(long)pct+50)/100);}
}
