package align2;

import java.util.Arrays;

public class Solver {
	
	
	public static final long bruteForce(int[] offsets, int[] lengths, int chunk, int minLists, int maxTotalLength){
		
		int bits=offsets.length;
		int max=(1<<bits)-1;
		
		for(long i=0; i<=max; i++){
			long x=evaluate(offsets, lengths, chunk, i);
		}
		
		assert(false);
		return 0;
	}
	
	
	public static final void findWorstGreedy(final int[] offsets, final int[] lengths,
			final int chunk, final int[] lists, int[] r){
		assert(r!=null && r.length==2);
		
		long min=Long.MAX_VALUE;
		int worstIndex=-1;
		for(int i=0; i<lists.length; i++){
			long value=valueOfElement(offsets, lengths, 1f, chunk, lists, i);
			if(value<min){
				if(min<EARLY_TERMINATION_SCORE){//Can speed up greedy algo
					r[0]=i;
					r[1]=(value<Integer.MIN_VALUE ? Integer.MIN_VALUE : value>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)value);
					return;
				}
				min=value;
				worstIndex=i;
			}
		}
//		if(min>0){worstIndex=-1;}
		r[0]=worstIndex;
		r[1]=(min<Integer.MIN_VALUE ? Integer.MIN_VALUE : min>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)min);
	}
	
	
	public static final void findWorstGreedy(final int[] offsets, final int[] lengths,
			final float[] weights, final int chunk, final int[] lists, int[] r){
		assert(r!=null && r.length==2);
		
		long min=Long.MAX_VALUE;
		int worstIndex=-1;
		for(int i=0; i<lists.length; i++){
//		for(int i=lists.length-1; i>=0; i--){
			long value=valueOfElement(offsets, lengths, weights[i], chunk, lists, i);
			if(value<min){
				if(min<EARLY_TERMINATION_SCORE && i!=0){//Can speed up greedy algo
					r[0]=i;
					r[1]=(value<Integer.MIN_VALUE ? Integer.MIN_VALUE : value>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)value);
//					System.out.print(".");
					return;
				}
				min=value;
				worstIndex=i;
			}
		}
//		if(min>0){worstIndex=-1;}
		r[0]=worstIndex;
		r[1]=(min<Integer.MIN_VALUE ? Integer.MIN_VALUE : min>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)min);
	}
	
	
	public static long valueOfElement(final int[] offsets, final int[] lengths, float keyWeight,
			final int chunk, final int[] lists, int index){
		
		final int numlists=lists.length;
		if(numlists<1){return 0;}
		
		final int prospect=lists[index];
		if(lengths[prospect]==0){return -999999;}
		
		long valuep=POINTS_PER_LIST+(POINTS_PER_LIST*2/lists.length)+((POINTS_PER_LIST*10)/lengths[prospect]);
		long valuem=POINTS_PER_SITE*lengths[prospect];
		
		if(prospect==0 || (prospect==offsets.length-1)){
			valuep+=BONUS_POINTS_FOR_END_LIST;
		}
		
		if(numlists==1){
			valuep+=(POINTS_FOR_TOTAL_LIST_WIDTH+POINTS_PER_BASE1)*chunk;
			return ((long)(valuep*keyWeight))+valuem;
		}
		
		
		final int first=lists[0];
		final int last=lists[numlists-1];
		
		//Offsets of elements to the left and right of the prospect
//		final int offL=(prospect==first ? - : offsets[lists[index-1]]);
//		final int offP=offsets[prospect];
//		final int offR=(prospect==last ? offsets[offsets.length-1] : offsets[lists[index+1]]);
//		assert(offL<=offP);
//		assert(offP<=offR);
//		assert(offL<offR) : "\noffsets.length="+offsets.length+", lengths.length="+lengths.length+"\n"+
//			", chunk="+chunk+", lists.length="+lists.length+", index="+index+"\n"+
//			", offL="+offL+", offR="+offR+", prospect="+prospect+", first="+first+", last="+last+"\n"+
//			", valuep="+valuep+", valuem="+valuem+", weight="+keyWeight+"\n"+
//			"offsets = "+Arrays.toString(offsets)+"\tlengths = "+Arrays.toString(lengths)+"\nlists = "+Arrays.toString(lists)+"\n";
		
		final int offL=(prospect==first ? -1 : offsets[lists[index-1]]);
		final int offP=offsets[prospect];
		final int offR=(prospect==last ? offsets[offsets.length-1]+1 : offsets[lists[index+1]]);
		assert(offL<=offP);
		assert(offP<=offR);
		assert(offL<offR) : "\noffsets.length="+offsets.length+", lengths.length="+lengths.length+"\n"+
			", chunk="+chunk+", lists.length="+lists.length+", index="+index+"\n"+
			", offL="+offL+", offR="+offR+", prospect="+prospect+", first="+first+", last="+last+"\n"+
			", valuep="+valuep+", valuem="+valuem+", weight="+keyWeight+"\n"+
			"offsets = "+Arrays.toString(offsets)+"\tlengths = "+Arrays.toString(lengths)+"\nlists = "+Arrays.toString(lists)+"\n";

		int oldLeftSpace=offP-offL;
		int oldRightSpace=offR-offP;
		int newSpace=offR-offL;
		
//		int oldLeftSpace=Tools.max((offP-offL)-1, 0)+1;
//		int oldRightSpace=Tools.max((offR-offP)-1, 0)+1;
//		int newSpace=Tools.max(offR-offL;
		
		long spaceScore=((oldLeftSpace*oldLeftSpace+oldRightSpace*oldRightSpace)-(newSpace*newSpace))*MULT_FOR_SPACING_PENALTY;
		assert(spaceScore>0) : "\n"+spaceScore+", "+oldLeftSpace+", "+oldRightSpace+", "+newSpace+"\n"+
			Arrays.toString(offsets)+"\nprospect="+prospect+"\n";
		valuep+=spaceScore;
		
		int uniquelyCovered;
		if(prospect==first){
			uniquelyCovered=offR-offP; //Technically, -1 should be added
		}else if(prospect==last){
			uniquelyCovered=offP-offL; //Technically, -1 should be added
		}else{
			int a=offL+chunk;
			int b=offR-a;
			uniquelyCovered=(b>0 ? b : 0);
		}
		
		if(prospect==first || prospect==last){
			valuep+=(POINTS_PER_BASE1+POINTS_FOR_TOTAL_LIST_WIDTH)*uniquelyCovered;
		}else{
			valuep+=POINTS_PER_BASE1*uniquelyCovered;
		}
		
		return ((long)(valuep*keyWeight))+valuem;
	}
	
	public static int[] toBitList(final int key){
		final int numlists=Integer.bitCount(key);
		final int[] lists=new int[numlists];
		for(int i=0, ptr=0; ptr<numlists; i++){
			if((masks32[i]&key)!=0){
				lists[ptr]=i;
				ptr++;
			}
		}
		return lists;
	}
	
	public static int[] toBitList(final long key){
		final int numlists=Long.bitCount(key);
		assert(numlists>0);
		final int[] lists=new int[numlists];
		for(int i=0, ptr=0; ptr<numlists; i++){
			if((masks[i]&key)!=0){
				lists[ptr]=i;
				ptr++;
			}
		}
		return lists;
	}
	
	public static long evaluate(int[] offsets, int[] lengths, final int chunk, final long key){
		
		long score=0;
		
		final int[] lists=toBitList(key);
		final int numlists=lists.length;
		
		final int first=lists[0];
		final int last=lists[numlists-1];
		
		score+=numlists*POINTS_PER_LIST;
		for(int i=0; i<numlists; i++){
			int list=lists[i];
			score+=POINTS_PER_SITE*lengths[list];
		}
		if(first==0){score+=BONUS_POINTS_FOR_END_LIST;}
		if(last==offsets.length-1){score+=BONUS_POINTS_FOR_END_LIST;}
		
		score+=(POINTS_FOR_TOTAL_LIST_WIDTH*(offsets[last]-offsets[first]+chunk));
		
		//TODO: Special case both ends
		for(int i=1; i<numlists; i++){
			int list1=lists[i-1];
			int list2=lists[i];
			int space=offsets[list2]-offsets[list1];
			int uncovered=space>chunk ? space-chunk : 0;
			
			score+=MULT_FOR_SPACING_PENALTY*(space*space);
			score-=POINTS_PER_BASE1*uncovered;
		}
		
		if(first>0){
			long x=offsets[first];
			score+=MULT_FOR_SPACING_PENALTY*(x*x);
			score-=POINTS_PER_BASE1*x;
		}
		
		if(last<(offsets.length-1)){
			long x=offsets[offsets.length-1]-offsets[last];
			score+=MULT_FOR_SPACING_PENALTY*(x*x);
			score-=POINTS_PER_BASE1*x;
		}
		
		return score;
	}
	
	public static final int BASE_POINTS_PER_SITE=-50; //Used to set POINTS_PER_SITE
	public static long POINTS_PER_SITE=-50; //TODO:  Make private with a get() and set() function
	
	public static final long MULT_FOR_SPACING_PENALTY=-30;
	
	public static long EARLY_TERMINATION_SCORE=(POINTS_PER_SITE*2000); //TODO: Should be set dynamically
	
	public static final long POINTS_PER_LIST=30000;
	public static final long POINTS_PER_BASE1=6000; //Points for a base covered once
	public static final long POINTS_PER_BASE2=1000;//POINTS_PER_BASE1/4; //Points for a base covered twice
	public static final long BONUS_POINTS_FOR_END_LIST=40000; //Extra points for the first and last list
	public static final long POINTS_FOR_TOTAL_LIST_WIDTH=5500; //multiplier for distance between first and last list

	public static final long[] masks=new long[64];
	public static final int[] masks32=new int[32];
	static{
		for(int i=0; i<masks.length; i++){masks[i]=(1L<<i);}
		for(int i=0; i<masks32.length; i++){masks32[i]=(1<<i);}
	}
}
