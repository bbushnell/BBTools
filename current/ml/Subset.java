package ml;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import shared.Shared;
import shared.Tools;
import structures.FloatList;

public class Subset {

	public Subset(ArrayList<Sample> sampleList) {
		final Sample[] x=new Sample[0];
		final ArrayList<Sample> pos=new ArrayList<Sample>();
		final ArrayList<Sample> neg=new ArrayList<Sample>();
		samples=sampleList.toArray(x);
		for(Sample s : samples) {
			
			if(s.goal[0]>=0.5f) {
				pos.add(s);
			}else{
				neg.add(s);
			}
		}
		positive=pos.toArray(x);
		negative=neg.toArray(x);
	}
	
	//Sorts the samples by error magnitude, then interleaves positive and negative errors
	void sortSamples(float fraction, boolean allowMultithreadedSort) {
		fraction=(Tools.min(fraction, 1f));
//		for(Sample s : samples) {
//			assert(s.calcPivot()==s.pivot); //Should be done in subthread
//			s.setPivot();
//		}
		if(allowMultithreadedSort) {
			Shared.sort(positive, 0, (int)Math.ceil(fraction*positive.length));
			Shared.sort(negative, 0, (int)Math.ceil(fraction*negative.length));
		}else {
			Arrays.sort(positive, 0, (int)Math.ceil(fraction*positive.length));
			Arrays.sort(negative, 0, (int)Math.ceil(fraction*negative.length));
		}
		assert(positive.length<2 || positive[0].pivot>=positive[1].pivot) : positive[0].pivot+", "+positive[1].pivot;
//		{
//			final int plim=(int)(fraction*positive.length);
//			final int nlim=(int)(fraction*negative.length);
//			for(int i=0; i<plim; i++) {positive[i].pivot=positive[i].calcPivot();}
//			for(int i=0; i<nlim; i++) {negative[i].pivot=negative[i].calcPivot();}
//			Shared.sort(positive, 0, plim-1);
//			Shared.sort(negative, 0, nlim-1);
//			assert(plim<2 || positive[0].pivot>=positive[1].pivot) : positive[0].pivot+", "+positive[1].pivot;
//		}
		
		int apos=0, ppos=0, npos=0;
		while(apos<samples.length) {
			if(npos<negative.length) {
				samples[apos]=negative[npos];
				apos++;
				npos++;
			}
			if(ppos<positive.length) {
				samples[apos]=positive[ppos];
				apos++;
				ppos++;
			}
		}
		assert(apos==samples.length);
		assert(ppos==positive.length);
		assert(npos==negative.length);
	}
	
	//Sorts the samples by error magnitude, then interleaves positive and negative errors
	void sortSamples2(float sortFraction, int useSamples, boolean allowMultithreadedSort, FloatList pivots) {
		sortFraction=(Tools.min(sortFraction, 1f));
//		for(Sample s : samples) {
//			assert(s.calcPivot()==s.pivot); //Should be done in subthread
//			s.setPivot();
//		}

//		System.err.println(sortFraction+", "+useSamples);
		final int pSortLimit=(int)Math.ceil(sortFraction*positive.length);
		final int nSortLimit=(int)Math.ceil(sortFraction*negative.length);
		final int pSwapLimit=Tools.min((useSamples+1)/2, positive.length);
		final int nSwapLimit=Tools.min((useSamples+1)/2, negative.length);
		

//		assert(pSwapLimit<=pSortLimit && pSwapLimit>0) : 
//			"\n"+pSwapLimit+", "+pSortLimit+", "+positive.length+
//			"\n"+nSwapLimit+", "+nSortLimit+", "+negative.length+
//			"\n"+useSamples+", "+sortFraction;
//		assert(nSwapLimit<=nSortLimit && nSwapLimit>0) : 
//			"\n"+pSwapLimit+", "+pSortLimit+", "+positive.length+
//			"\n"+nSwapLimit+", "+nSortLimit+", "+negative.length+
//			"\n"+useSamples+", "+sortFraction;
		
		//****TODO: fraction should be different for positive and negative.
		//Basically, the same number from each should get sorted.
		if(Trainer.PIVOT_SORT) {
			pivotSort(positive, pivots, pSortLimit, pSwapLimit, allowMultithreadedSort);
			pivotSort(negative, pivots, nSortLimit, nSwapLimit, allowMultithreadedSort);
		}else {
			sort(positive, pivots, pSortLimit, allowMultithreadedSort);
			sort(negative, pivots, nSortLimit, allowMultithreadedSort);
		}
		
		int apos=0, ppos=0, npos=0;
		while(apos<samples.length) {//Honestly, they don't all need to be inserted most of the time...
			if(npos<negative.length) {
				samples[apos]=negative[npos];
				apos++;
				npos++;
			}
			if(ppos<positive.length) {
				samples[apos]=positive[ppos];
				apos++;
				ppos++;
			}
		}
		assert(apos==samples.length);
		assert(ppos==positive.length);
		assert(npos==negative.length);
	}
	
	private static final void sort(final Sample[] samples, final FloatList pivots, 
			final int lim, final boolean mt) {
		if(mt) {
			Shared.sort(samples, 0, lim);
		}else {
			Arrays.sort(samples, 0, lim);
		}
		assert(samples.length<2 || samples[0].pivot>=samples[1].pivot) : samples[0].pivot+", "+samples[1].pivot;
	}
	
	private static final void pivotSort(final Sample[] samples, final FloatList pivots, 
			int sortLim, final int swapLim, final boolean mt) {
		assert(sortLim>0 && sortLim<=samples.length);
		assert(swapLim>0 && swapLim<=samples.length);
		sortLim=Tools.mid(sortLim, swapLim*2, samples.length);
//		if(swapLim>=samples.length) {return;}
		pivots.clear();
		for(int i=0; i<sortLim; i++) {
//			assert(s.calcPivot()==s.pivot);
//			s.setPivot();//Should be done in subthread
			pivots.add(samples[i].pivot);
		}
		if(mt) {
			Shared.sort(pivots.array, 0, sortLim);
		}else {
			Arrays.sort(pivots.array, 0, sortLim);
		}
//		assert(swapLim<=sortLim && swapLim>0) : swapLim+", "+sortLim+", "+samples.length;
		//Here, we ensure the pivot is within the sorted range; otherwise there was no point in sorting
		
		final int cutoffLoc=Tools.max(sortLim-swapLim, 0);
		final float cutoff=pivots.get(cutoffLoc);
		for(int src=0, dst=0; src<samples.length && dst<swapLim; src++) {
//			System.err.println(src+", "+dst+", "+swapLim);
			final float p=samples[src].pivot;
			if(p>=cutoff) {
				Sample temp=samples[dst];
				samples[dst]=samples[src];
				samples[src]=temp;
				dst++;
			}
		}
//		assert(samples.length<2 || lim>=samples.length-1 || samples[0].pivot>=samples[lim].pivot) : 
//			samples[0].pivot+", "+samples[lim].pivot+", "+samples[samples.length-1].pivot;
	}
	
	void triage(long currentEpoch, long startTriage, float positiveTriage, float negativeTriage) {
		assert(currentEpoch>=startTriage);
		if(positiveTriage<=0 && negativeTriage<=0) {return;}
		final int distance=500;
		{
			final int max=Math.round(positiveTriage*positive.length);
			if(max>0){
				for(Sample s : positive) {s.setPivot();;}
				Shared.sort(positive);
//				Tools.reverseInPlace(positive);
				for(int i=0; i<max; i++) {positive[i].setEpoch(currentEpoch+distance);}//send it to the future
				for(int i=Tools.max(max, positive.length-max); i<positive.length; i++) {
//					assert(positive[i].epoch>=currentEpoch) : positive[i].epoch+", "+currentEpoch+", "+i;
					positive[i].setEpoch(currentEpoch);//Not necessary, but resets the old triage victims.
				}
			}
		}
		{
			final int max=Math.round(negativeTriage*negative.length);
			if(max>0){
				for(Sample s : negative) {s.setPivot();}
				Shared.sort(negative);
//				Tools.reverseInPlace(negative);
				for(int i=0; i<max; i++) {negative[i].setEpoch(currentEpoch+distance);}//send it to the future
				for(int i=Tools.max(max, negative.length-max); i<negative.length; i++) {
//					assert(negative[i].epoch>=currentEpoch) : negative[i].epoch+", "+currentEpoch+", "+i;
					negative[i].setEpoch(currentEpoch);//Not necessary, but resets the old triage victims.
				}
			}
		}
	}
	
	void shuffle() {
//		Random randy=new Random(numShuffles);
		for(int i=0; i<samples.length; i++) {
			int idx=randy.nextInt(samples.length);
			Sample s=samples[idx];
			samples[idx]=samples[i];
			samples[i]=s;
		}
		numShuffles++;
	}
	public void reset() {
		nextFullPassEpoch=-1;
		numShuffles=0;
		randy.setSeed(0);
	}
	
	final Sample[] samples;
	final Sample[] positive;
	final Sample[] negative;
	private final Random randy=new Random(0);
	
	int nextFullPassEpoch=-1;
	int numShuffles=0;
	
}
