package align2;

import java.util.Arrays;

import dna.AminoAcid;
import dna.ChromosomeArray;
import shared.Tools;

public final class KeyRing {
	
	public static final void main(String[] args){
		int len=Integer.parseInt(args[0]);
		float density=(float) Double.parseDouble(args[1]);
		int chunksize=13;
		if(args.length>2){chunksize=Integer.parseInt(args[2]);}
		
		byte[] qual=new byte[len];
		Arrays.fill(qual, (byte)20);
		
		int[] offsets=KeyRing.makeOffsets(qual, chunksize, density, 2);
		System.out.println(Arrays.toString(offsets));
	}
	
	public static int[] makeKeys(byte[] s, int[] offsets, int chunksize){
		if(offsets==null){return null;}
		assert(chunksize>0 && chunksize<16);
		assert(offsets!=null) : s.length+", "+new String(s);
		int[] keys=new int[offsets.length];
		
//		System.out.println(Arrays.toString(offsets));
		
		for(int i=0; i<offsets.length; i++){
//			System.out.println(s.length()+", "+offsets.length+", "+chunksize+", "+keys.length+", "+i);
			keys[i]=ChromosomeArray.toNumber(offsets[i], offsets[i]+chunksize-1, s);
		}
		return keys;
	}
	
	public static int[] reverseComplementKeys(int[] keys, int k){
//		assert(!cs);
		int[] r=new int[keys.length];
		for(int i=0, x=keys.length-1; i<r.length; i++){
			r[i]=AminoAcid.reverseComplementBinaryFast(keys[x-i], k);
		}
		return r;
	}
	
	public static int reverseComplementKey(int key, int k){
//		return cs ? reverseComplementKey_old(key, k, cs) : AminoAcid.reverseComplementBinaryFast(key, k);
		return AminoAcid.reverseComplementBinaryFast(key, k);
	}
	
	
	public static final String decode(int key, int chunksize){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<chunksize; i++){
			int temp=(key>>(2*(chunksize-i-1)));
			temp=(temp&3);
			sb.append((char)AminoAcid.numberToBase[temp]);
		}
		
		String s=sb.toString();
		
		assert(key==ChromosomeArray.toNumber(0, s.length()-1, s)) :
			Integer.toHexString(key)+" -> "+s+" != "+Integer.toHexString(ChromosomeArray.toNumber(0, s.length()-1, s));
		
		return sb.toString();
	}

	/*
	public static final int[] makeOffsets(int readlen, int blocksize, int overlap, int minKeysDesired){
		assert(blocksize>0);
		assert(overlap<blocksize);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+overlap+", "+minKeysDesired;
		
		int slots=readlen-blocksize+1;
		int midslots=slots-2;
		int spacing=blocksize-overlap;

		if(slots==1){return new int[] {0};}
		if(slots<=spacing+1){return new int[] {0, slots-1};}
		
//		int middles=(midslots/spacing);
//
//		if(middles+2<minKeysDesired && midslots+2>=minKeysDesired){
//			while(middles+2<minKeysDesired){
//				spacing--;
//				assert(spacing>0);
//				middles=(midslots/spacing);
//			}
//		}
		
		int middles=(midslots/spacing);
		if(middles<minKeysDesired-2){
			middles=Tools.max(minKeysDesired-2, midslots);
		}
		
		assert(middles>0); //due to the escape conditions

//		float fspacing=midslots/(float)(middles+1);
		float fspacing=midslots/(float)(middles);
		assert(fspacing>=1);
		
		int[] offsets=new int[middles+2];
		offsets[0]=0;
		offsets[offsets.length-1]=slots-1;
		
		for(int i=1; i<=middles; i++){
			offsets[i]=Math.round(fspacing*i);
		}

//		System.out.println("readlen = \t"+readlen);
//		System.out.println("blocksize = \t"+blocksize);
//		System.out.println("overlap = \t"+overlap);
//		System.out.println("slots = \t"+slots);
//		System.out.println("midslots = \t"+midslots);
//		System.out.println("spacing = \t"+spacing);
//		System.out.println("middles = \t"+middles);
//		System.out.println("fspacing = \t"+fspacing);
//		System.out.println("Offsets = \t"+Arrays.toString(offsets));
		return offsets;
		
	}*/

	/** This is only useful for low-quality reads, with no-calls.  Otherwise it just wastes time... */
	public static final int[] reverseOffsets(final int[] offsetsP, final int k, final int readlen){
		int[] offsetsM=new int[offsetsP.length];
		for(int i=0; i<offsetsP.length; i++){
			int x=offsetsP[offsetsP.length-i-1];
			assert(x>=0);
			assert(x+k<=readlen);
			x=readlen-(x+k);
			assert(x>=0);
			assert(x+k<=readlen) : "\n"+Arrays.toString(offsetsP)+"\n"+Arrays.toString(offsetsM)+"\n"+i+"\n"+x+"\n"+readlen;
			offsetsM[i]=x;
		}
		return offsetsM;
	}
	
	public static final int[] makeOffsetsWithDensity(int readlen, int blocksize, float density, int minKeysDesired){
		assert(blocksize>0);
		assert(density<blocksize);
		assert(density>0);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+density+", "+minKeysDesired;
		
		int slots=readlen-blocksize+1;
		int midslots=slots-2;
		
		int desired=(int)Math.ceil((readlen*density)/blocksize);
		assert(desired>=0);
		desired=Tools.max(minKeysDesired, desired);
		desired=Tools.min(slots, desired);

		if(slots==1 || desired==1){return new int[] {0};}
		if(desired==2){return new int[] {0, slots-1};}
		
		int middles=desired-2;
		
		assert(middles>0); //due to the escape conditions

//		float fspacing=midslots/(float)(middles+1);
		float fspacing=midslots/(float)(middles);
		assert(fspacing>=1);
		
		int[] offsets=new int[desired];
		offsets[0]=0;
		offsets[offsets.length-1]=slots-1;
		
		for(int i=1; i<=middles; i++){
			offsets[i]=Math.round(fspacing*i);
		}

//		System.out.println("readlen = \t"+readlen);
//		System.out.println("blocksize = \t"+blocksize);
//		System.out.println("overlap = \t"+overlap);
//		System.out.println("slots = \t"+slots);
//		System.out.println("midslots = \t"+midslots);
//		System.out.println("spacing = \t"+spacing);
//		System.out.println("middles = \t"+middles);
//		System.out.println("fspacing = \t"+fspacing);
//		System.out.println("Offsets = \t"+Arrays.toString(offsets));
		return offsets;
		
	}
	
	
	public static final int[] makeOffsetsWithNumberOfKeys(int readlen, int blocksize, int maxKeys){
		assert(maxKeys>0);
//		System.err.println("readlen, blocksize, maxKeys = "+readlen+","+blocksize+","+maxKeys);
		if(blocksize>readlen){return null;}
		int slots=readlen-blocksize+1;
//		System.err.println("slots = "+slots);
		if(slots==1 || maxKeys==1){return new int[] {slots/2};}
		if(slots==2 || maxKeys==2){return new int[] {0, slots-1};}
		if(slots==3 || maxKeys==3){return new int[] {0, slots/2, slots-1};}
		
		int midslots=slots-2;
		maxKeys=Tools.min(maxKeys, slots);
		int middles=Tools.min(maxKeys-2, midslots);
//		System.err.println("midslots = "+midslots);
//		System.err.println("middles = "+middles);
		
		assert(middles>0); //due to the escape conditions
		
//		float fspacing=midslots/(float)(middles+0); //Bad - leaves 2 adjacent keys at the end.
		float fspacing=midslots/(float)(middles+1f);
		fspacing=Tools.max(1f, fspacing);
		assert(fspacing>=1);
		
		int[] offsets=new int[middles+2];
		offsets[0]=0;
		offsets[offsets.length-1]=slots-1;
		

//		for(int i=1; i<=middles; i++){
//			offsets[i]=Math.round(fspacing*i);
//		}
		
		
		
		for(int i=1; i<=middles; i++){
			offsets[i]=Math.round(fspacing*i);
		}
		if(middles>2){
			offsets[1]=(int)fspacing;
			offsets[middles]=(int) Math.ceil(fspacing*middles);
		}

//		System.out.println("readlen = \t"+readlen);
//		System.out.println("blocksize = \t"+blocksize);
////		System.out.println("overlap = \t"+overlap);
//		System.out.println("slots = \t"+slots);
//		System.out.println("midslots = \t"+midslots);
////		System.out.println("spacing = \t"+spacing);
//		System.out.println("middles = \t"+middles);
//		System.out.println("fspacing = \t"+fspacing);
//		System.out.println("Offsets = \t"+Arrays.toString(offsets));
		
		for(int i=1; i<offsets.length; i++){
			if(offsets[i]<=offsets[i-1]){assert(false) : "fspacing "+fspacing+"\nmidslots "+midslots+"\nmiddles "+middles+
				"\nmaxKeys "+maxKeys+"\nslots "+slots+"\noffsets "+Arrays.toString(offsets);}
		}
		
		return offsets;
		
	}
	
//	public static final int desiredKeys(int readlen, int blocksize, int overlap, int minKeysDesired){
//		assert(blocksize>0);
//		assert(overlap<blocksize);
//		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+overlap+", "+minKeysDesired;
//		assert(minKeysDesired>=2);
//
//		int slots=readlen-blocksize+1;
//		int midslots=slots-2;
//		int spacing=blocksize-overlap;
//
//		if(slots<=minKeysDesired){return slots;}
//		if(slots<=spacing+1){return Tools.min(3, slots);}
//
//		int middles=(midslots/spacing);
//		if(middles<minKeysDesired-2){
//			middles=Tools.max(minKeysDesired-2, midslots);
//		}
//
//		assert(middles>0); //due to the escape conditions
//		return middles+2;
//	}
	
	public static final int desiredKeysFromDensity(int readlen, int blocksize, float density, int minKeysDesired){
		assert(blocksize>0);
		assert(density<=blocksize) : density+", "+blocksize;
		assert(density>0);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+density+", "+minKeysDesired;
		
		int slots=readlen-blocksize+1;
		
		int desired=(int)Math.ceil((readlen*density)/blocksize);
		assert(desired>=0);
		desired=Tools.max(minKeysDesired, desired);
		desired=Tools.min(slots, desired);
		return desired;
	}
	
	public static final int[] makeOffsets(final int readlen, int blocksize, float density, int minKeysDesired){
		assert(blocksize>0);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+density+", "+minKeysDesired;
		
		if(readlen<blocksize){return null;}

		int desiredKeys=desiredKeysFromDensity(readlen, blocksize, density, minKeysDesired);
		assert(desiredKeys>0) : readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;
		
		int[] offsets=makeOffsetsWithNumberOfKeys(readlen, blocksize, desiredKeys);
//		System.out.println("desiredKeys="+desiredKeys+", actual="+(offsets==null ? 0 : offsets.length));
		assert(offsets!=null) :readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;
		return offsets;
	}
	
	public static final int[] makeOffsets(byte[] qual, int blocksize, float density, int minKeysDesired){
		int readlen=qual.length;
		assert(blocksize>0);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+density+", "+minKeysDesired;
		
		int left=0, right=readlen-1;

		for(int i=left, cntr=0; i<readlen && cntr<blocksize; i++, cntr++){
			if(qual[i]<1){
				left=i+1;
				cntr=0;
			}
		}
		for(int i=right, cntr=0; i>=0 && cntr<blocksize; i--, cntr++){
			if(qual[i]<1){
				right=i-1;
				cntr=0;
			}
		}
		
//		System.out.println("left="+left+", right="+right+", readlen="+readlen+", " +
//				"blocksize="+blocksize+", density="+density+", minKeysDesired="+minKeysDesired);
		
		readlen=right-left+1;
		assert(readlen<=qual.length);
		if(readlen<blocksize){return null;}

		int desiredKeys=desiredKeysFromDensity(qual.length, blocksize, density, minKeysDesired);
		assert(desiredKeys>0) : qual.length+","+readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;

//		System.out.println("desiredKeys="+desiredKeys);
//		System.out.println("Resulting density = "+(desiredKeys*blocksize)/(float)qual.length);
		
		int[] offsets=makeOffsetsWithNumberOfKeys(readlen, blocksize, desiredKeys);
//		System.out.println("desiredKeys="+desiredKeys+", actual="+(offsets==null ? 0 : offsets.length));
		assert(offsets!=null) : qual.length+","+readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;
		if(left>0){
			for(int i=0; i<offsets.length; i++){offsets[i]+=left;}
		}
		return offsets;
	}
	
//	public static final int[] makeOffsets2(float[] keyErrorProb,
//			final int readlenOriginal, int blocksize, float density, int minKeysDesired){
//		return makeOffsets2(keyErrorProb, readlenOriginal, blocksize, density, 2*density, minKeysDesired);
//	}
	
	public static final int[] makeOffsets2(float[] keyErrorProb,
			final int readlenOriginal, int blocksize, float density, float maxDensity, int minKeysDesired){
		int readlen=readlenOriginal;
		assert(maxDensity>=density);
		assert(blocksize>0);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+density+", "+minKeysDesired;
		
		int left=0, right=readlen-blocksize;
		
		//This can be set as low as .90 for long reads, if qualities are accurate.
		final float errorLimit=KEEP_BAD_KEYS ? 2f : 0.94f; //Default: .95f

		while(left<=right && keyErrorProb[left]>errorLimit){left++;}
		while(right>=left && keyErrorProb[right]>errorLimit){right--;}
		
//		System.out.println("left="+left+", right="+right+", readlen="+readlen+", " +
//				"blocksize="+blocksize+", density="+density+", minKeysDesired="+minKeysDesired);
		
		if(right<left){return null;}
		readlen=right-left+blocksize;
		assert(readlen<=readlenOriginal);
		if(readlen<blocksize){
			assert(false);
			return null;
		}
		
//		System.out.println("Left="+left+", right="+right);

		int desiredKeys=desiredKeysFromDensity(readlenOriginal, blocksize, density, minKeysDesired);
		if(readlen<readlenOriginal){
			int desiredKeys2=desiredKeysFromDensity(readlen, blocksize, maxDensity, minKeysDesired);
			desiredKeys=Tools.min(desiredKeys, desiredKeys2);
		}
		assert(desiredKeys>0) : readlenOriginal+","+readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;

//		System.out.println("desiredKeys="+desiredKeys);
//		System.out.println("Resulting density = "+(desiredKeys*blocksize)/(float)qual.length);
		
		int[] offsets=makeOffsetsWithNumberOfKeys(readlen, blocksize, desiredKeys);
		
//		System.out.println("offsets initial = "+Arrays.toString(offsets));
		
//		System.out.println("desiredKeys="+desiredKeys+", actual="+(offsets==null ? 0 : offsets.length));
		assert(offsets!=null) : readlenOriginal+","+readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;
		if(left>0){
			for(int i=0; i<offsets.length; i++){offsets[i]+=left;}
		}
		return offsets;
	}
	
	public static final int[] makeOffsets3(float[] keyErrorProb,
			final int readlenOriginal, int blocksize, float density, float maxDensity, int minKeysDesired, boolean semiperfectmode){
		int readlen=readlenOriginal;
		assert(maxDensity>=density);
		assert(blocksize>0);
		assert(blocksize<=readlen) : readlen+", "+blocksize+", "+density+", "+minKeysDesired;
		
		final int maxProbIndex=readlen-blocksize;
//		assert(maxProbIndex==keyErrorProb.length-1);
		assert(maxProbIndex<=keyErrorProb.length-1) : maxProbIndex+", "+keyErrorProb.length;
		int left=0, right=maxProbIndex;
		
		final float errorLimit2=KEEP_BAD_KEYS ? 2f : 0.9999f; //Default: .95f
		
		//This can be set as low as .90 for long reads, if qualities are accurate.
		final float errorLimit1=KEEP_BAD_KEYS ? 2f : (semiperfectmode ? 0.99f : 0.94f); //Default: .95f

		while(left<=right && keyErrorProb[left]>=errorLimit1){left++;}
		while(right>=left && keyErrorProb[right]>=errorLimit1){right--;}

//		System.out.println("Left="+left+", right="+right);
		
		int potentialKeys=0;
		for(int i=left; i<=right; i++){
			if(keyErrorProb[i]<errorLimit2){potentialKeys++;}
		}
		if(potentialKeys==0){return null;}
		
//		System.out.println("left="+left+", right="+right+", readlen="+readlen+", " +
//				"blocksize="+blocksize+", density="+density+", minKeysDesired="+minKeysDesired);
		
		if(right<left){return null;}
		readlen=right-left+blocksize;
		assert(readlen<=readlenOriginal);
		if(readlen<blocksize){
			assert(false);
			return null;
		}
		
		int desiredKeys=desiredKeysFromDensity(readlenOriginal, blocksize, density, minKeysDesired);
		if(readlen<readlenOriginal){
			int desiredKeys2=desiredKeysFromDensity(readlen, blocksize, maxDensity, minKeysDesired);
			desiredKeys=Tools.min(desiredKeys, desiredKeys2);
		}
		desiredKeys=Tools.min(desiredKeys, potentialKeys);
		assert(desiredKeys>0) : readlenOriginal+","+readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;

//		System.out.println("desiredKeys="+desiredKeys);
//		System.out.println("Resulting density = "+(desiredKeys*blocksize)/(float)qual.length);
		
		int[] offsets=new int[desiredKeys];
		float interval=(right-left)/(float)(Tools.max(desiredKeys-1, 1));
		
		int intervalInt=((int)interval)+1;
		
		float f=left;
		int prev=-1;
		int misses=0;
		for(int i=0, j=left; i<offsets.length; i++){
			int x=-1;
			
//			System.out.println("prev="+prev+", j="+j+", intervalInt="+intervalInt);
			
			if(prev<j){
				if(keyErrorProb[j]<errorLimit2 && (prev<0 || j-prev>0)){
					x=j;
//					System.out.println("A: x="+x);
				}else{
					for(int k=j-1, lim=prev+2; k>lim; k--){
						if(keyErrorProb[k]<errorLimit2){x=k;break;}
					}
//					System.out.println("B: x="+x);
					if(x<0){
						for(int k=j+1, lim=Tools.min(j+intervalInt, right); k<lim; k++){
							if(keyErrorProb[k]<errorLimit2){x=k;break;}
						}
					}
//					System.out.println("C: x="+x);
				}
			}
			
			offsets[i]=x;
			if(x>-1){
				assert(keyErrorProb[x]<errorLimit2);
				prev=x;
			}else{
				misses++;
				prev=Tools.max(prev, j-2);
			}
			
			f+=interval;
			j=Tools.min(maxProbIndex, (Tools.max(j+1, (int)Math.round(f))));
		}
//		System.out.println("offsets initial = "+Arrays.toString(offsets));
		
		if(misses>0){
			int[] offsets2=new int[offsets.length-misses];
			for(int i=0, j=0; i<offsets.length; i++){
				if(offsets[i]>=0){
					offsets2[j]=offsets[i];
					j++;
				}
			}
			offsets=offsets2;
		}
//		System.out.println("offsets shrunk = "+Arrays.toString(offsets));
		
//		System.out.println("desiredKeys="+desiredKeys+", actual="+(offsets==null ? 0 : offsets.length));
		assert(offsets!=null) : readlenOriginal+","+readlen+","+blocksize+","+density+","+minKeysDesired+","+desiredKeys;
		return offsets;
	}
	
	public static boolean KEEP_BAD_KEYS=false;
	
}
