package dna;

public abstract class Motif {
	
	public Motif(String name_, int length_, int center_){
		center=center_;
		length=length_;
		suffix=length-center-1;
		name=name_;
		
//		assert(center>=0 && center<length);
	}
	
	
	public final int countExact(String s){
		return countExact(s, 0, s.length());
	}
	
	
	public final int countExact(String s, int a, int b){
		return countExact(s.getBytes(), a, b);
	}
	
	
	public final int countExtended(String s){
		return countExtended(s, 0, s.length());
	}
	
	
	public final int countExtended(String s, int a, int b){
		return countExtended(s.getBytes(), a, b);
	}
	
	
	public final int countExact(byte[] source, int a, int b){
		
		int max=min(b, source.length-1)-length+1;
		
		int count=0;
		
		for(int i=a; i<=max; i++){
			if(matchesExactly(source, i)){count++;}
		}
		
		return count;
		
	}
	
	
	public final int countExtended(byte[] source, int a, int b){
		
		int max=min(b, source.length-1)-length+1;
		
		int count=0;
		
		for(int i=a; i<=max; i++){
			if(matchesExtended(source, i)){count++;}
		}
		
		return count;
		
	}
	
	
	public boolean matchesExactly(byte[] source, int a){
		throw new RuntimeException();
	}
	
	
	public boolean matchesExtended(byte[] source, int a){
		throw new RuntimeException();
	}
	
	public float normalize(double strength){
		return (float)strength;
	}
	
	
	public float matchStrength(byte[] source, int a){
		return(matchesExactly(source, a) ? 1 : 0);
	}
	
	
	public static final int minPos(float[] array){
		int pos=0;
		for(int i=1; i<array.length; i++){
			if(array[i]<array[pos]){pos=i;}
		}
		return pos;
	}
	
	
	public static final int maxPos(float[] array){
		int pos=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>array[pos]){pos=i;}
		}
		return pos;
	}
	
	@Override
	public String toString(){
		return name+", "+length+", "+center;
	}
	

	public final String name;
	public String commonLetters;
	public final int center;
	public final int length;
	public final int suffix;
	

	static final int min(int x, int y){return x<y ? x : y;}
	static final int max(int x, int y){return x>y ? x : y;}
	static final float min(float x, float y){return x<y ? x : y;}
	static final float max(float x, float y){return x>y ? x : y;}
	
	static final byte[] numberToBase=AminoAcid.numberToBase;
	static final byte[] numberToBaseExtended=AminoAcid.numberToBaseExtended;
	static final byte[] baseToNumber=AminoAcid.baseToNumberACGTN;
	static final byte[] baseToNumberExtended=AminoAcid.baseToNumberExtended;
	
	static final float[] baseProb1={0.256614f, 0.226617f, 0.238012f, 0.278756f};
	
	//Within 200 of exon and gene ends only
	static final float[] baseProb2={
		0.076019f, 0.046405f, 0.071754f, 0.062437f, 0.067143f, 0.066057f, 0.020333f, 0.073085f,
		0.060553f, 0.054897f, 0.068741f, 0.053822f, 0.052896f, 0.059260f, 0.077188f, 0.089412f
	};
	
	//name: Overall Frequency MP3
	static final float[] baseProb3={
		0.027343f, 0.011857f, 0.018295f, 0.018524f, 0.015942f, 0.012337f, 0.003792f, 0.014333f,
		0.019988f, 0.015837f, 0.020411f, 0.015518f, 0.014382f, 0.011355f, 0.016466f, 0.020234f,
		0.014364f, 0.014299f, 0.022875f, 0.015605f, 0.018893f, 0.019412f, 0.006677f, 0.021076f,
		0.003629f, 0.005854f, 0.006783f, 0.004067f, 0.010491f, 0.018413f, 0.024257f, 0.019924f,
		0.018029f, 0.010640f, 0.019427f, 0.012458f, 0.015158f, 0.017025f, 0.006167f, 0.016547f,
		0.018098f, 0.016891f, 0.020042f, 0.013710f, 0.010580f, 0.010773f, 0.018026f, 0.014443f,
		0.016281f, 0.009609f, 0.011157f, 0.015849f, 0.017150f, 0.017284f, 0.003696f, 0.021130f,
		0.018839f, 0.016316f, 0.021506f, 0.020527f, 0.017442f, 0.018720f, 0.018440f, 0.034811f
	};
	
//	protected static final Hashtable<String, float[]> percentTable=makePercentTable();
//
//	private static final Hashtable<String, float[]> makePercentTable(){
//
//		String[] keys={
//				"Exon Stops MP3",
//		};
//
//		float[][] values={
//				{
//					0.00234f, 0.01071f, 0.02476f, 0.05155f, 0.08682f, 0.1453f, 0.22434f, 0.29615f, 0.36233f, 0.41034f,
//					0.46028f, 0.52224f, 0.58198f, 0.63879f, 0.68356f, 0.70622f, 0.7268f, 0.75131f, 0.77065f, 0.79546f,
//					0.82445f, 0.85279f, 0.86899f, 0.88287f, 0.89197f, 0.90166f, 0.91405f, 0.93129f, 0.94708f, 0.95521f,
//					0.96106f, 0.96293f, 0.9663f, 0.97242f, 0.97662f, 0.97866f, 0.98017f, 0.98242f, 0.98459f, 0.98703f,
//					0.98957f, 0.99064f, 0.99157f, 0.99286f, 0.9952f, 0.99721f, 0.99858f, 0.99914f, 0.99967f, 0.9999f, 0.99998f
//				},
//		};
//
//		Hashtable<String, float[]> r= new Hashtable<String, float[]>();
//		for(int i=0; i<keys.length; i++){
//			r.put(keys[i], values[i]);
//		}
//
//		return r;
//	}

	static final float[] invBaseProb1=invert(baseProb1);
	
	static final float[] invBaseProb2=invert(baseProb2);
	
	static final float[] invBaseProb3=invert(baseProb3);
	
	static final float[][] baseProbN={
		null,
		baseProb1,
		baseProb2,
		baseProb3
	};
	
	static final float[][] invBaseProbN={
		null,
		invBaseProb1,
		invBaseProb2,
		invBaseProb3
	};
	
	private static final float[] invert(float[] in){
		float[] out=new float[in.length];
		for(int i=0; i<in.length; i++){
			out[i]=1f/in[i];
		}
		return out;
	}
	
	protected float[] percentile;

	public abstract int numBases();
	
	public float percentile(float strength){
//		float[] array=percentiles[numBases()];
		
		if(percentile==null){
			throw new RuntimeException("Can't find percentile array for "+this);
		}
		
		float[] array=percentile;
		
		int index=(int)(strength*array.length);
		
//		System.out.print(" *** index = "+index+" -> "+array[index]+" -> "+array[index+1]+" *** ");
		
		if(index>=array.length-1){return 1;}
		
		float a, b;
		if(index==0){
			a=0;
			b=array[0];
		}else{
			a=array[index];
			b=array[index+1];
		}
		
		float ratio=strength-(index/((float)array.length));
		
		return ratio*b+(1-ratio)*a;
		
	}
	
}
