package dna;

import shared.KillSwitch;
import shared.Tools;

public class MotifProbsN extends Motif {
	
	public static void main(String args[]){
		
		String s1="ATN";
		String s2="CTATGCCCATCTGATGGCATGAGGATGAA";
		
//		if(args.length>0){s1=args[0];}
//		if(args.length>1){s2=args[1];}
		
		MotifProbsN m=makeMotif("Exon Stops MP3", 10, 3, 3);
		
		System.out.println("Made motif "+m.name);
		
		String source=s2;
		
		
		int x=m.countExact(source);
		System.out.println(x+" matches.");
		
		byte[] sbytes=source.getBytes();
		
		for(int i=0; i<s2.length(); i++){
			String sub=s2.substring(i, min(i+m.probs.length, s2.length()));
			float p=m.matchStrength(sbytes, i);
			System.out.println(sub+Tools.format(": \t%.5f ->\t%.5f", p, m.normalize(p)));
		}
		
	}
	
	public static MotifProbsN makeMotif(String name_, int length_, int center_, int n_){
		Matrix mat=Matrix.get(name_);
		assert(mat!=null) : "\nCan't find '"+name_+"' in:\n"+Matrix.keys()+"\n\n";
		float[][] sub=mat.subGrid(center_, length_);
		
//		System.out.println("Found "+name+":\n"+Arrays.toString(sub[preLen]));
		
		assert(sub[0].length==(1<<(2*n_)));
		
		MotifProbsN r=new MotifProbsN(name_, sub, center_, n_);
		
		Matrix percentMatrix=null;
		
		
		try {
			percentMatrix=Matrix.get(name_+", "+r.length+", "+r.center);
		} catch (Exception e) {
			// TODO Auto-generated catch block
//			System.out.println("\nIgnoring missing percentMatrix:\n"+e);
		}
		
		if(percentMatrix!=null){
			r.percentile=percentMatrix.grid[0];
		}
//		r.percentile=percentTable.get(name);
		
		return r;
	}
	
	public MotifProbsN(String name_, float[][] p, int cen, int n){
		super(name_, p.length, cen);
		
		N=n;
		chunk=KillSwitch.allocByte1D(N);
		baseProb=Motif.baseProbN[N];
		
		probs=p;
		importance=positionImportance(probs);
		
		adjustForBaseProb(probs, baseProb);
		
		double pmin=1, pmax=1;
		
		double sum=0;
		for(int i=0; i<p.length; i++){
			for(int j=0; j<p[i].length; j++){
				sum+=p[i][j];
			}
		}
		matrixAvg=(float)(sum/(p.length*p[0].length));
		
		
		//Adjusts for importance
		for(int i=0; i<probs.length; i++){
			for(int j=0; j<probs[i].length; j++){
				probs[i][j]=(float)Math.pow(probs[i][j], 1+(importance[i]*.8));
			}
		}
		
		
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<probs.length; i++){
			int x=maxPos(probs[i]);
			int y=minPos(probs[i]);
			sb.append((char)numberToBase[x>>(2*(N-1))]);

//			pmax*=probs[i][x]*4; //TODO Note the .25; could be an empirical inverse probability, but that causes complications
//			pmin*=probs[i][y]*4;

			pmax*=probs[i][x];
			pmin*=probs[i][y];
			
//			pmax*=Math.pow(probs[i][x], 1+importance[i]);
//			pmin*=Math.pow(probs[i][y], 1+importance[i]);
			
//			pmax*=(probs[i][x]+(matrixAvg*importance[i]*.1f));
//			pmin*=(probs[i][y]+(matrixAvg*importance[i]*.1f));
		}
		

		maxProb=(float)pmax;
		minProb=(float)pmin;

		invProbDif=1f/(maxProb-minProb);
		invLength=1f/(length);
		
		commonLetters=sb.toString();
		
		lettersUpper=commonLetters.toUpperCase().getBytes();
		lettersLower=commonLetters.toLowerCase().getBytes();
		
		numbers=new byte[commonLetters.length()];
		numbersExtended=new byte[commonLetters.length()];
		
		for(int i=0; i<lettersUpper.length; i++){
			byte b=lettersUpper[i];
			numbers[i]=baseToNumber[b];
			numbersExtended[i]=baseToNumberExtended[b];
		}
		
	}
	
	
	public void adjustForBaseProb(float[][] grid, float[] base){
		for(int i=0; i<grid.length; i++){
			for(int j=0; j<grid[i].length; j++){
				grid[i][j]/=base[j];
			}
		}
	}
	
	
	@Override
	public float normalize(double strength){
		double r=strength-minProb;
//		r=r/(maxProb-minProb);
//		r=Math.pow(r, 1d/length);
		r=r*invProbDif;
		r=Math.pow(r, invLength);
		return (float)r;
	}
	
	
	public float normalize2(double strength){
		double r=Math.log(strength)-Math.log(minProb);
		
		double r2=Math.log(maxProb)-Math.log(minProb);
		
		r=r/r2;
		return (float)r;
	}
	
	
	@Override
	public boolean matchesExactly(byte[] source, int a){
		
		a=a-center;
		if(a<0 || a+length>source.length){return false;}
		
		for(int i=0; i<lettersUpper.length; i++){
			int x=i+a;
			if(source[x]!=lettersUpper[i] && source[x]!=lettersLower[i]){
				return false;
			}
		}
		return true;
	}
	
	
	@Override
	public float matchStrength(byte[] source, int a){
		
		a=a-center;
		if(a<0 || a+length+1>source.length){return minProb;}
		
		float r=1;
		
		for(int i=0; i<probs.length; i++){
			int x=i+a;
			
			for(int c=0; c<N; c++){
				chunk[c]=source[x+c];
			}
			
			int n=AminoAcid.baseTupleToNumber(chunk);
			if(n<0 || n>baseProb.length){return minProb;}
			
//			float p1=(probs[i][n]+(matrixAvg*importance[i]*.1f));
			
//			float p1=(float)Math.pow(probs[i][n], 1+importance[i]); //Note:  Assumes (A,C,G,T) only.
			float p1=probs[i][n]; //Note:  Assumes (A,C,G,T) only.
			
//			float p2=invBaseProb2[n];
//			float p2=4; //TODO
//
//			r=r*p1*p2;
			
			r=r*p1;
		}
		return r;
	}
	
	
	public float[] positionImportance(float[][] rawProbs){
		float[] base=baseProb;
		float[] out=new float[rawProbs.length];
		
		double maxSum=0;
		
		for(int i=0; i<out.length; i++){
			float[] array=rawProbs[i];
			double sum=0;
			for(int code=0; code<array.length; code++){
				double dif=Math.abs(array[code]-base[code]);
				sum+=Math.pow(dif,1.5); //Raise to a power to increase the effect
			}
			sum=Math.pow(sum, 0.75);
			out[i]=(float)sum;
			if(sum>maxSum){
				maxSum=sum;
			}
		}
		
		for(int i=0; i<out.length; i++){
			out[i]=(float)(out[i]/maxSum);
//			out[i]=out[i]*.9f+.1f; //Weakens the effect
//			out[i]=out[i]*.5f; //makes the scale 0 to .5
		}
		
		return out;
	}

	@Override
	public int numBases() {
		return N;
	}
	
	public final int N;
	
	public final float[][] probs;
	public final float[] importance;
	public final float matrixAvg;

	public final byte[] lettersUpper;
	public final byte[] lettersLower;
	public final byte[] numbers;
	public final byte[] numbersExtended;
	
	private final byte[] chunk;
	private final float[] baseProb;
	
	public float maxProb;
	public float minProb;
	
	public final float invProbDif;
	public final float invLength;
	
}
