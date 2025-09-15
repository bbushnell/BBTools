package icecream;

public class PolymerTrimmer {
	
	public static boolean parse(String arg, String a, String b){
		if(a.equalsIgnoreCase("minPolymer")){
			minPolymer=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("minFraction")){
			float f=Float.parseFloat(b);
			setMinFraction(f);
		}else if(a.equalsIgnoreCase("polyerror")){
			float f=Float.parseFloat(b);
			setMinFraction(1-f);
		}else{
			return false;
		}
		return true;
	}
	
	public static int testLeft(byte[] bases, char symbol){return testLeft(bases, (byte)symbol);}
	
	public static int testLeft(byte[] bases, byte symbol){
		float score=0;
		float max=0;
		int maxPos=-1;
		for(int i=0; i<bases.length && score>=minScore; i++){
			byte b=bases[i];
			if(b==symbol){
				score++;
				if(score>max){
					max=score;
					maxPos=i;
				}
			}else{
				score-=penalty;
			}
		}
		int trim=maxPos+1;
		return (trim<minPolymer ? 0 : trim);
	}
	
	public static int testRight(byte[] bases, char symbol){return testRight(bases, (byte)symbol);}
	
	public static int testRight(byte[] bases, byte symbol){
		float score=0;
		float max=0;
		int maxPos=bases.length;
		for(int i=bases.length-1; i>=0 && score>=minScore; i--){
			byte b=bases[i];
			if(b==symbol){
				score++;
				if(score>max){
					max=score;
					maxPos=i;
				}
			}else{
				score-=penalty;
			}
		}
		int trim=bases.length-maxPos;
		return (trim<minPolymer ? 0 : trim);
	}
	
	public static void setMinFraction(float f){
		assert(f>=0 && f<=1) : f;
		minFraction=f;
		penalty=(f>=1 ? 99 : ((1f/(1-minFraction))-1));
		minScore=(f>=1 ? 0 : -4*penalty);
	}
	
	static int minPolymer=5;
	private static float minFraction=0.8f;
	private static float penalty=(1f/(1-minFraction))-1;
	private static float minScore=-4*penalty;
	
}
