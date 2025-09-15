package dna;

public class MotifSimple extends Motif {
	
	public static void main(String args[]){
		
		String s1="ATN";
		String s2="ATGCCCATCTGATG";

		if(args.length>0){s1=args[0];}
		if(args.length>1){s2=args[1];}
		
		MotifSimple m=new MotifSimple(s1, 0);
		String source=s2;
		
		
		int x=m.countExtended(source);
		System.out.println(x+" matches.");
	}
	
	public MotifSimple(String s, int cen){
		super(s, s.length(), cen);
		
		commonLetters=s;
		lettersUpper=commonLetters.toUpperCase().getBytes();
		lettersLower=commonLetters.toLowerCase().getBytes();
		
		boolean x=false;
		for(int i=0; i<lettersUpper.length; i++){
			if(lettersUpper[i]!='A' && lettersUpper[i]!='C' && lettersUpper[i]!='G' && lettersUpper[i]!='T'){
				x=true;
			}
		}
		extended=x;

		numbers=new byte[s.length()];
		numbersExtended=new byte[s.length()];
		
		for(int i=0; i<lettersUpper.length; i++){
			byte b=lettersUpper[i];
			numbers[i]=baseToNumber[b];
			numbersExtended[i]=baseToNumberExtended[b];
		}
	}
	
	
	@Override
	public boolean matchesExactly(byte[] source, int a){
		assert(!extended);
		
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
	public boolean matchesExtended(byte[] source, int a){
		
		a=a-center;
		if(a<0 || a+length>source.length){return false;}
		
		for(int i=0; i<lettersUpper.length; i++){
			int x=i+a;
			
			byte s=source[x];
			byte n=baseToNumberExtended[s];
			
			if((n&numbersExtended[i])!=n){
				return false;
			}
		}
		return true;
	}

	@Override
	public int numBases() {
		return numbers.length;
	}
	

	public final byte[] lettersUpper;
	public final byte[] lettersLower;
	public final byte[] numbers;
	public final byte[] numbersExtended;
	
	public final boolean extended;
	
}
