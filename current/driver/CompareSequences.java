package driver;

import dna.ChromosomeArray;
import shared.Tools;

public class CompareSequences {

	public static void main(String[] args){

		ChromosomeArray cha1=ChromosomeArray.read(args[0]);
		ChromosomeArray cha2=ChromosomeArray.read(args[1]);
		
		long different=0;
		long same=0;
		long nToBase=0;
		long baseToN=0;
		long caseDifferent=0;
		long toUpper=0;
		long toLower=0;
		long difLen=cha2.maxIndex-cha1.maxIndex;
		
		int lim=cha2.maxIndex>cha1.maxIndex ? cha1.maxIndex : cha2.maxIndex;
		
		for(int i=0; i<lim; i++){
			char a=(char) cha1.get(i);
			char b=(char) cha2.get(i);
			if(a==b){
				same++;
			}else{
				different++;
				if(a=='N' && b!='N'){
					nToBase++;
				}else if(a!='N' && b=='N'){
					baseToN++;
				}
				
				if(Tools.toLowerCase(a)==Tools.toLowerCase(b)){
					caseDifferent++;
					if(a==Tools.toLowerCase(a)){
						toUpper++;
					}else{
						toLower++;
					}
				}
				
			}
		}
		
		same+=caseDifferent;
		different-=caseDifferent;
		
		System.out.println("Length Difference: "+difLen);
		System.out.println("Same bases:        "+same);
		System.out.println("Different bases:   "+different+" ("+(100f*different/(float)(different+same))+"%)");
		System.out.println("Base-to-N:         "+baseToN);
		System.out.println("N-To-Base:         "+nToBase);
		System.out.println("Changed case:      "+caseDifferent);
		System.out.println("toUpperCase:       "+toUpper);
		System.out.println("toLowerCase:       "+toLower);
		
	}
	
}
