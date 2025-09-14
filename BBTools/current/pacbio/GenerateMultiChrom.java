package pacbio;

import java.io.File;
import java.util.Random;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 16, 2012
 *
 */
public class GenerateMultiChrom {
	
	public static void main(String[] args){
		
		ChromosomeArray cha=null;
		
		try {
			int genomeIn=Integer.parseInt(args[0]);
			Data.setGenome(genomeIn);
			cha=Data.getChromosome(1);
			Data.unload(1, true);
		} catch (NumberFormatException e) {
			//ignore
		}
		
		if(cha==null){
			String inname=args[0];
			cha=ChromosomeArray.read(inname, 1);
		}
		
		assert(cha!=null);
		
		int copies=Integer.parseInt(args[1]);
		int build=Integer.parseInt(args[2]);
		
		int mincontig=-1;
		int maxcontig=-1;
		int buffer=-1;
		if(args.length>3){
			mincontig=Integer.parseInt(args[3]);
			maxcontig=Integer.parseInt(args[4]);
			buffer=Integer.parseInt(args[5]);
			System.out.println("Multichrom will be overlayed with blocks of "+buffer+" 'N'");
		}
		
		
//		String pattern=ROOT_GENOME+GENOME_BUILD+"/chr"+chrom+".chromC";
		
		File f=new File(Data.ROOT_GENOME+build);
		if(!f.exists()){f.mkdirs();}
		
		for(int i=1; i<=copies; i++){
			ChromosomeArray chb=makeSynthetic(cha, i);
			if(buffer>0){
				addN(chb, mincontig, maxcontig, buffer);
			}
			ReadWrite.write(chb, Data.ROOT_GENOME+build+"/chr"+i+Data.chromExtension(), false);
		}
		FastaToChromArrays2.writeInfo(build, copies, Data.name, "multiple_"+Data.GENOME_BUILD, false, false);
		
	}
	
	private static void addN(ChromosomeArray cha, int minContig, int maxContig, int buffer){
		
		final int spread=maxContig-minContig+1;
		final Random randy=new Random(cha.chromosome);
		final int lim=cha.maxIndex-Tools.max(maxContig, minContig+buffer);

		int contig=0;
		int nextContig=minContig+randy.nextInt(spread);
		
		for(int i=0; i<lim; i++){
			byte b=cha.get(i);
			if(b=='N'){contig=0;}
			else{
				contig++;
				if(contig>=nextContig){
					contig=0;
					int lim2=i+buffer;
					while(i<lim2){
						cha.set(i, 'N');
						i++;
					}
					nextContig=minContig+(randy.nextInt(spread)+randy.nextInt(spread))/2;
				}
			}
		}
	}

	/**
	 * @param cha
	 * @param i
	 * @return
	 */
	private static ChromosomeArray makeSynthetic(ChromosomeArray cha, int chrom) {
//		assert(false) : cha.array.length+", "+cha.maxIndex;
		ChromosomeArray chb=new ChromosomeArray(chrom, Shared.PLUS, cha.minIndex, cha.array.length+40);
		chb.maxIndex=-1;
		
		int dif=0;
		final int MIN_DIF=-12;
		final int MAX_DIF=12;
		final int INDEL_PERCENT=10;
		final int SUB_PERCENT=1;
		final int ERROR_PERCENT=INDEL_PERCENT+SUB_PERCENT;
		final int ERROR_LENGTH=3;
		
		Random randy=new Random(chrom);
		
		int a=cha.minIndex;
		int b=chb.minIndex;
		
		while(a<=cha.array.length){
			byte c=cha.get(a);
			int x=(c=='N' ? 100 : randy.nextInt(100));
			if(x>=ERROR_PERCENT){ //No error
				chb.set(b, c);
				a++;
				b++;
			}else if(x>=INDEL_PERCENT){//sub
				byte e=c;
				while(e==c){
					e=AminoAcid.numberToBase[randy.nextInt(4)];
				}
				chb.set(b, e);
				a++;
				b++;
			}else{//indel
				boolean ins=Tools.nextBoolean(randy);
				int len=Tools.min(randy.nextInt(ERROR_LENGTH), randy.nextInt(ERROR_LENGTH), randy.nextInt(ERROR_LENGTH+1))+1;
				if(ins && dif+len>MAX_DIF){
					ins=false;
				}else if(!ins && dif-len<MIN_DIF){
					ins=true;
				}
				
				if(ins){
					for(int i=0; i<len; i++){
						boolean same=randy.nextFloat()<0.6f; //Additional 60% chance that inserted base will be a duplicate of an existing base
						int n=randy.nextInt(4);
						byte e=(same ? c : AminoAcid.numberToBase[n]);
						chb.set(b, e);
						b++;
						dif++;
					}
				}else{
					a+=len;
					dif-=len;
				}
			}
		}
		
		return chb;
	}
	
	
	
}
