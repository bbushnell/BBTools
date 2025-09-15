package fileIO;

import dna.Data;
import dna.Gene;
import shared.Shared;

public class ChainLine implements Comparable<ChainLine> {
	
	
	public static void main(String[] args){

		int chrom=Gene.toChromosome(args[0]);

		ChainLine[][] lines=ChainBlock.loadChainLines(Data.ROOT_CHAIN+"hg18ToHg19.over.chain");

		for(int i=1; i<args.length; i++){
			int loc=Integer.parseInt(args[i]);
			int[] result=translate(loc, lines[chrom]);
			System.out.print(chrom+"\t+\t"+loc+"\t->\t");
			System.out.println(result==null ? "null" : result[0]+"\t"+Shared.strandCodes[result[1]]+"\t"+result[2]);
		}
		
	}
	
	
	public ChainLine(int chromT, byte strandT, int startT, int stopT, int chromQ, byte strandQ, int startQ, int stopQ){
		tChrom=chromT;
		tStrand=strandT;
		tStart=startT;
		tStop=stopT;
		
		qChrom=chromQ;
		qStrand=strandQ;
		qStart=startQ;
		qStop=stopQ;
	}
	
	
	@Override
	public String toString(){
		return tChrom+"\t"+Shared.strandCodes[tStrand]+"\t"+tStart+"\t"+tStop+"\t"+
		qChrom+"\t"+Shared.strandCodes[qStrand]+"\t"+qStart+"\t"+qStop;
	}
	
	
	public static int binarySearch(int loc, ChainLine[] array){
		return binarySearch(loc, array, 0, array.length-1);
	}
	
	
	public static int binarySearch(int loc, ChainLine[] array, int first, int last){
//		if(first>=last){
//			if(first>last){return -1;}
//			assert(first==last && first<array.length);
//			return (array[first].tStart<=loc && array[first].tStop>=loc) ? first : -1;
//		}
//		System.out.println("BinarySearch "+loc+", "+first+", "+last);
		if(first>last){return -1;}
		int mid=(first+last)/2;
		ChainLine midcl=array[mid];
//		System.out.println("mid = "+midcl);
		if(loc<midcl.tStart){return binarySearch(loc, array, first, mid-1);}
		else if(loc>midcl.tStop){return binarySearch(loc, array, mid+1, last);}
		return mid;
	}
	
	/** Returns {chrom, strand, loc} */
	public static int[] translate(int loc, ChainLine[] array){
		int index=binarySearch(loc, array);
		if(index<0){return null;}
		ChainLine cl=array[index];
		return cl.translate(loc);
	}
		
	public int[] translate(int loc){
		if(loc<tStart || loc>tStop){return null;}
//		assert(loc>=tStart && loc<=tStop);
		if(qChrom<1 || qChrom>25){return null;}
		if(qStrand==Shared.PLUS){
			return new int[] {qChrom, qStrand, qStart+loc-tStart};
		}else{
			assert(qStart>=qStop) : this;
			return new int[] {qChrom, qStrand, qStart-(loc-tStart)};
		}
	}
	
	
	public boolean contains(int a, int b){
		assert(b>=a);
		return a>=tStart && b<=tStop;
	}
	
	
	public boolean contains(int a){
		return a>=tStart && a<=tStop;
	}


	@Override
	public int compareTo(ChainLine other) {
		int temp;
		
		temp=tChrom-other.tChrom;
		if(temp!=0){return temp;}
		
		assert(tStrand==other.tStrand);
		
		temp=tStart-other.tStart;
		if(temp!=0){return temp;}

		temp=tStop-other.tStop;
		return temp;
	}
	
	public int tChrom;
	public byte tStrand;
	public int tStart;
	public int tStop;
	
	public int qChrom;
	public byte qStrand;
	public int qStart;
	public int qStop;
	
}
