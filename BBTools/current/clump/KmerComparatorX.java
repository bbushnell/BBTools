package clump;

public class KmerComparatorX extends KmerComparator2 {
	
	private KmerComparatorX(){}
	
	@Override
	public int compare(ReadKey a, ReadKey b){
//		assert(FlowcellCoordinate.spanTiles || Clump.forceSortXY);
		if(a.kmer!=b.kmer){return a.kmer>b.kmer ? -1 : 1;} //Bigger kmers first...
		if(a.kmerMinusStrand!=b.kmerMinusStrand){return a.kmerMinusStrand ? 1 : -1;}
		if(a.position!=b.position){return a.position<b.position ? 1 : -1;}
		if(Clump.opticalOnly){
			if(a.lane!=b.lane){return a.lane-b.lane;}
			if(!ReadKey.spanTilesX){
				if(a.tile!=b.tile){return a.tile-b.tile;}
			}
			if(a.x!=b.x){return a.x-b.x;}
		}
		return 0;
	}
	
	static final KmerComparatorX comparator=new KmerComparatorX();
	
}
