package clump;

public class KmerComparatorY extends KmerComparator2 {
	
	private KmerComparatorY(){}
	
	@Override
	public int compare(ReadKey a, ReadKey b){
//		assert(FlowcellCoordinate.spanTiles || Clump.forceSortXY);
		if(a.kmer!=b.kmer){return a.kmer>b.kmer ? -1 : 1;} //Bigger kmers first...
		if(a.kmerMinusStrand!=b.kmerMinusStrand){return a.kmerMinusStrand ? 1 : -1;}
		if(a.position!=b.position){return a.position<b.position ? 1 : -1;}
		if(Clump.opticalOnly){
			if(a.lane!=b.lane){return a.lane-b.lane;}
			if(!ReadKey.spanTilesY){
				if(a.tile!=b.tile){return a.tile-b.tile;}
			}
			if(a.y!=b.y){return a.y-b.y;}
		}
		return 0;
	}
	
	static final KmerComparatorY comparator=new KmerComparatorY();
	
}
