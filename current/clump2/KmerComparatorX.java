package clump2;

/**
 * Orders duplicate candidates by shared k-mer, read orientation, k-mer
 * position, and X-axis flowcell coordinate.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class KmerComparatorX extends KmerComparator2{
	
	/** Prevents construction outside the singleton instance. */
	private KmerComparatorX(){}
	
	//Within-clump ordering for optical-duplicate detection (Clump.sort,
	//opticalOnly pass). Shared pivot prefix (kmer desc, plus-strand first,
	//position desc) matches ReadKey/KmerComparator. Optical tail: lane, then
	//tile (skipped when spanTilesX - an X-run spans tiles), then x.
	//The lane/tile/x subtractions are safe for normal bounded flowcell ints
	//(ReadKey<-FlowcellCoordinate, Illumina header). opticalOnly=false skips
	//only this coordinate tail. Deliberate X-axis twin of KmerComparatorY
	//(which uses spanTilesY/a.y) - two passes, one per axis.
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
	
	/** Singleton instance for X-axis optical sorting. */
	static final KmerComparatorX comparator=new KmerComparatorX();
	
}
