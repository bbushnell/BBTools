package clump2;

/**
 * Orders duplicate candidates by shared k-mer, read orientation, k-mer
 * position, and Y-axis flowcell coordinate.
 *
 * @author Brian Bushnell
 */
public class KmerComparatorY extends KmerComparator2{
	
	/** Prevents construction outside the singleton instance. */
	private KmerComparatorY(){}
	
	//Y-axis twin of KmerComparatorX: identical pivot prefix + optical logic,
	//but the tail uses spanTilesY and a.y instead of spanTilesX/a.x
	//(two-pass optical sort, one comparator per flowcell axis).
	//opticalOnly=false skips only the coordinate tail; the shared
	//kmer/strand/position prefix remains active.
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
	
	/** Singleton instance for Y-axis optical sorting. */
	static final KmerComparatorY comparator=new KmerComparatorY();
	
}
