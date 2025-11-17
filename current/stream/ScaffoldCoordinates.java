package stream;

import dna.Data;

/**
 * Transforms BBMap index coordinates into scaffold-relative coordinates.
 * @author Brian Bushnell
 * @date Aug 26, 2014
 *
 */
public class ScaffoldCoordinates {

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates an empty ScaffoldCoordinates object with default values */
	public ScaffoldCoordinates(){}
	
	/** Creates ScaffoldCoordinates from a mapped Read.
	 * @param r The Read to extract coordinate information from */
	public ScaffoldCoordinates(Read r){set(r);}
	
	/** Creates ScaffoldCoordinates from a SiteScore alignment.
	 * @param ss The SiteScore to extract coordinate information from */
	public ScaffoldCoordinates(SiteScore ss){set(ss);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets coordinates from a mapped Read's alignment information.
	 * Only processes reads that are successfully mapped to a reference.
	 * @param r The Read to extract coordinates from
	 * @return true if coordinates were successfully set, false otherwise
	 */
	public boolean set(Read r){
		valid=false;
		if(r.mapped()){setFromIndex(r.chrom, r.start, r.stop, r.strand(), r);}
		return valid;
	}
	
	/**
	 * Sets coordinates from a SiteScore alignment.
	 * @param ss The SiteScore containing alignment coordinates
	 * @return true if coordinates were successfully set, false otherwise
	 */
	public boolean set(SiteScore ss){
		return setFromIndex(ss.chrom, ss.start, ss.stop, ss.strand, ss);
	}
	
	/**
	 * Sets coordinates from BBMap index positions to scaffold-relative coordinates.
	 * Validates that the alignment falls within a single scaffold and calculates
	 * the relative position within that scaffold. Only succeeds for single-scaffold alignments.
	 *
	 * @param iChrom_ Index chromosome identifier
	 * @param iStart_ Start position in index coordinates
	 * @param iStop_ Stop position in index coordinates
	 * @param strand_ Alignment strand (0=forward, 1=reverse)
	 * @param o Object used for assertion error context (typically Read or SiteScore)
	 * @return true if conversion succeeded and alignment is within single scaffold, false otherwise
	 */
	public boolean setFromIndex(int iChrom_, int iStart_, int iStop_, int strand_, Object o){
		valid=false;
		if(iChrom_>=0){
			iChrom=iChrom_;
			iStart=iStart_;
			iStop=iStop_;
			if(Data.isSingleScaffold(iChrom, iStart, iStop)){
				assert(Data.scaffoldLocs!=null) : "\n\n"+o+"\n\n";
				scafIndex=Data.scaffoldIndex(iChrom, (iStart+iStop)/2);
				name=Data.scaffoldNames[iChrom][scafIndex];
				scafLength=Data.scaffoldLengths[iChrom][scafIndex];
				start=Data.scaffoldRelativeLoc(iChrom, iStart, scafIndex);
				stop=start-iStart+iStop;
				strand=(byte)strand_;
				valid=true;
			}
		}
		if(!valid){clear();}
		return valid;
	}
	
	/** Resets all coordinate fields to invalid/default values.
	 * Sets valid flag to false and clears all position and identifier fields. */
	public void clear(){
		valid=false;
		scafIndex=-1;
		iChrom=-1;
		iStart=-1;
		start=-1;
		iStop=-1;
		stop=-1;
		strand=-1;
		scafLength=0;
		name=null;
		valid=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Index of the scaffold within the chromosome's scaffold array */
	public int scafIndex=-1;
	/** Index chromosome identifier from BBMap indexing */
	public int iChrom=-1;
	public int iStart=-1, iStop=-1;
	public int start=-1, stop=-1;
	/** Alignment strand (0=forward, 1=reverse, -1=unset) */
	public byte strand=-1;
	/** Length of the scaffold containing this alignment */
	public int scafLength=0;
	/** Name of the scaffold containing this alignment */
	public byte[] name=null;
	/** Flag indicating whether coordinate conversion was successful */
	public boolean valid=false;
	
}
