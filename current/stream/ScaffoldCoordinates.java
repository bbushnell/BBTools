package stream;

import dna.Data;

/**
 * Transforms BBMap index coordinates into scaffold-relative coordinates.
 * Handles coordinate conversion from global index positions to local scaffold positions,
 * validating that alignments fall within single scaffolds and calculating relative positions.
 *
 * @author Brian Bushnell
 * @date Aug 26, 2014
 */
public class ScaffoldCoordinates {

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates an empty ScaffoldCoordinates with default/invalid fields. */
	public ScaffoldCoordinates(){}
	
	/** Initializes coordinates from a mapped Read by calling set(r).
	 * @param r Read to extract coordinates from */
	public ScaffoldCoordinates(Read r){set(r);}
	
	/** Initializes coordinates from a SiteScore by calling set(ss).
	 * @param ss SiteScore alignment */
	public ScaffoldCoordinates(SiteScore ss){set(ss);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets coordinates from a mapped Read; returns true if successful.
	 * @param r Read to extract coordinates from
	 * @return true if set
	 */
	public boolean set(Read r){
		valid=false;
		if(r.mapped()){setFromIndex(r.chrom, r.start, r.stop, r.strand(), r);}
		return valid;
	}
	
	/**
	 * Sets coordinates from a SiteScore alignment.
	 * @param ss SiteScore containing alignment info
	 * @return true if set
	 */
	public boolean set(SiteScore ss){
		return setFromIndex(ss.chrom, ss.start, ss.stop, ss.strand, ss);
	}
	
	/**
	 * Converts BBMap index coordinates to scaffold-relative coordinates, validating single-scaffold alignment.
	 * @param iChrom_ Index chromosome
	 * @param iStart_ Index start
	 * @param iStop_ Index stop
	 * @param strand_ Strand (0/1)
	 * @param o Context object for assertions
	 * @return true if conversion succeeds
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
	
	/** Resets all coordinate fields and validity flag. */
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
	
	public int scafIndex=-1;
	public int iChrom=-1;
	public int iStart=-1, iStop=-1;
	public int start=-1, stop=-1;
	public byte strand=-1;
	public int scafLength=0;
	public byte[] name=null;
	public boolean valid=false;
	
}
