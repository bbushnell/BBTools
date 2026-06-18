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
		//ASYMMETRY with set(SiteScore): an UNMAPPED read skips setFromIndex entirely, so valid is left false but the coordinate fields are NOT cleared - stale values from a prior set() persist. set(SiteScore) always calls setFromIndex, which clears on failure. Callers must check the returned valid before trusting any field.
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
				//Midpoint is safe as the lookup point because isSingleScaffold just guaranteed the WHOLE span [iStart,iStop] lies in one scaffold, so any interior point (incl. the midpoint) maps to that same scaffold.
				//TODO: Possible bug [stream/ScaffoldCoordinates#002] - (iStart+iStop) can overflow int when index coords approach Integer.MAX_VALUE/2, yielding a negative/wrong midpoint; overflow-safe idiom is iStart+(iStop-iStart)/2. Reachability depends on the max index-coordinate size BBMap produces - verify before fixing.
				scafIndex=Data.scaffoldIndex(iChrom, (iStart+iStop)/2);
				name=Data.scaffoldNames[iChrom][scafIndex];
				scafLength=Data.scaffoldLengths[iChrom][scafIndex];
				start=Data.scaffoldRelativeLoc(iChrom, iStart, scafIndex);
				//stop is derived as relativeStart + span (iStop-iStart) rather than a second scaffoldRelativeLoc lookup; valid because the global->scaffold transform is a constant shift, so the span is preserved.
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
