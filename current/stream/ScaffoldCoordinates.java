package stream;

import dna.Data;

/**
 * Transforms BBMap index coordinates into scaffold-relative coordinates.
 * Uses loaded Data scaffold metadata and its padding-aware single-scaffold policy,
 * then chooses a scaffold by midpoint and subtracts its offset from both endpoints.
 * This is not a strict check that the resulting span lies within scaffold bases.
 *
 * <p>Instances are mutable scratch objects. Check each setter's result before reading
 * the fields; an unmapped Read invalidates the result without clearing old coordinates.
 * Callers must coordinate access and keep the shared Data metadata stable during use.
 * The name array is borrowed from Data, not copied.
 *
 * @author Brian Bushnell
 * @date Aug 26, 2014
 */
public class ScaffoldCoordinates{

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Creates an empty ScaffoldCoordinates with default/invalid fields. */
	public ScaffoldCoordinates(){}
	
	/** Initializes coordinates from a mapped Read by calling set(r).
	 * @param r Nonnull read; an unmapped read leaves this new instance invalid */
	public ScaffoldCoordinates(Read r){set(r);}
	
	/** Initializes coordinates from a SiteScore by calling set(ss).
	 * @param ss Nonnull SiteScore supplying index coordinates */
	public ScaffoldCoordinates(SiteScore ss){set(ss);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets coordinates from a mapped Read using the loaded scaffold metadata.
	 * An unmapped read sets valid=false but leaves every other field untouched.
	 * @param r Nonnull read to extract coordinates from
	 * @return true if conversion succeeds; false if unmapped or the index-span check rejects it
	 */
	public boolean set(Read r){
		valid=false;
		//ASYMMETRY with set(SiteScore): an UNMAPPED read skips setFromIndex entirely, so valid is left false but the coordinate fields are NOT cleared - stale values from a prior set() persist. set(SiteScore) always calls setFromIndex, which clears on a normal false return. Callers must check the returned valid before trusting any field.
		if(r.mapped()){setFromIndex(r.chrom, r.start, r.stop, r.strand(), r);}
		return valid;
	}
	
	/**
	 * Passes a SiteScore's index coordinates to setFromIndex, without a separate mapped-status check.
	 * @param ss Nonnull SiteScore containing the index span and strand
	 * @return true on conversion; a normal false return clears the fields
	 */
	public boolean set(SiteScore ss){
		return setFromIndex(ss.chrom, ss.start, ss.stop, ss.strand, ss);
	}
	
	/**
	 * Converts an inclusive index span under Data.isSingleScaffold's padding-aware policy.
	 * Requires compatible loaded scaffold locations, names and lengths for a nonnegative chromosome.
	 * The midpoint selects the scaffold; subtracting its offset preserves the supplied span.
	 * Relative coordinates are not clipped or checked against scafLength. A normal false
	 * return clears the fields; errors from invalid metadata are not handled here.
	 * @param iChrom_ Index chromosome; a negative value yields false after clearing
	 * @param iStart_ Zero-based inclusive index start, no greater than iStop_
	 * @param iStop_ Zero-based inclusive index stop
	 * @param strand_ Strand value copied by byte cast, normally 0/1
	 * @param o Context for the metadata assertion message; may be null
	 * @return true after conversion, false for a negative chromosome or rejected index span
	 */
	public boolean setFromIndex(int iChrom_, int iStart_, int iStop_, int strand_, Object o){
		valid=false;
		if(iChrom_>=0){
			iChrom=iChrom_;
			iStart=iStart_;
			iStop=iStop_;
			if(Data.isSingleScaffold(iChrom, iStart, iStop)){
				assert(Data.scaffoldLocs!=null) : "\n\n"+o+"\n\n";
				//Data applies its padding-aware span policy; midpoint lookup selects the scaffold. Strict scaffold-base bounds are not checked here.
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
	
	/** Clears validity and name, sets coordinate/index/strand fields to -1 and scaffold length to zero. */
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
	
	/** Selected scaffold index within the index chromosome. */
	public int scafIndex=-1;
	/** BBMap index chromosome for the supplied span. */
	public int iChrom=-1;
	/** Zero-based inclusive index endpoints as supplied. */
	public int iStart=-1, iStop=-1;
	/** Inclusive scaffold-relative endpoints; not clamped to the scaffold length. */
	public int start=-1, stop=-1;
	/** Copied strand value, normally 0/1, or -1 after clearing. */
	public byte strand=-1;
	/** Selected scaffold's length from Data, or zero after clearing. */
	public int scafLength=0;
	/** Shared name bytes from Data; treat contents as read-only. Null after clearing. */
	public byte[] name=null;
	/** Last setter's success flag; does not validate subsequent manual field changes. */
	public boolean valid=false;
	
}
