package bin;

import json.JsonObject;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;

/**
 * Interface defining sketch generation capabilities for genomic objects.
 * Provides contract for objects that can generate MinHash sketches for similarity comparison.
 * Enables polymorphic sketch creation and supports taxonomic metadata access.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public interface Sketchable extends Comparable<Sketchable> {

	/** Initializes this object from JSON data.
	 * @param jo JSON object containing initialization data */
	public void setFrom(JsonObject jo);
	/**
	 * Generates a MinHash sketch from this object and associated read data.
	 * @param smm Sketch maker for generating the MinHash sketch
	 * @param r Associated read data for sketch generation
	 * @return Generated MinHash sketch for similarity comparison
	 */
	public Sketch toSketch(SketchMakerMini smm, Read r);
	/** Sets the numeric identifier for this object.
	 * @param id Numeric identifier to assign */
	public void setID(int id);
	/** Returns the numeric identifier for this object */
	public int id();
	/** Returns the GC content as a fraction between 0 and 1 */
	public float gc();
	/** Returns the size in base pairs of the genomic sequence */
	public long size();
	/** Returns the taxonomic identifier for this object */
	public int taxid();
	/** Returns the number of contigs contained in this object */
	public int numContigs();
	/** Returns the total size in bases that was used for sketch generation */
	public long sketchedSize();
	/** Clears the taxonomic assignment for this object */
	public void clearTax();
	
}
