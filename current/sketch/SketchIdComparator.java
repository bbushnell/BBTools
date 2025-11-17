package sketch;

import java.util.Comparator;

/**
 * Comparator for sorting Sketch objects by their sketch ID.
 * Provides integer-based ordering of sketches using their sketchID field.
 * @author Brian Bushnell
 */
public class SketchIdComparator implements Comparator<Sketch> {

	/** Private constructor to enforce singleton pattern usage */
	private SketchIdComparator(){};
	
	@Override
	public int compare(Sketch a, Sketch b) {
		return a.sketchID-b.sketchID;
	}

	/** Singleton instance of the comparator for reuse across the application */
	public static final SketchIdComparator comparator=new SketchIdComparator();
	
}
