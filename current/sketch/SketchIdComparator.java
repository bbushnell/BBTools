package sketch;

import java.util.Comparator;

public class SketchIdComparator implements Comparator<Sketch> {

	private SketchIdComparator(){};
	
	@Override
	public int compare(Sketch a, Sketch b) {
		return a.sketchID-b.sketchID;
	}

	public static final SketchIdComparator comparator=new SketchIdComparator();
	
}
