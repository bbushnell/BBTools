package bin;

import json.JsonObject;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;

public interface Sketchable extends Comparable<Sketchable> {

	public void setFrom(JsonObject jo);
	public Sketch toSketch(SketchMakerMini smm, Read r);
	public void setID(int id);
	public int id();
	public float gc();
	public long size();
	public int taxid();
	public int numContigs();
	public long sketchedSize();
	public void clearTax();
	
}
