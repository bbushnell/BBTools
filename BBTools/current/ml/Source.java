package ml;

public abstract class Source {
	
	public Source() {}
	
	public float value(){return value;}
	public int id() {return 0;}
	public float bias() {return 0;}
	
	protected float value=0f;
	
	void setValue(float v) {
		value=v;
	}
	
	public abstract boolean check();
	
	/** True only if this has no inputs */
	public abstract boolean terminal();

}
