package clade;


public class SeqIndexResult implements Comparable<SeqIndexResult> {

	SeqIndexResult(Sequence query_, Sequence ref_, float ani_){
		query=query_;
		ref=ref_;
		ani=ani_;
	}
	
	@Override
	public int compareTo(SeqIndexResult o){
		if(ani!=o.ani) {return ani>o.ani ? 1 : 0;}
		return ref.compareTo(o.ref);
	}
	
	Sequence query;
	Sequence ref;
	float ani;
	sketch.Comparison c;
	
}
