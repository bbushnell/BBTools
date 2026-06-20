package clade;


public class SeqIndexResult implements Comparable<SeqIndexResult> {

	SeqIndexResult(Sequence query_, Sequence ref_, float ani_){
		query=query_;
		ref=ref_;
		ani=ani_;
	}
	
	@Override
	public int compareTo(SeqIndexResult o){
		//[clade/SeqIndexResult#001] FIXED - was `ani>o.ani ? 1 : 0`: returned 0 (==) when this.ani<o.ani, violating the comparator contract (a<b gave 0 while b>a gave 1 -> asymmetric -> Collections.sort at SeqIndex:170 misorders / can throw). Now -1. Latent (SeqIndex subsystem is dead). ani is a fraction so no overflow concern.
		if(ani!=o.ani) {return ani>o.ani ? 1 : -1;}
		return ref.compareTo(o.ref);
	}
	
	Sequence query;
	Sequence ref;
	float ani;
	sketch.Comparison c;
	
}
