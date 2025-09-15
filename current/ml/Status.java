package ml;

class Status implements Comparable<Status>{

	Status(CellNet net_, int epoch_, float alpha_, float anneal_){
		net=net_;
		epoch=epoch_;
		alpha=alpha_;
		anneal=anneal_;
	}
	
	@Override
	public int compareTo(Status b) {
		if(b==null) {
			return 1;
		}
		return net.compareTo(b.net);
	}
	
	@Override
	public int hashCode(){
		return epoch;
	}
	
	@Override
	public String toString() {return s;}
	
	@Override
	public boolean equals(Object b) {return equals((Status)b);}
	public boolean equals(Status b) {return epoch==b.epoch;} //For hashing only
	
	final CellNet net;
	String s;
	final int epoch;
	final float alpha;
	final float anneal;
	int count=1;
	
}
