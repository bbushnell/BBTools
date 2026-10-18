package ml;

/**
 * Network-status holder retained for older priority queues.
 * TODO: This class has no live constructor references.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
class Status implements Comparable<Status>{

	/**
	 * Constructs a status holder.
	 *
	 * @param net_ neural network being trained
	 * @param epoch_ current training epoch number
	 * @param alpha_ learning rate for training
	 * @param anneal_ annealing factor for learning rate adjustment
	 */
	Status(CellNet net_, int epoch_, float alpha_, float anneal_){
		net=net_;
		epoch=epoch_;
		alpha=alpha_;
		anneal=anneal_;
	}
	
	/** Returns +1 for a null argument; otherwise delegates to net.compareTo. */
	@Override
	public int compareTo(Status b){
		if(b==null){
			return 1;
		}
		return net.compareTo(b.net);
	}
	
	/** Returns the epoch as the hash code. */
	@Override
	public int hashCode(){
		return epoch;
	}
	
	/** Returns the cached status string, which starts null and is not updated here. */
	@Override
	public String toString(){return s;}
	
	@Override
	public boolean equals(Object b){return equals((Status)b);}
	//TODO: equals(Object) throws for null and non-Status objects.
	/** Checks equality by epoch only, matching hashCode but not compareTo. */
	public boolean equals(Status b){return epoch==b.epoch;} //For hashing only
	
	/** Neural network being tracked. */
	final CellNet net;
	/** String representation of the status. */
	String s;
	/** Current training epoch number. */
	final int epoch;
	/** Learning rate parameter. */
	final float alpha;
	/** Annealing factor. */
	final float anneal;
	/** Count field retained with the status holder. */
	int count=1;
	
}
