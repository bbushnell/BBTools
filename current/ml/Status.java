package ml;

/**
 * Represents a tracking object for neural network training state with training parameters.
 * Stores key training hyperparameters and network reference for comparison during training.
 * Enables sorting and identification of network instances based on training progress.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
class Status implements Comparable<Status>{

	/**
	 * Constructs a Status object with neural network training parameters.
	 *
	 * @param net_ The neural network being trained
	 * @param epoch_ Current training epoch number
	 * @param alpha_ Learning rate for training
	 * @param anneal_ Annealing factor for learning rate adjustment
	 */
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
	
	/** Neural network being tracked during training */
	final CellNet net;
	/** String representation of the status for display purposes */
	String s;
	/** Current training epoch number */
	final int epoch;
	/** Learning rate parameter for neural network training */
	final float alpha;
	/** Annealing factor for learning rate adjustment during training */
	final float anneal;
	/** Counter for tracking status occurrence or usage frequency */
	int count=1;
	
}
