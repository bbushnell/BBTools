package aligner;

import stream.Read;

public interface MicroAligner {
	
	/** Returns identity */
	public float map(Read r);
	
	/** Returns identity */
	public float map(Read r, float minid);
	
}
