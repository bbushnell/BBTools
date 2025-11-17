package bin;

import java.util.Comparator;

/**
 * Comparator for Bin objects to enable sorting based on their ID.
 * Provides a singleton comparator that compares Bin objects by their integer ID.
 * @author Brian Bushnell
 * @date 2025
 */
public class IDComparator implements Comparator<Bin>{

	private IDComparator() {}
	
	@Override
	public int compare(Bin a, Bin b) {
		return a.id()-b.id();
	}
	
	/** Singleton instance of the IDComparator for reuse */
	static final IDComparator comparator=new IDComparator();
	
}
