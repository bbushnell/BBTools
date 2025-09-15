package bin;

import java.util.Comparator;

public class IDComparator implements Comparator<Bin>{

	private IDComparator() {}
	
	@Override
	public int compare(Bin a, Bin b) {
		return a.id()-b.id();
	}
	
	static final IDComparator comparator=new IDComparator();
	
}
