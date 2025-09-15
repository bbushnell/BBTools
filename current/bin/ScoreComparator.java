package bin;

import java.util.Comparator;

class ScoreComparator implements Comparator<Bin> {
	
	private ScoreComparator() {}
	
	@Override
	public int compare(Bin a, Bin b) {
		if(a.score!=b.score) {return a.score<b.score ? -1 : 1;}
		if(a.size()!=b.size()) {return a.size()<b.size() ? -1 : 1;}
		return a.id()-b.id();
	}
	
	public static final ScoreComparator comparator=new ScoreComparator();
	
}