package align2;

public class Quad64 implements Comparable<Quad64>{
	
	public Quad64(int col_, int row_, int val_){
		column=col_;
		row=row_;
		site=val_;
	}
	
	@Override
	public boolean equals(Object other){
		assert(false);
		return site==((Quad64)other).site;
	}
	
	@Override
	public int hashCode(){return (int)site;}
	
	@Override
	public int compareTo(Quad64 other) {
		return site>other.site ? 1 : site<other.site ? -1 : column-other.column;
//		int x=site-other.site;
//		return(x>0 ? 1 : x<0 ? -1 : column-other.column);
	}
	
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}
	
	public final int column;
	public int row;
	public long site;
	public int list[];
	
}
