package align2;

public class Quad implements Comparable<Quad>{
	
	public Quad(int col_, int row_, int val_){
		column=col_;
		row=row_;
		site=val_;
	}
	
	@Override
	public boolean equals(Object other){
		return site==((Quad)other).site;
	}
	
	@Override
	public int hashCode(){return site;}
	
	@Override
	public int compareTo(Quad other) {
		int x=site-other.site;
		return(x==0 ? column-other.column : x);
	}
	
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}
	
	public final int column;
	public int row;
	public int site;
	public int list[];
	
}
