package structures;

public class Point implements Comparable<Point> {
	
	public Point(double x_, double y_) {
		x=x_;
		y=y_;
	}
	
	@Override
	public int compareTo(Point p) {
		if(x!=p.x) {return x>p.x ? 1 : -1;}
		return y>p.y ? 1 : y<p.y ? -1 : 0;
	}
	
	@Override
	public boolean equals(Object p) {
		return p!=null && getClass()==p.getClass() && equals((Point)p);
	}
	
	public boolean equals(Point p) {
		return p!=null && x==p.x && y==p.y;
	}
	
	public String toString() {
		return x+","+y;
	}
	
	public double x;
	public double y;
	
}
