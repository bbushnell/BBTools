package structures;

/**
 * Represents a 2D coordinate point with double precision.
 * Implements Comparable for natural ordering by x-coordinate first, then y-coordinate.
 * Used throughout BBTools for spatial positioning and geometric calculations in genomic analysis.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class Point implements Comparable<Point> {
	
	/**
	 * Constructs a Point with specified x and y coordinates.
	 * @param x_ The x-coordinate value
	 * @param y_ The y-coordinate value
	 */
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
	
	/**
	 * Tests equality with another Point using coordinate comparison.
	 * @param p The Point to compare against
	 * @return true if both Points have identical x and y coordinates
	 */
	public boolean equals(Point p) {
		return p!=null && x==p.x && y==p.y;
	}
	
	/** Returns a string representation of this Point in "x,y" format.
	 * @return String representation as "x,y" */
	public String toString() {
		return x+","+y;
	}
	
	/** The x-coordinate of this Point */
	public double x;
	/** The y-coordinate of this Point */
	public double y;
	
}
