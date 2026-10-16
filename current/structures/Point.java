package structures;

/**
 * Mutable two-dimensional point ordered by x coordinate, then y coordinate.
 * Object equality requires matching runtime classes and exact coordinate equality.
 * Coordinates are public and unsynchronized; callers sharing or sorting points
 * must coordinate externally. Hash-key callers must not mutate coordinates
 * while a Point is indexed.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class Point implements Comparable<Point>{

	/**
	 * Constructs a point with the specified coordinates.
	 * @param x_ X coordinate.
	 * @param y_ Y coordinate.
	 */
	public Point(double x_, double y_){
		x=x_;
		y=y_;
	}

	/** Orders points by x coordinate, then y coordinate. */
	@Override
	public int compareTo(Point p){
		if(x!=p.x){return x>p.x ? 1 : -1;} //TODO: Possible Bug - x=NaN makes compareTo(this) return -1.
		return y>p.y ? 1 : y<p.y ? -1 : 0;
	}

	@Override
	public boolean equals(Object p){
		return p!=null && getClass()==p.getClass() && equals((Point)p);
	}

	/**
	 * Tests whether another Point has the same exact coordinates.
	 * @param p Point to compare; may be null.
	 * @return True if p has identical x and y coordinates.
	 */
	public boolean equals(Point p){
		return p!=null && x==p.x && y==p.y;
	}

	/**
	 * Returns a hash based on exact coordinate equality.
	 * @return Hash code for this point.
	 */
	@Override
	public int hashCode(){
		return 31*coordinateHash(x)+coordinateHash(y);
	}

	/** Returns a hash for one coordinate, normalizing signed zeroes. */
	private static int coordinateHash(double value){
		final long bits=Double.doubleToLongBits(value==0 ? 0 : value);
		return (int)(bits^(bits>>>32));
	}

	/**
	 * Formats the coordinates separated by one comma.
	 * @return x and y in {@code x,y} form.
	 */
	public String toString(){
		return x+","+y;
	}

	/** X coordinate. */
	public double x;
	/** Y coordinate. */
	public double y;

}
