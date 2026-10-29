package tax;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.PriorityQueue;

import shared.Tools;

/**
 * Mutable cluster node for greedy maximum-identity agglomeration.
 * Leaves borrow lower-triangular identity rows; internal nodes copy and combine
 * those rows by elementwise maximum, clearing entries for leaves already merged.
 * Nodes and their public arrays/bitsets are caller-owned and not thread-safe.
 * This is an identity clustering tree, not an NCBI taxonomy node.
 *
 * @author Brian Bushnell
 * @date July 1, 2016
 */
public class IDNode implements Comparable<IDNode>{
	
	/** Greedily merges the highest-priority unparented clusters.
	 * Mutates parent links of the supplied nodes and newly created internal nodes.
	 * When a zero-similarity cluster's row maximum points within its own
	 * membership, the next active heap root is used to connect disconnected
	 * components.
	 * @param nodes Nonempty leaves in row order, with numbers matching array indices
	 * @return Remaining cluster root
	 * @throws IllegalArgumentException If the leaf array is empty */
	public static IDNode makeTree(IDNode[] nodes){
		
		PriorityQueue<IDNode> heap=new PriorityQueue<IDNode>(nodes.length);
		ArrayList<IDNode> list=new ArrayList<IDNode>(2*nodes.length);
		for(IDNode n : nodes){
			list.add(n);
			heap.add(n);
		}
		
		while(heap.size()>1){
			IDNode a=heap.poll();
//			System.err.println("Found A node "+a);
			if(a.parent==null){
				IDNode b=nodes[a.maxPos];
				if(b.parent!=null){
//					System.err.println("Skipped node "+b);
				}
				while(b.parent!=null){
					b=b.parent;
//					System.err.println("to parent    "+b);
				}
				if(b==a){
					b=nextActiveRoot(heap);
					if(b==null){return a;}
				}
//				System.err.println("Found B node "+b);
				IDNode c=new IDNode(a, b, list.size());
				list.add(c);
				heap.add(c);
//				System.err.println("Made C node  "+c);
//				System.err.println(c.toNewick());
			}

//			System.err.println();
		}
		
		return heap.poll();
	}

	/** Polls stale parented entries until the next active root is found.
	 * @param heap Remaining candidate roots ordered by maximum similarity
	 * @return The next unparented root, or null if none remain */
	private static IDNode nextActiveRoot(final PriorityQueue<IDNode> heap){
		while(!heap.isEmpty()){
			IDNode b=heap.poll();
			if(b.parent==null){return b;}
		}
		return null;
	}

	/** Orders larger maximum similarities first, breaking ties by node number.
	 * Similarities must not be NaN for comparator ordering to remain meaningful.
	 * @param idn Other cluster
	 * @return Negative for higher priority, zero when both maximum and number match */
	@Override
	public int compareTo(IDNode idn){
		if(max==idn.max){return number-idn.number;}
		return max<idn.max ? 1 : -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates one leaf and borrows its identity row without copying.
	 * @param array_ Identities to earlier leaves, indexed by leaf number
	 * @param number_ Nonnegative row/leaf number
	 * @param name_ Newick label, or null for an unlabeled leaf */
	public IDNode(double[] array_, int number_, String name_){
		array=array_;
		number=number_;
		name=name_;
		left=right=null;
		maxPos=(array.length>0 ? Tools.maxIndex(array) : 0);
		max=(maxPos>=array.length ? 0 : array[maxPos]);
		bs=new BitSet(number+1);
		bs.set(number);
	}
	
	/** Selects the shorter array; ties return the second argument.
	 * @param a First array
	 * @param b Second array
	 * @return Borrowed shorter array */
	double[] shorter(double[] a, double[] b){
		return a.length<b.length ? a : b;
	}
	
	/** Selects the longer array; ties return the first argument.
	 * @param a First array
	 * @param b Second array
	 * @return Borrowed longer array */
	double[] longer(double[] a, double[] b){
		return a.length<b.length ? b : a;
	}
	
	/** Merges two distinct unparented clusters and assigns their parent to this node.
	 * Copies the longer row, takes elementwise maxima and clears member-leaf slots.
	 * @param a Higher-priority cluster, with maximum at least that of b
	 * @param b Other cluster, with a distinct identity-row array
	 * @param number_ Number assigned to this internal node */
	public IDNode(IDNode a, IDNode b, int number_){
		assert(a!=b) : a+"; "+a.parent+"; "+a.left+"; "+a.right;
		assert(a.parent==null);
		assert(b.parent==null);
		assert(a.max>=b.max);
		
		number=number_;

		double[] array1=longer(a.array, b.array);
		double[] array2=shorter(a.array, b.array);
		assert(array1!=array2) : a.array.length+", "+b.array.length+", "+a.array+", "+b.array;

		bs=new BitSet();
		bs.or(a.bs);
		bs.or(b.bs);
		array=array1.clone();
		for(int i=0; i<array2.length; i++){
			array[i]=Tools.max(array[i], array2[i]);
		}
		array[a.maxPos]=0;
		for(int i=0; i<array.length; i++){
			if(bs.get(i)){array[i]=0;}
		}
		maxPos=Tools.maxIndex(array);
		max=array[maxPos];
		left=a;
		right=b;
		a.parent=b.parent=this;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Serializes this subtree, with four-decimal branch lengths and no final semicolon.
	 * @return Newly allocated Newick builder */
	public StringBuilder toNewick(){
		StringBuilder sb=new StringBuilder();
		toNewick(sb);
		return sb;
	}
	
	/** Appends recursive topology, sanitized labels and non-root branch lengths.
	 * Labels replace parentheses, colons, commas, semicolons and whitespace with
	 * underscores. Leaves use one minus their parent's maximum child similarity;
	 * internal nodes use maximum child similarity minus their own maximum.
	 * @param sb Destination builder */
	private void toNewick(StringBuilder sb){
		if(left!=null){
			sb.append('(');
			left.toNewick(sb);
			sb.append(',');
			right.toNewick(sb);
			sb.append(')');
		}
		if(name!=null){
			for(int i=0; i<name.length(); i++){
				char c=name.charAt(i);
				if(c=='(' || c==')' || c==':' || c==',' || c==';' || Character.isWhitespace(c)){c='_';}
				sb.append(c);
			}
		}
//		sb.append(Tools.format(":%.4f", max));
		if(parent!=null){
			sb.append(':');
			
			double len;
			if(left==null){
				len=1-Tools.max(parent.left.max, parent.right.max);
			}else{
				len=Tools.max(left.max, right.max)-max;
			}
			
			sb.append(Tools.format("%.4f", len));
//			assert(Tools.max(parent.left.max, parent.right.max)-parent.max<0.4) : parent+"\n"+parent.left+"\n"+parent.right+"\n";
		}
	}
	
	/** Returns string representation with node number, max similarity, and array values.
	 * @return Formatted string showing node details */
	@Override
	public String toString(){
		return "("+number+/*" "+name+*/" "+Tools.format("%.4f", max)+" "+toString(array)+")";
	}
	
	/** Formats an identity row for diagnostics.
	 * @param array Row values
	 * @return Bracketed four-decimal values */
	private static String toString(double[] array){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		sb.append(' ');
		for(double d : array){
			sb.append(Tools.format("%.4f ", d));
		}
		sb.append(']');
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Optional leaf label; internal nodes start unlabeled. */
	public String name;
	/** Borrowed leaf row or independently copied internal-cluster row. */
	public double[] array;
	/** Leaf row index or monotonically assigned internal-cluster number. */
	public final int number;
	/** Index of the largest row value when this node was constructed. */
	public final int maxPos;
	/** Largest row value at construction, or zero for an empty leaf row. */
	public final double max;
	/** Membership bits indexed by original leaf number. */
	public final BitSet bs;
	
	/** Assigned enclosing cluster, or null while this node remains a root. */
	public IDNode parent;
	/** First merged child, or null for a leaf. */
	public final IDNode left;
	/** Second merged child, or null for a leaf. */
	public final IDNode right;

}
