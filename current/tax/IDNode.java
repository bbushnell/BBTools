package tax;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.PriorityQueue;

import shared.Tools;

/**
 * Support class for IDTree.
 * @author Brian Bushnell
 * @date July 1, 2016
 *
 */
public class IDNode implements Comparable<IDNode>{
	
	/**
	 * Builds a hierarchical tree from an array of leaf nodes using similarity-based clustering.
	 * Uses priority queue to repeatedly merge nodes with highest similarity until single root remains.
	 * @param nodes Array of leaf nodes to cluster into tree
	 * @return Root node of the constructed tree
	 */
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

	@Override
	public int compareTo(IDNode idn) {
		if(max==idn.max){return number-idn.number;}
		return max<idn.max ? 1 : -1;
	}
	
	/**
	 * Creates a leaf node with similarity array and identifier.
	 * Finds maximum similarity value and position within the array.
	 *
	 * @param array_ Similarity values to other entities
	 * @param number_ Unique identifier number for this node
	 * @param name_ Optional name for this entity
	 */
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
	
	/**
	 * Returns the shorter of two arrays.
	 * @param a First array
	 * @param b Second array
	 * @return The array with fewer elements
	 */
	double[] shorter(double[] a, double[] b){
		return a.length<b.length ? a : b;
	}
	
	/**
	 * Returns the longer of two arrays.
	 * @param a First array
	 * @param b Second array
	 * @return The array with more elements
	 */
	double[] longer(double[] a, double[] b){
		return a.length<b.length ? b : a;
	}
	
	/**
	 * Creates internal node by merging two child nodes with updated similarity array.
	 * Takes maximum similarities from both children, zeros out positions for merged entities.
	 * Updates BitSet to track all entities contained within this subtree.
	 *
	 * @param a First child node (must have max >= b.max)
	 * @param b Second child node
	 * @param number_ Unique identifier for this internal node
	 */
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
	
	/** Generates Newick format tree representation.
	 * @return StringBuilder containing complete Newick format string */
	public StringBuilder toNewick(){
		StringBuilder sb=new StringBuilder();
		toNewick(sb);
		return sb;
	}
	
	/**
	 * Recursively builds Newick format string with branch lengths.
	 * Branch lengths calculated from similarity differences between nodes.
	 * Escapes special Newick characters in node names.
	 * @param sb StringBuilder to append Newick format data
	 */
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
	
	@Override
	public String toString(){
		return "("+number+/*" "+name+*/" "+Tools.format("%.4f", max)+" "+toString(array)+")";
	}
	
	/**
	 * Formats double array as bracketed string with fixed decimal precision.
	 * @param array Array of double values to format
	 * @return Formatted string representation of array
	 */
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

	/** Optional name identifier for this node */
	public String name;
	/** Similarity values to other entities in the dataset */
	public double[] array;
	/** Unique identifier number for this node */
	public final int number;
	/** Index position of maximum similarity value in array */
	public final int maxPos;
	/** Maximum similarity value from this node's array */
	public final double max;
	/** BitSet tracking all entity indices contained within this subtree */
	public final BitSet bs;
	
	/** Parent node in the tree hierarchy; null for root node */
	public IDNode parent;
	/** Left child node; null for leaf nodes */
	public final IDNode left;
	/** Right child node; null for leaf nodes */
	public final IDNode right;

}
