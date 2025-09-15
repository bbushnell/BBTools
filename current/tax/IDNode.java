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
	
	double[] shorter(double[] a, double[] b){
		return a.length<b.length ? a : b;
	}
	
	double[] longer(double[] a, double[] b){
		return a.length<b.length ? b : a;
	}
	
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
	
	public StringBuilder toNewick(){
		StringBuilder sb=new StringBuilder();
		toNewick(sb);
		return sb;
	}
	
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

	public String name;
	public double[] array;
	public final int number;
	public final int maxPos;
	public final double max;
	public final BitSet bs;
	
	public IDNode parent;
	public final IDNode left;
	public final IDNode right;

}
