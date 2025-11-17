package fun;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;

import fileIO.TextFile;
import shared.Timer;

/**
 * Graph pathfinding utility that finds shortest paths between nodes using
 * Dijkstra's algorithm. Reads weighted edges from tab-delimited files and
 * constructs bidirectional graphs for optimal path computation.
 *
 * Input format: source_node TAB target_node TAB weight
 *
 * @author Brian Bushnell
 */
public class FindPath {

	/**
	 * Main entry point for pathfinding operations. Expects three command-line
	 * arguments: start node name, stop node name, and graph file path.
	 * Loads graph from file, computes shortest path, and outputs result with timing.
	 * @param args Command-line arguments [start_node, stop_node, graph_file]
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		String start=args[0];
		String stop=args[1];
		String fname=args[2];
		
		makeGraph(fname);
		Path path;
		if(!start.equals(stop)){
			path=findPath(map.get(start), map.get(stop));
		}else{
			path=new Path(new Node(start));
		}
		printPath(path);
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	/**
	 * Implements Dijkstra's algorithm to find shortest path between two nodes.
	 * Uses breadth-first expansion with distance tracking to guarantee optimal paths.
	 * Maintains path mapping and processes nodes in distance-sorted order.
	 *
	 * @param start Starting node for pathfinding
	 * @param stop Target destination node
	 * @return Path object containing node sequence and total distance, null if unreachable
	 */
	private static Path findPath(Node start, Node stop) {
		HashMap<Node, Path> pmap=new HashMap<Node, Path>();
		pmap.put(start, new Path(start));
		LinkedHashSet<Node> seen=new LinkedHashSet<Node>();
		seen.add(start);
		
		while(seen.size()>0){
			LinkedHashSet<Node> seen2=new LinkedHashSet<Node>();
			for(Node n : seen){
				Path current=pmap.get(n);
				for(Edge e : n.edges){
					Path p=pmap.get(e.b);
					if(p==null || p.dist>current.dist+e.dist){
						p=current.copy();
						p.add(e);
						pmap.put(e.b, p);
						seen2.add(e.b);
					}
				}
			}
			seen=seen2;
		}
		return pmap.get(stop);
	}

	/**
	 * Outputs path results in comma-separated format with total distance.
	 * Prints "Unreachable." for null paths indicating no valid route exists.
	 * @param path Path to display, may be null for unreachable destinations
	 */
	private static void printPath(Path path) {
		if(path==null){
			System.out.println("Unreachable.");
			return;
		}
		String comma="";
		for(Node n : path.list){
			System.out.print(comma+n.name);
			comma=",";
		}
		System.out.println("  \t"+path.dist);
	}

	/**
	 * Constructs bidirectional graph from tab-delimited input file.
	 * Each line format: source TAB destination TAB weight
	 * Creates symmetric edges for undirected graph representation.
	 * @param fname Path to graph definition file
	 */
	static void makeGraph(String fname){
		map=new HashMap<String, Node>();
		TextFile tf=new TextFile(fname);
		String line=tf.nextLine();
		while(line!=null){
			String[] split=line.split("\t");
			Node a=fetch(split[0]), b=fetch(split[1]);
			int dist=Integer.parseInt(split[2]);
			a.edges.add(new Edge(a, b, dist));
			b.edges.add(new Edge(b, a, dist));
			line=tf.nextLine();
		}
	}
	
	/**
	 * Retrieves existing node or creates new one for given name.
	 * Maintains global node registry to ensure unique instances per name.
	 * @param s Node name identifier
	 * @return Node instance for the given name
	 */
	static Node fetch(String s){
		Node n=map.get(s);
		if(n==null){
			n=new Node(s);
			map.put(s, n);
		}
		return n;
	}
	
	/** Global registry mapping node names to Node instances */
	static HashMap<String, Node> map;
	
	/** Represents a graph node with name identifier and adjacency list.
	 * Stores outgoing edges for pathfinding traversal. */
	static class Node{
		
		/** Creates node with specified name identifier.
		 * @param s Unique name for this node */
		Node(String s){
			name=s;
		}
		/** Unique identifier for this node */
		String name;
		/** List of outgoing edges from this node */
		ArrayList<Edge> edges=new ArrayList<Edge>();
		
	}
	
	/** Represents weighted edge connecting two nodes in the graph.
	 * Stores source node, destination node, and traversal cost. */
	static class Edge{
		/**
		 * Creates weighted edge between two nodes.
		 * @param a_ Source node
		 * @param b_ Destination node
		 * @param dist_ Edge traversal cost
		 */
		Edge(Node a_, Node b_, int dist_){
			a=a_;
			b=b_;
			dist=dist_;
		}
		Node a, b;
		/** Weight/cost for traversing this edge */
		int dist;
	}
	
	/** Represents a sequence of nodes forming a path through the graph.
	 * Tracks cumulative distance and provides path manipulation methods. */
	static class Path{
		/** Creates path starting with specified node and zero distance.
		 * @param start Initial node for the path */
		Path(Node start){
			list.add(start);
		}
		/** Default constructor for internal path copying operations. */
		private Path() {
			// TODO Auto-generated constructor stub
		}
		/** Extends path by adding edge destination and updating total distance.
		 * @param e Edge to append to current path */
		public void add(Edge e){
			list.add(e.b);
			dist+=e.dist;
		}
		/** Creates deep copy of current path with same node sequence and distance.
		 * @return Independent copy of this path */
		public Path copy(){
			Path p=new Path();
			p.list.addAll(list);
			p.dist=dist;
			return p;
		}
		/** Ordered sequence of nodes forming the path */
		public ArrayList<Node> list=new ArrayList<Node>();
		/** Cumulative distance/cost for traversing entire path */
		public int dist=0;
	}
	
}
