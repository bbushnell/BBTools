package fun;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;

import fileIO.TextFile;
import shared.Timer;

public class FindPath {

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
	
	static Node fetch(String s){
		Node n=map.get(s);
		if(n==null){
			n=new Node(s);
			map.put(s, n);
		}
		return n;
	}
	
	static HashMap<String, Node> map;
	
	static class Node{
		
		Node(String s){
			name=s;
		}
		String name;
		ArrayList<Edge> edges=new ArrayList<Edge>();
		
	}
	
	static class Edge{
		Edge(Node a_, Node b_, int dist_){
			a=a_;
			b=b_;
			dist=dist_;
		}
		Node a, b;
		int dist;
	}
	
	static class Path{
		Path(Node start){
			list.add(start);
		}
		private Path() {
			// TODO Auto-generated constructor stub
		}
		public void add(Edge e){
			list.add(e.b);
			dist+=e.dist;
		}
		public Path copy(){
			Path p=new Path();
			p.list.addAll(list);
			p.dist=dist;
			return p;
		}
		public ArrayList<Node> list=new ArrayList<Node>();
		public int dist=0;
	}
	
}
