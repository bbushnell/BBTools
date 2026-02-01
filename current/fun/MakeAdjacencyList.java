package fun;

import java.util.Arrays;
import shared.Random;

import fileIO.TextStreamWriter;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Generates a randomized adjacency matrix representing a graph with configurable
 * parameters, writing the result to a text file.
 * Creates probabilistic graphs with randomly generated edge weights between nodes.
 * @author Brian Bushnell
 */
public class MakeAdjacencyList {
	
	/**
	 * Main entry point for generating and writing an adjacency matrix.
	 * Parses command-line arguments, generates matrix, and writes to output file.
	 * @param args Command-line arguments for configuration
	 */
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
		}
		
		parse(args);
		int[][] matrix=genMatrix();
		writeMatrix(matrix);
	}
	
	/**
	 * Parses command-line arguments to configure matrix generation parameters.
	 * Supports: out/out1, nodes/n, minlen/min, maxlen/max, prob, seed.
	 * @param args Array of command-line arguments in key=value format
	 * @throws RuntimeException if unknown parameter is encountered
	 */
	public static void parse(String[] args){
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("out") || a.equals("out1")){
				out=b;
			}else if(a.equals("nodes") || a.equals("n")){
				nodes=Integer.parseInt(b);
			}else if(a.equals("minlen") || a.equals("min")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("maxlen") || a.equals("max")){
				maxlen=Integer.parseInt(b);
			}else if(a.equals("prob")){
				prob=Float.parseFloat(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
	}
	
	/**
	 * Generates a random symmetric adjacency matrix for an undirected graph.
	 * Creates edges probabilistically between node pairs with random weights.
	 * @return Symmetric integer matrix where -1 indicates no edge,
	 * positive values indicate edge weights
	 */
	public static int[][] genMatrix(){
		
		final Random randy=Shared.threadLocalRandom(seed);
		final int[][] matrix=new int[nodes][nodes];
		final int range=maxlen-minlen+1;
		for(int[] array : matrix){
			Arrays.fill(array, -1);
		}
		
		for(int i=0; i<nodes; i++){
			for(int j=i+1; j<nodes; j++){
				if(randy.nextFloat()<prob){
					int dist=minlen+(range<1 ? 0 : randy.nextInt(range));
					matrix[i][j]=matrix[j][i]=dist;
				}
			}
		}
		return matrix;
	}
	
	/**
	 * Writes the adjacency matrix to a tab-delimited text file.
	 * Outputs only upper triangle plus diagonal to avoid duplicate edges.
	 * Format: nodeA\tnodeB\tweight for each edge.
	 * @param matrix The adjacency matrix to write
	 */
	public static void writeMatrix(int[][] matrix){
		TextStreamWriter tsw=new TextStreamWriter(out, false, false, false);
		tsw.start();
		for(int i=0; i<nodes; i++){
			for(int j=i; j<nodes; j++){
				int dist=matrix[i][j];
				if(dist>=0){
					tsw.print(toString(i)+"\t"+toString(j)+"\t"+dist+"\n");
				}
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Converts an integer to a base-26 alphabetic string representation.
	 * Uses A-Z characters where A=0, B=1, etc. Handles zero as "A".
	 * @param number Integer to convert (non-negative)
	 * @return Base-26 string representation using letters A-Z
	 */
	public static String toString(int number){
		ByteBuilder sb=new ByteBuilder();
		while(number>0){
			int x='A'+number%26;
			sb.append((char)x);
			number=number/26;
		}
		return (sb.length()<1 ? "A" : sb.reverseInPlace().toString());
	}
	
	/** Number of nodes in the generated graph (default 10) */
	public static int nodes=10;
	/** Minimum edge weight for randomly generated edges (default 5) */
	public static int minlen=5;
	/** Maximum edge weight for randomly generated edges (default 25) */
	public static int maxlen=25;
	/** Probability of creating an edge between any two nodes (default 0.3) */
	public static float prob=0.3f;
	/**
	 * Random seed for reproducible graph generation (default -1 for random seed)
	 */
	public static long seed=-1;
	/** Output file path for the generated adjacency list (default "stdout.txt") */
	public static String out="stdout.txt";
	
}
