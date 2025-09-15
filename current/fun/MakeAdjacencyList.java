package fun;

import java.util.Arrays;
import java.util.Random;

import fileIO.TextStreamWriter;
import shared.PreParser;
import shared.Shared;
import structures.ByteBuilder;

public class MakeAdjacencyList {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
		}
		
		parse(args);
		int[][] matrix=genMatrix();
		writeMatrix(matrix);
	}
	
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
	
	public static String toString(int number){
		ByteBuilder sb=new ByteBuilder();
		while(number>0){
			int x='A'+number%26;
			sb.append((char)x);
			number=number/26;
		}
		return (sb.length()<1 ? "A" : sb.reverseInPlace().toString());
	}
	
	public static int nodes=10;
	public static int minlen=5;
	public static int maxlen=25;
	public static float prob=0.3f;
	public static long seed=-1;
	public static String out="stdout.txt";
	
}
