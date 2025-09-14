package fun;

import java.util.LinkedList;
import java.util.Random;

public class Life {
	
	public static void main(String[] args){
		Life life=new Life(args);
		life.simulate();
	}
	public Life(String[] args){
		xdim=Integer.parseInt(args[0]);
		ydim=Integer.parseInt(args[1]);
		rounds=Integer.parseInt(args[2]);
		prob=Float.parseFloat(args[3]);
	}
	
	void simulate(){
		grid=new int[xdim][ydim];
		int[][] nextGrid=new int[xdim][ydim];
		initialize();
		
		LinkedList<int[][]> queue=new LinkedList<int[][]>();
		
		for(int i=0; i<rounds; i++){

			print(i);
			int count=fill(nextGrid);
			int[][] temp=grid;
			grid=nextGrid;
			nextGrid=temp;
//			if(count<1){break;}
//			if(equals(grid, nextGrid)){break;}

			for(int[][] x : queue){
				if(equals(grid, x)){return;}
			}
			queue.add(copy(grid));
			if(queue.size()>10){queue.poll();}
			
//			long time=System.nanoTime();
//			long next=time+50000000;
//			while(System.nanoTime()<next);
		}
	}
	
	int[][] copy(int[][] a){
		int[][] b=new int[xdim][ydim];
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				b[x][y]=a[x][y];
			}
		}
		return b;
	}
	
	boolean equals(int[][] a, int[][] b){
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				if(a[x][y]!=b[x][y]){return false;}
			}
		}
		return true;
	}
	
	void initialize(){
		Random randy=new Random();
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				grid[x][y]=(randy.nextFloat()<prob ? 1 : 0);
			}
		}
	}
	
	int fill(int[][] nextGrid){
		int count=0;
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				int z=next(x, y);
				nextGrid[x][y]=z;
				count+=z;
			}
		}
		return count;
	}
	
	int next(int x, int y){
		int sum=neighbors(x, y);
		return (sum==3 || (sum==2 && grid[x][y]==1)) ? 1 : 0;
	}
	
	int neighbors(int x, int y){
//		int minX=Tools.max(x-1, 0);
//		int minY=Tools.max(y-1, 0);
//		int maxX=Tools.min(x+1, xdim-1);
//		int maxY=Tools.min(y+1, ydim-1);
		
		int sum=-grid[x][y];
//		for(int i=minX; i<=maxX; i++){
//			for(int j=minY; j<=maxY; j++){
//				sum+=grid[i][j];
//			}
//		}
		for(int i=-1; i<=1; i++){
			for(int j=-1; j<=1; j++){
				sum+=grid[(i+x+xdim)%xdim][(j+y+ydim)%ydim];
			}
		}
		return sum;
	}
	
	void print(int round){
		
		StringBuilder sb=new StringBuilder();
		System.out.print("\033[H\033[2J");
		sb.append("\nRound "+round+"\n");
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				sb.append(grid[x][y]==0 ? ' ' : '@');
			}
			sb.append('\n');
		}
//		System.out.print("\033[H\033[2J");
//		System.out.flush();
		System.out.println(sb);
		System.out.flush();
	}
	
	int[][] grid;
	int xdim, ydim, rounds;
	float prob;
	
}
