package fun;

import java.util.LinkedList;
import java.util.Random;

/**
 * Implements Conway's Game of Life cellular automata simulation.
 * Creates a 2D grid of cells that evolve through discrete generations based
 * on specific neighborhood rules. Uses a toroidal (wrap-around) grid topology
 * and detects repeating states to prevent infinite loops.
 *
 * @author Brian Bushnell
 */
public class Life {
	
	/**
	 * Program entry point for Conway's Game of Life simulation.
	 * Expects 4 command-line arguments: width, height, rounds, and initial
	 * live cell probability.
	 *
	 * @param args [xdim, ydim, rounds, prob] where xdim=grid width,
	 * ydim=grid height, rounds=max generations, prob=initial
	 * probability of live cells (0.0-1.0)
	 */
	public static void main(String[] args){
		Life life=new Life(args);
		life.simulate();
	}
	/**
	 * Constructs a Life simulation with parameters from command line arguments.
	 * Parses grid dimensions, maximum rounds, and initial cell probability.
	 * @param args Array containing [xdim, ydim, rounds, prob] as strings
	 */
	public Life(String[] args){
		xdim=Integer.parseInt(args[0]);
		ydim=Integer.parseInt(args[1]);
		rounds=Integer.parseInt(args[2]);
		prob=Float.parseFloat(args[3]);
	}
	
	/**
	 * Runs the complete Game of Life simulation.
	 * Initializes the grid, then iterates through generations applying Conway's
	 * rules. Tracks previous states in a queue to detect cycles and terminate
	 * early if a repeating pattern is found.
	 */
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
	
	/**
	 * Creates a deep copy of a 2D integer grid.
	 * @param a The grid to copy
	 * @return New grid with identical cell values
	 */
	int[][] copy(int[][] a){
		int[][] b=new int[xdim][ydim];
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				b[x][y]=a[x][y];
			}
		}
		return b;
	}
	
	/**
	 * Compares two 2D grids for equality by checking all cell values.
	 * @param a First grid to compare
	 * @param b Second grid to compare
	 * @return true if all corresponding cells have identical values
	 */
	boolean equals(int[][] a, int[][] b){
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				if(a[x][y]!=b[x][y]){return false;}
			}
		}
		return true;
	}
	
	/**
	 * Randomly initializes the grid with live and dead cells.
	 * Each cell has a probability 'prob' of being alive (1) and
	 * '1-prob' of being dead (0).
	 */
	void initialize(){
		Random randy=new Random();
		for(int x=0; x<xdim; x++){
			for(int y=0; y<ydim; y++){
				grid[x][y]=(randy.nextFloat()<prob ? 1 : 0);
			}
		}
	}
	
	/**
	 * Calculates the next generation state for all cells and fills the next grid.
	 * Applies Conway's Game of Life rules to determine each cell's fate.
	 * @param nextGrid Grid to populate with next generation cell states
	 * @return Total count of live cells in the next generation
	 */
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
	
	/**
	 * Determines the next state of a cell using Conway's Game of Life rules.
	 * Live cells with 2-3 neighbors survive; dead cells with exactly 3 neighbors
	 * become alive; all other cells die or remain dead.
	 *
	 * @param x X coordinate of the cell
	 * @param y Y coordinate of the cell
	 * @return 1 if cell will be alive next generation, 0 if dead
	 */
	int next(int x, int y){
		int sum=neighbors(x, y);
		return (sum==3 || (sum==2 && grid[x][y]==1)) ? 1 : 0;
	}
	
	/**
	 * Counts live neighbors of a cell using toroidal topology.
	 * The grid wraps around at edges using modular arithmetic, creating
	 * a torus where edge cells have neighbors on the opposite side.
	 *
	 * @param x X coordinate of the cell
	 * @param y Y coordinate of the cell
	 * @return Number of live neighboring cells (0-8)
	 */
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
	
	/**
	 * Displays the current grid state to the terminal.
	 * Clears the screen and prints the grid using '@' for live cells and
	 * space for dead cells, with the current round number as a header.
	 * @param round Current generation number to display
	 */
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
	
	/** Current state grid where each cell is 0 (dead) or 1 (alive) */
	int[][] grid;
	/** Maximum number of generations to simulate */
	/** Height of the simulation grid (number of rows) */
	/** Width of the simulation grid (number of columns) */
	int xdim, ydim, rounds;
	/** Initial probability of a cell being alive during initialization */
	float prob;
	
}
