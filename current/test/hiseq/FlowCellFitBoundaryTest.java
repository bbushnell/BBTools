package test.hiseq;

import hiseq.FlowCell;
import hiseq.MicroTile;
import hiseq.Tile;

/** Literal insufficient-variation fit cases and an independent two-point regression control for HIL-110.
 * @author Brian Bushnell, Sucrose */
public class FlowCellFitBoundaryTest{

	/** Runs phase=baseline observations or phase=fixed loud-rejection checks in a dedicated JVM. */
	public static void main(String[] args){
		boolean fixed=false;
		for(String arg : args){
			if(arg.equals("phase=fixed")){fixed=true;}
			else if(!arg.equals("phase=baseline")){throw new IllegalArgumentException(arg);}
		}
		MicroTile.TRACK_CYCLES=false;
		Tile.xSize=Tile.ySize=500;
		System.out.println("case\tread_coefficients\troot_base_coefficients\tstatus");
		control();
		degenerate("one-tile", 1, fixed);
		degenerate("constant-uniqueness", 2, fixed);
	}

	/** Each tile has1200 aligned reads and4800 aligned bases, so the fit gate and eligibility threshold are met. */
	private static void populate(FlowCell flowcell, int index, long hits, long misses){
		final MicroTile tile=flowcell.getMicroTile(1, 1101, index*500, 0);
		tile.readCount=tile.alignedReadCount=1200;
		tile.baseCount=tile.alignedBaseCount=4800;
		tile.hits=hits; tile.misses=misses;
		tile.readErrorCount=index==0 ? 12 : 24;
		tile.baseErrorCount=index==0 ? 48 : 96;
		tile.acgtn[0]=4800;
	}

	/** Two x values25/75 determine the line directly, independently of Tools.linearRegression. */
	private static void control(){
		final FlowCell flowcell=new FlowCell(31);
		populate(flowcell, 0, 3, 1);
		populate(flowcell, 1, 1, 3);
		require(flowcell.calcStats().size()==2, "Both literal control microtiles must be retained");
		final double[] read=flowcell.uniqueToReadErrorRateFormula, base=flowcell.uniqueToBaseErrorRateFormula;
		print("distinct-uniqueness", read, base);
		require(finite(read) && finite(base), "Two distinct uniqueness values must give finite coefficients");
		checkLine(read, 12.2/1201, 24.2/1201, "read error");
		checkLine(base, Math.sqrt(48.002/4801), Math.sqrt(96.002/4801), "root base error");
	}

	/** Baseline exposes nonfinite fits; the repair must reject rather than guess coefficients or return null. */
	private static void degenerate(String name, int count, boolean fixed){
		final FlowCell flowcell=new FlowCell(31);
		for(int i=0; i<count; i++){populate(flowcell, i, 1, 1);}
		try{flowcell.calcStats();}
		catch(IllegalStateException e){
			if(!fixed){throw new AssertionError("Baseline unexpectedly rejected "+name, e);}
			System.out.println(name+"\t-\t-\tPASS: rejected ("+e.getMessage()+")");
			return;
		}
		final double[] read=flowcell.uniqueToReadErrorRateFormula, base=flowcell.uniqueToBaseErrorRateFormula;
		print(name, read, base);
		require(!fixed, "Non-identifiable regression must throw IllegalStateException: "+name);
		require(read!=null && base!=null && read.length==2 && base.length==2,
			"The baseline aligned-read gate must attempt both fits: "+name);
		require(!finite(read) && !finite(base), "Zero x variation must expose nonfinite coefficients: "+name);
	}

	/** A two-point line must interpolate the two literal values at25 and75 percent misses. */
	private static void checkLine(double[] fit, double first, double second, String name){
		final double slope=(second-first)/50, intercept=first-25*slope;
		require(Math.abs(fit[0]-intercept)<1e-10 && Math.abs(fit[1]-slope)<1e-10,
			name+" line differs: observed="+fit[0]+","+fit[1]+" expected="+intercept+","+slope);
	}

	/** Tests all expected coefficients, without accepting an absent fit as a finite fallback. */
	private static boolean finite(double[] fit){
		return fit!=null && fit.length==2 && Double.isFinite(fit[0]) && Double.isFinite(fit[1]);
	}

	/** Prints raw coefficient arrays before the case assertions. */
	private static void print(String name, double[] read, double[] base){
		System.out.println(name+"\t"+java.util.Arrays.toString(read)+"\t"+java.util.Arrays.toString(base)+"\tOBSERVED");
	}

	/** Checks literal expectations independently of the fixture's JVM assertion setting. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
