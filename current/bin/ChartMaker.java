package bin;

import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import shared.Tools;

/**
 * Utility class for generating chart data files from bin statistics.
 * Creates tab-separated output files suitable for plotting bin size distributions,
 * contamination vs completeness plots, and contamination histograms.
 * @author Brian Bushnell
 */
public class ChartMaker {
	
	/**
	 * Generates cumulative bin size chart data sorted by contamination.
	 * Creates a tab-separated file with cumulative size, clean bases, and dirty bases.
	 * Bins are sorted by contamination level using BinStatsComparator.
	 *
	 * @param fname Output filename for the chart data
	 * @param list List of BinStats to process
	 */
	static void makeChartFromBinStats(String fname, ArrayList<BinStats> list) {
		Collections.sort(list, new BinStatsComparator());
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		bsw.print("#Bin\tSize\tClean\tDirty\n");
		double size=0, clean=0, dirty=0;
		int i=0;
		for(BinStats b : list) {
			double contam=b.size*Tools.max(0, b.contam);
			size+=b.size;
			clean+=(b.size-contam);
			dirty+=contam;
			bsw.print(i).tab().print((long)size).tab().print((long)clean).tab().print((long)dirty).println();
			i++;
		}
		bsw.poisonAndWait();
	}
	
	/**
	 * Generates completeness vs contamination plot data.
	 * Creates a tab-separated file with bin index, completeness, contamination, and size.
	 * Bins are sorted by contamination level for consistent ordering.
	 *
	 * @param fname Output filename for the plot data
	 * @param list List of BinStats to process
	 */
	static void writeCCPlot(String fname, ArrayList<BinStats> list) {
		Collections.sort(list, new BinStatsComparator());
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		bsw.print("#Bin\tComplt\tContam\tSize\n");
		int i=0;
		for(BinStats b : list) {
			bsw.print(i).tab().print(b.complt, 4).tab().print(b.contam, 4).tab().print(b.size).println();
			i++;
		}
		bsw.poisonAndWait();
	}
	
	/**
	 * Generates contamination histogram data.
	 * Creates a tab-separated file with contamination level (0.1% increments),
	 * bin count, and total size for each contamination level.
	 *
	 * @param fname Output filename for the histogram data
	 * @param list List of BinStats to process
	 */
	static void writeContamHist(String fname, ArrayList<BinStats> list) {
		Collections.sort(list, new BinStatsComparator());
		int[] count=new int[1001];
		long[] size=new long[1001];
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		bsw.print("#Contam\tCount\tSize\n");
		int max=0;
		for(BinStats b : list) {
			int contam=(int)(b.contam*1000);
			count[contam]++;
			size[contam]+=b.size;
			max=Math.max(contam, max);
		}
		for(int contam=0; contam<=max; contam++) {
			bsw.print(contam*0.1f, 1).tab().print(count[contam]).tab().print(size[contam]).println();
		}
		bsw.poisonAndWait();
	}
	
}
