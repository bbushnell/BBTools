package bin;

import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteStreamWriter;
import shared.Tools;

public class ChartMaker {
	
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
