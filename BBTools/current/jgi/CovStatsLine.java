package jgi;

import java.util.Arrays;
import java.util.HashMap;

import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Oct 10, 2014
 *
 */
public class CovStatsLine {
	
	public CovStatsLine(String s){
		this(s.split("\t"));
	}

	/**
	 * ID	Avg_fold	Length	Ref_GC	Covered_percent	Covered_bases	Plus_reads	Minus_reads	(optional.... Read_GC)
	 * @param split
	 */
	public CovStatsLine(String[] split) {
		assert(split.length>=8) : Arrays.toString(split);
		assert(!split[0].startsWith("#")) : Arrays.toString(split);
		assert(id_FNUM==0) : "Not initialized with header: "+id_FNUM;

		if(id_FNUM>=0){id=split[id_FNUM];}
		if(length_FNUM>=0){length=Integer.parseInt(split[length_FNUM]);}
		if(coveredBases_FNUM>=0){coveredBases=Integer.parseInt(split[coveredBases_FNUM]);}
		if(plusReads_FNUM>=0){plusReads=Long.parseLong(split[plusReads_FNUM]);}
		if(minusReads_FNUM>=0){minusReads=Long.parseLong(split[minusReads_FNUM]);}
		if(avgFold_FNUM>=0){avgFold=Double.parseDouble(split[avgFold_FNUM]);}
		if(refGC_FNUM>=0){refGC=Double.parseDouble(split[refGC_FNUM]);}
		if(median_FNUM>=0){median=Long.parseLong(split[median_FNUM]);}
		if(underMin_FNUM>=0){underMin=Integer.parseInt(split[underMin_FNUM]);}
		if(readGC_FNUM>=0){readGC=Double.parseDouble(split[readGC_FNUM]);}
		if(stdDev_FNUM>=0){stdDev=Double.parseDouble(split[stdDev_FNUM]);}
		
//		assert(split.length>=8) : Arrays.toString(split);
//		assert(!split[0].startsWith("#")) : Arrays.toString(split);
//		id=split[0];
//		avgFold=Double.parseDouble(split[1]);
//		length=Integer.parseInt(split[2]);
//		refGC=Double.parseDouble(split[3]);
////		coveredPercent=Double.parseDouble(split[4]);
//		coveredBases=Integer.parseInt(split[5]);
//		plusReads=Long.parseLong(split[6]);
//		minusReads=Long.parseLong(split[7]);
//		if(split.length==11){
//			median=Integer.parseInt(split[8]);
//			underMin=Integer.parseInt(split[9]);
//			readGC=Double.parseDouble(split[10]);
//		}else if(split.length==10){
//			median=Integer.parseInt(split[8]);
//			if(CoveragePileup.USE_BITSETS && CoveragePileup.USE_WINDOW){
//				underMin=Integer.parseInt(split[9]);
//			}else{
//				readGC=Double.parseDouble(split[9]);
//			}
//		}else if(split.length==9){
//			readGC=Double.parseDouble(split[8]);
//		}else if(split.length<9){
//			//do nothing
//		}
	}
	
	public final double coveredPercent(){
		return (100.0*coveredBases)/Tools.max(1, length);
	}
	
	public final long reads(){return plusReads+minusReads;}
	
	/**
	 * @param csl
	 */
	public void add(CovStatsLine csl) {
		double invlen2=1d/Tools.max(1, length+csl.length);
		avgFold=((avgFold*length)+(csl.avgFold*csl.length))*invlen2;
		refGC=((refGC*length)+(csl.refGC*csl.length))*invlen2;
		readGC=((readGC*reads())+(csl.readGC*csl.reads()))*1.0/(Tools.max(1, reads()+csl.reads()));
		
		length+=csl.length;
		coveredBases+=csl.coveredBases;
		plusReads+=csl.plusReads;
		minusReads+=csl.minusReads;
		median=median+csl.median;
		underMin=underMin+csl.underMin;
	}
	
	@Override
	public String toString(){
		return Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f", id, avgFold, length,
				refGC, coveredPercent(), coveredBases, plusReads, minusReads, median, underMin, readGC, stdDev);
	}
	
	public static void initializeHeader(String header){
		while(header.startsWith("#")){header=header.substring(1);}
		String[] split=header.split("\t");
		HashMap<String, Integer> map=new HashMap<String, Integer>(23);
		for(int i=0; i<split.length; i++){
			String s=split[i].toLowerCase();
			if(s.startsWith("under_")){s="under_min";}
			map.put(split[i].toLowerCase(), i);
		}
		id_FNUM=map.containsKey("id") ? map.get("id") : -1;
		avgFold_FNUM=map.containsKey("avg_fold") ? map.get("avg_fold") : -1;
		length_FNUM=map.containsKey("length") ? map.get("length") : -1;
		refGC_FNUM=map.containsKey("ref_gc") ? map.get("ref_gc") : -1;
		coveredBases_FNUM=map.containsKey("covered_bases") ? map.get("covered_bases") : -1;
		plusReads_FNUM=map.containsKey("plus_reads") ? map.get("plus_reads") : -1;
		minusReads_FNUM=map.containsKey("minus_reads") ? map.get("minus_reads") : -1;
		median_FNUM=map.containsKey("median_fold") ? map.get("median_fold") : -1;
		underMin_FNUM=map.containsKey("under_min") ? map.get("under_min") : -1;
		readGC_FNUM=map.containsKey("read_gc") ? map.get("read_gc") : -1;
		stdDev_FNUM=map.containsKey("std_dev") ? map.get("std_dev") : -1;
		
		assert(id_FNUM==0) : "Bad header: "+id_FNUM+"\n"+header+"\n"+map;
	}
	
//	public static final String header1="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tMedian_fold\tUnder_min\tRead_GC";
//	public static final String header2="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC";
	
	public String id;
	public int length;
	public int coveredBases;
	public long plusReads;
	public long minusReads;
	public double avgFold;
	public double refGC;
	public long median;
	public int underMin;
	public double readGC;
	public double stdDev;
	
	private static int id_FNUM=-1;
	private static int length_FNUM=-1;
	private static int coveredBases_FNUM=-1;
	private static int plusReads_FNUM=-1;
	private static int minusReads_FNUM=-1;
	private static int avgFold_FNUM=-1;
	private static int refGC_FNUM=-1;
	private static int median_FNUM=-1;
	private static int underMin_FNUM=-1;
	private static int readGC_FNUM=-1;
	private static int stdDev_FNUM=-1;
}
