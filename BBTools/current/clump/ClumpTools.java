package clump;

import java.util.ArrayList;

import bloom.KCountArray;
import fileIO.ReadWrite;
import shared.Shared;
import stream.ConcurrentCollectionReadInputStream;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 12, 2015
 *
 */
public class ClumpTools {
	
	public static KCountArray table(){
		return table;
	}
	
	public static synchronized KCountArray getTable(ArrayList<Read> reads, int k, int minCount){
		fname1=fname2=null;
		table=null;
		ConcurrentCollectionReadInputStream cris=new ConcurrentCollectionReadInputStream(reads, null, -1);
		cris.start();
		table=PivotSet.makeKcaStatic(cris, k, minCount, Shared.AMINO_IN);
		ReadWrite.closeStream(cris);
		return table;
	}
	
	public static synchronized KCountArray getTable(String fname1_, String fname2_, int k_, int minCount_){
		if(fname1==null || !fname1.equals(fname1_) || table==null){
			fname1=fname1_;
			fname2=fname2_;
			String[] args=new String[] {"in1="+fname1, "in2="+fname2, "k="+k_, "minCount="+minCount_};
			table=PivotSet.makeSet(args);
		}
		return table;
	}
	
	public static synchronized void clearTable() {
		fname1=fname2=null;
		table=null;
	}
	
	private static String fname1=null, fname2=null;
	private static KCountArray table=null;
	
}
