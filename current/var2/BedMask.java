package var2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;

/**
 * An immutable set of genomic intervals loaded from a BED file, supporting fast
 * "is this position covered?" queries.  Used to restrict VCF processing to a set
 * of regions (e.g. a GIAB high-confidence benchmark BED).
 *
 * BED coordinates are 0-based half-open: an interval [start, end) covers the
 * 1-based positions start+1 through end.  A 1-based VCF POS is therefore inside
 * the interval iff start &lt; POS &lt;= end.
 *
 * Intervals are sorted and merged per scaffold on load, so {@link #contains} is a
 * single binary search and is correct even if the input BED is unsorted or has
 * overlapping intervals.  After construction the instance is read-only and its
 * {@link #contains} method is safe to call from multiple threads.
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class BedMask {

	/**
	 * Loads and indexes intervals from a BED file.
	 * @param fname Path to a BED file (may be gzipped); columns: scaffold, start, end (tab-delimited)
	 */
	public BedMask(String fname){
		//Accumulate raw intervals per scaffold
		HashMap<String, ArrayList<int[]>> tmp=new HashMap<String, ArrayList<int[]>>();
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(fname, FileFormat.TXT, null, true, false));
		long raw=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}//Comment line
			if(Tools.startsWith(line, "track") || Tools.startsWith(line, "browser")){continue;}//BED header line
			String[] f=new String(line).split("\t");
			if(f.length<3){continue;}//Not a valid BED record
			String scaf=f[0];
			int start=Integer.parseInt(f[1]);
			int stop=Integer.parseInt(f[2]);
			ArrayList<int[]> list=tmp.get(scaf);
			if(list==null){list=new ArrayList<int[]>(); tmp.put(scaf, list);}
			list.add(new int[] {start, stop});
			raw++;
		}
		bf.close();

		//Sort and merge each scaffold's intervals into parallel start/end arrays
		long merged=0;
		for(Map.Entry<String, ArrayList<int[]>> e : tmp.entrySet()){
			ArrayList<int[]> list=e.getValue();
			Collections.sort(list, new Comparator<int[]>(){
				@Override
				public int compare(int[] x, int[] y){return Integer.compare(x[0], y[0]);}
			});
			ArrayList<int[]> mlist=new ArrayList<int[]>(list.size());
			int[] cur=null;
			for(int[] iv : list){
				if(cur==null){cur=new int[] {iv[0], iv[1]};}
				else if(iv[0]<=cur[1]){cur[1]=Tools.max(cur[1], iv[1]);}//Overlapping or abutting: extend
				else{mlist.add(cur); cur=new int[] {iv[0], iv[1]};}
			}
			if(cur!=null){mlist.add(cur);}
			int n=mlist.size();
			int[] starts=new int[n], stops=new int[n];
			for(int i=0; i<n; i++){starts[i]=mlist.get(i)[0]; stops[i]=mlist.get(i)[1];}
			startMap.put(e.getKey(), starts);
			stopMap.put(e.getKey(), stops);
			merged+=n;
		}
		intervalsLoaded=merged;
		rawIntervals=raw;
		scaffolds=tmp.size();
	}

	/**
	 * Tests whether a 1-based position on a scaffold falls inside any interval.
	 * @param scaf Scaffold/chromosome name (must match the BED's naming)
	 * @param pos 1-based position
	 * @return true if the position is covered by an interval
	 */
	public boolean contains(String scaf, int pos){
		final int[] starts=startMap.get(scaf);
		if(starts==null){return false;}//No intervals on this scaffold
		final int[] stops=stopMap.get(scaf);
		//Binary search for the rightmost interval whose start < pos
		int lo=0, hi=starts.length-1, ans=-1;
		while(lo<=hi){
			int mid=(lo+hi)>>>1;
			if(starts[mid]<pos){ans=mid; lo=mid+1;}
			else{hi=mid-1;}
		}
		//Inside iff that interval also ends at or after pos (start < pos <= end)
		return ans>=0 && pos<=stops[ans];
	}

	/** Number of merged intervals loaded across all scaffolds. */
	public long intervalsLoaded(){return intervalsLoaded;}
	/** Number of raw BED records read before merging. */
	public long rawIntervals(){return rawIntervals;}
	/** Number of distinct scaffolds with at least one interval. */
	public int scaffolds(){return scaffolds;}

	/** Sorted, merged interval starts (0-based) per scaffold. */
	private final HashMap<String, int[]> startMap=new HashMap<String, int[]>();
	/** Sorted, merged interval ends (0-based, exclusive) per scaffold. */
	private final HashMap<String, int[]> stopMap=new HashMap<String, int[]>();

	private long intervalsLoaded=0;
	private long rawIntervals=0;
	private int scaffolds=0;

}
