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
 * A read-only index of genomic intervals loaded from a BED file, supporting fast
 * "is this position covered?" queries.  Used to restrict VCF processing to a set
 * of regions (e.g. a GIAB high-confidence benchmark BED).
 *
 * BED coordinates are 0-based half-open: an interval [start, end) covers the
 * 1-based positions start+1 through end.  A 1-based VCF POS is therefore inside
 * the interval iff start &lt; POS &lt;= end.
 *
 * Intervals are sorted and overlapping or abutting intervals are merged per
 * scaffold. For valid intervals, {@link #contains} performs one binary search
 * even when the input was unsorted. Coordinate bounds and interval ordering are
 * not validated during parsing. Only the first three tab-delimited columns are
 * used; later columns, including strand and block structure, are ignored.
 * After construction the index is not modified. Publish the completed instance
 * before sharing it with query threads; queries do not mutate its arrays or maps.
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class BedMask{

	/**
	 * Loads and indexes intervals from a BED file, matching scaffold names exactly.
	 * Skips empty lines, lines beginning with '#', and rows with fewer than three
	 * tab-delimited fields. Header keywords must be complete words; tab-delimited
	 * records with integer coordinates take precedence over the reserved names
	 * "track" and "browser". Numeric fields use Integer.parseInt without trimming.
	 * The reader is closed before indexing, including when parsing fails. Reported
	 * read/close errors prevent construction; a parsing failure retains precedence.
	 * @param fname Path to a BED file (may be gzipped); columns: scaffold, start, end (tab-delimited)
	 * @throws NumberFormatException If a retained record's coordinates are not parseable ints
	 * @throws RuntimeException If the reader reports an input or closure error
	 */
	public BedMask(String fname){
		//Accumulate raw intervals per scaffold
		HashMap<String, ArrayList<int[]>> tmp=new HashMap<String, ArrayList<int[]>>();
		ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(fname, FileFormat.TXT, null, true, false));
		long raw=0;
		Throwable failure=null;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}//Comment line
				String[] f=new String(line).split("\t");
				if(isHeader(line, f)){continue;}//VAR2-039: do not discard scaffold-name prefixes.
				if(f.length<3){continue;}//Not a valid BED record
				String scaf=f[0];
				int start=Integer.parseInt(f[1]);
				int stop=Integer.parseInt(f[2]);
				ArrayList<int[]> list=tmp.get(scaf);
				if(list==null){list=new ArrayList<int[]>(); tmp.put(scaf, list);}
				list.add(new int[] {start, stop});
				raw++;
			}
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			//VAR2-040: close on rejection and never publish a mask after a latched read failure.
			try{
				if(bf.close()){throw new RuntimeException("Error reading or closing BED file: "+fname);}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}

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

	/** Recognizes complete header keywords while giving parseable BED rows precedence. */
	private static boolean isHeader(final byte[] line, final String[] fields){
		assert(line.length>0) : "Constructor removes empty lines before BED header classification";
		final int length=Tools.startsWith(line, "track") ? 5 : Tools.startsWith(line, "browser") ? 7 : 0;
		if(length==0 || (line.length>length && line[length]!=' ' && line[length]!='\t')){return false;}
		if(fields.length>=3 && (fields[0].equals("track") || fields[0].equals("browser"))){
			try{
				Integer.parseInt(fields[1]);
				Integer.parseInt(fields[2]);
				return false;
			}catch(NumberFormatException e){/* Header metadata is not an interval. */}
		}
		return true;
	}

	/**
	 * Tests whether a 1-based position on a scaffold falls inside any interval.
	 * @param scaf Scaffold/chromosome name (must match the BED's naming)
	 * @param pos 1-based position
	 * @return true if start &lt; pos &lt;= end for an indexed interval; false for an absent scaffold
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
	/** Number of parsed intervals before merging, excluding skipped input rows. */
	public long rawIntervals(){return rawIntervals;}
	/** Number of distinct scaffolds with at least one interval. */
	public int scaffolds(){return scaffolds;}

	/** Sorted, merged interval starts (0-based) per scaffold. */
	private final HashMap<String, int[]> startMap=new HashMap<String, int[]>();
	/** Sorted, merged interval ends (0-based, exclusive) per scaffold. */
	private final HashMap<String, int[]> stopMap=new HashMap<String, int[]>();

	/** Number of indexed intervals after merging, including retained empty intervals. */
	private long intervalsLoaded=0;
	/** Number of parsed interval rows, including duplicates. */
	private long rawIntervals=0;
	/** Number of exact scaffold names with at least one parsed interval. */
	private int scaffolds=0;

}
