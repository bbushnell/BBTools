package stream;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import fileIO.TextFile;
import parse.LineParserS1;
import shared.Tools;

/**
 * Filters SAM/BAM records by overlap with a BED file, using a continuous
 * minimum-overlap-fraction threshold rather than a fixed mode.
 *
 * <p>A read's continuous reference span, excluding clips, is intersected
 * with the merged BED intervals on its scaffold. Let {@code covered} be the number of
 * the read's reference bases that fall inside the BED, and {@code span} its reference
 * length. Deletions and skipped regions contribute to this span; insertions and
 * clipping do not. The read <i>matches</i> the BED when {@code covered>0 && covered>=mof*span}:
 * <ul>
 * <li>{@code mof=0.0} (SamStreamerWrapper default) &rarr; any overlap (the {@code covered>0} guard means
 *     &ge;1bp, since fraction&ge;0 alone is trivially true).</li>
 * <li>{@code mof=0.5} &rarr; at least half of the reference span in the BED.</li>
 * <li>{@code mof=1.0} &rarr; full containment ({@code covered==span}).</li>
 * </ul>
 *
 * <p>{@code include=true} (SamStreamerWrapper default) keeps reads that match the BED; {@code include=false}
 * keeps reads that do not (reverse/exclude mode). {@link #passes(SamLine)} returns the
 * keep/discard decision for this filter alone. Callers combine it with any other
 * filters and decide how to route records; this class does not write output.
 *
 * <p>BED is parsed as 0-based half-open [start,stop); intervals are sorted and merged per
 * scaffold for binary-search overlap, without strand filtering. Scaffold names match
 * exactly. Null or unmapped records, absent scaffolds and nonpositive calculated spans
 * do not match. The private interval storage is unchanged after construction and
 * matching uses local scratch; callers remain responsible for safe publication and
 * avoiding concurrent changes to the supplied SAM records.
 *
 * @author UMP45
 */
public class BedReadFilter{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads BED intervals, sorts them by start and merges overlapping or abutting intervals.
	 * Uses the first three tab-separated columns and ignores later columns. Coordinates
	 * have surrounding whitespace trimmed; scaffold names are retained verbatim.
	 * Empty lines and lines starting with #, track or browser are skipped.
	 * Expects valid BED coordinates; this constructor is not a complete BED validator.
	 * @param bedPath Path to a BED file (&ge;3 columns: scaffold, 0-based start, exclusive stop)
	 * @param minOverlapFraction Value in [0,1]: fraction of the reference span that must lie in the BED
	 *        to count as a match; 0.0=any overlap, 1.0=full containment
	 * @param include true keeps reads matching the BED; false keeps reads that do not
	 */
	public BedReadFilter(String bedPath, double minOverlapFraction, boolean include){
		this.mof=minOverlapFraction;
		this.include=include;
		assert(mof>=0 && mof<=1) : "minoverlapfraction must be in [0,1]: "+mof;

		final HashMap<String, ArrayList<int[]>> tmp=new HashMap<String, ArrayList<int[]>>();
		final LineParserS1 lp=new LineParserS1('\t');
		TextFile tf=new TextFile(bedPath);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.length()==0 || line.charAt(0)=='#' || line.startsWith("track") || line.startsWith("browser")){continue;}
			lp.set(line);
			assert(lp.terms()>=3) : "Expected >=3-column BED (scaffold start stop), got: "+line;
			final String scaffold=lp.parseString(0);
			ArrayList<int[]> list=tmp.get(scaffold);
			if(list==null){list=new ArrayList<int[]>(); tmp.put(scaffold, list);}
			list.add(new int[]{Integer.parseInt(lp.parseString(1).trim()), Integer.parseInt(lp.parseString(2).trim())});
		}
		tf.close();

		map=new HashMap<String, int[][]>(tmp.size()*2);
		for(Map.Entry<String, ArrayList<int[]>> e : tmp.entrySet()){
			map.put(e.getKey(), mergeSorted(e.getValue()));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Sorts the supplied list in place and returns merged intervals as {starts[], stops[]}.
	 * Overlapping and abutting intervals are combined into new arrays. */
	private static int[][] mergeSorted(ArrayList<int[]> list){
		list.sort((a, b)->Integer.compare(a[0], b[0]));
		final ArrayList<int[]> merged=new ArrayList<int[]>(list.size());
		int[] cur=null;
		for(int[] iv : list){
			if(cur==null || iv[0]>cur[1]){cur=new int[]{iv[0], iv[1]}; merged.add(cur);}
			else if(iv[1]>cur[1]){cur[1]=iv[1];}
		}
		final int n=merged.size();
		final int[] starts=new int[n], stops=new int[n];
		for(int i=0; i<n; i++){starts[i]=merged.get(i)[0]; stops[i]=merged.get(i)[1];}
		return new int[][]{starts, stops};
	}

	/**
	 * Applies include/exclude mode to the overlap result.
	 * Null and unmapped records are nonmatches, so exclude mode keeps them.
	 * @param sl SAM record to filter; may be null
	 * @return true when include equals matchesBed(sl): keep matches in include mode, nonmatches otherwise
	 */
	public boolean passes(SamLine sl){return include==matchesBed(sl);}

	/**
	 * Tests positive BED coverage against the configured fraction of the continuous reference span.
	 * The span starts at pos-1 and counts CIGAR M, =, X, D and N operations via SamLine.
	 * Overlapping BED intervals have already been merged, so coverage is counted once.
	 * @param sl SAM record with a CIGAR supported by SamLine; may be null
	 * @return true for positive coverage meeting the threshold; false for null/unmapped
	 * records, an absent scaffold, a nonpositive calculated span or insufficient coverage
	 */
	public boolean matchesBed(SamLine sl){
		if(sl==null || !sl.mapped()){return false;}
		final int[][] iv=map.get(sl.rnameS());
		if(iv==null){return false;}
		final int rs=sl.start(false, false);//0-based inclusive ref start
		final int span=sl.calcCigarLength(false, false);//ref-consuming length
		if(span<=0){return false;}
		final long covered=coveredBases(iv, rs, rs+span);//rs+span = 0-based exclusive end
		return covered>0 && covered>=mof*span;
	}

	/** Counts bases of [rs,re) covered by sorted, disjoint intervals in {starts[], stops[]} form.
	 * Binary search locates the first interval ending after rs; only intersecting intervals are scanned. */
	private static long coveredBases(int[][] iv, int rs, int re){
		final int[] starts=iv[0], stops=iv[1];
		int a=0, b=starts.length;
		while(a<b){int m=(a+b)>>>1; if(stops[m]<=rs){a=m+1;}else{b=m;}}//first interval with stop>rs
		long covered=0;
		for(int i=a; i<starts.length && starts[i]<re; i++){
			final int s=Tools.max(starts[i], rs);
			final int e=Tools.min(stops[i], re);
			if(e>s){covered+=(e-s);}
		}
		return covered;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Merged, sorted intervals per scaffold as {starts[], stops[]}, 0-based half-open. */
	private final HashMap<String, int[][]> map;
	/** Minimum fraction of the read's reference span that must lie in the BED to match. */
	private final double mof;
	/** true=keep reads matching the BED; false=keep reads that do not. */
	private final boolean include;

}
