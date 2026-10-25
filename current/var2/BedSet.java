package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Performs set operations on BED interval files: union, intersection, and
 * subtraction (the first file minus the union of all the rest). Inputs may also
 * be VCF files, converted to padded reference spans. Reports each input's coverage,
 * the union and all-input intersection, and result coverage. Shared and unique
 * coverage per input are additionally reported when there are exactly two inputs.
 * This is the BED analog of {@link CompareVCF}.
 *
 * BED coordinates are 0-based half-open: an interval [start, end) has length
 * end-start and covers no base when start==end.  Intervals are sorted and merged
 * per scaffold on load, so each input contributes coverage 0 or 1 at any point;
 * the combined event depth at a coordinate therefore equals the number of inputs
 * covering it. Each operation sorts the combined boundary events per scaffold,
 * then performs a linear coordinate sweep over that depth:
 * union keeps depth&gt;=1, intersection keeps depth==N, and subtraction keeps the
 * points where the first file covers but none of the others do.  This is correct
 * even when valid inputs are unsorted or self-overlapping. The current container
 * asserts that merged intervals have positive length; isolated zero-length input
 * intervals therefore fail that check (see the retained BedSet#001 note).
 *
 * @author UMP45
 * @contributor Brian Bushnell
 * @date June 19, 2026
 */
public class BedSet{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the command and closes a redirected diagnostic stream after normal completion. */
	public static void main(String[] args){
		Timer t=new Timer();
		BedSet x=new BedSet(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructs a BedSet instance from command-line arguments.
	 * Parses mode and comma-separated inputs; one input is sufficient. Mode keyword
	 * presence selects the operation regardless of its assigned value, with the last
	 * keyword taking precedence. Padding and multiallelic selection apply only to VCF.
	 * Also enables shared pigz/unpigz settings and configures their thread count.
	 * @param args Command-line arguments containing file paths and options
	 */
	public BedSet(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		int mode_=SUBTRACT;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("subtract") || a.equals("subtraction") || a.equals("minus") || a.equals("difference") || a.equals("dif") || a.equals("diff")){
				mode_=SUBTRACT;
			}else if(a.equals("union") || a.equals("plus") || a.equals("merge")){
				mode_=UNION;
			}else if(a.equals("intersection") || a.equals("intersect") || a.equals("shared")){
				mode_=INTERSECTION;
			}else if(a.equals("pad") || a.equals("padding")){
				pad=Integer.parseInt(b);
			}else if(a.equals("multiallelic") || a.equals("ma") || a.equals("multi")){
				multiallelicOnly=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		mode=mode_;
		{//Process parser fields
			in1=(parser.in1==null ? null : parser.in1.split(","));
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(in1==null || in1.length<1){throw new RuntimeException("Error - at least one input file is required (two or more for set operations).");}

		ffin1=new FileFormat[in1.length];
		for(int i=0; i<in1.length; i++){
			ffin1[i]=FileFormat.testInput(in1[i], FileFormat.TXT, null, true, false);
		}
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads every input, applies the configured set operation, writes the result
	 * BED (if out= was given), and prints coverage statistics. Input processing is
	 * sequential and fully materialized before output opens. With one input, every
	 * operation returns its merged intervals. The processed-line counter counts
	 * retained input intervals, not all physical rows; byte counts are not accumulated.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){

		final Bed[] beds=new Bed[ffin1.length];
		for(int i=0; i<beds.length; i++){
			if(verbose){outstream.println("Loading intervals from "+ffin1[i].name());}
			beds[i]=new Bed(ffin1[i], pad, multiallelicOnly);
			linesProcessed+=beds[i].rawIntervals;
		}

		if(verbose){outstream.println("Applying "+modeName(mode)+" to "+beds.length+" interval sets");}
		final Bed result=op(beds, mode);

		if(ffout1!=null){writeBed(result, ffout1);}

		//Statistics block (always computed; the investigation payload)
		final Bed unionAll=op(beds, UNION);
		final Bed interAll=op(beds, INTERSECTION);

		t.stop();
		//TODO: Possible bug [var2/BedSet#002] (LOW) - bytesProcessed is never incremented anywhere, so the
		//reported bytes/sec is always 0. Wire it from the input byte counts, or drop it from this call.
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Operation:         \t"+modeName(mode));
		outstream.println("Inputs:            \t"+beds.length);
		for(int i=0; i<beds.length; i++){
			outstream.println("Input "+i+" ("+beds[i].name+"):");
			outstream.println("  bp covered:      \t"+beds[i].bp);
			outstream.println("  intervals:       \t"+beds[i].intervals);
			outstream.println("  scaffolds:       \t"+beds[i].scaffolds());
		}
		outstream.println("Union bp:          \t"+unionAll.bp);
		outstream.println("Intersection bp:   \t"+interAll.bp);

		if(beds.length==2){
			final long shared=interAll.bp;
			final long uniqueA=beds[0].bp-shared;
			final long uniqueB=beds[1].bp-shared;
			//Runtime self-checks that the sweep algebra is internally consistent
			final Bed aMinusB=op(new Bed[]{beds[0], beds[1]}, SUBTRACT);
			final Bed bMinusA=op(new Bed[]{beds[1], beds[0]}, SUBTRACT);
			assert(shared>=0 && shared<=Tools.min(beds[0].bp, beds[1].bp)) : "shared="+shared;
			assert(uniqueA==aMinusB.bp) : "uniqueA "+uniqueA+" != A-B "+aMinusB.bp;
			assert(uniqueB==bMinusA.bp) : "uniqueB "+uniqueB+" != B-A "+bMinusA.bp;
			assert(unionAll.bp==beds[0].bp+beds[1].bp-shared) : "union "+unionAll.bp+" != "+(beds[0].bp+beds[1].bp-shared);
			outstream.println("Shared bp:         \t"+shared);
			outstream.println("Unique to input 0: \t"+uniqueA);
			outstream.println("Unique to input 1: \t"+uniqueB);
		}
		outstream.println("Output ("+modeName(mode)+") bp:");
		outstream.println("  bp covered:      \t"+result.bp);
		outstream.println("  intervals:       \t"+result.intervals);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Set Operations         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Applies a set operation across all inputs via a per-scaffold coordinate sweep.
	 * @param beds Loaded inputs; the first is "A" for subtraction
	 * @param opMode SUBTRACT, UNION, or INTERSECTION
	 * @return A new merged Bed holding the result intervals
	 */
	private Bed op(Bed[] beds, int opMode){
		final Bed out=new Bed(modeName(opMode));
		final HashSet<String> scafs=new HashSet<String>();
		for(Bed bed : beds){scafs.addAll(bed.starts.keySet());}
		final int nInputs=beds.length;
		for(String scaf : scafs){
			ArrayList<int[]> segs=sweepScaf(beds, scaf, opMode, nInputs);
			if(segs!=null && !segs.isEmpty()){out.putScaf(scaf, segs);}
		}
		out.finish();
		return out;
	}

	/**
	 * Sweeps one scaffold's combined boundary events and emits the output intervals
	 * that satisfy the operation's depth predicate.
	 * @param beds Inputs (beds[0] is A, the rest are subtracted/compared against)
	 * @param scaf Scaffold name
	 * @param opMode SUBTRACT, UNION, or INTERSECTION
	 * @param nInputs Number of inputs
	 * @return Merged output intervals on this scaffold, or null if none
	 */
	private static ArrayList<int[]> sweepScaf(Bed[] beds, String scaf, int opMode, int nInputs){
		//Each event is {coord, deltaA, deltaRest}: A is input 0, rest are inputs 1..N-1
		final ArrayList<int[]> events=new ArrayList<int[]>();
		for(int b=0; b<nInputs; b++){
			final int[] s=beds[b].starts.get(scaf), t=beds[b].stops.get(scaf);
			if(s==null){continue;}
			final boolean isA=(b==0);
			for(int i=0; i<s.length; i++){
				if(isA){
					events.add(new int[]{s[i], 1, 0});
					events.add(new int[]{t[i], -1, 0});
				}else{
					events.add(new int[]{s[i], 0, 1});
					events.add(new int[]{t[i], 0, -1});
				}
			}
		}
		if(events.isEmpty()){return null;}
		Collections.sort(events, EVENT_COMP);

		final ArrayList<int[]> out=new ArrayList<int[]>();
		int depthA=0, depthRest=0;
		int idx=0;
		final int n=events.size();
		int lastCoord=events.get(0)[0];
		while(idx<n){
			final int coord=events.get(idx)[0];
			//The segment [lastCoord, coord) carries the depth accumulated BEFORE this coordinate
			if(coord>lastCoord && qualifies(opMode, depthA, depthRest, nInputs)){
				out.add(new int[]{lastCoord, coord});
			}
			//Apply every event at this coordinate before moving on
			while(idx<n && events.get(idx)[0]==coord){
				depthA+=events.get(idx)[1];
				depthRest+=events.get(idx)[2];
				idx++;
			}
			assert(depthA>=0 && depthA<=1) : "depthA out of range (input 0 should be merged): "+depthA+" @"+scaf+":"+coord;
			assert(depthRest>=0 && depthRest<=nInputs-1) : "depthRest out of range: "+depthRest+" @"+scaf+":"+coord;
			lastCoord=coord;
		}
		assert(depthA==0 && depthRest==0) : "Unbalanced events on "+scaf+": "+depthA+","+depthRest;
		return out.isEmpty() ? null : mergeList(out);
	}

	/**
	 * The depth predicate for each operation.  Because inputs are merged, depthA is
	 * 0 or 1 and depthRest is the count of the remaining inputs (0..N-1) covering the point.
	 * @return true if the current point belongs in the output
	 */
	private static boolean qualifies(int opMode, int depthA, int depthRest, int nInputs){
		if(opMode==UNION){return depthA+depthRest>=1;}
		if(opMode==INTERSECTION){return depthA>=1 && depthRest>=nInputs-1;}
		//SUBTRACT: covered by A but by none of the rest
		return depthA>=1 && depthRest==0;
	}

	/** Merges a coordinate-sorted list, copying intervals and joining overlaps or adjacency.
	 * Does not modify the supplied list or its interval arrays; bounds are not validated. */
	private static ArrayList<int[]> mergeList(ArrayList<int[]> list){
		final ArrayList<int[]> out=new ArrayList<int[]>(list.size());
		int[] cur=null;
		for(int[] iv : list){
			if(cur==null){cur=new int[]{iv[0], iv[1]};}
			else if(iv[0]<=cur[1]){cur[1]=Tools.max(cur[1], iv[1]);}//Overlapping or abutting: extend
			else{out.add(cur); cur=new int[]{iv[0], iv[1]};}
		}
		if(cur!=null){out.add(cur);}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Writes a Bed's intervals as a 3-column BED file, scaffolds sorted by name.
	 * Waits for normal writer completion and accumulates its error flag. The reusable
	 * formatting buffer is copied by print before it is cleared.
	 * @param bed Result intervals to write
	 * @param ff Output file format
	 */
	private void writeBed(Bed bed, FileFormat ff){
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ArrayList<String> scafs=new ArrayList<String>(bed.starts.keySet());
		Collections.sort(scafs);
		final ByteBuilder bb=new ByteBuilder(33000);
		for(String scaf : scafs){
			final int[] s=bed.starts.get(scaf), t=bed.stops.get(scaf);
			for(int i=0; i<s.length; i++){
				bb.append(scaf).tab().append(s[i]).tab().append(t[i]).nl();
				if(bb.length>=32000){bsw.print(bb); bb.clear();}
			}
		}
		if(bb.length>0){bsw.print(bb); bb.clear();}
		errorState|=bsw.poisonAndWait();
	}

	/** Human-readable name for an operation mode. */
	private static String modeName(int m){
		return m==UNION ? "union" : m==INTERSECTION ? "intersection" : "subtract";
	}

	/*--------------------------------------------------------------*/
	/*----------------         Bed Container        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * A loaded BED: sorted, merged intervals per scaffold (0-based half-open),
	 * stored as parallel start/stop int arrays, plus cached coverage totals.
	 */
	private static class Bed{

		/** Creates an empty Bed (for holding set-operation results). */
		Bed(String name_){name=name_;}

		/** Loads BED intervals or padded VCF reference spans, then sorts, merges and counts.
		 * The filename-based parser choice is made by isVcf; pad and selection do not affect BED. */
		Bed(FileFormat ff, int pad, boolean multiallelicOnly){
			name=ff.name();
			final HashMap<String, ArrayList<int[]>> tmp=new HashMap<String, ArrayList<int[]>>();
			if(isVcf(ff)){loadVcf(ff, pad, multiallelicOnly, tmp);}
			else{loadBed(ff, tmp);}
			for(Map.Entry<String, ArrayList<int[]>> e : tmp.entrySet()){
				final ArrayList<int[]> list=e.getValue();
				Collections.sort(list, IV_COMP);
				putScaf(e.getKey(), mergeList(list));
			}
			finish();
		}

		/** Selects VCF when the lowercased filename contains .vcf anywhere; ignores ff's format flag. */
		private static boolean isVcf(FileFormat ff){
			return ff.name().toLowerCase().contains(".vcf");
		}

		/** Appends three-column BED intervals using exact scaffold names; ignores later columns.
		 * Skips empty/comment/header/short rows; parseable BED rows take precedence over
		 * reserved header names. Coordinates are untrimmed ints; start&lt;=stop is asserted,
		 * but nonnegative bounds are not checked. Reader errors reject the load. */
		private void loadBed(FileFormat ff, HashMap<String, ArrayList<int[]>> tmp){
			final ByteFile bf=ByteFile.makeByteFile(ff);
			long raw=0;
			Throwable failure=null;
			try{
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(line.length==0 || line[0]=='#'){continue;}//Comment line
					final String[] f=new String(line).split("\t");
					if(isHeader(line, f)){continue;}//VAR2-041: preserve scaffold names sharing a header prefix.
					if(f.length<3){continue;}//Not a valid BED record
					final int start=Integer.parseInt(f[1]);
					final int stop=Integer.parseInt(f[2]);
					assert(start<=stop) : "Invalid BED interval (start>stop): "+new String(line);
					ArrayList<int[]> list=tmp.get(f[0]);
					if(list==null){list=new ArrayList<int[]>(); tmp.put(f[0], list);}
					list.add(new int[]{start, stop});
					raw++;
				}
			}catch(RuntimeException | Error t){
				failure=t;
				throw t;
			}finally{closeInput(bf, failure);}
			rawIntervals=raw;
		}

		/**
		 * Reads a VCF and emits each variant's reference span as a BED interval, padded by
		 * 'pad' bp on each side so that deletions are fully covered with margin.  When
		 * multiallelicOnly is true, only sites whose first-sample GT carries an allele index
		 * &gt;=2 are emitted. GT is located through FORMAT; a missing key or sample value
		 * does not satisfy that filter. Uses POS and literal REF
		 * length (zero for '.'); INFO/END and symbolic allele semantics are not decoded.
		 * Clamps only the padded start to zero; padding affects VCF inputs only.
		 */
		private void loadVcf(FileFormat ff, int pad, boolean multiallelicOnly, HashMap<String, ArrayList<int[]>> tmp){
			final ByteFile bf=ByteFile.makeByteFile(ff);
			long raw=0;
			Throwable failure=null;
			try{
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					if(line.length==0 || line[0]=='#'){continue;}//Header/comment line
					final String[] f=new String(line).split("\t");
					if(f.length<8){continue;}//Not a complete variant record
					if(multiallelicOnly){
						if(f.length<10){continue;}//No genotype column to test
						final String gt=genotype(f[8], f[9]);//VAR2-043: DP and other fields are not allele indices.
						if(gt==null || !hasMultiallelicAllele(gt)){continue;}
					}
					final int pos=Integer.parseInt(f[1]);//1-based POS
					final int reflen=(f[3].equals(".") ? 0 : f[3].length());
					final int start=Tools.max(0, (pos-1)-pad);
					final int stop=(pos-1)+reflen+pad;
					ArrayList<int[]> list=tmp.get(f[0]);
					if(list==null){list=new ArrayList<int[]>(); tmp.put(f[0], list);}
					list.add(new int[]{start, stop});
					raw++;
				}
			}catch(RuntimeException | Error t){
				failure=t;
				throw t;
			}finally{closeInput(bf, failure);}
			rawIntervals=raw;
		}

		/** Recognizes whole header keywords while preserving tab-delimited integer BED rows. */
		private static boolean isHeader(final byte[] line, final String[] fields){
			assert(line.length>0) : "loadBed removes empty lines before header classification";
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

		/** Closes an owned input and rejects latched errors without replacing an earlier failure. */
		private static void closeInput(final ByteFile bf, final Throwable failure){
			//VAR2-042: callers invoke this from finally, before publishing loaded interval totals.
			try{
				if(bf.close()){throw new RuntimeException("Error reading or closing interval file: "+bf.name());}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}

		/** Returns the first sample's GT field, or null if the key/value is absent.
		 * Scans colon-delimited fields without allocating arrays for unrelated values. */
		private static String genotype(final String format, final String sample){
			assert(format!=null && sample!=null) : "loadVcf supplies retained FORMAT and first-sample columns";
			int keyStart=0, valueStart=0;
			while(keyStart<format.length()){
				int keyEnd=format.indexOf(':', keyStart), valueEnd=sample.indexOf(':', valueStart);
				if(keyEnd<0){keyEnd=format.length();}
				if(valueEnd<0){valueEnd=sample.length();}
				if(keyEnd-keyStart==2 && format.charAt(keyStart)=='G' && format.charAt(keyStart+1)=='T'){
					return sample.substring(valueStart, valueEnd);
				}
				if(keyEnd==format.length() || valueEnd==sample.length()){return null;}
				keyStart=keyEnd+1;
				valueStart=valueEnd+1;
			}
			return null;
		}

		/** Tests slash/pipe-delimited integer tokens for an index >=2; skips missing or malformed tokens.
		 * Does not check the site's ALT count or confirm that the caller supplied a GT field. */
		private static boolean hasMultiallelicAllele(String gt){
			final String[] alleles=gt.split("[/|]");
			for(String al : alleles){
				if(al.isEmpty() || al.equals(".")){continue;}
				try{if(Integer.parseInt(al)>=2){return true;}}catch(NumberFormatException e){}
			}
			return false;
		}

		/** Copies one scaffold's intervals into parallel arrays; replaces the key and leaves totals stale. */
		void putScaf(String scaf, ArrayList<int[]> merged){
			final int n=merged.size();
			final int[] s=new int[n], t=new int[n];
			for(int i=0; i<n; i++){s[i]=merged.get(i)[0]; t[i]=merged.get(i)[1];}
			starts.put(scaf, s);
			stops.put(scaf, t);
		}

		/** Recomputes total bp and interval counts from the stored intervals. */
		void finish(){
			long b=0, iv=0;
			for(Map.Entry<String, int[]> e : starts.entrySet()){
				final int[] s=e.getValue(), t=stops.get(e.getKey());
				for(int i=0; i<s.length; i++){
					//TODO: Possible bug [var2/BedSet#001] (MEDIUM, QUESTION-to-Brian) - loadBed admits start==stop
					//(assert start<=stop, a valid zero-length BED interval) but this rejects it (assert s[i]<t[i]),
					//so an ISOLATED zero-length interval crashes here under -ea. Crash-loud not silent-wrong, so
					//intent question: relax to s[i]<=t[i] (0-bp interval is harmless) OR reject at load with <stop.
					assert(s[i]<t[i]) : "Zero/negative-length interval "+s[i]+","+t[i]+" on "+e.getKey();
					assert(i==0 || s[i]>t[i-1]) : "Unmerged/unsorted intervals on "+e.getKey();
					b+=(t[i]-s[i]);
				}
				iv+=s.length;
			}
			bp=b; intervals=iv;
		}

		/** Number of distinct scaffolds with at least one interval. */
		int scaffolds(){return starts.size();}

		/** Sorted, merged interval starts (0-based) per scaffold. */
		final HashMap<String, int[]> starts=new HashMap<String, int[]>();
		/** Sorted, merged interval ends (0-based, exclusive) per scaffold. */
		final HashMap<String, int[]> stops=new HashMap<String, int[]>();
		/** Total base pairs covered. */
		long bp=0;
		/** Total merged intervals. */
		long intervals=0;
		/** Retained BED intervals or selected VCF spans before merging; operation results retain zero. */
		long rawIntervals=0;
		/** Display label (input file name or operation name). */
		String name="";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained input intervals, summed after each input is loaded. */
	private long linesProcessed=0;
	/** Currently never incremented; see the retained BedSet#002 diagnostic note. */
	private long bytesProcessed=0;

	/** Ordered comma-separated input names; input zero is the subtraction minuend. */
	private String in1[]=null;
	/** Optional result BED destination; null requests statistics only. */
	private String out1=null;

	/** Descriptors in input order, opened when process loads each Bed. */
	private final FileFormat ffin1[];
	/** Optional output descriptor, opened only after all inputs and the result are loaded. */
	private final FileFormat ffout1;

	/** Last selected operation keyword, or SUBTRACT by default. */
	public final int mode;
	/** Signed VCF span padding; input parsing does not validate its range. */
	private int pad=0;
	/** Whether the VCF loader applies its first-sample allele-index test. */
	private boolean multiallelicOnly=false;

	/** Diagnostic destination chosen by PreParser. */
	private PrintStream outstream=System.err;
	/** Enables input-loading and operation progress messages on the diagnostic stream. */
	public static boolean verbose=false;
	/** Accumulated normal writer-close error; reader errors are not folded into it. */
	public boolean errorState=false;
	/** Output file policy passed to FileFormat. */
	private boolean overwrite=true;
	/** Whether the output descriptor requests appending to an existing file. */
	private boolean append=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Constant for the subtraction operation mode (first file minus the rest). */
	public static final int SUBTRACT=0;
	/** Constant for the union operation mode. */
	public static final int UNION=1;
	/** Constant for the intersection operation mode. */
	public static final int INTERSECTION=2;

	/** Orders raw intervals by start, then end. */
	private static final Comparator<int[]> IV_COMP=new Comparator<int[]>(){
		@Override
		public int compare(int[] x, int[] y){
			return x[0]!=y[0] ? Integer.compare(x[0], y[0]) : Integer.compare(x[1], y[1]);
		}
	};

	/** Orders sweep events by coordinate. */
	private static final Comparator<int[]> EVENT_COMP=new Comparator<int[]>(){
		@Override
		public int compare(int[] x, int[] y){return Integer.compare(x[0], y[0]);}
	};

}
