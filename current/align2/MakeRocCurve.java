package align2;

import java.io.File;
import java.util.Arrays;
import java.util.BitSet;

import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Timer;
import shared.Tools;
import stream.CustomHeader;
import stream.Read;
import stream.SamLine;

/**
 * Emits cumulative mapping/accuracy percentages at observed quality-score thresholds.
 * This is a threshold table, not a conventional FPR/TPR curve with separate
 * negative-class denominators. All percentages use the expected read-end count.
 * Strict/loose truth means both/either endpoint of one read, not its two mates.
 * Quality bins are clamped to 0..999. Current SYN names supply optional truth;
 * parsecustom=f reports mapping totals without truth or numeric-ID deduplication.
 * Configuration and accumulation arrays are global; serialize calls. main resets
 * run state, while direct process calls accumulate into the current arrays.
 * @author Brian Bushnell
 * @date 2013
 */
public class MakeRocCurve{

	public static void main(String[] args){
		resetRunState();

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}

		Timer t=new Timer();
		String in=null;
		long reads=-1;

		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("in") || a.equals("in1")){
				in=b;
			}else if(a.equals("reads")){
				reads=Parse.parseKMG(b);
			}else if(a.equals("parsecustom")){
				parsecustom=Parse.parseBoolean(b);
			}else if(a.equals("blasr")){
				BLASR=Parse.parseBoolean(b);
			}else if(a.equals("bitset")){
				USE_BITSET=Parse.parseBoolean(b);
			}else if(a.equals("thresh")){
				THRESH2=Integer.parseInt(b);
			}else if(a.equals("allowspaceslash")){
				allowSpaceslash=Parse.parseBoolean(b);
			}else if(a.equals("outputerrors")){
			}else if(i==0 && args[i].indexOf('=')<0 && (a.startsWith("stdin") || new File(args[0]).exists())){
				in=args[0];
			}else if(i==1 && args[i].indexOf('=')<0 && Tools.isDigit(a.charAt(0))){
				reads=Parse.parseKMG(a);
			}
		}

		if(in==null){throw new IllegalArgumentException("MakeRocCurve requires in=<SAM input>");}
		if(reads<1){throw new IllegalArgumentException("A positive reads=<expected read-end count> is required for ROC percentages");}
		if(USE_BITSET && parsecustom){
			int x=400000;
			if(reads>0 && reads<=Integer.MAX_VALUE){x=(int)reads;}
			try{
				seen=new BitSet(x);
			}catch(Exception e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.out.println("Did not have enough memory to allocate bitset; duplicate mappings will not be detected.");
			}
		}

		process(in);

		System.out.println("ROC Curve for "+in);
		System.out.println(header());
		gradeList(reads);
		t.stop();
		System.err.println("Time: \t"+t);

	}

	/**
	 * Processes a SAM file to collect alignment statistics for ROC analysis.
	 * Reads each SAM line, converts to Read objects, and calculates statistics for primary alignments while avoiding duplicate counting.
	 * @param samfile Path to the SAM format alignment file
	 */
	public static void process(String samfile){
		ByteFile tf=ByteFile.makeByteFile(samfile, false);
		try{
			LineParser1 lp=new LineParser1('\t');
			for(byte[] s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.length==0){continue;}
				byte c=s[0];
				if(c!='@'/* && c!=' ' && c!='\t'*/){
					SamLine sl=new SamLine(lp.set(s));
					//TODO: Probable bug - narrowing/shifting a long ID can alias reads or produce a negative BitSet index.
					final int id=(parsecustom && seen!=null ? ((((int)sl.parseNumericId())<<1)|sl.pairnum()) : -1);
					assert(sl!=null);
					Read r=sl.toRead(parsecustom);
					if(r!=null){
						r.samline=sl;
						if(sl.nonSecondary() && (!parsecustom || seen==null || !seen.get(id))){
							if(parsecustom && seen!=null){seen.set(id);}
							calcStatistics1(r, sl);
						}
					}else{
						assert(false) : "'"+"'";
						System.err.println("Bad read from line '"+s+"'");
					}
				}
			}
		}finally{
			if(tf.close()){throw new RuntimeException("I/O failure while reading ROC input: "+samfile);}
		}
	}

	public static String header(){
		return "minScore\tmapped\tretained\ttruePositiveStrict\tfalsePositiveStrict\ttruePositiveLoose" +
				"\tfalsePositiveLoose\tfalseNegative\tdiscarded\tambiguous";
	}

	/**
	 * Generates and prints the ROC curve data by iterating through quality scores from highest to lowest.
	 * Calculates cumulative statistics including true/false positives and outputs mapped/retained/ambiguous percentages.
	 * @param reads Total number of reads for percentage calculations
	 */
	public static void gradeList(long reads){
		if(reads<1){throw new IllegalArgumentException("ROC percentage denominator must be positive: "+reads);}

		int truePositiveStrict=0;
		int falsePositiveStrict=0;

		int truePositiveLoose=0;
		int falsePositiveLoose=0;

		int mapped=0;
		int mappedRetained=0;
		int unmapped=0;

		int discarded=0;
		int ambiguous=0;

		int primary=0;

		for(int q=truePositiveStrictA.length-1; q>=0; q--){
			//TODO: Probable bug - an unmapped ambiguous read increments neither mappedA nor unmappedA;
			//a bin containing only those reads is skipped, dropping its ambiguity count.
			if(mappedA[q]>0 || unmappedA[q]>0){
				truePositiveStrict+=truePositiveStrictA[q];
				falsePositiveStrict+=falsePositiveStrictA[q];
				truePositiveLoose+=truePositiveLooseA[q];
				falsePositiveLoose+=falsePositiveLooseA[q];
				mapped+=mappedA[q];
				mappedRetained+=mappedRetainedA[q];
				unmapped+=unmappedA[q];
				discarded+=discardedA[q];
				ambiguous+=ambiguousA[q];
				primary+=primaryA[q];

				double tmult=100d/reads;

				double mappedB=mapped*tmult;
				double retainedB=mappedRetained*tmult;
				double truePositiveStrictB=truePositiveStrict*tmult;
				double falsePositiveStrictB=falsePositiveStrict*tmult;
				double truePositiveLooseB=truePositiveLoose*tmult;
				double falsePositiveLooseB=falsePositiveLoose*tmult;
				double falseNegativeB=(reads-mapped)*tmult;
				double discardedB=discarded*tmult;
				double ambiguousB=ambiguous*tmult;

				StringBuilder sb=new StringBuilder();
				sb.append(q);
				sb.append('\t');
				sb.append(Tools.format("%.4f", mappedB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", retainedB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", truePositiveStrictB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", falsePositiveStrictB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", truePositiveLooseB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", falsePositiveLooseB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", falseNegativeB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", discardedB));
				sb.append('\t');
				sb.append(Tools.format("%.4f", ambiguousB));

				System.out.println(sb);
			}else{
				assert(truePositiveStrictA[q]==0) : q;
				assert(falsePositiveStrictA[q]==0) : q;
				assert(truePositiveLooseA[q]==0) : q;
				assert(falsePositiveLooseA[q]==0) : q;
			}

		}
	}

	public static void calcStatistics1(final Read r, SamLine sl){

		int q=r.mapScore;

		int THRESH=0;
		//[align2/MakeRocCurve#002 FIXED 2026-07-03] primaryA[q]++ was BEFORE the clamp below, so it indexed with the raw
		//mapScore -> AIOOBE when mapScore<0 or >=discardedA.length(1000), reachable for long/high-scoring reads (samtoroc.sh).
		//Every other array access uses the clamped q; moved the increment after the clamp so primaryA is bucketed consistently.
		if(q<0){q=0;}
		if(q>=discardedA.length){q=discardedA.length-1;}
		primaryA[q]++;

		if(r.discarded()/* || r.mapScore==0*/){
			discardedA[q]++;
			unmappedA[q]++;
		}else if(r.ambiguous()){
			if(r.mapped()){mappedA[q]++;}
			ambiguousA[q]++;
		}else if(r.mapScore<1){
			unmappedA[q]++;
		}else if(!r.mapped()){
			unmappedA[q]++;
		}
		else{

			mappedA[q]++;
			mappedRetainedA[q]++;

			if(parsecustom){
				CustomHeader h=new CustomHeader(sl.qname, sl.pairnum());
				boolean strict=isCorrectHit(sl, h);
				boolean loose=isCorrectHitLoose(sl, h);

				if(loose){
					truePositiveLooseA[q]++;
				}else{
					falsePositiveLooseA[q]++;
				}

				if(strict){
					truePositiveStrictA[q]++;
				}else{
					falsePositiveStrictA[q]++;
				}
			}
		}
	}

	/**
	 * Determines if an alignment is a strict true positive by comparing the mapped position to the true position in custom headers.
	 * Requires exact match of reference name, strand, start, and stop positions.
	 * @param sl SAM line with alignment information
	 * @param h Custom header containing true position information
	 * @return true if alignment exactly matches the true position
	 */
	public static boolean isCorrectHit(SamLine sl, CustomHeader h){
		if(!sl.mapped()){return false;}
		if(h.strand!=sl.strand()){return false;}
		int start=sl.start(true, true);
		int stop=sl.stop(start, true, true);
		if(h.start!=start){return false;}
		if(h.stop!=stop){return false;}
		if(!h.rname.equals(sl.rnameS())){return false;}
		return true;
	}

	/**
	 * Determines if an alignment is a loose true positive by allowing positional tolerance defined by THRESH2.
	 * More permissive than strict matching for evaluating alignment accuracy while still requiring the correct reference and strand.
	 * @param sl SAM line with alignment information
	 * @param h Custom header containing true position information
	 * @return true if alignment is within acceptable distance of true position
	 */
	public static boolean isCorrectHitLoose(SamLine sl, CustomHeader h){
		if(!sl.mapped()){return false;}
		if(h.strand!=sl.strand()){return false;}
		int start=sl.start(true, true);
		int stop=sl.stop(start, true, true);
		if(!h.rname.equals(sl.rnameS())){return false;}

		//[align2/MakeRocCurve#001 FIXED 2026-07-03] Removed two exact-match early returns (if(h.start!=start)return false; and
		//if(h.stop!=stop)return false;) that ran BEFORE the tolerance check below - they made this "loose" check require an EXACT
		//start AND stop, i.e. identical to isCorrectHit (strict), so the loose ROC columns silently equaled the strict ones.
		//The commented-out original isCorrectHitLoose went straight to the absdif tolerance; restored that. THRESH2=20 default.
		return(absdif(h.start, start)<=THRESH2 || absdif(h.stop, stop)<=THRESH2);
	}

	private static final long absdif(int a, int b){
		return a>b ? (long)a-b : (long)b-a;
	}

	/** Clears accumulated results and duplicate tracking without changing caller configuration. */
	private static void resetRunState(){
		for(int[] counts : new int[][]{truePositiveStrictA, falsePositiveStrictA,
				truePositiveLooseA, falsePositiveLooseA, mappedA, mappedRetainedA,
				unmappedA, discardedA, ambiguousA, primaryA}){
			Arrays.fill(counts, 0);
		}
		seen=null;
	}

	//TODO: Probable bug - per-bin and cumulative int counts can overflow on huge read sets.
	public static int truePositiveStrictA[]=new int[1000];
	public static int falsePositiveStrictA[]=new int[1000];

	public static int truePositiveLooseA[]=new int[1000];
	public static int falsePositiveLooseA[]=new int[1000];

	public static int mappedA[]=new int[1000];
	public static int mappedRetainedA[]=new int[1000];
	public static int unmappedA[]=new int[1000];

	public static int discardedA[]=new int[1000];
	public static int ambiguousA[]=new int[1000];

	public static int primaryA[]=new int[1000];

	public static boolean parsecustom=true;

	public static int THRESH2=20;
	public static boolean BLASR=false;
	public static boolean USE_BITSET=true;
	public static BitSet seen=null;
	public static boolean allowSpaceslash=true;

}
