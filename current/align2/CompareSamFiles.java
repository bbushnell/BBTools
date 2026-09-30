package align2;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.BitSet;

import dna.Data;
import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import stream.CustomHeader;
import stream.SamLine;
import stream.SiteScore;

/**
 * Compares strict-error membership for synthetic-read SAM files.
 * Emits headerless SAM records from either input when an accepted mapping is
 * not strictly correct there and the same numeric ID has no accepted strict
 * error in the other input. Absence, low quality or unmapped status in the other
 * file qualifies; correctness in the other file is not required. With one input,
 * emits its strict errors. Loose correctness is classified but also counts as a
 * strict error. Current SYN truth names are decoded through CustomHeader.
 * Inputs are read twice and must be reopenable; this is not a stdin pipeline.
 * State/configuration counters are process-global and not reset between calls.
 * @author Brian Bushnell
 * @date 2013
 */
public class CompareSamFiles{

	/**
	 * Program entry point that compares two SAM files and outputs reads with differing mapping status.
	 * Processes command-line arguments to specify input files, quality thresholds, and comparison parameters, and uses BitSets to track true/false positives.
	 * @param args Command-line arguments including in1, in2, reads, thresh, quality, parsecustom, blasr
	 */
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}

		String in1=null;
		String in2=null;
		long reads=-1;

		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("path") || a.equals("root")){
				Data.setPath(b);
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("parsecustom")){
				parsecustom=Parse.parseBoolean(b);
			}else if(a.equals("thresh")){
				THRESH2=Integer.parseInt(b);
			}else if(a.equals("printerr")){
				printerr=Parse.parseBoolean(b);
			}else if(a.equals("blasr")){
				BLASR=Parse.parseBoolean(b);
			}else if(a.equals("q") || a.equals("quality") || a.startsWith("minq")){
				minQuality=Integer.parseInt(b);
			}else if(in1==null && i==0 && args[i].indexOf('=')<0 && (a.startsWith("stdin") || new File(args[i]).exists())){
				in1=args[i];
			}else if(in2==null && i==1 && args[i].indexOf('=')<0 && (a.startsWith("stdin") || new File(args[i]).exists())){
				in2=args[i];
			}else if(a.equals("reads")){
				reads=Parse.parseKMG(b);
			}else if(i==2 && args[i].indexOf('=')<0 && Tools.isDigit(a.charAt(0))){
				reads=Parse.parseKMG(a);
			}
		}

		assert(in1!=null) : "CompareSamFiles requires in1=<synthetic-read SAM>";

		if(reads<1){
			reads=100000;
			System.err.println("Warning - number of expected reads was not specified.");
		}

		//TODO: Probable bug - an exception during parsing bypasses close; both input lifetimes need finally cleanup.
		ByteFile tf1=ByteFile.makeByteFile(in1, false);
		//TODO: Probable bug - stdin is accepted above but the output pass requires reset()/a second read.
		ByteFile tf2=null;
		if(in2!=null){tf2=ByteFile.makeByteFile(in2, false);}

		//TODO: Probable bug - long read counts/IDs are narrowed to int; large IDs can alias or become negative.
		//Paired ends sharing a numericID are also collapsed, so this is not an end-specific comparison.
		BitSet truePos1=new BitSet((int)reads);
		BitSet falsePos1=new BitSet((int)reads);
		BitSet truePos2=new BitSet((int)reads);
		BitSet falsePos2=new BitSet((int)reads);

		byte[] s=null;

		ByteFile tf;
		LineParser1 lp=new LineParser1('\t');
		{
			tf=tf1;
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.length==0){continue;}
				byte c=s[0];
				if(c!='@'){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.nonSecondary()){
						Read r=sl.toRead(parsecustom);
						if(parsecustom && r.originalSite==null){
							assert(false);
							System.err.println("Turned off custom parsing.");
							parsecustom=false;
						}
						//System.out.println(r);
						int type=type(r, sl);
						int id=(int)r.numericID;
						if(type==2){truePos1.set(id);}
						else if(type>2){falsePos1.set(id);}
					}
				}
			}
			tf.close();
		}
		if(tf2!=null){
			tf=tf2;
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.length==0){continue;}
				byte c=s[0];
				if(c!='@'){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.nonSecondary()){
						Read r=sl.toRead(parsecustom);
						if(parsecustom && r.originalSite==null){
							assert(false);
							System.err.println("Turned off custom parsing.");
							parsecustom=false;
						}
						//System.out.println(r);
						int type=type(r, sl);
						int id=(int)r.numericID;
						if(type==2){truePos2.set(id);}
						else if(type>2){falsePos2.set(id);}
					}
				}
			}
			tf.close();
		}

		BitSet added=new BitSet((int)reads);
		{
			tf=tf1;
			tf.reset();
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.length==0){continue;}
				byte c=s[0];
				if(c!='@'){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.nonSecondary()){
						Read r=sl.toRead(parsecustom);
						int id=(int)r.numericID;
						if(!added.get(id)){
							if(falsePos1.get(id) && !falsePos2.get(id)){
								System.out.write(s, 0, s.length); System.out.println();
								added.set(id);
							}
						}
					}
				}
			}
			tf.close();
		}
		if(tf2!=null){
			tf=tf2;
			tf.reset();
			for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(s.length==0){continue;}
				byte c=s[0];
				if(c!='@'){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.nonSecondary()){
						Read r=sl.toRead(parsecustom);
						int id=(int)r.numericID;
						if(!added.get(id)){
							if(falsePos2.get(id) && !falsePos1.get(id)){
								System.out.write(s, 0, s.length); System.out.println();
								added.set(id);
							}
						}
					}
				}
			}
			tf.close();
		}
	}

	/** Adds one read to legacy diagnostic counters; main uses type() instead. */
	public static void calcStatistics1(final Read r, SamLine sl){

		int THRESH=0;
		primary++;

		if(r.discarded()/* || r.mapScore==0*/){
			discarded++;
			unmapped++;
		}else if(r.ambiguous()){
			if(r.mapped()){mapped++;}
			ambiguous++;
		}else if(r.mapScore<1){
			unmapped++;
		}else if(r.mapScore<=minQuality){
			if(r.mapped()){mapped++;}
			ambiguous++;
		}else{
			if(!r.mapped()){
				unmapped++;
			}else{
				mapped++;
				mappedRetained++;

				if(parsecustom){
					SiteScore os=r.originalSite;
					assert(os!=null);
					if(os!=null){
						int trueChrom=os.chrom;
						byte trueStrand=os.strand;
						int trueStart=os.start;
						int trueStop=os.stop;
						SiteScore ss=new SiteScore(r.chrom, r.strand(), r.start, r.stop, 0, 0);
						CustomHeader truth=new CustomHeader(sl.qname, sl.pairnum());
						if(truth.rname==null){throw new IllegalArgumentException("Missing synthetic reference name: "+sl.qname);}
						byte[] originalContig=truth.rname.getBytes(StandardCharsets.UTF_8);
						if(BLASR){
							originalContig=(originalContig==null || Tools.indexOf(originalContig, (byte)'/')<0 ? originalContig :
								KillSwitch.copyOfRange(originalContig, 0, Tools.lastIndexOf(originalContig, (byte)'/')));
						}
						int cstart=truth.start;

						boolean strict=isCorrectHit(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH, originalContig, sl.rname(), cstart);
						boolean loose=isCorrectHitLoose(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH+THRESH2, originalContig, sl.rname(), cstart);

						if(loose){
							truePositiveLoose++;
						}else{
							falsePositiveLoose++;
						}

						if(strict){
							truePositiveStrict++;
						}else{
							falsePositiveStrict++;
						}
					}
				}
			}
		}
	}

	/** Returns 0=rejected/unmapped/ungraded, 1=ambiguous, 2=strict, 3=loose-only, 4=incorrect; increments primary. */
	public static int type(final Read r, SamLine sl){

		int THRESH=0;
		primary++;

		if(r.discarded()/* || r.mapScore==0*/){
			return 0;
		}else if(r.ambiguous()){
			return 1;
		}else if(r.mapScore<1){
			return 0;
		}else if(r.mapScore<=minQuality){
			return 1;
		}else{
			if(!r.mapped()){
				return 0;
			}else{

				if(parsecustom){
					SiteScore os=r.originalSite;
					assert(os!=null);
					if(os!=null){
						int trueChrom=os.chrom;
						byte trueStrand=os.strand;
						int trueStart=os.start;
						int trueStop=os.stop;
						SiteScore ss=new SiteScore(r.chrom, r.strand(), r.start, r.stop, 0, 0);
						CustomHeader truth=new CustomHeader(sl.qname, sl.pairnum());
						if(truth.rname==null){throw new IllegalArgumentException("Missing synthetic reference name: "+sl.qname);}
						byte[] originalContig=truth.rname.getBytes(StandardCharsets.UTF_8);
						if(BLASR){
							originalContig=(originalContig==null || Tools.indexOf(originalContig, (byte)'/')<0 ? originalContig :
								KillSwitch.copyOfRange(originalContig, 0, Tools.lastIndexOf(originalContig, (byte)'/')));
						}
						int cstart=truth.start;

						boolean strict=isCorrectHit(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH, originalContig, sl.rname(), cstart);
						boolean loose=isCorrectHitLoose(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH+THRESH2, originalContig, sl.rname(), cstart);

						if(strict){return 2;}
						if(loose){return 3;}
						return 4;
					}
				}
			}
		}
		return 0;
	}

	/**
	 * Evaluates whether an alignment represents a correct hit using strict criteria.
	 * Compares strand, contig/chromosome, and position coordinates within a threshold on contig-relative start and stop positions.
	 * @param ss The alignment site score to evaluate
	 * @param trueChrom True chromosome number
	 * @param trueStrand True strand orientation
	 * @param trueStart True start position
	 * @param trueStop True stop position
	 * @param thresh Maximum allowed position deviation
	 * @param originalContig Original contig name from alignment
	 * @param contig Current contig name
	 * @param cstart Contig-relative start position
	 * @return true if alignment meets strict correctness criteria
	 */
	public static boolean isCorrectHit(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh,
			byte[] originalContig, byte[] contig, int cstart){
		if(ss.strand!=trueStrand){return false;}
		if(originalContig!=null){
			if(!Arrays.equals(originalContig, contig)){return false;}
		}else{
			if(ss.chrom!=trueChrom){return false;}
		}

		assert(ss.stop>=ss.start) : "Inclusive mapped endpoints must be ordered: "+ss.toText();
		assert(trueStop>=trueStart) : "Inclusive truth endpoints must be ordered: "+trueStart+", "+trueStop;
		long cstop=(long)cstart+trueStop-trueStart;
		return (absdif(ss.start, cstart)<=thresh && absdif(ss.stop, cstop)<=thresh);
	}

	/**
	 * Evaluates whether an alignment represents a correct hit using loose criteria.
	 * Similar to isCorrectHit but uses OR logic for contig-relative start/stop positions instead of AND.
	 * @param ss The alignment site score to evaluate
	 * @param trueChrom True chromosome number
	 * @param trueStrand True strand orientation
	 * @param trueStart True start position
	 * @param trueStop True stop position
	 * @param thresh Maximum allowed position deviation
	 * @param originalContig Original contig name from alignment
	 * @param contig Current contig name
	 * @param cstart Contig-relative start position
	 * @return true if alignment meets loose correctness criteria
	 */
	public static boolean isCorrectHitLoose(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh,
			byte[] originalContig, byte[] contig, int cstart){
		if(ss.strand!=trueStrand){return false;}
		if(originalContig!=null){
			if(!Arrays.equals(originalContig, contig)){return false;}
		}else{
			if(ss.chrom!=trueChrom){return false;}
		}

		assert(ss.stop>=ss.start) : "Inclusive mapped endpoints must be ordered: "+ss.toText();
		assert(trueStop>=trueStart) : "Inclusive truth endpoints must be ordered: "+trueStart+", "+trueStop;
		long cstop=(long)cstart+trueStop-trueStart;
		return (absdif(ss.start, cstart)<=thresh || absdif(ss.stop, cstop)<=thresh);
	}

	private static final long absdif(long a, long b){
		return a>b ? a-b : b-a;
	}

	public static int truePositiveStrict=0;
	public static int falsePositiveStrict=0;

	public static int truePositiveLoose=0;
	public static int falsePositiveLoose=0;

	public static int mapped=0;
	public static int mappedRetained=0;
	public static int unmapped=0;

	public static int discarded=0;
	public static int ambiguous=0;

	public static long lines=0;
	public static long primary=0;
	public static long secondary=0;

	public static int minQuality=3;

	public static boolean parsecustom=true;
	public static boolean printerr=false;

	public static int THRESH2=20;
	public static boolean BLASR=false;

}
