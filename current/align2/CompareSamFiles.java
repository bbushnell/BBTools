package align2;

import java.io.File;
import java.util.Arrays;
import java.util.BitSet;

import dna.Data;
import fileIO.ByteFile;
import shared.KillSwitch;
import shared.LineParser1;
import shared.Parse;
import shared.PreParser;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/** Generate a file containing reads mapped correctly in one file and incorrectly in another file. */
public class CompareSamFiles {
	
	
	/**
	 * Program entry point that compares two SAM files and outputs reads with differing mapping status.
	 * Processes command-line arguments to specify input files, quality thresholds, and comparison parameters.
	 * Uses BitSet structures to track true positives and false positives across both input files.
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
//			}else if(a.equals("ssaha2") || a.equals("subtractleadingclip")){
//				SamLine.SUBTRACT_LEADING_SOFT_CLIP=Parse.parseBoolean(b);
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

		assert(in1!=null) : args[0]+".exists() ? "+new File(args[0]).exists();
//		assert(in2!=null) : args[1]+".exists() ? "+new File(args[1]).exists();
		
		if(reads<1){
//			assert(false) : "Number of expected reads was not specified.  Please add a parameter reads=<number> or disable assertions.";
			reads=100000;
			System.err.println("Warning - number of expected reads was not specified.");
		}

		ByteFile tf1=ByteFile.makeByteFile(in1, false);
		ByteFile tf2=null;
		if(in2!=null){tf2=ByteFile.makeByteFile(in2, false);}

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
				byte c=s[0];
				if(c!='@'/* && c!=' ' && c!='\t'*/){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.primary()){
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
				byte c=s[0];
				if(c!='@'/* && c!=' ' && c!='\t'*/){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.primary()){
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
				byte c=s[0];
				if(c!='@'/* && c!=' ' && c!='\t'*/){
					SamLine sl=new SamLine(lp.set(s));
//					assert(false) : s+", "+truePos1.cardinality()+", "+truePos2.cardinality()+", "+falsePos1.cardinality()+", "+falsePos2.cardinality()+", ";
					if(sl.primary()){
						Read r=sl.toRead(parsecustom);
						int id=(int)r.numericID;
						if(!added.get(id)){
//							if(truePos1.get(id)!=truePos2.get(id) || falsePos1.get(id)!=falsePos2.get(id)){
//								System.out.println(s);
//								added.set(id);
//							}
//							if(falsePos1.get(id) && truePos2.get(id)){
//								System.out.println(s);
//								added.set(id);
//							}
							if(falsePos1.get(id) && !falsePos2.get(id)){
								System.out.println(s);
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
				byte c=s[0];
				if(c!='@'/* && c!=' ' && c!='\t'*/){
					SamLine sl=new SamLine(lp.set(s));
					if(sl.primary()){
						Read r=sl.toRead(parsecustom);
						int id=(int)r.numericID;
						if(!added.get(id)){
//							if(truePos1.get(id)!=truePos2.get(id) || falsePos1.get(id)!=falsePos2.get(id)){
//								System.out.println(s);
//								added.set(id);
//							}
//							if(falsePos2.get(id) && truePos1.get(id)){
//								System.out.println(s);
//								added.set(id);
//							}
							if(falsePos2.get(id) && !falsePos1.get(id)){
								System.out.println(s);
								added.set(id);
							}
						}
					}
				}
			}
			tf.close();
		}
	}
	

	/**
	 * Calculates alignment statistics for a single read and updates global counters.
	 * Classifies reads as discarded, unmapped, ambiguous, or properly mapped based on quality scores.
	 * For custom parsing mode, evaluates alignment correctness using strict and loose criteria.
	 *
	 * @param r The read to evaluate
	 * @param sl Corresponding SAM line containing alignment details
	 */
	public static void calcStatistics1(final Read r, SamLine sl){
		
		int THRESH=0;
		primary++;
		
		if(r.discarded()/* || r.mapScore==0*/){
			discarded++;
			unmapped++;
		}else if(r.ambiguous()){
//			assert(r.mapped()) : "\n"+r+"\n"+sl+"\n";
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
						byte[] originalContig=sl.originalContig();
						if(BLASR){
							originalContig=(originalContig==null || Tools.indexOf(originalContig, (byte)'/')<0 ? originalContig :
								KillSwitch.copyOfRange(originalContig, 0, Tools.lastIndexOf(originalContig, (byte)'/')));
						}
						int cstart=sl.originalContigStart();

						boolean strict=isCorrectHit(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH, originalContig, sl.rname(), cstart);
						boolean loose=isCorrectHitLoose(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH+THRESH2, originalContig, sl.rname(), cstart);

						//				if(!strict){
						//					System.out.println(ss+", "+new String(originalContig)+", "+new String(sl.rname()));
						//					assert(false);
						//				}

						//				System.out.println("loose = "+loose+" for "+r.toText());

						if(loose){
							//					System.err.println("TPL\t"+trueChrom+", "+trueStrand+", "+trueStart+", "+trueStop+"\tvs\t"
							//							+ss.chrom+", "+ss.strand+", "+ss.start+", "+ss.stop);
							truePositiveLoose++;
						}else{
							//					System.err.println("FPL\t"+trueChrom+", "+trueStrand+", "+trueStart+", "+trueStop+"\tvs\t"
							//							+ss.chrom+", "+ss.strand+", "+ss.start+", "+ss.stop);
							falsePositiveLoose++;
						}

						if(strict){
							//					System.err.println("TPS\t"+trueStart+", "+trueStop+"\tvs\t"+ss.start+", "+ss.stop);
							truePositiveStrict++;
						}else{
							//					System.err.println("FPS\t"+trueStart+", "+trueStop+"\tvs\t"+ss.start+", "+ss.stop);
							falsePositiveStrict++;
						}
					}
				}
			}
		}
	}
	
	

	/**
	 * Determines the classification type for a read based on mapping quality and correctness.
	 * Returns integer codes: 0=discarded/unmapped, 1=ambiguous, 2=correct strict, 3=correct loose, 4=incorrect.
	 * Uses original site information when parsecustom is enabled to evaluate mapping accuracy.
	 *
	 * @param r The read to classify
	 * @param sl Corresponding SAM line with alignment information
	 * @return Integer type code indicating read classification
	 */
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
						byte[] originalContig=sl.originalContig();
						if(BLASR){
							originalContig=(originalContig==null || Tools.indexOf(originalContig, (byte)'/')<0 ? originalContig :
								KillSwitch.copyOfRange(originalContig, 0, Tools.lastIndexOf(originalContig, (byte)'/')));
						}
						int cstart=sl.originalContigStart();

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
	 * Compares strand, contig/chromosome, and position coordinates within specified threshold.
	 * Uses absolute position differences for both start and stop coordinates.
	 *
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

		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
		int cstop=cstart+trueStop-trueStart;
//		return (absdif(ss.start, trueStart)<=thresh && absdif(ss.stop, trueStop)<=thresh);
		return (absdif(ss.start, cstart)<=thresh && absdif(ss.stop, cstop)<=thresh);
	}
	
	
	/**
	 * Evaluates whether an alignment represents a correct hit using loose criteria.
	 * Similar to isCorrectHit but uses OR logic for position matching instead of AND.
	 * Allows hits where either start OR stop position is within threshold.
	 *
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

		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
		int cstop=cstart+trueStop-trueStart;
//		return (absdif(ss.start, trueStart)<=thresh || absdif(ss.stop, trueStop)<=thresh);
		return (absdif(ss.start, cstart)<=thresh || absdif(ss.stop, cstop)<=thresh);
	}
	
	/**
	 * Calculates absolute difference between two integers.
	 * @param a First integer
	 * @param b Second integer
	 * @return Absolute difference |a - b|
	 */
	private static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}

	/** Count of reads correctly mapped using strict criteria */
	public static int truePositiveStrict=0;
	/** Count of reads incorrectly mapped using strict criteria */
	public static int falsePositiveStrict=0;
	
	/** Count of reads correctly mapped using loose criteria */
	public static int truePositiveLoose=0;
	/** Count of reads incorrectly mapped using loose criteria */
	public static int falsePositiveLoose=0;

	/** Total count of reads with mapping assignments */
	public static int mapped=0;
	/** Count of mapped reads retained after quality filtering */
	public static int mappedRetained=0;
	/** Count of reads without mapping assignments */
	public static int unmapped=0;
	
	/** Count of reads discarded during processing */
	public static int discarded=0;
	/** Count of reads with ambiguous mapping assignments */
	public static int ambiguous=0;

	/** Total count of lines processed from input files */
	public static long lines=0;
	/** Count of primary alignment records processed */
	public static long primary=0;
	/** Count of secondary alignment records processed */
	public static long secondary=0;
	
	/** Minimum mapping quality threshold for read retention */
	public static int minQuality=3;

	/**
	 * Enable custom parsing to extract original site information from read names
	 */
	public static boolean parsecustom=true;
	/** Enable error message printing during processing */
	public static boolean printerr=false;

	/**
	 * Additional threshold value added to base threshold for loose matching criteria
	 */
	public static int THRESH2=20;
	/** Enable BLASR-specific processing for contig name parsing */
	public static boolean BLASR=false;
	
}
