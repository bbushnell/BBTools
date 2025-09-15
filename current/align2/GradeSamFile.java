package align2;

import java.io.File;
import java.util.BitSet;

import fileIO.FileFormat;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Tools;
import stream.CustomHeader;
import stream.Read;
import stream.SamLine;

public class GradeSamFile {
	
	
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		String in=null, outl=null, outs=null;
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
			}else if(a.equals("bitset")){
				USE_BITSET=Parse.parseBoolean(b);
			}else if(a.equals("outloose") || a.equals("outl")){
				outl=b;
			}else if(a.equals("outstrict") || a.equals("outs")){
				outs=b;
			}else if(i==0 && args[i].indexOf('=')<0 && (a.startsWith("stdin") || new File(args[0]).exists())){
				in=args[0];
			}else if(i==1 && args[i].indexOf('=')<0 && Tools.isDigit(a.charAt(0))){
				reads=Parse.parseKMG(a);
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		if(outl!=null){
			ffLoose=FileFormat.testOutput(outl, FileFormat.SAM, null, false, true, false, false);
			tswLoose=new TextStreamWriter(ffLoose);
			tswLoose.start();
		}
		
		if(outs!=null){
			ffStrict=FileFormat.testOutput(outs, FileFormat.SAM, null, false, true, false, false);
			tswStrict=new TextStreamWriter(ffStrict);
			tswStrict.start();
		}
		
		if(USE_BITSET){
			int x=400000;
			if(reads>0 && reads<=Integer.MAX_VALUE){x=(int)reads;}
			try {
				seen=new BitSet(x);
			} catch (Throwable e) {
				seen=null;
				e.printStackTrace();
				System.out.println("Did not have enough memory to allocate bitset; duplicate mappings will not be detected.");
			}
		}
		
		assert(in!=null) : args[0]+".exists() ? "+new File(args[0]).exists();
		
		if(reads<1){
			assert(false) : "Number of expected reads was not specified.  Please add a parameter reads=<number> or disable assertions.";
			System.err.println("Warning - number of expected reads was not specified.");
		}
		
		TextFile tf=new TextFile(in, false);
		
		String s=null;
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			char c=s.charAt(0);
//			System.out.println(s);
			if(c!='@'/* && c!=' ' && c!='\t'*/){
				SamLine sl=new SamLine(s);
				lines++;
				
				int id=(parsecustom && seen!=null ? ((((int)sl.parseNumericId())<<1)|sl.pairnum()) : (int)lines);
//				System.out.println(sl.parseNumericId()+", "+sl.pairnum()+", "+id+"");
//				if(id%500==10){assert(false);}
				if(sl.primary() && (!parsecustom || seen==null || !seen.get(id))){
					Read r=sl.toRead(parsecustom);
					if(seen!=null){seen.set(id);}
					if(parsecustom && r.originalSite==null){
						assert(false);
						System.err.println("Turned off custom parsing.");
						parsecustom=false;
					}
					//System.out.println(r);
					calcStatistics1(r, sl);
				}else{
					secondary++;
				}
			}
		}

		if(tswLoose!=null){tswLoose.poisonAndWait();}
		if(tswStrict!=null){tswStrict.poisonAndWait();}
		
		if(reads<-1){reads=primary;}
		
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
		
		System.out.println();
		System.out.println("Mapping Statistics for "+args[0]+":");
		System.out.println("primary alignments:    \t"+primary+" found of "+reads+" expected");
		System.out.println("secondary alignments:  \t"+secondary+" found");
		System.out.println(Tools.format("mapped:                \t"+(mappedB<10?" ":"")+"%.3f", mappedB)+"%");
		System.out.println(Tools.format("retained:              \t"+(retainedB<10?" ":"")+"%.3f", retainedB)+"%");
		System.out.println(Tools.format("discarded:             \t"+(discardedB<10?" ":"")+"%.3f", discardedB)+"%");
		System.out.println(Tools.format("ambiguous:             \t"+(ambiguousB<10?" ":"")+"%.3f", ambiguousB)+"%");
		if(parsecustom){
			System.out.println();
			System.out.println("Strict correctness (both ends exactly correct):");
			System.out.println(Tools.format("true positive:         \t"+(truePositiveStrictB<10?" ":"")+"%.3f", truePositiveStrictB)+"%");
			System.out.println(Tools.format("false positive:        \t"+(falsePositiveStrictB<10?" ":"")+"%.3f", falsePositiveStrictB)+"%");
			System.out.println();
			System.out.println("Loose correctness (one end approximately correct):");
			System.out.println(Tools.format("true positive:         \t"+(truePositiveLooseB<10?" ":"")+"%.3f", truePositiveLooseB)+"%");
			System.out.println(Tools.format("false positive:        \t"+(falsePositiveLooseB<10?" ":"")+"%.3f", falsePositiveLooseB)+"%");
		}
		System.out.println();
		System.out.println(Tools.format("false negative:        \t"+(falseNegativeB<10?" ":"")+"%.3f", falseNegativeB)+"%");
		
		if(printerr){
			System.err.println();
			System.err.println("Mapping Statistics for "+args[0]+":");
			System.err.println("primary alignments:    \t"+primary+" found of "+reads+" expected");
			System.err.println("secondary alignments:  \t"+secondary+" found");
			System.err.println(Tools.format("mapped:                \t"+(mappedB<10?" ":"")+"%.3f", mappedB)+"%");
			System.err.println(Tools.format("retained:              \t"+(retainedB<10?" ":"")+"%.3f", retainedB)+"%");
			System.err.println(Tools.format("discarded:             \t"+(discardedB<10?" ":"")+"%.3f", discardedB)+"%");
			System.err.println(Tools.format("ambiguous:             \t"+(ambiguousB<10?" ":"")+"%.3f", ambiguousB)+"%");
			if(parsecustom){
				System.err.println();
				System.err.println("Strict correctness (both ends exactly correct):");
				System.err.println(Tools.format("true positive:         \t"+(truePositiveStrictB<10?" ":"")+"%.3f", truePositiveStrictB)+"%");
				System.err.println(Tools.format("false positive:        \t"+(falsePositiveStrictB<10?" ":"")+"%.3f", falsePositiveStrictB)+"%");
				System.err.println();
				System.err.println("Loose correctness (one end approximately correct):");
				System.err.println(Tools.format("true positive:         \t"+(truePositiveLooseB<10?" ":"")+"%.3f", truePositiveLooseB)+"%");
				System.err.println(Tools.format("false positive:        \t"+(falsePositiveLooseB<10?" ":"")+"%.3f", falsePositiveLooseB)+"%");
			}
			System.err.println();
			System.err.println(Tools.format("false negative:        \t"+(falseNegativeB<10?" ":"")+"%.3f", falseNegativeB)+"%");
		}
		
		
	}
	
	
	public static void calcStatistics1(final Read r, SamLine sl){
		
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
					CustomHeader h=new CustomHeader(sl.qname, sl.pairnum());
					boolean strict=isCorrectHit(sl, h);
					boolean loose=isCorrectHitLoose(sl, h);

//					System.err.println(sl.strand()+", "+sl.start(true, true)+", "+sl.stop(sl.start(true, true), true, true)+", "+sl.rnameS());
//					System.err.println(h.strand+", "+h.start+", "+h.stop+", "+h.rname);
//					System.err.println(strict+", "+loose);
//					System.err.println();
//					assert(false);
					
//					SiteScore os=r.originalSite;
//					System.out.println("A1: "+os);
//							assert(os!=null);
//					final int trueChrom=os.chrom;
//					final byte trueStrand=os.strand;
//					final int trueStart=os.start;
//					final int trueStop=os.stop;
//					//						System.err.println();
//					//						System.err.println(sl);
//					//						System.err.println();
//					//						System.err.println(r);
//					//						System.err.println();
//					SiteScore ss=new SiteScore(r.chrom, r.strand(), r.start, r.stop, 0, 0);
//					byte[] originalContig=sl.originalContig();
//					if(BLASR){
//						originalContig=(originalContig==null || Tools.indexOf(originalContig, (byte)'/')<0 ? originalContig :
//							KillSwitch.copyOfRange(originalContig, 0, Tools.lastIndexOf(originalContig, (byte)'/')));
//					}
//					int cstart=sl.originalContigStart();
//
//					//						System.out.println("A2: "+trueStart+", "+cstart);
//					boolean strict=isCorrectHit(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH, originalContig, sl.rname(), cstart, r);
//					boolean loose=isCorrectHitLoose(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH+THRESH2, originalContig, sl.rname(), cstart);


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
						if(tswLoose!=null){
							if(ffLoose.samOrBam()){
								tswLoose.println(sl.toText());
							}else{
								tswLoose.println(r);
							}
						}
					}

					if(strict){
						//					System.err.println("TPS\t"+trueStart+", "+trueStop+"\tvs\t"+ss.start+", "+ss.stop);
						truePositiveStrict++;
					}else{
						//					System.err.println("FPS\t"+trueStart+", "+trueStop+"\tvs\t"+ss.start+", "+ss.stop);
						falsePositiveStrict++;
						if(tswStrict!=null){
							if(ffStrict.samOrBam()){
								tswStrict.println(sl.toText());
							}else{
								tswStrict.println(r);
							}
						}
					}
				}
			}
		}
		
	}
	
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
	
	public static boolean isCorrectHitLoose(SamLine sl, CustomHeader h){
		if(!sl.mapped()){return false;}
		if(h.strand!=sl.strand()){return false;}
		int start=sl.start(true, true);
		int stop=sl.stop(start, true, true);
		if(!h.rname.equals(sl.rnameS())){return false;}

		if(h.start!=start){return false;}
		if(h.stop!=stop){return false;}
		return(absdif(h.start, start)<=THRESH2 || absdif(h.stop, stop)<=THRESH2);
	}
	
//	public static boolean isCorrectHit(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh,
//			byte[] originalContig, byte[] contig, int cstart, Read r){
//		
//		final int cstop=cstart+trueStop-trueStart;
//		
////		System.out.println("\n"+r.id);
////		System.out.println("         \tstrand"+/*"\tchrom"+*/"\tstart\tstop\t");//+"scaf");
////		System.out.println("Original:\t"+trueStrand+/*"\t"+trueChrom+*/"\t"+trueStart+"\t"+trueStop+"\t");//+new String(originalContig));
////		System.out.println("Mapped:  \t"+ss.strand+/*"\t"+ss.chrom+*/"\t"+ss.start+"\t"+ss.stop+"\t");//+new String(contig));
////		System.out.println("cs:      \t"+trueStrand+/*"\t"+trueChrom+*/"\t"+cstart+"\t"+cstop+"\t");//+new String(contig));
//		
//		if(ss.strand!=trueStrand){return false;}
//		if(originalContig!=null){
//			if(!Arrays.equals(originalContig, contig)){return false;}
//		}else{
//			if(ss.chrom!=trueChrom){return false;}
//		}
//
//		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
//		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
////		return (absdif(ss.start, trueStart)<=thresh && absdif(ss.stop, trueStop)<=thresh);
//		return (absdif(ss.start, cstart)<=thresh && absdif(ss.stop, cstop)<=thresh);
//	}
//	
//	
//	public static boolean isCorrectHitLoose(SiteScore ss, int trueChrom, byte trueStrand, int trueStart, int trueStop, int thresh,
//			byte[] originalContig, byte[] contig, int cstart){
//		if(ss.strand!=trueStrand){return false;}
//		if(originalContig!=null){
//			if(!Arrays.equals(originalContig, contig)){return false;}
//		}else{
//			if(ss.chrom!=trueChrom){return false;}
//		}
//
//		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
//		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
//		int cstop=cstart+trueStop-trueStart;
////		return (absdif(ss.start, trueStart)<=thresh || absdif(ss.stop, trueStop)<=thresh);
//		return (absdif(ss.start, cstart)<=thresh || absdif(ss.stop, cstop)<=thresh);
//	}
	
	private static final int absdif(int a, int b){
		return a>b ? a-b : b-a;
	}

	public static FileFormat ffLoose=null;
	public static FileFormat ffStrict=null;
	public static TextStreamWriter tswLoose=null;
	public static TextStreamWriter tswStrict=null;

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
	public static boolean USE_BITSET=true;
	public static BitSet seen=null;
	
}
