package align2;

import java.io.File;
import java.util.BitSet;

import fileIO.TextFile;
import shared.Parse;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import stream.CustomHeader;
import stream.Read;
import stream.SamLine;

public class MakeRocCurve {
	
	
	public static void main(String[] args){

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
//			}else if(a.equals("ssaha2") || a.equals("subtractleadingclip")){
//				SamLine.SUBTRACT_LEADING_SOFT_CLIP=Parse.parseBoolean(b);
			}else if(a.equals("blasr")){
				BLASR=Parse.parseBoolean(b);
			}else if(a.equals("bitset")){
				USE_BITSET=Parse.parseBoolean(b);
			}else if(a.equals("thresh")){
				THRESH2=Integer.parseInt(b);
			}else if(a.equals("allowspaceslash")){
				allowSpaceslash=Parse.parseBoolean(b);
			}else if(a.equals("outputerrors")){
//				OUTPUT_ERRORS=true;
			}else if(i==0 && args[i].indexOf('=')<0 && (a.startsWith("stdin") || new File(args[0]).exists())){
				in=args[0];
			}else if(i==1 && args[i].indexOf('=')<0 && Tools.isDigit(a.charAt(0))){
				reads=Parse.parseKMG(a);
			}
		}
		
		if(USE_BITSET){
			int x=400000;
			if(reads>0 && reads<=Integer.MAX_VALUE){x=(int)reads;}
			try {
				seen=new BitSet(x);
			} catch (Exception e) {
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
	
	public static void process(String samfile){
		TextFile tf=new TextFile(samfile, false);
		
		String s=null;
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			char c=s.charAt(0);
			if(c!='@'/* && c!=' ' && c!='\t'*/){
				SamLine sl=new SamLine(s);
				final int id=((((int)sl.parseNumericId())<<1)|sl.pairnum());
				assert(sl!=null);
				Read r=sl.toRead(true);
				if(r!=null){
					r.samline=sl;
					if(sl.primary() && (seen==null || !seen.get(id))){
						if(seen!=null){seen.set(id);}
						calcStatistics1(r, sl);
					}
				}else{
					assert(false) : "'"+"'";
					System.err.println("Bad read from line '"+s+"'");
				}
//				calcStatistics1(r);
			}
		}
		tf.close();
	}
	
	public static String header(){
		return "minScore\tmapped\tretained\ttruePositiveStrict\tfalsePositiveStrict\ttruePositiveLoose" +
				"\tfalsePositiveLoose\tfalseNegative\tdiscarded\tambiguous";
	}
	
	public static void gradeList(long reads){

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
		primaryA[q]++;
		if(q<0){q=0;}
		if(q>=discardedA.length){q=discardedA.length-1;}
		
		if(r.discarded()/* || r.mapScore==0*/){
			discardedA[q]++;
			unmappedA[q]++;
		}else if(r.ambiguous()){
//			assert(r.mapped()) : "\n"+r+"\n"+sl+"\n";
			if(r.mapped()){mappedA[q]++;}
			ambiguousA[q]++;
		}else if(r.mapScore<1){
			unmappedA[q]++;
		}else if(!r.mapped()){
			unmappedA[q]++;
		}
//		else if(r.mapScore<=minQuality){
//			if(r.mapped()){mappedA[q]++;}
//			ambiguousA[q]++;
//		}
		else{

			mappedA[q]++;
			mappedRetainedA[q]++;

			if(parsecustom){
				CustomHeader h=new CustomHeader(sl.qname, sl.pairnum());
				boolean strict=isCorrectHit(sl, h);
				boolean loose=isCorrectHitLoose(sl, h);

//				SiteScore os=r.originalSite;
//				int trueChrom=os.chrom;
//				byte trueStrand=os.strand;
//				int trueStart=os.start;
//				int trueStop=os.stop;
//				SiteScore ss=new SiteScore(r.chrom, r.strand(), r.start, r.stop, 0, 0);
//				byte[] originalContig=sl.originalContig();
//				if(BLASR){
//					originalContig=(originalContig==null || Tools.indexOf(originalContig, (byte)'/')<0 ? originalContig :
//						KillSwitch.copyOfRange(originalContig, 0, Tools.lastIndexOf(originalContig, (byte)'/')));
//				}
//				int cstart=sl.originalContigStart();
//
//				boolean strict=isCorrectHit(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH, originalContig, sl.rname(), cstart);
//				boolean loose=isCorrectHitLoose(ss, trueChrom, trueStrand, trueStart, trueStop, THRESH+THRESH2, originalContig, sl.rname(), cstart);
//
//				//				if(!strict){
//				//					System.out.println(ss+", "+new String(originalContig)+", "+new String(sl.rname()));
//				//					assert(false);
//				//				}
//
//				//				System.out.println("loose = "+loose+" for "+r.toText());

				if(loose){
					//					System.err.println("TPL\t"+trueChrom+", "+trueStrand+", "+trueStart+", "+trueStop+"\tvs\t"
					//							+ss.chrom+", "+ss.strand+", "+ss.start+", "+ss.stop);
					truePositiveLooseA[q]++;
				}else{
					//					System.err.println("FPL\t"+trueChrom+", "+trueStrand+", "+trueStart+", "+trueStop+"\tvs\t"
					//							+ss.chrom+", "+ss.strand+", "+ss.start+", "+ss.stop);
					falsePositiveLooseA[q]++;
				}

				if(strict){
					//					System.err.println("TPS\t"+trueStart+", "+trueStop+"\tvs\t"+ss.start+", "+ss.stop);
					truePositiveStrictA[q]++;
				}else{
					//					System.err.println("FPS\t"+trueStart+", "+trueStop+"\tvs\t"+ss.start+", "+ss.stop);
					falsePositiveStrictA[q]++;
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
//			byte[] originalContig, byte[] contig, int cstart){
//		if(ss.strand!=trueStrand){return false;}
//		if(originalContig!=null){
//			if(!Arrays.equals(originalContig, contig)){
//				if(allowSpaceslash && originalContig.length==contig.length+3 && Tools.startsWith(originalContig, contig) &&
//						(Character.isWhitespace(originalContig[originalContig.length-3]))){
//					//do nothing
//				}else{
//					return false;
//				}
//			}
//		}else{
//			if(ss.chrom!=trueChrom){return false;}
//		}
//
//		assert(ss.stop>ss.start) : ss.toText()+", "+trueStart+", "+trueStop;
//		assert(trueStop>trueStart) : ss.toText()+", "+trueStart+", "+trueStop;
//		int cstop=cstart+trueStop-trueStart;
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
