package pacbio;

import java.util.ArrayList;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.SiteScoreR;

/**
 * @author Brian Bushnell
 * @date Jul 18, 2012
 *
 */
public class ProcessStackedSitesNormalized {
	
	/**
	 * Program entry point. Parses command-line parameters and executes site processing.
	 * Accepts input/output files and various filtering parameters.
	 * @param args Command-line arguments including input file, output file, and parameters
	 */
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		String infile=args[0];
		String outfile=args[1];
		
		for(int i=2; i<args.length; i++){
			String[] split=args[i].toLowerCase().split("=");
			String a=split[0];
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("scorethresh")){
				SCORE_THRESH=Float.parseFloat(b);
			}else if(a.equals("interval")){
				INTERVAL=Integer.parseInt(b);
			}else if(a.equals("minsitestodiscard")){
				MIN_SITES_TO_DISCARD=Integer.parseInt(b);
			}else if(a.equals("minlength")){
				MIN_LENGTH_TO_RETAIN=Integer.parseInt(b);
			}else if(a.equals("retainall")){
				RETAIN_ALL=Parse.parseBoolean(b);
				if(RETAIN_ALL){MIN_VOTES_TO_RETAIN=0;}
			}else if(a.equals("fractiontoretain1")){
				FRACTION_TO_RETAIN1=Float.parseFloat(b);
			}else if(a.equals("fractiontoretain2")){
				FRACTION_TO_RETAIN2=Float.parseFloat(b);
			}else if(a.equals("centerweight")){
				CENTER_WEIGHT=Float.parseFloat(b);
			}else if(a.equals("sitestoretain1")){
				SITES_TO_RETAIN1=Integer.parseInt(b);
			}else if(a.equals("sitestoretain2")){
				SITES_TO_RETAIN2=Integer.parseInt(b);
			}else if(a.equals("minvotestoretain")){
				MIN_VOTES_TO_RETAIN=Integer.parseInt(b);
			}else if(a.equals("mindistfromreadends")){
//				MIN_DIST_FROM_READ_ENDS=Integer.parseInt(b);
//				throw new RuntimeException("Deprecated - use minfractionfromreadends instead.");
				int x=Integer.parseInt(b);
				float f=x/((150-INTERVAL)*.5f);
				System.err.println("Warning - mindistfromreadends is deprecated.  Setting minfractionfromreadends = "+Tools.format("%.3f",f));
				MIN_FRACTION_FROM_READ_ENDS=f;
			}else if(a.equals("minfractionfromreadends")){
				MIN_FRACTION_FROM_READ_ENDS=Float.parseFloat(b);
			}else{
				assert(false) : "Unknown parameter "+a;
			}
		}
		
		process(infile, outfile);
		
		System.out.println("Sites In:\t"+sitesIn+"    \t"+Tools.format("%.3f%% correct",correctIn*100d/sitesIn));
		System.out.println("Sites Out:\t"+sitesOut+"    \t"+Tools.format("%.3f%% correct",correctOut*100d/sitesOut));
		t.stop();
		System.out.println("Time: \t"+t);
	}

	/**
	 * @param infile
	 * @param outfile
	 */
	public static void process(String infile, String outfile) {
		
		Buffer buffer=new Buffer(3, infile, outfile);
		
		int chrom=buffer.chrom;
		int start=buffer.min;
		int stop=buffer.min+INTERVAL-1;
		
		assert(buffer.array[0]!=null);
		while(buffer.array[0]!=null){
			
			processInterval(buffer, chrom, start, stop);
			
			start+=INTERVAL;
			stop+=INTERVAL;
			boolean success=buffer.advanceToInterval(start, stop, chrom);
			if(!success){
				chrom=buffer.chrom;
				start=buffer.min;
				stop=start+INTERVAL-1;
			}
		}
		buffer.close();
	}
	
	/**
	 * Processes alignment sites within a specific genomic interval.
	 * Separates sites by strand, applies distance filtering from read ends,
	 * and calculates normalized scores for qualifying sites.
	 *
	 * @param buffer Buffer containing stacked site data
	 * @param chrom Chromosome identifier
	 * @param start Start position of the interval
	 * @param stop End position of the interval
	 */
	private static void processInterval(Buffer buffer, int chrom, int start, int stop){

		ArrayList<SiteScoreR> plus=new ArrayList<SiteScoreR>();
		ArrayList<SiteScoreR> minus=new ArrayList<SiteScoreR>();

		for(Ssra ssra : buffer.array){
//			if(Tools.isWithin(start-MIN_DIST_FROM_READ_ENDS, stop+MIN_DIST_FROM_READ_ENDS,  ssra.min, ssra.max)){
			if(Tools.isWithin(start, stop,  ssra.min, ssra.max)){
				for(SiteScoreR ssr : ssra.array){
					
					int x=(int)((((ssr.stop-ssr.start+1)-INTERVAL)/2)*MIN_FRACTION_FROM_READ_ENDS);
					if(x<0){x=0;}
					
					if(ssr.readlen>=MIN_LENGTH_TO_RETAIN){
						if(Tools.isWithin(start, stop, ssr.start+x, ssr.stop-x)){
							ssr.normalizedScore=normalizedScore(ssr, Tools.min(start-ssr.start, ssr.stop-stop));
							if(ssr.strand==Shared.PLUS){
								plus.add(ssr);
							}else{
								minus.add(ssr);
							}
						}
					}

				}
			}
		}
		markRetain(plus);
		markRetain(minus);
		
	}
	
//	private static final int markRetain_old(ArrayList<SiteScoreR> list){
////		Shared.sort(list, SiteScoreR.NCOMP);
//		assert(list.size()<2 || list.get(0).normalizedScore>=list.get(1).normalizedScore) : list.get(0)+"\t"+list.get(1);
//
//		int sites=list.size()-MIN_SITES_TO_DISCARD; //Always ignore worst site(s).
//
//		int retain=(int)(sites*FRACTION_TO_RETAIN1);
//		if(retain>SITES_TO_RETAIN1){
//			int temp=(int)((retain-SITES_TO_RETAIN1)*FRACTION_TO_RETAIN2);
////			System.out.println("sites="+sites+", retain="+retain+", temp="+temp);
//			retain=SITES_TO_RETAIN1+temp;
//		}
//		retain=Tools.min(retain, SITES_TO_RETAIN2);
////		System.out.println("retain2="+retain);
//
////		for(int i=0; i<retain; i++){
////			list.get(i).retainVotes++;
////		}
//		Shared.sort(list);
//
//		final SiteScoreR best=(list!=null && list.size()>0 ? list.get(0) : null);
//		for(int i=0; i<retain; i++){
//			final SiteScoreR b=list.get(i);
//			if(i>0){
////				SiteScoreR a=list.get(i-1);
////				if(a.score-b.score>a.score*0.03f){break;}
//				if(best.score-b.score>best.score*0.034f){break;}
//			}
//
//			if(i==0){
//				b.retainVotes+=5;
//			}else if(i<3){
//				b.retainVotes+=3;
//			}else if(i<6){
//				b.retainVotes+=2;
//			}else{
//				b.retainVotes++;
//			}
//		}
//
//		return retain;
//	}
	
	/**
	 * Marks alignment sites for retention based on scoring hierarchy.
	 * Applies tiered retention fractions and assigns retention votes
	 * with highest weight for best-scoring sites.
	 *
	 * @param list List of alignment sites to evaluate
	 * @return Number of sites marked for retention
	 */
	private static final int markRetain(ArrayList<SiteScoreR> list){
//		Shared.sort(list, SiteScoreR.NCOMP);
//		assert(list.size()<2 || list.get(0).normalizedScore>=list.get(1).normalizedScore) : list.get(0)+"\t"+list.get(1);
		
		int sites=list.size()-MIN_SITES_TO_DISCARD; //Always ignore worst site(s).
		
		int retain=(int)(sites*FRACTION_TO_RETAIN1);
		if(retain>SITES_TO_RETAIN1){
			int temp=(int)((retain-SITES_TO_RETAIN1)*FRACTION_TO_RETAIN2);
//			System.out.println("sites="+sites+", retain="+retain+", temp="+temp);
			retain=SITES_TO_RETAIN1+temp;
		}
		retain=Tools.min(retain, SITES_TO_RETAIN2);
		
		if(RETAIN_ALL){retain=sites;}
		
//		System.out.println("retain2="+retain);
		
//		for(int i=0; i<retain; i++){
//			list.get(i).retainVotes++;
//		}
		Shared.sort(list, SiteScoreR.NCOMP);
//		assert(false) : SCORE_THRESH;
		final SiteScoreR best=(list!=null && list.size()>0 ? list.get(0) : null);
		for(int i=0; i<retain; i++){
			final SiteScoreR b=list.get(i);
			if(i>0){
//				SiteScoreR a=list.get(i-1);
//				if(a.score-b.score>a.score*0.03f){break;}
				if(!RETAIN_ALL && best.score-b.score>best.score*SCORE_THRESH){break;}
			}
			
			if(i==0){
				b.retainVotes+=5;
			}else if(i<4){
				b.retainVotes+=3;
			}else if(i<8){
				b.retainVotes+=2;
			}else{
				b.retainVotes++;
			}
		}
		
		return retain;
	}
	
	/**
	 * Parses a tab-delimited string into a SiteScoreR array wrapper.
	 * Extracts alignment scores, determines genomic bounds, and tracks
	 * best/worst scores across all sites in the line.
	 *
	 * @param s Tab-delimited string containing alignment site data
	 * @return Ssra object containing parsed alignment sites and metadata
	 */
	public static Ssra toSrar(String s){
		String[] split=s.split("\t");
		SiteScoreR[] scores=new SiteScoreR[split.length];
		int min=Integer.MAX_VALUE;
		int max=Integer.MIN_VALUE;
		int worst=Integer.MAX_VALUE;
		int best=Integer.MIN_VALUE;
		int chrom=-1;
		
		for(int i=0; i<split.length; i++){
			SiteScoreR ssr=scores[i]=SiteScoreR.fromText(split[i]);
			
//			int dif=ssr.readlen-ssr.reflen(); //Positive for insertions, negative for deletions
//			float modifier=dif/(float)(ssr.readlen*4);
//			if(modifier<lim2){modifier=lim2;}
//			if(modifier>lim1){modifier=lim1;}
//			ssr.normalizedScore=(int)ssr.score*(1+modifier);
			
			
			min=Tools.min(min, ssr.start);
			max=Tools.max(max, ssr.stop);
			worst=Tools.min(worst, ssr.score);
			best=Tools.max(best, ssr.score);
			assert(chrom==-1 || chrom==ssr.chrom);
			chrom=ssr.chrom;
		}
		Ssra ssra=new Ssra(scores, chrom, min, max, best, worst);
		return ssra;
	}
	
	/**
	 * Calculates normalized score for an alignment site.
	 * Applies insertion/deletion bias correction and center-weighting
	 * to prioritize sites near interval centers.
	 *
	 * @param ssr Alignment site to score
	 * @param endDist Distance from interval center
	 * @return Normalized score value
	 */
	public static float normalizedScore(SiteScoreR ssr, int endDist){
		final float lim1=0.008f;
		final float lim2=-lim1;
		
		
		int dif=ssr.readlen-ssr.reflen(); //Positive for insertions, negative for deletions
		float modifier=dif/(float)(ssr.readlen*4); //Prioritize reads with insertions over deletions, to correct for scoring bias
		if(modifier<lim2){modifier=lim2;}
		if(modifier>lim1){modifier=lim1;}
		
		int maxEndDist=(ssr.reflen()-INTERVAL)/2;
//		float modifier2=(0.03f*endDist)/maxEndDist;
		float modifier2=CENTER_WEIGHT*endDist/(float)maxEndDist; //Prioritize reads centered on this interval

		float f=ssr.score*(1+modifier+modifier2);
		return f;
	}
	
	/** Finds highest score of ssr's fully covering this site */
	public static int maxScore(Ssra ssra, final int min, final int max){
		assert(Tools.overlap(min, max, ssra.min, ssra.max));
		assert(Tools.isWithin(min, max, ssra.min, ssra.max));
		
		int best=-1;
		for(SiteScoreR ssr : ssra.array){
			if(ssr.start>min){break;}
			if(max>=ssr.stop){
				best=Tools.max(best, ssr.score);
				if(best>=ssra.best){break;}
			}
		}
		return best;
	}
	
	/** Container for SiteScoreR array with genomic bounds and score statistics.
	 * Tracks chromosome, position range, and best/worst scores for efficient processing. */
	public static class Ssra{

		public Ssra(){}
		
		/**
		 * Constructs Ssra with alignment sites and precomputed metadata.
		 *
		 * @param array_ Array of alignment sites
		 * @param chrom_ Chromosome identifier
		 * @param min_ Minimum genomic position
		 * @param max_ Maximum genomic position
		 * @param best_ Highest alignment score
		 * @param worst_ Lowest alignment score
		 */
		public Ssra(SiteScoreR[] array_, int chrom_, int min_, int max_, int best_, int worst_){
			array=array_;
			chrom=chrom_;
			min=min_;
			max=max_;
			best=best_;
			worst=worst_;
		}
		
		/** SiteScoreR array sorted by start loc, ascending */
		SiteScoreR[] array;
		/** All contents must have same chromosome / contig */
		int chrom;
		/** Minimum location in array */
		int min;
		/** Maximum location in array */
		int max;
		/** Top score in array */
		int best;
		/** Bottom score in array */
		int worst;
		
	}
	
	/**
	 * Sliding window buffer for processing alignment site files.
	 * Maintains array of Ssra objects for efficient interval-based processing
	 * with automatic file I/O and genomic coordinate tracking.
	 */
	public static class Buffer{
		
		/**
		 * Constructs buffer with specified size and file paths.
		 * Initializes file readers/writers and fills initial buffer contents.
		 *
		 * @param size Number of Ssra objects to buffer simultaneously
		 * @param infname_ Input file path
		 * @param outfname_ Output file path
		 */
		public Buffer(int size, String infname_, String outfname_){
			assert(!infname_.equalsIgnoreCase(outfname_)) : infname_+" == "+outfname_; //Not a complete test
			array=new Ssra[size];
			infname=infname_;
			outfname=outfname_;
			tf=new TextFile(infname, true);
			tsw=new TextStreamWriter(outfname, true, false, true);
			tsw.start();
			nextSsra=read();
			fill();
			
		}
		
		/**
		 * Reads next line from input file and parses into Ssra object.
		 * Updates global site count statistics.
		 * @return Parsed Ssra object, or null if end of file reached
		 */
		public Ssra read(){
			String s=tf.nextLine();
			if(s==null){
				tf.close();
				return null;
			}
			Ssra ssra=toSrar(s);
			sitesIn+=ssra.array.length;
			return ssra;
		}
		
		/**
		 * Advances buffer by one position, adding next Ssra and writing oldest.
		 * Maintains sliding window behavior for continuous processing.
		 * @return true if advancement succeeded, false if no more data
		 */
		private boolean advance(){
			if(nextSsra==null){return false;}
			
			Ssra old=add(nextSsra);
			nextSsra=read();
			if(old!=null){write(old);}
			return true;
		}
		
		/** Starting with an empty array, fill with next chrom */
		private boolean fill(){
			assert(array[0]==null);
			if(nextSsra==null){return false;}
			int c=nextSsra.chrom;
			for(int i=0; i<array.length && nextSsra!=null && c==nextSsra.chrom; i++){
				array[i]=nextSsra;
				nextSsra=read();
			}
			setLimits();
			return true;
		}
		
		/** Lowest alignment score among all sites in array */
		public boolean advanceToInterval(final int a, final int b, final int c){
			
			while(chrom<c || (chrom==c && max<a)){
				purge();
				boolean success=fill();
				if(!success){return false;}
			}
			
			assert(array[0]!=null && chrom>=c);
//			if(chrom>c || min>b){return false;} //Went past target
			
			while(array[0].max<a && nextSsra!=null && nextSsra.chrom==c){
				advance();
			}
			
			return chrom==c && Tools.overlap(a, b, min, max);
		}
		
		/** Writes all buffered Ssra objects to output and clears the buffer.
		 * Used when switching chromosomes or finishing processing. */
		private void purge() {
			for(int i=0; i<array.length; i++){
				Ssra ssra=array[i];
				if(ssra!=null){write(ssra);}
				array[i]=null;
			}
		}
		
		/**
		 * Writes Ssra object to output file after applying retention filtering.
		 * Only outputs sites with sufficient retention votes, updating statistics.
		 * @param ssra Ssra object containing sites to potentially write
		 */
		private void write(Ssra ssra) {
			String tab="";
			StringBuilder sb=new StringBuilder(ssra.array.length*48);
			
			final long sitesOut_0=sitesOut;
			for(SiteScoreR ssr : ssra.array){
				
//				if(ssr.weight>0){
//					ssr.normalizedScore/=ssr.weight;
//				}
				
				if(ssr.correct){correctIn++;}
				if(ssr.retainVotes>=MIN_VOTES_TO_RETAIN){
					sitesOut++;
					if(ssr.correct){correctOut++;}
					sb.append(tab);
					sb.append(ssr.toText());
					tab="\t";
				}
			}
			if(sitesOut_0==sitesOut){return;}
			sb.append('\n');
			tsw.print(sb);
		}

		/**
		 * Adds new Ssra to buffer using sliding window mechanism.
		 * Maintains buffer size by returning oldest Ssra when buffer is full.
		 * @param s New Ssra object to add to buffer
		 * @return Oldest Ssra that was displaced, or null if buffer had space
		 */
		public Ssra add(Ssra s){
			
			assert(array[0]==null || array[0].chrom==s.chrom);
			
			Ssra r=null;
			if(array[array.length-1]==null){
				//insert in first null loc
				for(int i=0; i<array.length; i++){
					if(array[i]==null){
						array[i]=s;
						break;
					}
				}
			}else{
				r=array[0];
				for(int i=1; i<array.length; i++){
					array[i-1]=array[i];
				}
				array[array.length-1]=s;
			}
			
			setLimits();
			
			return r;
		}
		
		/** Recalculates genomic bounds across all buffered Ssra objects.
		 * Updates min, max, and chrom fields based on current buffer contents. */
		private void setLimits(){
			max=Integer.MIN_VALUE;
			min=Integer.MAX_VALUE;
			chrom=array[0].chrom;
			for(int i=0; i<array.length; i++){
				if(array[i]!=null){
					min=Tools.min(min, array[i].min);
					max=Tools.max(max, array[i].max);
				}
			}
		}
		
		/** Closes buffer and associated file streams.
		 * Writes any remaining buffered data and properly terminates I/O. */
		public void close(){
			purge();
			while(fill()){purge();}
			tf.close();
			tsw.poison();
		}
		
		/** Maximum genomic position across all currently buffered Ssra objects */
		public int max=-1;
		/** Minimum genomic position across all currently buffered Ssra objects */
		public int min=-1;
		/** Current chromosome identifier for all buffered data */
		public int chrom=-1;
		
		/**
		 * Fixed-size array holding buffered Ssra objects for sliding window processing
		 */
		public final Ssra[] array;
		/** Next Ssra object read from file, waiting to be added to buffer */
		private Ssra nextSsra;
		/** Input file path for reading alignment site data */
		public final String infname;
		/** Output file path for writing filtered alignment results */
		public final String outfname;
		/** Text file reader for input alignment data */
		private TextFile tf;
		/** Text stream writer for output filtered results */
		private TextStreamWriter tsw;
		
	}
	
	/** Minimum read length required for a site to be considered for retention */
	public static int MIN_LENGTH_TO_RETAIN=0;
	/** If true, bypasses normal filtering and retains all qualifying sites */
	public static boolean RETAIN_ALL=false;
	
	/** Total count of input alignment sites processed */
	public static long sitesIn=0;
	/** Count of correct input alignment sites (for accuracy tracking) */
	public static long correctIn=0;
	/** Total count of output alignment sites after filtering */
	public static long sitesOut=0;
	/** Count of correct output alignment sites (for accuracy tracking) */
	public static long correctOut=0;
	/** Primary fraction of sites to retain in first tier of filtering */
	public static float FRACTION_TO_RETAIN1=0.75f;
	/** Secondary fraction applied to excess sites beyond SITES_TO_RETAIN1 */
	public static float FRACTION_TO_RETAIN2=0.3f;
	/** Minimum number of lowest-scoring sites to always discard */
	public static int MIN_SITES_TO_DISCARD=0;
	/** Maximum sites retained at primary fraction before secondary filtering */
	public static int SITES_TO_RETAIN1=8;
	/** Absolute maximum number of sites to retain regardless of fractions */
	public static int SITES_TO_RETAIN2=16;
	/** Minimum retention votes required for a site to be written to output */
	public static int MIN_VOTES_TO_RETAIN=5;
//	public static int MIN_DIST_FROM_READ_ENDS=25;
	/** Minimum distance from read ends as fraction of read length for inclusion */
	public static float MIN_FRACTION_FROM_READ_ENDS=0.35f;
	/** Score difference threshold as fraction of best score for retention cutoff */
	public static float SCORE_THRESH=0.034f;
	/** Weighting factor for prioritizing sites centered within intervals */
	public static float CENTER_WEIGHT=0.015f;
	/** Size of genomic intervals for processing alignment sites */
	public static int INTERVAL=12;
	
}
