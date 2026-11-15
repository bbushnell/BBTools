package pacbio;

import java.util.ArrayList;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.SiteScoreR;

/**
 * Votes on already stacked {@link SiteScoreR} arrays by sweeping fixed genomic
 * intervals, scoring sites relative to each interval's center, and retaining
 * the sites that earn enough accumulated votes before writing the original
 * tab-separated arrays back out.
 * @author Brian Bushnell
 * @date Jul 18, 2012
 */
public class ProcessStackedSitesNormalized{
	
	/**
	 * Parses the input and output paths plus filtering knobs, processes the
	 * stacked-site table, and prints retained/correct site counters.
	 * @param args Input path, output path, and optional flag=value settings
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
		
		try{
			process(infile, outfile);
		}catch(RuntimeException e){
			System.err.println("ProcessStackedSitesNormalized failed: "+e.getMessage());
			throw e;
		}

		//n [PSSN#003 note] LOW/cosmetic: args[0..1] remain unguarded (AIOOBE<2).
		//n Buffer.close() now waits for the writer and reports close errors; process() aborts it
		//n on interval scoring exceptions so rejected scores cannot leave the JVM alive.
		System.out.println("Sites In:\t"+sitesIn+"    \t"+correctPercent(correctIn, sitesIn));
		System.out.println("Sites Out:\t"+sitesOut+"    \t"+correctPercent(correctOut, sitesOut));
		t.stop();
		System.out.println("Time: \t"+t);
	}

	/** Formats a correct-site percentage, treating an empty denominator as 0%. */
	private static String correctPercent(long correct, long sites){
		double percent=(sites<1 ? 0 : correct*100d/sites);
		return Tools.format("%.3f%% correct", percent);
	}

	/**
	 * Reads stacked site records, scans each chromosome in {@link #INTERVAL}
	 * windows, and writes arrays after all overlapping windows have voted.
	 * @param infile Input tab-delimited SiteScoreR arrays
	 * @param outfile Output arrays containing sites with enough retain votes
	 */
	public static void process(String infile, String outfile){
		validateSettings();
		
		Buffer buffer=new Buffer(3, infile, outfile);
		boolean success=false;
		try{
			int chrom=buffer.chrom;
			int start=buffer.min;
			int stop=buffer.min+INTERVAL-1;
			
			while(buffer.array[0]!=null){

				processInterval(buffer, chrom, start, stop);

				start+=INTERVAL;
				stop+=INTERVAL;
				boolean advanced=buffer.advanceToInterval(start, stop, chrom);
				if(!advanced){
					chrom=buffer.chrom;
					start=buffer.min;
					stop=start+INTERVAL-1;
				}
			}
			success=true;
		}finally{
			if(success){buffer.close();}
			else{buffer.abort();}
		}
	}

	/** Validates finite public knobs before output construction where possible. */
	private static void validateSettings(){
		if(INTERVAL<1){
			throw new IllegalArgumentException("interval must be positive: "+INTERVAL);
		}
		checkFinite("fractiontoretain1", FRACTION_TO_RETAIN1);
		checkFinite("fractiontoretain2", FRACTION_TO_RETAIN2);
		checkFinite("minfractionfromreadends", MIN_FRACTION_FROM_READ_ENDS);
		checkFinite("scorethresh", SCORE_THRESH);
		checkFinite("centerweight", CENTER_WEIGHT);
	}

	/** Rejects NaN and infinity settings that would poison interval scores. */
	private static void checkFinite(String name, float value){
		if(!Float.isFinite(value)){
			throw new IllegalArgumentException(name+" must be finite: "+value);
		}
	}
	
	/**
	 * Scores plus- and minus-strand candidates contained in one interval.
	 * Sites too short or too near their own ends for the configured fractional
	 * end trim do not receive votes for this interval.
	 * @param buffer Sliding chromosome buffer
	 * @param chrom Chromosome being scanned
	 * @param start Inclusive interval start
	 * @param stop Inclusive interval stop
	 */
	private static void processInterval(Buffer buffer, int chrom, int start, int stop){

		ArrayList<SiteScoreR> plus=new ArrayList<SiteScoreR>();
		ArrayList<SiteScoreR> minus=new ArrayList<SiteScoreR>();

		for(Ssra ssra : buffer.array){
			if(ssra==null){continue;}
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
	 * Sorts one strand's candidates by normalized score and adds weighted retain
	 * votes to the retained prefix, subject to count and score-drop thresholds.
	 * @param list Same-strand candidates for the current interval
	 * @return Planned retained count before a possible score-drop break
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
		//n [ProcessStackedSitesNormalized#002 comprehension] LOW/note: list is sorted by NCOMP (normalizedScore) but the break
		//n test below uses best.SCORE (`best.score-b.score>best.score*SCORE_THRESH`) — a DIFFERENT key — so b.score can exceed
		//n best.score and the drop-off cutoff is approximate w.r.t. the sort order. Not a crash; intended-ish (score gates while
		//n normalizedScore orders). Also a single-site interval → sites=1 → retain=(int)(1*0.75)=0 → 0 votes → that lone site is
		//n dropped unless it earns votes from OTHER overlapping intervals (vote-based design, retainVotes accumulate across the
		//n sliding window). best==null only when list is empty, and then retain is also 0 so the loop body (which derefs best) never runs — safe.
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
	 * Parses one tab-delimited stacked group and summarizes its coordinate and
	 * score bounds.
	 * @param s One input line containing tab-delimited SiteScoreR records
	 * @return Parsed group summary
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
	 * Returns a score adjusted for read/reference length bias and distance from
	 * the current interval's nearest end.
	 * @param ssr Candidate site
	 * @param endDist Distance from the site to the interval edge being scored
	 * @return Adjusted score used only for interval-local ranking
	 */
	public static float normalizedScore(SiteScoreR ssr, int endDist){
		final float lim1=0.008f;
		final float lim2=-lim1;
		final int reflen=ssr.reflen();
		if(ssr.readlen<=0){
			throw new IllegalArgumentException("read length must be positive: "+ssr.readlen);
		}
		
		int dif=ssr.readlen-reflen; //Positive for insertions, negative for deletions
		float modifier=dif/(ssr.readlen*4f); //Prioritize reads with insertions over deletions, to correct for scoring bias
		if(modifier<lim2){modifier=lim2;}
		if(modifier>lim1){modifier=lim1;}
		
		int maxEndDist=(reflen-INTERVAL)/2;
		if(maxEndDist<=0){
			throw new IllegalArgumentException("score normalization requires reflen greater than INTERVAL+1: reflen="+reflen+", interval="+INTERVAL);
		}
		float modifier2=CENTER_WEIGHT*endDist/(float)maxEndDist; //Prioritize reads centered on this interval

		float f=ssr.score*(1+modifier+modifier2);
		if(!Float.isFinite(f)){
			throw new IllegalArgumentException("normalized score must be finite: "+f);
		}
		return f;
	}
	
	/**
	 * Finds the maximum alignment score among sites in an Ssra entry that overlap
	 * the specified genomic interval.
	 * @param ssra Container with alignment sites
	 * @param min Start position of target region
	 * @param max End position of target region
	 * @return Highest score among overlapping sites, or -1 if none qualify
	 */
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
	
	public static class Ssra{

		/** Empty constructor for legacy callers that fill fields directly. */
		public Ssra(){}
		
		/**
		 * Stores one stacked input line and its precomputed bounds.
		 * @param array_ Parsed SiteScoreR records
		 * @param chrom_ Chromosome shared by the array
		 * @param min_ Minimum start across the array
		 * @param max_ Maximum stop across the array
		 * @param best_ Highest raw score in the array
		 * @param worst_ Lowest raw score in the array
		 */
		public Ssra(SiteScoreR[] array_, int chrom_, int min_, int max_, int best_, int worst_){
			array=array_;
			chrom=chrom_;
			min=min_;
			max=max_;
			best=best_;
			worst=worst_;
		}
		
		/** Array of SiteScoreR objects representing stacked sites at the current
		 * chromosome. */
		SiteScoreR[] array;
		/** Chromosome or contig identifier for the buffered sites. */
		int chrom;
		/** Minimum genomic position across all buffered sites. */
		int min;
		/** Maximum genomic position across all buffered sites. */
		int max;
		/** Highest alignment score observed across all buffered sites. */
		int best;
		int worst;
		
	}
	
	public static class Buffer{
		
		/**
		 * Opens the input and output streams and primes a chromosome-local
		 * sliding window with up to {@code size} stacked groups.
		 * @param size Number of stacked input groups kept in memory
		 * @param infname_ Input path
		 * @param outfname_ Output path
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
		
		/** Reads one stacked group, updating input counters. */
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
		
		/** Advances by one group, flushing the evicted group if the buffer was full. */
		private boolean advance(){
			if(nextSsra==null){return false;}
			
			Ssra old=add(nextSsra);
			nextSsra=read();
			if(old!=null){write(old);}
			return true;
		}
		
		/** Fills the buffer with groups from the next chromosome. */
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
		
		/**
		 * Advances the buffer so that it contains data overlapping the specified
		 * genomic interval, skipping ahead in chromosome or position as needed.
		 * @param a Start position of target interval
		 * @param b End position of target interval
		 * @param c Target chromosome identifier
		 * @return true if the buffer contains overlapping data, false otherwise
		 */
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
		
		/** Writes all buffered groups and clears the active window. */
		private void purge(){
			for(int i=0; i<array.length; i++){
				Ssra ssra=array[i];
				if(ssra!=null){write(ssra);}
				array[i]=null;
			}
		}
		
		/**
		 * Writes retained sites from one original stacked group, preserving the
		 * tab-delimited group layout when at least one site survives.
		 * @param ssra Original stacked group
		 */
		private void write(Ssra ssra){
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
		 * Appends a same-chromosome group to the window, evicting the oldest
		 * group when the fixed-size buffer is already full.
		 * @param s Group to add
		 * @return Evicted group, or null if no group was evicted
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
		
		/** Recomputes the active chromosome and coordinate range. */
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
		
		/** Flushes all remaining groups, closes input, and joins output. */
		public void close(){
			purge();
			while(fill()){purge();}
			boolean error=tf.close();
			error|=tsw.poisonAndWait();
			if(error){throw new RuntimeException("Error closing "+outfname);}
		}

		/** Closes both streams after a processing failure without writing more groups. */
		public void abort(){
			tf.close();
			tsw.poisonAndWait();
		}
		
		public int max=-1;
		public int min=-1;
		public int chrom=-1;
		
		public final Ssra[] array;
		private Ssra nextSsra;
		public final String infname;
		public final String outfname;
		private TextFile tf;
		private TextStreamWriter tsw;
		
	}
	
	public static int MIN_LENGTH_TO_RETAIN=0;
	public static boolean RETAIN_ALL=false;
	
	public static long sitesIn=0;
	public static long correctIn=0;
	public static long sitesOut=0;
	public static long correctOut=0;
	public static float FRACTION_TO_RETAIN1=0.75f;
	public static float FRACTION_TO_RETAIN2=0.3f;
	public static int MIN_SITES_TO_DISCARD=0;
	public static int SITES_TO_RETAIN1=8;
	public static int SITES_TO_RETAIN2=16;
	public static int MIN_VOTES_TO_RETAIN=5;
//	public static int MIN_DIST_FROM_READ_ENDS=25;
	public static float MIN_FRACTION_FROM_READ_ENDS=0.35f;
	public static float SCORE_THRESH=0.034f;
	public static float CENTER_WEIGHT=0.015f;
	public static int INTERVAL=12;
	
}
