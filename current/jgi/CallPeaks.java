package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Oct 28, 2014
 *
 */
public class CallPeaks {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point for peak calling analysis.
	 * Creates CallPeaks instance and processes histogram data.
	 * @param args Command-line arguments for configuration
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		CallPeaks x=new CallPeaks(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs CallPeaks instance with command-line argument parsing.
	 * Configures all parameters for peak detection and genome analysis.
	 * Processes input/output file formats and filtering thresholds.
	 * @param args Command-line arguments containing parameters and file paths
	 */
	public CallPeaks(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, printClass ? getClass() : null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				in=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("minheight") || a.equals("h")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume") || a.equals("v")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth") || a.equals("w")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak") || a.equals("minp")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak") || a.equals("maxp")){
				maxPeak=(int)Tools.min(Integer.MAX_VALUE, Parse.parseKMG(b));
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(a.equals("smoothradius")){
				smoothRadius=Integer.parseInt(b);
			}else if(a.equals("smoothprogressive")){
				smoothProgressiveFlag=Parse.parseBoolean(b);
			}else if(a.equals("maxradius")){
				maxRadius=Integer.parseInt(b);
			}else if(a.equals("progressivemult")){
				progressiveMult=Float.parseFloat(b);
			}else if(a.equals("ploidy")){
				ploidyClaimed=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("ploidyLogic")){
				ploidyLogic=Integer.parseInt(b);
			}else if(a.equals("column") || a.equals("col") || a.equals("countcolumn")){
				countColumn=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("logscale")){
				logScale=Parse.parseBoolean(b);
			}else if(a.equals("logwidth")){
				logWidth=Double.parseDouble(b);
			}else if(a.equals("logpasses")){
				logPasses=Integer.parseInt(b);
			}
			
			else if(a.equals("byheight")){
				CALL_MODE=BY_HEIGHT;
			}else if(a.equals("byvolume")){
				CALL_MODE=BY_VOLUME;
			}else if(a.equalsIgnoreCase("weightByRelief")){
				weightByRelief=Parse.parseBoolean(b);
			}
			
			else if(in==null && i==0 && Tools.looksLikeInputStream(arg)){
				in=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
//		assert(logScale);
		if(out==null){out="stdout.txt";}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TEXT, null, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that loads histogram data and calls peaks.
	 * Analyzes k-mer frequency distribution to identify coverage peaks
	 * and calculate genome statistics including size and ploidy.
	 * @param t Timer for tracking execution time
	 */
	public void process(Timer t){
		LongList lists[]=loadHistogram(ffin, k);//TODO: Add support for GC here
		LongList hist=lists[0];
		LongList gchist=lists[1];
		ArrayList<Peak> peaks=callPeaks(hist, gchist);
		long sum=shared.Vector.sum(hist.array);
		printPeaks(peaks, k, sum, hist.array, (gchist==null ? null : gchist.array));
		hist=null;
		t.stop();
		System.err.println("Found "+peaks.size()+" peaks in "+t);
		
		peaks=null;
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Static method to create CallPeaks instance and analyze provided histogram arrays.
	 * Convenience method for programmatic peak calling with explicit parameters.
	 *
	 * @param array Main k-mer frequency histogram
	 * @param gcArray GC content histogram (may be null)
	 * @param fname Output file name
	 * @param ow Overwrite existing output files
	 * @param minHeight Minimum peak height threshold
	 * @param minVolume Minimum peak volume threshold
	 * @param minWidth Minimum peak width in bases
	 * @param minPeak Minimum peak position to consider
	 * @param maxPeak Maximum peak position to consider
	 * @param maxPeakCount Maximum number of peaks to report
	 * @param k K-mer size used for histogram generation
	 * @param ploidy Expected genome ploidy
	 * @param logScale Whether to apply log scaling
	 * @param logWidth Log scaling window width
	 * @param list Parameter list for instance creation
	 * @return true if analysis completed without errors
	 */
	public static boolean printPeaks(long[] array, long[] gcArray, String fname, boolean ow, long minHeight, long minVolume, int minWidth,
			int minPeak, int maxPeak, int maxPeakCount, int k, int ploidy, boolean logScale, double logWidth, ArrayList<String> list){
		if(list==null){list=new ArrayList<String>();}
		list.add("out="+fname);
		list.add("ow="+ow);
		list.add("minheight="+minHeight);
		list.add("minvolume="+minVolume);
		list.add("minwidth="+minWidth);
		list.add("minpeak="+minPeak);
		list.add("maxpeak="+maxPeak);
		list.add("maxpeaks="+maxPeakCount);
		list.add("k="+(k<1 ? 31 : k));
		if(logScale){
			list.add("logscale=t");
			list.add("logwidth="+logWidth);
		}
		if(ploidy>0){list.add("ploidy="+ploidy);}
		CallPeaks cp=new CallPeaks(list.toArray(new String[0]));
		ArrayList<Peak> peaks=cp.callPeaks(array, gcArray, array.length);
		cp.printPeaks(peaks, k, shared.Vector.sum(array), array, gcArray);
		return cp.errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads k-mer frequency histogram from input file.
	 * Parses tab-delimited histogram data with optional GC content column.
	 *
	 * @param ff Input file format specification
	 * @param k K-mer size for GC content calculation
	 * @return Array containing [main histogram, GC histogram] as LongList objects
	 */
	public static LongList[] loadHistogram(FileFormat ff, int k){//TODO: Add GC support
		LongList list=new LongList(8000);
		LongList gcList=new LongList(8000);
		boolean hasGC=false;
		TextFile tf=new TextFile(ff);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith("#")){
				//ignore
			}else{
				String[] split=line.split("\\s+");
				if(split.length==1){
					list.add(Long.parseLong(split[0]));
				}else{
					final int x=Integer.parseInt(split[0]);
					final long y=Long.parseLong(split[countColumn]);
					list.set(x, y);
					String last=split[split.length-1];
					if(split.length>2 && last.indexOf('.')>=0 && Tools.isNumeric(last)){
						hasGC=true;
						double gc=Double.parseDouble(last)*0.01;
						gcList.set(x, Math.round(y*gc*k));
					}
				}
			}
		}
		boolean errorState_=tf.close();
		assert(!errorState_) : "Encountered an error when reading "+ff.name()+".\n" +
				"To skip this error message, run with the '-da' flag.";
		
		return new LongList[] {list, hasGC ? gcList : null};
	}
	
	/**
	 * Reduces peak count by condensing smaller peaks into larger ones.
	 * Merges compatible peaks when total count exceeds maximum limit.
	 * Uses height and volume thresholds to determine which peaks to preserve.
	 *
	 * @param in Input list of peaks to condense
	 * @param maxCount Maximum number of peaks to retain
	 * @return Condensed list of peaks within the specified limit
	 */
	private static ArrayList<Peak> condense(ArrayList<Peak> in, int maxCount){
		if(in==null || in.isEmpty()){return in;}
		maxCount=Tools.max(Tools.min(in.size(), maxCount), 1);
		ArrayList<Peak> out=new ArrayList<Peak>(Tools.min(maxCount, in.size()));
		
		final long hlimit, vlimit;
		{
			{
				long[] heights=new long[in.size()];
				for(int i=0; i<in.size(); i++){
					Peak p=in.get(i);
					heights[i]=(callByRawCount ? p.maxHeight2() : p.maxHeight);
				}
				Arrays.sort(heights);
				hlimit=heights[heights.length-maxCount];
			}

			{
				int mc2=(maxCount+1)/2;
				long[] volumes=new long[in.size()];
				for(int i=0; i<in.size(); i++){
					Peak p=in.get(i);
					volumes[i]=(callByRawCount ? p.volume2 : p.volume);
				}
				Arrays.sort(volumes);
				vlimit=volumes[volumes.length-mc2];
			}
		}
		
		for(Peak p : in){
			final long height=(callByRawCount ? p.maxHeight2() : p.maxHeight);
			final long volume=(callByRawCount ? p.volume2 : p.volume);
			if(volume>=vlimit || height>=hlimit){
//				assert(!out.contains(p)) : height+", "+volume+", "+hlimit+", "+vlimit+", "+p;
				out.add(p);
//				System.err.println("Adding "+height+", "+volume+", "+hlimit+", "+vlimit+", "+p);
			}
		}
		
//		for(Peak p : in){
//			final long height=(callByRawCount ? p.maxHeight2() : p.maxHeight);
//			final long volume=(callByRawCount ? p.volume2 : p.volume);
//			if(volume>=vlimit || height>=hlimit){
//				assert(out.contains(p));
//			}else{
//				assert(!out.contains(p));
//				assert(volume<vlimit && height<hlimit);
//			}
//		}
		
		for(Peak p : in){
			final long height=(callByRawCount ? p.maxHeight2() : p.maxHeight);
			final long volume=(callByRawCount ? p.volume2 : p.volume);
			if(volume<vlimit && height<hlimit){
//				assert(!out.contains(p)) : height+", "+volume+", "+hlimit+", "+vlimit+", "+p;
				Peak p2=out.get(0);
				for(Peak temp : out){
					if(Tools.absdif(p.center, temp.center)<Tools.absdif(p.center, p2.center)){
						p2=temp;
					}
				}
				if(p2.compatibleWith(p)){
//					assert(!out.contains(p)) : height+", "+volume+", "+hlimit+", "+vlimit+", "+p;
//					assert((callByRawCount ? p2.maxHeight2() : p2.maxHeight)>=hlimit || (callByRawCount ? p2.volume2 : p2.volume)>=vlimit);
					p2.absorb(p);
//					assert(!out.contains(p)) : height+", "+volume+", "+hlimit+", "+vlimit+", "+p;
//					assert((callByRawCount ? p2.maxHeight2() : p2.maxHeight)>=hlimit || (callByRawCount ? p2.volume2 : p2.volume)>=vlimit);
				}//else discard
			}
		}
		
		return out;
	}
	
	/**
	 * Limits peak width by applying maximum width multiplier constraint.
	 * Adjusts peak start/stop positions relative to center position.
	 * Recalculates peak statistics after width adjustment.
	 *
	 * @param peaks List of peaks to adjust
	 * @param maxWidthMult Maximum width as multiplier of center position
	 * @param counts Histogram array for recalculation
	 */
	private static void capWidth(ArrayList<Peak> peaks, float maxWidthMult, long[] counts){
		float mult=1/maxWidthMult;
		for(Peak p : peaks){
			p.start=(int)Math.round(Tools.max(p.start, p.center*mult));
			p.stop=(int)Math.round(Tools.min(p.stop, p.center*maxWidthMult));
			p.recalculate(counts);
		}
		
//		for(int i=0; i<peaks.)
	}
	
	/**
	 * Outputs peak analysis results including genome statistics.
	 * Calculates and reports ploidy, genome size, heterozygosity rate,
	 * repeat content, and GC content from identified peaks.
	 *
	 * @param peaks List of called peaks
	 * @param k K-mer size used for histogram generation
	 * @param uniqueKmers Total unique k-mers in histogram
	 * @param hist Main frequency histogram array
	 * @param gcHist GC content histogram array (may be null)
	 */
	private void printPeaks(ArrayList<Peak> peaks, int k, long uniqueKmers, long[] hist, long[] gcHist){
		if(ffout==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		double gcmult=1.0/k;
		if(peaks.size()>0){
			try {
//				final Peak p0=peaks.get(0);
//				final int center0=p0.center;
				float minHetRate=0.0003f; //Less than the human rate of ~0.001, but this will miss the first peak for inbred species.
				final float minVolumeFraction=Tools.min(1, singleCopyKmerFraction(minHetRate, k, 2));
				final int ploidyEstimate=calcPloidy(peaks, minVolumeFraction);
				final int ploidy=ploidyClaimed>0 ? ploidyClaimed : ploidyEstimate;
				final float haploidPeakCenter=haploidPeakCenter(peaks, ploidy);
				final long errorKmers=errorKmers(peaks, hist, minVolumeFraction);
				final long genomeSizeInPeaks=genomeSizeInPeaks(peaks, ploidy, haploidPeakCenter);
				final long totalGenomeSize=genomeSize2(peaks, ploidy, haploidPeakCenter, hist);
				final double gcContent=(gcHist==null ? -1 : gcContent(peaks, k));
				final double gcContent2=(gcHist==null ? -1 : gcContent2(peaks, k, hist, gcHist));
				final long repeatSize=repeatSize(peaks, ploidy, haploidPeakCenter);
				final long repeatSize2=hist==null ? -1 : repeatSize2(peaks, ploidy, haploidPeakCenter, hist);
				final long haploidSize=totalGenomeSize/ploidy;
				final long hetLocs=calcHetLocations(peaks, ploidy, haploidPeakCenter, k);
				final double hetRate=(hetLocs/(double)haploidSize)/2;//Final /2 is because a SNP will make BOTH paired chromosomes single-copy
				final double repeatRate=repeatSize*1.0/genomeSizeInPeaks;
				final double repeatRate2=repeatSize2*1.0/totalGenomeSize;
				
				if(verbose) {
					System.err.println("genomeSize="+genomeSizeInPeaks+", \trepeatSize="+repeatSize);
					System.err.println("genomeSize2="+totalGenomeSize+",\trepeatSize2="+repeatSize2);
				}
				
				Peak p0=peaks.get(0);
				Peak ploidyPeak=p0, mainPeak=p0;
				float target=haploidPeakCenter*ploidy;
				int haploidCov;
				for(Peak p : peaks){
					if(p.volume>mainPeak.volume){
						mainPeak=p;
					}
					if(Tools.absdif(p.center, target)<Tools.absdif(ploidyPeak.center, target)){
						ploidyPeak=p;
					}
				}
				if(Tools.max(target,ploidyPeak.center)/(float)Tools.min(target,ploidyPeak.center)<1.3f){
					haploidCov=ploidyPeak.center;
				}else{
					haploidCov=(int)target;
				}
				
//				System.err.println("ploidyEstimate="+ploidyEstimate);
//				System.err.println("genomeSize="+genomeSize);
//				System.err.println("repeatSize="+repeatSize);
//				System.err.println("haploidSize="+haploidSize);
//				System.err.println("hetLocs="+hetLocs);
//				System.err.println("biggestPeak="+biggestPeak(peaks));
//				System.err.println("homozygousPeak("+ploidy+")="+homozygousPeak(peaks, ploidy));
//				System.err.println("homozygousPeak("+ploidyEstimate+")="+homozygousPeak(peaks, ploidyEstimate));
				
				if(ploidy!=ploidyEstimate){
					System.err.println("Warning - ploidy detected at "+ploidyEstimate+" differs from stated ploidy of "+ploidyClaimed);
				}

				if(k>0){bsw.println("#k\t"+k);}
				bsw.println("#unique_kmers\t"+uniqueKmers);
				bsw.println("#error_kmers\t"+errorKmers);
				bsw.println("#genomic_kmers\t"+(uniqueKmers-errorKmers));
				bsw.println("#main_peak\t"+mainPeak.center);
				bsw.println("#genome_size_in_peaks\t"+genomeSizeInPeaks);
				bsw.println("#genome_size\t"+totalGenomeSize);
				if(mainPeak.gc>0){
					bsw.println("#main_peak_gc\t"+Tools.format("%.3f", Tools.mid(0, 1, mainPeak.gc*gcmult/mainPeak.volume)));
				}
				if(gcHist!=null){
					bsw.println("#gc_content_in_peaks\t"+Tools.format("%.3f", Tools.mid(0, 1, gcContent)));
					bsw.println("#gc_content\t"+Tools.format("%.3f", Tools.mid(0, 1, gcContent2)));
				}
				if(ploidy>1 || true){bsw.println("#haploid_genome_size\t"+haploidSize);}
				bsw.println("#fold_coverage\t"+Math.round(haploidPeakCenter));
				if(ploidy>1 || true){bsw.println("#haploid_fold_coverage\t"+haploidCov);}
				bsw.println("#ploidy\t"+ploidy);
				if(ploidy!=ploidyEstimate){bsw.println("#ploidy_detected\t"+ploidyEstimate);}
				if(ploidy>1){bsw.println("#het_rate\t"+Tools.format("%.5f", hetRate));}
				bsw.println("#percent_repeat_in_peaks\t"+Tools.format("%.3f", (100*repeatRate)));
				if(repeatSize2>=0){bsw.println("#percent_repeat\t"+Tools.format("%.3f", (100*repeatRate2)));}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		bsw.println("#start\tcenter\tstop\tmax\tvolume"+(gcHist!=null ? "\tgc":""));
		ByteBuilder bb=new ByteBuilder(200);
		for(Peak p : peaks){
			if(p.volume>=minVolume){
				p.toBytes(bb);
				if(gcHist!=null){bb.tab().append(Tools.mid(0, 1, p.gc*gcmult/p.volume), 3);}
				bb.nl();
				bsw.print(bb);
				bb.setLength(0);
			}
		}
		errorState|=bsw.poisonAndWait();
	}
	
	/**
	 * Identifies the first genomic peak based on volume threshold.
	 * Returns the first peak meeting minimum volume fraction of largest peak.
	 *
	 * @param peaks List of peaks to search
	 * @param minFraction Minimum volume as fraction of largest peak
	 * @return First genomic peak or null if none found
	 */
	private static Peak firstGenomicPeak(ArrayList<Peak> peaks, float minFraction) {
		if(peaks.size()<1) {return null;}
		assert(minFraction<=1);
		final Peak biggest=peaks.get(biggestPeak(peaks));
		final long minVolume=(long)(biggest.volume*minFraction);
		for(Peak p : peaks) {
			if(p.volume>=minVolume) {return p;}
		}
		assert(false) : peaks.size()+", "+biggest;
		return null;
	}
	
	//single-copy kmers as a fraction of haploid genome size (result can range from 0 to ploidy).
	//The math here is all totally wrong and terrible but might be within a factor of 4 of correct
	/**
	 * Calculates expected fraction of single-copy k-mers based on heterozygosity.
	 * Uses mathematical model to estimate single-copy k-mer abundance
	 * as function of heterozygosity rate, k-mer size, and ploidy.
	 *
	 * @param hetRate Heterozygosity rate (fraction of heterozygous sites)
	 * @param k K-mer length
	 * @param ploidy Genome ploidy level
	 * @return Fraction of k-mers expected to be single-copy
	 */
	private static float singleCopyKmerFraction(float hetRate, int k, int ploidy) {
		if(ploidy<2) {return 1;}
		float kmersPerSnp=k;//*(ploidy==2 ? 2 : 1);//Could be 1 or 2 in higher ploidies depending on the structure
		float singleCopyKmers=hetRate*kmersPerSnp; //Fraction of whole genome size
		float asymptote=singleCopyKmers/(1+singleCopyKmers);//Prevents the number from going over 1
		return asymptote*2;
	}
	
	/**
	 * Estimates number of erroneous k-mers before first genomic peak.
	 * Sums histogram counts below the start position of first significant peak.
	 *
	 * @param peaks List of called peaks
	 * @param hist Frequency histogram array
	 * @param minVolumeFraction Minimum volume fraction for genomic peak detection
	 * @return Estimated count of error k-mers
	 */
	private static long errorKmers(ArrayList<Peak> peaks, long[] hist, float minVolumeFraction){
		if(peaks.size()<1) {return 0;}
		final Peak first=firstGenomicPeak(peaks, minVolumeFraction);
		if(first==null || hist==null){return -1;}
		long sum=0;
		for(int i=0; i<first.start; i++){
			sum+=hist[i];
		}
		return sum;
	}
	
	/**
	 * Calculates genome size based on peak volumes and estimated copy numbers.
	 * Multiplies each peak volume by its inferred copy number relative
	 * to haploid coverage level.
	 *
	 * @param peaks List of called peaks
	 * @param ploidy Genome ploidy level
	 * @param haploidPeakCenter Coverage level of haploid peak
	 * @return Estimated genome size in base pairs
	 */
	private static long genomeSizeInPeaks(ArrayList<Peak> peaks, final int ploidy, final float haploidPeakCenter){
		if(peaks.size()<1){return 0;}
		
		long sizeSum=0;
		final double mult=1.0/(Tools.max(1, haploidPeakCenter));
		for(Peak p : peaks){
			long size=p.volume*(Math.round(p.center*mult));
			sizeSum+=size;
		}
		return sizeSum;
	}
	
	/**
	 * Alternative genome size calculation using complete histogram data.
	 * Analyzes all histogram positions starting from first peak
	 * to estimate total genome size including repetitive regions.
	 *
	 * @param peaks List of called peaks for reference
	 * @param ploidy Genome ploidy level
	 * @param haploidPeakCenter Coverage level of haploid peak
	 * @param hist Complete frequency histogram array
	 * @return Estimated genome size including repetitive content
	 */
	private static long genomeSize2(ArrayList<Peak> peaks, final int ploidy, final float haploidPeakCenter, long[] hist){
		if(peaks.size()<1){return 0;}
		
		long sizeSum=0;
		final double mult=1.0/(Tools.max(1, haploidPeakCenter));

		final Peak p0=peaks.get(0);
		for(int i=p0.start; i<hist.length; i++){
			long size=hist[i]*(Tools.max(1, (long)Math.round(i*mult)));
			sizeSum+=size;
		}
		return sizeSum;
	}
	
	/**
	 * Calculates GC content from peak-based analysis.
	 * Weighs each peak's GC content by its volume and copy number.
	 *
	 * @param peaks List of peaks with GC information
	 * @param k K-mer size for normalization
	 * @return Estimated GC content as fraction (0.0 to 1.0)
	 */
	private static double gcContent(ArrayList<Peak> peaks, int k){
		if(peaks.size()<1){return 0;}
		
		long sizeSum=0, gcSum=0;
		final Peak p0=peaks.get(0);
		final int center0=p0.center;
		final double mult=1.0/(Tools.max(1, center0));
		for(Peak p : peaks){
			final long copies=Math.round(p.center*mult);
			long size=p.volume*copies;
			long gc=p.gc*copies;
			sizeSum+=size;
			gcSum+=gc;
		}
		return (gcSum*1.0)/(sizeSum*k);
	}
	
	/**
	 * Alternative GC content calculation using complete histogram data.
	 * Analyzes GC content across all coverage levels starting from first peak.
	 *
	 * @param peaks List of peaks for reference
	 * @param k K-mer size for normalization
	 * @param hist Main frequency histogram
	 * @param gcHist GC content histogram
	 * @return Estimated GC content as fraction (0.0 to 1.0)
	 */
	private static double gcContent2(ArrayList<Peak> peaks, int k, long[] hist, long[] gcHist){
		if(peaks.size()<1){return 0;}
		
		long sizeSum=0, gcSum=0;
		final Peak p0=peaks.get(0);
		final int center0=p0.center;
		final double mult=1.0/(Tools.max(1, center0));
		
		for(int i=p0.start; i<hist.length; i++){
			final long copies=(Tools.max(1, Math.round(i*mult)));
			long size=hist[i]*copies;
			long gc=gcHist[i]*copies;
			sizeSum+=size;
			gcSum+=gc;
		}
		return (gcSum*1.0)/(sizeSum*k);
	}
	
	/** Based on all peaks */
	private static long repeatSize(ArrayList<Peak> peaks, final int ploidy, final float haploidPeakCenter){
		if(peaks.size()<2){return 0;}
		assert(ploidy>0) : ploidy;
		final int homozygousLoc=homozygousPeak(peaks, ploidy, haploidPeakCenter);
		final double mult=1.0/(Tools.max(1, haploidPeakCenter));
		
		long sizeSum=0;
		for(int i=homozygousLoc+1; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			long size=p.volume*((Math.round(p.center*mult)-1));
			sizeSum+=size;
		}
		return sizeSum;
	}
	
	/** Based on ploidy peak only */
	private static long repeatSize2(ArrayList<Peak> peaks, final int ploidy, final float haploidPeakCenter, long[] hist){
		assert(ploidy>0) : ploidy;
		final double mult=1.0/(Tools.max(1, haploidPeakCenter));
		int valley=(int)Math.ceil(haploidPeakCenter*ploidy*(1.2f+1.0f/Tools.max(2,ploidy)));
		
		//Optional
		final int homozygousLoc=homozygousPeak(peaks, ploidy, haploidPeakCenter);
		if(ploidy>1 && homozygousLoc>=0) {
			Peak p=peaks.get(homozygousLoc);
			valley=p.stop+1;
		}
		
		long sizeSum=0;
		for(int i=valley; i<hist.length; i++){
			long size=hist[i]*(Math.round(i*mult)-1);
			sizeSum+=size;
		}
		return sizeSum;
	}
	
//	/** Based on primary peak only */
//	private static long repeatSize2(ArrayList<Peak> peaks, final float haploidPeakCenter, final int ploidy, long[] hist){
//		assert(ploidy>0) : ploidy;
////		final int homozygousLoc=homozygousPeak(peaks, haploidPeakCenter, ploidy);
//		final Peak p0=peaks.get(0);
//		final int center0=p0.center;
//		final double mult=1.0/(Tools.max(1, center0));
//		
//		
//		final int valley=(int)Math.ceil(center0*1.7f);
//		long sizeSum=0;
//		for(int i=valley; i<hist.length; i++){
//			long size=hist[i]*(Math.round(i*mult)-1);
//			sizeSum+=size;
//		}
//		return sizeSum;
//	}
	
	/**
	 * Finds index of peak with largest volume.
	 * @param peaks List of peaks to search
	 * @return Index of peak with maximum volume
	 */
	private static int biggestPeak(ArrayList<Peak> peaks){
		if(peaks.size()<2){return peaks.size()-1;}

		final Peak p0=peaks.get(0);
		Peak biggest=p0;
		int loc=0;
		for(int i=1; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			if(p.volume>biggest.volume){
				loc=i;
				biggest=p;
			}
		}
		return loc;
	}
	
	/**
	 * Finds index of peak with second-largest volume.
	 * @param peaks List of peaks to search
	 * @return Index of peak with second-largest volume
	 */
	private static int secondBiggestPeak(ArrayList<Peak> peaks){
		if(peaks.size()<2){return peaks.size()-1;}
		
		Peak biggest=peaks.get(0);
		Peak second=peaks.get(1);
		int bloc=0;
		int sloc=1;
		if(second.volume>biggest.volume){
			Peak temp=second;
			second=biggest;
			biggest=temp;
			bloc=1;
			sloc=0;
		}
		for(int i=2; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			if(p.volume>second.volume){
				sloc=i;
				second=p;
				if(second.volume>biggest.volume){
					Peak temp=second;
					second=biggest;
					biggest=temp;
					sloc=bloc;
					bloc=i;
				}
			}
		}
		return sloc;
	}
	
	/**
	 * Identifies peak corresponding to homozygous coverage level.
	 * Finds peak closest to expected homozygous coverage (ploidy * haploid).
	 *
	 * @param peaks List of peaks to search
	 * @param ploidy Genome ploidy level
	 * @param haploidPeakCenter Coverage level of haploid peak
	 * @return Index of peak closest to homozygous coverage
	 */
	private static int homozygousPeak(ArrayList<Peak> peaks, final int ploidy, final float haploidPeakCenter){
		if(peaks.size()<2){return peaks.size()-1;}
		assert(ploidy>0) : ploidy;

		final float target=haploidPeakCenter*ploidy;
//		System.err.println("target="+target);
		float bestDif=Integer.MAX_VALUE;
		int loc=0;
		for(int i=0; i<peaks.size(); i++){
			Peak p=peaks.get(i);
			float dif=Tools.absdif(target, p.center);
//			System.err.println("dif="+dif+" for peak "+i+", center "+p.center);
			if(dif<bestDif){
				bestDif=dif;
				loc=i;
//				System.err.println("New best at loc "+i);
			}
		}
		return loc;
	}
	
	/**
	 * Determines haploid peak coverage level from peak analysis.
	 * Uses peak volume ratios to distinguish between haploid and ploidy peaks.
	 *
	 * @param peaks List of called peaks
	 * @param ploidy Expected genome ploidy
	 * @return Estimated coverage level of haploid peak
	 */
	private static float haploidPeakCenter(ArrayList<Peak> peaks, int ploidy) {
		assert(ploidy>0);
		final Peak biggest=peaks.get(biggestPeak(peaks));
		final Peak second=peaks.get(secondBiggestPeak(peaks));
		if(second.volume*4>=biggest.volume) {//Similar volume; lowest is haploid peak
			return Tools.min(biggest.center, second.center);
		}
		//Assume biggest is ploidy peak
		return biggest.center/(float)ploidy;
	}
	
	/**
	 * Estimates genome ploidy from peak coverage pattern.
	 * Analyzes peak positions and volumes to determine most likely ploidy level.
	 * Supports different logic modes for ploidy estimation.
	 *
	 * @param peaks List of called peaks
	 * @param minVolumeFraction Minimum volume fraction for significant peaks
	 * @return Estimated ploidy level
	 */
	private static int calcPloidy(ArrayList<Peak> peaks, float minVolumeFraction){
		if(peaks.size()<2){return 1;}

		final Peak p0=peaks.get(0);
		final Peak biggest=peaks.get(biggestPeak(peaks));
		final Peak second=peaks.get(secondBiggestPeak(peaks));
		
		if(ploidyLogic==1) {//Old code path, assumes no error/contamination peaks
			if(biggest==p0){//p0 is biggest.
				//			System.err.println("a: "+biggest+" > "+second);
				if(second.volume*4<biggest.volume){return 1;} //Probably second is a repeat peak.
				int ratio=Math.round((second.center/(float)biggest.center));
				//			System.err.println("b: "+ratio);
				return Tools.max(1, ratio);
			}else {//p0 is not biggest.
				//I wonder if this should be biggest/second...
				int ratio=Math.round((biggest.center/(float)p0.center));
				//			System.err.println("c: "+biggest+", "+((biggest.center/(float)p0.center)));
				return ratio;
			}
		}else {

			if(second==biggest) {
				return 1;
			}else if(second.center<biggest.center) {//Second volume could be low, and is based on het rate
				if(second.volume<biggest.volume*minVolumeFraction) {
					//Second peak is tiny and likely contamination.
					return 1;
				}
			}else {//second volume should be substantial since it is the heterozygous fraction
				if(second.volume*4<biggest.volume) {
					//Second peak is small and likely repeat.
					return 1;
				}
			}
			final int max=Tools.max(biggest.center, second.center);
			final int min=Tools.min(biggest.center, second.center);
			return Math.round(max/(float)min);
		}
	}
	
	/**
	 * Calculates number of heterozygous genomic positions.
	 * Sums volumes of sub-ploidy peaks to estimate heterozygous sites.
	 *
	 * @param peaks List of called peaks
	 * @param ploidy Genome ploidy level
	 * @param haploidPeakCenter Coverage level of haploid peak
	 * @param k K-mer size for position calculation
	 * @return Estimated number of heterozygous positions
	 */
	private static long calcHetLocations(ArrayList<Peak> peaks, final int ploidy, final float haploidPeakCenter, final int k){
		if(peaks.size()<2){return 0;}
		assert(ploidy>0) : ploidy;
		final int homozygousLoc=homozygousPeak(peaks, ploidy, haploidPeakCenter);
		final Peak homoPeak=peaks.get(homozygousLoc);
		long sum=0;
		final int lim=ploidy/2;
		for(int i=0; i<homozygousLoc; i++){
			final Peak p=peaks.get(i);
			final int copyCount=Math.round((p.center*ploidy)/(float)homoPeak.center);
//			System.err.println("lim="+lim+". For peak "+i+", copyCount="+copyCount+", volume="+p.volume);
			if(copyCount>lim){break;}
//			double peakLocs=(p.volume/(double)k);
			sum=sum+p.volume;
		}
		return sum/k;
	}
	
	/**
	 * Calls peaks from LongList histogram data.
	 * Convenience wrapper for callPeaks method with arrays.
	 *
	 * @param list Main k-mer frequency histogram
	 * @param gcList GC content histogram (may be null)
	 * @return List of identified peaks
	 */
	public ArrayList<Peak> callPeaks(LongList list, LongList gcList){
		return callPeaks(list.array, gcList==null ? null : gcList.array, list.size);
	}
	
	/**
	 * Main peak calling algorithm for k-mer frequency histograms.
	 * Applies optional log scaling and smoothing, then identifies peaks
	 * using state machine approach. Filters peaks by size and volume thresholds.
	 *
	 * @param original Input histogram array
	 * @param gcArray GC content histogram (may be null)
	 * @param length Active length of histogram arrays
	 * @return List of identified peaks meeting filtering criteria
	 */
	public ArrayList<Peak> callPeaks(final long[] original, final long gcArray[], final int length){
		
		long[] array=original;
		
		if(logScale){
			array=logScale(array, logWidth, 1, logPasses);
		}
		
		if(smoothRadius>0){
			if(smoothProgressiveFlag){
				array=smoothProgressive(array, smoothRadius);
			}else{
				array=smooth(array, smoothRadius);
			}
		}
		
		ArrayList<Peak> peaks=new ArrayList<Peak>();
		
		int dip0=-1;
		for(int i=1; i<length; i++){
			if(array[i-1]<array[i]){
				dip0=i-1;
				break;
			}
		}
		if(dip0<0){return peaks;}
//		assert(false) : dip0+", "+array[dip0);
		
		final int UP=0, DOWN=1;
		int mode=UP;
		int start=dip0, center=-1;
		long prev=array[dip0];
		long sum=prev;
		long gcSum=gcArray==null ? -1 : gcArray[dip0];
		long sum2=prev*dip0;
		for(int i=dip0+1; i<length; i++){
			final long x=array[i];
			
//			if(i<16){System.err.println("i="+i+", x="+x+", mode="+mode+", center="+center+", start="+start+", dip0="+dip0);}
			
			if(mode==UP){
				if(x<prev){
					mode=DOWN;
					center=i-1;
				}
			}else{
				if(x>prev){
					mode=UP;
					int stop=i-1;
					long max=array[center];
					if(center>=minPeak && center<=maxPeak && max>=minHeight && (stop-start)>=minWidth && sum>=minVolume){
						for(int j=center-1; j>=0; j--){//find middle of mesas
							if(array[j]!=max){
								center=(center+j+2)/2;
								break;
							}
						}
						{
							long valley=array[stop];
							for(int j=stop; j>=0; j--){//find middle of valleys
								if(array[j]!=valley){
									if(valley==0){stop=j+1;}
									else{stop=(stop+j+2)/2;}
									break;
								}
							}
						}
						
						long height1=array[start], height2=array[stop];
						Peak p=new Peak(center, start, stop, center, max, height1, height2, height1, height2, sum, sum2, gcSum);
						peaks.add(p);
					}else{
//						Peak p=new Peak(center, start, stop, max, sum);
//						System.err.println("*"+p);
					}
					start=stop;
					stop=-1;
					sum=sum2=0;
					gcSum=0;
					center=-1;
					if(i>maxPeak){break;}
					while(i<array.length && array[i]==0){i++;}//Skip zero regions
				}
			}
			gcSum=(gcArray==null ? -1 : gcSum+gcArray[i]);
			sum+=x;
			sum2+=(x*i);
			prev=x;
		}
		
		if(mode==DOWN){
			int stop=length;
			long max=array[center];
			for(int j=center-1; j>=0; j--){//find middle of mesas
				if(array[j]!=max){
					center=(center+j+2)/2;
					break;
				}
			}
			{
				long valley=array[stop-1];
				for(int j=stop-1; j>=0; j--){//find middle of valleys
					if(array[j]!=valley){
						if(valley==0){stop=j+1;}
						else{stop=(stop+j+2)/2;}
						break;
					}
				}
			}
			if(center>=minPeak && center<=maxPeak && max>=minHeight && (stop-start)>=minWidth && sum>=minVolume){
				long height1=array[start], height2=array[Tools.min(stop, length-1)];
				Peak p=new Peak(center, start, Tools.min(stop, length-1), center, max, height1, height2, height1, height2, sum, sum2, gcSum);
				peaks.add(p);
			}else{
//				Peak p=new Peak(center, start, stop, max, sum);
//				System.err.println("*"+p);
			}
		}
		
		//Ensure peaks don't violate maxWidthMult by shrinking those that do, for each side, relative to the center.
		capWidth(peaks, maxWidthMult, array);
		
		//Have peaks absorb adjacent smaller peaks
		if(maxPeakCount<peaks.size()){
			peaks=condense(peaks, maxPeakCount);
		}
		
//		System.err.println("start\tcenter\tstop\tmax\tvolume");
//		for(Peak p : peaks){
//			System.err.println(p.toString());
//		}
//		assert(false);
		
		capWidth(peaks, maxWidthMult, array);
		
		if(peaks.size()>1){
			Peak biggest=peaks.get(biggestPeak(peaks));
			while(peaks.size()>1 && peaks.get(0).volume<0.0001*biggest.volume){
				peaks.remove(0);
			}
		}
		
		ArrayList<Peak> peaksOut=peaks;
		if(array!=original){
			peaksOut=new ArrayList<Peak>();
			recalculate(peaks, original);
			for(Peak p : peaks){
				if(p.volume>=minVolume){peaksOut.add(p);}
			}
		}
		
		return peaksOut;
	}
	
	/**
	 * Recalculates peak statistics using provided histogram array.
	 * Updates peak volumes, heights, and positions based on new data.
	 * @param peaks List of peaks to recalculate
	 * @param array Histogram array for recalculation
	 */
	private static void recalculate(ArrayList<Peak> peaks, long[] array){
		for(Peak p : peaks){
			p.recalculate(array);
		}
	}
	
	/**
	 * Applies logarithmic scaling to smooth histogram data.
	 * Uses variable-width sliding window that increases with position.
	 * Multiple passes can be applied for stronger smoothing effect.
	 *
	 * @param array Input histogram to scale
	 * @param width Scaling window width parameter
	 * @param scale Scaling factor for output values
	 * @param passes Number of scaling passes to apply
	 * @return Log-scaled histogram array
	 */
	public static long[] logScale(long[] array, double width, double scale, int passes){
		assert(passes>0);
		long[] log=array;
		for(int pass=0; pass<passes; pass++){
			final double halfWidth=width/2;
			final double limit=array.length-0.00001;
			log=new long[array.length];
			for(int pos=1; pos<array.length; pos++){
				final double center=pos+0.5;
				final double min=Tools.max(0, center-halfWidth*pos);
				final double max=Tools.min(limit, center+halfWidth*pos);
				final int mini=(int)min;
				final int maxi=(int)max;
				if(mini==maxi){
					log[pos]=Math.round((max-min)*array[mini]*scale);
				}else{
					double sum=0;
					sum+=(array[mini]*(mini+1-min));
					sum+=(array[maxi]*(max-maxi));
					for(int i=mini+1; i<maxi; i++){
						sum+=array[i];
					}
					log[pos]=Math.round(sum*scale);
				}
			}
			array=log;
		}
		return log;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Smoothing          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Applies progressive smoothing with increasing radius.
	 * Smoothing radius grows larger at higher positions to handle
	 * varying peak widths across coverage spectrum.
	 *
	 * @param data Input histogram array
	 * @param radius0 Initial smoothing radius
	 * @return Progressively smoothed histogram
	 */
	public static long[] smoothProgressive(final long[] data, int radius0){
		int radius=radius0;
		long div=radius*radius;
		double mult=1.0/div;
		long[] smoothed=new long[data.length];
		for(int i=0, next=5; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(long)Math.ceil(product);}
//			else{smoothed[i]=(long)product;}
			smoothed[i]=Math.round(product);
			if(i>next){
				next=(int)Math.ceil(1+next*progressiveMult);
				radius+=1;
				div=radius*radius;
				mult=1.0/div;
				if(radius>maxRadius){next=Integer.MAX_VALUE;}
//				System.err.println(radius+", "+div);
			}
		}
		return smoothed;
	}
	
	/**
	 * Applies uniform smoothing with fixed radius.
	 * Uses weighted averaging within specified radius around each position.
	 *
	 * @param data Input histogram array
	 * @param radius Smoothing radius in array positions
	 * @return Smoothed histogram array
	 */
	public static long[] smooth(final long[] data, int radius){
		final long div=radius*radius;
		final double mult=1.0/div;
		long[] smoothed=new long[data.length];
		for(int i=0; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(long)Math.ceil(product);}
//			else{smoothed[i]=(long)product;}
			smoothed[i]=Math.round(product);
		}
		return smoothed;
	}
	
	/**
	 * Calculates weighted sum around specified position for smoothing.
	 * Uses triangular weighting scheme with maximum weight at center.
	 *
	 * @param data Histogram array
	 * @param loc Center position for calculation
	 * @param radius Radius for weighted sum
	 * @return Weighted sum for smoothing
	 */
	private static long sumPoint(long[] data, int loc, int radius){
		long sum=0;
		int start=loc-radius+1;
		int stop=loc+radius-1;
		for(int i=start, x=1; i<loc; i++, x++){
			int i2=Tools.max(i, 0);
			sum+=data[i2]*x;
		}
		for(int i=loc, x=radius, max=data.length-1; i<=stop; i++, x--){
			int i2=Tools.min(i, max);
			sum+=data[i2]*x;
		}
		return sum;
	}
	
	/**
	 * Progressive smoothing for integer arrays.
	 * Integer version of progressive smoothing with increasing radius.
	 *
	 * @param data Input integer histogram
	 * @param radius0 Initial smoothing radius
	 * @return Progressively smoothed integer histogram
	 */
	public static int[] smoothProgressive(final int[] data, int radius0){
		int radius=radius0;
		long div=radius*radius;
		double mult=1.0/div;
		int[] smoothed=new int[data.length];
		for(int i=0, next=5; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(long)Math.ceil(product);}
//			else{smoothed[i]=(long)product;}
			smoothed[i]=(int)Math.round(product);
			if(i>next){
				next=(int)Math.ceil(next*2);
				radius+=1;
				div=radius*radius;
				mult=1.0/div;
				if(radius>10){next=Integer.MAX_VALUE;}
//				System.err.println(radius+", "+div);
			}
		}
		return smoothed;
	}
	
	/**
	 * Uniform smoothing for integer arrays.
	 * Integer version of fixed-radius smoothing algorithm.
	 *
	 * @param data Input integer histogram
	 * @param radius Smoothing radius
	 * @return Smoothed integer histogram
	 */
	public static int[] smooth(final int[] data, int radius){
		final long div=radius*radius;
		final double mult=1.0/div;
		int[] smoothed=new int[data.length];
		for(int i=0; i<data.length; i++){
			long sum=sumPoint(data, i, radius);
			double product=sum*mult;
//			if(data[i]>=product){smoothed[i]=(int)Math.ceil(product);}
//			else{smoothed[i]=(int)product;}
			smoothed[i]=(int)Math.round(product);
		}
		return smoothed;
	}
	
	/**
	 * Weighted sum calculation for integer arrays.
	 * Integer array version of weighted sum for smoothing.
	 *
	 * @param data Integer histogram array
	 * @param loc Center position
	 * @param radius Calculation radius
	 * @return Weighted sum for smoothing
	 */
	private static long sumPoint(int[] data, int loc, int radius){
		long sum=0;
		int start=loc-radius+1;
		int stop=loc+radius-1;
		for(int i=start, x=1; i<loc; i++, x++){
			int i2=Tools.max(i, 0);
			sum+=data[i2]*x;
		}
		for(int i=loc, x=radius, max=data.length-1; i<=stop; i++, x--){
			int i2=Tools.min(i, max);
			sum+=data[i2]*x;
		}
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Represents a coverage peak in k-mer frequency histogram.
	 * Stores peak boundaries, height, volume, and GC content statistics.
	 * Provides methods for peak manipulation and output formatting.
	 */
	private class Peak{
		
		/**
		 * Constructs Peak with complete statistics.
		 * Initializes all peak parameters including boundaries, heights, and volumes.
		 *
		 * @param center_ Peak center position
		 * @param start_ Peak start position (inclusive)
		 * @param stop_ Peak stop position (exclusive)
		 * @param maxPos_ Position of maximum height within peak
		 * @param maxHeight_ Maximum height value in peak
		 * @param startHeight_ Height at start position
		 * @param stopHeight_ Height at stop position
		 * @param leftMin_ Minimum height in left half
		 * @param rightMin_ Minimum height in right half
		 * @param volume_ Total peak volume (sum of heights)
		 * @param volume2_ Weighted peak volume (sum of heights * positions)
		 * @param gc_ GC content sum for this peak
		 */
		Peak(int center_, int start_, int stop_, int maxPos_, long maxHeight_, long startHeight_, long stopHeight_,
				long leftMin_, long rightMin_, long volume_, long volume2_, long gc_){
			
			center=center_;
			start=Tools.max(0, start_);
			stop=stop_;
			maxPos=maxPos_;
			
			maxHeight=maxHeight_;
			startHeight=startHeight_;
			stopHeight=stopHeight_;
			volume=volume_;
			volume2=volume2_;
			
			leftMin=leftMin_;
			rightMin=rightMin_;
			
			gc=gc_;
			
			assert(center>=0) : this;
			assert(start<center) : this;
			assert(stop>center) : this;
		}
		
		/**
		 * Determines if this peak is compatible for merging with another.
		 * Uses maximum width multiplier to check if peaks can be combined.
		 * @param p Peak to test compatibility with
		 * @return true if peaks are compatible for merging
		 */
		public boolean compatibleWith(Peak p) {
			int min=Tools.min(center, p.stop);
			int max=Tools.max(stop, p.center);
//			assert(min*maxWidthMult>=max) : this+", "+p+", "+(min*maxWidthMult)+", "+max;
			return min*maxWidthMult>=max;
		}

		/**
		 * @param array
		 */
		public void recalculate(long[] array) {
			maxHeight=array[center];
			startHeight=array[start];
			stopHeight=array[stop];
			leftMin=startHeight;
			rightMin=stopHeight;
			maxPos=center;
			volume=0;
			volume2=0;
			for(int i=start; i<stop; i++){
				long x=array[i];
				if(x>maxHeight){
					maxPos=i;
					maxHeight=x;
				}
				if(i<center){leftMin=Tools.min(leftMin,  x);}
				else if(i>center){rightMin=Tools.min(rightMin, x);}
				volume+=x;
				volume2+=(x*i);
			}
		}

		/**
		 * @param p
		 */
		public void absorb(Peak p) {
			assert(this!=p);
			
			if(center>p.center){
				assert(p.stop<stop) : "\n"+this+"\n"+p+"\n";
				if(start>p.start){
					start=p.start;
					startHeight=p.startHeight;
				}
				leftMin=Tools.min(leftMin, p.leftMin);
			}else{
				assert(p.start>start) : "\n"+this+"\n"+p+"\n";
				if(stop<p.stop){
					stop=p.stop;
					stopHeight=p.stopHeight;
				}
				rightMin=Tools.min(rightMin, p.rightMin);
			}
			
			//Potentially shift the center
//			{
//			long c1=callMetric();
//			long c2=p.callMetric();
//			System.err.print(this+" absorbed "+p);
//			if(c1<c2){
//				assert(false);
//				center=p.center;
//			}else{
//				
//			}
//			}
			
			if(maxHeight<p.maxHeight){
				maxHeight=p.maxHeight;
				maxPos=p.maxPos;
			}
			
			volume+=p.volume;
			volume2+=p.volume2;
			gc+=p.gc;
//			System.err.println(" -> "+this);
		}

		/** Returns peak width as difference between stop and start positions */
		int width(){return stop-start;}
		
		@Override
		public String toString(){
			return start+"\t"+center+"\t"+stop+"\t"+maxHeight+"\t"+volume;
		}
		
		/**
		 * Formats peak data for output as tab-delimited text.
		 * Appends start, center, stop, max height, and volume to ByteBuilder.
		 * @param bb ByteBuilder for output (created if null)
		 * @return ByteBuilder with formatted peak data
		 */
		public ByteBuilder toBytes(ByteBuilder bb){
			if(bb==null){bb=new ByteBuilder();}
			bb.append(start);
			bb.tab();
			bb.append(center);
			bb.tab();
			bb.append(stop);
			bb.tab();
			bb.append(maxHeight);
			bb.tab();
			bb.append(volume);
			return bb;
		}

		/** Inclusive */
		public int start;
		/** Peak center position */
		public int center;
		/** Exclusive */
		public int stop;
		/** Position of maximum height within peak */
		public int maxPos;
		
		//Unique counts
		/** Height at peak start position */
		public long startHeight;
//		public long centerHeight;
		/** Height at peak stop position */
		public long stopHeight;
		/** Total peak volume (sum of all heights in peak region) */
		public long volume;
		/** Weighted peak volume (sum of heights multiplied by positions) */
		public long volume2;
		
		/** Minimum height in left half of peak */
		public long leftMin;
		/** Minimum height in right half of peak */
		public long rightMin;
		/** Maximum height within peak region */
		public long maxHeight;
		/** Sum of GC content values for this peak */
		public long gc;

		//Raw counts
		/** Returns raw start height multiplied by start position */
		public long startHeight2(){return startHeight*start;}
//		public long centerHeight2(){return centerHeight*center;}
		/** Returns raw max height multiplied by center position */
		public long maxHeight2(){return maxHeight*center;}
		/** Returns raw stop height multiplied by stop position */
		public long stopHeight2(){return stopHeight*stop;}
		/**
		 * Returns the metric used for peak calling decisions.
		 * Choice between volume/height and raw/normalized based on call mode settings.
		 * @return Peak metric value for comparison
		 */
		public long callMetric(){
			if(CALL_MODE==BY_VOLUME){return callByRawCount ? volume2 : volume;}
			return callByRawCount ? maxHeight2() : maxHeight;
		}
		
		
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum peak height threshold for peak detection */
	private long minHeight=2;
	/** Minimum peak volume threshold for peak detection */
	private long minVolume=5;
	/** Minimum peak width in positions */
	private int minWidth=3;
	/** Minimum position to consider for peak detection */
	private int minPeak=2;
	/** Maximum position to consider for peak detection */
	private int maxPeak=Integer.MAX_VALUE;
	/** Maximum number of peaks to report */
	private int maxPeakCount=10;
	/** Maximum peak width as multiplier of center position */
	private float maxWidthMult=2.5f;
	/** Radius for histogram smoothing operations */
	private int smoothRadius=0;
	/** Whether to use progressive smoothing with increasing radius */
	private boolean smoothProgressiveFlag=true;
	/** K-mer size used for histogram generation */
	private int k=31;
	
	/** User-specified ploidy level (overrides automatic detection) */
	private int ploidyClaimed=-1;
	
	/** Whether to apply logarithmic scaling to histogram */
	private boolean logScale=false;
	/** Window width parameter for logarithmic scaling */
	private double logWidth=0.1;
	/** Number of logarithmic scaling passes to apply */
	private int logPasses=1;
	
	/** Input file path for histogram data */
	private String in;
	/** Output file path for peak analysis results */
	private String out;
	
	/** Input file format specification */
	private final FileFormat ffin;
	/** Output file format specification */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maximum radius allowed for progressive smoothing */
	public static int maxRadius=10;
	/** Multiplier for radius growth in progressive smoothing */
	public static float progressiveMult=2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Column index containing count values in input histogram */
	private static int countColumn=1;
	
	/** Output stream for error messages and diagnostics */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and detailed logging */
	public static boolean verbose=false;
	/** Flag indicating whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	/** Whether to print class name in help and error messages */
	public static boolean printClass=true;

	/** Whether to use raw counts instead of normalized counts for peak calling */
	public static boolean callByRawCount=false;
	/** Whether to weight peaks by relief (difference from surrounding valleys) */
	public static boolean weightByRelief=false;
	/** Constant for volume-based peak calling mode */
	public static final int BY_VOLUME=0, BY_HEIGHT=1;
	/** Current peak calling mode (BY_VOLUME or BY_HEIGHT) */
	public static int CALL_MODE=BY_VOLUME;
	/** Logic mode for ploidy detection algorithm */
	public static int ploidyLogic=2;
	
}
