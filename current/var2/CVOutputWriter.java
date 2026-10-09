package var2;

import fileIO.FileFormat;
import fileIO.TextStreamWriter;
import shared.Timer;
import shared.Tools;

/**
 * Emits CallVariants variant files and histogram summaries with timing output.
 * Variant outputs share one sorted VcfWriter, in VAR/VCF/GFF order; histogram
 * writers are created per file. Inputs are borrowed and must remain stable while
 * the methods run. Variant outputs throw after the affected writer reports
 * errorState; histogram helpers accumulate returned writer error flags across
 * requested files and throw after the last helper returns.
 * @author Brian Bushnell
 */
public class CVOutputWriter{

	/**
	 * Writes selected variant outputs with composite scoring (no supplied network).
	 * @param varMap Borrowed quiescent variants and dataset statistics
	 * @param varFilter Nonnull output thresholds
	 * @param ffout Optional ordered VAR descriptor
	 * @param vcf Optional VCF output path
	 * @param gffout Optional GFF output path
	 * @param readsProcessed Reported read count for headers
	 * @param pairedInSequencingReadsProcessed Reported paired-read count
	 * @param properlyPairedReadsProcessed Reported properly paired-read count
	 * @param trimmedBasesProcessed Reported base count for headers
	 * @param ref Optional reference path printed in headers
	 * @param trimWhitespace Trim VCF scaffold names
	 * @param sampleName VCF sample label
	 */
	public static void writeOutput(VarMap varMap, VarFilter varFilter, FileFormat ffout, String vcf, String gffout,
			long readsProcessed, long pairedInSequencingReadsProcessed, long properlyPairedReadsProcessed,
			long trimmedBasesProcessed, String ref, boolean trimWhitespace, String sampleName){
		writeOutput(varMap, varFilter, ffout, vcf, gffout, readsProcessed,
				pairedInSequencingReadsProcessed, properlyPairedReadsProcessed,
				trimmedBasesProcessed, ref, trimWhitespace, sampleName, null);
	}

	/**
	 * Creates one sorted writer and sequentially emits each requested format.
	 * No output means no writer allocation. The borrowed network is copied by each
	 * VcfWriter formatter. Recorded variant-output errors stop this method with an exception.
	 * @param net Optional scoring network; other parameters match the no-network overload
	 */
	public static void writeOutput(VarMap varMap, VarFilter varFilter, FileFormat ffout, String vcf, String gffout,
			long readsProcessed, long pairedInSequencingReadsProcessed, long properlyPairedReadsProcessed,
			long trimmedBasesProcessed, String ref, boolean trimWhitespace, String sampleName, ml.CellNet net){

		if(ffout!=null || vcf!=null || gffout!=null){
			Timer t3=new Timer("Sorting variants.");
			VcfWriter vw=new VcfWriter(varMap, varFilter, readsProcessed,
					pairedInSequencingReadsProcessed, properlyPairedReadsProcessed,
						trimmedBasesProcessed, ref, trimWhitespace, sampleName);
			vw.setNet(net);
			t3.stop("Time: ");

			if(ffout!=null){
				t3.start("Writing Var file.");
				vw.writeVarFile(ffout);
				if(vw.errorState){throw new RuntimeException("Failed to write VAR output: "+ffout.name());}
				t3.stop("Time: ");
			}
			if(vcf!=null){
				t3.start("Writing VCF file.");
				vw.writeVcfFile(vcf);
				if(vw.errorState){throw new RuntimeException("Failed to write VCF output: "+vcf);}
				t3.stop("Time: ");
			}
			if(gffout!=null){
				t3.start("Writing GFF file.");
				vw.writeGffFile(gffout);
				if(vw.errorState){throw new RuntimeException("Failed to write GFF output: "+gffout);}
				t3.stop("Time: ");
			}
		}
	}

	/**
	 * Writes only requested histograms. Score and average-quality counts use row 0;
	 * zygosity and maximum-quality counts are one-dimensional. Arrays for skipped
	 * outputs may be null; requested outputs require the helper's array preconditions.
	 * Returned writer error flags are accumulated and reported after the last
	 * requested helper call returns.
	 * @param scoreHistFile Optional score output path
	 * @param zygosityHistFile Optional zygosity output path
	 * @param qualityHistFile Optional quality output path
	 * @param scoreArray Score histogram with a nonnull row 0 when requested
	 * @param ploidyArray Nonempty zygosity histogram when requested
	 * @param avgQualityArray Average-quality histogram with a nonnull row 0 when requested
	 * @param maxQualityArray Maximum-quality histogram at least as long as average row 0
	 */
	public static void writeHistograms(String scoreHistFile, String zygosityHistFile, String qualityHistFile,
			long[][] scoreArray, long[] ploidyArray, long[][] avgQualityArray, long[] maxQualityArray){

		if(scoreHistFile!=null || zygosityHistFile!=null || qualityHistFile!=null){
			Timer t3=new Timer("Writing histograms.");
			boolean errorState=false;
			if(scoreHistFile!=null){
				errorState=writeScoreHist(scoreHistFile, scoreArray[0], true)|errorState;
			}
			if(zygosityHistFile!=null){
				errorState=writeZygosityHist(zygosityHistFile, ploidyArray)|errorState;
			}
			if(qualityHistFile!=null){
				errorState=writeQualityHist(qualityHistFile, avgQualityArray[0], maxQualityArray)|errorState;
			}
			t3.stop("Time: ");
			if(errorState){throw new RuntimeException("Failed to write histogram output.");}
		}
	}

	/** Writes a score histogram whose terminal bin reports overflow when populated. */
	static boolean writeScoreHist(String fname, long[] array){
		return writeScoreHist(fname, array, true);
	}

	/**
	 * Writes count-weighted score statistics and rows through the last nonzero bin.
	 * Mean is NaN for zero total count; median uses Tools histogram conventions.
	 * Overflow makes mean and median lower bounds; the mode is a histogram bucket.
	 * Opens an overwriting, non-appending TextStreamWriter and waits for completion.
	 * @param fname Output path
	 * @param array Nonnull score histogram; index is score and value is count
	 * @param overflow True if the terminal array index is an overflow bin
	 * @return True if the writer recorded an error, false otherwise
	 */
	static boolean writeScoreHist(String fname, long[] array, boolean overflow){
		int max=array.length-1;
		for(; max>=0; max--){
			if(array[max]!=0){break;}
		}
		overflow&=array[array.length-1]>0;
		long sum=0, sum2=0;
		for(int i=0; i<=max; i++){
			sum+=array[i];
			sum2+=(i*array[i]);
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, true, false, false);
		tsw.start();
		tsw.println("#ScoreHist");
		tsw.println("#Vars\t"+sum);
		tsw.println((overflow ? "#MeanLowerBound\t" : "#Mean\t")+Tools.format("%.2f", sum2*1.0/sum));
		tsw.println((overflow ? "#MedianLowerBound\t" : "#Median\t")+Tools.medianHistogram(array));
		tsw.println((overflow ? "#ModeBucket\t" : "#Mode\t")+scoreBucketName(Tools.calcModeHistogram(array), array.length, overflow));
		tsw.println("#Quality\tCount");
		for(int i=0; i<=max; i++){
			tsw.println(scoreBucketName(i, array.length, overflow)+"\t"+array[i]);
		}
		tsw.poisonAndWait();
		return tsw.errorState;
	}

	/** Returns the displayed name for score histogram bucket {@code i}. */
	private static String scoreBucketName(int i, int length, boolean overflow){
		return overflow && i==length-1 ? ">="+i : Integer.toString(i);
	}

	/**
	 * Writes all zygosity bins and their count-weighted mean. The final bin supplies
	 * the reported homozygous fraction. Zero total count produces NaN ratios.
	 * @param fname Output path, overwritten without append
	 * @param array Nonempty histogram indexed by zygosity
	 * @return True if the writer recorded an error, false otherwise
	 */
	static boolean writeZygosityHist(String fname, long[] array){
		int max=array.length-1;
		long sum=0, sum2=0;
		for(int i=0; i<=max; i++){
			sum+=array[i];
			sum2+=(i*array[i]);
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, true, false, false);
		tsw.start();
		tsw.println("#ZygoHist");
		tsw.println("#Vars\t"+sum);
		tsw.println("#Mean\t"+Tools.format("%.3f", sum2*1.0/sum));
		tsw.println("#HomozygousFraction\t"+Tools.format("%.3f", array[max]*1.0/sum));
		tsw.println("#Zygosity\tCount");
		for(int i=0; i<=max; i++){
			tsw.println(i+"\t"+array[i]);
		}
		tsw.poisonAndWait();
		return tsw.errorState;
	}

	/**
	 * Writes paired average/maximum-quality rows through their last nonzero bin.
	 * Summary count, mean, median and mode use only the average-quality histogram.
	 * Zero average total count gives a NaN mean. Waits for writer completion.
	 * @param fname Output path, overwritten without append
	 * @param avgQualArray Nonnull average-quality histogram
	 * @param maxQualArray Nonnull maximum-quality histogram at least as long as avgQualArray
	 * @return True if the writer recorded an error, false otherwise
	 */
	static boolean writeQualityHist(String fname, long[] avgQualArray, long[] maxQualArray){
		int max=avgQualArray.length-1;
		for(; max>=0; max--){
			if(avgQualArray[max]!=0 || maxQualArray[max]!=0){break;}
		}
		long avgsum=0, avgsum2=0;
		for(int i=0; i<=max; i++){
			avgsum+=avgQualArray[i];
			avgsum2+=(i*avgQualArray[i]);
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, true, false, false);
		tsw.start();
		tsw.println("#BaseQualityHist");
		tsw.println("#Vars\t"+avgsum);
		tsw.println("#Mean\t"+Tools.format("%.2f", avgsum2*1.0/avgsum));
		tsw.println("#Median\t"+Tools.medianHistogram(avgQualArray));
		tsw.println("#Mode\t"+Tools.calcModeHistogram(avgQualArray));
		tsw.println("#Quality\tAvgCount\tMaxCount");
		for(int i=0; i<=max; i++){
			tsw.println(i+"\t"+avgQualArray[i]+"\t"+maxQualArray[i]);
		}
		tsw.poisonAndWait();
		return tsw.errorState;
	}
}
