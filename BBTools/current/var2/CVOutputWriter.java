package var2;

import fileIO.FileFormat;
import fileIO.TextStreamWriter;
import shared.Timer;
import shared.Tools;

/**
 * Utility class for writing CallVariants output files.
 * Contains static methods for generating various output formats.
 * 
 * @author Brian Bushnell
 * @author Isla Winglet
 * @date June 2025
 */
public class CVOutputWriter {
	
	/**
	 * Writes output files (Var, VCF, GFF) based on provided parameters.
	 * @param varMap Variant map containing all variants
	 * @param varFilter Variant filtering parameters
	 * @param ffout FileFormat for var file output
	 * @param vcf VCF output file path
	 * @param gffout GFF output file path
	 * @param readsProcessed Total reads processed (minus discarded)
	 * @param pairedInSequencingReadsProcessed Paired reads processed
	 * @param properlyPairedReadsProcessed Properly paired reads processed
	 * @param trimmedBasesProcessed Total trimmed bases processed
	 * @param ref Reference file path
	 * @param trimWhitespace Whether to trim whitespace in output
	 * @param sampleName Sample name for output files
	 */
	public static void writeOutput(VarMap varMap, VarFilter varFilter, FileFormat ffout, String vcf, String gffout,
			long readsProcessed, long pairedInSequencingReadsProcessed, long properlyPairedReadsProcessed,
			long trimmedBasesProcessed, String ref, boolean trimWhitespace, String sampleName) {
		
		if(ffout!=null || vcf!=null || gffout!=null){
			Timer t3=new Timer("Sorting variants.");
			VcfWriter vw=new VcfWriter(varMap, varFilter, readsProcessed, 
					pairedInSequencingReadsProcessed, properlyPairedReadsProcessed,
						trimmedBasesProcessed, ref, trimWhitespace, sampleName);
			t3.stop("Time: ");
			
			if(ffout!=null){
				t3.start("Writing Var file.");
				vw.writeVarFile(ffout);
				t3.stop("Time: ");
			}
			if(vcf!=null){
				t3.start("Writing VCF file.");
				vw.writeVcfFile(vcf);
				t3.stop("Time: ");
			}
			if(gffout!=null){
				t3.start("Writing GFF file.");
				vw.writeGffFile(gffout);
				t3.stop("Time: ");
			}
		}
	}
	
	/**
	 * Writes histogram files for score, zygosity, and quality distributions.
	 * @param scoreHistFile Score histogram output file path
	 * @param zygosityHistFile Zygosity histogram output file path
	 * @param qualityHistFile Quality histogram output file path
	 * @param scoreArray Score distribution array
	 * @param ploidyArray Ploidy distribution array
	 * @param avgQualityArray Average quality distribution array
	 * @param maxQualityArray Maximum quality distribution array
	 */
	public static void writeHistograms(String scoreHistFile, String zygosityHistFile, String qualityHistFile,
			long[][] scoreArray, long[] ploidyArray, long[][] avgQualityArray, long[] maxQualityArray) {
		
		if(scoreHistFile!=null || zygosityHistFile!=null || qualityHistFile!=null){
			Timer t3=new Timer("Writing histograms.");
			if(scoreHistFile!=null){
				writeScoreHist(scoreHistFile, scoreArray[0]);
			}
			if(zygosityHistFile!=null){
				writeZygosityHist(zygosityHistFile, ploidyArray);
			}
			if(qualityHistFile!=null){
				writeQualityHist(qualityHistFile, avgQualityArray[0], maxQualityArray);
			}
			t3.stop("Time: ");
		}
	}
	
	/**
	 * Writes score histogram to specified file.
	 * @param fname Output file name
	 * @param array Score histogram array
	 * @return True if successful
	 */
	static boolean writeScoreHist(String fname, long[] array){
		int max=array.length-1;
		for(; max>=0; max--){
			if(array[max]!=0){break;}
		}
		long sum=0, sum2=0;
		for(int i=0; i<=max; i++){
			sum+=array[i];
			sum2+=(i*array[i]);
		}
		TextStreamWriter tsw=new TextStreamWriter(fname, true, false, false);
		tsw.start();
		tsw.println("#ScoreHist");
		tsw.println("#Vars\t"+sum);
		tsw.println("#Mean\t"+Tools.format("%.2f", sum2*1.0/sum));
		tsw.println("#Median\t"+Tools.medianHistogram(array));
		tsw.println("#Mode\t"+Tools.calcModeHistogram(array));
		tsw.println("#Quality\tCount");
		for(int i=0; i<=max; i++){
			tsw.println(i+"\t"+array[i]);
		}
		tsw.poisonAndWait();
		return tsw.errorState;
	}
	
	/**
	 * Writes zygosity histogram to specified file.
	 * @param fname Output file name
	 * @param array Zygosity histogram array
	 * @return True if successful
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
	 * Writes quality histogram to specified file.
	 * @param fname Output file name
	 * @param avgQualArray Average quality histogram array
	 * @param maxQualArray Maximum quality histogram array
	 * @return True if successful
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