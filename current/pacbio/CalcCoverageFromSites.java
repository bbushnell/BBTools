package pacbio;

import java.util.BitSet;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.PreParser;
import shared.Timer;
import shared.Tools;
import stream.SiteScoreR;
import structures.CoverageArray3;
import var.GenerateVarlets;

/**
 * Calculates coverage and correctness summaries from tab-delimited SiteScoreR
 * groups.  Each site contributes after {@link #MIN_END_DIST} bases are
 * trimmed from both ends, reports scan complete loaded chromosomes, and total
 * plus correct coverage are accumulated in 32-bit arrays so output and
 * no-output modes report the same depths.
 * @author Brian Bushnell
 * @date Jul 19, 2012
 */
public class CalcCoverageFromSites{

	/**
	 * Parses the input table, optional output pattern, genome build and
	 * mincoverage threshold, then runs coverage accumulation.
	 * @param args Input path, coverage output path or null, genome build, and
	 * optional mincoverage flag
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
		if(outfile.equalsIgnoreCase("null")){outfile=null;}
		assert(outfile==null || outfile.contains("#"));
		int genome=Integer.parseInt(args[2]);
		int mincoverage=1;
		for(int i=3; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split[1];

			if(a.equals("mincoverage")){
				mincoverage=Integer.parseInt(b);
			}
		}
		processAndWrite(infile, genome, mincoverage, outfile);
		t.stop();
		System.out.println("Time: \t"+t);
	}

	/**
	 * Computes coverage, optionally writes per-chromosome coverage objects, and
	 * prints aggregate coverage/correctness statistics. Each site is trimmed by
	 * {@link #MIN_END_DIST} bases on both ends before coverage is accumulated;
	 * site and read counters still include sites that become empty after trimming.
	 * Total and correct coverage use 8 bytes per reference position for two
	 * 32-bit arrays and throw before overflowing.
	 * Legacy read IDs above Integer.MAX_VALUE are excluded from Reads Represented.
	 * @param fname Input tab-delimited SiteScoreR groups
	 * @param genome Genome build to load
	 * @param mincoverage Minimum depth for a base to count as covered
	 * @param outpattern Output coverage pattern containing {@code #}, or null
	 */
	public static void processAndWrite(final String fname, final int genome,
			final int mincoverage, final String outpattern){
		Data.setGenome(genome);

		BitSet bs=new BitSet();

		CoverageArray3[] coverage=new CoverageArray3[Data.numChroms+1];
		CoverageArray3[] correct=new CoverageArray3[Data.numChroms+1];

		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			coverage[chrom]=new CoverageArray3(chrom, Data.chromLengths[chrom]);
			correct[chrom]=new CoverageArray3(chrom, Data.chromLengths[chrom]);
		}

		TextFile tf=new TextFile(fname, true);
		String s=tf.nextLine();

		long totalSites=0;
		long correctSites=0;
		long totalSiteLen=0;
		long correctSiteLen=0;

		while(s!=null){
			SiteScoreR[] sites=toSites(s);
			for(SiteScoreR ssr : sites){

				if(ssr.numericID<=Integer.MAX_VALUE){
					bs.set((int)ssr.numericID);
				}

				long len=(long)ssr.stop-ssr.start+1;
				totalSites++;
				totalSiteLen+=len;
				if(ssr.correct){
					correctSites++;
					correctSiteLen+=len;
				}

				int chrom=ssr.chrom;
				long minLong=Math.max((long)ssr.start+MIN_END_DIST, 0);
				long maxLong=Math.min((long)ssr.stop-MIN_END_DIST,
						(long)Data.chromLengths[chrom]-1);

				if(minLong<=maxLong){
					incrementCoverage(coverage[chrom], correct[chrom], ssr.correct,
							(int)minLong, (int)maxLong, chrom);
				}
			}

			s=tf.nextLine();
		}
		tf.close();

		if(outpattern!=null){
			for(int i=1; i<coverage.length; i++){
				CoverageArray3 cov=coverage[i];
				if(cov.maxIndex>=0 && cov.arrayLength()-cov.maxIndex>2000){
					cov.resize(cov.maxIndex+1);
				}
				ReadWrite.writeObjectInThread(cov, outpattern.replaceFirst("#", ""+i), false);
			}
		}

		long totalCoverage=0;
		long totalCoverageBase=0;
		long totalCoverageN=0;
		long correctCoverage=0;
		long correctCoverageBase=0;
		long correctCoverageN=0;

		long onlyCorrectBase=0;
		long onlyIncorrectBase=0;
		long onlyCorrectN=0;
		long onlyIncorrectN=0;
		long mostlyCorrectBase=0;
		long mostlyIncorrectBase=0;
		long mostlyCorrectN=0;
		long mostlyIncorrectN=0;
		long anyCorrectBase=0;
		long anyCorrectN=0;
		long noCorrectBase=0;
		long noCoverageBase=0;
		long noCoverageN=0;

		long baseCount=0;
		long nCount=0;
		long nCountCovered=0;

		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			ChromosomeArray cha=Data.getChromosome(chrom);
			CoverageArray3 cov=coverage[chrom];
			CoverageArray3 cor=correct[chrom];
			for(int i=0; i<Data.chromLengths[chrom]; i++){
				char b=Tools.toUpperCase((char)cha.get(i));
				if(!AminoAcid.isFullyDefined(b)){b='N';}

				int total=cov.get(i);
				int good=cor.get(i);
				int bad=total-good;

				totalCoverage+=total;
				correctCoverage+=good;

				if(b=='N'){
					nCount++;
					if(total>=mincoverage){
						totalCoverageN+=total;
						correctCoverageN+=good;
						nCountCovered++;
						if(total==good){
							onlyCorrectN++;
							mostlyCorrectN++;
						}else if(good>bad){
							mostlyCorrectN++;
						}else if(good==0){
							onlyIncorrectN++;
							mostlyIncorrectN++;
						}else if(bad>good){
							mostlyIncorrectN++;
						}
						if(good>0){anyCorrectN++;}
					}else{
						noCoverageN++;
					}
				}else{
					baseCount++;
					if(total>=mincoverage){
						totalCoverageBase+=total;
						correctCoverageBase+=good;
						if(total==good){
							onlyCorrectBase++;
							mostlyCorrectBase++;
						}else if(good>bad){
							mostlyCorrectBase++;
						}else if(good==0){
							onlyIncorrectBase++;
							mostlyIncorrectBase++;
							noCorrectBase++;
						}else if(bad>good){
							mostlyIncorrectBase++;
						}
						if(good>0){anyCorrectBase++;}
					}else{
						noCoverageBase++;
						noCorrectBase++;
					}
				}
			}
			Data.unload(chrom, true);
			coverage[chrom]=null;
			correct[chrom]=null;
		}

		long length=nCount+baseCount;

		double totalCoverageB=ratio(totalCoverage, length);
		double totalCoverageBaseB=ratio(totalCoverageBase, baseCount);
		double totalCoverageNB=ratio(totalCoverageN, nCount);
		double correctCoverageB=ratio(correctCoverage, length);
		double correctCoverageBaseB=ratio(correctCoverageBase, baseCount);
		double correctCoverageNB=ratio(correctCoverageN, nCount);

		double onlyCorrectBaseB=percent(onlyCorrectBase, baseCount);
		double onlyIncorrectBaseB=percent(onlyIncorrectBase, baseCount);
		double onlyCorrectNB=percent(onlyCorrectN, nCountCovered);
		double onlyIncorrectNB=percent(onlyIncorrectN, nCountCovered);
		double mostlyCorrectBaseB=percent(mostlyCorrectBase, baseCount);
		double mostlyIncorrectBaseB=percent(mostlyIncorrectBase, baseCount);
		double mostlyCorrectNB=percent(mostlyCorrectN, nCountCovered);
		double mostlyIncorrectNB=percent(mostlyIncorrectN, nCountCovered);
		double anyCorrectBaseB=percent(anyCorrectBase, baseCount);
		double anyCorrectNB=percent(anyCorrectN, nCountCovered);
		double noCorrectBaseB=percent(noCorrectBase, baseCount);
		double noCoverageBaseB=percent(noCoverageBase, baseCount);
		double noCoverageNB=percent(noCoverageN, nCount);

		double correctSitesB=percent(correctSites, totalSites);
		double correctSiteLenB=percent(correctSiteLen, totalSiteLen);

		System.out.println("\nOverall Statistics");

		if(bs!=null){
			System.out.println("Reads Represented:       \t"+bs.cardinality());
		}
		System.out.println(Tools.format("Total Correct Sites:     \t"+(correctSitesB<10?" ":"")+"%.3f%%  ", correctSitesB)+" \t"+correctSites);
		System.out.println(Tools.format("Total Correct Site Length:\t"+(correctSiteLenB<10?" ":"")+"%.3f%%  ", correctSiteLenB)+" \t"+correctSiteLen);

		System.out.println("\nCoverage Statistics");

		System.out.println(Tools.format("Avg Coverage:            \t"+(totalCoverageB<10?" ":"")+"%.3f", totalCoverageB)+"  \t"+totalCoverage);

		System.out.println(Tools.format("Avg Coverage Base:       \t"+(totalCoverageBaseB<10?" ":"")+"%.3f", totalCoverageBaseB)+"  \t"+totalCoverageBase);

		System.out.println(Tools.format("Avg Coverage N:          \t"+(totalCoverageNB<10?" ":"")+"%.3f", totalCoverageNB)+"  \t"+totalCoverageN);

		System.out.println(Tools.format("Correct Coverage:        \t"+(correctCoverageB<10?" ":"")+"%.3f", correctCoverageB)+"  \t"+correctCoverage);

		System.out.println(Tools.format("Correct Coverage Base:   \t"+(correctCoverageBaseB<10?" ":"")+"%.3f", correctCoverageBaseB)+"  \t"+correctCoverageBase);

		System.out.println(Tools.format("Correct Coverage N:      \t"+(correctCoverageNB<10?" ":"")+"%.3f", correctCoverageNB)+"  \t"+correctCoverageN);

		System.out.println("\nStatistics over Defined Bases");

		System.out.println(Tools.format("onlyCorrect:             \t"+(onlyCorrectBaseB<10?" ":"")+"%.3f", onlyCorrectBaseB)+"%");
		System.out.println(Tools.format("mostlyCorrect:           \t"+(mostlyCorrectBaseB<10?" ":"")+"%.3f", mostlyCorrectBaseB)+"%");
		System.out.println(Tools.format("anyCorrect:              \t"+(anyCorrectBaseB<10?" ":"")+"%.3f", anyCorrectBaseB)+"%");
		System.out.println(Tools.format("noCorrect:               \t"+(noCorrectBaseB<10?" ":"")+"%.3f", noCorrectBaseB)+"%");
		System.out.println(Tools.format("mostlyIncorrect:         \t"+(mostlyIncorrectBaseB<10?" ":"")+"%.3f", mostlyIncorrectBaseB)+"%");
		System.out.println(Tools.format("onlyIncorrect:           \t"+(onlyIncorrectBaseB<10?" ":"")+"%.3f", onlyIncorrectBaseB)+"%");
		System.out.println(Tools.format("noCoverage:              \t"+(noCoverageBaseB<10?" ":"")+"%.3f", noCoverageBaseB)+"%");

		System.out.println("\nStatistics over N (for covered locations)");

		System.out.println(Tools.format("onlyCorrect:             \t"+(onlyCorrectNB<10?" ":"")+"%.3f", onlyCorrectNB)+"%");
		System.out.println(Tools.format("mostlyCorrect:           \t"+(mostlyCorrectNB<10?" ":"")+"%.3f", mostlyCorrectNB)+"%");
		System.out.println(Tools.format("anyCorrect:              \t"+(anyCorrectNB<10?" ":"")+"%.3f", anyCorrectNB)+"%");
		System.out.println(Tools.format("mostlyIncorrect:         \t"+(mostlyIncorrectNB<10?" ":"")+"%.3f", mostlyIncorrectNB)+"%");
		System.out.println(Tools.format("onlyIncorrect:           \t"+(onlyIncorrectNB<10?" ":"")+"%.3f", onlyIncorrectNB)+"%");
		System.out.println(Tools.format("noCoverage (over all N): \t"+(noCoverageNB<10?" ":"")+"%.3f", noCoverageNB)+"%");


	}

	/**
	 * Computes coverage without writing coverage objects.
	 * @param fname Input tab-delimited SiteScoreR groups
	 * @param genome Genome build to load
	 * @param mincoverage Minimum depth for a base to count as covered
	 */
	public static void process(final String fname, final int genome, final int mincoverage){
		processAndWrite(fname, genome, mincoverage, null);
	}

	/**
	 * Increments total and optional correct coverage over one inclusive range.
	 * @param total Total coverage array to update
	 * @param correct Correct coverage array to update when isCorrect is true
	 * @param isCorrect True if this range should increment correct coverage
	 * @param min First position
	 * @param max Last position
	 * @param chrom Chromosome identifier for diagnostics
	 */
	private static void incrementCoverage(CoverageArray3 total, CoverageArray3 correct,
			boolean isCorrect, int min, int max, int chrom){
		int[] totalArray=total.array;
		int[] correctArray=(isCorrect ? correct.array : null);
		for(int i=min; i<=max; i++){
			if(totalArray[i]==Integer.MAX_VALUE){
				throw new ArithmeticException("Coverage depth overflow on chrom "
						+chrom+" at position "+i);
			}
			totalArray[i]++;
			if(isCorrect){correctArray[i]++;}
		}
		total.minIndex=Tools.min(min, total.minIndex);
		total.maxIndex=Tools.max(max, total.maxIndex);
		if(isCorrect){
			correct.minIndex=Tools.min(min, correct.minIndex);
			correct.maxIndex=Tools.max(max, correct.maxIndex);
		}
	}

	/**
	 * Divides two integers with an explicit zero-denominator guard.
	 */
	private static double ratio(long numerator, long denominator){
		return denominator<1 ? 0 : numerator/(double)denominator;
	}

	/**
	 * Converts a numerator and denominator to a guarded percentage.
	 */
	private static double percent(long numerator, long denominator){
		return denominator<1 ? 0 : numerator*100d/denominator;
	}

	/**
	 * Parses one tab-delimited line of serialized SiteScoreR records.
	 * @param s Tab-delimited serialized SiteScoreR records
	 * @return Parsed sites in input order
	 */
	public static SiteScoreR[] toSites(String s){
		String[] split=s.split("\t");
		SiteScoreR[] scores=new SiteScoreR[split.length];
		for(int i=0; i<split.length; i++){
			scores[i]=SiteScoreR.fromText(split[i]);
		}
		return scores;
	}

	/** End distance shared with GenerateVarlets site filtering */
	public static int MIN_END_DIST=GenerateVarlets.MIN_END_DIST;

}
