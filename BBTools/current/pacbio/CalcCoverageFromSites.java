package pacbio;

import java.util.BitSet;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import stream.SiteScoreR;
import structures.CoverageArray;
import structures.CoverageArray2;
import var.GenerateVarlets;

/**
 * @author Brian Bushnell
 * @date Jul 19, 2012
 *
 */
public class CalcCoverageFromSites {
	
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
		if(outfile==null){
			process(infile, genome, mincoverage);
		}else{
			processAndWrite(infile, genome, mincoverage, outfile);
		}
		t.stop();
		System.out.println("Time: \t"+t);
	}
	
	
	public static void processAndWrite(final String fname, final int genome, final int mincoverage, final String outpattern){
		Data.setGenome(genome);
		
		BitSet bs=new BitSet();
		
		CoverageArray[] coverage=new CoverageArray[Data.numChroms+1];
		byte[][] correct=new byte[Data.numChroms+1][];
		
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			coverage[chrom]=new CoverageArray2(chrom, Data.chromLengths[chrom]);
			correct[chrom]=new byte[Data.chromLengths[chrom]];
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
				
				if(bs!=null && ssr.numericID<=Integer.MAX_VALUE){
					bs.set((int)ssr.numericID);
				}
				
				int len=ssr.stop-ssr.start+1;
				totalSites++;
				totalSiteLen+=len;
				if(ssr.correct){
					correctSites++;
					correctSiteLen+=len;
				}
				
				
				int chrom=ssr.chrom;
				int min=Tools.max(ssr.start+MIN_END_DIST,  0);
				int max=Tools.min(ssr.stop-MIN_END_DIST,  Data.chromLengths[chrom]-1);
				
				CoverageArray ca=coverage[chrom];
				for(int i=min; i<=max; i++){
					ca.increment(i);
				}
				
				if(ssr.correct){
					byte[] array=correct[chrom];
					for(int i=min; i<=max; i++){
						if(array[i]<Byte.MAX_VALUE){array[i]++;}
					}
				}
			}
			
			s=tf.nextLine();
		}
		tf.close();
		
		for(int i=1; i<coverage.length; i++){
			if(coverage[i].arrayLength()-coverage[i].maxIndex>2000){coverage[i].resize(coverage[i].maxIndex+1);}
			ReadWrite.writeObjectInThread(coverage[i], outpattern.replaceFirst("#", ""+i), false);
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
			CoverageArray cov=coverage[chrom];
			byte[] cor=correct[chrom];
			for(int i=0; i<cov.maxIndex; i++){
				char b=Tools.toUpperCase((char)cha.get(i));
				if(!AminoAcid.isFullyDefined(b)){b='N';}
				
				int total=cov.get(i);
				int good=cor[i];
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
		double invlen=1.0/length;
		double invbase=1.0/baseCount;
		double invn=1.0/nCount;
		double invnc=1.0/nCountCovered; //covered N's
		
		double totalCoverageB=totalCoverage*invlen;
		double totalCoverageBaseB=totalCoverageBase*invbase;
		double totalCoverageNB=totalCoverageN*invn;
		double correctCoverageB=correctCoverage*invlen;
		double correctCoverageBaseB=correctCoverageBase*invbase;
		double correctCoverageNB=correctCoverageN*invn;
		
		double onlyCorrectBaseB=onlyCorrectBase*invbase*100;
		double onlyIncorrectBaseB=onlyIncorrectBase*invbase*100;
		double onlyCorrectNB=onlyCorrectN*invnc*100;
		double onlyIncorrectNB=onlyIncorrectN*invnc*100;
		double mostlyCorrectBaseB=mostlyCorrectBase*invbase*100;
		double mostlyIncorrectBaseB=mostlyIncorrectBase*invbase*100;
		double mostlyCorrectNB=mostlyCorrectN*invnc*100;
		double mostlyIncorrectNB=mostlyIncorrectN*invnc*100;
		double anyCorrectBaseB=anyCorrectBase*invbase*100;
		double anyCorrectNB=anyCorrectN*invnc*100;
		double noCorrectBaseB=noCorrectBase*invbase*100;
		double noCoverageBaseB=noCoverageBase*invbase*100;
		double noCoverageNB=noCoverageN*invn*100;
		


		double correctSitesB=correctSites*100d/totalSites;
		double correctSiteLenB=correctSiteLen*100d/totalSiteLen;
		
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
	
	
	public static void process(final String fname, final int genome, final int mincoverage){
		Data.setGenome(genome);
		
		BitSet bs=new BitSet();
		
		byte[][] coverage=new byte[Data.numChroms+1][];
		byte[][] correct=new byte[Data.numChroms+1][];
		
		for(int chrom=1; chrom<=Data.numChroms; chrom++){
			coverage[chrom]=new byte[Data.chromLengths[chrom]];
			correct[chrom]=new byte[Data.chromLengths[chrom]];
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
				
				if(bs!=null){
					bs.set((int)ssr.numericID);
				}
				
				int len=ssr.stop-ssr.start+1;
				totalSites++;
				totalSiteLen+=len;
				if(ssr.correct){
					correctSites++;
					correctSiteLen+=len;
				}
				
				
				int chrom=ssr.chrom;
				int min=Tools.max(ssr.start,  0);
				int max=Tools.min(ssr.stop,  Data.chromLengths[chrom]-1);
				byte[] array=coverage[chrom];
				for(int i=min; i<=max; i++){
					if(array[i]<Byte.MAX_VALUE){array[i]++;}
				}
				if(ssr.correct){
					array=correct[chrom];
					for(int i=min; i<=max; i++){
						if(array[i]<Byte.MAX_VALUE){array[i]++;}
					}
				}
			}
			
			s=tf.nextLine();
		}
		tf.close();
		
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
			byte[] cov=coverage[chrom];
			byte[] cor=correct[chrom];
			for(int i=0; i<cov.length; i++){
				char b=Tools.toUpperCase((char)cha.get(i));
				if(!AminoAcid.isFullyDefined(b)){b='N';}
				
				int total=cov[i];
				int good=cor[i];
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
		double invlen=1.0/length;
		double invbase=1.0/baseCount;
		double invn=1.0/nCount;
		double invnc=1.0/nCountCovered; //covered N's
		
		double totalCoverageB=totalCoverage*invlen;
		double totalCoverageBaseB=totalCoverageBase*invbase;
		double totalCoverageNB=totalCoverageN*invn;
		double correctCoverageB=correctCoverage*invlen;
		double correctCoverageBaseB=correctCoverageBase*invbase;
		double correctCoverageNB=correctCoverageN*invn;
		
		double onlyCorrectBaseB=onlyCorrectBase*invbase*100;
		double onlyIncorrectBaseB=onlyIncorrectBase*invbase*100;
		double onlyCorrectNB=onlyCorrectN*invnc*100;
		double onlyIncorrectNB=onlyIncorrectN*invnc*100;
		double mostlyCorrectBaseB=mostlyCorrectBase*invbase*100;
		double mostlyIncorrectBaseB=mostlyIncorrectBase*invbase*100;
		double mostlyCorrectNB=mostlyCorrectN*invnc*100;
		double mostlyIncorrectNB=mostlyIncorrectN*invnc*100;
		double anyCorrectBaseB=anyCorrectBase*invbase*100;
		double anyCorrectNB=anyCorrectN*invnc*100;
		double noCorrectBaseB=noCorrectBase*invbase*100;
		double noCoverageBaseB=noCoverageBase*invbase*100;
		double noCoverageNB=noCoverageN*invn*100;
		


		double correctSitesB=correctSites*100d/totalSites;
		double correctSiteLenB=correctSiteLen*100d/totalSiteLen;
		
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
	
	
	
	public static SiteScoreR[] toSites(String s){
		String[] split=s.split("\t");
		SiteScoreR[] scores=new SiteScoreR[split.length];
		for(int i=0; i<split.length; i++){
			scores[i]=SiteScoreR.fromText(split[i]);
		}
		return scores;
	}
	
	public static int MIN_END_DIST=GenerateVarlets.MIN_END_DIST; //These must be the same.
	
}
