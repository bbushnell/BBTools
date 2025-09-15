package jgi;

import java.io.File;
import java.io.PrintStream;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import bloom.KCountArray;
import bloom.KmerCountAbstract;
import bloom.ReadCounter;
import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;



/**
 * This class is designed to visualize the distribution of kmer depths across individual reads.
 * @author Brian Bushnell
 * @date May 15, 2013
 *
 */
public class ReadKmerDepthDistribution {

	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		String reads1=(args[0].indexOf("=")>0 ? null : args[0]);
		String reads2=(reads1!=null && args.length>1 ? args[1] : null);
		if(reads2!=null && "null".equalsIgnoreCase(reads2)){reads2=null;}
		
		{
			if(reads1!=null && !reads1.contains(",")){
				File f=new File(reads1);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(reads1+" does not exist.");}
			}
			if(reads2!=null && !reads2.contains(",")){
				File f=new File(reads2);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(reads2+" does not exist.");}
				if(reads2.equalsIgnoreCase(reads1)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		KmerCountAbstract.minQuality=4;
		KmerCountAbstract.minProb=0.4f;
		KmerCountAbstract.CANONICAL=true;
		
		int k=31;
		int cbits=32;
		int hashes=3;
//		int matrixbits=-1;
		long cells=-1;
		long maxReads=-1;
		int buildpasses=1;
		long tablereads=-1; //How many reads to process when building the hashtable
		int buildStepsize=4;
		String outKeep=null;
		int prehashes=-1;
		long precells=-1;
		String histFile=null;
		int threads=-1;
		ReadWrite.ZIPLEVEL=2;
		int minq=KmerCountAbstract.minQuality;
		
		boolean auto=true;
		boolean deterministic=true;
		
		List<String> extra=null;
		
		long memory=Runtime.getRuntime().maxMemory();
		
		Parser parser=new Parser();
		for(int i=(reads1==null ? 0 : 1); i<args.length; i++){
			if(args[i]==null){args[i]="null";}
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.equals("in") || a.equals("in1")){
				reads1=b;
			}else if(a.equals("in2")){
				reads2=b;
			}else if(a.startsWith("bits") ||a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("histlen") ||a.startsWith("histogramlen")){
				HIST_LEN_PRINT=Tools.min(Integer.MAX_VALUE, Long.parseLong(b)+1);
			}else if(a.startsWith("matrixbits")){
				int matrixbits=Integer.parseInt(b);
				assert(matrixbits<63);
				cells=1L<<matrixbits;
			}else if(a.startsWith("cells")){
				cells=Parse.parseKMG(b);
			}else if(a.startsWith("precells") || a.startsWith("prefiltercells")){
				precells=Parse.parseKMG(b);
				prefilter=prefilter || precells!=0;
			}else if(a.startsWith("minq")){
				minq=Byte.parseByte(b);
			}else if(a.equals("zerobin")){
				ZERO_BIN=Parse.parseBoolean(b);
			}else if(a.equals("deterministic") || a.equals("dr")){
				boolean x=Parse.parseBoolean(b);
				deterministic=x;
			}else if(a.startsWith("minprob")){
				KmerCountAbstract.minProb=Float.parseFloat(b);
			}else if(a.startsWith("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.startsWith("prehashes") || a.startsWith("prefilterhashes")){
				prehashes=Integer.parseInt(b);
				prefilter=prefilter || prehashes!=0;
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}else if(a.startsWith("stepsize") || a.startsWith("buildstepsize")){
				buildStepsize=Integer.parseInt(b);
			}else if(a.startsWith("passes") || a.startsWith("buildpasses")){
				buildpasses=Integer.parseInt(b);
			}else if(a.equals("printcoverage")){
				assert(false) : "This is not the program you are looking for.  Try KmerCoverage.";
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("tablereads") || a.startsWith("buildreads")){
				tablereads=Parse.parseKMG(b);
			}else if(a.equals("out") || a.equals("outk") || a.equals("outkeep") || a.equals("outgood")){
				outKeep=b;
			}else if(a.startsWith("hist")){
				histFile=b;
			}else if(a.startsWith("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("auto") || a.equals("automatic")){
				auto=Parse.parseBoolean(b);
			}else if(a.equals("canonical")){
				CANONICAL=KmerCountAbstract.CANONICAL=Parse.parseBoolean(b);
			}else if(a.equals("fixspikes")){
				FIX_SPIKES=Parse.parseBoolean(b);
			}else if(a.equals("printzerocoverage") || a.equals("pzc")){
				PRINT_ZERO_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("removeduplicatekmers") || a.equals("rdk")){
				KmerCountAbstract.KEEP_DUPLICATE_KMERS=!Parse.parseBoolean(b);
			}else if(a.equals("target") || a.equals("targetdepth")){
				TARGET_DEPTH=Integer.parseInt(b);
			}else if(a.equals("max") || a.equals("maxdepth")){
				MAX_DEPTH=Integer.parseInt(b);
			}else if(a.equals("min") || a.equals("mindepth")){
				MIN_DEPTH=Integer.parseInt(b);
			}else if(a.equals("minkmers") || a.equals("minkmersovermindepth") || a.equals("mingoodkmersperread") || a.equals("mgkpr")){
				MIN_KMERS_OVER_MIN_DEPTH=Tools.max(1, Integer.parseInt(b));
			}else if(a.equals("percentile") || a.equals("depthpercentile") || a.equals("dp")){
				DEPTH_PERCENTILE=Float.parseFloat(b);
				if(DEPTH_PERCENTILE>1 && DEPTH_PERCENTILE<=100){DEPTH_PERCENTILE/=100;}
				assert(DEPTH_PERCENTILE>=0 && DEPTH_PERCENTILE<=1) : "Depth percentile must be between 0 and 100.";
			}else if(a.equals("extra")){
				if(b!=null && !b.equalsIgnoreCase("null")){
					if(new File(b).exists()){
						extra=new ArrayList<String>();
						extra.add(b);
					}else{
						extra=Arrays.asList(b.split(","));
					}
				}
			}else{
				throw new RuntimeException("Unknown parameter "+arg);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		MAX_DEPTH=Tools.max(MAX_DEPTH, TARGET_DEPTH);
		assert(TARGET_DEPTH>0);
		
		assert(FastaReadInputStream.settingsOK());
		if(k>31){CANONICAL=KmerCountAbstract.CANONICAL=false;}
		assert(CANONICAL==KmerCountAbstract.CANONICAL);
		
//		if(output!=null && reads1.contains(",")){
//			throw new RuntimeException("\nLists of input files can only be used with histogram output, not full output.\n" +
//					"Please set output=null or move the extra input files to 'extra=file1,file2,...fileN'");
//		}
		
		{
			if(histFile==null){
//				HIST_LEN=Tools.min(20000, HIST_LEN);
//				HIST_LEN_PRINT=Tools.min(20000, HIST_LEN_PRINT);
			}else{
				USE_HISTOGRAM=true;
			}
			
			final int maxCount=(int)(cbits>16 ? Integer.MAX_VALUE : (1L<<cbits)-1);
			assert(maxCount>0);
			HIST_LEN_PRINT=Tools.max(1, Tools.min(HIST_LEN_PRINT, maxCount+1));
			assert(HIST_LEN_PRINT<=Integer.MAX_VALUE) : HIST_LEN_PRINT+", "+Integer.MAX_VALUE;
			HIST_LEN=(int)Tools.min(maxCount+1, Tools.max(HIST_LEN_PRINT, HIST_LEN));
			THREAD_HIST_LEN=Tools.min(THREAD_HIST_LEN, HIST_LEN);

			histogram_total=new AtomicLongArray(HIST_LEN);
		}
		
		if(extra!=null){
			for(String s : extra){
				File f=new File(s);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(s+" does not exist.");}
				assert(!s.equalsIgnoreCase(reads1) && (reads2==null || !s.equalsIgnoreCase(reads2))) : "\nInput file "+s+" should not be included as an extra file.\n";
			}
		}
		
//		outstream.println("ForceInterleaved = "+FASTQ.FORCE_INTERLEAVED);
		
//		assert(false) : reads1+", "+reads2+", "+output;
//		if(FASTQ.FORCE_INTERLEAVED && in2==null){
//			outstream.println()
//		}
		
		if(threads<=0){
			if(auto){THREADS=Shared.LOGICAL_PROCESSORS;}
			else{THREADS=8;}
		}else{
			THREADS=threads;
		}
		
//		System.err.println("THREADS="+THREADS+", KmerCountAbstract.THREADS="+KmerCountAbstract.THREADS);
		
		if(auto && cells==-1){
			final long usable=(long)Tools.max(((memory-96000000)*.73), memory*0.45);
			long mem=usable-(USE_HISTOGRAM ? (HIST_LEN*8*(1)) : 0);
			if(buildpasses>1){mem/=2;}
			cells=(mem*8)/cbits;
//
//			long tablebytes=((1L<<matrixbits)*cbits)/8;
//			if(tablebytes*3<usable){matrixbits++;}
//			outstream.println(tablebytes/1000000+", "+usable/1000000+", "+(tablebytes*3)/1000000);
			
		}else if(cells==-1){
			cells=1L<<34;
		}
		
		if(prefilter){
			if(precells<1){
				long totalbits=cells*cbits;
				long prebits=(long)(totalbits*0.35);
				precells=prebits/2;
				cells=(totalbits-prebits+cbits-1)/cbits; //Steal memory from cell allocation
			}
			if(prehashes<1){
				prehashes=(hashes+1)/2;
			}
		}
		
		{
			outstream.println("\nSettings:");
			outstream.println("threads:          \t"+THREADS);
			outstream.println("k:                \t"+k);
			outstream.println("deterministic:    \t"+deterministic);
			outstream.println("passes:           \t"+buildpasses);
			outstream.println("bits per cell:    \t"+cbits);
//			outstream.println("matrixbits: \t"+matrixbits);
			outstream.println("cells:            \t"+Tools.toKMG(cells));
			outstream.println("hashes:           \t"+hashes);
			if(prefilter){
				outstream.println("prefilter bits:   \t"+2);
//				outstream.println("matrixbits: \t"+matrixbits);
				outstream.println("prefilter cells:  \t"+(precells>0 && prehashes>0 ? Tools.toKMG(precells) : "?"));
				outstream.println("prefilter hashes: \t"+(precells>0 && prehashes>0 ? ""+prehashes : "?"));
			}
			outstream.println("base min quality: \t"+KmerCountAbstract.minQuality);
			outstream.println("kmer min prob:    \t"+KmerCountAbstract.minProb);
			
			outstream.println();
			outstream.println("target depth:     \t"+TARGET_DEPTH);
			outstream.println("min depth:        \t"+MIN_DEPTH);
			outstream.println("max depth:        \t"+MAX_DEPTH);
			outstream.println("min good kmers:   \t"+MIN_KMERS_OVER_MIN_DEPTH);
			outstream.println("depth percentile: \t"+Tools.format("%.1f", 100*DEPTH_PERCENTILE));
			outstream.println("remove duplicates:\t"+!KmerCountAbstract.KEEP_DUPLICATE_KMERS);
			outstream.println("fix spikes:       \t"+FIX_SPIKES);
			if(USE_HISTOGRAM && HIST_LEN>0){
				outstream.println("histogram length: \t"+(USE_HISTOGRAM ? HIST_LEN : 0));
			}
			if(histFile!=null){
				outstream.println("print zero cov:   \t"+PRINT_ZERO_COVERAGE);
			}
			
			outstream.println();
		}
		
		if(!prefilter && k<32 && cells>(1L<<(2*k))){cells=(1L<<(2*k));}
		assert(cells>0);
		
//		KmerCountAbstract.THREADS=Tools.max(THREADS/2, KmerCountAbstract.THREADS);  //Seems like 4 is actually optimal...
		
		FastaReadInputStream.MIN_READ_LEN=k;
		
		Timer t=new Timer();
		Timer ht=new Timer();
		t.start();
		ht.start();
		KCountArray kca;
		KCountArray prefilterArray=null;
//		outstream.println();
		ReadCounter rc=new ReadCounter(k, true, false, false, Shared.AMINO_IN);
		rc.detectStepsize=buildStepsize;
		if(prefilter){
			prefilterArray=rc.makeKca(reads1, reads2, extra, 2, precells, prehashes, minq,
					tablereads, 1, 1, 1, null, 0);
			outstream.println("Made prefilter:   \t"+prefilterArray.toShortString(prehashes));
			double uf=prefilterArray.usedFraction();
			if(uf>0.6){
				outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
					uf>0.7 ? "fairly" : "somewhat")+" full, which may reduce accuracy for kmers of depth under 3.  Ideal load is under 60% used." +
						"\nFor better accuracy, run on a node with more memory; quality-trim or error-correct reads; " +
						"or increase the values of the minprob flag to reduce spurious kmers.");
			}
		}
		kca=rc.makeKca(reads1, reads2, extra, cbits, cells, hashes, minq, 
				tablereads, buildpasses, 2, 2, prefilterArray, (prefilterArray==null ? 0 : prefilterArray.maxValue));
		ht.stop();
		
		outstream.println("Made hash table:  \t"+kca.toShortString(hashes));
		double uf=kca.usedFraction();
		if(uf>0.6){
			outstream.println("Warning:  This table is "+(uf>0.995 ? "totally" : uf>0.99 ? "crazy" : uf>0.95 ? "incredibly" : uf>0.9 ? "extremely" : uf>0.8 ? "very" :
				uf>0.7 ? "fairly" : "somewhat")+" full, which may reduce accuracy.  Ideal load is under 60% used." +
				"\nFor better accuracy, use the 'prefilter' flag; run on a node with more memory; quality-trim or error-correct reads; " +
					"or increase the values of the minprob flag to reduce spurious kmers.  In practice you should still get good normalization results " +
					"even with loads over 90%, but the histogram and statistics will be off.");
		}
		
		long estUnique;
		outstream.println();
		if(prefilterArray!=null){
			int lim1=prefilterArray.maxValue, lim2=prefilterArray.maxValue+1;
			double a=prefilterArray.estimateUniqueKmers(prehashes);
			double b=kca.estimateUniqueKmers(hashes, lim2);
			a=a-b;
			if(CANONICAL){
//				a=(a*KCountArray.canonMask)/(KCountArray.canonMask+1);
//				b=(b*KCountArray.canonMask)/(KCountArray.canonMask+1);
			}else{
				a/=2;
				b/=2;
			}
			estUnique=((long)((a+b)));
			outstream.println("Estimated kmers of depth 1-"+lim1+": \t"+(long)a);
			outstream.println("Estimated kmers of depth "+lim2+"+ : \t"+(long)b);
		}else{
//			double est=kca.cells*(1-Math.pow(1-Math.sqrt(kca.usedFraction()), 1.0/hashes));
//			double est=kca.cells*(1-Math.pow(1-kca.usedFraction(), 1.0/hashes));
			double est=kca.estimateUniqueKmers(hashes);
//			outstream.println("Used cells: "+kca.cellsUsed(1));
			if(CANONICAL){
//				est=(est*KCountArray.canonMask)/(KCountArray.canonMask+1);
			}else{
				est/=2;
			}
			estUnique=((long)((est)));
			
		}
		outstream.println("Estimated unique kmers:     \t"+estUnique);//+", or "+estUnique+" counting forward kmers only.");
//		outstream.println("(Includes forward and reverse kmers)");
		outstream.println();
		outstream.println("Table creation time:\t\t"+ht);//+"   \t"+Tools.format("%.2f", totalBases*1000000.0/(ht.elapsed))+" kb/sec");
		
		long bases=0;
		
		ListNum.setDeterministicRandom(deterministic);
		
		if(reads1!=null && reads1.contains(",") && !new File(reads1).exists()){
			throw new RuntimeException("This class is not designed to deal with lists of input files.");
		}else{
			bases=count(reads1, reads2, kca, k, maxReads, outKeep, overwrite, histFile, estUnique);
		}
		printTopology();
		
		t.stop();
		outstream.println("\nTotal time:      \t\t"+t+"   \t"+Tools.format("%.2f", bases*1000000.0/(t.elapsed))+" kb/sec");
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	
	public static void printTopology(){
		long total=peaks.get()+spikes.get()+flats.get()+valleys.get()+slopes.get();
		double mult=100.0/total;
		
		long sp=spikes.get();
		long pe=peaks.get();
		long va=valleys.get();
		long sl=slopes.get();
		long fl=flats.get();
		double dsp=mult*sp;
		double dpe=mult*pe;
		double dva=mult*va;
		double dsl=mult*sl;
		double dfl=mult*fl;
		
		System.err.println("\nDepth Topology:\t");
		System.err.println("Spikes:     \t\t\t"+(dsp<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dsp,sp));
		System.err.println("Peaks:      \t\t\t"+(dpe<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dpe,pe));
		System.err.println("Valleys:    \t\t\t"+(dva<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dva,va));
		System.err.println("Slopes:     \t\t\t"+(dsl<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dsl,sl));
		System.err.println("Flats:      \t\t\t"+(dfl<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dfl,fl));
	}


	public static long count(String in1, String in2, KCountArray kca, int k, long maxReads,
			String outKeep, boolean overwrite, String histFile, long estUnique) {
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		ConcurrentReadOutputStream rosKeep=null;
		if(outKeep!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			String out1=outKeep.replaceFirst("#", "1");
			String out2=null;
			
			if(cris.paired()){
				if(outKeep.contains("#")){
					out2=outKeep.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1));
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2)));
			
//			assert(false) : out1+", "+out2;
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, "attachment", true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, "attachment", true, overwrite, append, ordered);
			rosKeep=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
		}
		
		if(rosKeep!=null){
			rosKeep.start();
			outstream.println("Started output threads.");
		}
		
		long bases=downsample(cris, kca, k, maxReads, rosKeep, histFile, overwrite, estUnique);
		
		ReadWrite.closeStreams(cris, rosKeep);
		if(verbose){System.err.println("Closed streams");}
		
		return bases;
	}
	
	
	
	public static long downsample(ConcurrentReadInputStream cris, KCountArray kca, int k, long maxReads, ConcurrentReadOutputStream rosKeep,
			String histFile, boolean overwrite, long estUnique) {
		Timer tdetect=new Timer();
		tdetect.start();

		long totalBases=0;
		long totalReads=0;
		long basesKept=0;
		long readsKept=0;
		long basesTossed=0;
		long readsTossed=0;
		
//		assert(false) : THREADS;
		ProcessThread[] pta=new ProcessThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new ProcessThread(cris, kca, k, rosKeep);
			pta[i].start();
		}
		
		for(int i=0; i<pta.length; i++){
			ProcessThread ct=pta[i];
			synchronized(ct){
				while(ct.getState()!=State.TERMINATED){
					try {
						ct.join(1000);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				totalBases+=ct.totalBases;
				totalReads+=ct.totalReads;
				basesKept+=ct.basesKept;
				readsKept+=ct.readsKept;
				basesTossed+=ct.basesTossed;
				readsTossed+=ct.readsTossed;
				
				for(int j=0; j<ct.hist.length; j++){
					histogram_total.addAndGet(j, ct.hist[j]);
				}
			}
		}
		
		if(!ZERO_BIN && histogram_total!=null && histogram_total.length()>1){
			histogram_total.addAndGet(1, histogram_total.get(0));
			histogram_total.set(0, 0);
		}
		
//		outstream.println();
		tdetect.stop();
		outstream.println("Table read time: \t\t"+tdetect+"   \t"+Tools.format("%.2f", totalBases*1000000.0/(tdetect.elapsed))+" kb/sec");
		
		{
			String pad="";
			String s=""+totalReads;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Total reads in:  \t\t"+totalReads+pad+Tools.format("\t(%.3f%% Kept)", (readsKept*100.0/totalReads)));
			s=""+totalBases;
			while(pad.length()+s.length()<9){pad+=" ";}
			outstream.println("Total bases in:  \t\t"+totalBases+pad+Tools.format("\t(%.3f%% Kept)", (basesKept*100.0/totalBases)));
		}
//		outstream.println();
		if(histogram_total!=null){
			TextStreamWriter tswh=null;
			StringBuilder sb=new StringBuilder(100);
			if(USE_HISTOGRAM){
				tswh=new TextStreamWriter(histFile, overwrite, false, false);
				tswh.start();
				tswh.print("#Depth\tRaw_Count\tUnique_Kmers\n");
			}
			int lim=(int)(HIST_LEN_PRINT-1);
			long remaining=Tools.sum(histogram_total);
			long sumRaw1=0;
			long sumRaw2=0;
			long sum1=0;
			long sum2=0;
			long sumsquare=0;
			for(int i=0; i<lim; i++){
				long x=histogram_total.get(i);
				long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
//				long y=((x)/(i<1 ? 1 : i));
				sumRaw1+=x;
				sum1+=y;
				sumsquare+=(x*Tools.max(1, i));
				if(tswh!=null){
					if(PRINT_ZERO_COVERAGE /*|| x>0*/ || y>0){
						sb.append(i).append('\t');
						sb.append(x).append('\t');
						sb.append(y).append('\n');
					}
					tswh.print(sb.toString());
					sb.setLength(0);
				}
				if(sumRaw1>=remaining){break;} //Stop once there is no more coverage, even if PRINT_ZERO_COVERAGE is not set.
			}
			for(int i=lim; i<histogram_total.length(); i++){
				long x=histogram_total.get(i);
				sumRaw2+=x;
				long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
//				long y=((x)/(i<1 ? 1 : i));
				sum2+=y;
			}
			if(tswh!=null){
				if(sumRaw2>0 || sum2>0){
					sb.append(lim).append('\t');
					sb.append(sumRaw2).append('\t');
					sb.append(sum2).append('\n');
				}
				tswh.print(sb.toString());
				tswh.poison();
				tswh.waitForFinish();
				outstream.println("Wrote histogram to "+histFile);
			}
			
			long histCount=Tools.sum(histogram_total); //Total number of kmers counted
			long halfCount=(histCount+1)/2;
			double histCountU=0; //Unique kmers counted
			long temp1=0;
			double temp2=0;
			int median_all=-1;
			int median_unique=-1;
			for(int i=0; i<histogram_total.length(); i++){
				long x=histogram_total.get(i);
				temp1+=x;
				if(temp1>=halfCount && median_all<0){median_all=i;}
//				histSum+=(x*(double)i);
				histCountU+=(x/(double)Tools.max(1, i));
			}
			double halfCount2=(histCountU)/2;
			for(int i=0; i<histogram_total.length(); i++){
				long x=histogram_total.get(i);
				temp2+=(x/Tools.max(i, 1.0));
				if(temp2>=halfCount2 && median_unique<0){
					median_unique=i;
					break;
				}
			}
			if(median_all<0){median_all=0;}
			double avg_all=sumsquare/(double)histCount;
			double avg_unique=histCount/histCountU;
			double stdev_unique=Tools.standardDeviationHistogramKmer(histogram_total);
			double stdev_all=Tools.standardDeviationHistogram(histogram_total);
			outstream.println("Total kmers counted:          \t"+(sumRaw1+sumRaw2));
			
			double uniqueC=((sum1+sum2)*100.0/(sumRaw1+sumRaw2));
			double uniqueE=((estUnique)*100.0/(sumRaw1+sumRaw2));
			double uniqueM=Tools.max(uniqueC, uniqueE);
			outstream.println("Total unique kmer count:      \t"+(sum1+sum2));
			if(CANONICAL){outstream.println("Includes forward kmers only.");}
			outstream.println("The unique kmer estimate can be more accurate than the unique count, if the tables are very full.");
			outstream.println("The most accurate value is the greater of the two.");
			outstream.println();
			
			outstream.println("Percent unique:               \t"+(uniqueM<10 ? " " : "")+Tools.format("%.2f%%", uniqueM));

			outstream.println("Depth average:                \t"+Tools.format("%.2f\t(unique kmers)", avg_unique));
			outstream.println("Depth median:                 \t"+Tools.format("%d\t(unique kmers)", median_unique));
			outstream.println("Depth standard deviation:     \t"+Tools.format("%.2f\t(unique kmers)", stdev_unique));
			
			outstream.println("\nDepth average:                \t"+Tools.format("%.2f\t(all kmers)", avg_all));
			outstream.println("Depth median:                 \t"+Tools.format("%d\t(all kmers)", median_all));
			outstream.println("Depth standard deviation:     \t"+Tools.format("%.2f\t(all kmers)", stdev_all));
		}
		
		return totalBases;
	}
	
	
	
	/**
	 * Locates and fixes spikes in a coverage profile (potentially) caused by false positives in a bloom filter.
	 * Theory:  If a high-count kmer is adjacent on both sides to low-count kmers, it may be a false positive.
	 * It could either be reduced to the max of the two flanking points or examined in more detail.
	 * @param array An array of kmer counts for adjacent kmers in a read.
	 */
	private static void fixSpikes(int[] array){
		
		for(int i=1; i<array.length-1; i++){
			long a=Tools.max(1, array[i-1]);
			int b=array[i];
			long c=Tools.max(1, array[i+1]);
			if(b>1 && b>a && b>c){
				//peak
				if((b>=2*a || b>a+2) && (b>=2*c || b>c+2)){
					//spike
					array[i]=(int)Tools.max(a, c);
				}
			}
		}
	}
	private static void fixSpikes(int[] array, long[] kmers, KCountArray kca, int k){
		if(array.length<3){return;}
		if(array[1]-array[0]>1){
			array[0]=kca.readPrecise(kmers[0], k, CANONICAL);
		}
		if(array[array.length-1]-array[array.length-2]>1){
			array[array.length-1]=kca.readPrecise(kmers[array.length-1], k, CANONICAL);
		}
		
		for(int i=1; i<array.length-1; i++){
			int b=array[i];
			if(b>1){
				long a=Tools.max(1, array[i-1]);
				long c=Tools.max(1, array[i+1]);
				long key=kmers[i];

				if(b>a && b>c){
					//peak
					if(b<6 || b>a+1 || b>c+1){
						array[i]=kca.readPreciseMin(key, k, CANONICAL);
					}
					//				if((b>=2*a || b>a+2) && (b>=2*c || b>c+2)){
					//					//spike
					//					int b1=(int)((a+c)/2);
					//					int b2=kca.readLeft(key, k, CANONICAL);
					//					int b3=kca.readRight(key, k, CANONICAL);
					//					array[i]=Tools.min(b, b1, b2, b3);
					//				}
					//				else
					//				{
					////					array[i]=kca.readPreciseMin(key, k, CANONICAL);
					//				}
				}
				//			else
				//				if(Tools.max(ada, adc)>=Tools.max(2, Tools.min((int)a, b, (int)c)/4)){
				//					array[i]=kca.readPrecise(key, k, CANONICAL);
				//				}
				//			else
				//				if(b>a+1 || b>c+1){
				//					//steep
				//					array[i]=kca.readPrecise(key, k, CANONICAL);
				//				}
			}
		}
	}
	
	
	private static void analyzeSpikes(int[] array, int width){
		if(array.length<3){return;}
		int peakcount=0, valleycount=0, spikecount=0, flatcount=0, slopecount=0;
		for(int i=1; i<array.length-1; i++){
			long a=array[i-1];
			int b=array[i];
			long c=array[i+1];
			if(b>a && b>c){
				peakcount++;
				if((b>=2*a || b>a+2) && (b>=2*c || b>c+2)){
					spikecount++;
				}
			}else if(b<a && b<c){
				valleycount++;
			}else if(b==a && b==c){
				flatcount++;
			}else{
				slopecount++;
			}
		}
		if(peakcount>0){peaks.addAndGet(peakcount);}
		if(valleycount>0){valleys.addAndGet(valleycount);}
		if(spikecount>0){spikes.addAndGet(spikecount);}
		if(flatcount>0){flats.addAndGet(flatcount);}
		if(slopecount>0){slopes.addAndGet(slopecount);}
	}
	
	public static int[] generateCoverage(Read r, KCountArray kca, int k, int[] out, long[] kmers) {
		if(k>31){return generateCoverageLong(r, kca, k, out);}
		if(r==null || r.bases==null || r.length()<k){return new int[] {0};}
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		
		int len=0;
		long kmer=0;
		final byte[] bases=r.bases;
		final int arraylen=r.length()-k+1;
		if(out==null || out.length!=arraylen){out=new int[arraylen];}
		Arrays.fill(out, -1);
		if(FIX_SPIKES){
			if(kmers==null || kmers.length!=arraylen){kmers=new long[arraylen];}
			Arrays.fill(kmers, -1);
		}
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;

				if(len>=k){
					//						int count=kca.readPrecise(kmer, k, CANONICAL);
					int count=kca.read(kmer, k, CANONICAL);
					out[i-k+1]=count;
					if(kmers!=null){kmers[i-k+1]=kmer;}
				}
			}
		}
		
		if(FIX_SPIKES){fixSpikes(out, kmers, kca, k);}
//		fixSpikes(out, 1);
		
		analyzeSpikes(out, 1);
		return out;
	}
	
	public static int[] generateCoverageLong(Read r, KCountArray kca, int k, int[] out) {
		assert(k>31);
		if(r==null || r.bases==null || r.length()<k){return new int[] {0};}
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		
		int len=0;
		long kmer=0;
		final byte[] bases=r.bases;
		
		final int arraylen=r.length()-k+1;
		if(out==null || out.length!=arraylen){out=new int[arraylen];}
		Arrays.fill(out, -1);
		
		int tailshift=k%32;
		int tailshiftbits=tailshift*2;
		
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=Long.rotateLeft(kmer, 2);
				kmer=kmer^x;
				len++;
				if(len>k){
					long x2=AminoAcid.baseToNumber[bases[i-k]];
					kmer=kmer^(x2<<tailshiftbits);
				}

				if(len>=k){
					int count=kca.read(kmer);
					out[i-k+1]=count;
				}
			}
		}
		
		fixSpikes(out);
		
		analyzeSpikes(out, 1);
		return out;
	}
	
	
	private static class ProcessThread extends Thread{
		
		ProcessThread(ConcurrentReadInputStream cris_, KCountArray kca_, int k_, ConcurrentReadOutputStream rosk_){
			cris=cris_;
			kca=kca_;
			k=k_;
			rosk=rosk_;
		}
		
		@Override
		public void run(){
			countInThread();
		}
		
		void countInThread() {
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			final ArrayList<Read> keep=new ArrayList<Read>(Shared.bufferLen());
			
			int[] cov1=null;
			long[] kmers1=null;
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(int rnum=0; rnum<reads.size(); rnum++){
					Read r=reads.get(rnum);
					Read r2=r.mate;
					assert(r!=r2);
					
					int depth=-1;
					
					int readcount=0;
					int basecount=0;
					
					int min=0;
					int max=0;
					int[] cov=null;
					long[] kmers=null;
					
					if(r!=null && r.bases!=null){
						readcount++;
						basecount+=r.length();
						if(r.length()>=k){
							if(verbose){outstream.println();}
							if(FIX_SPIKES && k<32){
								final int arraylen=r.length()-k+1;
								if(kmers1==null || kmers1.length!=arraylen){kmers1=new long[arraylen];}
								kmers=kmers1;
							}
							cov=getSortedCoverageAndIncrementHistogram(r, cov1, kmers1);
							if(cov!=null){;
								int i=cov.length-1;
								while(i>=0 && cov[i]<MIN_DEPTH){i--;}
								if(i+1>=MIN_KMERS_OVER_MIN_DEPTH){depth=cov[(int)(i*(1-DEPTH_PERCENTILE))];}
								cov1=cov;
								min=cov[cov.length-1];
								max=cov[(int)(cov.length*0.05f)];
							}
						}
					}
					

					totalReads+=readcount;
					totalBases+=basecount;
					if(max>TARGET_DEPTH && max>2*min){
						readsKept+=readcount;
						basesKept+=basecount;
						StringBuilder sb=new StringBuilder();
						sb.append(cov[0]);
						for(int i=1; i<cov.length; i++){
							sb.append('\t');
							sb.append(cov[i]);
						}
						r.obj=sb.toString();
						keep.add(r);
					}else{
						readsTossed+=readcount;
						basesTossed+=basecount;
					}
				}
				

				if(rosk!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					rosk.add(keep, ln.id);
				}
				keep.clear();
				
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		private final int[] getSortedCoverageAndIncrementHistogram(Read r, int[] cov, long[] kmers){
			assert(r!=null && r.bases!=null && r.length()>=k) : r;
			cov=generateCoverage(r, kca, k, cov, kmers);
			if(cov!=null){
				Arrays.sort(cov);
				Tools.reverseInPlace(cov);
				incrementHistogramSorted(cov);
			}
			return cov;
		}
		
		private final void incrementHistogramSorted(int[] cov){
			if(hist==null || cov==null || cov.length==0){return;}
			
//			outstream.println(Arrays.toString(cov));
			
			int last=cov[0];
			long sum=0;
//			long sum2=0;
			for(int x : cov){
//				outstream.println("Processing "+x);
				if(x<0){break;}
				int y=Tools.min(x, HIST_LEN-1);
				if(y==last){sum++;}
				else if(sum>0){
//					outstream.println("Incrementing "+last+" by "+sum);
//					sum2+=sum;
					if(last<hist.length){hist[last]+=sum;}
					else{histogram_total.addAndGet(last, sum);}
					sum=1;
				}
				last=y;
			}
//			outstream.println("Ended loop");
			if(sum>0){
//				outstream.println("Incrementing "+last+" by "+sum);
//				sum2+=sum;
				if(last<hist.length){hist[last]+=sum;}
				else{histogram_total.addAndGet(last, sum);}
			}
//			assert(sum2==cov.length) : sum2+", "+cov.length+", "+last+", "+sum;
		}
		
		private final ConcurrentReadInputStream cris;
		private final KCountArray kca;
		private final int k;
		/** Stream for kept reads */
		private final ConcurrentReadOutputStream rosk;
		public final long[] hist=new long[THREAD_HIST_LEN];//(USE_HISTOGRAM ? new long[HIST_LEN] : null);
		
		private long totalBases=0;
		private long totalReads=0;

		public long readsKept=0;
		public long readsTossed=0;
		public long basesKept=0;
		public long basesTossed=0;
	}
	
	public static PrintStream outstream=System.err;

	public static int THREAD_HIST_LEN=1<<12;
	public static int HIST_LEN=1<<20;
	public static long HIST_LEN_PRINT=HIST_LEN;
	public static boolean USE_HISTOGRAM=false;
	public static boolean PRINT_ZERO_COVERAGE=false;
	public static AtomicLongArray histogram_total;
	
	private static int THREADS=8;
	private static boolean verbose=false;
	
	
	private static int TARGET_DEPTH=50;
	private static int MAX_DEPTH=-1;
	private static int MIN_DEPTH=3;
	private static int MIN_KMERS_OVER_MIN_DEPTH=10;
	private static float DEPTH_PERCENTILE=0.5f;
	
	
	public static boolean CANONICAL=true;
	public static boolean ZERO_BIN=false;
	public static boolean FIX_SPIKES=true;
	public static boolean ordered=false;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean prefilter=false;

	public static AtomicLong peaks=new AtomicLong();
	public static AtomicLong spikes=new AtomicLong();
	public static AtomicLong flats=new AtomicLong();
	public static AtomicLong valleys=new AtomicLong();
	public static AtomicLong slopes=new AtomicLong();
}
