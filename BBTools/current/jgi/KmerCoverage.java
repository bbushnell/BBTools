package jgi;

import java.io.File;
import java.io.PrintStream;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

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
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;


/**
 * @author Brian Bushnell
 * @date Oct 11, 2012
 *
 */
public class KmerCoverage {

	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		String in1=(args[0].indexOf("=")>0 ? null : args[0]);
		String in2=(in1!=null && args.length>1 ? args[1] : null);
		if(in2!=null && "null".equalsIgnoreCase(in2)){in2=null;}
		
		{
			if(in1!=null && !in1.contains(",")){
				File f=new File(in1);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(in1+" does not exist.");}
			}
			if(in2!=null && !in2.contains(",")){
				File f=new File(in2);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(in2+" does not exist.");}
				if(in1.equalsIgnoreCase(in2)){
					throw new RuntimeException("Both input files are the same.");
				}
			}
		}
		
		KmerCountAbstract.minQuality=4;
		KmerCountAbstract.minProb=0.1f;
		KmerCountAbstract.CANONICAL=true;
		
		int k=31;
		int cbits=16;
//		int gap=0;
		int hashes=4;
		long cells=-1;
		long maxReads=-1;
		int buildpasses=1;
		long tablereads=-1; //How many reads to process when building the hashtable
		int buildStepsize=4;
		String output=null;
		int prehashes=-1;
		long precells=-1;
		String histFile=null;
		int threads=-1;
		int minq=KmerCountAbstract.minQuality;
		boolean auto=true;
		List<String> extra=null;
		Parser parser=new Parser();
		
		long memory=Runtime.getRuntime().maxMemory();
		
		for(int i=(in1==null ? 0 : 1); i<args.length; i++){
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
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
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
			}else if(a.startsWith("minmedian")){
				MIN_MEDIAN=Integer.parseInt(b);
			}else if(a.startsWith("minaverage")){
				MIN_AVERAGE=Integer.parseInt(b);
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
				OUTPUT_ATTACHMENT=Parse.parseBoolean(b);
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("tablereads") || a.startsWith("buildreads")){
				tablereads=Parse.parseKMG(b);
			}else if(a.startsWith("out")){
				output=b;
			}else if(a.startsWith("hist")){
				histFile=b;
			}else if(a.startsWith("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.startsWith("datainheader") || a.startsWith("datainname") || a.startsWith("useheader")){
				USE_HEADER=Parse.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("auto") || a.equals("automatic")){
				auto=Parse.parseBoolean(b);
			}else if(a.equals("samplewhenreadingtable") || a.equals("sampleoutput")){
				DONT_SAMPLE_OUTPUT=!Parse.parseBoolean(b);
			}else if(a.equals("dontsamplewhenreadingtable") || a.equals("dontsampleoutput")){
				DONT_SAMPLE_OUTPUT=Parse.parseBoolean(b);
			}else if(a.startsWith("kmersample")){
				kmersamplerate=Integer.parseInt(b);
//				KmerCountAbstract.kmersamplerate=kmersamplerate;
			}else if(a.startsWith("sample") || a.startsWith("readsample")){
				readsamplerate=Integer.parseInt(b);
//				KmerCountAbstract.readsamplerate=readsamplerate;
			}else if(a.startsWith("canonical")){
				CANONICAL=KmerCountAbstract.CANONICAL=Parse.parseBoolean(b);
			}else if(a.startsWith("fixspikes")){
				FIX_SPIKES=Parse.parseBoolean(b);
			}else if(a.equals("printzerocoverage") || a.equals("pzc")){
				PRINT_ZERO_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("removeduplicatekmers") || a.equals("rdk")){
				KmerCountAbstract.KEEP_DUPLICATE_KMERS=!Parse.parseBoolean(b);
			}else if(a.startsWith("extra")){
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
		
		assert(FastaReadInputStream.settingsOK());
		if(k>31){CANONICAL=KmerCountAbstract.CANONICAL=false;}
		assert(CANONICAL==KmerCountAbstract.CANONICAL);
		
//		if(output!=null && reads1.contains(",")){
//			throw new RuntimeException("\nLists of input files can only be used with histogram output, not full output.\n" +
//					"Please set output=null or move the extra input files to 'extra=file1,file2,...fileN'");
//		}
		
		{
			if(histFile==null){
				
				
			}else{
				USE_HISTOGRAM=true;
			}
			
			final int maxCount=(int)(cbits>16 ? Integer.MAX_VALUE : (1L<<cbits)-1);
			assert(maxCount>0);
			HIST_LEN_PRINT=Tools.max(1, Tools.min(HIST_LEN_PRINT, maxCount+1));
			assert(HIST_LEN_PRINT<=Integer.MAX_VALUE) : HIST_LEN_PRINT+", "+Integer.MAX_VALUE;
			HIST_LEN=(int)Tools.min(maxCount+1, Tools.max(HIST_LEN_PRINT, HIST_LEN));
			
			histogram_total=new long[HIST_LEN];
		}
		
		if(extra!=null){
			for(String s : extra){
				File f=new File(s);
				if(!f.exists() || !f.isFile()){throw new RuntimeException(s+" does not exist.");}
				assert(!s.equalsIgnoreCase(in1) && (in2==null || !s.equalsIgnoreCase(in2))) : "\nInput file "+s+" should not be included as an extra file.\n";
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
			long mem=usable-(USE_HISTOGRAM ? (HIST_LEN*8*(THREADS+1)) : 0);
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
		outstream.println();
		ReadCounter rc=new ReadCounter(k, true, false, false, Shared.AMINO_IN);
		rc.detectStepsize=buildStepsize;
		if(prefilter){
			prefilterArray=rc.makeKca(in1, in2, extra, 2, precells, prehashes, minq,
					tablereads, 1, 1, 1, null, 0);
			outstream.println("Made prefilter:   \t"+prefilterArray.toShortString(prehashes));
		}
		kca=rc.makeKca(in1, in2, extra, cbits, cells, hashes, minq,
				tablereads, buildpasses, 2, 2, prefilterArray, (prefilterArray==null ? 0 : prefilterArray.maxValue));
		ht.stop();
		
		outstream.println("Made hash table:  \t"+kca.toShortString(hashes));
		
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
//			System.out.println("Used cells: "+kca.cellsUsed(1));
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
		
		if(in1!=null && in1.contains(",") && !new File(in1).exists()){
			String[] list1=in1.split(",");
			String[] list2=(in2==null ? null : in2.split(","));
			bases=count(list1, list2, kca, k, maxReads, output, ordered, overwrite, histFile, estUnique);
		}else{
			bases=count(in1, in2, kca, k, maxReads, output, ordered, overwrite, histFile, estUnique);
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
		
		System.err.println("\nDepth Topology\t");
		System.err.println("Spikes:     \t\t\t"+(dsp<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dsp,sp));
		System.err.println("Peaks:      \t\t\t"+(dpe<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dpe,pe));
		System.err.println("Valleys:    \t\t\t"+(dva<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dva,va));
		System.err.println("Slopes:     \t\t\t"+(dsl<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dsl,sl));
		System.err.println("Flats:      \t\t\t"+(dfl<10 ? " " : "")+Tools.format("%.3f%%  \t%d",dfl,fl));
	}


	public static long count(String reads1, String reads2, KCountArray kca, int k, long maxReads,
			String output, boolean ordered, boolean overwrite, String histFile, long estUnique) {
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		assert(cris!=null) : reads1;
		
		if(fileIO.FileFormat.hasFastaExtension(reads1)){
			ADD_CARROT=false;
		}
		
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		final ConcurrentReadOutputStream ros;
		if(output!=null){
			final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
			
			String out1=output.replaceFirst("#", "1");
			String out2=null;
			
			if(cris.paired()){
				if(output.contains("#")){
					out2=output.replaceFirst("#", "2");
				}else{
					outstream.println("Writing interleaved.");
				}
			}

			assert(!out1.equalsIgnoreCase(reads1) && !out1.equalsIgnoreCase(reads1));
			assert(out2==null || (!out2.equalsIgnoreCase(reads1) && !out2.equalsIgnoreCase(reads2)));
			
			FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, OUTPUT_ATTACHMENT ? "attachment" : null, true, overwrite, append, ordered);
			FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, OUTPUT_ATTACHMENT ? "attachment" : null, true, overwrite, append, ordered);
			ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			
			ros.start();
			outstream.println("Started output threads.");
		}else{
			ros=null;
		}
		
		long bases=calcCoverage(cris, kca, k, maxReads, ros, histFile, overwrite, estUnique);
		
		ReadWrite.closeStreams(cris, ros);
		if(verbose){System.err.println("Closed stream");}
		return bases;
	}
	
	
	public static long count(String[] list1, String[] list2, KCountArray kca, int k, long maxReads,
			String output, boolean ordered, boolean overwrite, String histFile, long estUnique) {
		
		ConcurrentReadOutputStream ros=null;
		String[] out1=null, out2=null;
		

		final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
		if(output!=null){
			if(!new File(output).exists()){
				out1=output.split(",");
			}else{
				out1=new String[] {output};
			}
			out2=new String[out1.length];
			for(int i=0; i<out1.length; i++){
				if(out1[i].contains("#")){
					out2[i]=out1[i].replaceFirst("#", "2");
					out1[i]=out1[i].replaceFirst("#", "1");
				}
			}
		}
		
		long bases=0;
		
		for(int x=0; x<list1.length; x++){
			
			if(out1!=null){
				if(x==0 || out1.length>1){
					if(ros!=null){
						ReadWrite.closeStream(ros);
					}
					
					FileFormat ff1=FileFormat.testOutput(out1[x], FileFormat.FASTQ, OUTPUT_ATTACHMENT ? "attachment" : null, true, overwrite, append, ordered);
					FileFormat ff2=out2==null ? null : FileFormat.testOutput(out2[x], FileFormat.FASTQ, OUTPUT_ATTACHMENT ? "attachment" : null, true, overwrite, append, ordered);
					ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
					
					ros.start();
					outstream.println("Started output threads.");
				}else if(ros!=null){
					ros.resetNextListID();
				}
			}
				
			String reads1=list1[x];
			String reads2=(list2==null || list2.length<=x ? null : list2[x]);

			final ConcurrentReadInputStream cris;
			{
				FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
				FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
				if(verbose){System.err.println("Started cris");}
				cris.start(); //4567
				if(ff1.fasta()){ADD_CARROT=false;}
			}
			boolean paired=cris.paired();
			if(verbose){System.err.println("Paired: "+paired);}

			bases+=calcCoverage(cris, kca, k, maxReads, ros, histFile, overwrite, estUnique);

			ReadWrite.closeStream(cris);
			if(verbose){System.err.println("Closed stream");}
			
		}

		//Wait until threads finish!

		ReadWrite.closeStream(ros);

		return bases;
	}
	

	
	public static long calcCoverage(ConcurrentReadInputStream cris, KCountArray kca, int k, long maxReads, ConcurrentReadOutputStream ros,
			String histFile, boolean overwrite, long estUnique) {
		Timer tdetect=new Timer();
		tdetect.start();
		
		long totalBases=0;
		long totalReads=0;
		
//		assert(false) : THREADS;
		ProcessThread[] pta=new ProcessThread[THREADS];
		for(int i=0; i<pta.length; i++){
			pta[i]=new ProcessThread(cris, kca, k, ros);
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
				
				for(int j=0; j<histogram_total.length; j++){
					histogram_total[j]+=ct.hist[j];
				}
			}
		}
		
		if(!ZERO_BIN && histogram_total!=null && histogram_total.length>1){
			histogram_total[1]+=histogram_total[0];
			histogram_total[0]=0;
		}
		
//		outstream.println();
		tdetect.stop();
		outstream.println("Table read time: \t\t"+tdetect+"   \t"+Tools.format("%.2f", totalBases*1000000.0/(tdetect.elapsed))+" kb/sec");
		outstream.println("Total reads:     \t\t"+totalReads);
		outstream.println("Total bases:     \t\t"+totalBases);
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
			long remaining=shared.Vector.sum(histogram_total);
			long sumRaw1=0;
			long sumRaw2=0;
			long sum1=0;
			long sum2=0;
			long sumsquare=0;
			for(int i=0; i<lim; i++){
				long x=histogram_total[i];
				long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
//				long y=((x)/(i<1 ? 1 : i));
				sumRaw1+=x;
				sum1+=y;
				sumsquare+=(x*Tools.max(1, i));
				if(tswh!=null){
					if(PRINT_ZERO_COVERAGE || x>0 || y>0){
						sb.append(i).append('\t');
						sb.append(x).append('\t');
						sb.append(y).append('\n');
					}
					tswh.print(sb.toString());
					sb.setLength(0);
				}
				if(sumRaw1>=remaining){break;} //Stop once there is no more coverage, even if PRINT_ZERO_COVERAGE is not set.
			}
			for(int i=lim; i<histogram_total.length; i++){
				long x=histogram_total[i];
				sumRaw2+=x;
				long y=((x+i/2)/(i<1 ? 1 : i)); //x+i/2 rounds to compensate for colliding kmers being put in an overly high bin
//				long y=((x)/(i<1 ? 1 : i));
				sum2+=y;
			}
			if(tswh!=null){
				if(PRINT_ZERO_COVERAGE || sumRaw2>0 || sum2>0){
					sb.append(lim).append('\t');
					sb.append(sumRaw2).append('\t');
					sb.append(sum2).append('\n');
				}
				tswh.print(sb.toString());
				tswh.poison();
				tswh.waitForFinish();
				outstream.println("Wrote histogram to "+histFile);
			}
			
			long histCount=shared.Vector.sum(histogram_total); //Total number of kmers counted
			long halfCount=(histCount+1)/2;
			double histCountU=0; //Unique kmers counted
			long temp1=0;
			double temp2=0;
			int median_all=-1;
			int median_unique=-1;
			for(int i=0; i<histogram_total.length; i++){
				long x=histogram_total[i];
				temp1+=x;
				if(temp1>=halfCount && median_all<0){median_all=i;}
//				histSum+=(x*(double)i);
				histCountU+=(x/(double)Tools.max(1, i));
			}
			double halfCount2=(histCountU)/2;
			for(int i=0; i<histogram_total.length; i++){
				long x=histogram_total[i];
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
	
	public static int[] generateCoverage(Read r, KCountArray kca, int k) {
		if(k>31){return generateCoverageLong(r, kca, k);}
		if(r==null || r.bases==null || r.length()<k){return new int[] {0};}
		
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		
		int len=0;
		long kmer=0;
		final byte[] bases=r.bases;
		final int[] out;
		final long[] kmers=(FIX_SPIKES ? new long[r.length()-k+1] : null);
		
		if(kmersamplerate<2 || DONT_SAMPLE_OUTPUT){
			out=new int[r.length()-k+1];
			Arrays.fill(out, -1);
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
		}else{
			out=new int[(r.length()-k+1+(kmersamplerate-1))/kmersamplerate];
			Arrays.fill(out, -1);
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				int x=AminoAcid.baseToNumber[b];
				if(x<0){
					len=0;
					kmer=0;
				}else{
					kmer=((kmer<<2)|x)&mask;
					len++;

					if(len>=k && i%kmersamplerate==0){
//						int count=kca.readPrecise(kmer, k, CANONICAL);
						int count=kca.read(kmer, k, CANONICAL);
						out[(i-k+1)/kmersamplerate]=count;
						if(kmers!=null){kmers[(i-k+1)/kmersamplerate]=kmer;}
					}
				}
			}
		}
		if(FIX_SPIKES){fixSpikes(out, kmers, kca, k);}
//		fixSpikes(out, 1);
		
		analyzeSpikes(out, 1);
		return out;
	}
	
	public static int[] generateCoverageLong(Read r, KCountArray kca, int k) {
		if(r==null || r.bases==null || r.length()<k){return new int[] {0};}
		
		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		
		int len=0;
		long kmer=0;
		final byte[] bases=r.bases;
		final int[] out;
		
		int tailshift=k%32;
		int tailshiftbits=tailshift*2;
		
		if(kmersamplerate<2 || DONT_SAMPLE_OUTPUT){
			out=new int[r.length()-k+1];
			Arrays.fill(out, -1);
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
		}else{
			out=new int[(r.length()-k+1+(kmersamplerate-1))/kmersamplerate];
			Arrays.fill(out, -1);
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
					
					if(len>=k && i%kmersamplerate==0){
						int count=kca.read(kmer);
						out[(i-k+1)/kmersamplerate]=count;
					}
				}
			}
		}
		fixSpikes(out);
		
		analyzeSpikes(out, 1);
		return out;
	}
	
	
	private static class ProcessThread extends Thread{
		
		ProcessThread(ConcurrentReadInputStream cris_, KCountArray kca_, int k_, ConcurrentReadOutputStream ros_){
			cris=cris_;
			kca=kca_;
			k=k_;
			ros=ros_;
		}
		
		@Override
		public void run(){
			countInThread();
		}
		
		void countInThread() {
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(int rnum=0; rnum<reads.size(); rnum++){
					Read r=reads.get(rnum);
					Read r2=r.mate;
					
					if(DONT_SAMPLE_OUTPUT || r.numericID%readsamplerate==0){
						boolean toss1=false;
						boolean toss2=false;
						{
							totalReads++;
							if(verbose){outstream.println();}
							if(OUTPUT_ATTACHMENT && ros!=null){
//								assert(false) : ros.FASTA+", "+ros.FASTQ+", "+ros.ATTACHMENT;
								r.obj=(ros.ff1.fastq() ? toFastqString(r) : toFastaString(r));
								toss1=r.discarded();
							}else{
								 int[] cov=getCoverageAndIncrementHistogram(r);
								 if(cov==null){toss1=true;}
								 else{
									 Arrays.sort(cov);
									 toss1=(cov[cov.length/2]<MIN_MEDIAN && Tools.averageInt(cov)<MIN_AVERAGE);
								 }
							}
						}
						if(r2!=null){
							totalReads++;
							if(verbose){outstream.println();}
							if(OUTPUT_ATTACHMENT && ros!=null){
								r2.obj=(ros.ff1.fastq() ? toFastqString(r2) : toFastaString(r2));
								toss2=r.discarded();
							}else{
								 int[] cov=getCoverageAndIncrementHistogram(r2);
								 if(cov==null){toss2=true;}
								 else{
									 Arrays.sort(cov);
									 toss2=(cov[cov.length/2]<MIN_MEDIAN && Tools.averageInt(cov)<MIN_AVERAGE);
								 }
							}
						}
						if(toss1 && (toss2 || r2==null)){reads.set(rnum, null);}
					}
				}
				
				if(ros!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
//					System.err.println("Adding list "+ln.id+" of length "+reads.size());
					ros.add(reads, ln.id);
				}
				
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(verbose){System.err.println("Finished reading");}
			cris.returnList(ln);
			if(verbose){System.err.println("Returned list");}
		}
		
		private int[] getCoverageAndIncrementHistogram(Read r){
			if(r.bases==null || r.length()<k){
				return null;
			}else{
				totalBases+=r.length();

				int[] cov=generateCoverage(r, kca, k);
				
				if(hist!=null){
					for(int i=0; i<cov.length; i++){
						int x=Tools.min(cov[i], HIST_LEN-1);
						if(x>=0){hist[x]++;}
					}
				}
				return cov;
			}
		}
		
		private String toFastaString(Read r){
			if(r.bases==null || r.length()<k){
				if(MIN_MEDIAN>0 || MIN_AVERAGE>0){r.setDiscarded(true);}
				if(USE_HEADER){
					return (ADD_CARROT ? ">" : "")+r.id+";0;0 0 0 0 0\n"+r.bases==null ? "" : new String(r.bases);
				}else{
					return (ADD_CARROT ? ">" : "")+r.id+"\n"+(r.bases==null ? "" : new String(r.bases))+"\n0\n0 0 0 0 0";
				}
			}else{
				totalBases+=r.length();

				int[] cov=generateCoverage(r, kca, k);
				
				if(hist!=null){
					for(int i=0; i<cov.length; i++){
						int x=Tools.max(0, Tools.min(cov[i], HIST_LEN-1));
						hist[x]++;
					}
				}
				
				StringBuilder sb=new StringBuilder(cov.length*4+r.length()+(r.id==null ? 4 : r.id.length())+10);
				
				if(USE_HEADER){
					if(ADD_CARROT || r.id.charAt(0)!='>'){sb.append('>');}
					sb.append(r.id).append(';');
					
					int min=cov[0], max=cov[0], sum=0;
					for(int i=0; i<cov.length; i++){
						sb.append(cov[i]+" ");
						min=Tools.min(min, cov[i]);
						max=Tools.max(max, cov[i]);
						sum+=cov[i];
					}

					sb.append(';');
					Arrays.sort(cov);
					int median=cov[cov.length/2];
					sb.append(median).append(' ');
					sb.append(Tools.format("%.3f ", sum/(float)cov.length));
					sb.append(Tools.format("%.3f ", Tools.standardDeviation(cov)));
					sb.append(min).append(' ');
					sb.append(max).append('\n');
					
					sb.append(new String(r.bases));
					
					if(median<MIN_MEDIAN || sum/cov.length<MIN_AVERAGE){r.setDiscarded(true);}
				}else{

					if(ADD_CARROT || r.id.charAt(0)!='>'){sb.append('>');}
					sb.append(r.id).append('\n');
					sb.append(new String(r.bases)).append('\n');

					int min=cov[0], max=cov[0], sum=0;
					for(int i=0; i<cov.length; i++){
						sb.append(cov[i]+" ");
						min=Tools.min(min, cov[i]);
						max=Tools.max(max, cov[i]);
						sum+=cov[i];
					}

					sb.append('\n');
					Arrays.sort(cov);
					int median=cov[cov.length/2];
					sb.append(median).append(' ');
					sb.append(Tools.format("%.3f ", sum/(float)cov.length));
					sb.append(Tools.format("%.3f ", Tools.standardDeviation(cov)));
					sb.append(min).append(' ');
					sb.append(max);
					
					if(median<MIN_MEDIAN || sum/cov.length<MIN_AVERAGE){r.setDiscarded(true);}
				}
				return sb.toString();
			}
		}
		
		private ByteBuilder toFastqString(Read r){
			ByteBuilder sb=r.toFastq();
			if(r.bases==null || r.length()<k){
				if(MIN_MEDIAN>0 || MIN_AVERAGE>0){r.setDiscarded(true);}
				sb.append("\n0\n0 0 0 0 0");
				return sb;
			}else{
				totalBases+=r.length();

				int[] cov=generateCoverage(r, kca, k);
				
				if(hist!=null){
					for(int i=0; i<cov.length; i++){
						int x=Tools.max(0, Tools.min(cov[i], HIST_LEN-1));
						assert(x>=0) : i+", "+cov[i]+", "+HIST_LEN;
						hist[x]++;
					}
				}
				sb.append('\n');

				int min=cov[0], max=cov[0], sum=0;
				for(int i=0; i<cov.length; i++){
					sb.append(cov[i]+" ");
					min=Tools.min(min, cov[i]);
					max=Tools.max(max, cov[i]);
					sum+=cov[i];
				}

				sb.append('\n');
				Arrays.sort(cov);
				int median=cov[cov.length/2];
				sb.append(median).append(' ');
				sb.append(Tools.format("%.3f ", sum/(float)cov.length));
				sb.append(Tools.format("%.3f ", Tools.standardDeviation(cov)));
				sb.append(min).append(' ');
				sb.append(max);

				if(median<MIN_MEDIAN || sum/cov.length<MIN_AVERAGE){r.setDiscarded(true);}
				return sb;
			}
		}
		
		private final ConcurrentReadInputStream cris;
		private final KCountArray kca;
		private final int k;
		private final ConcurrentReadOutputStream ros;
		public final long[] hist=new long[HIST_LEN];//(USE_HISTOGRAM ? new long[HIST_LEN] : null);
		
		private long totalBases=0;
		private long totalReads=0;
		
	}
	
	public static PrintStream outstream=System.err;

	
	public static int HIST_LEN=1<<14;
	public static long HIST_LEN_PRINT=HIST_LEN;
	public static boolean USE_HISTOGRAM=false;
	public static boolean PRINT_ZERO_COVERAGE=false;
	public static long[] histogram_total;
	
	private static int THREADS=8;
	private static boolean verbose=false;
	private static boolean USE_HEADER=false;

	private static boolean ADD_CARROT=false;
	private static boolean OUTPUT_ATTACHMENT=true;
	private static int MIN_MEDIAN=0;
	private static int MIN_AVERAGE=0;

	public static int kmersamplerate=1;
	public static int readsamplerate=1;
	public static boolean DONT_SAMPLE_OUTPUT=false;
	public static boolean CANONICAL=true;
	public static boolean ZERO_BIN=false;
	public static boolean FIX_SPIKES=true;
	public static boolean ordered=true;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean prefilter=false;

	public static AtomicLong peaks=new AtomicLong();
	public static AtomicLong spikes=new AtomicLong();
	public static AtomicLong flats=new AtomicLong();
	public static AtomicLong valleys=new AtomicLong();
	public static AtomicLong slopes=new AtomicLong();
}
