package cardinality;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import tracker.ReadStats;

class LogLogWrapper {
	
	public static final void main(String[] args){
		LogLogWrapper llw=new LogLogWrapper(args);
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		Timer t=new Timer();
		LongList results=new LongList();
		for(int i=0; i<llw.trials; i++){
			long card=llw.process();
			results.add(card);
		}
		
		if(llw.trials>1){
			results.sort();

			long kmers=llw.kmersProcessed;
			double mean=results.mean();
			double hmean=results.harmonicMean();
			double gmean=results.geometricMean();
			final double div=hmean;
//			long expected=(150-31+1)*llw.maxReads;
//			if(!llw.synth){expected=(long)mean;}
			
			long median=results.median();
			
			
			
			long min=results.min();
			long max=results.max();
			long p05=results.percentile(0.05);
			long p95=results.percentile(0.95);
			double stdev=results.stdev();
			double avgDif=results.avgDif(mean);
//			double rmsDif=results.rmsDif(mean);

			double range=(p95-p05)/div;
			
//			System.err.println(stdev);
//			System.err.println(rmsDif);
			
			stdev=stdev/div;
			avgDif=avgDif/div;
//			rmsDif=rmsDif/mean;
			
			/* Testing indicates that more buckets and fewer trials is more accurate. */
			/* For example, 8 buckets 256 trials < 32 buckets 64 trials < 256 buckets 8 trials, */
			/* from the standpoint of minimizing the standard deviation of the hmean of the estimates over 45 reps. */
			
			t.stopAndPrint();
			System.err.println("#Trials\tkmers\thmean\tmedian\tmin\tmax\t5%ile\t95%ile\trange\tstdev\tavgDif");
			System.err.println(llw.trials+"\t"+kmers+"\t"+String.format("%.2f", hmean)+"\t"+median+"\t"+min+"\t"+max
					+"\t"+p05+"\t"+p95+"\t"+String.format("%.5f", range)+"\t"+String.format("%.5f", stdev)+"\t"+String.format("%.5f", avgDif));
			
			if(CardinalityTracker.trackCounts){
				System.err.println("Avg Count:\t"+String.format("%.4f", llw.countSum/(llw.trials*(double)llw.buckets)));
			}
//			System.err.println("#Trials\tkmers\tmean\thmean\tgmean\tmedian\tmin\tmax\t5%ile\t95%ile\trange\tstdev\tavgDif");
//			System.err.println(llw.trials+"\t"+kmers+"\t"+String.format("%.2f", mean)+"\t"+String.format("%.2f", hmean)+"\t"+String.format("%.2f", gmean)+"\t"+median+"\t"+min+"\t"+max
//					+"\t"+p05+"\t"+p95+"\t"+String.format("%.5f", range)+"\t"+String.format("%.5f", stdev)+"\t"+String.format("%.5f", avgDif));
			
		}
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
	}
	
	public LogLogWrapper(String[] args){

		Shared.capBufferLen(200);
		Shared.capBuffers(8);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("buckets") || a.equals("loglogbuckets")){
				buckets=Parse.parseIntKMG(b);
			}else if(a.equals("k") || a.equals("loglogk")){
				k=Integer.parseInt(b);
			}else if(a.equals("seed") || a.equals("loglogseed")){
				seed=Long.parseLong(b);
			}else if(a.equals("seed2")){
				seed2=Long.parseLong(b);
			}else if(a.equals("minprob") || a.equals("loglogminprob")){
				minProb=Float.parseFloat(b);
			}else if(a.equals("synth")){
				synth=Parse.parseBoolean(b);
			}else if(a.equals("trials")){
				trials=Parse.parseIntKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("loglogcounts") || a.equals("loglogcount") || 
					a.equals("count") || a.equals("counts") || a.equals("trackcounts")){
				CardinalityTracker.trackCounts=Parse.parseBoolean(b);
			}else if(a.equals("atomic")){
				assert(false) : "Atomic flag disabled.";
//				CardinalityTracker.atomic=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=b;
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;

			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=(parser.in1==null ? null : parser.in1.split(","));
			in2=(parser.in2==null ? null : parser.in2.split(","));
			out=parser.out1;
		}
		
		assert(synth || (in1!=null && in1.length>0)) : "No primary input file specified.";
		if(synth){
			ffin1=ffin2=null;
		}else{
			ffin1=new FileFormat[in1.length];
			ffin2=new FileFormat[in1.length];
			
			for(int i=0; i<in1.length; i++){
				String a=in1[i];
				String b=(in2!=null && in2.length>i ? in2[i] : null);
				assert(a!=null) : "Null input filename.";
				if(b==null && a.indexOf('#')>-1 && !new File(a).exists()){
					b=a.replace("#", "2");
					a=a.replace("#", "1");
				}

				ffin1[i]=FileFormat.testInput(a, FileFormat.FASTQ, null, true, true);
				ffin2[i]=FileFormat.testInput(b, FileFormat.FASTQ, null, true, true);
			}
		}

		threads=Shared.threads();
		assert(FastaReadInputStream.settingsOK());
	}
	
	
	long process(){
		Timer t=new Timer();
		readsProcessed=basesProcessed=kmersProcessed=0;
		
		CardinalityTracker log=CardinalityTracker.makeTracker(buckets, k, seed, minProb);
		
		for(int ffnum=0, max=(synth ? 1 : ffin1.length); ffnum<max; ffnum++){
			ConcurrentReadInputStream cris=null;
			if(!synth){
				cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, ffin1[ffnum], ffin2[ffnum]);
				cris.start();
			}

			LogLogThread[] threadArray=new LogLogThread[threads];
			for(int tid=0; tid<threadArray.length; tid++){
				threadArray[tid]=new LogLogThread((CardinalityTracker.atomic ? log : CardinalityTracker.makeTracker(buckets, k, seed, minProb)), cris, tid);
			}
			for(LogLogThread llt : threadArray){
				llt.start();
			}
			for(LogLogThread llt : threadArray){
				while(llt.getState()!=Thread.State.TERMINATED){
					try {
						llt.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				readsProcessed+=llt.readsProcessedT;
				basesProcessed+=llt.basesProcessedT;
				kmersProcessed+=llt.kmersProcessedT;
				if(!CardinalityTracker.atomic){log.add(llt.log);}
			}

			if(cris!=null){errorState|=ReadWrite.closeStreams(cris);}
		}
		
//		final int[] max=new int[buckets];
//		if(CardinalityTracker.atomic){
//			for(int i=0; i<log.maxArray.length(); i++){
//				//				System.err.println(log.maxArray.get(i));
//				max[i]=log.maxArray.get(i);
//			}
//		}
		
		t.stop();
		
		
		long cardinality=log.cardinality();
		countSum+=(CardinalityTracker.trackCounts ? log.countSum() : 0);
		
		if(out!=null){
			ReadWrite.writeString(cardinality+"\n", out);
		}
		
		if(!Parser.silent) {
//		Arrays.sort(copy);
//		System.err.println("Median:        "+copy[Tools.median(copy)]);
		
//		System.err.println("Mean:          "+Tools.mean(copy));
//		System.err.println("Harmonic Mean: "+Tools.harmonicMean(copy));
		System.err.println("Cardinality:   "+cardinality);
//		System.err.println("CardinalityH:  "+log.cardinalityH());
		
//		for(long i : log.counts){System.err.println(i);}
		
		System.err.println("Time: \t"+t);
		}
		
		return cardinality;
	}
	
	private static Read makeRead(int len, Random randy, Read r){
		if(r==null || r.bases==null || r.bases.length!=len){
			r=new Read(null, null, 0);
			r.bases=new byte[len];
		}
		byte[] bases=r.bases;
		
		int pos=0;
		final int basesPerRand=4;//Fewer calls to rand should be faster
		for(int max=bases.length-(bases.length%basesPerRand); pos<max;){
			int x=randy.nextInt()%prime;
			for(int i=0; i<basesPerRand; i++){
				int num=x&3;
				byte b=AminoAcid.numberToBase[num];
				bases[pos]=b;
				pos++;
				x>>=2;
			}
		}
		for(; pos<bases.length; pos++){
			int x=randy.nextInt()%prime;
			int num=x&3;
			byte b=AminoAcid.numberToBase[num];
			bases[pos]=b;
		}
		return r;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class LogLogThread extends Thread{
		
		LogLogThread(CardinalityTracker log_, ConcurrentReadInputStream cris_, int tid_){
			log=log_;
			cris=cris_;
			tid=tid_;
		}
		
		@Override
		public void run(){
			if(cris!=null){runCris();}
			else{runSynth();}
		}
		
		public void runCris(){
			final int kt=k;
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(Read r : reads){
//					if(!r.validated()){r.validate(true);}
//					if(r.mate!=null && !r.mate.validated()){r.mate.validate(true);}
					log.hash(r);
					readsProcessedT+=r.pairCount();
					basesProcessedT+=r.pairLength();
					kmersProcessedT+=r.numPairKmers(kt);
				}
				
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}
		
		public void runSynth(){
			final int kt=k;
			assert(maxReads>0 && maxReads<Long.MAX_VALUE);
			long readsLeft=maxReads/threads;
			readsLeft+=(maxReads%threads>tid ? 1 : 0);
			Random randy=Shared.threadLocalRandom(seed2<0 ? seed2 : seed2+999999L);
			
			Read r=null;
			while(readsLeft>0){
				r=makeRead(150, randy, r);
				log.hash(r);
				readsProcessedT+=r.pairCount();
				basesProcessedT+=r.pairLength();
				kmersProcessedT+=r.numPairKmers(kt);
				readsLeft--;
			}
		}
		
		private final CardinalityTracker log;
		private final ConcurrentReadInputStream cris;
		private final int tid;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long kmersProcessedT=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int buckets=2048;
	private int k=31;
	private long seed=-1;
	private long seed2=-1;
	private float minProb=0;
	private long countSum=0;
	
	private String[] in1=null;
	private String[] in2=null;
	private String out=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long kmersProcessed=0;
	
	private long maxReads=-1;
	
	boolean overwrite=true;
	boolean append=false;
	boolean errorState=false;

	private int trials=1;
	private boolean synth=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat[] ffin1;
	private final FileFormat[] ffin2;
	
	final int threads;
	private static final int prime=32452843; //A prime number
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
}
