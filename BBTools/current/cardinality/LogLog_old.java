package cardinality;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.Dedupe;
import shared.Parse;
import shared.Parser;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;
import ukmer.Kmer;

/**
 * @author Brian Bushnell
 * @date Sep 30, 2015
 *
 */
public class LogLog_old {
	
	/** Create a LogLog with default parameters */
	public LogLog_old(){
		this(1999, 8, 31, -1, 0);
	}
	
	/** Create a LogLog with parsed parameters */
	public LogLog_old(Parser p){
		this(p.loglogbuckets, p.loglogbits, p.loglogk, p.loglogseed, p.loglogMinprob);
	}
	
	/**
	 * Create a LogLog with specified parameters
	 * @param buckets_ Number of buckets (counters)
	 * @param bits_ Bits hashed per cycle
	 * @param k_ Kmer length
	 * @param seed Random number generator seed; -1 for a random seed
	 * @param minProb_ Ignore kmers with under this probability of being correct
	 */
	public LogLog_old(int buckets_, int bits_, int k_, long seed, float minProb_){
//		hashes=hashes_;
//		if((buckets_&1)==0){buckets_=(int)Primes.primeAtLeast(buckets_);}
		buckets=buckets_;
		assert(Integer.bitCount(buckets)==1) : "Buckets must be a power of 2: "+buckets;
		bucketMask=buckets-1;
		bits=bits_;
		k=Kmer.getKbig(k_);
		minProb=minProb_;
		//assert(atomic);
		maxArrayA=(atomic ? new AtomicIntegerArray(buckets) : null);
		maxArray=(atomic ? null : new int[buckets]);
		steps=(63+bits)/bits;
		tables=new long[numTables][][];
		for(int i=0; i<numTables; i++){
			tables[i]=makeCodes(steps, bits, (seed<0 ? -1 : seed+i));
		}
		
//		assert(false) : "steps="+steps+", "+tables.length+", "+tables[0].length+", "+tables[0][0].length;
	}
	
	public static void main(String[] args){
		LogLogWrapper llw=new LogLogWrapper(args);
		
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		llw.process();
		
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
	}
	
//	public final long cardinality(boolean weighted){
//		double mult=0.7947388;
//		if(weighted){mult=0.7600300;}
//		return cardinality(mult);
//	}
	
	public final long cardinality(){
		return cardinality(0.7947388);
	}
	
	public final long cardinality(double mult){
		long sum=0;
		//assert(atomic);
		if(atomic){
			for(int i=0; i<maxArrayA.length(); i++){
				sum+=maxArrayA.get(i);
			}
		}else{
			for(int i=0; i<maxArray.length; i++){
				sum+=maxArray[i];
			}
		}
		double mean=sum/(double)buckets;
		long cardinality=(long)((((Math.pow(2, mean)-1)*buckets*SKIPMOD))/1.258275);
		lastCardinality=cardinality;
		return cardinality;
	}
	
	public final long cardinalityH(){
		double sum=0;
		for(int i=0; i<maxArrayA.length(); i++){
			int x=Tools.max(1, maxArrayA.get(i));
			sum+=1.0/x;
		}
		double mean=buckets/sum;
		return (long)((Math.pow(2, mean)*buckets*SKIPMOD));
	}
	
//	public long hashOld(final long value0, final long[][] table){
//		long value=value0, code=value0;
//		long mask=(bits>63 ? -1L : ~((-1L)<<bits));
//		
//		for(int i=0; i<steps; i++){
//			int x=(int)(value&mask);
//			value>>=bits;
//			code=Long.rotateLeft(code^table[i][x], 3);
//		}
//		return Long.rotateLeft(code, (int)(value0&31));
//	}
	
	public long hash(final long value0, final long[][] table){
		long value=value0, code=0;
		long mask=(bits>63 ? -1L : ~((-1L)<<bits));

		for(int i=0; i<steps; i++){//I could also do while value!=0
			int x=(int)(value&mask);
			value>>=bits;
			code=code^table[i][x];
		}
		return code;
	}
	
	public void add(long number){
		hash(number);
	}
	
	public void hash(Read r){
		if(r==null){return;}
		if(r.length()>=k){hash(r.bases, r.quality);}
		if(r.mateLength()>=k){hash(r.mate.bases, r.mate.quality);}
	}
	
	public void hash(byte[] bases, byte[] quals){
		if(k<32){hashSmall(bases, quals);}
		else{hashBig(bases, quals);}
	}
	
	public void hashSmall(byte[] bases, byte[] quals){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		int len=0;
		
		long kmer=0, rkmer=0;
		
		if(minProb>0 && quals!=null){//Debranched loop
			assert(quals.length==bases.length) : quals.length+", "+bases.length;
			float prob=1;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				
				{//Update probability
					byte q=quals[i];
					prob=prob*PROB_CORRECT[q];
					if(len>k){
						byte oldq=quals[i-k];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}
				if(x>=0){
					len++;
				}else{
					len=0;
					kmer=rkmer=0;
					prob=1;
				}
				if(len>=k && prob>=minProb){
					add(Tools.max(kmer, rkmer));
				}
			}
		}else{

			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				
				if(x>=0){
					len++;
				}else{
					len=0;
					kmer=rkmer=0;
				}
				if(len>=k){
					add(Tools.max(kmer, rkmer));
				}
			}
		}
	}
	
	public void hashBig(byte[] bases, byte[] quals){
		
		Kmer kmer=getLocalKmer();
		int len=0;
		float prob=1;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			kmer.addRightNumeric(x);
			if(minProb>0 && quals!=null){//Update probability
				prob=prob*PROB_CORRECT[quals[i]];
				if(len>k){
					byte oldq=quals[i-k];
					prob=prob*PROB_CORRECT_INVERSE[oldq];
				}
			}
			if(AminoAcid.isFullyDefined(b)){
				len++;
			}else{
				len=0;
				prob=1;
			}
			if(len>=k && prob>=minProb){
				add(kmer.xor());
			}
		}
	}
	
	public void add(LogLog_old log){
		if(atomic && maxArrayA!=log.maxArrayA){
			for(int i=0; i<buckets; i++){
				maxArrayA.set(maxArrayA.get(i), log.maxArrayA.get(i));
			}
		}else{
			for(int i=0; i<buckets; i++){
				maxArray[i]=Tools.max(maxArray[i], log.maxArray[i]);
			}
		}
	}
	
	public void hash(final long number){
		if(number%SKIPMOD!=0){return;}
		long key=number;
		
//		int i=(int)(number%5);
//		key=Long.rotateRight(key, 1);
//		key=hash(key, tables[i%numTables]);
		key=hash(key, tables[((int)number)&numTablesMask]);
		int leading=Long.numberOfLeadingZeros(key);
//		counts[leading]++;
		
		if(leading<3){return;}
//		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		final int bucket=(int)(key&bucketMask);
		
		if(atomic){
			int x=maxArrayA.get(bucket);
			while(leading>x){
				boolean b=maxArrayA.compareAndSet(bucket, x, leading);
				if(b){x=leading;}
				else{x=maxArrayA.get(bucket);}
			}
		}else{
			maxArray[bucket]=Tools.max(leading, maxArray[bucket]);
		}
	}
	
	private static long[][] makeCodes(int length, int bits, long seed){
		Random randy=Shared.threadLocalRandom(seed);
		int modes=1<<bits;
		long[][] r=new long[length][modes];
		for(int i=0; i<length; i++){
			for(int j=0; j<modes; j++){
				long x=randy.nextLong();
				while(Long.bitCount(x)>33){
					x&=(~(1L<<randy.nextInt(64)));
				}
				while(Long.bitCount(x)<31){
					x|=(1L<<randy.nextInt(64));
				}
				r[i][j]=x;
				
			}
		}
		return r;
	}
	
	public final int k;
	public final int numTables=4;
	public final int numTablesMask=numTables-1;
	public final int bits;
	public final float minProb;
//	public final int hashes;
	public final int steps;
	private final long[][][] tables;
	public final AtomicIntegerArray maxArrayA;
	public final int[] maxArray;
	public final int buckets;
	public final int bucketMask;
	private final ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	
	protected Kmer getLocalKmer(){
		Kmer kmer=localKmer.get();
		if(kmer==null){
			localKmer.set(new Kmer(k));
			kmer=localKmer.get();
		}
		kmer.clearFast();
		return kmer;
	}
	
	private static class LogLogWrapper{
		
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
				
				if(parser.parse(arg, a, b)){
					//do nothing
				}else if(a.equals("buckets") || a.equals("loglogbuckets")){
					long x=Parse.parseKMG(b);
					buckets=(int)Primes.primeAtLeast(Tools.min(1000000, x));
				}else if(a.equals("bits") || a.equals("loglogbits")){
					bits=Integer.parseInt(b);
				}else if(a.equals("k") || a.equals("loglogk")){
					k=Integer.parseInt(b);
				}else if(a.equals("seed") || a.equals("loglogseed")){
					seed=Long.parseLong(b);
				}else if(a.equals("minprob") || a.equals("loglogminprob")){
					minProb=Float.parseFloat(b);
				}else if(a.equals("verbose")){
					verbose=Parse.parseBoolean(b);
				}else if(a.equals("atomic")){
					assert(false) : "Atomic flag disabled.";
//					atomic=Parse.parseBoolean(b);
				}else if(a.equals("parse_flag_goes_here")){
					//Set a variable here
				}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
					parser.in1=b;
				}else{
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
			
			assert(in1!=null && in1.length>0) : "No primary input file specified.";
			{
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

			assert(FastaReadInputStream.settingsOK());
		}
		
		
		void process(){
			Timer t=new Timer();
			
			LogLog_old log=new LogLog_old(buckets, bits, k, seed, minProb);
			
			for(int ffnum=0; ffnum<ffin1.length; ffnum++){
				ConcurrentReadInputStream cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, ffin1[ffnum], ffin2[ffnum]);
				cris.start();

				LogLogThread[] threads=new LogLogThread[Shared.threads()];
				for(int i=0; i<threads.length; i++){
					threads[i]=new LogLogThread((atomic ? log : new LogLog_old(buckets, bits, k, seed, minProb)), cris);
				}
				for(LogLogThread llt : threads){
					llt.start();
				}
				for(LogLogThread llt : threads){
					while(llt.getState()!=Thread.State.TERMINATED){
						try {
							llt.join();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					if(!atomic){log.add(llt.log);}
				}

				errorState|=ReadWrite.closeStreams(cris);
			}
			
			final int[] max=new int[buckets];
			if(atomic){
				for(int i=0; i<log.maxArrayA.length(); i++){
					//				System.err.println(log.maxArray.get(i));
					max[i]=log.maxArrayA.get(i);
				}
			}
			
			t.stop();
			
			
			long cardinality=log.cardinality();
			
			if(out!=null){
				ReadWrite.writeString(cardinality+"\n", out);
			}
			
//			Arrays.sort(copy);
//			System.err.println("Median:        "+copy[Tools.median(copy)]);
			
//			System.err.println("Mean:          "+Tools.mean(copy));
//			System.err.println("Harmonic Mean: "+Tools.harmonicMean(copy));
			System.err.println("Cardinality:   "+log.cardinality());
//			System.err.println("CardinalityH:  "+log.cardinalityH());
			
//			for(long i : log.counts){System.err.println(i);}
			
			System.err.println("Time: \t"+t);
		}
		
		private class LogLogThread extends Thread{
			
			LogLogThread(LogLog_old log_, ConcurrentReadInputStream cris_){
				log=log_;
				cris=cris_;
			}
			
			@Override
			public void run(){
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					
					for(Read r : reads){
//						if(!r.validated()){r.validate(true);}
//						if(r.mate!=null && !r.mate.validated()){r.mate.validate(true);}
						log.hash(r);
					}
					
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				cris.returnList(ln);
			}
			
			private final LogLog_old log;
			private final ConcurrentReadInputStream cris;
			
		}
		
		/*--------------------------------------------------------------*/
		/*----------------            Fields            ----------------*/
		/*--------------------------------------------------------------*/
		
		private int buckets=2048;//1999
		private int bits=8;
		private int k=31;
		private long seed=-1;
		private float minProb=0;
		
		
		private String[] in1=null;
		private String[] in2=null;
		private String out=null;
		
		/*--------------------------------------------------------------*/
		
		protected long readsProcessed=0;
		protected long basesProcessed=0;
		
		private long maxReads=-1;
		
		boolean overwrite=true;
		boolean append=false;
		boolean errorState=false;
		
		/*--------------------------------------------------------------*/
		/*----------------         Final Fields         ----------------*/
		/*--------------------------------------------------------------*/
		
		private final FileFormat[] ffin1;
		private final FileFormat[] ffin2;
		
		/*--------------------------------------------------------------*/
		/*----------------        Common Fields         ----------------*/
		/*--------------------------------------------------------------*/
	}
	
	public static final float[] PROB_CORRECT=Arrays.copyOf(align2.QualityTools.PROB_CORRECT, 128);
	public static final float[] PROB_CORRECT_INVERSE=Arrays.copyOf(align2.QualityTools.PROB_CORRECT_INVERSE, 128);
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static final boolean atomic=true;
	private static final long SKIPMOD=3;
	public static long lastCardinality=-1;
	
}
