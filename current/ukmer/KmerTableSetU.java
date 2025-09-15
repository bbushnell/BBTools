package ukmer;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicLong;

import assemble.Contig;
import bloom.KmerCountAbstract;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import kmer.AbstractKmerTableSet;
import kmer.DumpThread;
import kmer.ScheduleMaker;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import tracker.ReadStats;


/**
 * Loads and holds kmers for Tadpole2/KmerCountExact
 * @author Brian Bushnell
 * @date Jun 22, 2015
 *
 */
public class KmerTableSetU extends AbstractKmerTableSet {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		KmerTableSetU set=new KmerTableSetU(args);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		set.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	private KmerTableSetU(String[] args){
		this(args, 0);//+5 if using ownership and building contigs
	}
	
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerTableSetU(String[] args, int extraBytesPerKmer_){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}//TODO - no easy way to close outstream
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		boolean prealloc_=false;
		int kbig_=62;
		int ways_=-1;
		int filterMax_=2;
		boolean ecco_=false, merge_=false;
		boolean rcomp_=true;
		double minProb_=defaultMinprob;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
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
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in1.add(ss);
					}
				}
			}else if(a.equals("in2")){
				in2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in2.add(ss);
					}
				}
			}else if(a.equals("extra")){
				extra.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						extra.add(ss);
					}
				}
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=Parse.parseIntKMG(b);
			}else if(a.equals("showstats") || a.equals("stats")){
				showStats=Parse.parseBoolean(b);
			}else if(a.equals("ways")){
				ways_=Parse.parseIntKMG(b);
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
				buflen=Parse.parseIntKMG(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nk needs an integer value such as k=50.  Default is 62.\n";
				kbig_=Parse.parseIntKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco_=Parse.parseBoolean(b);
			}else if(a.equals("merge")){
				merge_=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
//				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbose2")){
//				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				verbose2=Parse.parseBoolean(b);
			}else if(a.equals("minprob")){
				minProb_=Double.parseDouble(b);
			}else if(a.equals("minprobprefilter") || a.equals("mpp")){
				minProbPrefilter=Parse.parseBoolean(b);
			}else if(a.equals("minprobmain") || a.equals("mpm")){
				minProbMain=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
					prealloc_=Parse.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc_=(preallocFraction>0);
				}
			}else if(a.equals("prefilter")){
				if(b==null || b.length()<1 || !Tools.isDigit(b.charAt(0))){
					prefilter=Parse.parseBoolean(b);
				}else{
					filterMax_=Parse.parseIntKMG(b);
					prefilter=filterMax_>0;
				}
			}else if(a.equals("prefiltersize") || a.equals("prefilterfraction") || a.equals("pff")){
				prefilterFraction=Tools.max(0, Double.parseDouble(b));
				assert(prefilterFraction<=1) : "prefiltersize must be 0-1, a fraction of total memory.";
				prefilter=prefilterFraction>0;
			}else if(a.equals("prehashes") || a.equals("hashes")){
				prehashes=Parse.parseIntKMG(b);
			}else if(a.equals("prefilterpasses") || a.equals("prepasses")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.equalsIgnoreCase("auto")){
					prepasses=-1;
				}else{
					prepasses=Parse.parseIntKMG(b);
				}
			}else if(a.equals("onepass")){
				onePass=Parse.parseBoolean(b);
			}else if(a.equals("passes")){
				int passes=Parse.parseIntKMG(b);
				onePass=(passes<2);
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("filterMemoryOverride") || a.equalsIgnoreCase("filterMemory") || 
					a.equalsIgnoreCase("prefilterMemory") || a.equalsIgnoreCase("filtermem")){
				filterMemoryOverride=Parse.parseKMG(b);
			}
			
			else if(IGNORE_UNKNOWN_ARGS){
				//Do nothing
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();
			
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
		}
		
		if(prepasses==0 || !prefilter){
			prepasses=0;
			prefilter=false;
		}
		
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
//			long tmemory=Runtime.getRuntime().totalMemory();
			usableMemory=(long)Tools.max(((memory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
			if(prepasses==0 || !prefilter){
				filterMemory0=filterMemory1=0;
			}else if(filterMemoryOverride>0){
				filterMemory0=filterMemory1=filterMemoryOverride;
			}else{
				double low=Tools.min(prefilterFraction, 1-prefilterFraction);
				double high=1-low;
				if(prepasses<0 || (prepasses&1)==1){//odd passes
					filterMemory0=(long)(usableMemory*low);
					filterMemory1=(long)(usableMemory*high);
				}else{//even passes
					filterMemory0=(long)(usableMemory*high);
					filterMemory1=(long)(usableMemory*low);
				}
			}
			tableMemory=(long)(usableMemory*.95-filterMemory0);
		}
		
		mult=Kmer.getMult(kbig_);
		k=Kmer.getK(kbig_);
		kbig=k*mult;
		kbig2=kbig-1;
		assert(k<=31);

		prealloc=prealloc_;
		bytesPerKmer=4+8*mult+extraBytesPerKmer_;
		assert(bytesPerKmer>=4+8*mult) : bytesPerKmer+", "+mult+", "+k+", "+kbig+", "+(4+8*mult)+", "+extraBytesPerKmer_;
		if(ways_<1){
			long maxKmers=(2*tableMemory)/bytesPerKmer;
			long minWays=Tools.min(10000, maxKmers/Integer.MAX_VALUE);
			ways_=(int)Tools.max(31, (int)(Tools.min(96, Shared.threads())*2.5), minWays);
			ways_=(int)Primes.primeAtLeast(ways_);
			assert(ways_>0);
//			System.err.println("ways="+ways_);
		}
//		assert(false) : extraBytesPerKmer_+bytesPerKmer+", "+mult+", "+k+", "+kbig+", "+(4+8*mult)+", "+ways_;
		
		/* Set final variables; post-process and validate argument combinations */
		
		onePass=onePass&prefilter;
		ways=ways_;
		filterMax=Tools.min(filterMax_, 0x7FFFFFFF);
		ecco=ecco_;
		merge=merge_;
		minProb=(float)minProb_;
		rcomp=rcomp_;
		Kmer.rcomp=rcomp;
//		System.err.println("rcomp="+rcomp);
		estimatedKmerCapacity=(long)((tableMemory*1.0/bytesPerKmer)*((prealloc ? preallocFraction : 0.81))*(HashArrayU.maxLoadFactorFinal*0.97));
		KmerCountAbstract.minProb=(minProbPrefilter ? minProb : 0);
		
		if(kbig!=kbig_){
			System.err.println("K was changed from "+kbig_+" to "+kbig);
		}
		
		if(k<1){throw new RuntimeException("\nk needs an integer value above 0, such as k=27.  Default is 62.\n");}
		
		if(initialSize<1){
			final long memOverWays=tableMemory/(bytesPerKmer*ways);
			final double mem2=(prealloc ? preallocFraction : 1)*tableMemory;
			initialSize=(prealloc || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(bytesPerKmer*ways))) : initialSizeDefault);
			if(initialSize!=initialSizeDefault){
				System.err.println("Initial size set to "+initialSize);
			}
		}
		
		/* Adjust I/O settings and filenames */
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
		
		for(int i=0; i<in1.size(); i++){
			String s=in1.get(i);
			if(s!=null && s.contains("#") && !new File(s).exists()){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				in1.set(i, a+1+b);
				in2.add(a+2+b);
			}
		}
		
		if(!extra.isEmpty()){
			ArrayList<String> temp=(ArrayList<String>) extra.clone();
			extra.clear();
			for(String s : temp){
				if(s!=null && s.contains("#") && !new File(s).exists()){
					int pound=s.lastIndexOf('#');
					String a=s.substring(0, pound);
					String b=s.substring(pound+1);
					extra.add(a);
					extra.add(b);
				}else{
					extra.add(s);
				}
			}
		}
		
		{
			boolean allowDuplicates=true;
			if(!Tools.testInputFiles(allowDuplicates, true, in1, in2, extra)){
				throw new RuntimeException("\nCan't read some input files.\n");  
			}
		}
		assert(THREADS>0);
		
		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			outstream.println("Ways="+ways+", initialSize="+initialSize+", prefilter="+(prefilter ? "t" : "f")+", prealloc="+(prealloc ? (""+preallocFraction) : "f"));
			Shared.printMemory();
			outstream.println();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void clear(){
		tables=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	protected void allocateTables(){
		assert(tables==null);
		tables=null;
		final int tableType=AbstractKmerTableU.ARRAY1D;
		
		ScheduleMaker scheduleMaker=new ScheduleMaker(ways, bytesPerKmer, prealloc, 
				(prealloc ? preallocFraction : 1.0), -1, (prefilter ? prepasses : 0), prefilterFraction, filterMemoryOverride);
		int[] schedule=scheduleMaker.makeSchedule();
		tables=AbstractKmerTableU.preallocate(ways, tableType, schedule, k, kbig);
	}
	
	/**
	 * Load reads into tables, using multiple LoadThread.
	 */
	@Override
	public long loadKmers(String fname1, String fname2){
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff1, ff2);
			cris.start(); //4567
		}
		
		/* Create ProcessThreads */
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(THREADS);
		for(int i=0; i<THREADS; i++){alpt.add(new LoadThread(cris));}
		for(LoadThread pt : alpt){pt.start();}
		
		long added=0;
		
		/* Wait for threads to die, and gather statistics */
		for(LoadThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			added+=pt.added;
			
			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
			readsTrimmed+=pt.readsTrimmedT;
			basesTrimmed+=pt.basesTrimmedT;
			kmersIn+=pt.kmersInT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads kmers.
	 */
	private class LoadThread extends Thread{
		
		/**
		 * Constructor
		 * @param cris_ Read input stream
		 */
		public LoadThread(ConcurrentReadInputStream cris_){
			cris=cris_;
			table=new HashBufferU(tables, buflen, kbig, false);
			kmer=new Kmer(k, mult);
		}
		
		@Override
		public void run(){
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					Read r1=reads.get(i);
					Read r2=r1.mate;
					
					if(!r1.validated()){r1.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}
					
					if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}
					
					readsInT++;
					basesInT+=r1.length();
					if(r2!=null){
						readsInT++;
						basesInT+=r2.length();
					}
					
					//Determine whether to discard the reads based on average quality
					if(minAvgQuality>0){
						if(r1!=null && r1.quality!=null && r1.avgQuality(false, minAvgQualityBases)<minAvgQuality){r1.setDiscarded(true);}
						if(r2!=null && r2.quality!=null && r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){r2.setDiscarded(true);}
					}
					
					if(r1!=null){
						if(qtrimLeft || qtrimRight){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r1.length()<kbig){r1.setDiscarded(true);}
					}
					if(r2!=null){
						if(qtrimLeft || qtrimRight){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1);
							basesTrimmedT+=x;
							readsTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2.length()<kbig){r2.setDiscarded(true);}
					}
					
					if((ecco || merge) && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){
						if(merge){
							final int insert=BBMerge.findOverlapStrict(r1, r2, false);
							if(insert>0){
								r2.reverseComplement();
								r1=r1.joinRead(insert);
								r2=null;
							}
						}else if(ecco){
							BBMerge.findOverlapStrict(r1, r2, true);
						}
					}

					if(r1!=null){
						if(r1.discarded()){
							lowqBasesT+=r1.length();
							lowqReadsT++;
						}else{
							long temp=addKmersToTable(r1, kmer);
							added+=temp;
							if(verbose){System.err.println("A: Added "+temp);}
						}
					}
					if(r2!=null){
						if(r2.discarded()){
							lowqBasesT+=r2.length();
							lowqReadsT++;
						}else{
							long temp=addKmersToTable(r2, kmer);
							added+=temp;
							if(verbose){System.err.println("B: Added "+temp);}
						}
					}
				}
				
				//Fetch a new read list
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			long temp=table.flush();
			if(verbose){System.err.println("Flush: Added "+temp);}
			added+=temp;
		}
		
		
		private final int addKmersToTable(final Read r, Kmer kmer){
			if(onePass){return addKmersToTable_onePass(r, kmer);}
			if(r==null || r.bases==null){return 0;}
			final float minProb2=(minProbMain ? minProb : 0);
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			int created=0;
			int len=0;

			if(bases==null || bases.length<kbig){return -1;}
			kmer.clear();
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
//				assert(x>=0) : ((char)b)+", "+x+", "+new String(bases);

				//Update kmers
//				kmer.addRight(b);
				kmer.addRightNumeric(x);

				if(minProb2>0 && quals!=null){//Update probability
					prob=prob*PROB_CORRECT[quals[i]];
					if(len>kbig){
						byte oldq=quals[i-kbig];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					prob=1;
				}else{len++;}

				assert(len==kmer.len);
				
//				if(verbose){System.err.println("A: Scanning i="+i+", len="+len+", kmer="+kmer.toString()+"\t"+new String(bases, Tools.max(0, i-kbig2), Tools.min(i+1, kbig)));}
				if(len>=kbig && prob>=minProb2){
					kmersInT++;
//					System.err.println("kmer="+kmer+"; xor()="+kmer.xor()+"; filterMax2="+filterMax2+"; prefilter="+prefilter);
//					System.err.println("prefilterArray.read(xor.key())="+prefilterArray.read(kmer.xor())+"");
//					System.err.println("prefilterArray.read(kmer.key())="+prefilterArray.read(kmer.key())+"");
					if(!prefilter || prefilterArray.read(kmer.xor())>filterMax2){
						int temp=table.incrementAndReturnNumCreated(kmer);
//						System.err.println("kmer="+kmer+"; xor()="+kmer.xor()+"; temp="+temp+" ");
						created+=temp;
						if(verbose){System.err.println("C: Added "+temp);}
					}
				}
			}
			
			return created;
		}
		
		
		private final int addKmersToTable_onePass(final Read r, Kmer kmer){
			assert(prefilter);
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			int created=0;
			int len=0;

			if(bases==null || bases.length<kbig){return -1;}
			kmer.clear();
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];

				//Update kmers
				kmer.addRight(b);

				if(minProb>0 && quals!=null){//Update probability
					prob=prob*PROB_CORRECT[quals[i]];
					if(len>kbig){
						byte oldq=quals[i-kbig];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					prob=1;
				}else{len++;}

				assert(len==kmer.len);
				
				if(verbose){System.err.println("B: Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-kbig2), Tools.min(i+1, kbig)));}
				if(len>=kbig && prob>=minProb){
					final long xor=kmer.xor();
					int count=prefilterArray.incrementAndReturnUnincremented(xor, 1);
					if(count>=filterMax2){
						int temp=table.incrementAndReturnNumCreated(kmer);
						created+=temp;
						if(verbose){System.err.println("D: Added "+temp);}
					}
				}
			}
			return created;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Input read stream */
		private final ConcurrentReadInputStream cris;
		
		private final HashBufferU table;
		
		public long added=0;
		
		private long readsInT=0;
		private long basesInT=0;
		private long lowqReadsT=0;
		private long lowqBasesT=0;
		private long readsTrimmedT=0;
		private long basesTrimmedT=0;
		private long kmersInT=0;
		private final Kmer kmer;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------          Convenience         ----------------*/
	/*--------------------------------------------------------------*/
	
	public int regenerateCounts(byte[] bases, IntList counts, final int ca, final Kmer kmer){
		final int b=ca+kbig; //first base changed
		final int lim=Tools.min(counts.size, ca+kbig+1); //count limit
//		System.err.println("ca="+ca+", b="+b+", lim="+lim);
//		System.err.println("Regen from count "+(ca+1)+"-"+lim);
		int len=0;
		int valid=0;
		kmer.clear();
//		System.err.println("ca="+ca+", b="+b+", lim="+lim+", "+counts);
		
		//Generate initial kmer
		for(int i=b-kbig; i<b; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			
			kmer.addRight(base);
			
			if(x<0){
				len=0;
			}else{len++;}
			assert(len==kmer.len);
		}
		assert(len==kbig || Tools.indexOf(bases, (byte)'N')>=ca) : new String(bases)+"\n"+ca+", "+len;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts.
		 * i is an index in the base array, j is an index in the count array. */
		for(int i=b, j=ca+1; j<lim; i++, j++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			
			kmer.addRight(base);
			
			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);
			
			if(len>=kbig){
				valid++;
				int count=getCount(kmer);
				counts.set(j, count);
			}else{
				counts.set(j, 0);
			}
		}
//		System.err.println("ca="+ca+", b="+b+", lim="+lim+", "+counts);
		return valid;
	}
	
	@Override
	public int regenerateCounts(byte[] bases, IntList counts, final Kmer kmer, BitSet changed){
		assert(!changed.isEmpty());
		final int firstBase=changed.nextSetBit(0), lastBase=changed.length()-1;
		final int ca=firstBase-kbig;
//		final int b=changed.nextSetBit(0);ca+kbig; //first base changed
		final int firstCount=Tools.max(firstBase-kbig+1, 0), lastCount=Tools.min(counts.size-1, lastBase); //count limit
//		System.err.println("ca="+ca+", b="+b+", lim="+lim);
//		System.err.println("Regen from count "+(ca+1)+"-"+lim);
		int len=0;
		int valid=0;
		kmer.clear();
//		System.err.println("ca="+ca+", b="+b+", lim="+lim+", "+counts);
		
		//Generate initial kmer
		for(int i=Tools.max(0, firstBase-kbig+1); i<firstBase; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			
			kmer.addRight(base);
			
			if(x<0){
				len=0;
			}else{len++;}
			assert(len==kmer.len);
		}
		assert(len==kbig-1 || Tools.indexOf(bases, (byte)'N')>=0 || firstBase<kbig) :
			new String(bases)+"\n"+ca+", "+len+", "+firstBase;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts.
		 * i is an index in the base array, j is an index in the count array. */
		for(int i=firstBase, lim=Tools.min(lastBase+k-1, bases.length-1); i<=lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			
			kmer.addRight(base);
			
			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);
			
			final int c=i-kbig+1;
			if(i>=firstBase){
				if(len>=kbig){
					valid++;
					int count=getCount(kmer);
					counts.set(c, count);
				}else if(c>=0){
					counts.set(c, 0);
				}
			}
		}
//		System.err.println("ca="+ca+", b="+b+", lim="+lim+", "+counts);
		return valid;
	}
	
	@Override
	public int fillSpecificCounts(byte[] bases, IntList counts, BitSet positions, final Kmer kmer){
		counts.clear();
		
		{
			Kmer x=leftmostKmer(bases, bases.length, kmer);
			assert((x!=null)==(kmer.len==kbig));
		}
		int len=kmer.len;
		int valid=0;
		if(len>=kbig){
			valid++;
			int count=getCount(kmer);
			counts.add(count);
		}else{
			counts.add(0);
		}
		assert(kmer.len==len);
		assert(len<=kbig) : len+", "+kbig;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=kbig, j=0; i<bases.length; i++, j++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			
			kmer.addRight(base);
			
			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);
			
			if(len>=kbig && (positions==null || positions.get(j))){
				valid++;
				int count=getCount(kmer);
				counts.add(Tools.max(count, 0));
			}else{
				counts.add(0);
			}
		}
		return valid;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public long regenerate(final int limit){
		long sum=0;
		for(AbstractKmerTableU akt : tables){
			sum+=akt.regenerate(limit);
		}
		return sum;
	}

	public HashArrayU1D getTable(Kmer kmer){
		return (HashArrayU1D) tables[kmer.mod(ways)];
	}
	
	@Override
	public HashArrayU1D getTable(int tnum){
		return (HashArrayU1D) tables[tnum];
	}
	
	@Override
	public long[] fillHistogram(int histMax) {
		return HistogramMakerU.fillHistogram(tables, histMax);
	}
	
	@Override
	public void countGC(long[] gcCounts, int max) {
		for(AbstractKmerTableU set : tables){
			set.countGC(gcCounts, max);
		}
	}
	
	@Override
	public void initializeOwnership(){
		OwnershipThread.initialize(tables);
	}
	
	@Override
	public void clearOwnership(){
		OwnershipThread.clear(tables);
	}
	
	public Kmer rightmostKmer(final ByteBuilder bb, Kmer kmer){
		return rightmostKmer(bb.array, bb.length(), kmer);
	}
	
	public Kmer rightmostKmer(final byte[] bases, final int blen, final Kmer kmer){
		kmer.clear();
		if(blen<kbig){return null;}
		int len=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			for(int i=blen-kbig; i<blen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				kmer.addRight(b);
				
				if(x<0){len=0;}
				else{len++;}
				assert(len==kmer.len);
				
				//if(verbose){outstream.println("C: Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-kbig2), Tools.min(i+1, kbig)));}
			}
		}
		
		if(len<kbig){return null;}
		else{assert(len==kbig);}
		return kmer;
	}
	
	public Kmer leftmostKmer(final ByteBuilder bb, final Kmer kmer){
		return leftmostKmer(bb.array, bb.length(), kmer);
	}
	
	public Kmer leftmostKmer(final byte[] bases, final int blen, final Kmer kmer){
		kmer.clear();
		if(blen<kbig){return null;}
		int len=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			for(int i=0; i<kbig; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				kmer.addRight(b);
				
				if(x<0){len=0;}
				else{len++;}
				assert(len==kmer.len);
				
				if(verbose){outstream.println("D: Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-kbig2), Tools.min(i+1, kbig)));}
			}
		}
		
		if(len<kbig){return null;}
		else{assert(len==kbig);}
		return kmer;
	}
	
	public boolean doubleClaim(final ByteBuilder bb, final int id, Kmer kmer){
		return doubleClaim(bb.array, bb.length(), id, kmer);
	}
	
	/** Ensures there can be only one owner. */
	public boolean doubleClaim(final byte[] bases, final int blength, final int id, Kmer kmer){
		boolean success=claim(bases, blength, id, true, kmer);
		if(verbose){outstream.println("success1="+success+", id="+id+", blength="+blength);}
		if(!success){return false;}
		success=claim(bases, blength, id+CLAIM_OFFSET, true, kmer);
		if(verbose){outstream.println("success2="+success+", id="+id+", blength="+blength);}
		return success;
	}
	
	public boolean claim(final ByteBuilder bb, final int id, final boolean exitEarly, Kmer kmer){
		return claim(bb.array, bb.length(), id, exitEarly, kmer);
	}
	
	public float calcCoverage(final byte[] bases, final int blen, final Kmer kmer){
		if(blen<kbig){return 0;}
		int len=0;
		kmer.clear();
		long sum=0;
		int kmers=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		for(int i=0; i<blen; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);

			if(len>=kbig){
				int count=getCount(kmer);
				sum+=count;
				kmers++;
			}
		}
		return sum==0 ? 0 : sum/(float)kmers;
	}
	
	public float calcCoverage(final Contig contig, final Kmer kmer){
		final byte[] bases=contig.bases;
		if(bases.length<kbig){return 0;}
		int len=0;
		kmer.clear();
		long sum=0;
		int max=0, min=Integer.MAX_VALUE;
		int kmers=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);

			if(len>=kbig){
				int count=getCount(kmer);
				sum+=count;
				max=Tools.max(count, max);
				min=Tools.min(count, min);
				kmers++;
			}
		}
		contig.coverage=sum==0 ? 0 : sum/(float)kmers;
		contig.maxCov=max;
		contig.minCov=sum==0 ? 0 : min;
		return contig.coverage;
	}
	
	public boolean claim(final byte[] bases, final int blen, final int id, boolean exitEarly, final Kmer kmer){
		if(blen<kbig){return false;}
		if(verbose){outstream.println("Thread "+id+" claim start.");}
		int len=0;
		kmer.clear();
		boolean success=true;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blen && success; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);

			if(len>=kbig){
				success=claim(kmer, id/*, rid, i*/);
				success=(success || !exitEarly);
			}
		}
		return success;
	}
	
	public boolean claim(Kmer kmer, final int id/*, final long rid, final int pos*/){
		//TODO: rid and pos are just for debugging.
		final int way=kmer.mod(ways);
		final AbstractKmerTableU table=tables[way];
		final int count=table.getValue(kmer);
		assert(count==-1 || count>0) : count;
//		if(verbose  /*|| true*/){outstream.println("Count="+count+".");}
		if(count<0){return true;}
		assert(count>0) : count;
		final int owner=table.setOwner(kmer, id);
		if(verbose){outstream.println("owner="+owner+".");}
//		assert(owner==id) : id+", "+owner+", "+rid+", "+pos;
		return owner==id;
	}
	
	public void release(ByteBuilder bb, final int id, final Kmer kmer){
		release(bb.array, bb.length(), id, kmer);
	}
	
	public void release(final byte[] bases, final int blen, final int id, final Kmer kmer){
		if(verbose  /*|| true*/){outstream.println("*Thread "+id+" release start.");}
		int len=0;
		kmer.clear();
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blen; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);

			if(len>=kbig){
				release(kmer, id);
			}
		}
	}
	
	public boolean release(Kmer kmer, final int id){
		final int way=kmer.mod(ways);
		final AbstractKmerTableU table=tables[way];
		final int count=table.getValue(kmer);
//		if(verbose  /*|| true*/){outstream.println("Count="+count+".");}
		if(count<1){return true;}
		return table.clearOwner(kmer, id);
	}
	
	public int findOwner(ByteBuilder bb, final int id, final Kmer kmer){
		return findOwner(bb.array, bb.length(), id, kmer);
	}
	
	public int findOwner(final byte[] bases, final int blen, final int id, final Kmer kmer){
		int len=0;
		kmer.clear();
		int maxOwner=-1;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blen; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);
			//if(verbose){System.err.println("E: Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-kbig2), Tools.min(i+1, kbig)));}
			if(len>=kbig){
				int owner=findOwner(kmer);
				maxOwner=Tools.max(owner, maxOwner);
				if(maxOwner>id){break;}
			}
		}
		return maxOwner;
	}
	
	public int findOwner(final Kmer kmer){
		final int way=kmer.mod(ways);
		final AbstractKmerTableU table=tables[way];
		final int count=table.getValue(kmer);
		if(count<0){return -1;}
		final int owner=table.getOwner(kmer);
		return owner;
	}

	public int getCount(Kmer kmer){
		int way=kmer.mod(ways);
		return tables[way].getValue(kmer);
	}

	public int getOwner(Kmer kmer){
		int way=kmer.mod(ways);
		return tables[way].getOwner(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Fill Counts         ----------------*/
	/*--------------------------------------------------------------*/
	
	public int fillRightCounts(Kmer kmer, int[] counts){
		if(FAST_FILL && MASK_CORE && k>2){
			return fillRightCounts_fast(kmer, counts);
		}else{
			return fillRightCounts_safe(kmer, counts);
		}
	}
	
	public int fillLeftCounts(Kmer kmer, int[] counts){
		if(FAST_FILL && MASK_CORE && k>2){
			return fillLeftCounts_fast(kmer, counts);
		}else{
			return fillLeftCounts_safe(kmer, counts);
		}
	}
	
	public int fillRightCounts_fast(Kmer kmer, int[] counts){
		assert(MASK_CORE);
		final long old=kmer.addRightNumeric(0);
		if(kmer.corePalindrome()){
			kmer.addLeftNumeric(old);
			return fillRightCounts_safe(kmer, counts);
		}
		final int way=kmerToWay(kmer);
		final HashArrayU table=(HashArrayU) tables[way];
		final int cell=table.kmerToCell(kmer);
		int max=-1, maxPos=0;
		
		{//Iteration 0
			final int count=Tools.max(0, table.getValue(kmer, cell));
			assert(count==NOT_PRESENT || count>=0);
			counts[0]=count;
			max=count;
			maxPos=0;
			kmer.addLeftNumeric(old);
		}

		for(int i=1; i<=3; i++){
			kmer.addRightNumeric(i);
			if(verbose){outstream.println("kmer:               "+kmer);}
			final int count=Tools.max(0, table.getValue(kmer, cell));
			assert(count==NOT_PRESENT || count>=0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
			kmer.addLeftNumeric(old);
		}
		return maxPos;
	}
	
	public int fillRightCounts_safe(Kmer kmer, int[] counts){
		assert(kmer.len>=kbig);
		if(verbose){outstream.println("fillRightCounts:   "+kmer);}
		int max=-1, maxPos=0;
		
//		final Kmer kmer2=new Kmer(kmer);//123 TODO: Slow, for an assertion only;
		
		for(int i=0; i<=3; i++){
			final long old=kmer.addRightNumeric(i);
			if(verbose){outstream.println("kmer:               "+kmer);}
			int way=kmer.mod(ways);
			int count=tables[way].getValue(kmer);
			assert(count==NOT_PRESENT || count>=0);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
			kmer.addLeftNumeric(old);
		}
		return maxPos;
	}
	
//	public int fillLeftCounts_fast(Kmer kmer, int[] counts){
//		assert(MASK_CORE);
//		final long old=kmer.addLeftNumeric(0);
//		if(kmer.corePalindrome()){
//			kmer.addRightNumeric(old);
//			return fillLeftCounts_safe(kmer, counts);
//		}
//		final int way=kmerToWay(kmer);
//		final HashArrayU table=(HashArrayU) tables[way];
//		final int cell=table.kmerToCell(kmer);
//		int max=-1, maxPos=0;
//		
//		for(int i=0; i<=3; i++){
//			kmer.addLeftNumeric(i);
//			if(verbose){outstream.println("kmer:               "+kmer);}
//			final int count=Tools.max(0, table.getValue(kmer, cell));
//			assert(count==NOT_PRESENT || count>=0);
//			counts[i]=count;
//			if(count>max){
//				max=count;
//				maxPos=i;
//			}
//			kmer.addRightNumeric(old);
//		}
//		return maxPos;
//	}
	
	public int fillLeftCounts_fast(Kmer kmer, int[] counts){
		assert(MASK_CORE);
		final long old=kmer.addLeftNumeric(0);
		if(kmer.corePalindrome()){
			kmer.addRightNumeric(old);
			return fillLeftCounts_safe(kmer, counts);
		}
		final int way=kmerToWay(kmer);
		final HashArrayU table=(HashArrayU) tables[way];
		final int cell=table.kmerToCell(kmer);
		int max=-1, maxPos=0;
		
		{//Iteration 0
			if(verbose){outstream.println("kmer:               "+kmer);}
			final int count=Tools.max(0, table.getValue(kmer, cell));
			assert(count==NOT_PRESENT || count>=0);
			counts[0]=count;
			max=count;
			maxPos=0;
			kmer.addRightNumeric(old);
		}
		
		for(int i=1; i<=3; i++){
			kmer.addLeftNumeric(i);
			if(verbose){outstream.println("kmer:               "+kmer);}
			final int count=Tools.max(0, table.getValue(kmer, cell));
			assert(count==NOT_PRESENT || count>=0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
			kmer.addRightNumeric(old);
		}
		return maxPos;
	}
	
	
	public int fillLeftCounts_safe(final Kmer kmer, int[] counts){
		assert(kmer.len>=kbig);
		if(verbose){outstream.println("fillLeftCounts:    "+kmer);}
		int max=-1, maxPos=0;
		
//		final Kmer kmer2=new Kmer(kmer);//123 TODO: Slow, for an assertion only;
//		assert(false) : kmer+", "+kmer2;
		
		for(int i=0; i<=3; i++){
			if(verbose){
				outstream.println("kmer:               "+kmer+" (key==array1 ? "+(kmer.key()==kmer.array1()));
//				outstream.println("kmer2:              "+kmer2);
			}
			final long old=kmer.addLeftNumeric(i);
			if(verbose){
				outstream.println("after:              "+kmer+" (key==array1 ? "+(kmer.key()==kmer.array1()));
				outstream.println("i="+i+", old="+old);
			}
			int way=kmer.mod(ways);
			int count=tables[way].getValue(kmer);
			assert(count==NOT_PRESENT || count>=0);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
			kmer.addRightNumeric(old);
			if(verbose){outstream.println("restored:           "+kmer);}
//			assert(kmer.equals(kmer2)) : kmer+", "+kmer2+", "+kmer.xor()+", "+kmer2.xor();
		}
		return maxPos;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Printing Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean dumpKmersAsBytes(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining){
		if(fname==null){return false;}
		Timer t=new Timer();
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, true);
		bsw.start();
		for(AbstractKmerTableU set : tables){
			set.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		bsw.poisonAndWait();
		
		t.stop();
		if(printTime){outstream.println("Kmer Dump Time:             \t"+t);}
		return bsw.errorState;
	}
	
	@Override
	public boolean dumpKmersAsBytes_MT(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining){
		
		final int threads=Tools.min(Shared.threads(), tables.length);
		if(threads<3 || DumpThread.NUM_THREADS==1){return dumpKmersAsBytes(fname, mincount, maxcount, printTime, remaining);}
		
		if(fname==null){return false;}
		Timer t=new Timer();
		
		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, true);
		bsw.start();
		DumpThreadU.dump(k, mincount, maxcount, tables, bsw, remaining);
		bsw.poisonAndWait();
		
		t.stop();
		if(printTime){outstream.println("Kmer Dump Time:             \t"+t);}
		return bsw.errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private final StringBuilder toText(long[] kmer){return AbstractKmerTableU.toText(kmer, k);}
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int kbig(){return kbig;}
	@Override
	public long filterMemory(int pass){return ((pass&1)==0) ? filterMemory0 : filterMemory1;}
	@Override
	public boolean ecco(){return ecco;}
	@Override
	public boolean qtrimLeft(){return qtrimLeft;}
	@Override
	public boolean qtrimRight(){return qtrimRight;}
	@Override
	public float minAvgQuality(){return minAvgQuality;}
	@Override
	public long tableMemory(){return tableMemory;}
	@Override
	public long estimatedKmerCapacity(){return estimatedKmerCapacity;}
	@Override
	public int ways(){return ways;}
	@Override
	public boolean rcomp(){return rcomp;}
	
	public final int kmerToWay(final Kmer kmer){
		final int way=(int)(kmer.xor()%ways);
		return way;
	}
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in tables[Y] */
	private AbstractKmerTableU[] tables;
	public AbstractKmerTableU[] tables(){return tables;}
	
	public long filterMemoryOverride=0;
	
	private final int bytesPerKmer;

	private final long usableMemory;
	private final long filterMemory0;
	private final long filterMemory1;
	private final long tableMemory;
	private final long estimatedKmerCapacity;
	
	/** Number of tables (and threads, during loading) */
	private final boolean prealloc;
	
	/** Number of tables (and threads, during loading) */
	public final int ways;
	
	/** Total kmer length */
	public final int kbig;
	/** Normal kmer length */
	public final int k;
	/** kbig-1; used in some expressions */
	public final int kbig2;
	/** Number of little kmers in a big kmer */
	public final int mult;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	
	/** Quality-trim the left side */
	public final boolean qtrimLeft;
	/** Quality-trim the right side */
	public final boolean qtrimRight;
	/** Trim bases at this quality or below.  Default: 4 */
	public final float trimq;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;
	
	/** Throw away reads below this average quality after trimming.  Default: 0 */
	public final float minAvgQuality;
	/** If positive, calculate average quality from the first X bases only.  Default: 0 */
	public final int minAvgQualityBases;
	
	/** Ignore kmers with probability of correctness less than this */
	public final float minProb;
	
	/** Correct via overlap */
	private final boolean ecco;
	
	/** Attempt to merge via overlap prior to counting kmers */
	private final boolean merge;
	
	/*--------------------------------------------------------------*/
	/*----------------            Walker            ----------------*/
	/*--------------------------------------------------------------*/
	
	public WalkerU walk(){
		return new WalkerKSTU();
	}
	
	public class WalkerKSTU extends WalkerU {
		
		WalkerKSTU(){
			w=tables[0].walk();
		}
		
		/** 
		 * Fills this object with the next key and value.
		 * @return True if successful.
		 */
		public boolean next(){
			if(w==null){return false;}
			if(w.next()){return true;}
			if(tnum<tables.length){tnum++;}
			w=(tnum<tables.length ? tables[tnum].walk() : null);
			return w==null ? false : w.next();
		}
		
		public Kmer kmer(){return w.kmer();}
		public int value(){return w.value();}
		
		private WalkerU w=null;

		/** current table number */
		private int tnum=0;
	}
	
}
