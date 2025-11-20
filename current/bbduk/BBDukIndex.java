package bbduk;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;

import aligner.SideChannel3;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import json.JsonObject;
import kmer.AbstractKmerTable;
import kmer.ScheduleMaker;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.IntList;
import structures.ListNum;
import var2.ScafMap;
import var2.VarMap;
import var2.VcfLoader;

/**
 * Index and loader for BBDuk
 * @author Brian Bushnell
 * @date November 19, 2025
 *
 */
public class BBDukIndex {
	
	/**
	 * Constructor.
	 * @param p Parser with command line arguments
	 */
	public BBDukIndex(BBDukParser p){
		
		scaffoldNames.add(""); //Necessary so that the first real scaffold gets an id of 1, not zero
		scaffoldLengths.add(0);
		
		silent=p.silent;
		json=p.json;
		boolean prealloc=p.prealloc;
		preallocFraction=p.preallocFraction;
		initialSize=p.initialSize;
		refNames=p.refNames;
		altRefNames=p.altRefNames;
		ALLOW_LOCAL_ARRAYS=p.ALLOW_LOCAL_ARRAYS;
		ref=p.ref;
		altref=p.altref;
		literal=p.literal;
		samref=p.samref;
		outrefstats=p.outrefstats;
		tossJunk=p.tossJunk;
		dump=p.dump;
		ordered=p.ordered;
		useShortKmers=p.useShortKmers;
		maskMiddle=p.maskMiddle;
		midMaskLen=p.midMaskLen;
		hammingDistance=p.hammingDistance;
		qHammingDistance=p.qHammingDistance;
		editDistance=p.editDistance;
		hammingDistance2=p.hammingDistance2;
		qHammingDistance2=p.qHammingDistance2;
		editDistance2=p.editDistance2;
		maxSkip=p.maxSkip;
		minSkip=p.minSkip;
		varFile=p.varFile;
		vcfFile=p.vcfFile;
		varMap=p.varMap;
		scafMap=p.scafMap;
		fixVariants=p.fixVariants;
		unfixVariants=p.unfixVariants;
		samFile=p.samFile;
		filterVars=p.filterVars;
		jsonStats=p.jsonStats;
		makeReadStats=p.makeReadStats;
		rcomp=p.rcomp;
		forbidNs=p.forbidNs;
		middleMask=p.middleMask;
		useForest=p.useForest;
		useTable=p.useTable;
		useArray=p.useArray;
		k=p.k;
		k2=p.k2;
		kbig=p.kbig;
		keff=p.keff;
		mink=p.mink;
		kfilter=p.kfilter;
		ktrimLeft=p.ktrimLeft;
		ktrimRight=p.ktrimRight;
		ktrimN=p.ktrimN;
		ksplit=p.ksplit;

		printNonZeroOnly=p.printNonZeroOnly;
		useRefNames=p.useRefNames;
		speed=p.speed;
		qSkip=p.qSkip;
		noAccel=p.noAccel;
		accel=p.accel;
		amino=p.amino;
		bitsPerBase=p.bitsPerBase;
		maxSymbol=p.maxSymbol;
		symbols=p.symbols;
		symbolMask=p.symbolMask;
		shift2=p.shift2;
		mask=p.mask;
		kmask=p.kmask;
		clearMasks=p.clearMasks;
		setMasks=p.setMasks;
		leftMasks=p.leftMasks;
		rightMasks=p.rightMasks;
		lengthMasks=p.lengthMasks;
		symbolToNumber=p.symbolToNumber;
		symbolToNumber0=p.symbolToNumber0;
		symbolToComplementNumber0=p.symbolToComplementNumber0;
		trimByOverlap=p.trimByOverlap;
		align=p.align;
		alignOut=p.alignOut;
		alignRef=p.alignRef;
		alignMinid1=p.alignMinid1;
		alignMinid2=p.alignMinid2;
		alignK1=p.alignK1;
		alignK2=p.alignK2;
		alignMM1=p.alignMM1;
		alignMM2=p.alignMM2;
		
		outstream=BBDukParser.outstream;
		overwrite=BBDukParser.overwrite;
		DISPLAY_PROGRESS=BBDukParser.DISPLAY_PROGRESS;
		THREADS=BBDukParser.THREADS;
		REPLICATE_AMBIGUOUS=BBDukParser.REPLICATE_AMBIGUOUS;
		RELEASE_TABLES=BBDukParser.RELEASE_TABLES;

		refScafCounts=new int[refNames.size()];
				
		if(ref!=null){
			for(String s0 : ref){
				assert(s0!=null) : "Specified a null reference.";
				String s=s0.toLowerCase();
				assert(s==null || s.startsWith("stdin") || s.startsWith("standardin") || new File(s0).exists()) : "Can't find "+s0;
			}
		}
		if(align) {
			sidechannel=new SideChannel3(alignRef, alignOut, null, alignK1, alignK2, 
					alignMinid1, alignMinid2, alignMM1, alignMM2, overwrite, ordered);
		}else {
			sidechannel=null;
		}
		
		//Initialize tables
		final int tableType=(useForest ? AbstractKmerTable.FOREST1D : useTable ? AbstractKmerTable.TABLE : useArray ? AbstractKmerTable.ARRAY1D : 0);
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, prealloc, (prealloc ? preallocFraction : 0.9));
		int[] schedule=scheduleMaker.makeSchedule();
		outstream.print("Allocating kmer table: \t");
		keySets=AbstractKmerTable.preallocate(WAYS, tableType, schedule, -1L);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void loadRef(String in1_for_header){
		
		if(samref!=null){
			scafMap=ScafMap.loadReference(samref, true);
		}
		
		if(varFile!=null || vcfFile!=null || filterVars){
			if(scafMap==null){scafMap=ScafMap.loadSamHeader(in1_for_header);}
			assert(scafMap!=null && scafMap.size()>0) : "No scaffold names were loaded.";
			if(varFile!=null){
				outstream.println("Loading variants.");
				varMap=VcfLoader.loadVarFile(varFile, scafMap);
			}else if(vcfFile!=null){
				outstream.println("Loading variants.");
				varMap=VcfLoader.loadVcfFile(vcfFile, scafMap, false, false);
			}
			fixVariants=(makeReadStats && varMap!=null && varMap.size()>0 && 
				scafMap!=null && scafMap.size()>0);
		}
		
		/* Start overall timer */
		Timer t=new Timer();
		
		loadRef(t.time1);
		
		/* Stop timer and calculate speed statistics */
		t.stop();
	}
	
	
	/**
	 * Core processing method that loads reference kmers and processes input reads.
	 * Fills kmer tables from reference sequences.
	 * @param startTime Start time in nanoseconds for timing calculations
	 */
	public void loadRef(long startTime){
		
		/* Start phase timer */
		Timer t=new Timer();
		
		/* Fill tables with reference kmers */
		if((ref!=null && ref.length>0) || (literal!=null && literal.length>0)){
			final boolean oldTI=FASTQ.TEST_INTERLEAVED; //TODO: This needs to be changed to a non-static field, or somehow 'read mode' and 'ref mode' need to be distinguished.
			final boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			final boolean oldPC=FASTQ.PARSE_CUSTOM;
			final boolean oldSplit=FastaReadInputStream.SPLIT_READS;
			final int oldML=FastaReadInputStream.MIN_READ_LEN;
			
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			FASTQ.PARSE_CUSTOM=false;
			FastaReadInputStream.SPLIT_READS=false;
			FastaReadInputStream.MIN_READ_LEN=1;
			
			
			storedKmers=spawnLoadThreads();
			if(storedKmers<1 && altRefNames.size()>0){
				outstream.println("Ref had no kmers; using alt ref.");
				ref=altref;
				refNames=altRefNames;
				refScafCounts=new int[refNames.size()];
				scaffoldNames.clear();
				scaffoldLengths.clear();
				storedKmers=spawnLoadThreads();
			}
			
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			FASTQ.PARSE_CUSTOM=oldPC;
			FastaReadInputStream.SPLIT_READS=oldSplit;
			FastaReadInputStream.MIN_READ_LEN=oldML;
			
			if(useRefNames){toRefNames();}
			t.stop();
		}
		
		{
			long ram=freeMemory();
			ALLOW_LOCAL_ARRAYS=(scaffoldNames!=null && Tools.max(THREADS, 1)*3*8*scaffoldNames.size()<ram*5);
		}
		
		/* Dump kmers to text */
		if(dump!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(dump, overwrite, false, true);
			bsw.start();
			for(AbstractKmerTable set : keySets){
				set.dumpKmersAsBytes(bsw, k, 0, Integer.MAX_VALUE, null);
			}
			bsw.poisonAndWait();
		}
		
		if(storedKmers<1 && (ktrimRight || ktrimLeft || ktrimN || ksplit)){
			outstream.println("******  WARNING! A KMER OPERATION WAS CHOSEN BUT NO KMERS WERE LOADED.  ******");
			if(ref==null && literal==null){
				outstream.println("******  YOU NEED TO SPECIFY A REFERENCE FILE OR LITERAL SEQUENCE.       ******\n");
			}else{
				outstream.println("******  PLEASE ENSURE K IS LESS THAN OR EQUAL TO REF SEQUENCE LENGTHS.  ******\n");
			}
			if(ktrimRight && trimByOverlap){
				ktrimRight=false;
			}else{
				assert(false) : "You can bypass this assertion with the -da flag.";
			}
		}
	}
	
	public void cleanup() {
		
		/* Unload kmers to save memory */
		if(RELEASE_TABLES){unloadKmers();}
		
		/* Unload sequence data to save memory */
		if(RELEASE_TABLES){unloadScaffolds();}
	}
	
	/**
	 * Clear stored kmers.
	 */
	public void unloadKmers(){
		if(keySets!=null){
			for(int i=0; i<keySets.length; i++){keySets[i]=null;}
		}
	}
	
	/**
	 * Clear stored sequence data.
	 */
	public void unloadScaffolds(){
		if(scaffoldNames!=null && !scaffoldNames.isEmpty()){
			scaffoldNames.clear();
			scaffoldNames.trimToSize();
		}
		scaffoldReadCounts=null;
		scaffoldBaseCounts=null;
		scaffoldLengths=null;
	}
	
	/**
	 * Write statistics on a per-reference basis.
	 */
	void writeRefStats(String in1, String in2, long readsIn){
		if(outrefstats==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(outrefstats, overwrite, false, false);
		tsw.start();
		
		/* Count mapped reads */
		long mapped=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mapped+=scaffoldReadCounts.get(i);
		}
		
		final int numRefs=refNames.size();
		long[] refReadCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		long[] refLengths=new long[numRefs];
		
		for(int r=0, s=1; r<numRefs; r++){
			final int lim=s+refScafCounts[r];
			while(s<lim){
				refReadCounts[r]+=scaffoldReadCounts.get(s);
				refBaseCounts[r]+=scaffoldBaseCounts.get(s);
				refLengths[r]+=scaffoldLengths.get(s);
				s++;
			}
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mapped));
		tsw.print(Tools.format("#References\t%d\n",Tools.max(0, refNames.size())));
		tsw.print("#Name\tLength\tScaffolds\tBases\tCoverage\tReads\tRPKM\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=0; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long bases=refBaseCounts[i];
			final long len=refLengths[i];
			final int scafs=refScafCounts[i];
			final String name=ReadWrite.stripToCore(refNames.get(i));
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\n",name,len,scafs,bases,bases*invlen,reads,reads*mult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Fills the scaffold names array with reference names.
	 */
	private void toRefNames(){
		final int numRefs=refNames.size();
		for(int r=0, s=1; r<numRefs; r++){
			final int scafs=refScafCounts[r];
			final int lim=s+scafs;
			final String name=ReadWrite.stripToCore(refNames.get(r));
			while(s<lim){
				scaffoldNames.set(s, name);
				s++;
			}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long spawnLoadThreads(){
		Timer t=new Timer();
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		/* For each reference file... */
		int refNum=0;
		if(ref!=null){
			for(String refname : ref){

				/* Start an input stream */
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null, null, null, Shared.USE_MPI, true);
				cris.start(); //4567
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				/* Iterate through read lists from the input stream */
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					{
						/* Assign a unique ID number to each scaffold */
						ArrayList<Read> reads2=new ArrayList<Read>(reads);
						for(Read r1 : reads2){
							final Read r2=r1.mate;
							final Integer id=scaffoldNames.size();
							refScafCounts[refNum]++;
							scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
							int len=r1.length();
							r1.obj=id;
							if(r2!=null){
								r2.obj=id;
								len+=r2.length();
							}
							scaffoldLengths.add(len);
						}
						
						if(REPLICATE_AMBIGUOUS){
							reads2=Tools.replicateAmbiguous(reads2, Tools.min(k, mink));
						}

						/* Send a pointer to the read list to each LoadThread */
						for(LoadThread lt : loaders){
							boolean b=true;
							while(b){
								try {
									lt.queue.put(reads2);
									b=false;
								} catch (InterruptedException e) {
									//TODO:  This will hang due to still-running threads.
									throw new RuntimeException(e);
								}
							}
						}
					}

					/* Dispose of the old list and fetch a new one */
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				/* Cleanup */
				cris.returnList(ln);
				errorState|=ReadWrite.closeStream(cris);
				refNum++;
			}
		}

		/* If there are literal sequences to use as references */
		if(literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(literal.length);
			if(verbose){outstream.println("Adding literals "+Arrays.toString(literal));}

			/* Assign a unique ID number to each literal sequence */
			for(int i=0; i<literal.length; i++){
				final Integer id=scaffoldNames.size();
				final Read r=new Read(literal[i].getBytes(), null, id);
				refScafCounts[refNum]++;
				scaffoldNames.add(id.toString());
				scaffoldLengths.add(r.length());
				r.obj=id;
				list.add(r);
			}
			
			if(REPLICATE_AMBIGUOUS){
				list=Tools.replicateAmbiguous(list, Tools.min(k, mink));
			}

			/* Send a pointer to the read list to each LoadThread */
			for(LoadThread lt : loaders){
				boolean b=true;
				while(b){
					try {
						lt.queue.put(list);
						b=false;
					} catch (InterruptedException e) {
						//TODO:  This will hang due to still-running threads.
						throw new RuntimeException(e);
					}
				}
			}
		}
		
		/* Signal loaders to terminate */
		for(LoadThread lt : loaders){
			boolean b=true;
			while(b){
				try {
					lt.queue.put(POISON);
					b=false;
				} catch (InterruptedException e) {
					//TODO:  This will hang due to still-running threads.
					throw new RuntimeException(e);
				}
			}
		}
		
		/* Wait for loaders to die, and gather statistics */
		boolean success=true;
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					lt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			added+=lt.addedT;
			refKmers+=lt.refKmersT;
			refBases+=lt.refBasesT;
			refReads+=lt.refReadsT;
			success&=lt.success;
		}
		if(!success){KillSwitch.kill("Failed loading ref kmers; aborting.");}
		
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=WAYS;
		refBases/=WAYS;
		refReads/=WAYS;
		
		scaffoldReadCounts=new AtomicLongArray(scaffoldNames.size());
		scaffoldBaseCounts=new AtomicLongArray(scaffoldNames.size());

		t.stop();
		if(DISPLAY_PROGRESS && !json){
			outstream.println("Added "+added+" kmers; time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		
		if(verbose){
			TextStreamWriter tsw=new TextStreamWriter("stdout", false, false, false, FileFormat.TEXT);
			tsw.start();
			for(AbstractKmerTable table : keySets){
				table.dumpKmersAsText(tsw, k, 1, Integer.MAX_VALUE);
			}
			tsw.poisonAndWait();
		}
		
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads kmers into a table.  Each thread handles all kmers X such that X%WAYS==tnum.
	 */
	private class LoadThread extends Thread{
		
		public LoadThread(final int tnum_){
			tnum=tnum_;
			map=keySets[tnum];
		}
		
		/**
		 * Get the next list of reads (or scaffolds) from the queue.
		 * @return List of reads
		 */
		private ArrayList<Read> fetch(){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=queue.take();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			return list;
		}
		
		@Override
		public void run(){
			ArrayList<Read> reads=fetch();
			while(reads!=POISON){
				for(Read r1 : reads){
					assert(r1.pairnum()==0);
					final Read r2=r1.mate;

					final int rblen=(r1==null ? 0 : r1.length());
					final int rblen2=r1.mateLength();
					
					addedT+=addToMap(r1, rblen>20000000 ? k : rblen>5000000 ? 11 : rblen>500000 ? 2 : 0);
					if(r2!=null){
						addedT+=addToMap(r2, rblen2>20000000 ? k : rblen2>5000000 ? 11 : rblen2>500000 ? 2 : 0);
					}
				}
				reads=fetch();
			}
			
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
			success=true;
		}

		/**
		 * Store the read's kmers in a table.
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(Read r, int skip){
			skip=Tools.max(minSkip, Tools.min(maxSkip, skip));
			final byte[] bases=r.bases;
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;
			
			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return 0;}
			
			final int id=(Integer)r.obj;
			
			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
					if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
					if(verbose){outstream.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
							added+=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
							if(useShortKmers){
								if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
								if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
							}
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final long x=symbolToNumber0[b];
					final long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
					if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
					if(verbose){
						if(verbose){
							String fwd=new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k));
							String rev=AminoAcid.reverseComplementBases(fwd);
							String fwd2=kmerToString(kmer, Tools.min(len, k));
							outstream.println("fwd="+fwd+", fwd2="+fwd2+", rev="+rev+", kmer="+kmer+", rkmer="+rkmer);
							outstream.println("b="+(char)b+", x="+x+", x2="+x2+", bitsPerBase="+bitsPerBase+", shift2="+shift2);
							if(!amino){
								assert(AminoAcid.stringToKmer(fwd)==kmer) : fwd+", "+AminoAcid.stringToKmer(fwd)+", "+kmer+", "+len;
								if(len>=k){
									assert(rcomp(kmer, Tools.min(len, k))==rkmer);
									assert(rcomp(rkmer, Tools.min(len, k))==kmer);
									assert(AminoAcid.kmerToString(kmer, Tools.min(len, k)).equals(fwd));
									assert(AminoAcid.kmerToString(rkmer, Tools.min(len, k)).equals(rev)) : AminoAcid.kmerToString(rkmer, Tools.min(len, k))+" != "+rev+" (rkmer)";
								}
								assert(fwd.equalsIgnoreCase(fwd2)) : fwd+", "+fwd2; //may be unsafe
							}
							outstream.println("Scanning6 i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+", bases="+fwd+", rbases="+rev);
						}
					}
					if(len>=k){
						refKmersT++;
						final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
						final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
						added+=atm;
						if(useShortKmers){
							if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
							if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
						}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the left end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapLeftShift(long kmer, long rkmer, final long extraBase, final int id){
			if(verbose){outstream.println("addToMapLeftShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				kmer=kmer&rightMasks[i];
				rkmer=rkmer>>>bitsPerBase;
				long x=addToMap(kmer, rkmer, i, extraBase, id, lengthMasks[i], hammingDistance2, editDistance2);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, lengthMasks[i]))%WAYS==tnum){
						outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added left-shift kmer "+kmerToString(kmer&~lengthMasks[i], i)+"; value="+(toValue(kmer, rkmer, lengthMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+lengthMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						outstream.println("i="+i+"; tnum="+tnum+"; Looking for left-shift kmer "+kmerToString(kmer&~lengthMasks[i], i));
						final long value=toValue(kmer, rkmer, lengthMasks[i]);
						if(map.contains(value)){outstream.println("Found "+value);}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the right end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapRightShift(long kmer, long rkmer, final int id){
			if(verbose){outstream.println("addToMapRightShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				long extraBase=kmer&symbolMask;
				kmer=kmer>>>bitsPerBase;
				rkmer=rkmer&rightMasks[i];
				long x=addToMap(kmer, rkmer, i, extraBase, id, lengthMasks[i], hammingDistance2, editDistance2);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, lengthMasks[i]))%WAYS==tnum){
						outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added right-shift kmer "+kmerToString(kmer&~lengthMasks[i], i)+"; value="+(toValue(kmer, rkmer, lengthMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+lengthMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						outstream.println("i="+i+"; tnum="+tnum+"; Looking for right-shift kmer "+kmerToString(kmer&~lengthMasks[i], i));
						final long value=toValue(kmer, rkmer, lengthMasks[i]);
						if(map.contains(value)){outstream.println("Found "+value);}
					}
				}
			}
			return added;
		}
		
		
		/**
		 * Adds this kmer to the table, including any mutations implied by editDistance or hammingDistance.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param len Kmer length
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0
		 * @return Number of kmers stored
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0, final int hdist, final int edist){
			
			assert(kmask0==lengthMasks[len]) : kmask0+", "+len+", "+lengthMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(lengthMasks[len]);
			
			if(verbose){outstream.println("addToMap_A; len="+len+"; kMasks[len]="+lengthMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hdist==0){
				final long key=toValue(kmer, rkmer, kmask0);
				if(verbose){outstream.println("toValue ("+kmerToString(kmer, len)+", "+kmerToString(rkmer, len)+") = "+kmerToString(key, len)+" = "+key);}
				if(failsSpeed(key)){return 0;}
				if(key%WAYS!=tnum){return 0;}
				if(verbose){outstream.println("addToMap_B: "+kmerToString(key, len)+" ("+key+")");}
				added=map.setIfNotPresent(key, id);
			}else if(edist>0){
//				long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber2bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, edist, extraBase);
			}else{
				added=mutate(kmer, rkmer, len, id, hdist, -1);
			}
			if(verbose){outstream.println("addToMap added "+added+" keys.");}
			return added;
		}
		
		/**
		 * Mutate and store this kmer through 'dist' recursions.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @param dist Number of mutations
		 * @param extraBase Base added to end in case of deletions
		 * @return Number of kmers stored
		 */
		private long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase){
			long added=0;
			
			final long key=toValue(kmer, rkmer, lengthMasks[len]);
			
//			if(dist==1){System.err.println(".\t.\t"+kmerToString(kmer, k)+" initial.");}//123
			
			if(verbose){outstream.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+lengthMasks[len]);}
			if(key%WAYS==tnum){
				if(verbose){outstream.println("mutate_B: "+kmerToString(kmer&~lengthMasks[len], len)+" = "+key);}
				int x=map.setIfNotPresent(key, id);
//				if(x>0){System.err.println(".\t.\t"+kmerToString(kmer, k)+" Added!");}//123
				if(verbose){outstream.println("mutate_B added "+x+" keys.");}
				added+=x;
				assert(map.contains(key));
			}
			
			if(dist>0){
				final int dist2=dist-1;
				
				//Sub
				for(int j=0; j<symbols; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i]; //TODO:  6/14/23, fixed incorrect description of setMasks that swapped i and j; may need changing in Seal and etc
						if(temp!=kmer){
							long rtemp=rcomp(temp, len);
							added+=mutate(temp, rtemp, len, id, dist2, extraBase);
						}
					}
				}
				
				if(editDistance>0){
					//Del
					if(extraBase>=0 && extraBase<=maxSymbol){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<bitsPerBase)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=rcomp(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1);
							}
						}
					}

					//Ins
					final long eb2=kmer&symbolMask;
					for(int i=1; i<len; i++){
						final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>bitsPerBase);
						for(int j=0; j<symbols; j++){
							final long temp=temp0|setMasks[j][i-1];
							if(temp!=kmer){
								long rtemp=rcomp(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, eb2);
							}
						}
					}
				}

			}
			
			return added;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Number of kmers stored by this thread */
		public long addedT=0;
		/** Number of items encountered by this thread */
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		/** Thread number; used to determine which kmers to store */
		public final int tnum;
		/** Buffer of input read lists */
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);
		
		/** Destination for storing kmers */
		private final AbstractKmerTable map;
		
		/** Completed successfully */
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Index Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Transforms a kmer into all canonical values for a given Hamming distance.
	 * Returns the related id stored in the tables.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @param qPos Position of kmer in query
	 * @param len kmer length
	 * @param qHDist Hamming distance
	 * @param sets Kmer hash tables
	 * @return Value stored in table, or -1
	 */
	final int getValue(final long kmer, final long rkmer, final long lengthMask, final int qPos, final int len, final int qHDist){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		if(verbose){outstream.println("getValue()");}
		int id=getValueInner(kmer, rkmer, lengthMask, len, qPos);
		if(id<1 && qHDist>0){
			final int qHDist2=qHDist-1;
			
			//Sub
			for(int j=0; j<symbols && id<1; j++){
				for(int i=0; i<len && id<1; i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=rcomp(temp, len);
						id=getValue(temp, rtemp, lengthMask, qPos, len, qHDist2);
					}
				}
			}
		}
		return id;
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table and search.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @param qPos Position of kmer in query
	 * @param sets Kmer hash tables
	 * @return Value stored in table
	 */
	private final int getValueInner(final long kmer, final long rkmer, final long lengthMask, final int len, final int qPos){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		if(qSkip>1 && (qPos%qSkip!=0)){return -1;}

		if(verbose){
			outstream.println("getValueInner(kmer="+AminoAcid.kmerToString(kmer, len)+", rkmer="+AminoAcid.kmerToString(rkmer, len)+", len="+len+", mask="+lengthMask+")");
			outstream.println("getValueInner(kmer="+kmer+", rkmer="+rkmer+", len="+len+", mask="+lengthMask+")");
		}
		final long max=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		if(verbose){outstream.println("max="+AminoAcid.kmerToString(max, len)+" ("+max+")");}
		final long key=(max&middleMask)|lengthMask;
		if(verbose){outstream.println("key="+AminoAcid.kmerToString(key, len)+" ("+key+")");}
		if(passesSpeed(key)){
			if(verbose){outstream.println("Testing key "+kmerToString(key, len)+" ("+key+")");}
			AbstractKmerTable set=keySets[(int)(key%WAYS)];
			final int id=set.getValue(key);
			if(verbose){outstream.println("getValueInner("+kmerToString(kmer, len)+", "+kmerToString(rkmer, len)+") > "+kmerToString(key, len)+" ("+key+") = "+id);}
			return id;
		}
		if(verbose){outstream.println("Invalid key.");}
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Current available memory */
	private static final long freeMemory(){
		Runtime rt=Runtime.getRuntime();
		return rt.freeMemory();
	}
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : 
			"\n"+Long.toBinaryString(lengthMask)+
			"\n"+Long.toBinaryString(kmer)+
			"\n"+Long.toBinaryString(rkmer)+
			"\n"+Long.toBinaryString(rcomp(kmer, k));
		if(verbose){outstream.println("toValue("+AminoAcid.kmerToString(kmer, k)+", "+AminoAcid.kmerToString(rkmer, k)+", "+lengthMask+")");}
		final long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		if(verbose){outstream.println("value="+AminoAcid.kmerToString(value, k)+" = "+value);}
		final long ret=(value&middleMask)|lengthMask;
		if(verbose){outstream.println("ret="+AminoAcid.kmerToString(ret, k)+" = "+ret);}
		return ret;
	}
	
	/**
	 * Computes reverse complement of kmer.
	 * @param kmer Input kmer
	 * @param len Kmer length
	 * @return Reverse complement kmer
	 */
	final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
	/**
	 * Determines if kmer passes speed filtering threshold.
	 * @param key Kmer hash value
	 * @return true if kmer should be processed
	 */
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}
	
	/**
	 * Determines if kmer fails speed filtering threshold.
	 * @param key Kmer hash value
	 * @return true if kmer should be skipped
	 */
	final boolean failsSpeed(long key){
		return speed>0 && ((key&Long.MAX_VALUE)%17)<speed;
	}
	
	/** For verbose / debugging output */
	final String kmerToString(long kmer, int k){
		return amino ? AminoAcid.kmerToStringAA(kmer, k) : AminoAcid.kmerToString(kmer, k);
	}
	
	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Has this class encountered errors while processing? */
	public boolean errorState;
	
	/** Stores JSON output */
	private final JsonObject jsonStats;
	
	/** Number of reference reads processed for kmer loading */
	long refReads;
	/** Number of reference bases processed for kmer loading */
	long refBases;
	/** Number of reference kmers encountered during loading */
	long refKmers;
	/** Number of unique kmers actually stored in hash tables */
	long storedKmers;
	
	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */
	AtomicLongArray scaffoldReadCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */
	AtomicLongArray scaffoldBaseCounts;
	/** scaffoldLengths[id] stores the length of that scaffold */
	IntList scaffoldLengths=new IntList();
	/** hitCounts[x] stores the number of reads with exactly x kmer matches */
	private long[] hitCounts;
	/** Set to false to force threads to share atomic counter arrays. */
	boolean ALLOW_LOCAL_ARRAYS;
	
	/*--------------------------------------------------------------*/
	
	/** Fraction of available memory preallocated to arrays */
	private final double preallocFraction;
	/** Initial size of data structures */
	private final int initialSize;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	final AbstractKmerTable[] keySets;
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	ArrayList<String> refNames;
	final ArrayList<String> altRefNames;
	/** Number of scaffolds per reference. */
	int[] refScafCounts;
	/** Array of reference files from which to load kmers */
	String[] ref;
	/** Alternate reference to be used if main reference has no kmers */
	private final String[] altref;
	/** Array of literal strings from which to load kmers */
	private final String[] literal;
	/** Optional reference for sam file */
	private final String samref;
	
	/*--------------------------------------------------------------*/
	/*----------------          Immutable           ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean silent;
	private final boolean json;
	
	/** Statistics output files */
	private final String outrefstats;
	
	final boolean tossJunk;
	
	/** Dump kmers here. */
	private final String dump;
	
	/** Output reads in input order.  May reduce speed. */
	private final boolean ordered;
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	private final boolean useShortKmers;
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	private final boolean maskMiddle;
	private final int midMaskLen;
	
	/** Store reference kmers with up to this many substitutions */
	private final int hammingDistance;
	/** Search for query kmers with up to this many substitutions */
	private final int qHammingDistance;
	/** Store reference kmers with up to this many edits (including indels) */
	private final int editDistance;
	/** Store short reference kmers with up to this many substitutions */
	private final int hammingDistance2;
	/** Search for short query kmers with up to this many substitutions */
	private final int qHammingDistance2;
	/** Store short reference kmers with up to this many edits (including indels) */
	private final int editDistance2;
	/** Never skip more than this many consecutive kmers when hashing reference. */
	private final int maxSkip;
	/** Always skip at least this many consecutive kmers when hashing reference.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int minSkip;
	
	/*--------------------------------------------------------------*/
	/*----------------       Variant-Related        ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	private final String varFile;
	private final String vcfFile;
	VarMap varMap;
	ScafMap scafMap;
	boolean fixVariants;
	private final boolean unfixVariants;
	
	/** Optional file for quality score recalibration */
	private final String samFile;
	
	/** Filter reads with unsupported substitutions */
	private final boolean filterVars;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** True if a ReadStats object is being used for collecting data */
	private final boolean makeReadStats;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** Don't allow a read 'N' to match a reference 'A'.
	 * Reduces sensitivity when hdist>0 or edist>0.  Default: false. */
	private final boolean forbidNs;
	/** AND bitmask with 0's at the middle base */
	private final long middleMask;
	/** Use HashForest data structure */
	private final boolean useForest;
	/** Use KmerTable data structure */
	private final boolean useTable;
	/** Use HashArray data structure (default) */
	private final boolean useArray;
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Emulated kmer greater than k */
	private final int kbig;
	/** Effective kmer size */
	private final int keff;
	/** Shortest kmer to use for trimming */
	private final int mink;
	
	/** Filter reads by whether or not they have matching kmers */
	private final boolean kfilter;
	/** Trim matching kmers and all bases to the left */
	private final boolean ktrimLeft;
	/** Trim matching kmers and all bases to the right */
	boolean ktrimRight;
	/** Don't trim, but replace matching kmers with a symbol (default N) */
	private final boolean ktrimN;
	/** Split into two reads around the kmer */
	private final boolean ksplit;
	
	/** Print only statistics for scaffolds that matched at least one read
	 * Default: true. */
	private final boolean printNonZeroOnly;
	/** Use names of reference files instead of scaffolds.
	 * Default: false. */
	private final boolean useRefNames;
	
	/** Fraction of kmers to skip, 0 to 16 out of 17 */
	private final int speed;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int qSkip;
	
	/** noAccel is true if speed and qSkip are disabled, accel is the opposite. */
	private final boolean noAccel;
	private final boolean accel;
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	private final boolean amino;
	private final int bitsPerBase;
	private final int maxSymbol;
	private final int symbols;
	private final long symbolMask;
	private final int shift2;
	private final long mask;
	private final long kmask;
	
	/** x&clearMasks[i] will clear base i */
	private final long[] clearMasks;
	/** x|setMasks[j][i] will set position i to symbol j */
	private final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	private final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	private final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private final long[] lengthMasks;
	
	/** Symbol code; -1 for undefined */
	private final byte[] symbolToNumber;
	/** Symbol code; 0 for undefined */
	private final byte[] symbolToNumber0;
	/** Complementary symbol code; 0 for undefined */
	private final byte[] symbolToComplementNumber0;
	
	/*--------------------------------------------------------------*/
	/*----------------         BBMerge Flags        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Trim implied adapters based on overlap, for reads with insert size shorter than read length */
	private final boolean trimByOverlap;
	
	/*--------------------------------------------------------------*/
	/*----------------         Side Channel         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final boolean align;
	private final String alignOut;
	private final String alignRef;
	private final float alignMinid1;
	private final float alignMinid2;
	private final int alignK1;
	private final int alignK2;
	private final int alignMM1;
	private final int alignMM2;
	final SideChannel3 sidechannel;
	
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	private final PrintStream outstream;
	/** Permission to overwrite existing files */
	private final boolean overwrite;
	/** Display progress messages such as memory usage */
	private final boolean DISPLAY_PROGRESS;
	/** Number of ProcessThreads */
	private final int THREADS;
	/** Release memory used by kmer storage after processing reads */
	private final boolean RELEASE_TABLES;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	private final boolean REPLICATE_AMBIGUOUS;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of tables (and threads, during loading) */
	private static final int WAYS=7; //123
	/** Default initial size of data structures */
	private static final int initialSizeDefault=128000;
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Max value of hitCount array */
	private static final int HITCOUNT_LEN=1000;
	/** Indicates end of input stream */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	
}
