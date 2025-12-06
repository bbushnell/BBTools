package bbduk;

import java.io.File;
import java.io.PrintStream;
import aligner.SideChannel3;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import kmer.AbstractKmerTable;
import kmer.ScheduleMaker;
import shared.Tools;
import var2.ScafMap;
import var2.VcfLoader;

/**
 * Index and loader for BBDuk
 * @author Brian Bushnell
 * @date November 19, 2025
 *
 */
public class BBDukIndexMask extends BBDukIndex {

	/**
	 * Constructor.
	 * @param p Parser with command line arguments
	 */
	public BBDukIndexMask(BBDukParser p){
		
		scaffoldNames.add(""); //Necessary so that the first real scaffold gets an id of 1, not zero
		scaffoldLengths.add(0);
		
		silent=p.silent;
		ordered=p.ordered;
		boolean prealloc=p.prealloc;
		preallocFraction=p.preallocFraction;
		initialSize=p.initialSize;
		refNames=p.refNames;
		altRefNames=p.altRefNames;
		ref=p.ref;
		altref=p.altref;
		literal=p.literal;
		samref=p.samref;
		outrefstats=p.outrefstats;
		tossJunk=p.tossJunk;
		editDistance=p.editDistance;
		hammingDistance2=p.hammingDistance2;
		editDistance2=p.editDistance2;
		varFile=p.varFile;
		vcfFile=p.vcfFile;
		varMap=p.varMap;
		scafMap=p.scafMap;
		fixVariants=p.fixVariants;
		filterVars=p.filterVars;
		makeReadStats=p.makeReadStats;
		rcomp=p.rcomp;
		middleMask=p.middleMask;
		useTable=p.useTable;
		useArray=p.useArray;
		k=p.k;
		mink=p.mink;

		printNonZeroOnly=p.printNonZeroOnly;
		speed=p.speed;
		qSkip=p.qSkip;
		amino=p.amino;
		bitsPerBase=p.bitsPerBase;
		maxSymbol=p.maxSymbol;
		symbols=p.symbols;
		symbolMask=p.symbolMask;
		clearMasks=p.clearMasks;
		setMasks=p.setMasks;
		leftMasks=p.leftMasks;
		rightMasks=p.rightMasks;
		lengthMasks=p.lengthMasks;
		symbolToNumber=p.symbolToNumber;
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
		RELEASE_TABLES=BBDukParser.RELEASE_TABLES;
		int ways=p.WAYS;
		if(Integer.bitCount(ways)>1) {
			ways=Math.min(512, Integer.highestOneBit(ways<<1));
			System.err.println("ways changed to "+ways);
		}
		WAYS=ways;
		WAYMASK=ways-1;
		
		
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
		final int tableType=AbstractKmerTable.TABLE;
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, prealloc, (prealloc ? preallocFraction : 0.9));
		int[] schedule=scheduleMaker.makeSchedule();
		outstream.print("Allocating kmer table: \t");
		keySets=AbstractKmerTable.preallocate(WAYS, tableType, schedule, -1L);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean loaded() {
		return assortedLoaded && kmersLoaded;
	}
	
	synchronized void setKmersLoaded() {
		assert(!kmersLoaded);
		kmersLoaded=true;
	}
	
	@Override
	public synchronized void loadAssorted(String in1_for_header){
		assert(!assortedLoaded);
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
		assortedLoaded=true;
	}

	@Override
	public void cleanup() {
		
		/* Unload kmers to save memory */
		if(RELEASE_TABLES){unloadKmers();}
		
		/* Unload sequence data to save memory */
		if(RELEASE_TABLES){unloadScaffolds();}
	}
	
	/** Clear stored kmers. */
	@Override
	public void unloadKmers(){
		if(keySets!=null){
			for(int i=0; i<keySets.length; i++){keySets[i]=null;}
		}
	}
	
	/** Clear stored sequence data. */
	@Override
	public void unloadScaffolds(){
		if(scaffoldNames!=null && !scaffoldNames.isEmpty()){
			scaffoldNames.clear();
			scaffoldNames.trimToSize();
		}
		scaffoldReadCounts=null;
		scaffoldBaseCounts=null;
		scaffoldLengths=null;
	}
	
	/** Write statistics on a per-reference basis. */
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
	
	/** Fills the scaffold names array with reference names. */
	void toRefNames(){
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
	
	@Override
	void rebalance(int way) {//TODO:  Don't do this if near memory limit!
//		AbstractKmerTable map=keySets[way];
//		if(map.canRebalance() && map.size()>2L*map.arrayLength()){
//			map.rebalance();
//		}
	}
	
	void dump(ByteStreamWriter bsw, int minCount, int maxCount) {
		for(AbstractKmerTable set : keySets){
			set.dumpKmersAsBytes(bsw, k, minCount, maxCount, null);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds short kmers on the left end of the read.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param extraBase Base added to end in case of deletions
	 * @param id Scaffold number
	 * @return Number of kmers stored
	 */
	long addToMapLeftShift(long kmer, long rkmer, final long extraBase, final int id, final int tnum){
		if(verbose){outstream.println("addToMapLeftShift");}
		long added=0;
		for(int i=k-1; i>=mink; i--){
			kmer=kmer&rightMasks[i];
			rkmer=rkmer>>>bitsPerBase;
			long x=addToMap(kmer, rkmer, i, extraBase, id, lengthMasks[i], hammingDistance2, editDistance2, tnum);
			added+=x;
			if(verbose){
				final long value=toValue(kmer, rkmer, lengthMasks[i]);
				final long hash=hash(value);
				if((hash&WAYMASK)==tnum){
					outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added left-shift kmer "+kmerToString(kmer&~lengthMasks[i], i)+"; value="+(toValue(kmer, rkmer, lengthMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+lengthMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
					outstream.println("i="+i+"; tnum="+tnum+"; Looking for left-shift kmer "+kmerToString(kmer&~lengthMasks[i], i));
					final AbstractKmerTable map=keySets[tnum];
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
	long addToMapRightShift(long kmer, long rkmer, final int id, final int tnum){
		if(verbose){outstream.println("addToMapRightShift");}
		long added=0;
		for(int i=k-1; i>=mink; i--){
			long extraBase=kmer&symbolMask;
			kmer=kmer>>>bitsPerBase;
			rkmer=rkmer&rightMasks[i];
			long x=addToMap(kmer, rkmer, i, extraBase, id, lengthMasks[i], hammingDistance2, editDistance2, tnum);
			added+=x;
			if(verbose){
				final long value=toValue(kmer, rkmer, lengthMasks[i]);
				final long hash=hash(value);
				if((hash&WAYMASK)==tnum){
					outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added right-shift kmer "+kmerToString(kmer&~lengthMasks[i], i)+"; value="+(toValue(kmer, rkmer, lengthMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+lengthMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
					outstream.println("i="+i+"; tnum="+tnum+"; Looking for right-shift kmer "+kmerToString(kmer&~lengthMasks[i], i));
					final AbstractKmerTable map=keySets[tnum];
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
	long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0, final int hdist, final int edist, final int tnum){
		if(verbose){outstream.println("addToMap_A; len="+len+"; kMasks[len]="+lengthMasks[len]);}
		final AbstractKmerTable map=keySets[tnum];
		assert(kmask0==lengthMasks[len]) : kmask0+", "+len+", "+lengthMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(lengthMasks[len]);
		assert((kmer&kmask0)==0);
		final long added;
		if(hdist==0){
			final long key=toValue(kmer, rkmer, kmask0);
			final long hash=hash(key);
			if(verbose){outstream.println("toValue ("+kmerToString(kmer, len)+", "+kmerToString(rkmer, len)+") = "+kmerToString(key, len)+" = "+key);}
			if((hash&WAYMASK)!=tnum || failsSpeed(hash)){return 0;}
			if(verbose){outstream.println("addToMap_B: "+kmerToString(key, len)+" ("+key+")");}
			added=map.setIfNotPresent(key, id);
		}else if(edist>0){
			//				long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber2bases[i+1]]);
			added=mutate(kmer, rkmer, len, id, edist, extraBase, tnum);
		}else{
			added=mutate(kmer, rkmer, len, id, hdist, -1, tnum);
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
	long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase, final int tnum){
		final AbstractKmerTable map=keySets[tnum];
		long added=0;

		final long key=toValue(kmer, rkmer, lengthMasks[len]);
		final long hash=hash(key);

		//			if(dist==1){System.err.println(".\t.\t"+kmerToString(kmer, k)+" initial.");}//123

		if(verbose){outstream.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+lengthMasks[len]);}
		if((hash&WAYMASK)==tnum) {
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
						added+=mutate(temp, rtemp, len, id, dist2, extraBase, tnum);
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
							added+=mutate(temp, rtemp, len, id, dist2, -1, tnum);
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
							added+=mutate(temp, rtemp, len, id, dist2, eb2, tnum);
						}
					}
				}
			}

		}

		return added;
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
		
		final long key=toValue(kmer, rkmer, lengthMask);
		final long hash=hash(key);
		if(passesSpeed(hash)){
			if(verbose){outstream.println("Testing key "+kmerToString(key, len)+" ("+key+")");}
			AbstractKmerTable set=keySets[(int)(hash&WAYMASK)];
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
	private final long toValue(long kmer, long rkmer, long lengthMask){
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
	
	private static final long hash(long value) {
		return Tools.hash64shift(value)&Long.MAX_VALUE;
	}
	
	/**
	 * Computes reverse complement of kmer.
	 * @param kmer Input kmer
	 * @param len Kmer length
	 * @return Reverse complement kmer
	 */
	private final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
	/**
	 * Determines if kmer passes speed filtering threshold.
	 * @param hash Kmer hash value
	 * @return true if kmer should be processed
	 */
	final boolean passesSpeed(long hash){
		return speed<2 || (((hash>>16)&15)+1)>=speed;
	}
	
	/**
	 * Determines if kmer fails speed filtering threshold.
	 * @param hash Kmer hash value
	 * @return true if kmer should be skipped
	 */
	final boolean failsSpeed(long hash){//
		return speed>1 && (((hash>>16)&15)+1)<speed;
	}
	
	/** For verbose / debugging output */
	final String kmerToString(long kmer, int k){
		return amino ? AminoAcid.kmerToStringAA(kmer, k) : AminoAcid.kmerToString(kmer, k);
	}
	
	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
	@Override
	public int ways() {return WAYS;}
	
	@Override
	SideChannel3 sidechannel() {return sidechannel;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private boolean kmersLoaded=false;
	private boolean assortedLoaded=false;


	/** Number of tables (and threads, during loading) */
	private final int WAYS;
	private final long WAYMASK;

	/*--------------------------------------------------------------*/
	
	/** Fraction of available memory preallocated to arrays */
	private final double preallocFraction;
	/** Initial size of data structures */
	private final int initialSize;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	final AbstractKmerTable[] keySets;
	
	/*--------------------------------------------------------------*/
	/*----------------          Immutable           ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean silent;
	private final boolean ordered;
	
	/** Statistics output files */
	private final String outrefstats;
	
	final boolean tossJunk;
	
	/** Store reference kmers with up to this many edits (including indels) */
	private final int editDistance;
	/** Store short reference kmers with up to this many substitutions */
	private final int hammingDistance2;
	/** Store short reference kmers with up to this many edits (including indels) */
	private final int editDistance2;
	
	/*--------------------------------------------------------------*/
	/*----------------       Variant-Related        ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Document
	private final String varFile;
	private final String vcfFile;
	
	/** Filter reads with unsupported substitutions */
	private final boolean filterVars;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** True if a ReadStats object is being used for collecting data */
	private final boolean makeReadStats;
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** AND bitmask with 0's at the middle base */
	private final long middleMask;
	/** Use KmerTable data structure */
	private final boolean useTable;
	/** Use HashArray data structure (default) */
	private final boolean useArray;
	
	/** Normal kmer length */
	private final int k;
	/** Shortest kmer to use for trimming */
	private final int mink;
	
	/** Print only statistics for scaffolds that matched at least one read
	 * Default: true. */
	private final boolean printNonZeroOnly;
	
	/** Fraction of kmers to skip, 0 to 16 out of 17 */
	private final int speed;
	
	/** Skip this many kmers when examining the read.  Default 1.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int qSkip;
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	private final boolean amino;
	private final int bitsPerBase;
	private final int maxSymbol;
	private final int symbols;
	private final long symbolMask;
	
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
	/** Release memory used by kmer storage after processing reads */
	private final boolean RELEASE_TABLES;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Default initial size of data structures */
	private static final int initialSizeDefault=128000;
	/** Verbose messages */
	private static final boolean verbose=false;
	
}
