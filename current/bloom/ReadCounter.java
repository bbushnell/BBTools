package bloom;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.BitSet;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import jgi.Dedupe;
import kmer.KmerTableSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.SketchObject;
import stream.Streamer;
import stream.StreamerFactory;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import ukmer.Kmer;

/**
 * Builds packed counting tables from streamed sequence input using counting workers.
 * Nucleotide paths support primitive short keys or multiword hashes; amino-acid
 * counting uses a separate encoding. Short-key sampling, duplicate suppression,
 * quality filtering and long-key counting have different option coverage.
 * Paired-overlap processing requires nucleotide input.
 * Amino dispatch counts only the root, although processed-read totals include mates.
 * Short-key sampling bypasses quality, canonicalization and IDMASK rules. In
 * unsampled short-key counting, IDMASK selects sorted occurrence indices with
 * duplicate suppression and window-ending base indices otherwise. Long-key counting omits short-path sampling,
 * duplicate suppression and minQuality filtering.
 *
 * Each invocation owns its reader and workers but borrows any destination or
 * trusted table. Counting can mutate temporary reads through overlap processing
 * or trusted-base masking. Caught counting/flush failures propagate after worker
 * aggregation and reader closure; they do not roll back borrowed-table updates.
 * Reader/worker construction and startup are not protected by cleanup finally blocks.
 * Multiworker destinations must support concurrent increments; trusted tables and
 * shared counting settings, including BUFFERED, must remain stable during counting.
 * Inherited options/statistics are process-wide, so separate instances do not
 * isolate concurrent invocations with different settings.
 *
 * @author Brian Bushnell
 * @date Jul 5, 2012
 *
 */
public class ReadCounter extends KmerCountAbstract{
	
	/**
	 * Program entry point for standalone k-mer counting operations.
	 * Uses positional input paths followed by counting options. A second positional
	 * slot is reserved for the mate path; options are parsed starting at index two.
	 * A case-insensitive literal null in the mate slot selects unpaired input with options.
	 * @param args Primary path, optional/null mate placeholder, then named options
	 */
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object(){}.getClass().getEnclosingClass(), false);
			args=pp.args;
//			outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		String fname1=args[0];
		String fname2=(args.length>1 && !"null".equalsIgnoreCase(args[1]) ? args[1] : null);
		int k=14;
		int cbits=16;
		int matrixbits=-1;
		int hashes=1;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes")){
				hashes=Integer.parseInt(b);
			}else if(a.equals("canonical")){
				CANONICAL=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		int kbits=Tools.min(2*k, 62);
		if(matrixbits<0){
			matrixbits=kbits;
		}
		matrixbits=Tools.min(kbits, matrixbits);
		
		if(fileIO.FileFormat.hasFastaExtension(fname1)){
			assert(!FastaReadInputStream.SPLIT_READS);
			FastaReadInputStream.MIN_READ_LEN=k;
		}
		
		KCountArray counts=KCountArray.makeNew(1L<<matrixbits, cbits, hashes);
		ReadCounter rc=new ReadCounter(k, true, false, false, false);
		try{
			rc.count(fname1, fname2, counts);
		}catch(Exception e){
			throw new RuntimeException("Kmer counting failed; no complete count array is available.", e);
		}
		counts.shutdown();
		
		t.stop();
		System.out.println("Finished counting; time = "+t);
		
		rc.printStatistics(counts);
		
	}
	
	/** Selects nucleotide counting with reverse-complement matching and no overlap processing. */
	public ReadCounter(final int k_){this(k_, true, false, false, false);}
	
	/**
	 * Selects a representation and counting policy; does not open input or allocate a table.
	 *
	 * @param k_ K-mer length for counting
	 * @param rcomp_ Canonicalize short nucleotide keys; long-key counting uses Kmer's own representation
	 * @param ecco_ Correct paired overlaps when merge is not selected
	 * @param merge_ Count a joined sequence when a positive paired insert is found
	 * @param amino_ Whether to process amino acid sequences instead of nucleotides
	 */
	public ReadCounter(final int k_, final boolean rcomp_, 
			final boolean ecco_, final boolean merge_, final boolean amino_){
		this(k_, rcomp_, ecco_, merge_, amino_, false);
	}

	/** Creates a counter with optional independent secondary hashes for Bloom lanes.
	 * Requires positive k; amino-acid packed keys must fit below 64 bits and cannot
	 * request reverse complements. Short rolling paths also check inherited CANONICAL.
	 * Long keys use ukmer.Kmer geometry; constructor bit fields alone do not define it. */
	public ReadCounter(final int k_, final boolean rcomp_, final boolean ecco_, final boolean merge_,
			final boolean amino_, final boolean dualHash_){
		k=k_;
		k2=k-1;
		rcomp=rcomp_;
		ecco=ecco_;
		merge=merge_;
		amino=amino_;
		dualHash=dualHash_;

		final int bitsPerChar=(amino ? AminoAcid.AMINO_SHIFT : 2);
		aminoShift=AminoAcid.AMINO_SHIFT;
		shift=bitsPerChar*k;
		shift2=shift-bitsPerChar;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		assert(!amino || k*bitsPerChar<64);
		assert(!amino || !rcomp);
		assert(k>0);
	}

	/**
	 * Prints a table-cell histogram and ratios using shared keysCounted.
	 * The historical "Unique" label denotes occupied cells, not exact distinct keys
	 * for a hashed table. Shared counters need not describe just this supplied table.
	 * Empty occupied-cell statistics are reported as zero; histogram totals use longs.
	 * @param counts Readable table whose cell-frequency histogram will be summarized
	 */
	public void printStatistics(KCountArray counts){
		long[] freq=counts.transformToFrequency();
		
		long sum=sum(freq);
		System.out.println("Kmer fraction:");
		int lim1=8, lim2=16;
		for(int i=0; i<lim1; i++){
			String prefix=i+"";
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+Tools.format("%.3f%%   ",(100l*freq[i]/(double)sum))+"\t"+freq[i]);
		}
		while(lim1<=freq.length){
			long x=0;
			for(int i=lim1; i<lim2; i++){
				x+=freq[i];
			}
			String prefix=lim1+"-"+(lim2-1);
			if(lim2>=freq.length){prefix=lim1+"+";}
			while(prefix.length()<8){prefix=prefix+" ";}
			System.out.println(prefix+"\t"+Tools.format("%.3f%%   ",(100l*x/(double)sum))+"\t"+x);
			lim1*=2;
			lim2=min(lim2*2, freq.length);
		}
		
		long sum2=sum-freq[0];
		long x=freq[1];
		System.out.println();
		System.out.println("Keys Counted:  \t         \t"+keysCounted);
		System.out.println("Unique:        \t         \t"+sum2);
		System.out.println("Avg Sites/Key: \t         \t"+Tools.format("%.3f    ",(sum2==0 ? 0 : keysCounted*1d/sum2)));
		System.out.println();
		System.out.println("Singleton:     \t"+Tools.format("%.3f%%   ",(sum2==0 ? 0 : 100l*x/(double)sum2))+"\t"+x);
		x=sum2-x;
		System.out.println("Useful:        \t"+Tools.format("%.3f%%   ",(sum2==0 ? 0 : 100l*x/(double)sum2))+"\t"+x);
	}
	
	/**
	 * Pairs only the first entry from each primary/mate list. Remaining entries from
	 * both lists are appended as independent extra inputs, followed by extraFiles.
	 * Does not modify the supplied lists; delegates table creation and error handling.
	 *
	 * @param fname1 Primary input files as ArrayList
	 * @param fname2 Secondary input files as ArrayList (for paired reads)
	 * @param extraFiles Additional input files to process
	 * @param cbits Bits per count cell in the array
	 * @param cells Total number of cells in the count array
	 * @param hashes Number of hash functions to use
	 * @param minqual Minimum quality score for k-mer inclusion
	 * @param maxreads Reader-dependent input limit; see the delegated makeKca overload
	 * @param passes Number of processing passes for multi-pass counting
	 * @param thresh1 Threshold for first pass filtering
	 * @param thresh2 Threshold for final pass filtering
	 * @param prefilter Optional prefilter for two-stage counting
	 * @param prefilterLimit_ Prefilter threshold limit
	 * @return Populated KCountArray with k-mer counts
	 */
	public KCountArray makeKca_als(ArrayList<String> fname1, ArrayList<String> fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int thresh1, int thresh2, 
			KCountArray prefilter, int prefilterLimit_){
		String a=null, b=null;
		ArrayList<String> list=new ArrayList<String>();
		if(fname1!=null){
			for(int i=0; i<fname1.size(); i++){
				if(i==0){a=fname1.get(i);}else{list.add(fname1.get(i));}
			}
		}
		if(fname2!=null){
			for(int i=0; i<fname2.size(); i++){
				if(i==0){b=fname2.get(i);}else{list.add(fname2.get(i));}
			}
		}
		if(extraFiles!=null){
			for(String s : extraFiles){
				list.add(s);
			}
		}
		return makeKca(a, b, list.isEmpty() ? null : list, cbits, cells, hashes, minqual, 
				maxreads, passes, thresh1, thresh2,
				prefilter, prefilterLimit_);
	}
	
	/**
	 * Creates a table with minimum quality zero, unlimited input and one counting pass.
	 * Does not reset other inherited options such as minProb or duplicate policy.
	 *
	 * @param fname1 Primary input file
	 * @param fname2 Secondary input file (for paired reads)
	 * @param extraFiles Additional input files to process
	 * @param cbits Bits per count cell in the array
	 * @param cells Total number of cells in the count array
	 * @param hashes Number of hash functions to use
	 * @param prefilter Optional prefilter for two-stage counting
	 * @param prefilterLimit Prefilter threshold limit
	 * @return Populated KCountArray with k-mer counts
	 */
	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes,
			KCountArray prefilter, int prefilterLimit){
		return makeKca(fname1, fname2, extraFiles, cbits, cells, hashes, 0, -1, 1, 1, 1, prefilter, prefilterLimit);
	}
	
	/**
	 * Creates a KCountArray by processing input files with specified parameters.
	 * Handles comma-separated file lists, multi-pass counting and trusted-base masking.
	 * A null primary path is skipped in every pass; extra inputs may supply all reads.
	 * Saves/restores inherited minQuality, maxReads and FASTA split mode on every exit.
	 * Other global options/counters remain shared. Exceptions reject incomplete results;
	 * returned arrays are shut down after loading, except the all-null empty-table path.
	 * Extra inputs are read completely in every pass, independent of the primary limit.
	 *
	 * @param fname1 Primary input file
	 * @param fname2 Secondary input file (for paired reads)
	 * @param extraFiles Additional input files to process
	 * @param cbits Bits per count cell in the array
	 * @param cells Total number of cells in the count array
	 * @param hashes Number of hash functions to use
	 * @param minqual Minimum quality score for k-mer inclusion
	 * @param maxreads Reader-dependent primary limit (-1 for unlimited): FASTQ roots/pairs
	 * or SAM alignment records
	 * @param passes Number of processing passes (1 for single-pass, >1 for multi-pass)
	 * @param thresh1 Threshold for early pass filtering
	 * @param thresh2 Threshold for final pass filtering
	 * @param prefilter Optional prefilter for two-stage counting
	 * @param prefilterLimit_ Prefilter threshold limit
	 * @return Populated KCountArray with k-mer counts
	 */
	public KCountArray makeKca(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int thresh1, int thresh2,
			KCountArray prefilter, int prefilterLimit_){
		final boolean oldsplit=FastaReadInputStream.SPLIT_READS;
		final long oldmax=maxReads;
		final byte oldq=minQuality;
		try{
			return makeKcaInner(fname1, fname2, extraFiles, cbits, cells, hashes,
					minqual, maxreads, passes, thresh1, thresh2, prefilter, prefilterLimit_);
		}catch(Exception e){
			throw new RuntimeException("Kmer counting failed; refusing to return incomplete counts for "+fname1, e);
		}finally{
			minQuality=oldq;
			maxReads=oldmax;
			FastaReadInputStream.SPLIT_READS=oldsplit;
		}
	}

	/** Builds successive tables; each completed table becomes the next pass's trusted set.
	 * Extra inputs require a repeatable Iterable. The public wrapper restores its three
	 * saved counting settings even if loading or finalization fails. */
	private KCountArray makeKcaInner(String fname1, String fname2, Iterable<String> extraFiles,
			int cbits, long cells, int hashes, int minqual,
			long maxreads, int passes, int thresh1, int thresh2,
			KCountArray prefilter, int prefilterLimit_) throws Exception{
		assert(passes>0) : "At least one complete input pass is required before returning kmer counts: "+passes;
		if(verbose){System.err.println("Making kca from ("+fname1+", "+fname2+")\nk="+k+", cells="+Tools.toKMG(cells)+", cbits="+cbits);}
		
		if(fname1==null && fname2==null && extraFiles==null){
			return KCountArray.makeNew(cells, cbits, hashes, prefilter, prefilterLimit_);
		}
		
		maxReads=maxreads;
		minQuality=(byte)minqual;
		KCountArray kca=KCountArray.makeNew(cells, cbits, hashes, prefilter, prefilterLimit_);
		
		{//For processing input lists
			ArrayList<String> extra2=null;
			if(fname1!=null && fname1.contains(",")){
				String[] s=fname1.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname1=s[0];
			}
			if(fname2!=null && fname2.contains(",")){
				String[] s=fname2.split(",");
				if(extra2==null){extra2=new ArrayList<String>();}
				for(int i=1; i<s.length; i++){extra2.add(s[i]);}
				fname2=s[0];
			}
			if(extra2!=null){
				if(extraFiles!=null){
					for(String s : extraFiles){
						extra2.add(s);
					}
				}
				extraFiles=extra2;
			}
		}
		
		if(extraFiles!=null){
			for(String s : extraFiles){
				if(fileIO.FileFormat.hasFastaExtension(s)){
					assert(!FastaReadInputStream.SPLIT_READS);
				}
			}
		}
		
		if(passes==1){
			if(fname1!=null){
				count(fname1, fname2, kca);
			}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					count(s, null, kca);
				}
			}
			kca.shutdown();

		}else{
			assert(passes>1);
			KCountArray trusted=null;
			for(int i=1; i<passes; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
//				int step=(stepsize==1 ? 1 : stepsize+i%2);
//				//			if(!conservative){step=(step+3)/4;}
//				if(!conservative){step=Tools.min(3, (step+3)/4);}

				if(fname1!=null){count(fname1, fname2, cbits, kca, trusted, maxreads, thresh1, conservative);}
				if(extraFiles!=null){
					maxReads=-1;
					for(String s : extraFiles){
						count(s, null, cbits, kca, trusted, -1, thresh1, conservative);
					}
				}
				kca.shutdown();
				
				System.out.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(cells, cbits, hashes, prefilter, prefilterLimit_);

			}

			if(fname1!=null){count(fname1, fname2, cbits, kca, trusted, maxreads, thresh2, true);}
			if(extraFiles!=null){
				maxReads=-1;
				for(String s : extraFiles){
					count(s, null, cbits, kca, trusted, -1, thresh2, true);
				}
			}
			kca.shutdown();
		}
		return kca;
	}
	
	/**
	 * Opens/starts an owned Streamer, runs counting workers and closes it after joins.
	 * Uses the inherited reader-dependent limit (FASTQ roots/pairs, SAM alignment records)
	 * and at most MAX_COUNT_THREADS workers. Startup failures lack finally cleanup.
	 * An absent primary path containing # expands to mate1/2 paths when reads2 is null.
	 * Joins continue after interruption; the interrupted flag is not restored here.
	 * The borrowed destination is mutated and not shut down. On failure discard its
	 * partial contents; no successful return or rollback is promised.
	 *
	 * @param reads1 Non-null primary input file path
	 * @param reads2 Secondary input file path (null for single-end)
	 * @param counts Borrowed destination supporting concurrent increments when multiple
	 * workers run; shared counting settings must remain stable throughout the call
	 * @return The updated KCountArray
	 * @throws Exception If file I/O or processing errors occur
	 */
	public KCountArray count(String reads1, String reads2, KCountArray counts) throws Exception{
		assert(counts!=null);
		
		{
			int pound=reads1.lastIndexOf('#');
			if(pound>=0 && reads2==null && !new File(reads1).exists()){
				String a=reads1.substring(0, pound);
				String b=reads1.substring(pound+1);
				reads1=a+1+b;
				reads2=a+2+b;
			}
		}
		
		final Streamer cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			if(ff2==null){ff1.preferShreds=true;}
//			if(ff2!=null){ //TODO - interleaved flag
			cris=StreamerFactory.getReadInputStream(maxReads, true, ff1, ff2, 1);
			cris.start();
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}
		
		CountThread[] cta=new CountThread[Tools.min(Shared.threads(), MAX_COUNT_THREADS)];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, counts);
			cta[i].start();
		}
		boolean errorState=false;
		for(int i=0; i<cta.length; i++){
			CountThread ct=cta[i];
			synchronized(ct){
				while(ct.getState()!=State.TERMINATED){
					try{
						ct.join(2000);
					}catch(InterruptedException e){
						//Interruption is logged; joining continues without restoring interrupt status.
						e.printStackTrace();
					}
				}
			}
			errorState|=ct.errorState;
		}

		errorState|=ReadWrite.closeStream(cris);
		if(verbose){System.err.println("Closed stream");}
		if(verbose){System.err.println("Processed "+readsProcessed+" reads.");}

		if(errorState){
			throw new RuntimeException("Input or CountThread failure while counting "+reads1+"; kmer counts are incomplete.");
		}
		return counts;
	}
	
	/**
	 * Opens an owned reader and counts after masking untrusted bases on temporary reads.
	 * Uses the explicit reader-dependent limit and Shared.threads workers without the ordinary
	 * overload's MAX_COUNT_THREADS cap. Reader closure follows worker joins; the borrowed
	 * destination remains open for more increments. A failure leaves partial counts.
	 * Startup failures lack finally cleanup; interrupted joins continue without restoring
	 * interrupt status. Shared counting settings must remain stable throughout the call.
	 *
	 * @param reads1 Primary input file path
	 * @param reads2 Secondary input file path (null for single-end)
	 * @param cbits Bits per count cell (unused in this context)
	 * @param counts Borrowed destination supporting concurrent increments when multiple workers run
	 * @param trusted Optional completed trusted table, stable and safe for concurrent reads
	 * @param maxReads Reader-dependent input limit: FASTQ roots/pairs or SAM alignment records;
	 * negative means unlimited
	 * @param thresh Minimum count threshold in trusted set
	 * @param conservative Whether to use conservative error detection
	 * @return The updated KCountArray
	 * @throws Exception If file I/O or processing errors occur
	 */
	public KCountArray count(String reads1, String reads2, final int cbits, 
			KCountArray counts, final KCountArray trusted, final long maxReads,
			final int thresh, final boolean conservative) throws Exception{
		
		assert(counts!=null);
		
		{
			int pound=reads1.lastIndexOf('#');
			if(pound>=0 && reads2==null && !new File(reads1).exists()){
				String a=reads1.substring(0, pound);
				String b=reads1.substring(pound+1);
				reads1=a+1+b;
				reads2=a+2+b;
			}
		}
		
		final Streamer cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			if(ff2==null){ff1.preferShreds=true;}
			cris=StreamerFactory.getReadInputStream(maxReads, true, ff1, ff2, 1);
			cris.start();
		}
		
		assert(cris!=null) : reads1;
		if(verbose){System.err.println("Started cris");}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Paired: "+paired);}

		CountThread[] cta=new CountThread[Shared.threads()];
		for(int i=0; i<cta.length; i++){
			cta[i]=new CountThread(cris, counts, trusted, thresh, conservative);
			cta[i].start();
		}
		
		boolean errorState=false;
		for(int i=0; i<cta.length; i++){
			CountThread ct=cta[i];
			synchronized(ct){
				while(ct.isAlive()){
					try{
						ct.join(1000);
					}catch(InterruptedException e){
						//Interruption is logged; joining continues without restoring interrupt status.
						e.printStackTrace();
					}
				}
			}
			errorState|=ct.errorState;
		}

		errorState|=ReadWrite.closeStream(cris);
		if(verbose){System.err.println("Closed stream");}

		if(errorState){
			throw new RuntimeException("Input or CountThread failure while counting "+reads1+"; kmer counts are incomplete.");
		}
		
		return counts;
	}
	
	/**
	 * Delegates paired-overlap analysis to the selected strict or very-strict BBMerge path.
	 * With ecc enabled, the overlap routine may modify the paired reads.
	 *
	 * @param r1 First read of the pair
	 * @param r2 Second read of the pair
	 * @param ecc Whether this is for error correction (affects algorithm choice)
	 * @return Positive insert length on a suitable match, otherwise a nonpositive value
	 */
	private final int findOverlap(Read r1, Read r2, boolean ecc){
		if(vstrict){
			return BBMerge.findOverlapVStrict(r1, r2, ecc);
		}else{
			return BBMerge.findOverlapStrict(r1, r2, ecc);
		}
	}
	
	/**
	 * Packs a selected interval using Dedupe's permissive two-bit alphabet.
	 * A/C/G/T/U and lowercase forms have explicit codes; other ASCII symbols map to
	 * zero rather than an invalid-window sentinel. More than 32 bases lose high bits.
	 * @param bases Non-null ASCII base array, read without modification
	 * @param start Nonnegative start index
	 * @param klen Number of characters within the array bounds
	 * @return Low 64 bits of the packed interval, without canonicalization or quality filtering
	 */
	public static final long toKmer(final byte[] bases, final int start, final int klen){
		final int stop=start+klen;
		assert(stop<=bases.length) : klen+", "+bases.length;
		long kmer=0;
		
		for(int i=start; i<stop; i++){
			final byte b=bases[i];
			final long x=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|x);
		}
		return kmer;
	}
	
	/** Owns scratch buffers while borrowing the shared reader and destination/trusted tables. */
	private class CountThread extends Thread{
		
		/** Selects counting without a trusted table. */
		CountThread(final Streamer cris_, final KCountArray counts_){
			this(cris_, counts_, null, 2, true);
		}
		
		/** Captures shared tables/reader and this worker's trusted-mask policy; does not start. */
		CountThread(final Streamer cris_,
				final KCountArray counts_, final KCountArray trusted_, final int thresh_,
				final boolean conservative_){
			cris=cris_;
			counts=counts_;
			trusted=trusted_;
			thresh=thresh_;
			conservative=conservative_;
		}
		
		/** Counts, flushes the final buffer, records failures and aggregates observed totals.
		 * Totals may be partial on failure and do not certify a successful count table. */
		@Override
		public void run(){
			//2026-08-21 (Nepgear+Amber): count(cris) used to be able to throw an uncaught
			//AssertionError from PairStreamer (a genuine R1/R2 desync bug, now fixed at its
			//source -- see PairStreamer's carry-over buffer) and silently drop this thread's
			//contribution to the shared totals with no error visible to the caller. Same
			//try/catch+error-flag shape as BloomFilterCorrectorWrapper$ProcessThread and
			//FastqScanStreamer$ReaderRunnable: catch Throwable (assertions are Errors, not
			//Exceptions), record it, and let makeKca/count's caller-side check turn it into a
			//loud, reported failure instead of quietly-wrong counts.
			try{
				count(cris);
				//The final partial buffer contributes counts too; its failure must reach the caller.
				if(BUFFERED){dumpBufferT();}
			}catch(Throwable e){
				e.printStackTrace();
				errorState=true;
			}

			synchronized(getClass()){
				keysCounted+=keysCountedLocal;
				increments+=incrementsLocal;
				readsProcessed+=readsProcessedLocal;

				if(verbose){System.err.println(keysCounted+", "+keysCountedLocal);}
				if(verbose){System.err.println(readsProcessed+", "+readsProcessedLocal);}
			}
		}

		/** Set if count(cris) threw; checked by the caller after join() so a failure is reported, not silently dropped. */
		boolean errorState=false;
		
		/** Buffers one selected key or immediately applies its update to the table. */
		private void increment(final long key){
			if(BUFFERED){
				buffer.add(key);
				if(buffer.size>=BUFFERLEN){
					dumpBufferT();
				}
			}else{
				incrementByAmountT(key, 1);
			}
		}

		/** Applies a two-hash update directly; sketch mode instead uses the single-key path. */
		private void increment(final long key, final long key2){
			if(SKETCH_MODE){
				incrementByAmountT(key, 1);
			}else{
				counts.incrementDual(key, key2);
				keysCountedLocal++;
				incrementsLocal++;
			}
		}
		
		/** Samples a forward interval chosen from numericID and knum, without ordinary
		 * quality/canonicalization rules. Matches PlotFlowCell's raw sampled-key convention. */
		private void incrementOneKmer(Read r, int knum){
			if(r==null || r.length()<k){return;}
			final long kmer=toKmer(r.bases, (int)((r.numericID+knum*k)%(r.length()-k2)), k);
			if(kmer>=0){increment(kmer);}
		}
		
		/** Sorts the worker buffer, combines equal keys, emits updates and clears the buffer. */
		private void dumpBufferT(){
			final int lim=buffer.size;
			if(lim<1){return;}
			final long[] array=buffer.array;
			buffer.sort();//Can be disabled via parallelsort flag but only affects arrays>10k
			long kmer=array[0]-1;//Ensures a nonmatch
			int streak=0;
			for(int i=0; i<lim; i++){
				long x=array[i];
				if(x==kmer){streak++;}else{
					if(streak>0){incrementByAmountT(kmer, streak);}
					kmer=x;
					streak=1;
				}
			}
			assert(streak>0);
			incrementByAmountT(kmer, streak);
			buffer.clear();
		}
		
		/** Applies an aggregate single-key update, honoring sketch selection/hash storage.
		 * Counts accepted multiplicity in keysCountedLocal and one call in incrementsLocal. */
		private void incrementByAmountT(final long key, int amt){
			if(SKETCH_MODE){
				final long code=SketchObject.hash(key);
				if(code<sketchMinHashValue){return;}
				counts.increment(STORE_HASHED ? code : key, amt);
			}else if(dualHash){
				counts.increment(key, KCountArray.secondaryHash(key), amt);
			}else{
				counts.increment(key, amt);
			}
			keysCountedLocal+=amt;
			incrementsLocal++;
		}
		
		/** Consumes batches, optionally merges/corrects/masks temporary reads, and selects
		 * amino, short-key or long-key counting. The parent owns reader shutdown. */
		private final void count(final Streamer cris){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			final LongList buffer=new LongList(300);
			final Kmer kmer=(k>maxShortKmerLength ? new Kmer(k) : null);
			
			while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning
				for(Read r1 : reads){
					readsProcessedLocal+=r1.pairCount();
					
					Read r2=r1.mate;
					if((ecco || merge) && r2!=null){
						if(merge){
							final int insert=findOverlap(r1, r2, false);
							if(insert>0){
								//r2 is rcomp'd and not restored (unlike BloomFilterWrapper) - safe HERE: this path
								//only COUNTS the merged r1 and discards r2 (no output), so r2's orientation is moot.
								r2.reverseComplementFast();
								r1=r1.joinRead(insert);
								r2=null;
							}
						}else if(ecco){
							findOverlap(r1, r2, true);
						}
					}
					
					if(trusted!=null){
						clearUntrustedBases(r1);
						clearUntrustedBases(r2);
					}

					if(amino){
						addReadAmino(r1);
					}else if(k<=maxShortKmerLength){
						addRead_Advanced(r1, buffer);
					}else{
						addReadBig(r1, kmer);
						addReadBig(r1.mate, kmer);
					}
				}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			if(verbose){System.err.println("Finished reading");}
		}
		
		/** Replaces uncovered bases with N and existing qualities with zero; null masks leave
		 * the read unchanged. Returns the number of masked positions, not verified errors. */
		private int clearUntrustedBases(Read r){
			if(r==null || trusted==null){return 0;}
			BitSet bs=(conservative ? ErrorCorrect.detectErrorsBulk(r, trusted, k, thresh, detectStepsize) :
				ErrorCorrect.detectTrusted(r, trusted, k, thresh, detectStepsize));
			if(bs==null){return 0;}
			int ns=0;
			for(int i=bs.nextClearBit(0); i<r.length(); i=bs.nextClearBit(i+1)){
				r.bases[i]='N';
				if(r.quality!=null){r.quality[i]=0;}
				ns++;
			}
			return ns;
		}
		
		/**
		 * Hash a read's kmers into the KCountArray.
		 * Chooses sampled, per-pair deduplicated or all-window short-key counting.
		 * PREJOIN may reverse the mate in place and replace the local root with a join;
		 * these reads are discarded after counting. Only duplicate suppression sorts keys.
		 * @param r1 Temporary root with an optional linked mate
		 * @param buffer Reusable pair-local key scratch, cleared on entry
		 */
		private final void addRead_Advanced(Read r1, final LongList buffer){
			assert(k<=maxShortKmerLength);
			assert(!amino);
			
			buffer.clear();
			if(PREJOIN && r1.mate!=null && r1.insert()>0){
				r1.mate.reverseComplementFast();
				r1=r1.joinRead();
				assert(r1.mate==null);
			}
			Read r2=r1.mate;
			final int len1=Tools.max(0, r1.length()-k+1);
			final int len2=(r2==null || r2.bases==null) ? 0 : Tools.max(0, r2.length()-k+1);
			
			if(KMERS_PER_READ>0){
				final int lim1=Tools.min(KMERS_PER_READ, r1.length()-k2);
				final int lim2=Tools.min(KMERS_PER_READ, r1.mateLength()-k2);
				for(int i=0; i<lim1; i++){incrementOneKmer(r1, i);}
				for(int i=0; i<lim2; i++){incrementOneKmer(r2, i);}
			}else if(!KEEP_DUPLICATE_KMERS){
				fillKmerArray(r1, buffer);
				if(r2!=null){fillKmerArray(r2, buffer);}
				if(buffer.size()<1){return;}
				
				buffer.sort();
				long prev=-1;
				final long[] array=buffer.array;
				final int idmod=(int)(r1.numericID&IDMASK);
				for(int i=0, lim=buffer.size(); i<lim; i++){
					long kmer=array[i];
					if(kmer!=prev && (i&IDMASK)==idmod){
						increment(kmer);
						prev=kmer;
					}
				}
			}else{
				if(len1>0){addRead(r1);}
				if(len2>0){addRead(r2);}
			}
		}
		
		/** Counts multiword nucleotide windows with the table's one-/two-hash mode.
		 * Uses minProb when qualities exist; short-path sampling, minQuality, duplicate
		 * suppression and per-instance rcomp selection are not implemented here.
		 * A qualifying quality suffix prevents permanent product underflow. Double
		 * comparisons can differ from the historical float recurrence near the cutoff. */
		private final void addReadBig(Read r, Kmer kmer){
			if(r==null || r.bases==null){return;}
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			int len=0;
			
			if(bases==null || bases.length<k){return;}
			kmer.clear();
			
			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			final double cutoff=minProb;
			final boolean useQuality=(cutoff>0 && quals!=null);
			double prob=1;
			int qualityStart=0;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];

				//Update kmers
				kmer.addRight(b);
				
				if(useQuality){
					if(i-qualityStart>=k){prob/=KmerTableSet.PROB_CORRECT[quals[qualityStart++]];}
					prob*=KmerTableSet.PROB_CORRECT[quals[i]];
					//Factors are in (0,1]: adding bases cannot rescue a prefix already below cutoff.
					//Each index enters/leaves once. Keeping a qualifying suffix also bounds the product
					//above positive-float cutoff times the smallest quality factor, safely within double.
					while(prob<cutoff && qualityStart<=i){
						prob/=KmerTableSet.PROB_CORRECT[quals[qualityStart++]];
					}
					if(qualityStart>i){prob=1;}
				}

				//Handle Ns
				if(x<0){
					len=0;
					prob=1;
					qualityStart=i+1;
				}else{len++;}
				
				assert(len==kmer.len);
				
				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
				if(len>=k && prob>=cutoff && (!useQuality || i-qualityStart+1>=k)){
					if(dualHash){increment(kmer.xor(), kmer.xor2());}else{increment(kmer.xor());}
				}
			}
		}
		
		/** Appends eligible short keys for pair deduplication without clearing the list.
		 * Applies quality/probability and strand selection; long/amino requests delegate. */
		private final void fillKmerArray(Read r, final LongList list){
			if(amino){
				fillKmerArrayAmino(r, list);
				return;
			}
			if(k>maxShortKmerLength){
				fillKmerArrayLong(r, list);
				return;
			}
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(list!=null);
			
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			
			if(bases==null || bases.length<k){return;}
			
			long kmer=0;
			long rkmer=0;
			int len=0;
			float prob=1;
			
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				final byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(len>=k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}

				if(x<0 || q<minQuality){
					len=0;
					kmer=rkmer=0;
					prob=1;
				}else{
					len++;
					if(len>=k && prob>=minProb){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						list.add(key);
					}
				}
			}
		}
		
		/** Unimplemented amino-key collection: asserts under -ea and does nothing under -da. */
		private final void fillKmerArrayAmino(Read r, final LongList array){
			assert(false) : "TODO"; //TODO
		}
		
		/** Counts all eligible short windows, applying quality/probability, strand and IDMASK rules. */
		private final void addRead(Read r){
			if(amino){
				addReadAmino(r);
				return;
			}
			assert(k<=maxShortKmerLength);
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);

			final byte[] bases=r.bases;
			final byte[] quals=r.quality;

			if(bases==null || bases.length<k){return;}
			
			long kmer=0;
			long rkmer=0;
			int len=0;
			float prob=1;
			final int idmod=(int)(r.numericID&IDMASK);
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				final byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(len>=k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}

				if(x<0 || q<minQuality){
					len=0;
					kmer=rkmer=0;
					prob=1;
				}else{
					len++;
					if(len>=k && prob>=minProb && (i&IDMASK)==idmod){
						long key=(rcomp ? Tools.max(kmer, rkmer) : kmer);
						increment(key);
					}
				}
			}
		}
		
		/** Counts valid packed amino windows in this read only, ignoring quality and mate fields. */
		private final void addReadAmino(Read r){
			if(r==null){return;}
			assert(!PREJOIN || r.mate==null);

			final byte[] bases=r.bases;

			if(bases==null || bases.length<k){return;}
			
			long kmer=0;
			int len=0;

			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				long x=AminoAcid.acidToNumber[b];
				kmer=((kmer<<aminoShift)|x)&mask;

				if(x<0){
					len=0;
					kmer=0;
				}else{
					len++;
					if(len>=k){
						increment(kmer);
					}
				}
			}
		}
		
		/** Appends canonical multiword hashes with quality/probability filtering; creates
		 * one local Kmer object per call. This helper does not clear the destination. */
		private final void fillKmerArrayLong(Read r, final LongList list){
			assert(k>maxShortKmerLength) : k;
			assert(!PREJOIN || r.mate==null);
			assert(CANONICAL);
			assert(list!=null);
			Kmer kmer=new Kmer(k);
			
			float prob=1;
			byte[] bases=r.bases;
			byte[] quals=r.quality;
			
			kmer.clear();
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				kmer.addRight(b);
				
				byte q;
				if(quals==null){
					q=50;
				}else{
					q=quals[i];
					prob=prob*align2.QualityTools.PROB_CORRECT[q];
					if(kmer.len>k){
						byte oldq=quals[i-k];
						prob=prob*align2.QualityTools.PROB_CORRECT_INVERSE[oldq];
					}
				}
				
				if(!AminoAcid.isFullyDefined(b) || q<minQuality){
					kmer.clear();
					prob=1;
				}
				if(kmer.len>=k && prob>=minProb){
					list.add(kmer.xor());
				}
			}
		}
		
		/** Accepted update multiplicities, accumulated even if this worker later fails. */
		private long keysCountedLocal=0;
		/** Number of issued table-update calls after buffering/sketch selection. */
		private long incrementsLocal=0;
		/** Individual input reads, counting both mates before optional merging. */
		private long readsProcessedLocal=0;
		
		/** Shared started reader; this worker consumes but does not close it. */
		private final Streamer cris;
		
		/** Shared destination, mutated through its supported concurrent increment API. */
		private final KCountArray counts;
		/** Optional completed table used for trusted-base masks. */
		private final KCountArray trusted;
		/** Minimum trusted-table support. */
		private final int thresh;
		/** Selects bulk error detection instead of sampled trusted-window detection. */
		private final boolean conservative;
		/** Sketch hash cutoff captured when this worker is constructed. */
		private final long sketchMinHashValue=SketchObject.minHashValue;
		
		/** Optional single-key update buffer; BUFFERED must remain stable during this worker's life. */
		private final LongList buffer=(BUFFERED ? new LongList(BUFFERLEN) : null);
	}
	
	/**
	 * Number of bases to skip when detecting errors during trusted k-mer processing
	 */
	public int detectStepsize=1;//Jump this many bases when detecting errors
	
	/** K-mer length used for counting */
	private final int k;
	/** K-mer length minus 1, used for boundary calculations */
	private final int k2;
	/** Bit shift amount for amino acid encoding */
	private final int aminoShift;
	/** Total bit shift for k-mer encoding */
	private final int shift;
	/** Bit shift for reverse complement k-mer construction */
	private final int shift2;
	/** Bit mask for extracting k-mer bits from rolling hash */
	private final long mask;
	/** Whether to use reverse complement canonical k-mers */
	private final boolean rcomp;
	/** Whether error correction is enabled */
	private final boolean ecco;
	/** Whether to merge overlapping paired reads */
	private final boolean merge;
	/** Whether processing amino acid sequences instead of nucleotides */
	private final boolean amino;
	/** Whether kmers use independent secondary hashes for later Bloom lanes. */
	private final boolean dualHash;
	
	/**
	 * Whether to use very strict overlap detection instead of regular strict mode
	 */
	public static boolean vstrict=false;
	
}
