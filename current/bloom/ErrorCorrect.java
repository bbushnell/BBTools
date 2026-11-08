package bloom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Legacy packed-kmer error detection and substitution correction.
 * The correction loop runs on the calling thread; input, output and count-array
 * implementations may use their own workers. This is separate from BBCMS's
 * BloomFilterCorrector pipeline. KmerCount5/6 and ReadCounter also use the static
 * coverage detectors without invoking this command-line driver.
 *
 * Detection queries canonical keys and returns either window-start or base-index
 * masks, as documented per method. Correction edits the supplied Read and uses
 * forward-only candidate queries, so callers must supply compatible counts.
 * Shared flags and KmerCount6 settings are not restored between invocations;
 * concurrent pipelines with differing settings are not isolated. BLM205 records
 * the bounded owned-stream cleanup repair and its remaining coverage limits.
 *
 * @author Brian Bushnell
 * @date Aug 20, 2012
 */
public class ErrorCorrect extends Thread{
	
	/**
	 * Builds counts and processes reads using the legacy positional command line.
	 * The first two entries are input paths; option scanning starts at index 2,
	 * so a flag in the second position is treated as a filename. PreParser handles
	 * help/configuration first. Changes shared parsing and correction settings.
	 * @param args First input path, optional second input path, then flag=value options
	 */
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		String reads1=args[0];
		String reads2=(args.length>1 ? args[1] : null);
		
		int k=23;
		int cbits=4;
		int hashes=1;
		int thresh1=1;
		int thresh2=2;
		int matrixbits=34;
		long maxReads=-1;
		int buildpasses=1;
		long tablereads=-1; //How many reads to process when building the hashtable
		int buildStepsize=4;
		String output=null;
		boolean ordered=true;
		boolean overwrite=true;
		boolean append=false;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(a.startsWith("cbits") || a.startsWith("cellbits")){
				cbits=Integer.parseInt(b);
			}else if(a.equals("initialthresh") || a.equals("thresh1")){
				thresh1=Integer.parseInt(b);
			}else if(a.equals("thresh") || a.equals("thresh2")){
				thresh2=Integer.parseInt(b);
			}else if(a.startsWith("matrixbits")){
				matrixbits=Integer.parseInt(b);
			}else if(a.startsWith("hashes") || a.startsWith("multihash")){
				hashes=Integer.parseInt(b);
			}else if(a.startsWith("maxerrors")){
				ERROR_CORRECTION_LIMIT=Integer.parseInt(b);
			}else if(a.startsWith("passes")){
				buildpasses=Integer.parseInt(b);
			}else if(a.startsWith("stepsize") || a.startsWith("buildstepsize")){
				buildStepsize=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				System.err.println("Can't change threadcount for this class."); //THREADS=Integer.parseInt(b);
			}else if(a.startsWith("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.startsWith("tablereads")){
				tablereads=Parse.parseKMG(b);
			}else if(a.startsWith("build") || a.startsWith("genome")){
				Data.setGenome(Integer.parseInt(b));
				Data.sysout.println("Set genome to "+Data.GENOME_BUILD);
			}else if(a.equals("outputinfo") || a.startsWith("info")){
				OUTPUT_INFO=Parse.parseBoolean(b);
			}else if(a.startsWith("out")){
				output=b;
			}else if(a.startsWith("verbose")){
				KCountArray.verbose=Parse.parseBoolean(b);
//				verbose=KCountArray.verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered") || a.equals("ord")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		KCountArray kca=makeTable(reads1, reads2, k, cbits, hashes, buildpasses, matrixbits, tablereads, buildStepsize, thresh1, thresh2);
		
		detect(reads1, reads2, kca, k, thresh2, maxReads, output, ordered, append, overwrite);
	}
	
	/**
	 * Builds a count table, optionally using earlier passes as trusted-count filters.
	 * Calls shutdown on each table after its normal counting pass, then returns
	 * the final table. Writes KmerCount6.maxReads and does not restore that global.
	 * Intermediate passes adjust sampling; the final pass uses stepsize directly.
	 *
	 * @param reads1 First input reads file
	 * @param reads2 Second input reads file (may be null)
	 * @param k Packed DNA k-mer length in bases, from 1 through 31
	 * @param cbits Bits per count cell in bloom filter
	 * @param hashes Number of hash functions in bloom filter
	 * @param buildpasses Positive number of counting passes
	 * @param matrixbits Base-two logarithm of requested cell count, capped at 2*k
	 * @param maxreads Input-stream limit for table building; -1 means unlimited
	 * @param stepsize Positive detection-sampling step for the final trusted-table pass
	 * @param thresh1 Trusted-count threshold for intermediate passes
	 * @param thresh2 Trusted-count threshold for the final pass
	 * @return Final populated array after shutdown, owned by the caller
	 */
	public static KCountArray makeTable(String reads1, String reads2, int k, int cbits, int hashes, int buildpasses, int matrixbits,
			long maxreads, int stepsize, int thresh1, int thresh2){
		
		Timer thash=new Timer();
		
		KmerCount6.maxReads=maxreads;
		int kbits=2*k;
		matrixbits=Tools.min(kbits, matrixbits);
		
		thash.start();
//		Data.sysout.println("kbits="+(kbits)+" -> "+(1L<<kbits)+", matrixbits="+(matrixbits)+" -> "+(1L<<matrixbits)+", cbits="+cbits+", gap="+gap+", hashes="+hashes);
		KCountArray kca=KCountArray.makeNew(1L<<matrixbits, cbits, hashes);
		
		if(buildpasses==1){
			
			KmerCount6.count(reads1, reads2, k, cbits, true, kca);
			kca.shutdown();
			
		}else{
			assert(buildpasses>1);
			KCountArray trusted=null;
			for(int i=1; i<buildpasses; i++){
				boolean conservative=i>2;// /*or, alternately, (trusted==null || trusted.capacity()>0.3)
				int step=(stepsize==1 ? 1 : stepsize+i%2);
//				if(!conservative){step=(step+3)/4;}
				if(!conservative){step=Tools.min(3, (step+3)/4);}
				
				KmerCount6.count(reads1, reads2, k, cbits, true, kca, trusted, maxreads, thresh1, step, conservative);
				
				kca.shutdown();
				Data.sysout.println("Trusted:   \t"+kca.toShortString());
				trusted=kca;
				kca=KCountArray.makeNew(1L<<matrixbits, cbits, hashes);
				
			}
			
			KmerCount6.count(reads1, reads2, k, cbits, true, kca, trusted, maxreads, thresh2, stepsize, true);
			
			kca.shutdown();
		}
		
		
		thash.stop();
//		Data.sysout.println(kca.toString());
		Data.sysout.println("Table :    \t"+kca.toShortString());
		Data.sysout.println("Hash time:      \t"+thash);
		return kca;
	}
	
	/**
	 * Opens and starts owned read streams, then runs detection/correction on the calling thread.
	 * With paired input, a # in output selects separate mate files; otherwise
	 * output is interleaved. No output stream is created for null output. Attempts
	 * closure of each acquired stream in finally, including construction/processing failures.
	 * A recorded stream error prevents successful return; cleanup failures are suppressed
	 * on an existing primary exception rather than replacing it.
	 * Failure to start an acquired output halts loudly: its API cannot close NEW backends,
	 * and joining them would spin forever. This fatal-start path does not return to the caller.
	 *
	 * @param reads1 First input reads file
	 * @param reads2 Second input reads file (may be null)
	 * @param kca Borrowed count table compatible with detection and correction queries
	 * @param k Packed DNA k-mer length in bases
	 * @param thresh Minimum k-mer count threshold
	 * @param maxReads Limit passed to the input stream; -1 means unlimited
	 * @param output Output file path (may be null)
	 * @param ordered Whether to maintain read order in output
	 * @param append Whether to append to existing output file
	 * @param overwrite Whether to overwrite existing output file
	 */
	public static void detect(String reads1, String reads2, KCountArray kca, int k, int thresh, long maxReads, String output, boolean ordered, boolean append, boolean overwrite){

		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(reads1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(reads2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
		}
		assert(cris!=null) : reads1;
		
		ConcurrentReadOutputStream ros=null;
		boolean outputStarted=false;
		Throwable primaryFailure=null;
		try{
			cris.start();
			if(verbose){System.err.println("Started cris");}
			boolean paired=cris.paired();
			if(verbose){System.err.println("Paired: "+paired);}
			if(output!=null){
				String out1=output.replaceFirst("#", "1"), out2=null;
				if(cris.paired()){
					if(output.contains("#")){
						out2=output.replaceFirst("#", "2");
					}else{
						System.err.println("Writing interleaved.");
					}
				}
				final int buff=(!ordered ? 8 : Tools.max(16, 2*THREADS));
				FileFormat ff1=FileFormat.testOutput(out1, FileFormat.FASTQ, OUTPUT_INFO ? ".info" : null, true, overwrite, append, ordered);
				FileFormat ff2=FileFormat.testOutput(out2, FileFormat.FASTQ, OUTPUT_INFO ? ".info" : null, true, overwrite, append, ordered);
				assert(!ff1.sam()) : "Sam files need reference info for the header.";
				ros=ConcurrentReadOutputStream.getStream(ff1, ff2, buff, null, true);
			}
			if(ros!=null){
				ros.start();
				outputStarted=true;
				Data.sysout.println("Started output threads.");
			}
			detect(cris, kca, k, thresh, maxReads, ros);
		}catch(RuntimeException | Error failure){
			primaryFailure=failure;
			throw failure;
		}finally{
			if(ros!=null && !outputStarted){
				//A partial start can leave NEW backends: abort does not close them and join spins.
				if(primaryFailure!=null){primaryFailure.printStackTrace(System.err);}
				KillSwitch.kill("ErrorCorrect output failed to start completely; cannot safely join unstarted writers. Output is incomplete.");
			}
			//BLM205: close independently so failure from one owner cannot skip the other.
			boolean streamError=false;
			Throwable cleanupFailure=null;
			try{streamError|=ReadWrite.closeStream(cris);}
			catch(RuntimeException | Error failure){cleanupFailure=failure;}
			try{streamError|=ReadWrite.closeStream(ros);}
			catch(RuntimeException | Error failure){
				if(cleanupFailure==null){cleanupFailure=failure;}
				else if(cleanupFailure!=failure){cleanupFailure.addSuppressed(failure);}
			}
			if(streamError){
				final RuntimeException failure=new RuntimeException("ErrorCorrect input/output stream failed; output may be incomplete.");
				if(cleanupFailure==null){cleanupFailure=failure;}
				else{cleanupFailure.addSuppressed(failure);}
			}
			if(cleanupFailure!=null){
				if(primaryFailure!=null){
					if(primaryFailure!=cleanupFailure){primaryFailure.addSuppressed(cleanupFailure);}
				}else if(cleanupFailure instanceof RuntimeException){throw (RuntimeException)cleanupFailure;}
				else{throw (Error)cleanupFailure;}
			}
		}
		if(verbose){System.err.println("Closed stream");}
	}
	
	/**
	 * Consumes started input batches, corrects reads in place and optionally submits them to output.
	 * Processes both mates and reports individual-read/base totals. Returns batches
	 * to cris, but does not start or close either borrowed stream. Every processed batch ID
	 * is sent to ros when present, including lists whose roots were filtered out.
	 * Short reads retain bases, qualities and error metadata; the other mate is still processed.
	 *
	 * @param cris Started input stream supplying non-null read roots
	 * @param kca Borrowed count table compatible with detection and correction queries
	 * @param k Packed DNA k-mer length in bases
	 * @param thresh Minimum k-mer count threshold
	 * @param maxReads Unused here; the caller must configure any limit on cris
	 * @param ros Started output stream, or null to omit writing
	 */
	public static void detect(ConcurrentReadInputStream cris, KCountArray kca, int k, int thresh, long maxReads, ConcurrentReadOutputStream ros){
		Timer tdetect=new Timer();
		tdetect.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		long covered=0;
		long uncovered=0;
		
		long coveredFinal=0;
		long uncoveredFinal=0;
		
		long fullyCorrected=0;
		long failed=0;
		
		long totalBases=0;
		long totalReads=0;
		
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r : reads){
				Read r2=r.mate;
				{
					
//					if(r.numericID==23){verbose=true;}
					
					totalReads++;
					if(verbose){System.err.println();}
					totalBases+=r.length();
//					BitSet bs=detectErrors(r, kca, k, thresh);
					BitSet bs=detectErrorsBulk(r, kca, k, thresh, 1);
					if(verbose && bs!=null){System.err.println(toString(bs, r.length()));}
//					Data.sysout.println(toString(detectErrorsTips(r, kca, k, thresh), r.length()));
					if(verbose && bs!=null){System.err.println(toString(detectErrors(r, kca, k, thresh), r.length()-k+1));}
					if(bs==null){//No complete window; this does not establish an error-free read.
						//Preserve existing metadata and continue to the other mate.
					}else{
						int x=bs.cardinality();
						covered+=x;
						uncovered+=(r.length()-x);
						if(x<r.length()){
							bs=correctErrors(r, kca, k, thresh, bs, ERROR_CORRECTION_LIMIT, MAX_ERROR_BURST);
						}
						int y=bs.cardinality();
						coveredFinal+=y;
						uncoveredFinal+=(r.length()-y);
						if(x<r.length()){
							if(y==r.length()){
								fullyCorrected++;
							}else{
								failed++;
							}
						}
					}
				}
				if(r2!=null){
					totalReads++;
					totalBases+=r2.length();
//					BitSet bs=detectErrors(r2, kca, k, thresh);
					BitSet bs=detectErrorsBulk(r2, kca, k, thresh, 1);
					if(verbose && bs!=null){System.err.println(toString(bs, r2.length()));}
//					Data.sysout.println(toString(detectErrorsTips(r2, kca, k, thresh), r2.length()));
					if(verbose && bs!=null){System.err.println(toString(detectErrors(r2, kca, k, thresh), r2.length()-k+1));}
					if(bs==null){//No complete window; preserve metadata without counting a correction.
					}else{
						int x=bs.cardinality();
						covered+=x;
						uncovered+=(r2.length()-x);
						if(x<r2.length()){
							bs=correctErrors(r2, kca, k, thresh, bs, ERROR_CORRECTION_LIMIT, MAX_ERROR_BURST);
						}
						int y=bs.cardinality();
						coveredFinal+=y;
						uncoveredFinal+=(r2.length()-y);
						if(x<r2.length()){
							if(y==r2.length()){
								fullyCorrected++;
							}else{
								failed++;
							}
						}
					}
				}
			}
			
			if(ros!=null){ //Important to send all lists to output, even empty ones, to keep list IDs straight.
				if(DONT_OUTPUT_BAD_READS){removeBad(reads);}
				for(Read r : reads){
					if(r!=null){
						r.samline=null;
						assert(r.bases!=null);
						if(r.sites!=null && r.sites.isEmpty()){r.sites=null;}
					}
				}
//				System.err.println("Adding list of length "+readlist.size());
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
		
		tdetect.stop();
		Data.sysout.println("Detect time:    \t"+tdetect);
		Data.sysout.println("Total reads:    \t"+totalReads);
		Data.sysout.println("Total bases:    \t"+totalBases);
		Data.sysout.println("Reads Corrected:\t"+fullyCorrected);
		Data.sysout.println("Reads Failed:   \t"+failed);
		
		Data.sysout.println("\n - before correction - ");
		Data.sysout.println("Covered:        \t"+covered);
		Data.sysout.println("Uncovered:      \t"+uncovered);
		
		Data.sysout.println("\n -  after correction - ");
		Data.sysout.println("Covered:        \t"+coveredFinal);
		Data.sysout.println("Uncovered:      \t"+uncoveredFinal);
	}
	
	/**
	 * Marks the start index of every valid k-mer whose canonical count reaches thresh.
	 * This is a window mask, not a base-coverage mask. Undefined bases break the
	 * rolling window; no input arrays or Read fields are modified.
	 * @param r Non-null read to inspect
	 * @param kca Table supporting canonical packed-key queries
	 * @param k Packed DNA k-mer length in bases, from 1 through 31
	 * @param thresh Minimum supporting count
	 * @return New window-start BitSet, or null when the read is shorter than k
	 */
	public static BitSet detectErrors(final Read r, final KCountArray kca, final int k, final int thresh){
		
		final int kbits=2*k;
		final int shift2=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;

		int bslen=r.length()-k+1;
		if(bslen<1){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(bslen);

		int len=0;
		long kmer=0, rkmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				rkmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|((long)AminoAcid.baseToComplementNumber[b]<<shift2))&mask;
				len++;

				if(len>=k){
					int count=kca.read(Tools.max(kmer, rkmer));//canonical read - see detectErrorsBulk [bloom/ErrorCorrect#001]
					if(count>=thresh){
						bs.set(i+1-k);
					}
				}
			}
		}

		return bs;
	}
	
	/**
	 * Marks every base covered by a sampled valid window whose canonical count reaches thresh.
	 * Sampling starts at the first complete window in each defined-base run;
	 * a valid window ending at the final base is always sampled. A step of one
	 * examines every valid window. Bases and qualities are not modified. Also
	 * writes r.errors as the nonnegative number of bases without sampled support.
	 * @param r Non-null read to inspect
	 * @param kca Table supporting canonical packed-key queries
	 * @param k Packed DNA k-mer length in bases, from 1 through 31
	 * @param thresh Minimum supporting count
	 * @param stepsize Positive step between sampled window starts within each valid run
	 * @return New base-coverage BitSet, or null for absent bases or a read shorter than k
	 */
	public static BitSet detectErrorsBulk(final Read r, final KCountArray kca, final int k, final int thresh, final int stepsize){

		final int kbits=2*k;
		final int shift2=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;

		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.length());
		final int setlen=k;

		int len=0;
		long kmer=0, rkmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				rkmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|((long)AminoAcid.baseToComplementNumber[b]<<shift2))&mask;
				len++;

				//FIXED 2026-07-05 [bloom/ErrorCorrect#001]: read the CANONICAL kmer max(kmer,rkmer),
				//not forward-only. The count arrays are canonical (ReadCounter rcomp=true stores at
				//max(kmer,rkmer); KmerCount5/6 store both strands), so a forward read missed ~half the
				//kmers on canonical arrays (bbnorm/ecc multipass), N-ing correct bases. max() is right
				//for both array types (both-strands populate both cells with the same count).
				if(len>=k && ((len-k)%stepsize==0 || i==bases.length-1)){
					int count=kca.read(Tools.max(kmer, rkmer));
					if(count>=thresh){
						bs.set(i+1-setlen, i+1);
					}
				}
			}
		}
		//BLM205: errors counts unsupported bases, consistently with correctErrors.
		r.errors=r.length()-bs.cardinality();

		return bs;
	}
	
	/**
	 * Starts with all bases retained, then clears bases covered by sampled low-count valid windows.
	 * Samples by absolute window-end index modulo detectStepsize, also sampling a
	 * valid final-base window. Undefined windows are not queried, so a retained bit
	 * may simply be unexamined; it does not prove coverage by a high-count window.
	 * Does not modify the read or table.
	 * @param r Non-null read to inspect
	 * @param kca Table supporting canonical packed-key queries
	 * @param k Packed DNA k-mer length in bases, from 1 through 31
	 * @param thresh Count below which a sampled window clears its covered bases
	 * @param detectStepsize Positive sampling step on absolute window-end indices
	 * @return New retained-base BitSet, or null for absent bases or a read shorter than k
	 */
	public static BitSet detectTrusted(final Read r, final KCountArray kca, final int k, final int thresh, final int detectStepsize){

		final int kbits=2*k;
		final int shift2=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;

		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.length());
		bs.set(0, r.length());
		final int setlen=k;

		int len=0;
		long kmer=0, rkmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				rkmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|((long)AminoAcid.baseToComplementNumber[b]<<shift2))&mask;
				len++;

				//FIXED 2026-07-05 [bloom/ErrorCorrect#001]: canonical read max(kmer,rkmer), not forward-only.
				//Arrays are canonical (ReadCounter) or both-strands (KmerCount5/6); forward-only misread ~half
				//the kmers on canonical arrays, wrongly clearing trusted bases. See detectErrorsBulk.
				if(len>=k && (i%detectStepsize==0 || i==bases.length-1)){
					int count=kca.read(Tools.max(kmer, rkmer));
					if(count<thresh){
						bs.clear(i+1-setlen, i+1);
//						bs.clear(i+1-setlen+detectStepsize, i+1-detectStepsize);
//						bs.clear(i+k/2-detectStepsize, i+k/2+detectStepsize);
//						bs.clear(i+k/2);
					}
				}
			}
		}
//		assert(bases.length==r.length());
		return bs;
	}
	
	/**
	 * Marks the first and last base of each valid window with canonical count at least thresh.
	 * Interior bases are not marked merely because they lie inside that window.
	 * Does not modify the read or table.
	 *
	 * @param r Non-null read to inspect
	 * @param kca Table supporting canonical packed-key queries
	 * @param k Packed DNA k-mer length in bases, from 1 through 31
	 * @param thresh Minimum k-mer count threshold
	 * @return New base-index BitSet, or null for absent bases or a read shorter than k
	 */
	public static BitSet detectErrorsTips(final Read r, final KCountArray kca, final int k, final int thresh){
		
		final int kbits=2*k;
		final int shift2=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final byte[] bases=r.bases;

		if(r.bases==null || r.length()<k){return null;} //Read is too short to detect errors
		BitSet bs=new BitSet(r.length());
		final int setlen=k;

		int len=0;
		long kmer=0, rkmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				rkmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|((long)AminoAcid.baseToComplementNumber[b]<<shift2))&mask;
				len++;

				if(len>=k){
					int count=kca.read(Tools.max(kmer, rkmer));//canonical read - see detectErrorsBulk [bloom/ErrorCorrect#001]
					if(count>=thresh){
						bs.set(i+1-setlen);
						bs.set(i);
					}
				}
			}
		}
		return bs;
	}
	

	/**
	 * Attempts substitutions next to supported regions, using a full bulk-coverage mask.
	 * Mutates bases, writes existing qualities on edits/failures and updates r.errors. After
	 * each successful substitution, recomputes coverage with sampling step one;
	 * callers must use the returned mask rather than assume the original was updated.
	 * A failed attempt can mark its base N and return early. Absent qualities stay
	 * absent. Nearby correction attempts honor the passed maxBurst limit. Correction with
	 * OUTPUT_INFO is explicitly unsupported by an assertion.
	 * @param r Read to modify; bases and coverage must describe the same sequence
	 * @param kca Borrowed table compatible with canonical detection and forward-only candidate queries
	 * @param k Packed DNA k-mer length in bases, from 1 through 31
	 * @param thresh Minimum supporting count
	 * @param bs Non-null bulk-coverage mask computed with every valid window examined
	 * @param maxCorrections Maximum successful substitutions in this invocation
	 * @param maxBurst Maximum allowed nearby correction attempts in a burst
	 * @return Current coverage mask, possibly a new object after successful correction
	 */
	public static BitSet correctErrors(final Read r, final KCountArray kca, final int k, final int thresh, BitSet bs, final int maxCorrections, final int maxBurst){
		//Production callers pass MAX_ERROR_BURST; direct callers may choose another limit.
		
		assert(!OUTPUT_INFO) : "TODO: Outputting correction data is not yet supported.";
		
		int corrections=0; //Alternately, corrections=r.errorsCorrected
		r.errors=0;
		
		if(bs.cardinality()==0){//Cannot be corrected
			r.errors=r.length();
			return bs;
		}

//		verbose=!bs.get(0);
		if(verbose){
			Data.sysout.println();
			Data.sysout.println(toString(bs, r.length()));
			Data.sysout.println(toString(detectErrorsTips(r, kca, k, thresh), r.length()));
			Data.sysout.println(toString(detectErrors(r, kca, k, thresh), r.length()-k+1));
		}
		
		
		int lastloc=-99;
		int burst=1;
		while(!bs.get(0) && corrections<maxCorrections){//While the read starts with a '0', correct from the right.
//			Data.sysout.println("Could not correct.");
//			return bs;
			int errorLoc=bs.nextSetBit(0)-1;//Location to left of first '1'
			if(Tools.absdif(errorLoc,lastloc)<=BURST_THRESH){burst++;}
			else{burst=1;}
			lastloc=errorLoc;
			boolean success=(burst<=maxBurst) && correctFromRight(r, kca, k, thresh, bs, errorLoc);
			if(success){
				corrections++;
				bs=detectErrorsBulk(r, kca, k, thresh, 1);
				if(verbose){System.err.println(">\n"+toString(bs, r.length()));}
			}else{
				r.errors=r.length()-bs.cardinality();
//				r.errorsCorrected+=corrections;
				if(verbose){System.err.println("Could not correct.");}
				r.bases[errorLoc]='N';
				if(r.quality!=null){r.quality[errorLoc]=0;}
				return bs;
			}
		}
		
		burst=1;
		while(bs.cardinality()<r.length() && corrections<maxCorrections){
			if(bs.get(0)){//First bit is a "1", can correct from the left
				int errorLoc=bs.nextClearBit(0);//Location to left of first '0'
				if(Tools.absdif(errorLoc,lastloc)<=BURST_THRESH){burst++;}
				else{burst=1;}
				lastloc=errorLoc;
				boolean success=(burst<=maxBurst) && correctFromLeft(r, kca, k, thresh, bs, errorLoc);
				if(success){
					corrections++;
					bs=detectErrorsBulk(r, kca, k, thresh, 1);
					if(verbose){System.err.println(">\n"+toString(bs, r.length()));}
				}else{
					r.errors=r.length()-bs.cardinality();
//					r.errorsCorrected+=corrections;
					r.bases[errorLoc]='N';
					if(r.quality!=null){r.quality[errorLoc]=0;}
					if(verbose){System.err.println("Could not correct.");}
					return bs;
				}
			}
		}

		r.errors=r.length()-bs.cardinality();
//		r.errorsCorrected+=corrections;
		assert(corrections<=maxCorrections);
		return bs;
	}
	
	
	/**
	 * Tries each replacement base after a supported prefix and compares rightward support.
	 * Queries forward-only packed keys. Restores the original base if no candidate wins;
	 * a winning candidate changes the base and any existing quality, without updating bs.
	 * @param r Read whose bases and optional quality array may be modified
	 * @param kca Counts used to evaluate candidate k-mers
	 * @param k K-mer length in bases
	 * @param thresh Minimum supporting k-mer count
	 * @param bs Base-coverage mask from bulk detection
	 * @param error Index of the first unsupported base after the supported prefix
	 * @return true if a different base wins the minimum-advance comparison
	 */
	private static boolean correctFromLeft(Read r, KCountArray kca, int k, int thresh, BitSet bs, int error){
		final int kbits=2*k;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final int setlen=k;
		final int startLoc=error-(setlen)+1;
		final byte oldBase=r.bases[error];
		final byte[] bases=r.bases;
		
		final int minAdvance=Tools.min(MIN_ADVANCE, bases.length-error);
		
		long kmer=0;
		int len=0;
		for(int i=startLoc; i<error; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from left!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases)+"\nreadID: "+r.numericID);
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+startLoc;
		
		int[] counts=KillSwitch.allocInt1D(4);
		int[] dists=KillSwitch.allocInt1D(4);
		int maxLoc=Tools.min(bases.length-1, error+setlen-1);
		if(!bs.get(error+1)){maxLoc=Tools.min(maxLoc, error+9);}
		else{
			for(int i=error+2; i<=maxLoc; i++){
				if(!bs.get(i)){
					maxLoc=i-1;
					break;
				}
			}
		}
		
		if(verbose){System.err.println("correctFromLeft.  Error = "+error+", maxloc="+maxLoc);}
		for(int bnum=0; bnum<4; bnum++){
			byte c=AminoAcid.numberToBase[bnum];
			bases[error]=c;
			if(verbose){System.err.println("Considering "+(char)c);}
			long key=kmer;
			for(int loc=error; loc<=maxLoc; loc++){
				c=bases[loc];
				int x=AminoAcid.baseToNumber[c];
				if(x<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key<<2)|x)&mask;
				int count=kca.read(key);
				if(count<thresh){
					if(verbose){System.err.println("break: count="+count);}
					break;
				}
				dists[bnum]++;
				counts[bnum]+=count;
			}
		}
		bases[error]=oldBase;
		
		//Note:  I could require both to be the same, to decrease false-positives
		
		final int muid=maxUniqueIndex(dists);
		Arrays.sort(dists);
		final int advance=dists[3];
		final int delta=dists[3]-dists[2];
//		if(advance<minAdvance){return false;}
		if(delta<minAdvance){return false;}
		
		int best=(muid<0 ? maxUniqueIndex(counts) : muid);
		
		if(verbose){System.err.println("Best="+best+": "+Arrays.toString(dists)+"  \t"+Arrays.toString(counts));}
		if(best<0){return false;}
		byte bestC=AminoAcid.numberToBase[best];
		if(bestC==oldBase){return false;}
		bases[error]=bestC;
		
		if(r.quality!=null){r.quality[error]=(byte)Tools.min(10, 3+delta);}
		
		return true;
	}
	
	
	
	/**
	 * Tries each replacement base before a supported suffix and compares leftward support.
	 * Queries forward-only packed keys. Restores the original base if no candidate wins;
	 * a winning candidate changes the base and any existing quality, without updating bs.
	 * @param r Read whose bases and optional quality array may be modified
	 * @param kca Counts used to evaluate candidate k-mers
	 * @param k K-mer length in bases
	 * @param thresh Minimum supporting k-mer count
	 * @param bs Base-coverage mask from bulk detection
	 * @param error Index of the last unsupported base before the supported suffix
	 * @return true if a different base wins the minimum-advance comparison
	 */
	private static boolean correctFromRight(Read r, KCountArray kca, int k, int thresh, BitSet bs, int error){
		final int kbits=2*k;
		final int shift=kbits-2;
		final long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		final int setlen=k;
		final int stopLoc=error+(setlen)-1;
		final byte oldBase=r.bases[error];
		final byte[] bases=r.bases;
		
		final int minAdvance=Tools.min(MIN_ADVANCE, error+1);
		
		long kmer=0;
		int len=0;
		
		for(int i=error+1; i<=stopLoc; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
				throw new RuntimeException("Can't correct from right!\nerror="+error+"\n"+toString(bs, bases.length)+"\n"+new String(bases));
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
			}
//			Data.sysout.print((char)b);
		}
		kmer<<=2;
		
		if(verbose){
			Data.sysout.println();
			String s=Long.toBinaryString(kmer);
			while(s.length()<kbits){s="0"+s;}
			Data.sysout.println("kmer = \t"+s);
		}
		assert(len==setlen-1) : setlen+", "+len+", "+error+", "+stopLoc;
		
		int[] counts=KillSwitch.allocInt1D(4);
		int[] dists=KillSwitch.allocInt1D(4);
		int minLoc=Tools.max(0, error-setlen+1);
		if(error==0 || !bs.get(error-1)){minLoc=Tools.max(minLoc, error-9);}
		else{
			for(int i=error-2; i>=minLoc; i--){
				if(!bs.get(i)){
					minLoc=i+1;
					break;
				}
			}
		}
		
		if(verbose){
			Data.sysout.println("correctFromRight.  Error = "+error+", minloc="+minLoc);
			Data.sysout.println(new String(r.bases));
		}
		for(int bnum=0; bnum<4; bnum++){
			byte c=AminoAcid.numberToBase[bnum];
			bases[error]=c;
			if(verbose){System.err.println("Considering "+(char)c);}
			long key=kmer;
			for(int loc=error; loc>=minLoc; loc--){
				c=bases[loc];
				int x=AminoAcid.baseToNumber[c];
				if(x<0){
					if(verbose){System.err.println("break: N");}
					break;
				}
				key=((key>>2)|(((long)x)<<shift))&mask;
//				{
//					String s=Long.toBinaryString(key);
//					while(s.length()<kbits){s="0"+s;}
//					Data.sysout.println("mask="+Long.toBinaryString(mask)+", shift="+shift+", c="+c+", x="+x+", key  = \t"+s);
//				}
				int count=kca.read(key);
				if(count<thresh){
					if(verbose){System.err.println("break: count="+count);}
					break;
				}
				dists[bnum]++;
				counts[bnum]+=count;
			}
		}
		bases[error]=oldBase;
		
		//Note:  I could require both to be the same, to decrease false-positives
		
		final int muid=maxUniqueIndex(dists);
		Arrays.sort(dists);
		final int advance=dists[3];
		final int delta=dists[3]-dists[2];
//		if(advance<minAdvance){return false;}
		if(delta<minAdvance){return false;}
		
		int best=(muid<0 ? maxUniqueIndex(counts) : muid);
		
		if(verbose){System.err.println("Best="+best+": "+Arrays.toString(dists)+"  \t"+Arrays.toString(counts));}
		if(best<0){return false;}
		byte bestC=AminoAcid.numberToBase[best];
		if(bestC==oldBase){return false;}
		bases[error]=bestC;
		
		if(r.quality!=null){r.quality[error]=(byte)Tools.min(10, 3+delta);}
		
		return true;
	}
	

	/** Returns the unique maximum's index, or -1 for a tie; array must be nonempty and is not modified. */
	private static int maxUniqueIndex(int[] array){
		int max=array[0];
		int maxIndex=0;
		for(int i=1; i<array.length; i++){
			if(array[i]>max){
				max=array[i];
				maxIndex=i;
			}else if(max==array[i]){
				maxIndex=-1;
			}
		}
		return maxIndex;
	}

	/**
	 * Formats exactly len bit positions as '1'/'0' characters for diagnostics.
	 * Does not alter the mask; unset positions beyond its current logical end print zero.
	 *
	 * @param bs Non-null BitSet to inspect
	 * @param len Nonnegative number of positions to print
	 * @return String representation of BitSet
	 */
	public static final String toString(BitSet bs, int len){
//		assert(verbose);
		StringBuilder sb=new StringBuilder(len);
		for(int i=0; i<len; i++){sb.append(bs.get(i) ? '1' : '0');}
		return sb.toString();
	}
	
	/**
	 * Nulls a root only when it has errors and its mate is absent or also has errors.
	 * A pair with either error-free mate is retained together. Does not resize the
	 * list or mutate either read, preserving output list positions and IDs.
	 * @param list Mutable list of non-null read roots
	 */
	private static void removeBad(ArrayList<Read> list){
		
		for(int i=0; i<list.size(); i++){
			Read r=list.get(i);
			if(r.errors>0){
				if(r.mate==null || r.mate.errors>0){
					list.set(i, null);
				}
			}
		}
		
	}
	
	/** Controls verbose debug output during error correction */
	public static boolean verbose=false;
	/** Bails out if a read still has errors after correcting this many. */
	public static int ERROR_CORRECTION_LIMIT=6;
	/** Maximum allowed nearby correction attempts; limits correction bursts that may reflect low coverage. */
	public static int MAX_ERROR_BURST=3;
	/** Bursts have at most this distance between errors. E.G. '1' means errors are adjacent. */
	public static int BURST_THRESH=2;
	/** Filter errored single reads and pairs with errors in both mates; retain pairs with either good mate. */
	public static boolean DONT_OUTPUT_BAD_READS=false;
	/** Required advantage in supported-window advance over the second-best candidate, clipped near read ends. */
	public static int MIN_ADVANCE=1;

	/** Used in output-buffer sizing, not to start correction workers; correction itself runs on the caller thread. */
	public static final int THREADS=1;
	
	/** Select .info output; actual correction in this mode remains unsupported by correctErrors. */
	public static boolean OUTPUT_INFO=false;
	
	
}
