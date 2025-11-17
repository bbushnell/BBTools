package aligner;

import java.util.ArrayList;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBDuk;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamByteWriter;
import stream.SamLine;

/**
 * Side channel mapper for performing lightweight alignment of reads to reference sequences.
 * Uses dual k-mer based indices for fast mapping with configurable identity thresholds.
 * Provides concurrent output streams for mapped and unmapped reads with statistics tracking.
 * @author Brian Bushnell
 */
public class SideChannel3 {
	
	/*--------------------------------------------------------------*/
	/*----------------          Constructor         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a SideChannel3 mapper with single k-mer configuration.
	 * Creates indices and aligners using the specified k-mer size and identity threshold.
	 *
	 * @param ref_ Reference sequence file path
	 * @param out_ Output file for mapped reads
	 * @param outu_ Output file for unmapped reads
	 * @param k1_ K-mer size for indexing
	 * @param minid1 Minimum identity threshold for alignment
	 * @param midMaskLen1 Length of middle masking for fuzzy k-mer matching
	 * @param overwrite_ Whether to overwrite existing output files
	 * @param ordered_ Whether to maintain input order in output
	 */
	public SideChannel3(String ref_, String out_, String outu_, int k1_, float minid1, 
			int midMaskLen1, boolean overwrite_, boolean ordered_) {
		this(ref_, out_, outu_, k1_, -1, minid1, 1, midMaskLen1, 0, overwrite_, ordered_);
	}
	
	/**
	 * Constructs a SideChannel3 mapper with dual k-mer configuration.
	 * Creates two indices with different k-mer sizes for improved sensitivity.
	 * The larger k-mer is used for initial mapping, smaller for unmapped reads.
	 *
	 * @param ref_ Reference sequence file path
	 * @param out_ Output file for mapped reads
	 * @param outu_ Output file for unmapped reads
	 * @param k1_ First k-mer size (will be set to max of k1_, k2_)
	 * @param k2_ Second k-mer size (will be set to min of k1_, k2_)
	 * @param minid1 Minimum identity threshold for first mapper
	 * @param minid2 Minimum identity threshold for second mapper
	 * @param midMaskLen1 Middle masking length for first k-mer size
	 * @param midMaskLen2 Middle masking length for second k-mer size
	 * @param overwrite_ Whether to overwrite existing output files
	 * @param ordered_ Whether to maintain input order in output
	 */
	public SideChannel3(String ref_, String out_, String outu_, int k1_, int k2_, float minid1, float minid2, 
			int midMaskLen1, int midMaskLen2, boolean overwrite_, boolean ordered_) {
		Timer t=new Timer();
		ref=fixRefPath(ref_);
		out=out_;
		outu=outu_;
		k1=Tools.max(k1_, k2_);
		k2=Tools.min(k1_, k2_);
		minIdentity1=fixID(minid1);
		minIdentity2=fixID(minid2);
		overwrite=overwrite_;
		ordered=ordered_;
		assert(k1>0);

		ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, overwrite, false, ordered);
		ffoutu=FileFormat.testOutput(outu, FileFormat.FASTQ, null, true, overwrite, false, ordered);
		samOut=((ffout!=null && ffout.samOrBam()) || ((ffoutu!=null && ffoutu.samOrBam())));
		final Read r=MicroIndex3.loadRef(ref, samOut);
		index1=new MicroIndex3(k1, midMaskLen1, r);
		index2=(k2<1 ? null : new MicroIndex3(k2, midMaskLen2, r));
		mapper1=new MicroAligner3(index1, minIdentity1, true);
		mapper2=(k2<1 ? null : new MicroAligner3(index2, minIdentity2, true));

		if(samOut) {ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;}
		final int buff=(!ordered ? 12 : Tools.max(32, 2*Shared.threads()));
		if(ffout!=null) {
			cros=ConcurrentReadOutputStream.getStream(ffout, null, buff, null, false);
		}else {cros=null;}
		if(ffoutu!=null) {
			crosu=ConcurrentReadOutputStream.getStream(ffoutu, null, buff, null, false);
		}else {crosu=null;}
		t.stop("Created side channel"+(out==null ? "" : (" for "+out))+": ");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Maps a read pair using the configured aligners.
	 * Convenience method that delegates to the main mapping method.
	 *
	 * @param r1 First read of the pair
	 * @param r2 Second read of the pair (may be null for single-end)
	 * @return true if at least one read mapped successfully
	 */
	public boolean map(Read r1, Read r2) {
		return map(r1, r2, mapper1, mapper2);
	}
	
	/**
	 * Maps a read pair using specified aligners with dual k-mer strategy.
	 * Uses the first mapper for initial alignment, then tries the second mapper
	 * on unmapped reads. Sets proper pair flags for paired reads on same chromosome.
	 *
	 * @param r1 First read of the pair
	 * @param r2 Second read of the pair (may be null for single-end)
	 * @param mapper1 Primary aligner (larger k-mer)
	 * @param mapper2 Secondary aligner (smaller k-mer, may be null)
	 * @return true if at least one read mapped with sufficient identity
	 */
	public boolean map(Read r1, Read r2, MicroAligner3 mapper1, MicroAligner3 mapper2) {
		float id1=mapper1.map(r1);
		float id2=mapper1.map(r2);
		if(id1+id2<=0) {return false;}//Common case
		
		if(r2!=null) {
			if(mapper2!=null) {
				if(r1.mapped() && !r2.mapped()) {id2=mapper2.map(r2);}
				else if(r2.mapped() && !r1.mapped()) {id1=mapper2.map(r1);}
			}
			boolean properPair=(r1.mapped() && r2.mapped() && r1.chrom==r2.chrom && 
					r1.strand()!=r2.strand() && Tools.absdif(r1.start, r2.start)<=1000);
			r1.setPaired(properPair);
			r2.setPaired(properPair);
		}
		
		if(!r1.mapped()) {id1=0;}
		if(r2==null || !r2.mapped()) {id2=0;}
		long idsum=(long)((id1+id2)*10000);
		if(idsum<=0) {return false;}
		
		return true;
	}
	
	/**
	 * Generates alignment statistics summary string.
	 * Reports mapped read counts, percentages, and average identity.
	 *
	 * @param readsIn Total number of input reads processed
	 * @param basesIn Total number of input bases processed
	 * @return Formatted statistics string with counts and percentages
	 */
	public String stats(long readsIn, long basesIn) {
		long ro, bo, idsum, rm;
		ro=readsOut; bo=basesOut; idsum=identitySum; rm=readsMapped;
		String s=("Aligned reads:          \t"+ro+" reads ("+BBDuk.toPercent(ro, readsIn)+") \t"+
				+bo+" bases ("+BBDuk.toPercent(bo, basesIn)+") \tavgID="+Tools.format("%.4f", idsum/(100.0*rm)));
		return s;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Starts the concurrent output streams for writing results.
	 * Must be called before writing mapped or unmapped reads. */
	public void start() {
		if(cros!=null) {cros.start();}
		if(crosu!=null) {crosu.start();}
	}
	
	/**
	 * Shuts down output streams and finalizes output files.
	 * Should be called when all reads have been processed.
	 * @return true if any errors occurred during shutdown
	 */
	public boolean shutdown() {
		errorState=ReadWrite.closeOutputStreams(cros, crosu)|errorState;
		return errorState;
	}

	/** Everything in this list gets written */
	public int writeToMapped(ArrayList<Read> reads, long num) {
		if(reads==null || (reads.isEmpty() && !ordered) || (cros==null && !TRACK_STATS)) {
			assert(!ordered || reads!=null);
			return 0;
		}
		
		int rm=0, ro=0, bo=0;
		double idsum=0;
		for(Read r1 : reads) {
			Read r2=r1.mate;
			ro+=r1.pairCount();
			bo+=r1.pairLength();
			rm+=r1.pairMappedCount();

			if(TRACK_STATS) {
				if(r1.mapped()) {idsum+=r1.identity();}
				if(r2!=null && r2.mapped()) {idsum+=r2.identity();}
			}
		}


		if(TRACK_STATS && ro>0) {
			synchronized(this) {
				readsMapped+=rm;
				readsOut+=ro;
				basesOut+=bo;
				identitySum+=(long)(idsum*10000);
			}
		}
		if(cros==null) {return ro;}
		
		final boolean makeSamLine=(samOut && ReadStreamByteWriter.USE_ATTACHED_SAMLINE);
		if(makeSamLine) {
			for(Read r1 : reads) {
				Read r2=r1.mate;
				r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
				if(r2!=null) {r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));}
			}
		}
		cros.add(reads, num);
		return ro;
	}

	/** Expects a list of mixed mapped and unmapped reads */
	public int writeByStatus(ArrayList<Read> reads, long num) {
		int ro=writeMappedOnly(reads, num);
		writeUnmappedOnly(reads, num);
		return ro;
	}

	/** Expects a list of mixed mapped and unmapped reads; 
	 * only writes mapped reads (plus mates) */
	private int writeMappedOnly(ArrayList<Read> reads, long num) {
		if(reads==null || (reads.isEmpty() && !ordered)) {
			assert(!ordered || reads!=null);
			return 0;
		}
		
		int listSize=0;
		int rm=0, ro=0, bo=0, rou=0, bou=0;
		double idsum=0;
		for(Read r1 : reads) {
			boolean mapped=r1.eitherMapped();
			if(mapped) {
				Read r2=r1.mate;
				listSize++;
				ro+=r1.pairCount();
				bo+=r1.pairLength();
				rm+=r1.pairMappedCount();

				if(TRACK_STATS) {
					if(r1.mapped()) {idsum+=r1.identity();}
					if(r2!=null && r2.mapped()) {idsum+=r2.identity();}
				}
			}else{
				rou+=r1.pairCount();
				bou+=r1.pairLength();
			}
		}


		if(TRACK_STATS && ro>0) {
			synchronized(this) {
				readsMapped+=rm;
				readsOut+=ro;
				basesOut+=bo;
				identitySum+=(long)(idsum*10000);
			}
		}

		if(cros==null || (!ordered && listSize<1)) {return ro;}

		ArrayList<Read> list=new ArrayList<Read>(Tools.max(1, listSize));
		final boolean makeSamline=(samOut && ReadStreamByteWriter.USE_ATTACHED_SAMLINE);
		for(Read r1 : reads) {
			boolean mapped=r1.eitherMapped();
			if(mapped && list!=null) {
				Read r2=r1.mate;
				if(makeSamline) {
					r1.samline=(r1.samline!=null ? r1.samline : new SamLine(r1, 0));
					if(r2!=null) {r2.samline=(r2.samline!=null ? r2.samline : new SamLine(r2, 1));}
				}
				list.add(r1);
			}
		}
		cros.add(list, num);
		return ro;
	}

	/** Expects a list of mixed mapped and unmapped reads; 
	 * only writes unmapped reads with unmapped mates */
	private int writeUnmappedOnly(ArrayList<Read> reads, long num) {
		if(reads==null || crosu==null || (reads.isEmpty() && !ordered)) {
			assert(!ordered || reads!=null);
			return 0;
		}

		ArrayList<Read> list=new ArrayList<Read>(8);
		int rou=0;
		for(Read r1 : reads) {
			if(!r1.eitherMapped()) {
				rou+=r1.pairCount();
				list.add(r1);
			}
		}
		crosu.add(list, num);
		return rou;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses comma-separated k-mer values from command line argument.
	 * Returns array with up to 2 k-mer values for dual k-mer configuration.
	 *
	 * @param arg Original argument name (unused)
	 * @param a Argument key (unused)
	 * @param b Comma-separated k-mer values
	 * @return Array containing parsed k-mer integers
	 */
	public static int[] parseK(String arg, String a, String b) {
		int[] ret=new int[2];
		String[] terms=b.split(",");
		for(int i=0; i<terms.length; i++) {
			ret[i]=Integer.parseInt(terms[i]);
		}
		return ret;
	}
	
	/**
	 * Normalizes identity values to 0-1 range.
	 * Converts percentage values (>1) to decimal format.
	 * @param id Identity value to normalize
	 * @return Normalized identity in 0-1 range
	 */
	static float fixID(float id) {
		if(id>1) {id=id/100;}
		assert(id<=1);
		return id;
	}
	
	/**
	 * Resolves reference file paths, handling special cases.
	 * Converts "phix" shorthand to full PhiX reference path.
	 * @param refPath Original reference path
	 * @return Resolved reference file path
	 */
	static String fixRefPath(String refPath) {
		if(refPath==null || Tools.isReadableFile(refPath)) {return refPath;}
		if("phix".equalsIgnoreCase(refPath)){return Data.findPath("?phix2.fa.gz");}
		return refPath;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary k-mer index for first alignment pass */
	public final MicroIndex3 index1;
	/** Secondary k-mer index for second alignment pass (may be null) */
	public final MicroIndex3 index2;
	/** Primary aligner using index1 */
	public final MicroAligner3 mapper1;
	/** Secondary aligner using index2 (may be null) */
	public final MicroAligner3 mapper2;
	
	/** Primary k-mer size (larger of the two k values) */
	public final int k1;
	/** Secondary k-mer size (smaller of the two k values) */
	public final int k2;
	/** Minimum identity threshold for primary aligner */
	public final float minIdentity1;
	/** Minimum identity threshold for secondary aligner */
	public final float minIdentity2;
	
	/** Reference sequence file path */
	public final String ref;
	/** Output file path for mapped reads */
	public final String out;
	/** Output file path for unmapped reads */
	public final String outu;
	/** Whether output format is SAM/BAM */
	public final boolean samOut;
	/** File format specification for mapped output */
	public final FileFormat ffout;
	/** File format specification for unmapped output */
	public final FileFormat ffoutu;
	/** Concurrent output stream for mapped reads */
	private final ConcurrentReadOutputStream cros;
	/** Concurrent output stream for unmapped reads */
	private final ConcurrentReadOutputStream crosu;
	
	/** Count of reads that successfully mapped */
	public long readsMapped=0;
	/** Total count of reads written to output */
	public long readsOut=0;
	/** Total count of bases written to output */
	public long basesOut=0;
	/** Cumulative identity sum scaled by 10000 (0-10000 scale per read) */
	public long identitySum=0;//x100%; 0-10000 scale. 
	
	/** Whether to overwrite existing output files */
	public final boolean overwrite;
	/** Whether to maintain input order in output files */
	public final boolean ordered;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Global flag controlling whether to track alignment statistics */
	public static boolean TRACK_STATS=true;

}
