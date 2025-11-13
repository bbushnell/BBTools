package stream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;

/**
 * Factory for creating appropriate Streamer implementations based on file format.
 * Handles single files, paired files, and interleaved formats.
 * Supports FASTQ, SAM/BAM, and FASTA.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class StreamerFactory {

	/**
	 * Creates a Streamer for one or two input files with default settings.
	 * 
	 * @param ff1 Primary input file (R1 for paired data)
	 * @param ff2 Secondary input file (R2 for paired data), or null
	 * @param ordered True to maintain input order in output
	 * @param maxReads Maximum reads to process, or -1 for unlimited
	 * @return Appropriate Streamer implementation
	 */
	public static Streamer makeStreamer(FileFormat ff1, FileFormat ff2, 
			boolean ordered, long maxReads){
		return makeStreamer(ff1, ff2, ordered, maxReads, false, false);
	}
	
	/**
	 * Creates a Streamer for one or two input files.
	 * If ff2 is null, returns a single-file streamer.
	 * If ff2 is non-null, returns a PairStreamer wrapping both files.
	 * Forces ordering when pairing to ensure mate synchronization.
	 * 
	 * @param ff1 Primary input file (R1 for paired data)
	 * @param ff2 Secondary input file (R2 for paired data), or null
	 * @param ordered True to maintain input order in output (forced true if ff2!=null)
	 * @param maxReads Maximum reads to process, or -1 for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @return Appropriate Streamer implementation
	 */
	public static Streamer makeStreamer(FileFormat ff1, FileFormat ff2, 
		boolean ordered, long maxReads, boolean saveHeader, boolean makeReads){
		Streamer s1=makeStreamer(ff1, 0, ordered || ff2!=null, maxReads, saveHeader, makeReads);
		Streamer s2=makeStreamer(ff2, 1, true, maxReads, saveHeader, makeReads);
		return s2==null ? s1 : new PairStreamer(s1, s2);
	}
	
	public static Streamer makeStreamer(FileFormat ff1, FileFormat ff2, 
		boolean ordered, long maxReads, boolean saveHeader, boolean makeReads, int threads){
		Streamer s1=makeStreamer(ff1, 0, ordered || ff2!=null, maxReads, saveHeader, makeReads, threads);
		Streamer s2=makeStreamer(ff2, 1, true, maxReads, saveHeader, makeReads, threads);
		return s2==null ? s1 : new PairStreamer(s1, s2);
	}

	/**
	 * Creates a Streamer for a single input file with default settings.
	 * 
	 * @param ff Input file format
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Maximum reads to process, or -1 for unlimited
	 * @return Appropriate Streamer implementation, or null if ff is null
	 */
	public static Streamer makeStreamer(FileFormat ff, int pairnum, 
			boolean ordered, long maxReads){
		return makeStreamer(ff, pairnum, ordered, maxReads, false, false);
	}
	
	/**
	 * Creates a Streamer for a single input file with full configuration.
	 * 
	 * @param ff Input file format
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Maximum reads to process, or -1 for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @return Appropriate Streamer implementation, or null if ff is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Streamer makeStreamer(FileFormat ff, int pairnum, 
			boolean ordered, long maxReads, boolean saveHeader, boolean makeReads){
		return makeStreamer(ff, pairnum, ordered, maxReads, saveHeader, makeReads, -1);
	}
	
	/**
	 * Creates a Streamer for a single input file with full configuration.
	 * 
	 * @param ff Input file format
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Maximum reads to process, or -1 for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @param threads Worker threads; -1 for auto, 0 for singlethreaded (zero workers)
	 * @return Appropriate Streamer implementation, or null if ff is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Streamer makeStreamer(FileFormat ff, int pairnum, 
			boolean ordered, long maxReads, boolean saveHeader, boolean makeReads, int threads){
		if(ff==null){
			return null;
			
		}else if(ff.fastq()){
			if(Shared.threads()>=4 && threads!=0 && (threads>1 || FastqStreamer.DEFAULT_THREADS>1)) {
				return new FastqStreamer(ff, threads, pairnum, maxReads);
			}else {
				return new FastqStreamerST(ff, pairnum, maxReads);
			}
			
		}else if(ff.fasta()){
			if(threads==0 || (threads<0 && FastaStreamer.DEFAULT_THREADS<1) || 
					Shared.threads()<4 || Shared.LOW_MEMORY) {
				if(FASTA_STREAMER_2 && Shared.SIMD && !ff.interleaved()) {
					return new FastaStreamer2ZT(ff, pairnum, maxReads);
				}else {
					return new FastaStreamerZT(ff, pairnum, maxReads);
				}
			}else if(threads==1 || (threads<0 && FastaStreamer.DEFAULT_THREADS<2) || Shared.threads()<8) {
				if(FASTA_STREAMER_2 && Shared.SIMD && !ff.interleaved()) {
					return new FastaStreamer2ST(ff, pairnum, maxReads);
				}else {
					return new FastaStreamerST(ff, pairnum, maxReads);
				}
			}else {
				return new FastaStreamer(ff, threads, pairnum, maxReads);
			}
			
		}else if(ff.bam() && ReadWrite.nativeBamIn()){
			return new BamStreamer(ff, threads, saveHeader, ordered, maxReads, makeReads);
			
		}else if(ff.samOrBam()){
			if(Shared.threads()>=4 && threads!=0 && (threads>1 || SamStreamer.DEFAULT_THREADS>0)) {
				return new SamStreamer(ff, threads, saveHeader, ordered, maxReads, makeReads);
			}else {
				return new SamStreamerST(ff, saveHeader, maxReads, makeReads);
			}
			
		}
		
		throw new RuntimeException("Unsupported file format: "+ff);
	}
	
	//Legacy
	
	public static Streamer makeSamOrBamStreamer(String fname, int threads, boolean saveHeader, boolean ordered, long maxReads, boolean makeReads) {
		return makeSamOrBamStreamer(FileFormat.testInput(fname, FileFormat.SAM, null, true, false), threads, saveHeader, ordered, maxReads, makeReads);
	}
	
	public static Streamer makeSamOrBamStreamer(FileFormat ffin, int threads, boolean saveHeader, boolean ordered, long maxReads, boolean makeReads) {
		return makeStreamer(ffin, 0, ordered, maxReads, saveHeader, makeReads, threads);
	}
	
	//Generally faster and makes fewer objects
	//It does scan for newlines twice though
	public static boolean FASTA_STREAMER_2=true;
	
}