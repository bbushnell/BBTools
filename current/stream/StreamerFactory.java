package stream;

import fileIO.FileFormat;

/**
 * Factory for creating appropriate Streamer implementations based on file format.
 * Handles single files, paired files, and interleaved formats.
 * Supports FASTQ, SAM/BAM, and FASTA (when implemented).
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
		if(ff==null){
			return null;
		}else if(ff.fastq()){
			return new FastqStreamer(ff, FastqStreamer.DEFAULT_THREADS, pairnum, maxReads);
		}else if(ff.samOrBam()){
			return SamStreamer.makeStreamer(ff, SamStreamer.DEFAULT_THREADS, saveHeader, 
				ordered, maxReads, makeReads);
		}else if(ff.fasta()){//TODO
			//return new FastaStreamer(ff, FastaStreamer.DEFAULT_THREADS, pairnum, maxReads);
			throw new RuntimeException("FASTA streaming not yet implemented");
		}
		
		throw new RuntimeException("Unsupported file format: "+ff);
	}
	
}