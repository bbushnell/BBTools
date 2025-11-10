package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;

/**
 * Factory for creating appropriate Writer implementations based on file format.
 * Handles single files (interleaved or unpaired) and paired files.
 * Supports FASTQ and SAM/BAM output.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class WriterFactory {
	
	/*--------------------------------------------------------------*/
	/*----------------           Legacy             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** For drop-in ConcurrentReadOutputStream support */
	public static Writer getStream(FileFormat ffout1, FileFormat ffout2, 
			int buffersUnused, ArrayList<byte[]> header, boolean useSharedHeader) {
		return makeWriter(ffout1, ffout2, header, useSharedHeader);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Twin Files           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a Writer for one or two output files with default settings.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @return Appropriate Writer implementation
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2){
		return makeWriter(ffout1, ffout2, null, false);
	}
	
	/**
	 * Creates a Writer for one or two output files with specified thread count.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param threads Number of compression/formatting threads per file
	 * @return Appropriate Writer implementation
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2, int threads){
		return makeWriter(ffout1, ffout2, threads, null, false);
	}
	
	/**
	 * Creates a Writer for one or two output files.
	 * If ffout2 is null, returns a single-file writer (interleaved or unpaired).
	 * If ffout2 is non-null, returns a PairedWriter for separate R1/R2 files.
	 * Both files must be ordered when paired to ensure mate synchronization.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param header SAM/BAM header lines, or null
	 * @param useSharedHeader True to share header reference across threads (SAM/BAM only)
	 * @return Appropriate Writer implementation
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2,
			ArrayList<byte[]> header, boolean useSharedHeader){
		if(ffout2==null){
			// Single file - interleaved or unpaired
			return makeWriter(ffout1, true, true, header, useSharedHeader);
		}else{
			// Paired files
			assert(ffout1.ordered());
			assert(ffout2.ordered());
			Writer w1=makeWriter(ffout1, true, false, header, false);  // R1 only
			Writer w2=makeWriter(ffout2, false, true, header, false);  // R2 only
			return new PairedWriter(w1, w2);
		}
	}
	
	/**
	 * Creates a Writer for one or two output files with full configuration.
	 * If ffout2 is null, returns a single-file writer (interleaved or unpaired).
	 * If ffout2 is non-null, returns a PairedWriter for separate R1/R2 files.
	 * Both files must be ordered when paired to ensure mate synchronization.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param threads Number of compression/formatting threads per file
	 * @param header SAM/BAM header lines, or null
	 * @param useSharedHeader True to share header reference across threads (SAM/BAM only)
	 * @return Appropriate Writer implementation
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2, int threads,
			ArrayList<byte[]> header, boolean useSharedHeader){
		if(ffout2==null){
			// Single file - interleaved or unpaired
			return makeWriter(ffout1, true, true, threads, header, useSharedHeader);
		}else{
			// Paired files
			assert(ffout1.ordered());
			assert(ffout2.ordered());
			Writer w1=makeWriter(ffout1, true, false, threads, header, false);  // R1 only
			Writer w2=makeWriter(ffout2, false, true, threads, header, false);  // R2 only
			return new PairedWriter(w1, w2);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Single File          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a Writer for a single output file with default settings.
	 * Writes both R1 and R2 (interleaved) by default.
	 * 
	 * @param ffout Output file format
	 * @return Appropriate Writer implementation, or null if ffout is null
	 */
	public static Writer makeWriter(FileFormat ffout){
		return makeWriter(ffout, true, true, null, false);
	}

	/**
	 * Creates a Writer for a single output file with specified read selection.
	 * 
	 * @param ffout Output file format
	 * @param writeR1 True to write R1 reads (pairnum 0)
	 * @param writeR2 True to write R2 reads (pairnum 1)
	 * @return Appropriate Writer implementation, or null if ffout is null
	 */
	public static Writer makeWriter(FileFormat ffout, boolean writeR1, boolean writeR2){
		return makeWriter(ffout, writeR1, writeR2, null, false);
	}

	/**
	 * Creates a Writer for a single output file.
	 * 
	 * @param ffout Output file format
	 * @param writeR1 True to write R1 reads (pairnum 0)
	 * @param writeR2 True to write R2 reads (pairnum 1)
	 * @param header SAM/BAM header lines, or null
	 * @param useSharedHeader True to share header reference across threads (SAM/BAM only)
	 * @return Appropriate Writer implementation, or null if ffout is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Writer makeWriter(FileFormat ffout, boolean writeR1, boolean writeR2, 
			ArrayList<byte[]> header, boolean useSharedHeader){
		if(ffout==null){
			return null;
		}else if(ffout.fastq() || ffout.fasta()){
			if(Shared.threads()>=4 && FastqWriter.DEFAULT_THREADS>0) {
				return new FastqWriter(ffout, -1, writeR1, writeR2);
			}else {
				return new FastqWriterST(ffout, writeR1, writeR2);
			}
		}else if(ffout.header()){
			return new FastqWriterST(ffout, writeR1, writeR2);
		}else if(ffout.bam() && ReadWrite.nativeBamOut()){
			return new BamWriter(ffout, -1, header, useSharedHeader);
		}else if(ffout.samOrBam()){
			if(Shared.threads()>=4 && SamWriter.DEFAULT_THREADS>0) {
				return new SamWriter(ffout, -1, header, useSharedHeader);
			}else {
				return new SamWriterST(ffout, header, useSharedHeader);
			}
		}

		throw new RuntimeException("Unsupported file format: "+ffout);
	}

	/**
	 * Creates a Writer for a single output file with full configuration.
	 * 
	 * @param ffout Output file format
	 * @param writeR1 True to write R1 reads (pairnum 0)
	 * @param writeR2 True to write R2 reads (pairnum 1)
	 * @param threads Number of compression/formatting threads
	 * @param header SAM/BAM header lines, or null
	 * @param useSharedHeader True to share header reference across threads (SAM/BAM only)
	 * @return Appropriate Writer implementation, or null if ffout is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Writer makeWriter(FileFormat ffout, boolean writeR1, boolean writeR2, int threads,
			ArrayList<byte[]> header, boolean useSharedHeader){
		if(ffout==null){
			return null;
		}else if(ffout.fastq() || ffout.fasta()){
			if(Shared.threads()>=4 && threads!=0 && (threads>1 || FastqWriter.DEFAULT_THREADS>0)) {
				return new FastqWriter(ffout, -1, writeR1, writeR2);
			}else {
				return new FastqWriterST(ffout, writeR1, writeR2);
			}
		}else if(ffout.header()){
			return new FastqWriterST(ffout, writeR1, writeR2);
		}else if(ffout.bam() && ReadWrite.nativeBamOut()){
			return new BamWriter(ffout, threads, header, useSharedHeader);
		}else if(ffout.samOrBam()){
			if(Shared.threads()>=4 && threads!=0 && (threads>1 || SamWriter.DEFAULT_THREADS>0)) {
				return new SamWriter(ffout, threads, header, useSharedHeader);
			}else {
				return new SamWriterST(ffout, header, useSharedHeader);
			}
		}

		throw new RuntimeException("Unsupported file format: "+ffout);
	}

}