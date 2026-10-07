package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;

/**
 * Creates format-specific writers for one output or two separate mate outputs.
 * Supports FASTQ, FASTA, FASTA+QUAL, HEADER, SCARF and SAM/BAM. Each output
 * descriptor is routed independently; paired descriptors need not match formats.
 * Constructors can open files, but this factory does not call start(); the caller
 * owns startup, submission and completion under the selected writer's contract.
 * Explicit byte-array header lists are passed through without copying the list
 * or its arrays; the String-header adapter creates new arrays.
 *
 * Thread requests select a backend, not an exact total thread or compression budget.
 * FASTQ/FASTA and SAM use worker pools only above one requested worker, with at
 * least eight shared threads and LOW_MEMORY disabled. Otherwise their ST2 backend
 * may use one worker or caller-driven output. Native BAM, SCARF, HEADER and
 * FASTA+QUAL have separate policies; see the fully configured single-file method.
 * Use {@link LightweightWriterFactory} when lightweight leaf selection is required.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class WriterFactory{
	
	/*--------------------------------------------------------------*/
	/*----------------           Legacy             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Adapts a newline-delimited header and legacy buffer argument to the writer API.
	 * Header lines use String.split semantics and the platform-default byte encoding.
	 * @param ffout1 Primary output; required when ffout2 is nonnull
	 * @param ffout2 Optional separate mate output
	 * @param qf1 Optional QUAL filename for a FASTA primary output
	 * @param qf2 Optional QUAL filename for a separate FASTA mate output
	 * @param buffersUnused Ignored legacy buffer count
	 * @param header Optional newline-delimited SAM/BAM header
	 * @param useSharedHeader Request shared input-header lookup for single-file SAM/BAM
	 * @param threads Backend thread request; not a total thread limit
	 * @return New writer, or null when both output descriptors are null */
	public static Writer getStreamS(FileFormat ffout1, FileFormat ffout2, String qf1, String qf2,
			int buffersUnused, String header, boolean useSharedHeader, int threads){
		ArrayList<byte[]> headers=null;
		if(header!=null){
			headers=new ArrayList<byte[]>();
			for(String s : header.split("\n")){
				headers.add(s.getBytes());
			}
		}
		return makeWriter(ffout1, ffout2, qf1, qf2, threads, headers, useSharedHeader);
	}
	
	/** Adapts legacy arguments without separate QUAL output; does not start the writer.
	 * @param ffout1 Primary output; required when ffout2 is nonnull
	 * @param ffout2 Optional separate mate output
	 * @param buffersUnused Ignored legacy buffer count
	 * @param header Optional borrowed SAM/BAM header lines
	 * @param useSharedHeader Request shared input-header lookup for single-file SAM/BAM
	 * @param threads Backend thread request; not a total thread limit
	 * @return New writer, or null when both output descriptors are null */
	public static Writer getStream(FileFormat ffout1, FileFormat ffout2,
			int buffersUnused, ArrayList<byte[]> header, boolean useSharedHeader, int threads){
		return makeWriter(ffout1, ffout2, threads, header, useSharedHeader);
	}
	
	/** Adapts legacy arguments including optional QUAL outputs; ignores buffer count.
	 * @param ffout1 Primary output; required when ffout2 is nonnull
	 * @param ffout2 Optional separate mate output
	 * @param qf1 Optional QUAL filename for a FASTA primary output
	 * @param qf2 Optional QUAL filename for a separate FASTA mate output
	 * @param buffersUnused Ignored legacy buffer count
	 * @param header Optional borrowed SAM/BAM header lines
	 * @param useSharedHeader Request shared input-header lookup for single-file SAM/BAM
	 * @param threads Backend thread request; not a total thread limit
	 * @return New writer, or null when both output descriptors are null */
	public static Writer getStream(FileFormat ffout1, FileFormat ffout2, String qf1, String qf2,
			int buffersUnused, ArrayList<byte[]> header, boolean useSharedHeader, int threads){
		return makeWriter(ffout1, ffout2, qf1, qf2, threads, header, useSharedHeader);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Twin Files           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates one or two outputs with backend defaults and no supplied/shared header.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @return New writer, or null when both output descriptors are null
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2){
		return makeWriter(ffout1, ffout2, null, false);
	}
	
	/**
	 * Creates a Writer for one or two output files with specified thread count.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param threads Per-file backend thread request; not a total thread limit
	 * @return New writer, or null when both output descriptors are null
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2, int threads){
		return makeWriter(ffout1, ffout2, threads, null, false);
	}
	
	/**
	 * Creates a Writer for one or two output files.
	 * If ffout2 is null, returns a single-file writer (interleaved or unpaired).
	 * If ffout2 is non-null, returns a PairedWriter for separate R1/R2 files.
	 * Both paired descriptors must request ordering (checked with assertions enabled).
	 * Separate mate files disable shared-header lookup for both children and pass
	 * the same supplied header to each. A primary descriptor is required when paired.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param header Borrowed SAM/BAM header lines, or null for backend selection
	 * @param useSharedHeader Request shared input-header lookup; ignored for paired files
	 * @return New writer, or null when both output descriptors are null
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
			Writer w1=makeWriter(ffout1, true, false, header, false);// R1 only
			Writer w2=makeWriter(ffout2, false, true, header, false);// R2 only
			return new PairedWriter(w1, w2);
		}
	}
	
	/**
	 * Creates one or two outputs without separate QUAL files.
	 * Paired output requires a primary descriptor and ordering on both descriptors;
	 * it disables shared-header lookup for both children.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param threads Per-file backend thread request; not a total thread limit
	 * @param header Borrowed SAM/BAM header lines, or null for backend selection
	 * @param useSharedHeader Request shared input-header lookup; ignored for paired files
	 * @return New writer, or null when both output descriptors are null
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2, int threads,
			ArrayList<byte[]> header, boolean useSharedHeader){
		return makeWriter(ffout1, ffout2, null, null, threads, header, useSharedHeader);
	}
	
	/**
	 * Creates a Writer for one or two output files with full configuration.
	 * If ffout2 is null, returns a single-file writer (interleaved or unpaired).
	 * If ffout2 is non-null, returns a PairedWriter for separate R1/R2 files.
	 * Both paired descriptors must request ordering (checked with assertions enabled).
	 * Separate mate files require a primary descriptor, disable shared-header lookup
	 * and pass the same supplied header to each child. qf2 is ignored for one output;
	 * QUAL filenames are used only for FASTA outputs. Children are created in order.
	 * 
	 * @param ffout1 Primary output file (R1 for paired data, or interleaved/unpaired)
	 * @param ffout2 Secondary output file (R2 for paired data), or null
	 * @param qf1 Optional QUAL filename for a FASTA primary output
	 * @param qf2 Optional QUAL filename for a separate FASTA mate output
	 * @param threads Per-file backend thread request; not a total thread limit
	 * @param header Borrowed SAM/BAM header lines, or null for backend selection
	 * @param useSharedHeader Request shared input-header lookup; ignored for paired files
	 * @return New writer, or null when both output descriptors are null
	 */
	public static Writer makeWriter(FileFormat ffout1, FileFormat ffout2, String qf1, String qf2,
			int threads, ArrayList<byte[]> header, boolean useSharedHeader){
		if(ffout2==null){
			// Single file - interleaved or unpaired
			return makeWriter(ffout1, qf1, true, true, threads, header, useSharedHeader);
		}else{
			// Paired files
			assert(ffout1.ordered());
			assert(ffout2.ordered());
			Writer w1=makeWriter(ffout1, qf1, true, false, threads, header, false);// R1 only
			Writer w2=makeWriter(ffout2, qf2, false, true, threads, header, false);// R2 only
			return new PairedWriter(w1, w2);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Single File          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a Writer for a single output file with default settings.
	 * Selects both mates, backend thread defaults and no supplied/shared header.
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
	 * @param writeR1 Select first entries according to the backend's Read conversion
	 * @param writeR2 Select second mates according to the backend's Read conversion
	 * @return Appropriate Writer implementation, or null if ffout is null
	 */
	public static Writer makeWriter(FileFormat ffout, boolean writeR1, boolean writeR2){
		return makeWriter(ffout, writeR1, writeR2, null, false);
	}

	/**
	 * Creates a Writer for a single output file.
	 * 
	 * @param ffout Output file format
	 * @param writeR1 Select first entries according to the backend's Read conversion
	 * @param writeR2 Select second mates according to the backend's Read conversion
	 * @param header Borrowed SAM/BAM header lines, or null for backend selection
	 * @param useSharedHeader Request shared input-header lookup for SAM/BAM
	 * @return Appropriate Writer implementation, or null if ffout is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Writer makeWriter(FileFormat ffout, boolean writeR1, boolean writeR2,
			ArrayList<byte[]> header, boolean useSharedHeader){
		return makeWriter(ffout, writeR1, writeR2, -1, header, useSharedHeader);
	}

	/**
	 * Creates one output without a separate QUAL file; delegates backend selection.
	 * 
	 * @param ffout Output file format
	 * @param writeR1 Select first entries according to the backend's Read conversion
	 * @param writeR2 Select second mates according to the backend's Read conversion
	 * @param threads Backend thread request; negative selects backend defaults
	 * @param header Borrowed SAM/BAM header lines, or null for backend selection
	 * @param useSharedHeader Request shared input-header lookup for SAM/BAM
	 * @return Appropriate Writer implementation, or null if ffout is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Writer makeWriter(FileFormat ffout, boolean writeR1,
			boolean writeR2, int threads, ArrayList<byte[]> header, boolean useSharedHeader){
		return makeWriter(ffout, null, writeR1, writeR2, threads, header, useSharedHeader);
	}

	/**
	 * Creates a writer using the descriptor's format and current shared settings.
	 * FASTQ and FASTA without QUAL select FastqWriter or FastqWriterST2. Their ST2
	 * path uses a worker only for a positive thread request and at least four shared
	 * threads; its queue capacity is three for FASTA and five for FASTQ.
	 * HEADER always selects caller-driven FastqWriterST2 with capacity three.
	 * Native BAM output always selects BamWriter and leaves thread normalization to it.
	 * Other SAM/BAM selects SamWriter or SamWriterST2; that ST2 path uses a worker
	 * for any positive thread request, with capacity five. SCARF always selects
	 * FastqWriter with at least one requested worker. FASTA+QUAL always selects
	 * FastaQualWriterST, ignoring the thread request. Negative FASTQ/FASTA/SCARF
	 * requests use FastqWriter.DEFAULT_THREADS; nonnative SAM/BAM uses SamWriter's.
	 * The worker-selection predicates use these normalized requests.
	 *
	 * Header lookup/fallback is delegated to the selected SAM/BAM backend, whose
	 * behavior may differ when the shared input header is absent. No header copying
	 * or mate conversion is performed here; read selection follows the backend.
	 * 
	 * @param ffout Output file format
	 * @param qf Optional QUAL filename, used only when ffout is FASTA
	 * @param writeR1 Select first entries according to the backend's Read conversion
	 * @param writeR2 Select second mates according to the backend's Read conversion
	 * @param threads Backend thread request; not a total thread or compression limit
	 * @param header Borrowed SAM/BAM header lines, or null for backend selection
	 * @param useSharedHeader Request shared input-header lookup for SAM/BAM
	 * @return Appropriate Writer implementation, or null if ffout is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Writer makeWriter(FileFormat ffout, String qf, boolean writeR1,
		boolean writeR2, int threads, ArrayList<byte[]> header, boolean useSharedHeader){
//		System.err.println("makeWriter "+ffout+", "+qf);
		if(ffout==null){
			return null;
		}else if(ffout.fastq() || (ffout.fasta() && qf==null)){
			threads=(threads<0 ? FastqWriter.DEFAULT_THREADS : threads);
			boolean fa=ffout.fasta();
			if(threads>1 && Shared.threads()>=8 && !Shared.LOW_MEMORY){
				return new FastqWriter(ffout, threads, writeR1, writeR2);
			}else{
				return new FastqWriterST2(ffout, writeR1, writeR2, threads>0 && Shared.threads()>=4, fa ? 3 : 5);
			}
		}else if(ffout.header()){
			return new FastqWriterST2(ffout, writeR1, writeR2);
		}else if(ffout.bam() && ReadWrite.nativeBamOut()){
			return new BamWriter(ffout, threads, header, useSharedHeader, writeR1, writeR2);
		}else if(ffout.samOrBam()){
			threads=(threads<0 ? SamWriter.DEFAULT_THREADS : threads);
			if(threads>1 && Shared.threads()>=8 && !Shared.LOW_MEMORY){
				return new SamWriter(ffout, threads, header, useSharedHeader, writeR1, writeR2);
			}else{
				return new SamWriterST2(ffout, header, useSharedHeader, threads>0, 5, writeR1, writeR2);
			}
		}else if(ffout.scarf()){
			threads=(threads<0 ? FastqWriter.DEFAULT_THREADS : threads);
			return new FastqWriter(ffout, Math.max(1, threads), writeR1, writeR2);
		}else if(ffout.fasta() && qf!=null){
			return new FastaQualWriterST(ffout, qf, writeR1, writeR2);
		}

		throw new RuntimeException("Unsupported file format: "+ffout);
	}

}
