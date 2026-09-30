package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import structures.ListNum;

/**
 * Selects readers for FASTQ, FASTA with optional QUAL, SAM/BAM, GFA and SCARF.
 * Factory methods return unstarted streams; callers start, consume and close them.
 * The loadSharedHeader and getReads convenience methods perform the documented
 * consumption themselves. Format, requested threads and global settings determine
 * the concrete reader; a thread hint is not a total background-thread limit.
 * Two-file factories force ordering and wrap nonnull secondary inputs in PairStreamer.
 * Their inputs must satisfy that adapter's non-interleaved R1/R2 contracts.
 * Limits are forwarded unchanged: selected readers define their counting units
 * (for example, interleaved FASTQ counts pairs). Negative limits mean unlimited.
 *
 * @author Isla
 * @date October 31, 2025
 */
public class StreamerFactory{

	/*--------------------------------------------------------------*/
	/*----------------           Legacy             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Legacy entry point for an ordered, unstarted reader with Read conversion enabled.
	 * @param maxReads Limit forwarded unchanged to each input; units depend on the reader, negative for unlimited
	 * @param keepSamHeader Whether to publish the SAM/BAM header in shared state
	 * @param ff1 Primary input; required when ff2 is present
	 * @param ff2 Optional second input
	 * @param threads Reader-selection hint; negative selects format defaults
	 * @return Unstarted reader, or null when both inputs are null
	 */
	public static Streamer getReadInputStream(final long maxReads, final boolean keepSamHeader,
			final FileFormat ff1, final FileFormat ff2, final int threads){
		return makeStreamer(ff1, ff2, null, null, true, maxReads, keepSamHeader, true, threads);
	}

	/**
	 * Legacy ordered-reader entry point with optional FASTA quality files.
	 * @param maxReads Limit forwarded unchanged to each input; units depend on the reader, negative for unlimited
	 * @param keepSamHeader Whether to publish the SAM/BAM header in shared state
	 * @param ff1 Primary input; required when ff2 is present
	 * @param ff2 Optional second input
	 * @param qf1 Optional QUAL path for primary FASTA input
	 * @param qf2 Optional QUAL path for secondary FASTA input
	 * @param threads Reader-selection hint; negative selects format defaults
	 * @return Unstarted reader with Read conversion enabled, or null for no inputs
	 */
	public static Streamer getReadInputStream(final long maxReads, final boolean keepSamHeader,
			final FileFormat ff1, final FileFormat ff2, final String qf1, final String qf2, final int threads){
		return makeStreamer(ff1, ff2, qf1, qf2, true, maxReads, keepSamHeader, true, threads);
	}

	/**
	 * Resolves an input path using SAM as the default format, then delegates selection.
	 * @param fname Input path
	 * @param threads Reader-selection hint; negative selects format defaults
	 * @param saveHeader Whether to publish the SAM/BAM header in shared state
	 * @param ordered Whether to request ordered output batches
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @param makeReads Whether SAM/BAM readers also construct Read objects
	 * @return Unstarted reader selected from the resolved descriptor
	 */
	public static Streamer makeSamOrBamStreamer(final String fname, final int threads, final boolean saveHeader,
			final boolean ordered, final long maxReads, final boolean makeReads){
		return makeSamOrBamStreamer(FileFormat.testInput(fname, FileFormat.SAM, null, true, false), threads, saveHeader, ordered, maxReads, makeReads);
	}

	/**
	 * Delegates reader selection with pair number zero; does not enforce a SAM/BAM format.
	 * @param ffin Input descriptor, or null
	 * @param threads Reader-selection hint; negative selects format defaults
	 * @param saveHeader Whether to publish the SAM/BAM header in shared state
	 * @param ordered Whether to request ordered output batches
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @param makeReads Whether SAM/BAM readers also construct Read objects
	 * @return Unstarted reader, or null for a null descriptor
	 */
	public static Streamer makeSamOrBamStreamer(final FileFormat ffin, final int threads, final boolean saveHeader,
			final boolean ordered, final long maxReads, final boolean makeReads){
		return makeStreamer(ffin, 0, ordered, maxReads, saveHeader, makeReads, threads);
	}

	/**
	 * Resolves a SAM/BAM path and loads its shared header through the descriptor overload.
	 * @param s SAM/BAM input path
	 * @return Shared header list, without a defensive copy
	 */
	public static synchronized ArrayList<byte[]> loadSharedHeader(final String s){
		final FileFormat ff=FileFormat.testInput(s, FileFormat.SAM, null, false, false);
		return loadSharedHeader(ff);
	}

	/**
	 * Starts a SAM/BAM reader with header retention and a one-record limit, drains its
	 * SamLine output, closes it, then obtains the shared header. Reader buffering and
	 * decompression can read ahead; the limit does not bound physical input to one record.
	 * Synchronization serializes these helpers, not all users of global header state.
	 * @param ff Nonnull SAM/BAM input descriptor
	 * @return Shared header list, without a defensive copy
	 */
	public static synchronized ArrayList<byte[]> loadSharedHeader(final FileFormat ff){
		//Limit record consumption while allowing the selected reader to load and publish its header.
		final Streamer st=makeSamOrBamStreamer(ff, -1, true, true, 1, false);
		st.start();
		while(st.nextLines()!=null){}
		ReadWrite.closeStream(st);
		final ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Twin Files           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an unstarted reader with format-default threads; SAM/BAM header retention and Read conversion are disabled.
	 *
	 * @param ff1 Primary input file (R1 for paired data)
	 * @param ff2 Secondary input file (R2 for paired data), or null
	 * @param ordered True to maintain input order in output
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @return Appropriate Streamer implementation
	 */
	public static Streamer makeStreamer(final FileFormat ff1, final FileFormat ff2,
			final boolean ordered, final long maxReads){
		return makeStreamer(ff1, ff2, ordered, maxReads, false, false);
	}

	/**
	 * Creates an unstarted reader for one or two files with format-default threads.
	 * If ff2 is null, returns a single-file streamer.
	 * If ff2 is non-null, returns a PairStreamer wrapping both files.
	 * Forces ordering when pairing to ensure mate synchronization.
	 *
	 * @param ff1 Primary input file (R1 for paired data)
	 * @param ff2 Secondary input file (R2 for paired data), or null
	 * @param ordered True to maintain input order in output (forced true if ff2!=null)
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @return Appropriate Streamer implementation
	 */
	public static Streamer makeStreamer(final FileFormat ff1, final FileFormat ff2,
		final boolean ordered, final long maxReads, final boolean saveHeader, final boolean makeReads){
		//[stream/StreamerFactory#002 LOW/latent] if ff1==null but ff2!=null (R2 without R1 — invalid, normally unreachable),
		//s1=null and this builds PairStreamer(null, s2). The reachable contract is ff1 present, ff2 optional; latent on bad input.
		final Streamer s1=makeStreamer(ff1, 0, ordered || ff2!=null, maxReads, saveHeader, makeReads);
		final Streamer s2=makeStreamer(ff2, 1, true, maxReads, saveHeader, makeReads);
		return s2==null ? s1 : new PairStreamer(s1, s2);
	}

	/**
	 * Creates unstarted readers with a common thread hint, forcing order for two inputs.
	 * @param ff1 Primary input; required when ff2 is present
	 * @param ff2 Optional non-interleaved R2 input compatible with PairStreamer
	 * @param ordered Whether to request ordered batches; forced true for two inputs
	 * @param maxReads Limit forwarded unchanged to each input; units depend on the reader, negative for unlimited
	 * @param saveHeader Whether to publish SAM/BAM headers in shared state
	 * @param makeReads Whether SAM/BAM readers also construct Read objects
	 * @param threads Reader-selection hint; negative selects format defaults
	 * @return Primary reader, PairStreamer, or null when both inputs are null
	 */
	public static Streamer makeStreamer(final FileFormat ff1, final FileFormat ff2,
		final boolean ordered, final long maxReads, final boolean saveHeader, final boolean makeReads, final int threads){
		final Streamer s1=makeStreamer(ff1, 0, ordered || ff2!=null, maxReads, saveHeader, makeReads, threads);
		final Streamer s2=makeStreamer(ff2, 1, true, maxReads, saveHeader, makeReads, threads);
		return s2==null ? s1 : new PairStreamer(s1, s2);
	}

	/**
	 * Creates unstarted readers with optional FASTA quality files and a common thread hint.
	 * @param ff1 Primary input; required when ff2 is present
	 * @param ff2 Optional non-interleaved R2 input compatible with PairStreamer
	 * @param qf1 Optional QUAL path for primary FASTA input
	 * @param qf2 Optional QUAL path for secondary FASTA input
	 * @param ordered Whether to request ordered batches; forced true for two inputs
	 * @param maxReads Limit forwarded unchanged to each input; units depend on the reader, negative for unlimited
	 * @param saveHeader Whether to publish SAM/BAM headers in shared state
	 * @param makeReads Whether SAM/BAM readers also construct Read objects
	 * @param threads Reader-selection hint; negative selects format defaults
	 * @return Primary reader, PairStreamer, or null when both inputs are null
	 */
	public static Streamer makeStreamer(final FileFormat ff1, final FileFormat ff2, final String qf1, final String qf2,
		final boolean ordered, final long maxReads, final boolean saveHeader, final boolean makeReads, final int threads){
		final Streamer s1=makeStreamer(ff1, qf1, 0, ordered || ff2!=null, maxReads, saveHeader, makeReads, threads);
		final Streamer s2=makeStreamer(ff2, qf2, 1, true, maxReads, saveHeader, makeReads, threads);
		return s2==null ? s1 : new PairStreamer(s1, s2);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Single File          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an unstarted reader with format-default threads; SAM/BAM header retention and Read conversion are disabled.
	 *
	 * @param ff Input file format
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @return Appropriate Streamer implementation, or null if ff is null
	 */
	public static Streamer makeStreamer(final FileFormat ff, final int pairnum,
			final boolean ordered, final long maxReads){
		return makeStreamer(ff, pairnum, ordered, maxReads, false, false);
	}

	/**
	 * Creates an unstarted reader with format-default threads and no separate quality file.
	 *
	 * @param ff Input file format
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @return Appropriate Streamer implementation, or null if ff is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Streamer makeStreamer(final FileFormat ff, final int pairnum,
			final boolean ordered, final long maxReads, final boolean saveHeader, final boolean makeReads){
		return makeStreamer(ff, pairnum, ordered, maxReads, saveHeader, makeReads, -1);
	}

	/**
	 * Creates an unstarted reader with a thread hint and no separate quality file.
	 *
	 * @param ff Input file format
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @param threads Reader-selection hint; negative selects defaults, zero requests a simpler path where supported
	 * @return Appropriate Streamer implementation, or null if ff is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Streamer makeStreamer(final FileFormat ff, final int pairnum,
			final boolean ordered, final long maxReads, final boolean saveHeader, final boolean makeReads, final int threads){
		return makeStreamer(ff, null, pairnum, ordered, maxReads, saveHeader, makeReads, threads);
	}

	/**
	 * Selects an unstarted reader from the format, thread hint and current global settings.
	 * FASTQ uses its parallel reader only with more than eight global threads and a
	 * hint above one. FASTA selection also considers low-memory mode, SIMD and
	 * interleaving; a QUAL path selects the legacy quality-aware reader. Native BAM
	 * uses BamStreamer; remaining SAM/BAM inputs use parallel SAM parsing only with
	 * at least four global threads and a hint above one. GFA and SCARF use their own
	 * readers. Negative hints select format defaults; native BAM resolves its hint
	 * internally. Header retention marks native SAM/BAM input before startup so
	 * shared-header consumers can wait for publication.
	 *
	 * @param ff Input file format
	 * @param qf Optional QUAL path for FASTA; ignored for other formats
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Limit forwarded unchanged; units depend on the reader, negative for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @param threads Reader-selection hint; zero is not a guarantee of zero background threads
	 * @return Appropriate Streamer implementation, or null if ff is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Streamer makeStreamer(final FileFormat ff, final String qf, final int pairnum,
			final boolean ordered, final long maxReads, final boolean saveHeader, final boolean makeReads, int threads){
		if(ff==null){
			return null;

		}else if(ff.fastq()){
			threads=(threads<0 ? FastqStreamer.DEFAULT_THREADS : threads);
			if(Shared.threads()>8 && threads>1){
				return new FastqStreamer(ff, threads, pairnum, maxReads);
			}else{
				return new FastqScanStreamer(ff, pairnum, maxReads);
			}

		}else if(ff.fasta() && qf==null){
			threads=(threads<0 ? FastaStreamer.DEFAULT_THREADS : threads);
			if(threads==0 || Shared.threads()<4 || Shared.LOW_MEMORY){
				if(FASTA_STREAMER_2 && Shared.SIMD && !ff.interleaved()){
					return new FastaStreamer2ZT(ff, pairnum, maxReads);
				}else{
					return new FastaStreamerZT(ff, pairnum, maxReads);
				}
			}else if(threads==1 || Shared.threads()<8){
				if(FASTA_STREAMER_2 && Shared.SIMD && !ff.interleaved()){
					return new FastaStreamer2ST(ff, pairnum, maxReads);
				}else{
					return new FastaStreamerST(ff, pairnum, maxReads);
				}
			}else{
				return new FastaStreamer(ff, threads, pairnum, maxReads);
			}

		}else if(ff.bam() && ReadWrite.nativeBamIn()){
			//A native SAM/BAM input with saveHeader will call setSharedHeader from its worker thread; flag it
			//now (synchronously, on this thread) so getSharedHeader(true) blocks for it rather than racing to
			//null via the #002 gate.  Fixes var2/ScafMap.loadSamHeader race. See SamReadInputStream.markSamInputPresent.
			if(saveHeader){SamReadInputStream.markSamInputPresent();}
			return new BamStreamer(ff, threads, saveHeader, ordered, maxReads, makeReads);

		}else if(ff.samOrBam()){
			if(saveHeader){SamReadInputStream.markSamInputPresent();}
			threads=(threads<0 ? SamStreamer.DEFAULT_THREADS : threads);
			if(Shared.threads()>=4 && threads>1){
				return new SamStreamer(ff, threads, saveHeader, ordered, maxReads, makeReads);
			}else{
				return new SamStreamerST(ff, saveHeader, maxReads, makeReads);
			}

		}else if(ff.gfa()){
			return new GfaStreamerST(ff, pairnum, maxReads);

		}else if(ff.scarf()){
			return new ScarfStreamer(ff, pairnum, maxReads);

		//FASTA with a QUAL path bypasses the qf==null branch and uses the quality-aware reader.
		//Keep the distinction explicit so a general FASTA branch does not discard the quality path.
		}else if(ff.fasta() && qf!=null){
			threads=(threads<0 ? FastqStreamer.DEFAULT_THREADS : threads);
			return new FastaQualStreamerZT(ff, qf, pairnum, maxReads);

		}

		throw new RuntimeException("Unsupported file format: "+ff);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Convenience          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates and starts an ordered reader, materializes all returned reads and closes it.
	 * Uses a thread hint of one and enables SAM/BAM Read conversion.
	 * @param maxReads Limit forwarded unchanged to each input; units depend on the reader, negative for unlimited
	 * @param keepSamHeader Whether to publish the SAM/BAM header in shared state
	 * @param ff1 Nonnull primary input
	 * @param ff2 Optional non-interleaved R2 input compatible with PairStreamer
	 * @param qf1 Optional QUAL path for primary FASTA input
	 * @param qf2 Optional QUAL path for secondary FASTA input
	 * @return Collected read entries; mates remain attached rather than added separately
	 */
	public static ArrayList<Read> getReads(final long maxReads, final boolean keepSamHeader,
		final FileFormat ff1, final FileFormat ff2, final String qf1, final String qf2){
		final Streamer st=makeStreamer(ff1, ff2, qf1, qf2, true, maxReads, keepSamHeader, true, 1);
		st.start();
		return getReads(st);
	}

	/**
	 * Drains and closes an already-started, nonnull reader supporting nextList().
	 * Copies list entries into one result list without copying the Read objects.
	 * After normal exhaustion and close, a reported read error terminates the process
	 * rather than returning partial reference data. Exceptions during iteration
	 * propagate; this method does not provide a finally-based close guarantee.
	 * @param st Started reader; SAM/BAM readers must have Read conversion enabled
	 * @return Collected entries in the order supplied by the reader
	 */
	public static ArrayList<Read> getReads(final Streamer st){
		final ArrayList<Read> out=new ArrayList<Read>();
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			out.addAll(ln.list);
		}
		st.close();
		//[stream/StreamerFactory#001 FIXED 2026-06-21] crash LOUD on a read error instead of returning PARTIAL reads with only a stderr
		//warning: a truncated/corrupt input silently yields fewer reads, and st is already closed so the caller can't re-check errorState().
		//The ~6 callers (ddl SSU loaders, ifa aligners) load this as REFERENCE data -> partial reference = wrong results. KillSwitch.kill
		//aborts non-zero (BBTools contract: crash, never silently wrong). Sibling: ConcurrentReadInputStream.getReads (same fix; AllToAll:223 resolved).
		if(st.errorState()){
			KillSwitch.kill("Error: a read error (corrupt or truncated input) occurred reading "+st.fname()
				+"; aborting rather than returning partial/incorrect read data.");
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Enables version-two FASTA readers when SIMD is enabled and input is not interleaved.
	 * Consulted for the zero/one-worker selection paths; configure before creating readers.
	 * The original implementation rationale was fewer objects and generally higher speed,
	 * at the cost of scanning newlines twice; this is not a universal benchmark guarantee.
	 */
	public static boolean FASTA_STREAMER_2=true;

}
