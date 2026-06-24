package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import structures.ListNum;

/**
 * Factory for creating appropriate Streamer implementations based on file format.
 * Handles single files, paired files, and interleaved formats.
 * Supports FASTQ, SAM/BAM, and FASTA.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class StreamerFactory {
	
	/*--------------------------------------------------------------*/
	/*----------------           Legacy             ----------------*/
	/*--------------------------------------------------------------*/

	public static Streamer getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, int threads){
		return makeStreamer(ff1, ff2, null, null, true, maxReads, keepSamHeader, true, threads);
	}

	public static Streamer getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2, int threads){
		return makeStreamer(ff1, ff2, qf1, qf2, true, maxReads, keepSamHeader, true, threads);
	}
	
	public static Streamer makeSamOrBamStreamer(String fname, int threads, boolean saveHeader, 
			boolean ordered, long maxReads, boolean makeReads) {
		return makeSamOrBamStreamer(FileFormat.testInput(fname, FileFormat.SAM, null, true, false), threads, saveHeader, ordered, maxReads, makeReads);
	}
	
	public static Streamer makeSamOrBamStreamer(FileFormat ffin, int threads, boolean saveHeader, 
			boolean ordered, long maxReads, boolean makeReads) {
		return makeStreamer(ffin, 0, ordered, maxReads, saveHeader, makeReads, threads);
	}
	
	public static synchronized ArrayList<byte[]> loadSharedHeader(String s){
		FileFormat ff=FileFormat.testInput(s, FileFormat.SAM, null, false, false);
		return loadSharedHeader(ff);
	}
	
	public static synchronized ArrayList<byte[]> loadSharedHeader(FileFormat ff){
		//comprehension: maxReads=1 caps the scan — start() loads the SAM header (SamReadInputStream.getSharedHeader), and the
		//drain consumes only ~1 read's lines, not the whole file, so loading the shared header is O(header), not O(file).
		Streamer st=makeSamOrBamStreamer(ff, -1, true, true, 1, false);
		st.start();
		while(st.nextLines()!=null) {}
		ReadWrite.closeStream(st);
		ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Twin Files           ----------------*/
	/*--------------------------------------------------------------*/

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
		//[stream/StreamerFactory#002 LOW/latent] if ff1==null but ff2!=null (R2 without R1 — invalid, normally unreachable),
		//s1=null and this builds PairStreamer(null, s2). The reachable contract is ff1 present, ff2 optional; latent on bad input.
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
	
	public static Streamer makeStreamer(FileFormat ff1, FileFormat ff2, String qf1, String qf2,
		boolean ordered, long maxReads, boolean saveHeader, boolean makeReads, int threads){
		Streamer s1=makeStreamer(ff1, qf1, 0, ordered || ff2!=null, maxReads, saveHeader, makeReads, threads);
		Streamer s2=makeStreamer(ff2, qf2, 1, true, maxReads, saveHeader, makeReads, threads);
		return s2==null ? s1 : new PairStreamer(s1, s2);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Single File          ----------------*/
	/*--------------------------------------------------------------*/

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
		return makeStreamer(ff, null, pairnum, ordered, maxReads, saveHeader, makeReads, threads);
	}
	
	/**
	 * Creates a Streamer for a single input file with full configuration.
	 * 
	 * @param ff Input file format
	 * @param qf Qual file, for legacy support
	 * @param pairnum 0 for R1 or unpaired, 1 for R2
	 * @param ordered True to maintain input order in output
	 * @param maxReads Maximum reads to process, or -1 for unlimited
	 * @param saveHeader True to preserve SAM/BAM header information
	 * @param makeReads True to convert SamLines to Read objects (SAM/BAM only)
	 * @param threads Worker threads; -1 for auto, 0 for singlethreaded (zero workers)
	 * @return Appropriate Streamer implementation, or null if ff is null
	 * @throws RuntimeException if file format is unsupported
	 */
	public static Streamer makeStreamer(FileFormat ff, String qf, int pairnum, 
			boolean ordered, long maxReads, boolean saveHeader, boolean makeReads, int threads){
		if(ff==null){
			return null;
			
		}else if(ff.fastq()){
			threads=(threads<0 ? FastqStreamer.DEFAULT_THREADS : threads);
			if(Shared.threads()>8 && threads>1) {
				return new FastqStreamer(ff, threads, pairnum, maxReads);
			}else {
//				return new FastqStreamerST(ff, pairnum, maxReads);
				return new FastqScanStreamer(ff, pairnum, maxReads);
			}
			
		}else if(ff.fasta() && qf==null){
			threads=(threads<0 ? FastaStreamer.DEFAULT_THREADS : threads);
			if(threads==0 || Shared.threads()<4 || Shared.LOW_MEMORY) {
				if(FASTA_STREAMER_2 && Shared.SIMD && !ff.interleaved()) {
					return new FastaStreamer2ZT(ff, pairnum, maxReads);
				}else {
					return new FastaStreamerZT(ff, pairnum, maxReads);
				}
			}else if(threads==1 || Shared.threads()<8) {
				if(FASTA_STREAMER_2 && Shared.SIMD && !ff.interleaved()) {
					return new FastaStreamer2ST(ff, pairnum, maxReads);
				}else {
					return new FastaStreamerST(ff, pairnum, maxReads);
				}
			}else {
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
			if(Shared.threads()>=4 && threads>1) {
				return new SamStreamer(ff, threads, saveHeader, ordered, maxReads, makeReads);
			}else {
				return new SamStreamerST(ff, saveHeader, maxReads, makeReads);
			}
			
		}else if(ff.gfa()){
			return new GfaStreamerST(ff, pairnum, maxReads);
			
		}else if(ff.scarf()){
			return new ScarfStreamer(ff, pairnum, maxReads);
			
		//comprehension: reached only by a fasta WITH a qual file — it falls through the qf==null fasta branch (above) and the
		//intervening bam/sam/gfa/scarf branches (a fasta matches none), so the order is load-bearing: qf==null fasta MUST precede this.
		}else if(ff.fasta() && qf!=null) {
			threads=(threads<0 ? FastqStreamer.DEFAULT_THREADS : threads);
//			return new FastaQualStreamer(ff, qf, pairnum, maxReads);
			return new FastaQualStreamerZT(ff, qf, pairnum, maxReads);
			
		}
		
		throw new RuntimeException("Unsupported file format: "+ff);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Convenience          ----------------*/
	/*--------------------------------------------------------------*/

	public static ArrayList<Read> getReads(long maxReads, boolean keepSamHeader,
		FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		Streamer st=makeStreamer(ff1, ff2, qf1, qf2, true, maxReads, keepSamHeader, true, 1);
		st.start();
		return getReads(st);
	}

	public static ArrayList<Read> getReads(Streamer st){
		ArrayList<Read> out=new ArrayList<Read>();
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
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
	
	//Generally faster and makes fewer objects
	//It does scan for newlines twice though
	public static boolean FASTA_STREAMER_2=true;
	
}