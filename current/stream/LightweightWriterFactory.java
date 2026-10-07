package stream;

import java.util.ArrayList;

import fileIO.FileFormat;

/**
 * Creates ordinary sequence Writers with a fixed ST or ZT execution policy.
 * Format selection belongs here, independently for each destination. This factory
 * owns no routing, rotation, or per-name buffers. ST uses at most one writer thread
 * per physical output; FASTA+QUAL shares one thread for its two outputs. ZT starts
 * no writer threads. Compression runs on the same thread as the leaf's output.
 *
 * Supports FASTQ, FASTA, HEADER, SCARF, SAM, BAM and separate FASTA+QUAL. Native BAM,
 * raw output, gzip/BGZF and ZIP need no compressor subprocess. Other compression
 * types currently lack a supported lightweight backend and are rejected explicitly.
 * ZIP cannot be reopened for append. An append caller must retain the BAM reference
 * dictionary and order; this factory does not inspect the existing BAM dictionary.
 * Lightweight BAM always receives dictionary/header lines, even when global SAM
 * header-suppression flags are set. The BAM backend owns binary header emission.
 *
 * Instances are immutable; create one after parsing startup options. Output writers
 * have their own ownership and ordering contracts. In particular, supply dense IDs
 * from zero to each new ordered leaf and retain input payloads until consumed.
 * Output paths must name distinct physical files; literal duplicate names are
 * rejected, but path aliases and symbolic links are the caller's responsibility.
 *
 * @author Shinobu
 * @date October 1, 2026
 */
public final class LightweightWriterFactory{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves the static startup policy once, using current configured resources.
	 * @param maximumOutputs Upper bound on simultaneously open physical outputs;
	 * include mates and QUAL files, or use the handle cap when fan-out is unknown
	 */
	public LightweightWriterFactory(final int maximumOutputs){
		this(MultiWriterPolicy.snapshot().select(maximumOutputs), 2);
	}

	/**
	 * Creates an explicit ST or ZT factory with a bounded leaf queue capacity.
	 * Queue capacity bounds batches, not bytes; the caller must bound batch sizes.
	 * @param mode_ Resolved ST or ZT, never AUTO
	 * @param queueCapacity_ At least two queued batches for leaves that use a queue
	 */
	public LightweightWriterFactory(final int mode_, final int queueCapacity_){
		if(mode_!=MultiWriterPolicy.ST && mode_!=MultiWriterPolicy.ZT){
			throw new IllegalArgumentException("Resolve AUTO before creating a writer factory: "+mode_);
		}
		//JobQueue requires capacity>1, including ordered leaves in host-driven mode.
		if(queueCapacity_<2){throw new IllegalArgumentException("Writer queue capacity must be at least two: "+queueCapacity_);}
		mode=mode_;
		queueCapacity=queueCapacity_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates one output or a paired wrapper, without starting output threads.
	 * Each descriptor selects its own format. The first writer selects R1 and the
	 * second selects R2; a single output selects both. Optional QUAL is FASTA-only.
	 * @param first Primary output, or null when no output is requested
	 * @param second Optional separate mate output
	 * @param qual1 Optional numeric QUAL file for the primary FASTA output
	 * @param qual2 Optional numeric QUAL file for the mate FASTA output
	 * @param header Optional SAM/BAM header, with the same dictionary on each reopen
	 * @param useSharedHeader Whether SAM/BAM leaves should use the shared input header
	 * @return Unstarted Writer, or null if all output names are absent
	 */
	public Writer getStream(final FileFormat first, final FileFormat second,
			final String qual1, final String qual2, final ArrayList<byte[]> header,
			final boolean useSharedHeader){
		if(first==null){
			if(second!=null || qual1!=null || qual2!=null){
				throw new IllegalArgumentException("Mate/QUAL outputs require a primary output");
			}
			return null;
		}
		validate(first, qual1);
		if(second!=null){validate(second, qual2);}else if(qual2!=null){
			throw new IllegalArgumentException("Mate QUAL requires a mate FASTA output");
		}
		final String[] names={first.name(), second==null ? null : second.name(), qual1, qual2};
		for(int i=0; i<names.length; i++){
			for(int j=0; j<i; j++){
				if(names[i]!=null && names[i].equals(names[j])){
					throw new IllegalArgumentException("Output filenames must be distinct: "+names[i]);
				}
			}
		}
		final Writer w1=makeWriter(first, qual1, true, second==null, header, useSharedHeader);
		if(second==null){return w1;}
		try{
			final Writer w2=makeWriter(second, qual2, false, true, header, useSharedHeader);
			return new PairedWriter(w1, w2);
		}catch(RuntimeException|Error e){
			try{w1.poisonAndWait();}catch(RuntimeException|Error closeError){e.addSuppressed(closeError);}
			throw e;
		}
	}

	/**
	 * Creates an independently formatted leaf with explicit mate selection.
	 * @param ff Output descriptor
	 * @param qual Optional numeric QUAL filename for FASTA
	 * @param writeR1 Select first mates
	 * @param writeR2 Select second mates
	 * @param header Optional SAM/BAM header
	 * @param useSharedHeader Select the shared input header
	 * @return Unstarted ST or ZT leaf
	 */
	public Writer makeWriter(final FileFormat ff, final String qual, final boolean writeR1,
			final boolean writeR2, final ArrayList<byte[]> header, final boolean useSharedHeader){
		validate(ff, qual);
		if(!writeR1 && !writeR2){throw new IllegalArgumentException("A writer must select at least one mate");}
		final boolean threaded=(mode==MultiWriterPolicy.ST);
		if(ff.samOrBam()){
			return new SamWriterST2(ff, header, useSharedHeader, threaded, queueCapacity, writeR1, writeR2, true);
		}else if(qual!=null){
			return threaded ? new FastaQualWriterST(ff, qual, writeR1, writeR2, true)
				: new FastaQualWriterZT(ff, qual, writeR1, writeR2, true);
		}
		return new FastqWriterST2(ff, writeR1, writeR2, threaded, queueCapacity, true);
	}

	/** Rejects unsupported formats/compression before opening a destination. */
	private static void validate(final FileFormat ff, final String qual){
		if(ff==null || !(ff.fastq() || ff.fasta() || ff.header() || ff.scarf() || ff.samOrBam())){
			throw new IllegalArgumentException("Unsupported lightweight sequence output: "+ff);
		}
		LightweightOutputStream.validate(ff);
		if(qual!=null){
			if(!ff.fasta()){throw new IllegalArgumentException("Separate QUAL output requires FASTA: "+ff);}
			if(qual.equals(ff.name())){throw new IllegalArgumentException("FASTA and QUAL require distinct outputs: "+qual);}
			LightweightOutputStream.validate(qualFormat(ff, qual));
		}
	}

	/** Creates a QUAL descriptor with the sequence output's file options. */
	static FileFormat qualFormat(final FileFormat ff, final String qual){
		return FileFormat.testOutput(qual, FileFormat.QUAL, null, false, ff.overwrite(), ff.append(), false);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolved ST or ZT mode, unchanged by later static configuration updates. */
	public final int mode;
	/** Requested leaf queue capacity, in batches; FASTA+QUAL ST uses its fixed queue. */
	public final int queueCapacity;
}
