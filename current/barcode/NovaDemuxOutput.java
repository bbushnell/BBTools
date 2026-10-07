package barcode;

import java.util.ArrayList;
import java.util.Set;

import fileIO.FileFormat;
import shared.Tools;
import stream.FileWriterTarget;
import stream.LightweightWriterFactory;
import stream.MultiFileWriter;
import stream.MultiWriterPolicy;
import stream.Read;
import stream.Writer;
import structures.ByteBuilder;

/**
 * Opt-in NovaDemux output integration using ordinary ST/ZT Writers.
 * The multi-file container owns matched routing; its targets own serialization.
 * Unmatched output has its own dense batch IDs, continued by residual draining.
 * The physical output cap reserves unmatched files before budgeting matched leaves.
 * This helper starts no dispatcher or compressor threads of its own.
 *
 * @author Shinobu
 * @date October 1, 2026
 */
final class NovaDemuxOutput{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates matched routing and starts the optional unmatched Writer.
	 * @param out1 Primary matched pattern, or null to disable matched output
	 * @param out2 Optional matched mate pattern
	 * @param unmatched1 Optional unmatched primary filename
	 * @param unmatched2 Optional unmatched mate filename
	 * @param extout Explicit unmatched format override, as in the legacy caller
	 * @param overwrite Permit replacing output files
	 * @param append Append to existing output files
	 * @param sharedHeader Use the shared input SAM header for matched outputs
	 * @param maxFiles Total physical cap including unmatched outputs
	 * @param bufferBytes Estimated matched routing budget; leaf buffers and trackers are additional
	 * @param minReads_ Individual-read threshold for matched destinations
	 * @param cardinality Track per-destination paired sequence cardinality */
	NovaDemuxOutput(final String out1, final String out2, final String unmatched1, final String unmatched2,
			final String extout, final boolean overwrite, final boolean append, final boolean sharedHeader,
			final int maxFiles, final long bufferBytes, final long minReads_, final boolean cardinality){
		final int unmatchedFiles=unmatched1==null ? 0 : (unmatched2==null ? 1 : 2);
		final int matchedFiles=out1==null ? 0 : (out2==null ? 1 : 2);
		if(maxFiles<1 || bufferBytes<1 || minReads_<0 || matchedFiles+unmatchedFiles>maxFiles){
			throw new IllegalArgumentException("Multi-writer file cap must fit one matched destination plus unmatched files; budgets must be positive");
		}
		minReads=minReads_;
		leaves=new LightweightWriterFactory(maxFiles);
		targets=out1==null ? null : FileWriterTarget.patternFactory(leaves, out1, out2, null, null,
			FileFormat.FASTQ, overwrite, append, null, sharedHeader);
		matched=targets==null ? null : new MultiFileWriter(targets, maxFiles-unmatchedFiles, bufferBytes,
			256, Math.min(1_000_000L, bufferBytes), minReads, -1, cardinality);
		final boolean ordered=leaves.mode==MultiWriterPolicy.ST;
		final FileFormat first=FileFormat.testOutput(unmatched1, FileFormat.FASTQ, extout, false, overwrite, append, ordered);
		final FileFormat second=FileFormat.testOutput(unmatched2, FileFormat.FASTQ, extout, false, overwrite, append, ordered);
		unmatched=leaves.getStream(first, second, null, null, null, true);
		if(unmatched!=null){unmatched.start();}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Submits one input batch; payloads remain borrowed until output completes.
	 * @param reads All processed reads with their matched destination in Read.obj
	 * @param rejected Ordinary unmatched reads, or all reads in nosplit mode */
	void add(final ArrayList<Read> reads, final ArrayList<Read> rejected){
		if(closed){throw new IllegalStateException("Cannot submit after output completion");}
		assert(reads!=null && rejected!=null) : "NovaDemux always supplies complete input and unmatched lists";
		if(unmatched!=null){unmatched.add(new ArrayList<Read>(rejected), nextUnmatchedId++);}
		if(matched!=null){matched.add(reads);}
	}

	/** Completes matched data, residuals, absent expected files and unmatched output.
	 * Empty outputs are opened through their ordinary format Writer, including BAM metadata.
	 * @param expected Expected assignment names
	 * @param remap Optional filename symbol mapping
	 * @param writeEmpty Create files for expected names never observed
	 * @return Sticky output error indication */
	boolean close(final Set<String> expected, final byte[] remap, final boolean writeEmpty){
		if(closed){return error;}
		assert(!writeEmpty || expected!=null) : "Creating absent expected outputs requires the assignment set";
		final long start=System.nanoTime();
		try{
			if(matched!=null){
				error|=matched.close();
				if(!error && minReads>0){matched.dumpResidual(unmatched, nextUnmatchedId);}
				if(!error && writeEmpty){writeEmpty(expected, remap);}
			}
		}finally{
			if(unmatched!=null){error|=unmatched.poisonAndWait();}
			closed=true;
			closeNanos=System.nanoTime()-start;
		}
		return error;
	}

	/** Creates only never-observed expected destinations; residual names remain observed. */
	private void writeEmpty(final Set<String> expected, final byte[] remap){
		assert(expected!=null && matched!=null && targets!=null) : "Empty matched outputs require expected names and a matched resolver";
		final Set<String> observed=matched.getKeys();
		for(String name : expected){
			final String key=remap==null ? name : Tools.remap(remap, name);
			if(!observed.contains(key)){
				final Writer empty=targets.create(key).open(false);
				empty.start();
				error|=empty.poisonAndWait();
				observed.add(key);
			}
		}
	}

	/** @return Individual matched reads redirected as residuals */
	long residualReads(){return matched==null ? 0 : matched.residualReads();}
	/** @return Matched bases redirected as residuals */
	long residualBases(){return matched==null ? 0 : matched.residualBases();}

	/** Returns the legacy row layout with optional cardinality on emitted destinations. */
	ByteBuilder report(){
		assert(closed) : "NovaDemux reports final destination and residual totals after output completion";
		final ByteBuilder bb=new ByteBuilder();
		if(matched!=null){
			if(minReads>0){bb.append("Residual").tab().append(residualReads()).tab().append(residualBases()).nl();}
			bb.append(matched.report());
		}
		return bb;
	}

	/** Reports modern output counters and finalization time, not legacy retirement phases. */
	String diagnostics(){
		assert(closed) : "Finalization timing is available only after output completion";
		final ByteBuilder bb=new ByteBuilder();
		bb.append("Multi-writer mode:\t").append(leaves.mode==MultiWriterPolicy.ST ? "ST" : "ZT").nl();
		if(matched!=null){
			bb.append("Matched writer opens:\t").append(matched.writersOpened()).nl();
			bb.append("Matched rotations:\t").append(matched.rotations()).nl();
			bb.append("Peak matched files:\t").append(matched.peakOpenOutputs()).nl();
		}
		bb.append("Output finalization:\t").append(closeNanos*1e-9, 3).append(" s").nl();
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Stable ordinary Writer policy and matched target resolver. */
	private final LightweightWriterFactory leaves;
	private final MultiFileWriter.TargetFactory targets;
	/** Matched router and caller-owned unmatched leaf, independently optional. */
	private final MultiFileWriter matched;
	private final Writer unmatched;
	/** Matched eligibility threshold and unmatched batch continuation. */
	private final long minReads;
	private long nextUnmatchedId;
	/** Normal completion state and observed finalization duration. */
	private boolean closed, error;
	private long closeNanos;
}
