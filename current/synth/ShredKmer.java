package synth;

import java.io.File;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Partitions original contigs after branch k-mers without collapsing repeat copies.
 * A branch has multiple present predecessors or successors in the native canonical
 * k-mer tables; occurrence multiplicity alone does not create a branch. Exact or
 * hash-only large-k storage follows the requested native-table representation.
 * Output batches are bounded even when consecutive branch positions yield1-base pieces.
 *
 * Use each helper once, from one caller. The input must remain unchanged between
 * counting and the second reading pass; matching read/base totals do not compare
 * sequence contents. Configure shared input and k-mer settings before use.
 * @author Brian Bushnell, Nilou */
final class ShredKmer{

	/** Configures two-pass branch shredding and checks format/encoding preconditions.
	 * Requires one regular replayable DNA FASTA/FASTQ input, whole-contig FASTA
	 * settings and unmasked k-mers. Paired reads are rejected during processing.
	 * @param k_ Positive k-mer length, represented exactly by the native layout
	 * @param hashOnly_ Use hash-only keys for k greater than31; ignored for smaller k
	 * @param input_ Nonnull input descriptor; stdin and comma-separated names are unsupported
	 * @param output_ FASTA/FASTQ output descriptor, or null to count generated pieces without writing
	 * @param maxReads_ Input read limit applied to both passes; negative means unlimited
	 * @param prefix_ Optional header prefix followed by the input numeric ID; null retains the input ID
	 * @param fileTid_ Parse the filename's taxonomic ID; header-derived IDs take precedence
	 * @param headerTid_ Parse IDs from each header instead of using the filename ID
	 */
	ShredKmer(int k_, boolean hashOnly_, FileFormat input_, FileFormat output_, long maxReads_, String prefix_, boolean fileTid_, boolean headerTid_){
		k=k_; hashOnly=hashOnly_; input=input_; output=output_; maxReads=maxReads_;
		prefix=prefix_; fileTid=fileTid_; headerTid=headerTid_;
		if(k<1 || Shared.AMINO_IN || input==null || input.stdio() || !new File(input.name()).isFile()
				|| (!input.fasta() && !input.fastq())){
			throw new IllegalArgumentException("k= branch shredding requires a regular replayable FASTA/FASTQ input file and positive k");
		}
		if(input.name().indexOf(',')>=0){throw new IllegalArgumentException("k= branch mode requires one input filename without commas");}
		if(output!=null && !output.fasta() && !output.fastq()){
			throw new IllegalArgumentException("k= branch shredding writes FASTA or FASTQ");
		}
		if(FastaReadInputStream.SPLIT_READS || FastaReadInputStream.MIN_READ_LEN!=1){
			throw new IllegalArgumentException("Branch mode requires whole contigs: fastasplit=f and fastaminlen=1");
		}
		if(AbstractKmerTableSet.MASK_CORE || AbstractKmerTableSet.MASK_MIDDLE || Kmer.MASK_CORE){
			throw new IllegalArgumentException("Branch mode requires exact unmasked k-mers");
		}
		if(k>31 && Kmer.getKbig(k)!=k){throw new IllegalArgumentException("Requested k requires exact packed k-mers; enable packed layout");}
	}

	/** Counts k-mers, then streams the original contigs in order through bounded batches.
	 * Checks normal input/output status and read/base conservation after finalizing output.
	 * On an incomplete streaming pass, requests writer abandonment before closing the
	 * reader and clearing the tables. Abandonment does not guarantee drained output or
	 * closed output handles; a thrown cleanup failure can replace an earlier failure.
	 * Counts describe generated pieces even when no output writer was requested.
	 */
	void process(){
		final String[] args={"in="+input.name(), "k="+k, "rcomp=t", "prefilter=f", "prealloc=f", "minprob=0", "qtrim=f",
			"minavgquality=0", "reads="+maxReads, "hashonly="+(k>31 && hashOnly)};
		final KmerTableSet small=k<=31 ? new KmerTableSet(args, 12) : null;
		final KmerTableSetU large=k>31 ? new KmerTableSetU(args, 0) : null;
		final AbstractKmerTableSet tables=small==null ? large : small;
		if(large!=null && large.hashOnly()!=hashOnly){throw new IllegalStateException("Native table did not select the requested key representation");}
		tables.process(new Timer());
		if(tables.kbig()!=k){throw new IllegalStateException("Native table changed the requested k");}
		final Cursor cursor=new Cursor(k, small, large);
		final Streamer reader=StreamerFactory.makeStreamer(input, null, true, maxReads, false, true, 1);
		writer=WriterFactory.makeWriter(output, null, 1, null, false);
		boolean complete=false;
		try{
			reader.start(); if(writer!=null){writer.start();}
			final int tid=fileTid ? bin.BinObject.parseTaxID(input.name()) : -1;
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read:batch){
					if(read.mate!=null){throw new IllegalArgumentException("Branch shredding requires unpaired contigs");}
					readsIn++; basesIn+=read.length();
					cut(read, cursor, headerTid ? bin.BinObject.parseTaxID(read.id) : tid);
				}
			}
			flush();
			if(writer!=null && writer.poisonAndWait()){throw new IllegalStateException("Branch shred output failed");}
			if(reader.errorState()){throw new IllegalStateException("Branch shred input failed");}
			if(readsIn!=tables.readsIn || basesIn!=tables.basesIn || basesOut!=basesIn){
				throw new IllegalStateException("Branch shredding must preserve the same input on both passes and every base: counted="
					+tables.basesIn+", read="+basesIn+", emitted="+basesOut);
			}
			if(writer!=null && (writer.basesWritten()!=basesOut || writer.readsWritten()!=readsOut)){
				throw new IllegalStateException("Writer output disagrees with emitted branch pieces");
			}
			complete=true;
		}finally{
			if(!complete && writer!=null){writer.finishError();}
			final boolean inputError=ReadWrite.closeStream(reader); tables.clear();
			if(complete && inputError){throw new IllegalStateException("Branch shred input close failed");}
		}
	}

	/** Partitions one contig into half-open intervals, including ambiguous bases and short tails.
	 * A branch increments the branch-site counter even at the final base, where no
	 * additional cut is needed. Adjacent branch sites deliberately produce short pieces.
	 */
	private void cut(Read read, Cursor cursor, int tid){
		cursor.clear();
		final long before=basesOut;
		int start=0;
		for(int i=0; i<read.length(); i++){
			if(cursor.add(read.bases[i])){
				branches++;
				if(i+1<read.length()){emit(read, start, i+1, tid); start=i+1;}
			}
		}
		if(start<read.length()){emit(read, start, read.length(), tid);}
		assert(basesOut-before==read.length()) : "Half-open branch fragments must partition each original contig without gaps or overlaps";
	}

	/** Copies one nonempty interval and its qualities into a result record.
	 * Coordinates are half-open internally and inclusive in the output header.
	 * Positive taxonomic IDs are appended; no tag is added for zero or negative IDs.
	 * The generated-record counter supplies the output numeric ID before incrementing.
	 */
	private void emit(Read read, int start, int end, int tid){
		assert(start>=0 && start<end && end<=read.length()) : "Every branch piece must be a nonempty original-contig interval";
		final String name=prefix==null ? read.id : prefix+read.numericID;
		final String id=(name==null ? "" : name+"_")+start+"-"+(end-1)+(tid>0 ? "_tid_"+tid : "");
		pending.add(new Read(KillSwitch.copyOfRange(read.bases, start, end),
			read.quality==null ? null : KillSwitch.copyOfRange(read.quality, start, end), id, readsOut++, read.flags));
		basesOut+=end-start;
		if(pending.size()>=BATCH){flush();}
	}

	/** Transfers the pending batch to the writer, or discards it when no output was requested.
	 * Installs a fresh list so later pieces cannot mutate a queued batch.
	 */
	private void flush(){
		if(pending.isEmpty()){return;}
		if(writer!=null){writer.addReads(new ListNum<Read>(pending, outputBatch++));}
		pending=new ArrayList<Read>(BATCH);
	}

	/** Tracks input orientation and valid-prefix length for native four-way neighbor queries.
	 * Small k uses paired forward/reverse-complement longs; large k uses a mutable
	 * Kmer. Undefined bases clear the rolling state. Query return values are ignored;
	 * branching depends on how many of the four populated counts are positive.
	 */
	private static final class Cursor{
		/** Binds exactly one table representation and prepares its rolling encoding. */
		Cursor(int k_, KmerTableSet small_, KmerTableSetU large_){
			k=k_; small=small_; large=large_; word=large==null ? null : new Kmer(k);
			shift=2*(Math.min(k,31)-1); mask=small==null ? 0 : (1L<<(2*k))-1;
			assert((small==null)!=(large==null)) : "Exactly one native table representation serves this k";
		}
		/** Resets the run of valid bases at a contig boundary or undefined base. */
		void clear(){valid=0; forward=reverse=0; if(word!=null){word.clear();}}
		/** Adds an ASCII DNA symbol; returns true only for a complete k-mer with a branch. */
		boolean add(byte base){
			final int x=AminoAcid.baseToNumber[base];
			if(x<0){clear(); return false;}
			valid=Math.min(k, valid+1);
			if(small!=null){
				forward=((forward<<2)|x)&mask;
				reverse=(reverse>>>2)|((long)(3-x)<<shift);
			}else{word.addRightNumeric(x);}
			if(valid<k){return false;}
			if(small!=null){
				small.fillRightCounts(forward, reverse, counts, mask, shift);
				if(branched(counts)){return true;}
				small.fillLeftCounts(forward, reverse, counts, mask, shift);
			}else{
				// Neighbor queries restore the word but advance its len counter while rolling.
				// Track valid bases separately and bound len before each full-word query.
				word.len=k;
				large.fillRightCounts(word, counts);
				if(branched(counts)){return true;}
				large.fillLeftCounts(word, counts);
			}
			return branched(counts);
		}
		/** Tests distinct present extensions, irrespective of their occurrence multiplicity. */
		private static boolean branched(int[] counts){
			assert(counts.length==4) : "Native DNA neighbor enumeration returns one count per A/C/G/T extension";
			int present=0;
			for(int count:counts){if(count>0 && ++present>1){return true;}}
			return false;
		}
		final int k, shift;
		final long mask;
		final KmerTableSet small;
		final KmerTableSetU large;
		final Kmer word;
		final int[] counts=new int[4];
		int valid;
		long forward, reverse;
	}

	/** Second-pass input totals, generated piece/base totals, and branch-site occurrences. */
	long readsIn, basesIn, readsOut, basesOut, branches;
	/** Requested native k-mer length. */
	private final int k;
	/** Caller-supplied descriptors; output may be null. */
	private final FileFormat input, output;
	/** Same read limit for counting and original-contig streaming. */
	private final long maxReads;
	/** Optional replacement for the original header before coordinate suffixes. */
	private final String prefix;
	/** Header ID takes precedence over filename ID; hashOnly applies only to large k. */
	private final boolean fileTid, headerTid, hashOnly;
	/** Writer created for process(), or null when only generating counts. */
	private Writer writer;
	/** Current batch, replaced after each handoff or discard. */
	private ArrayList<Read> pending=new ArrayList<Read>(BATCH);
	/** Dense zero-based batch ID advanced only for actual writer handoffs. */
	private long outputBatch;
	/** Maximum number of generated records retained in one pending batch. */
	private static final int BATCH=256;
}
