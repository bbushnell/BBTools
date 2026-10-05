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

/** Lossless branch cuts on the original contigs, retaining repeat copies.
 * Counting and neighbor lookup use the native canonical k-mer tables. Output
 * batches are bounded even when consecutive branch positions yield1-base pieces.
 * @author Brian Bushnell, Nilou */
final class ShredKmer {

	/** Two passes require a regular replayable DNA FASTA/FASTQ file. */
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

	/** Counts once, then streams original contigs in order through bounded output batches. */
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

	/** Stop coordinates are exclusive; adjacent branches deliberately produce short pieces. */
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

	/** Copies each base/quality exactly once into a result record, preserving existing header conventions. */
	private void emit(Read read, int start, int end, int tid){
		assert(start>=0 && start<end && end<=read.length()) : "Every branch piece must be a nonempty original-contig interval";
		final String name=prefix==null ? read.id : prefix+read.numericID;
		final String id=(name==null ? "" : name+"_")+start+"-"+(end-1)+(tid>0 ? "_tid_"+tid : "");
		pending.add(new Read(KillSwitch.copyOfRange(read.bases, start, end),
			read.quality==null ? null : KillSwitch.copyOfRange(read.quality, start, end), id, readsOut++, read.flags));
		basesOut+=end-start;
		if(pending.size()>=BATCH){flush();}
	}

	/** Hand off a fresh native batch; no retained list aliases the writer's queue. */
	private void flush(){
		if(pending.isEmpty()){return;}
		if(writer!=null){writer.addReads(new ListNum<Read>(pending, outputBatch++));}
		pending=new ArrayList<Read>(BATCH);
	}

	/** Rolling input orientation plus native four-way predecessor/successor queries. */
	private static final class Cursor {
		Cursor(int k_, KmerTableSet small_, KmerTableSetU large_){
			k=k_; small=small_; large=large_; word=large==null ? null : new Kmer(k);
			shift=2*(Math.min(k,31)-1); mask=small==null ? 0 : (1L<<(2*k))-1;
			assert((small==null)!=(large==null)) : "Exactly one native table representation serves this k";
		}
		void clear(){valid=0; forward=reverse=0; if(word!=null){word.clear();}}
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

	long readsIn, basesIn, readsOut, basesOut, branches;
	private final int k;
	private final FileFormat input, output;
	private final long maxReads;
	private final String prefix;
	private final boolean fileTid, headerTid, hashOnly;
	private Writer writer;
	private ArrayList<Read> pending=new ArrayList<Read>(BATCH);
	private long outputBatch;
	private static final int BATCH=256;
}
