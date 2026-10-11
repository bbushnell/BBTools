package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import gff.GffLine;
import ml.CellNet;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Formats a sorted snapshot of borrowed variants as ordered VCF, VAR or GFF output.
 * The calling thread produces batches; the captured formatter workers submit dense
 * batch IDs to one ByteStreamWriter. Header ID 0 precedes data IDs 1..N, including
 * empty output batches. Snapshotting does not copy the mutable variants or registry.
 *
 * The outer omission test only applies the nearby-count rule, with forced bypass.
 * VCF formatting separately computes its filter label; VAR/GFF formatting scores
 * records without applying the complete VarFilter. A nonnull filter is required
 * by all headers, despite the worker's null-filter branch. GFF supports the variant
 * types handled by GffLine.toText. VCF/GFF require an appropriate scaffold registry.
 * Formatting can update variant scoring caches;
 * VAR output also enables global Var.useIdentity. Do not mutate inputs concurrently.
 *
 * Recorded errors accumulate across sequential writes. Headers are built before
 * opening output. On a formatter failure, workers drain queued batches and submit
 * empty jobs for their IDs before the first failure is rethrown after closure.
 * Output from a failed write is incomplete and must not be used. Queue-handoff
 * failures are not covered by this draining protocol.
 * @author Brian Bushnell
 * @contributor Isla
 */
public class VcfWriter{

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Captures dataset statistics and allocates a sorted array of the map's Var references.
	 * Does not start workers or open output. The filter and scaffold registry are borrowed;
	 * headers also derive read length from the supplied read/base totals.
	 * @param varMap_ Quiescent source map, with dataset statistics and scaffold registry
	 * @param filter_ Nonnull output configuration; retained by reference
	 * @param reads_ Reported read count for headers
	 * @param pairs_ Reported paired-read count, not a count of read pairs
	 * @param properPairs_ Reported properly paired-read count
	 * @param bases_ Reported base count used in header statistics
	 * @param ref_ Optional reference path printed in headers
	 * @param trimWhitespace_ Trim scaffold names in VCF header/records
	 * @param sampleName_ VCF sample name
	 */
	public VcfWriter(VarMap varMap_, VarFilter filter_, long reads_,
			long pairs_, long properPairs_, long bases_, String ref_,
			boolean trimWhitespace_, String sampleName_){

		varMap=varMap_;
		filter=filter_;
		trimWhitespace=trimWhitespace_;

		reads=reads_;
		pairs=pairs_;
		properPairs=properPairs_;
		bases=bases_;
		ref=ref_;
		sampleName=sampleName_;

		ploidy=varMap.ploidy;
		properPairRate=varMap.properPairRate;
		pairedInSequencingRate=varMap.pairedInSequencingRate;
		totalQualityAvg=varMap.totalQualityAvg;
		totalMapqAvg=varMap.totalMapqAvg;
		readLengthAvg=varMap.readLengthAvg;
		scafMap=varMap.scafMap;

		threads=Tools.max(1, Shared.threads());

		array=varMap.toArray(true);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a writer and emits VCF, using filename-overwrite defaults.
	 * @return True if normal completion recorded an error; failures may also throw or block
	 */
	public static boolean writeVcf(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeVcfFile(fname);
		return vw.errorState;
	}

	/**
	 * Creates a writer and emits VAR, using filename-overwrite defaults.
	 * @return True if normal completion recorded an error; failures may also throw or block
	 */
	public static boolean writeVar(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeVarFile(fname);
		return vw.errorState;
	}

	/**
	 * Creates a writer and emits GFF, using filename-overwrite defaults.
	 * @return True if normal completion recorded an error; failures may also throw or block
	 */
	public static boolean writeGff(String fname, VarMap varMap, VarFilter filter, boolean trimWhitespace,
			long reads, long pairs, long properPairs, long bases, String ref, String sampleName){
		VcfWriter vw=new VcfWriter(varMap, filter, reads, pairs, properPairs, bases, ref, trimWhitespace, sampleName);
		vw.writeGffFile(fname);
		return vw.errorState;
	}

	/**
	 * Writes VCF, permitting overwrite and subprocesses, without append.
	 * @param fname Output path; format is forced to VCF
	 */
	public void writeVcfFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.VCF, "vcf", true, true, false, true);
		writeFile(ff, VCFMODE);
	}

	/**
	 * Writes VAR, permitting overwrite and subprocesses, without append.
	 * @param fname Output path; format is forced to VAR
	 */
	public void writeVarFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.VAR, "var", true, true, false, true);
		writeFile(ff, VARMODE);
	}

	/**
	 * Writes GFF, permitting overwrite and subprocesses, without append.
	 * @param fname Output path; format is forced to GFF
	 */
	public void writeGffFile(final String fname){
		final FileFormat ff=FileFormat.testOutput(fname, FileFormat.GFF, "gff", true, true, false, true);
		writeFile(ff, GFFMODE);
	}

	/**
	 * Writes VCF using the supplied ordered descriptor and its output policies.
	 * @param ff Descriptor with VCF format and ordered output enabled
	 */
	public void writeVcfFile(final FileFormat ff){
		assert(ff.vcf());
		writeFile(ff, VCFMODE);
	}

	/**
	 * Writes VAR using the supplied ordered descriptor and its output policies.
	 * @param ff Descriptor with VAR format and ordered output enabled
	 */
	public void writeVarFile(final FileFormat ff){
		assert(ff.var()) : "Incorrect file extension: "+ff;
		writeFile(ff, VARMODE);
	}

	/**
	 * Writes GFF using the supplied ordered descriptor and its output policies.
	 * @param ff Descriptor with GFF format and ordered output enabled
	 */
	public void writeGffFile(final FileFormat ff){
		assert(ff.gff());
		writeFile(ff, GFFMODE);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds the header before opening output, then starts the writer and formatters.
	 * After formatter joins, folds the writer's close/error flag into sticky errorState.
	 * Formatter failures are propagated after draining and closure (VAR2-034).
	 * A prior formatter failure prevents reusing this writer for another output.
	 * @param ff Ordered output descriptor
	 * @param mode One of VARMODE, VCFMODE or GFFMODE
	 */
	private void writeFile(final FileFormat ff, final int mode){
		assert(ff.ordered());
		throwFailure();
		final ByteBuilder header=makeHeader(mode);
		final ArrayBlockingQueue<ListNum<Var>> inq=new ArrayBlockingQueue<ListNum<Var>>(threads+1);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		final ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		boolean poisoned=false;
		bsw.start();
		try{
			bsw.add(header, 0);
			spawnThreads(alpt, bsw, inq, mode);
			makeLists(inq);
			waitForFinish(alpt);
			poisoned=poisonWriter(bsw);
		}catch(Throwable t){
			recordFailure(t);
			shutdownFormatters(alpt, inq);
			if(!poisoned){poisonWriter(bsw);}
		}
		throwFailure();
	}

	/**
	 * Builds the selected metadata header before any output resource is opened.
	 * All branches read filter thresholds; null filter is not supported.
	 * @param mode Selected output format
	 * @return Complete header, to be submitted later as ordered job 0
	 */
	private ByteBuilder makeHeader(int mode){
		ByteBuilder bb=new ByteBuilder(1000);
		if(mode==VCFMODE){
			bb.append(VarHelper.toVcfHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					filter.minAlleleDepth, ploidy, reads, pairs, properPairs, bases, ref, scafMap, sampleName, trimWhitespace, net0)).append('\n');
		}else if(mode==VARMODE){
			bb.append(VarHelper.toVarHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					filter.minAlleleDepth, ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
		}else if(mode==GFFMODE){
			bb.append(GffLine.toHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
					ploidy, reads, pairs, properPairs, bases, ref)).append('\n');
		}else{assert(false);}
		return bb;
	}

	/** Keeps the first failure authoritative across workers and later output closure. */
	private synchronized void recordFailure(Throwable t){
		assert(t!=null) : "A formatter/closure failure must retain its original cause";
		if(failure==null){failure=t;}
		errorState=true;
	}

	/** Propagates a recorded failure; successful draining never makes a failed write succeed. */
	private void throwFailure(){
		final Throwable cause=failure;
		if(cause!=null){throw new RuntimeException("VcfWriter failed; output is incomplete", cause);}
	}

	/**
	 * Starts the captured number of formatters; each receives its own net0.copy(false)
	 * when a network is supplied. The shared input queue and output writer are borrowed.
	 * The caller retains the thread list so it can release already-started workers if
	 * construction or startup throws before all formatters are live.
	 * @param alpt Destination for constructed formatters
	 * @param bsw Started ordered byte writer
	 * @param inq Bounded input queue populated by the calling thread
	 * @param mode Selected output format
	 */
	private void spawnThreads(ArrayList<ProcessThread> alpt, ByteStreamWriter bsw, final ArrayBlockingQueue<ListNum<Var>> inq, int mode){

		//Do anything necessary prior to processing

		//Fill a list with ProcessThreads
		for(int i=0; i<threads; i++){
			ProcessThread pt=new ProcessThread(i, bsw, inq, mode);
			pt.net=(net0==null ? null : net0.copy(false));
			alpt.add(pt);
		}
		if(verbose){outstream.println("Spawned threads.");}

		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}

		//Do anything necessary after processing
	}

	/** Finishes the output writer and folds close status into sticky errorState. */
	private boolean poisonWriter(ByteStreamWriter bsw){
		try{
			errorState=bsw.poisonAndWait()|errorState;
		}catch(Throwable t){
			recordFailure(t);
		}
		return true;
	}

	/** Releases any started formatters after producer-side failure. */
	private void shutdownFormatters(ArrayList<ProcessThread> alpt, ArrayBlockingQueue<ListNum<Var>> inq){
		for(ProcessThread pt : alpt){
			if(pt.getState()!=Thread.State.NEW){
				putVars(POISON_VARS, inq);
				break;
			}
		}
		waitForFinish(alpt);
	}

	/**
	 * Joins each formatter and sets errorState if any lacks its success flag.
	 * Interrupted joins print and retry without restoring interruption. There is no
	 * failure cancellation, so a blocked worker can prevent this method from returning.
	 */
	private void waitForFinish(ArrayList<ProcessThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			if(verbose){outstream.println("Waiting for thread "+pt.tid);}
			while(pt.getState()!=Thread.State.TERMINATED && pt.getState()!=Thread.State.NEW){
				try{
					//Attempt a join operation
					pt.join();
				}catch(InterruptedException e){
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}

			allSuccess&=pt.success;
		}

		//Track whether any threads failed
		if(!allSuccess){errorState=true;}
	}

	/**
	 * Partitions the sorted array into batches of at most 200 variants with IDs starting at 1.
	 * Transfers Var references without copying and finishes with the shared terminal.
	 * Header ID 0 is supplied separately. Even an empty array receives the terminal.
	 * @param inq Formatter input queue; a failed consumer can leave puts blocked
	 */
	void makeLists(ArrayBlockingQueue<ListNum<Var>> inq){
		ArrayList<Var> list=new ArrayList<Var>(LIST_SIZE);
		long nextJobID=1;
		for(Var v : array){
			list.add(v);
			if(list.size()>=LIST_SIZE){
				putVars(new ListNum<Var>(list, nextJobID), inq);
				nextJobID++;
				list=new ArrayList<Var>(LIST_SIZE);
			}
		}
		if(list.size()>0){
			putVars(new ListNum<Var>(list, nextJobID), inq);
			nextJobID++;
			list=null;
		}
		if(verbose){outstream.println("tid "+0+" done making var lists.");}
		putVars(POISON_VARS, inq);
		if(verbose){outstream.println("tid "+0+" done poisoning.");}
	}

	/**
	 * Transfers a nonnull batch or terminal by blocking put, retrying interruptions.
	 * The queue receives the original list and variants; the producer must not mutate them.
	 * @param list Data batch or POISON_VARS
	 * @param inq Destination input queue
	 */
	final void putVars(ListNum<Var> list, final ArrayBlockingQueue<ListNum<Var>> inq){
		if(verbose){outstream.println("tid "+0+" putting vlist "+list.id+", size "+list.size());}
		while(list!=null){
			try{
				inq.put(list);
				list=null;
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		if(verbose){outstream.println("tid "+0+" done putting vlist");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats batches and sends copied byte arrays under their original ordered IDs. */
	private class ProcessThread extends Thread{

		/**
		 * Stores a formatter role, borrowed queues/writer and selected output format.
		 */
		public ProcessThread(int tid_, ByteStreamWriter bsw_, ArrayBlockingQueue<ListNum<Var>> inq_, int mode_){
			tid=tid_;
			bsw=bsw_;
			inq=inq_;
			mode=mode_;
		}

		/** Marks success only after formatting and terminal forwarding without a shared failure. */
		@Override
		public void run(){
			makeBytes();
			success=(failure==null);
		}

		/**
		 * Formats each batch into reusable scratch, transferring a copied buffer for every ID,
		 * even when nearby filtering emits no records. VCF computes PASS/filter labels;
		 * VAR/GFF scoring does not perform complete emission filtering. Forced records bypass
		 * the outer nearby-count omission. After transfer, clears and caps scratch at 200000
		 * bytes. After a formatter exception, drops that partial batch and drains remaining
		 * batches as empty ordered jobs. The terminal is still forwarded; failures in the
		 * queue or byte-writer handoff itself are outside this formatter-recovery scope.
		 */
		void makeBytes(){
			if(verbose){outstream.println("tid "+tid+" started makeBytes.");}

			ListNum<Var> list=takeVars();
			final ByteBuilder bb=new ByteBuilder();
			while(list!=null && list.size()>0 && list!=POISON_VARS){

				if(failure==null){
					try{
						formatBatch(list, bb);
					}catch(Throwable t){
						recordFailure(t);
						bb.clear();
					}
				}

				//Every assigned ID is required even when formatting failed (VAR2-034).
				bsw.add(new ByteBuilder(bb.toBytes()), list.id);
				list=takeVars();

				bb.clear();
				bb.shrinkTo(SHRINK_SIZE);
			}
			if(list==POISON_VARS){
				putVars(POISON_VARS, inq);
			}
			if(verbose){outstream.println("tid "+tid+" done making bytes.");}
		}

		/** Appends one data batch; the caller discards partial bytes if formatting throws. */
		private void formatBatch(ListNum<Var> list, ByteBuilder bb){
			assert(list.id>0) : "makeLists reserves ID 0 for the header and assigns positive data IDs";
			for(Var v : list){
				if(v.forced() || filter==null || !filter.failNearby || v.nearbyVarCount<=filter.maxNearbyCount){
					if(mode==VCFMODE){
						v.toVCF(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, filter, net, trimWhitespace);
					}else if(mode==VARMODE){
						v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap, net);
					}else if(mode==GFFMODE){
						GffLine.toText(bb, v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap, net);
					}
					bb.nl();
				}
			}
		}

		/**
		 * Blocks for an input batch or terminal; interruptions print and retry.
		 * @return Original queued batch, never null on normal return
		 */
		final ListNum<Var> takeVars(){
			if(verbose){outstream.println("tid "+tid+" taking vlist");}
			ListNum<Var> list=null;
			while(list==null){
				try{
					list=inq.take();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			if(verbose){outstream.println("tid "+tid+" took vlist "+list.id+", size "+list.size());}
			return list;
		}

		/** Borrowed input queue; data batches have dense positive IDs. */
		final ArrayBlockingQueue<ListNum<Var>> inq;
		/** Formatter identifier used only for diagnostic output. */
		final int tid;
		/** Borrowed started output writer; ownership remains with writeFile. */
		final ByteStreamWriter bsw;
		/** Private formatter network copy, or null for composite scoring. */
		CellNet net;
		/** Selected output format. */
		final int mode;
		/** Set only after normal makeBytes completion; read after join. */
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed network prototype copied per formatter before each write. */
	CellNet net0;
	/**
	 * Sets the borrowed network prototype for subsequent writes; null selects composite
	 * scoring. Must not race a write that is creating worker copies.
	 */
	public void setNet(CellNet n){net0=n;}
	/** Sorted array snapshot containing borrowed mutable variants. */
	final Var[] array;
	/** Ploidy captured from the source map. */
	final int ploidy;
	/** Proper-pair rate captured for scoring and headers. */
	final double properPairRate;
	/** Captured pairing statistic, currently unused by this writer. */
	final double pairedInSequencingRate;
	/** Captured average base quality. */
	final double totalQualityAvg;
	/** Captured average mapping quality. */
	final double totalMapqAvg;
	/** Captured scoring read length; header read length derives from supplied counts. */
	final double readLengthAvg;
	/** Borrowed registry and reference bases used by formatters. */
	final ScafMap scafMap;
	/** Original borrowed map; writing iterates the captured array instead. */
	final VarMap varMap;
	/** Borrowed nonnull thresholds and VCF filter-label configuration. */
	final VarFilter filter;
	/** Whether VCF formatting trims scaffold names. */
	final boolean trimWhitespace;

	/** VCF sample label. */
	final String sampleName;
	/** Reported read count for headers. */
	final long reads;
	/** Reported paired-read count for headers. */
	final long pairs;
	/** Reported properly paired-read count for headers. */
	final long properPairs;
	/** Reported base count for headers. */
	final long bases;
	/** Optional reference path included in headers. */
	final String ref;

	/** Formatter count captured as max(1, Shared.threads()) at construction. */
	final int threads;

	/** Sticky recorded error status across sequential output files. */
	boolean errorState=false;
	/** First formatter/closure failure, retained across writes and read between batches. */
	private volatile Throwable failure;

	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Internal output format selectors. */
	public static final int VARMODE=0, VCFMODE=1, GFFMODE=2;
	/** Shared terminal identity forwarded on normal formatter completion. */
	private static final ListNum<Var> POISON_VARS=new ListNum<Var>(null, Long.MAX_VALUE, true, false);

	/** Enable diagnostic output. */
	private static boolean verbose=false;

	/** Maximum variants per formatting batch. */
	private static final int LIST_SIZE=200;
	/** Maximum retained scratch capacity after each emitted batch. */
	private static final int SHRINK_SIZE=1000*LIST_SIZE;

	/** Borrowed diagnostic stream, never closed here. */
	protected static PrintStream outstream=System.err;
}
