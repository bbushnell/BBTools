package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import map.ObjectIntMap;
import prok.NcrnaTextIO.SequenceInput;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.FloatList;
import structures.IntList;
import structures.ListNum;

/** Exhaustive Quantum support scan shared by model reduction and boundary training.
 * Adapted from template.A_SampleStreamerMT: workers consume shared Read batches,
 * retain primitive result columns locally, then accumulate once after joining.
 * The library is read-only during scanning. Quantum owns its per-call alignment
 * buffers; this engine adds no alignment cache or shared mutable query scratch.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSupport{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaSupport(final Library library_, final float cutoff_){
		library=library_;
		cutoff=cutoff_;
		if(library==null || library.names.length<1 || !(cutoff>=0 && cutoff<=1)){
			throw new IllegalArgumentException("Support scanning requires a nonempty library and cutoff in [0,1]");
		}
	}

	/** Scans every input gene against every model without shortlist filtering.
	 * Caller must not modify the library until this method returns.
	 * @param input FASTA descriptor; paired input is not supported
	 * @param library Consensus names and sense-oriented bases
	 * @param cutoff Inclusive native identity threshold
	 * @param workers Number of alignment workers
	 * @param threadsIn Reader thread hint passed to StreamerFactory
	 * @return Original-library graph including unsupported genes
	 */
	public static Graph scan(final FileFormat input, final Library library,
			final float cutoff, final int workers, final int threadsIn){
		return scanWithStats(input, library, cutoff, workers, threadsIn).graph;
	}

	/** Runs the same scan and retains coarse worker participation/work counters. */
	public static ScanResult scanWithStats(final FileFormat input, final Library library,
			final float cutoff, final int workers, final int threadsIn){
		if(input==null || !input.fasta() || workers<1){
			throw new IllegalArgumentException("Support scan requires FASTA and at least one worker");
		}
		final NcrnaSupport scan=new NcrnaSupport(library, cutoff);
		final Graph graph=scan.process(input, workers, threadsIn);
		return new ScanResult(graph, scan.batches, scan.activeWorkers);
	}

	/** Scores borrowed uppercase in-memory genes through the same alignment worker routine.
	 * Null members means every gene; otherwise indices select an ordered subset.
	 * Gene/model arrays and the index list must remain unchanged until return.
	 * No Read wrappers, copied sequences or temporary FASTA round trip are needed.
	 */
	public static ScanResult scanWithStats(final Library genes, final IntList members,
			final Library library, final float cutoff, final int workers){
		if(genes==null || workers<1 || (members!=null && members.isEmpty())){
			throw new IllegalArgumentException("In-memory support requires genes and a positive worker budget");
		}
		final int count=members==null ? genes.names.length : members.size;
		if(members!=null){
			for(int i=0; i<count; i++){
				final int gene=members.get(i);
				if(gene<0 || gene>=genes.names.length){throw new IllegalArgumentException("Member index is outside the retained corpus: "+gene);}
			}
		}
		final NcrnaSupport scan=new NcrnaSupport(library, cutoff);
		final AtomicInteger next=new AtomicInteger();
		final int threads=Math.min(workers, (count-1)/MEMORY_BATCH+1);
		final ArrayList<ProcessThread> cohort=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			cohort.add(new ProcessThread(genes, members, next, library, cutoff));
		}
		try(NcrnaWorkers owner=new NcrnaWorkers()){
			owner.startAndJoin(cohort, null);
			for(final ProcessThread worker:cohort){scan.accumulate(worker);}
		}
		if(scan.errorState){throw new IllegalStateException("In-memory support scan did not complete");}
		return new ScanResult(scan.columns.graph(library.names, cutoff, false), scan.batches, scan.activeWorkers);
	}

	/** Loads consensus records through the native sequence reader, once per library. */
	public static Library loadLibrary(final String path){
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		if(!ff.fasta()){throw new IllegalArgumentException("Consensus input must be FASTA: "+path);}
		final Streamer st=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<String> headers=new ArrayList<String>();
		final ArrayList<byte[]> bases=new ArrayList<byte[]>();
		try{
			st.start();
			for(ListNum<Read> batch=st.nextList(); batch!=null; batch=st.nextList()){
				for(final Read read:batch){
					validateRead(read);
					names.add(recordName(read.id));
					headers.add(read.id);
					bases.add(read.bases);
				}
			}
		}catch(Throwable failure){
			// Native FASTA close is not cancellation: a failed consumer can leave
			// a producer blocked on queued unread batches. End the process loudly.
			KillSwitch.exceptionKill(failure);
		}
		// exceptionKill halts the VM. Only normal EOF reaches native close, which
		// is safe after all batches have been consumed and cannot mask a read failure.
		if(ReadWrite.closeStream(st)){throw new IllegalStateException("Error reading consensus FASTA: "+path);}
		return new Library(names.toArray(new String[names.size()]), bases.toArray(new byte[bases.size()][]),
			headers.toArray(new String[headers.size()]));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts one shared reader and joins template-style alignment workers. */
	private Graph process(final FileFormat input, final int workers, final int threadsIn){
		assert(workers>0) : "At least one consumer is needed to drain the sequence reader";
		try(NcrnaWorkers cohort=new NcrnaWorkers();
				SequenceInput owner=new SequenceInput(StreamerFactory.makeStreamer(input, null, false, -1, false, true, threadsIn))){
			final Streamer st=owner.stream;
			final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
			for(int i=0; i<workers; i++){
				threads.add(new ProcessThread(st, library, cutoff));
			}
			st.start();
			cohort.startAndJoin(threads, st);
			for(final ProcessThread worker:threads){accumulate(worker);}
		}
		if(errorState){throw new IllegalStateException("Support scan did not complete successfully");}
		return columns.graph(library.names, cutoff, false);
	}

	/** Merges each worker once after the shared cohort helper joins every producer. */
	private void accumulate(final ProcessThread worker){
		errorState|=!worker.success;
		columns.add(worker.columns);
		batches+=worker.batches;
		if(!worker.columns.names.isEmpty()){activeWorkers++;}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Workers           ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns its result columns and each fetched batch; library bytes are borrowed. */
	static final class ProcessThread extends Thread{
		ProcessThread(final Streamer input_, final Library library_, final float cutoff_){
			this(input_, null, null, null, library_, cutoff_);
		}

		ProcessThread(final Library genes_, final IntList members_, final AtomicInteger next_,
				final Library library_, final float cutoff_){
			this(null, genes_, members_, next_, library_, cutoff_);
		}

		private ProcessThread(final Streamer input_, final Library genes_, final IntList members_,
				final AtomicInteger next_, final Library library_, final float cutoff_){
			input=input_;
			genes=genes_;
			members=members_;
			next=next_;
			library=library_;
			cutoff=cutoff_;
		}

		@Override
		public void run(){
			try{
				if(input==null){
					processMemory();
				}else{
					processInner();
				}
				success=true;
			}catch(Throwable failure){
				// A dead consumer can strand the shared producer on its bounded queue.
				KillSwitch.exceptionKill(failure);
			}
		}

		private void processInner(){
			assert(library.bases.length==library.names.length) : "Each aligned consensus needs the same model index in the output axis";
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				batches++;
				for(final Read read:batch){processRead(read);}
			}
		}

		/** Stores only accepted edges and one overall winner, including below cutoff. */
		private void processRead(final Read read){
			validateRead(read);
			processSequence(recordName(read.id), read.bases);
		}

		/** Claims coarse index blocks; all sequence and name storage remains borrowed. */
		private void processMemory(){
			final int count=members==null ? genes.names.length : members.size;
			for(int from=next.getAndAdd(MEMORY_BATCH); from>=0 && from<count; from=next.getAndAdd(MEMORY_BATCH)){
				batches++;
				final int to=(int)Math.min(count, (long)from+MEMORY_BATCH);
				for(int i=from; i<to; i++){
					final int gene=members==null ? i : members.get(i);
					processSequence(genes.names[gene], genes.bases[gene]);
				}
			}
		}

		/** One exhaustive scoring implementation for native streams and retained genes. */
		private void processSequence(final String name, final byte[] bases){
			int best=-1;
			float bestIdentity=-1;
			for(int model=0; model<library.bases.length; model++){
				final float identity=QuantumAligner.alignStatic(bases, library.bases[model], null);
				if(!(identity>=0 && identity<=1)){
					throw new IllegalStateException("Quantum returned an invalid identity for "+name+" and "+library.names[model]);
				}
				if(identity>bestIdentity){
					best=model;
					bestIdentity=identity;
				}
				if(identity>=cutoff){
					columns.models.add(model);
					columns.identities.add(identity);
				}
			}
			columns.finish(name, best, bestIdentity);
		}

		final Streamer input;
		/** Optional retained input and selection, sharing one coarse-block dispatcher. */
		final Library genes;
		final IntList members;
		final AtomicInteger next;
		final Library library;
		final float cutoff;
		final Columns columns=new Columns();
		boolean success=false;
		long batches=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Retained Data         ----------------*/
	/*--------------------------------------------------------------*/

	/** Scientific graph plus counters from the completed scan, not serialized graph state. */
	public static final class ScanResult{
		private ScanResult(final Graph graph_, final long batches_, final int activeWorkers_){
			graph=graph_;
			batches=batches_;
			activeWorkers=activeWorkers_;
			alignments=(long)graph.queryNames.length*graph.modelNames.length;
		}
		/** Complete scientific result, including unsupported queries. */
		public final Graph graph;
		/** Number of native input batches consumed across all workers. */
		public final long batches;
		/** Number of query/model alignments performed by the exhaustive scan. */
		public final long alignments;
		/** Workers that processed at least one query. */
		public final int activeWorkers;
	}

	/** Consensus axes owned by the caller; treat their arrays as read-only. */
	public static final class Library{
		/** Borrows named consensus arrays, using each identifier as its FASTA header. */
		public Library(final String[] names_, final byte[][] bases_){
			this(names_, bases_, names_);
		}

		private Library(final String[] names_, final byte[][] bases_, final String[] headers_){
			if(names_==null || bases_==null || names_.length<1 || names_.length!=bases_.length){
				throw new IllegalArgumentException("A consensus library requires aligned, nonempty name and sequence arrays");
			}
			checkNames(names_);
			for(final byte[] sequence:bases_){
				if(sequence==null || sequence.length==0){throw new IllegalArgumentException("Empty consensus sequence");}
			}
			names=names_;
			bases=bases_;
			headers=headers_;
		}
		/** Unique FASTA identifiers indexed in the same order as bases. */
		public final String[] names;
		/** Borrowed, nonempty sense-oriented sequences; callers preserve them during scans. */
		public final byte[][] bases;
		/** Original FASTA headers, retained once so subset export need not reread input. */
		public final String[] headers;
	}

	/** Sparse support graph with named axes and explicit original/reduced semantics.
	 * Arrays are owned by this graph and must be treated as read-only by consumers.
	 * Original best is over all models; reduced best is over retained supporters.
	 */
	public static final class Graph{
		Graph(final String[] modelNames_, final String[] queryNames_, final int[] starts_,
				final int[] models_, final float[] identities_, final int[] bestModels_,
				final float[] bestIdentities_, final float cutoff_, final boolean reduced_){
			modelNames=modelNames_;
			queryNames=queryNames_;
			starts=starts_;
			models=models_;
			identities=identities_;
			bestModels=bestModels_;
			bestIdentities=bestIdentities_;
			cutoff=cutoff_;
			reduced=reduced_;
			validate();
		}

		/** Filters scored edges and renumbers models into the exported subset's order.
		 * No alignments are performed. Unsupported reduced queries use -1/NaN.
		 */
		public Graph subset(final boolean[] kept){
			if(kept==null || kept.length!=modelNames.length){throw new IllegalArgumentException("Subset flags must match the graph's model axis");}
			final int[] remap=new int[kept.length];
			Arrays.fill(remap, -1);
			final ArrayList<String> names=new ArrayList<String>();
			for(int model=0; model<kept.length; model++){
				if(kept[model]){
					remap[model]=names.size();
					names.add(modelNames[model]);
				}
			}
			if(names.isEmpty()){throw new IllegalArgumentException("A reduced library must retain at least one model");}
			final Columns data=new Columns();
			for(int gene=0; gene<queryNames.length; gene++){
				int best=-1;
				float identity=Float.NaN;
				for(int edge=starts[gene]; edge<starts[gene+1]; edge++){
					final int model=remap[models[edge]];
					if(model<0){continue;}
					data.models.add(model);
					data.identities.add(identities[edge]);
					if(best<0 || identities[edge]>identity){
						best=model;
						identity=identities[edge];
					}
				}
				data.finish(queryNames[gene], best, identity);
			}
			return data.graph(names.toArray(new String[names.size()]), cutoff, true);
		}

		/** Checks the graph's framing, unique axes and best/edge consistency once. */
		private void validate(){
			if(!(cutoff>=0 && cutoff<=1) || modelNames.length<1 || queryNames.length<1 ||
					starts.length!=queryNames.length+1 || models.length!=identities.length ||
					bestModels.length!=queryNames.length || bestIdentities.length!=queryNames.length ||
					starts[0]!=0 || starts[queryNames.length]!=models.length){
				throw new IllegalArgumentException("Malformed support graph dimensions or cutoff");
			}
			checkNames(modelNames);
			checkNames(queryNames);
			final int[] seen=new int[modelNames.length];
			Arrays.fill(seen, -1);
			for(int gene=0; gene<queryNames.length; gene++){
				final int from=starts[gene], to=starts[gene+1];
				if(from<0 || to<from || to>models.length){throw new IllegalArgumentException("Invalid support offsets for "+queryNames[gene]);}
				final int best=bestModels[gene];
				final float identity=bestIdentities[gene];
				if(reduced && from==to){
					if(best!=-1 || !Float.isNaN(identity)){throw new IllegalArgumentException("Unsupported reduced genes require NONE/NA: "+queryNames[gene]);}
					continue;
				}
				if(best<0 || best>=modelNames.length || !(identity>=0 && identity<=1)){
					throw new IllegalArgumentException("Invalid best model or identity for "+queryNames[gene]);
				}
				boolean foundBest=false;
				for(int edge=from; edge<to; edge++){
					final int model=models[edge];
					final float score=identities[edge];
					if(model<0 || model>=modelNames.length || seen[model]==gene || !(score>=cutoff && score<=identity)){
						throw new IllegalArgumentException("Invalid, duplicate or inconsistent support edge for "+queryNames[gene]);
					}
					seen[model]=gene;
					if(model==best && score==identity){foundBest=true;}
				}
				if((from==to)!=(identity<cutoff) || (from<to && !foundBest)){
					throw new IllegalArgumentException("Best model and support list disagree for "+queryNames[gene]);
				}
			}
		}

		/** Unique identifiers for the consensus and input-gene axes. */
		public final String[] modelNames, queryNames;
		/** Half-open edge offsets per query, with a final sentinel equal to models.length. */
		public final int[] starts;
		/** Model index for each retained support edge. */
		public final int[] models;
		/** Winning model per query; unsupported reduced queries use -1. */
		public final int[] bestModels;
		/** Native identity per retained support edge, at least cutoff. */
		public final float[] identities;
		/** Winning identity per query; unsupported reduced queries use NaN. */
		public final float[] bestIdentities;
		/** Inclusive identity threshold for support edges. */
		public final float cutoff;
		/** True when winners were chosen only among retained supporters. */
		public final boolean reduced;
	}

	/** Append-only primitive columns, local to a worker or serial file reader. */
	static final class Columns{
		Columns(){starts.add(0);}

		void finish(final String name, final int best, final float identity){
			assert(models.size==identities.size) : "Every stored support edge needs its own native identity";
			names.add(name);
			bestModels.add(best);
			bestIdentities.add(identity);
			starts.add(models.size);
		}

		void add(final Columns other){
			assert(other.starts.size==other.names.size()+1) : "Worker batches must finish every query before accumulation";
			final int offset=models.size;
			for(int i=1; i<other.starts.size; i++){
				starts.add(Math.addExact(offset, other.starts.get(i)));
			}
			names.addAll(other.names);
			models.addAll(other.models);
			identities.addAll(other.identities);
			bestModels.addAll(other.bestModels);
			bestIdentities.addAll(other.bestIdentities);
		}

		Graph graph(final String[] modelNames, final float cutoff, final boolean reduced){
			assert(names.size()==bestModels.size) : "Query names and winner columns must share one index space";
			return new Graph(modelNames.clone(), names.toArray(new String[names.size()]), starts.toArray(),
				models.toArray(), identities.toArray(), bestModels.toArray(), bestIdentities.toArray(), cutoff, reduced);
		}

		final ArrayList<String> names=new ArrayList<String>();
		final IntList starts=new IntList();
		final IntList models=new IntList();
		final FloatList identities=new FloatList();
		final IntList bestModels=new IntList();
		final FloatList bestIdentities=new FloatList();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the FASTA identifier; full headers must use printable ASCII. */
	static String recordName(final String header){
		if(header==null || header.isEmpty()){throw new IllegalArgumentException("FASTA records require a name");}
		int end=header.length();
		for(int i=0; i<header.length(); i++){
			final char c=header.charAt(i);
			if(c<32 || c>126){throw new IllegalArgumentException("FASTA headers must be printable ASCII");}
			if(c==' ' && end==header.length()){end=i;}
		}
		if(end==0){throw new IllegalArgumentException("FASTA identifier is empty");}
		return end==header.length() ? header : header.substring(0, end);
	}

	/** Detects ambiguous joins at load time; retained names are the map keys. */
	static void checkNames(final String[] names){
		final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);
		for(int i=0; i<names.length; i++){
			final String name=names[i];
			if(!recordName(name).equals(name) || index.contains(name)){
				throw new IllegalArgumentException("Duplicate or invalid graph identifier: "+name);
			}
			index.put(name, i);
		}
	}

	/** Normalizes owned input bases before alignment and rejects missing sequences. */
	private static void validateRead(final Read read){
		if(read==null || read.bases==null || read.bases.length==0 || read.mate!=null){
			throw new IllegalArgumentException("Support input must contain nonempty, unpaired gene sequences");
		}
		if(!read.validated()){read.validate(true);}
		Tools.toUpperCase(read.bases);
	}

	private final Library library;
	private final float cutoff;
	private final Columns columns=new Columns();
	private boolean errorState=false;
	private long batches=0;
	private int activeWorkers=0;
	private static final int MEMORY_BATCH=64;
}
