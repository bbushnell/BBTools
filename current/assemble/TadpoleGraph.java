package assemble;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectIntMap;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import template.ThreadWaiter;
import ukmer.Kmer;

/**
 * Explores read-supported connections between unchanged assembly endpoints.
 * Each source end uses bounded breadth-first search; an encountered contig end
 * terminates that route. Branches are retained, not resolved. One shortest
 * representative per target is recorded, with A/C/G/T tie order. These are
 * possible kmer paths, not evidence that a read spans the whole connection.
 * @author Fischl
 */
final class TadpoleGraph {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Supplies canonical DNA counts while search keys retain exact orientation. */
	interface Counts {
		int count(Kmer kmer);
	}

	/** Creates an explorer without modifying contig sequences, names, or ordering. */
	TadpoleGraph(final ArrayList<Contig> contigs_, final Counts counts_, final int k_,
			final int minDepth_, final int maxDistance_, final int maxStates_){
		assert(contigs_!=null && counts_!=null) : "Graph exploration needs contigs and read-kmer evidence.";
		if(k_<1 || minDepth_<1 || maxDistance_<0 || maxStates_<1){
			throw new IllegalArgumentException("Graph K, depth, and state bound must be positive; distance must be nonnegative.");
		}
		contigs=contigs_;
		counts=counts_;
		k=k_;
		minDepth=minDepth_;
		maxDistance=maxDistance_;
		maxStates=maxStates_;
		indexEnds();
	}

	/** Recognizes the new mode before the existing multi-K assembly dispatcher. */
	static boolean requested(final String[] args){
		assert(args!=null) : "Dispatch requires expanded command-line arguments.";
		for(String arg : args){
			final String s=arg.toLowerCase(java.util.Locale.ROOT);
			if(s.equals("mode=graph") || s.equals("mode=5") || s.startsWith("contigs=")){return true;}
		}
		return false;
	}

	/** Checks physical input/output aliases before any graph output may be opened. */
	static void checkPaths(final String assembly, final ArrayList<String> inputs, final String... outputs){
		assert(assembly!=null) : "Graph input must be specified before checking destinations.";
		final ArrayList<String> used=new ArrayList<String>(inputs);
		used.add(assembly);
		try{
			for(String output : outputs){
				if(output==null){continue;}
				if(FileFormat.isStdio(output)){
					for(String input : used){
						if(input!=null && ((FileFormat.isStdout(output) && FileFormat.isStdout(input))
								|| (FileFormat.isStderr(output) && FileFormat.isStderr(input)))){
							throw new IllegalArgumentException("Repeated graph stream destination: "+output);
						}
					}
				}else{
					final File dst=new File(output).getCanonicalFile();
					for(String input : used){
						if(input==null || FileFormat.isStdio(input)){continue;}
						final File src=new File(input).getCanonicalFile();
						if(src.equals(dst) || (src.exists() && dst.exists() && Files.isSameFile(src.toPath(), dst.toPath()))){
							throw new IllegalArgumentException("Graph output would overwrite an input or another output: "+output);
						}
					}
				}
				used.add(output);
			}
		}catch(IOException e){throw new RuntimeException("Could not verify graph input/output paths.", e);}
	}

	/** Reads complete FASTA records with header trimming and sequence validation disabled. */
	static ArrayList<Contig> readContigs(final String path){
		assert(path!=null) : "A supplied assembly is required for graph-only mode.";
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, validate=Read.VALIDATE_IN_CONSTRUCTOR;
		final boolean force=FASTQ.FORCE_INTERLEAVED, test=FASTQ.TEST_INTERLEAVED;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		final ArrayList<Contig> result=new ArrayList<Contig>();
		Streamer reader=null;
		try{
			final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
			if(!ff.fasta()){throw new IllegalArgumentException("contigs= must be a FASTA assembly.");}
			reader=StreamerFactory.makeStreamer(ff, 0, true, -1, false, true, 0);
			reader.start();
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read r : batch.list){
					if(r.bases==null){throw new IllegalArgumentException("Assembly record has no sequence: "+r.id);}
					result.add(new Contig(r.bases, r.id, result.size()));
				}
			}
		}finally{
			if(reader!=null){reader.close();}
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.VALIDATE_IN_CONSTRUCTOR=validate;
			FASTQ.FORCE_INTERLEAVED=force;
			FASTQ.TEST_INTERLEAVED=test;
		}
		if(reader.errorState()){throw new RuntimeException("Error reading assembly: "+path);}
		return result;
	}

	/** Indexes every incoming end, retaining all contigs that share an anchor. */
	private void indexEnds(){
		final Kmer key=new Kmer(k);
		assert(key.kbig==k) : "Graph index and count table must use the same effective K.";
		for(int i=0; i<contigs.size(); i++){
			final Contig c=contigs.get(i);
			assert(c.id==i) : "Graph contig IDs must preserve assembly input order.";
			for(int side=0; side<2; side++){
				if(!loadEnd(c, side, false, key)){continue;}
				final Key saved=new Key(key.array1().clone());
				IntList list=ends.get(saved);
				if(list==null){list=new IntList(1); ends.put(saved, list);}
				list.add(2*i+side);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Exploration         ----------------*/
	/*--------------------------------------------------------------*/

	/** Processes contigs independently; failures propagate before any graph is written. */
	void explore(final int threads){
		assert(threads>0) : "Graph worker count must be positive.";
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(threads, contigs.size()); i++){workers.add(new Worker(this));}
		ThreadWaiter.startAndWait(workers);
		for(Worker worker : workers){
			if(worker.failure!=null){throw new RuntimeException("Contig graph worker failed.", worker.failure);}
		}
	}

	/** Loads a complete terminal kmer in incoming or outgoing orientation. */
	private boolean loadEnd(final Contig c, final int side, final boolean outward, final Kmer key){
		assert(side==0 || side==1) : "A contig endpoint is left (0) or right (1).";
		key.clear();
		if(c.length()<k){return false;}
		final int start=(side==0 ? 0 : c.length()-k);
		for(int i=start; i<start+k; i++){
			final int base=AminoAcid.baseToNumber[c.bases[i]];
			if(base<0){key.clear(); return false;}
			key.addRightNumeric(base);
		}
		if(outward ? side==0 : side==1){key.rcomp();}
		return true;
	}

	/** One worker owns its reusable BFS storage and all metadata for each claimed contig. */
	private static final class Worker extends Thread {

		/** Allocates scratch only; the endpoint index and counts remain shared/read-only. */
		Worker(final TadpoleGraph graph_){
			graph=graph_;
			key=new Kmer(graph.k);
			next=new Kmer(graph.k);
			probe=new Key(next.array1());
		}

		/** Captures failures so a worker assertion cannot produce a partial successful graph. */
		@Override
		public void run(){
			try{
				for(int id=graph.nextContig.getAndIncrement(); id<graph.contigs.size(); id=graph.nextContig.getAndIncrement()){
					final Contig c=graph.contigs.get(id);
					measureDepth(c);
					c.leftEdges=walk(c, 0);
					c.rightEdges=walk(c, 1);
				}
			}catch(Throwable t){failure=t;}
		}

		/** Counts all valid assembly kmer positions, including those absent from the reads. */
		private void measureDepth(final Contig c){
			assert(c.bases!=null) : "Coverage annotation must not replace missing sequence with a guess.";
			key.clear();
			long sum=0;
			int positions=0, min=Integer.MAX_VALUE, max=0;
			for(byte b : c.bases){
				final int x=AminoAcid.baseToNumber[b];
				if(x<0){key.clear(); continue;}
				key.addRightNumeric(x);
				if(key.len()<graph.k){continue;}
				final int depth=Math.max(0, graph.counts.count(key));
				sum+=depth;
				positions++;
				min=Math.min(min, depth);
				max=Math.max(max, depth);
			}
			c.coverage=(positions==0 ? 0 : (float)((double)sum/positions));
			c.minCov=(positions==0 ? 0 : min);
			c.maxCov=max;
		}

		/** Finds shortest first-end connections, retaining branches instead of selecting one. */
		private ArrayList<Edge> walk(final Contig c, final int side){
			assert(c.id>=0) : "Reachability requires a stable source contig ID.";
			visited.clear();
			used=0;
			limited=false;
			final ArrayList<Edge> edges=new ArrayList<Edge>();
			if(!graph.loadEnd(c, side, true, key)){return edges;}
			final int seedDepth=Math.max(0, graph.counts.count(key));
			if(seedDepth<graph.minDepth){return edges;}
			add(key, 0, seedDepth, seedDepth, seedDepth, seedDepth);
			for(int head=0; head<used; head++){
				final WalkState state=states.get(head);
				final IntList targets=graph.ends.get(state.key);
				boolean terminal=false;
				if(targets!=null){
					for(int i=0; i<targets.size; i++){
						final int target=targets.get(i);
						if(state.distance==0 && target==2*c.id+side){continue;}
						final Edge edge=new Edge(c.id, target/2, state.distance, side+2*(target&1),
								state.firstDepth, null);
						edge.pathMinDepth=state.minDepth;
						edge.pathMaxDepth=state.maxDepth;
						edge.pathMeanDepth=state.sumDepth/(double)(state.distance+1);
						edges.add(edge);
						terminal=true;
					}
				}
				if(terminal){continue;}
				key.setFrom(state.key.words);
				key.len=graph.k;
				for(int base=0; base<4; base++){
					next.setFrom(key);
					next.addRightNumeric(base);
					final int depth=Math.max(0, graph.counts.count(next));
					if(depth<graph.minDepth){continue;}
					probe.words=next.array1();
					final int old=visited.get(probe);
					//A positive walk may return to its own endpoint (including a palindromic seed).
					if(old==0 && graph.ends.containsKey(probe)){
						if(state.distance<graph.maxDistance){
							addCycleEdges(c, side, state, depth, graph.ends.get(probe), edges);
						}else{limited=true;}
						continue;
					}
					if(old>=0){continue;}
					if(state.distance>=graph.maxDistance || used>=graph.maxStates){limited=true; continue;}
					add(next, state.distance+1, Math.min(state.minDepth, depth), Math.max(state.maxDepth, depth),
							state.sumDepth+depth, state.distance==0 ? depth : state.firstDepth);
				}
			}
			if(side==0){c.graphLeftLimited=limited;}
			else{c.graphRightLimited=limited;}
			return edges;
		}

		/** Records one shortest positive return to a seed endpoint rather than dropping its cycle. */
		private void addCycleEdges(final Contig c, final int side, final WalkState state, final int depth,
				final IntList targets, final ArrayList<Edge> edges){
			assert(targets!=null) : "Only an indexed seed endpoint can terminate a return path.";
			for(int i=0; i<targets.size; i++){
				final int target=targets.get(i);
				boolean present=false;
				for(Edge e : edges){
					if(e.destination==target/2 && e.orientation==side+2*(target&1)){present=true; break;}
				}
				if(present){continue;}
				final Edge e=new Edge(c.id, target/2, state.distance+1, side+2*(target&1),
						state.distance==0 ? depth : state.firstDepth, null);
				e.pathMinDepth=Math.min(state.minDepth, depth);
				e.pathMaxDepth=Math.max(state.maxDepth, depth);
				e.pathMeanDepth=(state.sumDepth+depth)/(double)(state.distance+2);
				edges.add(e);
			}
		}

		/** Reuses retained state objects only after the previous endpoint's map was cleared. */
		private void add(final Kmer kmer, final int distance, final int min, final int max, final long sum, final int first){
			assert(used<graph.maxStates) : "BFS may not exceed its explicitly reported state budget.";
			if(used==states.size()){states.add(new WalkState(kmer.array1().length));}
			final WalkState state=states.get(used);
			System.arraycopy(kmer.array1(), 0, state.key.words, 0, state.key.words.length);
			state.distance=distance;
			state.minDepth=min;
			state.maxDepth=max;
			state.sumDepth=sum;
			state.firstDepth=first;
			visited.put(state.key, used++);
		}

		final TadpoleGraph graph;
		final Kmer key, next;
		final Key probe;
		final ObjectIntMap<Key> visited=new ObjectIntMap<Key>(Key.class);
		final ArrayList<WalkState> states=new ArrayList<WalkState>();
		int used;
		boolean limited;
		Throwable failure;
	}

	/** Exact oriented packed words; canonical hashes alone are insufficient for endpoints. */
	private static final class Key {
		/** Retains the supplied array; callers copy it unless this is the lookup-only probe. */
		Key(final long[] words_){words=words_;}
		@Override
		public int hashCode(){return Arrays.hashCode(words);}
		@Override
		public boolean equals(final Object other){
			return other instanceof Key && Arrays.equals(words, ((Key)other).words);
		}
		long[] words;
	}

	/** Reusable queue slot; depths include the source anchor and each subsequent kmer. */
	private static final class WalkState {
		/** Allocates one exact key buffer for this retained queue slot. */
		WalkState(final int words){key=new Key(new long[words]);}
		final Key key;
		int distance, minDepth, maxDepth, firstDepth;
		long sumDepth;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes graph-only DOT; never invokes graph simplification or contig output. */
	void writeDot(final String path, final boolean pretty){
		assert(path!=null) : "DOT output requires a destination.";
		final ByteStreamWriter writer=TadpoleDot.open(path, pretty, k);
		final ByteBuilder bb=new ByteBuilder();
		final double logLengthPerInch=TadpoleDot.logLengthPerInch(contigs);
		for(Contig c : contigs){
			TadpoleDot.node(c, 0, k, pretty, true, logLengthPerInch, bb);
			for(Edge e : c.leftEdges){bb.tab(); TadpoleDot.edge(e, k, pretty, bb);}
			for(Edge e : c.rightEdges){bb.tab(); TadpoleDot.edge(e, k, pretty, bb);}
			writer.print(bb);
			bb.clear();
		}
		writer.print("}\n");
		if(writer.poisonAndWait()){throw new RuntimeException("Error writing graph: "+path);}
	}

	/** Writes input sequences unchanged; names are tags since GFA IDs cannot contain spaces. */
	void writeGfa(final String path){
		assert(path!=null) : "GFA output requires a destination.";
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.GFA, null, true,
				Tadpole.overwrite, Tadpole.append, false);
		final ByteStreamWriter writer=new ByteStreamWriter(ff);
		writer.start();
		final ByteBuilder bb=new ByteBuilder();
		bb.append("H\tVN:Z:1.0\tTS:Z:Tadpole-graph\tKS:i:").append(k).nl();
		for(Contig c : contigs){
			bb.append("S\tcontig_").append(c.id).tab().append(c.bases);
			bb.append("\tLN:i:").append(c.length()).append("\tDP:f:").append(c.coverage, 2);
			bb.append("\tSN:Z:");
			TadpoleDot.escape(c.name, bb);
			bb.append("\tSL:i:").append(c.graphLeftLimited ? 1 : 0);
			bb.append("\tSR:i:").append(c.graphRightLimited ? 1 : 0).nl();
			for(Edge e : c.leftEdges){gfaEdge(e, bb);}
			for(Edge e : c.rightEdges){gfaEdge(e, bb);}
			writer.print(bb);
			bb.clear();
		}
		//An empty assembly still has a valid GFA header.
		if(bb.length()>0){writer.print(bb);}
		if(writer.poisonAndWait()){throw new RuntimeException("Error writing graph: "+path);}
	}

	/** Labels arbitrary traversed gaps with unknown overlap, not a fabricated overlap CIGAR. */
	private void gfaEdge(final Edge e, final ByteBuilder bb){
		assert(e.pathMinDepth>=0) : "Graph-only links require measured representative path depth.";
		bb.append("L\tcontig_").append(e.origin).tab().append(e.sourceRight() ? '+' : '-');
		bb.append("\tcontig_").append(e.destination).tab().append(e.destRight() ? '-' : '+');
		bb.append("\t*\tKS:i:").append(k).append("\tEC:i:").append(e.depth);
		bb.append("\tEL:i:").append(e.length).append("\tMN:i:").append(e.pathMinDepth);
		bb.append("\tMD:f:").append(e.pathMeanDepth, 2).append("\tMX:i:").append(e.pathMaxDepth).nl();
	}

	/** Counts directed connections and bounded endpoints for an honest summary. */
	void report(){
		long edges=0, limited=0;
		for(Contig c : contigs){
			edges+=c.leftEdges.size()+c.rightEdges.size();
			limited+=(c.graphLeftLimited ? 1 : 0)+(c.graphRightLimited ? 1 : 0);
		}
		Tadpole.outstream.println("Graph-only: "+contigs.size()+" unchanged contigs; "+edges+
				" directed connections; "+limited+" search-limited endpoints; k="+k+".");
		Tadpole.outstream.println("Connections are first-end reachability, not spanning-read support. Depths describe one shortest representative path.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final ArrayList<Contig> contigs;
	final Counts counts;
	final int k, minDepth, maxDistance, maxStates;
	private final HashMap<Key, IntList> ends=new HashMap<Key, IntList>();
	private final AtomicInteger nextContig=new AtomicInteger();
}
