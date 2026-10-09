package assemble;

import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Opt-in streaming census of reciprocal fusion proposals before safety vetoes.
 * Counts are borrowed from exactly one phase table; reference truth is deferred
 * until after assembly. No selection or merging policy depends on these outputs.
 * @author Fischl
 */
final class FusionJoinCollector implements AutoCloseable {

	/** Opens new sidecars only after checking every input/output physical alias. */
	FusionJoinCollector(final String prefix, final ArrayList<String> inputs,
			final String assembly, final String gfa){
		assert(prefix!=null && assembly!=null) : "Corpus output needs an explicit prefix and assembly path.";
		final String[] paths={prefix+".diagnostic.tsv", prefix+".vectors.tsv",
				prefix+".histograms.tsv", prefix+".trace.tsv"};
		final ArrayList<String> protectedPaths=new ArrayList<String>(inputs);
		protectedPaths.add(gfa);
		TadpoleGraph.checkPaths(assembly, protectedPaths, paths);
		if(!Tools.testOutputFiles(false, false, false, paths)){
			throw new IllegalArgumentException("Fusion sidecars must be new writable files: "+prefix);
		}
		try{
			for(int i=0; i<paths.length; i++){
				writers[i]=new ByteStreamWriter(paths[i], false, false, false);
				writers[i].start();
			}
			FusionJoinDiagnostic.writeDiagnosticHeader(writers[0]);
			FusionJoinDiagnostic.writeVectorHeader(writers[1]);
			writers[2].print("event\tphase_k\tquery_k\tregion\tdepth\tpositions\n");
			writers[3].print("#cohort=reciprocal_pre_guard; labels are deferred; one event per frozen pair\n");
		}catch(RuntimeException | Error e){
			try{close();}catch(RuntimeException | Error closing){e.addSuppressed(closing);}
			throw e;
		}
	}

	/** Binds the current phase's provider; the preceding provider must be released. */
	void beginPhase(final int k_, final int maxOverlap, final boolean graphCounts, final TadpoleGraph.Counts counts){
		assert(measurer==null && counts!=null) : "End the previous phase before borrowing another count table.";
		assert(maxOverlap>=k_) : "FusionPass must describe the same overlap range as its overlapper.";
		k=k_;
		measurer=new FusionJoinDiagnostic(k, counts);
		row.clear().append("FusionPass\t").append(k).tab().append(maxOverlap).tab().append(graphCounts).nl();
		writers[3].print(row);
	}

	/** Drops the only borrowed count reference before its owner clears or replaces the table. */
	void endPhase(){measurer=null;}

	/** Records a frozen candidate without querying truth, modifying contigs or short-circuiting features. */
	void record(final Contig source, final boolean sr, final int st,
			final Contig dest, final boolean dr, final int dt, final int overlap){
		assert(measurer!=null && !closed) : "Feature extraction must be inside an open, live-count phase.";
		final int start=Math.max(0, source.length()-st-overlap-200);
		final int stop=Math.min(dest.length(), dt+overlap+200);
		final byte[] a=orientedSlice(source, sr, start, source.length());
		final byte[] b=orientedSlice(dest, dr, 0, stop);
		final FusionJoinDiagnostic.Join join=new FusionJoinDiagnostic.Join(++events, k,
				source.id, dest.id, overlap, st, dt, a, b);
		row.clear().append("OverlapPair\t").append(source.id).tab().append(dest.id);
		row.tab().append(source.length()).tab().append(dest.length()).tab().append(overlap);
		row.tab().append(st).tab().append(dt).tab().append(start).tab().append(a).tab().append(b).nl();
		writers[3].print(row);
		row.clear().append(events).tab().append(k).tab().append(k).tab().append(source.id).tab().append(dest.id);
		row.tab().append(overlap).tab().append(st).tab().append(dt);
		row.tab().append(join.left.length-overlap).tab().append(join.right.length-overlap);
		// Sentinel metadata is not truth; the post-assembly labeler supplies binary labels.
		row.append("\tpending\t-1\t-1\t-1\t-1");
		measurer.features(join, row, writers[1], writers[2]);
		writers[0].print(row.nl());
	}

	/** Copies only the bounded oriented context; cached contig metadata is untouched. */
	private static byte[] orientedSlice(final Contig contig, final boolean reverse, final int start, final int stop){
		assert(start>=0 && start<stop && stop<=contig.length()) : "Fusion context must remain inside the supplied contig.";
		final byte[] bases=new byte[stop-start];
		for(int p=start; p<stop; p++){
			bases[p-start]=reverse ? AminoAcid.baseToComplementExtended[contig.bases[contig.length()-1-p]] : contig.bases[p];
		}
		return bases;
	}

	/** Closes every started writer, including failure paths and legitimate zero-candidate assemblies. */
	@Override
	public void close(){
		if(closed){return;}
		closed=true;
		endPhase();
		boolean error=false;
		for(ByteStreamWriter writer : writers){if(writer!=null){error=writer.poisonAndWait() | error;}}
		if(error){throw new RuntimeException("Fusion vector sidecar write failed.");}
		System.err.println("FUSION_VECTOR_OUTPUT_CLOSED rows="+events);
	}

	/** Proves that rejected cycle proposals are collected without changing the graph. */
	static void selfTest(final String directory){
		final String prefix=directory+"/cycle";
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		for(String bases : new String[]{"GGATTAAC", "AACGGCCT", "CCTCCGGA"}){
			final Contig contig=new Contig(bases.getBytes(java.nio.charset.StandardCharsets.US_ASCII), contigs.size());
			contig.coverage=20;
			contig.leftCode=contig.rightCode=Tadpole.DEAD_END;
			contig.leftBridgeEndpoint=contig.rightBridgeEndpoint=true;
			contigs.add(contig);
		}
		try(FusionJoinCollector collector=new FusionJoinCollector(prefix,
				new ArrayList<String>(), directory+"/unused.fa", null)){
			collector.beginPhase(3, 3, false, new TadpoleGraph.Counts(){
				@Override
				public int count(final ukmer.Kmer key){return 1;}
			});
			final CrossKTipOverlapper overlapper=new CrossKTipOverlapper(contigs, 3, 3);
			overlapper.collector=collector;
			if(overlapper.addEdges()!=0){throw new AssertionError("Collector changed cycle rejection.");}
			for(Contig contig : contigs){
				if(contig.leftEdgeCount()!=0 || contig.rightEdgeCount()!=0){
					throw new AssertionError("Rejected cycle left graph edges.");
				}
			}
		}
		if(FusionJoinDiagnostic.readJoins(prefix+".trace.tsv", true).size()!=3){
			throw new AssertionError("Rejected cycle proposals were lost from the pre-guard census.");
		}
		System.err.println("FUSION_COLLECTOR_CYCLE_PASS rejected_proposals=3");
	}

	private final ByteStreamWriter[] writers=new ByteStreamWriter[4];
	private final ByteBuilder row=new ByteBuilder();
	private FusionJoinDiagnostic measurer;
	private int k, events;
	private boolean closed;
}
