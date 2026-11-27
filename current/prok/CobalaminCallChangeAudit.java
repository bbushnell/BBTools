package prok;

import java.util.ArrayList;

import prok.CompareCobalaminBenchmarkGrades.Genome;
import prok.CompareCobalaminBenchmarkGrades.Query;
import prok.CompareCobalaminBenchmarkGrades.Reference;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;

/** Explains saved call geometry without rematching or rerunning a caller.
 * The existing template-based normalized reader owns parsing and validation.
 * Each candidate unmatched call has one disposition. Every same-input,
 * same-contig, same-strand positive overlap is retained as an edge, including
 * many-to-many edges. Overlap is evidence of proximity, never candidate ancestry.
 * Lost scored references form a separate population. No historical CM score,
 * contact flag or unavailable-anchor cause can be recovered from these tables.
 * @author Raiden
 */
final class CobalaminCallChangeAudit{

	CobalaminCallChangeAudit(final String out_){out=out_;}

	static void checkOutputs(final String out){
		for(final String suffix:SUFFIXES){
			if(!Tools.testOutputFiles(false, false, false, out+suffix)){
				throw new IllegalArgumentException("Saved call audit requires fresh output: "+out+suffix);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	void write(final ArrayList<Genome> genomes){
		assert(genomes!=null && !genomes.isEmpty()) : "Normalized reader validates the full input population before an audit";
		try(Output calls=new Output(out+SUFFIXES[0]); Output links=new Output(out+SUFFIXES[1]);
			Output losses=new Output(out+SUFFIXES[2]); Output lossLinks=new Output(out+SUFFIXES[3])){
			calls.writer.println("input_id\tafter_call\tseqid\tstrand\tstart\tstop\tmember\tmodel\toutcome\tdisposition\told_links\told_lost_tp\told_retained_tp\told_masked\told_unmatched\tmutual_unique_overlap\tgeometry\told_reference_overlap_bp\tqualification");
			links.writer.println("input_id\tafter_call\tbefore_call\tseqid\tstrand\tbefore_start\tbefore_stop\tafter_start\tafter_stop\tbefore_member\tafter_member\tbefore_model\tafter_model\tbefore_outcome\tbefore_reference\tbefore_reference_retained\toverlap_bp\tbefore_bp\tafter_bp\tlength_delta_bp\td5_bp\td3_bp\tgeometry\told_reference_overlap_bp\told_call_after_degree\tqualification");
			losses.writer.println("input_id\treference_id\tseqid\tstrand\tref_start\tref_stop\tbefore_start\tbefore_stop\tbefore_member\tbefore_model\tafter_links\tafter_unmatched\tafter_masked\tafter_tp_elsewhere\tgeometry\tqualification");
			lossLinks.writer.println("input_id\treference_id\tafter_call\tseqid\tstrand\tbefore_start\tbefore_stop\tafter_start\tafter_stop\tbefore_member\tafter_member\tbefore_model\tafter_model\tafter_outcome\tafter_reference\told_call_overlap_bp\told_reference_overlap_bp\tbefore_bp\tafter_bp\tlength_delta_bp\td5_bp\td3_bp\tgeometry\tqualification");
			for(final Genome genome:genomes){
				for(int side=0; side<2; side++){
					for(int type=0; type<4; type++){totals[side][type]+=genome.actual[side][type];}
				}
				for(int i=0; i<genome.after.size(); i++){
					final Query after=genome.after.get(i);
					if(unmatched(after)){unmatched(genome, i, calls, links);}
				}
				for(final Reference ref:genome.references){
					if(ref.scored && ref.winners[0]!=null && ref.winners[1]==null){lost(ref, losses, lossLinks);}
				}
			}
		}
		long n=0;
		for(final long[] row:callCounts){for(final long count:row){n+=count;}}
		if(n!=totals[1][3]){throw new IllegalStateException("Every candidate unmatched call must have exactly one audit disposition");}
		try(Output summary=new Output(out+SUFFIXES[4])){
			summary.writer.println("population\tdisposition\tgeometry\tcount\tqualification");
			for(int side=0; side<2; side++){
				for(int type=0; type<4; type++){count(summary, side==0 ? "before" : "after", TOTALS[type], "ALL", totals[side][type]);}
			}
			for(int state=0; state<STATES.length; state++){
				for(int geometry=0; geometry<GEOMETRIES.length; geometry++){
					if(callCounts[state][geometry]>0){count(summary, "candidate_unmatched", STATES[state], GEOMETRIES[geometry], callCounts[state][geometry]);}
				}
			}
			for(int geometry=0; geometry<GEOMETRIES.length; geometry++){
				if(lossCounts[geometry]>0){count(summary, "lost_reference", "LOST_SCORED_MATCH", GEOMETRIES[geometry], lossCounts[geometry]);}
			}
		}
	}

	private void unmatched(final Genome genome, final int afterIndex, final Output calls, final Output links){
		final Query after=genome.after.get(afterIndex);
		assert(unmatched(after)) : "The primary audit population is exactly the saved unmatched queries";
		int n=0, mask=0, only=-1, referenceOverlap=0;
		final int[] outcomes=new int[4];
		for(int i=0; i<genome.before.size(); i++){
			final Query before=genome.before.get(i);
			final int bp=overlap(before.seqid, before.strand, before.start, before.stop, after);
			if(bp==0){continue;}
			final Reference ref=assigned(genome, before);
			final int state=ref!=null ? ref.winners[1]==null ? 0 : 1 : before.outcome.equals("MASKED") ? 2 : 3;
			outcomes[state]++;
			mask|=1<<state;
			n++;
			only=i;
			final int oldRefOverlap=ref==null ? 0 : overlap(ref.seqid, ref.strand, ref.start, ref.stop, after);
			referenceOverlap=Math.max(referenceOverlap, oldRefOverlap);
			row.clear().append(genome.id).tab().append(afterIndex).tab().append(i).tab().append(after.seqid).tab().append(after.strand)
				.tab().append(before.start).tab().append(before.stop).tab().append(after.start).tab().append(after.stop)
				.tab().append(before.member).tab().append(after.member).tab().append(before.model).tab().append(after.model)
				.tab().append(before.outcome).tab().append(before.reference).tab().append(ref==null ? "NA" : ref.winners[1]!=null ? "true" : "false")
				.tab().append(bp);
			geometryFields(before, after);
			links.writer.print(row.tab().append(oldRefOverlap).tab().append(degree(genome, before)).tab().append(QUALIFICATION).nl());
		}
		final int state=mask==0 ? 5 : Integer.bitCount(mask)>1 ? 4 : Integer.numberOfTrailingZeros(mask);
		final Query before=n==1 ? genome.before.get(only) : null;
		final int shape=n==0 ? 5 : n>1 ? 4 : geometry(before, after);
		callCounts[state][shape]++;
		row.clear().append(genome.id).tab().append(afterIndex).tab().append(after.seqid).tab().append(after.strand)
			.tab().append(after.start).tab().append(after.stop).tab().append(after.member).tab().append(after.model).tab().append(after.outcome)
			.tab().append(STATES[state]).tab().append(n);
		for(final int count:outcomes){row.tab().append(count);}
		calls.writer.print(row.tab().append(n==1 && degree(genome, before)==1 ? "true" : "false")
			.tab().append(GEOMETRIES[shape]).tab().append(referenceOverlap).tab().append(QUALIFICATION).nl());
	}

	private void lost(final Reference ref, final Output losses, final Output links){
		final Query before=ref.winners[0];
		assert(before!=null && ref.winners[1]==null) : "Lost references are validated old winners without a candidate winner, independently of unmatched-call links";
		int n=0, unmatched=0, masked=0, tp=0, shape=5;
		for(int i=0; i<ref.genome.after.size(); i++){
			final Query after=ref.genome.after.get(i);
			final int callBp=overlap(before.seqid, before.strand, before.start, before.stop, after);
			final int refBp=overlap(ref.seqid, ref.strand, ref.start, ref.stop, after);
			if(callBp==0 && refBp==0){continue;}
			shape=n==0 ? geometry(before, after) : 4;
			n++;
			if(unmatched(after)){unmatched++;}
			else if(after.outcome.equals("MASKED")){masked++;}
			else{tp++;}
			row.clear().append(ref.genome.id).tab().append(ref.id).tab().append(i).tab().append(ref.seqid).tab().append(ref.strand)
				.tab().append(before.start).tab().append(before.stop).tab().append(after.start).tab().append(after.stop)
				.tab().append(before.member).tab().append(after.member).tab().append(before.model).tab().append(after.model)
				.tab().append(after.outcome).tab().append(after.reference).tab().append(callBp).tab().append(refBp);
			geometryFields(before, after);
			links.writer.print(row.tab().append(QUALIFICATION).nl());
		}
		assert(n==unmatched+masked+tp) : "Each linked candidate retains one saved outcome; lost-reference links are never extra unmatched calls";
		lossCounts[shape]++;
		losses.writer.print(row.clear().append(ref.genome.id).tab().append(ref.id).tab().append(ref.seqid).tab().append(ref.strand)
			.tab().append(ref.start).tab().append(ref.stop).tab().append(before.start).tab().append(before.stop)
			.tab().append(before.member).tab().append(before.model).tab().append(n).tab().append(unmatched).tab().append(masked).tab().append(tp)
			.tab().append(GEOMETRIES[shape]).tab().append(QUALIFICATION).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Geometry          ----------------*/
	/*--------------------------------------------------------------*/

	private void geometryFields(final Query before, final Query after){
		assert(before.strand==after.strand && before.seqid.equals(after.seqid)) : "Signed endpoint changes compare only the same contig and strand";
		final int beforeLength=before.stop-before.start+1, afterLength=after.stop-after.start+1;
		final long d5=after.strand=='-' ? (long)before.stop-after.stop : (long)after.start-before.start;
		final long d3=after.strand=='-' ? (long)before.start-after.start : (long)after.stop-before.stop;
		row.tab().append(beforeLength).tab().append(afterLength).tab().append((long)afterLength-beforeLength)
			.tab().append(d5).tab().append(d3).tab().append(GEOMETRIES[geometry(before, after)]);
	}

	private static int geometry(final Query before, final Query after){
		assert(before!=null && after!=null) : "A geometric change needs two saved spans, never a missing-call sentinel";
		final int a=before.stop-before.start, b=after.stop-after.start;
		return b<a ? 0 : b>a ? 1 : before.start!=after.start ? 2 : 3;
	}
	private static int overlap(final String seqid, final char strand, final int start, final int stop, final Query query){
		assert(start>0 && stop>=start) : "Normalized positive inclusive intervals prevent negative lengths and overflow";
		if(strand!=query.strand || !seqid.equals(query.seqid)){return 0;}
		return Math.max(0, Math.min(stop, query.stop)-Math.max(start, query.start)+1);
	}
	private static int degree(final Genome genome, final Query before){
		assert(before!=null) : "Reverse overlap degree is defined only for an existing old call";
		int n=0;
		for(final Query after:genome.after){if(overlap(before.seqid, before.strand, before.start, before.stop, after)>0){n++;}}
		return n;
	}
	private static Reference assigned(final Genome genome, final Query query){
		if(!query.outcome.equals("TP")){return null;}
		for(final Reference ref:genome.references){if(ref.winners[0]==query){return ref;}}
		throw new IllegalStateException("A saved TP must retain its validated original reference binding");
	}
	private static boolean unmatched(final Query query){return query.outcome.startsWith("UNMATCHED_");}
	private void count(final Output output, final String population, final String state, final String geometry, final long count){
		assert(count>=0) : "Audit counts are conserved populations, never signed differences";
		output.writer.print(row.clear().append(population).tab().append(state).tab().append(geometry).tab().append(count).tab().append(QUALIFICATION).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------             State            ----------------*/
	/*--------------------------------------------------------------*/

	private final String out;
	private final ByteBuilder row=new ByteBuilder();
	private final long[][] totals=new long[2][4], callCounts=new long[6][6];
	private final long[] lossCounts=new long[6];
	private static final String[] SUFFIXES={".unmatched.tsv", ".unmatched_links.tsv", ".lost.tsv", ".lost_links.tsv", ".change_summary.tsv"};
	private static final String[] TOTALS={"CALLS", "MATCHED", "MASKED", "UNMATCHED"};
	private static final String[] STATES={"ONLY_OLD_LOST_TP", "ONLY_OLD_RETAINED_TP", "ONLY_OLD_MASKED", "ONLY_OLD_UNMATCHED", "MIXED_OLD_OUTCOMES", "NO_OLD_CALL_OVERLAP"};
	private static final String[] GEOMETRIES={"SHRUNK", "GROWN", "SHIFTED_EQUAL_LENGTH", "UNCHANGED_SPAN", "MULTIPLE_LINKS", "NO_LINK"};
	private static final String QUALIFICATION="SAVED_SPATIAL_LINKS_NOT_ANCESTRY_CM_SCORE_AND_FLAGS_UNAVAILABLE";
}
