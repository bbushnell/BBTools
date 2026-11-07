package prok;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaNativeTrainer.Method;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Detection and endpoint comparisons on one complete, immutable reference axis.
 * Uses the public half-shorter-span one-to-one matcher. Candidate matching gives
 * baseline-detected references priority, avoiding preventable apparent losses.
 * Endpoint offsets follow transcription: negative is upstream at either end.
 * Results borrow their input lists; callers must not mutate lists or loci afterward.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryComparison{

	/** Grades calls once. A nonnull baseline must use this exact reference-list object. */
	public static Result grade(final ArrayList<Locus> truth, final ArrayList<Locus> calls, final Result baseline){
		if(truth==null || truth.isEmpty() || calls==null || baseline!=null && baseline.truth!=truth){
			throw new IllegalArgumentException("Boundary grading requires calls and one shared nonempty reference axis");
		}
		final boolean[] protectedReferences=baseline==null ? null : new boolean[truth.size()];
		if(baseline!=null){
			for(int i=0; i<protectedReferences.length; i++){protectedReferences[i]=baseline.matches[i]>=0;}
		}
		return new Result(truth, calls, NcrnaLocusComparison.match(truth, calls, protectedReferences));
	}

	/** Returns null for OFF when neither candidate meets the declared eligibility rule. */
	public static Method choose(final Result off, final Result classifier, final Result regression,
			final int maxLost, final double minPrecision, final Method tiePreference){
		validateLimits(maxLost, minPrecision);
		if(tiePreference==null){throw new IllegalArgumentException("Specify a classifier or regression tie preference");}
		final Change c=new Change(off, classifier);
		final Change r=new Change(off, regression);
		final boolean classifierPass=c.eligible(maxLost, minPrecision);
		final boolean regressionPass=r.eligible(maxLost, minPrecision);
		if(!classifierPass){return regressionPass ? Method.REGRESSION : null;}
		if(!regressionPass){return Method.CLASSIFIER;}
		if(c.exactGain!=r.exactGain){return c.exactGain>r.exactGain ? Method.CLASSIFIER : Method.REGRESSION;}
		if(c.lost!=r.lost){return c.lost<r.lost ? Method.CLASSIFIER : Method.REGRESSION;}
		if(classifier.matched!=regression.matched){return classifier.matched>regression.matched ? Method.CLASSIFIER : Method.REGRESSION;}
		if(classifier.unmatched!=regression.unmatched){return classifier.unmatched<regression.unmatched ? Method.CLASSIFIER : Method.REGRESSION;}
		return tiePreference;
	}

	/** Writes all detection and boundary counts, with rates explicitly NA when undefined. */
	public static void writeSummary(final String path, final Result off, final Result classifier, final Result regression){
		try(Output output=new Output(path)){
			output.writer.print("method\ttp\tfp\tfn\tprecision\trecall\texact_5prime\texact_3prime\tboth_exact\tlost\tgained\tboth_exact_gain\n");
			final ByteBuilder row=new ByteBuilder();
			writeRow(output.writer, row, "off", off, off);
			writeRow(output.writer, row, "classifier", off, classifier);
			writeRow(output.writer, row, "regression", off, regression);
		}
	}

	/** Writes two saved-call results with the overlap rule and unadjudicated-count limitation explicit. */
	public static void writeSavedSummary(final String path, final Result before, final Result after){
		requireSameTruth(before, after);
		try(Output output=new Output(path)){
			output.writer.print("# overlap_rule=HALF_SHORTER; one-to-one, same contig/strand, overlap>=half the shorter span; baseline references prioritized\n");
			output.writer.print("# Saved calls only; fp means unmatched calls, not independently adjudicated false positives\n");
			output.writer.print("method\ttp\tfp\tfn\tprecision\trecall\texact_5prime\texact_3prime\tboth_exact\tlost\tgained\tboth_exact_gain\n");
			final ByteBuilder row=new ByteBuilder();
			writeRow(output.writer, row, "baseline", before, before);
			writeRow(output.writer, row, "candidate", before, after);
		}
	}

	/** Writes one row per reference, including missed, gained and lost loci and signed end shifts. */
	public static void writeTransitions(final String path, final Result before, final Result after){
		requireSameTruth(before, after);
		try(Output output=new Output(path)){
			output.writer.print("reference_index\tseqid\ttruth_start\ttruth_stop\tstrand\tbefore_start\tbefore_stop\tafter_start\tafter_stop\tbefore_5prime_offset\tbefore_3prime_offset\tafter_5prime_offset\tafter_3prime_offset\tdetection\tboth_exact\n");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<before.truth.size(); i++){
				final Locus reference=before.truth.get(i);
				final Locus a=before.call(i);
				final Locus b=after.call(i);
				row.clear();
				row.append(i).tab().append(reference.seqid).tab().append(reference.start).tab().append(reference.stop).tab();
				row.append(reference.strand==0 ? '+' : '-').tab();
				coordinates(row, a);
				row.tab();
				coordinates(row, b);
				row.tab();
				offsets(row, a, reference);
				row.tab();
				offsets(row, b, reference);
				row.tab().append(a==null ? (b==null ? "missed" : "gained") : (b==null ? "lost" : "retained")).tab();
				final boolean aExact=exact(a, reference);
				final boolean bExact=exact(b, reference);
				row.append(aExact ? (bExact ? "retained" : "lost") : (bExact ? "gained" : "neither")).nl();
				output.writer.print(row);
			}
		}
	}

	private static void writeRow(final ByteStreamWriter output, final ByteBuilder row, final String method,
			final Result off, final Result result){
		final Change change=new Change(off, result);
		row.clear();
		row.append(method).tab().append(result.matched).tab().append(result.unmatched).tab().append(result.missed).tab();
		if(result.calls.isEmpty()){row.append("NA");}else{row.append(result.precision(), 8, true);}
		row.tab().append(result.recall(), 8, true).tab().append(result.exactStart).tab().append(result.exactStop).tab();
		row.append(result.exactBoth).tab().append(change.lost).tab().append(change.gained).tab().append(change.exactGain).nl();
		output.print(row);
	}

	private static void coordinates(final ByteBuilder row, final Locus call){
		if(call==null){row.append("NA\tNA");}else{row.append(call.start).tab().append(call.stop);}
	}

	private static void offsets(final ByteBuilder row, final Locus call, final Locus reference){
		if(call==null){row.append("NA\tNA");}else{row.append(offset(call, reference, true)).tab().append(offset(call, reference, false));}
	}

	/** Positive coordinates bound these differences within signed int, even at the maximum coordinate. */
	public static int offset(final Locus call, final Locus reference, final boolean fivePrime){
		if(call==null || reference==null || call.strand!=reference.strand){
			throw new IllegalArgumentException("Endpoint offsets require a matched same-strand call and reference");
		}
		if(reference.strand==0){return fivePrime ? call.start-reference.start : call.stop-reference.stop;}
		return fivePrime ? reference.stop-call.stop : reference.start-call.start;
	}

	private static boolean exact(final Locus call, final Locus reference){
		return call!=null && call.start==reference.start && call.stop==reference.stop;
	}

	private static void requireSameTruth(final Result before, final Result after){
		if(before==null || after==null || before.truth!=after.truth){
			throw new IllegalArgumentException("Boundary comparisons require results on the same immutable reference axis");
		}
	}

	private static void validateLimits(final int maxLost, final double minPrecision){
		if(maxLost<0 || !Double.isFinite(minPrecision) || minPrecision<0 || minPrecision>1){
			throw new IllegalArgumentException("Boundary selection requires maxlost>=0 and minprecision in [0,1]");
		}
	}

	/** A retained matching with primitive counts; input lists and matching array are read-only. */
	public static final class Result{
		private Result(final ArrayList<Locus> truth_, final ArrayList<Locus> calls_, final int[] matches_){
			truth=truth_;
			calls=calls_;
			matches=matches_;
			int tp=0, start=0, stop=0, both=0;
			for(int i=0; i<matches.length; i++){
				final Locus call=call(i);
				if(call==null){continue;}
				final Locus reference=truth.get(i);
				tp++;
				final boolean a=offset(call, reference, true)==0;
				final boolean b=offset(call, reference, false)==0;
				if(a){start++;}
				if(b){stop++;}
				if(a && b){both++;}
			}
			matched=tp;
			unmatched=calls.size()-tp;
			missed=truth.size()-tp;
			exactStart=start;
			exactStop=stop;
			exactBoth=both;
			assert(unmatched>=0 && missed>=0 && both<=start && both<=stop && start<=tp && stop<=tp)
				: "One-to-one detection and exact-end subsets must conserve reference and call counts";
		}
		/** Matched call for one reference, or null when that reference was missed. */
		public Locus call(final int reference){return matches[reference]<0 ? null : calls.get(matches[reference]);}
		public double precision(){return calls.isEmpty() ? Double.NaN : matched/(double)calls.size();}
		public double recall(){return matched/(double)truth.size();}
		public final ArrayList<Locus> truth, calls;
		public final int[] matches;
		public final int matched, unmatched, missed, exactStart, exactStop, exactBoth;
	}

	/** Reference-level detection changes and the raw change in both-ends-exact count. */
	public static final class Change{
		public Change(final Result before, final Result after_){
			requireSameTruth(before, after_);
			after=after_;
			int losses=0, gains=0;
			for(int i=0; i<before.matches.length; i++){
				if(before.matches[i]>=0 && after.matches[i]<0){losses++;}
				if(before.matches[i]<0 && after.matches[i]>=0){gains++;}
			}
			lost=losses;
			gained=gains;
			exactGain=after.exactBoth-before.exactBoth;
			assert(after.matched-before.matched==gained-lost) : "Reference gains minus losses must equal the TP change";
		}
		/** The same eligibility rule applies to tuning candidates and the all-data refit. */
		public boolean eligible(final int maxLost, final double minPrecision){
			validateLimits(maxLost, minPrecision);
			return lost<=maxLost && after.precision()>=minPrecision && exactGain>0;
		}
		public final Result after;
		public final int lost, gained, exactGain;
	}

	private NcrnaBoundaryComparison(){}
}
