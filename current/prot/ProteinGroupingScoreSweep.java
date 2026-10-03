package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import shared.Timer;

/**
 * Exhaustive protein-grouping scoring experiment requested for MAG-QC.
 *
 * <p>For every unordered sequence pair, this tool computes one local BLOSUM62
 * alignment under a caller-selected <em>linear</em> gap cost (no distinct opening
 * charge). The identical alignment is then rendered as two acceptance statistics:</p>
 * <ul>
 * <li>AAI: identical standard-residue columns / all alignment columns, including gaps.</li>
 * <li>BSR: symmetric aligned-span BLOSUM score ratio, multiplied by a gentle terminal
 * coverage penalty C^gamma, where C=(qSpan+tSpan)/(qLen+tLen).</li>
 * </ul>
 *
 * <p>Both statistics are crossed with score and bidirectional-coverage thresholds.
 * Each point forms deterministic single-linkage groups in memory; the compact output
 * is a group-count/edge-count surface rather than another multi-terabyte edge list.
 * This is an assay, not the production grouping implementation.</p>
 *
 * <p>Usage: {@code java prot.ProteinGroupingScoreSweep in=family.faa out=surface.tsv
 * gapmult=1.0 threads=64}. With the nominal BLOSUM scale fixed at +4 for a common
 * perfect match, gap multipliers 0/0.5/1/1.5 become linear residue costs 0/2/4/6.</p>
 *
 * @author Elly
 */
public final class ProteinGroupingScoreSweep {

	private static final double[] DEFAULT_AAI_THRESHOLDS={85,87,88,89,89.5,90,90.5,91,92,95};
	private static final double[] DEFAULT_BSR_THRESHOLDS={70,75,80,82,84,85,86,87,88,89,
		90,91,92,93,94,95,96,97,98,99};
	private static final double[] DEFAULT_COVERAGE_THRESHOLDS={60,70,75,80,85,90};
	private static final int NOMINAL_MATCH_SCORE=4;

	private ProteinGroupingScoreSweep(){}

	public static void main(final String[] args) throws Exception {
		String in=null, out=null;
		double gapMultiplier=1.0, gamma=0.25;
		int threads=Runtime.getRuntime().availableProcessors();
		boolean overwrite=true;
		double[] aaiThresholds=DEFAULT_AAI_THRESHOLDS;
		double[] bsrThresholds=DEFAULT_BSR_THRESHOLDS;
		double[] coverageThresholds=DEFAULT_COVERAGE_THRESHOLDS;

		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0,eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("in") || a.equals("input")){in=b;}
			else if(a.equals("out") || a.equals("output")){out=b;}
			else if(a.equals("gapmult") || a.equals("gapmultiplier")){gapMultiplier=Double.parseDouble(b);}
			else if(a.equals("gamma")){gamma=Double.parseDouble(b);}
			else if(a.equals("threads") || a.equals("t")){threads=Integer.parseInt(b);}
			else if(a.equals("aaithresholds") || a.equals("aai")){aaiThresholds=parseList(b);}
			else if(a.equals("bsrthresholds") || a.equals("bsr")){bsrThresholds=parseList(b);}
			else if(a.equals("covthresholds") || a.equals("cov")){coverageThresholds=parseList(b);}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(in==null || out==null){throw new RuntimeException("Required: in=<proteins.faa> out=<surface.tsv>");}
		if(gapMultiplier<0 || gamma<0 || threads<1){
			throw new IllegalArgumentException("Require gapmult>=0, gamma>=0, threads>=1.");
		}
		final double gapExact=NOMINAL_MATCH_SCORE*gapMultiplier;
		final int gapCost=(int)Math.round(gapExact);
		final double gammaFinal=gamma;
		if(Math.abs(gapExact-gapCost)>1e-9){
			throw new IllegalArgumentException("gapmult*"+NOMINAL_MATCH_SCORE+
				" must be an integer score in this assay: "+gapMultiplier+" -> "+gapExact);
		}

		final ArrayList<ProteinSequence> seqs=new ArrayList<ProteinSequence>(ProteinSearch.readFasta(in));
		Collections.sort(seqs, Comparator.comparing(s -> s.id));
		checkIds(seqs);
		final int n=seqs.size();
		final long pairCountLong=((long)n*(n-1))/2;
		if(pairCountLong>Integer.MAX_VALUE){
			throw new RuntimeException("This exhaustive assay uses int-indexed metric arrays; "+n+
				" sequences imply "+pairCountLong+" pairs. Stratify or sample this family.");
		}
		final int pairCount=(int)pairCountLong;
		System.err.println("Loaded "+n+" sorted unique sequences; unordered pairs="+pairCount+
			", gapmult="+gapMultiplier+", linear_gap_cost="+gapCost+", gamma="+gamma+".");

		final float[] aai=new float[pairCount];
		final float[] bsr=new float[pairCount];
		final float[] qcov=new float[pairCount];
		final float[] tcov=new float[pairCount];
		final AtomicInteger nextQuery=new AtomicInteger(0);
		final AtomicLong positiveAlignments=new AtomicLong(0);
		final ExecutorService pool=Executors.newFixedThreadPool(threads);
		final Timer alignTimer=new Timer();
		final ArrayList<Future<Void>> workers=new ArrayList<Future<Void>>();
		for(int tid=0; tid<threads; tid++){
			workers.add(pool.submit(new Callable<Void>(){
				@Override
				public Void call(){
					final Metrics m=new Metrics();
					long localPositive=0;
					for(int i=nextQuery.getAndIncrement(); i<n-1; i=nextQuery.getAndIncrement()){
						final ProteinSequence q=seqs.get(i);
						int idx=pairIndex(i, i+1, n);
						for(int j=i+1; j<n; j++, idx++){
							if(!scorePair(q.enc, seqs.get(j).enc, gapCost, gammaFinal, m)){continue;}
							aai[idx]=(float)(100*m.aai);
							bsr[idx]=(float)(100*m.adjustedScoreRatio);
							qcov[idx]=(float)(100*m.qCoverage);
							tcov[idx]=(float)(100*m.tCoverage);
							localPositive++;
						}
					}
					positiveAlignments.addAndGet(localPositive);
					return null;
				}
			}));
		}
		for(Future<Void> f : workers){f.get();}
		alignTimer.stop();
		System.err.println("Aligned all pairs in "+alignTimer+"; positive alignments="+
			positiveAlignments.get()+".");

		final ArrayList<Point> points=new ArrayList<Point>();
		addPoints(points, "AAI", aaiThresholds, coverageThresholds);
		addPoints(points, "BSR", bsrThresholds, coverageThresholds);
		final AtomicInteger nextPoint=new AtomicInteger(0);
		final ArrayList<Future<Void>> sweepWorkers=new ArrayList<Future<Void>>();
		final Timer sweepTimer=new Timer();
		for(int tid=0; tid<Math.min(threads, points.size()); tid++){
			sweepWorkers.add(pool.submit(new Callable<Void>(){
				@Override
				public Void call(){
					for(int pi=nextPoint.getAndIncrement(); pi<points.size(); pi=nextPoint.getAndIncrement()){
						final Point p=points.get(pi);
						p.evaluate(p.method.equals("AAI") ? aai : bsr, qcov, tcov, n);
					}
					return null;
				}
			}));
		}
		for(Future<Void> f : sweepWorkers){f.get();}
		pool.shutdown();
		sweepTimer.stop();
		Collections.sort(points);
		writeResults(out, overwrite, points, in, n, pairCount, positiveAlignments.get(),
			gapMultiplier, gapCost, gamma, threads, alignTimer.elapsed, sweepTimer.elapsed);
		System.err.println("Evaluated "+points.size()+" threshold points in "+sweepTimer+
			"; wrote "+out+".");
	}

	/** Metrics for one pair, expressed as fractions in [0,1]. */
	static final class Metrics {
		double aai, overlapScoreRatio, lengthCoverage, adjustedScoreRatio;
		double qCoverage, tCoverage;
		int rawScore;
		Metrics(){}
		Metrics(double aai_, double overlap_, double lengthCoverage_, double adjusted_,
				double qcov_, double tcov_, int rawScore_){
			set(aai_, overlap_, lengthCoverage_, adjusted_, qcov_, tcov_, rawScore_);
		}
		void set(double aai_, double overlap_, double lengthCoverage_, double adjusted_,
				double qcov_, double tcov_, int rawScore_){
			aai=aai_; overlapScoreRatio=overlap_; lengthCoverage=lengthCoverage_;
			adjustedScoreRatio=adjusted_; qCoverage=qcov_; tCoverage=tcov_; rawScore=rawScore_;
		}
	}

	/** Convenience allocating form used by focused tests. */
	static Metrics scorePair(final byte[] q, final byte[] t, final int gapCost, final double gamma){
		final Metrics out=new Metrics();
		return scorePair(q,t,gapCost,gamma,out) ? out : null;
	}

	/**
	 * Scores one pair under BLOSUM62 plus a linear per-gap-residue cost, writing into
	 * caller-owned metric scratch. AAAligner reuses its per-thread DP matrices, though
	 * each successful alignment still returns one AAAlignment result object.
	 */
	static boolean scorePair(final byte[] q, final byte[] t, final int gapCost,
			final double gamma, final Metrics out){
		final AAAlignment aln=AAAligner.alignLinear(q, t, gapCost, false);
		if(aln==null){return false;}
		final int qSpan=aln.qStop-aln.qStart+1;
		final int tSpan=aln.tStop-aln.tStart+1;
		final int qMax=selfScore(q, aln.qStart, aln.qStop);
		final int tMax=selfScore(t, aln.tStart, aln.tStop);
		final int denom=qMax+tMax;
		double overlap=(denom<=0 ? 0 : (2.0*aln.rawScore)/denom);
		// Standard BLOSUM62 rows have their maximum on the diagonal. X contributes zero
		// attainable self score below, so numerical noise/ambiguity must not create >1.
		if(overlap<0){overlap=0;}
		if(overlap>1){overlap=1;}
		final double qCov=qSpan/(double)q.length;
		final double tCov=tSpan/(double)t.length;
		final double lengthCoverage=(qSpan+tSpan)/(double)(q.length+t.length);
		final double adjusted=overlap*Math.pow(lengthCoverage, gamma);
		out.set(aln.length==0 ? 0 : aln.identities/(double)aln.length,
			overlap, lengthCoverage, adjusted, qCov, tCov, aln.rawScore);
		return true;
	}

	/** Sum of attainable diagonal self scores over an inclusive aligned span; X contributes zero. */
	static int selfScore(final byte[] seq, final int start, final int stop){
		int sum=0;
		for(int i=start; i<=stop; i++){
			final int x=Blosum62.score(seq[i], seq[i]);
			if(x>0){sum+=x;}
		}
		return sum;
	}

	static int pairIndex(final int i, final int j, final int n){
		if(i<0 || j<=i || j>=n){throw new IllegalArgumentException("Bad pair: "+i+","+j+" n="+n);}
		final long offset=((long)i*(2L*n-i-1))/2;
		return (int)(offset+j-i-1);
	}

	private static void addPoints(final ArrayList<Point> out, final String method,
			final double[] scores, final double[] coverages){
		for(double score : scores){
			for(double cov : coverages){out.add(new Point(method, score, cov));}
		}
	}

	static final class Point implements Comparable<Point> {
		final String method;
		final double threshold, coverage;
		long passed;
		int groups, largest, singletons;
		Point(String method_, double threshold_, double coverage_){
			method=method_; threshold=threshold_; coverage=coverage_;
		}
		void evaluate(final float[] score, final float[] qcov, final float[] tcov, final int n){
			final int[] parent=new int[n], rank=new int[n];
			for(int i=0; i<n; i++){parent[i]=i;}
			long pass=0;
			for(int i=0; i<n-1; i++){
				int idx=pairIndex(i, i+1, n);
				for(int j=i+1; j<n; j++, idx++){
					if(score[idx]>=threshold && qcov[idx]>=coverage && tcov[idx]>=coverage){
						IdentityGroupBuilder.union(parent, rank, i, j);
						pass++;
					}
				}
			}
			final int[] sizes=new int[n];
			for(int i=0; i<n; i++){sizes[IdentityGroupBuilder.find(parent,i)]++;}
			int g=0, l=0, s=0;
			for(int x : sizes){if(x>0){g++; if(x>l){l=x;} if(x==1){s++;}}}
			passed=pass; groups=g; largest=l; singletons=s;
		}
		@Override
		public int compareTo(Point b){
			int x=method.compareTo(b.method);
			if(x!=0){return x;}
			x=Double.compare(threshold,b.threshold);
			return x!=0 ? x : Double.compare(coverage,b.coverage);
		}
	}

	private static void writeResults(final String out, final boolean overwrite,
			final List<Point> points, final String in, final int n, final int pairCount,
			final long positive, final double gapMult, final int gapCost, final double gamma,
			final int threads, final long alignMillis, final long sweepMillis){
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.TXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("#method\tscore_threshold_pct\tqcov_threshold_pct\ttcov_threshold_pct"+
			"\tedges_passed\tgroups\tlargest_group\tsingletons");
		for(Point p : points){
			bsw.println(p.method+'\t'+Double.toString(p.threshold)+'\t'+Double.toString(p.coverage)+
				'\t'+Double.toString(p.coverage)+'\t'+p.passed+'\t'+p.groups+'\t'+p.largest+'\t'+p.singletons);
		}
		bsw.poisonAndWait();

		final FileFormat mf=FileFormat.testOutput(out+".meta", FileFormat.TXT, null, false, overwrite, false, false);
		final ByteStreamWriter meta=new ByteStreamWriter(mf);
		meta.start();
		meta.println("tool\tprot.ProteinGroupingScoreSweep");
		meta.println("input\t"+in);
		meta.println("sequences\t"+n);
		meta.println("unordered_pairs\t"+pairCount);
		meta.println("positive_alignments\t"+positive);
		meta.println("substitution_matrix\tBLOSUM62");
		meta.println("alignment\tlocal");
		meta.println("gap_model\tlinear_per_residue_no_opening_differential");
		meta.println("nominal_match_score\t"+NOMINAL_MATCH_SCORE);
		meta.println("gap_multiplier\t"+gapMult);
		meta.println("gap_cost_per_residue\t"+gapCost);
		meta.println("bsr_formula\t2*rawScore/(qAlignedSelfScore+tAlignedSelfScore)*lengthCoverage^gamma");
		meta.println("length_coverage_formula\t(qAlignedSpan+tAlignedSpan)/(qLength+tLength)");
		meta.println("gamma\t"+gamma);
		meta.println("aai_formula\tidentical_standard_columns/alignment_columns_including_gaps");
		meta.println("metric_storage\tfloat32");
		meta.println("threshold_comparison\tfloat_metric_greater_than_or_equal_to_double_threshold_boundary_approximate");
		meta.println("threads\t"+threads);
		meta.println("alignment_nanos\t"+alignMillis);
		meta.println("sweep_nanos\t"+sweepMillis);
		meta.poisonAndWait();
	}

	private static double[] parseList(final String s){
		final String[] split=s.split(",");
		final double[] out=new double[split.length];
		for(int i=0; i<out.length; i++){
			out[i]=Double.parseDouble(split[i]);
			if(out[i]<0 || out[i]>100){throw new IllegalArgumentException("Threshold outside [0,100]: "+out[i]);}
		}
		return out;
	}

	private static void checkIds(final List<ProteinSequence> seqs){
		final HashSet<String> seen=new HashSet<String>();
		for(ProteinSequence s : seqs){
			if(!seen.add(s.id)){throw new RuntimeException("Duplicate sequence ID: "+s.id);}
		}
	}
}
