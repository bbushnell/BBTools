package prot.cluster;

import java.io.File;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import prot.Blosum62;
import prot.ConsensusScalars;
import structures.ByteBuilder;

/** Decision-time marginal diagnostics, never an admission gate. Each worker owns
 * one collector; merge only after all queries finish. Family inputs must describe
 * the candidate at decision time, before any online update. No pair payloads are
 * retained. Float quantiles are intervals, and wrong-side counts have bounds where
 * an inclusive comparison cuts through a bin. They cannot predict alternate winners.
 * @author Sayu
 */
public final class AssignmentHistograms {

	public AssignmentHistograms(int passes,int width,int shortlist){
		if(passes<0 || passes==Integer.MAX_VALUE || width<1 || shortlist<1){throw new IllegalArgumentException("Positive sketch width/shortlist and nonnegative passes required");}
		this.passes=passes; this.width=width; this.shortlist=shortlist;
		last=new int[]{width,1000,900,2000,1000,1000,2000,shortlist};
		bins=new long[passes+1][COHORTS.length][METRICS][];
		missing=new long[passes+1][COHORTS.length][METRICS];
		counts=new long[passes+1][COUNT_NAMES.length];
		joint=new long[passes+1][shortlist+1][shortlist+1];
		alignedCounts=new long[passes+1][shortlist+1];
		for(int p=0;p<=passes;p++){for(int c=0;c<COHORTS.length;c++){for(int m=0;m<METRICS;m++){bins[p][c][m]=new long[Math.addExact(last[m],1)];}}}
	}

	/** Compute full-query scalar values once, using the native ASCII scalar API. */
	public void beginQuery(String stage,int pass,byte[] encodedTail){
		if(active || encodedTail==null || encodedTail.length==0){throw new IllegalStateException("Begin one nonempty query only after finishing its predecessor");}
		phase=phase(stage,pass);
		if(raw.length<encodedTail.length){raw=new byte[encodedTail.length];}
		for(int i=0;i<encodedTail.length;i++){
			final int aa=encodedTail[i];
			if(aa==Blosum62.X_CODE){raw[i]='X';}
			else if(aa>=0 && aa<20){raw[i]=AminoAcid.numberToAcid[aa];}
			else{throw new IllegalArgumentException("Invalid encoded query residue at "+i+": "+aa);}
		}
		queryMask=ConsensusScalars.computeScalarMask(raw,encodedTail.length,queryScalars,0);
		queryRejected=0; queryRank=0; active=true;
	}

	public void rejected(int hits,double cosine,double[] familyScalars,int mask){
		requireActive(); record(REJECTED,hits,cosine,familyScalars,mask);
		queryRejected=Math.addExact(queryRejected,1);
	}

	/** referenceSize is the candidate's size BEFORE admission. Frozen reassignment
	 * must pass the previous-pass reference size, never cleared/incoming membership. */
	public void chosen(int hits,double cosine,double[] familyScalars,int mask,int winningRank,int referenceSize){
		requireActive();
		if(referenceSize<1){throw new IllegalArgumentException("Positive pre-admission reference size required");}
		if(queryRank!=0 || winningRank<1 || winningRank>shortlist){throw new IllegalArgumentException("One chosen candidate with rank in 1..shortlist required");}
		record(CHOSEN,hits,cosine,familyScalars,mask);
		increment(bins[phase][CHOSEN][RANK],winningRank,1); queryRank=winningRank;
		if(referenceSize>=5){recordSizeCohort(SIZE5,hits,cosine); increment(counts[phase],CHOSEN_SIZE5,1);}
		if(referenceSize>=10){recordSizeCohort(SIZE10,hits,cosine); increment(counts[phase],CHOSEN_SIZE10,1);}
	}

	/** No-filter compatibility: every ranked candidate was aligned. */
	public void finishQuery(int candidates,int passingCandidates,boolean hasWinner){
		finishQuery(candidates,passingCandidates,hasWinner,0,0);
	}
	/** Prefiltered candidates never enter rejected(). Winning rank refers to the
	 * original ranked shortlist, so it may exceed the number actually aligned.
	 * Filter counts are disjoint: test hits first, then hydro for hit survivors. */
	public void finishQuery(int candidates,int passingCandidates,boolean hasWinner,int hitFiltered,int hydroFiltered){
		requireActive();
		if(candidates<0 || candidates>shortlist || hitFiltered<0 || hitFiltered>candidates || hydroFiltered<0 || hydroFiltered>candidates-hitFiltered){
			throw new IllegalStateException("Candidate conservation failed: ranked="+candidates+" hitFiltered="+hitFiltered+" hydroFiltered="+hydroFiltered);
		}
		final int aligned=candidates-hitFiltered-hydroFiltered;
		if(passingCandidates<0 || passingCandidates>aligned ||
			queryRejected!=aligned-passingCandidates || hasWinner!=(queryRank>0) ||
			hasWinner!=(passingCandidates>0) || queryRank>candidates){
			throw new IllegalStateException("Candidate conservation failed: ranked="+candidates+" aligned="+aligned+" passing="+passingCandidates+" rejected="+queryRejected+" rank="+queryRank+" winner="+hasWinner);
		}
		increment(counts[phase],QUERIES,1); increment(counts[phase],CHOSEN_COUNT,hasWinner ? 1 : 0);
		increment(counts[phase],REJECTED_COUNT,queryRejected);
		increment(counts[phase],NONWINNERS,passingCandidates-(hasWinner ? 1 : 0));
		increment(counts[phase],ALIGNED,aligned);
		increment(counts[phase],RANKED,candidates); increment(counts[phase],FILTERED,candidates-aligned);
		increment(counts[phase],HIT_FILTERED,hitFiltered); increment(counts[phase],HYDRO_FILTERED,hydroFiltered);
		increment(alignedCounts[phase],aligned,1);
		increment(joint[phase][candidates],queryRank,1); active=false;
	}

	public void merge(AssignmentHistograms other){
		if(other==null || other==this || active || other.active || passes!=other.passes || width!=other.width || shortlist!=other.shortlist){
			throw new IllegalArgumentException("Merge distinct idle collectors with identical schemas");
		}
		for(int p=0;p<=passes;p++){
			for(int k=0;k<counts[p].length;k++){increment(counts[p],k,other.counts[p][k]);}
			for(int n=0;n<=shortlist;n++){for(int r=0;r<=n;r++){increment(joint[p][n],r,other.joint[p][n][r]);}}
			for(int n=0;n<=shortlist;n++){increment(alignedCounts[p],n,other.alignedCounts[p][n]);}
			for(int c=0;c<COHORTS.length;c++){for(int m=0;m<METRICS;m++){
				increment(missing[p][c],m,other.missing[p][c][m]);
				for(int b=0;b<=last[m];b++){increment(bins[p][c][m],b,other.bins[p][c][m][b]);}
			}}
		}
	}

	/** Prevent an externally supplied collector from silently using another schema. */
	public void validateSchema(int expectedPasses,int expectedWidth,int expectedShortlist){
		if(passes!=expectedPasses || width!=expectedWidth || shortlist!=expectedShortlist){
			throw new IllegalArgumentException("Histogram schema mismatch: actual passes/width/shortlist="+passes+"/"+width+"/"+shortlist+
				" expected="+expectedPasses+"/"+expectedWidth+"/"+expectedShortlist);
		}
	}

	/** Four deterministic sidecars, with no timing or worker-count fields. */
	public static String[] outputPaths(String prefix){
		if(prefix==null || prefix.isEmpty()){throw new IllegalArgumentException("Histogram output prefix required");}
		final String[] paths=new String[SUFFIXES.length];
		for(int i=0;i<paths.length;i++){paths[i]=prefix+SUFFIXES[i];}
		return paths;
	}

	public void write(String prefix,boolean overwrite){
		if(active || prefix==null || prefix.isEmpty()){throw new IllegalStateException("Writing requires a prefix and no unfinished query");}
		validate();
		for(String suffix:SUFFIXES){if(!overwrite && new File(prefix+suffix).exists()){throw new IllegalArgumentException("Output exists: "+prefix+suffix);}}
		for(int file=0;file<SUFFIXES.length;file++){
			final ByteStreamWriter out=new ByteStreamWriter(prefix+SUFFIXES[file],overwrite,false,true);
			out.start();
			try{
				out.println(HEADERS[file]);
				if(file==0){writeBins(out);}else if(file==1){writeQuantiles(out);}else if(file==2){writeCounts(out);}else{writeWrongSide(out);}
			}finally{if(out.poisonAndWait()){throw new IllegalStateException("Assignment histogram output failed: "+prefix+SUFFIXES[file]);}}
		}
	}

	private void record(int cohort,int hits,double cosine,double[] familyScalars,int mask){
		if(hits<0 || hits>width || familyScalars==null || familyScalars.length!=5){throw new IllegalArgumentException("Pair requires hits in 0..width and five family scalars");}
		final int cb=bin(COSINE,cosine);
		// A malformed/unset validity mask is unavailable; nonfinite numbers still fail.
		final int defined=(mask>=0 && mask<=ConsensusScalars.ALL_DEFINED) ? mask&queryMask : 0;
		for(int s=0;s<5;s++){
			if(!Double.isFinite(familyScalars[s])){throw new IllegalArgumentException("Nonfinite family scalar "+s+": "+familyScalars[s]);}
			scalarBins[s]=(defined&(1<<s))==0 ? -1 : bin(s+2,Math.abs(queryScalars[s]-familyScalars[s]));
		}
		increment(bins[phase][cohort][HITS],hits,1); increment(bins[phase][cohort][COSINE],cb,1);
		for(int s=0;s<5;s++){
			if(scalarBins[s]<0){increment(missing[phase][cohort],s+2,1);}
			else{increment(bins[phase][cohort][s+2],scalarBins[s],1);}
		}
	}
	/** Reuse the full pair validation and scalar bins just computed by record(CHOSEN). */
	private void recordSizeCohort(int cohort,int hits,double cosine){
		assert(cohort==SIZE5 || cohort==SIZE10) : "Conditional cohorts contain only chosen hits, cosine and hydrophobicity";
		increment(bins[phase][cohort][HITS],hits,1); increment(bins[phase][cohort][COSINE],bin(COSINE,cosine),1);
		if(scalarBins[0]<0){increment(missing[phase][cohort],2,1);}else{increment(bins[phase][cohort][2],scalarBins[0],1);}
	}
	private static boolean hasMetric(int cohort,int metric){return cohort>=SIZE5 ? metric<=2 : !(cohort==REJECTED && metric==RANK);}
	private static int countIndex(int cohort){return cohort==CHOSEN ? CHOSEN_COUNT : cohort==REJECTED ? REJECTED_COUNT : cohort==SIZE5 ? CHOSEN_SIZE5 : CHOSEN_SIZE10;}

	private int bin(int metric,double value){
		final double max=edge(metric,last[metric]);
		if(!Double.isFinite(value) || value<0 || value>max){throw new IllegalArgumentException("Invalid "+NAMES[metric]+" value "+value+"; range [0,"+max+"]");}
		int b=Math.min(last[metric],(int)(value*SCALES[metric]));
		// Correct multiplication rounding against the actual declared double edges.
		if(b>0 && value<edge(metric,b)){b--;}
		if(b<last[metric] && value>=edge(metric,b+1)){b++;}
		return b;
	}
	private static double edge(int metric,int b){return b/(double)SCALES[metric];}
	private static boolean exact(int metric){return metric==HITS || metric==RANK;}
	private int phase(String stage,int pass){
		if("greedy".equals(stage) && pass==0){return 0;}
		if("reassign".equals(stage) && pass>=1 && pass<=passes){return pass;}
		throw new IllegalArgumentException("Only greedy0 and reassignment1..passes are assignment phases: "+stage+"/"+pass);
	}
	private void requireActive(){if(!active){throw new IllegalStateException("A pair or finish requires an active query");}}
	private static void increment(long[] array,int index,long amount){array[index]=Math.addExact(array[index],amount);}
	private static long sum(long[] array){long n=0; for(long x:array){n=Math.addExact(n,x);} return n;}
	private static long rank(long n,int permille){
		assert(n>0 && permille>0 && permille<=1000) : "Nearest-rank quantiles require nonempty data and a positive probability";
		return Math.addExact(Math.multiplyExact(n/1000,permille),(n%1000*permille+999)/1000);
	}
	private static int quantile(long[] values,long rank){
		long n=0; for(int b=0;b<values.length;b++){n=Math.addExact(n,values[b]); if(n>=rank){return b;}}
		throw new IllegalStateException("Quantile rank exceeds histogram population");
	}
	private static ByteBuilder row(int p,String cohort,String metric){return new ByteBuilder().append(p==0 ? "greedy" : "reassign").tab().append(p).tab().append(cohort).tab().append(metric);}
	private void bounds(ByteBuilder bb,int metric,int b){
		bb.tab().append(edge(metric,b),6).tab().append(edge(metric,exact(metric) || b==last[metric] ? b : b+1),6).tab().append(exact(metric) || b==last[metric] ? "true" : "false");
	}

	private void validate(){
		for(int p=0;p<=passes;p++){
			final long[] c=counts[p];
			if(c[ALIGNED]!=Math.addExact(Math.addExact(c[CHOSEN_COUNT],c[REJECTED_COUNT]),c[NONWINNERS])){throw new IllegalStateException("Aligned candidate partition failed in phase "+p);}
			if(c[RANKED]!=Math.addExact(c[FILTERED],c[ALIGNED])){throw new IllegalStateException("Ranked/filter/aligned partition failed in phase "+p);}
			if(c[FILTERED]!=Math.addExact(c[HIT_FILTERED],c[HYDRO_FILTERED])){throw new IllegalStateException("Disjoint filter reasons must conserve prefiltered candidates in phase "+p);}
			if(c[CHOSEN_SIZE10]>c[CHOSEN_SIZE5] || c[CHOSEN_SIZE5]>c[CHOSEN_COUNT]){throw new IllegalStateException("Conditional chosen cohorts must be nested in phase "+p);}
			long q=0,w=0,ranked=0,aligned=0;
			for(int n=0;n<=shortlist;n++){
				aligned=Math.addExact(aligned,Math.multiplyExact((long)n,alignedCounts[p][n]));
				for(int r=0;r<=n;r++){q=Math.addExact(q,joint[p][n][r]); ranked=Math.addExact(ranked,Math.multiplyExact((long)n,joint[p][n][r])); if(r>0){w=Math.addExact(w,joint[p][n][r]);}}
			}
			if(q!=c[QUERIES] || w!=c[CHOSEN_COUNT] || sum(alignedCounts[p])!=q || ranked!=c[RANKED] || aligned!=c[ALIGNED]){throw new IllegalStateException("Rank/aligned candidate context does not conserve queries in phase "+p);}
			for(int cohort=0;cohort<COHORTS.length;cohort++){for(int m=0;m<METRICS;m++){
				if(!hasMetric(cohort,m)){continue;}
				if(Math.addExact(sum(bins[p][cohort][m]),missing[p][cohort][m])!=c[countIndex(cohort)]){throw new IllegalStateException("Metric denominator failed: "+p+"/"+cohort+"/"+NAMES[m]);}
			}}
		}
	}

	private void writeBins(ByteStreamWriter out){
		for(int p=0;p<=passes;p++){for(int c=0;c<COHORTS.length;c++){for(int m=0;m<METRICS;m++){
			if(!hasMetric(c,m)){continue;}
			for(int b=0;b<=last[m];b++){if(bins[p][c][m][b]>0){
				final ByteBuilder bb=row(p,COHORTS[c],NAMES[m]); bounds(bb,m,b); out.println(bb.tab().append(bins[p][c][m][b]));
			}}
			if(missing[p][c][m]>0){out.println(row(p,COHORTS[c],NAMES[m]).append("\tNA\tNA\tNA\t").append(missing[p][c][m]));}
		}}}
	}
	private void writeQuantiles(ByteStreamWriter out){
		for(int p=0;p<=passes;p++){for(int c=0;c<COHORTS.length;c++){for(int m=0;m<METRICS;m++){
			if(!hasMetric(c,m)){continue;}
			final long n=sum(bins[p][c][m]);
			for(int permille:QUANTILES){
				final ByteBuilder bb=row(p,COHORTS[c],NAMES[m]).tab().append(permille).tab().append(n).tab().append(missing[p][c][m]);
				if(n==0){bb.append("\tNA\tNA\tNA");}else{bounds(bb,m,quantile(bins[p][c][m],rank(n,permille)));}
				out.println(bb);
			}
		}}}
	}
	private void writeCounts(ByteStreamWriter out){
		for(int p=0;p<=passes;p++){
			for(int k=0;k<COUNT_NAMES.length;k++){out.println(row(p,COUNT_NAMES[k],"NA").append("\tNA\t").append(counts[p][k]));}
			out.println(row(p,"queries_without_winner","NA").append("\tNA\t").append(counts[p][QUERIES]-counts[p][CHOSEN_COUNT]));
			for(int n=0;n<=shortlist;n++){
				final long total=sum(joint[p][n]);
				if(total>0){out.println(row(p,"candidate_count",Integer.toString(n)).append("\tNA\t").append(total));}
				if(alignedCounts[p][n]>0){out.println(row(p,"aligned_candidate_count",Integer.toString(n)).append("\tNA\t").append(alignedCounts[p][n]));}
				for(int r=0;r<=n;r++){if(joint[p][n][r]>0){out.println(row(p,"rank_candidate_joint",Integer.toString(n)).tab().append(r).tab().append(joint[p][n][r]));}}
			}
		}
	}

	/** At an inclusive <= edge, a floating bin may contain equality plus values
	 * above it. Report both possibilities; >= edges are exact because bins are left closed. */
	private long surviving(int p,int cohort,int metric,int thresholdBin,boolean upperBound){
		long n=0;
		for(int b=0;b<=last[metric];b++){
			final boolean include=metric<=COSINE ? b>=thresholdBin :
				b<thresholdBin || (b==thresholdBin && (upperBound || exact(metric) || b==last[metric]));
			if(include){n=Math.addExact(n,bins[p][cohort][metric][b]);}
		}
		return n;
	}
	private void writeWrongSide(ByteStreamWriter out){
		for(int p=0;p<=passes;p++){for(int m=0;m<RANK;m++){
			final long chosen=sum(bins[p][CHOSEN][m]), rejected=sum(bins[p][REJECTED][m]);
			if(chosen==0){continue;}
			for(int permille:QUANTILES){
				final int b=quantile(bins[p][CHOSEN][m],rank(chosen,permille));
				final int edges=exact(m) || b==last[m] ? 1 : 2;
				for(int e=0;e<edges;e++){
					final int cut=b+e;
					final ByteBuilder bb=row(p,"chosen",NAMES[m]).tab().append(permille);
					bounds(bb,m,b);
					bb.tab().append(e==0 ? "lower" : "upper").tab().append(edge(m,cut),6).tab().append(m<=COSINE ? ">=" : "<=");
					bb.tab().append(chosen).tab().append(missing[p][CHOSEN][m]).tab().append(rejected).tab().append(missing[p][REJECTED][m]);
					bb.tab().append(chosen-surviving(p,CHOSEN,m,cut,true)).tab().append(chosen-surviving(p,CHOSEN,m,cut,false));
					bb.tab().append(surviving(p,REJECTED,m,cut,false)).tab().append(surviving(p,REJECTED,m,cut,true));
					out.println(bb);
				}
			}
		}}
	}

	public static final String[] SUFFIXES={".bins.tsv",".quantiles.tsv",".counts.tsv",".wrong_side.tsv"};
	private static final String[] HEADERS={
		"stage\tpass\tcohort\tmetric\tlower\tupper\tupper_inclusive\tcount",
		"stage\tpass\tcohort\tmetric\tquantile_permille\tavailable\tunavailable\tlower\tupper\tupper_inclusive",
		"stage\tpass\tcount_type\tcandidates\twinning_rank\tcount",
		"stage\tpass\tcutoff_cohort\tmetric\tquantile_permille\tquantile_lower\tquantile_upper\tquantile_upper_inclusive\tedge\tthreshold\toperator\tchosen_available\tchosen_unavailable\trejected_available\trejected_unavailable\tchosen_losses_min\tchosen_losses_max\trejected_surviving_min\trejected_surviving_max"};
	private static final String[] NAMES={"ddl_hits","monomer_cosine","abs_hydrophobicity_delta","abs_charge_delta","abs_aromaticity_delta","abs_small_delta","abs_pi_delta","winning_rank"};
	private static final String[] COHORTS={"chosen","rejected","chosen_refsize_ge5","chosen_refsize_ge10"};
	private static final String[] COUNT_NAMES={"queries","chosen_tails","rejected_pairs","passing_nonwinner_pairs","aligned_pairs","ranked_candidates","prefiltered_candidates","chosen_refsize_ge5_tails","chosen_refsize_ge10_tails","hit_filtered_candidates","hydro_filtered_candidates"};
	private static final int[] SCALES={1,1000,100,1000,1000,1000,100,1}, QUANTILES={1,10,50,100,900,950,990,999};
	private static final int CHOSEN=0,REJECTED=1,SIZE5=2,SIZE10=3,HITS=0,COSINE=1,RANK=7,METRICS=8;
	private static final int QUERIES=0,CHOSEN_COUNT=1,REJECTED_COUNT=2,NONWINNERS=3,ALIGNED=4,RANKED=5,FILTERED=6,CHOSEN_SIZE5=7,CHOSEN_SIZE10=8,HIT_FILTERED=9,HYDRO_FILTERED=10;
	private final int passes,width,shortlist;
	private final int[] last,scalarBins=new int[5];
	private final long[][][][] bins;
	private final long[][][] missing,joint;
	private final long[][] counts,alignedCounts;
	private final double[] queryScalars=new double[5];
	private byte[] raw=new byte[0];
	private int phase,queryMask,queryRejected,queryRank;
	private boolean active;
}
