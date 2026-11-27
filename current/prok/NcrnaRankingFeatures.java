package prok;

import java.util.Arrays;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Read-only training profiles and worker-owned window composition for ranking.
 * Profile ordinals must equal the current caller's full model axis. Missing
 * support is independent for GC, gene dimers and source-contig dimers.
 * No weights are fitted here. WindowRanking is an explicitly installed frozen
 * window-only order; the ordinary observation lifecycle remains independent.
 * @author Brian, Ganyu, Raiden
 */
public final class NcrnaRankingFeatures{

	/*--------------------------------------------------------------*/
	/*----------------       Training Profiles      ----------------*/
	/*--------------------------------------------------------------*/

	public static Profiles load(final String path, final String[] modelNames){
		return load(path, modelNames, false);
	}

	/** Current-window vectors have their own exact schema. Never silently use
	 * historical gene-core vectors as a protected-window reference profile. */
	public static Profiles loadWindows(final String path, final String[] modelNames){
		return load(path, modelNames, true);
	}

	private static Profiles load(final String path, final String[] modelNames, final boolean windows){
		require(path!=null && modelNames!=null && modelNames.length>0,
			"Ranking profiles require a file and the complete current caller model axis");
		final Profiles result=new Profiles(modelNames.length, windows);
		final LineParser1 p=new LineParser1('\t');
		int model=0;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, windows ? WINDOW_HEADER : HEADER, path);
			final int shift=windows ? 1 : 0;
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==41+shift && model<modelNames.length && p.termEquals(modelNames[model], 0),
					"Ranking profiles must preserve every model exactly once in its caller ordinal");
				final long gcSupport=p.parseLong(1);
				final long coreSupport=p.parseLong(3+shift);
				final long coreUnknown=p.parseLong(4+shift);
				final long contigSupport=p.parseLong(5+shift);
				final long contigUnknown=p.parseLong(6+shift);
				final long contigs=p.parseLong(39+shift);
				final long knownContigs=p.parseLong(40+shift);
				require(gcSupport>=0 && coreSupport>=0 && coreUnknown>=0 && contigSupport>=0
					&& contigUnknown>=0 && contigs>=0 && knownContigs>=0 && knownContigs<=contigs
					&& knownContigs<=contigSupport && contigs<=Math.addExact(contigSupport, contigUnknown)
					&& contigs-knownContigs<=contigUnknown
					&& (contigs==0)==(contigSupport==0 && contigUnknown==0)
					&& (knownContigs==0)==(contigSupport==0)
					&& Math.addExact(coreSupport, coreUnknown)==Math.addExact(contigSupport, contigUnknown),
					"Feature supports must preserve equal original/model incidences and actual distinct-contig counts");
				result.gcKnown[model]=gcSupport>0;
				if(gcSupport>0){result.gc[model]=fraction(p, 2);}
				else{require(p.termEquals("NA", 2), "Unsupported GC is explicitly NA, not a measured zero");}
				if(windows){
					require(gcSupport<=Math.addExact(coreSupport, coreUnknown), "Current GC and dimer support must describe the same original/window incidences");
					if(gcSupport==0){require(p.termEquals("NA", 3), "Unsupported GC standard deviation is NA");}
					else{
						result.gcSd[model]=fraction(p, 3);
						require(result.gcSd[model]<=.5 && (gcSupport>1 || result.gcSd[model]==0), "Population GC deviation is at most one half and zero for one observation");
					}
				}
				result.coreNorm[model]=readVector(p, 7+shift, coreSupport, result.core[model]);
				result.contigNorm[model]=readVector(p, 23+shift, contigSupport, result.contig[model]);
				model++;
			}
		}
		require(model==modelNames.length, "A profile cannot omit zero-support models from the caller axis");
		return result;
	}

	private static double readVector(final LineParser1 p, final int first, final long support, final double[] values){
		double sum=0;
		double squared=0;
		for(int bin=0; bin<16; bin++){
			if(support==0){require(p.termEquals("NA", first+bin), "Every unsupported vector coordinate must be NA");}
			else{
				final double value=fraction(p, first+bin);
				values[bin]=value;
				sum+=value;
				squared+=value*value;
			}
		}
		if(support==0){return 0;}
		require(Math.abs(sum-1)<1e-8 && squared>0, "Supported mean dimer frequencies must sum to one, not raw occurrence counts");
		return 1/Math.sqrt(squared);
	}
	private static double fraction(final LineParser1 p, final int column){
		final double value=p.parseDouble(column);
		require(Double.isFinite(value) && value>=0 && value<=1, "Measured composition coordinates must be finite fractions");
		return value;
	}

	/** Immutable after loading; share across workers, never expose mutable arrays. */
	public static final class Profiles{
		private Profiles(final int models, final boolean currentWindows_){
			currentWindows=currentWindows_;
			gc=new double[models];
			gcSd=new double[models];
			gcKnown=new boolean[models];
			core=new double[models][16];
			contig=new double[models][16];
			coreNorm=new double[models];
			contigNorm=new double[models];
		}
		public int models(){return gc.length;}
		private final boolean currentWindows;
		private final double[] gc, gcSd, coreNorm, contigNorm;
		private final boolean[] gcKnown;
		private final double[][] core, contig;
	}

	/** Frozen hits/1024 * exp(-(1-cosine)) order over the unchanged eligible set.
	 * Missing model/window support is neutral. One window scan, double scores and
	 * primitive heap scratch are reused; no GC or source-contig feature is read.
	 * The returned array is borrowed and only its supplied count-prefix is live. */
	public static final class WindowRanking{
		public WindowRanking(final Profiles profiles_){
			require(profiles_!=null && profiles_.models()>0 && profiles_.currentWindows,
				"Window ranking requires the complete current-window profile axis, never a historical gene-core profile");
			profiles=profiles_;
			weights=new double[profiles.models()];
			order=new int[profiles.models()];
		}
		public int[] build(final byte[] bases, final float[] raw, final int[] eligible, final int count){
			require(bases!=null && raw!=null && raw.length==order.length && eligible!=null
				&& count>=0 && count<=eligible.length && count<=order.length,
				"Window ranking borrows complete raw DDL counts and one unchanged eligible prefix");
			if(count==0){return order;}
			composition.scan(bases, 0, bases.length);
			final boolean known=composition.frequencies(window);
			final double queryNorm=known ? norm(window) : 0;
			for(int i=0; i<count; i++){
				final int model=eligible[i];
				require(model>=0 && model<order.length, "Window ranks retain the original complete model ordinals");
				final float hits=raw[model];
				require(Float.isFinite(hits) && hits>=0 && hits<=1024 && hits==(int)hits,
					"Frozen window weighting consumes raw integer hits from the 1024-bin pair, never an already weighted score");
				final double h=hits==0 ? 0 : hits/1024.0;
				weights[model]=known && profiles.coreNorm[model]>0
					? h*Math.exp(-distance(window, queryNorm, profiles.core[model], profiles.coreNorm[model])) : h;
				order[i]=model;
			}
			for(int i=count/2-1; i>=0; i--){sift(i, count);}
			for(int end=count-1; end>0; end--){
				final int model=order[end];
				order[end]=order[0];
				order[0]=model;
				sift(0, end);
			}
			return order;
		}
		private void sift(int root, final int count){
			assert(root>=0 && root<count && count<=order.length)
				: "Only the live eligible prefix forms the primitive window-ranking heap";
			while(2*root+1<count){
				int child=2*root+1;
				if(child+1<count && compare(order[child], order[child+1])<0){child++;}
				if(compare(order[root], order[child])>=0){return;}
				final int model=order[root];
				order[root]=order[child];
				order[child]=model;
				root=child;
			}
		}
		private int compare(final int a, final int b){
			assert(a>=0 && b>=0 && a<weights.length && b<weights.length)
				: "Heap entries are validated original model ordinals, not positions in a capped list";
			final int score=Double.compare(weights[b], weights[a]);
			return score==0 ? Integer.compare(a, b) : score;
		}
		private final Profiles profiles;
		private final NcrnaDimerComposition composition=new NcrnaDimerComposition();
		private final double[] window=new double[16];
		private final double[] weights;
		private final int[] order;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Worker Features        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaRankingFeatures(final Profiles profiles_){
		require(profiles_!=null, "Each ranking-feature worker needs an explicitly bound training profile");
		profiles=profiles_;
	}

	/** Call once per new source contig in its original forward orientation, even
	 * when the caller reuses a bases array. Identity-caching a mutable byte array
	 * would silently reuse a previous contig or the wrong orientation.
	 */
	public void beginContig(final byte[] forwardBases){
		contigReady=false;
		windowReady=false;
		require(forwardBases!=null, "A new contig requires its actual forward bases");
		contigComposition.scan(forwardBases, 0, forwardBases.length);
		contigKnown=contigComposition.frequencies(forward);
		// Reverse-complement permutation preserves the norm. Compute it once
		// per source contig, not again for each window on that contig.
		contextNorm=norm(forward);
		Arrays.fill(reverse, 0);
		for(int bin=0; bin<16; bin++){
			final int rc=((3-(bin&3))<<2)|(3-(bin>>>2));
			reverse[rc]=forward[bin];
		}
		contigReady=true;
	}

	/** Window bases are the actual oriented query borrowed from the verifier.
	 * This performs one composition scan per window, not one scan per model.
	 */
	public void beginWindow(final byte[] orientedWindow, final int strand){
		windowReady=false;
		require(contigReady && orientedWindow!=null && (strand==0 || strand==1),
			"Window features require the current source contig and explicit orientation");
		windowComposition.scan(orientedWindow, 0, orientedWindow.length);
		windowKnown=windowComposition.frequencies(window);
		windowGcKnown=windowComposition.knownBases()>0;
		windowGc=windowGcKnown ? windowComposition.gcBases()/(double)windowComposition.knownBases() : 0;
		windowNorm=norm(window);
		context=strand==0 ? forward : reverse;
		windowReady=true;
	}

	public boolean hasGc(final int model){check(model);return windowGcKnown && profiles.gcKnown[model];}
	public boolean hasWindowDimers(final int model){check(model);return windowKnown && profiles.coreNorm[model]>0;}
	public boolean hasContigDimers(final int model){check(model);return contigKnown && profiles.contigNorm[model]>0;}
	public double gcDistance(final int model){
		return Math.abs(gcDifference(model));
	}
	/** Signed window minus reference GC, in fractions rather than percent. */
	public double gcDifference(final int model){
		require(hasGc(model), "Missing window/model GC must be reported as unavailable, not a zero distance");
		return windowGc-profiles.gc[model];
	}
	public double windowDimerDistance(final int model){
		require(hasWindowDimers(model), "Missing gene/window dimers must retain their support mask");
		return distance(window, windowNorm, profiles.core[model], profiles.coreNorm[model]);
	}
	public double contigDimerDistance(final int model){
		require(hasContigDimers(model), "Missing source-contig dimers must retain their support mask");
		return distance(context, contextNorm, profiles.contig[model], profiles.contigNorm[model]);
	}
	private void check(final int model){
		require(windowReady && model>=0 && model<profiles.models(), "Feature access requires the current window and its original model ordinal");
	}
	private static double norm(final double[] vector){
		double sum=0;
		for(final double value:vector){sum+=value*value;}
		return sum==0 ? 0 : 1/Math.sqrt(sum);
	}
	private static double distance(final double[] a, final double aNorm, final double[] b, final double bNorm){
		double dot=0;
		for(int i=0; i<16; i++){dot+=a[i]*b[i];}
		final double cosine=dot*aNorm*bNorm;
		assert(cosine>=0 && cosine<=1+1e-12) : "Finite nonnegative normalized vectors bound cosine within roundoff of [0,1]";
		return 1-Math.min(1, cosine);
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Profiles profiles;
	private final NcrnaDimerComposition contigComposition=new NcrnaDimerComposition();
	private final NcrnaDimerComposition windowComposition=new NcrnaDimerComposition();
	private final double[] forward=new double[16], reverse=new double[16], window=new double[16];
	private double[] context;
	private double windowGc, windowNorm, contextNorm;
	private boolean contigReady, windowReady, windowKnown, windowGcKnown, contigKnown;
	public static final String HEADER=header();
	public static final String WINDOW_HEADER=windowHeader();
	private static String windowHeader(){
		final StringBuilder b=new StringBuilder("model\twindow_gc_support\tmean_window_gc\tsd_window_gc\twindow_support\twindow_unknown\tcontig_support\tcontig_unknown");
		for(int i=0; i<16; i++){b.append("\twindow_").append(i);}
		for(int i=0; i<16; i++){b.append("\tcontig_").append(i);}
		return b.append("\tdistinct_source_contigs\tknown_source_contigs").toString();
	}
	private static String header(){
		final StringBuilder b=new StringBuilder("model\twindow_gc_support\tmean_window_gc\tcore_support\tcore_unknown\tcontig_support\tcontig_unknown");
		for(int i=0; i<16; i++){b.append("\tcore_").append(i);}
		for(int i=0; i<16; i++){b.append("\tcontig_").append(i);}
		return b.append("\tdistinct_source_contigs\tknown_source_contigs").toString();
	}
}
