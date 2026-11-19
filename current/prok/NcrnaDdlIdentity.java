package prok;

import java.util.Arrays;

import idaligner.QuantumAligner;
import idaligner.Tracer;
import parse.LineParser1;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;

/** Identity to one frozen native-selected accepted consensus, using saved genes.
 * Model choice never examines new sketch ranks or maximizes measured identity.
 * The owner calls this before scoring workers and before deriving identity bins.
 * @author Brian, Raiden
 */
public final class NcrnaDdlIdentity{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaDdlIdentity(){}

	/** One retained result per original parent, including explicit missing values. */
	public static Locus[] measure(final NcrnaDdlSavedData data, final byte[][] references){
		return extend(data, references, null, null);
	}

	/** Reuses accepted identities only for identical source/parent/witness bindings.
	 * Previous values must have been measured or read against these references.
	 * New parents alone reach the aligner; old missing values remain explicit. */
	public static Locus[] extend(final NcrnaDdlSavedData data, final byte[][] references,
			final NcrnaDdlSavedData previousData, final Locus[] previous){
		require(data!=null && references!=null && references.length==data.models,
			"Identity must use the complete original model order bound by the saved-cache owner");
		require(previousData==null ? previous==null : previous!=null && previous.length==previousData.parents.size()
			&& previousData.models==data.models, "Reused identities require the complete previously validated parent/model roster");
		final Locus[] loci=witnesses(data);
		final byte[][] normalized=new byte[references.length][];
		for(int model=0; model<references.length; model++){
			require(references[model]!=null && references[model].length>0,
				"Every original consensus must retain its complete nonempty sequence for identity binding");
			normalized[model]=normalize(references[model]);
		}
		final int[] positions=new int[4];
		for(final Parent parent:data.parents){
			final Locus locus=loci[parent.index];
			require(parent.nativeRank<=0 || locus.model>=0,
				"A positive native rank needs its accepted original event/model before identity measurement or reuse");
			final Locus old=previousData==null ? null : previous(parent, previousData, previous);
			if(old!=null){
				require(old.model==locus.model && old.modelEvent==locus.modelEvent && old.geneEvent==locus.geneEvent
					&& old.containingWindows==locus.containingWindows && old.definedWindows==locus.definedWindows,
					"Reused identity must retain its exact native model/event and complete-gene witnesses");
				final String status=parent.nativeRank<=0 ? NOT_NATIVE_RANKABLE : locus.gene==null
					? locus.containingWindows>0 ? UNDEFINED_GENE : NO_COMPLETE_WINDOW : normalized[locus.model]==null
					? UNDEFINED_CONSENSUS : old.status.equals(MEASURED) ? MEASURED : INVALID_ALIGNMENT;
				require(old.status.equals(status), "Reused identity availability must agree with the unchanged saved gene and consensus");
				locus.identity=old.identity;
				locus.status=old.status;
				locus.reused=true;
				locus.gene=null;
				continue;
			}
			if(parent.nativeRank<=0){continue;}
			require(locus.model>=0, "A positive saved native locus minimum must have its original accepted event/model witness");
			if(locus.gene==null){
				locus.status=locus.containingWindows>0 ? UNDEFINED_GENE : NO_COMPLETE_WINDOW;
			}else if(normalized[locus.model]==null){
				locus.status=UNDEFINED_CONSENSUS;
			}else{
				locus.identity=align(normalized[locus.model], locus.gene, positions);
				locus.status=Float.isNaN(locus.identity) ? INVALID_ALIGNMENT : MEASURED;
			}
			locus.gene=null;
		}
		return loci;
	}

	/** Cold metadata join, using the existing source-panel parent maps. */
	private static Locus previous(final Parent parent, final NcrnaDdlSavedData data, final Locus[] values){
		for(final NcrnaDdlSavedData.Panel panel:data.panels){
			if(!panel.prefix.equals(parent.input.panel.prefix)){continue;}
			final Parent old=panel.parents.get(new String(parent.idBytes, java.nio.charset.StandardCharsets.UTF_8));
			if(old==null){return null;}
			require(old.input.slot==parent.input.slot && old.input.group==parent.input.group
				&& old.input.id.equals(parent.input.id) && old.input.path.equals(parent.input.path)
				&& Arrays.equals(old.seqidBytes, parent.seqidBytes) && old.start==parent.start && old.stop==parent.stop
				&& old.strand==parent.strand && old.stratum==parent.stratum && old.nativeRank==parent.nativeRank,
				"Identity reuse cannot exchange the source assembly, parent coordinates, seen label or native minimum");
			require(values[old.index]!=null, "A prior identity roster must retain an explicit value or missing disposition for every bound parent");
			return values[old.index];
		}
		return null;
	}

	/** Reconstructs saved-cache provenance and gene availability without alignment.
	 * Both the original measurement and frozen-value reader use this same routing;
	 * only measure() may call the aligner after these witnesses are established. */
	private static Locus[] witnesses(final NcrnaDdlSavedData data){
		require(data!=null && data.models>0, "Identity witnesses require the validated complete saved-cache model roster");
		final Locus[] loci=new Locus[data.parents.size()];
		for(final Parent parent:data.parents){
			require(parent.index>=0 && parent.index<loci.length && loci[parent.index]==null,
				"Identity rows must preserve every original selected parent exactly once");
			loci[parent.index]=new Locus(parent.nativeRank>0 ? NO_COMPLETE_WINDOW : NOT_NATIVE_RANKABLE);
		}
		for(final Window window:data.windows){
			for(final Edge edge:window.edges){
				final Parent parent=edge.parent;
				if(parent.nativeRank<=0){continue;}
				final Locus locus=loci[parent.index];
				if(edge.nativeRank==parent.nativeRank){
					final int model=edge.nativeTuple[3];
					require(model>=0 && model<data.models && (edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0,
						"A frozen native-selected identity model must be uniquely accepted for that original event/parent");
					if(locus.model<0 || window.event<locus.modelEvent
						|| window.event==locus.modelEvent && model<locus.model){
						locus.model=model;
						locus.modelEvent=window.event;
					}
				}
				gene(window, parent, locus);
			}
		}
		return loci;
	}

	/** Reads accepted identity values once, binding every original parent and witness.
	 * Only saved-window extraction is repeated to establish provenance/availability;
	 * no model is selected from new scores and no alignment is performed.
	 * Bin/side arrays carry the stored labels to the cohort reader's cross-check. */
	static Locus[] read(final String path, final NcrnaDdlSavedData data, final byte[][] names,
			final byte[][] references, final int[] bins, final int[] sides){
		require(path!=null && data!=null && names!=null && names.length==data.models
			&& references!=null && references.length==data.models && bins!=null && sides!=null
			&& bins.length==data.parents.size() && sides.length==bins.length,
			"Saved identities require the same complete parent/model roster and cohort output arrays");
		final Locus[] loci=witnesses(data);
		final boolean[] defined=new boolean[references.length];
		for(int model=0; model<references.length; model++){
			require(names[model]!=null && names[model].length>0 && references[model]!=null && references[model].length>0,
				"Every frozen identity model needs its original label and consensus bytes");
			defined[model]=true;
			for(final byte base:references[model]){defined[model]&=RrnaResourceIO.base(base)>=0;}
		}
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaDdlCohorts.IDENTITY_HEADER, path);
			for(final Parent parent:data.parents){
				final byte[] row=in.reader.nextLine();
				check(row!=null, "Saved identities omit an original parent", path);
				p.set(row);
				check(p.terms()==18 && field(p, 0, parent.input.panel.nameBytes) && p.parseInt(1)==parent.input.slot
					&& field(p, 2, parent.input.idBytes) && field(p, 3, parent.idBytes)
					&& p.termEquals(NcrnaDdlSavedData.STRATA[parent.stratum], 4) && p.parseInt(5)==parent.input.group
					&& p.termEquals(parent.nativeRank>0 ? "true" : "false", 6),
					"Saved identity parent binding differs from the accepted panel/input/parent/stratum/group", path);
				final Locus locus=loci[parent.index];
				check(parent.nativeRank<=0 || locus.model>=0,
					"A positive native minimum needs an accepted original model/event before restoring identity", path);
				check(p.parseInt(7)==locus.model && (locus.model<0 ? p.termEquals(".", 8) : field(p, 8, names[locus.model]))
					&& p.parseLong(9)==locus.modelEvent && p.parseLong(10)==locus.geneEvent
					&& p.parseInt(11)==locus.containingWindows && p.parseInt(12)==locus.definedWindows,
					"Saved identity model/event/gene witnesses differ from the original saved windows", path);
				locus.identity=Float.intBitsToFloat(p.parseInt(14));
				final String status;
				if(parent.nativeRank<=0){status=NOT_NATIVE_RANKABLE;}
				else if(locus.gene==null){status=locus.containingWindows>0 ? UNDEFINED_GENE : NO_COMPLETE_WINDOW;}
				else if(!defined[locus.model]){status=UNDEFINED_CONSENSUS;}
				else{status=p.termEquals(MEASURED, 17) ? MEASURED : INVALID_ALIGNMENT;}
				check(p.termEquals(status, 17), "Saved identity status disagrees with original gene/consensus availability", path);
				if(status.equals(MEASURED)){
					check(!Float.isNaN(locus.identity) && !Float.isInfinite(locus.identity) && locus.identity>=0 && locus.identity<=1
						&& !p.termEquals("NA", 13), "Measured saved identity needs finite original float bits in[0,1]", path);
					final double displayed=p.parseDouble(13);
					check(!Double.isNaN(displayed) && !Double.isInfinite(displayed)
						&& Math.abs(displayed-locus.identity)<=0.000000000500001,
						"Saved identity decimal disagrees with the original float bits at nine-place rounding", path);
				}else{
					check(p.termEquals("NA", 13) && Float.isNaN(locus.identity), "Unavailable saved identity must remain explicit NA/NaN", path);
				}
				bins[parent.index]=label(p, 15, NcrnaDdlCohorts.BINS, 5, path);
				sides[parent.index]=p.termEquals("OUTSIDE_NOT_SEEN_ROSTER", 16) ? -1 : p.termEquals("A", 16) ? 0 : p.termEquals("B", 16) ? 1 : -2;
				check(sides[parent.index]>=-1, "Saved identity group side is unavailable or malformed", path);
				locus.status=status;
				locus.gene=null;
			}
			check(in.reader.nextLine()==null, "Saved identities append an extra or duplicate original parent", path);
		}
		return loci;
	}

	private static int label(final LineParser1 p, final int field, final String[] labels, final int count, final String path){
		assert(count<=labels.length) : "Only per-locus labels, never the aggregate ALL bin, may be restored";
		for(int i=0; i<count; i++){if(p.termEquals(labels[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown saved identity cohort label: "+path);
	}
	/** Path context is constructed only for a failure, never once per valid field. */
	private static void check(final boolean valid, final String reason, final String path){
		if(!valid){throw new IllegalArgumentException(reason+": "+path);}
	}
	private static boolean field(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static void gene(final Window window, final Parent parent, final Locus locus){
		assert(window.input==parent.input && window.strand==parent.strand) :
			"NcrnaDdlSavedData.edge binds each identity gene to the original input and sense strand";
		if(parent.start-1<window.start || parent.stop-1>window.stop){return;}
		locus.containingWindows++;
		final byte[] gene=NcrnaDdlLocusOverlap.extract(window.bases, window.start, window.stop,
			window.strand, parent.start, parent.stop, parent.strand);
		if(gene==null){return;}
		locus.definedWindows++;
		if(locus.gene!=null && !Arrays.equals(locus.gene, gene)){
			throw new IllegalArgumentException("Complete defined captured genes for one parent must agree before assigning its identity; parent="
				+parent.index+", event="+window.event);
		}
		if(locus.gene==null || window.event<locus.geneEvent){
			locus.gene=gene;
			locus.geneEvent=window.event;
		}
	}

	/** Normalize once per retained consensus; undefined bases produce explicit NA. */
	private static byte[] normalize(final byte[] bases){
		assert(bases!=null && bases.length>0) : "The model roster was validated before consensus normalization";
		final byte[] result=new byte[bases.length];
		for(int i=0; i<bases.length; i++){
			final int base=RrnaResourceIO.base(bases[i]);
			if(base<0){return null;}
			result[i]=DNA[base];
		}
		return result;
	}

	/** Native semiglobal consensus-to-gene identity, not ANI or member identity.
	 * QuantumAligner:151 retains query direction when positions is nonnull.
	 * IDAlignerStatics:36 reads positions[0] as an input hint: reset on every call.
	 * Tracer:295 defines identity as M/(M+S+I+D) from its packed native score;
	 * use that returned float, not the truncated integer count output slots. */
	public static float align(final byte[] consensus, final byte[] gene, final int[] positions){
		require(consensus!=null && consensus.length>0 && gene!=null && gene.length>0
			&& positions!=null && positions.length>=4,
			"Consensus identity requires nonempty normalized sequences and reusable native endpoint/score scratch");
		require(!Tracer.GLOBAL, "Saved consensus identity uses Tracer's ordinary semiglobal denominator, never a process-global global-alignment override");
		Arrays.fill(positions, 0);
		positions[0]=Integer.MIN_VALUE;
		final float identity=QuantumAligner.alignStatic(consensus, gene, positions);
		return !Float.isNaN(identity) && !Float.isInfinite(identity) && identity>=0 && identity<=1
			&& positions[0]>=0 && positions[1]>=positions[0] && positions[1]<gene.length ? identity : Float.NaN;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Locus{
		private Locus(final String status_){status=status_;}
		public int model(){return model;}
		public long modelEvent(){return modelEvent;}
		public long geneEvent(){return geneEvent;}
		public int containingWindows(){return containingWindows;}
		public int definedWindows(){return definedWindows;}
		public float identity(){return identity;}
		public String status(){return status;}
		public boolean reused(){return reused;}
		private int model=-1, containingWindows, definedWindows;
		private long modelEvent=-1, geneEvent=-1;
		private float identity=Float.NaN;
		private String status;
		private byte[] gene;
		private boolean reused;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final String MEASURED="MEASURED", NOT_NATIVE_RANKABLE="NOT_NATIVE_RANKABLE";
	public static final String NO_COMPLETE_WINDOW="NO_COMPLETE_WINDOW", UNDEFINED_GENE="UNDEFINED_GENE";
	public static final String UNDEFINED_CONSENSUS="UNDEFINED_CONSENSUS", INVALID_ALIGNMENT="INVALID_ALIGNMENT";
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
}
