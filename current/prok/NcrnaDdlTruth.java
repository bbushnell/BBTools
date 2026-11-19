package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import map.ObjectMap;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaSubtypeEvidence.Support;
import structures.ByteBuilder;

/** Binds the DDL experiment to the registered reference and subtype witnesses.
 * One owner loads the selected inputs before caller workers start. Each physical
 * cobalamin parent is qualified once through NcrnaSubtypeEvidence; no caller or
 * sketch score can add truth. Masks, protected priors and multiple scored subtype
 * supports retain their original disposition. RF00174 scored parents define the
 * target denominator; all scored cobalamin parents participate in conservative
 * unique-span attribution, including other subtypes, and masks never consume it.
 * Returned rows/arrays are borrowed read-only. This loads metadata, not genomes.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlTruth{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Joins exact selected FASTA paths in their caller order, including empty inputs.
	 * Unselected reference GFFs are never opened. Subtype support is derived once
	 * per selected cobalamin parent and the complete witness population is checked
	 * once after loading. Repeated or absent selected inputs fail loudly.
	 */
	public static Input[] load(final String manifest, final String[] selected,
			final String subtypes, final String witnesses, final String events, final Rule[] rules){
		require(selected!=null && selected.length>0, "Truth binding requires the complete selected input population");
		final ObjectMap<String, Binding> byFasta=new ObjectMap<String, Binding>(String.class, Binding.class);
		for(final Binding binding:GradeNcrnaReference.loadManifest(manifest)){byFasta.put(binding.fasta, binding);}
		final Binding[] bindings=new Binding[selected.length];
		final ArrayList<String> ids=new ArrayList<String>(selected.length);
		for(int i=0; i<selected.length; i++){
			require(selected[i]!=null && !selected[i].isEmpty(), "A selected input needs its exact original FASTA path");
			final Binding binding=byFasta.get(selected[i]);
			require(binding!=null, "Each exact selected FASTA must occur in the original reference manifest");
			bindings[i]=binding;
			ids.add(binding.id);
			byFasta.remove(selected[i]);
		}
		final NcrnaSubtypeEvidence evidence=new NcrnaSubtypeEvidence(subtypes, witnesses, events, rules, ids);
		final byte[][] labels=new byte[evidence.size()][];
		int target=-1;
		for(int type=0; type<evidence.size(); type++){
			labels[type]=evidence.referenceModel(type).getBytes(StandardCharsets.UTF_8);
			if("RF00174".equals(evidence.referenceModel(type))){
				require(target<0, "RF00174 must identify exactly one qualified target subtype");
				target=type;
			}
		}
		require(target>=0, "The explicit subtype map must contain RF00174");
		final Input[] result=new Input[bindings.length];
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ByteBuilder names=new ByteBuilder();
		for(int input=0; input<bindings.length; input++){
			result[input]=bind(bindings[input], codec, evidence, target, labels, names);
		}
		evidence.finish();
		return result;
	}

	private NcrnaDdlTruth(){}

	/*--------------------------------------------------------------*/
	/*----------------        Parent Binding        ----------------*/
	/*--------------------------------------------------------------*/

	private static Input bind(final Binding binding, final NcrnaReferenceGff codec,
			final NcrnaSubtypeEvidence evidence, final int target, final byte[][] labels, final ByteBuilder names){
		assert(binding!=null && target>=0) : "Each reference read must bind an exact selected input and the explicit target subtype";
		final ArrayList<Parent> parents=new ArrayList<Parent>();
		final ObjectMap<String, GroupBuilder> groups=new ObjectMap<String, GroupBuilder>(String.class, GroupBuilder.class);
		final ArrayList<GroupBuilder> order=new ArrayList<GroupBuilder>();
		for(final Row row:codec.read(binding.gff, binding.id)){
			if(!"cobalamin".equals(row.locus.family)){continue;}
			final Support support=evidence.support(row);
			GroupBuilder group=groups.get(row.locus.seqid);
			if(group==null){
				group=new GroupBuilder(row.locus.seqid);
				groups.put(row.locus.seqid, group);
				order.add(group);
			}
			final boolean selected=!row.mask && (support.scored&(1<<target))!=0;
			final Parent parent=new Parent(row, selected, support.priorType>=0, support.multi(),
				names(labels, support.scored, names), names(labels, support.marginal, names),
				row.mask ? -1 : group.scored.size());
			parents.add(parent);
			group.parents.add(parent);
			if(!row.mask){group.scored.add(row.locus);}
			if(selected){group.targets.add(parent);}
		}
		final Group[] finished=new Group[order.size()];
		for(int i=0; i<finished.length; i++){finished[i]=new Group(order.get(i));}
		return new Input(binding, parents.toArray(new Parent[parents.size()]), finished);
	}

	/** Reference-model IDs describe support without exposing unexplained bit masks. */
	private static String names(final byte[][] labels, final int bits, final ByteBuilder out){
		assert(labels!=null && out!=null) : "Support names use the once-encoded explicit mapping and one loader-owned output scratch";
		out.clear();
		for(int type=0; type<labels.length; type++){
			if((bits&(1<<type))==0){continue;}
			if(out.length()>0){out.append(',');}
			out.append(labels[type]);
		}
		return out.length()==0 ? "." : new String(out.array, 0, out.length(), StandardCharsets.UTF_8);
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Bound Metadata        ----------------*/
	/*--------------------------------------------------------------*/

	/** One selected input, retaining zero-parent and zero-target populations. */
	public static final class Input{
		private Input(final Binding binding_, final Parent[] parents_, final Group[] ordered){
			assert(binding_!=null && parents_!=null && ordered!=null)
				: "Truth inputs are published only after the exact selected manifest and complete parent groups bind";
			binding=binding_;
			parents=parents_;
			int count=0, masks=0;
			for(final Parent parent:parents){
				if(parent.target){count++;}
				if(parent.row.mask){masks++;}
			}
			targets=count;
			masked=masks;
			for(final Group group:ordered){groups.put(group.contig, group);}
		}
		/** Exact native GFF seqid, scoped to this input; null means no cobalamin parent. */
		public Group group(final String contig){return groups.get(contig);}
		public final Binding binding;
		public final Parent[] parents;
		public final int targets, masked;
		private final ObjectMap<String, Group> groups=new ObjectMap<String, Group>(String.class, Group.class);
	}

	/** One physical parent, with immutable qualification in the native reference view. */
	public static final class Parent{
		private Parent(final Row row_, final boolean target_, final boolean prior_, final boolean multi_,
				final String scored_, final String marginal_, final int scoredIndex_){
			row=row_;
			target=target_;
			prior=prior_;
			multi=multi_;
			scoredModels=scored_;
			marginalModels=marginal_;
			scoredIndex=scoredIndex_;
			assert((scoredIndex>=0)!=row.mask && (!target || !row.mask))
				: "Only scored physical parents enter span attribution; target recovery cannot promote a mask";
		}
		public final Row row;
		public final boolean target, prior, multi;
		public final String scoredModels, marginalModels;
		/** Index in this contig's all-subtype scored array; masks retain -1. */
		public final int scoredIndex;
	}

	/** Both strands share one contig group; the oracle checks the original strand. */
	public static final class Group{
		private Group(final GroupBuilder builder){
			contig=builder.contig;
			parents=builder.parents.toArray(new Parent[builder.parents.size()]);
			targets=builder.targets.toArray(new Parent[builder.targets.size()]);
			scored=builder.scored.toArray(new Locus[builder.scored.size()]);
			for(final Parent parent:parents){
				assert(parent.scoredIndex<0 || scored[parent.scoredIndex]==parent.row.locus)
					: "Unique attribution must point back to the same physical parent, not a subtype witness";
			}
		}
		public final String contig;
		public final Parent[] parents, targets;
		public final Locus[] scored;
	}

	private static final class GroupBuilder{
		GroupBuilder(final String contig_){contig=contig_;}
		final String contig;
		final ArrayList<Parent> parents=new ArrayList<Parent>(), targets=new ArrayList<Parent>();
		final ArrayList<Locus> scored=new ArrayList<Locus>();
	}
}
