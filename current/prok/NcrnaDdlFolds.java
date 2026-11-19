package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Selected;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntHashMap;
import structures.IntList;

/** Assembly-group folds and exact-duplicate exclusions of checked member voters.
 * Borrows saved sequences and the existing exact-sequence set; never aligns or
 * changes consensus sequences, accepted-model bits, or original model ordinals.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlFolds implements NcrnaDdlMembers.KeptMemberSink{

	NcrnaDdlFolds(final NcrnaDdlSavedData data_, final int models, final String groups, final String genomes){
		data=data_;
		memberGroups=new IntList[models];
		for(int model=0; model<models; model++){memberGroups[model]=new IntList();}
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(groups)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GROUP_HEADER, groups);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==14 && p.parseInt(10)>=0, "Fold aliases require the accepted resolved assembly-group table");
				final int group=p.parseInt(10);
				for(int column=2; column<=3; column++){
					final String path=text(p, column);
					require(!pathGroups.contains(path) || pathGroups.get(path)==group, "A literal/canonical path cannot change assembly group");
					pathGroups.put(path, group);
				}
			}
		}
		try(Input in=new Input(genomes)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GENOME_HEADER, genomes);
			int ordinal=0;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==11 && p.parseInt(0)==ordinal++, "Fold donor binding retains the complete ordered genome roster");
				final String source=text(p, 1), path=text(p, 7);
				require(!sourceGroups.contains(source), "Each donor source must have one original genome binding");
				sourceGroups.put(source, pathGroups.get(path));
			}
		}
		for(final Selected input:data.selected){
			require(pathGroups.get(input.path)==input.group,
				"Saved query and donor assembly-group namespaces must agree before assigning folds: "+input.path);
		}
	}

	@Override
	public void accept(final String header, final int model, final byte[] bases){
		assert(model>=0 && model<memberGroups.length && bases!=null) : "Only checked retained members reach the fold sink";
		final String source=RrnaResourceIO.field(header, "source");
		require(sourceGroups.contains(source) && sourceGroups.get(source)>=0,
			"Every voting member needs a resolved donor assembly before fold exclusion: "+source);
		final int group=sourceGroups.get(source);
		memberGroups[model].add(group);
		contributingGroups.put(group, 1);
	}

	/** Original contribution precedes, and is independent of, judging exclusions. */
	boolean originalContributor(final Parent parent){return contributingGroups.contains(parent.input.group);}

	/** Extract once from saved windows; read a source genome only for missing genes. */
	void prepare(final NcrnaDdlMembers members, final NcrnaDdlCohorts cohorts){
		genes=new byte[data.parents.size()][];
		low=new boolean[genes.length];
		for(int model=0; model<memberGroups.length; model++){
			require(memberGroups[model].size==members.members(model).length,
				"Member group order must exactly follow the checked per-model voter arrays");
		}
		for(final Window window:data.windows){
			for(final Edge edge:window.edges){
				final Parent parent=edge.parent;
				if(!comparable(parent)){continue;}
				final byte[] gene=NcrnaDdlLocusOverlap.extract(window.bases, window.start, window.stop,
					window.strand, parent.start, parent.stop, parent.strand);
				retain(parent, gene);
			}
		}
		final FloatList finite=new FloatList();
		for(final Selected input:data.selected){
			boolean missing=false;
			for(final Parent parent:data.parents){missing|=parent.input==input && comparable(parent) && genes[parent.index]==null;}
			if(!missing){continue;}
			NcrnaDdlSavedData.sequences(input.path, read->{
				final byte[] seqid=RrnaResourceIO.first(read.id).getBytes(StandardCharsets.UTF_8);
				for(final Parent parent:data.parents){
					if(parent.input!=input || !comparable(parent) || genes[parent.index]!=null || !Arrays.equals(seqid, parent.seqidBytes)){continue;}
					retain(parent, NcrnaDdlLocusOverlap.extract(read.bases, 0, read.length()-1, 0,
						parent.start, parent.stop, parent.strand));
				}
			});
		}
		for(final Parent parent:data.parents){
			if(!comparable(parent)){continue;}
			require(genes[parent.index]!=null, "A comparable held-out locus needs complete defined DNA for exact-duplicate exclusion; parent="+parent.index);
			final float identity=cohorts.identity[parent.index].identity();
			if(Float.isFinite(identity)){finite.add(identity);}
		}
		final float[] cuts=NcrnaDdlCohorts.quantiles(finite.toArray());
		finiteIdentities=finite.size;
		require(finiteIdentities>0, "The all-comparable low-identity quartile needs at least one finite frozen identity");
		lowCut=cuts[0];
		for(final Parent parent:data.parents){low[parent.index]=NcrnaDdlCohorts.bin(cohorts.identity[parent.index].identity(), cuts)==0;}
	}

	private void retain(final Parent parent, final byte[] gene){
		assert(genes!=null && parent.index<genes.length) : "Gene witnesses bind the original complete parent roster";
		if(gene==null){return;}
		require(genes[parent.index]==null || Arrays.equals(genes[parent.index], gene), "Complete saved gene witnesses must agree across windows");
		if(genes[parent.index]==null){genes[parent.index]=gene;}
	}

	/** Balanced judging-group assignment, with background-only groups assigned separately. */
	int[] assign(final long seed, final int count, final Output out){
		require(genes!=null && seed>=0 && count>=2 && count<=5, "Prepare complete gene witnesses and declare two through five folds before splitting");
		foldCount=count;
		final IntHashMap target=new IntHashMap();
		for(final Parent parent:data.parents){if(comparable(parent)){target.put(parent.input.group, 1);}}
		final IntList groups=new IntList(), background=new IntList();
		for(final Selected input:data.selected){(target.contains(input.group) ? groups : background).add(input.group);}
		groups.sort();
		background.sort();
		final IntHashMap byGroup=new IntHashMap();
		shuffle(groups.array, groups.size, new Random(seed));
		shuffle(background.array, background.size, new Random(seed^0x5deece66dL));
		for(int i=0; i<groups.size; i++){byGroup.put(groups.get(i), i%count);}
		for(int i=0; i<background.size; i++){byGroup.put(background.get(i), i%count);}
		final int[] folds=new int[data.selected.size()];
		final ByteBuilder b=new ByteBuilder();
		for(final Selected input:data.selected){
			folds[input.index]=byGroup.get(input.group);
			out.writer.print(b.clear().append(input.index).tab().append(input.group).tab().append(folds[input.index])
				.tab().append(seed).tab().append(target.contains(input.group)).nl());
		}
		return folds;
	}

	/** Restores the fixed not-seen halves and their separate background split. */
	int[] fixed(final NcrnaDdlCohorts cohorts, final Output out){
		final IntHashMap byGroup=new IntHashMap();
		for(final Parent parent:data.parents){if(parent.stratum==0){byGroup.put(parent.input.group, cohorts.sides[parent.index]);}}
		final IntList extra=new IntList();
		for(final Selected input:data.selected){if(!byGroup.contains(input.group)){extra.add(input.group);}}
		extra.sort();
		final int[] sides=NcrnaDdlCohorts.split(extra.toArray());
		for(int i=0; i<sides.length; i++){byGroup.put(extra.get(i), sides[i]);}
		final int[] result=new int[data.selected.size()];
		final ByteBuilder b=new ByteBuilder();
		for(final Selected input:data.selected){
			result[input.index]=byGroup.get(input.group);
			out.writer.print(b.clear().append(input.index).tab().append(input.group).tab().append(result[input.index])
				.append("\t20261113\t").append(Arrays.binarySearch(extra.array, 0, extra.size, input.group)<0).nl());
		}
		return result;
	}

	/** Bind an already recorded assignment before reusing the exclusion routine. */
	void useSavedAssignment(final int[] assignment, final int count){
		require(genes!=null && assignment!=null && assignment.length==data.selected.size() && count>=2 && count<=5,
			"Saved fold reuse requires the prepared original gene witnesses and one recorded fold per selected input");
		for(final int fold:assignment){require(fold>=0 && fold<count, "A saved assignment cannot introduce an undeclared judging fold");}
		foldCount=count;
	}

	/** Remove entire judging assemblies plus exact sense-oriented duplicates anywhere. */
	boolean[][] exclusions(final int[] folds, final int judge, final NcrnaDdlMembers members){
		require(folds!=null && folds.length==data.selected.size() && judge>=0 && judge<foldCount,
			"Exclusions require one complete query-group assignment and an explicit judging fold");
		final IntHashMap heldGroups=new IntHashMap();
		for(final Selected input:data.selected){if(folds[input.index]==judge){heldGroups.put(input.group, 1);}}
		final NcrnaExactSequenceSet exact=new NcrnaExactSequenceSet();
		for(final Parent parent:data.parents){
			if(comparable(parent) && folds[parent.input.index]==judge){exact.add(genes[parent.index]);}
		}
		require(exact.size()>0, "Each judging fold needs at least one comparable complete locus");
		final boolean[][] excluded=new boolean[memberGroups.length][];
		for(int model=0; model<excluded.length; model++){
			final byte[][] sequences=members.members(model);
			excluded[model]=exclude(memberGroups[model].toArray(), sequences, heldGroups, exact);
		}
		return excluded;
	}

	/** Membership policy exposed for a literal cross-assembly exclusion fixture. */
	public static boolean[] exclude(final int[] groups, final byte[][] members,
			final IntHashMap heldGroups, final NcrnaExactSequenceSet exact){
		require(groups!=null && members!=null && groups.length==members.length && heldGroups!=null && exact!=null,
			"Every kept voter needs its aligned original assembly group and the complete held-out sequence set");
		final boolean[] excluded=new boolean[members.length];
		for(int i=0; i<members.length; i++){excluded[i]=heldGroups.contains(groups[i]) || exact.contains(members[i]);}
		return excluded;
	}

	static boolean comparable(final Parent parent){return parent.nativeRank>0 && parent.stratum<2;}

	private static void shuffle(final int[] values, final int length, final Random random){
		assert(values!=null && length<=values.length) : "Only the populated native-list prefix belongs to the fold shuffle";
		for(int i=length-1; i>0; i--){
			final int other=random.nextInt(i+1), value=values[i];
			values[i]=values[other];
			values[other]=value;
		}
	}
	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 cannot change fold source identity");
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private final NcrnaDdlSavedData data;
	private final IntList[] memberGroups;
	private final IntHashMap contributingGroups=new IntHashMap();
	private final ObjectIntMap<String> pathGroups=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> sourceGroups=new ObjectIntMap<String>(String.class);
	private byte[][] genes;
	boolean[] low;
	float lowCut;
	int finiteIdentities;
	private int foldCount;
}
