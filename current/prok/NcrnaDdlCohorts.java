package prok;

import java.util.Arrays;
import java.util.Random;

import map.IntHashSet;
import parse.LineParser1;
import prok.NcrnaDdlIdentity.Locus;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.Input;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntHashMap;
import structures.IntList;

/** Frozen identity bins and group cross-fit roles, before new sketch ranks.
 * Split all resolved groups with a selected NOT_SEEN target, including targets
 * with no native rank or finite identity. Quantiles instead use finite original
 * NOT_SEEN native-common loci. No model scores enter either decision.
 * @author Brian, Raiden
 */
public final class NcrnaDdlCohorts{

	NcrnaDdlCohorts(final NcrnaDdlSavedData data, final Locus[] identities){
		require(identities!=null && identities.length==data.parents.size(),
			"Frozen identity bins must retain every original selected parent");
		identity=identities;
		bins=new int[identities.length];
		sides=new int[identities.length];
		final FloatList finite=new FloatList();
		final IntHashSet seen=new IntHashSet(256);
		final IntList roster=new IntList();
		for(final Parent parent:data.parents){
			if(parent.stratum==0){
				if(seen.add(parent.input.group)){roster.add(parent.input.group);}
				if(parent.nativeRank>0 && valid(identity[parent.index].identity())){finite.add(identity[parent.index].identity());}
			}
		}
		cuts=quantiles(finite.toArray());
		finiteNotSeen=finite.size;
		roster.sort();
		groups=roster.toArray();
		final int[] assignments=split(groups);
		for(int i=0; i<groups.length; i++){byGroup.put(groups[i], assignments[i]);}
		for(final Parent parent:data.parents){
			bins[parent.index]=bin(identity[parent.index].identity(), cuts);
			sides[parent.index]=byGroup.get(parent.input.group);
		}
	}

	/** Restores the accepted identity/cohort files without measuring identities.
	 * Recomputing metadata from saved float values checks the frozen bins/group
	 * assignment; it does not refit identities or consult any new sketch score. */
	public static NcrnaDdlCohorts read(final String prefix, final NcrnaDdlSavedData data,
			final byte[][] names, final byte[][] references){
		require(prefix!=null && data!=null, "Restored cohorts require an explicit accepted prefix and complete saved parent roster");
		final int[] bins=new int[data.parents.size()], sides=new int[bins.length];
		final Locus[] identities=NcrnaDdlIdentity.read(prefix+".identity.tsv", data, names, references, bins, sides);
		final NcrnaDdlCohorts result=new NcrnaDdlCohorts(data, identities);
		require(Arrays.equals(bins, result.bins) && Arrays.equals(sides, result.sides),
			"Saved identity bins/group sides differ from the original float population and declared deterministic group split");
		result.checkGroups(prefix+".groups.tsv");
		result.checkDistribution(prefix+".identity_bins.tsv", data);
		return result;
	}

	private void checkGroups(final String path){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, GROUP_HEADER, path);
			for(final int group:groups){
				final byte[] row=in.reader.nextLine();
				check(row!=null, "Saved group roster is truncated", path);
				p.set(row);
				check(p.terms()==5 && p.parseInt(0)==group && p.termEquals(byGroup.get(group)==0 ? "A" : "B", 1)
					&& p.parseLong(2)==20261113 && p.termEquals("ALL_SELECTED_NOT_SEEN_TARGET_GROUPS", 3)
					&& p.termEquals("FLOOR_HALF_IN_A", 4),
					"Saved group roster/side/seed differs from the frozen identity cohort", path);
			}
			check(in.reader.nextLine()==null, "Saved group roster appends an extra or duplicate assembly group", path);
		}
	}

	private void checkDistribution(final String path, final NcrnaDdlSavedData data){
		final LineParser1 p=new LineParser1('\t');
		final int[] counts=new int[9];
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, DISTRIBUTION_HEADER, path);
			for(int state=0; state<3; state++){
				for(int bin=0; bin<5; bin++){
					final byte[] row=in.reader.nextLine();
					check(row!=null, "Saved identity distribution is truncated", path);
					p.set(row);
					check(p.terms()==17 && p.termEquals(NcrnaDdlSavedData.STRATA[state], 0) && p.termEquals(BINS[bin], 1),
						"Saved identity distribution must retain every stratum/bin in the original order", path);
					Arrays.fill(counts, 0);
					for(final Parent parent:data.parents){
						if(parent.stratum!=state || bins[parent.index]!=bin){continue;}
						counts[0]++;
						if(parent.nativeRank>0){counts[1]++;}
						final Locus value=identity[parent.index];
						if(valid(value.identity())){counts[2]++;}
						else{
							counts[3]++;
							for(int reason=0; reason<MISSING.length; reason++){
								if(value.status().equals(MISSING[reason])){counts[4+reason]++;}
							}
						}
					}
					for(int field=0; field<counts.length; field++){
						check(p.parseInt(field+2)==counts[field], "Saved identity distribution counts differ from its complete parent rows", path);
					}
					check(p.parseInt(11)==finiteNotSeen && p.termEquals(valid(cuts[0]) ? "true" : "false", 12)
						&& p.parseInt(16)==Float.floatToRawIntBits(cuts[0]),
						"Saved identity cutoff population/availability/float bits differ from the frozen values", path);
					for(int cut=0; cut<cuts.length; cut++){
						if(!valid(cuts[cut])){check(p.termEquals("NA", 13+cut), "Missing saved cutoff must remain NA", path);}
						else{
							final double displayed=p.parseDouble(13+cut);
							check(!Double.isNaN(displayed) && !Double.isInfinite(displayed)
								&& Math.abs(displayed-cuts[cut])<=0.000000000500001,
								"Saved identity cutoff decimal differs from its frozen float value", path);
						}
					}
				}
			}
			check(in.reader.nextLine()==null, "Saved identity distribution appends an extra stratum/bin", path);
		}
	}

	/** Nearest-rank Q1/Q2/Q3; finite identity population and ties are explicit. */
	public static float[] quantiles(final float[] values){
		require(values!=null, "An identity distribution must explicitly supply its finite population");
		final float[] sorted=values.clone();
		for(final float value:sorted){require(valid(value), "Only finite native identities in [0,1] enter the declared quartiles");}
		Arrays.sort(sorted);
		final float[] result={Float.NaN, Float.NaN, Float.NaN};
		for(int q=1; q<=3 && sorted.length>0; q++){result[q-1]=sorted[(int)(((long)sorted.length*q+3)/4)-1];}
		return result;
	}

	/** Returns sides in sorted-group order; seed and rounding are owner-declared. */
	public static int[] split(final int[] sortedGroups){
		require(sortedGroups!=null, "Group splitting requires the explicit sorted original group roster");
		final int[] order=new int[sortedGroups.length], sides=new int[sortedGroups.length];
		for(int i=0; i<order.length; i++){
			require(sortedGroups[i]>=0 && (i==0 || sortedGroups[i]>sortedGroups[i-1]),
				"Each resolved assembly group must occur exactly once in ascending order before shuffling");
			order[i]=i;
		}
		final Random random=new Random(20261113);
		for(int i=order.length-1; i>0; i--){
			final int j=random.nextInt(i+1), old=order[i];
			order[i]=order[j];
			order[j]=old;
		}
		for(int i=0; i<order.length; i++){sides[order[i]]=i<order.length/2 ? 0 : 1;}
		return sides;
	}

	public static int bin(final float identity, final float[] cuts){
		require(cuts!=null && cuts.length==3, "Identity bins use the three frozen nearest-rank quartiles");
		if(!valid(identity) || !valid(cuts[0])){return 4;}
		return identity<=cuts[0] ? 0 : identity<=cuts[1] ? 1 : identity<=cuts[2] ? 2 : 3;
	}

	/** direction0 fits A and judges B; direction1 fits B and judges A.
	 * Even SEEN loci from a judging group are excluded from fitting. */
	public static int role(final int stratum, final int side, final int direction){
		require(stratum>=0 && stratum<3 && side>=-1 && side<=1 && direction>=-1 && direction<=1,
			"Cross-fit roles require an original seen label, optional group side and declared direction");
		if(direction<0){return ALL;}
		if(stratum==2){return UNKNOWN;}
		final int judge=1-direction;
		if(stratum==0){
			require(side>=0, "Every NOT_SEEN target group must belong to exactly one frozen half");
			return side==judge ? JUDGE : TRAIN;
		}
		return side==judge ? EXCLUDED_SEEN : AUXILIARY_SEEN;
	}

	int role(final Parent parent, final int direction){return role(parent.stratum, sides[parent.index], direction);}
	boolean fit(final Parent parent, final int direction){
		final int role=role(parent, direction);
		return parent.nativeRank>0 && bins[parent.index]==0 && (role==TRAIN || role==AUXILIARY_SEEN);
	}

	void write(final String prefix, final NcrnaDdlSavedData data, final byte[][] names){
		final ByteBuilder b=new ByteBuilder();
		try(Output loci=new Output(prefix+".identity.tsv"); Output distribution=new Output(prefix+".identity_bins.tsv");
				Output split=new Output(prefix+".groups.tsv")){
			loci.writer.println(IDENTITY_HEADER);
			for(final Parent parent:data.parents){
				final Locus value=identity[parent.index];
				b.clear().append(parent.input.panel.nameBytes).tab().append(parent.input.slot).tab().append(parent.input.idBytes)
					.tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum]).tab().append(parent.input.group)
					.tab().append(parent.nativeRank>0).tab().append(value.model()).tab();
				if(value.model()<0){b.append('.');}else{b.append(names[value.model()]);}
				b.tab().append(value.modelEvent()).tab().append(value.geneEvent()).tab().append(value.containingWindows()).tab().append(value.definedWindows()).tab();
				if(valid(value.identity())){b.append(value.identity(), 9);}else{b.append("NA");}
				loci.writer.print(b.tab().append(Float.floatToRawIntBits(value.identity())).tab().append(BINS[bins[parent.index]])
					.tab().append(sides[parent.index]<0 ? "OUTSIDE_NOT_SEEN_ROSTER" : sides[parent.index]==0 ? "A" : "B").tab().append(value.status()).nl());
			}
			distribution.writer.println(DISTRIBUTION_HEADER);
			final int[] missing=new int[MISSING.length];
			for(int state=0; state<3; state++){
				for(int bin=0; bin<5; bin++){
					int targets=0, common=0, finite=0;
					Arrays.fill(missing, 0);
					for(final Parent parent:data.parents){
						if(parent.stratum!=state || bins[parent.index]!=bin){continue;}
						targets++;
						if(parent.nativeRank>0){common++;}
						if(valid(identity[parent.index].identity())){finite++;}
						else{
							boolean known=false;
							for(int reason=0; reason<MISSING.length; reason++){
								if(identity[parent.index].status().equals(MISSING[reason])){
									missing[reason]++;
									known=true;
								}
							}
							require(known, "Missing identities require one of the native helper's explicit dispositions");
						}
					}
					b.clear().append(NcrnaDdlSavedData.STRATA[state]).tab().append(BINS[bin]).tab().append(targets).tab().append(common)
						.tab().append(finite).tab().append(targets-finite);
					for(final int count:missing){b.tab().append(count);}
					b.tab().append(finiteNotSeen).tab().append(valid(cuts[0]));
					for(final float cut:cuts){
						b.tab();
						if(valid(cut)){b.append(cut, 9);}else{b.append("NA");}
					}
					distribution.writer.print(b.tab().append(Float.floatToRawIntBits(cuts[0])).nl());
				}
			}
			split.writer.println(GROUP_HEADER);
			for(final int group:groups){split.writer.print(b.clear().append(group).tab().append(byGroup.get(group)==0 ? "A" : "B")
				.append("\t20261113\tALL_SELECTED_NOT_SEEN_TARGET_GROUPS\tFLOOR_HALF_IN_A\n"));}
		}
	}

	private static boolean valid(final float value){return !Float.isNaN(value) && !Float.isInfinite(value) && value>=0 && value<=1;}
	private static void check(final boolean valid, final String reason, final String path){
		if(!valid){throw new IllegalArgumentException(reason+": "+path);}
	}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	final Locus[] identity;
	final int[] bins, sides, groups;
	final float[] cuts;
	final int finiteNotSeen;
	private final IntHashMap byGroup=new IntHashMap();
	static final String IDENTITY_HEADER="panel\tinput_slot\tinput_id\tparent\tseen_stratum\tassembly_group\tnative_common\tmodel_index\tmodel\tmodel_event\tgene_event\tcontaining_windows\tdefined_windows\tidentity\tidentity_float_bits\tidentity_bin\tgroup_side\tstatus";
	private static final String DISTRIBUTION_HEADER="seen_stratum\tidentity_bin\ttarget_loci\tnative_common_loci\tfinite_identities\tmissing_identities\tnot_native_rankable\tno_complete_window\tundefined_gene\tundefined_consensus\tinvalid_alignment\tcut_population\tcut_available\tq1\tq2\tq3\tq1_float_bits";
	private static final String GROUP_HEADER="assembly_group\tside\tseed\troster\trounding";
	public static final int ALL=0, JUDGE=1, TRAIN=2, EXCLUDED_SEEN=3, UNKNOWN=4, AUXILIARY_SEEN=5;
	public static final String[] ROLES={"ALL", "JUDGING", "TRAIN_SIDE", "EXCLUDED_SEEN_JUDGE_GROUP", "UNKNOWN", "SEEN_AUXILIARY"};
	public static final String[] BINS={"Q1_LOW", "Q2", "Q3", "Q4", "NA", "ALL"};
	private static final String[] MISSING={NcrnaDdlIdentity.NOT_NATIVE_RANKABLE, NcrnaDdlIdentity.NO_COMPLETE_WINDOW,
		NcrnaDdlIdentity.UNDEFINED_GENE, NcrnaDdlIdentity.UNDEFINED_CONSENSUS, NcrnaDdlIdentity.INVALID_ALIGNMENT};
}
