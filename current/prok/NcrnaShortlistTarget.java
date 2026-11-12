package prok;

import java.util.List;

import prok.NcrnaFamilyRegistry.Entry;
import prok.NcrnaReferenceCalls.Rule;

/** Resolves one generic caller index and its separate reference-grading scope.
 * Aggregate resource names never select a member implicitly. Group members may
 * be varied independently while all members remain enabled and union-graded.
 * This does not construct subtype truth or claim subtype recovery/exclusivity.
 * Instances borrow immutable registry metadata and are safe to share with workers.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaShortlistTarget{

	public NcrnaShortlistTarget(final String name, final String view){
		this(name, view, false);
	}

	/** Explicit joint mode selects every group member; legacy mode remains member-only. */
	public NcrnaShortlistTarget(final String name, final String view, final boolean joint_){
		joint=joint_;
		require(name!=null && !name.isEmpty(), "A shortlist target must name one registered caller");
		require(joint || (!name.equalsIgnoreCase("cobalamin") && !name.equalsIgnoreCase("6S")
			&& !name.equalsIgnoreCase("sixs") && !name.equalsIgnoreCase("SRP")),
			"Aggregate targets have independent indices; select an explicit cobalamin_rf, sixs_rf or srp member: "+name);
		final String resolved=joint && (name.equalsIgnoreCase("6S") || name.equalsIgnoreCase("sixs"))
			? "sixs_rf00013" : joint && name.equalsIgnoreCase("SRP") ? "srp_small" : name;
		entry=NcrnaFamilyRegistry.get(resolved);
		referenceFamily=referenceFamily(entry);
		require(entry.outputType==ProkObject.RNA && entry.boundary!=null && !entry.boundary.trna,
			"This native reference projection grades generic RNA; independent tRNA and specialized rRNA need their own integration: "+name);
		require("target".equals(view) || "union".equals(view), "grading must be target or union: "+view);
		union="union".equals(view);
		require(union==(entry.group!=null),
			"A grouped caller member requires explicit grading=union; a standalone caller requires grading=target: "+entry.id);
	}

	/** Find the actual enabled caller, then check every included model projection.
	 * Target-only maps must not accidentally grade another caller. A group union
	 * includes all of its loaded members, never just the member being varied.
	 */
	public int bind(final List<NcrnaFamily> families, final Rule[] rules, final String scoringFamily){
		require(families!=null && rules!=null && scoringFamily!=null, "Target binding needs loaded callers and an explicit reference family");
		require(referenceFamily.equals(scoringFamily),
			"Selected caller "+entry.id+" requires reference family "+referenceFamily+", not "+scoringFamily);
		final NcrnaReferenceCalls projection=new NcrnaReferenceCalls(rules);
		final boolean[] used=new boolean[rules.length];
		int selected=-1;
		int members=0;
		for(int i=0; i<families.size(); i++){
			final NcrnaFamily family=families.get(i);
			final Entry registration=NcrnaFamilyRegistry.get(family.name);
			final boolean include=union ? registration.group==entry.group : registration==entry;
			if(registration==entry){
				require(selected<0, "The target must be enabled exactly once: "+entry.id);
				selected=i;
			}
			if(family.outputType!=ProkObject.RNA){continue;}
			if(include){members++;}
			for(final String model:family.modelNames){
				final int index=projection.projectModel(model);
				final Rule rule=rules[index];
				require(rule.include==include && (!include || rule.family.equals(scoringFamily)),
					"Caller target and grading scope disagree for "+family.name+" model "+model);
				require(rule.displayLabel.equals(NcrnaFamilyRegistry.displayName(family.name)),
					"Model projection must retain the caller's display family: "+family.name);
				if(include){used[index]=true;}
			}
		}
		require(selected>=0, "The native caller configuration must enable target "+entry.id);
		require(members==(union ? entry.group.members.length : 1), "Union grading requires every registered member to remain enabled: "+entry.id);
		for(int i=0; i<rules.length; i++){
			require(!rules[i].include || used[i], "Included grading rule has no model in the selected caller scope: "+rules[i].modelPrefix);
		}
		return selected;
	}

	/** Ordered loaded indices for all varied members, after ordinary projection validation. */
	public int[] bindTargets(final List<NcrnaFamily> families, final Rule[] rules, final String scoringFamily){
		final int selected=bind(families, rules, scoringFamily);
		if(!joint || !union){return new int[]{selected};}
		final int[] indices=new int[entry.group.members.length];
		int count=0;
		for(int i=0; i<families.size(); i++){
			if(NcrnaFamilyRegistry.get(families.get(i).name).group==entry.group){indices[count++]=i;}
		}
		assert(count==indices.length) : "bind verified the complete union before assigning joint policy indices";
		return indices;
	}

	/** Freeze each family's packaged mapping mode, fractions and length/window
	 * settings. The sweep changes only cap and the matching absolute hit cutoff.
	 * Resource/model identity is supplied by the run's ordinary packaged loader.
	 */
	public void validate(final NcrnaFamily family, final boolean adaptiveColumn, final int baselineCap, final float baselineFloor, final int baselineMinHits){
		require(family!=null && entry.id.equals(family.name), "Policy validation must bind the selected loaded family");
		final boolean adaptive=entry.boundary.adaptive;
		final int expectedCap=entry.topN<0 || entry.dynamicModels && entry.models.length>0 ? family.library.length : entry.topN;
		require(family.indexK==entry.indexK && family.adaptive==adaptive && adaptiveColumn==adaptive
			&& family.adaptFloor==entry.boundary.floor && family.adaptTopFrac==entry.boundary.topFrac
			&& family.adaptQFrac==entry.boundary.qFrac && family.fixedMinHits==entry.minHits
			&& family.modelPatience==entry.modelPatience,
			"The target must retain its packaged mapping mode, verifier patience and constants; use floor for adaptive or minhits for fixed: "+entry.id);
		require(family.minLen==entry.minLen && family.windowPad==entry.windowPad && family.kLong==entry.seedK,
			"The shortlist sweep must not change packaged gene length, window padding or discovery seed K: "+entry.id);
		require(family.indexTopN==expectedCap && baselineCap==expectedCap
			&& (adaptive ? baselineFloor==family.adaptFloor : baselineMinHits==family.fixedMinHits),
			"First policy must equal the selected family's packaged caller baseline: "+entry.id);
		final boolean bounds=family.modelThresholds!=null && family.modelThresholds.hasLengths();
		final boolean packagedBounds=new FamilyResourceResolver().optional(entry.id, "length_bounds", null)!=null;
		require(bounds==packagedBounds, "The sweep must retain the packaged model-length gate: "+entry.id);
	}

	/** Canonical family identifiers in the native adjudicated-reference schema.
	 * An unknown binding needs an explicit schema decision, not a label inferred
	 * from the caller's own projection map. Specialized/eukaryotic truth is outside
	 * this prokaryotic reference contract until its canonical IDs are registered.
	 */
	private static String referenceFamily(final Entry registration){
		if(registration.group!=null){return registration.group.id;}
		switch(registration.id){
			case "rnasep": return "RNaseP";
			case "tmrna": return "tmRNA";
			case "fmn": return "FMN";
			case "sam1": return "SAM_I";
			case "tpp": return "TPP";
			case "glycine": return "glycine";
			default: throw new IllegalArgumentException("No native reference-family binding for shortlist target: "+registration.id);
		}
	}

	/** Parses a literal nonnegative integer without floating-point rounding.
	 * The parser borrows a byte range and allocates nothing on valid input.
	 */
	public static int fixedHits(final byte[] line, final int from, final int to){
		require(line!=null && from>=0 && from<to && to<=line.length, "Fixed minhits requires a nonempty byte field");
		int value=0;
		for(int i=from; i<to; i++){
			final int digit=line[i]-'0';
			require(digit>=0 && digit<=9 && value<=(Integer.MAX_VALUE-digit)/10,
				"Fixed minhits must contain decimal integer digits in [0,2147483647]; fractional or overflowing text is invalid");
			value=value*10+digit;
		}
		return value;
	}

	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	public final Entry entry;
	public final String referenceFamily;
	public final boolean union;
	public final boolean joint;
}
