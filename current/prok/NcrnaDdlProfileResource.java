package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import ddl.DDLIndexSmall;
import cardinality.DynamicDemiLog;
import map.LongHashSet;
import kmer.AbstractKmerTable;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Model-bound member-voted keys for cobalamin caller experiments.
 * V1 fixes forward middle-masked k10/2048 for native-first hybrid ordering.
 * V2 fixes partitioned plain9+masked10/4096 and long+.5*short, sketch-only.
 * V3 fixes plain7/masked10/shared-consensus quarters at512/1024, sketch-only.
 * V4 fixes forward8/masked10 native4+12 quarter-short keys, real member
 * plurality plus100perdistinctconsensuskey, equal-bin queries and saved padding.
 * Both use seed12345, exponent4 and raw-query-key bucket matches. Every row retains its
 * original ordinal, first-token ID, uppercase consensus bytes and voter count.
 * Loading compares the actual model identities and sequences, without hashes.
 * Ordering-only minor resources declare their member explicitly; the old
 * RF00174 ordering header and historical profile formats remain readable.
 * The distinct ordering-bonus header records maximum-only or all-key pseudo-votes
 * without changing real voter counts or the member-only zero-voter contract.
 * Keys are private and read-only after construction; query scratch is per worker.
 * No caller enables or loads this resource by default.
 * @author Brian, Raiden, Ganyu
 */
public final class NcrnaDdlProfileResource{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaDdlProfileResource(final char[][] keys_, final int[] members_, final boolean dual_, final boolean compact_,
			final boolean selected_, final LongHashSet blacklist_, final boolean orderOnly_, final long seed_,
			final byte[][] consensuses, final boolean suppressed_, final int bonus_, final int bonusRule_){
		this(keys_, members_, dual_, compact_, selected_, blacklist_, orderOnly_, seed_,
			consensuses, suppressed_, bonus_, bonusRule_, false);
	}

	private NcrnaDdlProfileResource(final char[][] keys_, final int[] members_, final boolean dual_, final boolean compact_,
			final boolean selected_, final LongHashSet blacklist_, final boolean orderOnly_, final long seed_,
			final byte[][] consensuses, final boolean suppressed_, final int bonus_, final int bonusRule_,
			final boolean consensusReplaced_){
		assert(keys_!=null && members_!=null && keys_.length==members_.length && keys_.length>0)
			: "Resource rows must bind one profile and voter count per current model";
		orderOnly=orderOnly_;
		suppressed=suppressed_;
		validateOrderingBonus(bonus_, bonusRule_);
		require(bonus_==0 || orderOnly && !suppressed && consensuses==null,
			"An explicit in-profile ordering bonus cannot also request a separate consensus index or suppression");
		orderingBonus=bonus_;
		orderingBonusRule=bonusRule_;
		consensusReplaced=consensusReplaced_;
		require(!consensusReplaced || orderOnly && !suppressed && bonus_==0 && consensuses==null,
			"Consensus replacement uses a second ordinary index, never an additive index or stored pseudo-votes");
		require(!suppressed || orderOnly && consensuses!=null,
			"Suppressed member profiles must always serve with their complementary consensus index");
		seed=seed_;
		keys=orderOnly ? null : keys_;
		members=members_;
		dual=dual_;
		compact=compact_;
		selected=selected_;
		blacklist=blacklist_;
		nativeProfiles=orderOnly ? new DynamicDemiLog[keys_.length] : null;
		if(orderOnly){
			for(int model=0; model<keys_.length; model++){
				nativeProfiles[model]=DynamicDemiLog.fromArray(keys_[model], model, null, 10, -1);
			}
		}
		index=orderOnly ? new DDLIndexSmall(nativeProfiles) : new DDLIndexSmall(keys);
		consensusIndex=consensuses==null ? null : buildConsensusIndex(consensuses);
		gc=null;
	}

	/** Builds a separate frozen index once; member keys and global DDL settings stay untouched. */
	private DDLIndexSmall buildConsensusIndex(final byte[][] consensuses){
		assert(orderOnly && consensuses.length==members.length)
			: "Runtime consensus scores must bind the same validated model ordinals as the member index";
		final NcrnaDdlSketchBuilder builder=NcrnaDdlSketchBuilder.forKeyArrays(buckets(), seed, blacklist);
		final DynamicDemiLog[] profiles=new DynamicDemiLog[consensuses.length];
		for(int model=0; model<profiles.length; model++){
			final DynamicDemiLog profile=DynamicDemiLog.create(buckets(), 10, seed, 0, false, false);
			builder.fillKeys(consensuses[model], profile.maxArray());
			if(suppressed){
				for(int bin=0; bin<buckets(); bin++){
					final char consensusKey=profile.maxArray()[bin];
					require(consensusKey==0 || consensusKey!=key(model, bin),
						"Suppressed member and consensus keys must be disjoint in every bin: model="+model+", bin="+bin);
				}
			}
			profiles[model]=profile;
		}
		return new DDLIndexSmall(profiles);
	}

	/** Shares immutable keys/index; adding a cold GC policy never rebuilds sketches. */
	private NcrnaDdlProfileResource(final NcrnaDdlProfileResource source, final NcrnaWindowGcProfile gc_){
		require(source.selected && source.gc==null && gc_!=null && gc_.models()==source.models(),
			"Window-GC scoring requires an unmodified selected-format resource and its exact model roster");
		keys=source.keys;
		members=source.members;
		dual=source.dual;
		compact=source.compact;
		selected=source.selected;
		orderOnly=source.orderOnly;
		suppressed=source.suppressed;
		orderingBonus=source.orderingBonus;
		orderingBonusRule=source.orderingBonusRule;
		consensusReplaced=source.consensusReplaced;
		seed=source.seed;
		nativeProfiles=source.nativeProfiles;
		blacklist=source.blacklist;
		index=source.index;
		consensusIndex=source.consensusIndex;
		gc=gc_;
	}

	/** The model-GC loader binds names to current ordinals, not to row positions. */
	public NcrnaDdlProfileResource withGc(final NcrnaWindowGcProfile profile){
		return new NcrnaDdlProfileResource(this, profile);
	}

	/** Called only by the paired loader with its validated model sequence axis.
	 * A nonempty consensus bin replaces its member key, even when the new key is
	 * smaller. Empty consensus bins retain the member key. No file or vote count
	 * is changed, and neither input is borrowed as mutable construction scratch. */
	private NcrnaDdlProfileResource withConsensusReplacement(final byte[][] consensuses){
		require(orderOnly && !suppressed && orderingBonus==0 && consensusIndex==null
			&& !consensusReplaced && gc==null && consensuses!=null && consensuses.length==models(),
			"Consensus replacement requires the unmodified member-only profiles and their bound model sequence axis");
		assert(nativeProfiles!=null && nativeProfiles.length==members.length)
			: "The ordering-only constructor binds one standard DDL per member-profile row before replacement";
		final int bins=buckets();
		final char[][] replaced=new char[models()][bins];
		final char[] consensusKeys=new char[bins];
		final NcrnaDdlSketchBuilder builder=NcrnaDdlSketchBuilder.forKeyArrays(bins, seed, blacklist);
		for(int model=0; model<models(); model++){
			require(consensuses[model]!=null && consensuses[model].length>0,
				"The bound consensus must exist even for a zero-member profile: model="+model);
			System.arraycopy(nativeProfiles[model].maxArray(), 0, replaced[model], 0, bins);
			builder.fillKeys(consensuses[model], consensusKeys);
			for(int bin=0; bin<bins; bin++){
				if(consensusKeys[bin]!=0){replaced[model][bin]=consensusKeys[bin];}
			}
		}
		return new NcrnaDdlProfileResource(replaced, members, dual, compact, selected, blacklist,
			true, seed, null, false, 0, BONUS_NONE, true);
	}

	/** Reads exactly the expected ordered models; caller-owned arrays are not modified. */
	public static NcrnaDdlProfileResource load(final String path, final String[] names, final byte[][] sequences){
		return load(path, names, sequences, false);
	}

	/** Optional consensus scores use a separate load-time sketch, never member pseudo-votes. */
	public static NcrnaDdlProfileResource load(final String path, final String[] names,
			final byte[][] sequences, final boolean consensus){
		final Rows rows=readRows(path, names, sequences, consensus);
		require(!rows.suppressed, "Suppressed profiles require the shared clan loader, not separate or historical consensus scoring");
		return new NcrnaDdlProfileResource(rows.keys, rows.members, rows.dual, rows.compact,
			rows.selected, rows.blacklist, rows.orderOnly, rows.seed, consensus ? sequences : null, false, rows.bonus, rows.bonusRule);
	}

	/** Reads member rows without constructing unused per-member indexes. Global
	 * ordinals follow the supplied main/second/third roster and each local order.
	 * All members must use the same member-only or suppressed-dual format. */
	public static NcrnaDdlProfileResource loadClan(final String[] paths, final String[][] names,
			final byte[][][] sequences){
		require(paths!=null && names!=null && sequences!=null && paths.length==3
			&& names.length==3 && sequences.length==3, "Joint profiles require exactly the three ordered clan members");
		final String[] families={"RF00174", "RF01482", "RF01689"};
		final Rows[] rows=new Rows[3];
		int count=0;
		for(int member=0; member<3; member++){
			require(memberFamily(names[member]).equals(families[member]),
				"Joint profile ordinals must follow main, second and third member order");
			rows[member]=readRows(paths[member], names[member], sequences[member], false);
			final Rows row=rows[member];
			require(row.orderOnly && (row.keys[0].length==512 || row.keys[0].length==1024),
				"Joint ranking requires ordering-only 512/1024-bin profiles");
			if(member>0){
				final long[] expected=rows[0].blacklist.toArray(), observed=row.blacklist.toArray();
				Arrays.sort(expected);
				Arrays.sort(observed);
				require(row.keys[0].length==rows[0].keys[0].length && row.seed==rows[0].seed
					&& row.suppressed==rows[0].suppressed && row.bonus==rows[0].bonus && row.bonusRule==rows[0].bonusRule
					&& Arrays.equals(expected, observed),
					"One joint query requires equal widths, suppression/bonus rules, seeds and actual frozen padding words for every member");
			}
			count=Math.addExact(count, row.keys.length);
		}
		final char[][] keys=new char[count][];
		final int[] members=new int[count];
		final byte[][] consensus=rows[0].suppressed ? new byte[count][] : null;
		int next=0;
		for(int member=0; member<rows.length; member++){
			final Rows row=rows[member];
			System.arraycopy(row.keys, 0, keys, next, row.keys.length);
			System.arraycopy(row.members, 0, members, next, row.members.length);
			if(consensus!=null){System.arraycopy(sequences[member], 0, consensus, next, row.keys.length);}
			next+=row.keys.length;
		}
		assert(next==count) : "The joint index must conserve every validated local model ordinal";
		return new NcrnaDdlProfileResource(keys, members, true, false, false, rows[0].blacklist, true, rows[0].seed,
			consensus, rows[0].suppressed, rows[0].bonus, rows[0].bonusRule);
	}

	/** Loads each member file once, binding both indexes to the same validated
	 * ordered names and consensus bytes. No arbitrary replacement-axis API is
	 * exposed after the loaded resource has forgotten those identities. */
	public static ConsensusPair loadClanWithConsensusReplacement(final String[] paths, final String[][] names,
			final byte[][][] sequences){
		final NcrnaDdlProfileResource memberOnly=loadClan(paths, names, sequences);
		final byte[][] axis=new byte[memberOnly.models()][];
		int next=0;
		for(final byte[][] member:sequences){
			System.arraycopy(member, 0, axis, next, member.length);
			next+=member.length;
		}
		assert(next==axis.length) : "The clan loader concatenates main, second and third names and sequences in this same order";
		final long cpu=NcrnaCpuProfile.cpuNow();
		final long wall=System.nanoTime();
		final NcrnaDdlProfileResource replacement=memberOnly.withConsensusReplacement(axis);
		return new ConsensusPair(memberOnly, replacement, cpu<0 ? -1 : NcrnaCpuProfile.cpuNow()-cpu, System.nanoTime()-wall);
	}

	/** Two immutable query resources from one validated member-profile load. */
	public static final class ConsensusPair{
		/** The served small indexes store sorted posting runs rather than nullable
		 * per-key list objects. Presence is the same nonempty-list test, without
		 * allocating another index or a separate reference-key bitset.
		 */
		public boolean containsReferenceKey(final int bin, final int key){
			return memberOnly.index.hasPosting(bin, key) || replacement.index.hasPosting(bin, key);
		}
		private ConsensusPair(final NcrnaDdlProfileResource memberOnly_, final NcrnaDdlProfileResource replacement_,
				final long buildCpuNanos_, final long buildWallNanos_){
			assert(!memberOnly_.consensusReplaced && replacement_.consensusReplaced
				&& memberOnly_.models()==replacement_.models())
				: "Only the paired loader can publish both indexes on its one validated model axis";
			memberOnly=memberOnly_;
			replacement=replacement_;
			buildCpuNanos=buildCpuNanos_;
			buildWallNanos=buildWallNanos_;
		}
		public final NcrnaDdlProfileResource memberOnly, replacement;
		/** Replacement key/index construction only; member-profile I/O is excluded.
		 * CPU is unavailable (-1) unless the existing CPU profiler was enabled. */
		public final long buildCpuNanos, buildWallNanos;
	}

	/** Decoder shared by separate and joint loads; it owns no query or ranking index. */
	private static Rows readRows(final String path, final String[] names, final byte[][] sequences, final boolean consensus){
		validateModels(names, sequences);
		final String family=memberFamily(names);
		final int[] members=new int[names.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			final byte[] header=in.reader.nextLine();
			boolean orderOnly=false;
			boolean suppressed=false;
			int orderBins=0;
			int bonus=0, bonusRule=BONUS_NONE;
			long seed=12345;
			if(header!=null){
				p.set(header);
				suppressed=p.termEquals(SUPPRESSED_ORDER_PREFIX, 0);
				final boolean bonusHeader=p.termEquals(BONUS_ORDER_PREFIX, 0);
				final boolean memberHeader=p.termEquals(MEMBER_ORDER_PREFIX, 0) || suppressed || bonusHeader;
				orderOnly=p.termEquals(ORDER_PREFIX, 0) || memberHeader;
				if(orderOnly){
					final int shift=memberHeader ? 1 : 0;
					require(p.terms()==(bonusHeader ? 8 : 5+shift) && p.termEquals(ORDER_GEOMETRY, 3+shift)
						&& p.termEquals(suppressed ? SUPPRESSED_SCORING : bonusHeader ? BONUS_SCORING : ORDER_SCORING, 4+shift)
						&& (memberHeader ? p.termEquals(family, 1) : family.equals("RF00174")),
						"Ordering-only profiles require the declared member, native geometry and explicit vote rule");
					orderBins=p.parseInt(1+shift);
					seed=p.parseLong(2+shift);
					if(bonusHeader){
						bonus=p.parseInt(6);
						bonusRule=parseBonusRule(p.parseString(7));
						validateOrderingBonus(bonus, bonusRule);
						require(bonus>0, "A bonus header must declare positive consensus votes, not impersonate a member-only resource");
					}
					require(!suppressed || orderBins==1024, "Suppressed dual ordering is defined for 1024 bins");
					require((orderBins==512 || orderBins==1024) && seed>=0,
						"Ordering-only geometry needs power-of-two512/1024bins and an explicit nonnegative seed");
				}
			}
			final boolean selected512=Arrays.equals(header, SELECTED512_HEADER), selected1024=Arrays.equals(header, SELECTED1024_HEADER);
			final boolean selected=selected512 || selected1024;
			final boolean compact512=Arrays.equals(header, COMPACT512_HEADER);
			final boolean compact1024=Arrays.equals(header, COMPACT1024_HEADER);
			final boolean compact=compact512 || compact1024;
			final boolean dual=orderOnly || selected || compact || Arrays.equals(header, DUAL_HEADER);
			require(dual || Arrays.equals(header, HEADER), "Unsupported profile header: "+path);
			require(orderOnly || family.equals("RF00174"), "Historical profile formats bind only RF00174");
			final int buckets=orderOnly ? orderBins : selected512 || compact512 ? 512 : selected1024 || compact1024 ? 1024 : dual ? DUAL_BUCKETS : BUCKETS;
			require(!consensus || family.equals("RF00174") && orderOnly && buckets==1024 && bonus==0,
				"Separate consensus scoring requires ordering-only1024-bin member profiles, not a historical bonus format");
			final LongHashSet blacklist=new LongHashSet(4);
			if(selected || orderOnly){
				final byte[] metadata=in.reader.nextLine();
				require(metadata!=null, "Selected profile resource lacks padding-word metadata");
				p.set(metadata);
				require(p.terms()==2 && p.termEquals("#PADDING_BLACKLIST", 0), "Selected profiles must retain their query/model padding-word set");
				final int count=p.parseInt(1);
				require(count==0 || count==2304, "Selected padding words are explicitly empty or the complete256/2048prefixes");
				int shortWords=0;
				for(int i=0; i<count; i++){
					final byte[] row=in.reader.nextLine();
					require(row!=null, "Truncated selected padding-word metadata");
					p.set(row);
					require(p.terms()==2 && p.termEquals("#PADDING_WORD", 0) && blacklist.add(p.parseLong(1)),
						"Selected padding words must retain distinct tagged identities");
					if(shortPaddingWord(p.parseLong(1))){shortWords++;}
				}
				require(count==0 || shortWords==256 && count-shortWords==2048, "Selected resource must retain both complete padding-word prefixes");
			}
			final char[][] keys=new char[names.length][buckets];
			for(int model=0; model<names.length; model++){
				final byte[] row=in.reader.nextLine();
				require(row!=null, "Missing profile row for original model "+model);
				p.set(row);
				require(p.terms()==buckets+4 && p.parseInt(0)==model,
					"Profile rows must retain the original model order and complete bucket vector: model="+model);
				require(Arrays.equals(p.parseByteArray(1), id(names[model])),
					"Profile model ID differs from the actual caller library: model="+model);
				final byte[] sequence=p.parseByteArray(2);
				require(sequence.length==sequences[model].length,
					"Profile consensus length differs from the actual caller model: model="+model);
				for(int base=0; base<sequence.length; base++){
					if(sequence[base]!=upper(sequences[model][base])){
						throw new IllegalArgumentException("Profile consensus content differs from the actual caller model: model="+model+", base="+base);
					}
				}
				members[model]=p.parseInt(3);
				require(orderOnly ? members[model]>=0 : members[model]>0,
					"Each model profile requires actual kept member votes, with zero allowed only after explicit holdout: model="+model);
				for(int bucket=0; bucket<buckets; bucket++){
					final int key=p.parseInt(bucket+4);
					if(key<0 || key>Character.MAX_VALUE){
						throw new IllegalArgumentException("Absolute profile keys must fit unsigned16bits: model="+model+", bucket="+bucket);
					}
					keys[model][bucket]=(char)key;
					require(members[model]>0 || key==0 || bonus>0, "A zero-voter member-only profile cannot contain invented keys");
				}
				if(members[model]==0 && bonus>0){validateConsensusOnly(keys[model], sequences[model], buckets, seed, blacklist);}
			}
			require(in.reader.nextLine()==null, "Profile resource has rows beyond the actual caller model count");
			return new Rows(keys, members, dual, compact, selected, blacklist, orderOnly, seed, suppressed, bonus, bonusRule);
		}
	}

	/** Appends already-built ordering rows without voting or rebuilding old keys.
	 * Both complete resources are decoded before output creation. The original
	 * metadata and row bytes are copied; only new row ordinals are rewritten.
	 * Sequence/name arrays are borrowed read-only for this call. */
	public static void appendOrdering(final String firstPath, final String extraPath, final String output,
			final String[] names, final byte[][] sequences, final String[] extraNames, final byte[][] extraSequences){
		final Rows first=readRows(firstPath, names, sequences, false);
		final Rows extra=readRows(extraPath, extraNames, extraSequences, false);
		require(memberFamily(names).equals(memberFamily(extraNames)), "Appended profiles must belong to the same member");
		compatibleOrdering(first, extra);
		final int total=Math.addExact(names.length, extraNames.length);
		final String[] combinedNames=Arrays.copyOf(names, total);
		final byte[][] combinedSequences=Arrays.copyOf(sequences, total);
		System.arraycopy(extraNames, 0, combinedNames, names.length, extraNames.length);
		System.arraycopy(extraSequences, 0, combinedSequences, sequences.length, extraSequences.length);
		validateModels(combinedNames, combinedSequences);
		RrnaResourceIO.fresh(output, "");
		try(Output out=new Output(output); Input old=new Input(firstPath); Input added=new Input(extraPath)){
			final ByteBuilder row=new ByteBuilder(16384);
			for(byte[] line; (line=old.reader.nextLine())!=null;){out.writer.print(row.clear().append(line).nl());}
			// Header + padding declaration + the complete declared padding set.
			for(int i=0; i<2+extra.blacklist.size(); i++){
				require(added.reader.nextLine()!=null, "Appended profile metadata changed after validation");
			}
			final LineParser1 p=new LineParser1('\t');
			for(int i=0; i<extraNames.length; i++){
				final byte[] line=added.reader.nextLine();
				require(line!=null, "Appended profile rows changed after validation");
				p.set(line);
				require(p.terms()==extra.keys[i].length+4 && p.parseInt(0)==i, "Appended profile ordinal or row width changed after validation");
				int tab=0;
				while(tab<line.length && line[tab]!='\t'){tab++;}
				assert(tab<line.length) : "The validated complete profile row contains its ordinal delimiter";
				out.writer.print(row.clear().append(names.length+i).append(line, tab, line.length-tab).nl());
			}
			require(added.reader.nextLine()==null, "Appended profile acquired extra rows after validation");
		}
		final Rows actual=readRows(output, combinedNames, combinedSequences, false);
		compatibleOrdering(first, actual);
		for(int i=0; i<total; i++){
			final Rows source=i<names.length ? first : extra;
			final int index=i<names.length ? i : i-names.length;
			require(actual.members[i]==source.members[index] && Arrays.equals(actual.keys[i], source.keys[index]),
				"Appended profile changed a source voter count or key vector: model="+i);
		}
	}

	/** Ordering policies and query geometry must agree before combining axes. */
	private static void compatibleOrdering(final Rows first, final Rows second){
		require(first.orderOnly && second.orderOnly && first.keys[0].length==second.keys[0].length
			&& first.seed==second.seed && first.suppressed==second.suppressed
			&& first.bonus==second.bonus && first.bonusRule==second.bonusRule,
			"Profile append requires identical ordering geometry, seed and vote policy");
		final long[] before=first.blacklist.toArray(), after=second.blacklist.toArray();
		Arrays.sort(before);
		Arrays.sort(after);
		require(Arrays.equals(before, after), "Profile append requires the identical padding-word set");
	}

	/** Publishes one fresh resource after member voting, without changing its inputs.
	 * Caller precondition: profiles were built with K/BUCKETS,forward middle mask,
	 * seed12345/exponent4, as BuildNcrnaDdlProfiles does. Generic profile arrays
	 * carry no geometry provenance; serialization cannot infer it from their keys. */
	public static void write(final String path, final NcrnaDdlProfile[] profiles,
			final String[] names, final byte[][] sequences){
		write(path, profiles, names, sequences, false);
	}

	/** The dual version fixes PARTITION9+masked10/4096 and long+.5*short.
	 * Its producer must use buildDual with that exact geometry before serialization. */
	public static void write(final String path, final NcrnaDdlProfile[] profiles,
			final String[] names, final byte[][] sequences, final boolean dual){
		final int buckets=dual ? DUAL_BUCKETS : BUCKETS;
		write(path, profiles, names, sequences, buckets, dual ? DUAL_FORMAT : FORMAT);
	}

	/** Producer must use buildPartition(members, consensuses, 7, buckets).
	 * The header fixes both member partitions and their shared consensus quarter. */
	public static void writeCompact(final String path, final NcrnaDdlProfile[] profiles,
			final String[] names, final byte[][] sequences, final int buckets){
		require(buckets==512 || buckets==1024, "Selected compact caller formats require512or1024total bins");
		write(path, profiles, names, sequences, buckets, buckets==512 ? COMPACT512_FORMAT : COMPACT1024_FORMAT);
	}

	/** Selected8/masked10 quarters: one greatest real-member key per bin votes;
	 * every distinct eligible consensus bin/key receives100additional votes. */
	public static void writeSelected(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist){
		writeSelected(path, keys, members, names, sequences, bins, blacklist, false, 12345);
	}

	/** Member-only profiles retain their hash seed and permit real holdout-empty
	 * models. This format requests only native-eligible reordering in the caller. */
	public static void writeOrdering(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist, final long seed){
		writeSelected(path, keys, members, names, sequences, bins, blacklist, true, seed);
	}

	/** Explicit bonus metadata preserves true member counts, including zero. */
	public static void writeBonusOrdering(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist,
			final long seed, final int bonus, final int bonusRule){
		validateOrderingBonus(bonus, bonusRule);
		require(bonus>0, "Use the unchanged member-only writer when no consensus bonus is requested");
		writeSelected(path, keys, members, names, sequences, bins, blacklist, true, seed, false, bonus, bonusRule);
	}

	/** Stores only suppressed member keys. The runtime builds the complementary
	 * consensus index from these exact bound sequences and sketch metadata. */
	public static void writeSuppressedOrdering(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist, final long seed){
		require(bins==1024, "Suppressed dual ordering preserves the fixed 1024-bin geometry");
		writeSelected(path, keys, members, names, sequences, bins, blacklist, true, seed, true);
	}

	private static void writeSelected(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist,
			final boolean orderOnly, final long seed){
		writeSelected(path, keys, members, names, sequences, bins, blacklist, orderOnly, seed, false);
	}

	private static void writeSelected(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist,
			final boolean orderOnly, final long seed, final boolean suppressed){
		writeSelected(path, keys, members, names, sequences, bins, blacklist, orderOnly, seed, suppressed, 0, BONUS_NONE);
	}

	private static void writeSelected(final String path, final char[][] keys, final int[] members,
			final String[] names, final byte[][] sequences, final int bins, final LongHashSet blacklist,
			final boolean orderOnly, final long seed, final boolean suppressed, final int bonus, final int bonusRule){
		validateModels(names, sequences);
		validateOrderingBonus(bonus, bonusRule);
		require(bonus==0 || orderOnly && !suppressed, "Consensus bonuses require their distinct ordering format");
		final String family=memberFamily(names);
		require(orderOnly || family.equals("RF00174"), "Historical consensus-bonus profiles bind only RF00174");
		require(seed>=0, "Every ordering profile must record the same explicit seed used by its queries");
		require((bins==512 || bins==1024) && keys!=null && members!=null && keys.length==names.length
			&& members.length==names.length && blacklist!=null && (blacklist.size()==0 || blacklist.size()==2304),
			"Selected export binds every actual model/voter count to the explicit512/1024layout and padding set");
		for(int model=0; model<keys.length; model++){
			require(keys[model]!=null && keys[model].length==bins && (orderOnly ? members[model]>=0 : members[model]>0),
				"Every model requires its complete key row and actual voter count");
			if(members[model]==0 && bonus==0){
				for(final char key:keys[model]){require(key==0, "A zero-voter profile must remain empty");}
			}else if(members[model]==0){
				validateConsensusOnly(keys[model], sequences[model], bins, seed, blacklist);
			}
		}
		final long[] words=blacklist.toArray();
		int shortWords=0;
		for(final long word:words){if(shortPaddingWord(word)){shortWords++;}}
		require(words.length==0 || shortWords==256 && words.length-shortWords==2048, "Selected export requires the complete256short/2048long padding set");
		RrnaResourceIO.fresh(path, "");
		try(Output output=new Output(path)){
			final String prefix=suppressed ? SUPPRESSED_ORDER_PREFIX+"\t"+family
				: family.equals("RF00174") ? ORDER_PREFIX : MEMBER_ORDER_PREFIX+"\t"+family;
			if(bonus>0){
				output.writer.println(BONUS_ORDER_PREFIX+"\t"+family+"\t"+bins+"\t"+seed+"\t"+ORDER_GEOMETRY
					+"\t"+BONUS_SCORING+"\t"+bonus+"\t"+bonusRuleName(bonusRule));
			}else{
				output.writer.println(orderOnly ? prefix+"\t"+bins+"\t"+seed+"\t"+ORDER_GEOMETRY+"\t"+(suppressed ? SUPPRESSED_SCORING : ORDER_SCORING)
					: bins==512 ? SELECTED512_FORMAT : SELECTED1024_FORMAT);
			}
			final ByteBuilder row=new ByteBuilder(16384);
			output.writer.print(row.clear().append("#PADDING_BLACKLIST\t").append(blacklist.size()).nl());
			Arrays.sort(words);
			for(final long word:words){output.writer.print(row.clear().append("#PADDING_WORD\t").append(word).nl());}
			for(int model=0; model<keys.length; model++){
				row.clear().append(model).tab().append(id(names[model])).tab();
				for(final byte base:sequences[model]){row.append(upper(base));}
				row.tab().append(members[model]);
				for(final char key:keys[model]){row.tab().append((int)key);}
				output.writer.print(row.nl());
			}
		}
	}

	private static void write(final String path, final NcrnaDdlProfile[] profiles,
			final String[] names, final byte[][] sequences, final int buckets, final String format){
		validateModels(names, sequences);
		require(memberFamily(names).equals("RF00174"), "Historical profile exports bind only RF00174");
		require(profiles!=null && profiles.length==names.length,
			"Member profiles must retain the complete ordered caller model roster");
		for(int model=0; model<profiles.length; model++){
			require(profiles[model]!=null && profiles[model].buckets()==buckets && profiles[model].members>0,
				"Export requires one voted profile of the declared width per model: model="+model+", buckets="+buckets);
		}
		RrnaResourceIO.fresh(path, "");
		try(Output out=new Output(path)){
			out.writer.println(format);
			final ByteBuilder b=new ByteBuilder(16384);
			for(int model=0; model<profiles.length; model++){
				b.clear().append(model).tab().append(id(names[model])).tab();
				for(final byte base:sequences[model]){b.append(upper(base));}
				b.tab().append(profiles[model].members);
				for(final char key:profiles[model].keys()){b.tab().append((int)key);}
				out.writer.print(b.nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Scratch       ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocate once per caller worker; neither the resource nor its keys are copied. */
	public Scorer newScorer(){return new Scorer(this);}

	public static final class Scorer{
		private Scorer(final NcrnaDdlProfileResource resource_){
			resource=resource_;
			query=resource.dual ? null : new NcrnaDdlQueryKeys(BUCKETS, K);
			dualQuery=resource.dual && !resource.compact && !resource.selected && !resource.orderOnly ? new NcrnaDdlDualKeys(9, 10, DUAL_BUCKETS, true) : null;
			partitionQuery=resource.compact ? new NcrnaDdlPartitionKeys(7, resource.buckets()) : null;
			selectedBuilder=resource.selected || resource.orderOnly ? NcrnaDdlSketchBuilder.forKeyArrays(resource.buckets(), resource.seed, resource.blacklist) : null;
			selectedQuery=resource.selected ? new char[resource.buckets()] : null;
			nativeQuery=resource.orderOnly ? DynamicDemiLog.create(resource.buckets(), 10, resource.seed, 0, false, false) : null;
			matches=new int[resource.models()];
			halves=dualQuery!=null ? new int[][]{new int[resource.models()], matches} : null;
			eligible=resource.gc==null ? null : new boolean[resource.models()];
		}

		/** Replaces the complete score vector with one immutable-index traversal. */
		public void score(final byte[] sequence, final float[] scores){
			require(scores!=null && scores.length==resource.models(),
				"Worker scores must cover every bound model before interleaving");
			if(resource.orderOnly){
				selectedBuilder.fillKeys(sequence, nativeQuery.maxArray());
				if(specificity==null){resource.index.query(nativeQuery, matches);}
				else{resource.index.querySpecificity(nativeQuery, matches, specificity, postingLengths);}
				if(resource.consensusIndex!=null){resource.consensusIndex.accumulateCounts(nativeQuery, matches);}
				if(resource.suppressed){
					int occupied=0;
					for(final char key:nativeQuery.maxArray()){if(key!=0){occupied++;}}
					for(final int count:matches){
						assert(count>=0 && count<=occupied) : "Disjoint reference keys permit at most one match per occupied query bin";
					}
				}
			}else if(resource.selected){
				selectedBuilder.fillKeys(sequence, selectedQuery);
				resource.index.query(selectedQuery, matches);
			}else if(resource.compact){
				partitionQuery.build(sequence);
				partitionQuery.count(resource.index, matches);
			}else if(resource.dual){
				dualQuery.build(sequence);
				dualQuery.countHalves(resource.index, halves);
			}else{
				//TODO: Probable bug - this ACGTN-only codec receives caller windows
				//from NcrnaScavenger.copyRegionUpper without IUPAC normalization.
				//Symbols such as R can contribute undefined spanning words. Preserve
				//B/C here; record this limitation before interpreting their comparison.
				query.build(sequence, null, false, true);
				query.count(resource.index, matches);
			}
			if(halves==null){
				for(int model=0; model<scores.length; model++){scores[model]=(float)matches[model];}
			}else{
				for(int model=0; model<scores.length; model++){
					scores[model]=(float)matches[model]+.5f*(float)halves[0][model];
				}
			}
			if(resource.gc!=null){resource.gc.apply(sequence, scores, eligible);}
		}

		/** Borrowed until the next score call; null retains the original complete order. */
		public boolean[] eligibility(){return eligible;}

		/** Builds one query and probes the validated member/replacement indexes.
		 * Both destinations are complete independent vectors, never summed. The
		 * pair's identity binds model ordinals as well as the query recipe. */
		public void scorePair(final byte[] sequence, final ConsensusPair pair,
				final float[] memberScores, final float[] replacementScores){
			require(pair!=null && pair.memberOnly==resource && specificity==null
				&& memberScores!=null && replacementScores!=null && memberScores!=replacementScores
				&& memberScores.length==resource.models() && replacementScores.length==resource.models(),
				"Paired scoring requires its exact validated member resource and separate complete score buffers");
			score(sequence, memberScores);
			pair.replacement.index.query(nativeQuery, matches);
			for(int model=0; model<matches.length; model++){replacementScores[model]=matches[model];}
		}

		/** Explicit query-only filter experiment on the standard paired indexes.
		 * Disallowed words never enter the native query's bin maxima. This leaves
		 * scorePair and every reference profile/file unchanged.
		 */
		public void scorePairWithReferenceKeys(final byte[] sequence, final ConsensusPair pair,
				final float[] memberScores, final float[] replacementScores, final NcrnaDdlReferenceKeys allowed){
			require(allowed!=null, "Filtered scoring requires its explicit existing-index lookup");
			scorePairObserved(sequence, pair, memberScores, replacementScores, allowed, null);
		}
		/** Optional census outside timing; null allowed preserves the ordinary query. */
		public void scorePairObserved(final byte[] sequence, final ConsensusPair pair,
				final float[] memberScores, final float[] replacementScores, final NcrnaDdlReferenceKeys allowed,
				final NcrnaDdlSketchBuilder.QueryStats stats){
			require(pair!=null && pair.memberOnly==resource && specificity==null
				&& (allowed==null || allowed.matches(pair))
				&& memberScores!=null && replacementScores!=null && memberScores!=replacementScores
				&& memberScores.length==resource.models() && replacementScores.length==resource.models(),
				"Filtered paired scoring requires its exact loaded reference pair, posting lookup and separate complete score buffers");
			selectedBuilder.fillKeys(sequence, nativeQuery.maxArray(), allowed, stats);
			resource.index.query(nativeQuery, matches);
			for(int model=0; model<matches.length; model++){memberScores[model]=matches[model];}
			pair.replacement.index.query(nativeQuery, matches);
			for(int model=0; model<matches.length; model++){replacementScores[model]=matches[model];}
		}

		/** Explicit member-only tie experiment; default query path is unchanged. */
		public void enableSpecificity(){
			require(resource.orderOnly && resource.consensusIndex==null && resource.gc==null,
				"Tie specificity is qualified only for unchanged member-only ordering profiles");
			require(specificity==null, "Configure worker specificity exactly once before querying");
			specificity=new double[resource.models()];
			postingLengths=new int[resource.buckets()];
		}
		public double[] specificitySums(){return specificity;}
		public int[] postingLengths(){return postingLengths;}
		public char[] queryKeys(){return nativeQuery.maxArray();}

		private final NcrnaDdlProfileResource resource;
		private final NcrnaDdlQueryKeys query;
		private final NcrnaDdlDualKeys dualQuery;
		private final NcrnaDdlPartitionKeys partitionQuery;
		private final NcrnaDdlSketchBuilder selectedBuilder;
		private final char[] selectedQuery;
		private final DynamicDemiLog nativeQuery;
		private final int[] matches;
		private final int[][] halves;
		private final boolean[] eligible;
		private double[] specificity;
		private int[] postingLengths;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private static void validateModels(final String[] names, final byte[][] sequences){
		require(names!=null && sequences!=null && names.length>0 && names.length==sequences.length,
			"Model binding requires nonempty aligned names and consensus sequences");
		final String prefix="cobalamin_"+memberFamily(names)+"_";
		final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		for(int model=0; model<names.length; model++){
			require(names[model]!=null && sequences[model]!=null && sequences[model].length>0,
				"Every profile row needs an intact named consensus: model="+model);
			final String name=RrnaResourceIO.first(names[model]);
			require(name.startsWith(prefix) && seen.add(name),
				"Member profile IDs must be unique and family-qualified: model="+model);
			for(final byte base:sequences[model]){
				final byte normalized=upper(base);
				if(normalized<'A' || normalized>'Z'){
					throw new IllegalArgumentException("Consensus bytes must be ASCII letters: model="+model);
				}
			}
		}
	}
	private static String memberFamily(final String[] names){return NcrnaDdlMembers.memberFamily(names);}
	private static byte[] id(final String name){return RrnaResourceIO.first(name).getBytes(StandardCharsets.UTF_8);}
	private static boolean shortPaddingWord(final long tagged){
		final long shortRaw=tagged^BuildNcrnaDdlLayoutCounts.wordTag(8, true, true);
		if(shortRaw>=0 && shortRaw<(1L<<16)){return true;}
		final long longRaw=tagged^BuildNcrnaDdlLayoutCounts.wordTag(10, false, true);
		require(longRaw>=0 && longRaw<(1L<<20) && (longRaw&LONG_PADDING_MASK)==longRaw,
			"Padding resource words must preserve the native forward8 or middle-masked10 tag/span");
		return false;
	}
	private static byte upper(final byte base){return base>='a' && base<='z' ? (byte)(base-32) : base;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public int models(){return members.length;}
	/** Dual-length, compact and selected resources use sketch-only caller ordering. */
	public boolean dual(){return dual;}
	public boolean compact(){return compact;}
	public boolean selected(){return selected;}
	public boolean orderOnly(){return orderOnly;}
	/** Cold compatibility check for separate arms using one frozen query recipe. */
	public boolean sameOrderingQuery(final NcrnaDdlProfileResource other){
		if(other==null || !orderOnly || !other.orderOnly || consensus() || other.consensus()
			|| buckets()!=other.buckets() || seed!=other.seed){return false;}
		final long[] expected=blacklist.toArray(), actual=other.blacklist.toArray();
		Arrays.sort(expected);
		Arrays.sort(actual);
		return Arrays.equals(expected, actual);
	}
	/** True only when a separate runtime consensus index contributes raw matches. */
	public boolean consensus(){return consensusIndex!=null;}
	/** Runtime bin replacement is distinct from additive scoring and stored votes. */
	public boolean consensusReplaced(){return consensusReplaced;}
	public boolean suppressed(){return suppressed;}
	/** Only the explicit ordering-bonus format; historical selected formats retain their own fixed header contract. */
	public int orderingBonus(){return orderingBonus;}
	public int orderingBonusRule(){return orderingBonusRule;}
	public long seed(){return seed;}
	public int buckets(){return index.buckets();}
	public int members(final int model){return members[model];}
	public char key(final int model, final int bucket){return orderOnly ? nativeProfiles[model].maxArray()[bucket] : keys[model][bucket];}
	private final char[][] keys;
	private final DynamicDemiLog[] nativeProfiles;
	private final int[] members;
	private final boolean dual, compact, selected, orderOnly, suppressed, consensusReplaced;
	private final int orderingBonus, orderingBonusRule;
	private final long seed;
	private final LongHashSet blacklist;
	private final DDLIndexSmall index, consensusIndex;
	private final NcrnaWindowGcProfile gc;

	/** Temporary decoded rows are consumed by exactly one immutable index build. */
	private static final class Rows{
		Rows(final char[][] keys_, final int[] members_, final boolean dual_, final boolean compact_,
				final boolean selected_, final LongHashSet blacklist_, final boolean orderOnly_, final long seed_, final boolean suppressed_,
				final int bonus_, final int bonusRule_){
			assert(keys_!=null && members_!=null && keys_.length==members_.length && keys_.length>0)
				: "Decoded profile rows must retain the validated model/voter axis before index construction";
			keys=keys_;
			members=members_;
			dual=dual_;
			compact=compact_;
			selected=selected_;
			blacklist=blacklist_;
			orderOnly=orderOnly_;
			suppressed=suppressed_;
			seed=seed_;
			bonus=bonus_;
			bonusRule=bonusRule_;
		}
		final char[][] keys;
		final int[] members;
		final boolean dual, compact, selected, orderOnly, suppressed;
		final LongHashSet blacklist;
		final long seed;
		final int bonus, bonusRule;
	}

	/** With no real voters, equal positive pseudo-votes plus greatest-key ties
	 * select exactly the native consensus maximum for BOTH declared bonus rules. */
	private static void validateConsensusOnly(final char[] keys, final byte[] sequence, final int bins,
			final long seed, final LongHashSet blacklist){
		final char[] expected=new char[bins];
		NcrnaDdlSketchBuilder.forKeyArrays(bins, seed, blacklist).fillKeys(sequence, expected);
		require(Arrays.equals(keys, expected), "Zero-real-voter bonus keys must be the actual native consensus maxima, never invented keys");
	}

	public static int parseBonusRule(final String name){
		require(name!=null, "An explicit consensus rule needs a value");
		if(name.equalsIgnoreCase("none")){return BONUS_NONE;}
		if(name.equalsIgnoreCase("maxima")){return BONUS_MAXIMA;}
		if(name.equalsIgnoreCase("allkeys")){return BONUS_ALL_KEYS;}
		throw new IllegalArgumentException("Consensus rule must be none, maxima or allkeys: "+name);
	}
	public static String bonusRuleName(final int rule){
		if(rule==BONUS_NONE){return "none";}
		if(rule==BONUS_MAXIMA){return "maxima";}
		if(rule==BONUS_ALL_KEYS){return "allkeys";}
		throw new IllegalArgumentException("Unknown consensus vote rule: "+rule);
	}
	public static void validateOrderingBonus(final int bonus, final int rule){
		require(bonus>=0 && (bonus==0 ? rule==BONUS_NONE : rule==BONUS_MAXIMA || rule==BONUS_ALL_KEYS),
			"Zero bonus pairs only with rule none; a positive bonus requires explicit maxima or allkeys");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final int K=10, BUCKETS=2048;
	public static final int BONUS_NONE=0, BONUS_MAXIMA=1, BONUS_ALL_KEYS=2;
	private static final String ORDER_PREFIX="#NCRNA_DDL_RF00174_ORDER_V1";
	private static final String MEMBER_ORDER_PREFIX="#NCRNA_DDL_MEMBER_ORDER_V1";
	private static final String SUPPRESSED_ORDER_PREFIX="#NCRNA_DDL_MEMBER_SUPPRESSED_ORDER_V1";
	private static final String BONUS_ORDER_PREFIX="#NCRNA_DDL_MEMBER_BONUS_ORDER_V1";
	private static final String BONUS_SCORING="MEMBER_PLURALITY_PLUS_CONSENSUS_BONUS_QUERY_EQUAL_BIN_MAX";
	private static final String SUPPRESSED_SCORING="MEMBER_BIN_MAX_WITHOUT_OWN_CONSENSUS_MAX_PLUS_RUNTIME_CONSENSUS";
	private static final String ORDER_GEOMETRY="FORWARD8_QUARTER_MASKED10_THREE_QUARTERS_NATIVE4_12";
	private static final String ORDER_SCORING="MEMBER_ONLY_QUERY_EQUAL_BIN_MAX";
	public static final int DUAL_BUCKETS=4096;
	private static final long LONG_PADDING_MASK=AbstractKmerTable.makeMiddleMask(10, false, 1);
	public static final String FORMAT="#NCRNA_DDL_RF00174_V1\t10\t2048\tFORWARD_MIDDLE_RIGHT\tRAW_ALL_QUERY_KEYS\t12345\t4\tCONSENSUS_UPPERCASE";
	public static final String DUAL_FORMAT="#NCRNA_DDL_RF00174_V2\tPARTITION\t9\t10\t4096\tLONG_PLUS_HALF_SHORT\tFORWARD_MIDDLE_RIGHT_LONG\tRAW_ALL_QUERY_KEYS\t12345\t4\tCONSENSUS_UPPERCASE";
	public static final String COMPACT512_FORMAT="#NCRNA_DDL_RF00174_V3\t512\tPLAIN7_QUARTER\tMASKED10_HALF\tTAGGED_SHARED_CONSENSUS_QUARTER\tMEMBER_PLURALITY_GREATEST_KEY\tCONSENSUS_MAX_KEY\tRAW_ALL_QUERY_KEYS\t12345\t4\tCONSENSUS_UPPERCASE";
	public static final String COMPACT1024_FORMAT=COMPACT512_FORMAT.replace("\t512\t", "\t1024\t");
	public static final String SELECTED512_FORMAT="#NCRNA_DDL_RF00174_V4\t512\tFORWARD8_QUARTER_MASKED10_THREE_QUARTERS\tNATIVE4_12\tMEMBER_BIN_MAX_PLURALITY_GREATEST_KEY\tCONSENSUS100_EACH_DISTINCT_ELIGIBLE_BIN_KEY\tQUERY_EQUAL_BIN_MAX\t12345\tCONSENSUS_UPPERCASE";
	public static final String SELECTED1024_FORMAT=SELECTED512_FORMAT.replace("\t512\t", "\t1024\t");
	private static final byte[] HEADER=FORMAT.getBytes(StandardCharsets.UTF_8);
	private static final byte[] DUAL_HEADER=DUAL_FORMAT.getBytes(StandardCharsets.UTF_8);
	private static final byte[] COMPACT512_HEADER=COMPACT512_FORMAT.getBytes(StandardCharsets.UTF_8);
	private static final byte[] COMPACT1024_HEADER=COMPACT1024_FORMAT.getBytes(StandardCharsets.UTF_8);
	private static final byte[] SELECTED512_HEADER=SELECTED512_FORMAT.getBytes(StandardCharsets.UTF_8);
	private static final byte[] SELECTED1024_HEADER=SELECTED1024_FORMAT.getBytes(StandardCharsets.UTF_8);
}
