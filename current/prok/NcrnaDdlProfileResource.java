package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import ddl.DDLIndexSmall;
import map.LongHashSet;
import kmer.AbstractKmerTable;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Model-bound member-voted keys for the fixed RF00174 caller experiment.
 * V1 fixes forward middle-masked k10/2048 for native-first hybrid ordering.
 * V2 fixes partitioned plain9+masked10/4096 and long+.5*short, sketch-only.
 * V3 fixes plain7/masked10/shared-consensus quarters at512/1024, sketch-only.
 * V4 fixes forward8/masked10 native4+12 quarter-short keys, real member
 * plurality plus100perdistinctconsensuskey, equal-bin queries and saved padding.
 * Both use seed12345, exponent4 and raw-query-key bucket matches. Every row retains its
 * original ordinal, first-token ID, uppercase consensus bytes and voter count.
 * Loading compares the actual model identities and sequences, without hashes.
 * Keys are private and read-only after construction; query scratch is per worker.
 * No caller enables or loads this resource by default.
 * @author Brian, Raiden
 */
public final class NcrnaDdlProfileResource{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaDdlProfileResource(final char[][] keys_, final int[] members_, final boolean dual_, final boolean compact_,
			final boolean selected_, final LongHashSet blacklist_){
		assert(keys_!=null && members_!=null && keys_.length==members_.length && keys_.length>0)
			: "Resource rows must bind one profile and voter count per current model";
		keys=keys_;
		members=members_;
		dual=dual_;
		compact=compact_;
		selected=selected_;
		blacklist=blacklist_;
		index=new DDLIndexSmall(keys);
		gc=null;
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
		blacklist=source.blacklist;
		index=source.index;
		gc=gc_;
	}

	/** The model-GC loader binds names to current ordinals, not to row positions. */
	public NcrnaDdlProfileResource withGc(final NcrnaWindowGcProfile profile){
		return new NcrnaDdlProfileResource(this, profile);
	}

	/** Reads exactly the expected ordered models; caller-owned arrays are not modified. */
	public static NcrnaDdlProfileResource load(final String path, final String[] names, final byte[][] sequences){
		validateModels(names, sequences);
		final int[] members=new int[names.length];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path)){
			final byte[] header=in.reader.nextLine();
			final boolean selected512=Arrays.equals(header, SELECTED512_HEADER), selected1024=Arrays.equals(header, SELECTED1024_HEADER);
			final boolean selected=selected512 || selected1024;
			final boolean compact512=Arrays.equals(header, COMPACT512_HEADER);
			final boolean compact1024=Arrays.equals(header, COMPACT1024_HEADER);
			final boolean compact=compact512 || compact1024;
			final boolean dual=selected || compact || Arrays.equals(header, DUAL_HEADER);
			require(dual || Arrays.equals(header, HEADER), "Unsupported profile header: "+path);
			final int buckets=selected512 || compact512 ? 512 : selected1024 || compact1024 ? 1024 : dual ? DUAL_BUCKETS : BUCKETS;
			final LongHashSet blacklist=new LongHashSet(4);
			if(selected){
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
				require(members[model]>0, "Each model profile requires actual kept member votes: model="+model);
				for(int bucket=0; bucket<buckets; bucket++){
					final int key=p.parseInt(bucket+4);
					if(key<0 || key>Character.MAX_VALUE){
						throw new IllegalArgumentException("Absolute profile keys must fit unsigned16bits: model="+model+", bucket="+bucket);
					}
					keys[model][bucket]=(char)key;
				}
			}
			require(in.reader.nextLine()==null, "Profile resource has rows beyond the actual caller model count");
			return new NcrnaDdlProfileResource(keys, members, dual, compact, selected, blacklist);
		}
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
		validateModels(names, sequences);
		require((bins==512 || bins==1024) && keys!=null && members!=null && keys.length==names.length
			&& members.length==names.length && blacklist!=null && (blacklist.size()==0 || blacklist.size()==2304),
			"Selected export binds every actual model/voter count to the explicit512/1024layout and padding set");
		for(int model=0; model<keys.length; model++){
			require(keys[model]!=null && keys[model].length==bins && members[model]>0, "Every selected model requires its complete key row and actual voters");
		}
		final long[] words=blacklist.toArray();
		int shortWords=0;
		for(final long word:words){if(shortPaddingWord(word)){shortWords++;}}
		require(words.length==0 || shortWords==256 && words.length-shortWords==2048, "Selected export requires the complete256short/2048long padding set");
		RrnaResourceIO.fresh(path, "");
		try(Output output=new Output(path)){
			output.writer.println(bins==512 ? SELECTED512_FORMAT : SELECTED1024_FORMAT);
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
			dualQuery=resource.dual && !resource.compact && !resource.selected ? new NcrnaDdlDualKeys(9, 10, DUAL_BUCKETS, true) : null;
			partitionQuery=resource.compact ? new NcrnaDdlPartitionKeys(7, resource.buckets()) : null;
			selectedBuilder=resource.selected ? NcrnaDdlSketchBuilder.forKeyArrays(resource.buckets(), 12345, resource.blacklist) : null;
			selectedQuery=resource.selected ? new char[resource.buckets()] : null;
			matches=new int[resource.models()];
			halves=dualQuery!=null ? new int[][]{new int[resource.models()], matches} : null;
			eligible=resource.gc==null ? null : new boolean[resource.models()];
		}

		/** Replaces the complete score vector with one immutable-index traversal. */
		public void score(final byte[] sequence, final float[] scores){
			require(scores!=null && scores.length==resource.models(),
				"Worker scores must cover every bound model before interleaving");
			if(resource.selected){
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

		private final NcrnaDdlProfileResource resource;
		private final NcrnaDdlQueryKeys query;
		private final NcrnaDdlDualKeys dualQuery;
		private final NcrnaDdlPartitionKeys partitionQuery;
		private final NcrnaDdlSketchBuilder selectedBuilder;
		private final char[] selectedQuery;
		private final int[] matches;
		private final int[][] halves;
		private final boolean[] eligible;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private static void validateModels(final String[] names, final byte[][] sequences){
		require(names!=null && sequences!=null && names.length>0 && names.length==sequences.length,
			"Model binding requires nonempty aligned names and consensus sequences");
		final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		for(int model=0; model<names.length; model++){
			require(names[model]!=null && sequences[model]!=null && sequences[model].length>0,
				"Every profile row needs an intact named consensus: model="+model);
			final String name=RrnaResourceIO.first(names[model]);
			require(name.startsWith("cobalamin_RF00174_") && seen.add(name),
				"RF00174 profile IDs must be unique and family-qualified: model="+model);
			for(final byte base:sequences[model]){
				final byte normalized=upper(base);
				if(normalized<'A' || normalized>'Z'){
					throw new IllegalArgumentException("Consensus bytes must be ASCII letters: model="+model);
				}
			}
		}
	}
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

	public int models(){return keys.length;}
	/** Dual-length, compact and selected resources use sketch-only caller ordering. */
	public boolean dual(){return dual;}
	public boolean compact(){return compact;}
	public boolean selected(){return selected;}
	public int buckets(){return keys[0].length;}
	public int members(final int model){return members[model];}
	public char key(final int model, final int bucket){return keys[model][bucket];}
	private final char[][] keys;
	private final int[] members;
	private final boolean dual, compact, selected;
	private final LongHashSet blacklist;
	private final DDLIndexSmall index;
	private final NcrnaWindowGcProfile gc;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final int K=10, BUCKETS=2048;
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
