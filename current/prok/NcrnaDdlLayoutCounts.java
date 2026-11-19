package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Validated plain-DDL geometry and immutable saved counts. Version2 uses only
 * key arrays and active score terms. Version1 remains explicitly historical
 * offline evidence; reading its counts does not restore an origin-aware scorer.
 * One sequential reader retains at most one layout's primitive window counts.
 * @author Ganyu, Raiden
 */
public final class NcrnaDdlLayoutCounts{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlLayoutCounts(final String path, final int models){
		require(path!=null && models>0 && models<=65536, "Layout metadata needs its complete bound model roster");
		final ArrayList<Layout> list=new ArrayList<Layout>();
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		String family=null;
		try(Input input=new Input(path)){
			final byte[] header=input.reader.nextLine();
			if(Arrays.equals(header, V1_LAYOUT_BYTES)){schemaVersion=1;}
			else{
				require(Arrays.equals(header, V2_LAYOUT_BYTES), "Unsupported layout metadata header: "+path);
				schemaVersion=2;
			}
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==(schemaVersion==1 ? 27 : 31), "Layout metadata must retain every declared versioned field");
				final Layout layout=new Layout(p, models, schemaVersion);
				require(ids.add(layout.id), "Duplicate layout identifier: "+layout.id);
				if(family==null){family=layout.family;}
				require(family.equals(layout.family), "One count invocation must retain one explicit family");
				list.add(layout);
			}
		}
		require(!list.isEmpty(), "A layout manifest cannot be empty");
		layouts=list.toArray(new Layout[list.size()]);
		historical=schemaVersion==1;
	}

	/** The accepted unit control is explicit and never inferred from its name. */
	public int control(final String id){
		for(int i=0; i<layouts.length; i++){
			final Layout layout=layouts[i];
			if(layout.id.equals(id)){
				require(layout.legacy && (layout.queryForm.equals("RAW_KEYS") || layout.queryForm.equals("QUERY_DDL"))
					&& layout.memberBins==384 && layout.consensusBins==128
					&& layout.memberShortBins==128 && layout.memberLongBins==256
					&& layout.terms==(layout.historical ? 4 : 3),
					"The explicit unit control must preserve accepted member384/consensus128 hashing and bin allocation");
				return i;
			}
		}
		throw new IllegalArgumentException("Explicit legacy control is absent: "+id);
	}

	public void write(final String path){
		try(Output output=new Output(path)){
			output.writer.println(schemaVersion==1 ? LAYOUT_HEADER : LAYOUT_HEADER_V2);
			for(final Layout layout:layouts){output.writer.println(layout.row);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Layout Metadata     ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Layout{
		private Layout(final LineParser1 p, final int expectedModels, final int version){
			schemaVersion=version;
			historical=version==1;
			final boolean lowBit;
			if(historical){
				id=p.parseString(0);
				family=p.parseString(1);
				memberMode=p.parseString(2);
				memberBins=p.parseInt(3);
				memberShortBins=p.parseInt(4);
				memberLongBins=p.parseInt(5);
				consensusMode=p.parseString(6);
				consensusBins=p.parseInt(7);
				exponentBits=p.parseInt(8);
				mantissaBits=p.parseInt(9);
				legacy=p.termEquals("LEGACY_PARTITION_ORIGINS", 10);
				lowBit=p.termEquals("LOW_KEY_BIT", 10);
				require(legacy || lowBit || p.termEquals("SIDECAR_ORIGIN", 10) || p.termEquals("PARTITION_AND_SIDECAR", 10),
					"Unknown historical length treatment; it must never be inferred from a display name");
				overflowPolicy=p.parseString(11);
				models=p.parseInt(12);
				shortK=7;
				longK=10;
				longMask="MIDDLE_RIGHT";
				consensusBonus=0;
				hashNamespace=legacy ? "LEGACY_PARTITIONS" : "LENGTH_TAGGED";
				queryForm="RAW_KEYS";
				terms=consensusMode.equals("BOTH") ? 4 : 3;
				require(terms==4 || consensusMode.equals("LONG_ONLY"), "Unknown historical consensus mode");
				termNames=terms==4 ? new String[]{"member_short", "member_long", "consensus_short", "consensus_long"}
					: new String[]{"member_short", "member_long", "consensus_long"};
				referenceTerm=1;
				memberKeyPayloadBytes=p.parseLong(13);
				memberOriginPayloadBytes=p.parseLong(14);
				consensusKeyPayloadBytes=p.parseLong(15);
				validateHistorical(p, lowBit);
			}else{
				require(p.parseInt(0)==2, "The version2 header requires schema_version=2 on every row");
				id=p.parseString(1);
				family=p.parseString(2);
				shortK=p.parseInt(3);
				longK=p.parseInt(4);
				longMask=p.parseString(5);
				memberMode=p.parseString(6);
				memberBins=p.parseInt(7);
				memberShortBins=p.parseInt(8);
				memberLongBins=p.parseInt(9);
				consensusMode=p.parseString(10);
				consensusBins=p.parseInt(11);
				consensusBonus=p.parseInt(12);
				exponentBits=p.parseInt(13);
				mantissaBits=p.parseInt(14);
				overflowPolicy=p.parseString(15);
				termNames=p.parseString(16).split(",", -1);
				terms=termNames.length;
				referenceTerm=p.parseInt(17);
				models=p.parseInt(18);
				memberKeyPayloadBytes=p.parseLong(19);
				memberOriginPayloadBytes=0;
				consensusKeyPayloadBytes=p.parseLong(20);
				hashNamespace=p.parseString(29);
				queryForm=p.parseString(30);
				legacy=hashNamespace.equals("LEGACY_PARTITIONS");
				lowBit=false;
				validateCurrent(p);
			}
			require(!id.isEmpty() && !family.isEmpty() && models==expectedModels,
				"Every layout needs explicit identifiers and the unchanged complete model count");
			require(memberMode.equals("MIXED") || memberMode.equals("PARTITIONED"), "Unknown member layout mode");
			require(memberBins>0 && memberBins<=65535 && consensusBins>=0 && consensusBins<=65535,
				"Saved unsigned16 counts require positive member width and representable sketch counts");
			require(memberMode.equals("MIXED") ? memberShortBins==0 && memberLongBins==0
				: memberShortBins>0 && memberLongBins>0 && memberShortBins+memberLongBins==memberBins,
				"Mixed bins are shared; explicit short/long partitions must conserve the complete member sketch");
			require(terms>=1 && terms<=4 && referenceTerm>=0 && referenceTerm<terms,
				"Active terms and their fixed reference must fit the four-slot count interchange");
			require(memberKeyPayloadBytes==2L*models*memberBins && consensusKeyPayloadBytes==2L*models*consensusBins,
				"Stored member and runtime consensus key payloads are exactly two bytes per model/bin");
			memberTerms=historical || memberMode.equals("PARTITIONED") ? 2 : 1;
			for(int term=0; term<4; term++){
				if(historical){limits[term]=term<2 ? memberBins : consensusBins;}
				else if(term<terms){
					final String name=termNames[term];
					limits[term]=name.equals("member_total") ? memberBins : name.equals("member_short")
						? memberShortBins : name.equals("member_long") ? memberLongBins : consensusBins;
				}
			}
			row=p.line().clone();
		}

		private void validateHistorical(final LineParser1 p, final boolean lowBit){
			require(consensusBins>0 && (memberMode.equals("MIXED")
				|| memberShortBins==memberBins/3 && memberLongBins==memberBins-memberShortBins),
				"Version1 preserves its original separate-consensus and one-third member geometry");
			require(exponentBits>0 && mantissaBits>0 && exponentBits+mantissaBits+(lowBit ? 1 : 0)<=16,
				"Historical key fields and any in-key length bit must fit sixteen bits");
			require(legacy ? memberMode.equals("PARTITIONED") && exponentBits==4 && mantissaBits==12 && overflowPolicy.equals("WRAP")
				: overflowPolicy.equals("CLAMP"), "Historical input must preserve its original overflow policy");
			for(int field=13; field<=24; field++){
				require(p.parseLong(field)>=0, "Historical payload accounting cannot be negative");
			}
			require(memberOriginPayloadBytes==(long)models*memberBins
				|| memberOriginPayloadBytes==0 && (!memberMode.equals("MIXED") || lowBit),
				"Historical mixed-sidecar byte costs must remain visible, never relabelled as plain DDL storage");
			require(p.parseLong(16)==(long)models*consensusBins, "Historical runtime origin accounting must retain its original units");
			final long slots=p.parseLong(17), populated=p.parseLong(18), empty=p.parseLong(19);
			final long reference=p.parseLong(20), posting=p.parseLong(22), extra=p.parseLong(24);
			require(p.parseLong(23)==slots*models && extra==2L*models,
				"Historical metadata must account for its then-owned length arrays; none are reconstructed by this reader");
			validateIndex(slots, populated, empty, reference, p.parseLong(21), posting, extra);
			seconds(p, 25);
			seconds(p, 26);
		}

		private void validateCurrent(final LineParser1 p){
			require(shortK>0 && shortK<longK && longK<=15 && (longMask.equals("MIDDLE_RIGHT") || longMask.equals("NONE")),
				"Current geometry must bind distinct short/long word lengths and the actual long-word mask");
			require(exponentBits==4 && mantissaBits==12 && overflowPolicy.equals("WRAP"),
				"The expanded recipe uses full4+12 keys with identical wrapping on profiles and queries");
			require(legacy || hashNamespace.equals("LENGTH_TAGGED") || hashNamespace.equals("NATIVE_DDL"),
				"Hash namespace must explicitly distinguish historical custom mapping from native DDL objects");
			if(hashNamespace.equals("NATIVE_DDL")){
				require(queryForm.equals("QUERY_DDL") && (memberShortBins==0 ? powerOfTwo(memberBins)
					: powerOfTwo(memberShortBins) && powerOfTwo(memberLongBins))
					&& (consensusBins==0 || powerOfTwo(consensusBins)),
					"Native DDL metadata requires object queries and power-of-two width for each actual part");
			}
			require(queryForm.equals("RAW_KEYS") || queryForm.equals("QUERY_DDL"), "Query form must explicitly name raw probes or greatest-key-per-bin sketching");
			final boolean separate=consensusMode.equals("BOTH") || consensusMode.equals("LONG_ONLY");
			final boolean pooled=consensusMode.equals("POOL_BOTH") || consensusMode.equals("POOL_LONG");
			require(separate ? consensusBins>0 && consensusBonus==0
				: pooled ? consensusBins==0 && consensusBonus>0 : consensusMode.equals("NONE") && consensusBins==0 && consensusBonus==0,
				"Separate consensus owns runtime bins; pooled bonuses and W0 member-only rows own no separate sketch");
			final boolean partitioned=memberMode.equals("PARTITIONED");
			final String[] expected;
			if(partitioned){
				expected=separate ? new String[]{"member_short", "member_long", "consensus_total"}
					: new String[]{"member_short", "member_long"};
			}else{
				expected=separate ? new String[]{"member_total", "consensus_total"} : new String[]{"member_total"};
			}
			require(Arrays.equals(termNames, expected) && referenceTerm==(partitioned ? 1 : 0),
				"Current score terms must follow bin partitions or whole mixed sketches without origin metadata");
			if(legacy){
				final int total=memberBins+consensusBins;
				require(shortK==7 && longK==10 && longMask.equals("MIDDLE_RIGHT") && partitioned
					&& consensusMode.equals("BOTH") && (total==512 || total==1024 || total==2048)
					&& memberShortBins==total/4 && memberLongBins==total/2 && consensusBins==total/4,
					"Legacy namespace must retain the actual accepted quarter/half/quarter key mapping");
			}
			if(queryForm.equals("QUERY_DDL") && (legacy || hashNamespace.equals("NATIVE_DDL"))){
				for(int field=21; field<=26; field++){
					require(p.parseLong(field)==0, "Direct native DDL comparison owns no index arrays or posting payload");
				}
			}else{
				validateIndex(p.parseLong(21), p.parseLong(22), p.parseLong(23), p.parseLong(24), p.parseLong(25), p.parseLong(26), 0);
			}
			seconds(p, 27);
			seconds(p, 28);
		}

		private void validateIndex(final long slots, final long populated, final long empty, final long reference,
				final long pointers, final long posting, final long historicalExtraPointers){
			require(slots==(long)memberBins+consensusBins && populated>=0 && populated<=slots && empty==slots-populated,
				"Posting references must conserve real nonempty and shared-empty bin slots");
			require((reference==4 || reference==8) && pointers==(slots+historicalExtraPointers)*reference,
				"Index pointer bytes must use the declared JVM reference width and actual arrays only");
			require(posting>=0 && posting%4==0 && posting/4>=populated && posting/4<=slots*models,
				"Key-derived postings contain one packed four-byte entry per nonempty model/bin");
		}

		private int limit(final int term){return limits[term];}

		public final String id, family, memberMode, consensusMode, longMask, overflowPolicy, hashNamespace, queryForm;
		public final int schemaVersion, shortK, longK, memberBins, memberShortBins, memberLongBins;
		public final int consensusBins, consensusBonus, terms, referenceTerm, models, exponentBits, mantissaBits;
		public final String[] termNames;
		public final boolean legacy, historical;
		public final long memberKeyPayloadBytes, memberOriginPayloadBytes, consensusKeyPayloadBytes;
		private final int memberTerms;
		private final int[] limits=new int[4];
		private final byte[] row;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Count Reader        ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Reader implements AutoCloseable{
		public Reader(final String path, final NcrnaDdlLayoutCounts manifest, final NcrnaDdlSavedData data_){
			require(path!=null && manifest!=null && data_!=null, "Counts need a fixed manifest and the original saved-window roster");
			input=new Input(path);
			layouts=manifest.layouts;
			data=data_;
			try{RrnaResourceIO.header(input.reader, manifest.schemaVersion==1 ? COUNT_HEADER : COUNT_HEADER_V2, path);}
			catch(RuntimeException | Error failure){
				try{input.close();}catch(RuntimeException | Error cleanup){failure.addSuppressed(cleanup);}
				throw failure;
			}
		}

		/** Returns four-slot model-major rows, with active terms FIRST in both versions. */
		public char[][] next(final boolean retain){
			require(next<layouts.length, "Counts cannot introduce an undeclared layout");
			final Layout layout=layouts[next++];
			final char[][] rows=retain ? new char[data.windows.size()][] : null;
			for(int w=0; w<data.windows.size(); w++){
				final byte[] text=input.reader.nextLine();
				require(text!=null, "Saved counts lost an original window: "+layout.id+", window="+w);
				p.set(text);
				final Window window=data.windows.get(w);
				require(p.terms()==5+4*layout.models && p.termEquals(layout.id, 0) && p.parseInt(1)==w
					&& same(p, 2, window.input.panel.nameBytes) && p.parseInt(3)==window.input.slot
					&& p.parseLong(4)==window.event, "Count rows must retain exact layout/window/panel/input/event/model order");
				final char[] row=retain ? new char[4*layout.models] : null;
				for(int m=0; m<layout.models; m++){
					if(layout.historical){readHistorical(layout, m, row);}
					else{readCurrent(layout, m, row);}
				}
				if(retain){rows[w]=row;}
			}
			return rows;
		}

		private void readHistorical(final Layout layout, final int model, final char[] row){
			int member=0, consensus=0;
			for(int term=0; term<4; term++){
				final int n=p.parseInt(5+4*model+term);
				require(n>=0 && n<=65535 && n<=layout.limit(term) && (layout.terms==4 || term!=2 || n==0),
					"Historical counts must preserve original bounds and the empty consensus-short slot in long-only rows");
				if(term<2){
					member+=n;
					if(layout.memberShortBins>0){
						require(n<=(term==0 ? layout.memberShortBins : layout.memberLongBins), "Historical member count exceeds its dedicated bins");
					}
				}else{consensus+=n;}
				if(row!=null && (layout.terms==4 || term!=2)){
					row[4*model+(layout.terms==3 && term==3 ? 2 : term)]=(char)n;
				}
			}
			require(member<=layout.memberBins && consensus<=layout.consensusBins, "Historical split counts cannot double-count a stored bin");
		}

		private void readCurrent(final Layout layout, final int model, final char[] row){
			int member=0, consensus=0;
			for(int term=0; term<4; term++){
				final int n=p.parseInt(5+4*model+term);
				require(n>=0 && n<=65535 && (term<layout.terms ? n<=layout.limit(term) : n==0),
					"Current counts need bounded active terms followed by zero padding; no origin-specific hidden term is allowed");
				if(term<layout.terms){
					if(term<layout.memberTerms){member+=n;}else{consensus+=n;}
				}
				if(row!=null){row[4*model+term]=(char)n;}
			}
			require(member<=layout.memberBins && consensus<=layout.consensusBins,
				"Disjoint member partitions and separate consensus totals cannot count a key twice");
		}

		public void finish(){require(next==layouts.length && input.reader.nextLine()==null, "Counts cannot omit layouts or append duplicate rows");}
		@Override
		public void close(){input.close();}
		private final Input input;
		private final Layout[] layouts;
		private final NcrnaDdlSavedData data;
		private final LineParser1 p=new LineParser1('\t');
		private int next;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Static Helpers      ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean powerOfTwo(final int value){return value>0 && (value&(value-1))==0;}

	private static void seconds(final LineParser1 p, final int field){
		final double value=p.parseDouble(field);
		require(Double.isFinite(value) && value>=0, "Construction time must be finite nonnegative seconds");
	}
	private static boolean same(final LineParser1 p, final int field, final byte[] expected){
		p.setBounds(field);
		if(p.b()-p.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	public final Layout[] layouts;
	public final int schemaVersion;
	public final boolean historical;
	public static final String LAYOUT_HEADER="layout_id\tfamily\tmember_mode\tmember_bins\tmember7_bins\tmember10_bins\tconsensus_mode\tconsensus_bins\texponent_bits\tmantissa_bits\tlength_treatment\toverflow_policy\tmodel_count\tmember_key_payload_bytes\tmember_origin_payload_bytes\tconsensus_key_payload_bytes\tconsensus_origin_payload_bytes\tindex_pointer_slots\tindex_populated_slots\tindex_empty_slots\tindex_reference_bytes\tindex_pointer_payload_bytes\tindex_posting_payload_bytes\tindex_origin_payload_bytes\tindex_origin_pointer_slots\tconsensus_build_cpu_seconds\tconsensus_build_elapsed_seconds";
	public static final String COUNT_HEADER="layout_id\twindow\tpanel\tinput_slot\tevent\tmodel_order_member7_member10_consensus7_consensus10";
	public static final String LAYOUT_HEADER_V2="schema_version\tlayout_id\tfamily\tshort_k\tlong_k\tlong_mask\tmember_mode\tmember_bins\tmember_short_bins\tmember_long_bins\tconsensus_mode\tconsensus_bins\tconsensus_bonus\texponent_bits\tmantissa_bits\toverflow_policy\tscore_terms\treference_term\tmodel_count\tmember_key_payload_bytes\tconsensus_key_payload_bytes\tindex_pointer_slots\tindex_populated_slots\tindex_empty_slots\tindex_reference_bytes\tindex_pointer_payload_bytes\tindex_posting_payload_bytes\tconsensus_build_cpu_seconds\tconsensus_build_elapsed_seconds\thash_namespace\tquery_form";
	public static final String COUNT_HEADER_V2="layout_id\twindow\tpanel\tinput_slot\tevent\tmodel_order_active_terms_zero_padded_to_four";
	private static final byte[] V1_LAYOUT_BYTES=LAYOUT_HEADER.getBytes(StandardCharsets.UTF_8);
	private static final byte[] V2_LAYOUT_BYTES=LAYOUT_HEADER_V2.getBytes(StandardCharsets.UTF_8);
}
