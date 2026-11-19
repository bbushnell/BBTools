package prot.cluster;

import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Arrays;

import fileIO.ByteFile;
import parse.LineParser1;
import prot.DigestSuffix;
import prot.ProteinSequence;
import structures.IntHashMap;

/** Hydrates selected seed descriptors from their pinned original assignment table.
 * The supplied population must be the original ordered input for this one source
 * codec. Every ordinal and ID is checked, including nonselected/unassigned rows.
 * Returns a reusable input cache and immutable descriptors; no current family
 * members or model payloads are fabricated, and no intermediate files are written.
 * @author Brian Bushnell, Keqing
 */
public final class SeedMemberRankings{

	public static SeedMemberRankings read(final Path members, final String pin, final ProteinSequence[] population,
			final int[] familyIds, final int[] counts, final int layer){
		require(members!=null && pin!=null && pin.matches("[0-9a-f]{20}") && population!=null &&
			familyIds!=null && counts!=null && familyIds.length>0 && familyIds.length<=population.length &&
			familyIds.length==counts.length && layer>=0,
			"Seed hydration needs pinned assignments, their original population, selected IDs/counts and a measurement layer");
		verify(members,pin);
		final IntHashMap offsets=new IntHashMap(Math.max(16,familyIds.length));
		final MemberRankStatistics[] statistics=new MemberRankStatistics[familyIds.length];
		long selectedMembers=0;
		for(int i=0; i<familyIds.length; i++){
			require(familyIds[i]>=0 && offsets.get(familyIds[i])<0 && counts[i]>0,
				"Initial seed descriptors require unique source-family IDs and positive original member counts");
			offsets.put(familyIds[i],i); statistics[i]=new MemberRankStatistics(); selectedMembers+=counts[i];
		}
		require(selectedMembers<=population.length,"Selected seed memberships cannot exceed the complete original input population");
		final MemberRankInputs inputs=new MemberRankInputs(population);
		final ByteFile file=ByteFile.makeByteFile1(members.toString(),false);
		final LineParser1 fields=new LineParser1('\t');
		int ordinal=0;
		try{
			require(Arrays.equals(file.nextLine(),HEADER),"Seed assignments require the native input_index/input_id/family_id/status header");
			for(byte[] row=file.nextLine(); row!=null; row=file.nextLine(),ordinal++){
				fields.set(row);
				require(fields.terms()==4 && ordinal<population.length && number(fields,0)==ordinal,
					"Seed assignments must contain one ordered row per original input ordinal");
				if(!fields.termEquals(population[ordinal].id,1)){
					throw new IllegalArgumentException("Seed assignment ID differs from the original protein at ordinal "+ordinal);
				}
				if(fields.termEquals("assigned",3)){
					final int at=offsets.get(number(fields,2));
					if(at>=0){
						require(statistics[at].count()<counts[at],"Seed assignment has more retained members than its selected codec count");
						inputs.addTo(statistics[at],ordinal);
					}
				}else{
					require(fields.termEquals("-1",2) && (fields.termEquals("unassigned",3) || fields.termEquals("filtered_length",3)),
						"Nonassigned seed rows require explicit native status and family ID -1");
				}
			}
			require(ordinal==population.length,"Seed assignment table omitted original input rows");
		}finally{if(file.close()){throw new IllegalStateException("Original seed assignment read failed: "+members);}}
		verify(members,pin);
		final MemberRankStatistics.Snapshot[] snapshots=new MemberRankStatistics.Snapshot[statistics.length];
		for(int i=0; i<snapshots.length; i++){
			require(statistics[i].count()==counts[i],"Hydrated member count differs from selected source-family count: "+familyIds[i]);
			snapshots[i]=statistics[i].freeze(layer,null);
		}
		return new SeedMemberRankings(inputs,offsets,snapshots);
	}

	/** Source IDs are mapped before selector renumbering; callers attach the result to the corresponding selected model. */
	public MemberRankStatistics.Snapshot family(final int sourceId){
		final int at=sourceId<0 ? -1 : offsets.get(sourceId);
		require(at>=0,"No selected original family descriptor for source ID "+sourceId);
		return snapshots[at];
	}
	public int families(){return snapshots.length;}

	private SeedMemberRankings(final MemberRankInputs inputs_, final IntHashMap offsets_, final MemberRankStatistics.Snapshot[] snapshots_){
		assert(inputs_!=null && offsets_!=null && snapshots_.length>0) : "Hydrated seeds retain a checked population and nonempty selected descriptor map";
		inputs=inputs_; offsets=offsets_; snapshots=snapshots_;
	}
	private static void verify(final Path path, final String pin){
		require(pin.equals(DigestSuffix.file(path.toString())),"Original seed membership fingerprint differs: "+path);
	}
	/** Checked decimal parsing prevents a wrapped input or family ordinal from aliasing a real row. */
	private static int number(final LineParser1 fields, final int term){
		fields.setBounds(term);
		final byte[] row=fields.line(); final int end=fields.b(); int value=0;
		require(fields.a()<end,"Empty seed ordinal or family ID");
		for(int i=fields.a(); i<end; i++){
			final int digit=row[i]-'0';
			require(digit>=0 && digit<=9 && value<=(Integer.MAX_VALUE-digit)/10,"Invalid or overflowing seed ordinal or family ID");
			value=value*10+digit;
		}
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	public final MemberRankInputs inputs;
	private final IntHashMap offsets;
	private final MemberRankStatistics.Snapshot[] snapshots;
	private static final byte[] HEADER="input_index\tinput_id\tfamily_id\tstatus".getBytes(StandardCharsets.US_ASCII);
}
