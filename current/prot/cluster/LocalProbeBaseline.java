package prot.cluster;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import parse.LineParser1;
import prot.ProteinSequence;
import structures.IntHashMap;
import structures.IntList;

/** Reconstructs the frozen selected-model baseline from one completed assignment.
 * Model payloads and original source counts come from the unchanged codec;
 * CURRENT memberships and the residual pool come solely from the ordered
 * assignment table. Selected zero-member models survive for neighbor evidence.
 * Neither source counts nor duplicate FASTA names invent original input rows.
 * All returned models/arrays are owned by this baseline and then borrowed
 * read-only by proposal and scoring workers.
 * @author Brian Bushnell, Keqing
 */
public final class LocalProbeBaseline{

	/** The driver verifies the input/codec/configuration binding before loading.
	 * expectedSelected is the actual declared roster size, not a fresh top-N rule:
	 * this loader consumes the setup's selection and never repeats its ranking.
	 */
	public static LocalProbeBaseline load(final ClusterSketchCodec.Loaded source,
			final String selectedPath,final String assignmentsPath,final ProteinSequence[] sequences,
			final int minimumLength,final int expectedSelected){
		return load(source,selectedPath,assignmentsPath,sequences,minimumLength,expectedSelected,false,-1,null);
	}

	/** Selected-original loader that remeasures current nonempty assignments and preserves
	 * serialized last-nonempty descriptors for empty selected templates.
	 */
	public static LocalProbeBaseline load(final ClusterSketchCodec.Loaded source,
			final String selectedPath,final String assignmentsPath,final ProteinSequence[] sequences,
			final int minimumLength,final int expectedSelected,final MemberRankInputs memberInputs){
		require(memberInputs!=null && memberInputs.belongsTo(sequences),
			"Ranked selected-baseline hydration requires the shared original-population cache");
		return load(source,selectedPath,assignmentsPath,sequences,minimumLength,expectedSelected,false,-1,memberInputs);
	}

	/** Imports every current model without selection or renumbering, including empty templates.
	 * The saved cursor prevents reuse of IDs retired by previous layers; codec counts must
	 * equal the supplied current memberships, rather than pre-assignment source counts.
	 */
	public static LocalProbeBaseline loadCurrent(final ClusterSketchCodec.Loaded source,
			final String assignmentsPath,final ProteinSequence[] sequences,final int minimumLength,
			final int nextFamilyId){
		require(source!=null && source.families!=null,"Current state requires a complete model codec");
		return load(source,null,assignmentsPath,sequences,minimumLength,source.families.length,true,nextFamilyId,null);
	}

	/** Current-state loader that verifies serialized nonempty descriptors exactly and
	 * keeps serialized last-nonempty descriptors on empty families.
	 */
	public static LocalProbeBaseline loadCurrent(final ClusterSketchCodec.Loaded source,
			final String assignmentsPath,final ProteinSequence[] sequences,final int minimumLength,
			final int nextFamilyId,final MemberRankInputs memberInputs){
		require(source!=null && source.families!=null,"Current state requires a complete model codec");
		require(memberInputs!=null && memberInputs.belongsTo(sequences),
			"Ranked current-state hydration requires the shared original-population cache");
		return load(source,null,assignmentsPath,sequences,minimumLength,source.families.length,true,nextFamilyId,memberInputs);
	}

	/** Takes the completed assignment and normalizes member order as the old ordinal TSV reader did.
	 * The producer has joined all workers; no concurrent reader may borrow these models yet.
	 * LiveProbeState performs the independent unique-owner partition validation.
	 */
	static LocalProbeBaseline fromAssignment(final GreedyClusterer.Result result,
			final ClusterSketchCodec.Header header,final prot.ProteinSequence[] sequences,
			final int minimumLength,final int nextFamilyId){
		require(result!=null && header!=null && sequences!=null && result.inputs()==sequences.length &&
			minimumLength>0 && header.families==result.families.length,"Resident assignment requires matching population and model census");
		final Family[] selected=result.families;
		final ArrayList<Family> active=new ArrayList<Family>();
		final IntHashMap index=new IntHashMap(Math.max(16,selected.length));
		final int[] counts=new int[selected.length];
		long assigned=0;int maximum=-1;
		for(int i=0;i<selected.length;i++){
			final Family family=selected[i];
			require(family!=null && index.get(family.familyId)<0,"Resident assignment needs unique model IDs");
			index.put(family.familyId,i);counts[i]=family.size();assigned+=family.size();maximum=Math.max(maximum,family.familyId);
			final int[] ordinals=family.memberIndexes();Arrays.sort(ordinals);
			family.clearMembers();
			for(int j=0;j<ordinals.length;j++){
				final int ordinal=ordinals[j];
				require(ordinal>=0 && ordinal<sequences.length && (j==0 || ordinal>ordinals[j-1]) &&
					result.familyFor(ordinal)==family.familyId,"Resident memberships must preserve the completed unique-owner assignment");
				family.addMember(ordinal);
			}
			if(family.size()>0){active.add(family);}
		}
		require(nextFamilyId>maximum && nextFamilyId<Integer.MAX_VALUE-8,"Resident fresh-ID cursor must exceed every retained model ID");
		final IntList residual=new IntList();int filtered=0;
		for(int i=0;i<sequences.length;i++){
			if(sequences[i].length()<minimumLength){require(result.familyFor(i)<0,"Filtered input cannot have an assigned family");filtered++;}
			else if(result.familyFor(i)<0){residual.add(i);}
		}
		require(assigned==result.assigned && assigned+filtered+residual.size==sequences.length,
			"Resident members, residuals and filtered inputs must conserve the assignment population");
		return new LocalProbeBaseline(selected,active.toArray(new Family[active.size()]),counts,residual.toArray(),
			index,nextFamilyId,Math.toIntExact(assigned),filtered,header);
	}

	private static LocalProbeBaseline load(final ClusterSketchCodec.Loaded source,
			final String selectedPath,final String assignmentsPath,final ProteinSequence[] sequences,
			final int minimumLength,final int expectedSelected,final boolean current,final int nextFamilyId,
			final MemberRankInputs memberInputs){
		require(source!=null && source.header!=null && source.families!=null && sequences!=null &&
			(current || selectedPath!=null) && assignmentsPath!=null && minimumLength>0 && expectedSelected>0 &&
			expectedSelected<=source.families.length,"Selected baseline requires a source codec, complete inputs and an explicit valid selection size");
		require(source.header.families==source.families.length,"Source codec header and roster counts differ");
		final IntHashMap sourceIndex=new IntHashMap(); int maxId=-1;
		for(int i=0;i<source.families.length;i++){
			final ClusterSketchCodec.Entry entry=source.families[i];
			require(entry!=null && entry.familyId>=0 && entry.familyId<Integer.MAX_VALUE-8 &&
				(current ? entry.members>=0 : entry.members>0) && sourceIndex.get(entry.familyId)<0,
				"Original codec needs unique bounded IDs with room for a fresh fission child");
			sourceIndex.put(entry.familyId,i); maxId=Math.max(maxId,entry.familyId);
		}
		require(!current || nextFamilyId>maxId && nextFamilyId<Integer.MAX_VALUE-8,
			"Current state requires its saved next-family cursor beyond every serialized ID, with room for new children");
		final Family[] selected=new Family[expectedSelected];
		final ClusterSketchCodec.Entry[] selectedEntries=memberInputs==null ? null : new ClusterSketchCodec.Entry[expectedSelected];
		final int[] sourceCounts=new int[expectedSelected];
		final IntHashMap selectedIndex=new IntHashMap();
		if(current){
			assert(expectedSelected==source.families.length) : "Current-state loading preserves the complete serialized roster";
			for(int i=0;i<selected.length;i++){
				final ClusterSketchCodec.Entry entry=source.families[i];
				final Family family=Family.fromSeed(entry,true);
				family.clearMembers();
				selected[i]=family; sourceCounts[i]=entry.members;
				if(selectedEntries!=null){selectedEntries[i]=entry;}
				selectedIndex.put(entry.familyId,i);
			}
		}else{
			final ByteFile roster=ByteFile.makeByteFile(selectedPath,false);
			int selectedCount=0;
			try{
				require(Arrays.equals(roster.nextLine(),SELECTED_HEADER),"Unexpected selected-model roster header: "+selectedPath);
				final LineParser1 fields=new LineParser1('\t');
				for(byte[] line=roster.nextLine();line!=null;line=roster.nextLine()){
					fields.set(line);
					require(fields.terms()==2 && selectedCount<selected.length,"Selected-model roster must have exactly the declared number of two-column rows");
					final int id=nonnegative(fields,0),count=nonnegative(fields,1),at=sourceIndex.get(id);
					requireAt(at>=0 && selectedIndex.get(id)<0,"Selected stable ID must occur once and exist in the original codec",id);
					final ClusterSketchCodec.Entry entry=source.families[at];
					requireAt(count==entry.members,"Selected source count differs from immutable codec for family",id);
					final Family family=Family.fromSeed(entry);
					// Remove the source-count fallback before attaching current memberships.
					family.clearMembers();
					selected[selectedCount]=family; sourceCounts[selectedCount]=count;
					if(selectedEntries!=null){selectedEntries[selectedCount]=entry;}
					selectedIndex.put(id,selectedCount++);
				}
				require(selectedCount==selected.length,"Selected-model roster omitted declared references");
			}finally{close(roster,selectedPath);}
		}

		final IntList residual=new IntList();
		int ordinal=0,assigned=0,filtered=0;
		final ByteFile assignments=ByteFile.makeByteFile(assignmentsPath,false);
		try{
			require(Arrays.equals(assignments.nextLine(),ASSIGNMENT_HEADER),"Unexpected one-pass assignment header: "+assignmentsPath);
			final LineParser1 fields=new LineParser1('\t');
			for(byte[] line=assignments.nextLine();line!=null;line=assignments.nextLine(),ordinal++){
				fields.set(line);
				requireAt(fields.terms()==4 && ordinal<sequences.length && nonnegative(fields,0)==ordinal,
					"One-pass assignments require exactly one ordered row per original input; expected ordinal",ordinal);
				final ProteinSequence sequence=sequences[ordinal];
				requireAt(sequence!=null && fields.termEquals(sequence.id,1),"Assignment identity differs from original input at ordinal",ordinal);
				final boolean shortInput=sequence.length()<minimumLength;
				if(fields.termEquals("assigned",3)){
					final int id=nonnegative(fields,2),at=selectedIndex.get(id);
					requireAt(!shortInput && at>=0,"Assigned row must be length-eligible and target a selected stable ID; ordinal",ordinal);
					selected[at].addMember(ordinal); assigned++;
				}else{
					requireAt(fields.termEquals("-1",2),"Nonassigned rows require family ID -1; ordinal",ordinal);
					if(fields.termEquals("unassigned",3)){
						requireAt(!shortInput,"A short input belongs to filtered_length, not the residual pool; ordinal",ordinal);
						residual.add(ordinal);
					}else{
						requireAt(shortInput && fields.termEquals("filtered_length",3),"Unknown or inconsistent assignment status; ordinal",ordinal);
						filtered++;
					}
				}
			}
			require(ordinal==sequences.length,"One-pass assignment table omitted original input rows");
		}finally{close(assignments,assignmentsPath);}
		final ArrayList<Family> active=new ArrayList<Family>();
		long counted=0;
		for(int i=0;i<selected.length;i++){
			final Family family=selected[i];
			requireAt(!current || family.size()==sourceCounts[i],
				"Current codec count differs from its assignment table for family",family.familyId);
			counted+=family.size(); if(family.size()>0){active.add(family);}
		}
		final int memberRankingLayer=memberInputs==null ? -1 : attachMemberRankings(selected,selectedEntries,memberInputs,current);
		require(counted==assigned && (long)assigned+filtered+residual.size==sequences.length,
			"Selected memberships, residual pool and length-filtered rows must partition the original inputs");
		return new LocalProbeBaseline(selected,active.toArray(new Family[active.size()]),sourceCounts,
			residual.toArray(),selectedIndex,current ? nextFamilyId : maxId+1,assigned,filtered,source.header,memberRankingLayer);
	}

	private static int attachMemberRankings(final Family[] selected,final ClusterSketchCodec.Entry[] entries,
			final MemberRankInputs inputs,final boolean current){
		int layer=0;
		for(int i=0;i<selected.length;i++){
			final ClusterSketchCodec.Entry entry=entries[i];
			requireAt(entry!=null && entry.memberRanking!=null,
				"Ranked hydration requires serialized member descriptors for family",selected[i].familyId);
			requireAt(entry.memberRanking.layer<Integer.MAX_VALUE,
				"Ranked hydration cannot advance beyond the maximum provenance layer for family",selected[i].familyId);
			layer=Math.max(layer,Math.addExact(entry.memberRanking.layer,1));
		}
		for(int i=0;i<selected.length;i++){
			final Family family=selected[i];
			final MemberRankStatistics.Snapshot serialized=entries[i].memberRanking;
			family.enableMemberRanking(inputs,current ? serialized.layer : layer,serialized);
			if(current && family.size()>0){
				requireAt(same(family.memberRankingSnapshot(),serialized),
					"Current ranked descriptor differs from its assignment table for family",family.familyId);
			}
		}
		return layer;
	}

	private static boolean same(final MemberRankStatistics.Snapshot a,final MemberRankStatistics.Snapshot b){
		return a!=null && b!=null && a.members==b.members && a.layer==b.layer &&
			Double.doubleToRawLongBits(a.medianLength)==Double.doubleToRawLongBits(b.medianLength) &&
			Double.doubleToRawLongBits(a.meanKd)==Double.doubleToRawLongBits(b.meanKd);
	}

	private LocalProbeBaseline(Family[] selected_,Family[] active_,int[] sourceCounts_,int[] residual_,
			IntHashMap selectedIndex_,int firstNewId_,int assigned_,int filtered_,ClusterSketchCodec.Header header_){
		this(selected_,active_,sourceCounts_,residual_,selectedIndex_,firstNewId_,assigned_,filtered_,header_,-1);
	}

	private LocalProbeBaseline(Family[] selected_,Family[] active_,int[] sourceCounts_,int[] residual_,
			IntHashMap selectedIndex_,int firstNewId_,int assigned_,int filtered_,ClusterSketchCodec.Header header_,
			int memberRankingLayer_){
		selected=selected_; active=active_; sourceCounts=sourceCounts_; residualOrdinals=residual_;
		selectedIndex=selectedIndex_; firstNewId=firstNewId_; assigned=assigned_; filtered=filtered_; sourceHeader=header_;
		memberRankingLayer=memberRankingLayer_;
	}

	/** Index in selected order; -1 means the stable ID was not selected. */
	public int selectedIndex(final int stableId){return stableId<0 ? -1 : selectedIndex.get(stableId);}

	/** Checked parsing prevents overflow from aliasing a different selected family. */
	private static int nonnegative(final LineParser1 fields,final int term){
		fields.setBounds(term); final byte[] line=fields.line();
		final int end=fields.b(); int value=0;
		requireAt(fields.a()<end,"Empty numeric baseline field",term);
		for(int i=fields.a();i<end;i++){
			final int digit=line[i]-'0';
			requireAt(digit>=0 && digit<=9 && value<=(Integer.MAX_VALUE-digit)/10,
				"Invalid or overflowing numeric baseline field",term);
			value=value*10+digit;
		}
		return value;
	}
	private static void close(final ByteFile file,final String path){
		if(file.close()){throw new IllegalStateException("Selected-baseline input failed: "+path);}
	}
	private static void require(final boolean valid,final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
	private static void requireAt(final boolean valid,final String message,final int position){
		if(!valid){throw new IllegalArgumentException(message+": "+position);}
	}
	public final Family[] selected,active;
	public final int[] sourceCounts,residualOrdinals;
	public final int firstNewId,assigned,filtered,memberRankingLayer;
	public final ClusterSketchCodec.Header sourceHeader;
	private final IntHashMap selectedIndex;
	private static final byte[] SELECTED_HEADER="family_id\tsource_members".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] ASSIGNMENT_HEADER="input_index\tinput_id\tfamily_id\tstatus".getBytes(StandardCharsets.US_ASCII);
}
