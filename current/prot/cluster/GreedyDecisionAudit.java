package prot.cluster;

import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prot.ConsensusScalars;
import prot.ProteinSequence;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Optional serial greedy-pass observations.  Old-reference family state is
 * updated once per actual admission; candidate rows observe existing screens.
 * @author Yelan
 */
public final class GreedyDecisionAudit implements AutoCloseable{

	public GreedyDecisionAudit(String prefix,String oldMembers,String lostMap,
			ProteinSequence[] sequences_,int shortlist_,int width_,boolean overwrite){
		if(prefix==null || prefix.length()==0 || oldMembers==null || lostMap==null ||
			sequences_==null || shortlist_<1 || shortlist_>500 || width_!=1024){
			throw new IllegalArgumentException("Greedy audit needs prefix, complete old members, lost map, inputs, 1-500 shortlist and 1024-wide sketches");
		}
		sequences=sequences_; shortlist=shortlist_; width=width_;
		counts=new long[5][width+1][3];
		selected=new boolean[sequences.length]; seen=new boolean[sequences.length];
		admitted=new boolean[sequences.length]; winnerByInput=new int[sequences.length];
		Arrays.fill(winnerByInput,Integer.MIN_VALUE);
		trackedOrdinalByInput=new int[sequences.length]; Arrays.fill(trackedOrdinalByInput,-1);
		oldFamilyByInput=new int[sequences.length]; Arrays.fill(oldFamilyByInput,-1);
		expectedOldFamilyByInput=new int[sequences.length]; Arrays.fill(expectedOldFamilyByInput,-1);
		readLostMap(lostMap);
		states=new OldFamilyState[oldFamilies];
		for(int i=0; i<states.length; i++){states[i]=new OldFamilyState();}
		readOldMembers(oldMembers);
		final String[] paths=outputPaths(prefix);
		try{
			queries=open(paths[0],overwrite); candidates=open(paths[1],overwrite); screens=open(paths[2],overwrite);
			queries.println(QUERY_HEADER);
			candidates.println("input_index\trank\tfamily_id\tcurrent_members\told_family_id\told_family_shared_members\thits\trequired_hits\toccupied_bins\thydro_available\thydro_delta\tscreen\talignment_attempted\tgate_reason\tscore\tquery_span\treference_span\tmetrics_available\tidentities\tcolumns\tdenominator_self_score\tidentity\tscore_ratio");
		}catch(RuntimeException failure){
			try{closeStreams();}catch(RuntimeException cleanup){failure.addSuppressed(cleanup);}
			throw failure;
		}
	}

	public static String[] outputPaths(String prefix){
		if(prefix==null || prefix.length()==0){throw new IllegalArgumentException("Empty greedy-audit prefix");}
		return new String[]{prefix+".greedy_queries.tsv",prefix+".greedy_candidates.tsv",
			prefix+".greedy_screens.tsv"};
	}

	public int trackedOldFamilyCount(){return oldFamilies;}
	public boolean selected(int input){return selected[input];}
	public int oldFamily(int input){return oldFamilyByInput[input];}
	public int historicalFamily(){
		require(active,"Greedy historical family is only defined during an active query");
		return currentHistoricalFamily;
	}
	public void validateSerial(boolean concurrent){
		if(concurrent){throw new IllegalArgumentException("Greedy decision audit v2 is serial-only");}
	}

	/** Ranking has completed.  The historical family is the live family with
	 * the most already-admitted members from this query's old reference family. */
	public void beginQuery(int input,String id,int count,CandidateFilter.Context query){
		require(!active && input>=0 && input<sequences.length && id!=null &&
			sequences[input].id.equals(id) && count>=0 && count<=shortlist && query!=null,
			"Greedy audit query must describe one ranked input");
		require(!seen[input],"Greedy audit saw duplicate input: "+input);
		seen[input]=true; currentInput=input; currentOrdinal=trackedOrdinalByInput[input];
		currentOldFamily=oldFamilyByInput[input]; currentSelected=selected[input];
		currentHistoricalFamily=-1; currentHistoricalShared=0; currentHistoricalTies=0;
		currentHistoricalRank=0; currentHistoricalHits=-1;
		currentHistoricalRequired=-1; currentHistoricalOccupied=-1; currentHistoricalScreen=-1;
		currentTargetHits=currentTargetRequired=currentTargetOccupied=currentTargetRank=currentLiveFamilies=-1;
		currentHistoricalState=false;
		queryRow.clear(); candidateRows.clear();
		if(currentSelected){
			require(currentOrdinal>=0,"Selected greedy-audit query must have a tracked old family");
			final OldFamilyState state=states[currentOrdinal];
			if(state.bestCount>0){
				currentHistoricalFamily=state.bestFamily;
				currentHistoricalShared=state.bestCount;
				currentHistoricalTies=state.bestTies;
			}
			queryRow.append(input).tab().append(id).tab().append(currentOldFamily).tab()
				.append(oldFamilySizes[currentOrdinal]).tab().append(count);
			selectedSeen++;
		}
		active=true;
	}

	/** Record full-index target state for the live historical family, whether or not
	 * it appears in the emitted shortlist. */
	public void historicalState(int hits,int requiredHits,int occupiedBins,int censoredRank,int liveFamilies){
		require(active && currentSelected && currentHistoricalFamily>=0 && !currentHistoricalState &&
			hits>=0 && requiredHits>=0 && requiredHits<counts.length && occupiedBins>=0 &&
			occupiedBins<=width && liveFamilies>0 && ((hits<requiredHits)==(censoredRank==0)) &&
			(hits<requiredHits || (censoredRank>=1 && censoredRank<=501 &&
				censoredRank<=liveFamilies)),
			"Greedy historical target state must be available once for a present old-family proxy");
		currentTargetHits=hits; currentTargetRequired=requiredHits;
		currentTargetOccupied=occupiedBins; currentTargetRank=censoredRank;
		currentLiveFamilies=liveFamilies; currentHistoricalState=true;
	}

	/** Record the actual greedy screen and optional gate result without invoking any scorer. */
	public void candidate(int rank,Family family,int hits,int requiredHits,int occupiedBins,
			CandidateFilter.Context query,int screen,TailAlignmentDecision decision){
		require(active && rank>=1 && rank<=shortlist && family!=null && query!=null && hits>=0 &&
				requiredHits>=0 && requiredHits<counts.length && occupiedBins>=0 && occupiedBins<=width &&
			screen>=CandidateFilter.PASS && screen<=CandidateFilter.HYDRO &&
			((screen==CandidateFilter.LOW_HITS)==(hits<requiredHits)) &&
			((screen==CandidateFilter.PASS)==(decision!=null)),
			"Greedy candidate audit must follow an actual ranked screen");
		counts[requiredHits][occupiedBins][screen]++;
		if(!currentSelected){return;}
		if(family.familyId==currentHistoricalFamily){
			require(currentHistoricalRank==0,"Historical family appeared twice in the same ranked shortlist");
			currentHistoricalRank=rank; currentHistoricalHits=hits;
			currentHistoricalRequired=requiredHits; currentHistoricalOccupied=occupiedBins;
			currentHistoricalScreen=screen;
		}
		final boolean hydroAvailable=query.defined() && hydroDefined(family);
		final double delta=hydroAvailable ?
			Math.abs(query.hydrophobicity()-family.scalarsUnsafe()[ConsensusScalars.HYDRO]) : 0;
		finite(delta);
		candidateRows.append(currentInput).tab().append(rank).tab().append(family.familyId).tab()
			.append(family.size()).tab().append(currentOldFamily).tab()
			.append(currentOrdinal<0 ? 0 : states[currentOrdinal].shared(family.familyId)).tab()
			.append(hits).tab().append(requiredHits).tab().append(occupiedBins).tab()
			.append(hydroAvailable ? 1 : 0).tab();
		value(candidateRows,hydroAvailable,delta);
		candidateRows.tab().append(SCREENS[screen]).tab().append(decision==null ? 0 : 1).tab();
		if(decision==null){candidateRows.append("not_aligned\tNA\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA");}
		else{
			candidateRows.append(decision.reason.name()).tab();
			final boolean hasTrace=decision.reason!=TailAlignmentDecision.Reason.NO_ALIGNMENT;
			if(hasTrace){candidateRows.append(decision.score).tab().append(decision.querySpan).tab().append(decision.referenceSpan);}
			else{candidateRows.append("NA\tNA\tNA");}
			final boolean metrics=decision.columns>0 && decision.denominatorSelfScore>0;
			candidateRows.tab().append(metrics ? 1 : 0).tab();
			if(metrics){
				finite(decision.identity); finite(decision.scoreRatio);
				candidateRows.append(decision.identities).tab().append(decision.columns).tab()
					.append(decision.denominatorSelfScore).tab().append(decision.identity,12).tab().append(decision.scoreRatio,12);
			}else{candidateRows.append("NA\tNA\tNA\tNA\tNA");}
		}
		candidateRows.nl();
	}

	public void finishQuery(int winner,int rank){
		require(active && rank>=0 && rank<=shortlist && ((winner>=0)==(rank>0)),
			"Greedy query winner must match its rank");
		if(currentSelected){
			queryRow.tab().append(winner).tab().append(rank).tab()
				.append(currentHistoricalFamily>=0 ? 1 : 0).tab()
				.append(currentHistoricalFamily).tab()
				.append(currentHistoricalShared).tab().append(currentHistoricalTies).tab()
				.append(currentHistoricalRank).tab();
			value(queryRow,currentHistoricalRank>0,currentHistoricalHits); queryRow.tab();
			value(queryRow,currentHistoricalRank>0,currentHistoricalRequired); queryRow.tab();
			value(queryRow,currentHistoricalRank>0,currentHistoricalOccupied); queryRow.tab();
			if(currentHistoricalRank>0){queryRow.append(SCREENS[currentHistoricalScreen]);}
			else{queryRow.append("NA");}
			appendHistoricalState(queryRow);
			queryRow.nl(); queries.print(queryRow); candidates.print(candidateRows);
		}
		winnerByInput[currentInput]=winner;
		active=false;
	}

	public void admit(int input,int family){
		require(input>=0 && input<seen.length && seen[input] && !admitted[input] &&
			winnerByInput[input]!=Integer.MIN_VALUE && family>=0 &&
			(winnerByInput[input]<0 || winnerByInput[input]==family),
			"Greedy audit admission must follow exactly one matching query: "+input);
		admitted[input]=true;
		final int ordinal=trackedOrdinalByInput[input];
		if(ordinal>=0){states[ordinal].admit(family);}
	}

	private void readLostMap(String path){
		final ByteFile input=ByteFile.makeByteFile1(path,false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			byte[] line=input.nextLine();
			while(line!=null && line.length>0 && line[0]=='#'){line=input.nextLine();}
			require(line!=null && AssignmentDecisionAudit.MAP_HEADER.equals(ascii(line)),"Unexpected lost-map header");
			for(line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==10,"Lost map needs ten fields");
				final int index=lp.parseInt(0), oldFamily=lp.parseInt(2), oldSize=lp.parseInt(3);
				require(index>=0 && index<sequences.length && oldFamily>=0 && oldSize>0,
					"Lost map fields are out of range");
				require(!selected[index] && lp.parseString(1).equals(sequences[index].id),
					"Lost map inputs must be unique and match the input FASTA");
				int ordinal=oldFamilyToOrdinal.get(oldFamily);
				if(ordinal<0){
					ordinal=oldFamilies++;
					oldFamilyToOrdinal.set(oldFamily,ordinal);
					if(ordinal>=oldFamilySizes.length){
						oldFamilySizes=Arrays.copyOf(oldFamilySizes,oldFamilySizes.length*2);
					}
					oldFamilySizes[ordinal]=oldSize;
				}else{require(oldFamilySizes[ordinal]==oldSize,"Old-family size disagrees for "+oldFamily);}
				selected[index]=true; expectedOldFamilyByInput[index]=oldFamily; selectedCount++;
			}
		}finally{require(!input.close(),"Could not read lost map: "+path);}
		oldFamilySizes=Arrays.copyOf(oldFamilySizes,oldFamilies);
	}

	private void readOldMembers(String path){
		final ByteFile input=ByteFile.makeByteFile1(path,false);
		final LineParser1 lp=new LineParser1('\t');
		final int[] observedOldFamilySizes=new int[oldFamilies];
		try{
			byte[] line=input.nextLine();
			require(line!=null && MEMBERS_HEADER.equals(ascii(line)),"Unexpected old-members header");
			int previous=-1;
			for(line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==4,"Old members need four fields");
				final int index=lp.parseInt(0), oldFamily=lp.parseInt(2);
				require(index==previous+1 && index<sequences.length &&
					lp.parseString(1).equals(sequences[index].id),"Old members must be ordered and match the input FASTA");
				previous=index;
				if(lp.termEquals("assigned",3)){
					require(oldFamily>=0,"Assigned old member needs a family ID");
					require(!selected[index] || oldFamily==expectedOldFamilyByInput[index],
						"Selected lost-query old family disagrees with the complete old map: "+index);
					final int ordinal=oldFamilyToOrdinal.get(oldFamily);
					if(ordinal>=0){
						oldFamilyByInput[index]=oldFamily; trackedOrdinalByInput[index]=ordinal;
						observedOldFamilySizes[ordinal]++;
					}
				}else{require(oldFamily<0,"Unassigned old member must use a negative family ID");}
			}
			require(previous==sequences.length-1,"Complete old members must cover every input");
		}finally{require(!input.close(),"Could not read complete old members: "+path);}
		for(int i=0; i<oldFamilies; i++){require(observedOldFamilySizes[i]==oldFamilySizes[i],
			"Lost-map old-family size disagrees with the complete old map");}
		for(int i=0; i<selected.length; i++){require(!selected[i] || trackedOrdinalByInput[i]>=0,
			"Every selected lost query must be assigned in the complete old map: "+i);}
	}

	private void writeScreens(){
		final ByteBuilder bb=new ByteBuilder("required_hits\toccupied_bins\tranked_candidates\thit_filtered\thydro_filtered\tpassed_screen\n");
		for(int required=0; required<counts.length; required++){
			for(int occupied=0; occupied<counts[required].length; occupied++){
				final long pass=counts[required][occupied][CandidateFilter.PASS];
				final long hit=counts[required][occupied][CandidateFilter.LOW_HITS];
				final long hydro=counts[required][occupied][CandidateFilter.HYDRO];
				if(pass+hit+hydro>0){
					bb.append(required).tab().append(occupied).tab().append(pass+hit+hydro).tab()
						.append(hit).tab().append(hydro).tab().append(pass).nl();
				}
			}
		}
		screens.print(bb);
	}

	private void appendHistoricalState(ByteBuilder bb){
		bb.tab();
		if(currentHistoricalFamily<0){
			value(bb,false,currentTargetHits); bb.tab();
			value(bb,false,currentTargetRequired); bb.tab();
			value(bb,false,currentTargetOccupied); bb.tab();
			value(bb,false,currentTargetRank); bb.tab();
			value(bb,false,currentLiveFamilies); bb.tab().append(NO_LIVE);
			return;
		}
		require(currentHistoricalState,"Present historical family needs full-index target state");
		if(currentHistoricalRank>0){
			require(currentHistoricalHits==currentTargetHits &&
				currentHistoricalRequired==currentTargetRequired &&
				currentHistoricalOccupied==currentTargetOccupied,
				"Historical target state disagrees with its emitted shortlist row");
		}
		if(currentTargetRank>0){
			require(currentTargetRank<=shortlist ? currentHistoricalRank==currentTargetRank :
				currentHistoricalRank==0,
				"Historical target rank disagrees with the emitted shortlist");
		}
		bb.append(currentTargetHits).tab().append(currentTargetRequired).tab()
			.append(currentTargetOccupied).tab();
		value(bb,currentTargetRank>0,currentTargetRank);
		bb.tab().append(currentLiveFamilies).tab().append(rankBucket(currentTargetRank));
	}

	private static boolean hydroDefined(Family family){return (family.scalarMask()&(1<<ConsensusScalars.HYDRO))!=0;}
	private static void finite(double value){require(Double.isFinite(value),"Greedy audit measurements must be finite");}
	private static void value(ByteBuilder bb,boolean available,double number){if(available){bb.append(number,12);}else{bb.append("NA");}}
	private static void value(ByteBuilder bb,boolean available,int number){if(available){bb.append(number);}else{bb.append("NA");}}
	private static String ascii(byte[] line){return new String(line,java.nio.charset.StandardCharsets.US_ASCII);}
	private static String rankBucket(int rank){
		if(rank==0){return BELOW_BAR;}
		if(rank<=50){return "rank_1_50";}
		if(rank<=100){return "rank_51_100";}
		if(rank<=500){return "rank_101_500";}
		return "rank_gt_500";
	}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalStateException(message);}}
	private static ByteStreamWriter open(String path,boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false); writer.start(); return writer;
	}
	private void closeStreams(){
		boolean error=false;
		RuntimeException failure=null;
		final ByteStreamWriter q=queries,c=candidates,s=screens;
		queries=null; candidates=null; screens=null;
		if(q!=null){
			try{error|=q.poisonAndWait();}catch(RuntimeException e){failure=merge(failure,e);}
		}
		if(c!=null){
			try{error|=c.poisonAndWait();}catch(RuntimeException e){failure=merge(failure,e);}
		}
		if(s!=null){
			try{error|=s.poisonAndWait();}catch(RuntimeException e){failure=merge(failure,e);}
		}
		if(error){failure=merge(failure,new IllegalStateException("Greedy audit output failed"));}
		if(failure!=null){throw failure;}
	}
	@Override
	public void close(){
		RuntimeException failure=null;
		try{
			require(!active && selectedSeen==selectedCount,
				"Greedy audit saw "+selectedSeen+" of "+selectedCount+" selected lost queries");
			writeScreens();
		}catch(RuntimeException e){failure=e;}
		try{closeStreams();}catch(RuntimeException e){failure=merge(failure,e);}
		if(failure!=null){throw failure;}
	}
	private static RuntimeException merge(RuntimeException failure,RuntimeException next){
		if(failure==null){return next;}
		failure.addSuppressed(next);
		return failure;
	}

	private final ProteinSequence[] sequences;
	private final int shortlist,width;
	private final boolean[] selected,seen,admitted;
	private final int[] winnerByInput;
	private final int[] trackedOrdinalByInput,oldFamilyByInput,expectedOldFamilyByInput;
	private int[] oldFamilySizes=new int[16];
	private OldFamilyState[] states;
	private final IntHashMap oldFamilyToOrdinal=new IntHashMap();
	private final long[][][] counts;
	private ByteStreamWriter queries,candidates,screens;
	private final ByteBuilder queryRow=new ByteBuilder(),candidateRows=new ByteBuilder();
	private int oldFamilies,selectedCount,selectedSeen,currentInput,currentOrdinal,currentOldFamily;
	private int currentHistoricalFamily,currentHistoricalShared,currentHistoricalTies;
	private int currentHistoricalRank,currentHistoricalHits,currentHistoricalRequired,currentHistoricalOccupied,currentHistoricalScreen;
	private int currentTargetHits,currentTargetRequired,currentTargetOccupied,currentTargetRank,currentLiveFamilies;
	private boolean active,currentSelected,currentHistoricalState;
	private static final String[] SCREENS={"pass","hits","hydro"};
	private static final String NO_LIVE="no_live_historical_family";
	private static final String BELOW_BAR="below_hit_bar";
	private static final String MEMBERS_HEADER="input_index\tinput_id\tfamily_id\tstatus";
	private static final String QUERY_HEADER="input_index\tinput_id\told_family_id\told_family_size\tcandidates\twinner_family_id\twinner_rank\thistorical_family_present\thistorical_family_id\thistorical_family_shared_members\thistorical_family_tied_families\thistorical_family_rank\thistorical_family_hits\thistorical_family_required_hits\thistorical_family_occupied_bins\thistorical_family_screen\thistorical_target_hits\thistorical_target_required_hits\thistorical_target_occupied_bins\thistorical_target_censored_rank\thistorical_live_families\thistorical_target_rank_bucket";

	private static final class OldFamilyState{
		void admit(int liveFamily){
			final int oldCount=Tools.max(counts.get(liveFamily),0);
			final int newCount=counts.increment(liveFamily);
			if(newCount>bestCount){bestFamily=liveFamily; bestCount=newCount; bestTies=1;}
			else if(newCount==bestCount && oldCount<bestCount){
				bestTies++;
				if(liveFamily<bestFamily){bestFamily=liveFamily;}
			}
		}
		int shared(int liveFamily){return Tools.max(counts.get(liveFamily),0);}
		final IntHashMap counts=new IntHashMap();
		int bestFamily=-1,bestCount,bestTies;
	}
}
