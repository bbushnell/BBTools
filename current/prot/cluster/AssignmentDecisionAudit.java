package prot.cluster;

import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prot.ConsensusScalars;
import prot.ProteinSequence;
import structures.ByteBuilder;

/** Optional serial final-pass observations. These streams never select a candidate.
 * Raw TSVs are packaged by the run controller, without changing global compression
 * settings or ordinary scientific output bytes. Hydro bits preserve exact values;
 * identity and ratio can be recovered exactly from the recorded integer counts.
 * @author Keqing
 */
public final class AssignmentDecisionAudit implements AutoCloseable {

	public AssignmentDecisionAudit(String prefix,String mapping,ProteinSequence[] sequences,
			int finalPass,int shortlist,boolean overwrite,TailAlignmentGate gate,CandidateFilter filter){
		if(prefix==null || prefix.length()==0 || mapping==null || sequences==null ||
			finalPass<1 || shortlist<1 || gate==null || filter==null || filter.passHitsFlex){
			throw new IllegalArgumentException("Audit needs a prefix, mapping, inputs, final pass, gate and fixed reassignment hit floor");
		}
		this.finalPass=finalPass; this.shortlist=shortlist;
		candidateFamilies=new int[shortlist]; acceptedCandidates=new boolean[shortlist];
		targets=new int[sequences.length]; Arrays.fill(targets,-2);
		seen=new boolean[sequences.length]; readTargets(mapping,sequences);
		final String[] paths=outputPaths(prefix);
		try{
			queries=open(paths[0],overwrite); candidates=open(paths[1],overwrite); references=open(paths[2],overwrite);
			final ByteBuilder meta=new ByteBuilder().append("#final_pass\t").append(finalPass)
				.append("\n#minimum_hits\t").append(filter.minimumHits)
				.append("\n#maximum_hydro_delta\t").append(Double.toString(filter.maximumHydroDelta))
				.append("\n#minimum_identity\t").append(Double.toString(gate.minIdentity))
				.append("\n#minimum_score_ratio\t").append(Double.toString(gate.minScoreRatio))
				.append("\n#minimum_span\t").append(gate.minSpan)
				.append("\n#denominator\t").append(gate.denominator.name())
				.append("\n#decimals\t12; hydro bits and integer alignment counts preserve exact evidence\n");
			queries.print(meta); candidates.print(meta); references.print(meta);
			queries.println("input_index\tinput_id\tcandidates\twinner_family_id\twinner_rank\tquery_hydro_available\tquery_hydro\tquery_hydro_bits\ttarget_requested\ttarget_family_id\ttarget_present\ttarget_rank\ttarget_hits");
			candidates.println("input_index\trank\tfamily_id\treference_members\thits\thydro_available\thydro_delta\tscreen\talignment_attempted\tgate_reason\tscore\tquery_span\treference_span\tmetrics_available\tidentities\tcolumns\tdenominator_self_score\tidentity\tscore_ratio");
			references.println("family_id\treference_members\tconsensus_length\thydro_available\thydro\thydro_bits");
		}catch(RuntimeException failure){
			try{closeStreams();}catch(RuntimeException cleanup){failure.addSuppressed(cleanup);}
			throw failure;
		}
	}

	public static String[] outputPaths(String prefix){
		if(prefix==null || prefix.length()==0){throw new IllegalArgumentException("Empty audit prefix");}
		return new String[]{prefix+".queries.tsv",prefix+".candidates.tsv",prefix+".references.tsv"};
	}

	public void validateSchema(int passes,int limit,int inputs,int threads){
		require(passes==finalPass && limit==shortlist && inputs==targets.length && threads==1,
			"Decision audit supports one serial final pass with the configured input/shortlist schema");
	}

	public int targetFor(int input){return targets[input];}

	public void beginPass(int pass,int families,int eligible){
		require(!started && pass==finalPass && families>=0 && eligible>=0 && eligible<=targets.length,
			"Audit begins once, before final-pass membership is cleared");
		started=true; expectedReferences=families; expectedQueries=eligible;
	}

	/** Called before clearMembers; sizes and scalars belong to the frozen reference. */
	public void reference(Family family,int size){
		require(started && queryCount==0 && !active && family!=null && size>0,
			"Frozen reference records precede all final-pass queries");
		final boolean defined=hydroDefined(family);
		final double hydro=family.scalarsUnsafe()[ConsensusScalars.HYDRO];
		finite(hydro);
		row.clear().append(family.familyId).tab().append(size).tab()
			.append(family.consensusUnsafe().length).tab().append(defined ? 1 : 0).tab();
		value(row,defined,hydro); row.tab(); bits(row,defined,hydro); row.nl();
		references.print(row); referenceCount++;
	}

	/** Ranking has completed. A present unshortlisted target still has an exact hit count. */
	public void beginQuery(int input,String id,int count,Family target,int targetHits,CandidateFilter.Context query){
		require(started && !finished && !active && referenceCount==expectedReferences &&
			input>=0 && input<seen.length && !seen[input] && count>=0 && count<=shortlist && query!=null,
			"Audit query must occur once after the complete frozen reference roster");
		active=true; seen[input]=true; currentInput=input; expectedCandidates=count; writtenCandidates=0;
		currentTarget=targets[input]; targetRank=0;
		require(target==null || (currentTarget>=0 && target.familyId==currentTarget && targetHits>=0),
			"Audit target metadata must describe the requested live family");
		queryRow.clear().append(input).tab().append(id).tab().append(count);
		queryHydro=query.hydrophobicity(); queryHydroDefined=query.defined(); finite(queryHydro);
		targetPresent=target!=null; currentTargetHits=targetHits;
		candidateRows.clear();
	}

	/** Copies the actual screen and gate result without invoking any scorer or mutating it. */
	public void candidate(int rank,Family family,int referenceSize,int hits,
			CandidateFilter.Context query,int screen,TailAlignmentDecision decision){
		require(active && rank==writtenCandidates+1 && rank<=expectedCandidates && family!=null &&
			referenceSize>0 && hits>=0 && screen>=CandidateFilter.PASS && screen<=CandidateFilter.HYDRO &&
			((screen==CandidateFilter.PASS)==(decision!=null)),
			"Audit rows must follow the actual ranked candidates and alignment eligibility");
		final boolean hydroAvailable=query.defined() && hydroDefined(family);
		final double delta=Math.abs(query.hydrophobicity()-family.scalarsUnsafe()[ConsensusScalars.HYDRO]);
		finite(delta);
		candidateRows.append(currentInput).tab().append(rank).tab().append(family.familyId).tab()
			.append(referenceSize).tab().append(hits).tab().append(hydroAvailable ? 1 : 0).tab();
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
		candidateRows.nl(); writtenCandidates++;
		candidateFamilies[rank-1]=family.familyId; acceptedCandidates[rank-1]=decision!=null && decision.accepted();
		if(family.familyId==currentTarget){targetRank=rank;}
		if(screen==CandidateFilter.PASS){alignedCount++;}else{filteredCount++;}
	}

	public void finishQuery(int winner,int rank){
		require(active && writtenCandidates==expectedCandidates && rank>=0 && rank<=expectedCandidates &&
			((winner>=0)==(rank>0)),"Audit query must contain its full shortlist and actual winner");
		require(winner<0 || (candidateFamilies[rank-1]==winner && acceptedCandidates[rank-1]),
			"Audit winner must be the accepted candidate at its recorded rank");
		queryRow.tab().append(winner).tab().append(rank).tab().append(queryHydroDefined ? 1 : 0).tab();
		value(queryRow,queryHydroDefined,queryHydro); queryRow.tab(); bits(queryRow,queryHydroDefined,queryHydro);
		queryRow.tab().append(currentTarget>=-1 ? 1 : 0).tab().append(currentTarget).tab()
			.append(targetPresent ? 1 : 0).tab().append(targetRank).tab();
		if(targetPresent){queryRow.append(currentTargetHits);}else{queryRow.append("NA");}
		queryRow.nl(); queries.print(queryRow); candidates.print(candidateRows);
		queryCount++; candidateCount+=writtenCandidates; active=false;
	}

	public void endPass(){
		require(started && !active && !finished && queryCount==expectedQueries &&
			referenceCount==expectedReferences && alignedCount+filteredCount==candidateCount,
			"Final-pass audit must conserve queries, references and ranked/aligned/filtered pairs");
		for(int i=0;i<targets.length;i++){require(targets[i]<-1 || seen[i],"Requested audit input was not queried: "+i);}
		finished=true;
	}

	private void readTargets(String path,ProteinSequence[] sequences){
		final ByteFile input=ByteFile.makeByteFile1(path,false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			byte[] line=input.nextLine();
			while(line!=null && line.length>0 && line[0]=='#'){line=input.nextLine();}
			require(line!=null && MAP_HEADER.equals(new String(line,java.nio.charset.StandardCharsets.US_ASCII)),"Unexpected audit mapping header");
			int previous=-1;
			for(line=input.nextLine();line!=null;line=input.nextLine()){
				lp.set(line); require(lp.terms()==10,"Audit mapping needs ten fields");
				final int index=integer(lp,0),oldFamily=integer(lp,2),size=integer(lp,3);
				final int target=integer(lp,7),shared=integer(lp,8),ties=integer(lp,9);
				require(index>previous && index<targets.length && oldFamily>=0 && size>0,
					"Audit mapping requires ordered unique input indices and positive old family sizes");
				require(lp.parseString(1).equals(sequences[index].id),"Audit input_id differs at index "+index);
				require(lp.termEquals("unassigned",5) || lp.termEquals("filtered_length",5),"Audit map must identify lost old assignments");
				final String bin=size<5 ? "<5" : size<10 ? "5-9" : size<50 ? "10-49" : "50+";
				require(lp.termEquals(bin,4),"Audit map old-family size bin disagrees");
				if(target<0){require(target==-1 && shared==0 && ties==0 && lp.termEquals("no_proxy",6),"Invalid absent audit proxy");}
				else{
					require(shared>0 && shared<=size && ties>0 && ties<=size &&
						lp.termEquals(ties==1 ? "unique_best" : "tie_lowest_new_family_id",6),"Invalid audit proxy counts/status");
				}
				targets[index]=target; previous=index;
			}
		}finally{require(!input.close(),"Could not read audit mapping: "+path);}
	}

	private static int integer(LineParser1 lp,int field){
		final int length=lp.length(field); require(length>0,"Empty audit integer");
		final boolean negative=lp.parseByte(field,0)=='-'; int value=0,start=negative ? 1 : 0;
		require(length>start,"Missing audit integer digits");
		for(int i=start;i<length;i++){
			final int digit=lp.parseByte(field,i)-'0'; require(digit>=0 && digit<=9,"Invalid audit integer digit");
			value=Math.addExact(Math.multiplyExact(value,10),digit);
		}
		return negative ? -value : value;
	}

	private static boolean hydroDefined(Family family){return (family.scalarMask()&(1<<ConsensusScalars.HYDRO))!=0;}
	private static void finite(double value){require(Double.isFinite(value),"Audit measurements must be finite");}
	private static void value(ByteBuilder bb,boolean available,double number){if(available){bb.append(number,12);}else{bb.append("NA");}}
	private static void bits(ByteBuilder bb,boolean available,double number){if(available){bb.append(Double.doubleToRawLongBits(number));}else{bb.append("NA");}}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalStateException(message);}}
	private static ByteStreamWriter open(String path,boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false); writer.start(); return writer;
	}
	private void closeStreams(){
		boolean error=false;
		if(queries!=null){error|=queries.poisonAndWait(); queries=null;}
		if(candidates!=null){error|=candidates.poisonAndWait(); candidates=null;}
		if(references!=null){error|=references.poisonAndWait(); references=null;}
		require(!error,"Decision audit output failed");
	}
	@Override public void close(){closeStreams(); require(finished,"Decision audit did not complete its final pass");}

	private final int finalPass,shortlist;
	private final int[] targets;
	private final boolean[] seen;
	private final int[] candidateFamilies;
	private final boolean[] acceptedCandidates;
	private ByteStreamWriter queries,candidates,references;
	private final ByteBuilder row=new ByteBuilder(),queryRow=new ByteBuilder(),candidateRows=new ByteBuilder();
	private int expectedReferences,expectedQueries,referenceCount,currentInput,expectedCandidates,writtenCandidates,currentTarget,targetRank,currentTargetHits;
	private long queryCount,candidateCount,alignedCount,filteredCount;
	private double queryHydro;
	private boolean started,finished,active,queryHydroDefined,targetPresent;
	private static final String[] SCREENS={"pass","hits","hydro"};
	public static final String MAP_HEADER="input_index\tinput_id\told_family_id\told_family_size\told_family_size_bin\tnew_status\tproxy_status\tproxy_new_family_id\tproxy_shared_old_members\tproxy_tied_new_families";
}
