package prot.cluster;

import java.util.Arrays;

import prot.AAGraph;
import prot.AAGraphNode;
import prot.Blosum62;
import prot.HbmPositionModel;
import prot.ProteinSubstitutionMatrix;
import structures.ByteBuilder;
import structures.IntList;

/** Single-threaded consensus accumulator in a fixed representative frame.
 * Only accepted decisions contribute. Insertions are attached to reference boundaries,
 * including boundaries zero and length, so terminal insertions and insertions following
 * a deletion cannot disappear during accumulation. Every member has unit weight.
 * @author Keqing
 */
public final class FamilyConsensus {

	public FamilyConsensus(byte[] pivotEncoded, boolean pinRight, boolean trimLeft, float trimDepthFraction){
		this(pivotEncoded, pinRight ? TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT :
			TailAlignmentGate.Mode.FULL_QUERY_UNPINNED, trimLeft, trimDepthFraction);
	}

	public FamilyConsensus(byte[] pivotEncoded, TailAlignmentGate.Mode mode, boolean trimEnd, float trimDepthFraction){
		this(pivotEncoded, mode, trimEnd, trimDepthFraction, ProteinSubstitutionMatrix.BLOSUM62);
	}

	/** Trace validation uses the same score and gap units as the producing gate. */
	public FamilyConsensus(byte[] pivotEncoded, TailAlignmentGate.Mode mode, boolean trimEnd,
			float trimDepthFraction, ProteinSubstitutionMatrix matrix_){
		if(matrix_==null){throw new IllegalArgumentException("Consensus needs an explicit substitution matrix");}
		matrix=matrix_;
		if(pivotEncoded==null || pivotEncoded.length==0){throw new IllegalArgumentException("Consensus needs a nonempty pivot");}
		if(mode==null){throw new IllegalArgumentException("Consensus needs an alignment mode");}
		if(!Float.isFinite(trimDepthFraction) || trimDepthFraction<0 || trimDepthFraction>1){
			throw new IllegalArgumentException("Consensus trim depth fraction must be in [0,1]");
		}
		for(byte aa : pivotEncoded){validateResidue(aa);}
		graph=new AAGraph(pivotEncoded.clone(), 0);
		// The pivot defines coordinates; it is not an extra member observation.
		for(AAGraphNode node : graph.ref){
			Arrays.fill(node.count, 0); Arrays.fill(node.weight, 0);
			node.countSum=node.weightSum=0;
		}
		insertions=new AAGraphNode[pivotEncoded.length+1];
		boundaryDepth=new int[insertions.length];
		pinRight=mode==TailAlignmentGate.Mode.PIN_RIGHT_TRIM_LEFT ||
			mode==TailAlignmentGate.Mode.PIN_RIGHT_FULL_QUERY || mode==TailAlignmentGate.Mode.PIN_RIGHT_OVERHANG_ONLY;
		pinLeft=mode==TailAlignmentGate.Mode.PIN_LEFT_TRIM_RIGHT;
		trimLeft=!pinLeft && trimEnd;
		trimRight=pinLeft ? trimEnd : !pinRight;
		this.trimDepthFraction=trimDepthFraction;
	}

	/** Validates the whole accepted trace before changing any counts. Rejections are inert. */
	public boolean add(byte[] query, TailAlignmentDecision decision){
		if(decision==null){throw new IllegalArgumentException("Missing clustering decision");}
		if(!decision.accepted()){return false;}
		validateTrace(query, decision);
		if(members==Integer.MAX_VALUE){throw new IllegalStateException("Consensus observation counts exceed integer capacity");}
		final long newResidues=Math.addExact(residues, decision.querySpan);
		for(int r=decision.referenceStart; r<=decision.referenceEnd+1; r++){boundaryDepth[r]++;}
		int q=decision.queryStart, r=decision.referenceStart;
		AAGraphNode previousInsertion=null;
		for(byte op : decision.pathUnsafe()){
			if(op=='m'){
				graph.ref[r++].add(query[q++], 1); previousInsertion=null;
			}else if(op=='D'){
				graph.del[r++].add((byte)-1, 1); previousInsertion=null;
			}else{
				assert(op=='I') : "validateTrace permits only m/I/D before graph mutation";
				AAGraphNode node=previousInsertion==null ? insertions[r] : previousInsertion.insEdge;
				if(node==null){
					node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, r);
					if(previousInsertion==null){insertions[r]=node;}
					else{previousInsertion.insEdge=node;}
				}
				node.add(query[q++], 1); previousInsertion=node;
			}
		}
		assert(q==decision.queryEnd+1 && r==decision.referenceEnd+1) :
			"Accumulation must consume exactly the trace validated before mutation";
		members++; residues=newResidues;
		return true;
	}

	/** Selects plurality columns using AAGraph allele thresholds; pinned ends are not trimmed. */
	public byte[] consensus(){
		if(members==0){throw new IllegalStateException("Cannot derive consensus without member observations");}
		final ByteBuilder sequence=new ByteBuilder();
		final IntList depths=new IntList();
		int maximum=0;
		for(int r=0; r<=graph.ref.length; r++){
			final int coverage=boundaryDepth[r];
			for(AAGraphNode node=insertions[r]; node!=null; node=node.insEdge){
				if(coverage==0 || node.countSum/(double)coverage<graph.MAF_ins){break;}
				sequence.append(node.consensus(graph.MAF_ins, graph.minDepth));
				depths.add(node.countSum); maximum=Math.max(maximum, node.countSum);
			}
			if(r==graph.ref.length){break;}
			final AAGraphNode ref=graph.ref[r], del=graph.del[r];
			final int depth=Math.addExact(ref.countSum, del.countSum);
			if(depth==0){continue;}// An unobserved scaffold column is not evidence.
			if(ref.weightSum>=del.weightSum || del.countSum/(double)depth<graph.MAF_del){
				sequence.append(ref.consensus(graph.MAF_sub, graph.minDepth));
				depths.add(depth); maximum=Math.max(maximum, depth);
			}
		}
		final byte[] result=sequence.toBytes();
		assert(result.length==depths.size) : "Each emitted consensus column needs its own coverage for trimming";
		final int minimum=trimDepthFraction>0 ? Math.max(1, (int)(maximum*trimDepthFraction)) : 0;
		int left=0, right=result.length;
		if(trimLeft){
			while(left<right && (result[left]==Blosum62.X_CODE || depths.get(left)<minimum)){left++;}
		}
		if(trimRight){
			while(right>left && (result[right-1]==Blosum62.X_CODE || depths.get(right-1)<minimum)){right--;}
		}
		if(left==right){throw new EmptyConsensusException();}
		return left==0 && right==result.length ? result : Arrays.copyOfRange(result, left, right);
	}

	/** Ensures a decision cannot accidentally accumulate a different query or representative trace. */
	private void validateTrace(byte[] query, TailAlignmentDecision decision){
		final byte[] path=decision.pathUnsafe();
		if(query==null || path==null || path.length!=decision.columns || path.length==0 ||
			decision.queryStart<0 || decision.queryEnd<decision.queryStart ||
			decision.queryEnd>=query.length ||
			decision.referenceStart<0 || decision.referenceEnd>=graph.pivot.length ||
			decision.referenceEnd<decision.referenceStart ||
			(!pinLeft && decision.queryEnd!=query.length-1) ||
			(pinRight && (decision.queryEnd!=query.length-1 ||
				decision.referenceEnd!=graph.pivot.length-1)) ||
			(pinLeft && (decision.queryStart!=0 || decision.referenceStart!=0))){
			throw new IllegalArgumentException("Accepted consensus decision has inconsistent sequence bounds");
		}
		int q=decision.queryStart, r=decision.referenceStart, identities=0, paired=0;
		long score=0, self=0;
		for(byte op : path){
			if(op=='m'){
				if(q>decision.queryEnd || r>decision.referenceEnd){throw new IllegalArgumentException("Consensus paired trace overrun");}
				validateResidue(query[q]);
				identities+=(query[q]==graph.pivot[r] && Blosum62.isStandard(query[q]) ? 1 : 0);
				score+=HbmPositionModel.SCALE*matrix.score(query[q], graph.pivot[r]);
				self+=HbmPositionModel.SCALE*matrix.score(query[q], query[q]);
				paired++; q++; r++;
			}else if(op=='I'){
				if(q>decision.queryEnd){throw new IllegalArgumentException("Consensus insertion trace overrun");}
				validateResidue(query[q]);
				self+=HbmPositionModel.SCALE*matrix.score(query[q], query[q]);
				score-=HbmPositionModel.SCALE*matrix.linearGap; q++;
			}else if(op=='D'){
				if(r>decision.referenceEnd){throw new IllegalArgumentException("Consensus deletion trace overrun");}
				score-=HbmPositionModel.SCALE*matrix.linearGap; r++;
			}else{throw new IllegalArgumentException("Unknown consensus trace operation: "+(char)op);}
		}
		if(q!=decision.queryEnd+1 || r!=decision.referenceEnd+1 ||
			identities!=decision.identities || paired!=decision.paired || score!=decision.score ||
			self!=decision.alignedSelfScore || decision.querySpan!=q-decision.queryStart){
			throw new IllegalArgumentException("Consensus input does not match the recorded gate decision");
		}
	}

	private static void validateResidue(byte aa){
		if(aa<0 || aa>Blosum62.X_CODE || aa==20){throw new IllegalArgumentException("Invalid encoded consensus residue: "+aa);}
	}

	public int acceptedMembers(){return members;}
	public long acceptedResidues(){return residues;}
	/** A valid accumulation whose explicit column-selection policy emits no sequence. */
	public static final class EmptyConsensusException extends IllegalStateException {
		private EmptyConsensusException(){super("Consensus contains no retained columns after explicit trimming");}
		private static final long serialVersionUID=1L;
	}
	private final AAGraph graph;
	private final ProteinSubstitutionMatrix matrix;
	private final AAGraphNode[] insertions;
	private final int[] boundaryDepth;
	private final boolean pinRight, pinLeft, trimLeft, trimRight;
	private final float trimDepthFraction;
	private int members;
	private long residues;
}
