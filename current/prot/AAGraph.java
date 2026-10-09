package prot;

import java.util.Arrays;
import java.util.ArrayList;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/**
 * A protein consensus graph: a 20-symbol amino-acid port of
 * {@code consensus.BaseGraph}.
 *
 * <p>The consensus of a cluster is the highest-scoring trace through this graph,
 * so two design points from the nucleotide original are load-bearing:</p>
 *
 * <ul>
 * <li><b>Insertions are real nodes</b> (an {@link AAGraphNode#insEdge} chain off
 * each REF column), not dropped. A residue most members carry can therefore enter
 * the consensus even if the pivot lacked it, instead of being lost to the seed's
 * shape.</li>
 * <li><b>Identity-inverted weighting</b> (optional, {@link #weightByIdentity})
 * up-weights LOW-identity members: {@code weight = baseWeight + max(0,
 * identityCeiling - pident)} on a 0-100 percent scale. This lets one consensus be
 * pulled to fit divergent / under-represented members rather than only the
 * over-represented ones -- and because the consensus IS the winning trace, this
 * changes the consensus itself, not merely a score.</li>
 * </ul>
 *
 * <p>Read-pileup machinery (quality, mapq, SAM) is deliberately omitted; proteins
 * have none of it. This is a correctness-first MVP: standalone (no
 * BaseGraph/ConsensusObject inheritance), single-threaded, O(members) memory.</p>
 *
 * @author Eru
 */
public final class AAGraph {

	/*--------------------------------------------------------------*/
	/*----------------        Tunable knobs         ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether to up-weight low-identity members (ANI-inverted weighting). */
	public boolean weightByIdentity=false;
	/** Percent-identity (0-100) at/above which a member gets no identity bonus. */
	public float identityCeiling=40f;
	/** Minimum per-column weight of any member observation. */
	public int baseWeight=1;
	/** Min allele fraction (by count) for a substitution to override the pivot residue. */
	public float MAF_sub=0.25f;
	/** Min deletion allele fraction (by count) required to drop a reference column. */
	public float MAF_del=0.5f;
	/** Min insertion allele fraction (by count) required to include an inserted column. */
	public float MAF_ins=0.5f;
	/** Minimum depth for a residue to win a column. */
	public int minDepth=1;
	/** Trim consensus ends whose depth is below this fraction of the max depth. */
	public float trimDepthFraction=0f;

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds a graph scaffolded on a pivot sequence, padded with X on both ends so
	 * members that overhang the pivot can grow the consensus outward.
	 *
	 * @param pivotEnc Encoded pivot residues (the cluster centroid/seed).
	 * @param pad_ Number of X residues to pad on each end.
	 */
	public AAGraph(byte[] pivotEnc, int pad_){
		pad=pad_;
		pivot=padX(pivotEnc, pad);
		ref=new AAGraphNode[pivot.length];
		del=new AAGraphNode[pivot.length];
		for(int i=0; i<pivot.length; i++){
			ref[i]=new AAGraphNode(pivot[i], AAGraphNode.REF, i);
			del[i]=new AAGraphNode(pivot[i], AAGraphNode.DEL, i);
		}
		for(int i=0; i<pivot.length; i++){//Weight-1 scaffold so every column has a defined residue.
			ref[i].add(pivot[i], 1);
		}
	}

	/** Pads both ends of an encoded sequence with X (Blosum62.X_CODE). */
	private static byte[] padX(byte[] in, int pad){
		if(pad<1){return in;}
		final byte[] out=new byte[in.length+2*pad];
		Arrays.fill(out, Blosum62.X_CODE);
		for(int i=0; i<in.length; i++){out[i+pad]=in[i];}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds a cluster member using a path-recorded glocal alignment to this graph's
	 * padded pivot ({@link AAAligner#alignGlocal(byte[],byte[],boolean)} with
	 * recordPath=true). The whole member is placed; overhang maps into the padding
	 * and can grow the consensus in either direction.
	 *
	 * @param memberEnc The member's encoded residues (the alignment's query).
	 * @param aln The member-vs-pivot glocal alignment with a populated match path.
	 */
	public void add(byte[] memberEnc, AAAlignment aln){
		final byte[] match=aln.match;
		assert(match!=null) : "add() needs a path-recorded alignment (alignGlocal(...,true)).";
		int weight=baseWeight;
		if(weightByIdentity){weight+=Tools.max(0, (int)(identityCeiling-aln.pident()));}
		addTraceInternal(memberEnc, aln.qStart, aln.tStart, match, weight);
	}

	/**
	 * Adds a full-query m/I/D trace without assigning a BLOSUM meaning to its score.
	 * Validates the entire path before mutation and derives identity from actual residues.
	 * Returns residues excluded by the existing leading/trailing-overhang policy.
	 */
	public int addTrace(final byte[] memberEnc, final int targetStart, final byte[] path){
		if(path==null || path.length==0 || memberEnc.length==0 || targetStart<0 || targetStart>=ref.length){
			throw new IllegalArgumentException("Graph accumulation needs a nonempty full-query trace anchored inside its pivot");
		}
		for(byte aa : memberEnc){
			if(aa<0 || aa>Blosum62.X_CODE || aa==20){throw new IllegalArgumentException("Invalid graph member residue: "+aa);}
		}
		int q=0, r=targetStart, identities=0;
		for(byte op : path){
			if(op=='m'){
				if(q>=memberEnc.length || r>=ref.length){throw new IllegalArgumentException("Graph trace emission exceeds query or pivot");}
				if(memberEnc[q]==pivot[r] && memberEnc[q]<20){identities++;}
				q++; r++;
			}else if(op=='I'){q++;}
			else if(op=='D'){r++;}
			else{throw new IllegalArgumentException("Unknown graph trace operation: "+op);}
			if(q>memberEnc.length || r>ref.length){throw new IllegalArgumentException("Graph trace exceeds query or pivot");}
		}
		if(q!=memberEnc.length){throw new IllegalArgumentException("Graph trace must consume the complete member");}
		int weight=baseWeight;
		if(weightByIdentity){weight+=Tools.max(0, (int)(identityCeiling-100.0*identities/path.length));}
		return addTraceInternal(memberEnc, 0, targetStart, path, weight);
	}

	/** Shared legacy accumulation; callers select validation and weighting contracts. */
	private int addTraceInternal(byte[] memberEnc, int qpos, int rpos, byte[] match, int weight){
		int dropped=0;
		AAGraphNode prevNode=(rpos<=0 ? null : ref[rpos-1]);
		for(int mpos=0; mpos<match.length && rpos<ref.length; mpos++){
			final byte m=match[mpos];
			final AAGraphNode next;
			if(m=='m'){
				next=ref[rpos];
				next.add(memberEnc[qpos], weight);
				qpos++; rpos++;
			}else if(m=='D'){
				next=del[rpos];
				next.add((byte)-1, weight);//DEL: only the count/weight sums matter
				rpos++;
			}else{//'I'
				//A leading insertion (prevNode==null) means the member overhangs the pivot's left
				//edge beyond the padding: those residues have no node to attach to, so drop them
				//(a few N-terminal residues of one member, under-counted) rather than crash. Enlarge
				//pad if this happens often. prevNode stays null until the first anchored column.
				if(prevNode==null){qpos++; dropped++; continue;}
				if(prevNode.insEdge==null){
					prevNode.insEdge=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, rpos);
				}
				next=prevNode.insEdge;
				next.add(memberEnc[qpos], weight);
				qpos++;
			}
			prevNode=next;
		}
		return dropped+memberEnc.length-qpos;
	}

	/**
	 * Traverses the graph to the highest-weight trace: the plurality residue at each
	 * kept reference column, plus insertion columns whose weight and allele fraction
	 * clear the thresholds. Low-depth and residual-X ends are trimmed.
	 *
	 * @return The consensus as encoded residues (0-19 or X_CODE).
	 */
	public byte[] traverse(){
		return traverse(null);
	}

	/** Returns detached selected column counts in exactly the emitted, trimmed consensus frame. */
	public Traversal traverseProfile(){
		final ArrayList<AAGraphNode> selected=new ArrayList<AAGraphNode>();
		final byte[] consensus=traverse(selected);
		return new Traversal(consensus, selected);
	}

	private byte[] traverse(final ArrayList<AAGraphNode> selected){
		final ByteBuilder bb=new ByteBuilder();
		final IntList depthList=new IntList();
		int maxDepth=0;

		for(int i=0; i<ref.length; i++){
			final AAGraphNode dnode=del[i], rnode=ref[i];
			AAGraphNode inode=rnode.insEdge;

			final int dw=dnode.weightSum, rw=rnode.weightSum;
			final int dc=dnode.countSum, rc=rnode.countSum;
			final int depth=dc+rc;
			maxDepth=Tools.max(maxDepth, depth);

			final float afMult=1f/Tools.max(1, dc+rc);
			final float daf=dc*afMult;
			final long weightSum=dw+rw;

			if(rw>=dw || daf<MAF_del){//Keep this reference column.
				bb.append(rnode.consensus(MAF_sub, minDepth));
				depthList.add(depth);
				if(selected!=null){selected.add(rnode);}
				//Then walk the insertion chain while it is the plurality outgoing allele.
				while(inode!=null && inode.weightSum>=(weightSum-inode.weightSum)
						&& inode.countSum*afMult>=MAF_ins){
					bb.append(inode.consensus(MAF_ins, minDepth));
					depthList.add(depth);
					if(selected!=null){selected.add(inode);}
					inode=inode.insEdge;
				}
			}//else: deletion column, emit nothing
		}

		final byte[] cons=bb.toBytes();
		return trim(cons, depthList, maxDepth, selected);
	}

	/** Trims residual-X and (optionally) low-depth residues from both consensus ends. */
	private byte[] trim(byte[] cons, IntList depthList, int maxDepth, ArrayList<AAGraphNode> selected){
		final int trimDepth=(trimDepthFraction>0 ? Tools.max(1, (int)(trimDepthFraction*maxDepth)) : 0);
		int left=0, right=0;
		while(left<cons.length &&
				(cons[left]==Blosum62.X_CODE || depthList.get(left)<trimDepth)){left++;}
		while(right<cons.length-left &&
				(cons[cons.length-right-1]==Blosum62.X_CODE
				|| depthList.get(cons.length-right-1)<trimDepth)){right++;}
		if(selected!=null){
			selected.subList(cons.length-right, selected.size()).clear();
			selected.subList(0, left).clear();
		}
		return (left==0 && right==0) ? cons : Arrays.copyOfRange(cons, left, cons.length-right);
	}

	/** Immutable emitted-column snapshot; it does not claim to rebase the complete DEL/INS graph. */
	public static final class Traversal{
		private final byte[] consensus;
		private final int[][] counts;
		private final int[] types, positions;
		private Traversal(byte[] sequence, ArrayList<AAGraphNode> nodes){
			assert(sequence.length==nodes.size()) : "Selected columns must follow exactly the consensus trim";
			consensus=sequence.clone(); counts=new int[nodes.size()][];
			types=new int[nodes.size()]; positions=new int[nodes.size()];
			for(int i=0; i<nodes.size(); i++){
				final AAGraphNode node=nodes.get(i);
				counts[i]=node.count.clone(); types[i]=node.type; positions[i]=node.rpos;
			}
		}
		public byte[] consensus(){return consensus.clone();}
		public int length(){return consensus.length;}
		public int count(int column, int residue){return counts[column][residue];}
		public int sourceType(int column){return types[column];}
		public int sourcePosition(int column){return positions[column];}
	}

	/** Convenience: the consensus decoded to an amino-acid string. */
	public String consensusString(){return decode(traverse());}

	/**
	 * Decodes encoded residues to an amino-acid string (X_CODE renders as 'X').
	 * @param enc Encoded residues (0-19 or X_CODE).
	 * @return The amino-acid letters.
	 */
	public static String decode(byte[] enc){
		final ByteBuilder bb=new ByteBuilder(enc.length);
		for(byte e : enc){
			bb.append(e==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[e]);
		}
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Padding residues added to each end of the pivot. */
	public final int pad;
	/** The X-padded, encoded pivot the graph is scaffolded on. */
	public final byte[] pivot;
	/** Reference (pivot) column nodes, one per padded pivot position. */
	public final AAGraphNode[] ref;
	/** Deletion column nodes, one per padded pivot position. */
	public final AAGraphNode[] del;
}
