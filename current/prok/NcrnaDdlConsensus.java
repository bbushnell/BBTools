package prok;

import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import stream.Read;
import structures.ByteBuilder;

/** Quantum consensus refinement with one real placement per member. A proposal
 * is published only if every member still passes; its retained traces also
 * supply the final HBM, without another alignment pass or a synthetic pivot vote.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlConsensus{

	/*--------------------------------------------------------------*/
	/*----------------          Refinement          ----------------*/
	/*--------------------------------------------------------------*/

	public static Fit refine(final ArrayList<Read> members, final byte[] initial,
			final float minimumIdentity, final int passes, final boolean trim){
		if(members==null || members.isEmpty() || initial==null || initial.length==0 || passes<0
				|| !Float.isFinite(minimumIdentity) || minimumIdentity<0 || minimumIdentity>1){
			throw new IllegalArgumentException("Consensus refinement needs real members, an accepted model, finite gate and nonnegative pass limit");
		}
		byte[] bases=initial;
		Placement placed=place(members, bases, minimumIdentity);
		long alignments=placed.alignments;
		if(!placed.success){throw new IllegalStateException("An assigned cluster member failed Quantum placement against its accepted consensus");}
		for(int pass=0; pass<passes; pass++){
			final byte[] proposal=consensus(bases, placed.reads, trim);
			if(proposal.length==0 || Arrays.equals(proposal, bases)){break;}
			final Placement next=place(members, proposal, minimumIdentity);
			alignments+=next.alignments;
			if(!next.success){break;}
			bases=proposal;
			placed=next;
		}
		assert(placed.reads.length==members.size()) : "A published consensus retains one actual Quantum placement for every assigned member";
		return new Fit(bases, placed, alignments);
	}

	/** Both orientations align the shorter operand as the full glocal query,
	 * matching Quantum's ordinary identity path. Failed proposals retain no model. */
	private static Placement place(final ArrayList<Read> members, final byte[] reference, final float minimumIdentity){
		final Placement result=new Placement(members.size());
		final Read[] placed=result.reads;
		final AlignmentStats stats=new AlignmentStats(true);
		for(int i=0; i<placed.length; i++){
			final Read member=members.get(i);
			if(member==null || member.bases==null || member.length()==0){
				throw new IllegalArgumentException("Consensus members must be real nonempty nucleotide sequences");
			}
			stats.clear();
			stats.doTrace=true;
			final boolean direct=member.length()<=reference.length;
			result.alignments++;
			final float identity=direct ? QuantumAligner.alignAndTraceStatic(member.bases, reference, stats)
				: QuantumAligner.alignAndTraceStatic(reference, member.bases, stats);
			if(!Float.isFinite(identity) || identity<minimumIdentity || stats.matchString==null){return result;}
			result.identities[i]=identity;
			placed[i]=toModelFrame(member, stats, direct);
			final Read read=placed[i];
			if(read==null){return result;}
			assert(read.start>=0 && read.stop<reference.length && read.stop>=read.start
					&& TrnaConsensusBuilder.alignmentConsumption(read.match, true)==read.length()
					&& TrnaConsensusBuilder.alignmentConsumption(read.match, false)==read.stop-read.start+1)
				: "Quantum trace conversion must conserve member/model consumption in the final HBM coordinate frame";
		}
		result.success=true;
		return result;
	}

	/** Converts an existing Quantum trace, without alignment. Direct means the
	 * member was query; otherwise the model was query. The caller owns stats. */
	public static Read toModelFrame(final Read member, final AlignmentStats stats, final boolean direct){
		if(member==null || member.bases==null || stats==null || stats.matchString==null){
			throw new IllegalArgumentException("Trace conversion requires an actual member and its recorded alignment");
		}
		if(!direct){return TrnaConsensusBuilder.invertToModelFrame(member.bases, stats);}
		stats.matchString=TrnaConsensusBuilder.normalizeOpposingIndels(stats.matchString);
		//BaseGraph reads the query at every operation, including D. Floating
		//terminal reference gaps cannot be retained after the last query base.
		int from=0, to=stats.matchString.length;
		while(from<to && stats.matchString[from]=='D'){from++;stats.rStart++;}
		while(to>from && stats.matchString[to-1]=='D'){to--;stats.rStop--;}
		if(from==to){return null;}
		if(from>0 || to<stats.matchString.length){stats.matchString=Arrays.copyOfRange(stats.matchString, from, to);}
		return TrnaConsensusBuilder.toAlignedRead(member.bases, stats, member.id, member.numericID);
	}

	/** Same nucleotide plurality/deletion and modal-end policy as the existing
	 * builder; the old reference supplies coordinates, never an additional vote. */
	private static byte[] consensus(final byte[] reference, final Read[] placed, final boolean trim){
		assert(reference.length>0 && placed.length>0) : "Only a nonempty actual cluster can propose a consensus";
		final int[] counts=new int[reference.length*5];
		final int[] starts=new int[reference.length], stops=new int[reference.length];
		for(final Read read:placed){
			starts[read.start]++;
			stops[read.stop]++;
			int q=0, r=read.start;
			for(final byte op:read.match){
				if(op=='I'){q++;continue;}
				if(op=='D'){counts[r*5+4]++;r++;continue;}
				assert(op=='m' || op=='S' || op=='N') : "BaseGraph placement operations must remain nucleotide diagonal/insertion/deletion operations";
				final int base=AminoAcid.baseToNumber[read.bases[q]];
				if(base>=0 && base<4){counts[r*5+base]++;}
				q++;
				r++;
			}
		}
		int first=0, last=reference.length-1;
		if(trim){
			for(int i=0; i<reference.length; i++){
				if(starts[i]>starts[first]){first=i;}
				if(stops[i]>stops[last]){last=i;}
			}
			if(first>=last){first=0;last=reference.length-1;}
		}
		final ByteBuilder out=new ByteBuilder(reference.length);
		for(int r=first; r<=last; r++){
			final int offset=r*5;
			int total=0, best=0;
			for(int base=0; base<4; base++){
				total+=counts[offset+base];
				if(counts[offset+base]>counts[offset+best]){best=base;}
			}
			if(total==0 && counts[offset+4]==0){out.append(reference[r]);}
			else if(counts[offset+4]<=total){out.append(AminoAcid.numberToBase[best]);}
		}
		return out.toBytes();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Final Model         ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Fit{
		private Fit(final byte[] bases_, final Placement placement, final long alignments_){
			bases=bases_;
			placed=placement.reads;
			identities=placement.identities;
			alignments=alignments_;
		}

		/** Makes an HBM from the already retained final-reference traces. */
		public BaseGraph graph(final String name){
			final BaseGraph graph=new BaseGraph(name, bases, null, 0, 0);
			//BaseGraph's constructor initializes a reference observation. The pivot
			//is a coordinate frame here; only the actual members below may vote.
			for(int pos=0; pos<bases.length; pos++){BaseGraphHelper.setCounts(graph, pos, 0, 0, 0, 0, 0);}
			for(final Read read:placed){graph.add(read);}
			BaseGraphHelper.initForScoring(graph);
			assert(Arrays.equals(graph.original, bases)) : "The HBM and consensus FASTA must name the identical final nucleotide coordinates";
			return graph;
		}

		public int members(){return placed.length;}
		public final byte[] bases;
		/** Exact Quantum identities in the input member order; no reporting realignment. */
		public final float[] identities;
		/** Every executed placement, including rejected refinement proposals. */
		public final long alignments;
		private final Read[] placed;
	}

	private static final class Placement{
		Placement(final int members){
			reads=new Read[members];
			identities=new float[members];
		}
		final Read[] reads;
		final float[] identities;
		long alignments;
		boolean success;
	}
}
