package assemble;

import java.util.ArrayList;

import assemble.HollowTipConnections.Connection;
import dna.AminoAcid;
import map.LongObjectMap;
import structures.ByteBuilder;
import ukmer.HollowAnchorTable;

/**
 * Experimental second-pass gap verification and conservative path spelling.
 * Each candidate uses its selected exact anchor pair and read windows matching
 * all known contig flanks. Gap strings are voted as whole sequences, never as
 * a synthetic per-base mosaic. One read record contributes at most one vote;
 * repeated occurrences of its selected pair invalidate that record's vote.
 * Read records are not necessarily independent molecules. Cycles are not joined.
 * @author Brian Bushnell, Fischl
 */
public final class HollowBridgeJoiner {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Binds the unchanged discovery snapshot and explicit sequence-vote thresholds. */
	public HollowBridgeJoiner(final ArrayList<Contig> contigs_, final HollowTipConnections finder,
			final ArrayList<Connection> connections, final int minimum_, final int ratio_){
		this(contigs_, finder, connections, minimum_, ratio_, 0, 16);
	}

	/** Optional outer-flank voting reuses the selected-pair pass; zero flank bases preserves old behavior. */
	public HollowBridgeJoiner(final ArrayList<Contig> contigs_, final HollowTipConnections finder,
			final ArrayList<Connection> connections, final int minimum_, final int ratio_,
			final int flankBases_, final int flankMin_){
		this(contigs_, finder, connections, minimum_, ratio_, flankBases_, flankMin_, false);
	}

	/** Scaffold mode preserves complete contigs and marks each uncertain adjacency with exactly ten Ns. */
	public HollowBridgeJoiner(final ArrayList<Contig> contigs_, final HollowTipConnections finder,
			final ArrayList<Connection> connections, final int minimum_, final int ratio_,
			final int flankBases_, final int flankMin_, final boolean fillGaps_){
		if(contigs_==null || finder==null || connections==null || minimum_<1 || ratio_<2){
			throw new IllegalArgumentException("Gap verification requires discovery, positive support and vote ratio>=2.");
		}
		contigs=contigs_;
		span=finder.span;
		minimum=minimum_;
		ratio=ratio_;
		fillGaps=fillGaps_;
		if(flankBases_<0 || flankBases_>4096 || flankMin_<1 || (flankBases_>0 && flankMin_>flankBases_)){
			throw new IllegalArgumentException("Outer-flank verification requires 1<=flankmin<=flankbases<=4096, or flankbases=0.");
		}
		flankBases=flankBases_;
		flankMin=flankMin_;
		queryFlank=new byte[flankBases];
		targetFlank=new byte[flankBases+1];
		delayed=new long[span-32];
		for(Connection connection : connections){
			final Evidence e=new Evidence(connection, finder.anchors);
			evidence.add(e);
			if(connection.conflicted()){e.reason="end_conflict"; continue;}
			if(connection.slideReason()!=null){e.reason=connection.slideReason(); continue;}
			if(!overlapMatches(e)){e.reason="overlap_mismatch"; continue;}
			index(connection.firstWord(), connection.secondWord(), false, e);
			index(HollowAnchorTable.reverseComplement(connection.secondWord()),
					HollowAnchorTable.reverseComplement(connection.firstWord()), true, e);
		}
	}

	/** Adds a retained selector; collisions enumerate rather than replace candidates. */
	private void index(final long first, final long second, final boolean reverse, final Evidence e){
		assert(e!=null) : "Every selector must retain its geometry and vote owner.";
		selectors.put(first, new Selector(second, reverse, e, selectors.get(first)));
	}

	/*--------------------------------------------------------------*/
	/*----------------      Read verification       ----------------*/
	/*--------------------------------------------------------------*/

	/** Streams one complete record, preserving read boundaries and rolling without per-window allocation. */
	public void addRead(final byte[] bases){
		if(finished){throw new IllegalStateException("Gap evidence is frozen after joining.");}
		if(bases==null){throw new IllegalArgumentException("Gap evidence must contain sequence.");}
		if(++readSerial<=0){throw new IllegalStateException("Read serial overflow would corrupt vote deduplication.");}
		long word=0;
		int valid=0, slot=0;
		for(int pos=0; pos<bases.length; pos++){
			final int base=AminoAcid.baseToNumber[bases[pos]];
			if(base<0){word=0; valid=0;}
			else{word=(word<<2)|base; valid++;}
			final long first=delayed[slot];
			delayed[slot]=word;
			if(++slot==delayed.length){slot=0;}
			if(valid<span){continue;}
			for(Selector s=selectors.get(first); s!=null; s=s.next){
				if(s.second==word){observe(s, bases, pos-span+1);}
			}
		}
	}

	/** Removes duplicate-locus votes, verifies known sequence, and retains only distinct gap strings. */
	private void observe(final Selector selector, final byte[] bases, final int start){
		final Evidence e=selector.evidence;
		assert(start>=0 && start+span<=bases.length) : "A selected pair must identify one complete read window.";
		if(e.lastRead==readSerial){
			if(e.lastVote!=null){e.lastVote.count--; e.lastVote=null;}
			if(e.lastFlank<0){e.outerContradictory--;}
			else if(e.lastFlank>0){e.outerConfirmed--;}
			e.lastFlank=0;
			e.repeated++;
			return;
		}
		e.lastRead=readSerial;
		e.lastVote=null;
		e.lastFlank=0;
		if(e.overflow){return;}
		final Connection c=e.connection;
		final byte[] left=contigs.get(c.leftEnd/2).bases, right=contigs.get(c.rightEnd/2).bases;
		for(int i=0; i<span; i++){
			final int base=readBase(bases, start, i, selector.reverse);
			if((i<e.sourceEnd && base!=contigBase(left, left.length-e.sourceEnd+i, (c.leftEnd&1)==0))
					|| (i>=e.targetStart && base!=contigBase(right, i-e.targetStart, (c.rightEnd&1)!=0))){
				e.flankRejected++;
				return;
			}
		}
		if(flankBases>0){
			final int before=checkFlank(bases, start, selector.reverse, left,
					(c.leftEnd&1)==0, left.length-e.sourceEnd-1, -1);
			final int after=checkFlank(bases, start, selector.reverse, right,
					(c.rightEnd&1)!=0, span-e.targetStart, 1);
			if(before<0 || after<0){
				e.outerContradictory++;
				e.lastFlank=-1;
				return;
			}else if(before>0 || after>0){
				e.outerConfirmed++;
				e.lastFlank=1;
			}
		}
		final int length=(int)Math.max(0, c.gap);
		Variant match=null;
		for(Variant variant : e.variants){
			boolean equal=true;
			for(int i=0; i<length && equal; i++){
				equal=AminoAcid.baseToNumber[variant.bases[i]]==readBase(bases, start, e.sourceEnd+i, selector.reverse);
			}
			if(equal){match=variant; break;}
		}
		if(match==null){
			if(e.variants.size()>=MAX_VARIANTS){e.overflow=true; return;}
			final byte[] gap=new byte[length];
			for(int i=0; i<length; i++){gap[i]=ALPHABET[readBase(bases, start, e.sourceEnd+i, selector.reverse)];}
			match=new Variant(gap);
			e.variants.add(match);
		}
		match.count++;
		e.lastVote=match;
	}

	/** Checks one outward flank: -1 contradicts, 0 lacks usable sequence, 1 agrees within one edit. */
	private int checkFlank(final byte[] read, final int start, final boolean reverseRead,
			final byte[] contig, final boolean reverseContig, final int contigStart, final int direction){
		assert(direction==-1 || direction==1) : "Flanks walk away from the anchored window, never through its unknown gap.";
		final int readAvailable=direction<0 ? (reverseRead ? read.length-start-span : start)
				: (reverseRead ? start : read.length-start-span);
		final int contigAvailable=direction<0 ? contigStart+1 : contig.length-contigStart;
		// Reserve one reference base so a deletion at the far end cannot manufacture a contradiction.
		final int length=Math.min(flankBases, Math.min(readAvailable, contigAvailable-1));
		if(length<flankMin){return 0;}
		for(int i=0; i<=length; i++){
			final int base=contigBase(contig, contigStart+direction*i, reverseContig);
			if(base<0){return 0;}
			targetFlank[i]=(byte)base;
			if(i<length){
				final int observed=readBase(read, start, direction<0 ? -1-i : span+i, reverseRead);
				if(observed<0){return 0;}
				queryFlank[i]=(byte)observed;
			}
		}
		return flankAgrees(queryFlank, length, targetFlank, length+1) ? 1 : -1;
	}

	/** Anchored prefix comparison allowing one substitution, insertion or deletion and a free far reference end. */
	public static boolean flankAgrees(final byte[] query, final int length, final byte[] target, final int targetLength){
		assert(length>0 && length<=query.length && targetLength>=length && targetLength<=target.length) :
				"Outer-flank comparison consumes the full observed query against a sufficiently long anchored reference prefix.";
		int i=0;
		while(i<length && query[i]==target[i]){i++;}
		if(i==length){return true;}
		return suffixMatches(query, length, i+1, target, targetLength, i+1)
				|| suffixMatches(query, length, i+1, target, targetLength, i)
				|| suffixMatches(query, length, i, target, targetLength, i+1);
	}

	/** After spending the one edit, the remaining observed suffix must match exactly. */
	private static boolean suffixMatches(final byte[] query, final int length, final int qStart,
			final byte[] target, final int targetLength, final int tStart){
		assert(qStart>=0 && qStart<=length && tStart>=0) : "The single-edit alternatives advance only past the first disagreement.";
		if(targetLength-tStart<length-qStart){return false;}
		for(int q=qStart, t=tStart; q<length; q++, t++){
			if(query[q]!=target[t]){return false;}
		}
		return true;
	}

	/** Actual overlap bases must agree in both contig orientations, independently of read selection. */
	private boolean overlapMatches(final Evidence e){
		final Connection c=e.connection;
		assert(c.gap>-Math.min(contigs.get(c.leftEnd/2).length(), contigs.get(c.rightEnd/2).length())) :
				"Discovery must exclude containment before overlap comparison.";
		final byte[] left=contigs.get(c.leftEnd/2).bases, right=contigs.get(c.rightEnd/2).bases;
		for(int i=0; i<-c.gap; i++){
			final int a=contigBase(left, (int)(left.length+c.gap+i), (c.leftEnd&1)==0);
			final int b=contigBase(right, i, (c.rightEnd&1)!=0);
			if(a<0 || a!=b){return false;}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Path spelling        ----------------*/
	/*--------------------------------------------------------------*/

	/** Chooses whole-gap votes, joins only linear paths, and preserves every unused contig unchanged. */
	public ArrayList<Contig> join(){
		if(finished){throw new IllegalStateException("Joining consumes one evidence snapshot.");}
		finished=true;
		final Evidence[] ends=new Evidence[2*contigs.size()];
		for(Evidence e : evidence){
			if(e.reason!=null){continue;}
			Variant best=null;
			long runner=0;
			for(Variant variant : e.variants){
				if(best==null || variant.count>best.count){runner=best==null ? 0 : best.count; best=variant;}
				else{runner=Math.max(runner, variant.count);}
			}
			e.support=best==null ? 0 : best.count;
			e.runner=runner;
			if(e.overflow){e.reason="too_many_gap_sequences";}
			else if(e.outerContradictory>=minimum && e.outerConfirmed/e.outerContradictory<ratio){
				e.reason="outer_flank_conflict";
			}
			else if(e.support<minimum){e.reason="insufficient_verified_reads";}
			else if(runner>0 && e.support/runner<ratio){e.reason="gap_sequence_ambiguity";}
			else{
				e.gap=best.bases;
				final Connection c=e.connection;
				assert(ends[c.leftEnd]==null && ends[c.rightEnd]==null) : "Unconflicted discovery must leave at most one edge per physical end.";
				ends[c.leftEnd]=ends[c.rightEnd]=e;
			}
		}
		final boolean[] used=new boolean[contigs.size()];
		final ArrayList<Contig> result=new ArrayList<Contig>();
		for(int id=0; id<contigs.size(); id++){
			if(used[id]){continue;}
			final boolean left=ends[2*id]!=null, right=ends[2*id+1]!=null;
			if(left==right){continue;}//Isolated contigs and closed components are preserved below.
			result.add(spellPath(id, left, ends, used, result.size()));
		}
		for(Evidence e : evidence){if(e.reason==null){e.reason="cycle_not_joined";}}
		for(int id=0; id<contigs.size(); id++){if(!used[id]){result.add(contigs.get(id));}}
		assert(result.size()==contigs.size()-joins) : "Every linear join must reduce the number of emitted contigs by exactly one.";
		return result;
	}

	/** Traverses a degree-one-ended component; orientation follows entry and exit sides. */
	private Contig spellPath(final int start, final boolean reverse, final Evidence[] ends,
			final boolean[] used, final int outputId){
		assert(!used[start]) : "Linear paths must own disjoint input contigs.";
		final ByteBuilder sequence=new ByteBuilder();
		final ByteBuilder name=new ByteBuilder();
		name.append(fillGaps ? "hollow_join_" : "hollow_scaffold_").append(outputId)
				.append(" members=").append(start).append(reverse ? '-' : '+');
		append(sequence, contigs.get(start).bases, reverse, 0);
		used[start]=true;
		int exit=2*start+(reverse ? 0 : 1);
		for(Evidence e=ends[exit]; e!=null; e=ends[exit]){
			final Connection c=e.connection;
			final boolean backward=exit==c.rightEnd;
			final int entry=backward ? c.leftEnd : c.rightEnd, next=entry/2;
			if(used[next]){throw new IllegalStateException("A degree-one-ended path unexpectedly revisited contig "+next);}
			final boolean nextReverse=(entry&1)!=0;
			if(fillGaps){
				if(c.gap>0){append(sequence, e.gap, backward, 0);}
				append(sequence, contigs.get(next).bases, nextReverse, (int)Math.max(0, -c.gap));
			}else{
				for(int i=0; i<SCAFFOLD_GAP; i++){sequence.append('N');}
				append(sequence, contigs.get(next).bases, nextReverse, 0);
			}
			name.append(',').append(next).append(nextReverse ? '-' : '+');
			used[next]=true;
			e.reason="joined";
			joins++;
			exit=entry^1;
		}
		return new Contig(sequence.toBytes(), name.toString(), outputId);
	}

	/** Appends an oriented suffix, preserving undefined bases outside the verified connection. */
	private static void append(final ByteBuilder target, final byte[] bases, final boolean reverse, final int skip){
		assert(skip>=0 && skip<bases.length) : "Overlap spelling must retain nonoverlapping contig sequence.";
		if((long)target.length()+bases.length-skip>Integer.MAX_VALUE-8){throw new IllegalArgumentException("Joined sequence exceeds array capacity.");}
		for(int i=skip; i<bases.length; i++){
			target.append(reverse ? AminoAcid.baseToComplementExtended[bases[bases.length-1-i]] : bases[i]);
		}
	}

	/** Returns a normalized base in one selected read-window orientation. */
	private int readBase(final byte[] bases, final int start, final int offset, final boolean reverse){
		final int base=AminoAcid.baseToNumber[bases[start+(reverse ? span-1-offset : offset)]];
		return reverse && base>=0 ? 3-base : base;
	}

	/** Returns a normalized base in one contig orientation; undefined bases do not match evidence. */
	private static int contigBase(final byte[] bases, final int offset, final boolean reverse){
		assert(offset>=0 && offset<bases.length) : "Anchor geometry must keep verified flanks within their contigs.";
		final int base=AminoAcid.baseToNumber[bases[reverse ? bases.length-1-offset : offset]];
		return reverse && base>=0 ? 3-base : base;
	}

	public int joins(){return joins;}
	public ArrayList<Evidence> evidence(){return new ArrayList<Evidence>(evidence);}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One retained geometry and its bounded set of whole-gap read votes. */
	public static final class Evidence {
		private Evidence(final Connection connection_, final HollowTipAnchors anchors){
			connection=connection_;
			sourceEnd=32+anchors.exitDistance(connection.leftAnchor(), (connection.leftEnd&1)==0);
			targetStart=sourceEnd+(int)connection.gap;
			assert(sourceEnd>=32 && sourceEnd<=anchors.contigLength(connection.leftAnchor())) :
					"The source arm must lie inside its contig; a deep overlap may place the target start before the read window.";
		}
		public String reason(){return reason;}
		public long support(){return support;}
		public long runner(){return runner;}
		public long repeated(){return repeated;}
		public long flankRejected(){return flankRejected;}
		public long outerConfirmed(){return outerConfirmed;}
		public long outerContradictory(){return outerContradictory;}
		public final Connection connection;
		private final int sourceEnd, targetStart;
		private final ArrayList<Variant> variants=new ArrayList<Variant>();
		private Variant lastVote;
		private long lastRead, support, runner, repeated, flankRejected;
		private long outerConfirmed, outerContradictory;
		private int lastFlank;
		private boolean overflow;
		private byte[] gap;
		private String reason;
	}

	/** One exact missing sequence and the number of qualifying read records containing it. */
	private static final class Variant {
		Variant(final byte[] bases_){bases=bases_;}
		final byte[] bases;
		long count;
	}

	/** One directional arm-pair selector; allocated only per retained candidate, not per read. */
	private static final class Selector {
		Selector(final long second_, final boolean reverse_, final Evidence evidence_, final Selector next_){
			second=second_;
			reverse=reverse_;
			evidence=evidence_;
			next=next_;
		}
		final long second;
		final boolean reverse;
		final Evidence evidence;
		final Selector next;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<Contig> contigs;
	private final ArrayList<Evidence> evidence=new ArrayList<Evidence>();
	private final LongObjectMap<Selector> selectors=new LongObjectMap<Selector>(Selector.class);
	private final int span, minimum, ratio;
	private final int flankBases, flankMin;
	private final boolean fillGaps;
	private final byte[] queryFlank, targetFlank;
	private final long[] delayed;
	private long readSerial;
	private int joins;
	private boolean finished;
	private static final int MAX_VARIANTS=16;
	public static final int SCAFFOLD_GAP=10;
	private static final byte[] ALPHABET={'A','C','G','T'};
}
