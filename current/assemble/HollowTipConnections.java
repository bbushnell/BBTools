package assemble;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import dna.AminoAcid;
import map.LongObjectMap;
import ukmer.HollowAnchorTable;

/**
 * Experimental read-supported connections between assembly-unique tip anchors.
 * Counts literal outer32bp pairs, retaining every partner of an indexed anchor,
 * including partners absent from the assembly. Discovery requires a sole
 * supported partner in both directions and records conflicting end geometries.
 * It does not merge contigs, reconstruct gaps, or change Tadpole defaults.
 * @author Brian Bushnell, Fischl
 */
public final class HollowTipConnections {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds immutable assembly anchors and an initially empty, single-span evidence table. */
	public HollowTipConnections(final ArrayList<Contig> contigs, final int span_,
			final int tipBases, final int minSupport_){
		this(contigs, span_, tipBases, minSupport_, 0, tipBases-32, 0);
	}

	/** Seed distances exclude the32bp arm; sliding evidence extends beyond the seed band. */
	public HollowTipConnections(final ArrayList<Contig> contigs, final int span_,
			final int tipBases, final int minSupport_, final int tipMin_, final int tipMax_, final int slide_){
		if(contigs==null || contigs.size()>Integer.MAX_VALUE/2 || span_<=64 || minSupport_<1){
			throw new IllegalArgumentException("Connections require representable contig ends, span>64 and positive support.");
		}
		if(tipBases<32 || tipMin_<0 || tipMax_<tipMin_ || tipMax_>tipBases-32 || slide_<0
				|| (long)tipBases+slide_>Integer.MAX_VALUE){
			throw new IllegalArgumentException("Require0<=tipmin<=tipmax<=tipbases-32 and a nonnegative representable slide window.");
		}
		span=span_;
		minSupport=minSupport_;
		tipMin=tipMin_;
		tipMax=tipMax_;
		slide=slide_;
		sequences=contigs;
		contigCount=contigs.size();
		anchors=new HollowTipAnchors(contigs, tipBases+slide, slide>0);
		table=new HollowAnchorTable(span, 1024);
		delayedWords=new long[span-32];
		delayedEligible=new boolean[span-32];
	}

	/*--------------------------------------------------------------*/
	/*----------------       Evidence loading       ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts one read once; only fully defined spans contribute, even in the masked interior. */
	public void addRead(final byte[] bases){
		if(discovered){throw new IllegalStateException("Evidence cannot change after connection discovery.");}
		if(bases==null){throw new IllegalArgumentException("Read bases must not be null.");}
		reads++;
		readBases+=bases.length;
		if(bases.length<span){return;}
		long forward=0, reverse=0;
		int valid=0, slot=0;
		for(int pos=0; pos<bases.length; pos++){
			final int base=AminoAcid.baseToNumber[bases[pos]];
			if(base<0){forward=reverse=0; valid=0;}
			else{
				forward=(forward<<2)|base;
				reverse=(reverse>>>2)|((3L-base)<<62);
				valid++;
			}
			final long first=delayedWords[slot];
			final boolean firstEligible=delayedEligible[slot];
			final long key=Long.compareUnsigned(forward, reverse)>=0 ? forward : reverse;
			final boolean lastEligible=valid>=32 && (slide>0 ? anchors.containsCanonical(key) : anchors.lookupCanonical(key)>=0);
			delayedWords[slot]=forward;
			delayedEligible[slot]=lastEligible;
			if(++slot==delayedWords.length){slot=0;}
			if(valid<span){continue;}
			windows++;
			if(firstEligible || lastEligible){
				table.increment(first, forward);
				retainedWindows++;
			}
		}
		assert(retainedWindows<=windows) : "Tip-restricted pair evidence must be a subset of complete observed windows.";
	}

	/*--------------------------------------------------------------*/
	/*----------------         Connections          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns deterministic candidate geometries. A conflicted candidate is
	 * evidence only, never permission to join. Support is the maximum individual
	 * pair count, not a sum of correlated overlapping windows or distinct reads.
	 */
	public ArrayList<Connection> discover(){
		if(discovered){throw new IllegalStateException("Connection discovery is a one-time snapshot operation.");}
		discovered=true;
		table.freeze();
		final LongObjectMap<Connection> byEnds=new LongObjectMap<Connection>(Connection.class);
		final ArrayList<Connection> result=new ArrayList<Connection>();
		final SolePartner outgoing=new SolePartner(minSupport), incoming=new SolePartner(minSupport);
		for(int left=0; left<anchors.size(); left++){
			if(!anchors.isUnique(left)){continue;}
			for(int direction=0; direction<2; direction++){
				final boolean leftReverse=direction==1;
				if(!seedDistance(anchors.exitDistance(left, leftReverse))){continue;}
				final long word=leftReverse ? HollowAnchorTable.reverseComplement(anchors.forwardWord(left)) : anchors.forwardWord(left);
				outgoing.reset();
				table.forEachPartner(word, outgoing);
				if(outgoing.supported!=1){continue;}
				final int right=anchors.lookup(outgoing.partner);
				if(right<0 || anchors.contigIndex(right)==anchors.contigIndex(left)){continue;}
				final boolean rightReverse=anchors.reverse(right, outgoing.partner);
				if(!seedDistance(anchors.entryDistance(right, rightReverse))){continue;}
				final int leftEnd=2*anchors.contigIndex(left)+(leftReverse ? 0 : 1);
				final int rightEnd=2*anchors.contigIndex(right)+(rightReverse ? 1 : 0);
				if(leftEnd>rightEnd){continue;}//The reverse traversal represents the same physical end pair.
				incoming.reset();
				table.forEachPartner(HollowAnchorTable.reverseComplement(outgoing.partner), incoming);
				if(incoming.supported!=1 || incoming.partner!=HollowAnchorTable.reverseComplement(word)){continue;}
				assert(incoming.count==outgoing.count) : "Reciprocal endpoint lookup must share the same canonical pair count.";
				final long gap=anchors.gap(left, leftReverse, right, rightReverse, span);
				if(gap<=-Math.min(anchors.contigLength(left), anchors.contigLength(right))){continue;}
				final long pair=((long)leftEnd<<32)|(rightEnd&0xffffffffL);
				final Connection head=byEnds.get(pair);
				Connection match=head;
				while(match!=null && match.gap!=gap){match=match.next;}
				if(match==null){
					match=new Connection(leftEnd, rightEnd, gap, left, right, word, outgoing.partner, outgoing.count);
					match.next=head;
					byEnds.put(pair, match);
					result.add(match);
				}else{match.addEvidence(left, right, word, outgoing.partner, outgoing.count);}
			}
		}
		final int[] degree=new int[2*contigCount];
		for(Connection connection : result){degree[connection.leftEnd]++; degree[connection.rightEnd]++;}
		for(Connection connection : result){
			connection.conflicted=degree[connection.leftEnd]>1 || degree[connection.rightEnd]>1;
			connection.next=null;//Lookup-only chains must not leak into the result graph.
			if(slide>0){checkSliding(connection, outgoing, incoming);}
		}
		Collections.sort(result, ORDER);
		return result;
	}

	public long reads(){return reads;}
	public long readBases(){return readBases;}
	public long windows(){return windows;}
	public long retainedWindows(){return retainedWindows;}
	public int distinctPairs(){return table.size();}
	public long pairArrayBytes(){return table.arrayBytes();}

	/** Seed exclusion is not evidence exclusion; terminal/repeated probe words remain countable. */
	private boolean seedDistance(final int distance){return distance>=tipMin && distance<=tipMax;}

	/**
	 * Audits the strongest deterministic seed, without retrying a different seed
	 * after failure. Both known contig arms slide together at every integer offset.
	 * Any supported alternative partner vetoes the join, including unmapped or
	 * nonunique partners. Missing expected support is reported separately.
	 */
	private void checkSliding(final Connection c, final SolePartner outgoing, final SolePartner incoming){
		assert(slide>0 && table.frozen()) : "Sliding verification consumes the complete immutable read evidence table.";
		final byte[] left=sequences.get(c.leftEnd/2).bases, right=sequences.get(c.rightEnd/2).bases;
		final boolean lr=(c.leftEnd&1)==0, rr=(c.rightEnd&1)!=0;
		final int ls=left.length-32-anchors.exitDistance(c.leftAnchor, lr);
		final int rs=anchors.entryDistance(c.rightAnchor, rr);
		if(ls<slide || rs<slide || left.length-32-ls<slide || right.length-32-rs<slide){
			c.slideReason="sliding_margin";
			return;
		}
		c.slideMinCount=Integer.MAX_VALUE;
		for(int offset=-slide; offset<=slide; offset++){
			long a=0, b=0;
			boolean defined=true;
			for(int i=0; i<32; i++){
				final int x=orientedBase(left, ls+offset+i, lr), y=orientedBase(right, rs+offset+i, rr);
				if(x<0 || y<0){defined=false; break;}
				a=(a<<2)|x;
				b=(b<<2)|y;
			}
			if(!defined){
				if(c.slideReason==null){c.slideReason="sliding_undefined"; c.slideFailureOffset=offset;}
				c.slideMinCount=0;
				continue;
			}
			final int count=table.get(a, b);
			c.slideChecked++;
			c.slideMinCount=Math.min(c.slideMinCount, count);
			outgoing.reset(); incoming.reset();
			table.forEachPartner(a, outgoing);
			table.forEachPartner(HollowAnchorTable.reverseComplement(b), incoming);
			final boolean branch=outgoing.supported>1 || incoming.supported>1
					|| (outgoing.supported==1 && outgoing.partner!=b)
					|| (incoming.supported==1 && incoming.partner!=HollowAnchorTable.reverseComplement(a));
			if(branch){
				if(!"sliding_branch".equals(c.slideReason)){c.slideReason="sliding_branch"; c.slideFailureOffset=offset;}
			}else if(count<minSupport && c.slideReason==null){
				c.slideReason="sliding_low_support";
				c.slideFailureOffset=offset;
			}
		}
	}

	/** Reads an oriented contig coordinate; callers verify both sliding bounds before entry. */
	private static int orientedBase(final byte[] bases, final int position, final boolean reverse){
		assert(position>=0 && position<bases.length) : "Sliding arms must stay wholly inside their known contig sequence.";
		final int base=AminoAcid.baseToNumber[bases[reverse ? bases.length-1-position : position]];
		return reverse && base>=0 ? 3-base : base;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained result for one physical end pair and one signed gap/overlap. */
	public static final class Connection {

		private Connection(final int leftEnd_, final int rightEnd_, final long gap_,
				final int left, final int right, final long first, final long second, final int count){
			assert(leftEnd_<rightEnd_ && count>0) : "Candidate canonicalization requires ordered distinct ends and positive evidence.";
			leftEnd=leftEnd_;
			rightEnd=rightEnd_;
			gap=gap_;
			addEvidence(left, right, first, second, count);
		}

		/** Overlapping anchor observations are correlated; retain their maximum count, never their sum. */
		private void addEvidence(final int left, final int right, final long first, final long second, final int count){
			assert(count>0) : "Only supported observed pairs can contribute candidate evidence.";
			anchorPairs++;
			if(count>support || (count==support && (left<leftAnchor || (left==leftAnchor && right<rightAnchor)))){
				support=count;
				leftAnchor=left;
				rightAnchor=right;
				firstWord=first;
				secondWord=second;
			}
		}

		public int support(){return support;}
		public int anchorPairs(){return anchorPairs;}
		public int leftAnchor(){return leftAnchor;}
		public int rightAnchor(){return rightAnchor;}
		public long firstWord(){return firstWord;}
		public long secondWord(){return secondWord;}
		public boolean conflicted(){return conflicted;}
		public String slideReason(){return slideReason;}
		public int slideChecked(){return slideChecked;}
		public int slideMinCount(){return slideMinCount;}
		public int slideFailureOffset(){return slideFailureOffset;}

		public final int leftEnd, rightEnd;
		public final long gap;
		private int support, anchorPairs, leftAnchor, rightAnchor;
		private long firstWord, secondWord;
		private boolean conflicted;
		private String slideReason;
		private int slideChecked, slideMinCount, slideFailureOffset;
		private Connection next;
	}

	/** Reusable callback counts every sufficiently supported partner, including unmapped ones. */
	private static final class SolePartner implements HollowAnchorTable.PartnerConsumer {
		SolePartner(final int minimum_){minimum=minimum_;}
		void reset(){supported=0; count=0; partner=0;}
		@Override
		public void accept(final long partner_, final int count_){
			if(count_<minimum){return;}
			if(supported<2){supported++;}
			partner=partner_;
			count=count_;
		}
		final int minimum;
		int supported, count;
		long partner;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final int span, minSupport;
	public final int tipMin, tipMax, slide;
	public final HollowTipAnchors anchors;
	private final ArrayList<Contig> sequences;
	private final int contigCount;
	private final HollowAnchorTable table;
	private final long[] delayedWords;
	private final boolean[] delayedEligible;
	private long reads, readBases, windows, retainedWindows;
	private boolean discovered;
	private static final Comparator<Connection> ORDER=new Comparator<Connection>(){
		@Override
		public int compare(final Connection a, final Connection b){
			int cmp=Integer.compare(a.leftEnd, b.leftEnd);
			if(cmp==0){cmp=Integer.compare(a.rightEnd, b.rightEnd);}
			return cmp==0 ? Long.compare(a.gap, b.gap) : cmp;
		}
	};
}
