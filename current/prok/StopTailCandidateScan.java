package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import dna.GeneticCode;
import structures.IntList;

/** Enumerates real stops independently of native start admission in all six frames.
 * Default queries contain exactly the last 200 unambiguous sense codons before a stop.
 * Explicit short-query APIs emit the whole uninterrupted interval when it is shorter.
 * Coordinates are zero-based on the original contig; strand 0 is plus, 1 minus.
 * No artificial contig-edge or ambiguity-boundary stop is emitted.
 * @author Keqing
 */
public final class StopTailCandidateScan{

	private StopTailCandidateScan(){}

	/** Borrows, but never changes, bases/code. The caller must not mutate bases during scanning.
	 * The result owns its query arrays and may be shared after publication. Each scan is scoped
	 * to this contig and exact code object, including custom codes whose numeric ID is zero. */
	public static Scan scan(final String contig, final byte[] bases, final GeneticCode code){
		return scan(contig, bases, code, TAIL_LENGTH);
	}

	/** Explicit research scan with 40..200 minimum residues. Longer intervals still
	 * yield only their last 200 residues; no padding or annotated start is used.
	 * Emitting a short query does not authorize using a 200-aa classifier for it. */
	public static Scan scan(final String contig, final byte[] bases, final GeneticCode code, final int minimumTailLength){
		if(minimumTailLength<40 || minimumTailLength>TAIL_LENGTH){
			throw new IllegalArgumentException("Minimum stop-tail length must be in [40,200]");
		}
		return scanConfigured(contig, bases, code, minimumTailLength, minimumTailLength<100 ? 31 : 15, TAIL_LENGTH);
	}

	/** Select separately calibrated short bins; the original 200-aa bin is always on.
	 * Enabling 100–124 and 150–199 does not enable the intervening bin. */
	public static Scan scanBins(final String contig, final byte[] bases, final GeneticCode code,
			final boolean bin100, final boolean bin125, final boolean bin150){
		return scanBins(contig, bases, code, false, bin100, bin125, bin150);
	}

	/** The independently calibrated40..99 interval has its own bit; old bin selection is unchanged. */
	public static Scan scanBins(final String contig, final byte[] bases, final GeneticCode code,
			final boolean bin40, final boolean bin100, final boolean bin125, final boolean bin150){
		return scanBins(contig, bases, code, bin40, bin100, bin125, bin150, 40);
	}

	/** Raising the scanner minimum never changes the calibrated40..99 child's upper bound. */
	public static Scan scanBins(final String contig, final byte[] bases, final GeneticCode code,
			final boolean bin40, final boolean bin100, final boolean bin125, final boolean bin150, final int minimumShortLength){
		return scanBins(contig, bases, code, bin40, bin100, bin125, bin150, minimumShortLength, false);
	}

	/** Optional earliest-start boundary within each uninterrupted frame interval.
	 * Uses the supplied code's start set; no start means no query. Contig edges and
	 * ambiguous codons reset the interval. Long queries retain their last200 residues. */
	public static Scan scanBins(final String contig, final byte[] bases, final GeneticCode code,
			final boolean bin40, final boolean bin100, final boolean bin125, final boolean bin150,
			final int minimumShortLength, final boolean firstStart){
		if(minimumShortLength<40 || minimumShortLength>99){
			throw new IllegalArgumentException("Short-bin scanner minimum must be in [40,99]");
		}
		final int minimum=bin40 ? minimumShortLength : bin100 ? 100 : bin125 ? 125 : bin150 ? 150 : TAIL_LENGTH;
		return scanConfigured(contig, bases, code, minimum,
			8|(bin40 ? 16 : 0)|(bin100 ? 1 : 0)|(bin125 ? 2 : 0)|(bin150 ? 4 : 0), TAIL_LENGTH, firstStart);
	}

	/** Research-only fixed100 experiment: every eligible interval contributes its last100 residues. */
	public static Scan scanFixed100(final String contig, final byte[] bases, final GeneticCode code){
		return scanConfigured(contig, bases, code, 100, 1, 100);
	}

	/** Emit every selected interval of at least minimum residues, keeping its terminal
	 * min(selected length, maximum) residues. Raw upstreamCodons and the unclipped
	 * queryCodons remain distinct; firstStart uses the supplied code's earliest start.
	 * This API selects representation only, not a classifier or calibrated gate. */
	public static Scan scanCapped(final String contig, final byte[] bases, final GeneticCode code,
			final int minimumTailLength, final int maximumTailLength, final boolean firstStart){
		if(minimumTailLength<40 || minimumTailLength>maximumTailLength || maximumTailLength>TAIL_LENGTH){
			throw new IllegalArgumentException("Capped stop-tail scan requires 40<=minimum<=maximum<=200");
		}
		return scanConfigured(contig, bases, code, minimumTailLength, 31, maximumTailLength, firstStart);
	}

	private static Scan scanConfigured(final String contig, final byte[] bases, final GeneticCode code,
			final int minimumTailLength, final int lengthBins, final int maximumTailLength){
		return scanConfigured(contig, bases, code, minimumTailLength, lengthBins, maximumTailLength, false);
	}
	private static Scan scanConfigured(final String contig, final byte[] bases, final GeneticCode code,
			final int minimumTailLength, final int lengthBins, final int maximumTailLength, final boolean firstStart){
		return scanConfigured(contig, bases, code, minimumTailLength, lengthBins, maximumTailLength, firstStart, null);
	}

	/** Opt-in edge observation in the SAME six translation passes as ordinary capped stops.
	 * Edge queries never require a start and always retain40..100 residues, independently
	 * of ordinary stop policy. A null scope disables the extension exactly. */
	public static Scan scanCappedWithFragments(final String contig, final byte[] bases, final GeneticCode code,
			final int minimumTailLength, final int maximumTailLength, final boolean firstStart,
			final ContigEdgeFragment.Scope scope){
		if(minimumTailLength<40 || minimumTailLength>maximumTailLength || maximumTailLength>TAIL_LENGTH){
			throw new IllegalArgumentException("Capped stop-tail scan requires 40<=minimum<=maximum<=200");
		}
		return scanConfigured(contig, bases, code, minimumTailLength, 31, maximumTailLength, firstStart, scope);
	}

	private static Scan scanConfigured(final String contig, final byte[] bases, final GeneticCode code,
			final int minimumTailLength, final int lengthBins, final int maximumTailLength, final boolean firstStart,
			final ContigEdgeFragment.Scope scope){
		return scanConfigured(contig,bases,code,minimumTailLength,lengthBins,maximumTailLength,firstStart,scope,false);
	}

	/** Explicit caller GC view. Ordinary scanners retain their original work path. */
	public static Scan scanCappedWithGc(final String contig,final byte[] bases,final GeneticCode code,
			final int minimumTailLength,final boolean firstStart){
		if(minimumTailLength<40 || minimumTailLength>99){
			throw new IllegalArgumentException("GC caller minimum must be40..99");
		}
		return scanConfigured(contig,bases,code,minimumTailLength,31,100,firstStart,null,true);
	}

	private static Scan scanConfigured(final String contig, final byte[] bases, final GeneticCode code,
			final int minimumTailLength, final int lengthBins, final int maximumTailLength, final boolean firstStart,
			final ContigEdgeFragment.Scope scope,final boolean measureGc){
		if(contig==null || contig.isEmpty() || bases==null || code==null){
			throw new IllegalArgumentException("Stop-tail scan requires a contig name, bases and genetic code");
		}
		final ArrayList<Candidate> candidates=new ArrayList<Candidate>();
		final ArrayList<ContigEdgeFragment> fragments=scope==null ? null : new ArrayList<ContigEdgeFragment>(6);
		final boolean edgeEnabled=scope!=null && scope.topology!=ContigEdgeFragment.Topology.CIRCULAR;
		final IntList starts=edgeEnabled ? new IntList() : null;
		long edgeShort=0, edgeEmpty=0;
		final long[] stopCounts=new long[3];
		long ambiguous=0, trailing=0, shortQueries=0, noStartStops=0, edgeQueries=0, ambiguityQueries=0;
		long atLeastMaximum=0, cappedQueries=0;
		final int ringLength=edgeEnabled ? Math.max(maximumTailLength, ContigEdgeFragment.MAXIMUM) : maximumTailLength;
		final byte[] ring=new byte[ringLength];
		for(int strand=0; strand<2; strand++){
			for(int frame=0; frame<3; frame++){
				int sense=0, next=0, first=-1, boundary=BOUNDARY_EDGE;
				if(edgeEnabled){starts.clear();}
				for(int p=frame; p<=bases.length-3; p+=3){
					final int left=strand==0 ? p : bases.length-p-3;
					int codon=GeneticCode.codon(bases, left);
					if(codon<0){
						sense=next=0; first=-1; boundary=BOUNDARY_AMBIGUITY; ambiguous++;
						if(edgeEnabled){starts.clear();}
						continue;
					}
					if(strand==1){codon=reverseComplementCodon(codon);}
					final byte aa=code.aminoAcid(codon);
					if(aa=='*'){
						stopCounts[sense<100 ? 0 : sense<TAIL_LENGTH ? 1 : 2]++;
						if(sense>=maximumTailLength){atLeastMaximum++;}
						final int queryCodons=firstStart ? (first<0 ? 0 : sense-first) : sense;
						final int rawLength=Math.min(sense, maximumTailLength);
						if(firstStart && first<0 && rawLength>=minimumTailLength && (lengthBins&lengthBit(rawLength))!=0){noStartStops++;}
						final int length=Math.min(queryCodons, maximumTailLength);
						if(length>=minimumTailLength && (lengthBins&lengthBit(length))!=0){
							final byte[] tail=new byte[length];
							if(length<TAIL_LENGTH){shortQueries++;}
							if(length==maximumTailLength){cappedQueries++;}
							copySuffix(ring, next, tail);
							final int gc=measureGc ? senseGcCount(bases,left,strand,length) : 0;
							final int shortGc=measureGc && length>=100 ? senseGcCount(bases,left,strand,99) : gc;
							candidates.add(new Candidate(left, strand, frame, sense, queryCodons, tail, maximumTailLength, boundary,
								measureGc,gc,shortGc));
							if(boundary==BOUNDARY_EDGE){edgeQueries++;}
							else if(boundary==BOUNDARY_AMBIGUITY){ambiguityQueries++;}
						}
						sense=next=0; first=-1; boundary=BOUNDARY_STOP;
						if(edgeEnabled){starts.clear();}
					}else{
						if(firstStart && first<0 && code.isStart(codon)){first=sense;}
						if(edgeEnabled && code.isStart(codon)){starts.add(p);}
						ring[next++]=aa;
						if(next==ringLength){next=0;}
						sense++;
					}
				}
				if(sense>0){trailing++;}
				if(edgeEnabled){
					if(sense==0){edgeEmpty++;}
					else if(sense<ContigEdgeFragment.MINIMUM){edgeShort++;}
					else{
						final byte[] peptide=new byte[Math.min(sense, ContigEdgeFragment.MAXIMUM)];
						copySuffix(ring, next, peptide);
						fragments.add(new ContigEdgeFragment(scope, contig, code, bases.length, strand, frame,
							sense, boundary, peptide, starts.toArray()));
					}
				}
			}
		}
		assert(firstStart || cappedQueries==((lengthBins&lengthBit(maximumTailLength))==0 ? 0 : atLeastMaximum))
			: "Whole-interval scans must emit every real stop reaching the enabled cap; cap="+maximumTailLength+
			", physical="+atLeastMaximum+", emitted="+cappedQueries;
		return new Scan(contig, bases.length, code, candidates, stopCounts, ambiguous, trailing,
			minimumTailLength, lengthBins, shortQueries, maximumTailLength, firstStart, noStartStops,
			edgeQueries, ambiguityQueries,
			scope, fragments, edgeShort, edgeEmpty);
	}
	private static void copySuffix(final byte[] ring, final int next, final byte[] suffix){
		assert(suffix.length<=ring.length && next>=0 && next<ring.length) :
			"The shared translation ring must retain every requested suffix codon";
		final int begin=(next-suffix.length+ring.length)%ring.length, part=Math.min(suffix.length, ring.length-begin);
		System.arraycopy(ring, begin, suffix, 0, part);
		System.arraycopy(ring, 0, suffix, part, suffix.length-part);
	}
	static int lengthBit(final int length){return length<100 ? 16 : length<125 ? 1 : length<150 ? 2 : length<200 ? 4 : 8;}

	/** Count sense nucleotides upstream of a real stop; terminal stop is excluded.
	 * The minus interval is on the original contig; GC is complement-invariant. */
	public static int senseGcCount(final byte[] bases,final int stopLeft,final int strand,final int aminoAcids){
		if(bases==null || stopLeft<0 || stopLeft>bases.length-3 || strand<0 || strand>1
				|| aminoAcids<1 || aminoAcids>TAIL_LENGTH){
			throw new IllegalArgumentException("Sense GC requires a valid stop coordinate, strand and tail length");
		}
		final int start=strand==0 ? stopLeft-3*aminoAcids : stopLeft+3,end=start+3*aminoAcids;
		if(start<0 || end>bases.length){throw new IllegalArgumentException("Sense GC interval exceeds its contig");}
		int gc=0;
		for(int i=start;i<end;i++){
			final int b=bases[i]&0xDF;
			if(b=='G' || b=='C'){gc++;}
			else if(b!='A' && b!='T' && b!='U'){
				throw new IllegalArgumentException("Ambiguous nucleotide in measured sense tail");
			}
		}
		return gc;
	}

	/** Reverses the three packed A=0,C=1,G=2,T=3 bases and complements each. */
	private static int reverseComplementCodon(final int codon){
		assert(codon>=0 && codon<64) : "Only a complete GeneticCode codon may be reverse complemented";
		return ((3^(codon&3))<<4)|((3^((codon>>>2)&3))<<2)|(3^(codon>>>4));
	}

	/** A stop key is unique only within a Scan's contig and exact genetic code. */
	public static long stopKey(final int stopLeft, final int strand){
		if(stopLeft<0 || strand<0 || strand>1){
			throw new IllegalArgumentException("Stop key requires nonnegative genomic left coordinate and strand 0 or 1");
		}
		return ((long)stopLeft<<1)|strand;
	}

	/** Immutable metadata with a privately owned query; raw interval and selected query span are distinct. */
	public static final class Candidate{
		private Candidate(final int stopLeft_, final int strand_, final int frame_, final int upstreamCodons_,
				final int queryCodons_, final byte[] tail_, final int maximum_, final int boundary_,
				final boolean gcAvailable_,final int gc_,final int shortGc_){
			assert(queryCodons_<=upstreamCodons_ && tail_.length>=40 && tail_.length<=TAIL_LENGTH
					&& tail_.length==Math.min(queryCodons_, maximum_))
				: "A query is the selected sense suffix within the uninterrupted interval, capped but never padded";
			stopLeft=stopLeft_; strand=strand_; frame=frame_; upstreamCodons=upstreamCodons_;
			queryCodons=queryCodons_; tail=tail_;
			intervalBoundary=boundary_;
			gcAvailable=gcAvailable_;gc=gc_;shortGc=shortGc_;
		}
		public long key(){return stopKey(stopLeft, strand);}
		/** Copies for external callers; the package-local adapter may borrow tail read-only. */
		public byte[] copyTail(){return tail.clone();}
		public double tailGc(){return gcAvailable ? gc/(3.0*tail.length) : 0;}
		public double shortTailGc(){return gcAvailable ? shortGc/(3.0*Math.min(tail.length,99)) : 0;}
		public final boolean gcAvailable;
		private final int gc,shortGc;
		public final int stopLeft, strand, frame, upstreamCodons;
		/** Selected start-to-stop span before clipping; equals upstreamCodons in the default mode. */
		public final int queryCodons;
		/** Upstream reset that bounds the raw interval: edge, ambiguity or a real stop. */
		public final int intervalBoundary;
		final byte[] tail;
	}

	/** Stop counts include zero-length intervals; ambiguity counts are frame-codon counts,
	 * not distinct ambiguous nucleotides. Trailing intervals contain sense codons but no real stop. */
	public static final class Scan{
		private Scan(final String contig_, final int length_, final GeneticCode code_, final ArrayList<Candidate> candidates_,
				final long[] counts, final long ambiguous_, final long trailing_, final int minimum_,
				final int bins_, final long shortQueries_, final int maximum_,
				final boolean firstStart_, final long noStartStops_, final long edgeQueries_, final long ambiguityQueries_,
				final ContigEdgeFragment.Scope scope_, final ArrayList<ContigEdgeFragment> fragments_,
				final long edgeShort_, final long edgeEmpty_){
			contig=contig_; length=length_; code=code_;
			candidates=Collections.unmodifiableList(candidates_);
			under100=counts[0]; from100To199=counts[1]; atLeast200=counts[2];
			ambiguousCodons=ambiguous_; unterminatedIntervals=trailing_;
			minimumTailLength=minimum_; maximumTailLength=maximum_; lengthBins=bins_; shortQueries=shortQueries_;
			firstStart=firstStart_;
			noStartStops=noStartStops_; edgeBoundedQueries=edgeQueries_; ambiguityBoundedQueries=ambiguityQueries_;
			edgeScope=scope_;
			fragments=fragments_==null ? Collections.<ContigEdgeFragment>emptyList() : Collections.unmodifiableList(fragments_);
			edgeShortIntervals=edgeShort_; edgeEmptyFrames=edgeEmpty_;
			edgeCircularFrames=scope_!=null && scope_.topology==ContigEdgeFragment.Topology.CIRCULAR ? 6 : 0;
			assert(scope_==null || fragments.size()+edgeShortIntervals+edgeEmptyFrames+edgeCircularFrames==6) :
				"Every oriented frame must be retained or explicitly excluded once; no annotation filters edge discovery";
		}
		public boolean emitsLength(final int length){
			return length>=minimumTailLength && length<=maximumTailLength && (lengthBins&lengthBit(length))!=0;
		}
		public long realStops(){return under100+from100To199+atLeast200;}
		public final String contig;
		public final int length;
		public final int minimumTailLength, maximumTailLength, lengthBins;
		public final long shortQueries;
		public final boolean firstStart;
		/** Otherwise-eligible whole-window candidates lost solely because the interval has no start. */
		public final long noStartStops, edgeBoundedQueries, ambiguityBoundedQueries;
		public final GeneticCode code;
		public final List<Candidate> candidates;
		public final ContigEdgeFragment.Scope edgeScope;
		public final List<ContigEdgeFragment> fragments;
		public final long edgeShortIntervals, edgeEmptyFrames, edgeCircularFrames;
		public final long under100, from100To199, atLeast200, ambiguousCodons, unterminatedIntervals;
	}

	public static final int TAIL_LENGTH=200;
	public static final int BOUNDARY_EDGE=0, BOUNDARY_AMBIGUITY=1, BOUNDARY_STOP=2;
}
