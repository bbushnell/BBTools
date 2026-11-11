package assemble;

import java.util.ArrayList;

import dna.AminoAcid;
import map.LongIntMap;
import structures.IntList;
import structures.LongList;
import ukmer.HollowAnchorTable;

/**
 * Immutable snapshot of assembly-wide unique 32 bp anchors near contig ends.
 * Candidate discovery is tip-restricted; uniqueness verification scans every
 * contig position, including interiors and reverse-complement duplicates.
 * Anchor identity is canonical, while stored coordinates retain the first
 * literal tip occurrence. Palindromic anchors are excluded from uniqueness
 * because one word cannot orient a connection. Rebuild after changing contig
 * sequences or their order; indices refer to the supplied list order, not
 * mutable Contig.id values.
 * @author Brian Bushnell, Fischl
 */
public final class HollowTipAnchors{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Indexes whole 32 bp words contained in either terminal tipBases interval. */
	public HollowTipAnchors(final ArrayList<Contig> contigs, final int tipBases_){
		this(contigs, tipBases_, false);
	}

	/** Optionally retains palindromic candidates for branch evidence, never as unique oriented seeds. */
	public HollowTipAnchors(final ArrayList<Contig> contigs, final int tipBases_, final boolean retainPalindromes_){
		if(contigs==null || tipBases_<32){
			throw new IllegalArgumentException("Tip indexing requires contigs and tipBases>=32.");
		}
		tipBases=tipBases_;
		retainPalindromes=retainPalindromes_;
		for(int id=0; id<contigs.size(); id++){
			final Contig contig=contigs.get(id);
			if(contig==null || contig.bases==null){throw new IllegalArgumentException("Null contig/bases at index "+id);}
			final int last=contig.length()-32;
			if(last<0){continue;}
			final int leftLast=Math.min(last, tipBases-32);
			addCandidates(contig.bases, id, 0, leftLast);
			// Avoid visiting the same positions twice when the two tip bands overlap.
			final int rightFirst=Math.max(leftLast+1, contig.length()-tipBases);
			if(rightFirst<=last){addCandidates(contig.bases, id, rightFirst, last);}
		}
		occurrences=new byte[words.size];
		for(Contig contig : contigs){countOccurrences(contig.bases);}
		int unique=0;
		for(int i=0; i<occurrences.length; i++){if(isUnique(i)){unique++;}}
		uniqueCount=unique;
		assert(words.size==positions.size && words.size==owners.size && words.size==lengths.size) :
				"Each candidate anchor must retain one consistent position/owner/length record.";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Queries            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a snapshot anchor index only when the word occurs exactly once in the assembly. */
	public int lookup(final long literalWord){
		return lookupCanonical(canonical(literalWord));
	}

	/** Package-local rolling scanners already maintain both strands; key must be their unsigned max. */
	int lookupCanonical(final long key){
		final int index=map.get(key)-1;
		return index>=0 && isUnique(index) ? index : -1;
	}

	/** Evidence membership ignores uniqueness: a repeated shifted arm must not hide branches. */
	boolean containsCanonical(final long key){return map.get(key)>0;}

	/** Number of distinct candidate canonical words seen in a tip interval. */
	public int size(){return words.size;}

	/** Number of nonpalindromic candidate words that occur exactly once across all contig positions. */
	public int uniqueCount(){return uniqueCount;}

	/** Whether the indexed candidate survives global uniqueness and palindrome filters. */
	public boolean isUnique(final int index){
		checkIndex(index);
		return occurrences[index]==1 && words.get(index)!=HollowAnchorTable.reverseComplement(words.get(index));
	}

	/** Literal 32 bp word as it appears in its owning contig, before canonicalization. */
	public long forwardWord(final int index){checkIndex(index); return words.get(index);}

	/** Index of the owning contig in the list used to build this snapshot. */
	public int contigIndex(final int index){checkIndex(index); return owners.get(index);}

	/** Zero-based start coordinate of the literal word on its stored forward strand. */
	public int start(final int index){checkIndex(index); return positions.get(index);}

	/** Length of the owning contig when this immutable snapshot was built. */
	public int contigLength(final int index){checkIndex(index); return lengths.get(index);}

	/** Whether this observed arm is reverse-complemented relative to its stored contig orientation. */
	public boolean reverse(final int index, final long observed){
		checkIndex(index);
		final long forward=words.get(index);
		if(observed==forward){return false;}
		if(observed==HollowAnchorTable.reverseComplement(forward)){return true;}
		throw new IllegalArgumentException("Observed word does not identify anchor "+index);
	}

	/** Bases from the 32 bp arm to the exit end in the indicated contig strand. */
	public int exitDistance(final int index, final boolean reverse){
		checkIndex(index);
		return reverse ? positions.get(index) : lengths.get(index)-positions.get(index)-32;
	}

	/** Bases from the entry end to the 32 bp arm in the indicated contig strand. */
	public int entryDistance(final int index, final boolean reverse){
		checkIndex(index);
		return reverse ? lengths.get(index)-positions.get(index)-32 : positions.get(index);
	}

	/** Rejects an anchor located at the wrong end for an outward-directed connection. */
	public boolean canExit(final int index, final boolean reverse){return exitDistance(index, reverse)<=tipBases-32;}

	/** Rejects an anchor located at the wrong end for an inward-directed connection. */
	public boolean canEnter(final int index, final boolean reverse){return entryDistance(index, reverse)<=tipBases-32;}

	/** Implied signed gap: positive means missing sequence, negative means overlap. */
	public long gap(final int left, final boolean leftReverse, final int right,
			final boolean rightReverse, final int span){
		if(span<=64){throw new IllegalArgumentException("Hollow anchor geometry requires span>64.");}
		assert(isUnique(left) && isUnique(right)) : "Connection geometry must use assembly-unique anchors.";
		return span-64L-exitDistance(left, leftReverse)-entryDistance(right, rightReverse);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Internal methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds first tip positions for canonical identities in a bounded start interval. */
	private void addCandidates(final byte[] bases, final int owner, final int first, final int last){
		assert(first>=0 && first<=last && last<=bases.length-32) : "Tip start bounds must contain complete32bp words.";
		long forward=0, reverse=0;
		int valid=0;
		for(int pos=first; pos<=last+31; pos++){
			final int base=AminoAcid.baseToNumber[bases[pos]];
			if(base<0){forward=reverse=0; valid=0; continue;}
			forward=(forward<<2)|base;
			reverse=(reverse>>>2)|((3L-base)<<62);
			valid++;
			if(valid<32 || (!retainPalindromes && forward==reverse)){continue;}
			final long key=Long.compareUnsigned(forward, reverse)>=0 ? forward : reverse;
			if(map.get(key)>=1){continue;}
			// The stored value is an index, not a counter. get/put avoid increment's resize-return contract.
			map.put(key, words.size+1);
			words.add(forward);
			owners.add(owner);
			positions.add(pos-31);
			lengths.add(bases.length);
		}
	}

	/** Counts physical canonical occurrences across all positions, capping at nonunique. */
	private void countOccurrences(final byte[] bases){
		assert(occurrences.length==words.size) : "All candidates must be collected before global uniqueness counting.";
		long forward=0, reverse=0;
		int valid=0;
		for(int pos=0; pos<bases.length; pos++){
			final int base=AminoAcid.baseToNumber[bases[pos]];
			if(base<0){forward=reverse=0; valid=0; continue;}
			forward=(forward<<2)|base;
			reverse=(reverse>>>2)|((3L-base)<<62);
			valid++;
			if(valid<32){continue;}
			final long key=Long.compareUnsigned(forward, reverse)>=0 ? forward : reverse;
			final int index=map.get(key)-1;
			if(index>=0 && occurrences[index]<2){occurrences[index]++;}
		}
	}

	/** Canonical32bp anchor identity is independent of the direction used for a connection. */
	private static long canonical(final long word){
		final long reverse=HollowAnchorTable.reverseComplement(word);
		return Long.compareUnsigned(word, reverse)>=0 ? word : reverse;
	}

	/** Asserts the caller supplied an index from this snapshot when assertions are enabled. */
	private void checkIndex(final int index){
		assert(index>=0 && index<words.size) : "Anchor index must come from this snapshot's lookup or iteration: "+index;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Terminal interval size used for candidate discovery and endpoint gates. */
	public final int tipBases;
	/** True only for evidence snapshots that keep strand-ambiguous candidates. */
	private final boolean retainPalindromes;
	/** Cached number of candidates that pass global uniqueness and orientation filters. */
	private final int uniqueCount;
	/** Canonical 32 bp word to one-based candidate index. */
	private final LongIntMap map=new LongIntMap();
	/** Literal forward words from the first retained tip occurrence of each canonical key. */
	private final LongList words=new LongList();
	/** Owning list indices, zero-based starts and snapshot contig lengths for each word. */
	private final IntList owners=new IntList(), positions=new IntList(), lengths=new IntList();
	/** Physical occurrence counts for candidate canonical words, capped at two. */
	private final byte[] occurrences;
}
