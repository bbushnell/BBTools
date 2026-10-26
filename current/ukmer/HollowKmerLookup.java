package ukmer;

import dna.AminoAcid;
import shared.Tools;

/** Queries an already loaded Hollow table from two oriented 32-base anchors.
 * No full-span Kmer is allocated, reconstructed or changed by a query. The table
 * must remain immutable for the lifetime of this view. Observed anchor spacing
 * is chosen by the caller; fingerprints always use the loaded table's K salt.
 * @author Brian Bushnell, Fischl */
public final class HollowKmerLookup {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Captures a loaded Hollow table and its strand convention without copying counts. */
	public HollowKmerLookup(final KmerTableSetU source){
		this(requireHollow(source), source.kbig, source.rcomp());
	}

	/** Creates a view of explicitly supplied Hollow partitions, also useful for fixtures.
	 * The caller guarantees that counts were built with this K and strand convention. */
	public HollowKmerLookup(final AbstractKmerTableU[] source, final int k, final boolean rcomp_){
		if(k<=64 || source==null || source.length==0){
			throw new IllegalArgumentException("Hollow lookup requires K>64 and loaded table partitions.");
		}
		tables=new HashArrayH1D[source.length];
		for(int i=0; i<source.length; i++){
			if(!(source[i] instanceof HashArrayH1D) || ((HashArrayH1D)source[i]).fingerprintOnly()){
				throw new IllegalArgumentException("Hollow lookup requires resizable two-hash partitions; invalid way="+i);
			}
			tables[i]=(HashArrayH1D)source[i];
		}
		rcomp=rcomp_;
		placementSalt=Tools.hash64shift(k);
		fingerprintSalt=Tools.splitMix64(k);
		assert(tables.length==source.length) : "Kmer.mod routes by the original partition count; a view must preserve every way.";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the count for forward-oriented anchors, normalizing an absent key to zero.
	 * The upper bits contain the first bases, exactly as in a full packed Kmer word.
	 * All 64 bits of each anchor are valid data, including negative long values. */
	public int count(long left, long right){
		assert(tables.length>0) : "The placement hash needs the same nonempty partition set used while loading counts.";
		if(rcomp){
			final long reverseLeft=AminoAcid.reverseComplementBinaryFast(right, 32);
			final int first=Long.compareUnsigned(left, reverseLeft);
			if(first<=0){
				final long reverseRight=AminoAcid.reverseComplementBinaryFast(left, 32);
				if(first<0 || Long.compareUnsigned(right, reverseRight)<0){
					left=reverseLeft;
					right=reverseRight;
				}
			}
		}
		// Keep these formulas identical to Kmer.xor0()/xor2() in Hollow mode.
		final long placement=(Long.rotateLeft(Tools.hash64shift(left), 25)
				^Tools.hash64shift(right)^placementSalt)&Long.MAX_VALUE;
		final long fingerprint=Tools.splitMix64(Long.rotateRight(Tools.splitMix64(left), 23)
				^right^fingerprintSalt)&Long.MAX_VALUE;
		return Math.max(0, tables[(int)(placement%tables.length)].getValue(placement, fingerprint));
	}

	/** Counts all right-shifted anchor pairs without trusting the interior ingress base.
	 * Each old32-base anchor loses its first base and receives one A/C/G/T base.
	 * counts[(leftBase<<2)|rightBase] uses forward sequence orientation, even when
	 * the stored key is reverse-complement canonical. The caller supplies a reused
	 * 16-element buffer; this method allocates nothing and does not change counts.
	 * These are observations only, not permission to correct a base or extend a path. */
	public void fillRightCounts16(final long left, final long right, final int[] counts){
		if(counts==null || counts.length!=16){
			throw new IllegalArgumentException("Joint anchor successors require exactly16 count slots.");
		}
		final long shiftedLeft=left<<2, shiftedRight=right<<2;
		for(int leftBase=0; leftBase<4; leftBase++){
			for(int rightBase=0; rightBase<4; rightBase++){
				counts[(leftBase<<2)|rightBase]=count(shiftedLeft|leftBase, shiftedRight|rightBase);
			}
		}
		assert(counts[0]>=0 && counts[15]>=0) :
				"Endpoint count normalizes absent keys to zero; successor selection must never receive table sentinels.";
	}

	/** Rejects an unloaded or solid source before any endpoint can query incompatible counts. */
	private static AbstractKmerTableU[] requireHollow(final KmerTableSetU source){
		if(source==null || !source.hollow){
			throw new IllegalArgumentException("Endpoint-only lookup requires a Hollow table.");
		}
		return source.tables();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final HashArrayH1D[] tables;
	private final boolean rcomp;
	private final long placementSalt, fingerprintSalt;
}
