package assemble;

import dna.AminoAcid;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Protects an original seed base when every overlapping solid window is supported.
 * This worker-local view borrows an immutable original-read-only table. It does
 * not require solid support elsewhere on an extension, edit bases, or own counts.
 * @author Brian Bushnell, Fischl */
public final class HollowSeedGuard {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires loaded canonical exact Solid64/96 counts, excluding hashed tables. */
	public HollowSeedGuard(final KmerTableSetU source, final int minimum){
		if(source==null || source.hollow || !source.rcomp() || (source.kbig!=64 && source.kbig!=96)
				|| source.tables()==null || source.tables().length==0 || minimum<1){
			throw new IllegalArgumentException("Seed protection needs loaded exact Solid64/96 counts and positive minimum depth.");
		}
		for(AbstractKmerTableU partition : source.tables()){
			if(!(partition instanceof HashArrayU1D)){
				throw new IllegalArgumentException("Original seed protection requires explicit exact count partitions, not fingerprints.");
			}
		}
		tables=source;
		k=source.kbig;
		minDepth=minimum;
		word=new Kmer(k, false);
		assert(!word.hollow && word.kbig==k) : "Solid evidence must be queried with the same full-key geometry; Kmer's boolean flag means Hollow.";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests the unchanged oriented seed before a proposed substitution is applied.
	 * originalLength excludes appended bases. A short seed, an appended position,
	 * an undefined overlapping window or any window below minDepth abstains from
	 * protection; lack of support is not permission to correct. The caller must
	 * still apply every ordinary recovery and rollback condition. */
	public boolean protects(final byte[] bases, final int originalLength, final int position){
		return originalMaximum(bases, originalLength, position)>=minDepth;
	}

	/** Optionally releases a supported original allele only when every proposed
	 * solid window has at least four times the strongest original-window count.
	 * The numeric replacement is inspected without editing bases or table counts;
	 * returning false merely releases this veto, not the ordinary rescue gates. */
	public boolean protects(final byte[] bases, final int originalLength, final int position, final int replacement){
		if(replacement<0 || replacement>3){
			throw new IllegalArgumentException("Relative seed protection requires an A/C/G/T numeric replacement.");
		}
		final int originalMax=originalMaximum(bases, originalLength, position);
		if(originalMax<minDepth){return false;}
		return proposalBelow(bases, originalLength, position, replacement, 4L*originalMax);
	}

	/** Keeps full original support absolute, but tests partial support relatively.
	 * Any supported original window protects unless every proposed window reaches
	 * four times the original maximum. Weak/absent originals do not require solid
	 * confirmation here; the separate proposed-window veto remains authoritative. */
	public boolean protectsPartial(final byte[] bases, final int originalLength, final int position, final int replacement){
		if(bases==null || originalLength<0 || originalLength>bases.length
				|| position<0 || position>=bases.length || replacement<0 || replacement>3){
			throw new IllegalArgumentException("Partial seed protection needs valid seed bounds and an A/C/G/T replacement.");
		}
		if(originalLength<k || position>=originalLength){return false;}
		final int first=Math.max(0, position-k+1), last=Math.min(position, originalLength-k);
		assert(first<=last && last+k<=originalLength) :
				"Partial support must inspect all edit-overlapping windows while excluding appended sequence.";
		int minimum=Integer.MAX_VALUE, maximum=0;
		word.clear();
		for(int i=first; i<last+k; i++){
			final int base=AminoAcid.baseToNumber[bases[i]];
			if(base<0){return false;}
			word.addRightNumeric(base);
			if(i>=first+k-1){
				final int count=Math.max(0, tables.getCount(word));
				minimum=Math.min(minimum, count);
				maximum=Math.max(maximum, count);
			}
		}
		if(minimum>=minDepth){return true;}
		if(maximum<minDepth){return false;}
		return proposalBelow(bases, originalLength, position, replacement, 4L*maximum);
	}

	/** Queries a defined candidate without changing the seed or immutable evidence.
	 * Both callers first validated every original base in this seed-contained range. */
	private boolean proposalBelow(final byte[] bases, final int originalLength, final int position,
			final int replacement, final long required){
		final int first=Math.max(0, position-k+1), last=Math.min(position, originalLength-k);
		assert(required>0 && first<=last && last+k<=originalLength) :
				"A protected original allele supplies defined seed-contained windows for the identical proposed-coordinate scan.";
		word.clear();
		for(int i=first; i<last+k; i++){
			final int base=i==position ? replacement : AminoAcid.baseToNumber[bases[i]];
			assert(base>=0) : "Original-support scans must validate the candidate range before comparing depths.";
			word.addRightNumeric(base);
			if(i>=first+k-1 && tables.getCount(word)<required){return true;}
		}
		return false;
	}

	/** Vetoes a seed substitution if any proposed seed-contained window lacks support.
	 * This optional check is independent of original support and requires depth one,
	 * not minDepth. Appended positions and seeds shorter than k are outside its scope.
	 * The replacement is queried without changing bases or the shared count table. */
	public boolean rejectsProposal(final byte[] bases, final int originalLength, final int position, final int replacement){
		if(bases==null || originalLength<0 || originalLength>bases.length
				|| position<0 || position>=bases.length || replacement<0 || replacement>3){
			throw new IllegalArgumentException("Proposed seed support needs valid seed bounds and an A/C/G/T replacement.");
		}
		if(originalLength<k || position>=originalLength){return false;}
		final int first=Math.max(0, position-k+1), last=Math.min(position, originalLength-k);
		assert(first<=last && last+k<=originalLength) :
				"Proposed support must cover the edit but exclude every appended base, matching the original-seed guard geometry.";
		word.clear();
		for(int i=first; i<last+k; i++){
			final int base=i==position ? replacement : AminoAcid.baseToNumber[bases[i]];
			if(base<0){return true;}
			word.addRightNumeric(base);
			if(i>=first+k-1 && tables.getCount(word)<1){return true;}
		}
		return false;
	}

	/** Returns the maximum only if ALL original windows are supported; otherwise -1.
	 * Unsupported or undefined originals exit before any proposed-word scan. */
	private int originalMaximum(final byte[] bases, final int originalLength, final int position){
		if(bases==null || originalLength<0 || originalLength>bases.length
				|| position<0 || position>=bases.length){
			throw new IllegalArgumentException("Seed protection needs valid original-length and proposed-position coordinates.");
		}
		if(originalLength<k || position>=originalLength){return -1;}
		final int first=Math.max(0, position-k+1), last=Math.min(position, originalLength-k);
		assert(first<=last && last+k<=originalLength) : "Only windows containing the proposed edit and entirely inside the original seed may veto it.";
		word.clear();
		int maximum=0;
		for(int i=first; i<last+k; i++){
			if(AminoAcid.baseToNumber[bases[i]]<0){return -1;}
			word.addRight(bases[i]);
			if(i>=first+k-1){
				final int count=tables.getCount(word);
				if(count<minDepth){return -1;}
				maximum=Math.max(maximum, count);
			}
		}
		return maximum;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final KmerTableSetU tables;
	private final Kmer word;
	private final int k, minDepth;
}
