package assemble;

import stream.Read;

/** Transactional application of one arbitrary single-base edit.
 * This is deliberately separate from homopolymer classification and batch
 * acceptance: positions are explicit original-read coordinates.
 * @author Diona
 */
final class LocalSingleBaseEdit {

	/** The operation-specific overloads keep deletion from needing a base sentinel. */
	enum Operation { SUBSTITUTION, DELETION, INSERTION }

	/** Apply a deletion at the original base position. */
	void apply(final Read read, final Operation operation, final int position){
		if(operation!=Operation.DELETION){
			throw new IllegalArgumentException("This overload is only for DELETION.");
		}
		applyInternal(read, operation, position, (byte)0);
	}

	/** Apply a substitution at position, or insert {@code base} before position. */
	void apply(final Read read, final Operation operation, final int position, final byte base){
		if(operation==Operation.DELETION){
			throw new IllegalArgumentException("DELETION does not accept a replacement base.");
		}
		if(operation!=Operation.SUBSTITUTION && operation!=Operation.INSERTION){
			throw new IllegalArgumentException("Unknown local edit operation: "+operation);
		}
		if(!defined(base)){
			throw new IllegalArgumentException("Local edit base must be uppercase A/C/G/T: "+(char)base);
		}
		applyInternal(read, operation, position, base);
	}

	private void applyInternal(final Read read, final Operation operation, final int position, final byte base){
		// Reject attached coordinates/provenance before examining or replacing arrays.
		HomopolymerIndelEdit.requireEditable(read);
		if(operation==null){throw new IllegalArgumentException("Local edit operation is null.");}
		if(read.bases==null){throw new IllegalArgumentException("Local edit requires non-null bases.");}
		final byte[] oldBases=read.bases, oldQuality=read.quality;
		if(oldQuality!=null && oldQuality.length!=oldBases.length){
			throw new IllegalArgumentException("Local edit requires equal sequence/quality lengths, or null qualities.");
		}
		final int length=oldBases.length;
		final boolean insertion=operation==Operation.INSERTION;
		final boolean substitution=operation==Operation.SUBSTITUTION;
		if(insertion ? position<0 || position>length : position<0 || position>=length){
			throw new IllegalArgumentException("Local edit position is outside its operation domain: "+operation+" "+position);
		}
		if(insertion && length==Integer.MAX_VALUE){
			throw new IllegalArgumentException("Insertion would overflow the maximum array length.");
		}

		final long newLength=(long)length+(insertion ? 1 : substitution ? 0 : -1);
		if(newLength<0 || newLength>Integer.MAX_VALUE){
			throw new IllegalArgumentException("Local edit result length is out of range: "+newLength);
		}
		final byte[] newBases=new byte[(int)newLength];
		final byte[] newQuality=oldQuality==null ? null : new byte[newBases.length];
		if(substitution){
			System.arraycopy(oldBases,0,newBases,0,length);
			newBases[position]=base;
			if(newQuality!=null){System.arraycopy(oldQuality,0,newQuality,0,length);newQuality[position]=0;}
		}else if(insertion){
			System.arraycopy(oldBases,0,newBases,0,position);
			newBases[position]=base;
			System.arraycopy(oldBases,position,newBases,position+1,length-position);
			if(newQuality!=null){
				System.arraycopy(oldQuality,0,newQuality,0,position);
				System.arraycopy(oldQuality,position,newQuality,position+1,length-position);
			}
		}else{
			System.arraycopy(oldBases,0,newBases,0,position);
			System.arraycopy(oldBases,position+1,newBases,position,length-position-1);
			if(newQuality!=null){
				System.arraycopy(oldQuality,0,newQuality,0,position);
				System.arraycopy(oldQuality,position+1,newQuality,position,length-position-1);
			}
		}
		assert((long)newBases.length==newLength && (newQuality==null || newQuality.length==newBases.length)) :
			"Local edit reconstruction length mismatch: expected="+newLength+
			" bases="+newBases.length+" quality="+(newQuality==null ? -1 : newQuality.length);
		// No operation after this point can fail; install only complete replacements.
		read.bases=newBases;
		read.quality=newQuality;
	}

	private static boolean defined(final byte base){return base=='A' || base=='C' || base=='G' || base=='T';}
}
