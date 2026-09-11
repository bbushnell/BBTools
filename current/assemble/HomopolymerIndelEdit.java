package assemble;

import stream.Read;

/**
 * Mechanical application of ONE already approved run-length edit. This helper
 * neither classifies errors nor changes original arrays. Its guarded Read
 * overload rejects alignment metadata before replacing any Read arrays. Native
 * callers require experimental opt-in; singleton +1 requires separate opt-in.
 * One reusable result holder belongs to one worker.
 * @author Fischl
 */
final class HomopolymerIndelEdit {
	HomopolymerIndelEdit(){this(false);}
	HomopolymerIndelEdit(final boolean singletons_){singletons=singletons_;}

	/**
	 * Apply to a worker-owned, unpaired, alignment-free Read. Nonzero edits reject
	 * metadata rather than silently invalidating somebody else's coordinate system.
	 * Both new arrays are built before either is installed; this is exception-safe,
	 * not synchronization for a concurrently shared Read. Rescan after each edit;
	 * intervals from the old length must not be reused. No-op leaves all state intact.
	 */
	void apply(final Read read, final int start, final int end, final int delta){
		editedBases=editedQualities=null;
		position=-1;
		if(read==null){throw new IllegalArgumentException("Cannot edit a null Read.");}
		if(delta!=0){requireEditable(read);}
		apply(read.bases, read.quality, start, end, delta);
		if(delta==0){return;}
		assert(editedBases.length==read.bases.length+delta &&
			(editedQualities==null || editedQualities.length==editedBases.length)) :
			"Read array installation requires a complete consistent splice result.";
		read.bases=editedBases;
		read.quality=editedQualities;
	}

	/** Shared metadata guard; do not use clearMapping as a substitute for this contract. */
	static void requireEditable(final Read read){
		if(read==null){throw new IllegalArgumentException("Cannot edit a null Read.");}
		if(read.mate!=null || read.paired() || read.pairnum()!=0 || read.mapped() || read.perfect() || read.shortmatch() ||
				read.secondary() || read.supplementary() || read.rescued() || read.ambiguous() || read.insertvalid() || read.aminoacid() ||
				read.samline!=null || read.match!=null || read.gaps!=null || read.sites!=null || read.originalSite!=null ||
				read.obj!=null || read.chrom!=-1 || read.start!=-1 || read.stop!=-1 || read.mapScore!=0 || read.errors!=0 ||
				read.insert()>=0){
			throw new IllegalArgumentException("Homopolymer editing requires an unpaired alignment-free Read with no attached coordinate metadata: "+read.id);
		}
	}

	/**
	 * Coordinates refer to the supplied oriented sequence, not a reference.
	 * A/C runs use their left edge; G/T runs use their right edge, so placement
	 * and retained quality bytes commute with reverse complementation, including
	 * even-length runs and quality ties. Placement is bookkeeping, not a claim
	 * about which physical homopolymer base was erroneous. An inserted base has
	 * quality zero (unknown); existing qualities are never inflated. Null quality
	 * stays null. A zero delta returns the original array identities unchanged.
	 */
	void apply(final byte[] bases, final byte[] qualities, final int start,
			final int end, final int delta){
		// Clear the reusable result even when an invalid request throws.
		editedBases=editedQualities=null;
		position=-1;
		if(bases==null || (qualities!=null && qualities.length!=bases.length)){
			throw new IllegalArgumentException("Edit requires equal sequence/quality lengths, or null qualities.");
		}
		if(delta<-1 || delta>1 || start<0 || end<start || end>bases.length){
			throw new IllegalArgumentException("Invalid one-base edit interval/delta: "+start+".."+end+" delta="+delta);
		}
		if(delta==0){editedBases=bases; editedQualities=qualities; return;}
		if(end-start<2 && !(singletons && end-start==1 && delta==1)){
			throw new IllegalArgumentException("Approved edit requires run>=2, or explicitly enabled singleton +1.");
		}
		final byte b=bases[start];
		if(b!='A' && b!='C' && b!='G' && b!='T'){
			throw new IllegalArgumentException("Approved run must contain uppercase A/C/G/T.");
		}
		for(int i=start+1; i<end; i++){
			if(bases[i]!=b){throw new IllegalArgumentException("Edit interval is not a homopolymer.");}
		}
		if((start>0 && bases[start-1]==b) || (end<bases.length && bases[end]==b)){
			throw new IllegalArgumentException("Edit interval must describe the entire maximal run.");
		}
		if(delta>0 && bases.length==Integer.MAX_VALUE){
			throw new IllegalArgumentException("Edited sequence length would overflow.");
		}
		final boolean left=b=='A' || b=='C';
		final int pos=left ? start : delta<0 ? end-1 : end;
		final int sourceTail=pos+(delta<0 ? 1 : 0);
		final int targetTail=pos+(delta>0 ? 1 : 0);
		final byte[] out=new byte[bases.length+delta];
		final byte[] outQ=qualities==null ? null : new byte[out.length];
		System.arraycopy(bases, 0, out, 0, pos);
		System.arraycopy(bases, sourceTail, out, targetTail, bases.length-sourceTail);
		if(delta>0){out[pos]=b;}
		if(outQ!=null){
			System.arraycopy(qualities, 0, outQ, 0, pos);
			System.arraycopy(qualities, sourceTail, outQ, targetTail, qualities.length-sourceTail);
			// The new slot stays zero: kmer support is not a calibrated base quality.
		}
		assert(targetTail+bases.length-sourceTail==out.length && (outQ==null || outQ.length==out.length)) :
			"Splice suffix and quality coordinates must cover exactly the length-changing output.";
		editedBases=out;
		editedQualities=outQ;
		position=pos;
	}

	byte[] editedBases, editedQualities;
	/** Insertion boundary or deleted-base index in the original oriented read. */
	int position=-1;
	private final boolean singletons;
}
