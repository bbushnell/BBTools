package assemble;

import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/** Experimental integrated single-edit caller: original depths -> local probes ->
 * full winning-context verification -> transactional edit. No native CLI wiring.
 * One edit per call; repeated calls MUST recompute the profile (this method does).
 * Canonical whole-read orientation is a correctness-first reference implementation,
 * not a claim of optimal allocation or runtime. Worker-local; immutable count table.
 * @author Fischl */
final class LocalEditCorrector {
	LocalEditCorrector(final int k_,final HomopolymerIndelProposal.CountLookup lookup_){this(k_,lookup_,1);}
	LocalEditCorrector(final int k_,final HomopolymerIndelProposal.CountLookup lookup_,final int windows){
		this(k_,lookup_,windows,false);
	}
	LocalEditCorrector(final int k_,final HomopolymerIndelProposal.CountLookup lookup_,final int windows,final boolean checkIndelCompetition_){
		if(lookup_==null || k_<5){throw new IllegalArgumentException("Corrector requires K>=5 and immutable counts.");}
		k=k_;lookup=lookup_;key=new Kmer(k);probe=new LocalEditKmerProbe(k,lookup,windows);locator=new LocalEditTroughLocator(k,3);
		checkIndelCompetition=checkIndelCompetition_;
		if(key.kbig!=k){throw new IllegalArgumentException("Correction K must equal table K.");}
	}
	/** Return 0 or 1 applied edits. All original arrays remain unchanged. */
	int correctOne(final Read read){
		return correctOne(read,false);
	}
	/** Experimental opt-in: commit only the first edit of a verified pair witness. */
	int correctOne(final Read read,final boolean nearbyPairs){
		return correctOne(read,nearbyPairs,-1);
	}
	/** On the first read pass, estimate separated errors from depth geometry,
	 * without probing mutations. This is a rejection heuristic, not verified repairs.
	 * Negative initialLimit preserves the ordinary first-edit discovery path. */
	int correctOne(final Read read,final boolean nearbyPairs,final int initialLimit){
		HomopolymerIndelEdit.requireEditable(read);
		if(read.bases==null || (read.quality!=null && read.quality.length!=read.bases.length)){
			throw new IllegalArgumentException("Correction requires bases with matching or null qualities.");
		}
		profileQueries=probeQueries=verificationQueries=0;lastOperation=null;lastPosition=-1;
		pairQueries=0;usedPairLookahead=false;initialEstimatedEdits=0;verifiedLoci=0;
		if(pairExcludedStarts!=null){pairExcludedStarts.clear();}
		lowRegions=skippedEdge=skippedWide=skippedUndefined=acceptedTroughs=0;
		ambiguousLoci=noSupportedCandidateLoci=allCandidatesRejectedLoci=0;
		supportedCandidates=verificationRejectedCandidates=0;callStatus=null;
		if(read.length()<k+2){callStatus=CallStatus.SHORT_READ;return 0;}
		for(byte b:read.bases){if(b!='A' && b!='C' && b!='G' && b!='T' && b!='N'){callStatus=CallStatus.UNSUPPORTED_BASE;return 0;}}
		final int orientation=orientation(read.bases);
		// A self-RC sequence cannot have a strand-invariant single nonsymmetric edit.
		if(orientation==0){callStatus=CallStatus.SELF_RC;return 0;}
		final boolean reverse=orientation>0;
		final byte[] bases=reverse ? reverseComplement(read.bases) : read.bases;
		fillDepths(bases);locator.reset(bases,counts);
		if(initialLimit>=0){
			initialEstimatedEdits=estimateInitialEdits(initialLimit);
			// Reuse the same depths for ordinary discovery; preflight adds no lookups.
			locator.reset(bases,counts);
			if(initialEstimatedEdits>initialLimit){finishScan(CallStatus.INITIAL_LIMIT);return 0;}
		}
		probe.beginRead(bases);
		while(locator.next()){
			acceptedTroughs++;final long supportedBefore=supportedCandidates;
			chosenOperation=null;chosenPosition=-1;ambiguous=false;
			final int a=locator.depthStart,end=locator.depthEnd,w=a+(end-a-1)/2;
			int originalMax=0;for(int j=a;j<end;j++){originalMax=Math.max(originalMax,counts.get(j));}
			assert(originalMax<3) : "LocalEditTroughLocator(k,3) returns only depth<3 windows; contrast multiplication is bounded.";
			final int threshold=Math.max(4,originalMax*4+1);
			for(int p=locator.baseStart;p<locator.baseEnd;p++){
				probe.substitutions(bases,w,p);probeQueries+=probe.queries;
				for(int b=0;b<4;b++){if(probe.substitutionDepth[b]>=threshold){consider(bases,LocalSingleBaseEdit.Operation.SUBSTITUTION,p,ALPHABET[b],a-1,end+k,threshold);}}
			}
			if(ambiguous){recordAmbiguous(a,nearbyPairs);continue;}
			//Approximate counts can support a false S beside a verified indel.
			//Opt-in callers must not hide that ambiguity through S-first ordering.
			if(chosenOperation==null || checkIndelCompetition){
				for(int p=locator.baseStart;p<locator.baseEnd;p++){
					probe.indelsAt(bases,w,p);probeQueries+=probe.queries;
					if(probe.deletionDepth>=threshold){consider(bases,LocalSingleBaseEdit.Operation.DELETION,p,(byte)0,a-1,end+k,threshold);}
					if(p>=locator.gapStart && p<locator.gapEnd){for(int b=0;b<4;b++){
						if(probe.insertionDepth[b]>=threshold){consider(bases,LocalSingleBaseEdit.Operation.INSERTION,p,ALPHABET[b],a-1,end+k,threshold);}
					}}
				}
			}
			if(ambiguous){recordAmbiguous(a,nearbyPairs);continue;}
			if(chosenOperation==null){
				if(supportedCandidates==supportedBefore){noSupportedCandidateLoci++;}
				else{allCandidatesRejectedLoci++;}
				continue;
			}
			verifiedLoci++;
			applyChosen(read,bases,reverse);finishScan(CallStatus.APPLIED);return 1;
		}
		if(nearbyPairs){
			if(pairLookahead==null){pairLookahead=new LocalEditPairLookahead(k,lookup);}
			final boolean found=pairLookahead.propose(bases,counts,pairExcludedStarts);pairQueries=pairLookahead.queries;
			if(found){
				chosenOperation=pairLookahead.firstOperation;chosenPosition=pairLookahead.firstPosition;chosenBase=pairLookahead.firstBase;
				applyChosen(read,bases,reverse);usedPairLookahead=true;finishScan(CallStatus.APPLIED_PAIR);return 1;
			}
		}
		finishScan(CallStatus.EXHAUSTED);
		return 0;
	}
	/** Count only K-1/K-wide troughs with strong flanks and separated contexts.
	 * An isolated missing base disrupts K-1 windows; a substitution/extra base K.
	 * Short collision holes, edge/N regions and merged troughs are not counted.
	 * Coverage/variants can mimic errors: this estimates burden, not biological truth.
	 * Underestimation is handled by the engine's max+1 whole-read rollback. */
	private long estimateInitialEdits(final int limit){
		assert(limit>=0) : "A negative initial limit disables preflight in correctOne.";
		long estimated=0,lastEnd=-1;
		while(locator.next()){
			final int a=locator.depthStart,end=locator.depthEnd;
			if(end-a<k-1){continue;}
			int originalMax=0;for(int j=a;j<end;j++){originalMax=Math.max(originalMax,counts.get(j));}
			assert(originalMax<3) : "LocalEditTroughLocator(k,3) bounds low-window depth at 2, so 4*depth+1 cannot overflow.";
			final int threshold=Math.max(4,originalMax*4+1);
			if(counts.get(a-1)<threshold || counts.get(end)<threshold){continue;}
			// Disjoint full flank-window spans avoid treating neighboring troughs
			// as independent evidence. No proposed or normalized edit is needed.
			if(a-1L<lastEnd){continue;}
			lastEnd=(long)end+k;
			if(++estimated>limit){break;}
		}
		return estimated;
	}
	/** Do not let lookahead reinterpret a locus with competing verified single edits. */
	private void recordAmbiguous(final int start,final boolean nearbyPairs){
		assert(start>=0) : "Excluded pair loci use original canonical kmer-start coordinates from the ordinary locator.";
		ambiguousLoci++;
		if(nearbyPairs){
			if(pairExcludedStarts==null){pairExcludedStarts=new IntList();}
			assert(pairExcludedStarts.size==0 || pairExcludedStarts.get(pairExcludedStarts.size-1)<start) :
				"The ordinary locator visits disjoint troughs in order; pair exclusion uses an ordered merge against the same profile.";
			pairExcludedStarts.add(start);
		}
	}
	private void applyChosen(final Read read,final byte[] bases,final boolean reverse){
		assert(chosenOperation!=null && chosenPosition>=0) : "Only a verified single-edit or pair-first proposal may reach transactional application.";
		final byte[] quality=reverse && read.quality!=null ? reversed(read.quality) : read.quality;
		final Read working=new Read(bases,quality,read.id,read.numericID,false);
		if(chosenOperation==LocalSingleBaseEdit.Operation.DELETION){editor.apply(working,chosenOperation,chosenPosition);}
		else{editor.apply(working,chosenOperation,chosenPosition,chosenBase);}
		final byte[] output=reverse ? reverseComplement(working.bases) : working.bases;
		final byte[] outputQ=reverse && working.quality!=null ? reversed(working.quality) : working.quality;
		lastOperation=chosenOperation;
		lastPosition=reverse ? read.length()-(chosenOperation==LocalSingleBaseEdit.Operation.INSERTION ? 0 : 1)-chosenPosition : chosenPosition;
		assert(outputQ==null || outputQ.length==output.length) : "Canonical correction must restore one quality per output base before installation.";
		read.bases=output;read.quality=outputQ;
	}
	/** Snapshot only the portion scanned in this call; APPLIED may stop early. */
	private void finishScan(final CallStatus status){
		callStatus=status;lowRegions=locator.lowRegions;
		skippedEdge=locator.skippedEdge;skippedWide=locator.skippedWide;skippedUndefined=locator.skippedUndefined;
		assert(lowRegions==skippedEdge+skippedWide+skippedUndefined+acceptedTroughs) :
			"Each encountered low region is skipped once or returned by the locator; counters must describe this call only.";
		assert(acceptedTroughs==ambiguousLoci+noSupportedCandidateLoci+allCandidatesRejectedLoci+verifiedLoci) :
				"Every correction-profile trough is ambiguous, unsupported, rejected or verified; depth-only preflight resets locator counters before discovery.";
	}
	private void fillDepths(final byte[] bases){
		counts.clear();key.clearFast();
		for(int i=0;i<bases.length;i++){
			if(bases[i]=='N'){key.clearFast();}else{key.addRight(bases[i]);}
			if(i>=k-1){counts.add(key.len()>=k ? count(false) : 0);}
		}
		assert(counts.size==bases.length-k+1) : "Locator requires one original count per kmer-start position.";
	}
	private void consider(final byte[] bases,final LocalSingleBaseEdit.Operation op,final int p,final byte b,final int from,final int to,final int threshold){
		supportedCandidates++;
		if(!verify(bases,op,p,b,from,to,threshold)){verificationRejectedCandidates++;return;}
		final int normalized=normalize(bases,op,p,b);
		if(chosenOperation==null){chosenOperation=op;chosenPosition=normalized;chosenBase=b;}
		else if(chosenOperation!=op || chosenPosition!=normalized || chosenBase!=b){ambiguous=true;}
	}
	/** Verify every affected kmer as well as the trough and its supported flank windows.
	 * Only candidates whose single-word probe passed reach this streaming check. */
	private boolean verify(final byte[] bases,final LocalSingleBaseEdit.Operation op,final int p,final byte b,final int contextFrom,final int contextTo,final int threshold){
		final int delta=op==LocalSingleBaseEdit.Operation.INSERTION ? 1 : op==LocalSingleBaseEdit.Operation.DELETION ? -1 : 0;
		//A narrow trough does not span every window changed by its candidate edit.
		//S changes starts [p-K+1,p]; D creates [p-K+1,p-1]; I creates [p-K+1,p].
		//Convert the last affected edited window back to an original exclusive bound.
		final int from=Math.max(0,Math.min(contextFrom,p-k+1));
		final int to=Math.min(bases.length,Math.max(contextTo,p+k-(delta>0 ? 1 : 0)));
		assert(from<=p && p<to && to<=bases.length) : "Selected edit must lie inside its complete original verification context.";
		key.clearFast();int measured=0;
		for(int i=0;i<to-from+delta;i++){
			final int relative=p-from;
			final byte x=op==LocalSingleBaseEdit.Operation.SUBSTITUTION ? (i==relative ? b : bases[from+i]) :
				delta>0 ? (i==relative ? b : bases[from+i-(i>relative ? 1 : 0)]) : bases[from+i+(i>=relative ? 1 : 0)];
			if(x!='A' && x!='C' && x!='G' && x!='T'){return false;}key.addRight(x);
			if(key.len()>=k){measured++;if(count(true)<threshold){return false;}}
		}
		assert(measured>0) : "A winning candidate must have at least one complete verification kmer.";return true;
	}
	/** Equivalent 1bp indel placements collapse without an HP-specific error model. */
	private static int normalize(final byte[] bases,final LocalSingleBaseEdit.Operation op,int p,final byte inserted){
		if(op==LocalSingleBaseEdit.Operation.SUBSTITUTION){return p;}
		final byte b=op==LocalSingleBaseEdit.Operation.DELETION ? bases[p] : inserted;
		if(b=='A' || b=='C'){while(p>0 && bases[p-1]==b){p--;}}
		else if(op==LocalSingleBaseEdit.Operation.DELETION){while(p+1<bases.length && bases[p+1]==b){p++;}}
		else{while(p<bases.length && bases[p]==b){p++;}}
		return p;
	}
	private int count(final boolean verification){
		final int value=lookup.count(key);if(verification){verificationQueries++;}else{profileQueries++;}
		if(value<-1){throw new IllegalStateException("Invalid count depth: "+value);}return Math.max(0,value);
	}
	private static int orientation(final byte[] b){
		for(int i=0;i<b.length;i++){final int delta=b[i]-complement(b[b.length-1-i]);if(delta!=0){return delta;}}return 0;
	}
	private static byte complement(final byte b){return b=='A' ? (byte)'T' : b=='C' ? (byte)'G' : b=='G' ? (byte)'C' : b=='T' ? (byte)'A' : (byte)'N';}
	private static byte[] reverseComplement(final byte[] b){final byte[] out=new byte[b.length];for(int i=0;i<b.length;i++){out[i]=complement(b[b.length-1-i]);}return out;}
	private static byte[] reversed(final byte[] b){final byte[] out=new byte[b.length];for(int i=0;i<b.length;i++){out[i]=b[b.length-1-i];}return out;}
	long profileQueries,probeQueries,verificationQueries;
	/** Per-call diagnostics; EXHAUSTED is the terminal zero call, even after earlier edits.
	 * A runner cap has no zero call: do not label these last-success counters as terminal-zero. */
	enum CallStatus {SHORT_READ,UNSUPPORTED_BASE,SELF_RC,EXHAUSTED,APPLIED,APPLIED_PAIR,INITIAL_LIMIT}
	CallStatus callStatus;
	long initialEstimatedEdits;
	private int verifiedLoci;
	int lowRegions,skippedEdge,skippedWide,skippedUndefined,acceptedTroughs;
	int ambiguousLoci,noSupportedCandidateLoci,allCandidatesRejectedLoci;
	long supportedCandidates,verificationRejectedCandidates;
	long pairQueries;
	boolean usedPairLookahead;
	private LocalEditPairLookahead pairLookahead;
	private IntList pairExcludedStarts;
	LocalSingleBaseEdit.Operation lastOperation;
	int lastPosition;
	private LocalSingleBaseEdit.Operation chosenOperation;
	private int chosenPosition;private byte chosenBase;private boolean ambiguous;
	private final int k;private final HomopolymerIndelProposal.CountLookup lookup;
	private final boolean checkIndelCompetition;
	private final Kmer key;private final LocalEditKmerProbe probe;private final LocalEditTroughLocator locator;
	private final IntList counts=new IntList();private final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();
	private static final byte[] ALPHABET={'A','C','G','T'};
}
