package assemble;

import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/** Worker-local single-edit caller used by native Tadpole and BBCMS:
 * dense/sparse original depths -> local probes -> full winning-context verification
 * -> transactional edit. Sparse scanning may miss narrow unsampled troughs.
 * correctOne installs at most one edit per call and recomputes the profile each time;
 * collect instead delivers original-coordinate proposals without installing edits.
 * Canonical whole-read orientation is a correctness-first reference implementation,
 * not a claim of optimal allocation or runtime. Worker-local; immutable count table.
 * Borrowed counts, layout settings and input arrays must remain stable during a call.
 * A successful edit replaces the Read arrays; it does not mutate the earlier arrays.
 * Counters and scratch describe this call, not a whole-read transaction or concurrent snapshot.
 * Public diagnostics are maintained by this worker; callers must not assign them.
 * Read them after a call, with no concurrent use of the corrector. Initial metadata
 * validation failures retain prior diagnostics; later exceptions may leave partial values.
 * @author Fischl */
public final class LocalEditCorrector{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses one probe window, dense scanning and substitution-first candidate selection. */
	public LocalEditCorrector(final int k_, final HomopolymerIndelProposal.CountLookup lookup_){
		this(k_, lookup_, 1);
	}
	/** Selects probe-window count while retaining dense, substitution-first discovery. */
	LocalEditCorrector(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final int windows){
		this(k_, lookup_, windows, false);
	}
	/** Optionally compares indels even after finding a supported substitution. */
	public LocalEditCorrector(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final int windows, final boolean checkIndelCompetition_){
		this(k_, lookup_, windows, checkIndelCompetition_, 1);
	}
	/** Selects original-profile stride; pair lookahead overrides it with dense scanning. */
	LocalEditCorrector(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final int windows,
			final boolean checkIndelCompetition_, final int profileStride_){
		this(k_, lookup_, windows, checkIndelCompetition_, profileStride_, null);
	}
	/** Allocates reusable discovery scratch for exactly representable K>=5 and stable counts.
	 * Probe windows must be positive, odd and no greater than K; stride must be positive.
	 * A supplied neural gate requires one window and indel competition. The caller owns
	 * gate/table compatibility and beginRead setup for depth-aware scoring; this constructor
	 * neither creates that gate nor begins its read scope. Neural callers must exclude
	 * pair lookahead, whose fallback is not scored here; LocalEditEngine enforces that rule.
	 * Retains the supplied count provider and optional gate; use one corrector/gate per worker.
	 * @param k_ Exactly representable count-table word length, at least 5
	 * @param lookup_ Nonnull stable count provider with the table's K/hash/canonicalization
	 * @param windows Positive odd probe-window count no greater than K; neural scoring requires 1
	 * @param checkIndelCompetition_ Whether indels compete with substitutions; required with a gate
	 * @param profileStride_ Positive original-profile sampling stride; 1 is dense
	 * @param neuralGate_ Optional worker-local veto, or null for heuristic-only correction
	 */
	public LocalEditCorrector(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final int windows,
			final boolean checkIndelCompetition_, final int profileStride_, final LocalEditNeuralGate neuralGate_){
		if(lookup_==null || k_<5){throw new IllegalArgumentException("Corrector requires K>=5 and immutable counts.");}
		if(profileStride_<1){throw new IllegalArgumentException("Profile stride must be positive.");}
		if(neuralGate_!=null && (windows!=1 || !checkIndelCompetition_)){
			throw new IllegalArgumentException("Frozen neural local edits require windows=1 and complete indel competition.");
		}
		k=k_;
		lookup=lookup_;
		key=new Kmer(k);
		probe=new LocalEditKmerProbe(k, lookup, windows);
		locator=new LocalEditTroughLocator(k, 3);
		profileStride=profileStride_; profile=new LocalEditDepthProfile(k, lookup); counts=profile.depths;
		checkIndelCompetition=checkIndelCompetition_; neuralGate=neuralGate_;
		if(key.kbig!=k){throw new IllegalArgumentException("Correction K must equal table K.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the number of installed edits, zero or one; old arrays remain intact, but Read references may change. */
	public int correctOne(final Read read){
		return correctOne(read, false);
	}
	/** Experimental opt-in: installs only the first edit of a two-edit witness after single-edit
	 * discovery is exhausted. This concerns one read, not mates; neural callers must pass false. */
	public int correctOne(final Read read, final boolean nearbyPairs){
		return correctOne(read, nearbyPairs, -1);
	}
	/** On the first read pass, estimate separated errors from depth geometry,
	 * without probing mutations. This is a rejection heuristic, not verified repairs.
	 * Negative initialLimit preserves the ordinary first-edit discovery path. */
	public int correctOne(final Read read, final boolean nearbyPairs, final int initialLimit){
		return scan(read, nearbyPairs, initialLimit, null);
	}
	/** Propose against one immutable original profile without applying any edit.
	 * Positions and dependency bounds address the supplied canonical bases.
	 * The sink must not mutate bases or reenter this worker; false stops scanning.
	 * Returns the number of sink invocations, including a final false response;
	 * that is proposals delivered, not edits accepted or applied by the sink.
	 * This is independent-locus discovery, not replay of sequential correction.
	 * Pair witnesses are intentionally excluded from this experimental API. */
	int collect(final Read read, final int initialLimit, final ProposalSink sink){
		if(sink==null){throw new IllegalArgumentException("Non-mutating discovery requires a proposal sink.");}
		return scan(read, false, initialLimit, sink);
	}
	/** Serial callback borrowing canonical bases and a half-open dependency interval.
	 * Position/base describe the normalized edit in that orientation; reverse says whether
	 * it differs from the supplied Read. False stops after this delivered proposal. */
	interface ProposalSink{
		/** Consumes one verified proposal without mutating bases or reentering this corrector. */
		boolean accept(byte[] canonicalBases, boolean reverse, LocalSingleBaseEdit.Operation operation,
			int position, byte base, int contextFrom, int contextTo, int threshold);
	}
	/** Validates editable metadata/array shape before resetting call diagnostics.
	 * Then chooses canonical orientation, scans bounded troughs and applies or delivers
	 * proposals. Later exceptions can leave partial scratch/counters; initial validation
	 * failures retain the previous call's diagnostics. A nonnull sink disables installation. */
	private int scan(final Read read, final boolean nearbyPairs, final int initialLimit, final ProposalSink sink){
		assert(sink==null || !nearbyPairs) : "Independent original-coordinate proposals do not implement pair-witness semantics.";
		HomopolymerIndelEdit.requireEditable(read);
		if(read.bases==null || (read.quality!=null && read.quality.length!=read.bases.length)){
			throw new IllegalArgumentException("Correction requires bases with matching or null qualities.");
		}
		profileQueries=probeQueries=verificationQueries=0; lastOperation=null; lastPosition=-1;
		pairQueries=0; usedPairLookahead=false; initialEstimatedEdits=0; verifiedLoci=0;
		if(pairExcludedStarts!=null){pairExcludedStarts.clear();}
		lowRegions=skippedEdge=skippedWide=skippedUndefined=acceptedTroughs=0;
		ambiguousLoci=noSupportedCandidateLoci=allCandidatesRejectedLoci=0;
		supportedCandidates=verificationRejectedCandidates=0; callStatus=null;
		neuralEvaluations=neuralRejectedLoci=neuralAcceptedLoci=0;
		if(read.length()<k+2){callStatus=CallStatus.SHORT_READ; return 0;}
		for(byte b:read.bases){
			if(b!='A' && b!='C' && b!='G' && b!='T' && b!='N'){
				callStatus=CallStatus.UNSUPPORTED_BASE; return 0;
			}
		}
		final int orientation=orientation(read.bases);
		// A self-RC sequence cannot have a strand-invariant single nonsymmetric edit.
		if(orientation==0){callStatus=CallStatus.SELF_RC; return 0;}
		final boolean reverse=orientation>0;
		final byte[] bases=reverse ? reverseComplement(read.bases) : read.bases;
		// Pair-witness code consumes a complete dense profile; preserve that contract.
		profile.fill(bases, nearbyPairs ? 1 : profileStride); profileQueries=profile.queries;
		locator.reset(bases, counts, profile.sparse);
		if(initialLimit>=0){
			initialEstimatedEdits=estimateInitialEdits(initialLimit);
			// Reuse the same depths for ordinary discovery; preflight adds no lookups.
			locator.reset(bases, counts, profile.sparse);
			if(initialEstimatedEdits>initialLimit){finishScan(CallStatus.INITIAL_LIMIT); return 0;}
		}
		probe.beginRead(bases); int delivered=0;
		while(locator.next()){
			acceptedTroughs++; final long supportedBefore=supportedCandidates;
			chosenOperation=null; chosenPosition=chosenEvidencePosition=-1; ambiguous=false;
			final int a=locator.depthStart, end=locator.depthEnd, w=a+(end-a-1)/2;
			int originalMax=0; for(int j=a; j<end; j++){originalMax=Math.max(originalMax, counts.get(j));}
			assert(originalMax<3) : "LocalEditTroughLocator(k,3) returns only depth<3 windows; contrast multiplication is bounded.";
			final int threshold=Math.max(4, originalMax*4+1);
			for(int p=locator.baseStart; p<locator.baseEnd; p++){
				probe.substitutions(bases, w, p); probeQueries+=probe.queries;
				for(int b=0; b<4; b++){
					if(probe.substitutionDepth[b]>=threshold){
						consider(bases, LocalSingleBaseEdit.Operation.SUBSTITUTION, p,
								ALPHABET[b], a-1, end+k, threshold);
					}
				}
			}
			if(ambiguous){recordAmbiguous(a, nearbyPairs); continue;}
			//Exact as well as approximate counts can support S beside a verified indel.
			//Guarded callers must not hide that ambiguity through S-first ordering.
			if(chosenOperation==null || checkIndelCompetition){
				for(int p=locator.baseStart; p<locator.baseEnd; p++){
					probe.indelsAt(bases, w, p); probeQueries+=probe.queries;
					if(probe.deletionDepth>=threshold){consider(bases, LocalSingleBaseEdit.Operation.DELETION, p, (byte)0, a-1, end+k, threshold);}
					if(p>=locator.gapStart && p<locator.gapEnd){for(int b=0; b<4; b++){
						if(probe.insertionDepth[b]>=threshold){consider(bases, LocalSingleBaseEdit.Operation.INSERTION, p, ALPHABET[b], a-1, end+k, threshold);}
					}}
				}
			}
			if(ambiguous){recordAmbiguous(a, nearbyPairs); continue;}
			if(chosenOperation==null){
				if(supportedCandidates==supportedBefore){noSupportedCandidateLoci++;}else{allCandidatesRejectedLoci++;}
				continue;
			}
			if(neuralGate!=null){
				neuralEvaluations++;
				if(!neuralAccept(read, reverse)){
					neuralRejectedLoci++; allCandidatesRejectedLoci++; continue;
				}
				neuralAcceptedLoci++;
			}
			verifiedLoci++;
			if(sink!=null){
				delivered++;
				// Include the union of candidate verification footprints, not just
				// the winning edit; later edits must not change competing evidence.
				final int dependencyFrom=Math.max(0, Math.min(a-1, locator.baseStart-k+1));
				final int dependencyTo=(int)Math.min((long)bases.length, Math.max((long)end+k, (long)locator.baseEnd+k-1));
				if(!sink.accept(bases, reverse, chosenOperation, chosenPosition, chosenBase, dependencyFrom, dependencyTo, threshold)){
					finishScan(CallStatus.COLLECTED); return delivered;
				}
				continue;
			}
			applyChosen(read, bases, reverse); finishScan(CallStatus.APPLIED); return 1;
		}
		if(nearbyPairs){
			if(pairLookahead==null){pairLookahead=new LocalEditPairLookahead(k, lookup);}
			final boolean found=pairLookahead.propose(bases, counts, pairExcludedStarts); pairQueries=pairLookahead.queries;
			if(found){
				chosenOperation=pairLookahead.firstOperation;
				chosenPosition=pairLookahead.firstPosition;
				chosenBase=pairLookahead.firstBase;
				applyChosen(read, bases, reverse); usedPairLookahead=true; finishScan(CallStatus.APPLIED_PAIR); return 1;
			}
		}
		finishScan(sink==null ? CallStatus.EXHAUSTED : CallStatus.COLLECTED);
		return delivered;
	}
	/** Count only K-1/K-wide troughs with strong flanks and separated contexts.
	 * An isolated missing base disrupts K-1 windows; a substitution/extra base K.
	 * Short collision holes, edge/N regions and merged troughs are not counted.
	 * Coverage/variants can mimic errors: this estimates burden, not biological truth.
	 * Underestimation is handled by the engine's max+1 whole-read rollback. */
	private long estimateInitialEdits(final int limit){
		assert(limit>=0) : "A negative initial limit disables preflight in correctOne.";
		long estimated=0, lastEnd=-1;
		while(locator.next()){
			final int a=locator.depthStart, end=locator.depthEnd;
			if(end-a<k-1){continue;}
			int originalMax=0; for(int j=a; j<end; j++){originalMax=Math.max(originalMax, counts.get(j));}
			assert(originalMax<3) : "LocalEditTroughLocator(k,3) bounds low-window depth at 2, so 4*depth+1 cannot overflow.";
			final int threshold=Math.max(4, originalMax*4+1);
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
	private void recordAmbiguous(final int start, final boolean nearbyPairs){
		assert(start>=0) : "Excluded pair loci use original canonical kmer-start coordinates from the ordinary locator.";
		ambiguousLoci++;
		if(nearbyPairs){
			if(pairExcludedStarts==null){pairExcludedStarts=new IntList();}
			assert(pairExcludedStarts.size==0 || pairExcludedStarts.get(pairExcludedStarts.size-1)<start) :
				"The ordinary locator visits disjoint troughs in order; pair exclusion uses an ordered merge against the same profile.";
			pairExcludedStarts.add(start);
		}
	}
	/** Edits a temporary Read in canonical orientation, then installs complete arrays
	 * in the caller's original orientation. lastPosition uses the pre-edit original
	 * read coordinates: existing-base position or insertion boundary as appropriate. */
	private void applyChosen(final Read read, final byte[] bases, final boolean reverse){
		assert(chosenOperation!=null && chosenPosition>=0) : "Only a verified single-edit or pair-first proposal may reach transactional application.";
		final byte[] quality=reverse && read.quality!=null ? reversed(read.quality) : read.quality;
		final Read working=new Read(bases, quality, read.id, read.numericID, false);
		if(chosenOperation==LocalSingleBaseEdit.Operation.DELETION){
			editor.apply(working, chosenOperation, chosenPosition);
		}else{
			editor.apply(working, chosenOperation, chosenPosition, chosenBase);
		}
		final byte[] output=reverse ? reverseComplement(working.bases) : working.bases;
		final byte[] outputQ=reverse && working.quality!=null ? reversed(working.quality) : working.quality;
		lastOperation=chosenOperation;
		lastPosition=reverse ? read.length()-(chosenOperation==LocalSingleBaseEdit.Operation.INSERTION ? 0 : 1)-chosenPosition : chosenPosition;
		assert(outputQ==null || outputQ.length==output.length) : "Canonical correction must restore one quality per output base before installation.";
		read.bases=output; read.quality=outputQ;
	}
	/** Snapshot only the portion scanned in this call; APPLIED may stop early. */
	private void finishScan(final CallStatus status){
		callStatus=status; lowRegions=locator.lowRegions;
		skippedEdge=locator.skippedEdge; skippedWide=locator.skippedWide; skippedUndefined=locator.skippedUndefined;
		assert(lowRegions==skippedEdge+skippedWide+skippedUndefined+acceptedTroughs) :
			"Each encountered low region is skipped once or returned by the locator; counters must describe this call only.";
		assert(acceptedTroughs==ambiguousLoci+noSupportedCandidateLoci+allCandidatesRejectedLoci+verifiedLoci) :
				"Every correction-profile trough is ambiguous, unsupported, rejected or verified; depth-only preflight resets locator counters before discovery.";
	}
	/** Verifies a supported raw candidate and compares its normalized edit with the current winner.
	 * Distinct verified edits mark the trough ambiguous; equivalent placements retain the first evidence locus. */
	private void consider(final byte[] bases, final LocalSingleBaseEdit.Operation op, final int p,
			final byte b, final int from, final int to, final int threshold){
		supportedCandidates++;
		if(!verify(bases, op, p, b, from, to, threshold)){verificationRejectedCandidates++; return;}
		final int normalized=normalize(bases, op, p, b);
		if(chosenOperation==null){
			chosenOperation=op;
			chosenPosition=normalized;
			chosenEvidencePosition=p;
			chosenBase=b;
		}else if(chosenOperation!=op || chosenPosition!=normalized || chosenBase!=b){
			ambiguous=true;
		}
	}
	/** Score the raw proposal against unchanged original-orientation read bytes.
	 * The normalized canonical position remains reserved for transactional apply. */
	private boolean neuralAccept(final Read read, final boolean reverse){
		assert(neuralGate!=null && chosenOperation!=null && chosenEvidencePosition>=0) :
			"Only a unique verified raw proposal may reach the neural veto.";
		final int position=reverse ? read.length()-(chosenOperation==LocalSingleBaseEdit.Operation.INSERTION ? 0 : 1)-chosenEvidencePosition :
				chosenEvidencePosition;
		final byte base=reverse && chosenOperation!=LocalSingleBaseEdit.Operation.DELETION ? complement(chosenBase) : chosenBase;
		return neuralGate.accept(read.bases, position, chosenOperation, base);
	}
	/** Verify every affected kmer as well as the trough and its supported flank windows.
	 * Only candidates whose single-word probe passed reach this streaming check. */
	private boolean verify(final byte[] bases, final LocalSingleBaseEdit.Operation op, final int p,
			final byte b, final int contextFrom, final int contextTo, final int threshold){
		final int delta=op==LocalSingleBaseEdit.Operation.INSERTION ? 1 : op==LocalSingleBaseEdit.Operation.DELETION ? -1 : 0;
		//A narrow trough does not span every window changed by its candidate edit.
		//S changes starts [p-K+1,p]; D creates [p-K+1,p-1]; I creates [p-K+1,p].
		//Convert the last affected edited window back to an original exclusive bound.
		final int from=Math.max(0, Math.min(contextFrom, p-k+1));
		final int to=Math.min(bases.length, Math.max(contextTo, p+k-(delta>0 ? 1 : 0)));
		assert(from<=p && p<to && to<=bases.length) : "Selected edit must lie inside its complete original verification context.";
		key.clearFast(); int measured=0;
		for(int i=0; i<to-from+delta; i++){
			final int relative=p-from;
			final byte x=op==LocalSingleBaseEdit.Operation.SUBSTITUTION ? (i==relative ? b : bases[from+i]) :
				delta>0 ? (i==relative ? b : bases[from+i-(i>relative ? 1 : 0)]) : bases[from+i+(i>=relative ? 1 : 0)];
			if(x!='A' && x!='C' && x!='G' && x!='T'){return false;}key.addRight(x);
			if(key.len()>=k){measured++; if(count(true)<threshold){return false;}}
		}
		assert(measured>0) : "A winning candidate must have at least one complete verification kmer."; return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Collapses equivalent one-base indel placements without changing sequence.
	 * A/C normalize left; G/T normalize right in canonical orientation. Substitutions stay put. */
	private static int normalize(final byte[] bases, final LocalSingleBaseEdit.Operation op, int p, final byte inserted){
		if(op==LocalSingleBaseEdit.Operation.SUBSTITUTION){return p;}
		final byte b=op==LocalSingleBaseEdit.Operation.DELETION ? bases[p] : inserted;
		if(b=='A' || b=='C'){
			while(p>0 && bases[p-1]==b){p--;}
		}else if(op==LocalSingleBaseEdit.Operation.DELETION){
			while(p+1<bases.length && bases[p+1]==b){p++;}
		}else{
			while(p<bases.length && bases[p]==b){p++;}
		}
		return p;
	}
	/** Counts the scratch word, mapping absent -1 to zero and rejecting lower depths.
	 * Increments the selected counter after lookup returns; a throwing lookup is not counted. */
	private int count(final boolean verification){
		final int value=lookup.count(key); if(verification){verificationQueries++;}else{profileQueries++;}
		if(value<-1){throw new IllegalStateException("Invalid count depth: "+value);}return Math.max(0, value);
	}
	/** Compares uppercase A/C/G/T/N bytes lexicographically to their reverse complement.
	 * Negative selects the original orientation; zero identifies a self-complementary read. */
	private static int orientation(final byte[] b){
		for(int i=0; i<b.length; i++){
			final int delta=b[i]-complement(b[b.length-1-i]);
			if(delta!=0){return delta;}
		}
		return 0;
	}
	/** Complements defined uppercase bases and maps all other bytes to N. */
	private static byte complement(final byte b){return b=='A' ? (byte)'T' : b=='C' ? (byte)'G' : b=='G' ? (byte)'C' : b=='T' ? (byte)'A' : (byte)'N';}
	/** Allocates the reverse-complement sequence; does not mutate the source array. */
	private static byte[] reverseComplement(final byte[] b){
		final byte[] out=new byte[b.length];
		for(int i=0; i<b.length; i++){out[i]=complement(b[b.length-1-i]);}
		return out;
	}
	/** Allocates reversed qualities; callers pass a nonnull array. */
	private static byte[] reversed(final byte[] b){
		final byte[] out=new byte[b.length];
		for(int i=0; i<b.length; i++){out[i]=b[b.length-1-i];}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-call reported profile/probe work and completed verification lookups; failures may leave partial totals. */
	public long profileQueries, probeQueries, verificationQueries;
	/** Per-call diagnostics; EXHAUSTED is the terminal zero call, even after earlier edits.
	 * A runner cap has no zero call: do not label these last-success counters as terminal-zero. */
	public enum CallStatus{SHORT_READ, UNSUPPORTED_BASE, SELF_RC, EXHAUSTED, APPLIED, APPLIED_PAIR, INITIAL_LIMIT, COLLECTED}
	/** Last normal completion reason; null during a validated call that has not finished. */
	public CallStatus callStatus;
	/** Depth-geometry estimate, possibly stopped at initialLimit+1; not a count of proven errors. */
	public long initialEstimatedEdits;
	/** Troughs with a unique candidate surviving all configured vetoes. */
	private int verifiedLoci;
	/** Locator outcomes copied at finishScan; acceptedTroughs counts delivered geometric regions. */
	public int lowRegions, skippedEdge, skippedWide, skippedUndefined, acceptedTroughs;
	/** Mutually exclusive unsuccessful outcomes among encountered geometric regions. */
	public int ambiguousLoci;
	public int noSupportedCandidateLoci, allCandidatesRejectedLoci;
	/** Per-call gate attempts and outcomes for unique ordinary candidates, excluding pair fallback.
	 * Reset after initial read validation; earlier validation failures retain the previous values
	 * and later exceptions can leave partial totals. Observe after the worker's call, not concurrently.
	 * Diagnostic fields are maintained by this corrector; callers should not assign them. */
	public int neuralEvaluations, neuralRejectedLoci, neuralAcceptedLoci;
	/** Supported raw candidate attempts and those rejected by full-context verification. */
	public long supportedCandidates, verificationRejectedCandidates;
	/** Lookup count reported by the optional pair-witness search for this call. */
	public long pairQueries;
	/** Whether the installed edit came from the optional pair-witness fallback. */
	public boolean usedPairLookahead;
	/** Lazily allocated worker-local pair search, reused across reads. */
	private LocalEditPairLookahead pairLookahead;
	/** Increasing canonical depth starts where ambiguous single edits exclude pair reinterpretation. */
	private IntList pairExcludedStarts;
	/** Last installed edit, or null when this validated call installed none. */
	public LocalSingleBaseEdit.Operation lastOperation;
	/** Pre-edit original-orientation coordinate of that installed edit, or -1. */
	public int lastPosition;
	/** Current trough's retained operation and normalized/raw evidence positions. */
	private LocalSingleBaseEdit.Operation chosenOperation;
	private int chosenPosition, chosenEvidencePosition;
	/** Retained proposed allele in canonical orientation; unused by deletion application. */
	private byte chosenBase;
	/** Whether distinct verified edits competed in the current trough. */
	private boolean ambiguous;
	/** Exactly represented table span shared by all scratch helpers. */
	private final int k;
	/** Borrowed stable count source; lookup may refresh key caches but not its represented word. */
	private final HomopolymerIndelProposal.CountLookup lookup;
	/** Whether a substitution must also compete with supported indels. */
	private final boolean checkIndelCompetition;
	/** Borrowed optional gate; caller establishes its original-read scope. */
	private final LocalEditNeuralGate neuralGate;
	/** Reused full-context verification word. */
	private final Kmer key;
	/** Worker-local candidate counts and lazy probe cache. */
	private final LocalEditKmerProbe probe;
	/** Worker-local scan over the current original-depth profile. */
	private final LocalEditTroughLocator locator;
	/** Requested profile sampling stride; pair mode always uses one. */
	private final int profileStride;
	/** Owns reused scan depths; counts aliases its output list. */
	private final LocalEditDepthProfile profile;
	private final IntList counts;
	/** Stateless single-edit applicator, installing new arrays only after reconstruction. */
	private final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
