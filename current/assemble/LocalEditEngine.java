package assemble;

import ml.CellNet;
import stream.Read;
import ukmer.Kmer;

/** Worker-local whole-read transaction over repeated conservative single-base edits.
 * Reuse serially across reads; do not reenter or share mutable scratch/networks.
 * The caller supplies stable counts with the same K/hash/canonicalization used
 * during loading. No table updates or paired-read metadata adaptation occur.
 * Accepted edits replace arrays; exceeding the edit budget or an exception restores
 * the original bases/quality references. Work counters are not rolled back.
 * @author Fischl */
public final class LocalEditEngine{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Stable source of nonnegative depths, or -1 for absence. Borrowed Kmer scratch
	 * must neither be retained nor mutated by the count provider. */
	public interface CountLookup{int count(Kmer key);}

	/** Constructs a dense engine that checks indels only when no substitution is chosen. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows){
		this(k, lookup, windows, false);
	}
	/** Constructs a dense engine; true also checks indels against a chosen substitution. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition){
		this(k, lookup, windows, checkIndelCompetition, 1);
	}
	/** Stride1 is dense; larger strides refine sampled low-depth regions but may miss
	 * unsampled troughs. Pair lookahead always requests a dense profile. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition, final int profileStride){
		this(k, lookup, windows, checkIndelCompetition, profileStride, null, Float.NaN);
	}
	/** Optional model vetoes an already unique, fully verified candidate. The network
	 * is retained, not copied; frozen weights do not make inference scratch shareable. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition,
			final int profileStride, final CellNet neuralNet, final float neuralCutoff){
		this(k, lookup, windows, checkIndelCompetition, profileStride, neuralNet, neuralCutoff, null, Float.NaN);
	}
	/** Constructs a serial engine with optional low/high models routed by original-read depth.
	 * Loaded model dimensions/cutoffs are checked here, not model identity or provenance;
	 * callers that require a pinned resource must validate it before construction.
	 * Networks are retained and have simdFF enabled; supply worker-local instances.
	 * @param k Effective count-table K, at least5; neural scoring requires62
	 * @param lookup Stable count provider retained for the engine lifetime
	 * @param windows Positive odd probe-window count no greater than K; neural scoring requires1
	 * @param checkIndelCompetition Whether indels compete with a supported substitution; neural scoring requires true
	 * @param profileStride Positive original-profile sampling stride;1 is dense
	 * @param neuralNet Low-depth or single network, or null to disable neural scoring
	 * @param neuralCutoff Inclusive score cutoff in [0,1] when neuralNet is present
	 * @param neuralHighNet Optional39-input high-depth model; requires a39-input low model
	 * @param neuralHighCutoff Finite [0,1] cutoff with a high model, otherwise NaN when a gate is constructed
	 */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition,
			final int profileStride, final CellNet neuralNet, final float neuralCutoff,
			final CellNet neuralHighNet, final float neuralHighCutoff){
		if(lookup==null){throw new IllegalArgumentException("Local edits require immutable counts.");}
		if(neuralHighNet!=null && neuralNet==null){
			throw new IllegalArgumentException("High-depth neural local edits require the low-depth model.");
		}
		if((neuralNet!=null || neuralHighNet!=null) && (windows!=1 || !checkIndelCompetition)){
			throw new IllegalArgumentException("Frozen neural local edits require windows=1 and complete indel competition.");
		}
		final HomopolymerIndelProposal.CountLookup adapted=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return lookup.count(key);}
		};
		final LocalEditNeuralGate gate=neuralNet==null ? null :
			new LocalEditNeuralGate(k, adapted, neuralNet, neuralCutoff, neuralHighNet, neuralHighCutoff);
		neuralGate=gate;
		neuralEnabled=neuralNet!=null;
		corrector=new LocalEditCorrector(k, adapted, windows, checkIndelCompetition, profileStride, gate);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs a whole-read transaction without nearby-pair lookahead; see the overload. */
	public int correct(final Read read, final int maxEdits){
		return correct(read, maxEdits, false);
	}
	/** Recomputes discovery after each edit and commits only after discovery stops.
	 * An excessive initial depth-based burden estimate skips the read; that estimate
	 * is a rejection heuristic, not a count of verified or biological errors.
	 * Discovering and applying edit maxEdits+1 restores both original array references
	 * and returns0. Exceptions in discovery also restore them, then propagate.
	 * Pair witnesses concern two edits in ONE read; only their first edit is applied
	 * per discovery pass. Neural scoring cannot be combined with pair lookahead.
	 * @param read Unpaired, alignment-free nucleotide read without attached coordinate
	 * metadata; requires bases and matching or null qualities. No other owner may mutate
	 * it during this call; this engine intentionally replaces its bases/quality arrays
	 * @param maxEdits Positive maximum committed edits for this whole read
	 * @param nearbyPairs Enable experimental two-edit witnesses after single-edit discovery
	 * @return Number of committed edits;0 includes unchanged, initially skipped and rolled-back reads
	 */
	public int correct(final Read read, final int maxEdits, final boolean nearbyPairs){
		if(maxEdits<1){throw new IllegalArgumentException("fixindelsmax must be positive.");}
		if(neuralEnabled && nearbyPairs){throw new IllegalArgumentException("Frozen neural local edits do not support pair lookahead.");}
		HomopolymerIndelEdit.requireEditable(read);
		final byte[] originalBases=read.bases, originalQuality=read.quality;
		if(neuralGate!=null){neuralGate.beginRead(originalBases);}
		int applied=0; long s=0, i=0, d=0; boolean commit=false;
		try{
			while(true){
				final int changed=corrector.correctOne(read, nearbyPairs, applied==0 ? maxEdits : -1);
				profileQueries+=corrector.profileQueries; probeQueries+=corrector.probeQueries;
				verificationQueries+=corrector.verificationQueries; pairQueries+=corrector.pairQueries;
				neuralEvaluations+=corrector.neuralEvaluations; neuralAcceptedLoci+=corrector.neuralAcceptedLoci;
				neuralRejectedLoci+=corrector.neuralRejectedLoci;
				if(corrector.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT){initialSkippedReads++; return 0;}
				if(changed==0){
					substitutions+=s; insertions+=i; deletions+=d;
					if(applied>0){changedReads++;}
					commit=true; return applied;
				}
				assert(changed==1 && corrector.lastOperation!=null) : "Only one verified edit per discovery pass may enter the read transaction.";
				attemptedEdits++;
				if(applied==maxEdits){cappedReads++; rolledBackReads++; return 0;}
				applied++;
				switch(corrector.lastOperation){
					case SUBSTITUTION: s++; break;
					case INSERTION: i++; break;
					case DELETION: d++; break;
					default: throw new AssertionError("Unknown local-edit operation.");
				}
			}
		}finally{
			// Editors install new arrays; original references preserve all original
			// bases/qualities without an extra whole-read copy. Exceptions also restore.
			if(!commit){read.bases=originalBases; read.quality=originalQuality;}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Worker-lifetime committed edit totals/changed reads, plus budget-exceeded reads.
	 * These fields are cumulative; merge only after the worker has joined. */
	public long substitutions, insertions, deletions, changedReads, cappedReads;
	/** Lookup work accumulated after completed discovery calls, including rolled-back reads.
	 * Work from a discovery call that throws before returning is not accumulated here. */
	public long profileQueries, probeQueries, verificationQueries;
	/** Initial heuristic skips, budget rollbacks, attempted edits (including max+1),
	 * and pair lookup work. rolledBackReads excludes exception-triggered restoration. */
	public long initialSkippedReads, rolledBackReads, attemptedEdits, pairQueries;
	/** Completed-call neural decisions, including candidates in later rolled-back reads. */
	public long neuralEvaluations, neuralAcceptedLoci, neuralRejectedLoci;
	/** Whether construction supplied a low/single neural network. */
	private final boolean neuralEnabled;
	/** Optional retained network gate and original-read context, reused serially. */
	private final LocalEditNeuralGate neuralGate;
	/** Reused discovery scratch; installs at most one edit per completed call. */
	private final LocalEditCorrector corrector;
}
