package assemble;

import ml.CellNet;
import stream.Read;
import ukmer.Kmer;

/** Worker-local public adapter for the conservative single-base edit kernel.
 * The caller supplies immutable counts with the same K/hash/canonicalization
 * used during loading. No table updates or paired-read metadata adaptation.
 * @author Fischl */
public final class LocalEditEngine {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable kmer-count source used by the local edit engine. */
	public interface CountLookup {int count(Kmer key);}

	/** Constructs an engine with sequential single-edit competition disabled. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows){
		this(k, lookup, windows, false);
	}
	/** Constructs an engine with explicit indel-competition behavior. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition){
		this(k, lookup, windows, checkIndelCompetition, 1);
	}
	/** Stride1 is dense; larger strides expand sampled low-depth regions exactly. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition, final int profileStride){
		this(k, lookup, windows, checkIndelCompetition, profileStride, null, Float.NaN);
	}
	/** Optional frozen model is a veto after unique full-context verification. */
	public LocalEditEngine(final int k, final CountLookup lookup, final int windows, final boolean checkIndelCompetition,
			final int profileStride, final CellNet neuralNet, final float neuralCutoff){
		this(k, lookup, windows, checkIndelCompetition, profileStride, neuralNet, neuralCutoff, null, Float.NaN);
	}
	/** Optional frozen low/high-depth models route each original read by median kmer depth. */
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

	/** Whole-read transaction: skip excessive initial depth-estimated burden, or
	 * roll back if later discovery requires more than maxEdits. */
	public int correct(final Read read, final int maxEdits){
		return correct(read, maxEdits, false);
	}
	/** Pair witnesses concern two edits in ONE unpaired read, not paired reads. */
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

	/** Worker-lifetime counters; merge only after the worker has joined. */
	public long substitutions, insertions, deletions, changedReads, cappedReads;
	public long profileQueries, probeQueries, verificationQueries;
	public long initialSkippedReads, rolledBackReads, attemptedEdits, pairQueries;
	public long neuralEvaluations, neuralAcceptedLoci, neuralRejectedLoci;
	private final boolean neuralEnabled;
	private final LocalEditNeuralGate neuralGate;
	private final LocalEditCorrector corrector;
}
