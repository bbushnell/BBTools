package assemble;

import stream.Read;
import ukmer.Kmer;

/** Worker-local public adapter for the conservative single-base edit kernel.
 * The caller supplies immutable counts with the same K/hash/canonicalization
 * used during loading. No table updates or paired-read metadata adaptation.
 * @author Fischl */
public final class LocalEditEngine {

	public interface CountLookup {int count(Kmer key);}

	public LocalEditEngine(final int k,final CountLookup lookup,final int windows){
		this(k,lookup,windows,false);
	}
	public LocalEditEngine(final int k,final CountLookup lookup,final int windows,final boolean checkIndelCompetition){
		this(k,lookup,windows,checkIndelCompetition,1);
	}
	/** Stride1 is dense; larger strides expand sampled low-depth regions exactly. */
	public LocalEditEngine(final int k,final CountLookup lookup,final int windows,final boolean checkIndelCompetition,final int profileStride){
		this(k,lookup,windows,checkIndelCompetition,profileStride,0);
	}
	/** Optional anchor policy; legacy overloads retain disabled0. */
	public LocalEditEngine(final int k,final CountLookup lookup,final int windows,final boolean checkIndelCompetition,final int profileStride,final int anchorMinimum){
		this(k,lookup,windows,checkIndelCompetition,profileStride,anchorMinimum,0);
	}
	/** Optional magnitude ratio; 0 disables, otherwise requires >=2. */
	public LocalEditEngine(final int k,final CountLookup lookup,final int windows,final boolean checkIndelCompetition,final int profileStride,final int anchorMinimum,final int magnitudeFactor){
		if(lookup==null){throw new IllegalArgumentException("Local edits require immutable counts.");}
		corrector=new LocalEditCorrector(k,new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return lookup.count(key);}
		},windows,checkIndelCompetition,profileStride,anchorMinimum,magnitudeFactor);
	}

	/** Whole-read transaction: skip excessive initial depth-estimated burden, or
	 * roll back if later discovery requires more than maxEdits. */
	public int correct(final Read read,final int maxEdits){
		return correct(read,maxEdits,false);
	}
	/** Pair witnesses concern two edits in ONE unpaired read, not paired reads. */
	public int correct(final Read read,final int maxEdits,final boolean nearbyPairs){
		if(maxEdits<1){throw new IllegalArgumentException("fixindelsmax must be positive.");}
		HomopolymerIndelEdit.requireEditable(read);
		final byte[] originalBases=read.bases,originalQuality=read.quality;
		int applied=0;long s=0,i=0,d=0;boolean commit=false;
		try{
			while(true){
				final int changed=corrector.correctOne(read,nearbyPairs,applied==0 ? maxEdits : -1);
				profileQueries+=corrector.profileQueries;probeQueries+=corrector.probeQueries;
				verificationQueries+=corrector.verificationQueries;pairQueries+=corrector.pairQueries;
				anchorEvaluatedLoci+=corrector.anchorEvaluatedLoci;anchorRejectedLoci+=corrector.anchorRejectedLoci;
				magnitudeEvaluatedLoci+=corrector.magnitudeEvaluatedLoci;magnitudeRejectedLoci+=corrector.magnitudeRejectedLoci;
				if(corrector.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT){initialSkippedReads++;return 0;}
				if(changed==0){
					substitutions+=s;insertions+=i;deletions+=d;
					if(applied>0){changedReads++;}
					commit=true;return applied;
				}
				assert(changed==1 && corrector.lastOperation!=null) : "Only one verified edit per discovery pass may enter the read transaction.";
				attemptedEdits++;
				if(applied==maxEdits){cappedReads++;rolledBackReads++;return 0;}
				applied++;
				switch(corrector.lastOperation){
					case SUBSTITUTION: s++;break;
					case INSERTION: i++;break;
					case DELETION: d++;break;
					default: throw new AssertionError("Unknown local-edit operation.");
				}
			}
		}finally{
			// Editors install new arrays; original references preserve all original
			// bases/qualities without an extra whole-read copy. Exceptions also restore.
			if(!commit){read.bases=originalBases;read.quality=originalQuality;}
		}
	}

	/** Worker-lifetime counters; merge only after the worker has joined. */
	public long substitutions,insertions,deletions,changedReads,cappedReads;
	public long profileQueries,probeQueries,verificationQueries;
	public long initialSkippedReads,rolledBackReads,attemptedEdits,pairQueries;
	/** Work counters, including rescans and transactions later rolled back. */
	public long anchorEvaluatedLoci,anchorRejectedLoci;
	public long magnitudeEvaluatedLoci,magnitudeRejectedLoci;
	private final LocalEditCorrector corrector;
}
