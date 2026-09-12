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
		if(lookup==null){throw new IllegalArgumentException("Local edits require immutable counts.");}
		corrector=new LocalEditCorrector(k,new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return lookup.count(key);}
		},windows,checkIndelCompetition);
	}

	/** Rescan after each edit. Hitting the cap is NOT evidence of exhaustion. */
	public int correct(final Read read,final int maxEdits){
		if(maxEdits<1){throw new IllegalArgumentException("localeditmax must be positive.");}
		int applied=0;
		while(applied<maxEdits){
			final int changed=corrector.correctOne(read);
			profileQueries+=corrector.profileQueries;probeQueries+=corrector.probeQueries;
			verificationQueries+=corrector.verificationQueries;
			if(changed==0){break;}
			assert(changed==1 && corrector.lastOperation!=null) : "The single-edit kernel commits one verified edit before recomputing its depth profile.";
			applied++;
			switch(corrector.lastOperation){
				case SUBSTITUTION: substitutions++;break;
				case INSERTION: insertions++;break;
				case DELETION: deletions++;break;
				default: throw new AssertionError("Unknown local-edit operation.");
			}
		}
		if(applied>0){changedReads++;}
		if(applied==maxEdits){cappedReads++;}
		return applied;
	}

	/** Worker-lifetime counters; merge only after the worker has joined. */
	public long substitutions,insertions,deletions,changedReads,cappedReads;
	public long profileQueries,probeQueries,verificationQueries;
	private final LocalEditCorrector corrector;
}
