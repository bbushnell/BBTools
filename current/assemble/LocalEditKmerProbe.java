package assemble;

import java.util.Arrays;
import ukmer.Kmer;

/** Reusable, non-mutating 3+1+4 local-edit count probes at a supplied position.
 * This enumerates evidence only: no localization, acceptance, or read editing.
 * Caller may stop after substitutions; indel probes are lazy and cached.
 * @author Fischl */
final class LocalEditKmerProbe {
	LocalEditKmerProbe(final int k_,final HomopolymerIndelProposal.CountLookup lookup_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Probe requires K>=5 and an immutable count lookup.");}
		k=k_;lookup=lookup_;key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Probe K must equal effective table K.");}
	}
	/** Original k-base window [windowStart,windowStart+k), position inside it.
	 * Requires a called A/C/G/T at position (three, not four substitutions).
	 * Each other undefined base makes a particular candidate window unavailable.
	 * The input/table must not change until the next call to substitutions. */
	void substitutions(final byte[] bases_,final int windowStart_,final int position_){
		final int original=prepare(bases_,windowStart_,position_);
		key.clearFast();
		for(int i=windowStart;i<windowStart+k;i++){
			if(code(bases[i])<0){substitutionsComplete=true;return;}key.addRight(bases[i]);
		}
		substitutionsAvailable=true;
		for(int b=0;b<4;b++){
			if(b==original){continue;}
			final long previous=key.substituteBase(position-windowStart,b);
			assert(previous==original) : "Each candidate must start from the original built window; the prior probe must restore it.";
			substitutionDepth[b]=count();
			key.substituteBase(position-windowStart,original);
		}
		assert(queries==3) : "A fully defined called base has exactly three substitution alternatives.";
		substitutionsComplete=true;
	}
	/** After a caller's multi-position substitution stage found no repair, probe
	 * only the five indels here; do not repeat the three substitution lookups. */
	void indelsAt(final byte[] bases_,final int windowStart_,final int position_){
		prepare(bases_,windowStart_,position_);substitutionsComplete=true;indels();
	}
	private int prepare(final byte[] bases_,final int windowStart_,final int position_){
		substitutionsComplete=false;bases=null;
		if(bases_==null || windowStart_<0 || (long)windowStart_+k>bases_.length || position_<windowStart_ || position_>=windowStart_+k){
			throw new IllegalArgumentException("Substitution probes require an in-bounds kmer and a position inside it.");
		}
		final int original=code(bases_[position_]);
		if(original<0){throw new IllegalArgumentException("Three-alternative probe requires a called A/C/G/T base.");}
		bases=bases_;windowStart=windowStart_;position=position_;queries=0;indelsEvaluated=false;
		Arrays.fill(substitutionDepth,-1);Arrays.fill(insertionDepth,-1);deletionDepth=-1;
		substitutionsAvailable=deletionAvailable=insertionsAvailable=false;
		return original;
	}
	/** Deletion drops bases[position] and draws in bases[windowStart+k].
	 * Insertion inserts before position and trims the original last window base.
	 * Each resulting candidate is exactly k bases. No candidate sequence allocation. */
	void indels(){
		if(!substitutionsComplete){throw new IllegalStateException("Successfully probe substitutions first; failed probes require a fresh input call.");}
		if(indelsEvaluated){return;}
		try{probeIndels();indelsEvaluated=true;}
		catch(RuntimeException failure){substitutionsComplete=false;throw failure;}
		catch(Error failure){substitutionsComplete=false;throw failure;}
	}
	private void probeIndels(){
		if((long)windowStart+k<bases.length){
			key.clearFast();boolean valid=true;
			for(int i=windowStart;i<=windowStart+k;i++){
				if(i==position){continue;}
				if(code(bases[i])<0){valid=false;break;}key.addRight(bases[i]);
			}
			if(valid){deletionAvailable=true;deletionDepth=count();}
		}
		key.clearFast();
		for(int i=0;i<k;i++){
			final int relative=position-windowStart;
			final byte b=i==relative ? (byte)'A' : bases[windowStart+i-(i>relative ? 1 : 0)];
			if(code(b)<0){return;}key.addRight(b);
		}
		insertionsAvailable=true;
		for(int b=0;b<4;b++){
			key.substituteBase(position-windowStart,b);insertionDepth[b]=count();
		}
		assert(queries<=8) : "One locus requires at most three substitutions, one deletion and four inserted-base lookups.";
	}
	private int count(){
		assert(key.len()>=k) : "Only complete edited kmers may reach the count table.";
		final int depth=lookup.count(key);queries++;
		if(depth<-1){throw new IllegalStateException("CountLookup permits depth>=0 or absent=-1, got "+depth);}
		return Math.max(0,depth);
	}
	private static int code(final byte b){return b=='A' ? 0 : b=='C' ? 1 : b=='G' ? 2 : b=='T' ? 3 : -1;}
	/** -1 means unavailable/not probed; a queried absent kmer has depth zero. */
	final int[] substitutionDepth=new int[4],insertionDepth=new int[4];
	int deletionDepth,queries;
	boolean substitutionsAvailable,deletionAvailable,insertionsAvailable;
	private boolean indelsEvaluated,substitutionsComplete;
	private byte[] bases;
	private int windowStart,position;
	private final int k;
	private final Kmer key;
	private final HomopolymerIndelProposal.CountLookup lookup;
}
