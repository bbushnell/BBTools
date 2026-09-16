package assemble;

import java.util.Arrays;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/** Unrouted local sequential prototype with explicit ORIGINAL-coordinate
 * regions. Each quadruple is [copyFrom,copyTo,coreFrom,coreTo), with at least
 * K immutable context bases on either side of the mutable core. Automatic
 * region discovery/grouping is deliberately not implemented here.
 * @author Fischl */
final class LocalEditPatchCorrector {

	LocalEditPatchCorrector(final int k_,final HomopolymerIndelProposal.CountLookup lookup_,final int windows,
		final boolean competition,final int maxPatchBases_){
		if(k_<5 || lookup_==null || maxPatchBases_<2L*k_+1){throw new IllegalArgumentException("Local regions require K>=5, immutable counts, and bounded space for two full flanks.");}
		k=k_;lookup=lookup_;maxPatchBases=maxPatchBases_;
		windowCount=windows;checkCompetition=competition;
		preflight=new LocalEditCorrector(k,lookup,windows,competition,1);
		local=new LocalEditCorrector(k,lookup,windows,competition,1);key=new Kmer(k);
	}
	/** Owns a fresh original-read selection, preserving the initial guard while
	 * avoiding a second full preflight scan. There is no caller-exposed switch
	 * to skip guards, and an old selection cannot be supplied through this API. */
	int correctSelected(final Read input,final LocalEditPatchRegions selector,final int maxEdits){
		clearCounters();transaction.clear();
		if(selector==null || !selector.matchesDenseExecutor(k,lookup,windowCount,checkCompetition,maxPatchBases)){
			throw new IllegalArgumentException("Shared preflight requires identical K, immutable count pool, windows, competition, dense stride and patch bound.");
		}
		final IntList regions=selector.select(input,maxEdits);
		try{return correct(input,regions,maxEdits,false);}
		finally{
			initialRejected|=selector.initialRejected;
			profileQueries+=selector.profileQueries;probeQueries+=selector.probeQueries;verificationQueries+=selector.verificationQueries;
		}
	}
	int correct(final Read input,final IntList regions,final int maxEdits){
		return correct(input,regions,maxEdits,true);
	}
	private int correct(final Read input,final IntList regions,final int maxEdits,final boolean preflightRequired){
		clearCounters();transaction.clear();
		if(maxEdits<1){throw new IllegalArgumentException("The original read's edit limit must be positive.");}
		HomopolymerIndelEdit.requireEditable(input);validateRegions(input,regions);
		if(regions.size==0){return 0;}
		final byte[] originalBases=input.bases,originalQuality=input.quality;
		try{
			if(preflightRequired){
				// Explicit supplied regions have not necessarily passed an original
				// burden guard. Only correctSelected owns evidence to omit this pass.
				preflight.collect(input,maxEdits,(b,reverse,op,p,base,from,to,threshold)->false);
				profileQueries+=preflight.profileQueries;probeQueries+=preflight.probeQueries;verificationQueries+=preflight.verificationQueries;
				initialRejected=preflight.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT;
				if(initialRejected || preflight.callStatus==LocalEditCorrector.CallStatus.SHORT_READ ||
					preflight.callStatus==LocalEditCorrector.CallStatus.UNSUPPORTED_BASE || preflight.callStatus==LocalEditCorrector.CallStatus.SELF_RC){return 0;}
			}
			transaction.reset(input);
			for(int i=0;i<regions.size;i+=4){
				final int from=regions.array[i],to=regions.array[i+1],coreStart=regions.array[i+2]-from,originalCoreEnd=regions.array[i+3]-from;
				if(to-from>maxPatchBases){sizeRejected++;continue;}
				final byte[] source=Arrays.copyOfRange(originalBases,from,to);
				final byte[] sourceQuality=originalQuality==null ? null : Arrays.copyOfRange(originalQuality,from,to);
				copiedScratchBases+=source.length;
				working.bases=source;working.quality=sourceQuality;working.id=input.id;working.numericID=input.numericID;
				int coreEnd=originalCoreEnd,s=0,a=0,d=0;boolean rejected=false;
				while(true){
					final int changed=local.correctOne(working,false,-1);
					profileQueries+=local.profileQueries;probeQueries+=local.probeQueries;verificationQueries+=local.verificationQueries;
					if(changed==0){break;}
					assert(changed==1 && local.lastOperation!=null) : "The local transaction accounts for one native S/I/D per successful discovery pass.";
					attemptedEdits++;
					final LocalSingleBaseEdit.Operation op=local.lastOperation;final int p=local.lastPosition;
					if(p<coreStart || (op==LocalSingleBaseEdit.Operation.INSERTION ? p>coreEnd : p>=coreEnd)){
						boundaryRejected++;rejected=true;break;
					}
					if((long)transaction.editCount()+s+a+d>=maxEdits){budgetRejected=true;return 0;}
					if(op==LocalSingleBaseEdit.Operation.SUBSTITUTION){s++;}
					else if(op==LocalSingleBaseEdit.Operation.INSERTION){a++;coreEnd++;}
					else if(op==LocalSingleBaseEdit.Operation.DELETION){d++;coreEnd--;}
					else{throw new IllegalStateException("Unexpected local operation.");}
					if(working.length()>maxPatchBases){sizeRejected++;rejected=true;break;}
				}
				if(rejected || s+a+d==0){continue;}
				// These are bookkeeping invariants, not soft scientific filters:
				// an edit confined to the core cannot alter either context flank.
				if(!sameFlanks(source,working.bases,coreStart,originalCoreEnd,coreEnd) ||
					(sourceQuality!=null && !sameFlanks(sourceQuality,working.quality,coreStart,originalCoreEnd,coreEnd))){
					throw new IllegalStateException("Local edits escaped their owned core or corrupted flank bookkeeping.");
				}
				if(!verifyFinal(source,working.bases)){finalRejected++;continue;}
				transaction.add(from,to,working.bases,working.quality,s,a,d);stagedPatches++;
			}
			if(input.bases!=originalBases || input.quality!=originalQuality){throw new IllegalStateException("Original arrays changed while local patches were being corrected.");}
			final int edits=transaction.commit(maxEdits);
			if(edits<0){throw new IllegalStateException("The per-patch remaining-budget check failed to protect the whole-read transaction.");}
			substitutions=transaction.substitutions;insertions=transaction.insertions;deletions=transaction.deletions;
			materializedBases=transaction.materializedBases;copiedOriginalBases=transaction.copiedOriginalBases;copiedReplacementBases=transaction.copiedReplacementBases;
			return edits;
		}finally{transaction.clear();working.bases=null;working.quality=null;working.id=null;}
	}
	private void validateRegions(final Read input,final IntList regions){
		if(input.bases==null || (input.quality!=null && input.quality.length!=input.length()) || regions==null || regions.size%4!=0){throw new IllegalArgumentException("Expected intact input and four original-coordinate fields per region.");}
		int previousEnd=-1;
		for(int i=0;i<regions.size;i+=4){
			final int from=regions.array[i],to=regions.array[i+1],a=regions.array[i+2],b=regions.array[i+3];
			if(from<0 || to>input.length() || from>=to || from<previousEnd || a<from || b>to || a>b || a-from<k || to-b<k){
				throw new IllegalArgumentException("Regions must be disjoint original intervals with a contained core and K context bases on each side.");
			}
			previousEnd=to;
		}
	}
	private boolean verifyFinal(final byte[] before,final byte[] after){
		assert(before!=null && after!=null) : "Final validation only sees a completed local candidate and its immutable original scratch.";
		int prefix=0,suffix=0;
		while(prefix<before.length && prefix<after.length && before[prefix]==after[prefix]){prefix++;}
		while(suffix<before.length-prefix && suffix<after.length-prefix && before[before.length-1-suffix]==after[after.length-1-suffix]){suffix++;}
		// Outside this envelope every word is literally unchanged. Within it,
		// the prototype deliberately requires full support, even for an unchanged
		// internal island. Each native step already enforced its own >=4 threshold.
		final int first=Math.max(0,prefix-k+1),last=Math.min(after.length-k,after.length-suffix-1);
		if(first>last){return true;}
		key.clearFast();
		for(int p=first;p<last+k;p++){
			final byte b=after[p];if(b!='A' && b!='C' && b!='G' && b!='T'){return false;}
			key.addRight(b);
			if(key.len()>=k){
				final int count=lookup.count(key);verificationQueries++;
				if(count<-1){throw new IllegalStateException("Invalid immutable count in final patch validation.");}
				if(count<4){return false;}
			}
		}
		return true;
	}
	private static boolean sameFlanks(final byte[] before,final byte[] after,final int coreStart,final int oldEnd,final int newEnd){
		assert(before!=null && after!=null && coreStart>=0 && oldEnd>=coreStart && newEnd>=coreStart) : "Core endpoints evolve only by accepted local insertions/deletions.";
		if(after.length-newEnd!=before.length-oldEnd){return false;}
		for(int p=0;p<coreStart;p++){if(before[p]!=after[p]){return false;}}
		for(int p=0;p<before.length-oldEnd;p++){if(before[oldEnd+p]!=after[newEnd+p]){return false;}}
		return true;
	}
	private void clearCounters(){
		initialRejected=budgetRejected=false;boundaryRejected=sizeRejected=finalRejected=stagedPatches=0;
		substitutions=insertions=deletions=0;profileQueries=probeQueries=verificationQueries=attemptedEdits=0;
		copiedScratchBases=materializedBases=copiedOriginalBases=copiedReplacementBases=0;
	}
	boolean initialRejected,budgetRejected;
	int boundaryRejected,sizeRejected,finalRejected,stagedPatches,substitutions,insertions,deletions;
	long profileQueries,probeQueries,verificationQueries,attemptedEdits,copiedScratchBases,materializedBases,copiedOriginalBases,copiedReplacementBases;
	private final int k,maxPatchBases;private final HomopolymerIndelProposal.CountLookup lookup;
	private final int windowCount;private final boolean checkCompetition;
	private final LocalEditCorrector preflight,local;private final Kmer key;
	private final LocalEditPatchBatch transaction=new LocalEditPatchBatch();
	private final Read working=new Read(new byte[0],null,"local-working",0,false);
}
