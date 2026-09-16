package assemble;

import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import shared.Shared;

/** Worker-reused transaction for already validated original-coordinate patches.
 * Does NOT discover corrections or validate biological support. Replacements
 * include their final quality bytes; all local S/I/D operations count against
 * the read-wide budget, not merely the number of replacement intervals.
 * Original arrays remain immutable until one complete final installation.
 * @author Fischl */
final class LocalEditPatchBatch {

	void reset(final Read input){
		clear();clearCounters();HomopolymerIndelEdit.requireEditable(input);
		if(input.bases==null || (input.quality!=null && input.quality.length!=input.bases.length)){
			throw new IllegalArgumentException("Patch transaction requires intact bases and matching or null qualities.");
		}
		read=input;bases=input.bases;quality=input.quality;
	}
	/** Copy a validated replacement into worker-owned staging storage. Intervals
	 * [from,to) must be disjoint and starts strictly increase; from==to inserts.
	 * Same-start operations must be combined into one patch by the caller.
	 * Any append failure abandons the entire transaction, never a partial list. */
	void add(final int from,final int to,final byte[] replacement,final byte[] replacementQuality,
		final int subs,final int ins,final int dels){
		boolean success=false;
		try{
			if(read==null){throw new IllegalStateException("Reset before staging a patch.");}
			if(from<0 || to<from || to>bases.length || replacement==null){throw new IllegalArgumentException("Patch interval/replacement outside the original read: "+from+".."+to);}
			if((quality==null)!=(replacementQuality==null) || (replacementQuality!=null && replacementQuality.length!=replacement.length)){
				throw new IllegalArgumentException("A patch must preserve quality nullness and supply one final quality per replacement base.");
			}
			final long operations=(long)subs+ins+dels,total=(long)stagedSubs+stagedIns+stagedDels+operations;
			if(subs<0 || ins<0 || dels<0 || operations<1 || total>Integer.MAX_VALUE){throw new IllegalArgumentException("Patch operation counts must be positive in total and fit the read-level counter.");}
			if((long)replacement.length-(to-from)!=(long)ins-dels){throw new IllegalArgumentException("Replacement length must agree with the actual insertion/deletion counts.");}
			if(spans.size>0 && (from<=spans.array[spans.size-4] || from<spans.array[spans.size-3])){
				throw new IllegalArgumentException("Patches overlap or share an original start; combine dependent edits before staging.");
			}
			if((long)replacementBases.length()+replacement.length>Shared.MAX_ARRAY_LEN){throw new IllegalArgumentException("Replacement staging exceeds the supported primitive array size.");}
			final int offset=replacementBases.length();
			replacementBases.append(replacement);
			if(quality!=null){replacementQualities.append(replacementQuality);}
			spans.add(from);spans.add(to);spans.add(offset);spans.add(replacement.length);
			stagedSubs+=subs;stagedIns+=ins;stagedDels+=dels;
			assert(spans.size%4==0 && (quality==null || replacementBases.length()==replacementQualities.length())) :
				"Each staged interval has four primitive fields and matching base/quality pool offsets.";
			success=true;
		}finally{if(!success){clear();}}
	}
	/** Return committed operation count, -1 if over budget, or 0 for no patches.
	 * Always closes the transaction. Failure never overwrites original arrays. */
	int commit(final int maxEdits){
		if(read==null){throw new IllegalStateException("No active patch transaction.");}
		clearCounters();
		try{
			if(maxEdits<0){throw new IllegalArgumentException("Read-level edit limit must be nonnegative.");}
			HomopolymerIndelEdit.requireEditable(read);
			if(read.bases!=bases || read.quality!=quality){throw new IllegalStateException("Read arrays changed during patch staging.");}
			final int operations=stagedSubs+stagedIns+stagedDels;
			if(operations>maxEdits){return -1;}
			if(spans.size==0){return 0;}
			final int length=LocalEditBatch.resultLength(bases.length,(long)stagedIns-stagedDels);
			final byte[] output=new byte[length],outputQuality=quality==null ? null : new byte[length];
			int sourceEnd=bases.length,destinationEnd=length;
			long originalCopied=0,replacementCopied=0;
			for(int i=spans.size-4;i>=0;i-=4){
				final int from=spans.array[i],to=spans.array[i+1],offset=spans.array[i+2],len=spans.array[i+3];
				final int tailLength=sourceEnd-to;destinationEnd-=tailLength;
				assert(tailLength>=0 && destinationEnd>=len) : "Ordered disjoint spans and net indel accounting must leave room for each tail and replacement.";
				System.arraycopy(bases,to,output,destinationEnd,tailLength);
				if(quality!=null){System.arraycopy(quality,to,outputQuality,destinationEnd,tailLength);}
				destinationEnd-=len;
				System.arraycopy(replacementBases.array,offset,output,destinationEnd,len);
				if(quality!=null){System.arraycopy(replacementQualities.array,offset,outputQuality,destinationEnd,len);}
				originalCopied+=tailLength;replacementCopied+=len;sourceEnd=from;
			}
			destinationEnd-=sourceEnd;
			assert(destinationEnd==0) : "Backward interval splicing must consume exactly the final array length.";
			System.arraycopy(bases,0,output,0,sourceEnd);
			if(quality!=null){System.arraycopy(quality,0,outputQuality,0,sourceEnd);}
			originalCopied+=sourceEnd;
			assert(originalCopied+replacementCopied==length) : "Each final base must come from exactly one original span or staged replacement.";
			copiedOriginalBases=originalCopied;copiedReplacementBases=replacementCopied;materializedBases=length;
			substitutions=stagedSubs;insertions=stagedIns;deletions=stagedDels;
			read.bases=output;read.quality=outputQuality;
			return operations;
		}finally{clear();}
	}
	void clear(){
		spans.clear();replacementBases.clear();replacementQualities.clear();
		read=null;bases=null;quality=null;stagedSubs=stagedIns=stagedDels=0;
	}
	int size(){return spans.size/4;}
	int editCount(){return stagedSubs+stagedIns+stagedDels;}
	private void clearCounters(){copiedOriginalBases=copiedReplacementBases=materializedBases=0;substitutions=insertions=deletions=0;}
	long copiedOriginalBases,copiedReplacementBases,materializedBases;
	int substitutions,insertions,deletions;
	private int stagedSubs,stagedIns,stagedDels;
	private Read read;private byte[] bases,quality;
	private final IntList spans=new IntList(32);
	private final ByteBuilder replacementBases=new ByteBuilder(256),replacementQualities=new ByteBuilder(256);
}
