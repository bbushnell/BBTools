package bbdukGemini;

import java.util.ArrayList;
import kmer.AbstractKmerTable;
import shared.Tools;
import stream.Read;
import shared.Parser;

/**
 * Manages the k-mer index for filtering and trimming.
 * Encapsulates the k-mer tables, canonicalization, and mutation logic.
 * Inherits core k-mer configuration constants from BBDukObject.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 18, 2025
 */
public class BBDukIndex extends BBDukObject{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * @param keySets_ The array of k-mer hash tables.
	 * @param refScaffoldNames_ List of reference scaffold names.
	 * @param bbdp BBDuk-specific configuration parameters.
	 * @param parser General parsing object for stream settings.
	 */
	public BBDukIndex(AbstractKmerTable[] keySets_, ArrayList<String> refScaffoldNames_,BBDukParser bbdp,Parser parser){
		super(bbdp,parser);
		keySets=keySets_;
		scaffoldNames=refScaffoldNames_;
		WAYS=keySets.length;
		hdist=bbdp.hammingDistance;
		edist=bbdp.editDistance;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Loading Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Stores the read's kmers, using mutation logic if required. */
	public long addReadKmers(final Read r,final int skip,final int id,final boolean forbidNs,final boolean useShortKmers,final int minSkip_,final int maxSkip_){
		final byte[] bases=r.bases;
		if(bases==null||bases.length<k)return 0;
		long kmer=0;
		long rkmer=0;
		long added=0;
		int len=0;
		final int skip_=Tools.max(minSkip_,Tools.min(maxSkip_,skip)); 

		for(int i=0;i<bases.length;i++){
			final byte b=bases[i];
			final long x=symbolToNumber0[b];
			final long x2=symbolToComplementNumber0[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
			if(forbidNs&&!isFullyDefined(b)){len=0;rkmer=0;}else len++;
			
			if(len>=k){
				final long extraBase=(i>=bases.length-1?-1:symbolToNumber[bases[i+1]]);
				
				if(skip_==1||(len%skip_==0)){
					added+=addToMap(kmer,rkmer,k,extraBase,id,BBDukIndexConstants.kmask,hdist,edist);
					
					if(useShortKmers){
						if(i==k-1)added+=addToMapRightShift(kmer,rkmer,id,hdist,edist);
						if(i==bases.length-1)added+=addToMapLeftShift(kmer,rkmer,extraBase,id,hdist,edist);
					}
				}
			}
		}
		return added;
	}
	
	/** Core logic for kmer insertion, handling no-mutation, hamming, and edit distance. */
	private long addToMap(final long kmer,final long rkmer,final int len,final long extraBase,final int id,final long kmask0,final int hdist_,final int edist_){
		
		final long added;
		if(hdist_==0){
			final long key=BBDukIndexConstants.toValue(kmer,rkmer,kmask0,rcomp,middleMask);
			if(failsSpeed(key))return 0; 
			if(key%WAYS!=tnum)return 0;
			added=keySets[tnum].setIfNotPresent(key,id);
		}else if(edist_>0){
			added=mutate(kmer,rkmer,len,id,edist_,extraBase,lengthMasks[len],hdist_);
		}else{
			added=mutate(kmer,rkmer,len,id,hdist_,-1,lengthMasks[len],hdist_);
		}
		return added;
	}
	
	/** Mutate and store this kmer through 'dist' recursions. */
	private long mutate(final long kmer,final long rkmer,final int len,final int id,final int dist,final long extraBase,final long kmask0,final int maxHdist){
		
		long added=0;
		final long key=BBDukIndexConstants.toValue(kmer,rkmer,kmask0,rcomp,middleMask);
		
		if(key%WAYS==tnum){
			added+=keySets[tnum].setIfNotPresent(key,id);
		}
		
		if(dist>0){
			final int dist2=dist-1;
			
			// Substitution
			for(int j=0;j<symbols;j++){
				for(int i=0;i<len;i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=BBDukIndexConstants.rcomp(temp,len); 
						added+=mutate(temp,rtemp,len,id,dist2,extraBase,kmask0,maxHdist);
					}
				}
			}
			
			// If maxHdist != edist, we are only doing substitution (Hamming), so stop here
			if(maxHdist!=edist)return added; 

			// Deletion
			if(extraBase>=0&&extraBase<=maxSymbol){
				for(int i=1;i<len;i++){
					final long temp=(kmer&leftMasks[i])|((kmer<<bitsPerBase)&rightMasks[i])|extraBase;
					if(temp!=kmer){
						long rtemp=BBDukIndexConstants.rcomp(temp,len); 
						added+=mutate(temp,rtemp,len,id,dist2,-1,kmask0,maxHdist);
					}
				}
			}

			// Insertion
			final long eb2=kmer&symbolMask;
			for(int i=1;i<len;i++){
				final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>bitsPerBase);
				for(int j=0;j<symbols;j++){
					final long temp=temp0|setMasks[j][i-1];
					if(temp!=kmer){
						long rtemp=BBDukIndexConstants.rcomp(temp,len); 
						added+=mutate(temp,rtemp,len,id,dist2,eb2,kmask0,maxHdist);
					}
				}
			}
		}
		return added;
	}
	
	/** Adds short kmers on the left end of the read (mutation logic). */
	private long addToMapLeftShift(long kmer,long rkmer,final long extraBase,final int id,final int hdist_,final int edist_){
		long added=0;
		for(int i=k-1;i>=mink;i--){
			kmer=kmer&rightMasks[i];
			rkmer=rkmer>>>bitsPerBase;
			added+=addToMap(kmer,rkmer,i,extraBase,id,lengthMasks[i],hdist_,edist_);
		}
		return added;
	}
	
	/** Adds short kmers on the right end of the read (mutation logic). */
	private long addToMapRightShift(long kmer,long rkmer,final int id,final int hdist_,final int edist_){
		long added=0;
		for(int i=k-1;i>=mink;i--){
			long extraBase=kmer&symbolMask;
			kmer=kmer>>>bitsPerBase;
			rkmer=rkmer&rightMasks[i];
			added+=addToMap(kmer,rkmer,i,extraBase,id,lengthMasks[i],hdist_,edist_);
		}
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Query Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Queries tables for a kmer, including hamming distance mutations if needed. */
	public final int queryKmerID(final long kmer,final long rkmer,final int len,final int qHDist){
		final long lengthMask=lengthMasks[len];
		int id=queryKmerIDInner(kmer,rkmer,lengthMask,len);
		
		if(id<1&&qHDist>0){
			final int qHDist2=qHDist-1;
			
			//Substitution
			for(int j=0;j<symbols&&id<1;j++){
				for(int i=0;i<len&&id<1;i++){
					final long temp=(kmer&clearMasks[i])|setMasks[j][i];
					if(temp!=kmer){
						long rtemp=BBDukIndexConstants.rcomp(temp,len); 
						id=queryKmerID(temp,rtemp,len,qHDist2);
					}
				}
			}
		}
		return id;
	}
	
	/** Queries tables for a kmer (canonicalizes and looks up). No mutation performed. */
	private final int queryKmerIDInner(final long kmer,final long rkmer,final long lengthMask,final int len){
		final long key=BBDukIndexConstants.toValue(kmer,rkmer,lengthMask,rcomp,middleMask); 
		
		if(passesSpeed(key)){ 
			AbstractKmerTable set=keySets[(int)(key%WAYS)];
			return set.getValue(key);
		}
		return -1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Getter Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public int getWays(){return WAYS;}
	public AbstractKmerTable[] getKeySets(){return keySets;}
	public ArrayList<String> getScaffoldNames(){return scaffoldNames;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final AbstractKmerTable[] keySets;
	private final ArrayList<String> scaffoldNames;
	private final int WAYS;
	private final int hdist;
	private final int edist;
	public int tnum;
}