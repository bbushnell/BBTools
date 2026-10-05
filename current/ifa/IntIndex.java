package ifa;

import dna.AminoAcid;
import map.IntListHashMap;
import structures.IntList;

/**
 * Forward, unmasked k-mer position index with a dormant candidate-search API.
 * Every fully defined reference window is indexed; positions are zero-based
 * starts in ascending order. Ambiguous bases break consecutive valid windows.
 * The constructor consumes the reference without retaining or modifying it.
 * This class is not synchronized; its map and position lists are mutable.
 * @author Brian Bushnell
 */
public class IntIndex{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds the forward position map; no middle masking or stride is applied.
	 * @param ref Non-null ASCII reference bases; map capacity is twice its length
	 * @param k_ K-mer length, expected to be 1 through 15 for the int encoding
	 */
	public IntIndex(byte[] ref, int k_){
		k=k_;
		map=new IntListHashMap(ref.length*2);
		indexRef(ref);
	}

	/**
	 * Appends each fully defined window's start to its two-bit encoded key's list.
	 * References shorter than k add no positions. Uses AminoAcid's base table,
	 * including its lowercase and U handling; negative table entries reset the run.
	 * @param ref Reference bases, read only during construction
	 */
	private void indexRef(byte[] ref){
		if(ref.length<k){return;}

		int kmer=0;
		int mask=(1<<(2*k))-1;
		int len=0;

		for(int i=0; i<ref.length; i++){
			byte b=ref[i];
			int x=AminoAcid.baseToNumber[b];

			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				if(++len>=k){
					IntList list=map.get(kmer);
					if(list==null){
						list=new IntList(4);
						map.put(kmer, list);
					}
					list.add(i-k+1);
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Unimplemented candidate-search placeholder; does not inspect either argument.
	 * @param query Intended query bases; currently ignored
	 * @param maxHits Intended hit limit; currently ignored
	 * @return A new empty list on every call, regardless of indexed positions
	 */
	public IntList getCandidates(byte[] query, int maxHits){
		IntList candidates=new IntList();
		// Extract k-mers from query and lookup
		// Could use spaced seeds or multiple k values
		return candidates;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Length of each unmasked forward k-mer. */
	final int k;
	/** Owned mutable map from encoded k-mer to ascending reference start positions. */
	final IntListHashMap map;
}
