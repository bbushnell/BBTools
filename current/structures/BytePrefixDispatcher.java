package structures;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/**
 * Immutable, ZERO-ALLOCATION prefix dispatcher: maps a small fixed set of byte-string keys to int ids
 * (0..N-1, input order) by matching a key that PREFIXES a query byte[] and is terminated by a delimiter.
 * Built for the common parser pattern "a labeled line 'name&lt;delim&gt;value...' -&gt; route by name"
 * without materializing a String per line (which an ObjectIntMap&lt;String&gt;/HashMap lookup would).
 *
 * <p>Layout: a first-byte-bucketed flat table, chosen over a trie/double-array-trie/hash for SMALL key
 * sets (~15-30 short keys). The properties that hold unconditionally: ZERO query allocation, and no
 * per-byte dependent pointer chase (a trie serializes one dependent node-load per byte, each blocking on
 * the previous). The DESIGN HYPOTHESIS -- to be confirmed by benchmark, not asserted -- is that a
 * first-byte dispatch narrowing to a tiny contiguous candidate range (usually 1-3), then a single
 * SEQUENTIAL byte compare the CPU can prefetch (compare bytes in cache alongside the query line), incurs
 * fewer serialized memory hops than a trie/DAT, and beats a hash's scan-then-verify. The latency margin is
 * exactly what the benchmark measures; the structure is right even if a given JVM makes the margin small.
 * <ul>
 * <li>{@code bucketStart[257]} -- CSR row starts; candidates for first byte b are
 *     {@code [bucketStart[b], bucketStart[b+1])}.</li>
 * <li>{@code keyData[]} -- all key bytes concatenated, in bucket-sorted order.</li>
 * <li>{@code keyOffset[]}, {@code keyLength[]} -- the slice of keyData for each candidate slot.</li>
 * <li>{@code id[]} -- the input-order id (0..N-1) of the key at each slot.</li>
 * </ul>
 *
 * <p>Case-sensitive. ASCII keys only. Immutable after construction.
 *
 * @author Noire (design with Barbara)
 * @date July 30, 2026
 */
public final class BytePrefixDispatcher {

	/** The delimiter byte that must terminate a matched key in the query. */
	private final byte delimiter;
	/** Number of keys. */
	private final int n;
	/** CSR bucket starts by first byte (length 257); firstByte b -> slots [bucketStart[b], bucketStart[b+1]). */
	private final int[] bucketStart;
	/** All key bytes concatenated, in bucket-sorted slot order. */
	private final byte[] keyData;
	/** Per-slot start into keyData. */
	private final int[] keyOffset;
	/** Per-slot key length. */
	private final int[] keyLength;
	/** Per-slot input-order id. */
	private final int[] id;
	/** Original keys by id, for diagnostics only. */
	private final String[] keysById;

	/**
	 * Builds the dispatcher from keys (id = input-order index) and the terminating delimiter.
	 * @param keys Distinct, non-empty, ASCII keys not containing the delimiter.
	 * @param delimiter The byte that must immediately follow a matched key in a query (e.g. '\t').
	 */
	public BytePrefixDispatcher(String[] keys, byte delimiter){
		assert(keys!=null && keys.length>0) : "No keys provided";
		this.delimiter=delimiter;
		this.n=keys.length;
		this.keysById=keys.clone();

		final byte[][] kb=new byte[n][];
		for(int i=0; i<n; i++){
			final String k=keys[i];
			assert(k!=null && k.length()>0) : "Null/empty key at index "+i;
			for(int c=0; c<k.length(); c++){
				assert(k.charAt(c)<128) : "Non-ASCII key '"+k+"' (US-ASCII would substitute '?')";
				assert((byte)k.charAt(c)!=delimiter) : "Key '"+k+"' contains the delimiter byte "+delimiter;
			}
			kb[i]=k.getBytes(StandardCharsets.US_ASCII);
		}
		for(int i=0; i<n; i++){
			for(int j=i+1; j<n; j++){
				assert(!Arrays.equals(kb[i], kb[j])) : "Duplicate key '"+keys[i]+"'";
			}
		}

		//CSR bucket starts by first byte.
		final int[] cnt=new int[256];
		for(int i=0; i<n; i++){cnt[kb[i][0]&0xFF]++;}
		bucketStart=new int[257];
		for(int b=0; b<256; b++){bucketStart[b+1]=bucketStart[b]+cnt[b];}

		//Assign each key to a slot within its first-byte bucket.
		id=new int[n];
		keyOffset=new int[n];
		keyLength=new int[n];
		final int[] cursor=new int[256];
		final int[] slotOfKey=new int[n];
		for(int i=0; i<n; i++){
			final int fb=kb[i][0]&0xFF;
			final int slot=bucketStart[fb]+(cursor[fb]++);
			slotOfKey[i]=slot;
			id[slot]=i;
			keyLength[slot]=kb[i].length;
		}

		int total=0;
		for(byte[] b : kb){total+=b.length;}
		keyData=new byte[total];
		int pos=0;
		for(int slot=0; slot<n; slot++){
			final byte[] b=kb[id[slot]];
			keyOffset[slot]=pos;
			System.arraycopy(b, 0, keyData, pos, b.length);
			pos+=b.length;
		}
	}

	/**
	 * Resolves the key that prefixes {@code q} starting at {@code start} and is terminated by the delimiter
	 * (which must fall strictly before {@code limit}). ZERO allocation.
	 * @param q Query bytes (e.g. a line).
	 * @param start Index where the prefix begins.
	 * @param limit Exclusive end of valid bytes in q (usually the line length).
	 * @return The key's input-order id (0..N-1), or -1 if no key matches (unknown name, missing/mispositioned
	 *         delimiter, truncation, or empty input).
	 */
	public int lookup(byte[] q, int start, int limit){
		if(q==null || start>=limit || start<0){return -1;}
		final int fb=q[start]&0xFF;
		final int hi=bucketStart[fb+1];
		for(int s=bucketStart[fb]; s<hi; s++){
			final int len=keyLength[s];
			final int end=start+len;
			//Delimiter must sit exactly one past the key: this rejects a key that is a PREFIX of a longer
			//token (base vs bases) and a key that runs past the query.
			if(end<limit && q[end]==delimiter){
				final int off=keyOffset[s];
				int j=1;//byte 0 already matched via the first-byte bucket
				while(j<len && q[start+j]==keyData[off+j]){j++;}
				if(j==len){return id[s];}
			}
		}
		return -1;
	}

	/** Convenience: lookup over the whole array. */
	public int lookup(byte[] q){return lookup(q, 0, q==null ? 0 : q.length);}

	/** Number of keys. */
	public int size(){return n;}

	/** The original key for an id (diagnostics only). */
	public String keyForId(int keyId){return keysById[keyId];}

	/*--------------------------------------------------------------*/
	/*----------------          Self-test           ----------------*/
	/*--------------------------------------------------------------*/

	/** Correctness self-test over the clade key set + the nasty inputs Barbara flagged. */
	public static void main(String[] args){
		final String[] keys={"tid","level","name","lineage","domain","gc","entropy","strandedness",
				"bases","contigs","1mers","2mers","3mers","4mers","5mers","ddl","16S","18S"};
		final BytePrefixDispatcher d=new BytePrefixDispatcher(keys, (byte)'\t');
		int fails=0;
		fails+=check(d, "tid\t123", 0);
		fails+=check(d, "level\t2\tspecies", 1);
		fails+=check(d, "name\tEscherichia coli", 2);
		fails+=check(d, "names\tfoo", -1);      //name is a prefix of names -> must NOT match name
		fails+=check(d, "lineage\td__Bacteria", 3);
		fails+=check(d, "level", -1);           //no delimiter
		fails+=check(d, "domain\t0", 4);
		fails+=check(d, "1mers\t1 2 3 4", 10);
		fails+=check(d, "16S\tACGT", 16);
		fails+=check(d, "18S\tACGT", 17);
		fails+=check(d, "potato\tx", -1);       //unknown key
		fails+=check(d, "", -1);                //empty
		fails+=check(d, "\t", -1);              //delimiter first
		fails+=check(d, "bases\t99", 8);
		fails+=check(d, "base\t99", -1);        //base is not a key (bases is) -> no match
		//nonzero start (substring lookup)
		byte[] q="XXtid\t7".getBytes(StandardCharsets.US_ASCII);
		if(d.lookup(q, 2, q.length)!=0){System.err.println("FAIL nonzero-start"); fails++;}
		//limit < a.length that still contains the delimiter -> match
		byte[] q2="tid\t7GARBAGE".getBytes(StandardCharsets.US_ASCII);
		if(d.lookup(q2, 0, 5)!=0){System.err.println("FAIL limit<length (delim inside)"); fails++;}
		//limit cutting the delimiter off (delimiter at index 3, limit 3 -> excluded) -> no match
		if(d.lookup(q2, 0, 3)!=-1){System.err.println("FAIL limit excludes delim"); fails++;}
		//exact key at end with no delimiter -> no match
		byte[] q3="domain".getBytes(StandardCharsets.US_ASCII);
		if(d.lookup(q3, 0, q3.length)!=-1){System.err.println("FAIL exact-key-no-delim"); fails++;}
		//unknown high-bit first byte -> empty bucket -> no match (no AIOOBE)
		byte[] q4={(byte)0x80, (byte)'\t'};
		if(d.lookup(q4, 0, q4.length)!=-1){System.err.println("FAIL high-bit first byte"); fails++;}
		//duplicate input keys must be rejected at build (assertions on)
		boolean rejectedDup=false;
		try{new BytePrefixDispatcher(new String[]{"a","b","a"}, (byte)'\t');}
		catch(AssertionError e){rejectedDup=true;}
		if(!rejectedDup){System.err.println("FAIL duplicate-key not rejected"); fails++;}
		System.err.println(fails==0 ? "ALL SELF-TESTS PASS ("+keys.length+" keys)" : fails+" SELF-TEST FAILURES");
		if(fails>0){System.exit(1);}
	}

	private static int check(BytePrefixDispatcher d, String line, int expect){
		byte[] q=line.getBytes(StandardCharsets.US_ASCII);
		int got=d.lookup(q, 0, q.length);
		if(got!=expect){System.err.println("FAIL: '"+line+"' -> "+got+" expected "+expect); return 1;}
		return 0;
	}
}
