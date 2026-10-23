package clump2;

import shared.Random;

import shared.Tools;
import stream.Read;

/**
 * Sequence hashes used as duplicate-detection hints, plus exact byte comparisons.
 * Hashing normalizes ASCII letter case; equality does not. Hash matches therefore
 * require an equality check. Read names, qualities and other metadata are ignored.
 * The initialized lookup table is read-only; callers must not mutate read bases
 * or mate links concurrently with these operations.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class Hasher{

	/**
	 * Creates hash code lookup tables for 128 ASCII characters with case normalization.
	 * Automatically maps lowercase letters to their uppercase equivalents for
	 * case-insensitive sequence hashing.
	 * @param modes Number of hash modes to generate per character
	 * @return Hash code table with lowercase letters mapped to uppercase values
	 */
	private static synchronized long[][] makeCodes2(int modes){
		long[][] r=makeCodes(128, modes);
		
		for(int i=0; i<26; i++){
			char c=(char)('A'+i);
			r[Tools.toLowerCase(c)]=r[c];
		}
		return r;
	}
	
	/**
	 * Generates randomized hash code lookup table using seeded random number generation.
	 * Uses fixed seed (1) to ensure deterministic hash values across program runs.
	 * @param symbols Number of ASCII characters to generate codes for
	 * @param modes Number of different hash modes per character
	 * @return 2D array of random long values for hash computation
	 */
	//CLEVER: seed=1 fixes this table for a given Shared.random implementation. Shared.random currently
	//selects FastRandomXoshiro; changing that generator can change hashes despite the unchanged seed.
	private static synchronized long[][] makeCodes(int symbols, int modes){
		Random randy=shared.Shared.random(1);
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong();
			}
		}
		return r;
	}
	
	/**
	 * Computes deterministic hash code for genetic sequence using length-seeded approach.
	 * Each base contributes to hash using mode-specific lookup value with left rotation
	 * for bit mixing. Mode selection cycles based on current hash state.
	 * @param bases Non-null ASCII sequence bytes (0 through 127); not modified
	 * @return Hash code incorporating sequence length and base-specific randomized values
	 */
	//Tabulation hash seeded by length: per base, pick mode=code&31 (in [0,31]; hashcodes has 32 modes),
	//XOR the base's per-mode random long, rotate left 1 to diffuse position. The L65 assert is effectively
	//dead: makeCodes fills ALL 128 entries (no null for any in-range byte; an out-of-range byte b<0/>127
	//AIOOBEs first), so it never catches an "invalid character" - but valid sequences (ACGTN ASCII) never
	//trigger it, and even a stray in-range char just gets a hash (no crash) since this is only a dedup HINT
	//(KmerSort's equalsPaired compares actual bytes and is the arbiter). Harmless.
	public static long hash(byte[] bases){
		long code=bases.length;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int mode=(int)(code&31);
			assert(hashcodes[b]!=null) : "Invalid sequence character: '"+(char)b+"'";
			code=code^hashcodes[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	/**
	 * Computes hash code for a single sequencing read.
	 * @param r Non-null read with non-null ASCII bases
	 * @return Hash code of the read's base sequence
	 */
	public static final long hash(Read r){
		return hash(r.bases);
	}
	
	/**
	 * Computes combined hash code for paired-end reads using XOR and rotation.
	 * For single reads returns individual hash; for pairs combines mate hashes
	 * with left rotation to maintain hash distribution properties.
	 * @param r Non-null primary read; it and any mate require non-null ASCII bases
	 * @return Combined hash for read pair or single read hash if unpaired
	 */
	//a^rotateLeft(b,1) is intentionally ASYMMETRIC (read1 vs read2 have fixed roles), so the pair hash
	//distinguishes (r,mate) order - correct for paired dedup. Matches equalsPaired (primary+mate in order).
	public static final long hashPair(Read r){
		long a=hash(r);
		if(r.mate==null){return a;}
		long b=hash(r.mate);
		return a^Long.rotateLeft(b, 1);
	}
	
	/**
	 * Tests equality for paired reads by comparing both primary and mate sequences.
	 * @param a First read pair
	 * @param b Second read pair
	 * @return true if both primary reads and their mates are sequence-identical
	 * @throws NullPointerException If both primary read references are null
	 */
	public static final boolean equalsPaired(Read a, Read b){
		return equals(a, b) && equals(a.mate, b.mate);
	}
	
	/**
	 * Tests sequence equality between two reads using byte-level comparison.
	 * Two null references are equal; one null reference is unequal. Zero-length
	 * reads compare equal even if one has null base storage and the other an empty
	 * array. Nonempty sequences are compared byte-for-byte, including case.
	 * @param a First read to compare
	 * @param b Second read to compare
	 * @return true if reads have identical base sequences
	 */
	public static final boolean equals(Read a, Read b){
		if(a==b){return true;}
		if(a==null || b==null){
			assert(a!=null || b!=null);
			return false;
		}
		if(a.length()!=b.length()){return false;}
		if(a.length()==0){return true;}
		byte[] ab=a.bases, bb=b.bases;
		for(int i=0; i<ab.length; i++){
			if(ab[i]!=bb[i]){return false;}
		}
		return true;
	}
	
	/**
	 * Precomputed hash lookup table with 32 modes per character for sequence hashing
	 */
	private static final long[][] hashcodes=makeCodes2(32);
	
}
