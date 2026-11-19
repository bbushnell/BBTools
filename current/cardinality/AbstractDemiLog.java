package cardinality;

import parse.Parser;
import simd.Vector;
import structures.ByteBuilder;

/**
 * Width-aware API for mutable DemiLog sketches; instances are not thread-safe.
 * Typed accessors expose borrowed register arrays and never silently narrow a register.
 * Treat borrowed arrays as read-only: direct edits bypass cached statistics and floors.
 * Existing 16-bit callers retain their char[] access; wide callers use int[].
 * @author Sayu
 */
public abstract class AbstractDemiLog extends CardinalityTracker{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes the shared hashing recipe from the parser's loglog settings. */
	protected AbstractDemiLog(Parser parser){
		super(parser);
	}

	/** Initializes the shared recipe; seed -1 selects CardinalityTracker's default seed. */
	protected AbstractDemiLog(int buckets, int k, long seed, float minProb){
		super(buckets, k, seed, minProb);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Sketch API          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the stored register width in bits: 16 or 32. */
	public abstract int registerBits();
	/** Returns the number of register bits reserved for the exponent. */
	public abstract int exponentWidth();
	/** Returns the minimum register rank used by this implementation's early-exit filter. */
	public abstract int getGlobalNLZ();
	/** Returns borrowed literal-kmer witnesses in bucket order, or null when not retained. */
	public abstract long[] kmerArray();
	/** Returns whether this sketch retains a literal-kmer witness array. */
	public abstract boolean hasKmers();
	/** Appends tab-separated A48 register values to bb without adding a newline; returns bb. */
	public abstract ByteBuilder toBytes(ByteBuilder bb);
	/** Appends retained A48 kmer witnesses to bb; appends nothing without witnesses. */
	public abstract ByteBuilder toBytesKmers(ByteBuilder bb);
	/**
	 * Folds to fewer buckets using the implementation's blacklist preference when requested.
	 * Payload retention and aliasing are implementation-specific; this is not a deep-copy API.
	 */
	public abstract AbstractDemiLog condenseSketch(int newBuckets, boolean blacklist);

	/** Bit pattern as an int; use unsigned comparison for width32. */
	public abstract int registerValue(int bucket);
	/** Returns the effective XOR value used before hashing, with the default seed resolved. */
	public final long seed(){return hashXor;}

	/** Returns borrowed 16-bit registers; unsupported widths throw rather than truncate. */
	public char[] maxArray(){throw new IllegalArgumentException("A 32-bit DDL cannot supply a 16-bit register array");}
	/** Returns borrowed 32-bit registers; unsupported widths throw rather than widen silently. */
	public int[] maxArray32(){throw new IllegalArgumentException("A 16-bit DDL cannot supply a 32-bit register array");}
	/** Returns a detached copy of 16-bit absolute registers; unsupported widths throw. */
	public char[] toAbsoluteArray(){return maxArray().clone();}
	/** Returns this instance as a legacy 16-bit sketch, or throws for another implementation. */
	public final DynamicDemiLog as16(){
		if(!(this instanceof DynamicDemiLog)){throw new IllegalArgumentException("This legacy operation requires a 16-bit DDL");}
		return (DynamicDemiLog)this;
	}

	/**
	 * Requires equal register width, bucket count, k and exponent width before dispatch.
	 * The effective seed must also match for 32-bit sketches; legacy 16-bit seed behavior
	 * is preserved. This does not validate count, witness or blacklist configuration.
	 * @throws IllegalArgumentException if other is null or incompatible
	 */
	public final void requireCompatible(AbstractDemiLog other){
		if(other==null || registerBits()!=other.registerBits() || buckets!=other.buckets || k!=other.k ||
			exponentWidth()!=other.exponentWidth() || (registerBits()==32 && seed()!=other.seed())){
			throw new IllegalArgumentException("DDL comparison/merge requires equal width, buckets, k, exponent and wide-sketch seed");
		}
		assert(buckets==other.buckets) : "SIMD comparison reads corresponding registers without array bounds reconciliation";
	}

	/**
	 * Compares corresponding registers after validating their layout.
	 * @return New array containing lower, equal nonzero, higher and both-empty counts
	 */
	public int[] compareToDetailed(AbstractDemiLog other){
		requireCompatible(other);
		return registerBits()==16 ? Vector.compareDDL(maxArray(), other.maxArray()) :
			Vector.compareDDL(maxArray32(), other.maxArray32());
	}
}
