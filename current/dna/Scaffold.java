package dna;

import parse.LineParser1;
import parse.LineParserS1;

/**
 * Mutable scaffold metadata and caller-managed coverage counters.
 * Constructors can parse a tab-delimited SAM @SQ header. Natural ordering and
 * hashing use name, while Object.equals retains identity equality. Callers own
 * synchronization and must not change names while using name-based collections.
 *
 * @author Brian Bushnell
 * @date Jan 4, 2013
 */
public class Scaffold implements Comparable<Scaffold>{

	/** Stores caller-supplied metadata without validation or interning. */
	public Scaffold(String name_, String assembly_, int length_){
		name=name_;
		assembly=assembly_;
		length=length_;
	}

	/** Parses split @SQ fields; ignores tags other than SN, LN and AS.
	 * @deprecated Prefer a reusable LineParser1 or LineParserS1 to splitting each line. */
	@Deprecated
	public Scaffold(String[] split){
		assert(split.length>2 && split[0].equals("@SQ"));
		for(String s : split){
			if(s.equals("@SQ")){
				//Do nothing
			}else if(s.startsWith("SN:")){
				assert(name==null);
				name=new String(s.substring(3)); //Data.forceIntern(s.substring(3));
			}else if(s.startsWith("LN:")){
				length=Integer.parseInt(s.substring(3));
			}else if(s.startsWith("AS:")){
				assembly=Data.forceIntern(s.substring(3));
			}
		}
		assert(length>-1);
		assert(name!=null);
	}

	/**
	 * Parses SN, LN and optional AS from a tab-delimited @SQ line; other tags are ignored.
	 * Interns the parsed name through Data.forceIntern. Mutates the parser's current
	 * field bounds while reading tags. Assertions check required name/nonnegative length
	 * and selected duplicate tags; this is not complete SAM-header validation.
	 * @param lp Borrowed parser positioned on the complete header line
	 */
	public Scaffold(LineParser1 lp){
		assert(lp.startsWith("@SQ"));
		for(int i=1, terms=lp.terms(); i<terms; i++){
			if(lp.termStartsWith("SN:", i)){
				assert(name==null);
				lp.incrementA(3);
				name=lp.parseStringFromCurrentField();
				name=Data.forceIntern(name);
			}else if(lp.termStartsWith("LN:", i)){
				assert(length<=0);
				lp.incrementA(3);
				length=lp.parseIntFromCurrentField();
			}else if(lp.termStartsWith("AS:", i)){
				assert(assembly==null);
				lp.incrementA(3);
				assembly=lp.parseStringFromCurrentField();
			}
		}
		assert(length>-1);
		assert(name!=null);
	}

	/**
	 * String-backed equivalent of the LineParser1 constructor, with the same tag
	 * selection, assertions and current-field mutation.
	 * @param lp Borrowed parser positioned on the complete tab-delimited @SQ line
	 */
	public Scaffold(LineParserS1 lp){
		assert(lp.startsWith("@SQ"));
		for(int i=1; i<lp.terms(); i++){
			if(lp.termStartsWith("SN:", i)){
				assert(name==null);
				lp.incrementA(3);
				name=lp.parseStringFromCurrentField();
				name=Data.forceIntern(name);
			}else if(lp.termStartsWith("LN:", i)){
				assert(length<=0);
				lp.incrementA(3);
				length=lp.parseIntFromCurrentField();
			}else if(lp.termStartsWith("AS:", i)){
				assert(assembly==null);
				lp.incrementA(3);
				assembly=lp.parseStringFromCurrentField();
			}
		}
		assert(length>-1);
		assert(name!=null);
	}

	/** Stores a name and length with no assembly identifier or validation. */
	public Scaffold(String name_, int length_){
		name=name_;
		length=length_;
	}

	/** Hash of the current nonnull name; equality itself remains object identity. */
	@Override
	public int hashCode(){
		return name.hashCode();
	}

	/** Lexicographic name ordering; ignores length, assembly and coverage statistics. */
	@Override
	public int compareTo(Scaffold other){
		return name.compareTo(other.name);
	}

	/** Emits @SQ, SN, LN and optional AS fields without a trailing newline. */
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+(assembly==null ? "" : "\tAS:"+assembly);
	}

	/** Extracts the first SN tag without interning, leaving parser bounds at its value.
	 * Missing SN fails an assertion when enabled, otherwise returns null. */
	public static String name(LineParser1 lp){
		assert(lp.startsWith("@SQ"));
		for(int i=1; i<lp.terms(); i++){
			if(lp.termStartsWith("SN:", i)){
				lp.incrementA(3);
				String name=lp.parseStringFromCurrentField();
				return name;
			}
		}
		assert(false);
		return null;
	}

	/** Reference identifier, also used by ordering and hashing. */
	public String name;
	/** Optional assembly identifier from AS. */
	public String assembly;
	/** Reference length in bases; -1 until supplied or parsed. */
	public int length=-1;
	/** Accumulated mapped-base contributions under the caller's coverage policy. */
	public long basehits=0;
	/** Accumulated accepted read/alignment hits. */
	public long readhits=0;
	/** Fragment-weight accumulator; CoveragePileup divides this by two for FPKM counts. */
	public long fraghits=0;
	/** Accepted read/alignment hits assigned to the minus strand. */
	public long readhitsMinus=0;

	/** Optional caller-owned base-category counts accumulated from aligned read sequences. */
	public long[] basecount;
	/** Reference GC fraction, populated by callers that read the reference. */
	public float gc;
	/** Per-bin GC fraction, computed from the reference by CoveragePileup.processReference(); null if no ref= was given */
	public float[] gcBins;

	/** Caller-owned auxiliary data, such as unstranded or plus-strand coverage storage. */
	public Object obj0;

	/** Second caller-owned auxiliary slot, used for minus-strand coverage when stranded. */
	public Object obj1;

}
