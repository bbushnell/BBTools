package prot;

import java.util.Map;

/** Immutable original-member bounds, independent of consensus and query lengths.
 * An absent option pair preserves the historical100..200 contract. Short packages
 * must supply both bounds; no process-global setting relaxes other readers.
 * @author Keqing
 */
public final class MemberLengthRange{

	public MemberLengthRange(int minimum, int maximum){
		if(minimum<40 || maximum>200 || minimum>maximum){
			throw new IllegalArgumentException("Member range requires40<=minmemberaa<=maxmemberaa<=200: "+minimum+".."+maximum);
		}
		min=minimum; max=maximum;
	}

	/** Parse a complete pair without modifying the caller's recipe. */
	public static MemberLengthRange from(Map<String,String> options){
		if(options==null){throw new IllegalArgumentException("Member range requires an options map");}
		return from(options.get(MIN_KEY), options.get(MAX_KEY));
	}

	public static MemberLengthRange from(String minimum, String maximum){
		if(minimum==null && maximum==null){return LEGACY;}
		if(minimum==null || maximum==null){throw new IllegalArgumentException("minmemberaa and maxmemberaa are required together");}
		return new MemberLengthRange(Integer.parseInt(minimum), Integer.parseInt(maximum));
	}

	public boolean contains(int length){return length>=min && length<=max;}

	/** Tests the aggregate before allocation; overflow is an invalid census, never a wrapped bound. */
	public boolean containsTotal(long members, long residues){
		if(members<0 || residues<0){return false;}
		try{return residues>=Math.multiplyExact(members, (long)min) && residues<=Math.multiplyExact(members, (long)max);}
		catch(ArithmeticException overflow){return false;}
	}

	public boolean same(MemberLengthRange other){return other!=null && min==other.min && max==other.max;}

	/** Existing recipes remain byte-identical when they use the historical range. */
	public String recipe(String separator){
		if(separator==null){throw new IllegalArgumentException("Range metadata needs an explicit key/value separator");}
		return same(LEGACY) ? "" : MIN_KEY+separator+min+"\n"+MAX_KEY+separator+max+"\n";
	}

	/** Adds only nonlegacy provenance; these keys describe members, not query filtering. */
	public void put(Map<String,String> options){
		if(options==null){throw new IllegalArgumentException("Member range requires a recipe destination");}
		if(!same(LEGACY)){options.put(MIN_KEY, Integer.toString(min)); options.put(MAX_KEY, Integer.toString(max));}
	}

	@Override public String toString(){return min+".."+max;}
	public final int min, max;
	public static final String MIN_KEY="minmemberaa", MAX_KEY="maxmemberaa";
	public static final MemberLengthRange LEGACY=new MemberLengthRange(100, 200);
}
