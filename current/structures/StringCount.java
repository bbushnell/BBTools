package structures;

/**
 * Object holding a String and numbers, for tracking the number of read and base hits per scaffold.
 * Used by BBDuk and Seal.
 */
public class StringCount implements Comparable<StringCount>{

	/** Creates StringCount with only a name, initializing counts to zero.
	 * @param name_ Scaffold or reference sequence name */
	public StringCount(String name_){
		name=name_;
	}
	/**
	 * Creates StringCount with name, length, and basic counts.
	 * Sets ambiguous reads count to zero.
	 * @param name_ Scaffold or reference sequence name
	 * @param len_ Length of the reference sequence
	 * @param reads_ Number of reads aligned to this scaffold
	 * @param bases_ Number of bases aligned to this scaffold
	 */
	public StringCount(String name_, int len_, long reads_, long bases_){
		this(name_, len_, reads_, bases_, 0);
	}
	/**
	 * Creates StringCount with full statistics including ambiguous read count.
	 * @param name_ Scaffold or reference sequence name
	 * @param len_ Length of the reference sequence
	 * @param reads_ Number of reads aligned to this scaffold
	 * @param bases_ Number of bases aligned to this scaffold
	 * @param ambigReads_ Number of ambiguous reads (multi-mapping or uncertain alignment)
	 */
	public StringCount(String name_, int len_, long reads_, long bases_, long ambigReads_){
		name=name_;
		length=len_;
		reads=reads_;
		bases=bases_;
		ambigReads=ambigReads_;
	}
	@Override
	public final int compareTo(StringCount o){
		if(bases!=o.bases){return o.bases>bases ? 1 : -1;}
		if(reads!=o.reads){return o.reads>reads ? 1 : -1;}
		return name.compareTo(o.name);
	}
	/**
	 * Tests equality based on compareTo result.
	 * @param o Other StringCount to compare
	 * @return true if all comparison fields are equal
	 */
	public final boolean equals(StringCount o){
		return compareTo(o)==0;
	}
	@Override
	public final int hashCode(){
		return name.hashCode();
	}
	@Override
	public final String toString(){
		return name+"\t"+length+"\t"+reads+"\t"+bases;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Scaffold or reference sequence identifier */
	public final String name;
	/** Length of the reference sequence in bases */
	public int length;
	public long reads, bases, ambigReads;
}