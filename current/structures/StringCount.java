package structures;

/**
 * Data structure for tracking scaffold statistics in sequence alignment tools.
 * Stores a scaffold name along with counts of reads and bases that aligned to it.
 * Used primarily by BBDuk and Seal for contamination detection and reference matching.
 * Implements Comparable to enable sorting by base count (descending), then read count,
 * then lexicographically by name.
 *
 * @author Brian Bushnell
 */
public class StringCount implements Comparable<StringCount>{

	public StringCount(String name_){
		name=name_;
	}
	public StringCount(String name_, int len_, long reads_, long bases_){
		this(name_, len_, reads_, bases_, 0);
	}
	public StringCount(String name_, int len_, long reads_, long bases_, long ambigReads_){
		name=name_;
		length=len_;
		reads=reads_;
		bases=bases_;
		ambigReads=ambigReads_;
	}
	/**
	 * Compares StringCount objects for sorting by alignment statistics.
	 * Primary sort: base count (descending)
	 * Secondary sort: read count (descending)
	 * Tertiary sort: name (ascending lexicographic)
	 * @param o Other StringCount to compare against
	 * @return Negative if this < o, positive if this > o, zero if equal
	 */
	@Override
	public final int compareTo(StringCount o){
		if(bases!=o.bases){return o.bases>bases ? 1 : -1;}
		if(reads!=o.reads){return o.reads>reads ? 1 : -1;}
		return name.compareTo(o.name);
	}
	public final boolean equals(StringCount o){
		return compareTo(o)==0;
	}
	/** Hash code based solely on the scaffold name.
	 * @return Hash code of the name string */
	@Override
	public final int hashCode(){
		return name.hashCode();
	}
	/**
	 * Tab-delimited string representation for output formatting.
	 * Format: name\tlength\treads\tbases
	 * @return Tab-separated values string
	 */
	@Override
	public final String toString(){
		return name+"\t"+length+"\t"+reads+"\t"+bases;
	}
	
	/*--------------------------------------------------------------*/
	
	public final String name;
	public int length;
	/** Number of ambiguous reads that map to multiple locations. */
	/** Number of bases aligned to this scaffold. */
	public long reads, bases, ambigReads;
}