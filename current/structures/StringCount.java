package structures;

/**
 * Mutable per-reference hit-count row for BBDuk and Seal statistics.
 * <p>
 * Sorting puts higher base counts first, then higher read counts, then lower
 * lexicographic names.  Seal also stores ambiguous-read counts for output;
 * those do not affect ordering.  The typed equals overload tests compareTo
 * equality but leaves inherited Object.equals as identity equality.
 *
 * @author Brian Bushnell
 */
public class StringCount implements Comparable<StringCount>{

	/** Creates a name-only row; callers fill counters later. */
	public StringCount(String name_){
		name=name_;
	}

	/** Creates a row without ambiguous-read counts. */
	public StringCount(String name_, int len_, long reads_, long bases_){
		this(name_, len_, reads_, bases_, 0);
	}

	/** Creates a complete row including Seal's ambiguous-read count. */
	public StringCount(String name_, int len_, long reads_, long bases_, long ambigReads_){
		name=name_;
		length=len_;
		reads=reads_;
		bases=bases_;
		ambigReads=ambigReads_;
	}

	/** Orders rows by bases descending, reads descending, then name ascending. */
	@Override
	public final int compareTo(StringCount o){
		if(bases!=o.bases){return o.bases>bases ? 1 : -1;}
		if(reads!=o.reads){return o.reads>reads ? 1 : -1;}
		return name.compareTo(o.name);
	}

	/** Tests whether two rows compare as equal. */
	public final boolean equals(StringCount o){
		return compareTo(o)==0;
	}

	/** Returns a name-only hash. */
	@Override
	public final int hashCode(){
		return name.hashCode();
	}

	/** Returns {@code name\tlength\treads\tbases}. */
	@Override
	public final String toString(){
		return name+"\t"+length+"\t"+reads+"\t"+bases;
	}

	/*--------------------------------------------------------------*/

	/** Reference or set name. */
	public final String name;

	/** Reference or set length. */
	public int length;

	/** Read, base and ambiguous-read counters. */
	public long reads, bases, ambigReads;

}
