package ddl;

/**
 * Lightweight transport object for SSU sequences found by gene-calling.
 * Carries metadata from GeneCaller/Orf back to the caller without
 * forcing DDLRecord to know about gene-calling internals.
 *
 * @author Noire, Brian Bushnell
 * @date May 19, 2026
 */
public class SSURecord {

	public byte[] bases;
	public int riboType=DDLRecord.RIBO_NONE;
	public String contigName;
	public int start;
	public byte strand;
	public long numericID;
	public String fileName;

	public boolean is16S(){return riboType==DDLRecord.RIBO_16S;}
	public boolean is18S(){return riboType==DDLRecord.RIBO_18S;}

	@Override
	public String toString(){
		return DDLRecord.riboName(riboType)+"("+bases.length+"bp) "+contigName+":"+start+(char)strand;
	}
}
