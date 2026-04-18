package ddl;

import cardinality.DynamicDemiLog;

/**
 * Container for a DynamicDemiLog sketch and its associated metadata.
 * The DDL is the sketch; this holds everything else about the organism.
 *
 * @author Brian Bushnell, Ady
 * @date April 17, 2026
 */
public class DDLRecord implements Comparable<DDLRecord> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public DDLRecord(DynamicDemiLog ddl_){
		ddl=ddl_;
	}

	public DDLRecord(DynamicDemiLog ddl_, int id_, int taxID_, String name_){
		ddl=ddl_;
		id=id_;
		taxID=taxID_;
		name=name_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Natural ordering by taxID, then by id. */
	@Override
	public int compareTo(DDLRecord o){
		int x=taxID-o.taxID;
		return x!=0 ? x : id-o.id;
	}

	@Override
	public String toString(){
		return "DDLRecord(id="+id+", tid="+taxID+", name="+name+", bases="+bases+")";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final DynamicDemiLog ddl;

	/** Load order index — always present, always unique */
	public int id=-1;
	/** Real taxonomy ID from tid_NNNN or tid|NNNN; -1 when unknown */
	public int taxID=-1;
	public String name;
	public String filename;
	public long bases;
	public int contigs;
	public long cardinality=-1;
	public float gc=-1;
}
