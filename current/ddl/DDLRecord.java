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

	public DDLRecord(DynamicDemiLog ddl_, int taxID_, String name_){
		ddl=ddl_;
		taxID=taxID_;
		name=name_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Natural ordering by taxID. */
	@Override
	public int compareTo(DDLRecord o){return taxID-o.taxID;}

	@Override
	public String toString(){
		return "DDLRecord(tid="+taxID+", name="+name+", bases="+bases+")";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final DynamicDemiLog ddl;

	public int taxID=-1;
	public String name;
	public String filename;
	public long bases;
	public int contigs;
	public long cardinality=-1;
	public float gc=-1;
}
