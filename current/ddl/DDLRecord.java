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
		id=-1;
		name=null;
	}

	public DDLRecord(DynamicDemiLog ddl_, long id_, int taxID_, String name_){
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
		return x!=0 ? x : Long.compare(id, o.id);
	}

	@Override
	public String toString(){
		return "DDLRecord(id="+id+", tid="+taxID+", name="+name+", bases="+bases+")";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final DynamicDemiLog ddl;

	/** First-observed sequence ID — numericID from Streamer, or file index */
	public final long id;
	/** Real taxonomy ID from tid_NNNN or tid|NNNN; -1 when unknown */
	public int taxID=-1;
	public String name;
	public String filename;
	public long bases;
	public int contigs;
	public long cardinality=-1;
	public float gc=-1;
	public String origin;
	public String lineage;

	public byte[] r16S;
	public byte[] r18S;
	public byte[] r5S;
	public byte[] r23S;
	public byte[] r28S;
	public byte[] r5p8S;
	public byte[] rITS;
	public String contigName;
	public int ssuStart=-1;
	public byte ssuStrand=0;

	public boolean hasSSU(){return r16S!=null | r18S!=null;}
	public boolean hasRibo(){return r16S!=null | r18S!=null | r5S!=null | r23S!=null | r28S!=null | r5p8S!=null | rITS!=null;}

	public static final int RIBO_NONE=0, RIBO_16S=1, RIBO_18S=2, RIBO_5S=3,
		RIBO_23S=4, RIBO_28S=5, RIBO_5p8S=6, RIBO_ITS=7;

	public int riboType(){
		if(r16S!=null){return RIBO_16S;}
		if(r18S!=null){return RIBO_18S;}
		if(r5S!=null){return RIBO_5S;}
		if(r23S!=null){return RIBO_23S;}
		if(r28S!=null){return RIBO_28S;}
		if(r5p8S!=null){return RIBO_5p8S;}
		if(rITS!=null){return RIBO_ITS;}
		return RIBO_NONE;
	}

	public byte[] riboBytes(){
		if(r16S!=null){return r16S;}
		if(r18S!=null){return r18S;}
		if(r5S!=null){return r5S;}
		if(r23S!=null){return r23S;}
		if(r28S!=null){return r28S;}
		if(r5p8S!=null){return r5p8S;}
		if(rITS!=null){return rITS;}
		return null;
	}

	public static String riboName(int type){
		switch(type){
			case RIBO_16S: return "16S";
			case RIBO_18S: return "18S";
			case RIBO_5S: return "5S";
			case RIBO_23S: return "23S";
			case RIBO_28S: return "28S";
			case RIBO_5p8S: return "5.8S";
			case RIBO_ITS: return "ITS";
			default: return "-";
		}
	}
}
