package sketch;

import java.util.HashMap;

class Record implements Cloneable, Comparable<Record> {
	
	Record(ResultLineParser parser){
		qTaxID=parser.qTaxID;
		rTaxID=parser.rTaxID;
		qBases=parser.qBases;
		rBases=parser.rBases;
		qSize=parser.qSize;
		rSize=parser.rSize;
		ani=parser.ani;
		ssu=parser.ssu;
		taxLevelExtended=parser.taxLevelExtended;
		text=parser.text;
	}
	
	Record copy(){
		try {
			return (Record) this.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	public void processSSU(){
		if(verbose){System.err.println("Record.processSSU(): qTaxID="+qTaxID+", rTaxID="+rTaxID+", ssu="+(float)ssu);}
		final HashMap<Integer, byte[]> map=SSUMap.r16SMap;
		if(ssu<=0 && map!=null){
			byte[] qssu=map.get(qTaxID);
			byte[] rssu=map.get(rTaxID);
			if(qssu!=null && rssu!=null){
				ssu=SketchObject.align(qssu, rssu);
				if(verbose){System.err.println("Aligned; ssu="+(float)ssu);}
			}else{
				if(verbose){System.err.println("Missing: "+qTaxID+"="+qssu+", "+rTaxID+"="+rssu);}
			}
		}else{
			if(verbose){System.err.println("Skipping: ssu="+(float)ssu+", map="+(map==null ? "null" : map.size()));}
		}
	}
	
	@Override
	public int compareTo(Record o) {
		return ani>o.ani ? -1 : ani<o.ani ? 1 : 0;
	}
	
	public double ssu(){return ssu<=0 ? 0 : ssu;}
	
	final int qTaxID;
	final int rTaxID;
	final long qBases;
	final long rBases;
	final long qSize;
	final long rSize;
	final double ani;
	private double ssu=-1;
	int taxLevelExtended;
	int correctNCBI;//0 is false, 1 is true
	int correctSSU;//0 is false, 1 is true
	boolean missingSSU;
	byte[] text;
	//All records are correctNCBI there is no later record with a lower taxLevel.
	//Records are correctSSU if:
	//Their SSUID is >= that of all subsequent records, and the next record (if present) has a SSUID
	
	static boolean verbose;
}