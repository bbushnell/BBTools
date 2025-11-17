package dna;

import shared.LineParser1;
import shared.LineParserS1;

/**
 * @author Brian Bushnell
 * @date Jan 4, 2013
 *
 */
public class Scaffold implements Comparable<Scaffold> {
	
	/**
	 * Creates a scaffold with complete metadata.
	 * @param name_ Scaffold name/identifier
	 * @param assembly_ Assembly build information
	 * @param length_ Total scaffold sequence length
	 */
	public Scaffold(String name_, String assembly_, int length_){
		name=name_;
		assembly=assembly_;
		length=length_;
	}
	
	/** Assumes SAM format 
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	@Deprecated
	public Scaffold(String[] split) {
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
	
	/** Should be faster. Assumes SAM format.
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	public Scaffold(LineParser1 lp) {
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
	
	/** Should be faster. Assumes SAM format.
	 * e.g.<br> @SQ	SN:scaffold_0	LN:1785514	AS:build 9 */
	public Scaffold(LineParserS1 lp) {
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
	
	/**
	 * Creates a scaffold with name and length only.
	 * Assembly information will be null.
	 * @param name_ Scaffold name/identifier
	 * @param length_ Total scaffold sequence length
	 */
	public Scaffold(String name_, int length_) {
		name=name_;
		length=length_;
	}
	
	@Override
	public int hashCode(){
		return name.hashCode();
	}
	
	@Override
	public int compareTo(Scaffold other){
		return name.compareTo(other.name);
	}
	
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+(assembly==null ? "" : "\tAS:"+assembly);
	}
	
	/**
	 * Extracts scaffold name from SAM header line without creating full Scaffold object.
	 * @param lp LineParser1 positioned at a @SQ header line
	 * @return Scaffold name from SN field
	 */
	public static String name(LineParser1 lp) {
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
	
	/** Scaffold name/identifier */
	public String name;
	/** Assembly build information */
	public String assembly;
	/** Total scaffold sequence length in bases */
	public int length=-1;
	/** Number of bases that received alignment hits */
	public long basehits=0;
	/** Total number of reads aligned to this scaffold */
	public long readhits=0;
	/** For calculating FPKM */
	public long fraghits=0;
	/** Number of reads aligned to the minus strand of this scaffold */
	public long readhitsMinus=0;
	
	/** {A,C,G,T,N} */
	public long[] basecount;
	/** GC content as a fraction (0.0 to 1.0) */
	public float gc;
	
	/** For attaching things */
	public Object obj0;
	
	/** For attaching things for strand1 */
	public Object obj1;
	
}
