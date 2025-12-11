package dna;

import shared.LineParser1;
import shared.LineParserS1;

/**
 * Represents a genomic scaffold with metadata for sequence assembly and tracking.
 * Manages genomic scaffold information, parsing SAM (Sequence Alignment/Map) format
 * entries and providing essential metadata about genomic sequences.
 *
 * @author Brian Bushnell
 * @date Jan 4, 2013
 */
public class Scaffold implements Comparable<Scaffold> {
	
	public Scaffold(String name_, String assembly_, int length_){
		name=name_;
		assembly=assembly_;
		length=length_;
	}
	
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
	
	/**
	 * Parses scaffold from SAM header format using LineParser1 for improved performance.
	 * Assumes SAM format: @SQ SN:scaffold_0 LN:1785514 AS:build_9
	 * @param lp LineParser1 positioned at a @SQ header line
	 */
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
	
	/**
	 * Parses scaffold from SAM header format using LineParserS1 for improved performance.
	 * Assumes SAM format: @SQ SN:scaffold_0 LN:1785514 AS:build_9
	 * @param lp LineParserS1 positioned at a @SQ header line
	 */
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
	
	public String name;
	public String assembly;
	public int length=-1;
	public long basehits=0;
	public long readhits=0;
	/** Number of fragments aligned to this scaffold for FPKM calculation */
	public long fraghits=0;
	public long readhitsMinus=0;
	
	public long[] basecount;
	public float gc;
	
	public Object obj0;
	
	public Object obj1;
	
}
