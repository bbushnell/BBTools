package hmm;

import java.nio.charset.StandardCharsets;

import shared.Parse;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Precise parser for multi-field delimited Hidden Markov Model (HMM) search result lines.
 * Extracts 23 distinct fields from space-delimited byte arrays using robust parsing with
 * error assertions and type conversion. Handles complex multi-space delimited lines with
 * strict field extraction and supports multiple numeric types (int, float, double).
 *
 * @author Brian Bushnell
 */
public class HMMSearchLine {

	/**
	 * Constructs HMMSearchLine by parsing 23 fields from space-delimited byte array.
	 * Uses sequential parsing with assertion checks to extract protein name, HMM details,
	 * scores, and alignment coordinates from HMM search result line.
	 * @param line_ Raw byte array containing space-delimited HMM search result
	 */
	public HMMSearchLine(byte[] line_){
		line=line_;
		
		Tools.FORCE_JAVA_PARSE_DOUBLE=true; //TODO: put this somewhere else

		int a=0, b=0;
		

//		System.err.println("1: a="+a+", b="+b);
		
		//These are for parsing multi-space-delimited lines
		while(b<line.length && line[b]!=' '){b++;}
//		System.err.println("2: a="+a+", b="+b);
		assert(b>a) : "Missing field 0: "+new String(line);
		name=new String(line, a, b-a, StandardCharsets.US_ASCII);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
//		System.err.println("3: a="+a+", b="+b);
		
		while(b<line.length && line[b]!=' '){b++;}
//		System.err.println("4: a="+a+", b="+b);
		assert(b>a) : "Missing field 1: "+new String(line);
		field1=new String(line, a, b-a, StandardCharsets.US_ASCII);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
//		System.err.println("5: a="+a+", b="+b);
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		length=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		hmmName=new String(line, a, b-a, StandardCharsets.US_ASCII);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		accession=new String(line, a, b-a, StandardCharsets.US_ASCII);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		field5=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		field6=Parse.parseDouble(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		field7=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		field8=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 9: "+new String(line);
		field9=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 10: "+new String(line);
		field10=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 11: "+new String(line);
		field11=Parse.parseDouble(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 12: "+new String(line);
		field12=Parse.parseDouble(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 13: "+new String(line);
		field13=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 14: "+new String(line);
		field14=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 15: "+new String(line);
		field15=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 16: "+new String(line);
		field16=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 17: "+new String(line);
		field17=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 18: "+new String(line);
		field18=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 19: "+new String(line);
		field19=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 20: "+new String(line);
		field20=Parse.parseInt(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 21: "+new String(line);
		field21=Parse.parseFloat(line, a, b);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
		while(b<line.length && line[b]!=' '){b++;}
		assert(b>a) : "Missing field 22: "+new String(line);
		field22=new String(line, a, b-a, StandardCharsets.US_ASCII);
		while(b<line.length && line[b]==' '){b++;}
		a=b;
		
	}
	
	/**
	 * Creates compact tab-delimited text representation of key search result fields.
	 * Includes protein name, sequence length, and HMM model name for summary output.
	 * @return ByteBuilder containing formatted search result summary
	 */
	public ByteBuilder toText(){
		ByteBuilder bb=new ByteBuilder();
		bb.append(name).tab();
		bb.append(length).tab();
		bb.append(hmmName);
		return bb;
	}
	
	/** Returns string representation of the search result using compact text format.
	 * @return String containing tab-delimited summary of key search fields */
	public String toString(){
		return toText().toString();
	}
	
	/** Raw byte array containing the original HMM search result line */
	public final byte[] line;

	//protein_1            -            257 ATP-synt_A           PF00119.18   211   
	//1.9e-49  159.6  27.5   
	//1   1   7.3e-51   2.5e-49  159.2  27.5
	//3   210    41   250    38   251 0.87 -
	
	//field0 - protein_1
	/** Protein sequence identifier (field 0) */
	String name;

	//field1 - -
	/** Search result field 1 (typically dash separator) */
	String field1;

	//field2 - 257
	/** Protein sequence length (field 2) */
	int length; //length?

	//field3 - ATP-synt_A
	/** HMM model name (field 3) */
	String hmmName;

	//field4 - PF00119.18
	/** HMM model accession number (field 4) */
	String accession;
	
	//field5 - 211 
	/** HMM search result field 5 (integer value) */
	int field5;

	//field6 - 1.9e-49
	/** HMM search result field 6 (double precision value) */
	double field6;

	//field7 - 159.6
	/** HMM search result field 7 (floating point value) */
	float field7;

	//field8 - 27.5
	/** HMM search result field 8 (floating point value) */
	float field8;

	//field9 - 1
	/** HMM search result field 9 (floating point value) */
	float field9;

	//field10 - 1
	/** HMM search result field 10 (floating point value) */
	float field10;

	//field11 - 7.3e-51
	/** HMM search result field 11 (double precision value) */
	double field11;

	//field12 - 2.5e-49
	/** HMM search result field 12 (double precision value) */
	double field12;

	//field13 - 159.2
	/** HMM search result field 13 (floating point value) */
	float field13;

	//field14 - 27.5
	/** HMM search result field 14 (floating point value) */
	float field14;

	//field15 - 3
	/** HMM search result field 15 (integer value) */
	int field15;

	//field16 - 210
	/** HMM search result field 16 (integer value) */
	int field16;

	//field17 - 41
	/** HMM search result field 17 (integer value) */
	int field17;
	
	//field18 - 250
	/** HMM search result field 18 (integer value) */
	int field18;

	//field19 - 38
	/** HMM search result field 19 (integer value) */
	int field19;

	//field20 - 251
	/** HMM search result field 20 (integer value) */
	int field20;

	//field21 - 0.87
	/** HMM search result field 21 (floating point value) */
	float field21;

	//field22 - -
	/** HMM search result field 22 (string value, typically dash separator) */
	String field22;
	
	
}
