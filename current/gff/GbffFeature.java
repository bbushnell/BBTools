package gff;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Parses and represents GenBank feature annotations, converting complex feature
 * data into a standardized GFF (Generic Feature Format) representation.
 * Processes raw GenBank feature annotation lines, extracting biological feature
 * details like coordinates, strand information, and qualifiers.
 *
 * @author Brian Bushnell
 */
public class GbffFeature {

	/**
	 * Constructs a GenBank feature from raw annotation lines.
	 * Parses coordinate information, feature type, and qualifiers from GenBank
	 * feature format lines. Performs validation and error checking during parsing.
	 *
	 * @param lines0 Raw GenBank feature annotation lines
	 * @param typeString Feature type identifier (e.g., "gene", "CDS", "rRNA")
	 * @param accessionString Accession number for the sequence containing this feature
	 */
	public GbffFeature(final ArrayList<byte[]> lines0, final String typeString, final String accessionString){
		accession=accessionString;
		setType(typeString);
		parseSlow(lines0);
		if(type==rRNA){
			setSubtype();
		}
		if(stop<start){error=true;}
	}
	
	/**
	 * Main parsing method that processes GenBank feature lines.
	 * Fixes line formatting, extracts coordinates, and parses qualifiers
	 * like product names and locus tags. Sets rRNA subtypes when applicable.
	 * @param lines0 Raw GenBank feature annotation lines to parse
	 */
	private void parseSlow(final ArrayList<byte[]> lines0){
		ArrayList<byte[]> lines=fixLines(lines0);
		parseStartStop(lines.get(0));
		for(int i=1; i<lines.size(); i++){
			byte[] line=lines.get(i);
			if(Tools.startsWith(line, "product=")){
				product=parseLine(line);
			}else if(Tools.startsWith(line, "locus_tag=")){
				locus_tag=parseLine(line);
			}else if(Tools.equals(line, "pseudo")){
				pseudo=true;
			}
			
//			else if(Tools.startsWith(line, "ID=")){
//				id=parseLine(line);
//			}else if(Tools.startsWith(line, "Name=")){
//				name=parseLine(line);
//			}
		}
//		System.err.println("\nvvvvv");
//		for(byte[] line : lines0){
//			System.err.println("'"+new String(line)+"'");
//		}
//		for(byte[] line : lines){
//			System.err.println("'"+new String(line)+"'");
//		}
//		System.err.println("^^^^^");
	}
	
	/**
	 * Combines multi-line GenBank qualifiers into single logical lines.
	 * GenBank format allows qualifiers to span multiple lines with specific
	 * formatting. This method reconstructs the complete qualifier strings.
	 *
	 * @param lines Raw input lines that may be fragmented
	 * @return Fixed lines with complete qualifier statements
	 */
	ArrayList<byte[]> fixLines(ArrayList<byte[]> lines){
		ArrayList<byte[]> fixed=new ArrayList<byte[]>();
		ByteBuilder bb=new ByteBuilder();
		for(byte[] line : lines){
			if(bb.length()>0 && line[21]=='/'){
				fixed.add(bb.toBytes());
				bb.clear();
			}
			append(bb, line);
		}
		if(bb.length()>0){
			fixed.add(bb.toBytes());
			bb.clear();
		}
		return fixed;
	}
	
	/**
	 * Appends GenBank line content to a ByteBuilder, handling line continuations.
	 * Processes GenBank format conventions for multi-line entries, removing
	 * formatting characters and preserving qualifier content.
	 *
	 * @param bb ByteBuilder to append content to
	 * @param line GenBank format line to process
	 */
	void append(ByteBuilder bb, byte[] line){
		assert(line[20]==' ');
		assert(line.length>21);
//		assert(line[21]!=' ') : "'"+new String(line)+"'";
		if(line[21]=='/'){
			bb.append(line, 22, line.length-22);
		}else{
//			System.err.println(line.length+", "+21+", "+(line.length-21+1)+"\n'"+new String(line)+"'");
			if(bb.length>0){bb.append(' ');}
			bb.append(line, 21, line.length-21);
		}
	}
	
	/**
	 * Sets the feature type from a string identifier.
	 * Validates the type against known GenBank feature types and assigns
	 * the corresponding numeric type code.
	 * @param typeString String representation of the feature type
	 */
	void setType(String typeString){
		int x=Tools.find(typeString, typeStrings);
		assert(x>=0) : x+", "+typeString;
		type=x;
	}
	
	/**
	 * Parses coordinate information from GenBank location strings.
	 * Handles complex location formats including complement notation for
	 * reverse strand features and join operations for split features.
	 * Extracts start and stop positions with proper error handling.
	 *
	 * @param line0 GenBank location line containing coordinate information
	 */
	void parseStartStop(final byte[] line0){
		byte[] line=line0;
		
		if(line[0]=='c'){
			assert(Tools.startsWith(line, "complement("));
			line=Arrays.copyOfRange(line, 11, line.length-1);
			strand=Shared.MINUS;
		}
		if(line[0]=='j'){
			assert(Tools.startsWith(line, "join("));
			line=Arrays.copyOfRange(line, 5, line.length-1);
			strand=Shared.MINUS;
		}
		
		int i=0;
		for(start=0; i<line.length; i++){
			int x=line[i];
			if(x=='.'){break;}
			else if(x!='<'){
				if(Tools.isDigit(x)){
					start=start*10+(x-'0');
				}else{
					//if(!error){System.err.println(new String(line0)+"\n"+new String(line));}
					error=true;
				}
			}
		}
//		while(line[i]=='.'){i++;} //Not needed
		for(stop=0; i<line.length; i++){
			int x=line[i];
			if(x=='.' || x==','){
				stop=0;
			}else if(x==' '){
				//do nothing; line wrap
			}else if(x!='>'){
				if(Tools.isDigit(x)){
					stop=stop*10+(x-'0');
				}else{
					//if(!error){System.err.println(new String(line0)+"\n"+new String(line));}
					error=true;
				}
			}
		}
	}
	
	/**
	 * Extracts the value from a GenBank qualifier line.
	 * Parses lines in the format "qualifier=value" and returns the unquoted
	 * value string, removing surrounding quotation marks.
	 *
	 * @param line GenBank qualifier line to parse
	 * @return Extracted qualifier value without quotes
	 */
	String parseLine(byte[] line){
		String[] split=Tools.equalsPattern.split(new String(line));
		String s=split[1];
		return s.substring(1, s.length()-1);
	}
	
	/**
	 * Determines rRNA subtype from the product field.
	 * For rRNA features, examines the product description to identify
	 * specific rRNA subtypes like 5S, 16S, or 23S ribosomal RNA.
	 */
	void setSubtype(){
		subtype=-1;
		if(product==null){return;}
		String[] split=Tools.spacePattern.split(product);
		subtype=Tools.find(split[0], typeStrings);
//		assert(false) : type+", "+subtype+", "+split[0]+", "+this.toString()+"\n"+product;
	}
	
	/**
	 * Writes this feature in GFF3 format to a ByteStreamWriter.
	 * Converts the GenBank feature representation to standard GFF3 format
	 * with proper field formatting and attribute handling.
	 * @param bsw Writer to output the GFF3 formatted line
	 */
	public void toGff(ByteStreamWriter bsw) {
		ByteBuilder bb=bsw.getBuffer();
		appendGff(bb);
		bb.nl();
		bsw.flushBuffer(false);
	}
	
	/**
	 * Appends GFF3 representation of this feature to a ByteBuilder.
	 * Formats all feature information according to GFF3 specification including
	 * coordinates, strand, type, and attributes. Handles pseudogene conversion
	 * and attribute formatting.
	 *
	 * @param bb ByteBuilder to append the GFF3 line to
	 * @return The same ByteBuilder for method chaining
	 */
	public ByteBuilder appendGff(ByteBuilder bb) {
//		bsw.print("#seqid	source	type	start	end	score	strand	phase	attributes\n".getBytes());
		bb.append(accession).tab();
		bb.append('.').tab();
		bb.append((pseudo && type==GENE) ? "pseudogene" : typeStringsGff[type]).tab();
		bb.append(start).tab();
		bb.append(stop).tab();
		bb.append('.').tab();
		bb.append(Shared.strandCodes2[strand]).tab();
		bb.append('.').tab();
		
		boolean attributes=false;
//		if(id!=null){
//			bb.append("ID=").append(id);
//			attributes=true;
//		}
//		if(name!=null){
//			if(attributes){bb.append(';');}
//			bb.append("Name=").append(name);
//			attributes=true;
//		}
		if(product!=null){
			if(attributes){bb.append(';');}
			bb.append("product=").append(product);
			attributes=true;
		}
		if(locus_tag!=null){
			if(attributes){bb.append(';');}
			bb.append("locus_tag=").append(locus_tag);
			attributes=true;
		}
		if(subtype>-1){
			if(attributes){bb.append(';');}
			bb.append("subtype=").append(typeStringsGff[subtype]);
			attributes=true;
		}
		if(!attributes){bb.append('.');}
		return bb;
	}
	
	
	@Override
	public String toString(){
		return appendGff(new ByteBuilder()).toString();
	}

	/** Numeric feature type code (gene, CDS, rRNA, etc.) */
	public int type=-1;
	/** Numeric subtype code for specialized features like rRNA subtypes */
	public int subtype=-1;
	//TODO: could have coding amino, for tRNA
	/** Product description from GenBank product qualifier */
	public String product;
	/** Locus tag identifier from GenBank locus_tag qualifier */
	public String locus_tag;
//	public String id;
//	public String name;
	
	/** Start coordinate of the feature (1-based) */
	public int start;
	/** Stop coordinate of the feature (1-based, inclusive) */
	public int stop;
	/** Strand orientation (Shared.PLUS or Shared.MINUS) */
	public byte strand=Shared.PLUS;
	/** Accession number of the sequence containing this feature */
	public String accession;
	/** True if this feature is marked as a pseudogene */
	public boolean pseudo=false;
	/** True if parsing errors were encountered for this feature */
	public boolean error=false;

	/** String representations of GenBank feature types */
	public static final String[] typeStrings={"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", 
			"5'UTR", "3'UTR", "intron", "exon", "5S", "16S", "23S"};
	/** GFF3-compatible string representations of feature types */
	public static final String[] typeStringsGff={"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", 
			"five_prime_UTR", "three_prime_UTR", "intron", "exon", "5S", "16S", "23S"};
	
	//types
	/** Constant for exon feature type */
	/** Constant for intron feature type */
	/** Constant for 3' untranslated region feature type */
	/** Constant for 5' untranslated region feature type */
	/** Constant for repeat region feature type */
	/** Constant for non-coding RNA feature type */
	/** Constant for transfer RNA feature type */
	/** Constant for ribosomal RNA feature type */
	/** Constant for coding sequence feature type */
	/** Constant for gene feature type */
	public static final int GENE=0, CDS=1, rRNA=2, tRNA=3, ncRNA=4, repeat_region=5, UTR5=6, UTR3=7, intron=8, exon=9;
	//subtypes
	/** Constant for 23S ribosomal RNA subtype */
	/** Constant for 16S ribosomal RNA subtype */
	/** Constant for 5S ribosomal RNA subtype */
	public static final int r5S=10, r16S=11, r23S=12;
	
}
