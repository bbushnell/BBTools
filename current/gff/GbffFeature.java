package gff;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Parses a single GenBank (GBFF) feature record and converts it to GFF3.
 * Reassembles wrapped qualifier lines (fixLines), parses the location
 * (start/stop/strand, including complement() and join()), and emits a GFF line
 * via appendGff. Sets error=true on an unparseable location or stop&lt;start.
 *
 * @author Brian Bushnell
 */
public class GbffFeature {

	public GbffFeature(final ArrayList<byte[]> lines0, final String typeString, final String accessionString){
		accession=accessionString;
		setType(typeString);
		parseSlow(lines0);
		if(type==rRNA){
			setSubtype();
		}
		if(stop<start){error=true;}
	}
	
	/** Reassembles wrapped feature lines (fixLines), parses the location from the first reassembled
	 * line (parseStartStop), then scans the rest for product= / locus_tag= / pseudo qualifiers.
	 * @param lines0 The raw fixed-width GBFF lines for this feature. */
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
	
	/** Reassembles wrapped qualifier lines into whole qualifiers: a new '/qualifier' (marker at col 21)
	 * flushes the buffer, and continuation lines are appended to the current one. Returns the rebuilt list
	 * with the location first, then one entry per qualifier. @param lines Raw GBFF feature lines. */
	ArrayList<byte[]> fixLines(ArrayList<byte[]> lines){
		ArrayList<byte[]> fixed=new ArrayList<byte[]>();
		ByteBuilder bb=new ByteBuilder();
		for(byte[] line : lines){
			if(bb.length()>0 && line[21]=='/'){//new '/' qualifier at col 21 -> flush the previous reassembled qualifier (line[21] unguarded; a <22-byte line AIOOBEs = malformed GBFF)
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
	
	/** Appends one fixed-width GBFF line's content to bb: from col 22 for a '/qualifier' line, else from
	 * col 21 (space-separated from prior text). Asserts the col-20 indentation (malformed GBFF crashes loud).
	 * @param bb Accumulator for the current qualifier/location. @param line One raw feature line. */
	void append(ByteBuilder bb, byte[] line){
		//GBFF is fixed-width: location/qualifier text begins at col 21 ('/' marker at 21, value at 22). The asserts enforce that indentation (malformed GBFF -> crash-loud).
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
	
	/** Resolves the feature-type string to its index in typeStrings; asserts it is a known type.
	 * @param typeString GBFF feature key (e.g. "CDS", "gene", "rRNA"). */
	void setType(String typeString){
		int x=Tools.find(typeString, typeStrings);
		assert(x>=0) : x+", "+typeString;
		type=x;
	}
	
	/** Parses the feature location into start/stop/strand: complement(...) sets the minus strand, join(...)
	 * is unwrapped without changing strand, then start=first integer and stop=last integer across the
	 * comma-separated ranges. Sets error=true on any non-numeric character. @param line0 The location text. */
	void parseStartStop(final byte[] line0){
		byte[] line=line0;
		
		if(line[0]=='c'){
			assert(Tools.startsWith(line, "complement("));
			line=Arrays.copyOfRange(line, 11, line.length-1);
			strand=Shared.MINUS;//complement() is the strand operator -> reverse strand
		}
		if(line[0]=='j'){
			assert(Tools.startsWith(line, "join("));
			line=Arrays.copyOfRange(line, 5, line.length-1);
			//NO strand change here: join() is multi-segment, not reverse-strand. complement(join(...)) gets MINUS from the block above; a forward join() must keep the PLUS default. (was: strand=Shared.MINUS -- a copy-paste bug, gff/GbffFeature#001)
		}
		
		//Only complement() and join() are handled. Other INSDC location operators (order(...), bond(...), gap(...))
		//fall through here, hit a non-digit in the loop below, set error=true, and are silently dropped by
		//GbffLocus.parseFeature. QUESTION [gff/GbffFeature path]: likely OK given the tool's stated "features I
		//care about" scope; revisit if order() CDS must survive (then handle order() like join).
		//start = first integer (parse until '.'); the stop loop below resets on each '.'/',' so stop ends as the LAST integer (spans comma-separated join ranges)
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
	
	/** Extracts the value from a key=value qualifier line, stripping the surrounding double-quote
	 * characters. @param line A reassembled qualifier line. @return The unquoted value string. */
	String parseLine(byte[] line){
		String[] split=Tools.equalsPattern.split(new String(line));
		String s=split[1];
		return s.substring(1, s.length()-1);
	}
	
	/** For rRNA features, sets subtype from the first word of the product (e.g. "16S"/"23S") if it
	 * matches a known type; leaves subtype=-1 when product is null or unrecognized. */
	void setSubtype(){
		subtype=-1;
		if(product==null){return;}
		String[] split=Tools.spacePattern.split(product);
		subtype=Tools.find(split[0], typeStrings);
//		assert(false) : type+", "+subtype+", "+split[0]+", "+this.toString()+"\n"+product;
	}
	
	public void toGff(ByteStreamWriter bsw) {
		ByteBuilder bb=bsw.getBuffer();
		appendGff(bb);
		bb.nl();
		bsw.flushBuffer(false);
	}
	
	/** Appends the 9-column GFF3 line for this feature: seqid=accession, source='.', type (pseudogene for a
	 * pseudo gene), start, stop, score='.', strand (Shared.strandCodes2[strand]), phase='.', and attributes
	 * built from product/locus_tag/subtype (or '.' when none). @param bb Output buffer. @return bb, for chaining. */
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
	
	
	/** Returns a string representation in GFF3 format.
	 * @return GFF3 formatted string representation of this feature */
	@Override
	public String toString(){
		return appendGff(new ByteBuilder()).toString();
	}

	public int type=-1;
	public int subtype=-1;
	//TODO: could have coding amino, for tRNA
	public String product;
	public String locus_tag;
//	public String id;
//	public String name;
	
	public int start;
	public int stop;
	public byte strand=Shared.PLUS;
	public String accession;
	public boolean pseudo=false;
	public boolean error=false;

	public static final String[] typeStrings={"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", 
			"5'UTR", "3'UTR", "intron", "exon", "5S", "16S", "23S"};
	public static final String[] typeStringsGff={"gene", "CDS", "rRNA", "tRNA", "ncRNA", "repeat_region", 
			"five_prime_UTR", "three_prime_UTR", "intron", "exon", "5S", "16S", "23S"};
	
	/** Feature-type constants; index into typeStrings/typeStringsGff (e.g. GENE=0 -> "gene", CDS=1 -> "CDS"). */
	public static final int GENE=0, CDS=1, rRNA=2, tRNA=3, ncRNA=4, repeat_region=5, UTR5=6, UTR3=7, intron=8, exon=9;
	/** rRNA subtype constants (5S/16S/23S); also index into typeStrings/typeStringsGff. */
	public static final int r5S=10, r16S=11, r23S=12;
	
}
