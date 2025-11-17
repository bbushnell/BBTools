package gff;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import shared.Tools;

/**
 * Parser and processor for GenBank locus records, extracting genomic annotation
 * details from text-based file formats.
 * Systematically parses different sections of a GenBank locus record, including
 * metadata, sequence information, and genomic features. Processes records
 * line-by-line, identifying and extracting specific annotations such as
 * accession, organism, species, and feature details.
 *
 * @author Brian Bushnell
 */
public class GbffLocus {

	/**
	 * Constructs a GbffLocus by parsing the provided GenBank file lines.
	 * Processes all lines sequentially, parsing each block type until the
	 * entire locus record is processed.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 */
	public GbffLocus(ArrayList<byte[]> lines) {
		while(num<lines.size()){
			parseBlock(lines);
		}
	}

	/**
	 * Parses a single block from GenBank format based on the block type identifier.
	 * Dispatches to appropriate parsing methods for LOCUS, DEFINITION, ACCESSION,
	 * VERSION, DBLINK, KEYWORDS, SOURCE, REFERENCE, COMMENT, FEATURES, CONTIG,
	 * ORIGIN, and PRIMARY blocks.
	 *
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Current line number after processing the block
	 */
	int parseBlock(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(Tools.startsWith(line, " ")){
			assert(false) : line;
			num++;
		}else if(Tools.startsWith(line, "LOCUS ")){
			parseLocus(lines);
		}else if(Tools.startsWith(line, "DEFINITION ")){
			parseDefinition(lines);
		}else if(Tools.startsWith(line, "ACCESSION ")){
			parseAccession(lines);
		}else if(Tools.startsWith(line, "VERSION ")){
			parseVersion(lines);
		}else if(Tools.startsWith(line, "DBLINK ")){
			parseDBLink(lines);
		}else if(Tools.startsWith(line, "KEYWORDS ")){
			parseKeywords(lines);
		}else if(Tools.startsWith(line, "SOURCE ")){
			parseSource(lines);
		}else if(Tools.startsWith(line, "REFERENCE ")){
			parseReference(lines);
		}else if(Tools.startsWith(line, "COMMENT ")){
			parseComment(lines);
		}else if(Tools.startsWith(line, "FEATURES ")){
			parseFeatures(lines);
		}else if(Tools.startsWith(line, "CONTIG ")){
			parseContig(lines);
		}else if(Tools.startsWith(line, "ORIGIN ")){
			parseOrigin(lines);
		}else if(Tools.startsWith(line, "PRIMARY ")){
			parsePrimary(lines);
		}else{
			assert(false) : "Unhandled block type: "+new String(line);
		}
		return num;
	}
	
	/**
	 * Advances to the next non-empty line in the GenBank file.
	 * Skips empty lines and null entries while maintaining bounds checking.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Next non-empty line as byte array, or null if at end
	 */
	private byte[] nextLine(ArrayList<byte[]> lines){
		byte[] line=null;
		for(final int lim=lines.size()-1; num<lim && (line==null || line.length==0); ){
//			System.err.println(num+", "+lim);
			num++;
			line=lines.get(num);
		}
//		System.err.println(line);
//		assert(line!=null);
		return line;
	}
	
	/**
	 * Gets the current line without advancing the line pointer.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Current line as byte array, or null if beyond bounds
	 */
	private byte[] getLine(ArrayList<byte[]> lines){
		return num>=lines.size() ? null : lines.get(num);
	}
	
	/** Move pointer to next block start */
	private int advanceBlock(ArrayList<byte[]> lines){
		for(num++; num<lines.size(); num++){
			byte[] line=lines.get(num);
			if(line!=null && line.length>0 && line[0]!=' '){break;}
		}
		return num;
	}
	
	/** Move pointer to next block start */
	private int advanceFeature(ArrayList<byte[]> lines){
		for(num++; num<lines.size(); num++){
			byte[] line=lines.get(num);
			if(line!=null && line.length>0 && (line[0]!=' ' || line[5]!=' ')){break;}
		}
		return num;
	}
	
	/**
	 * Extracts the content portion of a GenBank block line by removing the
	 * block identifier and formatting spaces. Assumes standard GenBank format
	 * with block name ending at position 11 followed by a space.
	 *
	 * @param line GenBank block line as byte array
	 * @return Content string without block identifier prefix
	 */
	private String trimBlockName(byte[] line){
		assert(line.length>=12 && line[11]==' ') : new String(line);
		return new String(line, 12, line.length-12);
	}
	
	/**
	 * Extracts the feature type from a GenBank FEATURES section line.
	 * Parses the feature type located between positions 5-20 in the standard
	 * GenBank format.
	 *
	 * @param line Feature line as byte array from FEATURES section
	 * @return Feature type string (e.g., "CDS", "gene", "tRNA")
	 */
	private String toFeatureType(byte[] line){
		assert(line[4]==' ');
		assert(line[5]!=' ');
		assert(line[20]==' ');
		int start=5, stop=6;
		for(; stop<21 && line[stop]!=' '; stop++){}
		return new String(line, start, stop-start);
	}
	
	/**
	 * Parses the LOCUS block to extract the primary accession identifier.
	 * If no accession has been set, uses the first whitespace-separated token
	 * from the LOCUS line as the accession.
	 *
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseLocus(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
//		assert(Tools.startsWith(line, "LOCUS")) : new String(line);
		if(accession==null){
			String s=trimBlockName(line);
			String[] split=Tools.whitespacePlus.split(s);
			accession=split.length>0 ? split[0] : null;
		}
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the DEFINITION block to extract organism information.
	 * Uses the first comma-separated token from the definition line as the
	 * organism name if not already set.
	 *
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseDefinition(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(organism==null){
			String s=trimBlockName(line);
			String[] split=Tools.commaPattern.split(s);
			organism=split.length>0 ? split[0] : null;
		}
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the ACCESSION block to extract the primary accession number.
	 * Uses the first whitespace-separated token as the accession if not
	 * already set.
	 *
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseAccession(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(accession==null){
			String s=trimBlockName(line);
			String[] split=Tools.whitespacePlus.split(s);
			accession=split.length>0 ? split[0] : null;
		}
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the VERSION block to extract versioned accession information.
	 * Updates the accession with version information if the version is longer
	 * than the current accession or if no accession exists.
	 *
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseVersion(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		String s=trimBlockName(line);
		String[] split=Tools.whitespacePlus.split(s);
		s=split.length>0 ? split[0] : null;
		if(accession==null || (s!=null && s.length()>1)){
			accession=s;
		}
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the DBLINK block containing database cross-references.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseDBLink(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the KEYWORDS block containing classification keywords.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseKeywords(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the SOURCE block to extract species information.
	 * Sets the species field from the source line content if not already set.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseSource(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(species==null){
			species=trimBlockName(line);
		}
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the REFERENCE block containing publication citations.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseReference(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the COMMENT block containing additional annotations.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseComment(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the FEATURES block to extract genomic feature annotations.
	 * Processes each feature line, identifying feature types and creating
	 * GbffFeature objects for supported types. Only processes features
	 * matching the predefined feature types array.
	 *
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after processing all features
	 */
	private int parseFeatures(ArrayList<byte[]> lines){
		for(byte[] line=nextLine(lines); line!=null && line[0]==' '; line=getLine(lines)){
//			System.err.println(num+": "+new String(line));
			String type=toFeatureType(line);
			int idx=Tools.find(type, featureTypes);
//			System.err.println("idx="+idx+" for '"+type+"'");
			if(idx>=0){
//				System.err.println("parseFeature");
				parseFeature(lines, type);
//				System.err.println(features.get(features.size()-1));
			}else{
//				System.err.println("advanceFeature");
				advanceFeature(lines);
			}
		}
		return num;
	}
	
	/** Move pointer to next block start */
	private int parseFeature(ArrayList<byte[]> lines, String type){
		ArrayList<byte[]> flist=new ArrayList<byte[]>();
		flist.add(lines.get(num));
		for(num++; num<lines.size(); num++){
			byte[] line=lines.get(num);
			if(line!=null && line.length>0 && (line[0]!=' ' || line[5]!=' ')){
//				assert(false) : Character.toString(line[0])+", "+Character.toString(line[5])+", "+Character.toString(line[6])+"\n"+new String(line);
				break;
			}
			flist.add(line);
		}
		GbffFeature f=new GbffFeature(flist, type, accession);
		if(!f.error){
			features.add(f);
		}else{
//			System.err.println("Failed to parse feature "+f);
		}
		return num;
	}
	
	/**
	 * Parses the CONTIG block containing assembly information.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseContig(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the ORIGIN block containing sequence data.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parseOrigin(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Parses the PRIMARY block containing primary sequence information.
	 * Currently skips processing and advances to the next block.
	 * @param lines ArrayList of byte arrays containing GenBank file lines
	 * @return Line number after advancing past this block
	 */
	private int parsePrimary(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	/**
	 * Converts the parsed locus information to GFF3 format and writes to output.
	 * Outputs sequence region information and converts all valid features
	 * (CDS, tRNA, rRNA) that are not pseudo-genes to GFF format.
	 * @param bsw ByteStreamWriter for output writing
	 */
	public void toGff(ByteStreamWriter bsw) {
		final byte[] accessionB=accession.getBytes();
		bsw.print(seqRegB);
		bsw.print(accessionB);
		if(start>0 && stop>0){
			bsw.print(' ').print(start).print(' ').print(stop);
		}
		bsw.println();
		for(GbffFeature f : features){
			if(f.type==GbffFeature.CDS || f.type==GbffFeature.tRNA || f.type==GbffFeature.rRNA){
				if(!f.pseudo && !f.error){
					f.toGff(bsw);
				}
			}
		}
	}
	
	
	/** Line number */
	int num=0;
	
	/** Flag to control whether gene features are printed in output */
	boolean printGene=false;
	/** Flag to control whether repeat features are printed in output */
	boolean printRepeat=false; 
	
	/** Array of supported GenBank feature types for parsing */
	public static String[] featureTypes=GbffFeature.typeStrings;
	/** Byte array constant for GFF sequence region header */
	private static final byte[] seqRegB="##sequence-region ".getBytes();
	
	/** GenBank accession number for this locus */
	String accession;
	/** Organism name extracted from DEFINITION block */
	String organism;
	/** Species name extracted from SOURCE block */
	String species;
	/** Start coordinate for sequence region */
	int start;
	/** Stop coordinate for sequence region */
	int stop;
	/** List of parsed genomic features from this locus */
	ArrayList<GbffFeature> features=new ArrayList<GbffFeature>();
}
