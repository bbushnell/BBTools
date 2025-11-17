package gff;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;

/**
 * Parses GenBank flat file (.gbff) format and converts it to GFF3 annotation files.
 * Provides sequential file parsing for GenBank files, extracting genomic feature
 * annotations from LOCUS records and transforming them into standardized GFF3 output.
 * Supports command-line usage for direct file conversion with thread-safe processing.
 *
 * @author Brian Bushnell
 */
public class GbffFile {
	
	/**
	 * Program entry point for GenBank to GFF conversion.
	 * Accepts input GenBank file and optional output GFF file path.
	 * Default output is "stdout.gff" if not specified.
	 * @param args Command-line arguments: [input.gbff] [output.gff]
	 */
	public static void main(String[] args){
		String gbff=args[0];
		String gff=(args.length>1 ? args[1] : "stdout.gff");

		if(gbff.indexOf('=')>=0){gbff=gbff.split("=")[1];}
		if(gff.indexOf('=')>=0){gff=gff.split("=")[1];}
		
		FileFormat ffin=FileFormat.testInput(gbff, ".gbff", true);
		FileFormat ffout=FileFormat.testOutput(gff, FileFormat.GFF, null, true, true, false, false);
		GbffFile file=new GbffFile(ffin);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		file.toGff(bsw, true);
		bsw.poisonAndWait();
	}
	
	
//	##gff-version 3
//	#!gff-spec-version 1.21
//	#!processor NCBI annotwriter
//	#!genome-build IMG-taxon 2724679794 annotated assembly
//	#!genome-build-accession NCBI_Assembly:GCF_900182635.1
//	#!annotation-date 07/14/2019 01:52:19
//	#!annotation-source NCBI RefSeq 
//	##sequence-region NZ_FXTD01000001.1 1 528269
//	##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=413815
	
	/**
	 * Converts all GenBank locus records to GFF3 format and writes to stream.
	 * Iterates through all locus records in the file, converting each to GFF format.
	 * Optional header includes GFF version, BBTools version, and column specifications.
	 *
	 * @param bsw Output stream writer for GFF data
	 * @param printHeader Whether to include GFF3 header with version information
	 */
	public void toGff(ByteStreamWriter bsw, boolean printHeader){
		if(printHeader){
			bsw.println("##gff-version 3".getBytes());
			bsw.println(("#BBTools "+Shared.BBTOOLS_VERSION_STRING+" GbffToGff").getBytes());
			bsw.println("#seqid	source	type	start	end	score	strand	phase	attributes".getBytes());
		}
		for(GbffLocus locus=nextLocus(); locus!=null; locus=nextLocus()){
			locus.toGff(bsw);
		}
	}
	
	/**
	 * Constructs a GbffFile parser for the specified GenBank file format.
	 * Validates the input format is GBFF and initializes the file reader.
	 * Automatically resets to beginning of file for parsing.
	 * @param ff_ FileFormat object specifying the GenBank file to parse
	 */
	public GbffFile(FileFormat ff_) {
		ff=ff_;
		assert(ff.format()==FileFormat.GBFF) : ff;
		reset();
	}
	
	/**
	 * Resets the file reader to the beginning of the GenBank file.
	 * Thread-safe method that closes existing reader and creates new ByteFile reader.
	 * Reads first line to prepare for locus parsing, handling empty files gracefully.
	 */
	public synchronized void reset(){
		if(bf!=null){
			bf.close();
			bf=null;
		}
		bf=ByteFile.makeByteFile(ff, FileFormat.GBFF);
		line=bf.nextLine();
		if(line==null){bf.close();}//empty
	}
	
	/**
	 * Parses and returns the next LOCUS record from the GenBank file.
	 * Reads lines from current position until next LOCUS or end of file.
	 * Skips sequence data (ORIGIN sections) and handles multi-line records.
	 * Returns null when no more locus records are available.
	 *
	 * @return GbffLocus object containing parsed locus data, or null if end of file
	 */
	public GbffLocus nextLocus(){
		assert(bf!=null);
		if(line==null){return null;}
		assert(Tools.startsWith(line, "LOCUS ")) : "Expecting: 'LOCUS ...'\nGot: '"+new String(line)+"'";
		ArrayList<byte[]> lines=new ArrayList<byte[]>();
		lines.add(line);
		boolean sequence=false;
		for(line=bf.nextLine(); line!=null && (line.length==0 || line[0]!='L' || !Tools.startsWith(line, "LOCUS ")); line=bf.nextLine()){
			if(line.length>0){
				final byte b=line[0];
				if(b=='/'){
					//skip
				}else if(b=='O' && Tools.startsWith(line, "ORIGIN ")){
					sequence=true;
				}else if(b==' ' && sequence){
						//do nothing
				}else{
					sequence=false;
					lines.add(line);
				}
			}
		}
		if(line==null){bf.close();}
		return new GbffLocus(lines);
	}
	
	/** FileFormat specification for the input GenBank file */
	private final FileFormat ff;
	/** ByteFile reader for processing the GenBank file line-by-line */
	private ByteFile bf;
	/** Current line being processed from the GenBank file */
	private byte[] line=null;
	
}
