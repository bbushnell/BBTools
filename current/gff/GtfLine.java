package gff;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Parses and represents Gene Transfer Format (GTF) annotation lines.
 * Provides simple GTF parsing and conversion to GFF format for genomic feature analysis.
 * Supports tab-delimited GTF files with standard fields including genomic coordinates,
 * strand information, and feature attributes.
 *
 * @author Brian Bushnell
 */
public class GtfLine {
	
	/**
	 * Converts GTF file to GFF3 format.
	 * Reads input GTF file line by line, parses each GTF entry, converts to GFF format,
	 * and writes to output file. Automatically adds GFF3 header and column descriptions.
	 * @param args Command line arguments: [0]=input_file, [1]=output_file (optional)
	 */
	public static void main(String[] args) {
		String in=args[0];
		String out=args.length>1 ? args[1] : "stdout.gff";
		ByteFile bf=ByteFile.makeByteFile(in, true);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(out, true, false, true);
		
		ByteBuilder bb=new ByteBuilder();
		bb.append("##gff-version 3\n");
		bb.append("#seqid	source	type	start	end	score	strand	phase	attributes\n");
		bsw.print(bb);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(Tools.startsWith(line, '#')) {continue;}
			GtfLine gtf=new GtfLine(line);
			GffLine gff=new GffLine(gtf);
			gff.appendTo(bb.clear());
			bsw.print(bb.nl());
		}
		bsw.poisonAndWait();
	}
	
	/**
	 * Constructs a GTF line by parsing tab-delimited input data.
	 * Parses all standard GTF fields including genomic coordinates, strand information,
	 * and feature attributes. Handles missing values by converting '.' to appropriate
	 * default values (-1 for score and frame).
	 *
	 * @param line Tab-delimited byte array containing GTF data
	 */
	public GtfLine(byte[] line){
		
		LineParser1 lp=new LineParser1('\t');
		lp.set(line);
		
	    seqname=lp.parseString(0);
	    source=lp.parseString(1);
	    feature=lp.parseString(2);
	    start=lp.parseInt(3);
	    end=lp.parseInt(4);
	    score=(lp.termEquals('.', 5) ? -1 : lp.parseFloat(5));
	    strand=lp.parseByte(6, 0);
	    frame=(lp.termEquals('.', 7) ? -1 : lp.parseInt(7));
	    attribute=lp.parseString(8);
//		if(start==267921) {System.err.println(lp+" -> "+attribute);}
	}
	
    /** Sequence name/chromosome identifier from GTF record */
    String seqname;
    /** Source program or database that generated the annotation */
    String source;
    /** Feature type (e.g., gene, exon, CDS) */
    String feature;
    /** Start coordinate of the feature (1-based inclusive) */
    int start;
    /** End coordinate of the feature (1-based inclusive) */
    int end;
    /** Feature score or confidence value; -1 if not specified */
    float score;
    /** Strand orientation: '+' for forward, '-' for reverse, '.' for unstranded */
    byte strand;
    /** Reading frame for CDS features (0, 1, 2); -1 if not applicable */
    int frame;
    /** Semicolon-separated attribute field containing feature metadata */
    String attribute;
	
}
