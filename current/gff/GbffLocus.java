package gff;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import shared.Tools;

public class GbffLocus {

	public GbffLocus(ArrayList<byte[]> lines) {
		while(num<lines.size()){
			parseBlock(lines);
		}
	}

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
	
	private String trimBlockName(byte[] line){
		assert(line.length>=12 && line[11]==' ') : new String(line);
		return new String(line, 12, line.length-12);
	}
	
	private String toFeatureType(byte[] line){
		assert(line[4]==' ');
		assert(line[5]!=' ');
		assert(line[20]==' ');
		int start=5, stop=6;
		for(; stop<21 && line[stop]!=' '; stop++){}
		return new String(line, start, stop-start);
	}
	
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
	
	private int parseDefinition(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(organism==null){
			String s=trimBlockName(line);
			String[] split=Tools.commaPattern.split(s);
			organism=split.length>0 ? split[0] : null;
		}
		return advanceBlock(lines);
	}
	
	private int parseAccession(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(accession==null){
			String s=trimBlockName(line);
			String[] split=Tools.whitespacePlus.split(s);
			accession=split.length>0 ? split[0] : null;
		}
		return advanceBlock(lines);
	}
	
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
	
	private int parseDBLink(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	private int parseKeywords(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	private int parseSource(ArrayList<byte[]> lines){
		byte[] line=lines.get(num);
		if(species==null){
			species=trimBlockName(line);
		}
		return advanceBlock(lines);
	}
	
	private int parseReference(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	private int parseComment(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
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
	
	private int parseContig(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	private int parseOrigin(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
	private int parsePrimary(ArrayList<byte[]> lines){
		return advanceBlock(lines);
	}
	
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
	
	boolean printGene=false;
	boolean printRepeat=false; 
	
	public static String[] featureTypes=GbffFeature.typeStrings;
	private static final byte[] seqRegB="##sequence-region ".getBytes();
	
	String accession;
	String organism;
	String species;
	int start;
	int stop;
	ArrayList<GbffFeature> features=new ArrayList<GbffFeature>();
}
