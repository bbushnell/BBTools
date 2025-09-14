package gff;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import ml.CellNet;
import prok.Orf;
import prok.ProkObject;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.Feature;
import structures.Range;
import var2.ScafMap;
import var2.VCFLine;
import var2.Var;

/**
 * Parses and represents GFF3 format annotation lines.
 * Supports constructors from multiple input formats including GTF, VCF, Var objects, and ORF predictions.
 * Implements Feature interface for genomic coordinate operations.
 * Provides batch loading with type filtering and prokaryotic feature classification.
 * @author Brian Bushnell
 * @date Sep 12, 2018
 */
public class GffLine implements Comparable<GffLine>, Feature, Cloneable {
	
	/**
	 * Parses tab-delimited GFF line into 9 standard fields.
	 * @param line GFF format line bytes
	 */
	//#seqid	source	type	start	end	score	strand	phase	attributes
	public GffLine(byte[] line){
		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		seqid=parseSeqid ? intern(new String(line, a, b-a)) : null;
//		assert(seqid==null || seqid.equals(new String(line, a, b-a)));
//		assert(seqid!=null) : new String(line, a, b-a)+", "+a+", "+b+"\n"+line;
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		if(b==a+1 && line[a]=='.'){source=DOTS;}
		else{source=paseSource ? intern(new String(line, a, b-a)) : null;}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		if(b==a+1 && line[a]=='.'){type=DOTS;}
		else{
			try { //This was to catch a probably intermittent hardware error; can't replicate.
				type=(parseType ? intern(new String(line, a, b-a)) : null);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println("\n"+new String(line)+"\n"+a+", "+b+", "+(b-a));
				assert(false);
			}
		}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		start=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		stop=Parse.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		if(b<=a){
			//Badly formatted line; common in IMG
			return;
		}
		assert(b>a) : "Missing field 5: "+new String(line);
		if(b==a+1 && line[a]=='.'){score=-1;}
		else{score=Parse.parseFloat(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		assert(b==a+1);
		strand=find(line[a], STRANDS);
//		assert(strand>0) : line[a]+", "+Arrays.toString(STRANDS)+", "+(char)line[b];
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		assert(b==a+1);
		if(line[a]=='.'){phase=-1;}
		else{phase=Parse.parseInt(line, a, b);}
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		if(b==a+1 && line[a]=='.'){attributes=DOTS;}
		else{attributes=parseAttributes ? new String(line, a, b-a) : null;}
		b++;
		a=b;
		
//		assert(strand>=0) : "\n"+this.toString()+"\n"+new String(line);
	}
	
	/** Creates GffLine from GTF format, converting attribute format. */
	public GffLine(GtfLine gtf){
		seqid=gtf.seqname;
		source=gtf.source;
		type=gtf.feature;
		start=gtf.start;
		stop=gtf.end;
		score=gtf.score;
		strand=find(gtf.strand, STRANDS);
		phase=gtf.frame;
		String att=gtf.attribute;
//		if(att.indexOf(';')<0) {att=att.replace(',', ';');}
//		if(att.contains("; ")) {att.replaceAll("; ", ";");}
		StringBuilder sb=new StringBuilder();
		int tnum=0;
		for(String term : att.split(";")) {
			int space=term.indexOf(' ');
			int equals=term.indexOf('=');
			if(tnum>0) {sb.append(';');}
			if(space>0) {
				sb.append(term.substring(0, space)).append('=')
				.append(term.substring(space+1, term.length()));
			}else if(equals>0) {
				sb.append(term);
			}else {
				sb.append("attrib").append(tnum).append('=').append(term);
			}
			tnum++;
		}
		attributes=sb.toString();
	}
	
	/** Creates sequence_variant_obs GffLine from VCF variant call. */
	public GffLine(VCFLine vcf){
		seqid=vcf.scaf;
		source=DOTS;
		type="sequence_variant_obs";
		start=vcf.start()+1;
		stop=vcf.stop()+1;
		score=(float)vcf.qual;
		strand=PLUS;
		phase=-1;
		final int vtype=vcf.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]).append(' ');
		if(vtype==Var.SUB){
			bb.append(vcf.ref).append('>').append(vcf.alt);
		}else if(vtype==Var.DEL){
			bb.append("length ").append(vcf.reflen()-vcf.readlen());
		}else if(vtype==Var.INS){
			int offset=vcf.reflen();
			int length=vcf.readlen()-offset;
			bb.append(vcf.alt, offset, length);
		}else if(vtype==Var.NOCALL){
			bb.append("length ").append(vcf.reflen());
		}
		attributes=bb.toString();
		bb.clear();
	}
	
	/** Creates sequence_variant_obs GffLine from Var with calculated score. */
	public GffLine(Var v, double properPairRate, double totalQualityAvg, double totalMapqAvg, 
			double readLengthAvg, double rarity, int ploidy, ScafMap map, CellNet net){
		seqid=v.scafName();
		source=DOTS;
		type="sequence_variant_obs";
		start=v.start+1;
		stop=Tools.max(v.start+1, v.stop);
		score=(float)v.score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map, net);
		strand=PLUS;
		phase=-1;
		final int vtype=v.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
		attributes=bb.toString();
		bb.clear();
	}
	
	/** Creates sequence_variant_obs GffLine from Var without quality scoring. */
	public GffLine(Var v){
		seqid=v.scafName();
		source="BBTools";
		type="sequence_variant_obs";
		start=v.start+1;
		stop=Tools.max(v.start+1, v.stop);
		score=-1;
		strand=PLUS;
		phase=-1;
		final int vtype=v.type();
		ByteBuilder bb=new ByteBuilder(16);
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
		attributes=bb.toString();
		bb.clear();
	}
	
	@Override
	public GffLine clone() {
		try {
			return (GffLine)super.clone();
		} catch (CloneNotSupportedException e) {
			throw new RuntimeException(e);
		}
	}
	
	/** Creates GffLine from prokaryotic ORF prediction with feature type and scoring. */
	public GffLine(Orf o){
		seqid=o.scafName;
		source="BBTools";
		type=ProkObject.typeStrings2[o.type];
		start=o.start+1;
		stop=o.stop+1;
		score=Tools.max(0, o.orfScore);
		strand=o.strand;
		phase=0;
		ByteBuilder bb=new ByteBuilder(16);
		
		bb.append(ProkObject.typeStrings[o.type]).append(',');
		if(o.type==0){
			bb.append("fr").append(o.frame).append(',');
//			phase=o.frame;
		}
//		bb.append(startCodon).append(',');
//		bb.append(stopCodon).append(',');
		bb.append("startScr:").append(o.startScore, 3).append(',');
		bb.append("stopScr:").append(o.stopScore, 3).append(',');
		bb.append("innerScr:").append(o.averageKmerScore(), 3).append(',');
		bb.append("len:").append(length());
		if(o.type==0){
			bb.append(',');
			bb.append("start:").append(AminoAcid.codonToString(o.startCodon)).append(',');
			bb.append("stop:").append(AminoAcid.codonToString(o.stopCodon));
		}
		attributes=bb.toString();
		bb.clear();
	}
	
	/**
	 * Loads GFF file filtering by feature types.
	 * @param fname File path
	 * @param types Comma-separated feature types (null for all)
	 * @param banUnprocessed Skip unprocessed prokaryotic types
	 * @return List of matching GffLine objects
	 */
	public static ArrayList<GffLine> loadGffFile(String fname, String types, boolean banUnprocessed){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.GFF, null, false, false);
		return loadGffFile(ff, types, banUnprocessed);
	}
	
	/**
	 * Loads GFF file into separate arrays by feature type.
	 * @param gff File path
	 * @param types Comma-separated feature types
	 * @param banUnprocessed Skip unprocessed prokaryotic types
	 * @return Array of lists, one per feature type
	 */
	public static ArrayList<GffLine>[] loadGffFileByType(String gff, String types, boolean banUnprocessed){
		FileFormat ff=FileFormat.testInput(gff, FileFormat.GFF, null, false, false);
		return loadGffFileByType(ff, types, banUnprocessed);
	}
	
	public static ArrayList<GffLine>[] loadGffFileByType(FileFormat ff, String types, boolean banUnprocessed){
		ArrayList<GffLine> list=loadGffFile(ff, types, banUnprocessed);
		String[] typeArray=types.split(",");
		ArrayList<GffLine>[] lists=new ArrayList[typeArray.length];
		for(int i=0; i<typeArray.length; i++){
			String type=typeArray[i];
			lists[i]=new ArrayList<GffLine>();
			for(GffLine gline : list){
				if(gline.type.equalsIgnoreCase(type)){
					lists[i].add(gline);
				}
			}
		}
		return lists;
	}
	
	public static ArrayList<GffLine> loadGffFile(FileFormat ff, String types, boolean banUnprocessed){
		HashSet<String> set=null;
		if(types!=null){
			String[] split=types.split(",");
			set=new HashSet<String>(split.length*2);
			for(String s : split){
				set.add(s);
			}
		}
		
		ArrayList<GffLine> list=new ArrayList<GffLine>();
		ByteFile bf=ByteFile.makeByteFile(ff);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line[0]=='#'){
				//skip
			}else{
				GffLine gline=new GffLine(line);
				assert(gline.strand>=0) : "\n"+gline.toString()+"\n"+new String(line)+"\n";
				if(set==null || (gline.type!=null && set.contains(gline.type))){
					if(!banUnprocessed || ProkObject.processType(gline.prokType())){
						list.add(gline);
					}
				}
			}
		}
		
		boolean error=bf.close();
		assert(!error) : "Problem with file "+ff.name();
		return list;
	}
	
	/**
	 * Appends GFF representation of variant to ByteBuilder.
	 * @param bb Output buffer
	 * @param v Variant object
	 * @param properPairRate Fraction of properly paired reads
	 * @param totalQualityAvg Average base quality
	 * @param totalMapqAvg Average mapping quality
	 * @param readLengthAvg Average read length
	 * @param rarity Variant rarity score
	 * @param ploidy Expected ploidy
	 * @param map Scaffold mapping
	 * @param net Neural network for scoring
	 */
	public static void toText(ByteBuilder bb, Var v, double properPairRate, double totalQualityAvg, 
			double totalMapqAvg, double readLengthAvg, double rarity, int ploidy, ScafMap map, CellNet net){
//		assert(false);
		bb.append(v.scafName(map)).append('\t');
		bb.append('.').append('\t');
		bb.append("sequence_variant_obs").append('\t');
		bb.append(v.start+1).append('\t');
		bb.append(Tools.max(v.start+1, v.stop)).append('\t');
		bb.append(v.score(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, rarity, ploidy, map, net), 2).append('\t');
		bb.append('+').append('\t');
		bb.append('.').append('\t');
//		System.err.println(v.typeString()+", "+v.start+", "+v.stop);
		final int vtype=v.type();
		bb.append("ID=").append(Var.typeArray[vtype]);
		if(vtype==Var.SUB || vtype==Var.INS){
			bb.append(' ').append(v.allele);
		}else if(vtype==Var.DEL || vtype==Var.NOCALL){
			bb.append(" length ").append(v.reflen());
		}else{assert(false) : vtype+"\n"+v;}
	}
	
	/**
	 * Creates GFF header with variant calling statistics.
	 * @param properPairRate Fraction of properly paired reads
	 * @param totalQualityAvg Average base quality
	 * @param mapqAvg Average mapping quality
	 * @param rarity Variant rarity threshold
	 * @param minAlleleFraction Minimum allele fraction for variant calls
	 * @param ploidy Expected ploidy
	 * @param reads Total read count
	 * @param pairs Paired read count
	 * @param properPairs Properly paired read count
	 * @param bases Total base count
	 * @param ref Reference genome name
	 * @return GFF header string with metadata
	 */
	public static String toHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int ploidy, 
			long reads, long pairs, long properPairs, long bases, String ref){
		StringBuilder sb=new StringBuilder();
		
		final double readLengthAvg=bases/Tools.max(1.0, reads);
		sb.append("##gff-version 3\n");
		sb.append("#BBMapVersion\t"+Shared.BBTOOLS_VERSION_STRING+"\n");
		sb.append("#ploidy\t"+ploidy+"\n");
		sb.append(Tools.format("#rarity\t%.5f\n", rarity));
		sb.append(Tools.format("#minAlleleFraction\t%.4f\n", minAlleleFraction));
		sb.append("#reads\t"+reads+"\n");
		sb.append("#pairedReads\t"+pairs+"\n");
		sb.append("#properlyPairedReads\t"+properPairs+"\n");
		sb.append(Tools.format("#readLengthAvg\t%.2f\n", readLengthAvg));
		sb.append(Tools.format("#properPairRate\t%.4f\n", properPairRate));
		sb.append(Tools.format("#totalQualityAvg\t%.4f\n", totalQualityAvg));
		sb.append(Tools.format("#mapqAvg\t%.2f\n", mapqAvg));
		if(ref!=null){sb.append("#reference\t"+ref+"\n");}
		
		sb.append("#seqid	source	type	start	end	score	strand	phase	attributes");
		return sb.toString();
	}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		appendTo(bb);
		return bb.toString();
	}
	
	public ByteBuilder appendTo(ByteBuilder bb){
		bb.append(seqid==null ? "." : seqid).append('\t');
		bb.append(source==null ? "." : source).append('\t');
		bb.append(type==null ? "." : type).append('\t');
		bb.append(start).append('\t');
		bb.append(stop).append('\t');
		if(score<0){bb.append('.').append('\t');}
		else{bb.append(score, 2).append('\t');}
		
		bb.append((strand>=0 ? STRANDS[strand] : (byte)'.')).append('\t');
		
		if(phase<0){bb.append('.').append('\t');}
		else{bb.append(phase).append('\t');}
		
		bb.append(attributes==null ? "." : attributes);
		return bb;
	}
	
	/** Returns feature length in bases. */
	public int length() {
		return stop-start+1;
	}
	
	private static int find(byte a, byte[] array){
		for(int i=0; i<array.length; i++){
			if(array[i]==a){return i;}
		}
		return -1;
	}
	
	private static String intern(String s){
		return Data.forceIntern(s);
	}
	
	@Override
	public int hashCode(){
		return trueStop()^seqid.hashCode();
	}
	
	@Override
	public boolean equals(Object o){
		GffLine b=(GffLine)o;
		if(start!=b.start){return false;}
		if(stop!=b.stop){return false;}
		if(strand!=b.strand){return false;}
		if(!seqid.equals(b.seqid)){return false;}
		if(!type.equals(b.type)){return false;}
		return true;
	}
	
	@Override
	public int compareTo(GffLine g) {
		int x=seqid.compareTo(g.seqid);
		if(x!=0) {return x;}
		if(start!=g.start) {return start-g.start;}
		if(stop!=g.stop) {return stop-g.stop;}
		if(strand!=g.strand) {return strand-g.strand;}
		return type.compareTo(g.type);
	}
	
	@Override
	public int start() {return start-1;}

	@Override
	public int stop() {return stop-1;}

	@Override
	public int strand() {return strand;}

	@Override
	public String name() {return attributes;}

	@Override
	public String seqid() {return seqid;}
	
	/**
	 * Creates Range map from GFF lines grouped by sequence ID.
	 * @param gffLines List of GFF features
	 * @return Map of sequence ID to Range arrays for overlap queries
	 */
	public static HashMap<String, Range[]> makeRangeMap(ArrayList<GffLine> gffLines){
		HashMap<String, ArrayList<GffLine>> lineMap=new HashMap<String, ArrayList<GffLine>>();
		HashMap<String, Range[]> rangeMap=new HashMap<String, Range[]>();
		for(GffLine gline : gffLines) {
			String key=gline.seqid;
			ArrayList<GffLine> list=lineMap.get(key);
			if(list==null) {
				list=new ArrayList<GffLine>();
				lineMap.put(key, list);
			}
			list.add(gline);
		}
		for(Entry<String, ArrayList<GffLine>> e : lineMap.entrySet()) {
			String key=e.getKey();
			ArrayList<GffLine> list=e.getValue();
			Range[] ranges=makeRangesOneSequence(list);
			rangeMap.put(key, ranges);
		}
		return rangeMap;
	}
	
	public static Range[] makeRangesOneSequence(ArrayList<GffLine> listForOneSequence){
		if(listForOneSequence==null || listForOneSequence.isEmpty()) {return null;}
		String name=listForOneSequence.get(0).seqid;
		assert(name.equals(listForOneSequence.get(listForOneSequence.size()-1).seqid));
		Range[] ranges=Range.toRanges(listForOneSequence);
		Range.populateRanges(ranges, listForOneSequence);
		for(int i=0; i<ranges.length; i++) {
			@SuppressWarnings("unchecked")
			ArrayList<GffLine> x=(ArrayList<GffLine>) ranges[i].obj1;
			if(x==null || x.isEmpty()) {ranges[i]=null;}
		}
		ranges=Tools.condenseStrict(ranges);
		return ranges;
	}

//	//These should be sorted and all share the same seqid
//	//Note that this gives one big range if a single feature overlaps everything
//	//Instead, genes should be registered with a minimal number of ranges,
//	//such that each feature start or stop triggers a new range, and each range contains all
//	//features that fully contain it.
//	//To do this you would maintain a list of currently open features...
//	//or...
//	//more simply, make a list of ranges, and register each gene with each one it overlaps.
//	public static Range[] toRanges(ArrayList<GffLine> lines, final int mode){
//
//		ArrayList<Range> list=new ArrayList<Range>(8192);
//		ArrayList<GffLine> glist=new ArrayList<GffLine>(64);
//
//		Range current=null;
//		GffLine prev=null;
//		for(int i=0; i<lines.size(); i++){
//			GffLine g=lines.get(i);
//			assert(prev==null || prev.seqid.equals(g.seqid)) :
//				"Different seqid:\n"+prev+"\n"+g+"\n";
//			assert(prev==null || prev.compareTo(g)<=0 || prev.start==g.start) :
//				"Sort order problem:\n"+prev+"\n"+g+"\n";
//			Range r;
//			
//			int a=g.start, b=g.stop;
//			
//			if(b>=a){//If the size is nonzero
//				r=new Range(a, b);
//				
//				if(current==null){
//					current=r;
//					glist.add(g);
//				}else if(current.touches(r)){
//					current=current.merge(r);
//					glist.add(g);
//				}else{
//					current.obj1=glist.toArray(new GffLine[glist.size()]);
//					glist.clear();
//					glist.add(g);
//					list.add(current);
//					current=r;
//				}
//			}
//		}
//		if(current!=null){ //i.e., if there were any GffLines
//			current.obj1=glist.toArray(new GffLine[glist.size()]);
//			list.add(current);
//		}
//
//		return list.toArray(new Range[list.size()]);
//	}
	
	/** Returns 5' coordinate considering strand orientation. */
	public int trueStart(){
		return strand==0 ? start : stop;
	}
	
	/** Returns 3' coordinate considering strand orientation. */
	public int trueStop(){
		return strand==0 ? stop : start;
	}
	
	/**
	 * Returns prokaryotic feature type constant.
	 * Classifies rRNA by subunit (16S, 23S, 18S, 5S) based on attributes.
	 * @return ProkObject type constant or -1 if unrecognized
	 */
	public final int prokType(){
		if(type.equals("CDS")){
			return ProkObject.CDS;
		}else if(type.equals("tRNA")){
			return ProkObject.tRNA;
		}else if(type.equals("rRNA")){
			if(attributes.contains("16S")){
				return ProkObject.r16S;
			}else if(attributes.contains("23S")){
				return ProkObject.r23S;
			}else if(attributes.contains("18S")){
				return ProkObject.r18S;
			}else if(attributes.contains("5S") && length()<300){
				return ProkObject.r5S;
			}
		}
		return -1;
	}
	
	/** Returns true if feature has partial=true attribute. */
	public final boolean partial(){return attributes!=null && attributes.contains("partial=true");}
	
	/** Returns true if feature coordinates are within scaffold bounds. */
	public final boolean inbounds(int scaflen){return start>=0 && stop<scaflen;}
	
	/** Sequence/scaffold identifier */
	public String seqid;
	/** Feature source (program/database) */
	public String source;
	/** Feature type (CDS, rRNA, gene, etc) */
	public String type;
	/** 1-based start coordinate */
	public int start;
	/** 1-based end coordinate (inclusive) */
	public int stop;
	/** Feature score (-1 if missing) */
	public float score;
	/** Strand (0=+, 1=-, 2=?, 3=.) */
	public int strand;
	/** Coding phase (0, 1, 2, or -1 if missing) */
	public int phase;
	/** Semicolon-separated attribute string */
	public String attributes;
	
	/** Strand character array for byte lookup */
	private static final byte[] STRANDS=new byte[] {'+', '-', '?', '.'};
	/** Strand constants: PLUS=0, MINUS=1, QMARK=2, DOT=3 */
	public static final int PLUS=0, MINUS=1, QMARK=2, DOT=3;
	/** Constant for missing field placeholder */
	public static final String DOTS=".";
	
	/** Controls whether to parse and intern sequence IDs */
	public static boolean parseSeqid=true;
	/** Controls whether to parse and intern source field (typo in variable name) */
	public static boolean paseSource=false;
	/** Controls whether to parse and intern feature types */
	public static boolean parseType=true;
	/** Controls whether to parse score field (unused) */
	public static boolean parseScore=false;
	/** Controls whether to parse attributes field */
	public static boolean parseAttributes=true;
	
	
//	public static boolean parseSeqid=true;
//	public static boolean paseSource=true;
//	public static boolean parseType=true;
//	public static boolean parseScore=true;
//	public static boolean parseAttributes=true;
	
}
