package scalar;

import java.util.ArrayList;

import bin.DataLoader;
import clade.Clade;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectDoubleMap;
import parse.LineParser1;
import parse.LineParserS1;
import parse.LineParserS4;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchMakerMini;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;
import tracker.KmerTracker;

public class ScalarData implements Comparable<ScalarData>{

	ScalarData(boolean storeNames, long numericID_){
		if(storeNames) {names=new ArrayList<String>();}
		numericID=numericID_;
	}
	
	public FloatList[] data() {return new FloatList[] {gc, hh, caga, depth};}
	public FloatList[] reorder(int[] order) {
		FloatList[] list=new FloatList[] {gc, hh, caga, depth};
		if(order!=null) {
			list=new FloatList[] {list[order[0]], list[order[1]], list[order[2]]};
		}
		return list;
	}
	
	public ScalarData readTSV(FileFormat ffin1){

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		LineParser1 parser=new LineParser1('\t');

		byte[] line;
		String prev=null;
		boolean newFormat=false;  // Detect format from header
		boolean headerSeen=false;
		while((line=bf.nextLine())!=null){
			if(line.length>0 && line[0]=='#' && !headerSeen){
				// Detect format: new format has "Length" column
				String header=new String(line);
				newFormat=header.contains("Length");
				headerSeen=true;
			}else if(line.length>0 && line[0]!='#'){
				bytesProcessed+=(line.length+1);
				pointsProcessed++;
				parser.set(line);

				int col=0;
				if(newFormat){
					// New format: Name Length GC HH CAGA Depth Start TaxID TaxID2 [...]
					if(names!=null) {
						String name=parser.parseString(col++);
						names.add(name!=null && name.equals(prev) ? prev : name);
						prev=name;
					}else{col++;}
					length.add(parser.parseFloat(col++));
					gc.add(parser.parseFloat(col++));
					hh.add(parser.parseFloat(col++));
					caga.add(parser.parseFloat(col++));
					depth.add(parser.parseFloat(col++));
					start.add(parser.parseInt(col++));
					taxIDs.add(parser.parseInt(col++));
					if(parser.terms()>col){taxIDs2.add(parser.parseInt(col++));}
					else{taxIDs2.add(0);}
				}else{
					// Old format: GC HH CAGA [Depth] TaxID [Name]
					gc.add(parser.parseFloat(col++));
					hh.add(parser.parseFloat(col++));
					caga.add(parser.parseFloat(col++));
					// Check if Depth column exists
					if(parser.terms()>=5){
						depth.add(parser.parseFloat(col++));
					}else{
						depth.add(0);
					}
					taxIDs.add(parser.parseInt(col++));
					taxIDs2.add(0);  // No second taxonomy in old format
					length.add(0);   // No length in old format
					start.add(0);    // No start in old format
					if(parser.terms()>col && names!=null) {
						String name=parser.parseString(col++);
						names.add(name!=null && name.equals(prev) ? prev : name);
						prev=name;
					}
				}
			}
		}
		bf.close();
		return this;
	}
	
	public final void print(String outname, boolean header,
			boolean printName, boolean printPos, int interval) {
		FileFormat ffout=FileFormat.testOutput(outname, FileFormat.TXT, null, true, true, false, false);
		print(ffout, header, printName, printPos, interval);
	}
	
	public final void print(FileFormat ffout, boolean header, 
			boolean printName, boolean printPos, int interval) {
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		print(bsw, header, printName, printPos, interval);
		bsw.poison();
	}
	
	public static String header(boolean sideHeader, boolean printName, boolean printPos) {
		ByteBuilder bb=new ByteBuilder();
		bb.append("#");
		if(sideHeader) {bb.appendt("Header");}
		if(printName) {bb.append("Name\t");}
		bb.append("Length\tGC\tHH\tCAGA\tDepth");
		if(printPos) {bb.append("\tStart");}
		bb.append("\tTaxID\tTaxID2");
		return bb.nl().toString();
	}
	
	public ByteBuilder mean(boolean sideHeader, String name) {
		ByteBuilder bb=new ByteBuilder();
		if(sideHeader) {bb.appendt("Mean");}
		bb.appendt(gc.mean(),5);
		bb.appendt(hh.mean(),5);
		bb.appendt(caga.mean(),5);
		bb.appendt(depth.size()>0 ? depth.mean() : 0,5);
		int tid=(taxIDs==null ? 0 : taxIDs.modeUnsorted());
		bb.appendt(tid);
		if(name!=null) {bb.appendt(name);}
		bb.replaceLast('\n');
		return bb;
	}
	
	public ByteBuilder stdev(boolean sideHeader, String name) {
		ByteBuilder bb=new ByteBuilder();
		if(sideHeader) {bb.appendt("STDev");}
		bb.appendt(gc.stdev(),5);
		bb.appendt(hh.stdev(),5);
		bb.appendt(caga.stdev(),5);
		bb.appendt(depth.size()>0 ? depth.stdev() : 0,5);
		int tid=(taxIDs==null ? 0 : taxIDs.modeUnsorted());
		bb.appendt(tid);
		if(name!=null) {bb.appendt(name);}
		bb.replaceLast('\n');
		return bb;
	}
	
	public final void print(ByteStreamWriter bsw, boolean header,
			boolean printName, boolean printPos, int interval) {
		if(bsw==null) {return;}
		if(header) {bsw.print(header(false, printName, printPos));}
		ByteBuilder bb=new ByteBuilder();
		String prevName=null;
		int pos=0;
		for(int i=0, len=gc.size(); i<len; i++) {
			// Column order: Name, Length, GC, HH, CAGA, Depth, Start, TaxID, TaxID2
			if(names!=null && printName) {
				String name=names.get(i);
				boolean match=(name!=null && name.equals(prevName));
				prevName=name;
				pos=match ? pos+interval : 0;
				bb.appendt(name);
			}
			bb.appendt(length.size()>0 ? length.get(i) : 0, decimals);
			bb.appendt(gc.get(i), decimals);
			bb.appendt(hh.get(i), decimals);
			bb.appendt(caga.get(i), decimals);
			bb.appendt(depth.size()>0 ? depth.get(i) : 0, decimals);
			if(printPos) {bb.appendt(start.size()>0 ? start.get(i) : pos);}
			bb.appendt(taxIDs.size()>0 ? taxIDs.get(i) : 0);
			bb.appendt(taxIDs2.size()>0 ? taxIDs2.get(i) : 0);
			bb.replaceLast('\n');
			bsw.print(bb);
			bb.clear();
		}
	}

	public void add(Read r, KmerTracker dimers, SketchMakerMini smm,
			int interval, int minlen, int tid, boolean breakOnContig) {
		if(r==null) {return;}
		final byte[] bases=r.bases;
		if(makeClade && r.length()>=minCladeSize) {
			if(clade==null) {clade=new Clade(-1, -1, null);}
			clade.add(bases, null);
		}
		if(makeSketch && r.length()>=minSketchSize) {smm.processRead(r);}
		if(breakOnContig && r.length()<minlen) {return;}
		readsProcessed++;
		basesProcessed+=r.length();
		if(breakOnContig) {dimers.clearAll();}
		if(parseTID && tid<0) {tid=bin.BinObject.parseTaxID(r.name());}
		if(dimers.window>0) {
			for(byte b : bases) {
				boolean newValid=dimers.addWindowed(b);
				if(newValid && interval>0 && dimers.count()>=interval) {
					toInterval(dimers, tid, r.name());
				}
			}
		}else {dimers.add(bases);}
		if(verbose) {System.err.println("dimers.count()="+dimers.count()+", minlen="+minlen);}
		if(dimers.count()>=minlen) {toInterval(dimers, tid, r.name());}
	}

	private void toInterval(KmerTracker dimers, int tid, String name) {
		if(verbose) {System.err.println("calling toInterval(dimers)");}
		gc.add(dimers.GC());
		hh.add(dimers.HH());
		caga.add(dimers.CAGA());
		depth.add(getDepth(name));
		length.add(dimers.count());  // Interval length as proxy for now
		start.add(0);  // TODO: Track position within contig
		taxIDs.add(tid);
		taxIDs2.add(0);  // No second taxonomy yet
		if(names!=null) {names.add(name);}
		dimers.resetCount();
		pointsProcessed++;
	}
	
	public void add(ScalarData sd) {
		assert(sd!=this);
		gc.append(sd.gc);
		hh.append(sd.hh);
		caga.append(sd.caga);
		depth.append(sd.depth);
		length.append(sd.length);
		start.addAll(sd.start);
		taxIDs.addAll(sd.taxIDs);
		taxIDs2.addAll(sd.taxIDs2);
		if(names!=null) {names.addAll(sd.names);}
		readsProcessed+=sd.readsProcessed;
		basesProcessed+=sd.basesProcessed;
		bytesProcessed+=sd.bytesProcessed;
		pointsProcessed+=sd.pointsProcessed;
//		if(clade!=null) {clade.add(sd.clade);}
	}

	public int tid(int i){return taxIDs.get(i);}
	public String name(int i){return names==null ? null : names.get(i);}

	/**
	 * Get a row object for the specified index.
	 * @param index Row index
	 * @param si ScalarInterval to populate (or null to create new)
	 * @return ScalarInterval populated with data from this row
	 */
	public ScalarInterval getRow(int index, ScalarInterval si){
		if(si==null){si=new ScalarInterval();}
		si.set(this, index);
		return si;
	}

	/**
	 * Add a row to this ScalarData.
	 * @param row ScalarInterval to add
	 */
	public void addRow(ScalarInterval row){
		if(row.name!=null && names!=null){names.add(row.name);}
		length.add(row.length);
		gc.add(row.gc);
		hh.add(row.hh);
		caga.add(row.caga);
		depth.add(row.depth);
		start.add(row.start);
		taxIDs.add(row.taxID);
		taxIDs2.add(row.taxID2);
		pointsProcessed++;
	}

	@Override
	public int compareTo(ScalarData o){
		return numericID<o.numericID ? -1 : numericID>o.numericID ? 1 : o.gc.size-gc.size;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Depth Loading          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Load depth/coverage data from a file.
	 * @param covFile Coverage file (pileup.sh or covmaker.sh format)
	 * @param depthFile SAM/BAM file for depth calculation
	 */
	public void loadDepth(String covFile, String depthFile){
		if(covFile!=null){
			depthMap=loadCoverageFile(covFile);
		}else if(depthFile!=null){
			depthMap=loadDepthFromAlignment(depthFile);
		}
	}

	/**
	 * Load coverage file (pileup or covmaker format).
	 * @param fname Coverage file path
	 * @return Map of contig name to depth
	 */
	private ObjectDoubleMap<String> loadCoverageFile(String fname){
		if(verbose){System.err.println("Loading coverage from "+fname);}
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		byte[] firstLine=bf.nextLine();
		bf.close();

		String header=new String(firstLine);
		if(header.startsWith("#ID")){
			return loadPileupCoverage(fname);
		}else if(header.startsWith("#Contigs") || header.startsWith("#Depths")){
			return DataLoader.loadCovFile2(fname);
		}else{
			throw new RuntimeException("Unknown coverage file format. Expected pileup (#ID) or covmaker (#Contigs) format.");
		}
	}

	/**
	 * Load coverage from pileup.sh output format.
	 * @param fname Pileup coverage file path
	 * @return Map of contig name to depth
	 */
	private ObjectDoubleMap<String> loadPileupCoverage(String fname){
		LineParser1 lp=new LineParser1('\t');
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		ObjectDoubleMap<String> map=new ObjectDoubleMap<String>();

		byte[] line=bf.nextLine();
		// Skip header
		while(line!=null && Tools.startsWith(line, '#')){
			line=bf.nextLine();
		}

		// Parse data lines
		while(line!=null){
			lp.set(line);
			String name=Tools.trimToWhitespace(lp.parseString(0));
			float depth=lp.parseFloat(1);  // Avg_fold column
			map.put(name, depth);
			line=bf.nextLine();
		}
		bf.close();
		if(verbose){System.err.println("Loaded "+map.size()+" contigs from pileup file.");}
		return map;
	}

	/**
	 * Parse depth from read header in Spades or Tadpole format.
	 * @param name Read name/header
	 * @return Depth value, or -1 if not found
	 */
	public float parseDepthFromHeader(String name){
		if(name==null){return -1;}
		if(name.startsWith("NODE_") && name.contains("_cov_")){//Spades
			lps.set(name);
			return lps.parseFloat(5);
		}else if(name.startsWith("contig_") && name.contains(",cov=")){//Tadpole
			lpt.set(name);
			return lpt.parseFloat(5);
		}else if(name.contains("_cov_")){//Generic
			lps.set(name);
			for(int i=0; i<lps.terms()-1; i++){
				if(lps.termEquals("cov", i)){
					return lps.parseFloat(i+1);
				}
			}
		}
		return -1;
	}

	/**
	 * Get depth for a given contig name.
	 * Tries depthMap first, then parses from header.
	 * @param name Contig name
	 * @return Depth value, or 0 if not found
	 */
	public float getDepth(String name){
		float depth;
		if(depthMap!=null && name!=null){
			String trimmedName=Tools.trimToWhitespace(name);
			depth=(float)depthMap.get(trimmedName);
			if(depth<0){depth=parseDepthFromHeader(name);}
		}else{
			depth=parseDepthFromHeader(name);
		}
		return depth<0 ? 0 : depth;
	}

	/**
	 * Calculate depth from SAM/BAM alignment file.
	 * @param fname SAM/BAM file path
	 * @return Map of contig name to depth
	 */
	private ObjectDoubleMap<String> loadDepthFromAlignment(String fname){
		System.err.print("Loading depth from alignment file "+fname+": ");
		shared.Timer t=new shared.Timer(System.err, false);

		//First pass: get contig lengths from @SQ headers
		ObjectDoubleMap<String> lengthMap=new ObjectDoubleMap<String>();
		ObjectDoubleMap<String> bpMap=new ObjectDoubleMap<String>();

		ByteFile bf=ByteFile.makeByteFile(fname, true);
		byte[] line=bf.nextLine();

		//Parse @SQ headers for contig lengths
		while(line!=null && Tools.startsWith(line, '@')){
			String lineStr=new String(line);
			if(lineStr.startsWith("@SQ")){
				//Parse @SQ\tSN:name\tLN:length
				String[] parts=lineStr.split("\t");
				String name=null;
				int length=0;
				for(String part : parts){
					if(part.startsWith("SN:")){
						name=Tools.trimToWhitespace(part.substring(3));
					}else if(part.startsWith("LN:")){
						length=Integer.parseInt(part.substring(3));
					}
				}
				if(name!=null && length>0){
					lengthMap.put(name, length);
					bpMap.put(name, 0);  //Initialize bp count
				}
			}
			line=bf.nextLine();
		}

		//Second pass: count aligned bp per contig
		int readCount=0;
		while(line!=null){
			if(Tools.startsWith(line, '@')){
				line=bf.nextLine();
				continue;
			}

			//Parse SAM line: QNAME FLAG RNAME POS MAPQ CIGAR
			String lineStr=new String(line);
			String[] parts=lineStr.split("\t");
			if(parts.length>=6){
				int flag=Integer.parseInt(parts[1]);
				String rname=Tools.trimToWhitespace(parts[2]);
				String cigar=parts[5];

				//Skip unmapped reads (FLAG & 4)
				if((flag&4)==0){
					//Calculate aligned bases from CIGAR
					int alignedBases=calculateAlignedBases(cigar);
					if(alignedBases>0 && bpMap.contains(rname)){
						bpMap.put(rname, bpMap.get(rname)+alignedBases);
					}
				}
			}
			readCount++;
			if(readCount%100000==0){System.err.print(".");}
			line=bf.nextLine();
		}
		bf.close();

		//Calculate depth = bp/length
		ObjectDoubleMap<String> depthMap=new ObjectDoubleMap<String>();
		Object[] contigKeys=lengthMap.keys();
		for(Object obj : contigKeys){
			if(obj!=null){
				String contig=(String)obj;
				double length=lengthMap.get(contig);
				double bp=bpMap.get(contig);
				double depth=bp/length;
				depthMap.put(contig, depth);
			}
		}

		t.stopAndPrint();
		System.err.println("Loaded depth for "+depthMap.size()+" contigs from "+readCount+" alignments.");
		return depthMap;
	}

	/**
	 * Calculates aligned bases from CIGAR string.
	 * Counts M, =, X operations (matches/mismatches).
	 * @param cigar CIGAR string
	 * @return Number of aligned bases
	 */
	private static int calculateAlignedBases(String cigar){
		if(cigar.equals("*")){return 0;}
		int total=0;
		int num=0;
		for(int i=0; i<cigar.length(); i++){
			char c=cigar.charAt(i);
			if(Character.isDigit(c)){
				num=num*10+(c-'0');
			}else{
				if(c=='M' || c=='=' || c=='X'){
					total+=num;
				}
				num=0;
			}
		}
		return total;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public FloatList gc=new FloatList();
	public FloatList hh=new FloatList();
	public FloatList caga=new FloatList();
	public FloatList depth=new FloatList();
	public FloatList length=new FloatList();
	public IntList start=new IntList();
	public IntList taxIDs=new IntList();
	public IntList taxIDs2=new IntList();
	public ArrayList<String> names;

	long readsProcessed=0;
	long basesProcessed=0;
	long bytesProcessed=0;
	long pointsProcessed=0;
	
	Clade clade;
	Sketch sketch;

	/** Depth map for coverage lookup by contig name (static so all instances can share) */
	static ObjectDoubleMap<String> depthMap=null;

	/** Line parsers for depth extraction from headers */
	private final LineParserS1 lps=new LineParserS1('_');  // Spades format
	private final LineParserS4 lpt=new LineParserS4(",,=,");  // Tadpole format

	final long numericID;
	
	public static boolean parseTID=false;
	public static boolean makeClade=false;
	public static boolean makeSketch=false;
	public static int minCladeSize=2000;
	public static int minSketchSize=5000;
	public static int decimals=5;

	public static boolean verbose=false;
	
}
