package scalar;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import bin.DataLoader;
import clade.Clade;
import dna.AminoAcid;
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

/**
 * Mutable column-oriented compositional measurements and optional contig labels.
 * Core metric columns must contain one value per row. Optional columns may be empty
 * in legacy callers; a populated column must align with the core columns.
 * Instances, parsers and row mutations are not thread-safe. Configure the shared
 * depth map and process-wide flags before launching workers.
 *
 * @author Brian Bushnell, Neptune
 */
public class ScalarData implements Comparable<ScalarData>{

	/** Validates command-line files before an output can truncate input or existing data.
	 * Repeated input files remain allowed; aliases of an existing output are rejected.
	 * Special stream names keep the standard BBTools handling. */
	static void checkFiles(ArrayList<String> inputs, String output, boolean overwrite, boolean append){
		if(inputs==null || inputs.isEmpty()){throw new IllegalArgumentException("At least one scalar input is required");}
		Tools.testInputFiles(true, true, inputs);
		if(!Tools.testOutputFiles(overwrite, append, false, output)){
			throw new IllegalArgumentException("Cannot write scalar output with overwrite="+overwrite+": "+output);
		}
		if(Tools.isNullFileName(output) || Tools.isSpecialOutputName(output)){return;}
		final Path target=Paths.get(output);
		if(!Files.isRegularFile(target)){return;}
		try{
			for(String input : inputs){
				if(!Tools.isSpecialInputName(input) && Files.isSameFile(Paths.get(input), target)){
					throw new IllegalArgumentException("Scalar input and output refer to the same file: "+input);
				}
			}
		}catch(IOException e){throw new RuntimeException("Cannot check scalar input/output aliases", e);}
	}

	/** Creates empty columns; names are retained only when requested. */
	ScalarData(boolean storeNames, long numericID_){
		if(storeNames){names=new ArrayList<String>();}
		numericID=numericID_;
	}

	/** Returns a new array containing the four borrowed metric columns. */
	public FloatList[] data(){return new FloatList[] {gc, hh, caga, depth};}
	/** Returns the first three requested metric columns, or all four when order is null.
	 * Indices address GC, HH, CAGA and depth in that order; lists are not copied. */
	public FloatList[] reorder(int[] order){
		FloatList[] list=new FloatList[] {gc, hh, caga, depth};
		if(order!=null){
			list=new FloatList[] {list[order[0]], list[order[1]], list[order[2]]};
		}
		return list;
	}

	/** Appends TSV rows to this table and closes the reader it opens.
	 * Legacy data stores GC, HH, CAGA, optional depth, TaxID and optional name.
	 * Modern interval data includes length and may include name/start/end columns.
	 * A missing End is unknown (-1), not Start+Length. */
	public ScalarData readTSV(FileFormat ffin1){
		final ByteFile bf=ByteFile.makeByteFile(ffin1);
		final LineParser1 parser=new LineParser1('\t');
		int[] columns=null;
		int headerTerms=0;
		String prev=null;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				parser.set(line);
				if(line[0]=='#'){
					final int[] candidate=readColumns(parser);
					if(candidate!=null){columns=candidate; headerTerms=parser.terms();}
					continue;
				}
				final int terms=parser.terms();
				if(terms<(columns==null ? 4 : headerTerms)){
					throw new IllegalArgumentException("Incomplete scalar row in "+ffin1.name()+": "+parser);
				}
				//Without a header retain the historical four/five/six-column layout.
				final int ni=columns==null ? (terms>5 ? 5 : -1) : columns[0];
				final int li=columns==null ? -1 : columns[1];
				final int gi=columns==null ? 0 : columns[2];
				final int hi=columns==null ? 1 : columns[3];
				final int ci=columns==null ? 2 : columns[4];
				final int di=columns==null ? (terms>=5 ? 3 : -1) : columns[5];
				final int si=columns==null ? -1 : columns[6];
				final int ti=columns==null ? (terms>=5 ? 4 : 3) : columns[7];
				final int t2i=columns==null ? -1 : columns[8];
				final int ei=columns==null ? -1 : columns[9];
				final float g=parser.parseFloat(gi), h=parser.parseFloat(hi), c=parser.parseFloat(ci);
				final float len=li<0 ? 0 : parser.parseFloat(li), dep=di<0 ? 0 : parser.parseFloat(di);
				final int pos=si<0 ? 0 : parser.parseInt(si), tid=ti<0 ? 0 : parser.parseInt(ti);
				final int tid2=t2i<0 ? 0 : parser.parseInt(t2i);
				if(names!=null){
					String name=ni<0 ? null : parser.parseString(ni);
					if(name!=null && name.equals(prev)){name=prev;}
					names.add(name);
					prev=name;
				}
				gc.add(g); hh.add(h); caga.add(c); depth.add(dep); length.add(len);
				start.add(pos); taxIDs.add(tid); taxIDs2.add(tid2);
				addEnd(ei<0 ? -1 : parser.parseInt(ei));
				bytesProcessed+=line.length+1L;
				pointsProcessed++;
			}
		}finally{
			if(bf.close()){throw new RuntimeException("Error reading scalar table: "+ffin1.name());}
		}
		return this;
	}

	/** Maps a recognized header once; unrelated comment lines leave the schema unchanged.
	 * Core metric columns are mandatory. Missing End uses -1; other optional columns
	 * retain their legacy zero/null defaults. */
	private static int[] readColumns(LineParser1 parser){
		final String[] labels={"Name", "Length", "GC", "HH", "CAGA", "Depth", "Start", "TaxID", "TaxID2", "End"};
		final int[] columns={-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		for(int i=0; i<parser.terms(); i++){
			String label=parser.parseString(i);
			if(i==0){label=label.substring(1);}
			for(int j=0; j<labels.length; j++){
				if(labels[j].equalsIgnoreCase(label)){
					if(columns[j]>=0){throw new IllegalArgumentException("Duplicate scalar column: "+label);}
					columns[j]=i;
				}
			}
		}
		if(columns[2]<0 && columns[3]<0 && columns[4]<0){return null;}
		if(columns[2]<0 || columns[3]<0 || columns[4]<0){
			throw new IllegalArgumentException("Scalar header must contain GC, HH and CAGA: "+parser);
		}
		return columns;
	}

	/** Opens a text output and writes this table with optional name/start columns. */
	public final void print(String outname, boolean header,
			boolean printName, boolean printPos, int interval){
		FileFormat ffout=FileFormat.testOutput(outname, FileFormat.TXT, null, true, true, false, false);
		print(ffout, header, printName, printPos, interval);
	}

	/** Owns, joins and checks the writer opened from ffout; null output discards rows.
	 * The borrowed-writer overload does not close it. Output failures throw. */
	public final void print(FileFormat ffout, boolean header,
			boolean printName, boolean printPos, int interval){
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		if(bsw==null){return;}
		try{print(bsw, header, printName, printPos, interval);}
		finally{
			if(bsw.poisonAndWait()){throw new RuntimeException("Error writing scalar table: "+ffout.name());}
		}
	}

	/** Builds the interval-table header; optional End follows the original columns. */
	public static String header(boolean sideHeader, boolean printName, boolean printPos){
		ByteBuilder bb=new ByteBuilder();
		bb.append("#");
		if(sideHeader){bb.appendt("Header");}
		if(printName){bb.append("Name\t");}
		bb.append("Length\tGC\tHH\tCAGA\tDepth");
		if(printPos){bb.append("\tStart");}
		bb.append("\tTaxID\tTaxID2");
		if(printPos){bb.append("\tEnd");}
		return bb.nl().toString();
	}

	/** Builds a summary row: GC, HH, CAGA, depth and modal TaxID, then optional name.
	 * This shorter summary schema deliberately differs from the interval-table header. */
	public ByteBuilder mean(boolean sideHeader, String name){
		ByteBuilder bb=new ByteBuilder();
		if(sideHeader){bb.appendt("Mean");}
		bb.appendt(gc.mean(), 5);
		bb.appendt(hh.mean(), 5);
		bb.appendt(caga.mean(), 5);
		bb.appendt(depth.size()>0 ? depth.mean() : 0, 5);
		int tid=(taxIDs==null ? 0 : taxIDs.modeUnsorted());
		bb.appendt(tid);
		if(name!=null){bb.appendt(name);}
		bb.replaceLast('\n');
		return bb;
	}

	/** Builds population standard deviations using the same schema as mean(). */
	public ByteBuilder stdev(boolean sideHeader, String name){
		ByteBuilder bb=new ByteBuilder();
		if(sideHeader){bb.appendt("STDev");}
		bb.appendt(gc.stdev(), 5);
		bb.appendt(hh.stdev(), 5);
		bb.appendt(caga.stdev(), 5);
		bb.appendt(depth.size()>0 ? depth.stdev() : 0, 5);
		int tid=(taxIDs==null ? 0 : taxIDs.modeUnsorted());
		bb.appendt(tid);
		if(name!=null){bb.appendt(name);}
		bb.replaceLast('\n');
		return bb;
	}

	/** Writes rows to a borrowed writer; null output discards them.
	 * Empty optional columns use zeros, or inferred positions for missing starts.
	 * Missing ends are -1, never inferred from Length, which may be a sampling weight.
	 * A requested Name column is blank when name storage was disabled. */
	public final void print(ByteStreamWriter bsw, boolean header,
			boolean printName, boolean printPos, int interval){
		if(bsw==null){return;}
		if(header){bsw.print(header(false, printName, printPos));}
		ByteBuilder bb=new ByteBuilder();
		String prevName=null;
		int pos=0;
		for(int i=0, len=gc.size(); i<len; i++){
			//Optional Name/Start/End surround the original metric/taxonomy columns.
			if(printName){
				String name=(names==null ? null : names.get(i));
				boolean match=(name!=null && name.equals(prevName));
				prevName=name;
				pos=match ? pos+interval : 0;
				bb.appendt(names==null ? "" : name);
			}
			//scalar/ScalarData#001: populated optional columns must align with core rows.
			//Empty legacy columns use defaults. FloatList.mean/stdev return zero when empty.
			bb.appendt(length.size()>0 ? length.get(i) : 0, decimals);
			bb.appendt(gc.get(i), decimals);
			bb.appendt(hh.get(i), decimals);
			bb.appendt(caga.get(i), decimals);
			bb.appendt(depth.size()>0 ? depth.get(i) : 0, decimals);
			if(printPos){bb.appendt(start.size()>0 ? start.get(i) : pos);}
			bb.appendt(taxIDs.size()>0 ? taxIDs.get(i) : 0);
			bb.appendt(taxIDs2.size()>0 ? taxIDs2.get(i) : 0);
			if(printPos){bb.appendt(end.size()>0 ? end.get(i) : -1);}
			bb.replaceLast('\n');
			bsw.print(bb);
			bb.clear();
		}
	}

	/** Adds one sequence to caller-owned dimer state and emits sampled measurements.
	 * window counts retained valid dimers, while interval/minlen count valid input bases.
	 * breakOnContig resets tracking between retained sequences. A non-null sketch maker
	 * is required when makeSketch is enabled. The mate is the caller's responsibility.
	 * Start/End bound the retained dimers, zero-based and half-open; ambiguity gaps
	 * inside that span need not contribute counts. Empty dimer sets and break=false
	 * have unknown coordinates (-1). Length remains the newly represented valid bases. */
	public void add(Read r, KmerTracker dimers, SketchMakerMini smm,
			int interval, int minlen, int tid, boolean breakOnContig){
		if(r==null){return;}
		assert(dimers.k==2) : "Scalar GC/HH/CAGA and coordinate spans require a dinucleotide tracker";
		final byte[] bases=r.bases;
		if(makeClade && r.length()>=minCladeSize){
			if(clade==null){clade=new Clade(-1, -1, null);}
			clade.add(bases, null);
		}
		if(makeSketch && r.length()>=minSketchSize){smm.processRead(r);}
		if(breakOnContig && r.length()<minlen){return;}
		readsProcessed++;
		basesProcessed+=r.length();
		if(breakOnContig){dimers.clearAll();}
		if(parseTID && tid<0){tid=bin.BinObject.parseTaxID(r.name());}
		int spanStart=-1, spanEnd=-1, newestEnd=-1, pending=0;
		if(dimers.window>0){
			for(int i=0; i<bases.length; i++){
				boolean newValid=dimers.addWindowed(bases[i]);
				if(newValid && breakOnContig){
					pending++;
					newestEnd=i+1;
				}
				if(newValid && interval>0 && dimers.count()>=interval){
					if(breakOnContig){
						spanStart=advanceWindowStart(bases, spanStart, spanEnd, newestEnd, pending, dimers.window);
						spanEnd=newestEnd;
						pending=0;
					}
					toInterval(dimers, tid, r.name(), spanStart, spanEnd);
				}
			}
		}else{dimers.add(bases);}
		if(verbose){System.err.println("dimers.count()="+dimers.count()+", minlen="+minlen);}
		if(dimers.count()>=minlen){
			if(breakOnContig && newestEnd>=0){
				spanStart=advanceWindowStart(bases, spanStart, spanEnd, newestEnd, pending, dimers.window);
				spanEnd=newestEnd;
			}else if(breakOnContig && bases!=null){
				//A partial ring (or unlimited accumulation) has never evicted a dimer.
				spanStart=nextDimerStart(bases, 0, bases.length-2);
				if(spanStart>=0){
					int i=bases.length-2;
					while(AminoAcid.baseToNumber[bases[i]]<0 || AminoAcid.baseToNumber[bases[i+1]]<0){i--;}
					spanEnd=i+2;
				}
			}
			toInterval(dimers, tid, r.name(), spanStart, spanEnd);
		}
	}

	/** Advances the oldest retained pair only when a row is needed. A contiguous
	 * span moves arithmetically; gapped spans use one validity lookup per scanned
	 * base. No coordinates are recomputed for each intervening tracker insertion. */
	private static int advanceWindowStart(byte[] bases, int from, int previousEnd,
			int newestEnd, int pending, int window){
		if(pending==0){return from;}
		if(from<0){
			from=nextDimerStart(bases, 0, newestEnd-2);
			pending--;//The first full ring initializes rather than evicts its oldest pair.
		}
		assert(from>=0) : "A full KmerTracker ring must contain a valid dimer";
		if(pending==0){return from;}
		if(newestEnd-previousEnd==pending && previousEnd-from-1==window){return from+pending;}
		int pos=from+1;//The previous oldest pair guarantees this base is valid.
		int previous=AminoAcid.baseToNumber[bases[pos]];
		while(pending>0){
			//At most one pair per new base: scan that lower bound without a branch
			//per pair, then make up exactly the missing pairs in the next block.
			final int limit=pos+pending;
			int missing=0;
			while(pos<limit){
				final int current=AminoAcid.baseToNumber[bases[++pos]];
				missing+=(previous|current)>>>31;
				previous=current;
			}
			pending=missing;
		}
		assert(pos<newestEnd) : "Oldest retained pair cannot follow the newest pair";
		return pos-1;
	}

	/** Finds the first valid pair start in an inclusive range, or -1 if none.
	 * Initializes full, partial and unlimited metric spans. */
	private static int nextDimerStart(byte[] bases, int from, int last){
		assert(from>=0 && last<bases.length-1) : "Dimer starts must leave room for the second base";
		for(int i=from; i<=last; i++){
			if(AminoAcid.baseToNumber[bases[i+1]]<0){i++;}
			else if(AminoAcid.baseToNumber[bases[i]]>=0){return i;}
		}
		return -1;
	}

	/** Samples current counts and resets only the emission counter, not the rolling window. */
	private void toInterval(KmerTracker dimers, int tid, String name, int spanStart, int spanEnd){
		if(verbose){System.err.println("calling toInterval(dimers)");}
		gc.add(dimers.GC());
		hh.add(dimers.HH());
		caga.add(dimers.CAGA());
		depth.add(getDepth(name));
		length.add(dimers.count());//Represented-base weight, not the overlapping window span.
		start.add(spanStart);
		addEnd(spanEnd);
		taxIDs.add(tid);
		taxIDs2.add(0);  // No second taxonomy yet
		if(names!=null){names.add(name);}
		dimers.resetCount();
		pointsProcessed++;
	}

	/** Appends columns and counters from a distinct table; clade/sketch objects are not merged. */
	public void add(ScalarData sd){
		assert(sd!=this);
		final int oldRows=gc.size();
		assert((end.size()==0 || end.size()==oldRows) &&
				(sd.end.size()==0 || sd.end.size()==sd.gc.size())) : "End columns must be empty or row-aligned before merging";
		while(end.size()<oldRows){end.add(-1);}//Pad legacy tables without End.
		if(sd.end.size()>0){end.addAll(sd.end);}
		else{for(int i=0; i<sd.gc.size(); i++){end.add(-1);}}
		gc.append(sd.gc);
		hh.append(sd.hh);
		caga.append(sd.caga);
		depth.append(sd.depth);
		length.append(sd.length);
		start.addAll(sd.start);
		taxIDs.addAll(sd.taxIDs);
		taxIDs2.addAll(sd.taxIDs2);
		if(names!=null){names.addAll(sd.names);}
		readsProcessed+=sd.readsProcessed;
		basesProcessed+=sd.basesProcessed;
		bytesProcessed+=sd.bytesProcessed;
		pointsProcessed+=sd.pointsProcessed;
//		if(clade!=null){clade.add(sd.clade);}
	}

	/** Returns the primary taxonomic label at row i. */
	public int tid(int i){return taxIDs.get(i);}
	/** Returns the stored name, or null when names are disabled. */
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
		if(names!=null){names.add(row.name);}
		length.add(row.length);
		gc.add(row.gc);
		hh.add(row.hh);
		caga.add(row.caga);
		depth.add(row.depth);
		start.add(row.start);
		addEnd(row.end);
		taxIDs.add(row.taxID);
		taxIDs2.add(row.taxID2);
		pointsProcessed++;
	}

	/** Appends the newest row's end, padding pre-existing legacy rows as unknown. */
	private void addEnd(int value){
		assert(end.size()<gc.size()) : "Append exactly one End after adding the row's core metrics";
		while(end.size()<gc.size()-1){end.add(-1);}
		end.add(value);
	}

	/** Orders by input ID, then by descending row count for equal IDs. */
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
		final ByteFile bf=ByteFile.makeByteFile(fname, true);
		final byte[] firstLine;
		try{firstLine=bf.nextLine();}
		finally{
			if(bf.close()){throw new RuntimeException("Error reading coverage file: "+fname);}
		}
		if(firstLine==null){throw new IllegalArgumentException("Empty coverage file: "+fname);}

		String header=new String(firstLine);
		if(header.startsWith("#ID")){
			return loadPileupCoverage(fname);
		}else if(header.startsWith("#Contigs") || header.startsWith("#Depths")){
			//The covmaker loader owns its reader and rejects read or close failures.
			return DataLoader.loadCovFile2(fname);
		}else{
			throw new RuntimeException("Unknown coverage file format. Expected pileup (#ID) or covmaker (#Contigs) format.");
		}
	}

	/**
	 * Loads pileup coverage and closes the owned reader even when parsing fails.
	 * Read errors reject the partial map rather than returning incomplete coverage.
	 * @param fname Pileup coverage file path
	 * @return Map of contig name to depth
	 */
	private ObjectDoubleMap<String> loadPileupCoverage(String fname){
		LineParser1 lp=new LineParser1('\t');
		final ByteFile bf=ByteFile.makeByteFile(fname, true);
		ObjectDoubleMap<String> map=new ObjectDoubleMap<String>(String.class);
		try{
			byte[] line=bf.nextLine();
			//Skip header.
			while(line!=null && Tools.startsWith(line, '#')){line=bf.nextLine();}
			while(line!=null){
				lp.set(line);
				String name=Tools.trimToWhitespace(lp.parseString(0));
				float depth=lp.parseFloat(1);//Avg_fold column
				map.put(name, depth);
				line=bf.nextLine();
			}
		}finally{
			if(bf.close()){throw new RuntimeException("Error reading pileup coverage: "+fname);}
		}
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
			//Contig.toHeader emits len, cov, gc; the ordered =,=, boundaries put cov at3.
			return lpt.parseFloat(3);
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
	 * Calculates depth from SAM/BAM and closes the owned reader on all exits.
	 * Read errors reject partial counts; only M, = and X CIGAR operations contribute.
	 * @param fname SAM/BAM file path
	 * @return Map of contig name to depth
	 */
	private ObjectDoubleMap<String> loadDepthFromAlignment(String fname){
		System.err.print("Loading depth from alignment file "+fname+": ");
		shared.Timer t=new shared.Timer(System.err, false);

		//First pass: get contig lengths from @SQ headers
		ObjectDoubleMap<String> lengthMap=new ObjectDoubleMap<String>(String.class);
		ObjectDoubleMap<String> bpMap=new ObjectDoubleMap<String>(String.class);

		final ByteFile bf=ByteFile.makeByteFile(fname, true);
		int readCount=0;
		try{
			byte[] line=bf.nextLine();
			//Parse @SQ headers for contig lengths.
			while(line!=null && Tools.startsWith(line, '@')){
				String lineStr=new String(line);
				if(lineStr.startsWith("@SQ")){
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
						bpMap.put(name, 0);
					}
				}
				line=bf.nextLine();
			}
			//Continue through alignments, accumulating aligned bases per contig.
			while(line!=null){
				if(Tools.startsWith(line, '@')){
					line=bf.nextLine();
					continue;
				}
				String lineStr=new String(line);
				String[] parts=lineStr.split("\t");
				if(parts.length>=6){
					int flag=Integer.parseInt(parts[1]);
					String rname=Tools.trimToWhitespace(parts[2]);
					String cigar=parts[5];
					if((flag&4)==0){
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
		}finally{
			if(bf.close()){throw new RuntimeException("Error reading depth alignment: "+fname);}
		}

		//Calculate depth = bp/length
		ObjectDoubleMap<String> depthMap=new ObjectDoubleMap<String>(String.class);
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

	/** GC fraction derived from counted dimer terminal bases. */
	public FloatList gc=new FloatList();
	/** GC-normalized homopolymer dimer metric. */
	public FloatList hh=new FloatList();
	/** GC-normalized CAGA dimer metric. */
	public FloatList caga=new FloatList();
	/** Coverage from the shared map or sequence header, with zero for unavailable depth. */
	public FloatList depth=new FloatList();
	/** Producer-supplied lengths; sequence sampling uses newly represented valid bases. */
	public FloatList length=new FloatList();
	/** Zero-based inclusive metric-span starts; -1 when sequence coordinates are unavailable. */
	public IntList start=new IntList();
	/** Exclusive metric-span ends; -1 when unavailable or absent from an older table. */
	public IntList end=new IntList();
	/** Primary taxonomic assignment, with nonpositive values indicating no assignment. */
	public IntList taxIDs=new IntList();
	/** Secondary taxonomic assignment; zero when unavailable. */
	public IntList taxIDs2=new IntList();
	/** Optional row-aligned names; null disables name storage. */
	public ArrayList<String> names;

	/** Retained input sequences processed by add. */
	long readsProcessed=0;
	/** Bases in those sequences, including ambiguous symbols. */
	long basesProcessed=0;
	/** TSV data-row bytes plus one presumed terminator per row. */
	long bytesProcessed=0;
	/** Rows appended by TSV loading or sequence sampling. */
	long pointsProcessed=0;

	/** Optional per-sequence compositional classifier input. */
	Clade clade;
	/** Optional per-sequence sketch for remote labeling. */
	Sketch sketch;

	//NOTE [scalar/ScalarData#002] LOW: depthMap is a SHARED-MUTABLE STATIC ("so all instances can share"). Footgun: if two
	//ScalarData instances loadDepth() DIFFERENT files, the second overwrites the first's map for BOTH instances (and getDepth
	//reads whichever was loaded last). Fine for the single-run CLI use, but a latent cross-contamination trap for any multi-
	//dataset use. Also loadDepthFromAlignment's Integer.parseInt(parts[1]) NFEs on a non-numeric FLAG (guarded only by len>=6),
	//and parseDepthFromHeader's parseFloat(5) assumes term 5 exists. All LOW (recent, well-structured; these are edges).
	/** Depth map for coverage lookup by contig name (static so all instances can share) */
	static ObjectDoubleMap<String> depthMap=null;

	/** Line parsers for depth extraction from headers */
	private final LineParserS1 lps=new LineParserS1('_');  // Spades format
	private final LineParserS4 lpt=new LineParserS4("=,=,");  // Tadpole format

	/** Source ordering key, fixed at construction. */
	final long numericID;

	/** Whether sequence/file headers may supply taxonomic IDs. */
	public static boolean parseTID=false;
	/** Whether add builds per-sequence clade inputs. */
	public static boolean makeClade=false;
	/** Whether add feeds a caller-provided sketch maker. */
	public static boolean makeSketch=false;
	/** Minimum sequence length submitted to the clade accumulator. */
	public static int minCladeSize=2000;
	/** Minimum sequence length submitted to the sketch accumulator. */
	public static int minSketchSize=5000;
	/** Fractional digits in interval output. Summary rows retain five digits. */
	public static int decimals=5;

	/** Emits diagnostic progress to stderr when true. */
	public static boolean verbose=false;

}
