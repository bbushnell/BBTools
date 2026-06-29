package ddl;

import java.io.PrintStream;

import cardinality.DynamicDemiLog;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Reads and writes DDLRecord objects (DynamicDemiLog + metadata) in
 * A48-encoded TSV format.
 *
 * File format:
 * <pre>
 * #tid	504728
 * #name	Meiothermus ruber
 * #file	Mruber.fna.gz
 * #bases	3097457
 * #contigs	1
 * #gc	0.634
 * #len	2048
 * 0m:	$g2	bCT	...	(tab-separated A48 bucket values)
 * </pre>
 *
 * @author Brian Bushnell, Ady
 * @date April 17, 2026
 */
public class DDLLoader {

	/*--------------------------------------------------------------*/
	/*----------------        Loading Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses a tab-delimited A48 data line into a DynamicDemiLog (legacy absolute format).
	 * @param line Raw byte line (A48-encoded bucket values)
	 * @param lp LineParser1 (caller should NOT have called lp.set yet)
	 * @param k K-mer length
	 * @return Populated DynamicDemiLog */
	public static DynamicDemiLog parseDDL(byte[] line, LineParser1 lp, int k){
		return parseDDL(line, lp, k, -1);
	}

	/** Parses a tab-delimited A48 data line into a DynamicDemiLog.
	 * @param line Raw byte line (A48-encoded bucket values)
	 * @param lp LineParser1 (caller should NOT have called lp.set yet)
	 * @param k K-mer length
	 * @param offset GlobalNLZ offset; -1 for legacy absolute format, >=0 for relative
	 * @return Populated DynamicDemiLog */
	public static DynamicDemiLog parseDDL(byte[] line, LineParser1 lp, int k, int offset){
		return parseDDL(line, lp, k, offset, null);
	}

	public static DynamicDemiLog parseDDL(byte[] line, LineParser1 lp, int k, int offset, long[] kmers){
		lp.set(line);
		final int terms=lp.terms();
		final char[] loaded=new char[terms];
		for(int i=0; i<terms; i++){
			loaded[i]=(char)lp.parseLongA48(i);
		}
		return DynamicDemiLog.fromArray(loaded, -1, null, k, offset, kmers);
	}

	/** Parses a tab-delimited A48 kmer data line into a long array. */
	public static long[] parseKmers(byte[] line, LineParser1 lp){
		lp.set(line);
		final int terms=lp.terms();
		final long[] kmers=new long[terms];
		for(int i=0; i<terms; i++){
			kmers[i]=lp.parseLongA48(i);
		}
		return kmers;
	}

	/** Loads all DDLRecords from a TSV file.
	 * @param path Input file path (may be gzipped)
	 * @param k K-mer length
	 * @return List of loaded records */
	public static ArrayList<DDLRecord> loadFile(String path, int k){
		final ArrayList<DDLRecord> records=new ArrayList<DDLRecord>();
		final FileFormat ff=FileFormat.testInput(path, FileFormat.TEXT, null, false, true);
		final ByteFile bf=ByteFile.makeByteFile(ff);
		final LineParser1 lp=new LineParser1('\t');

		long currentId=-1;
		int currentTid=-1;
		String currentName=null, currentFile=null, currentOrigin=null, currentLineage=null;
		long currentBases=0;
		int currentContigs=0;
		float currentGC=-1;
		int currentOffset=-1; // -1 = legacy absolute format
		boolean currentHasKmers=false;

		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length<1){continue;}
			if(line[0]=='#'){
				lp.set(line);
				if(lp.terms()<2){continue;}
				if(lp.termEquals("#k", 0)){
					int fileK=(int)lp.parseLong(1);
					if(fileK!=k){
						outstream.println("WARNING: DDL file k="+fileK+" does not match expected k="+k);
					}
				}
				else if(lp.termEquals("#seed", 0)){/*file-level, informational*/}
				else if(lp.termEquals("#exponent", 0)){DynamicDemiLog.setExponent((int)lp.parseLong(1));}
				else if(lp.termEquals("#id", 0)){currentId=lp.parseLong(1);}
				else if(lp.termEquals("#tid", 0)){currentTid=(int)lp.parseLong(1);}
				else if(lp.termEquals("#name", 0)){currentName=lp.parseString(1);}
				else if(lp.termEquals("#file", 0)){currentFile=lp.parseString(1);}
				else if(lp.termEquals("#bases", 0)){currentBases=lp.parseLong(1);}
				else if(lp.termEquals("#contigs", 0)){currentContigs=(int)lp.parseLong(1);}
				else if(lp.termEquals("#gc", 0)){currentGC=lp.parseFloat(1);}
				else if(lp.termEquals("#origin", 0)){currentOrigin=lp.parseString(1);}
				else if(lp.termEquals("#lineage", 0)){currentLineage=lp.parseString(1);}
				else if(lp.termEquals("#offset", 0)){currentOffset=(int)lp.parseLong(1);}
				else if(lp.termEquals("#haskmers", 0)){currentHasKmers=lp.parseLong(1)>0;}
				else if(lp.termEquals("#blacklist", 0)){lastBlacklistHeader=lp.parseString(1);}
				else if(lp.termEquals("#date", 0)){/*file-level, informational*/}
				continue;
			}

			long[] kmers=null;
			if(currentHasKmers){
				//Safe to read kmerLine before parsing 'line' below: ByteFile.nextLine() returns a fresh
				//array per call (copyOfRange), so this does not alias or overwrite the pending data line.
				byte[] kmerLine=bf.nextLine();
				if(kmerLine!=null && kmerLine.length>0){
					kmers=parseKmers(kmerLine, lp);
				}
			}

			DynamicDemiLog ddl=parseDDL(line, lp, k, currentOffset, kmers);
			DDLRecord rec=new DDLRecord(ddl, currentId, currentTid, currentName);
			rec.filename=currentFile;
			rec.bases=currentBases;
			rec.contigs=currentContigs;
			rec.gc=currentGC;
			rec.origin=currentOrigin;
			rec.lineage=currentLineage;
			rec.cardinality=ddl.cardinality();
			records.add(rec);

			currentId=-1L; currentTid=-1; currentName=null; currentFile=null; currentOrigin=null; currentLineage=null;
			currentBases=0; currentContigs=0; currentGC=-1; currentOffset=-1; currentHasKmers=false;
		}
		bf.close();
		return records;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Writing Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes DDLRecords to a TSV file.
	 * @param records List of records to write
	 * @param path Output file path (may end in .gz)
	 * @param overwrite Overwrite existing files */
	public static void writeFile(ArrayList<DDLRecord> records, String path, boolean overwrite){
		writeFile(records, path, overwrite, -1, -1, null);
	}

	public static void writeFile(ArrayList<DDLRecord> records, String path,
			boolean overwrite, int k, long seed){
		writeFile(records, path, overwrite, k, seed, null);
	}

	public static void writeFile(ArrayList<DDLRecord> records, String path,
			boolean overwrite, int k, long seed, String blacklistFile){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder();
		//Write file-level header
		if(k>0){bb.append("#k").tab().append(k).nl();}
		if(seed!=0){bb.append("#seed").tab().append(seed).nl();}
		bb.append("#exponent").tab().append(DynamicDemiLog.exponentBits()).nl();
		if(blacklistFile!=null){bb.append("#blacklist").tab().append(blacklistFile).nl();}
		bb.append("#date").tab().append(new java.text.SimpleDateFormat("yyyy-MM-dd").format(new java.util.Date())).nl();
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
		for(int i=0; i<records.size(); i++){
			appendRecord(bb, records.get(i));
			bsw.print(bb);
			bb.clear();
		}
		bsw.poisonAndWait();
	}

	/** Appends a full record (headers + A48 data) to a ByteBuilder. */
	public static void appendRecord(ByteBuilder bb, DDLRecord rec){
		if(rec.id>=0){bb.append("#id").tab().append(rec.id).nl();}
		if(rec.taxID>=0){bb.append("#tid").tab().append(rec.taxID).nl();}
		if(rec.name!=null){bb.append("#name").tab().append(rec.name).nl();}
		if(rec.filename!=null){bb.append("#file").tab().append(rec.filename).nl();}
		if(rec.bases>0){bb.append("#bases").tab().append(rec.bases).nl();}
		if(rec.contigs>0){bb.append("#contigs").tab().append(rec.contigs).nl();}
		if(rec.gc>=0){bb.append("#gc").tab().append(rec.gc, 4).nl();}
		if(rec.origin!=null){bb.append("#origin").tab().append(rec.origin).nl();}
		if(rec.lineage!=null){bb.append("#lineage").tab().append(rec.lineage).nl();}
		bb.append("#len").tab().append(rec.ddl.buckets).nl();
		if(rec.ddl.hasKmers()){bb.append("#haskmers").tab().append(1).nl();}
		assert(rec.ddl.getGlobalNLZ()<57) : "Bad offset in appendRecord: "+rec.ddl.getGlobalNLZ();
		if(rec.ddl.getGlobalNLZ()>=0){
			bb.append("#offset").tab().append(rec.ddl.getGlobalNLZ()).nl();
			rec.ddl.toBytesRelative(bb);
		}else{
			rec.ddl.toBytes(bb);
		}
		bb.nl();
		if(rec.ddl.hasKmers()){
			rec.ddl.toBytesKmers(bb);
			bb.nl();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Last #blacklist header value seen during loadFile(). Set during parsing,
	 * readable by callers (e.g., DDLMerger) to propagate to output. */
	public static String lastBlacklistHeader;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final PrintStream outstream=System.err;
}
