package prot;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.LongList;

/**
 * Consumer-side loader for {@code reference_cds_shred_survival_v1} tables written by
 * {@link ReferenceCdsShredSurvivalTable}, and the ONLY place the per-bin label operands are
 * derived from them. A consumer (MagQCVectorMaker, refcdstable=) keeps one {@code int} table
 * index per cached shred ({@link #requireShred}) and per bin calls {@link #counts} with the
 * selected native and foreign indices; the result is the same three-integer
 * {@link ReferenceCdsSurvivalLabelReader.Counts} the accepted per-bin sidecar path produces.
 *
 * <p>Semantics preserved from the sidecar + label reader (D39): a gene is retained iff one selected
 * source span fully contains it; a gene counts ONCE however many selected spans (or repeated
 * instances of one span) contain it; native_total is the tid's full original CDS universe
 * including genes no shred can hold; foreign genes are counted per foreign tid and summed.</p>
 *
 * <p>Provenance (Yoimiya review 2026-09-09): a consumer loads through a hash-bound TABLE MANIFEST
 * ({@code reference_cds_shred_survival_manifest_v1}, written by {@link #writeManifest} after a full
 * validated load) and must pass the sha256 of the shred list it actually uses (the per-contig
 * cache); every table's bytes are hashed and compared with the manifest before parsing, the
 * tables' {@code #shreds_sha256} must equal the caller's, coverage must be complete, and every
 * declared total is recomputed. Nothing in a table is trusted because the table says so.</p>
 *
 * @author UMP45
 */
public final class ReferenceCdsShredSurvivalTableReader {

	static final String MANIFEST_SCHEMA="reference_cds_shred_survival_manifest_v1";
	static final String MANIFEST_COLUMNS="table_path\ttable_sha256\tshreds_in_assembly\ttids";

	private final String[] ids;
	private final int[] tids;
	private final int[] rangeOff;//per shred: offset into ranges (pairs), size n+1
	private final int[] ranges;//flattened (lo,hi) inclusive ordinal ranges
	private final IntHashMap cdsTotalByTid;
	private final StringIntMap index;
	public final long shredsListed;
	public final String shredsSha256;
	public final int shredCount;

	private ReferenceCdsShredSurvivalTableReader(String[] ids_, int[] tids_, int[] rangeOff_, int[] ranges_, IntHashMap cds_, long listed_, String shredsSha_){
		ids=ids_; tids=tids_; rangeOff=rangeOff_; ranges=ranges_; cdsTotalByTid=cds_; shredsListed=listed_; shredsSha256=shredsSha_; shredCount=ids.length;
		index=new StringIntMap(ids.length*2);
		for(int i=0; i<ids.length; i++){
			if(index.put(ids[i], i)!=-1){throw new IllegalArgumentException("Duplicate shred id across tables: "+ids[i]);}
		}
	}

	public static void main(String[] args){
		if(args.length>=1 && args[0].equalsIgnoreCase("manifest")){
			String tables=null, out=null;
			for(int i=1; i<args.length; i++){
				final String arg=args[i]; final int eq=arg.indexOf('=');
				if(eq<1){throw new RuntimeException("Arguments must be key=value: "+arg);}
				final String key=arg.substring(0, eq).toLowerCase(), value=arg.substring(eq+1);
				if(key.equals("tables")){tables=value;}
				else if(key.equals("out")){out=value;}
				else{throw new RuntimeException("Unknown argument: "+arg);}
			}
			if(tables==null || out==null){throw new RuntimeException("Usage: java -ea prot.ReferenceCdsShredSurvivalTableReader manifest tables=<a.tsv,b.tsv,...> out=<manifest.tsv>");}
			final ReferenceCdsShredSurvivalTableReader r=writeManifest(tables, out);
			System.err.println("ReferenceCdsShredSurvivalTableReader manifest PASS: tables="+tables.split(",").length+" shreds="+r.shredCount+" listed="+r.shredsListed+" tids="+r.cdsTotalByTid.size()+" out="+out);
			return;
		}
		throw new RuntimeException("Usage: java -ea prot.ReferenceCdsShredSurvivalTableReader manifest tables=<a.tsv,...> out=<manifest.tsv>");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Manifest            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Validates the complete table set (full parse, coverage, totals) and writes the hash-bound manifest a consumer loads
	 * through {@link #load}. Table paths are recorded as given; relative paths are resolved against the manifest's directory
	 * at load time, so pass paths relative to that directory (or absolute ones).
	 */
	public static ReferenceCdsShredSurvivalTableReader writeManifest(final String commaSeparatedPaths, final String outPath){
		final File outFile=new File(outPath);
		final File base=(outFile.getAbsoluteFile().getParentFile());
		final ArrayList<String> paths=new ArrayList<String>();
		for(String p : commaSeparatedPaths.split(",")){if(!p.isEmpty()){paths.add(p);}}
		if(paths.isEmpty()){throw new IllegalArgumentException("No table paths");}
		final ArrayList<String> resolved=new ArrayList<String>();
		for(String p : paths){
			final String rp=resolve(base, p);
			if(ReferenceCdsShredSurvivalTable.samePath(rp, outPath)){throw new IllegalArgumentException("Manifest out= names a table it would describe: "+outPath);}
			resolved.add(rp);
		}
		if(outFile.exists()){throw new IllegalArgumentException("Manifest out= already exists; remove it deliberately first: "+outPath);}
		final ArrayList<FileTotals> totals=new ArrayList<FileTotals>();
		final ReferenceCdsShredSurvivalTableReader r=loadTables(resolved, true, null, totals);
		// Failure-safe publication (as the producer): unique owned temp beside out, atomic rename, no REPLACE of a destination that appeared meanwhile.
		final File partial;
		try{partial=File.createTempFile(outFile.getName()+".", ".partial", base);}
		catch(java.io.IOException e){throw new RuntimeException("Could not create a temporary manifest beside "+outPath, e);}
		boolean done=false;
		try{
			final ByteStreamWriter bsw=new ByteStreamWriter(partial.getPath(), true, false, true); bsw.start();
			final ByteBuilder bb=new ByteBuilder(4096);
			bb.append("#schema_version\t").append(MANIFEST_SCHEMA).nl();
			bb.append("#tool\tprot.ReferenceCdsShredSurvivalTableReader manifest").nl();
			bb.append("#shreds_sha256\t").append(r.shredsSha256).nl();
			bb.append("#shreds_listed\t").append(r.shredsListed).nl();
			bb.append("#tables\t").append(paths.size()).nl();
			bb.append("#columns\t").append(MANIFEST_COLUMNS).nl();
			long sum=0;
			for(int i=0; i<paths.size(); i++){
				final FileTotals ft=totals.get(i);
				bb.append(paths.get(i)).tab().append(ReferenceCdsSurvivalLabelReader.sha256File(resolved.get(i))).tab().append(ft.inAssembly).tab().append(ft.aRows).nl();
				sum+=ft.inAssembly;
			}
			bb.append("#end\ttables=").append(paths.size()).append("\tsum_in_assembly=").append(sum).nl();
			bsw.print(bb);
			if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+partial);}
			ReferenceCdsShredSurvivalTable.publishCompleted(partial, outFile, false);//create-if-absent: a manifest never replaces a destination
			done=true;
		}finally{
			if(!done && partial.exists()){partial.delete();}
		}
		return r;
	}

	/**
	 * Consumer entry point. BOTH pins are mandatory 64-hex sha256 values (Yoimiya review 2026-09-09: an unpinned manifest lets a
	 * regenerated table+manifest pair change labels against the same cache): {@code expectedManifestSha256} pins the manifest bytes,
	 * {@code expectedShredsSha256} is the sha256 of the shred list the caller actually uses (the per-contig cache file, whose column 0
	 * is the shred id). Completeness is always required. Unpinned loading exists only package-private ({@link #loadTables}).
	 */
	public static ReferenceCdsShredSurvivalTableReader load(final String manifestPath, final String expectedManifestSha256, final String expectedShredsSha256){
		if(manifestPath==null || manifestPath.isEmpty()){throw new IllegalArgumentException("Missing table manifest path");}
		if(!isHex64(expectedShredsSha256)){throw new IllegalArgumentException("A 64-hex shred-list (cache) sha256 is required to bind the tables to the caller's input; got "+expectedShredsSha256);}
		if(!isHex64(expectedManifestSha256)){throw new IllegalArgumentException("A 64-hex table-manifest sha256 pin is required at the production entry point; got "+expectedManifestSha256);}
		{
			final String actual=ReferenceCdsSurvivalLabelReader.sha256File(manifestPath);
			if(!expectedManifestSha256.equalsIgnoreCase(actual)){throw new IllegalArgumentException("Table manifest SHA-256 mismatch for "+manifestPath+": expected "+expectedManifestSha256+" got "+actual);}
		}
		final File base=new File(manifestPath).getAbsoluteFile().getParentFile();
		final ArrayList<String> paths=new ArrayList<String>(), shas=new ArrayList<String>();
		final LongList inAssembly=new LongList(), tidCounts=new LongList();
		String shredsSha=null; long listed=-1; int declaredTables=-1; boolean schema=false, columns=false, end=false;
		final ByteFile bf=ByteFile.makeByteFile(manifestPath, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		long lineNo=0; boolean ok=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNo++;
				if(line.length==0){continue;}
				if(end){throw bad(manifestPath, lineNo, "content after #end");}
				lp.set(line);
				if(line[0]=='#'){
					if(lp.termEquals("#schema_version", 0)){if(schema || lp.terms()!=2 || !lp.termEquals(MANIFEST_SCHEMA, 1)){throw bad(manifestPath, lineNo, "invalid schema line");} schema=true;}
					else if(lp.termEquals("#shreds_sha256", 0)){if(shredsSha!=null || lp.terms()!=2 || lp.length(1)!=64){throw bad(manifestPath, lineNo, "invalid #shreds_sha256");} shredsSha=lp.parseString(1);}
					else if(lp.termEquals("#shreds_listed", 0)){if(listed>=0 || lp.terms()!=2){throw bad(manifestPath, lineNo, "invalid #shreds_listed");} listed=parseCount(lp, 1, manifestPath, lineNo);}
					else if(lp.termEquals("#tables", 0)){if(declaredTables>=0 || lp.terms()!=2){throw bad(manifestPath, lineNo, "invalid #tables");} declaredTables=parseIntCount(lp, 1, manifestPath, lineNo);}
					else if(lp.termEquals("#columns", 0)){if(columns || lp.terms()!=5 || !new String(line).equals("#columns\t"+MANIFEST_COLUMNS)){throw bad(manifestPath, lineNo, "invalid #columns");} columns=true;}
					else if(lp.termEquals("#end", 0)){
						if(lp.terms()!=3 || !lp.termStartsWith("tables=", 1) || !lp.termStartsWith("sum_in_assembly=", 2)){throw bad(manifestPath, lineNo, "malformed #end");}
						final long t=parseCount(lp, 1, 7, manifestPath, lineNo), s=parseCount(lp, 2, 16, manifestPath, lineNo);
						long sum=0; for(int i=0; i<inAssembly.size(); i++){sum+=inAssembly.get(i);}
						if(t!=paths.size() || s!=sum){throw bad(manifestPath, lineNo, "#end tables="+t+" sum_in_assembly="+s+" but rows give "+paths.size()+"/"+sum);}
						end=true;
					}
					continue;
				}
				if(!schema || !columns){throw bad(manifestPath, lineNo, "row before #schema_version/#columns");}
				if(lp.terms()!=4 || lp.length(0)==0 || lp.length(1)!=64){throw bad(manifestPath, lineNo, "manifest row must be path, 64-hex sha256, shreds_in_assembly, tids");}
				paths.add(lp.parseString(0)); shas.add(lp.parseString(1));
				inAssembly.add(parseCount(lp, 2, manifestPath, lineNo)); tidCounts.add(parseCount(lp, 3, manifestPath, lineNo));
			}
			ok=true;
		}finally{
			if(bf.close() && ok){throw new RuntimeException("I/O error reading "+manifestPath);}
		}
		if(!schema || !columns || !end || shredsSha==null || listed<0 || declaredTables<0){throw new IllegalArgumentException("Incomplete table manifest headers/trailer: "+manifestPath);}
		if(declaredTables!=paths.size()){throw new IllegalArgumentException(manifestPath+": #tables "+declaredTables+" != rows "+paths.size());}
		if(!shredsSha.equalsIgnoreCase(expectedShredsSha256)){throw new IllegalArgumentException("Table manifest shred-list sha256 "+shredsSha+" != the caller's shred list (cache) sha256 "+expectedShredsSha256+"; the tables were built from a different shred list");}
		// Bytes first: every table must hash to its manifest row before it is parsed.
		final ArrayList<String> resolved=new ArrayList<String>();
		for(int i=0; i<paths.size(); i++){
			final String p=resolve(base, paths.get(i));
			final String actual=ReferenceCdsSurvivalLabelReader.sha256File(p);
			if(!actual.equalsIgnoreCase(shas.get(i))){throw new IllegalArgumentException("Table SHA-256 mismatch for "+p+": manifest "+shas.get(i)+" actual "+actual);}
			resolved.add(p);
		}
		final ArrayList<FileTotals> totals=new ArrayList<FileTotals>();
		final ReferenceCdsShredSurvivalTableReader r=loadTables(resolved, true, shredsSha, totals);
		if(r.shredsListed!=listed){throw new IllegalArgumentException("Manifest #shreds_listed "+listed+" != tables' "+r.shredsListed);}
		for(int i=0; i<paths.size(); i++){
			final FileTotals ft=totals.get(i);
			if(ft.inAssembly!=inAssembly.get(i) || ft.aRows!=tidCounts.get(i)){throw new IllegalArgumentException("Manifest row for "+paths.get(i)+" declares in_assembly="+inAssembly.get(i)+" tids="+tidCounts.get(i)+" but the table has "+ft.inAssembly+"/"+ft.aRows);}
		}
		return r;
	}

	static boolean isHex64(final String s){
		if(s==null || s.length()!=64){return false;}
		for(int i=0; i<64; i++){final char c=s.charAt(i); if(!((c>='0' && c<='9') || (c>='a' && c<='f') || (c>='A' && c<='F'))){return false;}}
		return true;
	}

	private static String resolve(final File base, final String p){
		final File f=new File(p);
		return f.isAbsolute() ? p : new File(base, p).getPath();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Table loading         ----------------*/
	/*--------------------------------------------------------------*/

	/** Package-private direct loader (tests, manifest builder): no byte pinning; requireComplete demands full coverage. */
	static ReferenceCdsShredSurvivalTableReader loadTables(final String commaSeparatedPaths, final boolean requireComplete){
		final ArrayList<String> paths=new ArrayList<String>();
		for(String p : commaSeparatedPaths.split(",")){if(!p.isEmpty()){paths.add(p);}}
		return loadTables(paths, requireComplete, null, null);
	}

	private static ReferenceCdsShredSurvivalTableReader loadTables(final ArrayList<String> paths, final boolean requireComplete,
			final String expectedShredsSha, final ArrayList<FileTotals> totalsOut){
		if(paths.isEmpty()){throw new IllegalArgumentException("No table paths");}
		final ArrayList<String> ids=new ArrayList<String>();
		final IntList tids=new IntList(), rangeOff=new IntList(), ranges=new IntList();
		final IntHashMap cds=new IntHashMap(4096), minTotal=new IntHashMap(4096);
		long listed=-1, inAssemblySum=0; String shredsSha=null;
		rangeOff.add(0);
		for(String path : paths){
			final FileTotals ft=readOne(path, ids, tids, rangeOff, ranges, cds, minTotal);
			if(totalsOut!=null){totalsOut.add(ft);}
			if(listed<0){listed=ft.listed; shredsSha=ft.shredsSha;}
			else if(listed!=ft.listed || !shredsSha.equalsIgnoreCase(ft.shredsSha)){
				throw new IllegalArgumentException("Table "+path+" cites a different shred list (listed="+ft.listed+" sha="+ft.shredsSha+") than the first table (listed="+listed+" sha="+shredsSha+")");
			}
			inAssemblySum+=ft.inAssembly;
		}
		if(expectedShredsSha!=null && !expectedShredsSha.equalsIgnoreCase(shredsSha)){throw new IllegalArgumentException("Tables' #shreds_sha256 "+shredsSha+" != expected "+expectedShredsSha);}
		if(requireComplete && inAssemblySum!=listed){
			throw new IllegalArgumentException("Tables cover "+inAssemblySum+" shreds but the shred list has "+listed+"; a per-phylum table is missing or the list differs");
		}
		if(inAssemblySum!=ids.size()){throw new AssertionError("Row count "+ids.size()+" != summed shreds_in_assembly "+inAssemblySum);}
		// Every tid seen on ANY S row (contained=0 included) needs an A row, and every ordinal must lie below its cds_total.
		// IntHashMap.get returns -1 for an absent key (structures/IntHashMap.java:103); stored values are all >= 0.
		for(int tid : minTotal.toArray()){
			final int total=cds.get(tid);
			if(total<0){throw new IllegalArgumentException("S rows for tid "+tid+" but no A row");}
			final int need=minTotal.get(tid);
			if(need>total){throw new IllegalArgumentException("tid "+tid+": ordinal "+(need-1)+" >= cds_total "+total);}
		}
		return new ReferenceCdsShredSurvivalTableReader(ids.toArray(new String[0]), tids.toArray(), rangeOff.toArray(), ranges.toArray(), cds, listed, shredsSha);
	}

	private static final class FileTotals { long listed=-1, inAssembly=-1, aRows=0; String shredsSha=null; }

	private static FileTotals readOne(final String path, final ArrayList<String> ids, final IntList tids, final IntList rangeOff,
			final IntList ranges, final IntHashMap cds, final IntHashMap minTotal){
		final FileTotals ft=new FileTotals();
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		boolean schema=false, end=false, ok=false; long sRows=0, aRows=0, lineNo=0, genesSum=0, containedSum=0, partialSum=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNo++;
				if(line.length==0){continue;}
				if(end){throw bad(path, lineNo, "content after #end trailer");}
				lp.set(line);
				if(line[0]=='#'){
					if(lp.termEquals("#schema_version", 0)){
						if(schema || lp.terms()!=2 || !lp.termEquals(ReferenceCdsShredSurvivalTable.SCHEMA, 1)){throw bad(path, lineNo, "invalid schema line");}
						schema=true;
					}else if(lp.termEquals("#shreds_listed", 0)){
						if(ft.listed>=0 || lp.terms()!=2){throw bad(path, lineNo, "duplicate/invalid #shreds_listed");}
						ft.listed=parseCount(lp, 1, path, lineNo);
					}else if(lp.termEquals("#shreds_in_assembly", 0)){
						if(ft.inAssembly>=0 || lp.terms()!=2){throw bad(path, lineNo, "duplicate/invalid #shreds_in_assembly");}
						ft.inAssembly=parseCount(lp, 1, path, lineNo);
					}else if(lp.termEquals("#shreds_sha256", 0)){
						if(ft.shredsSha!=null || lp.terms()!=3 || lp.length(2)!=64){throw bad(path, lineNo, "duplicate/invalid #shreds_sha256");}
						ft.shredsSha=lp.parseString(2);
					}else if(lp.termEquals("#end", 0)){
						// #end S=<n> A=<m> genes=<g> contained_total=<c> partial_total=<p> ...: every declared total is recomputed.
						if(lp.terms()<6 || !lp.termStartsWith("S=", 1) || !lp.termStartsWith("A=", 2) || !lp.termStartsWith("genes=", 3)
							|| !lp.termStartsWith("contained_total=", 4) || !lp.termStartsWith("partial_total=", 5)){throw bad(path, lineNo, "malformed #end trailer");}
						final long s=parseCount(lp, 1, 2, path, lineNo), a=parseCount(lp, 2, 2, path, lineNo), g=parseCount(lp, 3, 6, path, lineNo);
						final long c=parseCount(lp, 4, 16, path, lineNo), p=parseCount(lp, 5, 14, path, lineNo);
						if(s!=sRows || a!=aRows){throw bad(path, lineNo, "#end trailer S="+s+" A="+a+" but read S="+sRows+" A="+aRows);}
						if(g!=genesSum){throw bad(path, lineNo, "#end genes="+g+" but A rows sum to "+genesSum);}
						if(c!=containedSum || p!=partialSum){throw bad(path, lineNo, "#end contained_total="+c+" partial_total="+p+" but S rows sum to "+containedSum+"/"+partialSum);}
						end=true;
					}
					continue;
				}
				if(!schema){throw bad(path, lineNo, "data before #schema_version");}
				if(lp.termEquals('S', 0)){
					if(lp.terms()!=6){throw bad(path, lineNo, "S row must have 6 fields, got "+lp.terms());}
					final String id=lp.parseString(1);
					final int tid=parseIntCount(lp, 2, path, lineNo), contained=parseIntCount(lp, 3, path, lineNo), partial=parseIntCount(lp, 4, path, lineNo);
					if(id.isEmpty()){throw bad(path, lineNo, "empty shred id");}
					int count=0, prevHi=-1;
					if(lp.termEquals('-', 5)){
						if(contained!=0){throw bad(path, lineNo, "contained_count "+contained+" with no ordinals");}
					}else{
						lp.setBounds(5); final int a=lp.a(), b=lp.b();
						int i=a;
						while(i<b){
							final long lo=parseOrdinal(line, i, b, path, lineNo); i=ordinalEnd(line, i, b);
							long hi=lo;
							if(i<b && line[i]=='-'){
								i++; hi=parseOrdinal(line, i, b, path, lineNo); i=ordinalEnd(line, i, b);
								if(hi<=lo){throw bad(path, lineNo, "malformed ordinal range "+lo+"-"+hi);}
							}
							if(lo<=prevHi){throw bad(path, lineNo, "ordinal ranges not strictly increasing");}
							if(i<b){if(line[i]!=','){throw bad(path, lineNo, "malformed ordinal list separator");} i++; if(i>=b){throw bad(path, lineNo, "trailing comma");}}
							ranges.add((int)lo); ranges.add((int)hi); count+=(int)(hi-lo+1); prevHi=(int)hi;
						}
						if(count!=contained){throw bad(path, lineNo, "contained_count "+contained+" != ordinals listed "+count);}
					}
					if(prevHi+1>minTotal.get(tid)){minTotal.put(tid, prevHi+1);}//absent = -1; contained=0 rows register the tid with 0
					ids.add(id); tids.add(tid); rangeOff.add(ranges.size()/2);
					sRows++; containedSum+=contained; partialSum+=partial;
				}else if(lp.termEquals('A', 0)){
					if(lp.terms()!=7){throw bad(path, lineNo, "A row must have 7 fields, got "+lp.terms());}
					final int tid=parseIntCount(lp, 1, path, lineNo), total=parseIntCount(lp, 2, path, lineNo);
					if(cds.get(tid)>=0){throw bad(path, lineNo, "duplicate A row for tid "+tid+" (within or across tables)");}
					cds.put(tid, total);
					aRows++; genesSum+=total;
				}else{throw bad(path, lineNo, "unknown row type");}
			}
			ok=true;
		}finally{
			if(bf.close() && ok){throw new RuntimeException("I/O error reading "+path);}
		}
		if(!schema){throw new IllegalArgumentException("Missing #schema_version in "+path);}
		if(!end){throw new IllegalArgumentException("Missing #end trailer (truncated table?) in "+path);}
		if(ft.listed<0 || ft.inAssembly<0 || ft.shredsSha==null){throw new IllegalArgumentException("Missing #shreds_listed/#shreds_in_assembly/#shreds_sha256 in "+path);}
		if(ft.inAssembly!=sRows){throw new IllegalArgumentException(path+": #shreds_in_assembly "+ft.inAssembly+" != S rows "+sRows);}
		ft.aRows=aRows;
		return ft;
	}

	/** Non-negative decimal in [0, 2^31): digits only, no overflow (a 10+-digit value is rejected). */
	private static long parseOrdinal(final byte[] line, final int from, final int to, final String path, final long lineNo){
		long v=0; int i=from;
		while(i<to && line[i]>='0' && line[i]<='9'){v=v*10+(line[i]-'0'); i++; if(v>Integer.MAX_VALUE){throw bad(path, lineNo, "ordinal overflow");}}
		if(i==from){throw bad(path, lineNo, "malformed ordinal list");}
		return v;
	}
	private static int ordinalEnd(final byte[] line, int i, final int to){while(i<to && line[i]>='0' && line[i]<='9'){i++;} return i;}

	/** A whole term as a non-negative count (digits only, <= 2^62, overflow checked BEFORE multiplying); offset skips a "key=" prefix. */
	private static long parseCount(final LineParser1 lp, final int term, final String path, final long lineNo){return parseCount(lp, term, 0, path, lineNo);}
	private static long parseCount(final LineParser1 lp, final int term, final int offset, final String path, final long lineNo){
		lp.setBounds(term); final byte[] line=lp.line(); final int a=lp.a()+offset, b=lp.b();
		if(a>=b){throw bad(path, lineNo, "empty count in term "+term);}
		long v=0;
		for(int i=a; i<b; i++){
			final byte c=line[i];
			if(c<'0' || c>'9'){throw bad(path, lineNo, "non-digit in count term "+term);}
			if(v>((1L<<62)-9)/10){throw bad(path, lineNo, "count overflow in term "+term);}
			v=v*10+(c-'0');
		}
		return v;
	}
	/** A count that must fit an int (tids, per-table counts): bounds checked before narrowing. */
	private static int parseIntCount(final LineParser1 lp, final int term, final String path, final long lineNo){
		final long v=parseCount(lp, term, 0, path, lineNo);
		if(v>Integer.MAX_VALUE){throw bad(path, lineNo, "count "+v+" exceeds int range in term "+term);}
		return (int)v;
	}

	private static IllegalArgumentException bad(String f, long n, String s){return new IllegalArgumentException(f+":"+n+": "+s);}

	/*--------------------------------------------------------------*/
	/*----------------            Queries           ----------------*/
	/*--------------------------------------------------------------*/

	/** Table index of a shred (source span id, no __instance_ suffix); throws naming the shred if absent. */
	public int requireShred(final String shredId){
		final int i=index.get(shredId);
		if(i<0){throw new IllegalArgumentException("Shred absent from the reference-CDS shred survival table: "+shredId);}
		return i;
	}

	/** Table index or -1. */
	public int indexOf(final String shredId){return index.get(shredId);}

	public String idOf(final int shredIndex){return ids[shredIndex];}

	public int tidOf(final int shredIndex){return tids[shredIndex];}

	/** Original reference CDS universe of a tid (its A row); throws on an unknown tid. */
	public int cdsTotal(final int tid){
		final int v=cdsTotalByTid.get(tid);//-1 when absent (IntHashMap.get)
		if(v<0){throw new IllegalArgumentException("tid absent from the reference-CDS shred survival table: "+tid);}
		return v;
	}

	public boolean hasTid(final int tid){return cdsTotalByTid.get(tid)>=0;}

	public int tidCount(){return cdsTotalByTid.size();}

	/** Genes fully contained by one shred (its own count; NOT a union). */
	public int containedCount(final int shredIndex){
		int n=0;
		for(int r=rangeOff[shredIndex]; r<rangeOff[shredIndex+1]; r++){n+=ranges[2*r+1]-ranges[2*r]+1;}
		return n;
	}

	/**
	 * The three label operands for one bin: native_total = cds_total(targetTid); native_retained =
	 * distinct genes of targetTid contained by the DISTINCT selected native shreds; foreign_retained =
	 * the same per foreign tid, summed. Lists may repeat an index (repeated instances) -- it adds nothing.
	 */
	public ReferenceCdsSurvivalLabelReader.Counts counts(final int targetTid, final IntList nativeIdx, final IntList foreignIdx){
		final int nativeTotal=cdsTotal(targetTid);
		for(int i=0; i<nativeIdx.size(); i++){
			final int si=nativeIdx.get(i);
			if(tids[si]!=targetTid){throw new IllegalArgumentException("Native shred "+ids[si]+" has tid "+tids[si]+" != target tid "+targetTid);}
		}
		final int nativeRetained=distinctContained(nativeIdx);
		// Foreign: group by tid (sort packed tid<<32|index), union within each tid run.
		final LongList packed=new LongList(Math.max(16, foreignIdx.size()));
		for(int i=0; i<foreignIdx.size(); i++){
			final int si=foreignIdx.get(i);
			if(tids[si]==targetTid){throw new IllegalArgumentException("Foreign shred "+ids[si]+" has the target tid "+targetTid);}
			packed.add((((long)tids[si])<<32) | si);
		}
		packed.sort();
		int foreignRetained=0;
		final IntList run=new IntList(16);
		for(int i=0; i<packed.size();){
			final int tid=(int)(packed.get(i)>>>32);
			run.clear();
			while(i<packed.size() && (int)(packed.get(i)>>>32)==tid){run.add((int)(packed.get(i) & 0xFFFFFFFFL)); i++;}
			foreignRetained+=distinctContained(run);
		}
		return new ReferenceCdsSurvivalLabelReader.Counts(nativeTotal, nativeRetained, foreignRetained);
	}

	/** Distinct ordinals over the union of the (deduplicated) shreds' ranges; all shreds must share one tid. */
	private int distinctContained(final IntList idx){
		if(idx.size()==0){return 0;}
		final int[] sorted=idx.toArray(); Arrays.sort(sorted);
		final LongList pairs=new LongList(16);
		int prev=-1;
		for(int si : sorted){
			if(si==prev){continue;}//repeated instance of the same span
			prev=si;
			assert(tids[si]==tids[sorted[0]]) : "distinctContained mixes tids "+tids[si]+" and "+tids[sorted[0]]+"; caller must group by tid (ordinals are per tid)";
			for(int r=rangeOff[si]; r<rangeOff[si+1]; r++){pairs.add((((long)ranges[2*r])<<32) | ranges[2*r+1]);}
		}
		if(pairs.size()==0){return 0;}
		pairs.sort();
		int count=0, curLo=(int)(pairs.get(0)>>>32), curHi=(int)(pairs.get(0) & 0xFFFFFFFFL);
		for(int i=1; i<pairs.size(); i++){
			final int lo=(int)(pairs.get(i)>>>32), hi=(int)(pairs.get(i) & 0xFFFFFFFFL);
			if(lo<=curHi+1){curHi=Math.max(curHi, hi);}
			else{count+=curHi-curLo+1; curLo=lo; curHi=hi;}
		}
		count+=curHi-curLo+1;
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------      String -> int map       ----------------*/
	/*--------------------------------------------------------------*/

	/** Open-addressing String->int map (no boxing); values are >=0, absent = -1. */
	static final class StringIntMap {
		private String[] keys; private int[] vals; private int mask, size;
		StringIntMap(int expected){
			int cap=16; while(cap<expected*2){cap<<=1;}
			keys=new String[cap]; vals=new int[cap]; mask=cap-1;
		}
		private static int spread(int h){h^=(h>>>16); h*=0x7feb352d; h^=(h>>>15); return h;}
		/** Returns the previous value or -1. */
		int put(String k, int v){
			if(size*2>=keys.length){grow();}
			int i=spread(k.hashCode()) & mask;
			while(keys[i]!=null){
				if(keys[i].equals(k)){final int old=vals[i]; vals[i]=v; return old;}
				i=(i+1) & mask;
			}
			keys[i]=k; vals[i]=v; size++;
			return -1;
		}
		int get(String k){
			int i=spread(k.hashCode()) & mask;
			while(keys[i]!=null){
				if(keys[i].equals(k)){return vals[i];}
				i=(i+1) & mask;
			}
			return -1;
		}
		private void grow(){
			final String[] ok=keys; final int[] ov=vals;
			keys=new String[ok.length*2]; vals=new int[ok.length*2]; mask=keys.length-1; size=0;
			for(int i=0; i<ok.length; i++){if(ok[i]!=null){put(ok[i], ov[i]);}}
		}
		int size(){return size;}
	}
}
