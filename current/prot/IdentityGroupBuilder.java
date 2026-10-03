package prot;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Builds deterministic connected-component identity groups from an explicit mmseqs
 * createdb-&gt;prefilter-&gt;align-&gt;convertalis edge list, NOT `easy-cluster`'s own clustering
 * module -- see mag-qc results/mmseqs_grouping_fragrich_verdict_v1.md for why easy-cluster's
 * internal workflow was found to drop real edges that survive prefilter+align on their own.
 *
 * <p>Reads the full member-ID universe (so singletons with zero surviving edges are still
 * reported) and a convertalis edge list with `--format-output
 * query,target,nident,alnlen,qstart,qend,qlen,tstart,tend,tlen` (raw integer fields only --
 * identity and coverage are decided by exact integer cross-multiplication, never by parsing
 * mmseqs' own rounded percentage columns), applies caller-specified thresholds, and unions
 * every pair that clears both. Output is cluster.tsv-compatible (rep\tmember, one row per
 * member, canonically sorted) plus a sidecar with a whole-partition SHA-256 and per-group
 * membership hashes/sizes.</p>
 *
 * <p>Determinism: the member-ID universe is loaded then sorted lexicographically before
 * indexing (index order never depends on input file order); edges are streamed (never
 * retained) and unioned as read (union is commutative and path-compressed, so processing
 * order never affects the final partition); output rows are sorted by (representative,
 * member) rather than by any input order. The representative of each component is its
 * lexicographically smallest member ID.</p>
 *
 * <p>Usage: {@code idgroupbuilder.sh ids=<members.txt> edges=<convertalis.m8> out=<cluster.tsv>
 * minid=90 minqcov=80 mintcov=80}</p>
 *
 * <p><b>Edge-list generation gotcha (found 2026-08-30):</b> {@code mmseqs align}'s default
 * {@code --alignment-mode 0} (automatic) silently reports {@code nident=0} for every row,
 * including perfect self-hits. {@code --alignment-mode 3} is required to compute identity at
 * all, but {@code nident} still reports 0 unless the {@code -a} flag (full backtrace/cigar) is
 * also passed -- {@code -a} is what actually populates {@code nident}/{@code qaln}/{@code taln}.
 * Also: mmseqs' own displayed {@code pident} column uses a DIFFERENT convention than the exact
 * {@code nident/alnlen} ratio this tool requires (observed case: pident=89.600 vs exact
 * nident/alnlen=339/378=89.6825% for the same row) -- never threshold on mmseqs' {@code pident}
 * column, always on the raw integer fields. Frozen generation command:
 * {@code mmseqs align queryDB targetDB prefDB alnDB -a --min-seq-id 0.0 -c 0.0 --cov-mode 0
 * --seq-id-mode 0 --alignment-mode 3 -e <liberal> --max-accept <high> --max-rejected <high>},
 * then {@code mmseqs convertalis ... --format-output
 * query,target,nident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits}.</p>
 *
 * <p><b>I/O error propagation</b> (Elly's review, 2026-08-31, ahead of running this tool at
 * 4432-family batch scale): every {@code ByteFile.close()}/{@code ByteStreamWriter.poisonAndWait()}
 * boolean is checked and thrown on -- a latched read or write error (e.g. a truncated write from a
 * disk-full or filesystem hiccup mid-batch) previously went undetected, which matters far more at
 * thousands of unattended invocations than it did during hand-run stress-family testing.
 *
 * @author Eru
 */
public final class IdentityGroupBuilder {

	public static void main(String[] args){
		String idsFile=null, edgesFile=null, outFile=null;
		String minIdStr="90", minQCovStr="80", minTCovStr="80";
		boolean overwrite=true;

		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("ids")){idsFile=b;}
			else if(a.equals("edges")){edgesFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("minid")){minIdStr=b;}
			else if(a.equals("minqcov")){minQCovStr=b;}
			else if(a.equals("mintcov")){minTCovStr=b;}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(idsFile==null || edgesFile==null || outFile==null){
			throw new RuntimeException("Required: ids=<member list> edges=<m8 file> out=<cluster.tsv>");
		}

		//Percent thresholds parsed EXACTLY from their decimal text (never via a double
		//intermediate, per Elly's review -- a double cannot exactly represent most decimal
		//fractions, so Math.round(x*1e6) on a double-parsed input is not provably exact for
		//arbitrary sweep values like "89.6825"). BigDecimal(String) parses the literal digits;
		//exactRatio() converts to an exact (numerator, scale) pair with scale a power of ten.
		final long[] idRatio=exactPercentRatio(minIdStr);
		final long[] qCovRatio=exactPercentRatio(minQCovStr);
		final long[] tCovRatio=exactPercentRatio(minTCovStr);

		final HashMap<String,Integer> idToIndex=new HashMap<String,Integer>();
		final ArrayList<String> indexToId=new ArrayList<String>();
		loadIds(idsFile, idToIndex, indexToId);
		final int n=indexToId.size();
		if(n==0){throw new RuntimeException("No IDs loaded from "+idsFile);}

		final int[] parent=new int[n];
		final int[] rank=new int[n];
		for(int i=0; i<n; i++){parent[i]=i;}

		final long[] counters=new long[3];//[0]=edges read, [1]=edges passing threshold, [2]=self-hits skipped
		processEdges(edgesFile, idToIndex, parent, rank, idRatio, qCovRatio, tCovRatio, counters);

		final int groupCount=writeClusters(outFile, overwrite, parent, indexToId);

		System.err.println("IdentityGroupBuilder: "+n+" members, "+groupCount+" groups, "+counters[0]
			+" edges read, "+counters[2]+" self-hits skipped, "+counters[1]+" edges passed threshold "
			+"(minid="+minIdStr+" minqcov="+minQCovStr+" mintcov="+minTCovStr+").");
	}

	/**
	 * Parses a percent threshold string EXACTLY (no double intermediate) into a
	 * {numerator, scale} pair such that the true value equals numerator/scale, scale a power
	 * of ten matching the input's own decimal places (e.g. "89.6825" -&gt; {896825, 10000}).
	 * Rejects negative values, values over 100, and any non-decimal text (BigDecimal's own
	 * NumberFormatException propagates unconditionally).
	 */
	static long[] exactPercentRatio(String pctStr){
		final BigDecimal bd=new BigDecimal(pctStr);
		if(bd.signum()<0 || bd.compareTo(new BigDecimal(100))>0){
			throw new IllegalArgumentException("Percent threshold must be in [0,100]: "+pctStr);
		}
		final BigInteger unscaled=bd.unscaledValue();
		final int scaleExp=bd.scale();
		if(scaleExp<0){
			//e.g. "9E1" (90 with a negative BigDecimal scale) -- normalize to a nonnegative scale.
			return new long[]{unscaled.multiply(BigInteger.TEN.pow(-scaleExp)).longValueExact(), 1L};
		}
		return new long[]{unscaled.longValueExact(), BigInteger.TEN.pow(scaleExp).longValueExact()};
	}

	/**
	 * Loads the full member-ID universe, one ID per line, then sorts lexicographically so
	 * index order is independent of the input file's line order (determinism requirement).
	 * Throws unconditionally (not assert) on a duplicate or empty ID -- both indicate a
	 * corrupt manifest, not a recoverable condition.
	 */
	static void loadIds(String fname, HashMap<String,Integer> idToIndex, ArrayList<String> indexToId){
		final ByteFile bf=ByteFile.makeByteFile(fname, false);
		final ArrayList<String> raw=new ArrayList<String>();
		final HashMap<String,Integer> seen=new HashMap<String,Integer>();
		int lineNo=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			lineNo++;
			if(line.length>0 && line[0]=='#'){continue;}//comment lines only -- a blank line is an error, not skipped
			final String id=new String(line);
			if(id.isEmpty()){throw new RuntimeException("Blank member-ID line "+lineNo+" of "+fname+" -- fail loud, not silently skipped.");}
			final Integer prior=seen.put(id, lineNo);
			if(prior!=null){
				throw new RuntimeException("Duplicate member ID '"+id+"' at lines "+prior+" and "+lineNo+" of "+fname);
			}
			raw.add(id);
		}
		if(bf.close()){throw new RuntimeException("ByteFile reported an I/O error reading "+fname);}
		Collections.sort(raw);
		for(String id : raw){
			idToIndex.put(id, indexToId.size());
			indexToId.add(id);
		}
	}

	/**
	 * Streams the convertalis edge list (never retained in memory) and unions every pair
	 * clearing both thresholds. Identity/coverage are decided by exact integer
	 * cross-multiplication against nident/alnlen and the inclusive aligned span over qlen/tlen
	 * -- never by parsing mmseqs' own rounded pident/qcov/tcov percentage columns.
	 * Malformed integer fields throw NumberFormatException unconditionally (no -ea dependency).
	 * An edge ID absent from the member universe throws unconditionally -- likely cross-family
	 * contamination or a stale ID list, never silently dropped.
	 *
	 * @param fname convertalis m8 file, columns: query,target,nident,alnlen,qstart,qend,qlen,tstart,tend,tlen[,...]
	 * @param idRatio {numerator,scale} for the minid percent threshold (see {@link #exactPercentRatio}).
	 * @param qCovRatio {numerator,scale} for the minqcov percent threshold.
	 * @param tCovRatio {numerator,scale} for the mintcov percent threshold.
	 */
	static void processEdges(String fname, HashMap<String,Integer> idToIndex,
			int[] parent, int[] rank, long[] idRatio, long[] qCovRatio, long[] tCovRatio, long[] counters){
		final ByteFile bf=ByteFile.makeByteFile(fname, false);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final String s=new String(line);
			final String[] f=s.split("\t");
			if(f.length<10){
				throw new RuntimeException("Malformed edge row (need query,target,nident,alnlen,qstart,"
					+"qend,qlen,tstart,tend,tlen): "+s);
			}
			final String q=f[0], t=f[1];
			counters[0]++;
			//Validate BOTH IDs against the member universe unconditionally, before any threshold
			//check -- a foreign ID must crash loud whether or not its row happens to pass the
			//threshold (a threshold-gated check would let contamination on a failing row escape).
			final Integer qi=idToIndex.get(q), ti=idToIndex.get(t);
			if(qi==null){throw new RuntimeException("Edge references ID '"+q+"' not present in the ID universe file.");}
			if(ti==null){throw new RuntimeException("Edge references ID '"+t+"' not present in the ID universe file.");}
			//Parse and validate EVERY field for EVERY row -- including self-hits -- before the
			//self-hit skip below. A malformed self-hit row must still crash loud; skipping self
			//rows before field validation would let a corrupt self-row escape undetected, making
			//"zero guard trips" a claim about non-self rows only, not "every row" (Elly's review).
			final long nident=Long.parseLong(f[2]);
			final long alnlen=Long.parseLong(f[3]);
			final long qstart=Long.parseLong(f[4]), qend=Long.parseLong(f[5]), qlen=Long.parseLong(f[6]);
			final long tstart=Long.parseLong(f[7]), tend=Long.parseLong(f[8]), tlen=Long.parseLong(f[9]);
			if(alnlen<=0 || qlen<=0 || tlen<=0){
				throw new RuntimeException("Nonpositive length field in edge row: "+s);
			}
			if(nident<0 || nident>alnlen){
				throw new RuntimeException("nident "+nident+" out of range [0,alnlen="+alnlen+"] in edge row: "+s);
			}
			if(!coordsValid(qstart, qend, qlen) || !coordsValid(tstart, tend, tlen)){
				throw new RuntimeException("Coordinate out of range (1-based inclusive, 1<=start<=end<=len) in edge row: "+s);
			}
			if(q.equals(t)){counters[2]++; continue;}
			final long qSpan=qend-qstart+1, tSpan=tend-tstart+1;
			//percentGE computes nident*100*idRatio.scale >= idRatio.num*alnlen (and the coverage
			//analogs) with every multiplication overflow-checked, never a value precomputed here.
			final boolean idPass = percentGE(nident, alnlen, idRatio[0], idRatio[1]);
			final boolean qCovPass = percentGE(qSpan, qlen, qCovRatio[0], qCovRatio[1]);
			final boolean tCovPass = percentGE(tSpan, tlen, tCovRatio[0], tCovRatio[1]);
			if(!(idPass && qCovPass && tCovPass)){continue;}
			counters[1]++;
			union(parent, rank, qi, ti);
		}
		if(bf.close()){throw new RuntimeException("ByteFile reported an I/O error reading "+fname);}
	}

	/** Incremented every time {@link #percentGE} takes the BigInteger fallback path -- package-private
	 * and observable ONLY so a test can prove the fallback actually executed, not merely that the
	 * returned answer happens to be correct (Elly's review: an earlier "overflow test" used
	 * denominators of 1, so the fast-path multiplication never actually overflowed and silently
	 * gave the right answer via the normal path -- the fallback was never exercised). */
	static long percentGEFallbackCount=0;

	/**
	 * Exact test of num/den &gt;= pctNum/(100*pctScale) -- i.e. whether the fraction num/den clears
	 * a percent threshold expressed as the exact rational pctNum/pctScale (see
	 * {@link #exactPercentRatio}). The fast path computes ALL THREE multiplications
	 * (100*pctScale, num*that, and pctNum*den) via {@link Math#multiplyExact}, so overflow at
	 * ANY step is caught -- an earlier version computed {@code 100L*scale} unchecked before
	 * calling a 2-argument cross-multiply helper, so that specific multiplication could still
	 * silently wrap even though the helper itself was overflow-safe (Elly's review). On overflow,
	 * falls back to {@link BigInteger}, which cannot overflow, and increments
	 * {@link #percentGEFallbackCount} so a test can confirm the fallback path actually ran.
	 */
	static boolean percentGE(long num, long den, long pctNum, long pctScale){
		try{
			final long rhsFactor=Math.multiplyExact(100L, pctScale);
			final long lhs=Math.multiplyExact(num, rhsFactor);
			final long rhs=Math.multiplyExact(pctNum, den);
			return lhs>=rhs;
		}catch(ArithmeticException overflow){
			percentGEFallbackCount++;
			final BigInteger lhs=BigInteger.valueOf(num).multiply(BigInteger.valueOf(100)).multiply(BigInteger.valueOf(pctScale));
			final BigInteger rhs=BigInteger.valueOf(pctNum).multiply(BigInteger.valueOf(den));
			return lhs.compareTo(rhs)>=0;
		}
	}

	/**
	 * True iff (start,end,len) is a valid 1-based inclusive span: 1&lt;=start&lt;=end&lt;=len.
	 * mmseqs convertalis reports coordinates this way -- a full-length 100%-coverage alignment
	 * on a length-L sequence is start=1,end=L,len=L (end==len is normal, not an overflow). An
	 * earlier 0-based-assuming version of this check wrongly rejected every full-coverage row;
	 * found on shortfam's edge list (dense with near-identical full-length members). Fixed
	 * 2026-08-30; see {@link IdentityGroupBuilderTest} for the regression cases.
	 */
	static boolean coordsValid(long start, long end, long len){
		return start>=1 && end>=start && end<=len;
	}

	static int find(int[] parent, int x){
		int r=x;
		while(parent[r]!=r){r=parent[r];}
		while(parent[x]!=r){final int next=parent[x]; parent[x]=r; x=next;}
		return r;
	}

	static void union(int[] parent, int[] rank, int a, int b){
		final int ra=find(parent, a), rb=find(parent, b);
		if(ra==rb){return;}
		if(rank[ra]<rank[rb]){parent[ra]=rb;}
		else if(rank[ra]>rank[rb]){parent[rb]=ra;}
		else{parent[rb]=ra; rank[ra]++;}
	}

	/**
	 * Writes cluster.tsv-format output (rep\tmember per row, exactly one row per input member
	 * including singleton self-rows), sorted canonically by (representative, member) so output
	 * is byte-identical regardless of member-list or edge-list input order. Also writes a
	 * sidecar (`<out>.sidecar.tsv`) with a whole-partition SHA-256 plus per-group
	 * rep/size/member-set-SHA-256, so a downstream consumer can verify integrity without
	 * re-deriving the union-find.
	 *
	 * @return number of distinct groups (components)
	 */
	static int writeClusters(String outFile, boolean overwrite, int[] parent, ArrayList<String> indexToId){
		final int n=indexToId.size();
		final String[] repOfRoot=new String[n];//indexed by root; lazily filled with the min ID seen
		final int[] root=new int[n];
		for(int i=0; i<n; i++){
			final int r=find(parent, i);
			root[i]=r;
			final String id=indexToId.get(i);
			if(repOfRoot[r]==null || id.compareTo(repOfRoot[r])<0){repOfRoot[r]=id;}
		}

		//Group members by final representative, sorted, for both the main output and the sidecar.
		final TreeMap<String,ArrayList<String>> groups=new TreeMap<String,ArrayList<String>>();
		for(int i=0; i<n; i++){
			final String rep=repOfRoot[root[i]];
			ArrayList<String> list=groups.get(rep);
			if(list==null){list=new ArrayList<String>(); groups.put(rep, list);}
			list.add(indexToId.get(i));
		}

		final FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder();
		final MessageDigest wholeDigest=sha256();
		int rowCount=0;
		for(java.util.Map.Entry<String,ArrayList<String>> e : groups.entrySet()){
			final String rep=e.getKey();
			final ArrayList<String> members=e.getValue();
			Collections.sort(members);
			for(String m : members){
				bb.clear();
				bb.append(rep).tab().append(m);
				bsw.println(bb);
				wholeDigest.update(bb.toBytes());
				wholeDigest.update((byte)'\n');
				rowCount++;
			}
		}
		if(bsw.poisonAndWait()){
			throw new RuntimeException("ByteStreamWriter reported an I/O error writing "+outFile);
		}
		if(rowCount!=n){
			throw new RuntimeException("Wrote "+rowCount+" rows but expected exactly "+n+" (one per member).");
		}

		writeSidecar(outFile+".sidecar.tsv", overwrite, groups, wholeDigest);
		return groups.size();
	}

	/** Per-group rep/size/member-set-SHA-256, plus the whole-partition SHA-256, for integrity verification. */
	static void writeSidecar(String sidecarFile, boolean overwrite, TreeMap<String,ArrayList<String>> groups,
			MessageDigest wholeDigest){
		final FileFormat ff=FileFormat.testOutput(sidecarFile, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println(("#whole_partition_sha256\t"+toHex(wholeDigest.digest())).getBytes());
		bsw.println("#rep\tsize\tmembers_sha256".getBytes());
		for(java.util.Map.Entry<String,ArrayList<String>> e : groups.entrySet()){
			final MessageDigest gd=sha256();
			for(String m : e.getValue()){gd.update(m.getBytes()); gd.update((byte)'\n');}
			bsw.println((e.getKey()+"\t"+e.getValue().size()+"\t"+toHex(gd.digest())).getBytes());
		}
		if(bsw.poisonAndWait()){
			throw new RuntimeException("ByteStreamWriter reported an I/O error writing "+sidecarFile);
		}
	}

	static MessageDigest sha256(){
		try{return MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}
	}

	static String toHex(byte[] b){
		final StringBuilder sb=new StringBuilder(b.length*2);
		for(byte x : b){sb.append(String.format("%02x", x));}
		return sb.toString();
	}
}
