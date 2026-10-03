package prot;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.LongList;

/**
 * Per-assembly precomputation of the reference-CDS survival rule for every shred of an
 * assembly: for each shred (source span), which original reference genes it FULLY contains.
 *
 * <p>This is not a new definition. {@link ReferenceCdsSurvivalSidecar} evaluates, per bin,
 * "RETAINED iff one selected source span contains the whole gene interval" (its L150-160); that
 * predicate depends only on (gene interval, span), so it is evaluated here ONCE per assembly and
 * the per-bin label operands become a union over the bin's selected spans
 * ({@link ReferenceCdsShredSurvivalTableReader#counts}). Same parsing code (the sidecar's
 * {@code SHRED} pattern, {@code normalizeHeader}, {@code geneId}), same coordinate conventions,
 * same gene universe (protein-coding CDS rows of the original GFF, attribute-less rows keyed by
 * their physical line), same identity-envelope rule for multi-row gene ids.</p>
 *
 * <p>Inputs: {@code gff=} one whole-genome GFF (may be .gz), {@code assembly=} the matching
 * FASTA (contig lengths; every GFF seqid and every shred source key must exist in it),
 * {@code shreds=} the shred id list (a FASTA, whose sequence lengths are validated against the
 * encoded span like the sidecar does, or a text/TSV file whose first tab-field per non-# line is
 * the shred id, e.g. the per-contig cache), {@code out=}, optional {@code tid=} fallback for
 * contigs whose key carries no {@code _tid_N} suffix (fixtures). Shred ids are SOURCE spans
 * ({@code <contig>_tid_<tid>_<start>-<stop>} or the legacy {@code <contig>_<start>-<stop>_tid_<tid>});
 * an {@code __instance_} suffix is rejected because instances are a per-bin notion.</p>
 *
 * <p>Output {@code reference_cds_shred_survival_v1} (TSV): provenance headers (input sha256s,
 * shred-list totals), then one {@code S} row per shred of this assembly in GFF-encounter order
 * followed by shreds of contigs with no CDS row, then one {@code A} row per tid, then an
 * {@code #end} trailer carrying both row counts so truncation is detectable. Gene identities are
 * emitted as per-tid ORDINALS (0-based, in order of first appearance in the GFF; the A row's
 * cds_total is the ordinal count), compressed as {@code lo-hi} ranges. The list may run over
 * several per-phylum tables that all cite the same shred list: the reader checks that the
 * per-file {@code shreds_in_assembly} counts sum to the common {@code shreds_listed}.</p>
 *
 * <p>Shreds of one contig are NOT assumed to tile it: overlapping or nested spans (found in the
 * real corpus, see plans/SUBNET_INPUT_DEPENDENCY_BRIEF_20260909.md section 2.3) simply list the
 * same ordinals twice and the reader's union counts each gene once.</p>
 *
 * @author UMP45, Yoimiya
 */
public final class ReferenceCdsShredSurvivalTable {

	static final String SCHEMA="reference_cds_shred_survival_v1";
	static final String COLUMNS_S="S\tshred_id\ttid\tcontained_count\tpartial_count\tcontained_gene_ordinals";
	static final String COLUMNS_A="A\ttid\tcds_total\tcontigs\tcontigs_with_shreds\tnoshred_max_len\tnoshred_genes";

	private ReferenceCdsShredSurvivalTable(){}

	public static void main(String[] args){
		String gff=null, assembly=null, shreds=null, out=null; int tidFallback=-1; boolean overwrite=false;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new RuntimeException("Arguments must be key=value: "+arg);}
			final String key=arg.substring(0, eq).toLowerCase(), value=arg.substring(eq+1);
			if(key.equals("gff")){gff=value;}
			else if(key.equals("assembly")){assembly=value;}
			else if(key.equals("shreds")){shreds=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("tid")){tidFallback=Integer.parseInt(value);}
			else if(key.equals("overwrite") || key.equals("ow")){overwrite=parse.Parse.parseBoolean(value);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(gff==null || assembly==null || shreds==null || out==null){
			throw new RuntimeException("Usage: java -ea prot.ReferenceCdsShredSurvivalTable gff=<whole-genome.gff[.gz]> assembly=<contigs.fa[.gz]> shreds=<shred id list (fasta or tsv col 0)> out=<table.tsv> [tid=<fallback tid>] [overwrite=f]");
		}
		final Summary s=build(gff, assembly, shreds, out, tidFallback, overwrite);
		System.err.println("ReferenceCdsShredSurvivalTable PASS: tids="+s.tids+" contigs="+s.contigs+" genes="+s.genes
			+" shreds_listed="+s.shredsListed+" shreds_in_assembly="+s.shredsInAssembly+" contained_total="+s.containedTotal
			+" partial_total="+s.partialTotal+" noshred_contigs="+s.noShredContigs+" noshred_max_len="+s.noShredMaxLen
			+" noshred_genes="+s.noShredGenes+" out="+out);
	}

	/** Totals returned to the caller (and printed by main); every field is also in the file. */
	public static final class Summary {
		public int tids, contigs, noShredContigs, noShredMaxLen;
		public long genes, cdsRows, shredsListed, shredsInAssembly, containedTotal, partialTotal, noShredGenes;
	}

	/**
	 * Receives each physical CDS row from the same GFF stream that assigns the
	 * per-tid gene ordinals used by the survival tables. Several CDS rows may
	 * share one ordinal when their identity attributes name one multipart gene.
	 */
	interface GeneContributionSink {
		void observe(String contig, int localGeneIndex, int tid, int geneOrdinal,
				int start0, int end0);
	}

	/** One assembly contig: its length, tid, its shreds, and (while its GFF group is open) its genes. */
	private static final class ContigRec {
		final String key; final int length; final int tid; final int order;
		final IntList shredIdx=new IntList(4);//indices into the shred arrays
		boolean seenInGff=false; int genes=0;
		ContigRec(String key_, int length_, int tid_, int order_){key=key_; length=length_; tid=tid_; order=order_;}
	}

	/** Per-tid accumulators; IntHashMap keyed by tid, zero-boxing. */
	private static final class TidStats {
		final IntHashMap genes=new IntHashMap(1024), contigs=new IntHashMap(1024), withShreds=new IntHashMap(1024),
			noShredMaxLen=new IntHashMap(1024), noShredGenes=new IntHashMap(1024);
		/** IntHashMap.get returns -1 for an absent key (structures/IntHashMap.java:103); all stored values here are >= 0. */
		static int get(IntHashMap m, int key){final int v=m.get(key); return v<0 ? 0 : v;}
		static void inc(IntHashMap m, int key, int by){m.put(key, get(m, key)+by);}
		static void max(IntHashMap m, int key, int v){if(v>get(m, key)){m.put(key, v);}}
	}

	public static Summary build(final String gffPath, final String assemblyPath, final String shredsPath,
			final String outPath, final int tidFallback){
		return build(gffPath, assemblyPath, shredsPath, outPath, tidFallback, false);
	}

	/**
	 * Output discipline (Yoimiya review 2026-09-09): the output may not name an input, may not already exist unless
	 * overwrite is set, and is written to {@code <out>.partial} then atomically renamed only after the writer closed
	 * clean -- so a reader never sees a truncated table under the final name, and a failed run leaves no artifact.
	 */
	public static Summary build(final String gffPath, final String assemblyPath, final String shredsPath,
			final String outPath, final int tidFallback, final boolean overwrite){
		final java.io.File outFile=new java.io.File(outPath);
		for(String in : new String[]{gffPath, assemblyPath, shredsPath}){
			if(samePath(in, outPath)){throw new IllegalArgumentException("out= names an input file: "+outPath);}
		}
		if(outFile.exists() && !overwrite){throw new IllegalArgumentException("out= already exists (pass overwrite=t to replace it): "+outPath);}
		// A unique temporary owned by this run (never a shared, pre-existing name), in the output's directory so the final
		// rename is atomic; on failure only this file is removed and nothing else is touched.
		final java.io.File partial;
		try{
			final java.io.File dir=outFile.getAbsoluteFile().getParentFile();
			if(dir!=null && !dir.isDirectory()){throw new IllegalArgumentException("Output directory does not exist: "+dir);}
			partial=java.io.File.createTempFile(outFile.getName()+".", ".partial", dir);
		}catch(java.io.IOException e){throw new RuntimeException("Could not create a temporary output beside "+outPath, e);}
		final String partialPath=partial.getPath();
		boolean done=false;
		try{
			final Summary s=buildTo(gffPath, assemblyPath, shredsPath, partialPath, tidFallback);
			publishCompleted(partial, outFile, overwrite);
			done=true;
			return s;
		}finally{
			if(!done && partial.exists()){partial.delete();}
		}
	}

	/**
	 * Publishes a completed, owned temporary as {@code dest}. overwrite=t: atomic replace. overwrite=f: atomic CREATE-IF-ABSENT via a
	 * hard link (Files.createLink fails with FileAlreadyExistsException if the destination exists), then the temporary is unlinked --
	 * because Files.move(ATOMIC_MOVE) may replace an existing destination whether or not REPLACE_EXISTING is given (JDK Files.java
	 * L1237-1250: "ATOMIC_MOVE ... other options are ignored ... implementation specific"; Yoimiya review 19:35Z). A destination that
	 * appeared during the run is therefore never clobbered under overwrite=f. Fails closed where hard links are unsupported.
	 */
	static void publishCompleted(final java.io.File temp, final java.io.File dest, final boolean overwrite){
		// Thin context wrapper: the single publication implementation is the root-authored CompletedFilePublisher (Yoimiya 19:41Z).
		try{CompletedFilePublisher.publish(temp.toPath(), dest.toPath(), overwrite);}
		catch(UnsupportedOperationException e){throw new RuntimeException("Filesystem does not support atomic create-if-absent publication (hard links) for "+dest+"; refusing to publish without overwrite=t", e);}
		catch(java.nio.file.FileAlreadyExistsException e){throw new IllegalArgumentException("Destination appeared during the run and overwrite=f: "+dest, e);}
		catch(java.io.IOException e){throw new RuntimeException("Could not publish "+temp+" as "+dest, e);}
	}

	static boolean samePath(final String a, final String b){
		try{
			final java.io.File fa=new java.io.File(a), fb=new java.io.File(b);
			if(fa.getCanonicalPath().equals(fb.getCanonicalPath())){return true;}
			return fa.exists() && fb.exists() && java.nio.file.Files.isSameFile(fa.toPath(), fb.toPath());
		}catch(java.io.IOException e){throw new RuntimeException("Could not compare paths "+a+" and "+b, e);}
	}

	private static Summary buildTo(final String gffPath, final String assemblyPath, final String shredsPath,
			final String outPath, final int tidFallback){
		final Summary sum=new Summary();
		final HashMap<String,ContigRec> contigs=new HashMap<String,ContigRec>(1<<16);
		final ArrayList<ContigRec> contigOrder=new ArrayList<ContigRec>();
		readAssembly(assemblyPath, contigs, contigOrder, tidFallback);
		// Shred arrays (only shreds whose source contig is in THIS assembly are kept).
		final ArrayList<String> shredId=new ArrayList<String>();
		final IntList shredStart=new IntList(), shredEnd=new IntList();
		final ArrayList<ContigRec> shredContig=new ArrayList<ContigRec>();
		sum.shredsListed=readShreds(shredsPath, contigs, shredId, shredStart, shredEnd, shredContig, tidFallback);
		sum.shredsInAssembly=shredId.size();

		final String gffSha=sha256(gffPath), asmSha=sha256(assemblyPath), shredsSha=sha256(shredsPath);

		final ByteStreamWriter bsw=new ByteStreamWriter(outPath, true, false, true); bsw.start();
		final ByteBuilder bb=new ByteBuilder(4096);
		bb.append("#schema_version\t").append(SCHEMA).nl();
		bb.append("#tool\tprot.ReferenceCdsShredSurvivalTable").nl();
		bb.append("#coordinate_convention\tGFF3=1-based-inclusive; shred_suffix=0-based-inclusive; internal=0-based-half-open").nl();
		bb.append("#header_normalization\tfirst-token-before-tab; spaces-to-underscore; no other rewriting").nl();
		bb.append("#gene_universe\tprotein-coding CDS rows from original reference GFF; RNA and fresh calls excluded").nl();
		bb.append("#gene_id_attribute_priority\tgene_id,locus_tag,Parent,ID; multiple CDS rows with one identity form one interval envelope").nl();
		bb.append("#gene_id_fallback\tattribute-less CDS uses its 1-based physical GFF line number as identity; multipart CDS cannot be inferred and remains one row per CDS feature").nl();
		bb.append("#survival_rule\tfull reference interval contained by one selected source span; split pieces never survive").nl();
		bb.append("#gene_ordinal\t0-based per tid in order of first appearance in the GFF; contained_gene_ordinals are lo-hi ranges, comma-joined, - if none; a consumer unions ordinals over its DISTINCT selected spans").nl();
		bb.append("#instance_identity\tshred_id is the source span; __instance_ suffixes are a per-bin notion and are rejected here; repeated instances of one span retain each gene once").nl();
		bb.append("#tiling_assumed\tfalse; overlapping or nested spans list shared ordinals more than once and the reader dedups").nl();
		bb.append("#gff_sha256\t").append(gffPath).append('\t').append(gffSha).nl();
		bb.append("#assembly_sha256\t").append(assemblyPath).append('\t').append(asmSha).nl();
		bb.append("#shreds_sha256\t").append(shredsPath).append('\t').append(shredsSha).nl();
		bb.append("#shreds_listed\t").append(sum.shredsListed).nl();
		bb.append("#shreds_in_assembly\t").append(sum.shredsInAssembly).nl();
		bb.append("#tid_fallback\t").append(tidFallback<0 ? "-" : Integer.toString(tidFallback)).nl();
		bb.append("#columns_S\t").append(COLUMNS_S).nl();
		bb.append("#columns_A\t").append(COLUMNS_A).nl();
		bsw.print(bb); bb.clear();

		final TidStats ts=new TidStats();
		for(ContigRec c : contigOrder){
			TidStats.inc(ts.contigs, c.tid, 1);
			if(c.shredIdx.size()>0){TidStats.inc(ts.withShreds, c.tid, 1);}
			else{sum.noShredContigs++; sum.noShredMaxLen=Math.max(sum.noShredMaxLen, c.length); TidStats.max(ts.noShredMaxLen, c.tid, c.length);}
		}

		final long[] sRows=new long[1];
		long aRows=0; boolean ok=false;
		try{
			readGffAndEmit(gffPath, contigs, shredId, shredStart, shredEnd, shredContig, ts, bsw, bb, sum, sRows, null);
			// Contigs that never appeared in the GFF: their shreds contain nothing.
			final IntList emptyOrdinals=new IntList(0);
			for(ContigRec c : contigOrder){
				if(c.seenInGff){continue;}
				for(int i=0; i<c.shredIdx.size(); i++){
					final int si=c.shredIdx.get(i);
					appendSRow(bb, shredId.get(si), c.tid, 0, 0, emptyOrdinals); sRows[0]++;
					if(bb.length()>=(1<<16)){bsw.print(bb); bb.clear();}
				}
			}
			if(bb.length()>0){bsw.print(bb); bb.clear();}

			// A rows, tids ascending.
			final int[] tids=ts.contigs.toArray(); java.util.Arrays.sort(tids);
			for(int tid : tids){
				bb.append('A').tab().append(tid).tab().append(TidStats.get(ts.genes, tid)).tab().append(TidStats.get(ts.contigs, tid)).tab()
					.append(TidStats.get(ts.withShreds, tid)).tab().append(TidStats.get(ts.noShredMaxLen, tid)).tab().append(TidStats.get(ts.noShredGenes, tid)).nl();
				aRows++;
				if(bb.length()>=(1<<16)){bsw.print(bb); bb.clear();}
			}
			sum.tids=tids.length; sum.contigs=contigOrder.size();
			bb.append("#end\tS=").append(sRows[0]).append("\tA=").append(aRows).append("\tgenes=").append(sum.genes)
				.append("\tcontained_total=").append(sum.containedTotal).append("\tpartial_total=").append(sum.partialTotal)
				.append("\tnoshred_contigs=").append(sum.noShredContigs).append("\tnoshred_max_len=").append(sum.noShredMaxLen)
				.append("\tnoshred_genes=").append(sum.noShredGenes).nl();
			bsw.print(bb); bb.clear();
			ok=true;
		}finally{
			// A rejection thrown mid-stream must not leave the (non-daemon) writer thread alive: the JVM would never exit and the
			// partial file would look like output. Close it here; the success path checks the writer's error flag below.
			if(!ok){bsw.poisonAndWait();}
		}
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+outPath);}
		assert(sRows[0]==sum.shredsInAssembly) : "Emitted "+sRows[0]+" S rows for "+sum.shredsInAssembly+" shreds in the assembly; every"
			+" shred must get exactly one row (GFF-seen contigs then unseen contigs) or the reader's completeness check is void";
		return sum;
	}

	/**
	 * Streams one assembly's genes through the canonical survival-table parser
	 * without constructing shred rows. This is the shared ordinal-assignment seam
	 * for per-gene breakpoint contribution tables; it deliberately does not
	 * duplicate GFF identity or coordinate logic.
	 */
	static Summary scanGenes(final String gffPath, final String assemblyPath,
			final int tidFallback, final GeneContributionSink sink){
		if(sink==null){throw new IllegalArgumentException("Gene contribution sink is null");}
		final Summary sum=new Summary();
		final HashMap<String,ContigRec> contigs=new HashMap<String,ContigRec>(1<<16);
		final ArrayList<ContigRec> contigOrder=new ArrayList<ContigRec>();
		readAssembly(assemblyPath, contigs, contigOrder, tidFallback);
		final TidStats ts=new TidStats();
		for(ContigRec c : contigOrder){TidStats.inc(ts.contigs, c.tid, 1);}
		final ArrayList<String> shredId=new ArrayList<String>();
		final IntList shredStart=new IntList(), shredEnd=new IntList();
		final ArrayList<ContigRec> shredContig=new ArrayList<ContigRec>();
		final long[] sRows=new long[1];
		readGffAndEmit(gffPath, contigs, shredId, shredStart, shredEnd, shredContig,
			ts, null, new ByteBuilder(1), sum, sRows, sink);
		assert(sRows[0]==0) : "Gene-only scan emitted shred rows without a shred input";
		sum.tids=ts.contigs.size();
		sum.contigs=contigOrder.size();
		return sum;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	/** FASTA: header key (normalized exactly as the sidecar does) -> length; tid from the key or the fallback. */
	private static void readAssembly(final String path, final HashMap<String,ContigRec> contigs,
			final ArrayList<ContigRec> order, final int tidFallback){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		String key=null; int len=0; long lineNo=0; boolean ok=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNo++;
				if(line.length>0 && line[0]=='>'){
					if(key!=null){putContig(contigs, order, key, len, tidFallback, path);}
					key=ReferenceCdsSurvivalSidecar.normalizeHeader(new String(line, 1, line.length-1));
					if(key.isEmpty()){throw new IllegalArgumentException("Empty assembly FASTA header at "+path+":"+lineNo);}
					len=0;
				}else if(key!=null){len+=nonWhitespace(line);}
				else if(nonWhitespace(line)>0){throw new IllegalArgumentException("Sequence before header at "+path+":"+lineNo);}
			}
			if(key!=null){putContig(contigs, order, key, len, tidFallback, path);}
			ok=true;
		}finally{
			if(bf.close() && ok){throw new RuntimeException("I/O error reading "+path);}//a malformed-input exception in flight wins
		}
		if(contigs.isEmpty()){throw new IllegalArgumentException("Assembly FASTA has no records: "+path);}
	}

	private static void putContig(HashMap<String,ContigRec> contigs, ArrayList<ContigRec> order, String key, int len, int tidFallback, String path){
		if(len<=0){throw new IllegalArgumentException("Empty assembly sequence for "+key);}
		int tid=tidFromKey(key);
		if(tid<0){tid=tidFallback;}
		if(tid<0){throw new IllegalArgumentException("Assembly contig "+key+" carries no _tid_<n> suffix and no tid= fallback was given ("+path+")");}
		final ContigRec c=new ContigRec(key, len, tid, order.size());
		if(contigs.put(key, c)!=null){throw new IllegalArgumentException("Duplicate assembly FASTA header: "+key);}
		order.add(c);
	}

	/**
	 * The tid carried by a normalized contig key, or -1. Two corpus formats coexist (both appear in the whole-genome FASTA,
	 * its GFF seqids and the cache shred ids, consistently per contig -- checked on Dori 2026-09-09, e.g. Acidobacteriota):
	 * the suffix form {@code <contig>_tid_<digits>} and the prefix form {@code tid|<digits>|<header>}.
	 */
	static int tidFromKey(final String key){
		if(key.startsWith("tid|")){
			final int bar=key.indexOf('|', 4);
			if(bar<=4){return -1;}
			return parseDigits(key, 4, bar);
		}
		final int at=key.lastIndexOf("_tid_");
		if(at<0 || at+5>=key.length()){return -1;}
		return parseDigits(key, at+5, key.length());
	}

	/** Decimal digits of key[from,to) as a tid, -1 if any non-digit; overflow rejected loudly. */
	private static int parseDigits(final String key, final int from, final int to){
		long v=0;
		for(int i=from; i<to; i++){
			final char ch=key.charAt(i);
			if(ch<'0' || ch>'9'){return -1;}
			v=v*10+(ch-'0');
			if(v>Integer.MAX_VALUE){throw new IllegalArgumentException("tid overflow in contig key "+key);}
		}
		return (int)v;
	}

	/**
	 * Shred list: FASTA (headers + validated sequence length) or text (first tab field per non-# line).
	 * Keeps only shreds whose source key is a contig of this assembly; returns the total listed.
	 */
	private static long readShreds(final String path, final HashMap<String,ContigRec> contigs, final ArrayList<String> shredId,
			final IntList shredStart, final IntList shredEnd, final ArrayList<ContigRec> shredContig, final int tidFallback){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		long listed=0, lineNo=0;
		String pendingHeader=null; int pendingLen=0; boolean fasta=false, text=false, ok=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNo++;
				if(line.length==0){continue;}
				if(line[0]=='>'){
					if(text){throw new IllegalArgumentException("Shred list mixes FASTA headers and text ids at "+path+":"+lineNo);}
					fasta=true;
					if(pendingHeader!=null){listed++; addShred(pendingHeader, pendingLen, true, contigs, shredId, shredStart, shredEnd, shredContig, tidFallback, path, lineNo);}
					pendingHeader=ReferenceCdsSurvivalSidecar.normalizeHeader(new String(line, 1, line.length-1)); pendingLen=0;
				}else if(fasta){
					if(pendingHeader==null){throw new IllegalArgumentException("Sequence before shred header at "+path+":"+lineNo);}
					pendingLen+=nonWhitespace(line);
				}else{
					if(line[0]=='#'){continue;}
					text=true;
					int tab=0; while(tab<line.length && line[tab]!='\t'){tab++;}
					final String id=new String(line, 0, tab).trim();
					if(id.isEmpty()){throw new IllegalArgumentException("Empty shred id at "+path+":"+lineNo);}
					listed++; addShred(id, -1, false, contigs, shredId, shredStart, shredEnd, shredContig, tidFallback, path, lineNo);
				}
			}
			if(pendingHeader!=null){listed++; addShred(pendingHeader, pendingLen, true, contigs, shredId, shredStart, shredEnd, shredContig, tidFallback, path, lineNo);}
			ok=true;
		}finally{
			if(bf.close() && ok){throw new RuntimeException("I/O error reading "+path);}
		}
		if(listed==0){throw new IllegalArgumentException("Shred list is empty: "+path);}
		return listed;
	}

	private static void addShred(final String header, final int seqLen, final boolean validateLen, final HashMap<String,ContigRec> contigs,
			final ArrayList<String> shredId, final IntList shredStart, final IntList shredEnd, final ArrayList<ContigRec> shredContig,
			final int tidFallback, final String path, final long lineNo){
		if(ReferenceCdsSurvivalSidecar.INSTANCE_SUFFIX.matcher(header).matches()){
			throw new IllegalArgumentException("Shred list must hold source spans, not instances (__instance_ suffix) at "+path+":"+lineNo+": "+header);
		}
		final Matcher m=ReferenceCdsSurvivalSidecar.SHRED.matcher(header);
		if(!m.matches()){throw new IllegalArgumentException("Cannot parse shred source span header at "+path+":"+lineNo+": "+header);}
		final String source=ReferenceCdsSurvivalSidecar.normalizeHeader(m.group(1));
		final int start=Integer.parseInt(m.group(2)), stop=Integer.parseInt(m.group(3));
		if(stop<start){throw new IllegalArgumentException("Shred stop before start: "+header);}
		final int end=stop+1;
		if(validateLen && seqLen!=end-start){throw new IllegalArgumentException("Shred sequence length "+seqLen+" != encoded span length "+(end-start)+" at "+path+":"+lineNo+": "+header);}
		final ContigRec c=contigs.get(source);
		if(c==null){return;}//belongs to another assembly (per-phylum tables); the reader sums shreds_in_assembly over files
		if(end>c.length){throw new IllegalArgumentException("Shred span outside source contig "+source+": "+start+"-"+stop+" length="+c.length+" ("+header+")");}
		// The tid encoded on the shred (legacy trailing _tid_N after the span, else the source key's own prefix/suffix) must agree with the contig's tid.
		int shredTid=-1;
		final int tidAt=header.lastIndexOf("_tid_");
		if(tidAt>=0 && tidAt>=source.length()){shredTid=parseDigits(header, tidAt+5, header.length());}
		else{shredTid=tidFromKey(source);}
		if(shredTid<0){shredTid=tidFallback;}
		if(shredTid!=c.tid){throw new IllegalArgumentException("Shred "+header+" encodes tid "+shredTid+" but its contig "+source+" has tid "+c.tid);}
		for(int i=0; i<c.shredIdx.size(); i++){
			final int si=c.shredIdx.get(i);
			if(shredStart.get(si)==start && shredEnd.get(si)==end){throw new IllegalArgumentException("Duplicate shred source span at "+path+":"+lineNo+": "+header+" (already "+shredId.get(si)+")");}
		}
		c.shredIdx.add(shredId.size());
		shredId.add(header); shredStart.add(start); shredEnd.add(end); shredContig.add(c);
	}

	/*--------------------------------------------------------------*/
	/*----------------      GFF stream + emission    ----------------*/
	/*--------------------------------------------------------------*/

	/** Streams the GFF grouped by seqid; on each group end, evaluates containment for that contig's shreds and emits its S rows. */
	private static void readGffAndEmit(final String path, final HashMap<String,ContigRec> contigs, final ArrayList<String> shredId,
			final IntList shredStart, final IntList shredEnd, final ArrayList<ContigRec> shredContig, final TidStats ts,
			final ByteStreamWriter bsw, final ByteBuilder bb, final Summary sum, final long[] sRows,
			final GeneContributionSink geneSink){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		// Open contig group state.
		ContigRec cur=null; byte[] curKeyBytes=null; int curKeyLen=0;
		final IntList gStart=new IntList(1024), gEnd=new IntList(1024), gOrd=new IntList(1024);
		final HashMap<String,Integer> idToLocal=new HashMap<String,Integer>();//only touched for rows that carry a gene identity attribute
		final LongList sortBuf=new LongList(1024); final IntList prefMax=new IntList(1024), ordBuf=new IntList(256);
		long lineNo=0; int localGeneIndex=0; boolean ok=false;
		try{
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			lineNo++;
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<9){throw new IllegalArgumentException("GFF line "+lineNo+" has fewer than 9 fields: "+path);}
			if(!lp.termEquals("CDS", 2)){continue;}
			// Group boundary: compare the raw seqid bytes with the open group's raw bytes.
			lp.setBounds(0); final int a0=lp.a(), b0=lp.b();//a()=start, b()=end (exclusive) of term 0
			final int kl=b0-a0;
			boolean same=(cur!=null && kl==curKeyLen);
			if(same){for(int i=0; i<kl; i++){if(line[a0+i]!=curKeyBytes[i]){same=false; break;}}}
			if(!same){
				if(cur!=null){finishContig(cur, gStart, gEnd, gOrd, shredId, shredStart, shredEnd, ts, bsw, bb, sortBuf, prefMax, ordBuf, sum, sRows);}
				final String key=ReferenceCdsSurvivalSidecar.normalizeHeader(new String(line, a0, kl));
				cur=contigs.get(key);
				if(cur==null){throw new IllegalArgumentException("GFF CDS seqid absent from assembly FASTA: "+key+" line "+lineNo+" ("+path+")");}
				if(cur.seenInGff){throw new IllegalArgumentException("GFF is not grouped by seqid: "+key+" reappears at line "+lineNo+" ("+path+"); sort the GFF by seqid first");}
				cur.seenInGff=true;
				curKeyBytes=new byte[kl]; System.arraycopy(line, a0, curKeyBytes, 0, kl); curKeyLen=kl;
				gStart.clear(); gEnd.clear(); gOrd.clear(); idToLocal.clear();
				localGeneIndex=0;
			}
			final int start=lp.parseInt(3), end=lp.parseInt(4);
			if(start<1 || end<start || end>cur.length){throw new IllegalArgumentException("GFF CDS interval outside assembly at line "+lineNo+": "+start+"-"+end+" length="+cur.length+" ("+cur.key+")");}
			final int sl=lp.setBounds(6); final int sa=lp.a();
			if(sl!=1 || !(line[sa]=='+' || line[sa]=='-' || line[sa]=='.')){throw new IllegalArgumentException("Unsupported GFF strand at line "+lineNo+": "+new String(line, sa, Math.max(sl,0)));}
			final char strand=(char)line[sa];
			// Gene identity exactly as the sidecar: attribute priority gene_id,locus_tag,Parent,ID; else the physical line.
			String id=null;
			final int al=lp.setBounds(8); final int aa=lp.a();
			if(al>0 && !(al==1 && line[aa]=='.')){id=ReferenceCdsSurvivalSidecar.geneId(new String(line, aa, al));}
			final int geneOrdinal;
			if(id==null){
				// CDS_row_<line>: unique by construction; no map lookup needed.
				geneOrdinal=TidStats.get(ts.genes, cur.tid);
				gStart.add(start-1); gEnd.add(end); gOrd.add(geneOrdinal); TidStats.inc(ts.genes, cur.tid, 1);
			}else{
				final String k=id+"|"+strand;
				final Integer local=idToLocal.get(k);
				if(local==null){
					idToLocal.put(k, gStart.size());
					geneOrdinal=TidStats.get(ts.genes, cur.tid);
					gStart.add(start-1); gEnd.add(end); gOrd.add(geneOrdinal); TidStats.inc(ts.genes, cur.tid, 1);
				}else{
					final int li=local.intValue();
					geneOrdinal=gOrd.get(li);
					gStart.set(li, Math.min(gStart.get(li), start-1)); gEnd.set(li, Math.max(gEnd.get(li), end));
				}
			}
			if(geneSink!=null){geneSink.observe(cur.key, localGeneIndex, cur.tid, geneOrdinal, start-1, end);}
			localGeneIndex++;
			sum.cdsRows++;
		}
		if(cur!=null){finishContig(cur, gStart, gEnd, gOrd, shredId, shredStart, shredEnd, ts, bsw, bb, sortBuf, prefMax, ordBuf, sum, sRows);}
		ok=true;
		}finally{
			if(bf.close() && ok){throw new RuntimeException("I/O error reading "+path);}
		}
	}

	/** Containment for one closed contig group: sorts genes by start, evaluates every shred of the contig, emits its S rows. */
	private static void finishContig(final ContigRec c, final IntList gStart, final IntList gEnd, final IntList gOrd,
			final ArrayList<String> shredId, final IntList shredStart, final IntList shredEnd, final TidStats ts,
			final ByteStreamWriter bsw, final ByteBuilder bb, final LongList sortBuf, final IntList prefMax, final IntList ordBuf,
			final Summary sum, final long[] sRows){
		final int n=gStart.size();
		c.genes=n; sum.genes+=n;
		if(c.shredIdx.size()==0){sum.noShredGenes+=n; TidStats.inc(ts.noShredGenes, c.tid, n); return;}
		// Order genes by start (ties by local index) without boxing: pack (start<<32 | local).
		sortBuf.clear();
		for(int i=0; i<n; i++){sortBuf.add((((long)gStart.get(i))<<32) | i);}
		sortBuf.sort();
		final int[] order=new int[n]; final int[] sStart=new int[n], sEnd=new int[n], sOrd=new int[n];
		for(int i=0; i<n; i++){final int li=(int)(sortBuf.get(i) & 0xFFFFFFFFL); order[i]=li; sStart[i]=gStart.get(li); sEnd[i]=gEnd.get(li); sOrd[i]=gOrd.get(li);}
		prefMax.clear();
		for(int i=0, m=0; i<n; i++){m=Math.max(m, sEnd[i]); prefMax.add(m);}
		for(int k=0; k<c.shredIdx.size(); k++){
			final int si=c.shredIdx.get(k);
			final int s=shredStart.get(si), e=shredEnd.get(si);
			// First gene with start >= s.
			int lo=0, hi=n;
			while(lo<hi){final int mid=(lo+hi)>>>1; if(sStart[mid]<s){lo=mid+1;}else{hi=mid;}}
			ordBuf.clear(); int partial=0;
			for(int i=lo; i<n && sStart[i]<e; i++){
				if(sEnd[i]<=e){ordBuf.add(sOrd[i]);}
				else{partial++;}//starts inside, ends beyond: overlaps but is not contained
			}
			// Genes starting before s that still reach into the span (partial): bounded by the prefix max end.
			for(int i=lo-1; i>=0 && prefMax.get(i)>s; i--){if(sEnd[i]>s){partial++;}}
			ordBuf.sort();
			assert(ordBuf.unique()) : "Duplicate ordinal within one shred "+shredId.get(si)+": gene ordinals are unique per gene by construction (gOrd assigned once per identity)";
			sum.containedTotal+=ordBuf.size(); sum.partialTotal+=partial;
			appendSRow(bb, shredId.get(si), c.tid, ordBuf.size(), partial, ordBuf); sRows[0]++;
			if(bb.length()>=(1<<16)){bsw.print(bb); bb.clear();}
		}
	}

	/** S row with ordinals compressed to lo-hi ranges (ordBuf must be sorted, unique). */
	static void appendSRow(final ByteBuilder bb, final String id, final int tid, final int contained, final int partial, final IntList ordBuf){
		assert(contained==ordBuf.size()) : "contained_count "+contained+" != ordinal count "+ordBuf.size()+" for "+id;
		bb.append('S').tab().append(id).tab().append(tid).tab().append(contained).tab().append(partial).tab();
		if(ordBuf.size()==0){bb.append('-');}
		else{
			int lo=ordBuf.get(0), prev=lo; boolean first=true;
			for(int i=1; i<=ordBuf.size(); i++){
				final int v=(i<ordBuf.size() ? ordBuf.get(i) : Integer.MIN_VALUE);
				if(i<ordBuf.size() && v==prev+1){prev=v; continue;}
				if(!first){bb.append(',');} first=false;
				bb.append(lo); if(prev!=lo){bb.append('-').append(prev);}
				lo=v; prev=v;
			}
		}
		bb.nl();
	}

	/** The sidecar's file hash (raw bytes, so a .gz input hashes as its compressed bytes), unchecked wrapper. */
	private static String sha256(final String path){
		try{return ReferenceCdsSurvivalSidecar.sha256(path);}
		catch(java.io.IOException e){throw new RuntimeException("Could not hash "+path, e);}
	}

	private static int nonWhitespace(final byte[] line){
		int n=0; for(int i=0; i<line.length; i++){if(line[i]>' '){n++;}} return n;
	}
}
