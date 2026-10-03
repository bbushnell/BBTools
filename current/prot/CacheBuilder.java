package prot;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prok.Orf;
import prok.ProkObject;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import tax.TaxTree;
import tracker.KmerTracker;

/**
 * Builds the per-shred MAG-QC feature cache consumed by {@link MagQCVectorMaker#loadCache}.
 * Emits one tab-separated row per shred with the 19-field contract that loadCache reads:
 *
 * <p>contig_id, tid, domain, length, gc, acgt, cds, mapped, glenSum, glenSq, coding,
 * r16, r23, r5, rother, trna, families, anticodon, dimers
 *
 * <p>FIELD SOURCES (three buckets):
 * <ul>
 * <li><b>shred sequence/name</b> (all_shreds.fa): contig_id (first token), tid (parsed from name),
 *     length (bp), gc (G+C count), acgt (A/C/G/T count, excl N), dimers (16 raw dinucleotide
 *     counts, dense CSV, KmerTracker(k=2) native index order 0=AA..15=TT — ADDITIVE components
 *     only; HH/CAGA are ratios and must be derived downstream from SUMMED counts, never averaged
 *     per-contig, same discipline gc/acgt already use for GC).
 * <li><b>archaea4/bacteria4 filename sets</b>: domain. NOT a taxtree lookup — the previous
 *     generic-Java CacheBuilder called {@code tree.getNodeAtLevel(tid, SUPERKINGDOM)}, which threw
 *     an uncaught AssertionError on a tid missing from tree.taxtree.gz (job 25081829 crash-hung 90+
 *     min: main thread died, the ByteStreamWriter writer thread was never poisoned, so the non-daemon
 *     writer thread kept the JVM alive with no progress). archaea4/bacteria4 ARE the literal
 *     source-of-truth genome partition these shreds came from (filenames {@code tid_<N>_...}), so
 *     membership is a zero-taxonomy, zero-crash-surface lookup. One documented exception: tid 29447
 *     (Xanthomonas albilineans plasmid trio, embedded per-contig by multishred with a species-level
 *     tid distinct from the genome's filename tid 380358 — verified directly against shred headers,
 *     not assumed; see records/CACHEBUILDER_PROVENANCE.md). Any OTHER tid in neither set is a real
 *     anomaly and crashes loud via {@link KillSwitch#assertDie(String)} rather than silently
 *     defaulting to a domain.
 * <li><b>callgenes GFF</b> (shred_gff/*.gff.gz, from step 02b): cds (# CDS features), glenSum/glenSq
 *     (Sum and Sum-of-squares of CDS length = end-start+1), coding (Sum CDS bp), r16/r23/r5/rother
 *     (rRNA by attributes[0] subtype), trna (# tRNA features), anticodon (sparse packed-code:count
 *     pairs, from the tRNA attributes' literal {@code anticodon:XXX} tag, mirroring the rRNA
 *     subtype scan — ~31% of tRNA lines carry none; code=(b0&lt;&lt;4)|(b1&lt;&lt;2)|b2, 2 bits/base
 *     via {@link AminoAcid#baseToNumber}, 0-63, absent/ambiguous -&gt; not recorded).
 * <li><b>re-searched hits.m8</b> (per-gene query IDs shredname_gN, step 05 @ --max-seqs 25): families
 *     (sparse "rank:count;..." gene-copies per family, BEST-hit-per-gene mapped to its familylist rank),
 *     mapped (# distinct genes of the shred with a top-N-family best hit).
 * </ul>
 *
 * <p>Because each gene contributes exactly ONE family-copy (its best hit), per shred
 * Sum(families counts) == mapped EXACTLY (the tie-together gate). Every shred is emitted, including
 * zero-gene / zero-hit shreds (empty families field, zero stats), so the row count == the shred count.
 *
 * <p>SINGLE-THREADED (UMP45's ruling): this is a one-time build, not the iteration bottleneck, and
 * MT is the exact axis that crash-hung the previous version — correctness first on this pass.
 *
 * <p>Usage: cachebuilder.sh shreds=all_shreds.fa gff=a.gff.gz,b.gff.gz hits=hits.m8 \
 *          familylist=familylist.tsv archaea4=/path/to/archaea4 bacteria4=/path/to/bacteria4 \
 *          out=percontig_cache.tsv [topn=8000]
 *
 * @author Eru
 */
public class CacheBuilder {

	public static void main(String[] args){
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		Timer t=new Timer();
		CacheBuilder x=new CacheBuilder(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CacheBuilder(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("shreds") || a.equals("in")){shredsFile=b;}
			else if(a.equals("gff") || a.equals("gffs")){gffArg=b;}
			else if(a.equals("hits") || a.equals("m8")){hitsFile=b;}
			else if(a.equals("hitsformat") || a.equals("format")){hitsFormat=b.toLowerCase();}
			else if(a.equals("familylist") || a.equals("family")){familyFile=b;}
			else if(a.equals("archaea4")){archaea4Dir=b;}
			else if(a.equals("bacteria4")){bacteria4Dir=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("topn") || a.equals("top")){topN=Integer.parseInt(b);}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		assert(shredsFile!=null) : "shreds= is required.";
		assert(gffArg!=null) : "gff= is required.";
		assert(hitsFile!=null) : "hits= is required.";
		assert(familyFile!=null) : "familylist= is required.";
		assert(archaea4Dir!=null) : "archaea4= is required.";
		assert(bacteria4Dir!=null) : "bacteria4= is required.";
		assert(outFile!=null) : "out= is required.";
		assert(hitsFormat.equals("blast6") || hitsFormat.equals("genehit")) :
			"hitsformat= must be 'blast6' (default) or 'genehit', got '"+hitsFormat+"'.";
	}

	/** Per-shred accumulator: GFF-derived gene stats + hits-derived family copies + mapped genes.
	 *  fam is rank-&gt;copy-count, primitive int keys/values (never more than a few dozen entries
	 *  per shred, but there are 2.3M shreds, so boxed Integer keys/values here would be a real cost). */
	static final class Acc {
		int cds, coding, r16, r23, r5, rother, trna, mapped;
		long glenSum, glenSq;
		IntHashMap fam=new IntHashMap(4);
		/** Packed-anticodon-code (0-63, 2 bits/base, A=0/C=1/G=2/T=3) -&gt; tRNA count with that
		 *  anticodon. At most 64 distinct keys, so the default capacity is never grown. */
		IntHashMap anticodon=new IntHashMap(4);

		/** Adds one CDS feature using the same inclusive coordinates as GFF. */
		void addCds(final long start, final long stop){
			long len=stop-start+1;
			if(len<0){len=-len;}
			cds++; glenSum+=len; glenSq+=len*len; coding+=(int)len;
		}

		/** Adds one rRNA feature using CacheBuilder's frozen subtype vocabulary. */
		void addRrna(final String subtype){
			if("16S".equals(subtype)){r16++;}
			else if("23S".equals(subtype)){r23++;}
			else if("5S".equals(subtype)){r5++;}
			else{rother++;}
		}

		/** Adds one tRNA feature and, when structurally recoverable, its anticodon. */
		void addTrna(final int anticodonCode){
			trna++;
			if(anticodonCode>=0){
				final int cur=anticodon.get(anticodonCode);
				anticodon.put(anticodonCode, (cur<0 ? 1 : cur+1));
			}
		}

		/** Adds a real in-process CallGenes feature, preserving loadGff's field semantics. */
		void addOrf(final Orf orf){
			if(orf.type==ProkObject.CDS){addCds(orf.start+1L, orf.stop+1L);}
			else if(orf.type==ProkObject.tRNA){addTrna(parseAnticodonCode(orf.trnaAnticodon));}
			else if(orf.type==ProkObject.r16S){addRrna("16S");}
			else if(orf.type==ProkObject.r23S){addRrna("23S");}
			else if(orf.type==ProkObject.r5S){addRrna("5S");}
			else if(orf.type==ProkObject.r18S){addRrna("18S");}
		}

		/** Converts a tRNA anticodon string to the same packed code used by loadGff, or -1. */
		private static int parseAnticodonCode(final String anticodon){
			if(anticodon==null || anticodon.length()!=3){return -1;}
			int code=0;
			for(int i=0; i<3; i++){
				final byte b=(byte)Character.toUpperCase(anticodon.charAt(i));
				final int n=dna.AminoAcid.baseToNumber[b];
				if(n<0){return -1;}
				code=(code<<2)|n;
			}
			return code;
		}
	}

	/** Xanthomonas albilineans plasmid trio (NC_017555/6/7): multishred embeds this per-contig
	 *  species-level tid, distinct from the genome's filename tid (380358, in bacteria4). Verified
	 *  directly against shred headers (not assumed) — see records/CACHEBUILDER_PROVENANCE.md.
	 *  Brian's 08b ruling: plasmids keep their own tid, a documented characteristic not an anomaly. */
	static final int XANTHOMONAS_PLASMID_TID=29447;

	/** Uniform taxid REVISIONS (STATUS.md "Tid Uniqueness"): the genome's FILENAME carries a stale
	 *  tid, but NCBI's current record — which multishred reads per-contig — carries a newer one,
	 *  uniformly across every contig of that genome (unlike the Xanthomonas case, this is NOT a
	 *  split; one tid per genome, just not the filename's). Found by a full shred-tid vs
	 *  archaea4/bacteria4-filename-tid reconciliation on the real cluster data (not assumed) —
	 *  exactly 5 divergent tids total (this table + the plasmid), matching STATUS.md exactly; see
	 *  records/CACHEBUILDER_PROVENANCE.md. All 4 organisms are documented bacterial genera. */
	static final int[] REVISED_TID_BACTERIA={
		2666081, //stale filename tid 2666083, Duganella rivi
		2936273, //stale filename tid 2935863, Alkalimarinus coralli
		3014752, //stale filename tid 2995136, Dellaglioa carnosa
		3025676, //stale filename tid 766894,  Shouchella hunanensis
	};

	static boolean isRevisedBacteriaTid(int tid){
		for(int t : REVISED_TID_BACTERIA){if(t==tid){return true;}}
		return false;
	}

	void process(Timer t){
		HashMap<String, Integer> repRank=loadFamilyList();
		IntHashMap archaeaSet=loadTidSet(archaea4Dir);
		IntHashMap bacteriaSet=loadTidSet(bacteria4Dir);
		outstream.println("Domain tid sets: "+archaeaSet.size()+" archaea, "+bacteriaSet.size()+" bacteria.");
		HashMap<String, Acc> accs=new HashMap<String, Acc>(1<<22);
		loadGff(accs);
		//Snapshot of every seqid the GFF ever named, BEFORE hits can add more keys (hits use real
		//full shred names, never truncated). emit() removes each real shred name it visits; any
		//name still present after the full shreds-file pass never matched a real shred and is a
		//silent GFF-attribution loss (records/V4B_CACHE_TRD_DEFECT_v1.md's exact defect class --
		//e.g. callgenes trd=t collapsing every shred of a contig to one coordinate-free seqid).
		final HashSet<String> gffOnlySeqids=new HashSet<String>(accs.keySet());
		if(hitsFormat.equals("genehit")){loadGeneHits(accs, repRank);}
		else{loadHits(accs, repRank);}
		emit(accs, archaeaSet, bacteriaSet, gffOnlySeqids);
		if(!gffOnlySeqids.isEmpty()){
			final String first=gffOnlySeqids.iterator().next();
			//Plain RuntimeException, not KillSwitch.assertDie: this runs on process()'s single
			//main thread (CacheBuilder is explicitly single-threaded), so there is no
			//worker-thread-hang risk assertDie exists to convert -- an uncaught RuntimeException
			//here already crashes the whole JVM loud, AND (unlike assertDie's hard VM halt) stays
			//testable via a normal try/catch selftest.
			throw new RuntimeException("GFF seqid never matched any real shred name from the "
				+"shreds file ("+gffOnlySeqids.size()+" orphaned seqid(s), first: '"+first+"'). "
				+"Every GFF feature attributed to an orphaned seqid was silently lost from the "
				+"real shred's cache row instead of crashing loud -- this is the exact "
				+"truncated-header defect class documented in records/V4B_CACHE_TRD_DEFECT_v1.md.");
		}
		t.stop();
		outstream.println("Time: \t"+t);
	}

	/** rep_id -> prevalence rank, top-N only (the same 0-based rank space MagQCVectorMaker loads). */
	HashMap<String, Integer> loadFamilyList(){
		HashMap<String, Integer> m=new HashMap<String, Integer>(1<<14);
		final ByteFile bf=ByteFile.makeByteFile(familyFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<2){continue;}
			int rank=lp.parseInt(0);
			if(rank>=topN){continue;}
			String rep=lp.parseString(1);
			m.put(rep, Integer.valueOf(rank));
		}
		bf.close();
		outstream.println("Loaded "+m.size()+" family reps (top "+topN+").");
		return m;
	}

	/** Lists dir/*.fna.gz, parses the tid out of each filename (form tid_&lt;N&gt;_...), and returns
	 *  the membership set. Ground truth: archaea4/bacteria4 ARE the source genome partition these
	 *  shreds were shredded from — zero taxonomy-tree dependency. */
	static IntHashMap loadTidSet(String dir){
		File d=new File(dir);
		File[] files=d.listFiles();
		assert(files!=null) : "Could not list directory: "+dir;
		IntHashMap set=new IntHashMap(Math.max(128, files.length*2));
		int found=0;
		for(File f : files){
			String name=f.getName();
			if(!name.endsWith(".fna.gz")){continue;}
			int tid=TaxTree.parseTaxID(name);
			assert(tid>0) : "Could not parse tid from filename: "+name;
			set.put(tid, 1);
			found++;
		}
		assert(found>0) : "No .fna.gz files with a parseable tid found in "+dir;
		return set;
	}

	/** GFF (comma-separated .gff.gz list) -> per-shred cds/glen/coding/rRNA/trna. seqid=col0, type=col2,
	 *  coords=col3/col4, rRNA subtype = first token of the attributes col8 (ProkObject typeStrings).
	 *  A shred's GFF lines are contiguous (one shred lives in one phylum file, callgenes emits its
	 *  features in a run), so the current-seqid Acc is cached and only re-looked-up on a real change —
	 *  avoids allocating/hashing a String for every one of ~9M+ feature lines, only for the ~2.2M
	 *  shred-with-genes transitions. */
	void loadGff(HashMap<String, Acc> accs){
		String[] files=gffArg.split(",");
		final LineParser1 lp=new LineParser1((byte)'\t');
		long lines=0;
		String curSeqid=null;
		Acc curAcc=null;
		for(String f : files){
			final ByteFile bf=ByteFile.makeByteFile(f, true);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.terms()<9){continue;}
				int seqLen=lp.length(0);
				if(curSeqid==null || !regionEqualsString(line, lp.a(), seqLen, curSeqid)){
					curSeqid=lp.parseString(0);
					curAcc=accs.get(curSeqid);
					if(curAcc==null){curAcc=new Acc(); accs.put(curSeqid, curAcc);}
				}
				final Acc a=curAcc;
				final int typeLen=lp.length(2);
				final int typeA=lp.a();
				if(regionEqualsString(line, typeA, typeLen, "CDS")){
					a.addCds(lp.parseLong(3), lp.parseLong(4));
				}else if(regionEqualsString(line, typeA, typeLen, "tRNA")){
					final int attrLen=lp.length(8);
					final int attrA=lp.a();
					a.addTrna(parseAnticodonCode(line, attrA, attrLen));
				}else if(regionEqualsString(line, typeA, typeLen, "rRNA")){
					final int attrLen=lp.length(8);
					final int attrA=lp.a();
					int subLen=attrLen;
					for(int i=0; i<attrLen; i++){if(line[attrA+i]==','){subLen=i; break;}}
					a.addRrna(new String(line, attrA, subLen));
				}
				lines++;
				if(lines%10000000==0){outstream.println("  gff "+lines/1000000+"M lines, "+accs.size()+" shreds");}
			}
			bf.close();
		}
		outstream.println("GFF: "+lines+" feature lines over "+files.length+" file(s); "+accs.size()+" shreds with genes.");
	}

	/** One gene's running best-hit resolution -- shared by BOTH {@link #loadHits} (blast6) and
	 *  {@link #loadGeneHits} (gene_hits.tsv), so the two input formats structurally cannot diverge in
	 *  their selection/tiebreak behavior (plans/GENE_HIT_CONTRACT_DESIGN_v1.md §5: "one shared
	 *  selection method, two readers", sealed 2026-09-02). */
	static final class BestHitState {
		int rank=-1; double bits=-1; String rep=null;
	}

	/** Applies the FROZEN tiebreak (bitscore DESC primary, rep_id lexicographic ASC) and the
	 *  top-N-family exclusion rule (a rep absent from {@code repRank} is silently excluded, same as no
	 *  hit at all for that candidate) to one (rep,bits) candidate against {@code state}'s current best.
	 *  The legacy blast6 reader uses that exclusion directly. The canonical gene_hits reader first
	 *  requires every rep to occur in CacheBuilder's bound familylist/top-N roster, so an unlisted
	 *  family crashes before reaching this method. CacheBuilder is not role-aware: the upstream
	 *  GeneHitsAssigner gate is what prevents a bait present in a combined family list from being
	 *  emitted. Both valid readers still funnel accepted candidates through this one tiebreak
	 *  implementation. */
	static void considerHit(final BestHitState state, final String rep, final double bits,
			final HashMap<String, Integer> repRank){
		final Integer rank=repRank.get(rep);
		if(rank==null){return;}
		assert(Double.isFinite(bits)) : "Non-finite score "+bits+" for rep '"+rep+"' -- a malformed "
			+"row must never silently pass into the best-hit reduction.";
		if(bits>state.bits || (bits==state.bits && (state.rep==null || rep.compareTo(state.rep)<0))){
			state.bits=bits; state.rank=rank.intValue(); state.rep=rep;
		}
	}

	/** hits.m8 (blast6, per-gene query IDs, GROUPED BY QUERY) -> per-shred family copies + mapped, via
	 *  the shared {@link #considerHit} selection. Query-grouping uses the same contiguous-run
	 *  change-detection as loadGff (mmseqs emits all hits for one query consecutively) — only the
	 *  target rep_id (genuinely different every row) needs a String per hit. This grouped-scalar shape
	 *  is UNCHANGED from before the gene_hits.tsv contract existed (design doc §3: mmseqs honors its
	 *  documented grouping guarantee in practice, so this path is not broken today) -- only the NEW
	 *  {@link #loadGeneHits} path adopts the order-independent accumulator. */
	void loadHits(HashMap<String, Acc> accs, HashMap<String, Integer> repRank){
		final ByteFile bf=ByteFile.makeByteFile(hitsFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		long lines=0, genesWithHit=0;
		String curGene=null;
		BestHitState state=new BestHitState();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			final int terms=lp.terms();
			if(terms<2){continue;}
			final int geneLen=lp.length(0);
			final boolean sameGene=curGene!=null && regionEqualsString(line, lp.a(), geneLen, curGene);
			if(!sameGene){
				if(curGene!=null && state.rank>=0){addGeneHit(accs, curGene, state.rank); genesWithHit++;}
				curGene=lp.parseString(0); state=new BestHitState();
			}
			final String rep=lp.parseString(1);
			final double bits=(terms>11 ? lp.parseDouble(11) : 0);
			considerHit(state, rep, bits, repRank);
			lines++;
			if(lines%10000000==0){outstream.println("  hits "+lines/1000000+"M rows, "+genesWithHit+" genes mapped");}
		}
		if(curGene!=null && state.rank>=0){addGeneHit(accs, curGene, state.rank); genesWithHit++;}
		bf.close();
		outstream.println("Hits: "+lines+" rows; "+genesWithHit+" genes with a top-"+topN+"-family best hit.");
	}

	/** gene_hits.tsv (plans/GENE_HIT_CONTRACT_DESIGN_v1.md, sealed 2026-09-02) -> per-shred family
	 *  copies + mapped, via the SAME {@link #considerHit} selection {@link #loadHits} uses. ORDER-
	 *  INDEPENDENT (design §3): rows may arrive in any order, including interleaved/shuffled queries
	 *  -- one {@link BestHitState} entry per distinct gene persists across the whole file, closed out
	 *  (calling {@link #addGeneHit}) once, at end-of-file, by iterating the map. Presized to the
	 *  caller-measured entry count (never a resize-inducing raw count, design §3b's presizing catch) --
	 *  here sized generically from the file's own line count on a first pass, since a fixed corpus
	 *  constant would silently go stale for any other corpus. */
	void loadGeneHits(HashMap<String, Acc> accs, HashMap<String, Integer> repRank){
		//Upper bound on distinct genes (each row names at most one gene); presizing to this bound is
		//always safe (never undersized) and only wastes a bounded, small excess of empty slots.
		final int rowCountUpperBound=countDataRowsUpperBound(hitsFile);
		final HashMap<String, BestHitState> byGene=new HashMap<String, BestHitState>(
			(int)(rowCountUpperBound/0.75f)+1);
		final ByteFile bf=ByteFile.makeByteFile(hitsFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		boolean sawSchemaVersion=false, sawMethod=false, sawScoreSemantics=false;
		long lines=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String h=new String(line, 1, line.length-1);
				final int tab=h.indexOf('\t');
				final String key=(tab<0 ? h : h.substring(0, tab)).trim();
				final String val=(tab<0 ? "" : h.substring(tab+1)).trim();
				if(key.equals("schema_version")){
					if(sawSchemaVersion){throw new RuntimeException("Duplicate #schema_version header in "+hitsFile);}
					if(!val.equals("1")){throw new RuntimeException(hitsFile+" declares schema_version="+val+", only '1' is defined.");}
					sawSchemaVersion=true;
				}else if(key.equals("method")){
					if(val.isEmpty()){throw new RuntimeException("Empty #method value in "+hitsFile);}
					sawMethod=true;
				}else if(key.equals("score_semantics")){
					if(val.isEmpty()){throw new RuntimeException("Empty #score_semantics value in "+hitsFile);}
					sawScoreSemantics=true;
				}else{
					throw new RuntimeException("Unrecognized metadata line '#"+key+"' in "+hitsFile+" -- "
						+"only #schema_version, #method, #score_semantics are defined.");
				}
				continue;
			}
			if(!(sawSchemaVersion && sawMethod && sawScoreSemantics)){
				throw new RuntimeException("Data row before all 3 required metadata lines "
					+"(#schema_version, #method, #score_semantics) in "+hitsFile+".");
			}
			lp.set(line);
			if(lp.terms()!=3){
				throw new RuntimeException("Malformed gene_hits.tsv row (expected exactly 3 fields): "+new String(line));
			}
			final String gene=lp.parseString(0);
			final String rep=lp.parseString(1);
			//Canonical gene_hits.tsv must contain only IDs in CacheBuilder's bound familylist/top-N
			//roster. Unlike legacy blast6, silently dropping an unlisted rep here would conceal a
			//producer/consumer binding error. This is not a role check: GeneHitsAssigner owns the
			//separate rule that a bait already present in a combined family list is never emitted.
			if(!repRank.containsKey(rep)){
				throw new RuntimeException("gene_hits.tsv rep_id '"+rep+"' for gene '"+gene
					+"' is absent from the bound familylist/top-N roster in "+hitsFile+".");
			}
			//Standard Double.parseDouble, not the zero-allocation Parse.parseDouble: the latter cannot
			//parse literal "NaN"/"Infinity" text (crashes with an unrelated internal AssertionError
			//instead of failing cleanly) -- a real scorer whose own math produced a non-finite value
			//would format it exactly this way (String.valueOf(Double.NaN)=="NaN"), so this format's
			//malformed-score gate case must parse it correctly to reject it with a clear diagnostic,
			//never crash on the parse itself. One small allocation per row is an accepted cost here
			//(a one-time batch validation path, not the iteration bottleneck).
			final String scoreText=lp.parseString(2);
			final double bits;
			try{
				bits=Double.parseDouble(scoreText);
			}catch(final NumberFormatException e){
				throw new RuntimeException("Unparseable score '"+scoreText+"' for gene '"+gene+"' rep '"
					+rep+"' in "+hitsFile+".", e);
			}
			if(!Double.isFinite(bits)){
				throw new RuntimeException("Non-finite score "+bits+" for gene '"+gene+"' rep '"+rep
					+"' in "+hitsFile+" -- a malformed row must never silently pass.");
			}
			BestHitState state=byGene.get(gene);
			if(state==null){state=new BestHitState(); byGene.put(gene, state);}
			considerHit(state, rep, bits, repRank);
			lines++;
			if(lines%10000000==0){outstream.println("  gene_hits "+lines/1000000+"M rows, "+byGene.size()+" distinct genes");}
		}
		bf.close();
		if(!(sawSchemaVersion && sawMethod && sawScoreSemantics)){
			throw new RuntimeException(hitsFile+" is missing one or more required metadata lines "
				+"(#schema_version, #method, #score_semantics).");
		}
		long genesWithHit=0;
		for(final java.util.Map.Entry<String, BestHitState> e : byGene.entrySet()){
			if(e.getValue().rank>=0){addGeneHit(accs, e.getKey(), e.getValue().rank); genesWithHit++;}
		}
		outstream.println("gene_hits.tsv: "+lines+" rows, "+byGene.size()+" distinct genes; "
			+genesWithHit+" genes with a top-"+topN+"-family best hit.");
	}

	/** ONE fast pass counting lines starting with neither '#' nor empty, as a presizing estimate for
	 *  {@link #loadGeneHits}'s accumulator -- an upper bound on distinct genes (each row is at most one
	 *  gene), cheap relative to the real load pass, and never a stale corpus-specific constant. */
	static int countDataRowsUpperBound(final String file){
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		long n=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]!='#'){n++;}
		}
		bf.close();
		assert(n<=Integer.MAX_VALUE) : "Row count "+n+" exceeds Integer.MAX_VALUE -- presizing estimate overflow.";
		return (int)n;
	}

	/** Strip the _gN gene suffix -> shred name; add one gene-copy to fam[rank] and increment mapped. */
	static void addGeneHit(HashMap<String, Acc> accs, String geneId, int rank){
		String shred=stripGene(geneId);
		Acc a=accs.get(shred);
		if(a==null){a=new Acc(); accs.put(shred, a);}
		int cur=a.fam.get(rank);
		a.fam.put(rank, (cur<0 ? 1 : cur+1));
		a.mapped++;
	}

	/** shredname_gN -> shredname (only if the suffix after _g is all digits). */
	static String stripGene(String geneId){
		int i=geneId.lastIndexOf("_g");
		if(i>=0 && i+2<geneId.length()){
			boolean digits=true;
			for(int k=i+2; k<geneId.length(); k++){
				if(!Character.isDigit(geneId.charAt(k))){digits=false; break;}
			}
			if(digits){return geneId.substring(0, i);}
		}
		return geneId;
	}

	/** True if line[a,a+len) byte-equals s, without allocating. */
	static boolean regionEqualsString(byte[] line, int a, int len, String s){
		if(s.length()!=len){return false;}
		for(int i=0; i<len; i++){if(line[a+i]!=(byte)s.charAt(i)){return false;}}
		return true;
	}

	/** Byte form of "anticodon:", scanned for literally within a tRNA feature's attributes column
	 *  (primary-byte-confirmed 2026-08-24: real shred_gff carries a literal {@code ,anticodon:XXX}
	 *  tag at the end of tRNA attributes, XXX a 3-letter DNA code). */
	static final byte[] ANTICODON_TAG={'a','n','t','i','c','o','d','o','n',':'};

	/** Scans attrs[attrA, attrA+attrLen) for a literal "anticodon:" tag and returns the packed
	 *  2-bit-per-base code (0-63) of the 3 letters that follow, or -1 if the tag is absent or any
	 *  of the 3 letters isn't A/C/G/T. ~31% of real tRNA lines carry no tag at all (TrnaCaller's
	 *  structural extraction can fail) — callers must treat -1 as "not recorded", not an error. */
	static int parseAnticodonCode(byte[] line, int attrA, int attrLen){
		final int tagLen=ANTICODON_TAG.length;
		for(int i=0; i+tagLen+3<=attrLen; i++){
			boolean match=true;
			for(int j=0; j<tagLen; j++){
				if(line[attrA+i+j]!=ANTICODON_TAG[j]){match=false; break;}
			}
			if(match){
				final int p=attrA+i+tagLen;
				final int b0=AminoAcid.baseToNumber[line[p]];
				final int b1=AminoAcid.baseToNumber[line[p+1]];
				final int b2=AminoAcid.baseToNumber[line[p+2]];
				if(b0<0 || b1<0 || b2<0){return -1;}
				return (b0<<4)|(b1<<2)|b2;
			}
		}
		return -1;
	}

	/** Streams all_shreds.fa; for each shred computes f0-f5 and the dimers field from the
	 *  sequence/name, resolves domain from the archaea/bacteria tid sets, pulls f6-f17 from the
	 *  accumulators (zeros if absent), and writes the 19-field row. Reports aggregate gates. The
	 *  ByteStreamWriter's writer thread is a
	 *  non-daemon Thread that runs regardless of whether the calling algorithm is threaded — an
	 *  uncaught throw anywhere in this method, single-threaded or not, would otherwise leave it
	 *  un-poisoned and keep the JVM alive doing nothing (the exact class that hung job 25081829 for
	 *  90+ min). The try/finally is the actual structural fix, not thread-count. */
	void emit(HashMap<String, Acc> accs, IntHashMap archaeaSet, IntHashMap bacteriaSet,
			HashSet<String> gffOnlySeqids){
		FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();

		long rows=0, withGenes=0;
		long aR16=0, aR23=0, aR5=0, aRother=0, aTrna=0, aMapped=0, aFamCopies=0, aCds=0;

		try{
			ByteBuilder bb=new ByteBuilder(1<<16);
			bb.append("#contig_id\ttid\tdomain\tlength\tgc\tacgt\tcds\tmapped\tglenSum\tglenSq\tcoding"
				+"\tr16\tr23\tr5\trother\ttrna\tfamilies\tanticodon\tdimers").nl();
			bsw.print(bb); bb.clear();

			final IntList rankBuf=new IntList(64), countBuf=new IntList(64);
			// Reused across every shred (cleared, not reallocated) - additive dimer counts for
			// bin-faithful HH/CAGA (IMPLEMENTATION CORRECTNESS NOTE #1, magqc_rebuild_20260824.plan):
			// GC/ACGT stay as the existing hand-tallied fields above; this is purely for HH/CAGA.
			final KmerTracker dimerTracker=new KmerTracker(2);
			final ByteFile bf=ByteFile.makeByteFile(shredsFile, true);
			String name=null; int length=0, gc=0, acgt=0;
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0 && line[0]=='>'){
					if(name!=null){
						long[] agg=writeRow(bsw, bb, accs, archaeaSet, bacteriaSet, name, length, gc, acgt,
							dimerTracker.counts, rankBuf, countBuf);
						gffOnlySeqids.remove(name);
						rows++; withGenes+=agg[0];
						aCds+=agg[1]; aMapped+=agg[2]; aFamCopies+=agg[3];
						aR16+=agg[4]; aR23+=agg[5]; aR5+=agg[6]; aRother+=agg[7]; aTrna+=agg[8];
						if(rows%500000==0){outstream.println("  emitted "+rows+" rows");}
					}
					// Full pre-tab header, spaces->underscores: matches the 02b-normalized GFF seqid
					// and the hits.m8 query (04 _gN, stripped). all_shreds headers have no tabs, but
					// guard anyway. This is the canonical per-shred join key across all three sources.
					int start=1, end=line.length;
					for(int i=1; i<line.length; i++){if(line[i]=='\t'){end=i; break;}}
					StringBuilder sb=new StringBuilder(end-start);
					for(int i=start; i<end; i++){
						byte ch=line[i];
						sb.append(ch==' ' ? '_' : (char)ch);
					}
					name=sb.toString();
					length=0; gc=0; acgt=0;
					dimerTracker.clearAll();//safe: the prior shred's counts were already read by writeRow above
				}else{
					for(int i=0; i<line.length; i++){
						byte ch=line[i];
						length++;
						switch(ch){
							case 'G': case 'g': case 'C': case 'c': gc++; acgt++; break;
							case 'A': case 'a': case 'T': case 't': acgt++; break;
							default: break;
						}
						dimerTracker.add(ch);
					}
				}
			}
			bf.close();
			if(name!=null){
				long[] agg=writeRow(bsw, bb, accs, archaeaSet, bacteriaSet, name, length, gc, acgt,
					dimerTracker.counts, rankBuf, countBuf);
				gffOnlySeqids.remove(name);
				rows++; withGenes+=agg[0];
				aCds+=agg[1]; aMapped+=agg[2]; aFamCopies+=agg[3];
				aR16+=agg[4]; aR23+=agg[5]; aR5+=agg[6]; aRother+=agg[7]; aTrna+=agg[8];
			}
		}finally{
			bsw.poisonAndWait();
		}

		outstream.println("==== CacheBuilder gates ====");
		outstream.println("rows (== shred count): "+rows);
		outstream.println("shreds with >=1 gene:  "+withGenes);
		outstream.println("CDS total:             "+aCds);
		outstream.println("mapped total:          "+aMapped);
		outstream.println("family-copies total:   "+aFamCopies+"   (must EQUAL mapped total: "+(aFamCopies==aMapped)+")");
		outstream.println("aggregate rRNA:        16S="+aR16+" 23S="+aR23+" 5S="+aR5+" rother="+aRother+"  (PAYOFF: nonzero)");
		outstream.println("aggregate tRNA:        "+aTrna);
	}

	/** Writes one 19-field row; returns aggregate counters
	 *  [withGene, cds, mapped, famCopies, r16, r23, r5, rother, trna]. rankBuf/countBuf are caller-owned
	 *  reused scratch buffers (cleared each call) — the only per-row allocation avoided this way is the
	 *  family list, which can otherwise run to dozens of entries per shred across 2.3M shreds. dimers is
	 *  the CALLER's reused KmerTracker.counts snapshot for exactly this shred (read here synchronously,
	 *  before the caller clears it for the next shred — never retained past this call). */
	long[] writeRow(ByteStreamWriter bsw, ByteBuilder bb, HashMap<String, Acc> accs,
			IntHashMap archaeaSet, IntHashMap bacteriaSet, String name, int length, int gc, int acgt,
			long[] dimers, IntList rankBuf, IntList countBuf){
		final long[] result=appendRow(bb, accs, archaeaSet, bacteriaSet, name, length, gc, acgt,
			dimers, rankBuf, countBuf);
		bsw.print(bb); bb.clear();
		return result;
	}

	/** Appends the exact CacheBuilder 19-field row representation without performing I/O. */
	static long[] appendRow(ByteBuilder bb, HashMap<String, Acc> accs,
			IntHashMap archaeaSet, IntHashMap bacteriaSet, String name, int length, int gc, int acgt,
			long[] dimers, IntList rankBuf, IntList countBuf){
		int tid=TaxTree.parseTaxID(name);
		assert(tid>0) : "Non-positive tid parsed from shred name (corpus should carry only "
			+"valid tids after the source fix): "+name;

		final String domain;
		if(archaeaSet.get(tid)>=0){domain="Archaea";}
		else if(bacteriaSet.get(tid)>=0){domain="Bacteria";}
		else if(tid==XANTHOMONAS_PLASMID_TID){domain="Bacteria";}
		else if(isRevisedBacteriaTid(tid)){domain="Bacteria";}
		else{
			domain=KillSwitch.assertDie("Shred tid "+tid+" ("+name+") is in neither the archaea4 nor "
				+"bacteria4 corpus set, and is not one of the 5 documented exceptions (Xanthomonas "
				+"plasmid tid "+XANTHOMONAS_PLASMID_TID+" or the 4 revised-tid genomes). This is a "
				+"real anomaly, not silently defaulted.");
		}

		Acc a=accs.get(name);
		int cds=0, mapped=0, coding=0, r16=0, r23=0, r5=0, rother=0, trna=0;
		long glenSum=0, glenSq=0, famCopies=0;
		int withGene=0;
		rankBuf.clear(); countBuf.clear();
		if(a!=null){
			cds=a.cds; mapped=a.mapped; coding=a.coding; r16=a.r16; r23=a.r23; r5=a.r5;
			rother=a.rother; trna=a.trna; glenSum=a.glenSum; glenSq=a.glenSq;
			if(cds>0 || a.fam.size()>0){withGene=1;}
			if(a.fam.size()>0){
				final int[] keys=a.fam.keys(), vals=a.fam.values();
				final int invalid=a.fam.invalid();
				for(int i=0; i<keys.length; i++){
					if(keys[i]!=invalid){
						rankBuf.add(keys[i]); countBuf.add(vals[i]); famCopies+=vals[i];
					}
				}
				sortParallel(rankBuf.array, countBuf.array, rankBuf.size());
			}
		}

		if(mapped>0 && cds==0){
			//Plain RuntimeException -- same reasoning as the orphaned-seqid check in process():
			//single main thread, no hang risk, must stay testable.
			throw new RuntimeException("Shred '"+name+"' (tid "+tid+") has mapped="+mapped
				+" genes with a top-N-family hit but cds=0 GFF features -- impossible for a "
				+"correct build (every mapped gene's protein originates from a CDS call on this "
				+"same shred). This is the exact truncated-GFF-seqid defect class documented in "
				+"records/V4B_CACHE_TRD_DEFECT_v1.md, where a truncated seqid silently orphans all "
				+"of a shred's real CDS features under the wrong cache key.");
		}

		bb.append(name).append('\t').append(tid).append('\t').append(domain).append('\t')
			.append(length).append('\t').append(gc).append('\t').append(acgt).append('\t')
			.append(cds).append('\t').append(mapped).append('\t').append(glenSum).append('\t')
			.append(glenSq).append('\t').append(coding).append('\t').append(r16).append('\t')
			.append(r23).append('\t').append(r5).append('\t').append(rother).append('\t').append(trna).append('\t');
		for(int i=0; i<rankBuf.size(); i++){
			if(i>0){bb.append(';');}
			bb.append(rankBuf.array[i]).append(':').append(countBuf.array[i]);
		}

		// anticodon (sparse, packed-code:count;..., sorted by code for reproducible regen)
		bb.append('\t');
		rankBuf.clear(); countBuf.clear();
		if(a!=null && a.anticodon.size()>0){
			final int[] keys=a.anticodon.keys(), vals=a.anticodon.values();
			final int invalid=a.anticodon.invalid();
			for(int i=0; i<keys.length; i++){
				if(keys[i]!=invalid){rankBuf.add(keys[i]); countBuf.add(vals[i]);}
			}
			sortParallel(rankBuf.array, countBuf.array, rankBuf.size());
		}
		for(int i=0; i<rankBuf.size(); i++){
			if(i>0){bb.append(';');}
			bb.append(rankBuf.array[i]).append(':').append(countBuf.array[i]);
		}

		// dimers (dense, 16 raw additive counts, KmerTracker native index order 0=AA..15=TT)
		bb.append('\t');
		for(int i=0; i<16; i++){
			if(i>0){bb.append(',');}
			bb.append(dimers[i]);
		}
		bb.nl();
		return new long[]{withGene, cds, mapped, famCopies, r16, r23, r5, rother, trna};
	}

	/** In-place insertion sort of a[0,n) ascending, permuting b in step (rank/count pairs). Zero
	 *  allocation; n is small (typical per-shred family count is a handful, rarely more than dozens). */
	static void sortParallel(int[] a, int[] b, int n){
		for(int i=1; i<n; i++){
			int ka=a[i], kb=b[i]; int j=i-1;
			while(j>=0 && a[j]>ka){a[j+1]=a[j]; b[j+1]=b[j]; j--;}
			a[j+1]=ka; b[j+1]=kb;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Self-test            ----------------*/
	/*--------------------------------------------------------------*/
	// plans/GENE_HIT_CONTRACT_DESIGN_v1.md's 9 required gate cases, plus the true blast6-vs-gene_hits
	// differential. Reuses Step1AssayDriver's package-visible writeFile/readWholeFile/expectThrows
	// rather than duplicating them (same project convention this whole session).

	static final class Fixture {
		String work, shredsFile, gffFile, familyFile, archaea4Dir, bacteria4Dir;
	}

	/** Minimal but REAL fixture: 2 shreds (one Bacteria tid_100, one Archaea tid_200), an empty GFF
	 *  (no CDS/rRNA/tRNA -- irrelevant to the families/mapped fields under test), and a 2-family
	 *  familylist (repA rank 0, repB rank 1). Shared by every case below; only the hits/gene_hits file
	 *  and its format vary per case. */
	static Fixture buildBaseFixture() throws Exception {
		final String work=System.getProperty("java.io.tmpdir")+"/cachebuilder_selftest_"+System.nanoTime();
		new File(work).mkdirs();
		final Fixture f=new Fixture();
		f.work=work;
		f.shredsFile=work+"/shreds.fa";
		Step1AssayDriver.writeFile(f.shredsFile, ">tid_100_shredA\nACGTACGTACGT\n>tid_200_shredB\nACGTACGTACGT\n");
		f.gffFile=work+"/shreds.gff";
		//One real CDS per shred -- keeps mapped>0 implying cds>0 (assertDie'd in appendRow)
		//consistent with the fixtures' own hits, since real hits always originate from a real CDS.
		Step1AssayDriver.writeFile(f.gffFile, "##gff-version 3\n"
			+"tid_100_shredA\tx\tCDS\t1\t12\t.\t+\t0\tID=1\n"
			+"tid_200_shredB\tx\tCDS\t1\t12\t.\t+\t0\tID=1\n");
		f.familyFile=work+"/familylist.tsv";
		Step1AssayDriver.writeFile(f.familyFile, "#rank\trep_id\tocc_total\n0\trepA\t100\n1\trepB\t50\n");
		f.archaea4Dir=work+"/archaea4"; new File(f.archaea4Dir).mkdirs();
		Step1AssayDriver.writeFile(f.archaea4Dir+"/tid_200_x.fna.gz", "");
		f.bacteria4Dir=work+"/bacteria4"; new File(f.bacteria4Dir).mkdirs();
		Step1AssayDriver.writeFile(f.bacteria4Dir+"/tid_100_x.fna.gz", "");
		return f;
	}

	static CacheBuilder buildCacheBuilder(final Fixture f, final String hitsFile, final String hitsFormat, final String outFile){
		final ArrayList<String> args=new ArrayList<String>();
		args.add("shreds="+f.shredsFile);
		args.add("gff="+f.gffFile);
		args.add("hits="+hitsFile);
		args.add("hitsformat="+hitsFormat);
		args.add("familylist="+f.familyFile);
		args.add("archaea4="+f.archaea4Dir);
		args.add("bacteria4="+f.bacteria4Dir);
		args.add("out="+outFile);
		return new CacheBuilder(args.toArray(new String[0]));
	}

	static String blast6Row(final String query, final String target, final double bits){
		return query+"\t"+target+"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t"+bits+"\n";
	}

	static String geneHitHeader(){
		return "#schema_version\t1\n#method\tA\n#score_semantics\ttest\n";
	}

	static String famString(final HashMap<String,Acc> accs, final String shred){
		final Acc a=accs.get(shred);
		if(a==null){return "";}
		final int[] keys=a.fam.keys(), vals=a.fam.values();
		final int invalid=a.fam.invalid();
		final ArrayList<int[]> pairs=new ArrayList<int[]>();
		for(int i=0; i<keys.length; i++){if(keys[i]!=invalid){pairs.add(new int[]{keys[i], vals[i]});}}
		pairs.sort((x,y)->Integer.compare(x[0], y[0]));
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<pairs.size(); i++){
			if(i>0){sb.append(';');}
			sb.append(pairs.get(i)[0]).append(':').append(pairs.get(i)[1]);
		}
		return sb.toString();
	}

	static int mappedCount(final HashMap<String,Acc> accs, final String shred){
		final Acc a=accs.get(shred);
		return a==null ? 0 : a.mapped;
	}

	static void selftest(){
		try{
			testDifferentialBlast6VsGeneHits();
			System.err.println("  selftest[differential: blast6 vs. converted gene_hits.tsv over an "
				+"equivalent real percontig_cache row -- byte-identical full output]: PASS");

			testEqualScoreTies();
			System.err.println("  selftest[equal-score target ties: lexicographic rep_id ASC wins, "
				+"never first-seen, identical under both readers]: PASS");

			testShuffledInterleavedRows();
			System.err.println("  selftest[shuffled/interleaved gene_hits.tsv rows resolve identically "
				+"to grouped rows -- proves the order-independent accumulator]: PASS");

			testUnknownFamilyRepId();
			System.err.println("  selftest[gene_hits.tsv crashes loud on an unlisted rep_id, while "
				+"the shared legacy blast candidate reducer retains silent out-of-roster exclusion]: PASS");

			testDuplicateRows();
			System.err.println("  selftest[an exact duplicate (gene,rep,score) row is harmless, "
				+"identical to the duplicate removed]: PASS");

			testNonFiniteScoreRejected();
			System.err.println("  selftest[a NaN/Infinity score crashes loud, never silently passes]: PASS");

			testDeterminismAcrossInsertionOrders();
			System.err.println("  selftest[the SAME row multiset in >=2 different insertion orders "
				+"produces byte-identical resolved state]: PASS");

			testMalformedHeaderSubCases();
			System.err.println("  selftest[all 8 malformed-header sub-cases crash loud before any "
				+"data row is processed]: PASS");

			testOrphanedGffSeqidCrashesLoud();
			System.err.println("  selftest[a GFF seqid that never matches any real shred name -- "
				+"the V4B_CACHE_TRD_DEFECT_v1 truncated-header class -- crashes loud, never silently "
				+"drops the feature]: PASS");

			testMappedWithZeroCdsCrashesLoud();
			System.err.println("  selftest[a real shred with mapped>0 genes but cds=0 GFF features "
				+"-- an impossible state for a correct build -- crashes loud]: PASS");

			System.err.println("CacheBuilder SELFTEST PASS (11/11 gate cases; corpus-scale spot check, "
				+"case 9, is a separate real-cluster verification, not a unit test).");
		}catch(final Exception e){
			throw new RuntimeException("CacheBuilder SELFTEST FAILED", e);
		}
	}

	static void testDifferentialBlast6VsGeneHits() throws Exception {
		final Fixture f=buildBaseFixture();
		final String blast6File=f.work+"/hits.m8";
		Step1AssayDriver.writeFile(blast6File,
			blast6Row("tid_100_shredA_g1", "repA", 50.0)
			+blast6Row("tid_100_shredA_g2", "repB", 30.0)
			+blast6Row("tid_200_shredB_g1", "repA", 10.0));
		final String geneHitFile=f.work+"/gene_hits.tsv";
		Step1AssayDriver.writeFile(geneHitFile, geneHitHeader()
			+"tid_100_shredA_g1\trepA\t50.0\n"
			+"tid_100_shredA_g2\trepB\t30.0\n"
			+"tid_200_shredB_g1\trepA\t10.0\n");

		final String outBlast6=f.work+"/out_blast6.tsv";
		buildCacheBuilder(f, blast6File, "blast6", outBlast6).process(new Timer());
		final String outGeneHit=f.work+"/out_genehit.tsv";
		buildCacheBuilder(f, geneHitFile, "genehit", outGeneHit).process(new Timer());

		final String blast6Text=Step1AssayDriver.readWholeFile(outBlast6);
		final String geneHitText=Step1AssayDriver.readWholeFile(outGeneHit);
		if(!blast6Text.equals(geneHitText)){
			throw new RuntimeException("SELFTEST FAILED: blast6 and gene_hits.tsv outputs differ.\n"
				+"blast6:\n"+blast6Text+"\ngenehit:\n"+geneHitText);
		}
		if(!blast6Text.contains("tid_100_shredA\t100\tBacteria\t") || !blast6Text.contains("\t0:1;1:1\t")){
			throw new RuntimeException("SELFTEST FAILED: unexpected output content -- got:\n"+blast6Text);
		}
	}

	static void testEqualScoreTies() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();

		//blast6: repB listed FIRST, repA SECOND, both bits=50 -- repA must still win (lexicographic ASC).
		final String blast6File=f.work+"/ties.m8";
		Step1AssayDriver.writeFile(blast6File, blast6Row("geneX_g1", "repB", 50.0)+blast6Row("geneX_g1", "repA", 50.0));
		final HashMap<String,Acc> accsBlast6=new HashMap<String,Acc>();
		buildCacheBuilder(f, blast6File, "blast6", f.work+"/o1.tsv").loadHits(accsBlast6, repRank);
		if(!famString(accsBlast6, "geneX").equals("0:1")){
			throw new RuntimeException("SELFTEST FAILED: blast6 tie -- expected repA (rank 0) to win, got "+famString(accsBlast6, "geneX"));
		}

		final String geneHitFile=f.work+"/ties_gh.tsv";
		Step1AssayDriver.writeFile(geneHitFile, geneHitHeader()+"geneX_g1\trepB\t50.0\ngeneX_g1\trepA\t50.0\n");
		final HashMap<String,Acc> accsGeneHit=new HashMap<String,Acc>();
		buildCacheBuilder(f, geneHitFile, "genehit", f.work+"/o2.tsv").loadGeneHits(accsGeneHit, repRank);
		if(!famString(accsGeneHit, "geneX").equals("0:1")){
			throw new RuntimeException("SELFTEST FAILED: gene_hits.tsv tie -- expected repA (rank 0) to win, got "+famString(accsGeneHit, "geneX"));
		}
	}

	static void testShuffledInterleavedRows() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();

		//Grouped (contiguous) order: gene1's 2 rows together, then gene2's 2 rows.
		final String grouped=f.work+"/grouped.tsv";
		Step1AssayDriver.writeFile(grouped, geneHitHeader()
			+"gene1_g1\trepB\t10.0\ngene1_g1\trepA\t20.0\ngene2_g1\trepA\t5.0\ngene2_g1\trepB\t15.0\n");
		final HashMap<String,Acc> accsGrouped=new HashMap<String,Acc>();
		buildCacheBuilder(f, grouped, "genehit", f.work+"/o1.tsv").loadGeneHits(accsGrouped, repRank);

		//Shuffled/interleaved: same 4 rows, gene1 and gene2 rows alternate.
		final String shuffled=f.work+"/shuffled.tsv";
		Step1AssayDriver.writeFile(shuffled, geneHitHeader()
			+"gene1_g1\trepB\t10.0\ngene2_g1\trepA\t5.0\ngene1_g1\trepA\t20.0\ngene2_g1\trepB\t15.0\n");
		final HashMap<String,Acc> accsShuffled=new HashMap<String,Acc>();
		buildCacheBuilder(f, shuffled, "genehit", f.work+"/o2.tsv").loadGeneHits(accsShuffled, repRank);

		if(!famString(accsGrouped, "gene1").equals(famString(accsShuffled, "gene1"))
				|| !famString(accsGrouped, "gene2").equals(famString(accsShuffled, "gene2"))
				|| mappedCount(accsGrouped, "gene1")!=mappedCount(accsShuffled, "gene1")
				|| mappedCount(accsGrouped, "gene2")!=mappedCount(accsShuffled, "gene2")){
			throw new RuntimeException("SELFTEST FAILED: shuffled input resolved differently from grouped input -- "
				+"grouped gene1="+famString(accsGrouped, "gene1")+" gene2="+famString(accsGrouped, "gene2")
				+"; shuffled gene1="+famString(accsShuffled, "gene1")+" gene2="+famString(accsShuffled, "gene2"));
		}
		if(!famString(accsGrouped, "gene1").equals("0:1") || !famString(accsGrouped, "gene2").equals("1:1")){
			throw new RuntimeException("SELFTEST FAILED: unexpected baseline resolution -- gene1="
				+famString(accsGrouped, "gene1")+" gene2="+famString(accsGrouped, "gene2"));
		}
	}

	static void testUnknownFamilyRepId() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();

		final String withUnknown=f.work+"/unknown.tsv";
		Step1AssayDriver.writeFile(withUnknown, geneHitHeader()+"gY_g1\trepUnknown\t999.0\ngY_g1\trepA\t5.0\n");
		boolean threw=false;
		try{buildCacheBuilder(f, withUnknown, "genehit", f.work+"/o1.tsv")
			.loadGeneHits(new HashMap<String,Acc>(), repRank);}
		catch(final RuntimeException e){
			threw=e.getMessage()!=null && e.getMessage().contains("repUnknown") && e.getMessage().contains("gY_g1")
				&& e.getMessage().contains("bound familylist/top-N roster");
			if(!threw){throw new RuntimeException("unknown gene_hits rep_id threw without gene/rep/bound-roster context", e);}
		}
		if(!threw){throw new RuntimeException("unknown gene_hits rep_id did not crash loud");}

		//The shared reducer's legacy blast behavior remains unchanged: out-of-roster targets are
		//ordinary search noise and are silently excluded before a later tracked candidate wins.
		final BestHitState legacy=new BestHitState();
		considerHit(legacy, "repUnknown", 999.0, repRank);
		if(legacy.rank!=-1 || legacy.rep!=null){throw new RuntimeException("legacy unknown rep was not excluded");}
		considerHit(legacy, "repA", 5.0, repRank);
		if(legacy.rank!=0 || !"repA".equals(legacy.rep)){throw new RuntimeException("legacy tracked candidate did not win after unknown exclusion");}
	}

	static void testDuplicateRows() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();

		final String withDup=f.work+"/dup.tsv";
		Step1AssayDriver.writeFile(withDup, geneHitHeader()+"gZ_g1\trepA\t7.0\ngZ_g1\trepA\t7.0\n");
		final HashMap<String,Acc> accsA=new HashMap<String,Acc>();
		buildCacheBuilder(f, withDup, "genehit", f.work+"/o1.tsv").loadGeneHits(accsA, repRank);

		final String withoutDup=f.work+"/nodup.tsv";
		Step1AssayDriver.writeFile(withoutDup, geneHitHeader()+"gZ_g1\trepA\t7.0\n");
		final HashMap<String,Acc> accsB=new HashMap<String,Acc>();
		buildCacheBuilder(f, withoutDup, "genehit", f.work+"/o2.tsv").loadGeneHits(accsB, repRank);

		if(!famString(accsA, "gZ").equals(famString(accsB, "gZ")) || !famString(accsA, "gZ").equals("0:1")){
			throw new RuntimeException("SELFTEST FAILED: duplicate row changed the result -- with="+famString(accsA, "gZ")
				+" without="+famString(accsB, "gZ"));
		}
	}

	static void testNonFiniteScoreRejected() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();

		final String nanFile=f.work+"/nan.tsv";
		Step1AssayDriver.writeFile(nanFile, geneHitHeader()+"gN_g1\trepA\tNaN\n");
		Step1AssayDriver.expectThrows(()->buildCacheBuilder(f, nanFile, "genehit", f.work+"/o1.tsv")
			.loadGeneHits(new HashMap<String,Acc>(), repRank), "NaN score");

		final String infFile=f.work+"/inf.tsv";
		Step1AssayDriver.writeFile(infFile, geneHitHeader()+"gN_g1\trepA\tInfinity\n");
		Step1AssayDriver.expectThrows(()->buildCacheBuilder(f, infFile, "genehit", f.work+"/o2.tsv")
			.loadGeneHits(new HashMap<String,Acc>(), repRank), "Infinity score");
	}

	static void testDeterminismAcrossInsertionOrders() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();

		final String forward=f.work+"/forward.tsv";
		Step1AssayDriver.writeFile(forward, geneHitHeader()
			+"g1_g1\trepA\t10.0\ng1_g1\trepB\t20.0\ng2_g1\trepB\t5.0\ng2_g1\trepA\t15.0\ng3_g1\trepA\t8.0\n");
		final HashMap<String,Acc> accsForward=new HashMap<String,Acc>();
		buildCacheBuilder(f, forward, "genehit", f.work+"/o1.tsv").loadGeneHits(accsForward, repRank);

		final String reverse=f.work+"/reverse.tsv";
		Step1AssayDriver.writeFile(reverse, geneHitHeader()
			+"g3_g1\trepA\t8.0\ng2_g1\trepA\t15.0\ng2_g1\trepB\t5.0\ng1_g1\trepB\t20.0\ng1_g1\trepA\t10.0\n");
		final HashMap<String,Acc> accsReverse=new HashMap<String,Acc>();
		buildCacheBuilder(f, reverse, "genehit", f.work+"/o2.tsv").loadGeneHits(accsReverse, repRank);

		for(final String gene : new String[]{"g1", "g2", "g3"}){
			if(!famString(accsForward, gene).equals(famString(accsReverse, gene))
					|| mappedCount(accsForward, gene)!=mappedCount(accsReverse, gene)){
				throw new RuntimeException("SELFTEST FAILED: determinism violated for '"+gene+"' -- forward="
					+famString(accsForward, gene)+" reverse="+famString(accsReverse, gene));
			}
		}
	}

	static void testMalformedHeaderSubCases() throws Exception {
		final Fixture f=buildBaseFixture();
		final HashMap<String,Integer> repRank=buildCacheBuilder(f, f.work+"/dummy.m8", "blast6", f.work+"/dummy.tsv").loadFamilyList();
		final String[] badHeaders={
			"#method\t1\n#score_semantics\ttest\ng_g1\trepA\t1.0\n",                    //missing schema_version
			"#schema_version\t1\n#schema_version\t1\n#method\t1\n#score_semantics\ttest\ng_g1\trepA\t1.0\n", //duplicate schema_version
			"#schema_version\t2\n#method\t1\n#score_semantics\ttest\ng_g1\trepA\t1.0\n", //wrong-version schema_version
			"#schema_version\t1\n#score_semantics\ttest\ng_g1\trepA\t1.0\n",             //missing method
			"#schema_version\t1\n#method\t1\ng_g1\trepA\t1.0\n",                         //missing score_semantics
			"#schema_version\t1\n#method\t1\n#score_semantics\t\ng_g1\trepA\t1.0\n",     //empty score_semantics
			"#schema_version\t1\n#method\t1\n#score_semantics\ttest\n#extra\tstuff\ng_g1\trepA\t1.0\n", //extra metadata line
			"g_g1\trepA\t1.0\n#schema_version\t1\n#method\t1\n#score_semantics\ttest\n", //data-before-header
		};
		final String[] labels={"missing schema_version", "duplicate schema_version", "wrong-version schema_version",
			"missing method", "missing score_semantics", "empty score_semantics", "extra metadata line",
			"data-before-header"};
		for(int i=0; i<badHeaders.length; i++){
			final String file=f.work+"/badheader_"+i+".tsv";
			Step1AssayDriver.writeFile(file, badHeaders[i]);
			final String outFile=f.work+"/badheader_out_"+i+".tsv";
			Step1AssayDriver.expectThrows(()->buildCacheBuilder(f, file, "genehit", outFile)
				.loadGeneHits(new HashMap<String,Acc>(), repRank), labels[i]);
		}
	}

	/** A GFF seqid ('tid_999_ghost') that never matches either real shred name in the base
	 *  fixture must crash the full pipeline loud -- the exact V4B_CACHE_TRD_DEFECT_v1 class,
	 *  where callgenes trd=t collapsed every shred of a contig to one coordinate-free seqid that
	 *  no real shred name (which keeps its coordinates) could ever match. */
	static void testOrphanedGffSeqidCrashesLoud() throws Exception {
		final Fixture f=buildBaseFixture();
		Step1AssayDriver.writeFile(f.gffFile, "##gff-version 3\n"
			+"tid_100_shredA\tx\tCDS\t1\t12\t.\t+\t0\tID=1\n"
			+"tid_200_shredB\tx\tCDS\t1\t12\t.\t+\t0\tID=1\n"
			+"tid_999_ghost\tx\tCDS\t1\t12\t.\t+\t0\tID=1\n");//never a real shred name
		final String emptyHits=f.work+"/empty.m8";
		Step1AssayDriver.writeFile(emptyHits, "");
		Step1AssayDriver.expectThrows(()->buildCacheBuilder(f, emptyHits, "blast6", f.work+"/orphan_out.tsv")
			.process(new Timer()), "orphaned GFF seqid");
	}

	/** A real shred with a gene hit (mapped>0) but zero GFF CDS features must crash the full
	 *  pipeline loud -- this state is impossible in a correct build (every mapped gene's protein
	 *  originates from a CDS call on that same shred), and is exactly how
	 *  V4B_CACHE_TRD_DEFECT_v1's 556,294 affected rows looked before this assertion existed. */
	static void testMappedWithZeroCdsCrashesLoud() throws Exception {
		final Fixture f=buildBaseFixture();
		//Only shredB gets a CDS; shredA gets none, but WILL get a gene hit below.
		Step1AssayDriver.writeFile(f.gffFile, "##gff-version 3\n"
			+"tid_200_shredB\tx\tCDS\t1\t12\t.\t+\t0\tID=1\n");
		final String hits=f.work+"/mapped_zero_cds.m8";
		Step1AssayDriver.writeFile(hits, blast6Row("tid_100_shredA_g1", "repA", 50.0));
		Step1AssayDriver.expectThrows(()->buildCacheBuilder(f, hits, "blast6", f.work+"/mzc_out.tsv")
			.process(new Timer()), "mapped>0 && cds==0");
	}

	private String shredsFile=null, gffArg=null, hitsFile=null, familyFile=null;
	private String archaea4Dir=null, bacteria4Dir=null, outFile=null;
	/** 'blast6' (default, current mmseqs pipeline) or 'genehit' (the method-agnostic gene_hits.tsv
	 *  contract, plans/GENE_HIT_CONTRACT_DESIGN_v1.md, sealed 2026-09-02). */
	private String hitsFormat="blast6";
	private int topN=8000;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
