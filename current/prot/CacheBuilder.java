package prot;

import java.io.File;
import java.util.HashMap;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
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
		loadHits(accs, repRank);
		emit(accs, archaeaSet, bacteriaSet);
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
					long len=lp.parseLong(4)-lp.parseLong(3)+1;
					if(len<0){len=-len;}
					a.cds++; a.glenSum+=len; a.glenSq+=len*len; a.coding+=(int)len;
				}else if(regionEqualsString(line, typeA, typeLen, "tRNA")){
					a.trna++;
					final int attrLen=lp.length(8);
					final int attrA=lp.a();
					final int code=parseAnticodonCode(line, attrA, attrLen);
					if(code>=0){
						final int cur=a.anticodon.get(code);
						a.anticodon.put(code, (cur<0 ? 1 : cur+1));
					}
				}else if(regionEqualsString(line, typeA, typeLen, "rRNA")){
					final int attrLen=lp.length(8);
					final int attrA=lp.a();
					int subLen=attrLen;
					for(int i=0; i<attrLen; i++){if(line[attrA+i]==','){subLen=i; break;}}
					if(regionEqualsString(line, attrA, subLen, "16S")){a.r16++;}
					else if(regionEqualsString(line, attrA, subLen, "23S")){a.r23++;}
					else if(regionEqualsString(line, attrA, subLen, "5S")){a.r5++;}
					else{a.rother++;}
				}
				lines++;
				if(lines%10000000==0){outstream.println("  gff "+lines/1000000+"M lines, "+accs.size()+" shreds");}
			}
			bf.close();
		}
		outstream.println("GFF: "+lines+" feature lines over "+files.length+" file(s); "+accs.size()+" shreds with genes.");
	}

	/** hits.m8 (blast6, per-gene query IDs, GROUPED BY QUERY) -> per-shred family copies + mapped.
	 *  For each gene (query) take the BEST hit whose target is a top-N rep: bitscore=col11 DESC primary,
	 *  target rep_id lexicographic ASC tiebreak (mmseqs row order is not guaranteed stable across runs,
	 *  so ties need an explicit deterministic break for the cache to regenerate identically). Add ONE
	 *  copy to that gene's shred under the rep's rank, and count the gene as mapped. Query-grouping
	 *  uses the same contiguous-run change-detection as loadGff (mmseqs emits all hits for one query
	 *  consecutively) — only the target rep_id (genuinely different every row) needs a String per hit. */
	void loadHits(HashMap<String, Acc> accs, HashMap<String, Integer> repRank){
		final ByteFile bf=ByteFile.makeByteFile(hitsFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		long lines=0, genesWithHit=0;
		String curGene=null; int bestRank=-1; double bestBits=-1; String bestRep=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			final int terms=lp.terms();
			if(terms<2){continue;}
			final int geneLen=lp.length(0);
			final boolean sameGene=curGene!=null && regionEqualsString(line, lp.a(), geneLen, curGene);
			if(!sameGene){
				if(curGene!=null && bestRank>=0){addGeneHit(accs, curGene, bestRank); genesWithHit++;}
				curGene=lp.parseString(0); bestRank=-1; bestBits=-1; bestRep=null;
			}
			String rep=lp.parseString(1);
			Integer rank=repRank.get(rep);
			if(rank!=null){
				double bits=(terms>11 ? lp.parseDouble(11) : 0);
				if(bits>bestBits || (bits==bestBits && (bestRep==null || rep.compareTo(bestRep)<0))){
					bestBits=bits; bestRank=rank.intValue(); bestRep=rep;
				}
			}
			lines++;
			if(lines%10000000==0){outstream.println("  hits "+lines/1000000+"M rows, "+genesWithHit+" genes mapped");}
		}
		if(curGene!=null && bestRank>=0){addGeneHit(accs, curGene, bestRank); genesWithHit++;}
		bf.close();
		outstream.println("Hits: "+lines+" rows; "+genesWithHit+" genes with a top-"+topN+"-family best hit.");
	}

	/** Strip the _gN gene suffix -> shred name; add one gene-copy to fam[rank] and increment mapped. */
	void addGeneHit(HashMap<String, Acc> accs, String geneId, int rank){
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
	void emit(HashMap<String, Acc> accs, IntHashMap archaeaSet, IntHashMap bacteriaSet){
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
		bsw.print(bb); bb.clear();

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

	private String shredsFile=null, gffArg=null, hitsFile=null, familyFile=null;
	private String archaea4Dir=null, bacteria4Dir=null, outFile=null;
	private int topN=8000;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
