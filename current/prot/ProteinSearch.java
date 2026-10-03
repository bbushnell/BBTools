package prot;

import java.util.ArrayList;
import java.util.List;

import java.io.FileInputStream;
import java.io.InputStream;
import java.security.DigestInputStream;
import java.security.MessageDigest;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Command-line front end for {@link ProteinSearcher}: a blastp-style protein
 * search that reads a query protein FASTA and a database protein FASTA and
 * writes BLAST-tab {@code outfmt 6} TSV results.
 *
 * <p>This is the thin CLI/testing wrapper around the in-memory API; the actual
 * search logic lives in {@link ProteinSearcher} so it can be called directly
 * from other BBTools code without any file I/O.</p>
 *
 * <p>Usage: {@code proteinsearch.sh query=<q.faa> db=<db.faa> out=<hits.tsv>}</p>
 *
 * @author Eru
 */
public final class ProteinSearch {

	/**
	 * Program entry point: parses arguments, loads both FASTA inputs, runs the
	 * search, writes the TSV plus a sidecar manifest.
	 * @param args Command-line arguments (flag=value).
	 */
	public static void main(String[] args){
		final Timer t=new Timer();
		if(args.length==0){printUsageAndExit();}

		String queryFile=null, dbFile=null, out=null, queryLoader="generic";
		String sidecarFile=null, kmerSetsFile=null;
		boolean overwrite=true;
		final ProteinSearcher searcher=new ProteinSearcher();
		FamilySearchResources.applyFrozenV1(searcher);

		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("query") || a.equals("in") || a.equals("q")){queryFile=b;}
			else if(a.equals("db") || a.equals("ref") || a.equals("database") || a.equals("d")){dbFile=b;}
			else if(a.equals("out") || a.equals("o")){out=b;}
			else if(a.equals("k")){searcher.k=Integer.parseInt(b);}
			else if(a.equals("minseedhits")){searcher.minSeedHits=Integer.parseInt(b);}
			else if(a.equals("evalue") || a.equals("e")){searcher.evalueCutoff=Double.parseDouble(b);}
			else if(a.equals("minid") || a.equals("minpident")){searcher.minPident=Double.parseDouble(b);}
			else if(a.equals("minscore")){searcher.minRawScore=Integer.parseInt(b);}
			else if(a.equals("mincoverage") || a.equals("mincov")){searcher.minCoverage=Double.parseDouble(b);}
			else if(a.equals("threads") || a.equals("t")){searcher.threads=Integer.parseInt(b);}
			else if(a.equals("maxtargetseqs") || a.equals("mts")){searcher.maxTargetSeqs=Integer.parseInt(b);}
			else if(a.equals("reducedseed") || a.equals("reduced")){searcher.reducedSeed=Parse.parseBoolean(b);}
			else if(a.equals("queryloader") || a.equals("qloader")){queryLoader=(b==null ? null : b.toLowerCase());}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("shortlist")){searcher.shortlist=Integer.parseInt(b);}
			else if(a.equals("sidecar")){sidecarFile=b;}
			else if(a.equals("kmersets")){kmerSetsFile=b;}
			else if(a.equals("aligner")){searcher.aligner=b.toLowerCase();}
			else if(a.equals("triage")){searcher.triage=Parse.parseBoolean(b);}
			else if(a.equals("kmernorm")){searcher.kmerNorm=b.toLowerCase();}
			else if(a.equals("covdef")){searcher.covDef=b.toLowerCase();}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}

		if(queryFile==null || dbFile==null){
			System.err.println("Error: both query= and db= are required.\n");
			printUsageAndExit();
		}
		if(!"generic".equals(queryLoader) && !"callgenes".equals(queryLoader)){
			throw new RuntimeException("queryloader= must be 'generic' or 'callgenes', got: "+queryLoader);
		}
		if(searcher.shortlist>0){
			if(sidecarFile==null || kmerSetsFile==null){
				throw new RuntimeException("shortlist>0 requires both sidecar= and kmersets= "+
					"(FamilyShortlistSidecar.load needs the live covering-set file to cross-verify "+
					"against the sidecar's recorded hash, even though the sidecar already bundles "+
					"the k-mers themselves -- catches a stale/rebuilt covering-set file).");
			}
			System.err.println("Loading shortlist sidecar "+sidecarFile+"...");
			//A3 (2026-09-04, plans/PER_FAMILY_THRESHOLDS_v1.md sec 0): map the lowercase flag
			//value onto the sidecar's recorded capitalized aligner name -- UMP45's own words
			//(results/family_triage_gate_v1.md): "map flag aligner=aaaligner -> expectedAligner
			//AAAligner". null for glocal/blosum: neither is sidecar-bound yet (Eru's aligners are
			//still a placeholder in ProteinSearcher.alignWith()), so the aligner-mismatch check is
			//skipped for them rather than spuriously refusing every real aaaligner run.
			final String expectedAligner=(searcher.aligner.equals("aaaligner") ? "AAAligner" : null);
			searcher.sidecar=FamilyShortlistSidecar.load(sidecarFile, dbFile, kmerSetsFile, expectedAligner);
			System.err.println("Sidecar loaded: "+searcher.sidecar.nFamilies+" families, "+
				searcher.sidecar.dims+" composition dims.");
		}

		//"callgenes" mode routes through CallGenesProteinLoader.load(), which mirrors
		//04_fix_headers_v3.sh's per-gene _gN unique-ification -- required because CallGenes'
		//outa= output puts every gene on one contig under an IDENTICAL header, which the generic
		//readFasta() (first-whitespace-token ID) collapses onto one id, correctly tripped by
		//ProteinSearcher.checkDuplicateIds() (real failure, documented in CallGenesProteinLoader's
		//own javadoc, found 2026-08-29 running this on real CallGenes output: 4117 genes all
		//shared one id).
		final List<ProteinSequence> queries="callgenes".equals(queryLoader) ?
			CallGenesProteinLoader.loadValidated(queryFile) : readFasta(queryFile);
		final List<ProteinSequence> targets=readFasta(dbFile);
		System.err.println("Loaded "+queries.size()+" queries ("+queryLoader+" loader) and "+
			targets.size()+" database sequences.");

		final List<ProteinHit> hits=searcher.search(queries, targets);

		writeResults(hits, out, overwrite);
		writeSidecar(out, overwrite, queryFile, dbFile, queryLoader, searcher, queries.size(),
			targets.size(), hits.size(), sidecarFile, kmerSetsFile);

		t.stop();
		System.err.println("Wrote "+hits.size()+" hits in "+t);
	}

	/**
	 * Reads a protein FASTA into validated in-memory sequences. The identifier is
	 * the first whitespace-delimited token of each header (frozen contract).
	 * @param fname FASTA filename.
	 * @return List of protein sequences.
	 */
	static List<ProteinSequence> readFasta(final String fname){
		final ArrayList<ProteinSequence> list=new ArrayList<ProteinSequence>();
		final ByteFile bf=ByteFile.makeByteFile(fname, false);
		String id=null;
		int skippedStops=0;
		final ByteBuilder seq=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				if(id!=null && addOrSkipMalformedStop(list, id, seq.toBytes())){skippedStops++;}
				id=parseHeaderId(line);
				seq.clear();
			}else{
				seq.append(line);
			}
		}
		if(id!=null && addOrSkipMalformedStop(list, id, seq.toBytes())){skippedStops++;}
		bf.close();
		if(skippedStops>0){
			System.err.println("WARNING: skipped "+skippedStops+" FASTA record(s) with a leading/internal '*' stop marker from "+fname+".");
		}
		if(list.isEmpty()){throw new RuntimeException("No sequences found in "+fname);}
		return list;
	}

	/** Adds one strict sequence, except that malformed nonterminal-stop records are explicitly skipped. */
	static boolean addOrSkipMalformedStop(final List<ProteinSequence> list, final String id, final byte[] raw){
		if(ProteinSequence.hasNonterminalStopMarker(raw)){return true;}
		list.add(new ProteinSequence(id, raw));
		return false;
	}

	/**
	 * Streaming SHA-256 of a file, for exact reproducible-metadata recording in the .meta
	 * sidecar (added 2026-09-02, Elly's assignment) -- lets two runs on the same inputs be
	 * verified identical by input content, not just by path string (a path can be reused with
	 * different bytes).
	 * @param fname File to hash.
	 * @return Lowercase hex SHA-256 digest.
	 */
	static String sha256(final String fname){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			final byte[] buf=new byte[65536];
			try(InputStream in=new DigestInputStream(new FileInputStream(fname), md)){
				while(in.read(buf)>=0){}
			}
			final byte[] digest=md.digest();
			final StringBuilder sb=new StringBuilder(digest.length*2);
			for(final byte b : digest){
				final String h=Integer.toHexString(b&0xff);
				if(h.length()<2){sb.append('0');}
				sb.append(h);
			}
			return sb.toString();
		}catch(Exception e){
			throw new RuntimeException("Failed to hash "+fname, e);
		}
	}

	/** Extracts the first whitespace-delimited token from a FASTA header line. */
	static String parseHeaderId(final byte[] header){
		int start=1;//skip '>'
		while(start<header.length && (header[start]==' ' || header[start]=='\t')){start++;}
		int stop=start;
		while(stop<header.length && header[stop]!=' ' && header[stop]!='\t'){stop++;}
		if(stop<=start){throw new RuntimeException("Empty FASTA identifier in header line.");}
		return new String(header, start, stop-start);
	}

	/** Writes the outfmt6 TSV to a file or stdout. */
	static void writeResults(final List<ProteinHit> hits, final String out, final boolean overwrite){
		if(out==null || out.equalsIgnoreCase("stdout")){
			for(ProteinHit h : hits){System.out.println(h.toTsv());}
			return;
		}
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		for(ProteinHit h : hits){bsw.println(h.toTsv());}
		bsw.poisonAndWait();
	}

	/**
	 * Writes a minimal run-manifest sidecar recording the parameters and counts.
	 * A full immutable manifest (checksums, versioned lambda/K, DB build ID) is
	 * deferred; this records what the MVP actually used, including that the
	 * E-value is approximate (edge correction omitted).
	 */
	static void writeSidecar(final String out, final boolean overwrite, final String queryFile,
			final String dbFile, final String queryLoader, final ProteinSearcher s, final int nQ,
			final int nT, final int nHits, final String sidecarFile, final String kmerSetsFile){
		if(out==null || out.equalsIgnoreCase("stdout")){return;}
		final String metaName=out+".meta";
		final FileFormat ff=FileFormat.testOutput(metaName, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("#schema_version\t1");
		bsw.println("#ProteinSearch run manifest (MVP)");
		bsw.println("query_fasta\t"+queryFile);
		bsw.println("query_sha256\t"+sha256(queryFile));
		bsw.println("query_loader\t"+queryLoader);
		bsw.println("database_fasta\t"+dbFile);
		bsw.println("database_sha256\t"+sha256(dbFile));
		bsw.println("matrix\tBLOSUM62");
		bsw.println("gap_open\t"+Blosum62.GAP_OPEN);
		bsw.println("gap_extend\t"+Blosum62.GAP_EXTEND);
		bsw.println("lambda\t"+Blosum62.LAMBDA);
		bsw.println("K\t"+Blosum62.K);
		bsw.println("k_seed\t"+s.k);
		bsw.println("reduced_seed\t"+s.reducedSeed);
		bsw.println("min_seed_hits\t"+s.minSeedHits);
		bsw.println("evalue_cutoff\t"+s.evalueCutoff);
		bsw.println("min_pident\t"+s.minPident);
		bsw.println("min_coverage\t"+s.minCoverage);
		bsw.println("min_raw_score\t"+s.minRawScore);
		bsw.println("max_target_seqs\t"+s.maxTargetSeqs);
		bsw.println("threads\t"+s.threads);
		bsw.println("shortlist\t"+s.shortlist);
		if(s.shortlist>0){
			bsw.println("aligner\t"+s.aligner);
			bsw.println("sidecar_file\t"+sidecarFile);
			bsw.println("sidecar_consensus_sha256\t"+s.sidecar.consensusSha256);
			bsw.println("sidecar_kmersets_sha256\t"+s.sidecar.kmerSetsSha256);
			bsw.println("kmersets_file\t"+kmerSetsFile);
			bsw.println("kmer_norm\t"+s.kmerNorm);
			bsw.println("triage\t"+(s.triage ? "t" : "f"));
			if(s.triage){
				bsw.println("triage_covdef\t"+s.covDef);
				bsw.println("triage_stage1_survivors\t"+s.triageStage1Survivors.get());
				bsw.println("triage_stage2_survivors\t"+s.triageStage2Survivors.get());
				bsw.println("triage_stage3_survivors\t"+s.triageStage3Survivors.get());
				bsw.println("triage_aligned\t"+s.triageAligned.get());
				bsw.println("triage_hits\t"+s.triageHits.get());
			}
		}
		bsw.println("evalue_calibrated\tfalse");
		bsw.println("evalue_approximate\ttrue");
		bsw.println("queries\t"+nQ);
		bsw.println("database_seqs\t"+nT);
		bsw.println("hits\t"+nHits);
		bsw.poisonAndWait();
	}

	/** Prints usage text and exits. */
	static void printUsageAndExit(){
		System.err.println(
			"ProteinSearch (BBTools prot package) — blastp-style protein search MVP\n"+
			"\n"+
			"Usage: proteinsearch.sh query=<query.faa> db=<database.faa> out=<hits.tsv>\n"+
			"\n"+
			"Required:\n"+
			"  query=   Query protein FASTA (aa).\n"+
			"  db=      Database protein FASTA (aa).\n"+
			"Optional:\n"+
			"  Defaults: k=5, minseedhits=1, evalue=0.001, minid=30, mincoverage=0.7, mts=25,\n"+
			"           reducedseed=t; these are frozen v1 defaults from FamilySearchResources.\n"+
			"  out=     Output TSV (outfmt 6). Default: stdout.\n"+
			"  evalue=  E-value cutoff (default 0.001; frozen v1 via FamilySearchResources).\n"+
			"  minid=   Minimum percent identity (default 30; frozen v1 via FamilySearchResources).\n"+
			"  minscore=Minimum raw BLOSUM62 score (default 0).\n"+
			"  mincoverage=Minimum bidirectional coverage, alnLen/qLen AND alnLen/tLen (default 0.7;\n"+
			"           frozen v1 via FamilySearchResources). mmseqs' --cov-mode 0 equivalent.\n"+
			"  k=       Seed k-mer length (default 5).\n"+
			"  reduced= Use amino8 reduced-alphabet seeds (t/f, default t; frozen v1 via FamilySearchResources).\n"+
			"  queryloader=generic|callgenes  How to assign query IDs (default generic: first\n"+
			"           whitespace token). Use 'callgenes' when query= is CallGenes' raw outa=\n"+
			"           output -- every gene on one contig shares an identical header there, which\n"+
			"           the generic loader collapses onto one id (a real, documented crash via\n"+
			"           ProteinSearcher's duplicate-ID check); 'callgenes' routes through\n"+
			"           CallGenesProteinLoader.load() instead, matching 04_fix_headers_v3.sh's\n"+
			"           per-gene _gN convention.\n"+
			"  mts=     max-target-seqs: cap distinct targets per query.\n"+
			"  threads=Worker threads, parallelized over queries (default 1). Output is\n"+
			"           byte-identical to any other thread count -- pure a speed knob.\n"+
			"  ow=      Overwrite output (t/f, default t).\n"+
			"  shortlist=N  Score every query against every family by F4 (composition+covering-set+\n"+
			"           length fusion) and align only the top N families, instead of every\n"+
			"           seed-passing candidate (default 0: unchanged legacy behavior). Requires\n"+
			"           sidecar= and kmersets=.\n"+
			"  sidecar= FamilyShortlistSidecarBuilder's output artifact (required when shortlist>0).\n"+
			"           Hash-verified against db= and kmersets= at load time; a stale/rebuilt\n"+
			"           sidecar or covering-set file is refused loud, not silently scored.\n"+
			"  kmersets=Covering-set file the sidecar was built from (required when shortlist>0,\n"+
			"           for the hash cross-check above).\n"+
			"  aligner= aaaligner (legacy default) | d55 | glocal | blosum; only when shortlist>0.\n"+
			"           d55 is canonical BLOSUM62 with linear gap 4. glocal/blosum are historical\n"+
			"           aliases for AAAligner.alignGlocal, not the archived prototype classes.\n"+
			"  kmernorm=none (default) | sub | cs300 -- F4's k-mer term: 'none' = raw covering-set hits\n"+
			"           (today's path, byte-identical); 'sub' = hits minus the query's expected coincidental\n"+
			"           count rho_q*|C_f|; 'cs300' = hits/(|C_f|+300) (results/shortlist_kmer_normalization_v1.md).\n"+
			"           Only when shortlist>0.\n"+
			"  triage=  t/f (default f) -- apply the v3 sidecar's per-family thresholds (length-ratio\n"+
			"           window, dimer-cosine floor, covering-kmer floor as a pre-topN survivor mask;\n"+
			"           identity/score-ratio/raw-score/coverage inside the hit filter). Requires a v3\n"+
			"           sidecar (schema_version>=2, built with FamilyShortlistSidecarBuilder's\n"+
			"           thresholds=) -- crashes loud on a v2 sidecar. Only consulted when shortlist>0.\n"+
			"  covdef=  column (default) | span -- triage=t's stage-4 coverage definition: 'column'\n"+
			"           matches mincoverage='s existing alnLen/qLen,alnLen/tLen formula (gap columns\n"+
			"           counted); 'span' uses the aligned span per side, (qStop-qStart+1)/qLen etc.\n"+
			"\n"+
			"Output columns (tab-separated):\n"+
			"  query target pident length mismatch gapopen qstart qend tstart tend evalue bitscore\n"+
			"Note: bitscore is rigorous (gapped BLOSUM62 11/1); E-value is approximate\n"+
			"(edge-length correction omitted) and flagged in the .meta sidecar.\n");
		System.exit(0);
	}
}
