package prot;

import java.util.ArrayList;
import java.util.List;

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

		String queryFile=null, dbFile=null, out=null;
		boolean overwrite=true;
		final ProteinSearcher searcher=new ProteinSearcher();

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
			else if(a.equals("maxtargetseqs") || a.equals("mts")){searcher.maxTargetSeqs=Integer.parseInt(b);}
			else if(a.equals("reducedseed") || a.equals("reduced")){searcher.reducedSeed=Parse.parseBoolean(b);}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}

		if(queryFile==null || dbFile==null){
			System.err.println("Error: both query= and db= are required.\n");
			printUsageAndExit();
		}

		final List<ProteinSequence> queries=readFasta(queryFile);
		final List<ProteinSequence> targets=readFasta(dbFile);
		System.err.println("Loaded "+queries.size()+" queries and "+targets.size()+" database sequences.");

		final List<ProteinHit> hits=searcher.search(queries, targets);

		writeResults(hits, out, overwrite);
		writeSidecar(out, overwrite, queryFile, dbFile, searcher, queries.size(),
			targets.size(), hits.size());

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
		final ByteBuilder seq=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				if(id!=null){list.add(new ProteinSequence(id, seq.toBytes()));}
				id=parseHeaderId(line);
				seq.clear();
			}else{
				seq.append(line);
			}
		}
		if(id!=null){list.add(new ProteinSequence(id, seq.toBytes()));}
		bf.close();
		if(list.isEmpty()){throw new RuntimeException("No sequences found in "+fname);}
		return list;
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
			final String dbFile, final ProteinSearcher s, final int nQ, final int nT, final int nHits){
		if(out==null || out.equalsIgnoreCase("stdout")){return;}
		final String metaName=out+".meta";
		final FileFormat ff=FileFormat.testOutput(metaName, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("#ProteinSearch run manifest (MVP)");
		bsw.println("query_fasta\t"+queryFile);
		bsw.println("database_fasta\t"+dbFile);
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
		bsw.println("min_raw_score\t"+s.minRawScore);
		bsw.println("max_target_seqs\t"+s.maxTargetSeqs);
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
			"  out=     Output TSV (outfmt 6). Default: stdout.\n"+
			"  evalue=  E-value cutoff (default 10).\n"+
			"  minid=   Minimum percent identity (default 0).\n"+
			"  minscore=Minimum raw BLOSUM62 score (default 0).\n"+
			"  k=       Seed k-mer length (default 5).\n"+
			"  reduced= Use amino8 reduced-alphabet seeds (t/f, default f).\n"+
			"  mts=     max-target-seqs: cap distinct targets per query.\n"+
			"  ow=      Overwrite output (t/f, default t).\n"+
			"\n"+
			"Output columns (tab-separated):\n"+
			"  query target pident length mismatch gapopen qstart qend tstart tend evalue bitscore\n"+
			"Note: bitscore is rigorous (gapped BLOSUM62 11/1); E-value is approximate\n"+
			"(edge-length correction omitted) and flagged in the .meta sidecar.\n");
		System.exit(0);
	}
}
