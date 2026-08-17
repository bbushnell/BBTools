package prot;

import java.util.List;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import shared.Timer;

/**
 * Command-line front end for {@link ProteinClusterer}: greedy identity-threshold
 * protein clustering. Reads a protein FASTA and writes a representative-to-member
 * TSV, one row per member.
 *
 * <p>This is the thin CLI/testing wrapper around the in-memory API; the actual
 * clustering logic lives in {@link ProteinClusterer} so it can be called directly
 * from other BBTools code without any file I/O.</p>
 *
 * <p>Usage: {@code clusterproteins.sh in=<proteins.faa> out=<clusters.tsv>}</p>
 *
 * @author Eru
 */
public final class ClusterProteins {

	/**
	 * Program entry point: parses arguments, loads the FASTA, clusters, and
	 * writes the TSV plus a sidecar manifest.
	 * @param args Command-line arguments (flag=value).
	 */
	public static void main(String[] args){
		final Timer t=new Timer();
		if(args.length==0){printUsageAndExit();}

		String inFile=null, out=null;
		boolean overwrite=true;
		final ProteinClusterer clusterer=new ProteinClusterer();

		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("in") || a.equals("input") || a.equals("i")){inFile=b;}
			else if(a.equals("out") || a.equals("o")){out=b;}
			else if(a.equals("k")){clusterer.k=Integer.parseInt(b);}
			else if(a.equals("minseedhits")){clusterer.minSeedHits=Integer.parseInt(b);}
			else if(a.equals("minid") || a.equals("minidentity") || a.equals("id")){clusterer.minIdentity=parseIdentity(b);}
			else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("cov")){clusterer.minCoverage=Double.parseDouble(b);}
			else if(a.equals("reducedseed") || a.equals("reduced")){clusterer.reducedSeed=Parse.parseBoolean(b);}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}

		if(inFile==null){
			System.err.println("Error: in= is required.\n");
			printUsageAndExit();
		}

		final List<ProteinSequence> seqs=ProteinSearch.readFasta(inFile);
		System.err.println("Loaded "+seqs.size()+" sequences.");

		final List<ProteinCluster> clusters=clusterer.cluster(seqs);

		writeResults(clusters, out, overwrite);
		writeSidecar(out, overwrite, inFile, clusterer, seqs.size(), clusters);

		t.stop();
		System.err.println("Formed "+clusters.size()+" clusters from "+seqs.size()+
			" sequences in "+t);
	}

	/**
	 * Parses an identity threshold, accepting either a percent (e.g. 90) or a
	 * fraction (e.g. 0.9, which is scaled to 90).
	 * @param s Raw value.
	 * @return Percent identity in [0,100].
	 */
	static double parseIdentity(final String s){
		final double v=Double.parseDouble(s);
		return v<=1.0 ? v*100.0 : v;
	}

	/** Writes the cluster TSV to a file or stdout. */
	static void writeResults(final List<ProteinCluster> clusters, final String out,
			final boolean overwrite){
		if(out==null || out.equalsIgnoreCase("stdout")){
			System.out.println(header());
			for(ProteinCluster c : clusters){
				for(ClusterMember m : c.members){System.out.println(row(c, m));}
			}
			return;
		}
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println(header());
		for(ProteinCluster c : clusters){
			for(ClusterMember m : c.members){bsw.println(row(c, m));}
		}
		bsw.poisonAndWait();
	}

	/** TSV header line. @return Column header. */
	static String header(){
		return "#cluster_id\trepresentative\tmember\tidentity\tcoverage\tis_representative";
	}

	/**
	 * Formats one member as a TSV row. Kept on StringBuilder/String.format
	 * deliberately: ByteBuilder.append(double,int) special-cases whole numbers
	 * (returns "100" not "100.000" when x0==(long)x0, e.g. every representative's
	 * identity/coverage), which would silently change this column's formatting --
	 * confirmed empirically (representative rows print "100"/"1" under the
	 * ByteBuilder append instead of the required always-3-decimals "100.000"/
	 * "1.000") before reverting. Output-preservation wins over the stylistic
	 * ByteBuilder-everywhere convention when the two conflict.
	 */
	static String row(final ProteinCluster c, final ClusterMember m){
		final boolean isRep=m.seq==c.representative;
		final StringBuilder sb=new StringBuilder();
		sb.append(c.id).append('\t');
		sb.append(c.representative.id).append('\t');
		sb.append(m.id()).append('\t');
		sb.append(String.format("%.3f", m.identity)).append('\t');
		sb.append(String.format("%.3f", m.coverage)).append('\t');
		sb.append(isRep ? "1" : "0");
		return sb.toString();
	}

	/**
	 * Writes a minimal run-manifest sidecar recording the parameters and counts.
	 * A full immutable manifest (checksums, versioned parameters) is deferred;
	 * this records what this run actually used.
	 */
	static void writeSidecar(final String out, final boolean overwrite, final String inFile,
			final ProteinClusterer s, final int nSeqs, final List<ProteinCluster> clusters){
		if(out==null || out.equalsIgnoreCase("stdout")){return;}
		int singletons=0, largest=0;
		for(ProteinCluster c : clusters){
			if(c.isSingleton()){singletons++;}
			if(c.size()>largest){largest=c.size();}
		}
		final String metaName=out+".meta";
		final FileFormat ff=FileFormat.testOutput(metaName, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("#ClusterProteins run manifest");
		bsw.println("input_fasta\t"+inFile);
		bsw.println("matrix\tBLOSUM62");
		bsw.println("gap_open\t"+Blosum62.GAP_OPEN);
		bsw.println("gap_extend\t"+Blosum62.GAP_EXTEND);
		bsw.println("k_seed\t"+s.k);
		bsw.println("reduced_seed\t"+s.reducedSeed);
		bsw.println("min_seed_hits\t"+s.minSeedHits);
		bsw.println("min_identity_pct\t"+s.minIdentity);
		bsw.println("min_coverage\t"+s.minCoverage);
		bsw.println("stable_ids_cross_run\tfalse");
		bsw.println("sequences\t"+nSeqs);
		bsw.println("clusters\t"+clusters.size());
		bsw.println("singletons\t"+singletons);
		bsw.println("largest_cluster\t"+largest);
		bsw.poisonAndWait();
	}

	/** Prints usage text and exits. */
	static void printUsageAndExit(){
		System.err.println(
			"ClusterProteins (BBTools prot package) — greedy identity-threshold protein clustering\n"+
			"\n"+
			"Usage: clusterproteins.sh in=<proteins.faa> out=<clusters.tsv>\n"+
			"\n"+
			"Required:\n"+
			"  in=      Protein FASTA (aa).\n"+
			"Optional:\n"+
			"  out=     Output TSV. Default: stdout. A .meta sidecar is written beside it.\n"+
			"  minid=   Minimum percent identity to join a cluster (default 90; 0.9 also ok).\n"+
			"  mincov=  Minimum aligned fraction of member and rep (default 0.8).\n"+
			"  k=       Seed k-mer length (default 5).\n"+
			"  reduced= Use amino8 reduced-alphabet seeds (t/f, default f).\n"+
			"  ow=      Overwrite output (t/f, default t).\n"+
			"\n"+
			"Output columns (tab-separated):\n"+
			"  cluster_id representative member identity coverage is_representative\n"+
			"One row per member; the representative is listed as a member of its own cluster.\n"+
			"Note: cluster ids are stable within a run only (cross-run stable ids not implemented).\n");
		System.exit(0);
	}
}
