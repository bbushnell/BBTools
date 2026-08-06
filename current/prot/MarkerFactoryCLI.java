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
 * Command-line front end for {@link MarkerFactory}: builds per-domain single-copy
 * marker sets from a manifest of per-genome protein FASTAs and writes the marker
 * set as a TSV plus a provenance sidecar.
 *
 * <p>This is the thin CLI/testing wrapper around the in-memory API; the actual
 * marker-building logic lives in {@link MarkerFactory} so it can be called
 * directly from other BBTools code without any file I/O.</p>
 *
 * <p>Manifest format: a tab-separated text file, one genome per line, columns
 * {@code genome_id  domain  lineage  fasta_path}. Blank lines and lines starting
 * with {@code #} are ignored. Each genome's protein ids are namespaced with the
 * genome id ({@code genome_id~orig_id}) so ids stay pool-unique across genomes.</p>
 *
 * <p>Usage: {@code markerfactory.sh manifest=<genomes.tsv> out=<markers.tsv>}</p>
 *
 * @author Eru
 */
public final class MarkerFactoryCLI {

	/**
	 * Program entry point: parses arguments, loads the manifest, builds the marker
	 * sets, and writes the TSV plus a provenance sidecar.
	 * @param args Command-line arguments (flag=value).
	 */
	public static void main(String[] args){
		final Timer t=new Timer();
		if(args.length==0){printUsageAndExit();}

		String manifest=null, out=null, version="v1", timestamp=null, taxonomy=null;
		String repsOut=null;
		boolean overwrite=true;
		final MarkerFactory factory=new MarkerFactory();

		for(final String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("manifest") || a.equals("in") || a.equals("i")){manifest=b;}
			else if(a.equals("out") || a.equals("o")){out=b;}
			else if(a.equals("repsout") || a.equals("repfasta") || a.equals("fasta")){repsOut=b;}
			else if(a.equals("version") || a.equals("v")){version=b;}
			else if(a.equals("timestamp") || a.equals("ts")){timestamp=b;}
			else if(a.equals("taxonomy") || a.equals("taxversion")){taxonomy=b;}
			else if(a.equals("threshold") || a.equals("minsingle") || a.equals("t")){
				factory.selectionThreshold=Double.parseDouble(b);
			}
			else if(a.equals("minid") || a.equals("minidentity") || a.equals("id")){
				factory.clusterer.minIdentity=ClusterProteins.parseIdentity(b);
			}
			else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("cov")){
				factory.clusterer.minCoverage=Double.parseDouble(b);
			}
			else if(a.equals("k")){factory.clusterer.k=Integer.parseInt(b);}
			else if(a.equals("minseedhits")){factory.clusterer.minSeedHits=Integer.parseInt(b);}
			else if(a.equals("reducedseed") || a.equals("reduced")){
				factory.clusterer.reducedSeed=Parse.parseBoolean(b);
			}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}

		if(manifest==null){
			System.err.println("Error: manifest= is required.\n");
			printUsageAndExit();
		}

		final List<GenomeProteins> genomes=loadManifest(manifest);
		int totalProteins=0;
		for(final GenomeProteins g : genomes){totalProteins+=g.size();}
		System.err.println("Loaded "+genomes.size()+" genomes, "+totalProteins+" proteins.");

		final List<MarkerSet> sets=factory.build(genomes, version, timestamp, taxonomy);

		writeResults(sets, out, overwrite);
		writeProvenance(sets, out, overwrite, manifest, factory);
		if(repsOut!=null){writeRepresentatives(sets, repsOut, overwrite);}

		t.stop();
		final StringBuilder sb=new StringBuilder();
		for(final MarkerSet ms : sets){
			sb.append(" ").append(ms.domain).append("=").append(ms.selectedCount())
				.append("/").append(ms.families.size());
		}
		System.err.println("Built "+sets.size()+" domain marker set(s) [selected/families]:"+sb+
			" in "+t);
	}

	/**
	 * Loads a genome manifest into labeled in-memory genomes.
	 * @param manifestFile Manifest path (genome_id, domain, lineage, fasta per line).
	 * @return Labeled genomes in manifest order.
	 */
	static List<GenomeProteins> loadManifest(final String manifestFile){
		final ArrayList<GenomeProteins> genomes=new ArrayList<GenomeProteins>();
		final ByteFile bf=ByteFile.makeByteFile(manifestFile, false);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final String s=new String(line);
			if(s.trim().length()==0){continue;}
			final String[] parts=s.split("\t");
			if(parts.length<4){
				throw new RuntimeException("Manifest line needs 4 tab-separated columns "+
					"(genome_id domain lineage fasta): '"+s+"'");
			}
			final String gid=parts[0].trim();
			final String domain=parts[1].trim();
			final String lineage=parts[2].trim();
			final String fasta=parts[3].trim();
			final List<ProteinSequence> prots=readGenomeFasta(fasta, gid);
			genomes.add(new GenomeProteins(gid, domain, lineage.length()==0 ? null : lineage, prots));
		}
		bf.close();
		if(genomes.isEmpty()){throw new RuntimeException("No genomes in manifest "+manifestFile);}
		return genomes;
	}

	/**
	 * Reads one genome's protein FASTA, namespacing each id with the genome id so
	 * ids stay unique across the pool.
	 * @param fname FASTA path.
	 * @param genomeId Genome id used as the id prefix.
	 * @return Namespaced protein sequences.
	 */
	static List<ProteinSequence> readGenomeFasta(final String fname, final String genomeId){
		final ArrayList<ProteinSequence> list=new ArrayList<ProteinSequence>();
		final ByteFile bf=ByteFile.makeByteFile(fname, false);
		String id=null;
		final ByteBuilder seq=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				if(id!=null){list.add(new ProteinSequence(id, seq.toBytes()));}
				id=genomeId+"~"+ProteinSearch.parseHeaderId(line);
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

	/** Writes the marker-set TSV (all domains) to a file or stdout. */
	static void writeResults(final List<MarkerSet> sets, final String out, final boolean overwrite){
		if(out==null || out.equalsIgnoreCase("stdout")){
			System.out.println(header());
			for(final MarkerSet ms : sets){
				for(final MarkerFamily f : ms.families){System.out.println(row(ms, f));}
			}
			return;
		}
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println(header());
		for(final MarkerSet ms : sets){
			for(final MarkerFamily f : ms.families){bsw.println(row(ms, f));}
		}
		bsw.poisonAndWait();
	}

	/** TSV header line. @return Column header. */
	static String header(){
		return "#domain\tmarker_set_version\tfamily_id\trepresentative\tprevalence\t"+
			"n_genomes\tcopies_0\tcopies_1\tcopies_2\tcopies_3\tcopies_4plus\t"+
			"fraction_exactly_once\tselected_single_copy";
	}

	/** Formats one marker family as a TSV row. */
	static String row(final MarkerSet ms, final MarkerFamily f){
		final int[] bins=f.dist.bins;
		final StringBuilder sb=new StringBuilder();
		sb.append(ms.domain).append('\t');
		sb.append(ms.version).append('\t');
		sb.append(f.familyId).append('\t');
		sb.append(f.representative.id).append('\t');
		sb.append(String.format("%.3f", f.prevalence())).append('\t');
		sb.append(ms.genomeCount).append('\t');
		sb.append(bins[0]).append('\t');
		sb.append(bins[1]).append('\t');
		sb.append(bins[2]).append('\t');
		sb.append(bins[3]).append('\t');
		sb.append(bins[4]).append('\t');
		sb.append(String.format("%.3f", f.fractionExactlyOnce())).append('\t');
		sb.append(f.selectedSingleCopy ? "1" : "0");
		return sb.toString();
	}

	/**
	 * Writes a provenance sidecar recording the parameters and, per domain, the
	 * source genome ids. A full immutable build manifest (input/output checksums,
	 * commit state, taxonomy snapshot) is deferred.
	 */
	static void writeProvenance(final List<MarkerSet> sets, final String out,
			final boolean overwrite, final String manifest, final MarkerFactory factory){
		if(out==null || out.equalsIgnoreCase("stdout")){return;}
		final String provName=out+".prov";
		final FileFormat ff=FileFormat.testOutput(provName, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("#MarkerFactory provenance");
		bsw.println("manifest\t"+manifest);
		bsw.println("matrix\tBLOSUM62");
		bsw.println("selection_threshold\t"+factory.selectionThreshold);
		bsw.println("min_identity_pct\t"+factory.clusterer.minIdentity);
		bsw.println("min_coverage\t"+factory.clusterer.minCoverage);
		bsw.println("k_seed\t"+factory.clusterer.k);
		bsw.println("reduced_seed\t"+factory.clusterer.reducedSeed);
		bsw.println("stable_family_ids_cross_run\tfalse");
		for(final MarkerSet ms : sets){
			bsw.println("#domain "+ms.domain);
			bsw.println("domain\t"+ms.domain);
			bsw.println("marker_set_version\t"+ms.version);
			bsw.println("build_timestamp\t"+ms.provenance.buildTimestamp);
			bsw.println("taxonomy_version\t"+ms.provenance.taxonomyVersion);
			bsw.println("genomes\t"+ms.genomeCount);
			bsw.println("families\t"+ms.families.size());
			bsw.println("selected_markers\t"+ms.selectedCount());
			final StringBuilder ids=new StringBuilder();
			for(int i=0; i<ms.provenance.sourceGenomeIds.size(); i++){
				if(i>0){ids.append(',');}
				ids.append(ms.provenance.sourceGenomeIds.get(i));
			}
			bsw.println("source_genome_ids\t"+ids);
		}
		bsw.poisonAndWait();
	}

	/**
	 * Writes a self-describing marker-representatives FASTA: one record per family
	 * per domain, the header carrying the family metadata needed to reconstruct a
	 * marker set (family id, domain, selected flag, copy bins, genome count) and
	 * the body being the representative's residues. This is the file
	 * {@code markervector.sh} reads to score a bin against the marker set (the
	 * marker TSV alone omits the representative sequences).
	 * @param sets Marker sets to serialize.
	 * @param repsOut Output FASTA path.
	 * @param overwrite Overwrite existing output.
	 */
	static void writeRepresentatives(final List<MarkerSet> sets, final String repsOut,
			final boolean overwrite){
		final FileFormat ff=FileFormat.testOutput(repsOut, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		for(final MarkerSet ms : sets){
			for(final MarkerFamily f : ms.families){
				final int[] b=f.dist.bins;
				final StringBuilder sb=new StringBuilder();
				sb.append('>').append(f.representative.id);
				sb.append(" family_id=").append(f.familyId);
				sb.append(" domain=").append(ms.domain);
				sb.append(" version=").append(ms.version);
				sb.append(" selected=").append(f.selectedSingleCopy ? "1" : "0");
				sb.append(" genomes=").append(ms.genomeCount);
				sb.append(" copies=").append(b[0]).append(',').append(b[1]).append(',')
					.append(b[2]).append(',').append(b[3]).append(',').append(b[4]);
				bsw.println(sb.toString());
				bsw.println(decode(f.representative.enc));
			}
		}
		bsw.poisonAndWait();
	}

	/**
	 * Decodes BBTools-encoded residues back to ASCII amino-acid letters. Standard
	 * residues (0-19) map via {@code AminoAcid.numberToAcid}; the ambiguous code
	 * {@code Blosum62.X_CODE} maps to 'X'. Re-encoding the result reproduces the
	 * same encoded array, so the round-trip is stable for search/alignment.
	 * @param enc Encoded residues (values 0-19 or X_CODE).
	 * @return ASCII residue string.
	 */
	static String decode(final byte[] enc){
		final StringBuilder sb=new StringBuilder(enc.length);
		for(final byte e : enc){
			if(e>=0 && e<=19){sb.append((char)dna.AminoAcid.numberToAcid[e]);}
			else if(e==Blosum62.X_CODE){sb.append('X');}
			else{throw new RuntimeException("Unencoded residue in decode: "+e);}
		}
		return sb.toString();
	}

	/** Prints usage text and exits. */
	static void printUsageAndExit(){
		System.err.println(
			"MarkerFactory (BBTools prot package) — build per-domain single-copy marker sets\n"+
			"\n"+
			"Usage: markerfactory.sh manifest=<genomes.tsv> out=<markers.tsv>\n"+
			"\n"+
			"Required:\n"+
			"  manifest= Tab-separated file, one genome per line:\n"+
			"            genome_id <tab> domain <tab> lineage <tab> fasta_path\n"+
			"            (blank lines and #-comments ignored)\n"+
			"Optional:\n"+
			"  out=       Output TSV. Default: stdout. A .prov sidecar is written beside it.\n"+
			"  repsout=   Optional marker-representatives FASTA (self-describing headers);\n"+
			"             this is the file markervector.sh reads to score a bin.\n"+
			"  threshold= Min fraction of a domain's genomes carrying a family EXACTLY ONCE\n"+
			"             to select it as a single-copy marker (default 0.97).\n"+
			"  minid=     Clustering min percent identity (default 90; 0.9 also ok).\n"+
			"  mincov=    Clustering min aligned fraction (default 0.8).\n"+
			"  k=         Clustering seed k-mer length (default 5).\n"+
			"  reduced=   Use amino8 reduced-alphabet seeds (t/f, default f).\n"+
			"  version=   Marker-set version id (default v1).\n"+
			"  timestamp= Build timestamp recorded in provenance (default NA; not read from clock).\n"+
			"  taxonomy=  Taxonomy snapshot id recorded in provenance (default NA).\n"+
			"  ow=        Overwrite output (t/f, default t).\n"+
			"\n"+
			"Output columns (tab-separated):\n"+
			"  domain marker_set_version family_id representative prevalence n_genomes\n"+
			"  copies_0 copies_1 copies_2 copies_3 copies_4plus fraction_exactly_once\n"+
			"  selected_single_copy\n"+
			"One row per family per domain; families with zero prevalence in a domain are omitted.\n"+
			"Note: family ids are stable within a run only (cross-run stable ids not implemented).\n");
		System.exit(0);
	}
}
