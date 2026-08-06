package prot;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Command-line front end for {@link MarkerVectorizer}: turns a genome bin's
 * protein FASTA into a gene presence/count vector against a marker set and writes
 * it as a TSV plus the derived completeness/contamination scalars.
 *
 * <p>This is the thin CLI/testing wrapper around the in-memory API; the actual
 * vectorization lives in {@link MarkerVectorizer} so it can be called directly
 * from other BBTools code without any file I/O.</p>
 *
 * <p>The marker set is read from a marker-representatives FASTA as written by
 * {@code markerfactory.sh repsout=<file>} — a self-describing FASTA whose headers
 * carry {@code family_id}, {@code domain}, {@code selected}, {@code genomes} and
 * {@code copies} and whose bodies are the representative residues. (The plain
 * marker TSV omits the sequences, so it cannot be scored against on its own.)</p>
 *
 * <p>Usage: {@code markervector.sh bin=<bin.faa> markers=<markers.faa> out=<vec.tsv>}</p>
 *
 * @author Eru
 */
public final class MarkerVectorCLI {

	/**
	 * Program entry point: parses arguments, loads the bin and marker set, builds
	 * the vector, and writes the TSV plus derived scalars.
	 * @param args Command-line arguments (flag=value).
	 */
	public static void main(String[] args){
		final Timer t=new Timer();
		if(args.length==0){printUsageAndExit();}

		String binFile=null, markersFile=null, out=null, domain=null;
		double minId=-1, minCov=-1;
		boolean overwrite=true;

		for(final String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("bin") || a.equals("in") || a.equals("i")){binFile=b;}
			else if(a.equals("markers") || a.equals("markerset") || a.equals("m")){markersFile=b;}
			else if(a.equals("out") || a.equals("o")){out=b;}
			else if(a.equals("domain") || a.equals("d")){domain=b;}
			else if(a.equals("minid") || a.equals("minidentity") || a.equals("id")){
				minId=ClusterProteins.parseIdentity(b);
			}
			else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("cov")){
				minCov=Double.parseDouble(b);
			}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else if(a.equals("-h") || a.equals("--help") || a.equals("help")){printUsageAndExit();}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}

		if(binFile==null || markersFile==null){
			System.err.println("Error: bin= and markers= are both required.\n");
			printUsageAndExit();
		}

		final List<ProteinSequence> bin=ProteinSearch.readFasta(binFile);
		System.err.println("Loaded "+bin.size()+" bin proteins.");

		final LinkedHashMap<String, MarkerSet> byDomain=loadMarkerSets(markersFile);
		final MarkerSet ms=selectDomain(byDomain, domain);
		System.err.println("Scoring against domain '"+ms.domain+"' ("+ms.selectedCount()+
			" selected markers of "+ms.families.size()+" families).");

		final MarkerVectorizer vec=new MarkerVectorizer();
		vec.minIdentity=minId;
		vec.minCoverage=minCov;
		final MarkerVector mv=vec.vectorize(bin, ms);

		writeResults(mv, out, overwrite);

		t.stop();
		System.err.println("Vector dim="+mv.dimension()+" present="+mv.familiesPresent()+
			" exactlyOnce="+mv.familiesExactlyOnce()+" multiCopy="+mv.familiesMultiCopy()+
			" matched="+mv.proteinsMatched+" unmatched="+mv.proteinsUnmatched+" in "+t);
	}

	/**
	 * Loads a marker-representatives FASTA (markerfactory {@code repsout=}) into one
	 * {@link MarkerSet} per domain.
	 * @param markersFile Marker-representatives FASTA path.
	 * @return Domain -&gt; reconstructed marker set (first-seen domain order).
	 */
	static LinkedHashMap<String, MarkerSet> loadMarkerSets(final String markersFile){
		//Per domain: families, version, genome count (from the first record seen).
		final LinkedHashMap<String, ArrayList<MarkerFamily>> fams=
			new LinkedHashMap<String, ArrayList<MarkerFamily>>();
		final LinkedHashMap<String, String> versions=new LinkedHashMap<String, String>();
		final LinkedHashMap<String, Integer> genomeCounts=new LinkedHashMap<String, Integer>();

		final ByteFile bf=ByteFile.makeByteFile(markersFile, false);
		String header=null;
		final ByteBuilder seq=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				if(header!=null){addFamily(header, seq, fams, versions, genomeCounts);}
				header=new String(line);
				seq.clear();
			}else{
				seq.append(line);
			}
		}
		if(header!=null){addFamily(header, seq, fams, versions, genomeCounts);}
		bf.close();
		if(fams.isEmpty()){throw new RuntimeException("No marker families in "+markersFile);}

		final LinkedHashMap<String, MarkerSet> out=new LinkedHashMap<String, MarkerSet>();
		for(final String d : fams.keySet()){
			final ArrayList<MarkerFamily> list=fams.get(d);
			final int gc=genomeCounts.get(d).intValue();
			//Provenance thresholds default to the markerfactory defaults; the CLI's
			//minid/mincov overrides (if any) are applied on the vectorizer instead.
			final MarkerSetProvenance prov=new MarkerSetProvenance("NA", 0.97,
				MarkerVectorizer.DEFAULT_MIN_IDENTITY, MarkerVectorizer.DEFAULT_MIN_COVERAGE,
				5, new ArrayList<String>(), "NA");
			out.put(d, new MarkerSet(d, versions.get(d), prov, gc, list));
		}
		return out;
	}

	/**
	 * Parses one marker-representatives record and appends a {@link MarkerFamily} to
	 * its domain bucket.
	 * @param header Full FASTA header line (including '&gt;').
	 * @param seq Accumulated residues for this record.
	 * @param fams Domain -&gt; family list (updated).
	 * @param versions Domain -&gt; version id (updated).
	 * @param genomeCounts Domain -&gt; genome count (updated).
	 */
	private static void addFamily(final String header, final ByteBuilder seq,
			final LinkedHashMap<String, ArrayList<MarkerFamily>> fams,
			final LinkedHashMap<String, String> versions,
			final LinkedHashMap<String, Integer> genomeCounts){
		final String[] tokens=header.substring(1).trim().split("\\s+");
		if(tokens.length==0 || tokens[0].length()==0){
			throw new RuntimeException("Marker record has no representative id: "+header);
		}
		final String repId=tokens[0];
		int familyId=-1, genomes=0;
		String domain=null, version="v1", copies=null;
		boolean selected=false, sawSelected=false;
		for(int i=1; i<tokens.length; i++){
			final int eq=tokens[i].indexOf('=');
			if(eq<0){continue;}
			final String key=tokens[i].substring(0, eq);
			final String val=tokens[i].substring(eq+1);
			if(key.equals("family_id")){familyId=Integer.parseInt(val);}
			else if(key.equals("domain")){domain=val;}
			else if(key.equals("version")){version=val;}
			else if(key.equals("selected")){selected=val.equals("1"); sawSelected=true;}
			else if(key.equals("genomes")){genomes=Integer.parseInt(val);}
			else if(key.equals("copies")){copies=val;}
		}
		if(familyId<0){throw new RuntimeException("Marker record missing family_id: "+header);}
		if(domain==null){throw new RuntimeException("Marker record missing domain: "+header);}
		if(!sawSelected){throw new RuntimeException("Marker record missing selected flag: "+header);}

		final CopyNumberDistribution dist=parseCopies(copies);
		final ProteinSequence rep=new ProteinSequence(repId, seq.toBytes());
		ArrayList<MarkerFamily> list=fams.get(domain);
		if(list==null){
			list=new ArrayList<MarkerFamily>();
			fams.put(domain, list);
			versions.put(domain, version);
			genomeCounts.put(domain, Integer.valueOf(genomes));
		}
		list.add(new MarkerFamily(familyId, rep, dist, selected));
	}

	/**
	 * Rebuilds a copy-number distribution from a "b0,b1,b2,b3,b4" bins string; an
	 * absent/blank value yields an empty distribution (unused by vectorization).
	 * @param copies Comma-separated bin counts, or null.
	 * @return Reconstructed distribution.
	 */
	private static CopyNumberDistribution parseCopies(final String copies){
		final CopyNumberDistribution dist=new CopyNumberDistribution();
		if(copies==null || copies.length()==0){return dist;}
		final String[] parts=copies.split(",");
		//bins[c] genomes carry the family c times; replay each as add(c).
		for(int c=0; c<parts.length && c<CopyNumberDistribution.BINS; c++){
			final int n=Integer.parseInt(parts[c].trim());
			for(int j=0; j<n; j++){dist.add(c);}
		}
		return dist;
	}

	/**
	 * Chooses which domain's marker set to score against: the requested domain, or
	 * the sole domain when only one is present.
	 * @param byDomain Reconstructed marker sets.
	 * @param requested Requested domain, or null.
	 * @return The selected marker set.
	 */
	static MarkerSet selectDomain(final LinkedHashMap<String, MarkerSet> byDomain,
			final String requested){
		if(requested!=null){
			final MarkerSet ms=byDomain.get(requested);
			if(ms==null){
				throw new RuntimeException("Requested domain '"+requested+"' not in marker file; "+
					"available: "+byDomain.keySet());
			}
			return ms;
		}
		if(byDomain.size()!=1){
			throw new RuntimeException("Marker file has "+byDomain.size()+" domains "+
				byDomain.keySet()+"; specify domain= to choose one.");
		}
		return byDomain.values().iterator().next();
	}

	/** Writes the count vector TSV plus derived-scalar comment lines. */
	static void writeResults(final MarkerVector mv, final String out, final boolean overwrite){
		if(out==null || out.equalsIgnoreCase("stdout")){
			System.out.println(header());
			for(int i=0; i<mv.dimension(); i++){System.out.println(row(mv, i));}
			for(final String s : scalarLines(mv)){System.out.println(s);}
			return;
		}
		final FileFormat ff=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println(header());
		for(int i=0; i<mv.dimension(); i++){bsw.println(row(mv, i));}
		for(final String s : scalarLines(mv)){bsw.println(s);}
		bsw.poisonAndWait();
	}

	/** TSV header line. @return Column header. */
	static String header(){return "#family_id\trepresentative\tcount";}

	/** Formats one family's vector entry as a TSV row. */
	static String row(final MarkerVector mv, final int i){
		return mv.familyIds[i]+"\t"+mv.representativeIds[i]+"\t"+mv.counts[i];
	}

	/** Builds the derived-scalar comment lines appended after the vector. */
	static String[] scalarLines(final MarkerVector mv){
		return new String[]{
			"#domain\t"+mv.domain,
			"#dimension\t"+mv.dimension(),
			"#families_present\t"+mv.familiesPresent(),
			"#families_exactly_once\t"+mv.familiesExactlyOnce(),
			"#families_multi_copy\t"+mv.familiesMultiCopy(),
			"#proteins_matched\t"+mv.proteinsMatched,
			"#proteins_unmatched\t"+mv.proteinsUnmatched
		};
	}

	/** Prints usage text and exits. */
	static void printUsageAndExit(){
		System.err.println(
			"MarkerVector (BBTools prot package) — bin proteins -> gene presence/count vector\n"+
			"\n"+
			"Turns a genome bin's proteins into a fixed-length copy-count vector against a\n"+
			"marker set: entry i = how many of the bin's proteins match selected marker\n"+
			"family i (0 = absent). Feeds the MAG-QC completeness/contamination net.\n"+
			"\n"+
			"Usage: markervector.sh bin=<bin.faa> markers=<markers.faa> out=<vec.tsv>\n"+
			"\n"+
			"Required:\n"+
			"  bin=      Genome bin protein FASTA (aa).\n"+
			"  markers=  Marker-representatives FASTA from 'markerfactory.sh repsout=<file>'\n"+
			"            (self-describing headers; the plain marker TSV lacks sequences).\n"+
			"Optional:\n"+
			"  out=      Output TSV. Default: stdout.\n"+
			"  domain=   Which domain's marker set to score against (required only if the\n"+
			"            marker file holds more than one domain).\n"+
			"  minid=    Override min percent identity to assign a protein to a family\n"+
			"            (default: the marker set's build identity, else 90; 0.9 also ok).\n"+
			"  mincov=   Override min aligned fraction of both protein and representative\n"+
			"            (default: the marker set's build coverage, else 0.8).\n"+
			"  ow=       Overwrite output (t/f, default t).\n"+
			"\n"+
			"Output: one '#family_id  representative  count' row per selected marker family,\n"+
			"followed by '#'-prefixed derived scalars (families present / exactly-once /\n"+
			"multi-copy, proteins matched / unmatched).\n"+
			"Note: real genome FASTA -> proteins (CallGenes) and the QC net itself are not\n"+
			"part of this tool; this only produces the net's input feature vector.\n");
		System.exit(0);
	}
}
