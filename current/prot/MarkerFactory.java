package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Builds per-domain single-copy marker sets from a set of labeled genomes.
 *
 * <p><b>This is the load-bearing API.</b> {@link #build} takes in-memory
 * {@link GenomeProteins} and returns one {@link MarkerSet} per domain — no disk
 * round-trip; callable directly from other BBTools code. Pipeline:</p>
 * <ol>
 * <li>Pool every genome's proteins and remember each protein's genome (by object
 * identity, so pool-unique ids are enough).</li>
 * <li>Cluster the pool into families with {@link ProteinClusterer} — the family
 * grouping engine (Capability 2), not reinvented here.</li>
 * <li>Per domain, per family, tally the {@link CopyNumberDistribution} across
 * that domain's genomes (how many carry it 0/1/2/3/4+ times).</li>
 * <li>Select single-copy markers: families carried <b>exactly once</b> in at
 * least {@link #selectionThreshold} of the domain's genomes. This single test
 * excludes both low-prevalence families (present in too few genomes -&gt; low
 * exactly-once fraction) and multi-copy families (present but rarely exactly
 * once).</li>
 * <li>Emit a versioned, provenance-stamped {@code MarkerSet} per domain.</li>
 * </ol>
 *
 * <p>Only families with prevalence &gt; 0 in a domain are recorded in that
 * domain's set; a family made entirely of another domain's proteins is not.</p>
 *
 * <p>Deferred (reported, not built): real genome corpus + {@code CallGenes}
 * hookup (genome FASTA -&gt; proteins), cross-run stable family ids, the lineage
 * hierarchy with nearest-ancestor fallback, threshold/bin-edge tuning,
 * serialization-format finalization, and any speed/scale work.</p>
 *
 * @author Eru
 */
public final class MarkerFactory {

	/** Clustering engine (family grouping); tune its thresholds before {@link #build}. */
	public final ProteinClusterer clusterer=new ProteinClusterer();
	/** Minimum exactly-once fraction of a domain's genomes to select as a marker. */
	public double selectionThreshold=0.97;

	/**
	 * Builds one marker set per domain from the given genomes.
	 * @param genomes Labeled genomes (non-empty; protein ids pool-unique).
	 * @param version Version id stamped on every produced marker set.
	 * @param buildTimestamp Caller-supplied build timestamp (never read from the
	 *   clock here); may be null (recorded as "NA").
	 * @param taxonomyVersion Taxonomy snapshot id slot; may be null ("NA").
	 * @return Marker sets, one per domain, ordered by domain name.
	 */
	public List<MarkerSet> build(final List<GenomeProteins> genomes, final String version,
			final String buildTimestamp, final String taxonomyVersion){
		if(genomes==null || genomes.isEmpty()){throw new RuntimeException("No genomes supplied.");}

		//Pool all proteins; map each protein to its genome by object identity.
		final ArrayList<ProteinSequence> pool=new ArrayList<ProteinSequence>();
		final IdentityHashMap<ProteinSequence, GenomeProteins> genomeOf=
			new IdentityHashMap<ProteinSequence, GenomeProteins>();
		//Domain -> its genomes, in input order (LinkedHashMap for deterministic grouping).
		final LinkedHashMap<String, ArrayList<GenomeProteins>> byDomain=
			new LinkedHashMap<String, ArrayList<GenomeProteins>>();
		for(final GenomeProteins g : genomes){
			ArrayList<GenomeProteins> dl=byDomain.get(g.domain);
			if(dl==null){dl=new ArrayList<GenomeProteins>(); byDomain.put(g.domain, dl);}
			dl.add(g);
			for(final ProteinSequence p : g.proteins){
				pool.add(p);
				genomeOf.put(p, g);
			}
		}
		if(pool.isEmpty()){throw new RuntimeException("Genomes contain no proteins.");}

		//Cluster the pool into families (the grouping engine).
		final List<ProteinCluster> clusters=clusterer.cluster(pool);

		//For each cluster, tally copies per genome (by object identity), once.
		final ArrayList<IdentityHashMap<GenomeProteins, int[]>> clusterCopies=
			new ArrayList<IdentityHashMap<GenomeProteins, int[]>>(clusters.size());
		for(final ProteinCluster c : clusters){
			final IdentityHashMap<GenomeProteins, int[]> counts=
				new IdentityHashMap<GenomeProteins, int[]>();
			for(final ClusterMember m : c.members){
				final GenomeProteins g=genomeOf.get(m.seq);
				assert(g!=null) : "Clustered protein '"+m.seq.id+"' has no source genome.";
				int[] box=counts.get(g);
				if(box==null){box=new int[1]; counts.put(g, box);}
				box[0]++;
			}
			clusterCopies.add(counts);
		}

		//Build a MarkerSet per domain (domains sorted for deterministic output).
		final ArrayList<MarkerSet> result=new ArrayList<MarkerSet>();
		final ArrayList<String> domains=new ArrayList<String>(byDomain.keySet());
		Collections.sort(domains);
		for(final String domain : domains){
			final ArrayList<GenomeProteins> domGenomes=byDomain.get(domain);
			final ArrayList<MarkerFamily> families=new ArrayList<MarkerFamily>();
			for(int ci=0; ci<clusters.size(); ci++){
				final ProteinCluster c=clusters.get(ci);
				final IdentityHashMap<GenomeProteins, int[]> counts=clusterCopies.get(ci);
				final CopyNumberDistribution dist=new CopyNumberDistribution();
				for(final GenomeProteins g : domGenomes){
					final int[] box=counts.get(g);
					dist.add(box==null ? 0 : box[0]);
				}
				//Only families actually present in this domain belong in its set.
				if(dist.present()==0){continue;}
				final boolean selected=dist.fractionExactlyOnce()>=selectionThreshold;
				families.add(new MarkerFamily(c.id, c.representative, dist, selected));
			}
			Collections.sort(families, FAMILY_BY_ID);
			final ArrayList<String> ids=new ArrayList<String>();
			for(final GenomeProteins g : domGenomes){ids.add(g.id);}
			final MarkerSetProvenance prov=new MarkerSetProvenance(buildTimestamp,
				selectionThreshold, clusterer.minIdentity, clusterer.minCoverage,
				clusterer.k, ids, taxonomyVersion);
			result.add(new MarkerSet(domain, version, prov, domGenomes.size(), families));
		}
		return result;
	}

	/** Orders families by ascending family id (deterministic output). */
	private static final Comparator<MarkerFamily> FAMILY_BY_ID=new Comparator<MarkerFamily>(){
		@Override
		public int compare(final MarkerFamily a, final MarkerFamily b){
			return a.familyId<b.familyId ? -1 : (a.familyId>b.familyId ? 1 : 0);
		}
	};
}
