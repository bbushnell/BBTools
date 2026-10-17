package prot;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/** Reusable, read-only family reference with one alignment workspace per worker.
 * Loads the existing expanded compact HBM package without unpacking its files.
 * Selection is the same full-top50 rule as HbmCompetitiveAssign, not a new
 * functional-name classifier. Callers own translation and annotation policy.
 * @author Keqing
 */
public final class ProteinFamilyAssigner {

	/** Validates stored package bytes and the loader's logical-content bindings.
	 * The manifest is the six-column artifacts.tsv shipped with the reference.
	 * No files are created or changed. Model files must remain immutable while loading.
	 */
	public static ProteinFamilyAssigner load(final Path directory) throws IOException{
		if(directory==null || !Files.isDirectory(directory)){
			throw new IllegalArgumentException("A family-package directory is required: "+directory);
		}
		final HashMap<String,Artifact> artifacts=readArtifacts(directory.resolve("artifacts.tsv"));
		final HashMap<String,String> options=new HashMap<String,String>();
		String modelPin=null;
		for(int i=0; i<KEYS.length; i++){
			final Artifact artifact=artifacts.get(NAMES[i]);
			final Path compressed=directory.resolve(NAMES[i]+".gz");
			final boolean gzip=Files.isRegularFile(compressed);
			final Path file=gzip ? compressed : directory.resolve(NAMES[i]);
			if(!Files.isRegularFile(file)){throw new IllegalArgumentException("Missing family-package member: "+NAMES[i]);}
			final long size=Files.size(file), expected=gzip ? artifact.compressedBytes : artifact.rawBytes;
			if(size!=expected){throw new IllegalArgumentException("Package size mismatch for "+file+": expected "+expected+", observed "+size);}
			final String observed=DigestSuffix.file(file.toString());
			final String pin=gzip ? artifact.compressedPin : artifact.rawPin;
			if(!observed.equals(pin)){throw new IllegalArgumentException("Package stored-byte pin differs for "+file);}
			options.put(KEYS[i], file.toString());
			//Library compares these to decoded scientific pins, not gzip transport bytes.
			options.put(KEYS[i]+"sha80", artifact.rawPin);
			if(KEYS[i].equals("hbm")){modelPin=observed;}
		}
		//The experimental CLI verified decoded metadata before constructing Library.
		//Keep that prerequisite for the three logical pins Library uses directly.
		for(String key : LOGICAL_PINS){
			if(!MagQCTextResource.sha80(options.get(key)).equals(options.get(key+"sha80"))){
				throw new IllegalArgumentException("Package decoded-byte pin differs for "+key);
			}
		}
		final HbmCompetitiveAssign.Library library=new HbmCompetitiveAssign.Library(options);
		validateRoster(options.get("roster"), library.ids, library.memberCounts);
		return new ProteinFamilyAssigner(library, modelPin, options.get("identitysha80"), options.get("consensussha80"), options.get("backgroundsha80"));
	}

	private ProteinFamilyAssigner(final HbmCompetitiveAssign.Library library_,
			final String modelPin_, final String identityPin_, final String consensusPin_, final String backgroundPin_){
		assert(library_!=null && library_.ids.length>=SHORTLIST) : "The frozen competitive rule requires a full top50 reference";
		library=library_;
		modelPin=modelPin_;
		identityPin=identityPin_;
		consensusPin=consensusPin_;
		backgroundPin=backgroundPin_;
	}

	/** Creates mutable query/selection scratch owned exclusively by one caller thread. */
	public Worker newWorker(){return new Worker();}
	public int familyCount(){return library.ids.length;}
	/** Stored HBM content pin, also usable to bind an optional label sidecar. */
	public String modelPin(){return modelPin;}
	/** Decoded stable-family identity table pin; labels must bind this ordering too. */
	public String identityPin(){return identityPin;}
	/** Decoded sequences against which public naming evidence was generated. */
	public String consensusPin(){return consensusPin;}
	/** Verified decoded positional-background pin, independent of transport compression. */
	public String backgroundPin(){return backgroundPin;}
	/** Package-internal read-only access for combined reference competition. */
	HbmCompetitiveAssign.Library library(){return library;}
	public int familyId(final int rank){checkRank(rank); return library.familyIds[rank];}
	public String representative(final int rank){checkRank(rank); return library.ids[rank];}
	/** Number of construction members, not a count of independent taxa. */
	public long memberCount(final int rank){checkRank(rank); return library.memberCounts[rank];}

	private void checkRank(final int rank){
		if(rank<0 || rank>=familyCount()){throw new IllegalArgumentException("Family rank outside the reference: "+rank);}
	}

	/** A Worker is not shareable across threads. The surrounding library is read-only. */
	public final class Worker {
		private Worker(){
			shortlist=new ProteinSearcher.ShortlistScratch(library.sidecar.dims, familyCount(), SHORTLIST);
			selection=new HbmCompetitiveGate.Scratch(familyCount(), SHORTLIST);
		}

		/** Raw amino-acid bytes, with at most one conventional terminal '*'.
		 * Invalid residues and internal stops are rejected by ProteinSequence; callers
		 * should classify unsupported translations before invoking the alignment path.
		 */
		public Hit assign(final String queryId, final byte[] aminoAcids){
			final ProteinSequence query=new ProteinSequence(queryId, aminoAcids);
			ProteinSearcher.scoreFamiliesF4(query, library.sidecar, shortlist);
			HbmCompetitiveAssign.top(shortlist.f4, shortlist.topIdx);
			final int rank=library.gate.select(query.enc, shortlist.topIdx, SHORTLIST, selection);
			assert(selection.scored+selection.lengthRejected==SHORTLIST) : "The competitive gate must account for every shortlisted family";
			return new Hit(rank, rank<0 ? -1 : library.familyIds[rank], rank<0 ? null : library.ids[rank],
				query.length(), selection);
		}

		private final ProteinSearcher.ShortlistScratch shortlist;
		private final HbmCompetitiveGate.Scratch selection;
	}

	/** Detached per-query result; later calls cannot mutate previously returned evidence. */
	public static final class Hit {
		private Hit(final int rank_, final int family_, final String representative_,
				final int length_, final HbmCompetitiveGate.Scratch work){
			this(0, rank_, family_, representative_, length_, work, work.scored, work.lengthRejected);
		}
		/** Combined workers retain local ranks but account for the entire global shortlist. */
		Hit(final int pack_, final int rank_, final int family_, final String representative_,
				final int length_, final HbmCompetitiveGate.Scratch work, final int scored_, final int rejected_){
			packIndex=pack_;
			rank=rank_; familyId=family_; representative=representative_; queryLength=length_;
			score64=work.score64; referenceStart=work.start; referenceEnd=work.end;
			identity=work.metrics.identity; queryCoverage=work.metrics.coverageQ;
			coreCoverage=work.metrics.coverageT; mutualCoverage=work.metrics.coverage;
			scored=scored_; lengthRejected=rejected_;
		}
		public boolean assigned(){return rank>=0;}
		public final int rank, familyId, queryLength, score64, referenceStart, referenceEnd;
		/** Canonical pack index; zero for the legacy single-library API, -1 for combined no-hit. */
		public final int packIndex;
		public final String representative;
		public final float identity, queryCoverage, coreCoverage, mutualCoverage;
		public final int scored, lengthRejected;
	}

	/** Package metadata is small; existing BBTools ByteFile-backed table parsing suffices. */
	private static HashMap<String,Artifact> readArtifacts(final Path path){
		final List<String[]> rows=HmmComparisonData.rows(path.toString());
		if(rows.size()!=NAMES.length+1 || !Arrays.equals(rows.get(0), HEADER)){
			throw new IllegalArgumentException("Unsupported family-package artifacts.tsv header or member count");
		}
		final HashMap<String,Artifact> out=new HashMap<String,Artifact>();
		for(int i=1; i<rows.size(); i++){
			final String[] row=rows.get(i);
			if(row.length!=HEADER.length || !Arrays.asList(NAMES).contains(row[0])){
				throw new IllegalArgumentException("Unknown or malformed family-package artifact at row "+(i+1));
			}
			final long raw=Long.parseLong(row[1]), compressed=Long.parseLong(row[3]);
			if(raw<1 || compressed<1 || !(row[5].equals("bgzf9") || row[5].equals("pigz11"))){
				throw new IllegalArgumentException("Invalid package sizes/codec for "+row[0]);
			}
			final Artifact value=new Artifact(raw, compressed,
				DigestSuffix.requireSuffix(row[2], row[0]), DigestSuffix.requireSuffix(row[4], row[0]));
			if(out.put(row[0], value)!=null){throw new IllegalArgumentException("Duplicate package artifact: "+row[0]);}
		}
		return out;
	}

	/** Roster is a declared package member even though the older CLI derived order from FASTA. */
	private static void validateRoster(final String path, final String[] representatives, final long[] members){
		final List<String[]> rows=HmmComparisonData.rows(path);
		if(rows.size()!=representatives.length+1){throw new IllegalArgumentException("Package roster count differs from consensus");}
		//The exact schema is checked here rather than treating arbitrary rows as an ID list.
		if(!Arrays.equals(rows.get(0), new String[]{"rank", "rep_id", "members"})){
			throw new IllegalArgumentException("Unsupported package roster header");
		}
		for(int i=0; i<representatives.length; i++){
			final String[] row=rows.get(i+1);
			if(row.length!=3 || !row[0].equals(Integer.toString(i)) || !row[1].equals(representatives[i]) || Long.parseLong(row[2])!=members[i]){
				throw new IllegalArgumentException("Package roster order differs at family "+i);
			}
		}
	}

	private static final class Artifact {
		Artifact(long raw, long compressed, String rawPin_, String compressedPin_){
			rawBytes=raw; compressedBytes=compressed; rawPin=rawPin_; compressedPin=compressedPin_;
		}
		final long rawBytes, compressedBytes;
		final String rawPin, compressedPin;
	}

	public static final String RULE="hbm8537_fixed_competitive_v1";
	private static final int SHORTLIST=50;
	private static final String[] HEADER={"artifact", "raw_bytes", "raw_sha80", "compressed_bytes", "compressed_sha80", "codec"};
	private static final String[] KEYS={"hbm", "consensus", "provenance", "background", "identity", "roster", "cores", "kmers", "sidecar"};
	private static final String[] NAMES={"model.hbmc", "consensus.faa", "provenance.tsv", "background.tsv", "family_identity.tsv", "roster.tsv", "cores.tsv", "kmers.tsv", "sidecar.tsv"};
	private static final String[] LOGICAL_PINS={"identity", "consensus", "background"};
	private final HbmCompetitiveAssign.Library library;
	private final String modelPin, identityPin, consensusPin, backgroundPin;
}
