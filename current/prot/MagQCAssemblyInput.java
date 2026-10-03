package prot;

import java.util.ArrayList;
import java.util.HashMap;

import cardinality.DynamicDemiLog;
import clade.Clade;
import clade.SendClade;
import fileIO.FileFormat;
import prok.CallGenes;
import prok.GeneModel;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import bin.AdjustEntropy;
import tracker.EntropyTracker;

/**
 * Converts one assembly into the same sufficient statistics used by prepared-input inference.
 * Taxonomy comes from one QuickClade request, never from contig identifiers.
 * @author Yoimiya
 */
final class MagQCAssemblyInput {

	private MagQCAssemblyInput(){}

	/** Loads the bound assignment resources, classifies once, and calls the predicted-phylum model. */
	static MagQCPreparedBin build(HashMap<String,String> options){
		final Taxonomy supplied=overrideTaxonomy(options);
		final String id=required(options, "binid");
		if(id.indexOf('\t')>=0 || id.indexOf('\n')>=0 || id.indexOf('\r')>=0){
			throw new IllegalArgumentException("binid must fit one TSV field");
		}
		final String mode=options.containsKey("pgmmode") ? required(options, "pgmmode") : "taxonomy";
		if(!mode.equals("taxonomy") && !mode.equals("default")){
			throw new IllegalArgumentException("pgmmode must be taxonomy or default");
		}
		final int passes=options.containsKey("passes") ? Integer.parseInt(required(options, "passes")) : 1;
		if(passes<1){throw new IllegalArgumentException("passes must be positive");}
		final ProteinSearcher.AssignPolicy policy=ProteinSearcher.AssignPolicy.valueOf(required(options, "policy"));
		final int lookahead=Integer.parseInt(required(options, "lookahead"));
		ProteinSearcher.validateLookahead(policy, lookahead);
		final ProteinSearcher.AssignmentBinding binding=loadBinding(options);
		final ArrayList<Read> contigs=readContigs(required(options, "fasta"));
		final Taxonomy taxonomy=supplied==null ? classify(contigs, options.get("taxaddress")) : supplied;
		final GeneModel pgm=mode.equals("default") ? GeneCallAdapter.defaultModel() :
			CallGenes.getPhylumPGM(taxonomy.phylum.equals("unknown") ? null : taxonomy.phylum);
		final ProteinSearcher searcher=new ProteinSearcher();
		searcher.aligner="d55";
		final MagQCPreparedBin result=FastaInCacheRowBuilder.buildPreparedBin(contigs, id,
			taxonomy.status, taxonomy.domain, taxonomy.phylum, pgm, passes, binding, searcher, policy, lookahead);
		options.put("pgmmode", mode); options.put("passes", Integer.toString(passes));
		options.put("qc_status", taxonomy.status); options.put("qc_domain", taxonomy.domain);
		options.put("qc_phylum", taxonomy.phylum);
		return result;
	}

	/** Loads one immutable assignment binding shared by all bins in a dedicated batch. */
	static ProteinSearcher.AssignmentBinding loadBinding(HashMap<String,String> options){
		final String profile=required(options, "profile");
		final String profilePin=required(options, "profilesha80");
		if(!profilePin.matches("[0-9a-f]{20}") || !profilePin.equals(MagQCTextResource.sha80(profile))){
			throw new IllegalArgumentException("Assignment profile sha80 mismatch");
		}
		final ProteinSearcher.AssignmentBinding binding=ProteinSearcher.AssignmentBinding.loadSchema7(
			profile, profilePin, required(options, "roster"), required(options, "ref"),
			required(options, "rolemanifest"), required(options, "core"),
			required(options, "coveringsets"), required(options, "sidecar"),
			required(options, "hbmbundle"), required(options, "hbmprovenance"));
		final HbmRosterLoader features=HbmRosterLoader.load(required(options, "familylist"));
		if(features.size()!=binding.familyCount()){
			throw new IllegalArgumentException("Assignment and subnet family counts differ");
		}
		for(int rank=0; rank<features.size(); rank++){
			if(!features.repIdByRank[rank].equals(binding.repId(rank))){
				throw new IllegalArgumentException("Assignment/subnet representative mismatch at rank "+rank);
			}
		}
		return binding;
	}

	/** Reads a complete assembly through the native single-reader FASTA path. */
	static ArrayList<Read> readContigs(String path){
		final FileFormat input=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		if(!input.fasta()){throw new IllegalArgumentException("fasta= requires an assembly FASTA");}
		final ArrayList<Read> contigs=StreamerFactory.getReads(-1, false, input, null, null, null);
		if(contigs.isEmpty()){throw new IllegalArgumentException("Assembly FASTA contains no contigs");}
		return contigs;
	}

	/**
	 * Resolves an explicit offline override without a server call or header-derived TID.
	 * A phylum requires its domain; domain-only overrides retain unknown phylum.
	 */
	static Taxonomy overrideTaxonomy(HashMap<String,String> options){
		assert(options!=null) : "Taxonomy selection requires the parsed assembly configuration";
		if(!options.containsKey("taxdomain") && !options.containsKey("taxphylum")){return null;}
		final String value=required(options, "taxdomain");
		final String domain=value.equalsIgnoreCase("Bacteria") ? "Bacteria" :
			value.equalsIgnoreCase("Archaea") ? "Archaea" : null;
		if(domain==null){throw new IllegalArgumentException("taxdomain requires Bacteria or Archaea");}
		final String phylum=options.containsKey("taxphylum") ? required(options, "taxphylum") : "unknown";
		for(int i=0; i<phylum.length(); i++){
			if(Character.isISOControl(phylum.charAt(i))){throw new IllegalArgumentException("taxphylum must fit one TSV field");}
		}
		if(!phylum.equals(phylum.trim())){throw new IllegalArgumentException("taxphylum has leading or trailing whitespace");}
		return new Taxonomy(domain, phylum);
	}

	/** Uses the C1 CallGenes request recipe; all classification precedes MAG-QC network loading. */
	static Taxonomy classify(ArrayList<Read> contigs, String address){
		synchronized(Clade.class){
			try(SketchSession session=openSketchSession()){return session.classify(contigs, address);}
		}
	}

	/**
	 * Fixes the C1 sketch recipe once, before batch taxonomy workers start. The
	 * dedicated client JVM must not run unrelated Clade configuration concurrently.
	 * Classification finishes before subnet initialization and gene-model loading.
	 */
	static SketchSession openSketchSession(){
		synchronized(Clade.class){
			if(activeSketchSession!=null){throw new IllegalStateException("A taxonomy sketch session is already active");}
			final SketchSettings previous=new SketchSettings();
			try{
				Clade.MAKE_DDLS=true; Clade.DDL_K=25; Clade.DDL_BUCKETS=32768; Clade.DDL_SEED=12345L;
				DynamicDemiLog.setExponent(5);
				if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150){AdjustEntropy.load(4, 150);}
				activeSketchSession=new SketchSession(previous);
				return activeSketchSession;
			}catch(RuntimeException e){previous.close(); throw e;}
			catch(Error e){previous.close(); throw e;}
		}
	}

	/** Owns setup/restoration; per-bin sketch arrays and entropy counters are private. */
	static final class SketchSession implements AutoCloseable{
		private SketchSession(SketchSettings previous_){previous=previous_;}

		/** Sends the same one-query C1 wire request as the legacy single-bin path. */
		Taxonomy classify(ArrayList<Read> contigs, String address){
			final String endpoint=address==null || address.equals("refseq") ? null : address;
			if(endpoint!=null && endpoint.isEmpty()){throw new IllegalArgumentException("Empty taxaddress");}
			final SketchRequest request=request(contigs);
			final String response=SendClade.sendMessage(request.bytes, endpoint, false);
			return parseResponse(response, request.bases, request.contigs);
		}

		/**
		 * Serializes while the recipe is fixed. Network I/O and response parsing no
		 * longer depend on sketch globals and need not hold their configuration lock.
		 */
		SketchRequest request(ArrayList<Read> contigs){
			synchronized(Clade.class){
				if(activeSketchSession!=this){throw new IllegalStateException("Taxonomy sketch session is closed");}
				active++; peak=Math.max(peak, active);
			}
			try{
				if(contigs==null || contigs.isEmpty()){throw new IllegalArgumentException("QuickClade requires contigs");}
				final Clade query=new Clade(0, 0, QUERY_NAME);
				final EntropyTracker entropy=new EntropyTracker(4, 150, false);
				for(Read r:contigs){
					if(r==null || r.bases==null || r.bases.length==0){
						throw new IllegalArgumentException("Empty/null assembly contig");
					}
					query.add(r.bases, entropy);
				}
				query.finish();
				final ArrayList<Clade> queries=new ArrayList<Clade>(1); queries.add(query);
				final byte[] bytes=SendClade.toMessage(queries, true, 1, false, false, 1, 1);
				// Native transport includes undefined monomers in the echoed Q_Bases.
				return new SketchRequest(bytes, query.monomerSum(), contigs.size());
			}finally{
				synchronized(Clade.class){active--;}
			}
		}

		/** Diagnostic overlap count; it is not a throughput measurement. */
		int peakWorkers(){synchronized(Clade.class){return peak;}}

		/** Restores prior settings only after the caller has joined sketch workers. */
		@Override public void close(){
			synchronized(Clade.class){
				if(activeSketchSession!=this){return;}
				if(active!=0){throw new IllegalStateException("Cannot close taxonomy session with active sketch workers");}
				try{previous.close();}finally{activeSketchSession=null;}
			}
		}

		private final SketchSettings previous;
		private int active, peak;
	}

	/** Fully serialized request: no mutable Clade instance escapes to a later phase. */
	static final class SketchRequest{
		SketchRequest(byte[] bytes_, long bases_, long contigs_){
			assert(bytes_!=null && bytes_.length>0 && bases_>0 && contigs_>0) :
				"A serialized nonempty assembly must carry bytes and echoed query counts";
			bytes=bytes_; bases=bases_; contigs=contigs_;
		}
		final byte[] bytes;
		final long bases, contigs;
	}

	/** The same five restored globals as the original single-bin classifier. */
	private static final class SketchSettings implements AutoCloseable{
		@Override public void close(){
			Clade.MAKE_DDLS=ddl; Clade.DDL_K=k; Clade.DDL_BUCKETS=buckets; Clade.DDL_SEED=seed;
			DynamicDemiLog.setExponent(exponent);
		}
		private final boolean ddl=Clade.MAKE_DDLS;
		private final int k=Clade.DDL_K, buckets=Clade.DDL_BUCKETS, exponent=DynamicDemiLog.exponentBits();
		private final long seed=Clade.DDL_SEED;
	}

	private static SketchSession activeSketchSession;

	/**
	 * Accepts exactly one machine query and at most one result. A bare #Query1 is a genuine no-hit;
	 * null, error text, extra results and malformed rows are failures, never missing taxonomy.
	 * The lineage field is found independently of optional SSU/sketch columns. Confidence is not
	 * a taxonomy filter: C1 used the returned reference lineage, not ConfLevel/Confidence.
	 */
	static Taxonomy parseResponse(String response, long bases, long contigs){
		if(response==null || response.isEmpty()){throw malformed("empty server response");}
		final String[] lines=response.split("\n", -1);
		boolean query=false, result=false;
		Taxonomy taxonomy=new Taxonomy("unknown", "unknown");
		for(String raw:lines){
			final String line=raw.endsWith("\r") ? raw.substring(0, raw.length()-1) : raw;
			if(line.isEmpty()){continue;}
			if(!query){
				if(!line.equals("#Query1")){throw malformed("expected single-query machine envelope");}
				query=true; continue;
			}
			if(result || line.startsWith("#")){throw malformed("extra query or result");}
			result=true;
			final String[] fields=line.split("\t", -1);
			if(fields.length<20 || !fields[0].equals(QUERY_NAME)){throw malformed("mismatched query or short row");}
			if(Long.parseLong(fields[2])!=bases || Long.parseLong(fields[3])!=contigs){
				throw malformed("query base/contig counts differ from submitted assembly");
			}
			Long.parseLong(fields[5]); Long.parseLong(fields[7]); Long.parseLong(fields[8]);
			finite(fields[1]); finite(fields[6]);
			for(int i=10; i<17; i++){finite(fields[i]);}
			int lineage=-1;
			for(int i=17; i<fields.length-2; i++){
				if(fields[i].contains("__") || fields[i].equals("NA") ||
					fields[i].startsWith("d:") || fields[i].startsWith("p:")){
					lineage=i; break;
				}
			}
			if(lineage<0){throw malformed("missing lineage field");}
			for(int i=17; i<lineage; i++){finite(fields[i]);}
			final Taxonomy ranks=parseLineage(fields[lineage]);
			// Comparison.appendResultMachine writes six DDL fields immediately before
			// lineage, optionally preceded by one SSU field. ANI is a fraction.
			final int optional=lineage-17;
			if(optional!=0 && optional!=1 && optional!=6 && optional!=7){
				throw malformed("unsupported optional metric layout");
			}
			final double ani=optional>=6 ? Double.parseDouble(fields[lineage-6]) : -1;
			if(ani < -1 || ani>1){throw malformed("ANI outside its fraction range");}
			taxonomy=new Taxonomy(ranks.domain, ranks.phylum, fields[4], Long.parseLong(fields[5]), ani);
		}
		if(!query){throw malformed("missing query envelope");}
		return taxonomy;
	}

	/** Parses named taxonomy ranks only; duplicate or malformed ranks cannot invent a label. */
	private static Taxonomy parseLineage(String lineage){
		if(lineage.equals("NA")){return new Taxonomy("unknown", "unknown");}
		String domain="unknown", phylum="unknown";
		boolean domainSeen=false, phylumSeen=false;
		for(String raw:lineage.split(";", -1)){
			final String rank=raw.trim();
			if(rank.isEmpty()){continue;}
			final int doubleUnderscore=rank.indexOf("__");
			final int separator=doubleUnderscore<0 ? rank.indexOf(':') : doubleUnderscore;
			if(separator<1 || separator>2){throw malformed("malformed lineage rank");}
			final int prefixLength=doubleUnderscore<0 ? 1 : 2;
			final String prefix=rank.substring(0, separator), name=rank.substring(separator+prefixLength).trim();
			// Clade.lineage can emit sk__, st__, and repeated ss__ ancestors. Preserve C1's
			// domain/phylum projection without rejecting those valid finer ranks.
			if(prefix.equals("d")){
				if(domainSeen){throw malformed("repeated domain");} domainSeen=true;
				if(!name.isEmpty()){domain=name;}
			}else if(prefix.equals("p")){
				if(phylumSeen){throw malformed("repeated phylum");} phylumSeen=true;
				if(!name.isEmpty()){phylum=name;}
			}
		}
		return new Taxonomy(domain, phylum);
	}

	/** Rejects error text or nonfinite metrics in the numeric portion of a machine row. */
	private static void finite(String text){
		if(!Double.isFinite(Double.parseDouble(text))){throw malformed("nonfinite comparison metric");}
	}

	/** Adds assembly-only provenance without changing existing prepared-input output bytes. */
	static void appendProvenance(ByteBuilder out, HashMap<String,String> options){
		if(!options.containsKey("fasta")){return;}
		out.append("#input_mode\tassembly_fasta\n#taxonomy_source\t")
			.append(options.containsKey("taxdomain") ? "user-supplied" : "QuickClade server").nl();
		for(String key:new String[]{"profilesha80", "policy", "lookahead", "pgmmode", "passes", "qc_status", "qc_domain", "qc_phylum"}){
			out.append('#').append(key).tab().append(options.get(key)).nl();
		}
		out.append("#caller_asymmetry\ttraining true-phylum PGM; deployment ")
			.append(options.containsKey("taxdomain") ? "user-supplied-phylum" : "predicted-phylum")
			.append(" PGM unless pgmmode=default\n");
	}

	/** Requires explicit assignment inputs; the binding validates all coupled resource identities. */
	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Required assembly argument: "+key);}
		return value;
	}

	/** Keeps a classifier failure distinguishable from an ordinary unknown classification. */
	private static IllegalArgumentException malformed(String reason){
		return new IllegalArgumentException("Invalid QuickClade response: "+reason);
	}

	/** Immutable deployment taxonomy; reference identity is retained only for reporting. */
	static final class Taxonomy{
		Taxonomy(String domain_, String phylum_){
			this(domain_, phylum_, "NA", -1, -1);
		}
		Taxonomy(String domain_, String phylum_, String name_, long taxId_, double ani_){
			domain=domain_; phylum=phylum_;
			if(domain.equals("unknown") && !phylum.equals("unknown")){throw malformed("phylum without domain");}
			status=domain.equals("unknown") ? "unknown" : phylum.equals("unknown") ? "partial" : "classified";
			name=name_; taxId=taxId_; ani=ani_;
		}
		final String status, domain, phylum, name;
		final long taxId;
		final double ani;
	}

	static final String QUERY_NAME="magqc_bin";
	static final String[] OPTIONS={"fasta", "binid", "pgmmode", "passes", "taxaddress", "taxdomain", "taxphylum", "profile", "profilesha80",
		"roster", "ref", "rolemanifest", "core", "coveringsets", "sidecar",
		"hbmbundle", "hbmprovenance", "policy", "lookahead"};
}
