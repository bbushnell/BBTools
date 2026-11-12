package prok;

import java.util.ArrayList;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaTextIO.Input;

/** Imports native CM observations without admitting biological reference loci.
 * Reuses RrnaCmTruthBuilder's byte tokenizer and EukPanelCmEvidence's inclusive
 * overlap arithmetic. Model roles come from the reviewed roster, never an old
 * caller-family whitelist. All hits and zero-hit search counts survive import.
 * The caller supplies source-bound tables and retains the original native files.
 * Current FASTA observations do not prove the bytes used by an earlier search.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceEvidence{

	/*--------------------------------------------------------------*/
	/*----------------         Model roster         ----------------*/
	/*--------------------------------------------------------------*/

	/** Joins reviewed roles, captured model metadata and the exact searched query roster.
	 * Declared model hashes are preserved, not recomputed or independently certified here.
	 */
	public static Roster loadModels(final String scope, final String assets, final String queries){
		final ObjectMap<String,Scope> roles=new ObjectMap<String,Scope>(String.class, Scope.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(scope)){
			RrnaResourceIO.header(input.reader, SCOPE_HEADER, scope);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Expected six reviewed model-scope fields");
				final String accession=p.parseString(0), label=p.parseString(1), searchRole=p.parseString(2);
				final String scoringScope=p.parseString(3), family=p.parseString(4), role=p.parseString(5);
				require(validAccession(accession) && !label.isEmpty() && !searchRole.isEmpty()
					&& !scoringScope.isEmpty() && family.matches("[A-Za-z0-9_.]+"), "Incomplete model scope: "+accession);
				require(role.equals("TARGET_CLASS") || role.equals("RELATED_CLASS"), "Unknown evidence role: "+role);
				require(!roles.contains(accession), "Repeated model scope: "+accession);
				roles.put(accession, new Scope(label, searchRole, scoringScope, family, role));
			}
		}
		final Roster roster=new Roster();
		try(Input input=new Input(assets)){
			RrnaResourceIO.header(input.reader, ASSET_HEADER, assets);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Expected six captured-model fields");
				final String accession=p.parseString(0), release=p.parseString(1);
				final int length=p.parseInt(2);
				final double ga=p.parseDouble(3);
				final String hash=p.parseString(4), path=p.parseString(5);
				final Scope role=roles.get(accession);
				require(role!=null && !roster.models.contains(accession), "Unknown or repeated captured model: "+accession);
				require(!release.isEmpty() && length>0 && Double.isFinite(ga)
					&& SeedOffsetTable.validSha80(hash) && path.startsWith("/"), "Invalid captured model metadata: "+accession);
				if(roster.release==null){roster.release=release;}
				require(roster.release.equals(release), "One evidence generation must not silently mix Rfam releases");
				roster.models.put(accession, new Model(accession, role, release, length, ga, hash, path));
			}
		}
		require(roles.size()>0 && roles.size()==roster.models.size(), "Every reviewed model needs captured metadata");
		try(Input input=new Input(queries)){
			RrnaResourceIO.header(input.reader, QUERY_HEADER, queries);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==3, "Expected accession, query name and CLEN");
				final Model model=roster.models.get(p.parseString(0));
				require(model!=null && model.queryIndex<0 && p.parseInt(2)==model.length,
					"Queried model must occur once and agree with captured CLEN");
				model.queryName=p.parseString(1);
				require(!model.queryName.isEmpty(), "Native queries require their original model names");
				model.queryIndex=roster.queries.size();
				roster.queries.add(model);
			}
		}
		require(!roster.queries.isEmpty(), "The searched roster cannot be empty");
		return roster;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Raw import          ----------------*/
	/*--------------------------------------------------------------*/

	/** Imports one completed whole-genome search, preserving unresolved target/range observations.
	 * No coordinate clipping, accession aliases, score filtering, locus merging or admission occurs.
	 */
	public static Evidence read(final String table, final String completion, final Observation observed, final Roster roster){
		require(observed!=null && roster!=null, "Raw evidence requires observed contigs and an explicit model roster");
		final Evidence evidence=new Evidence(table, observed, roster);
		readCompletion(completion, evidence);
		final RrnaCmTruthBuilder.Fields fields=new RrnaCmTruthBuilder.Fields();
		boolean program=false, version=false, finished=false;
		long number=0;
		Model previousModel=null;
		String previousTarget=null;
		try(Input input=new Input(table)){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				number++;
				if(row.length==0){continue;}
				require(!finished, "Content follows the native completion footer: "+table+":"+number);
				if(row[0]=='#'){
					final String text=Euk18sShredCmSummary.text(row);
					if(text.startsWith("# Program:")){
						require(text.substring(10).trim().equals("cmsearch"), "Expected native cmsearch output: "+table);
						program=true;
					}else if(text.startsWith("# Version:")){
						require(text.substring(10).trim().startsWith("1.1.5 "), "Expected the recorded Infernal 1.1.5 format: "+table);
						version=true;
					}else if(text.equals("# [ok]")){finished=true;}
					continue;
				}
				fields.set(row);
				require(fields.eq(4, "cm"), "HMM-only rows cannot be imported as CM evidence: "+table+":"+number);
				if(previousModel==null || !fields.eq(3, previousModel.accession)){
					previousModel=roster.models.get(fields.text(3));
				}
				require(previousModel!=null && previousModel.queryIndex>=0 && fields.eq(2, previousModel.queryName),
					"Raw hit uses an unknown, unqueried or differently named model: "+table+":"+number);
				if(previousTarget==null || !fields.eq(0, previousTarget)){previousTarget=fields.text(0);}
				final Hit hit=new Hit(table, number, row, previousTarget, previousModel, fields, observed);
				evidence.hits.add(hit);
				evidence.observedCounts[previousModel.queryIndex]++;
			}
		}
		require(program && version && finished, "Missing native program/version/completion footer: "+table);
		for(int i=0; i<roster.queries.size(); i++){
			require(evidence.observedCounts[i]==evidence.expectedCounts[i],
				"Native completion and retained hits disagree for "+observed.genome.id+" "+roster.queries.get(i).accession
				+": "+evidence.observedCounts[i]+"/"+evidence.expectedCounts[i]);
		}
		return evidence;
	}

	/** Completion accounting includes models with no hits and the exact two-strand search space. */
	private static void readCompletion(final String path, final Evidence evidence){
		assert(evidence.hits.isEmpty()) : "Completion denominators are loaded before any raw hit is retained";
		final LineParser1 p=new LineParser1('\t');
		final boolean[] seen=new boolean[evidence.roster.queries.size()];
		final NcrnaEvidenceCatalog.Genome genome=evidence.observed.genome;
		final long searched=Math.multiplyExact(genome.bases, 2);
		int count=0;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, COMPLETION_HEADER, path);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==7 && p.termEquals(genome.id, 0) && p.termEquals("COMPLETE", 6),
					"Completion rows must identify this input and a successful native search");
				final Model model=evidence.roster.models.get(p.parseString(1));
				require(model!=null && model.queryIndex>=0 && !seen[model.queryIndex], "Unknown or repeated completion model");
				final long hits=p.parseLong(5);
				require(p.parseLong(2)==genome.records && p.parseLong(3)==genome.bases
					&& p.parseLong(4)==searched && hits>=0, "Native search dimensions disagree with the input census");
				evidence.expectedCounts[model.queryIndex]=hits;
				seen[model.queryIndex]=true;
				count++;
			}
		}
		require(count==seen.length, "Missing populated or zero-hit model completion rows: "+path);
	}

	/** Deterministic support order for hits already associated with the same physical historical locus.
	 * It ranks evidence; it never chooses a biological winner or merges features by itself.
	 */
	public static int compareSupport(final Hit a, final Hit b, final long start, final long stop){
		require(start>0 && stop>=start && a.genome.id.equals(b.genome.id) && a.genome.fasta.equals(b.genome.fasta)
			&& a.target.equals(b.target) && a.strand==b.strand,
			"Support ranking requires one input, literal target, strand and valid historical interval");
		int order=Long.compare(EukPanelCmEvidence.overlap(b.start, b.stop, start, stop),
			EukPanelCmEvidence.overlap(a.start, a.stop, start, stop));
		if(order==0){order=Long.compare((long)b.modelSpan*a.model.length, (long)a.modelSpan*b.model.length);}
		if(order==0){order=Double.compare(b.bits, a.bits);}
		if(order==0){order=Double.compare(a.evalue, b.evalue);}
		if(order==0){order=a.model.accession.compareTo(b.model.accession);}
		if(order==0){order=Integer.compare(a.start, b.start);}
		if(order==0){order=Integer.compare(a.stop, b.stop);}
		if(order==0){order=a.table.compareTo(b.table);}
		return order==0 ? Long.compare(a.line, b.line) : order;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Evidence records      ----------------*/
	/*--------------------------------------------------------------*/

	/** One native hit. Raw bytes and model restrictions remain available for later adjudication. */
	public static final class Hit{
		private Hit(final String table_, final long line_, final byte[] raw_, final String target_, final Model model_,
			final RrnaCmTruthBuilder.Fields f, final Observation observed){
			this(table_, line_, raw_, target_, target_, model_, f, observed, false);
		}

		/** Batch callers return hits only after verifying the native original-GA option binding. */
		Hit(final String table_, final long line_, final byte[] raw_, final String target_, final String nativeTarget_, final Model model_,
			final RrnaCmTruthBuilder.Fields f, final Observation observed, final boolean boundInclusion_){
			table=table_;
			genome=observed.genome;
			line=line_;
			raw=raw_;
			target=target_;
			nativeTarget=nativeTarget_;
			boundInclusion=boundInclusion_;
			model=model_;
			modelFrom=f.number(5);
			modelTo=f.number(6);
			seqFrom=f.number(7);
			seqTo=f.number(8);
			start=Math.min(seqFrom, seqTo);
			stop=Math.max(seqFrom, seqTo);
			strand=f.eq(9, "+") ? '+' : f.eq(9, "-") ? '-' : '?';
			truncation=f.eq(10, "no") ? 0 : f.eq(10, "5'") ? 1 : f.eq(10, "3'") ? 2 : f.eq(10, "5'&3'") ? 3 : -1;
			included=f.eq(16, "!");
			bits=Parse.parseDouble(f.bytes, f.start[14], f.end[14]);
			evalueText=f.text(15);
			evalue=Double.parseDouble(evalueText);
			if(!(modelFrom>0 && modelTo>=modelFrom && modelTo<=model.length && start>0)){
				throw new IllegalArgumentException("Invalid CM or sequence coordinates at native row "+line);
			}
			if(!(strand=='+' && seqFrom<=seqTo || strand=='-' && seqFrom>=seqTo)){
				throw new IllegalArgumentException("Native strand and sequence coordinate order disagree at row "+line);
			}
			if(!(truncation>=0 && (included || f.eq(16, "?")))){
				throw new IllegalArgumentException("Unknown native truncation/inclusion flag at row "+line);
			}
			if(!(Double.isFinite(bits) && Double.isFinite(evalue) && evalue>=0)){
				throw new IllegalArgumentException("Invalid score/E-value at native row "+line);
			}
			modelSpan=modelTo-modelFrom+1;
			contigLength=observed.length(target);
			range=contigLength<0 ? "UNRESOLVED_LITERAL_TARGET" : stop>contigLength ? "OUT_OF_RANGE" : "VALID_INTERVAL";
			fullModelSpan=modelFrom==1 && modelTo==model.length;
			atLeftEdge=range.equals("VALID_INTERVAL") && start==1;
			atRightEdge=range.equals("VALID_INTERVAL") && stop==contigLength;
			// Legacy observations preserve their historical printed-score basis. Batch
			// observations use a verified original-GA inclusion binding because Infernal
			// prints only one decimal. Neither basis alone establishes biological truth.
			//TODO: The unbound legacy path must not be used for exact sub-GA boundary decisions.
			gaPass=boundInclusion ? included : bits>=model.ga;
		}

		public final long line;
		/** Retained native line; consumers must not mutate it. */
		public final byte[] raw;
		public final String table, target, nativeTarget, range, evalueText;
		public final NcrnaEvidenceCatalog.Genome genome;
		public final Model model;
		public final int modelFrom, modelTo, modelSpan, seqFrom, seqTo, start, stop, contigLength;
		public final char strand;
		/** Native flags: 0=no truncation, 1=5-prime, 2=3-prime, 3=both; not a biological fragment label. */
		public final int truncation;
		public final double bits, evalue;
		public final boolean included, boundInclusion, fullModelSpan, atLeftEdge, atRightEdge, gaPass;
	}

	/** Per-input observations only; no admitted-locus count or truth flag is available. */
	public static final class Evidence{
		Evidence(final String table_, final Observation observed_, final Roster roster_){
			table=table_;
			observed=observed_;
			roster=roster_;
			expectedCounts=new long[roster.queries.size()];
			observedCounts=new long[roster.queries.size()];
		}
		public int size(){return hits.size();}
		public Hit get(final int index){return hits.get(index);}
		public long hitsForQuery(final int index){return observedCounts[index];}
		/** Adds a derived per-input observation after its batch target has been resolved.
		 * Native expected counts are validated for the whole batch, never fabricated
		 * by copying these per-input observed counts into the legacy expected array.
		 */
		void addBatchHit(final Hit hit){
			require(table==null && hit.genome==observed.genome && hit.boundInclusion, "Batch evidence requires resolved input identity and bound inclusion");
			hits.add(hit);
			observedCounts[hit.model.queryIndex]++;
		}
		/** Legacy single-genome native table; null for explicitly batched evidence. */
		public final String table;
		public final Observation observed;
		public final Roster roster;
		private final ArrayList<Hit> hits=new ArrayList<Hit>();
		private final long[] expectedCounts, observedCounts;
	}

	/** Complete captured model universe, with a distinct exact list of searched queries. */
	public static final class Roster{
		private Roster(){}
		public int size(){return models.size();}
		public int queried(){return queries.size();}
		public Model query(final int index){return queries.get(index);}
		public Model model(final String accession){return models.get(accession);}
		private final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class, Model.class);
		private final ArrayList<Model> queries=new ArrayList<Model>();
		private String release;
	}

	/** Declared class/structure and scoring restrictions are never discarded after normalization. */
	public static final class Model{
		private Model(final String accession_, final Scope scope_, final String release_, final int length_,
			final double ga_, final String hash_, final String path_){
			accession=accession_;
			scope=scope_;
			release=release_;
			length=length_;
			ga=ga_;
			sha80=hash_;
			path=path_;
		}
		public boolean relatedClass(){return scope.role.equals("RELATED_CLASS");}
		String queryName(){return queryName;}
		int queryIndex(){return queryIndex;}
		public final String accession, release, sha80, path;
		public final int length;
		public final double ga;
		public final Scope scope;
		private String queryName;
		private int queryIndex=-1;
	}

	public static final class Scope{
		private Scope(final String label_, final String search_, final String scoring_, final String family_, final String role_){
			label=label_;
			searchRole=search_;
			scoringScope=scoring_;
			family=family_;
			role=role_;
		}
		public final String label, searchRole, scoringScope, family, role;
	}

	private static boolean validAccession(final String value){return value.matches("RF[0-9]{5}");}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private NcrnaReferenceEvidence(){}

	public static final String SCOPE_HEADER="CM_accession\tbiological_label\tsearch_role\tscoring_scope\treference_family\tevidence_role";
	public static final String ASSET_HEADER="accession\trelease\tCLEN\tGA_bits\tsha80\tpath";
	public static final String QUERY_HEADER="accession\tname\tCLEN";
	public static final String COMPLETION_HEADER="input_id\tmodel\tFASTA_records\tinput_bases\tsearched_residues\treported_hits\tstatus";
}
