package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaSupport.Graph;

/** Complete original-locus or admitted-occurrence population for boundary learning.
 * Joins the scored original graph, genome roster and locus ledger before
 * grouping TaxIDs. Tuning partitions groups; fixed refit withholds none and
 * requires no feasible validation quota. Partial/excluded loci retain their original
 * model bins and split membership; eligibility only controls subsequent learning.
 * This input-only helper uses the ByteFile/LineParser idiom from the text template.
 * It retains primitive locus columns, not the original TSV lines or gene objects.
 * Returned arrays are owned here and must be treated as read-only by consumers.
 * A distinct occurrence-ledger header supplies observed coordinates without
 * historical annotation flags. Its admitted records remain eligible unless an
 * explicit exclusions file names them; unavailable sequence features are reported
 * separately by the downstream tables/vector tools. Legacy annotation behavior
 * is unchanged when the original 21-column header is supplied.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningPopulation{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads and validates all metadata before deriving a split or eligibility.
	 * @param graph Complete ORIGINAL graph, including below-cutoff winners
	 * @param originalNames Names in the supplied original consensus FASTA order
	 * @param roster Genome roster, including genomes with zero target loci
	 * @param loci Complete original annotation ledger
	 * @param exclusions Optional file of one original locus ID per line
	 * @param seed Whole-TaxID split seed
	 */
	public static NcrnaLearningPopulation load(final Graph graph, final String[] originalNames,
			final String roster, final String loci, final String exclusions, final long seed){
		return load(graph, originalNames, roster, loci, exclusions, seed, false);
	}

	/** Loads all eligible records without requesting a validation partition.
	 * Source identity/count checks and explicit partial/exclusion handling remain
	 * identical to tuning, including corpora with fewer than ten TaxID groups.
	 */
	public static NcrnaLearningPopulation loadForRefit(final Graph graph, final String[] originalNames,
			final String roster, final String loci, final String exclusions){
		return load(graph, originalNames, roster, loci, exclusions, 0, true);
	}

	private static NcrnaLearningPopulation load(final Graph graph, final String[] originalNames,
			final String roster, final String loci, final String exclusions, final long seed, final boolean refit){
		if(graph==null || graph.reduced || !Arrays.equals(graph.modelNames, originalNames)){
			throw new IllegalArgumentException("Coverage bins require an ORIGINAL graph on the supplied original library axis");
		}
		final NcrnaLearningPopulation population=new NcrnaLearningPopulation(graph);
		population.readRoster(roster);
		population.readLoci(loci);
		population.partition(seed, refit);
		population.readExclusions(exclusions);
		return population;
	}

	private NcrnaLearningPopulation(final Graph graph){
		ids=graph.queryNames.clone();
		models=graph.modelNames.clone();
		originalModel=graph.bestModels.clone();
		final int count=ids.length;
		source=new int[count];
		Arrays.fill(source, -1);
		length=new int[count];
		leftFlank=new int[count];
		rightFlank=new int[count];
		partial=new boolean[count];
		excluded=new boolean[count];
		for(int gene=0; gene<count; gene++){
			assert(originalModel[gene]>=0 && originalModel[gene]<models.length)
				: "Original graphs retain an overall winner even below the support cutoff";
			byLocus.put(ids[gene], gene);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Population Join      ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps one source object per genome; phylum identifiers include domain. */
	private void readRoster(final String path){
		final ArrayList<Source> rows=new ArrayList<Source>();
		final ObjectIntMap<String> assemblyTaxid=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> phyla=new ObjectIntMap<String>(String.class);
		try(Input input=new Input(path, ROSTER_HEADER, 11)){
			while(input.next()){
				if(input.number(0)!=rows.size()){throw input.bad("Roster indices must be dense and ordered");}
				final String name=input.text(1), assembly=input.text(2);
				final int taxid=input.number(3), expected=input.number(8);
				if(taxid==0 || bySource.contains(name)){throw input.bad("Source needs a positive TaxID and unique genome name: "+name);}
				final int prior=assemblyTaxid.get(assembly);
				if(prior>=0 && prior!=taxid){throw input.bad("Conflicting TaxIDs for assembly "+assembly);}
				assemblyTaxid.put(assembly, taxid);
				final String domain=input.text(4), phylum=input.text(5);
				final String key=domain+"\t"+phylum;
				int bin=phyla.get(key);
				if(bin<0){
					bin=phyla.size();
					phyla.put(key, bin);
				}
				bySource.put(name, rows.size());
				rows.add(new Source(name, assembly, domain, phylum, taxid, expected, bin));
			}
		}
		if(rows.isEmpty()){throw new IllegalArgumentException("Learning requires a nonempty original genome roster");}
		sources=rows.toArray(new Source[rows.size()]);
		phylumCount=phyla.size();
	}

	/** Checks source identity and count conservation for every original graph query. */
	private void readLoci(final String path){
		final int[] observed=new int[sources.length];
		int count=0;
		try(Input input=new Input(path, LOCUS_HEADER, 21, true)){
			while(input.next()){
				final LineParser1 p=input.parser;
				final String id=input.text(0);
				final int shift=input.occurrenceRows ? 1 : 0;
				final int gene=byLocus.get(id), genome=bySource.get(input.text(1+shift));
				if(gene<0 || genome<0 || source[gene]>=0){throw input.bad("Unknown or duplicate original locus/source: "+id);}
				final Source s=sources[genome];
				if(!p.termEquals(s.assembly, 2+shift) || input.number(3+shift)!=s.taxid
						|| !p.termEquals(s.domain, 4+shift) || !p.termEquals(s.phylum, 5+shift)){
					throw input.bad("Locus identity disagrees with its original source: "+id);
				}
				if(input.occurrenceRows){
					readOccurrence(input, gene, id);
				}else{
					partial[gene]=p.termEquals("SOURCE_PARTIAL", 14);
					if(!partial[gene] && !p.termEquals("ANNOTATION_FULL_COPY", 14)){
						throw input.bad("Unknown annotation status for "+id);
					}
					if(!p.termEquals(partial[gene] ? "true" : "false", 12)){
						throw input.bad("Partial flag disagrees with annotation status for "+id);
					}
					if(!p.termEquals("false", 13)){throw input.bad("Resolve pseudogene status in the corpus before learning: "+id);}
					if(!p.termEquals("+", 7) && !p.termEquals("-", 7)){throw input.bad("Locus needs an explicit strand: "+id);}
					length[gene]=input.number(9);
					final int contig=input.number(15), ambiguous=input.number(17);
					leftFlank[gene]=input.number(18);
					rightFlank[gene]=input.number(19);
					if(length[gene]==0 || length[gene]>contig || ambiguous>length[gene]){
						throw input.bad("Gene length or ambiguity count does not fit the contig: "+id);
					}
					if(!p.termEquals("true", 16) && !p.termEquals("false", 16)){throw input.bad("Invalid physical-edge flag: "+id);}
				}
				source[gene]=genome;
				observed[genome]++;
				count++;
			}
		}
		if(count!=ids.length){throw new IllegalArgumentException("Original graph has "+ids.length+" loci but the ledger joins "+count);}
		for(int genome=0; genome<sources.length; genome++){
			if(observed[genome]!=sources[genome].expected){
				throw new IllegalArgumentException("Roster/locus count mismatch for "+sources[genome].name+
					": expected "+sources[genome].expected+", observed "+observed[genome]);
			}
		}
	}

	/** Admitted occurrences have their own coordinate contract, not invented annotation flags.
	 * Every admitted record stays eligible unless an explicit exclusions file names it.
	 * Original annotation partial/pseudogene evidence stays in the admission provenance.
	 */
	private void readOccurrence(final Input input, final int gene, final String id){
		assert(source[gene]<0) : "Occurrence identity is unique and is assigned to its source after this validation";
		final LineParser1 p=input.parser;
		final String occurrence=input.text(1);
		if(!id.equals("ID="+occurrence)){throw input.bad("Record ID must match the admitted occurrence's CutGff identifier: "+id);}
		if(!p.termEquals("+", 8) && !p.termEquals("-", 8)){throw input.bad("Occurrence needs an explicit strand: "+id);}
		final int start=input.number(9), stop=input.number(10), contig=input.number(12);
		length[gene]=input.number(11);
		leftFlank[gene]=input.number(13);
		rightFlank[gene]=input.number(14);
		if(start<1 || stop<start || stop>contig || (long)stop-start+1!=length[gene]){
			throw input.bad("Occurrence length and coordinates disagree: "+id);
		}
		final boolean minus=p.termEquals("-", 8);
		final int availableLeft=minus ? contig-stop : start-1, availableRight=minus ? start-1 : contig-stop;
		if(leftFlank[gene]>availableLeft || rightFlank[gene]>availableRight){
			throw input.bad("Sense-oriented occurrence flanks exceed the source contig: "+id);
		}
	}

	/** Establishes source groups before eligibility; refit deliberately withholds none. */
	private void partition(final long seed, final boolean refit){
		final int[] taxa=new int[ids.length], phyla=new int[ids.length];
		for(int gene=0; gene<ids.length; gene++){
			assert(source[gene]>=0) : "Every original query must join a source before forming coverage quotas";
			final Source s=sources[source[gene]];
			taxa[gene]=s.taxid;
			phyla[gene]=s.bin;
		}
		split=refit ? NcrnaLearningSplit.all(taxa, phyla, originalModel, phylumCount, models.length)
			: NcrnaLearningSplit.split(taxa, phyla, originalModel, phylumCount, models.length, seed);
	}

	/** Exclusions remove teaching examples, never original population members. */
	private void readExclusions(final String path){
		assert(split!=null) : "Learning exclusions follow complete source grouping and cannot alter its population";
		if(path==null){return;}
		try(Input input=new Input(path, null, 1)){
			while(input.next()){
				final String id=input.text(0);
				final int gene=byLocus.get(id);
				if(gene<0 || excluded[gene]){throw input.bad("Unknown or duplicate learning exclusion: "+id);}
				excluded[gene]=true;
			}
		}
	}

	/** Exact locus lookup; no truncation or guessing of a supplied identifier. */
	public int index(final String id){return byLocus.get(id);}
	public boolean eligible(final int gene){return !partial[gene] && !excluded[gene];}
	public NcrnaLearningSplit.Result split(){return split;}
	public String sourceName(final int gene){return sources[source[gene]].name;}

	/*--------------------------------------------------------------*/
	/*----------------          Input Owner         ----------------*/
	/*--------------------------------------------------------------*/

	/** Serial byte reader with reusable field offsets; close preserves parse failures. */
	private static final class Input implements AutoCloseable{
		Input(final String path_, final String header, final int terms_){
			this(path_, header, terms_, false);
		}
		Input(final String path_, final String header, final int terms_, final boolean allowOccurrences_){
			path=path_;
			terms=terms_;
			allowOccurrences=allowOccurrences_;
			owner=new NcrnaTextIO.Input(path);
			file=owner.reader;
			// Constructor cannot throw after opening: header parsing is deferred to next().
			pendingHeader=header;
		}

		boolean next(){
			if(pendingHeader!=null){
				final byte[] header=file.nextLine();
				line++;
				final String text=header==null ? null : new String(header, StandardCharsets.UTF_8);
				if(allowOccurrences && OCCURRENCE_HEADER.equals(text)){
					occurrenceRows=true;
					terms=15;
				}else if(!pendingHeader.equals(text)){
					throw bad("Unexpected metadata header");
				}
				pendingHeader=null;
			}
			final byte[] row=file.nextLine();
			if(row==null){return false;}
			line++;
			parser.set(row);
			if(parser.terms()!=terms){throw bad("Expected "+terms+" tab-separated fields");}
			return true;
		}

		/** Decimal count parser rejects sign, overflow and missing values before narrowing. */
		int number(final int term){
			parser.setBounds(term);
			try{
				return NcrnaNumbers.unsignedInt(parser.line(), parser.a(), parser.b());
			}catch(IllegalArgumentException failure){
				throw bad("Invalid integer at field "+term+": "+failure.getMessage());
			}
		}

		String text(final int term){
			final String value=parser.parseString(term);
			if(value.isEmpty()){throw bad("Empty identifier at field "+term);}
			return value;
		}

		IllegalArgumentException bad(final String message){return new IllegalArgumentException(path+":"+line+": "+message);}

		@Override
		public void close(){
			owner.close();
		}

		final String path;
		int terms;
		final boolean allowOccurrences;
		boolean occurrenceRows;
		final NcrnaTextIO.Input owner;
		final ByteFile file;
		final LineParser1 parser=new LineParser1('\t');
		String pendingHeader;
		long line;
	}

	/** Persisted source identity shared by every locus from that source. */
	private static final class Source{
		Source(final String name_, final String assembly_, final String domain_, final String phylum_, final int taxid_, final int expected_, final int bin_){
			name=name_;
			assembly=assembly_;
			domain=domain_;
			phylum=phylum_;
			taxid=taxid_;
			expected=expected_;
			bin=bin_;
		}
		final String name, assembly, domain, phylum;
		final int taxid, expected, bin;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final String[] ids, models;
	public final int[] originalModel, source, length, leftFlank, rightFlank;
	/** Historical SOURCE_PARTIAL exclusion flag; occurrence metadata makes no such annotation claim. */
	public final boolean[] partial;
	/** Explicit caller-supplied teaching exclusions, preserving the complete population. */
	public final boolean[] excluded;
	private final ObjectIntMap<String> byLocus=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> bySource=new ObjectIntMap<String>(String.class);
	private Source[] sources;
	private int phylumCount;
	private NcrnaLearningSplit.Result split;

	public static final String ROSTER_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs";
	public static final String LOCUS_HEADER="locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines";
	public static final String OCCURRENCE_HEADER="record_id\toccurrence_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tstart1\tstop1\tlength\tcontig_length\tlflank\trflank";
}
