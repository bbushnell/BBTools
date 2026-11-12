package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import map.ObjectMap;
import map.ObjectSet;
import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.DoubleList;
import shared.Tools;

/** Joins declared prior decisions to current saved model evidence without admission.
 * Prior accepted cores, accepted whole intervals and construction-eligible cores
 * are distinct evidence bases. Missing hits, multiple own-model hits, exact repeated
 * anchors and paths outside the declared input namespace remain explicit. Only a
 * single own-model hit passing the existing reciprocal-half overlap rule contributes
 * a score/coverage observation; weaker overlaps remain explicit. Other
 * model overlaps are retained separately, not treated as resolved alternatives.
 *
 * Quantiles describe observed hits, not the uncensored positive distribution.
 * In particular, thresholded input cannot establish the missing low-score tail.
 * No threshold is chosen, no prior decision is revoked, and no GFF is written.
 * Native parsing/output and primitive numeric lists; one thread owns all state.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceAnchors{

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads saved observations, joins declared anchors, then closes all four outputs. */
	static void process(final String input, final String manifest, final String prefix,
		final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceAnchors report=new NcrnaReferenceAnchors(input, catalog, roster);
		report.readManifest(manifest);
		report.write(prefix);
	}

	private NcrnaReferenceAnchors(final String input, final NcrnaEvidenceCatalog catalog_, final Roster roster_){
		catalog=catalog_;
		roster=roster_;
		geometry=NcrnaReferenceGeometry.load(input, catalog, roster);
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			final InputState state=new InputState(genome, i);
			byId.put(genome.id, state);
			addPath(genome.fasta, state);
			for(int j=0; j<genome.memberships(); j++){addPath(genome.membership(j).suppliedFasta, state);}
		}
		for(final Hit hit:geometry.hits()){byId.get(hit.input.genome.id).hits.add(hit);}
		for(int i=0; i<roster.queried(); i++){model(roster.query(i).accession);}
		readContigs(input+".contigs.tsv.gz");
	}
	private void addPath(final String path, final InputState state){
		final InputState previous=byPath.get(path);
		require(previous==null || previous==state, "Declared FASTA path maps to multiple inputs: "+path);
		byPath.put(path, state);
	}
	/** Validates literal contigs even when the current scan has no hit there. */
	private void readContigs(final String path){
		final LineParser1 p=new LineParser1('\t');
		InputState state=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, "input_id\tfasta\tliteral_target\tobserved_length", path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==4, "Saved contig table needs four fields");
				if(state==null || !p.termEquals(state.genome.id, 0)){state=byId.get(p.parseString(0));}
				require(state!=null && geometry.profiles().selected(state.genome.id) && same(p, 1, state.path),
					"Saved contig must bind its selected input and exact FASTA");
				probe.set(p, 2);
				final int length=p.parseInt(3);
				require(probe.to>probe.from && length>0 && !state.lengths.contains(probe), "Saved contigs must be unique and have positive lengths");
				state.lengths.put(probe.copy(), length);
				state.bases=Math.addExact(state.bases, length);
			}
		}
		for(final InputState item:byId.values()){
			if(item==null || !geometry.profiles().selected(item.genome.id)){continue;}
			require(item.lengths.size()==item.genome.records && item.bases==item.genome.bases,
				"Saved contig observations must conserve each selected input's dimensions");
		}
	}
	private void checkContig(final InputState input, final int declaredLength){
		probe.set(key.array, 0, key.length());
		require(input.lengths.get(probe)==declaredLength && declaredLength>0,
			"Prior literal target or contig length differs from the saved input observations");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Prior sources        ----------------*/
	/*--------------------------------------------------------------*/

	private void readManifest(final String path){
		final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Anchor manifest needs six fields");
				final String name=text(p, 0), format=text(p, 1), file=text(p, 2), bindings=text(p, 3), authority=text(p, 4);
				require(!name.isEmpty() && names.add(name) && file.startsWith("/") && !authority.isEmpty()
					&& p.parseLong(5)>=0, "Anchor sources need unique names, absolute files, authority and exact row counts");
				require(format.equals(BOUND) && bindings.equals("NA") || format.equals(TPP) && bindings.startsWith("/"),
					"Use bound_prior_v1 with bindings=NA, or tpp_construction_v1 with its original genome roster");
				require(Tools.testInputFiles(false, true, file), "Declared prior decisions must be readable before output");
				if(format.equals(TPP)){require(Tools.testInputFiles(false, true, bindings), "Declared TPP bindings must be readable before output");}
				sources.add(new Source(name, format, file, bindings, authority, p.parseLong(5)));
			}
		}
		require(!sources.isEmpty(), "Declare the prior sources explicitly");
	}

	/** Existing bound rows retain their original decision-file/line and authorities. */
	private void readBound(final Source source, final Output rows, final Output links){
		final LineParser1 p=new LineParser1('\t');
		InputState input=null;
		try(Input in=new Input(source.file)){
			RrnaResourceIO.header(in.reader, BOUND_HEADER, source.file);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				require(p.terms()==21, "Bound prior table needs twenty-one fields");
				require(p.parseLong(2)>=2 && !p.termEquals("NA", 18), "Bound rows require their original physical line and declared authority");
				if(input==null || !p.termEquals(input.genome.id, 6)){input=byId.get(p.parseString(6));}
				require(input!=null && geometry.profiles().selected(input.genome.id)
					&& same(p, 7, input.path), "Bound prior must retain its exact selected census input/path");
				final String accession=p.parseString(10);
				final ModelState model=model(accession);
				final int basis=p.termEquals(accession, 11) ? CORE : p.termEquals("ACCEPTED_NOVEL_FMN", 11) ? WHOLE : -1;
				require(basis>=0 || p.termEquals("UNRESOLVED", 11) || p.termEquals("UNREVIEWED", 11) || p.termEquals("B12clade", 11),
					"Unknown prior assignment cannot be silently classified nonpositive");
				require(basis!=WHOLE || accession.equals("RF00050"), "Accepted whole-interval format is specific to the prior FMN decision");
				final int start=basis==CORE ? p.parseInt(14) : basis==WHOLE ? p.parseInt(12) : 0;
				final int stop=basis==CORE ? p.parseInt(15) : basis==WHOLE ? p.parseInt(13) : 0;
				key.clear();
				copy(key, p, 8);
				checkContig(input, p.parseInt(9));
				final char strand=strand(p, 16);
				if(basis>=0){require(start>0 && stop>=start && stop<=p.parseInt(9), "Prior anchor leaves its recorded contig");}
				observe(source, line, input, model, basis, start, stop, strand, links);
				begin(b.clear(), source, line);
				copy(b, p, 1).tab().append(p.parseLong(2)).tab();
				field(b, source.authority).tab();
				copy(b, p, 18).tab();
				copy(b, p, 19).tab();
				copy(b, p, 5).tab();
				writeRow(rows, b, input, model, basis, start, stop, strand, p, 11, true);
			}
		}
	}

	/** Construction eligibility is retained as its own basis, never relabeled scored truth.
	 * Paths absent from the supplied census/memberships are reported, not guessed by
	 * assembly, TaxID or basename. All original source rows remain accounted for.
	 */
	private void readTpp(final Source source, final Output rows, final Output links){
		final ObjectMap<Target,Binding> bindings=new ObjectMap<Target,Binding>(Target.class, Binding.class);
		final Target probe=new Target();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(source.bindings)){
			RrnaResourceIO.header(in.reader, NcrnaCalibrationPrep.ROSTER_HEADER, source.bindings);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==11, "Original TPP genome roster needs eleven fields");
				probe.set(p, 1);
				require(probe.to>probe.from && !bindings.contains(probe), "Original TPP genome identifiers must be unique");
				final String fasta=text(p, 7);
				require(fasta.startsWith("/"), "TPP source FASTA must be absolute");
				bindings.put(probe.copy(), new Binding(fasta, byPath.get(fasta), p.parseByteArray(2)));
			}
		}
		final ModelState model=model("RF00059");
		try(Input in=new Input(source.file)){
			RrnaResourceIO.header(in.reader, TPP_HEADER, source.file);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				require(p.terms()==46, "Reviewed TPP disposition table needs forty-six fields");
				probe.set(p, 1);
				final Binding binding=bindings.get(probe);
				require(binding!=null, "TPP decision lacks its original genome-to-FASTA binding");
				require(same(p, 2, binding.assembly), "TPP decision assembly differs from its original genome roster");
				require(p.termEquals("true", 41) || p.termEquals("false", 41), "TPP eligibility must be literal true or false");
				final int basis=p.termEquals("true", 41) ? CONSTRUCTION : -1;
				final int start=basis<0 ? 0 : p.parseInt(39), stop=basis<0 ? 0 : p.parseInt(40);
				final char strand=strand(p, 7);
				if(basis>=0){require(start>0 && stop>=start && stop<=p.parseInt(15), "TPP selected interval leaves its original contig");}
				p.setBounds(6);
				BuildNcrnaReference.field(key.clear(), p.line(), p.a(), p.b());
				if(binding.input!=null && geometry.profiles().selected(binding.input.genome.id)){
					checkContig(binding.input, p.parseInt(15));
				}
				observe(source, line, binding.input, model, basis, start, stop, strand, links);
				begin(b.clear(), source, line);
				field(b, source.file).tab().append(line).tab();
				field(b, source.authority).append("\tNA\tNA\t");
				field(b, binding.fasta).tab();
				writeRow(rows, b, binding.input, model, basis, start, stop, strand, p, 42, false);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Current evidence       ----------------*/
	/*--------------------------------------------------------------*/

	private void observe(final Source source, final long line, final InputState input,
		final ModelState model, final int basis, final int start, final int stop,
		final char strand, final Output links){
		source.rows++;
		own=other=0;
		duplicate=null;
		evaluated=false;
		if(basis<0){
			source.notPositive++;
			status="NOT_DECLARED_POSITIVE";
			return;
		}
		final Stats stats=model.stats[basis];
		stats.rows++;
		if(input==null || !geometry.profiles().selected(input.genome.id)){
			source.outside++;
			stats.outside++;
			status="OUTSIDE_DECLARED_INPUT_PATHS";
			return;
		}
		for(final Anchor previous:input.anchors){
			if(previous.model==model && previous.basis==basis && previous.start==start && previous.stop==stop
				&& previous.strand==strand && same(previous.target, key)){
				duplicate=previous;
				break;
			}
		}
		if(duplicate!=null){
			source.duplicates++;
			stats.duplicates++;
			status="REPEATED_EXACT_ANCHOR";
			return;
		}
		input.anchors.add(new Anchor(source, line, model, basis, start, stop, strand, key.toBytes()));
		source.unique++;
		stats.unique++;
		if(!stats.inputs[input.index]){
			stats.inputs[input.index]=true;
			stats.genomes++;
		}
		if(model.query<0){
			stats.unsearched++;
			status="CURRENT_MODEL_NOT_SEARCHED";
			return;
		}
		Hit single=null;
		evaluated=true;
		for(final Hit hit:input.hits){
			if(!hit.valid() || hit.strand!=strand || hit.stop<start || hit.start>stop || !same(hit.target.bytes, key)){continue;}
			final boolean ownModel=hit.query==model.query;
			if(ownModel){
				own++;
				single=hit;
			}else{other++;}
			final Model observed=roster.query(hit.query);
			begin(b.clear(), source, line);
			b.append(input.genome.id).tab().append(hit.table).tab().append(hit.line).tab();
			b.append(observed.accession).tab().append(ownModel).tab().append(hit.start).tab().append(hit.stop).tab();
			b.append((long)Math.min(stop, hit.stop)-Math.max(start, hit.start)+1).tab();
			b.append(hit.bits, 6).tab().append(hit.bits-observed.ga, 6).tab();
			b.append((hit.to-hit.from+1)/(double)observed.length, 6).tab().append(hit.truncation).tab();
			links.writer.print(b.append(hit.left).tab().append(hit.right).nl());
		}
		if(other>0){stats.alternatives++;}
		if(own==0){
			stats.missing++;
			status="NO_CURRENT_OWN_MODEL_OVERLAP";
		}else if(own>1){
			stats.multiple++;
			status="MULTIPLE_CURRENT_OWN_MODEL_HITS";
		}else{
			assert(single!=null) : "One own-model overlap must retain its exact native witness";
			final long overlap=(long)Math.min(stop, single.stop)-Math.max(start, single.start)+1;
			if(overlap<NcrnaGffFamilyGrader.MIN_OVERLAP_FRAC*((long)stop-start+1)
				|| overlap<NcrnaGffFamilyGrader.MIN_OVERLAP_FRAC*single.length()){
				stats.weakOverlap++;
				status="SINGLE_OWN_HIT_WEAK_OVERLAP";
				return;
			}
			stats.margins.add(single.bits-model.model.ga);
			stats.coverage.add((single.to-single.from+1)/(double)model.model.length);
			if(single.full()){
				stats.full++;
				stats.fullMargins.add(single.bits-model.model.ga);
			}else{stats.partial++;}
			status="SINGLE_CURRENT_OWN_MODEL_HIT";
		}
	}

	private ModelState model(final String accession){
		ModelState result=byModel.get(accession);
		if(result==null){
			final Model captured=roster.model(accession);
			require(captured!=null, "A declared prior accession lacks current model metadata: "+accession);
			int query=-1;
			for(int i=0; i<roster.queried(); i++){if(roster.query(i)==captured){query=i;}}
			result=new ModelState(captured, query, catalog.size());
			byModel.put(accession, result);
			models.add(result);
		}
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		try(Output rows=new Output(prefix+SUFFIXES[0]); Output links=new Output(prefix+SUFFIXES[1])){
			rows.writer.println("source\tsource_file\tsource_line\toriginal_file\toriginal_line\tdeclared_authority\toriginal_authority\toriginal_row_authority\toriginal_fasta\tinput_id\ttarget\taccession\tassignment\tbasis\tanchor_start1\tanchor_stop1\tstrand\tduplicate_source\tduplicate_line\tcurrent_own_hits\tcurrent_other_hits\tstatus");
			links.writer.println("source\tsource_file\tsource_line\tinput_id\tnative_table\tnative_line\taccession\town_model\tstart1\tstop1\toverlap_bases\tbits\tGA_margin_bits\tmodel_coverage\tnative_truncation\tleft_edge\tright_edge");
			for(final Source source:sources){
				if(source.format.equals(BOUND)){readBound(source, rows, links);}
				else{readTpp(source, rows, links);}
				require(source.rows==source.expected && source.rows==source.notPositive+source.outside+source.duplicates+source.unique,
					"Every declared prior row must be accounted for exactly once: "+source.name);
			}
		}
		writeModels(prefix+SUFFIXES[2]);
		try(Output out=new Output(prefix+SUFFIXES[3])){
			out.writer.println("source\trows\tnot_positive_anchor\toutside_input_paths\trepeated_exact_anchor\tunique_scoped_anchors\tinputs_observed\tcatalog_inputs\tevidence_status\tstatus");
			for(final Source source:sources){
				field(b.clear(), source.name).tab().append(source.rows).tab().append(source.notPositive).tab();
				out.writer.print(b.append(source.outside).tab().append(source.duplicates).tab().append(source.unique)
					.tab().append(geometry.profiles().observed()).tab().append(catalog.size()).tab().append(geometry.profiles().evidenceState())
					.append("\tPRIOR_EVIDENCE_JOIN_NOT_ADMISSION\n"));
			}
		}
	}
	private void writeRow(final Output out, final ByteBuilder b, final InputState input, final ModelState model,
		final int basis, final int start, final int stop, final char strand, final LineParser1 p,
		final int assignment, final boolean encoded){
		b.append(input==null ? "NA" : input.genome.id).tab().append(key).tab().append(model.model.accession).tab();
		if(encoded){copy(b, p, assignment);}
		else{
			p.setBounds(assignment);
			BuildNcrnaReference.field(b, p.line(), p.a(), p.b());
		}
		b.tab().append(basis<0 ? "NOT_POSITIVE" : BASES[basis]).tab().append(start).tab().append(stop).tab().append(strand).tab();
		if(duplicate==null){b.append("NA\tNA");}
		else{field(b, duplicate.source.name).tab().append(duplicate.line);}
		if(evaluated){b.tab().append(own).tab().append(other);}
		else{b.append("\tNA\tNA");}
		out.writer.print(b.tab().append(status).nl());
	}
	private void writeModels(final String path){
		try(Output out=new Output(path)){
			out.writer.println("accession\trelease\tmodel_sha80\tGA_bits\tmodel_length\tqueried\tbasis\tdeclared_positive_rows\toutside_input_paths\trepeated_exact_anchors\tunique_scoped_anchors\tdistinct_inputs\tno_own_hit\tmultiple_own_hits\tunsearched\tsingle_full\tsingle_not_full\tother_model_overlaps\tmargin_min\tmargin_p05\tmargin_median\tmargin_p95\tmargin_max\tcoverage_min\tcoverage_p05\tcoverage_median\tcoverage_p95\tcoverage_max\tsingle_weak_overlap\tfull_margin_min\tfull_margin_p05\tfull_margin_median\tfull_margin_p95\tfull_margin_max\tstatus");
			for(final ModelState item:models){
				final Model m=item.model;
				for(int basis=0; basis<BASES.length; basis++){
					final Stats s=item.stats[basis];
					require(s.rows==s.outside+s.duplicates+s.unique && s.unique==s.missing+s.multiple+s.unsearched+s.full+s.partial+s.weakOverlap,
						"Per-model anchor outcome counts must conserve declared rows: "+m.accession);
					b.clear().append(m.accession).tab().append(m.release).tab().append(m.sha80).tab().append(m.ga, 6).tab();
					b.append(m.length).tab().append(item.query>=0).tab().append(BASES[basis]).tab().append(s.rows).tab();
					b.append(s.outside).tab().append(s.duplicates).tab().append(s.unique).tab().append(s.genomes).tab();
					b.append(s.missing).tab().append(s.multiple).tab().append(s.unsearched).tab().append(s.full).tab().append(s.partial).tab().append(s.alternatives);
					quantiles(b, s.margins);
					quantiles(b, s.coverage);
					b.tab().append(s.weakOverlap);
					quantiles(b, s.fullMargins);
					out.writer.print(b.append("\tOBSERVED_HITS_ONLY_NOT_ADMISSION_THRESHOLDS\n"));
				}
			}
		}
	}
	/** Appends nearest-rank min/p05/median/p95/max; sorts the caller-owned live list. */
	static void quantiles(final ByteBuilder b, final DoubleList values){
		values.sort();
		for(final int percentile:PERCENTILES){
			b.tab();
			if(values.size==0){b.append("NA");}
			else{
				final int rank=(int)Math.max(1L, ((long)values.size*percentile+99)/100);
				b.append(values.array[rank-1], 6);
			}
		}
	}
	private static ByteBuilder begin(final ByteBuilder b, final Source source, final long line){
		field(b, source.name).tab();
		return field(b, source.file).tab().append(line).tab();
	}
	private static ByteBuilder field(final ByteBuilder b, final String value){return BuildNcrnaReference.field(b, value);}
	/** The existing bound table is already escaped; copy its UTF-8 bytes exactly. */
	private static ByteBuilder copy(final ByteBuilder b, final LineParser1 p, final int column){
		p.setBounds(column);
		return b.append(p.line(), p.a(), p.b()-p.a());
	}
	private static String text(final LineParser1 p, final int column){
		p.setBounds(column);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static boolean same(final LineParser1 p, final int column, final byte[] bytes){
		p.setBounds(column);
		if(p.b()-p.a()!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(bytes[i]!=p.line()[p.a()+i]){return false;}}
		return true;
	}
	private static boolean same(final byte[] bytes, final ByteBuilder b){
		if(bytes.length!=b.length()){return false;}
		for(int i=0; i<bytes.length; i++){if(bytes[i]!=b.array[i]){return false;}}
		return true;
	}
	private static char strand(final LineParser1 p, final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Prior strand must remain explicit");
		return p.termEquals("+", column) ? '+' : '-';
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	/*--------------------------------------------------------------*/
	/*----------------         Retained state       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Source{
		Source(final String name_, final String format_, final String file_, final String bindings_, final String authority_, final long expected_){
			name=name_;
			format=format_;
			file=file_;
			bindings=bindings_;
			authority=authority_;
			expected=expected_;
		}
		final String name, format, file, bindings, authority;
		final long expected;
		long rows, notPositive, outside, duplicates, unique;
	}
	private static final class InputState{
		InputState(final Genome genome_, final int index_){
			genome=genome_;
			index=index_;
			path=field(new ByteBuilder(), genome.fasta).toBytes();
		}
		final Genome genome;
		final int index;
		final byte[] path;
		final ObjectIntMap<Target> lengths=new ObjectIntMap<Target>(16, Target.class);
		long bases;
		final ArrayList<Hit> hits=new ArrayList<Hit>();
		final ArrayList<Anchor> anchors=new ArrayList<Anchor>();
	}
	private static final class Binding{
		Binding(final String fasta_, final InputState input_, final byte[] assembly_){
			fasta=fasta_;
			input=input_;
			assembly=assembly_;
		}
		final String fasta;
		final InputState input;
		final byte[] assembly;
	}
	private static final class Anchor{
		Anchor(final Source source_, final long line_, final ModelState model_, final int basis_, final int start_, final int stop_, final char strand_, final byte[] target_){
			source=source_;
			line=line_;
			model=model_;
			basis=basis_;
			start=start_;
			stop=stop_;
			strand=strand_;
			target=target_;
		}
		final Source source;
		final long line;
		final ModelState model;
		final int basis, start, stop;
		final char strand;
		final byte[] target;
	}
	private static final class ModelState{
		ModelState(final Model model_, final int query_, final int inputs){
			model=model_;
			query=query_;
			for(int i=0; i<stats.length; i++){stats[i]=new Stats(inputs);}
		}
		final Model model;
		final int query;
		final Stats[] stats=new Stats[BASES.length];
	}
	private static final class Stats{
		Stats(final int count){inputs=new boolean[count];}
		final boolean[] inputs;
		final DoubleList margins=new DoubleList(1), coverage=new DoubleList(1), fullMargins=new DoubleList(1);
		long rows, outside, duplicates, unique, genomes, missing, multiple, unsearched, full, partial, alternatives, weakOverlap;
	}
	private final NcrnaEvidenceCatalog catalog;
	private final Roster roster;
	private final NcrnaReferenceGeometry geometry;
	private final ObjectMap<String,InputState> byId=new ObjectMap<String,InputState>(String.class, InputState.class);
	private final ObjectMap<String,InputState> byPath=new ObjectMap<String,InputState>(String.class, InputState.class);
	private final ObjectMap<String,ModelState> byModel=new ObjectMap<String,ModelState>(String.class, ModelState.class);
	private final ArrayList<ModelState> models=new ArrayList<ModelState>();
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ByteBuilder b=new ByteBuilder(), key=new ByteBuilder();
	private Anchor duplicate;
	private final Target probe=new Target();
	private boolean evaluated;
	private int own, other;
	private String status;
	private static final int CORE=0, WHOLE=1, CONSTRUCTION=2;
	private static final String[] BASES={"PRIOR_ACCEPTED_CORE", "PRIOR_ACCEPTED_WHOLE_INTERVAL", "REVIEWED_CONSTRUCTION_CORE"};
	private static final int[] PERCENTILES={0, 5, 50, 95, 100};
	private static final String BOUND="bound_prior_v1", TPP="tpp_construction_v1";
	static final String[] SUFFIXES={".anchor_rows.tsv.gz", ".anchor_links.tsv.gz", ".anchor_models.tsv", ".anchor_summary.tsv"};
	public static final String MANIFEST_HEADER="source_id\tformat\tdecisions\tbindings\tauthority\texpected_rows";
	public static final String BOUND_HEADER=NcrnaReferencePrior.ROW_HEADER;
	public static final String TPP_HEADER="locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines\tcm_state\traw_hits\toverlapping_sense_hits\tcm_model_from\tcm_model_to\tcm_seq_from\tcm_seq_to\tcm_trunc\tcm_bits\toverlap_bp\tleft_excess\tright_excess\tgeometry\tselected_phase\tselected_state\tselected_context_start\tselected_context_stop\tselected_bits\tselected_genomic_start\tselected_genomic_stop\tconstruction_candidate\tdisposition\tcoverage_critical_note\tinitial_construction_candidate\tinitial_disposition";
}
