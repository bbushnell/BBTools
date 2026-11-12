package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectIntMap;
import map.LongIntMap;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaReferenceGeometry.Pair;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Constructs two explicitly hypothetical reference scenarios from saved evidence.
 * Adapted from A_SampleByteFile. Native catalog/model/geometry readers retain the
 * original population. Only a declared draft margin differs between the scenarios;
 * reviewed prior positives and all structural/conflict rules stay fixed. Native
 * events must all pass their captured model's GA: this tool consumes a GA-filtered
 * search, not an unfiltered sensitive scan. It rejects rather than drops violations.
 *
 * Model-implied family resolution is an assumption for these diagnostics, never
 * biological adjudication. No reference GFF, approval claim, caller or scan runs.
 * One thread owns all mutable state. Exact model intervals collapse duplicate
 * events; different models remain separate until identical same-family MASK rows
 * are collapsed. No connected-component interval is used as a gene boundary.
 * Original source links remain a required companion: witnesses identify their
 * exact input/native-table/native-line join key rather than copying that table.
 * @author Brian Bushnell, Raiden
 */
public final class DraftNcrnaScenarios{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, DraftNcrnaScenarios.class, false);
		try{new DraftNcrnaScenarios(pp.args, pp.outstream).process(timer);}
		finally{Shared.closeStream(pp.outstream);}
	}
	private DraftNcrnaScenarios(final String[] args, final PrintStream stream){
		outstream=stream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("census")){census=b;}
			else if(a.equals("memberships")){memberships=b;}
			else if(a.equals("scope")){scope=b;}
			else if(a.equals("assets")){assets=b;}
			else if(a.equals("queries")){queries=b;}
			else if(a.equals("policy")){policyPath=b;}
			else if(a.equals("priors")){priorPath=b;}
			else if(a.equals("links")){linksPath=b;}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		prefix=common.in1;
		require(prefix!=null && census!=null && memberships!=null && scope!=null && assets!=null && queries!=null
			&& policyPath!=null && priorPath!=null && linksPath!=null && common.out1!=null,
			"Supply in=evidence_prefix census= memberships= scope= assets= queries= policy= priors= links= out=new_directory");
		require(!common.overwrite && !common.append, "Hypothetical construction requires new outputs");
		require(Tools.testInputFiles(false, true, census, memberships, scope, assets, queries, policyPath, priorPath, linksPath), "Scenario metadata must be readable");
		output=new File(common.out1).getAbsoluteFile().toPath();
		require(!Files.exists(output), "Scenario output directory must not exist");
		policyEncoded=encoded(policyPath);
		scopeEncoded=encoded(scope);
		linksEncoded=encoded(linksPath);
		priorEncoded=encoded(priorPath);
		witnessEncoded=encoded(path("witnesses.tsv.gz"));
		priorWitnessEncoded=encoded(path("prior_witnesses.tsv.gz"));
	}

	private void process(final Timer timer){
		catalog=NcrnaEvidenceCatalog.load(census, memberships);
		roster=NcrnaReferenceEvidence.loadModels(scope, assets, queries);
		geometry=NcrnaReferenceGeometry.load(prefix, catalog, roster);
		hits=geometry.hits();
		checkSearchCutoffs();
		flags=new int[hits.length];
		canonical=new int[hits.length];
		leader=new int[hits.length];
		before=new byte[hits.length];
		after=new byte[hits.length];
		maskModels=new int[hits.length];
		counts=new long[roster.queried()][2][STATES.length];
		policy=new Policy[roster.queried()];
		familyNumber=new int[roster.queried()];
		final ObjectIntMap<String> familyOrder=new ObjectIntMap<String>(String.class);
		for(int i=0; i<roster.queried(); i++){
			final String family=roster.query(i).scope.family;
			family(family);
			if(!familyOrder.contains(family)){familyOrder.put(family, familyOrder.size());}
			familyNumber[i]=familyOrder.get(family);
		}
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			inputs.put(genome.id, new InputState(genome));
		}
		readPolicy();
		readPriors();
		checkPriorContigs();
		readLinks();
		collapseEvents();
		addCompetition();
		classify();
		collapseMasks();
		try{Files.createDirectory(output);}
		catch(IOException failure){throw new IllegalStateException("Could not create diagnostic output directory", failure);}
		writeSources();
		writeScenarios();
		writeWitnesses();
		writeCounts();
		writeSummary();
		timer.stop();
		outstream.println("Constructed two hypothetical scenarios from "+hits.length+" saved events; time: "+timer);
		outstream.println(QUALIFICATION);
	}
	/** Every possible competitor is GA-passing; prior positives are not tested here. */
	private void checkSearchCutoffs(){
		for(final Hit hit:hits){
			final Model model=roster.query(hit.query);
			final boolean pass=hit.bits>=model.ga;
			if(hit.ga!=pass){
				throw new IllegalArgumentException("Saved GA flag disagrees with captured threshold for "+model.accession+" input="+hit.input.genome.id+" native_line="+hit.line);
			}
			if(!pass){
				throw new IllegalArgumentException("Scenario construction requires GA-passing native events; below-GA "+model.accession+" input="+hit.input.genome.id+" native_line="+hit.line);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Declared policies      ----------------*/
	/*--------------------------------------------------------------*/

	private void readPolicy(){
		try(Input in=new Input(policyPath)){
			RrnaResourceIO.header(in.reader, POLICY_HEADER, policyPath);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				parser.set(row);
				require(parser.terms()==23, "Draft policy requires its complete twenty-three-field schema");
				int query=-1;
				for(int i=0; i<policy.length; i++){if(parser.termEquals(roster.query(i).accession, 0)){query=i;}}
				require(query>=0 && policy[query]==null, "Policy names an unqueried or repeated model");
				final Model model=roster.query(query);
				require(parser.termEquals(model.scope.family, 1) && parser.termEquals(model.release, 2)
					&& parser.termEquals(model.sha80, 3) && parser.parseDouble(4)==model.ga && parser.parseInt(5)==model.length,
					"Policy family/release/model/GA/length must match captured metadata");
				require(parser.termEquals("DRAFT_NOT_APPLIED", 20), "This tool accepts illustrative draft policies, not approval assertions");
				final int action=parser.termEquals("BOUNDARY_UNSELECTED", 18) ? MARGIN
					: parser.termEquals("NO_AUTOMATIC_ADDITION", 18) ? NO_ADD
					: parser.termEquals("CONTEXT_PENDING_BRIAN", 18) ? CONTEXT_ONLY
					: parser.termEquals("MASK_PENDING_SPARSE_ANCHOR_REVIEW", 18)
					|| parser.termEquals("MASK_PENDING_STRUCTURE_REVIEW", 18)
					|| parser.termEquals("MASK_PENDING_SUBTYPE_REVIEW", 18) ? FIXED_MASK : -1;
				require(action>=0, "Unknown draft action cannot acquire scoring authority");
				final double margin=action==MARGIN ? parser.parseDouble(17) : 0;
				require(action==MARGIN ? Double.isFinite(margin) && margin>=0 && parser.parseDouble(14)==margin : parser.termEquals("NA", 17),
					"Only margin-comparison rows may supply the identical declared full-margin p95");
				policy[query]=new Policy(action, margin, line);
				policyRows++;
				if(action==MARGIN){variableRows++;}
			}
		}
		require(policyRows>0, "An explicit nonempty draft model policy is required");
	}

	/** Retains every original classification; only exact accepted assignments score. */
	private void readPriors(){
		try(Input in=new Input(priorPath)){
			RrnaResourceIO.header(in.reader, NcrnaReferencePrior.ROW_HEADER, priorPath);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				parser.set(row);
				require(parser.terms()==21, "Bound priors require all original classification/interval/authority fields");
				final InputState input=inputs.get(decoded(6));
				require(input!=null && decoded(7).equals(input.genome.fasta), "Prior must bind an exact catalog input and FASTA");
				final String accession=decoded(10), assignment=decoded(11);
				final Model model=roster.model(accession);
				require(model!=null, "Prior source family requires captured model metadata");
				family(model.scope.family);
				final boolean core=assignment.equals(accession) && (accession.equals("RF00174") || accession.equals("RF01482") || accession.equals("RF01689"));
				final boolean whole=assignment.equals("ACCEPTED_NOVEL_FMN") && accession.equals("RF00050");
				require(core || whole || assignment.equals("UNRESOLVED") || assignment.equals("UNREVIEWED") || assignment.equals("B12clade"), "Unknown prior assignment must remain unresolved by an explicit supported rule");
				final int start=parser.parseInt(core ? 14 : 12), stop=parser.parseInt(core ? 15 : 13), length=parser.parseInt(9);
				final int sourceStart=parser.parseInt(12), sourceStop=parser.parseInt(13), coreStart=parser.parseInt(14), coreStop=parser.parseInt(15);
				require(sourceStart>0 && sourceStop>=sourceStart && sourceStop<=length
					&& (coreStart==0 && coreStop==0 || coreStart>0 && coreStop>=coreStart && coreStop<=length), "Original prior source/core coordinates are invalid");
				require(!whole || parser.termEquals("PRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS", 17), "FMN accepted prior must retain reviewed whole-interval semantics");
				require(!core || parser.termEquals("SOURCE_AND_REVIEWED_CORE", 17), "Cobalamin accepted prior must retain reviewed core semantics");
				require(!decoded(0).isEmpty() && !decoded(18).isEmpty() && !decoded(19).isEmpty(), "Prior source and authority fields must remain explicit");
				final char strand=strand(16);
				require(start>0 && stop>=start && stop<=length, "Prior interval must retain valid original core/whole/source bounds");
				probe.set(parser, 8);
				PriorTarget target=input.targets.get(probe);
				if(target==null){
					target=new PriorTarget(probe.copy(), length);
					input.targets.put(target.name, target);
				}
				require(target.length==length, "Prior rows disagree about observed literal contig length");
				final Prior prior=new Prior(input, target, model, start, stop, strand, core || whole, line, row);
				if(prior.accepted){
					for(final Prior existing:target.positive){
						require(existing.start!=start || existing.stop!=stop || existing.strand!=strand || !existing.model.scope.family.equals(model.scope.family), "Duplicate accepted physical prior needs explicit adjudication");
					}
					target.positive.add(prior);
					acceptedPriors++;
				}
				priors.add(prior);
			}
		}
	}

	/** Saved observations validate prior coordinates; no FASTA payload is reopened. */
	private void checkPriorContigs(){
		InputState input=null;
		try(Input in=new Input(prefix+".contigs.tsv.gz")){
			RrnaResourceIO.header(in.reader, WriteNcrnaReference.CONTIG_HEADER, prefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==4, "Saved contig observations require four fields");
				if(input==null || !parser.termEquals(input.genome.id, 0)){input=inputs.get(decoded(0));}
				require(input!=null && same(1, input.fasta), "Observed contigs must bind exact catalog inputs");
				probe.set(parser, 2);
				final PriorTarget target=input.targets.get(probe);
				if(target==null){continue;}
				require(!target.observed && parser.parseInt(3)==target.length, "Prior contig observation is duplicate or inconsistent");
				target.observed=true;
			}
		}
		for(final Prior prior:priors){require(prior.target.observed, "Every prior needs its original observed contig");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Sources and geometry    ----------------*/
	/*--------------------------------------------------------------*/

	private void readLinks(){
		NcrnaReferenceGeometry.InputState input=null;
		Source source=null;
		final LongHashSet pairs=new LongHashSet(16, .7f);
		int sourceRows=0;
		long count=0;
		try(Input in=new Input(linksPath)){
			RrnaResourceIO.header(in.reader, LINKS_HEADER, linksPath);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==10, "Source links require their complete ten-field schema");
				if(input==null || !parser.termEquals(input.genome.id, 0)){input=geometry.input(decoded(0));}
				require(input!=null, "Source link must retain the exact native input/table binding");
				final int index=input.find(parser, 3, parser.parseLong(4));
				require(index>=0, "Source link has no saved native event");
				probe.set(parser, 1);
				if(source==null || source.input!=input.genome || !source.path.equals(probe)){
					final InputState owner=inputs.get(input.genome.id);
					source=owner.sources.get(probe);
					if(source==null){
						boolean declared=false;
						for(int i=0; i<owner.genome.memberships(); i++){
							BuildNcrnaReference.field(buffer.clear(), owner.genome.membership(i).gff);
							if(same(1, buffer.toBytes())){declared=true;}
						}
						require(declared, "Source links must name a GFF declared for their exact input");
						source=new Source(owner.genome, probe.copy());
						owner.sources.put(source.path, source);
					}
				}
				final long sourceLine=parser.parseLong(2);
				require(sourceLine>0, "Source links must retain positive physical GFF lines");
				int ordinal=source.rows.get(sourceLine);
				if(ordinal<0){
					require(sourceRows<Integer.MAX_VALUE, "Source-line identity index exceeded its native capacity");
					ordinal=sourceRows++;
					source.rows.put(sourceLine, ordinal);
				}
				require(pairs.add(((long)ordinal<<32)|(index&0xffffffffL)), "Repeated input/GFF/source-line/native-line link");
				require(parser.termEquals("RNA_ROW", 6) || parser.termEquals("COMPONENT_ROW", 6)
					|| parser.termEquals("GENE_OR_TRANSCRIPT_CONTEXT", 6) || parser.termEquals("OTHER_ROW", 6), "Unknown source kind");
				require(parser.termEquals("SAME_STRAND_OVERLAP", 7) || parser.termEquals("OPPOSITE_STRAND_OVERLAP", 7)
					|| parser.termEquals("UNASSESSABLE", 7), "Unknown source association");
				final int claims=NcrnaReferenceStrata.readFlags(parser, 8), warnings=NcrnaReferenceStrata.readFlags(parser, 9);
				count++;
				if(!parser.termEquals("SAME_STRAND_OVERLAP", 7)){continue;}
				// A direct RFAM claim remains a claim on regulatory/tagged-gene rows.
				// Attribute warnings below remain specific to an RNA row itself.
				if((claims&(NcrnaReferenceStrata.OTHER_FAMILY_SOURCE|NcrnaReferenceStrata.UNKNOWN_RFAM))!=0){flags[index]|=CROSS_SOURCE;}
				if((claims&NcrnaReferenceStrata.SAME_FAMILY_SOURCE)!=0){flags[index]|=SOURCE_WARNING;}
				if(!parser.termEquals("RNA_ROW", 6)){continue;}
				if((warnings&RNA_WARNINGS)!=0){flags[index]|=SOURCE_WARNING;}
			}
		}
		require(count==geometry.profiles().historicalPairs(), "Source links must conserve the saved historical/model pair population");
		sourceLinks=count;
	}

	private void collapseEvents(){
		ordered=hits.clone();
		Arrays.sort(ordered, this::comparePosition);
		Hit previous=null;
		for(final Hit hit:ordered){
			final boolean duplicate=previous!=null && samePhysical(previous, hit) && previous.query==hit.query;
			final int first=duplicate ? canonical[previous.index] : hit.index;
			canonical[hit.index]=first;
			leader[hit.index]=first;
			assert(canonical[first]==first) : "A duplicate maps directly to its first model interval, never to a chain that can change after mask collapse";
			if(duplicate){
				final Hit original=hits[first];
				flags[first]|=flags[hit.index];
				if(original.bits!=hit.bits || original.from!=hit.from || original.to!=hit.to || original.truncation!=hit.truncation
					|| original.range!=hit.range || original.left!=hit.left || original.right!=hit.right || original.ga!=hit.ga){flags[first]|=CONFLICTING_DUPLICATE;}
			}else{modelIntervals++;}
			previous=hit;
		}
	}

	/** Every saved queried model participates, including those with no addition rule. */
	private void addCompetition(){
		for(final Pair pair:geometry.pairs()){
			final Hit a=pair.a, b=pair.b;
			if(!a.valid() || !b.valid() || a.strand!=b.strand){continue;}
			final int first=canonical[a.index], second=canonical[b.index];
			if(first==second){continue;}
			final int bit=familyNumber[a.query]==familyNumber[b.query] ? SAME_COMPETITOR : CROSS_COMPETITOR;
			flags[first]|=bit;
			flags[second]|=bit;
		}
	}

	private void classify(){
		InputState input=null;
		for(final Hit hit:ordered){
			if(canonical[hit.index]!=hit.index){continue;}
			final Model model=roster.query(hit.query);
			final Policy rule=policy[hit.query];
			if(input==null || input.genome!=hit.input.genome){input=inputs.get(hit.input.genome.id);}
			final PriorTarget target=input.targets.get(hit.target);
			if(target!=null && hit.valid()){
				for(final Prior prior:target.positive){
					if(prior.strand==hit.strand && prior.start<=hit.stop && prior.stop>=hit.start){flags[hit.index]|=PRESERVED_PRIOR;}
				}
			}
			byte state;
			final int blockers=flags[hit.index];
			if(rule==null || rule.action==CONTEXT_ONLY){state=CONTEXT;}
			else if(rule.action==NO_ADD || !hit.valid() || (blockers&(CROSS_SOURCE|CROSS_COMPETITOR|CONFLICTING_DUPLICATE|PRESERVED_PRIOR))!=0){state=UNRESOLVED;}
			else if(rule.action==FIXED_MASK || !hit.full() || hit.left || hit.right || (blockers&(SOURCE_WARNING|SAME_COMPETITOR))!=0){state=MASK;}
			else{state=SCORED;}
			before[hit.index]=state;
			// Addition avoids cancellation at a decimal boundary; no epsilon lowers it.
			after[hit.index]=state==SCORED && hit.bits<model.ga+rule.margin ? MASK : state;
			counts[hit.query][0][state]++;
			counts[hit.query][1][after[hit.index]]++;
		}
	}

	/** A repeated physical mask is one interval with multiple retained witnesses. */
	private void collapseMasks(){
		Hit firstMask=null, previous=null;
		for(final Hit hit:ordered){
			if(canonical[hit.index]!=hit.index){continue;}
			if(previous==null || !samePhysical(previous, hit)){firstMask=null;}
			if(before[hit.index]==MASK || after[hit.index]==MASK){
				if(firstMask==null){
					firstMask=hit;
					maskModels[hit.index]=1;
				}else{
					require(before[firstMask.index]==MASK && after[firstMask.index]==MASK && before[hit.index]==MASK && after[hit.index]==MASK,
						"Coincident same-family alternatives must be blocked by the saved positive-overlap competition rule in both scenarios");
					leader[hit.index]=firstMask.index;
					maskModels[firstMask.index]++;
				}
			}
			previous=hit;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Reports          ----------------*/
	/*--------------------------------------------------------------*/

	private void writeSources(){
		try(Output out=new Output(path("sources.tsv"))){
			out.writer.println("role\tpath\tqualification");
			final String[] roles={"evidence_prefix", "census", "memberships", "scope", "assets", "queries", "draft_policy", "bound_priors", "required_source_links"};
			final String[] paths={prefix, census, memberships, scope, assets, queries, policyPath, priorPath, linksPath};
			for(int i=0; i<roles.length; i++){
				buffer.clear().append(roles[i]).tab();
				field(paths[i]).tab().append(QUALIFICATION).nl();
				out.writer.print(buffer);
			}
		}
	}
	private void writeScenarios(){
		try(Output low=new Output(path("ga.tsv.gz")); Output high=new Output(path("margin.tsv.gz")); Output priorRows=new Output(path("prior_witnesses.tsv.gz"))){
			low.writer.println(CompareNcrnaScenarios.HEADER);
			high.writer.println(CompareNcrnaScenarios.HEADER);
			priorRows.writer.println("decision_id\tbound_line\t"+NcrnaReferencePrior.ROW_HEADER);
			for(final Prior prior:priors){
				buffer.clear().append("prior_").append(prior.line).tab().append(prior.line).tab().append(prior.original).nl();
				priorRows.writer.print(buffer);
				buffer.clear().append("prior_").append(prior.line).tab().append(prior.input.genome.id).tab()
					.append(prior.target.name.bytes).tab().append(prior.start).tab().append(prior.stop).tab().append(prior.strand).tab()
					.append(prior.model.scope.family).tab().append(prior.accepted ? "RESOLVED\tscored" : "UNRESOLVED\tunresolved")
					.append("\tpreserved_prior\t").append(prior.model.accession).append("\toriginal_classification_preserved\t")
					.append(priorEncoded).append(':').append(prior.line).tab().append(priorWitnessEncoded).append(':').append(prior.line).nl();
				low.writer.print(buffer);
				high.writer.print(buffer);
				physicalCounts[0][prior.accepted ? SCORED : UNRESOLVED]++;
				physicalCounts[1][prior.accepted ? SCORED : UNRESOLVED]++;
				family(prior.model.scope.family).counts[0][prior.accepted ? SCORED : UNRESOLVED]++;
				family(prior.model.scope.family).counts[1][prior.accepted ? SCORED : UNRESOLVED]++;
			}
			for(final Hit hit:ordered){
				if(canonical[hit.index]!=hit.index || leader[hit.index]!=hit.index){continue;}
				writeCandidate(low, hit, 0);
				writeCandidate(high, hit, 1);
			}
		}
	}
	private void writeCandidate(final Output out, final Hit hit, final int side){
		final Model model=roster.query(hit.query);
		final Policy rule=policy[hit.query];
		final byte state=side==0 ? before[hit.index] : after[hit.index];
		assert(maskModels[hit.index]<=1 || state==MASK) : "collapseMasks merges only fixed masks; it cannot merge scored loci or change their boundaries";
		// Catalog IDs, family names, accessions and these status tokens are validated
		// ASCII; target/path bytes already carry the evidence TSV escaping. Avoid
		// per-event Path resolution, String concatenation and UTF-8 array allocation.
		buffer.clear().append("event_").append(hit.line).tab().append(hit.input.genome.id).tab().append(hit.target.bytes)
			.tab().append(hit.start).tab().append(hit.stop).tab().append(hit.strand).tab().append(model.scope.family)
			.tab().append(state==SCORED || state==MASK ? "RESOLVED" : "UNRESOLVED").tab().append(STATES[state])
			.append("\tdiagnostic_declared_policy\t").append(maskModels[hit.index]>1 ? "." : model.accession)
			.tab().append(reason(hit, side)).tab().append(rule==null ? scopeEncoded : policyEncoded);
		if(rule!=null){buffer.append(':').append(rule.line);}
		buffer.tab().append(witnessEncoded).append('#').append(hit.input.genome.id).append("/event_").append(hit.line).nl();
		out.writer.print(buffer);
		physicalCounts[side][state]++;
		family(model.scope.family).counts[side][state]++;
	}
	private void writeWitnesses(){
		try(Output out=new Output(path("witnesses.tsv.gz"))){
			out.writer.println("decision_id\tinput_id\tnative_table\tnative_line\taccession\tmodel_interval_first_line\tflags\tmodel_before\tmodel_after\tbefore_reason\tafter_reason\trequired_source_links\tqualification");
			for(final Hit hit:hits){
				final int representative=canonical[hit.index], outputLeader=leader[representative];
				assert(canonical[representative]==representative && leader[outputLeader]==outputLeader)
					: "Every original event must reach one emitted decision through at most one duplicate and one exact-mask link";
				buffer.clear().append("event_").append(hits[outputLeader].line).tab().append(hit.input.genome.id).tab()
					.append(hit.table).tab().append(hit.line).tab().append(roster.query(hit.query).accession).tab()
					.append(hits[representative].line).tab().append(flags[representative]).tab()
					.append(STATES[before[representative]]).tab().append(STATES[after[representative]]).tab()
					.append(reason(hits[representative], 0)).tab().append(reason(hits[representative], 1)).tab();
				buffer.append(linksEncoded).tab().append(QUALIFICATION).nl();
				out.writer.print(buffer);
			}
		}
	}
	private void writeCounts(){
		try(Output out=new Output(path("counts.tsv"))){
			out.writer.println("accession\tfamily\tscenario\tunique_model_intervals\tscored\tdont_care\tunresolved\tcontext\tqualification");
			for(int i=0; i<counts.length; i++){
				for(int side=0; side<2; side++){
					long total=0;
					for(final long count:counts[i][side]){total+=count;}
					buffer.clear().append(roster.query(i).accession).tab().append(roster.query(i).scope.family).tab().append(side==0 ? "GA" : "GA_PLUS_MARGIN").tab().append(total);
					for(final long count:counts[i][side]){buffer.tab().append(count);}
					out.writer.print(buffer.tab().append(QUALIFICATION).nl());
				}
			}
		}
		try(Output out=new Output(path("family_counts.tsv"))){
			out.writer.println("family\tscenario\tscored_rows\tmask_rows\tunresolved_rows\tcontext_rows\tqualification");
			for(final FamilyCounts family:families){
				for(int side=0; side<2; side++){
					buffer.clear().append(family.name).tab().append(side==0 ? "GA" : "GA_PLUS_MARGIN");
					for(final long count:family.counts[side]){buffer.tab().append(count);}
					out.writer.print(buffer.tab().append(QUALIFICATION).nl());
				}
			}
		}
	}
	private void writeSummary(){
		long low=0, high=0;
		for(final long[][] model:counts){
			for(final long count:model[0]){low+=count;}
			for(final long count:model[1]){high+=count;}
		}
		require(low==modelIntervals && high==modelIntervals, "Each unique model interval must receive exactly one disposition on each side");
		for(int side=0; side<2; side++){
			for(int state=0; state<STATES.length; state++){
				long total=0;
				for(final FamilyCounts family:families){total+=family.counts[side][state];}
				require(total==physicalCounts[side][state], "Per-family emitted rows must conserve every disposition in the scenario summary");
			}
		}
		try(Output out=new Output(path("summary.tsv"))){
			out.writer.println("catalog_inputs\tinputs_observed\traw_events\tunique_model_intervals\tduplicate_events\tsource_links\tprior_rows\taccepted_priors\tpolicy_rows\tvariable_margin_rows\tscenario\tscored_rows\tmask_rows\tunresolved_rows\tcontext_rows\tevidence_status\tqualification");
			for(int side=0; side<2; side++){
				buffer.clear().append(catalog.size()).tab().append(geometry.profiles().observed()).tab().append(hits.length).tab().append(modelIntervals).tab().append(hits.length-modelIntervals)
					.tab().append(sourceLinks).tab().append(priors.size()).tab().append(acceptedPriors).tab().append(policyRows).tab().append(variableRows)
					.tab().append(side==0 ? "GA" : "GA_PLUS_MARGIN");
				for(final long count:physicalCounts[side]){buffer.tab().append(count);}
				out.writer.print(buffer.tab().append(geometry.profiles().evidenceState()).tab().append(QUALIFICATION).nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	private String reason(final Hit hit, final int side){
		final Policy rule=policy[hit.query];
		final int bits=flags[hit.index];
		if(rule==null){return "outside_declared_addition_scope";}
		if(rule.action==CONTEXT_ONLY){return "model_context_not_mature_tRNA_truth";}
		if(rule.action==NO_ADD){return "no_automatic_addition";}
		if(!hit.valid()){return "invalid_observed_interval";}
		if((bits&PRESERVED_PRIOR)!=0){return "overlaps_preserved_prior_no_widening_or_exemption";}
		if((bits&CONFLICTING_DUPLICATE)!=0){return "conflicting_duplicate_observations";}
		if((bits&(CROSS_SOURCE|CROSS_COMPETITOR))!=0){return "unresolved_cross_family_or_unknown_source_claim";}
		if(rule.action==FIXED_MASK){return "fixed_pending_structure_or_subtype_review";}
		if((bits&SAME_COMPETITOR)!=0){return "same_family_competing_interval";}
		if((bits&SOURCE_WARNING)!=0){return "direct_family_claim_or_RNA_attribute_warning";}
		if(!hit.full()){return "partial_model_support";}
		if(hit.left || hit.right){return "physical_contig_edge";}
		if(side==1 && hit.bits<roster.query(hit.query).ga+rule.margin){return "marginal_model_support";}
		return "hypothetical_model_family_assignment";
	}
	private boolean samePhysical(final Hit a, final Hit b){return a.input==b.input && a.target==b.target && a.strand==b.strand && a.start==b.start && a.stop==b.stop && familyNumber[a.query]==familyNumber[b.query];}
	private String decoded(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private boolean same(final int column, final byte[] bytes){
		parser.setBounds(column);
		if(bytes==null || parser.b()-parser.a()!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(parser.line()[parser.a()+i]!=bytes[i]){return false;}}
		return true;
	}
	private char strand(final int field){
		require(parser.termEquals("+", field) || parser.termEquals("-", field), "Prior strand must be explicit");
		return parser.termEquals("+", field) ? '+' : '-';
	}
	private ByteBuilder field(final String text){return BuildNcrnaReference.field(buffer, text);}
	/** Frames each constant path once; hot row writers borrow the immutable bytes. */
	private static byte[] encoded(final String text){return BuildNcrnaReference.field(new ByteBuilder(), text).toBytes();}
	/** Retains one primitive counter block per family, including zero-hit families. */
	private FamilyCounts family(final String name){
		FamilyCounts result=byFamily.get(name);
		if(result==null){
			result=new FamilyCounts(name);
			byFamily.put(name, result);
			families.add(result);
		}
		return result;
	}
	private String path(final String file){return output.resolve(file).toString();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/** Sorting starts after process has populated the per-model family ordinals. */
	private int comparePosition(final Hit a, final Hit b){
		int c=Integer.compare(a.input.index, b.input.index);
		if(c==0){c=Tools.compare(a.target.bytes, b.target.bytes);}
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		if(c==0){c=Integer.compare(familyNumber[a.query], familyNumber[b.query]);}
		if(c==0){c=Integer.compare(a.query, b.query);}
		if(c==0){c=Tools.compare(a.table, b.table);}
		return c==0 ? Long.compare(a.line, b.line) : c;
	}
	private static final class Policy{
		Policy(final int action_, final double margin_, final long line_){
			action=action_;
			margin=margin_;
			line=line_;
		}
		final int action;
		final double margin;
		final long line;
	}
	private static final class FamilyCounts{
		FamilyCounts(final String name_){name=name_;}
		final String name;
		final long[][] counts=new long[2][STATES.length];
	}
	private static final class InputState{
		InputState(final Genome genome_){
			genome=genome_;
			fasta=BuildNcrnaReference.field(new ByteBuilder(), genome.fasta).toBytes();
		}
		final Genome genome;
		final byte[] fasta;
		final ObjectMap<Target,PriorTarget> targets=new ObjectMap<Target,PriorTarget>(Target.class, PriorTarget.class);
		final ObjectMap<Target,Source> sources=new ObjectMap<Target,Source>(Target.class, Source.class);
	}
	private static final class Source{
		Source(final Genome input_, final Target path_){
			input=input_;
			path=path_;
		}
		final Genome input;
		final Target path;
		final LongIntMap rows=new LongIntMap(16);
	}
	private static final class PriorTarget{
		PriorTarget(final Target name_, final int length_){
			name=name_;
			length=length_;
		}
		final Target name;
		final int length;
		boolean observed;
		final ArrayList<Prior> positive=new ArrayList<Prior>();
	}
	private static final class Prior{
		Prior(final InputState input_, final PriorTarget target_, final Model model_, final int start_, final int stop_, final char strand_, final boolean accepted_, final long line_, final byte[] original_){
			input=input_;
			target=target_;
			model=model_;
			start=start_;
			stop=stop_;
			strand=strand_;
			accepted=accepted_;
			line=line_;
			original=original_;
		}
		final InputState input;
		final PriorTarget target;
		final Model model;
		final int start, stop;
		final char strand;
		final boolean accepted;
		final long line;
		final byte[] original;
	}
	private String prefix, census, memberships, scope, assets, queries, policyPath, priorPath, linksPath;
	private NcrnaEvidenceCatalog catalog;
	private Roster roster;
	private NcrnaReferenceGeometry geometry;
	private Hit[] hits, ordered;
	private Policy[] policy;
	private int[] flags, canonical, leader, maskModels, familyNumber;
	private byte[] before, after;
	private long[][][] counts;
	private final long[][] physicalCounts=new long[2][STATES.length];
	private long modelIntervals, sourceLinks, acceptedPriors, policyRows, variableRows;
	private final Path output;
	private final byte[] policyEncoded, scopeEncoded, linksEncoded, priorEncoded, witnessEncoded, priorWitnessEncoded;
	private final PrintStream outstream;
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final Target probe=new Target();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectMap<String,InputState> inputs=new ObjectMap<String,InputState>(String.class, InputState.class);
	private final ArrayList<Prior> priors=new ArrayList<Prior>();
	private final ObjectMap<String,FamilyCounts> byFamily=new ObjectMap<String,FamilyCounts>(String.class, FamilyCounts.class);
	private final ArrayList<FamilyCounts> families=new ArrayList<FamilyCounts>();
	private static final byte SCORED=0, MASK=1, UNRESOLVED=2, CONTEXT=3;
	private static final int MARGIN=0, NO_ADD=1, CONTEXT_ONLY=2, FIXED_MASK=3;
	private static final int CROSS_SOURCE=1, SOURCE_WARNING=2, CONFLICTING_DUPLICATE=4, SAME_COMPETITOR=8, CROSS_COMPETITOR=16, PRESERVED_PRIOR=32;
	private static final int RNA_WARNINGS=NcrnaReferenceStrata.PSEUDO_TRUE|NcrnaReferenceStrata.PSEUDOGENE_TAG|NcrnaReferenceStrata.PSEUDO_BIOTYPE
		|NcrnaReferenceStrata.PARTIAL_TRUE|NcrnaReferenceStrata.EXCEPTION|NcrnaReferenceStrata.UNKNOWN_BOOLEAN|NcrnaReferenceStrata.INCOMPATIBLE_ID;
	private static final String[] STATES={"scored", "dont_care", "unresolved", "context"};
	public static final String QUALIFICATION="HYPOTHETICAL_MODEL_FAMILY_ASSIGNMENT_NOT_BIOLOGICAL_REFERENCE";
	public static final String POLICY_HEADER="accession\tfamily\trelease\tmodel_sha80\tGA_bits\tmodel_length\tprior_anchor_N\tprior_score_N\tdistribution_basis\tanchor_N\tscore_N\tfull_score_N\tfull_margin_p05\tfull_margin_median\tfull_margin_p95\tcoverage_p05\tcoverage_median\tillustrative_p95_margin_bits\tnew_candidate_disposition\tqualification\tstatus\tanchor_tier\tdistinct_input_N";
	public static final String LINKS_HEADER="input_id\tgff\tline\tnative_table\tnative_line\tsource_type\tsource_kind\tassociation\trfam_comparison\tsource_flags";
}
