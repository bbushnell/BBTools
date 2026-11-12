package prok;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaReferenceGeometry.Pair;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaReferenceGff.Scores;
import prok.NcrnaReferenceStrata.SourceRow;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;

/** Applies a declared reference policy to captured evidence, without a new search.
 * Native original-GA inclusion supplies only the confident score component.
 * Full model span, no native truncation and no qualifying competing observation
 * or contradictory source warning are independently required. Marginal masking
 * requires the selected printed-score floor and no touched rounding boundary.
 * Below-floor and boundary-uncertain events remain unmasked provenance.
 *
 * All captured events survive with their own measurements and native witnesses.
 * Same-family witnesses collapse only with all-member reciprocal compatibility,
 * retaining every measurement and one original representative interval. Neither
 * weak witnesses nor broad tandem bridges choose representatives or join groups.
 * Reviewed prior intervals and assignments remain unchanged; overlapping new
 * evidence cannot widen them or create a wrong-family masking exemption.
 *
 * One thread owns the shared geometry and primitive decision arrays. GFF rows
 * live for one input at a time. Decimal objects are model-level setup only;
 * event thresholds and percentage output use integer arithmetic.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceTiers{

	static void process(final String input, final String attributes, final String attributeEvidence,
		final String priors, final String output, final String percent, final String collapsePercent, final String census,
		final String memberships, final String scope, final String assets, final String queries,
		final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceTiers tiers=new NcrnaReferenceTiers(input, output, percent, collapsePercent, catalog, roster);
		tiers.readThresholds(assets);
		tiers.readCaptureBounds(input+".batch_sources.tsv");
		tiers.strata=NcrnaReferenceStrata.load(input, attributes, attributeEvidence, tiers.geometry, catalog, roster);
		tiers.sourceFlags=tiers.strata.flags();
		tiers.readPriors(priors);
		tiers.checkPriorContigs(input+".contigs.tsv.gz");
		tiers.classify();
		tiers.write(input, attributes, attributeEvidence, priors, census, memberships, scope, assets, queries);
	}

	private NcrnaReferenceTiers(final String input, final String output, final String percent, final String collapsePercent,
		final NcrnaEvidenceCatalog catalog_, final Roster roster_){
		floor=new BigDecimal(percent);
		floorText=floor.toPlainString();
		final BigDecimal collapse=new BigDecimal(collapsePercent);
		collapseText=collapse.toPlainString();
		collapseMillionths=collapse.movePointRight(4).longValueExact();
		require(collapseMillionths>0 && collapseMillionths<=1000000,
			"collapsepercent must be within (0,100], with at most four decimal places");
		require(floor.compareTo(BigDecimal.valueOf(75))>=0 && floor.compareTo(BigDecimal.valueOf(100))<=0,
			"marginalpercent must lie within the captured 75..100 percent range");
		catalog=catalog_;
		roster=roster_;
		root=new File(output).getAbsoluteFile().toPath();
		require(rawSafe(root.toString()), "Reference output paths must be representable in the literal manifest");
		witnessPath=root.resolve("annotation_witnesses.tsv.gz").toString();
		require(!Files.exists(root), "Tier output directory must be fresh");
		geometry=NcrnaReferenceGeometry.load(input, catalog, roster, true);
		require(geometry.profiles().batchEvidence() && geometry.profiles().observed()==catalog.size(),
			"Tier construction requires complete original-GA-bound captured batch evidence");
		hits=geometry.hits();
		ordered=hits.clone();
		bits=new long[hits.length];
		state=new byte[hits.length];
		flags=new int[hits.length];
		groupOf=new int[hits.length];
		Arrays.fill(groupOf, -1);
		qualifies=new boolean[hits.length];
		uncertain=new boolean[hits.length];
		thresholds=new Threshold[roster.queried()];
		eventCounts=new long[roster.queried()+1][STATES.length];
		inputCounts=new long[catalog.size()][3];
		modelFamilies=new Family[roster.queried()];
		for(int i=0; i<modelFamilies.length; i++){
			final String name=roster.query(i).scope.family;
			Family family=familyByName.get(name);
			if(family==null){
				family=new Family(name);
				familyByName.put(name, family);
				families.add(family);
			}
			modelFamilies[i]=family;
		}
		Arrays.sort(ordered, (a, b)->compare(a, b));
	}

	/*--------------------------------------------------------------*/
	/*----------------      Thresholds and Priors    ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains literal GA decimals from captured model metadata, not a new CM parse. */
	private void readThresholds(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaReferenceEvidence.ASSET_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==6, "Captured model metadata requires six fields");
				final Model model=roster.model(text(0));
				require(model!=null, "Captured threshold uses an unknown model");
				final int index=model.queryIndex();
				if(index<0){continue;}
				require(thresholds[index]==null && parser.parseDouble(3)==model.ga, "Repeated or inconsistent captured GA");
				thresholds[index]=new Threshold(text(3), floor);
			}
		}
		for(final Threshold threshold:thresholds){require(threshold!=null, "Every queried model requires its captured GA");}
	}

	/** Reads only accepted per-batch completion metadata, never native alignment
	 * reports or sequence payloads. The promised75% re-threshold range must be
	 * supported by every original reporting cutoff, including zero-hit models. */
	private void readCaptureBounds(final String path){
		final ObjectSet<String> batches=new ObjectSet<String>(String.class), directories=new ObjectSet<String>(String.class);
		final boolean[] seen=new boolean[roster.queried()];
		final long[] counts=new long[roster.queried()];
		final LineParser1 p=new LineParser1('\t');
		long bases=0, records=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BATCH_SOURCES_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==11, "Captured batch source metadata requires eleven fields");
				final String batch=text(0), directory=text(2);
				require(batches.add(batch) && directories.add(directory) && new File(directory).isAbsolute(), "Captured batch IDs and native directories must be unique");
				Arrays.fill(seen, false);
				long batchBases=-1, batchRecords=-1;
				final byte[] encoded=BuildNcrnaReference.field(new ByteBuilder(), directory).toBytes();
				final String completionPath=new File(directory, "completion.tsv").getPath();
				try(Input completion=new Input(completionPath)){
					RrnaResourceIO.header(completion.reader, NcrnaReferenceBatch.COMPLETION_HEADER, completionPath);
					for(byte[] entry=completion.reader.nextLine(); entry!=null; entry=completion.reader.nextLine()){
						p.set(entry);
						require(p.terms()==9 && p.termEquals(batch, 0) && p.termEquals("COMPLETE", 8), "Incomplete or foreign captured model completion");
						int query=-1;
						for(int i=0; i<roster.queried(); i++){if(p.termEquals(roster.query(i).accession, 1)){query=i;}}
						require(query>=0 && !seen[query], "Unknown or repeated captured model threshold");
						seen[query]=true;
						final long n=p.parseLong(2), size=p.parseLong(3), hits=p.parseLong(5);
						if(batchBases<0){
							batchBases=size;
							batchRecords=n;
						}
						require(n==batchRecords && size==batchBases && n>0 && size>0 && p.parseLong(4)==Math.multiplyExact(size, 2) && hits>=0,
							"Captured model completions must conserve one batch's two-strand population");
						final BigDecimal ga=new BigDecimal(thresholds[query].gaText);
						require(new BigDecimal(p.parseString(6)).compareTo(ga)==0
							&& new BigDecimal(p.parseString(7)).multiply(BigDecimal.valueOf(100)).compareTo(ga.multiply(BigDecimal.valueOf(75)))==0,
							"Tier provenance requires the captured original GA and75%-of-GA reporting floor");
						counts[query]=Math.addExact(counts[query], hits);
						captureRows.add(buffer.clear().append(encoded).tab().append(entry).toBytes());
					}
				}
				for(final boolean present:seen){require(present, "Missing captured threshold cannot imply zero-hit coverage");}
				bases=Math.addExact(bases, batchBases);
				records=Math.addExact(records, batchRecords);
			}
		}
		long expectedBases=0, expectedRecords=0;
		for(int i=0; i<catalog.size(); i++){
			expectedBases=Math.addExact(expectedBases, catalog.get(i).bases);
			expectedRecords=Math.addExact(expectedRecords, catalog.get(i).records);
		}
		require(batches.size()>0 && bases==expectedBases && records==expectedRecords, "Captured reporting thresholds must cover the complete catalog");
		for(int i=0; i<counts.length; i++){
			if(counts[i]!=geometry.profiles().modelEvents(i)){
				throw new IllegalArgumentException("Captured threshold count differs for "+roster.query(i).accession
					+": saved="+geometry.profiles().modelEvents(i)+", captured="+counts[i]);
			}
		}
	}

	/** Only the existing explicit accepted-core/whole classifications are protected. */
	private void readPriors(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaReferencePrior.ROW_HEADER, path);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				parser.set(row);
				require(parser.terms()==21, "Bound prior rows require all original identity and authority fields");
				final String id=text(6);
				final NcrnaReferenceGeometry.InputState input=geometry.input(id);
				require(input!=null && text(7).equals(input.genome.fasta), "A prior must bind its original catalog input and FASTA");
				PriorInput owner=priorInputs.get(id);
				if(owner==null){
					owner=new PriorInput(input);
					priorInputs.put(id, owner);
				}
				final String accession=text(10), assignment=text(11);
				final Model model=roster.model(accession);
				require(model!=null && model.queryIndex()>=0, "Prior family requires a captured queried model");
				final boolean core=assignment.equals(accession) && (accession.equals("RF00174") || accession.equals("RF01482") || accession.equals("RF01689"));
				final boolean whole=assignment.equals("ACCEPTED_NOVEL_FMN") && accession.equals("RF00050");
				require(core || whole || assignment.equals("UNRESOLVED") || assignment.equals("UNREVIEWED") || assignment.equals("B12clade"),
					"Unknown prior classification cannot acquire protection automatically");
				final int length=parser.parseInt(9), start=parser.parseInt(core ? 14 : 12), stop=parser.parseInt(core ? 15 : 13);
				final int sourceStart=parser.parseInt(12), sourceStop=parser.parseInt(13), coreStart=parser.parseInt(14), coreStop=parser.parseInt(15);
				require(sourceStart>0 && sourceStop>=sourceStart && sourceStop<=length
					&& (coreStart==0 && coreStop==0 || coreStart>0 && coreStop>=coreStart && coreStop<=length), "Invalid retained prior source/core intervals");
				require(start>0 && stop>=start && stop<=length, "Prior bounds must retain their accepted core or whole-interval meaning");
				require(!core || parser.termEquals("SOURCE_AND_REVIEWED_CORE", 17), "Accepted core must retain original interval semantics");
				require(!whole || parser.termEquals("PRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS", 17), "FMN prior cannot silently acquire new aptamer bounds");
				require(!text(0).isEmpty() && !text(18).isEmpty() && !text(19).isEmpty(), "Prior source and authority remain required");
				final char strand=parser.termEquals("+", 16) ? '+' : parser.termEquals("-", 16) ? '-' : '?';
				require(strand!='?', "Prior strand must be explicit");
				probe.set(parser, 8);
				PriorTarget target=owner.targets.get(probe);
				if(target==null){
					target=new PriorTarget(probe.copy(), length);
					owner.targets.put(target.key, target);
				}
				require(target.length==length, "Prior contig lengths disagree");
				final Prior prior=new Prior(owner, target, model, start, stop, strand, core || whole, line, row);
				if(prior.accepted){
					for(final Prior existing:target.accepted){
						require(existing.start!=start || existing.stop!=stop || existing.strand!=strand || !existing.model.scope.family.equals(model.scope.family),
							"Repeated protected physical prior needs explicit resolution, not double counting");
					}
					target.accepted.add(prior);
					owner.accepted.add(prior);
				}
				priors.add(prior);
			}
		}
	}

	private void checkPriorContigs(final String path){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, WriteNcrnaReference.CONTIG_HEADER, path);
			PriorInput owner=null;
			String previous=null;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==4, "Saved contig observations require four fields");
				if(previous==null || !parser.termEquals(previous, 0)){
					previous=text(0);
					owner=priorInputs.get(previous);
				}
				if(owner==null){continue;}
				require(same(1, owner.fasta), "Prior contig must retain its original FASTA binding");
				probe.set(parser, 2);
				final PriorTarget target=owner.targets.get(probe);
				if(target==null){continue;}
				require(!target.observed && parser.parseInt(3)==target.length, "Prior contig is duplicate or has different bounds");
				target.observed=true;
			}
		}
		for(final Prior prior:priors){require(prior.target.observed, "Every prior must retain an observed literal contig");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Tier Assignment      ----------------*/
	/*--------------------------------------------------------------*/

	private void classify(){
		for(final Hit hit:hits){
			final long score=Math.round(hit.bits*10);
			if(!(score>-(1L<<50) && score<(1L<<50) && hit.bits==score/10.0)){
				throw new IllegalArgumentException("Saved score lacks native one-decimal precision: input="+hit.input.genome.id
					+", model="+roster.query(hit.query).accession+", native_line="+hit.line+", score="+hit.bits);
			}
			bits[hit.index]=score;
			final Threshold threshold=thresholds[hit.query];
			final long center=score*10;
			uncertain[hit.index]=center-5<=threshold.floorLow && center+5>=threshold.floorHigh;
			qualifies[hit.index]=center>=threshold.floorHigh && !uncertain[hit.index];
			final PriorInput owner=priorInputs.get(hit.input.genome.id);
			final PriorTarget target=owner==null ? null : owner.targets.get(hit.target);
			if(target==null || !hit.valid()){continue;}
			for(final Prior prior:target.accepted){
				if(prior.strand!=hit.strand || prior.start>hit.stop || prior.stop<hit.start){continue;}
				if(prior.model.scope.family.equals(roster.query(hit.query).scope.family)){
					flags[hit.index]|=PRIOR_SAME;
					prior.hits.add(hit.index);
				}else{flags[hit.index]|=PRIOR_OTHER;}
			}
		}
		for(final Pair pair:geometry.pairs()){
			final Hit a=pair.a, b=pair.b;
			if(!a.valid() || !b.valid() || a.strand!=b.strand || !relevant(a) || !relevant(b)){continue;}
			if(!sameFamily(a, b)){
				flags[a.index]|=COMPETITOR;
				flags[b.index]|=COMPETITOR;
			}
		}
		for(final Hit hit:hits){
			final int index=hit.index;
			final byte tier;
			if(!hit.valid()){tier=INVALID;}
			else if((flags[index]&PRIOR_OTHER)!=0){tier=REJECTED;}
			else if((flags[index]&PRIOR_SAME)!=0){tier=COVERED;}
			else if(hit.included && hit.full() && (sourceFlags[index]&SOURCE_BLOCKERS)==0
				&& (flags[index]&COMPETITOR)==0){tier=CONFIDENT;}
			else if(uncertain[index]){tier=UNCERTAIN;}
			else if(qualifies[index]){tier=MARGINAL;}
			else{tier=BELOW;}
			state[index]=tier;
			eventCounts[hit.query][tier]++;
			eventCounts[roster.queried()][tier]++;
			modelFamilies[hit.query].events[tier]++;
		}
		for(final Hit hit:hits){if(tandemBridge(hit)){flags[hit.index]|=BRIDGE;}}
	}
	/** A bound GA pass proves threshold passage for every permitted marginal floor. */
	private boolean relevant(final Hit hit){return (flags[hit.index]&PRIOR_OTHER)==0 && (hit.included || qualifies[hit.index]);}
	private boolean primary(final Hit hit){return state[hit.index]==CONFIDENT || state[hit.index]==MARGINAL;}
	private boolean sameFamily(final Hit a, final Hit b){return modelFamilies[a.query]==modelFamilies[b.query];}
	private long overlap(final Hit a, final Hit b){return Math.max(0L, (long)Math.min(a.stop, b.stop)-Math.max(a.start, b.start)+1);}
	private boolean compatible(final Hit a, final Hit b){
		final long shared=overlap(a, b)*1000000;
		return shared>=a.length()*collapseMillionths && shared>=b.length()*collapseMillionths;
	}
	/** A broad hit is witness-only when it substantially covers two disjoint,
	 * shorter native-included/full-model supports of the same family. Reverse
	 * reciprocal failure distinguishes a bridge from a legitimate group member.
	 * Extremal bounds detect a disjoint pair without a quadratic neighbor walk. */
	private boolean tandemBridge(final Hit hit){
		if(!primary(hit) || hit.neighbors==null){return false;}
		int minimumStop=Integer.MAX_VALUE, maximumStart=-1;
		for(int i=0; i<hit.neighbors.size; i++){
			final Hit other=hits[hit.neighbors.get(i)];
			if(!primary(other) || !other.included || !other.full() || !sameFamily(hit, other)
				|| hit.strand!=other.strand || other.length()>=hit.length() || compatible(hit, other)){continue;}
			if(overlap(hit, other)*1000000<other.length()*collapseMillionths){continue;}
			minimumStop=Math.min(minimumStop, other.stop);
			maximumStart=Math.max(maximumStart, other.start);
		}
		return maximumStart>minimumStop;
	}
	/** Qualifying observations alone form groups. Every added member must agree
	 * reciprocally with every existing member. Stable ranking makes input order
	 * irrelevant; the original representative is never widened to a union. */
	private ArrayList<Group> groups(final int from, final int to){
		assert(from>=0 && from<=to && to<=ordered.length) : "Grouping receives one catalog input's contiguous saved-event slice";
		final ArrayList<Hit> candidates=new ArrayList<Hit>();
		for(int i=from; i<to; i++){
			final Hit hit=ordered[i];
			if(primary(hit) && (flags[hit.index]&BRIDGE)==0){candidates.add(hit);}
		}
		candidates.sort((a, b)->rank(a, b));
		final ArrayList<Group> groups=new ArrayList<Group>();
		final IntList nearby=new IntList();
		for(final Hit hit:candidates){
			neighborGroups(hit, nearby);
			Group chosen=null;
			int selected=-1;
			for(int i=0; i<nearby.size; i++){
				final int index=nearby.get(i);
				final Group group=groups.get(index);
				boolean agree=true;
				for(int j=0; agree && j<group.members.size; j++){agree=compatible(hit, hits[group.members.get(j)]);}
				if(agree){
					chosen=group;
					selected=index;
					break;
				}
			}
			if(chosen==null){
				chosen=new Group(hit);
				selected=groups.size();
				groups.add(chosen);
			}
			chosen.members.add(hit.index);
			chosen.witnesses.add(hit.index);
			groupOf[hit.index]=selected;
		}
		// Provenance-only witnesses attach after grouping and never alter membership.
		for(int i=from; i<to; i++){
			final Hit hit=ordered[i];
			if(groupOf[hit.index]>=0 || state[hit.index]==INVALID || state[hit.index]==REJECTED || state[hit.index]==COVERED){continue;}
			neighborGroups(hit, nearby);
			for(int j=0; j<nearby.size; j++){groups.get(nearby.get(j)).witnesses.add(hit.index);}
		}
		groups.sort((a, b)->compare(a.representative, b.representative));
		return groups;
	}
	/** Uses the existing direct-overlap graph, not an all-loci scan. */
	private void neighborGroups(final Hit hit, final IntList result){
		assert(hit!=null && result!=null) : "Direct graph traversal uses a retained event and a caller-owned reusable index list";
		result.clear();
		if(hit.neighbors==null){return;}
		for(int i=0; i<hit.neighbors.size; i++){
			final Hit other=hits[hit.neighbors.get(i)];
			final int group=groupOf[other.index];
			if(group>=0 && hit.strand==other.strand && sameFamily(hit, other)){result.add(group);}
		}
		result.sort();
		int kept=0, previous=-1;
		for(int i=0; i<result.size; i++){
			final int value=result.get(i);
			if(value!=previous){result.set(kept++, value);}
			previous=value;
		}
		result.size=kept;
	}
	private int rank(final Hit a, final Hit b){
		int c=Boolean.compare(state[b.index]==CONFIDENT, state[a.index]==CONFIDENT);
		if(c==0){c=Boolean.compare(b.included, a.included);}
		if(c==0){c=Boolean.compare(b.full(), a.full());}
		if(c==0){c=Double.compare(b.bits/roster.query(b.query).ga, a.bits/roster.query(a.query).ga);}
		return c==0 ? compare(a, b) : c;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Serialization      ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String input, final String attributes, final String attributeEvidence, final String priorPath,
		final String census, final String memberships, final String scope, final String assets, final String queries){
		try{Files.createDirectory(root);}
		catch(IOException failure){throw new IllegalStateException("Could not create reference output directory", failure);}
		writeParameters();
		try(Output capture=new Output(path("capture_thresholds.tsv"))){
			capture.writer.println("native_directory\t"+NcrnaReferenceBatch.COMPLETION_HEADER);
			for(final byte[] row:captureRows){capture.writer.println(row);}
		}
		final String[] inputs={input, attributes, attributeEvidence, priorPath, census, memberships, scope, assets, queries};
		writeSources(inputs);
		writeConfig(inputs);
		final ArrayList<Row> rows=new ArrayList<Row>();
		int next=0;
		try(Output witnesses=new Output(path("annotation_witnesses.tsv.gz"))){
			witnesses.writer.println("annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation");
			for(int inputIndex=0; inputIndex<catalog.size(); inputIndex++){
				final NcrnaEvidenceCatalog.Genome genome=catalog.get(inputIndex);
				rows.clear();
				decodedTargets.clear();
				final PriorInput owner=priorInputs.get(genome.id);
				if(owner!=null){
					for(final Prior prior:owner.accepted){
						final Row row=priorRow(prior);
						rows.add(row);
						inputCounts[inputIndex][0]++;
						inputCounts[inputIndex][2]++;
						modelFamilies[prior.model.queryIndex()].scored++;
						modelFamilies[prior.model.queryIndex()].priors++;
						priorWitness(witnesses, row, prior);
						for(int i=0; i<prior.hits.size; i++){eventWitness(witnesses, row, hits[prior.hits.get(i)], prior.line);}
					}
				}
				int end=next;
				while(end<ordered.length && ordered[end].input.index==inputIndex){end++;}
				for(final Group group:groups(next, end)){
					final Hit hit=group.representative;
					final boolean mask=state[hit.index]==MARGINAL;
					final Hit[] members=new Hit[group.witnesses.size];
					for(int i=0; i<members.length; i++){members[i]=hits[group.witnesses.get(i)];}
					Arrays.sort(members, (a, b)->compare(a, b));
					final Scores scores=measurements(members, 0, members.length);
					final Row row=codec.create(genome.id, new Locus(targetName(hit.target), roster.query(hit.query).scope.family,
						hit.start, hit.stop, hit.strand), mask, mask ? reason(hit) : null, roster.query(hit.query).accession,
						witnessPath, scores);
					rows.add(row);
					inputCounts[inputIndex][mask ? 1 : 0]++;
					if(mask){modelFamilies[hit.query].masked++;}
					else{modelFamilies[hit.query].scored++;}
					for(final Hit member:members){eventWitness(witnesses, row, member, -1);}
				}
				next=end;
				codec.write(path("input_"+(inputIndex+1)+".gff.gz"), genome.id, rows);
			}
		}
		require(next==ordered.length, "Every captured event must belong to an emitted input view or provenance-only group");
		writeEvents();
		writePriors();
		writeSourceProvenance();
		writeCounts();
		publishManifest();
	}

	private Scores measurements(final Hit[] entries, final int from, final int to){
		assert(from>=0 && from<to && to<=entries.length) : "A physical annotation must retain at least one original model witness";
		final long[] scores=new long[to-from];
		final String[] ga=new String[scores.length];
		buffer.clear();
		for(int i=from; i<to; i++){
			final Hit hit=entries[i];
			if(i>from){buffer.append(',');}
			buffer.append(roster.query(hit.query).accession);
			scores[i-from]=bits[hit.index];
			ga[i-from]=thresholds[hit.query].gaText;
		}
		return codec.scores(buffer.toString(), scores, ga);
	}
	private Row priorRow(final Prior prior){
		final Scores scores;
		if(prior.hits.size==0){
			scores=codec.scores(prior.model.accession, new long[]{Long.MIN_VALUE},
				new String[]{thresholds[prior.model.queryIndex()].gaText});
		}else{
			final Hit[] associated=new Hit[prior.hits.size];
			for(int i=0; i<associated.length; i++){associated[i]=hits[prior.hits.get(i)];}
			Arrays.sort(associated, (a, b)->compare(a, b));
			scores=measurements(associated, 0, associated.length);
		}
		return codec.create(prior.owner.input.genome.id, new Locus(targetName(prior.target.key), prior.model.scope.family,
			prior.start, prior.stop, prior.strand), false, null, prior.model.accession, witnessPath, scores);
	}
	private void eventWitness(final Output out, final Row row, final Hit hit, final long priorLine){
		BuildNcrnaReference.field(buffer.clear(), row.id).tab().append(hit.input.genome.id).append("\tCM_HIT\t").append(hit.index+1).tab();
		buffer.append(hit.table).tab().append(hit.line).tab().append(roster.query(hit.query).accession).tab()
			.append(hit.start).tab().append(hit.stop).tab().append(priorLine).tab()
			.append(priorLine>=0 ? "SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS"
				: (flags[hit.index]&BRIDGE)!=0 ? "TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE"
				: !primary(hit) ? "WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT"
				: "QUALIFYING_SAME_FAMILY_GROUP_MEMBER").nl();
		out.writer.print(buffer);
	}
	private void priorWitness(final Output out, final Row row, final Prior prior){
		BuildNcrnaReference.field(buffer.clear(), row.id).tab().append(row.input).append("\tPROTECTED_PRIOR\t.\t.\t.\t");
		out.writer.print(buffer.append(prior.model.accession).tab().append(prior.start).tab().append(prior.stop).tab()
			.append(prior.line).append("\tORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY\n"));
	}

	private void writeEvents(){
		try(Output out=new Output(path("event_provenance.tsv.gz"))){
			out.writer.println("event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision");
			for(final Hit hit:hits){
				final int i=hit.index;
				final Model model=roster.query(hit.query);
				final Threshold threshold=thresholds[hit.query];
				buffer.clear().append(i+1).tab().append(hit.input.genome.id).tab().append(hit.target.bytes).tab()
					.append(hit.start).tab().append(hit.stop).tab().append(hit.strand).tab().append(model.accession).tab().append(model.scope.family).tab()
					.append(hit.from).tab().append(hit.to).tab().append(model.length).tab().append(hit.truncation).tab().append(hit.left).tab().append(hit.right).tab();
				NcrnaReferenceGff.appendFixed(buffer, bits[i], 10);
				buffer.tab().append(threshold.gaText).tab();
				NcrnaReferenceGff.appendFixed(buffer, NcrnaReferenceGff.percentage(bits[i], threshold.gaUnits, threshold.gaScale), 1000000);
				buffer.tab().append(hit.included).tab().append(hit.included!=(bits[i]*10>=threshold.gaHigh)).tab().append(floorText).tab()
					.append(bits[i]*10>=threshold.floorHigh).tab().append(uncertain[i]).tab().append(STATES[state[i]]).tab().append(reason(hit)).tab();
				NcrnaReferenceStrata.names(buffer, sourceFlags[i]).tab().append(flags[i]).tab().append(hit.table).tab().append(hit.line)
					.append("\tPRINTED_0.1_BIT_NOT_INTERNAL_PRECISION\t6_DECIMAL_PLACES_SOURCE_PAIR_RETAINED\n");
				out.writer.print(buffer);
			}
		}
	}
	private void writePriors(){
		try(Output out=new Output(path("prior_provenance.tsv.gz"))){
			out.writer.println("prior_line\tretained_disposition\t"+NcrnaReferencePrior.ROW_HEADER);
			for(final Prior prior:priors){
				out.writer.print(buffer.clear().append(prior.line).tab().append(prior.accepted ? "PROTECTED_SCORED" : "ORIGINAL_UNRESOLVED_PROVENANCE")
					.tab().append(prior.original).nl());
			}
		}
	}

	/** Original source rows remain claims, not copied reference genes. Counts use
	 * source-row units, including repeated source descriptions of one physical locus.
	 * A known RFAM family (or literal tRNA/tmRNA type) partitions rows; unresolved
	 * family claims remain UNASSIGNED. Only direct same-strand model links contribute. */
	private void writeSourceProvenance(){
		final ArrayList<SourceRow> rows=strata.sourceRows();
		final String[] family=new String[rows.size()];
		final int[] support=new int[rows.size()];
		for(final SourceRow row:rows){family[row.index]=sourceFamily(row);}
		final IntList sources=strata.linkedSources(), linked=strata.linkedHits();
		try(Output out=new Output(path("source_event_links.tsv.gz"))){
			out.writer.println("input_id\tgff\tsource_line\tevent_id\tsame_strand_valid\tfamily_agreement\tevent_tier");
			for(int i=0; i<sources.size; i++){
				final SourceRow row=rows.get(sources.get(i));
				final Hit hit=hits[linked.get(i)];
				final boolean valid=row.assessable && hit.valid() && row.strand==hit.strand;
				final boolean unknown=family[row.index].equals("UNASSIGNED");
				final boolean same=family[row.index].equals(roster.query(hit.query).scope.family);
				if(valid){
					support[row.index]|=1;
					if(same || unknown){support[row.index]|=1<<(state[hit.index]+1);}
				}
				buffer.clear().append(hit.input.genome.id).tab().append(row.source.path).tab().append(row.line).tab()
					.append(hit.index+1).tab().append(valid).tab().append(unknown ? "UNASSIGNED_SOURCE" : same ? "SAME_FAMILY" : "OTHER_FAMILY").tab().append(STATES[state[hit.index]]).nl();
				out.writer.print(buffer);
			}
		}
		final ObjectMap<String,long[]> counts=new ObjectMap<String,long[]>(String.class, long[].class);
		final ArrayList<String> names=new ArrayList<String>();
		try(Output out=new Output(path("source_provenance.tsv.gz"))){
			out.writer.println("input_id\tgff\tsource_line\ttarget\tsource_type\tstart\tstop\tstrand\tdeclared_family\tselection_bits\tdisposition\tsource_flags\tqualification");
			for(final SourceRow row:rows){
				final String name=family[row.index];
				long[] totals=counts.get(name);
				if(totals==null){
					totals=new long[SOURCE_STATES.length];
					counts.put(name, totals);
					names.add(name);
				}
				final int flags=support[row.index];
				final int disposition;
				if((row.selection&3)==0){disposition=0;}
				else if(!row.assessable){disposition=1;}
				else if(protectedSource(row, name)){disposition=2;}
				else if((flags&(1<<(CONFIDENT+1)))!=0){disposition=3;}
				else if((flags&(1<<(MARGINAL+1)))!=0){disposition=4;}
				else if((flags&(1<<(UNCERTAIN+1)))!=0){disposition=5;}
				else if((flags&(1<<(REJECTED+1)))!=0){disposition=6;}
				else if((flags&(1<<(COVERED+1)))!=0){disposition=7;}
				else if((flags&(1<<(BELOW+1)))!=0){disposition=8;}
				else{disposition=flags==0 ? 9 : 10;}
				totals[disposition]++;
				buffer.clear().append(row.source.input.state.genome.id).tab().append(row.source.path).tab().append(row.line).tab()
					.append(row.target.bytes).tab().append(row.type.bytes).tab().append(row.start).tab().append(row.stop).tab().append(row.strand).tab()
					.append(name).tab().append(row.selection).tab().append(SOURCE_STATES[disposition]).tab();
				NcrnaReferenceStrata.names(buffer, row.flags).append("\tSOURCE_CLAIM_NOT_INDEPENDENT_BIOLOGICAL_TRUTH\n");
				out.writer.print(buffer);
			}
		}
		try(Output out=new Output(path("source_family_counts.tsv"))){
			buffer.clear().append("declared_family\tsource_rows");
			for(final String status:SOURCE_STATES){buffer.tab().append(status);}
			out.writer.print(buffer.nl());
			long all=0;
			for(final String name:names){
				final long[] totals=counts.get(name);
				long total=0;
				for(final long n:totals){total+=n;}
				all+=total;
				buffer.clear().append(name).tab().append(total);
				for(final long n:totals){buffer.tab().append(n);}
				out.writer.print(buffer.nl());
			}
			require(all==geometry.profiles().retained(), "Every original retained source row needs one provenance disposition");
		}
	}
	private String sourceFamily(final SourceRow row){
		String family=null;
		for(final Model model:row.models){
			if(model==null){return "UNASSIGNED";}
			if(family!=null && !family.equals(model.scope.family)){return "UNASSIGNED";}
			family=model.scope.family;
		}
		if(family!=null){return family;}
		if(shared.Tools.equals(row.type.bytes, TRNA)){return "tRNA";}
		if(shared.Tools.equals(row.type.bytes, TMRNA)){return "tmRNA";}
		return "UNASSIGNED";
	}
	private boolean protectedSource(final SourceRow row, final String family){
		final PriorInput input=priorInputs.get(row.source.input.state.genome.id);
		final PriorTarget target=input==null ? null : input.targets.get(row.target);
		if(target==null){return false;}
		for(final Prior prior:target.accepted){
			if(prior.strand==row.strand && prior.start<=row.stop && prior.stop>=row.start && prior.model.scope.family.equals(family)){return true;}
		}
		return false;
	}
	private void writeParameters(){
		try(Output out=new Output(path("parameters.tsv"))){
			out.writer.println("parameter\tvalue");
			out.writer.println("schema\tNCRNA_REFERENCE_TIERS_V1");
			out.writer.println("marginal_percent\t"+floorText);
			out.writer.println("collapse_percent\t"+collapseText);
			out.writer.println("collapse_parameter_basis\tCHOSEN_GEOMETRY_NOT_MEASURED_OPTIMUM");
			out.writer.println("collapse_rule\tALL_MEMBER_RECIPROCAL_OVERLAP_QUALIFYING_SUPPORT_ONLY");
			out.writer.println("representative_order\tCONFIDENT_NATIVE_INCLUDED_FULL_MODEL_PERCENT_GA_COORDINATES_MODEL_NATIVE_KEY");
			out.writer.println("raw_event_tiers\tEVIDENCE_ELIGIBILITY_BEFORE_LOCUS_COLLAPSE_BRIDGES_FLAGGED_SEPARATELY");
			out.writer.println("captured_reporting_floor_percent\t75");
			out.writer.println("confident_score_basis\tORIGINAL_GA_BOUND_NATIVE_INCLUSION");
			out.writer.println("marginal_score_basis\tPRINTED_SCORE_WITH_0.05_BIT_BOUNDS");
			out.writer.println("floor_boundary_action\tUNCERTAIN_UNSCORED_UNMASKED");
			out.writer.println("printed_score_precision_bits\t0.1");
			out.writer.println("percentage_display\t6_DECIMAL_PLACES_HALF_AWAY_FROM_ZERO");
			out.writer.println("percentage_exact_basis\tSTORED_PRINTED_BITS_AND_ORIGINAL_GA");
			out.writer.println("prior_rule\tPRESERVE_NAMED_ACCEPTED_CORE_OR_WHOLE_INTERVAL_NO_WIDENING");
			out.writer.println("source_flag_dictionary\tNAMED_FLAGS_INLINE_IN_PROVENANCE");
			out.writer.println("decision_flag_dictionary\t1:OTHER_FAMILY_COMPETITOR;4:SAME_FAMILY_PRIOR;8:OTHER_FAMILY_PRIOR;16:TANDEM_BRIDGE_WITNESS");
			out.writer.println("rethresholding\tREBUILD_FROM_RETAINED_SAVED_EVIDENCE_AND_PRIORS_NO_NATIVE_SEARCH");
		}
	}
	private void writeSources(final String[] paths){
		try(Output out=new Output(path("sources.tsv"))){
			out.writer.println("role\tpath\tretention");
			for(int i=0; i<paths.length; i++){
				buffer.clear().append(SOURCE_ROLES[i]).tab();
				BuildNcrnaReference.field(buffer, paths[i]==null ? "SAME_AS_EVIDENCE" : paths[i]);
				out.writer.print(buffer.append("\tDURABLE_RETHRESHOLD_DEPENDENCY\n"));
			}
		}
	}
	/** Standard BBTools config is a reproducible re-threshold recipe. A caller must
	 * override out= with a fresh directory; later command-line floor values win. */
	private void writeConfig(final String[] paths){
		try(Output out=new Output(path("rethreshold.config"))){
			out.writer.println("mode=tiers");
			for(int i=0; i<paths.length; i++){
				if(paths[i]==null){continue;}
				require(rawSafe(paths[i]), "Re-threshold config paths must fit one literal configuration line");
				out.writer.print((CLI_KEYS[i]+'='+new File(paths[i]).getAbsolutePath()+'\n').getBytes(StandardCharsets.UTF_8));
			}
			out.writer.print(("out="+root+'\n').getBytes(StandardCharsets.UTF_8));
			out.writer.println("marginalpercent="+floorText);
			out.writer.println("collapsepercent="+collapseText);
			out.writer.println("ow=f");
			out.writer.println("t=1");
		}
	}
	private void writeCounts(){
		try(Output out=new Output(path("event_counts.tsv"))){
			buffer.clear().append("model\traw_events");
			for(final String tier:STATES){buffer.tab().append(tier);}
			out.writer.print(buffer.nl());
			for(int i=0; i<eventCounts.length; i++){
				long total=0;
				for(final long n:eventCounts[i]){total+=n;}
				require(total==(i==roster.queried() ? hits.length : geometry.profiles().modelEvents(i)), "Tier counts must conserve every original raw event");
				buffer.clear().append(i==roster.queried() ? "ALL_QUERIED" : roster.query(i).accession).tab().append(total);
				for(final long n:eventCounts[i]){buffer.tab().append(n);}
				out.writer.print(buffer.nl());
			}
		}
		try(Output out=new Output(path("summary.tsv"))){
			out.writer.println("input_id\tscored_annotations\tmasked_annotations\tprotected_priors\tstatus");
			for(int i=0; i<inputCounts.length; i++){
				buffer.clear().append(catalog.get(i).id);
				for(final long n:inputCounts[i]){buffer.tab().append(n);}
				out.writer.print(buffer.append("\tDECLARED_POLICY_APPLIED_NOT_INDEPENDENT_BIOLOGICAL_VALIDATION\n"));
			}
		}
		try(Output out=new Output(path("family_counts.tsv"))){
			buffer.clear().append("family\traw_events");
			for(final String tier:STATES){buffer.tab().append(tier);}
			out.writer.print(buffer.append("\tscored_annotations\tmasked_annotations\tprotected_priors\n"));
			long total=0, scored=0, masked=0, protectedCount=0;
			for(final Family family:families){
				long count=0;
				for(final long n:family.events){count+=n;}
				total+=count;
				scored+=family.scored;
				masked+=family.masked;
				protectedCount+=family.priors;
				buffer.clear().append(family.name).tab().append(count);
				for(final long n:family.events){buffer.tab().append(n);}
				out.writer.print(buffer.tab().append(family.scored).tab().append(family.masked).tab().append(family.priors).nl());
			}
			for(final long[] counts:inputCounts){
				scored-=counts[0];
				masked-=counts[1];
				protectedCount-=counts[2];
			}
			require(total==hits.length && scored==0 && masked==0 && protectedCount==0,
				"Family raw-event counts and physical annotations must separately conserve their populations");
		}
	}
	private void publishManifest(){
		final Path pending=root.resolve("reference.partial.tsv");
		try(Output out=new Output(pending.toString())){
			out.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			for(int i=0; i<catalog.size(); i++){
				final NcrnaEvidenceCatalog.Genome genome=catalog.get(i);
				require(rawSafe(genome.id) && rawSafe(genome.fasta), "Reference manifest identities must be representable as literal TSV fields");
				out.writer.print(buffer.clear().append(genome.id.getBytes(StandardCharsets.UTF_8)).tab()
					.append(genome.fasta.getBytes(StandardCharsets.UTF_8)).tab()
					.append(path("input_"+(i+1)+".gff.gz").getBytes(StandardCharsets.UTF_8)).nl());
			}
		}
		try{Files.move(pending, root.resolve("reference.tsv"), StandardCopyOption.ATOMIC_MOVE);}
		catch(IOException failure){throw new IllegalStateException("Could not publish the completed reference manifest", failure);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers and State   ----------------*/
	/*--------------------------------------------------------------*/

	private String reason(final Hit hit){
		final int i=hit.index;
		if(state[i]==INVALID){return "invalid_observed_interval";}
		if(state[i]==REJECTED){return "conflicts_with_reviewed_family";}
		if(state[i]==COVERED){return "overlaps_protected_prior_no_widening";}
		if((flags[i]&BRIDGE)!=0){return "tandem_bridge_witness_not_representative";}
		if(state[i]==CONFIDENT){return "native_GA_full_model_clean_family";}
		if(state[i]==UNCERTAIN){return "printed_rounding_interval_touches_floor";}
		if(state[i]==BELOW){return "below_selected_marginal_floor";}
		if((flags[i]&COMPETITOR)!=0){return "qualifying_competing_CM";}
		if((sourceFlags[i]&SOURCE_BLOCKERS)!=0){return "discordant_source_evidence";}
		if(!hit.full()){return "partial_model_support";}
		return "below_GA_above_marginal_floor";
	}
	private int compare(final Hit a, final Hit b){
		int c=Integer.compare(a.input.index, b.input.index);
		if(c==0){c=shared.Tools.compare(a.target.bytes, b.target.bytes);}
		if(c==0){c=Character.compare(a.strand, b.strand);}
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		if(c==0){c=roster.query(a.query).scope.family.compareTo(roster.query(b.query).scope.family);}
		if(c==0){c=roster.query(a.query).accession.compareTo(roster.query(b.query).accession);}
		if(c==0){c=Double.compare(b.bits, a.bits);}
		if(c==0){c=Boolean.compare(b.included, a.included);}
		if(c==0){c=Integer.compare(a.from, b.from);}
		if(c==0){c=Integer.compare(a.to, b.to);}
		if(c==0){c=Integer.compare(a.truncation, b.truncation);}
		if(c==0){c=shared.Tools.compare(a.table, b.table);}
		return c==0 ? Long.compare(a.line, b.line) : c;
	}
	private String text(final int field){
		parser.setBounds(field);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private String decode(final byte[] bytes){return decoder.decode(bytes, 0, bytes.length);}
	private String targetName(final Target target){
		String name=decodedTargets.get(target);
		if(name==null){
			name=decode(target.bytes);
			decodedTargets.put(target, name);
		}
		return name;
	}
	private boolean same(final int field, final byte[] expected){
		parser.setBounds(field);
		if(parser.b()-parser.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(parser.line()[parser.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private static boolean rawSafe(final String value){return value.indexOf('\t')<0 && value.indexOf('\r')<0 && value.indexOf('\n')<0;}
	private String path(final String name){return root.resolve(name).toString();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Threshold{
		Threshold(final String text, final BigDecimal percent){
			gaText=text;
			final BigDecimal ga=new BigDecimal(text);
			require(ga.signum()>0 && ga.precision()<=15 && ga.scale()>=0, "Captured GA must be a positive supported plain decimal");
			gaUnits=ga.unscaledValue().longValueExact();
			gaScale=ga.scale();
			final BigDecimal selected=ga.multiply(percent);
			floorLow=selected.setScale(0, RoundingMode.FLOOR).longValueExact();
			floorHigh=selected.setScale(0, RoundingMode.CEILING).longValueExact();
			gaHigh=ga.multiply(BigDecimal.valueOf(100)).setScale(0, RoundingMode.CEILING).longValueExact();
		}
		final String gaText;
		final long floorLow, floorHigh, gaHigh, gaUnits;
		final int gaScale;
	}
	private static final class PriorInput{
		PriorInput(final NcrnaReferenceGeometry.InputState input_){
			input=input_;
			fasta=BuildNcrnaReference.field(new ByteBuilder(), input.genome.fasta).toBytes();
		}
		final NcrnaReferenceGeometry.InputState input;
		final byte[] fasta;
		final ObjectMap<Target,PriorTarget> targets=new ObjectMap<Target,PriorTarget>(Target.class, PriorTarget.class);
		final ArrayList<Prior> accepted=new ArrayList<Prior>();
	}
	private static final class PriorTarget{
		PriorTarget(final Target key_, final int length_){
			key=key_;
			length=length_;
		}
		final Target key;
		final int length;
		boolean observed;
		final ArrayList<Prior> accepted=new ArrayList<Prior>();
	}
	private static final class Prior{
		Prior(final PriorInput owner_, final PriorTarget target_, final Model model_, final int start_, final int stop_,
			final char strand_, final boolean accepted_, final long line_, final byte[] original_){
			owner=owner_;
			target=target_;
			model=model_;
			start=start_;
			stop=stop_;
			strand=strand_;
			accepted=accepted_;
			line=line_;
			original=original_;
		}
		final PriorInput owner;
		final PriorTarget target;
		final Model model;
		final int start, stop;
		final char strand;
		final boolean accepted;
		final long line;
		final byte[] original;
		final IntList hits=new IntList();
	}
	private static final class Family{
		Family(final String name_){name=name_;}
		final String name;
		final long[] events=new long[STATES.length];
		long scored, masked, priors;
	}
	private static final class Group{
		Group(final Hit representative_){representative=representative_;}
		final Hit representative;
		final IntList members=new IntList(), witnesses=new IntList();
	}
	private final NcrnaEvidenceCatalog catalog;
	private final Roster roster;
	private final NcrnaReferenceGeometry geometry;
	private final Hit[] hits, ordered;
	private final long[] bits;
	private final byte[] state;
	private final int[] flags, groupOf;
	private int[] sourceFlags;
	private NcrnaReferenceStrata strata;
	private final boolean[] qualifies, uncertain;
	private final Threshold[] thresholds;
	private final BigDecimal floor;
	private final String floorText, collapseText, witnessPath;
	private final long collapseMillionths;
	private final Path root;
	private final long[][] eventCounts, inputCounts;
	private final ObjectMap<String,PriorInput> priorInputs=new ObjectMap<String,PriorInput>(String.class, PriorInput.class);
	private final ArrayList<Prior> priors=new ArrayList<Prior>();
	private final ArrayList<byte[]> captureRows=new ArrayList<byte[]>();
	private final ArrayList<Family> families=new ArrayList<Family>();
	private final ObjectMap<String,Family> familyByName=new ObjectMap<String,Family>(String.class, Family.class);
	private final Family[] modelFamilies;
	private final ObjectMap<Target,String> decodedTargets=new ObjectMap<Target,String>(Target.class, String.class);
	private final Target probe=new Target();
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final NcrnaReferenceGff codec=new NcrnaReferenceGff();
	private final ByteBuilder buffer=new ByteBuilder();
	private static final byte CONFIDENT=0, MARGINAL=1, BELOW=2, UNCERTAIN=3, REJECTED=4, COVERED=5, INVALID=6;
	private static final int COMPETITOR=1, PRIOR_SAME=4, PRIOR_OTHER=8, BRIDGE=16;
	private static final int SOURCE_BLOCKERS=NcrnaReferenceStrata.OTHER_FAMILY_SOURCE|NcrnaReferenceStrata.UNKNOWN_RFAM
		|NcrnaReferenceStrata.PSEUDO_TRUE|NcrnaReferenceStrata.PSEUDOGENE_TAG|NcrnaReferenceStrata.PSEUDO_BIOTYPE
		|NcrnaReferenceStrata.PARTIAL_TRUE|NcrnaReferenceStrata.EXCEPTION|NcrnaReferenceStrata.UNKNOWN_BOOLEAN|NcrnaReferenceStrata.INCOMPATIBLE_ID;
	private static final String[] STATES={"CONFIDENT", "MARGINAL", "BELOW_MARGINAL", "UNCERTAIN_FLOOR", "REJECTED_WRONG_FAMILY", "COVERED_BY_PROTECTED_PRIOR", "INVALID_INTERVAL"};
	private static final String[] SOURCE_ROLES={"evidence_prefix", "captured_attributes", "attribute_evidence_prefix", "protected_priors",
		"census", "memberships", "scope", "captured_model_metadata", "queried_models"};
	private static final String[] CLI_KEYS={"in", "attributes", "attributesevidence", "priors", "census", "memberships", "scope", "assets", "queries"};
	private static final String[] SOURCE_STATES={"CONTEXT_PROVENANCE", "UNASSESSABLE_SOURCE", "PROTECTED_PRIOR_SOURCE",
		"SUPPORTED_CONFIDENT_EVENT", "SUPPORTED_MARGINAL_EVENT", "UNCERTAINTY_UNMASKED", "REMOVED_REVIEWED_FAMILY_CONFLICT",
		"PROVENANCE_NEAR_PROTECTED_PRIOR", "REMOVED_BELOW_MARGINAL", "REMOVED_NO_CM_SUPPORT", "REMOVED_OTHER_FAMILY_ONLY"};
	private static final byte[] TRNA="tRNA".getBytes(StandardCharsets.US_ASCII), TMRNA="tmRNA".getBytes(StandardCharsets.US_ASCII);
	private static final String BATCH_SOURCES_HEADER="batch_id\tprepared_directory\tnative_directory\tfasta\tfasta_sha80\tmapping\tmapping_sha80\tmembers\tmembers_sha80\tsearch_space_Mb\tobservation_basis";
}
