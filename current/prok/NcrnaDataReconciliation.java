package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectMap;
import map.ObjectSet;
import map.IntObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Reconciles saved RF00174 evidence, assembly metadata and member coordinates.
 * Derived from template.A_SampleByteFile; one metadata worker, native byte I/O.
 * Never opens genome sequences, clusters sequences, sketches or calls genes.
 * Exact-coordinate groups and same-frame overlap components remain separate
 * outputs; neither proves a biological identity across renamed/reversed frames.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDataReconciliation{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaDataReconciliation tool=new NcrnaDataReconciliation(args);
		try{
			tool.process();
			timer.stop();
			tool.outstream.println("NCRNA_DATA_RECONCILIATION_COMPLETE_METADATA_ONLY time="+timer);
		}finally{Shared.closeStream(tool.outstream);}
	}

	public NcrnaDataReconciliation(final String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("groups")){groups=b;}
			else if(a.equals("mapping")){mapping=b;}
			else if(a.equals("roots")){roots=b.split(",", -1);}
			else if(a.equals("genomes")){genomes=b;}
			else if(a.equals("assignments")){assignments=b;}
			else if(a.equals("ledger")){ledger=b;}
			else if(a.equals("proxy")){proxy=b;}
			else if(a.equals("recovered")){recovered=b;}
			else if(a.equals("overlays")){overlays=b.split(",", -1);}
			else if(a.equals("expectedgroups")){expectedGroups=Parse.parseIntKMG(b);}
			else if(a.equals("expectedcohorts")){expectedCohorts=Parse.parseIntKMG(b);}
			else if(a.equals("expectedreference")){expectedReference=Parse.parseIntKMG(b);}
			else if(a.equals("expectedconfident")){expectedConfident=Parse.parseIntKMG(b);}
			else if(a.equals("expectedkept")){expectedKept=Parse.parseIntKMG(b);}
			else if(a.equals("expectedassigned")){expectedAssigned=Parse.parseIntKMG(b);}
			else if(a.equals("expectedgenomes")){expectedGenomes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmodels")){expectedModels=Parse.parseIntKMG(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		require(in!=null && out!=null && groups!=null && mapping!=null && roots!=null && genomes!=null && assignments!=null
			&& ledger!=null && proxy!=null && recovered!=null && overlays!=null,
			"Supply in=cohorts.tsv groups= mapping= roots= genomes= assignments= ledger= proxy= recovered= overlays= out=fresh_directory");
		require(!parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1
			&& parser.in2==null && parser.out2==null, "Reconciliation requires complete metadata and fresh outputs");
	}

	private void process(){
		assert(loci.isEmpty() && members.isEmpty()) : "One reconciler owns one complete metadata generation";
		readGroups();
		require(!new File(out).exists() && new File(out).mkdir(), "Output directory must be fresh: "+out);
		try(Output parameters=new Output(path("parameters.tsv"))){
			parameters.writer.println("format_version\t2");
			parameters.writer.println("qualification\tNative RF00174 event has valid one-based coordinates and native_included=true (the saved native GA decision), independent of prior annotation, tier, truncation or representative model.");
			parameters.writer.println("duplicates\tExact groups retain their original bounds; components merge overlapping inclusive qualified intervals on the same assembly/frame/contig/strand, never adjacency or unproved aliases.");
			parameters.writer.println("inventory\tCurrent live collection from the named roots; historical 51,079-path roster not consulted; no genome payload opened.");
			parameters.writer.println("member_matches\tExact bounds and reciprocal80 overlap are separate; no proven match with unresolved mapping is not a gene-absence claim.");
			parameters.writer.println("scan_coverage\tDeclared census binding or accepted assembly-group alias; not exact-current-version sequence validation.");
			parameters.writer.println("frame_cost\tDimensions copied from saved scan censuses once per canonical frame; no sequence comparison performed.");
		}
		readCohorts();
		inventory();
		try(Output annotations=new Output(path("annotations.tsv.gz")); Output events=new Output(path("events.tsv.gz"));
			Output links=new Output(path("witnesses.tsv.gz")); Output counts=new Output(path("cohorts.tsv"))){
			annotations.writer.println("cohort\tannotation_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tsaved_scored_RF00174\tqualifying_RF00174_group_events\tconfident_RF00174_group_events\tnative_GA_RF00174_events\tprotected_prior\tdisposition");
			events.writer.println("cohort\tassembly_group\tnative_GA_qualified\tprior_annotation_overlap\tqualifying_annotation_links\t"+EVENT_HEADER);
			links.writer.println("cohort\t"+WITNESS_HEADER);
			counts.writer.println("cohort\tcensus_inputs\tsaved_scored_RF00174_annotations\traw_RF00174_events\tnative_GA_RF00174_events");
			for(final Cohort cohort:cohorts){
				final int beforeRaw=rawEvents, beforeGa=nativeGaEvents;
				final int scored=cohort(cohort, annotations, events, links);
				counts.writer.print(b.clear().append(utf8(cohort.label)).tab().append(cohort.inputs.size()).tab()
					.append(scored).tab().append(rawEvents-beforeRaw).tab().append(nativeGaEvents-beforeGa).nl());
			}
		}
		require(originalScored==expectedReference && importConfident==expectedConfident,
			"Saved population checks differ: reference="+originalScored+" expected="+expectedReference
			+", import confident events="+importConfident+" expected="+expectedConfident);
		indexLoci();
		readMembers();
		positions(proxy, 0);
		positions(recovered, 1);
		for(final String overlay:overlays){positions(overlay, 2);}
		compareMembers();
		writeLoci();
		writeComponents();
		writeFrames();
		writeSummary();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Assembly Metadata       ----------------*/
	/*--------------------------------------------------------------*/

	private void readGroups(){
		try(Input input=new Input(groups)){
			RrnaResourceIO.header(input.reader, NcrnaConstructionMembers.GROUP_HEADER, groups);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==14, "Accepted assembly audit requires fourteen fields");
				final int group=p.parseInt(10);
				require(group>=0, "Accepted assembly groups must be nonnegative");
				final Assembly assembly=new Assembly(group, text(2), text(3), p.termEquals("true", 12));
				bind(assembly.literal, assembly);
				bind(assembly.canonical, assembly);
				auditedGroups.put(group, 1);
				if(assembly.uncertain){uncertainGroups.put(group, 1);}
				nextGroup=Math.max(nextGroup, group+1);
			}
		}
		require(auditedGroups.size()==expectedGroups, "Original accepted alias-audit count differs: "+auditedGroups.size());
	}

	/** Bind all declared census paths before asking which current live paths were
	 * scanned. The historical alias audit is an independent, unchanged check. */
	private void readCohorts(){
		require(expectedCohorts>=2, "Retain original and prior-import cohorts before adding new scan increments");
		final ObjectSet<String> labels=new ObjectSet<String>(String.class);
		try(Input input=new Input(in)){
			RrnaResourceIO.header(input.reader, COHORT_HEADER, in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Every scan cohort requires seven explicit manifest/count fields");
				final Cohort cohort=new Cohort(text(0), text(1), text(2), text(3), text(4), p.parseInt(5), p.parseInt(6));
				require(!labels.contains(cohort.label), "Scan cohort names must be unique");
				labels.add(cohort.label);
				cohorts.add(cohort);
			}
		}
		require(cohorts.size()==expectedCohorts, "Declared scan cohort count differs: "+cohorts.size());
		require(labels.contains("original") && labels.contains("import"), "Name the original and audited prior-import cohorts explicitly; row order is not a cohort role");
		for(final Cohort cohort:cohorts){
			try(Input input=new Input(cohort.census)){
				RrnaResourceIO.header(input.reader, CENSUS_HEADER, cohort.census);
				for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
					p.set(row);
					require(p.terms()==5, "Saved scan census requires all five fields");
					final String id=text(0), fasta=text(1), canonical=sourceFrame(fasta);
					Assembly assembly=assemblies.get(fasta);
					if(assembly==null){assembly=assemblies.get(canonical);}
					if(assembly==null){assembly=new Assembly(nextGroup++, fasta, canonical, false);}
					bind(fasta, assembly);
					bind(canonical, assembly);
					final Assembly frame=new Assembly(assembly.group, fasta, canonical, assembly.uncertain);
					require(!cohort.inputs.contains(id), "A cohort input ID cannot be repeated: "+id);
					cohort.inputs.put(id, frame);
					cohort.inputIds.add(id);
					scannedGroups.put(assembly.group, 1);
					final long bytes=p.parseLong(2), bases=p.parseLong(3), records=p.parseLong(4);
					require(bytes>=0 && bases>0 && records>0, "Saved census dimensions must describe a real scanned input");
					SourceFrame source=frames.get(canonical);
					if(source==null){
						source=new SourceFrame(frame, bytes, bases, records);
						frames.put(canonical, source);
						frameList.add(source);
						frameCounts.increment(frame.group, 1);
					}else{
						require(source.assembly.group==frame.group && source.bytes==bytes && source.bases==bases && source.records==records,
							"One canonical scan frame cannot have conflicting saved census dimensions: "+canonical);
					}
					source.inputs++;
				}
			}
			require(cohort.inputs.size()==cohort.expectedInputs, "Cohort census input count differs: "+cohort.label);
		}
	}

	private void bind(final String key, final Assembly assembly){
		assert(key!=null && assembly!=null) : "Source path binding must retain its accepted assembly row";
		final Assembly old=assemblies.get(key);
		require(old==null || old.group==assembly.group, "One source path cannot name two accepted groups: "+key);
		if(old==null){assemblies.put(key, assembly);}
	}

	/** Lists only immediate .fna.gz members of the named live roots; opens no payload. */
	private void inventory(){
		final ObjectSet<String> literalPaths=new ObjectSet<String>(String.class);
		final IntHashMap liveGroups=new IntHashMap();
		final IntHashMap missingGroups=new IntHashMap();
		try(Output inventory=new Output(path("assemblies.tsv.gz")); Output missing=new Output(path("unscanned.tsv.gz"))){
			inventory.writer.println("fasta\tcanonical_fasta\tassembly_group\tcoverage\tidentity_basis");
			missing.writer.println("assembly_group\tfasta\tcanonical_fasta\tstatus");
			for(final String root:roots){
				final File directory=new File(root);
				require(directory.isAbsolute() && directory.isDirectory(), "Each live source root must exist: "+root);
				final String[] names=directory.list((dir, name)->name.endsWith(".fna.gz"));
				require(names!=null && names.length>0, "Cannot enumerate the declared source root: "+root);
				Arrays.sort(names);
				for(final String name:names){
					final File file=new File(directory, name);
					final String literal=file.getAbsolutePath(), canonical=sourceFrame(literal);
					require(file.isFile() && !literalPaths.contains(literal), "Live source paths must be existing unique files: "+literal);
					literalPaths.add(literal);
					Assembly assembly=assemblies.get(literal);
					if(assembly==null){assembly=assemblies.get(canonical);}
					final boolean scanned=assembly!=null && scannedGroups.contains(assembly.group);
					if(assembly==null){
						assembly=new Assembly(nextGroup++, literal, canonical, false);
						bind(literal, assembly);
						bind(canonical, assembly);
					}
					liveGroups.put(assembly.group, 1);
					inventory.writer.print(b.clear().append(utf8(literal)).tab().append(utf8(canonical)).tab().append(assembly.group)
						.tab().append(scanned ? "SCANNED" : "NO_ACCEPTED_SCAN_BINDING").tab()
						.append(scanned ? "DECLARED_CENSUS_OR_ACCEPTED_GROUP_BINDING" : "CANONICAL_PATH_ONLY_OTHER_ALIASES_UNRESOLVED").nl());
					if(!scanned){
						missingPaths++;
						missing.writer.print(b.clear().append(assembly.group).tab().append(utf8(literal)).tab().append(utf8(canonical))
							.append("\tNO_DECLARED_SCAN_CENSUS_OR_ACCEPTED_GROUP_BINDING\n"));
						missingGroups.put(assembly.group, 1);
					}
				}
			}
		}
		livePaths=literalPaths.size();
		liveAssemblyGroups=liveGroups.size();
		unscannedGroups=missingGroups.size();
	}

	private static String canonical(final File file){
		try{return file.getCanonicalPath();}
		catch(IOException failure){throw new IllegalArgumentException("Cannot resolve source path: "+file, failure);}
	}

	private String sourceFrame(final String path){
		String resolved=canonicalPaths.get(path);
		if(resolved==null){
			resolved=canonical(new File(path));
			canonicalPaths.put(path, resolved);
		}
		return resolved;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Saved Call Evidence     ----------------*/
	/*--------------------------------------------------------------*/

	private int cohort(final Cohort cohort,
			final Output annotationOutput, final Output eventOutput, final Output linkOutput){
		final String label=cohort.label, reference=cohort.reference, eventPath=cohort.events, witnesses=cohort.witnesses;
		final ObjectMap<String, Assembly> inputs=cohort.inputs;
		final NcrnaSubtypeEvidence subtype=new NcrnaSubtypeEvidence(mapping, witnesses, eventPath, cohort.inputIds);
		int selectedType=-1;
		for(int i=0; i<subtype.size(); i++){
			if(FAMILY.equals(subtype.referenceModel(i))){selectedType=i;}
		}
		require(selectedType>=0, "Explicit member mapping must contain RF00174");
		final int selectedBit=1<<selectedType;
		final ObjectMap<String, Annotation> annotations=new ObjectMap<String, Annotation>(String.class, Annotation.class);
		final ArrayList<Annotation> rows=new ArrayList<Annotation>();
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ObjectSet<String> seenInputs=new ObjectSet<String>(String.class);
		int scored=0;
		try(Input input=new Input(reference)){
			RrnaResourceIO.header(input.reader, "input_id\tfasta\tgff", reference);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				final String id=text(0), fasta=text(1), gff=text(2);
				final Assembly assembly=inputs.get(id);
				require(p.terms()==3 && assembly!=null && assembly.literal.equals(fasta) && !seenInputs.contains(id),
					"Reference manifest must bind each exact census input once: "+id);
				seenInputs.add(id);
				for(final Row row:codec.read(gff, id)){
					if(!NcrnaSubtypeEvidence.FAMILY.equals(row.locus.family)){continue;}
					final NcrnaSubtypeEvidence.Support support=subtype.support(row);
					final boolean savedSupport=!row.mask && (support.scored&selectedBit)!=0;
					final Annotation annotation=new Annotation(row, assembly, savedSupport, support.priorType>=0);
					require(!annotations.contains(row.id), "Reference annotation IDs must be unique within a cohort");
					annotations.put(row.id, annotation);
					rows.add(annotation);
					if(savedSupport){scored++;}
				}
			}
		}
		subtype.finish();
		require(seenInputs.size()==cohort.expectedInputs && (cohort.expectedScored<0 || scored==cohort.expectedScored),
			"Saved reference population differs in "+label+": inputs="+seenInputs.size()+", scored="+scored);
		if(label.equals("original")){originalScored=scored;}
		final ObjectMap<String, Event> eventIds=new ObjectMap<String, Event>(String.class, Event.class);
		final ArrayList<Event> events=new ArrayList<Event>();
		try(Input input=new Input(eventPath)){
			RrnaResourceIO.header(input.reader, EVENT_HEADER, eventPath);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==30, "Event provenance must retain all thirty fields");
				if(!p.termEquals(FAMILY, 6)){continue;}
				final String id=text(0);
				final Assembly assembly=inputs.get(text(1));
				require(assembly!=null && !eventIds.contains(id), "Event identity must bind one input and native record: "+id);
				final boolean valid=p.parseInt(3)>0 && p.parseInt(4)>=p.parseInt(3)
					&& (p.termEquals("+", 5) || p.termEquals("-", 5));
				final Event event=new Event(line, assembly, p.termEquals("CONFIDENT", 22),
					valid && p.termEquals("true", 17), text(1), text(26), p.parseLong(27));
				eventIds.put(id, event);
				events.add(event);
				rawEvents++;
				if(event.qualified){nativeGaEvents++;}
				if(p.termEquals("MARGINAL", 22)){
					marginalEvents++;
					if(event.qualified){marginalGaEvents++;}
				}
				if(p.termEquals("UNCERTAIN_FLOOR", 22)){
					uncertainEvents++;
					if(event.qualified){uncertainGaEvents++;}
				}
				if(label.equals("import") && event.confident){importConfident++;}
			}
		}
		try(Input input=new Input(witnesses)){
			RrnaResourceIO.header(input.reader, WITNESS_HEADER, witnesses);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==11, "Annotation witnesses must retain eleven fields");
				if(!p.termEquals(FAMILY, 6)){continue;}
				linkOutput.writer.print(b.clear().append(utf8(label)).tab().append(line).nl());
				if(!p.termEquals("CM_HIT", 2)){continue;}
				final Event event=eventIds.get(text(3));
				require(event!=null, "Every RF00174 witness must retain its source event");
				require(event.input.equals(text(1)) && event.nativeTable.equals(text(4)) && event.nativeLine==p.parseLong(5),
					"RF00174 witness must identify the exact cohort input, native table and physical line");
				final Annotation annotation=annotations.get(decoded(0));
				if(annotation==null){continue;}//Other subtype's annotation remains event evidence only.
				event.prior|=annotation.prior;
				if(event.qualified){annotation.nativeGa++;}
				require(annotation.row.input.equals(text(1)) && annotation.assembly.group==event.assembly.group,
					"Annotation and event must retain the same cohort input and accepted assembly");
				if(p.termEquals("QUALIFYING_SAME_FAMILY_GROUP_MEMBER", 10)){
					annotation.support++;
					if(event.confident){annotation.confident++;}
					event.annotationLinks++;
				}
			}
		}
		for(final Annotation annotation:rows){
			final Row row=annotation.row;
			if(annotation.savedSupport){qualifiedAnnotations++;}
			if(annotation.savedSupport && annotation.prior && annotation.nativeGa==0){priorOnlyAnnotations++;}
			final String disposition=annotation.prior ? annotation.nativeGa>0 ? "PRIOR_WITH_NATIVE_GA_HIT" : "PRIOR_WITHOUT_NATIVE_GA_HIT_OUTSIDE_CM_ROSTER"
				: annotation.nativeGa>0 ? "ANNOTATION_WITH_NATIVE_GA_HIT" : "NO_NATIVE_GA_RF00174_SUPPORT";
			annotationOutput.writer.print(b.clear().append(utf8(label)).tab().append(utf8(row.id)).tab().append(annotation.assembly.group)
				.tab().append(utf8(coordinateFrame(annotation.assembly, row.locus.seqid))).tab().append(utf8(row.locus.seqid)).tab().append(row.locus.start).tab().append(row.locus.stop)
				.tab().append(row.locus.strand).tab().append(annotation.savedSupport).tab().append(annotation.support).tab()
				.append(annotation.confident).tab().append(annotation.nativeGa).tab().append(annotation.prior).tab().append(disposition).nl());
		}
		for(final Event event:events){
			if(event.qualified && event.prior){nativeGaPriorEvents++;}
			p.set(event.row);
			if(p.parseInt(3)>0 && p.parseInt(4)>=p.parseInt(3) && (p.termEquals("+", 5) || p.termEquals("-", 5))){
				final Locus locus=addLocus(event.assembly, decoded(2), p.parseInt(3), p.parseInt(4), strand(5));
				locus.events++;
				locus.qualified|=event.qualified;
				if(event.qualified){locus.nativeEvents++;}
				locus.prior|=event.prior && event.qualified;
			}
			eventOutput.writer.print(b.clear().append(utf8(label)).tab().append(event.assembly.group).tab()
				.append(event.qualified).tab().append(event.prior).tab()
				.append(event.annotationLinks).tab().append(event.row).nl());
		}
		return scored;
	}

	private Locus addLocus(final Assembly assembly, final String contig, final int start, final int stop, final char strand){
		require(start>0 && stop>=start, "Saved loci require one-based inclusive coordinates");
		final String frame=coordinateFrame(assembly, contig);
		final String key=coordinateKey(assembly.group, frame, contig, start, stop, strand);
		Locus locus=coordinates.get(key);
		if(locus==null){
			locus=new Locus(loci.size(), assembly.group, frame, contig, start, stop, strand);
			coordinates.put(key, locus);
			loci.add(locus);
		}
		return locus;
	}

	/** Whole-assembly multiset equality permits contig reversal and renaming.
	 * No string that resembles an accession substitutes for an actual coordinate map. */
	public static String coordinateFrame(final Assembly assembly, final String contig){
		assert(assembly!=null && contig!=null) : "Coordinates require the retained assembly and literal contig identity";
		return assembly.canonical;
	}

	private static String coordinateKey(final int group, final String frame, final String contig,
			final int start, final int stop, final char strand){
		return group+"\t"+frame+"\t"+contig+"\t"+strand+"\t"+start+"\t"+stop;
	}

	private void indexLoci(){
		for(final Locus locus:loci){
			final String key=bucket(locus.group, locus.contig, locus.strand);
			final Locus previous=byContig.get(key);
			locus.next=previous;
			byContig.put(key, locus);
			locus.nextGroup=byGroup.get(locus.group);
			byGroup.put(locus.group, locus);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------     Member Coordinates       ----------------*/
	/*--------------------------------------------------------------*/

	private void readMembers(){
		final ObjectMap<String, Assembly> sources=new ObjectMap<String, Assembly>(String.class, Assembly.class);
		try(Input input=new Input(genomes)){
			RrnaResourceIO.header(input.reader, NcrnaConstructionMembers.GENOME_HEADER, genomes);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==11, "Original source genomes require eleven fields");
				final String id=text(1);
				require(!sources.contains(id), "Original donor names must be unique");
				final Assembly bound=assemblies.get(text(7));
				require(bound!=null, "Every PGAP donor must occur in the complete accepted scan audit");
				sources.put(id, bound);
			}
		}
		require(sources.size()==expectedGenomes, "Original source-genome roster count differs: "+sources.size());
		final ObjectMap<String, Member> headers=new ObjectMap<String, Member>(String.class, Member.class);
		final ObjectSet<String> modelNames=new ObjectSet<String>(String.class);
		final ObjectSet<String> assignedHeaders=new ObjectSet<String>(String.class);
		int assigned=0;
		try(Input input=new Input(assignments)){
			RrnaResourceIO.header(input.reader, ASSIGNMENT_HEADER, assignments);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Assignments require all six fields");
				final String header=text(0);
				require(!assignedHeaders.contains(header), "Each complete builder header must have exactly one assignment");
				assignedHeaders.add(header);
				assigned++;
				if(!p.termEquals("kept", 5)){
					require((p.termEquals("orphan", 5) || p.termEquals("small", 5)) && text(4).isEmpty(),
						"Nonretained assignments need an explicit disposition and no model");
					continue;
				}
				final String model=RrnaResourceIO.first(text(4));
				require(model.startsWith("cobalamin_RF00174_"), "Kept members must bind an RF00174 library model");
				modelNames.add(model);
				final String id=RrnaResourceIO.first(header);
				final Assembly owner=sources.get(RrnaResourceIO.field(header, "source"));
				require(owner!=null && !memberIds.contains(id) && !headers.contains(header), "Each kept member needs one original donor and unique header");
				final Member member=new Member(id, header, owner);
				members.add(member);
				memberIds.put(id, member);
				headers.put(header, member);
			}
		}
		require(assigned==expectedAssigned && members.size()==expectedKept && modelNames.size()==expectedModels,
			"Original assignment, retained member or model counts differ");
		int joined=0;
		try(Input input=new Input(ledger)){
			RrnaResourceIO.header(input.reader, LEDGER_HEADER, ledger);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==8, "Original routing ledger requires eight fields");
				if(!p.termEquals(FAMILY, 3)){continue;}
				final String source=text(1), biological=text(2), raw=text(7);
				final String header=raw+(source.equals(FAMILY) ? " member_family="+biological : "")+" build_group="+FAMILY;
				final Member member=headers.get(header);
				if(member==null){continue;}//Orphans do not enter retained-member coordinate counts.
				require(!member.bound && member.id.equals(text(0)), "Each kept assignment needs its exact original routing witness");
				member.bound=true;
				member.length=p.parseInt(6);
					require(member.length>0, "A retained member's recorded length must be positive");
				joined++;
			}
		}
		require(joined==expectedKept, "Every retained member must join the original routing ledger");
	}

	/** Uses the same correction precedence as NcrnaDdlLocusOverlap.positions,
	 * without its rankable-query selection or any sequence/holdout classifier. */
	private void positions(final String inputPath, final int mode){
		try(Input input=new Input(inputPath)){
			RrnaResourceIO.header(input.reader, mode==0 ? NcrnaDdlLocusOverlap.RAW_HEADER
				: mode==1 ? NcrnaDdlLocusOverlap.RECOVERED_HEADER : NcrnaDdlLocusOverlap.OVERLAY_HEADER, inputPath);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==(mode==0 ? 21 : mode==1 ? 38 : 33), "Coordinate witness schema differs: "+inputPath);
				final Member member=memberIds.get(text(0));
				if(member==null){continue;}
				final String contig=text(6);
				final char strand=strand(7);
				require(RrnaResourceIO.field(member.header, "source").equals(text(1)), "A coordinate correction cannot change the original donor");
				if(mode==0){
					require(!member.raw, "Duplicate original coordinate row");
					member.raw=true;
					final String intervals=text(8);
					final int dash=intervals.indexOf('-');
					if(dash>0 && intervals.indexOf(',', 0)<0 && intervals.indexOf(';', 0)<0 && intervals.indexOf('-', dash+1)<0){
						member.start=Integer.parseInt(intervals.substring(0, dash));
						member.stop=Integer.parseInt(intervals.substring(dash+1));
					}
				}else{
					require(!member.corrected && (mode!=1 || member.raw), "Corrected member coordinates require one appropriate original witness");
					if(mode==2){require(p.termEquals("true", 23), "A retained cross-labelled member needs an admitted overlay");}
					final int offset=mode==1 ? 35 : 27;
					member.start=p.parseInt(offset);
					member.stop=p.parseInt(offset+1);
					require(member.stop-member.start+1==member.length && p.parseInt(offset+2)==member.length,
						"Corrected coordinates must span the retained member's recorded length");
					member.corrected=true;
				}
				require(member.contig==null || member.contig.equals(contig) && member.strand==strand,
					"Coordinate correction cannot change contig or strand");
				member.contig=contig;
				member.strand=strand;
				final String headerContig=RrnaResourceIO.optionalField(member.header, "contig");
				final boolean needsCorrection=RrnaResourceIO.optionalField(member.header, "recovery")!=null
					|| RrnaResourceIO.optionalField(member.header, "coordinate_model")!=null
					|| member.id.startsWith("RF01482_") || member.id.startsWith("RF01689_");
				boolean headerAgrees=true;
				if(headerContig!=null){
					headerAgrees=headerContig.equals(contig) && RrnaResourceIO.field(member.header, "strand").charAt(0)==strand
						&& Integer.parseInt(RrnaResourceIO.field(member.header, "start"))==member.start
						&& Integer.parseInt(RrnaResourceIO.field(member.header, "stop"))==member.stop;
				}
				require(mode==0 || headerAgrees, "Final corrected coordinates must agree with the retained construction header");
				member.known=member.start>0 && member.stop-member.start+1==member.length
					&& headerAgrees && (!needsCorrection || member.corrected);
			}
		}
	}

	private void compareMembers(){
		try(Output output=new Output(path("members.tsv.gz")); Output pairs=new Output(path("member_overlaps.tsv.gz"))){
			output.writer.println("member\tassembly_group\tcontig\tstart1\tstop1\tstrand\texact_qualified_loci\treciprocal80_qualified_loci\tother_overlap_loci\tunmapped_comparisons\tdisposition");
			pairs.writer.println("member\tlocus_id\tshared_bases\tmember_bases\tlocus_bases\trelationship");
			for(final Member member:members){
				int exact=0, reciprocal=0, other=0, unmapped=0;
				if(member.known){
					final String frame=coordinateFrame(member.assembly, member.contig);
					for(Locus locus=byGroup.get(member.assembly.group); locus!=null; locus=locus.nextGroup){
						if(locus.qualified && !locus.frame.equals(frame)){unmapped++;}
					}
					for(Locus locus=byContig.get(bucket(member.assembly.group, member.contig, member.strand)); locus!=null; locus=locus.next){
						if(!locus.qualified){continue;}
						if(!locus.frame.equals(frame)){
							continue;
						}
						final int shared=Math.min(member.stop, locus.stop)-Math.max(member.start, locus.start)+1;
						if(shared<=0){continue;}
						final int length=locus.stop-locus.start+1;
						final boolean same=member.start==locus.start && member.stop==locus.stop;
						final boolean compatible=reciprocal80(shared, member.length, length);
						if(same){
							exact++;
							locus.exactMembers++;
						}
						if(compatible){
							reciprocal++;
							locus.compatibleMembers++;
						}else{other++;}
						pairs.writer.print(b.clear().append(utf8(member.id)).tab().append(locus.id).tab().append(shared).tab()
							.append(member.length).tab().append(length).tab().append(same ? "EXACT" : compatible ? "RECIPROCAL80_NOT_IDENTICAL" : "OVERLAP_ONLY").nl());
					}
				}
				final String state;
				if(exact>0){
					state="EXACT_COORDINATES";
					exactMembers++;
				}else if(!member.known || unmapped>0 || uncertainGroups.contains(member.assembly.group)){
					state="UNRESOLVED_MAPPING_OR_ASSEMBLY";
					unmappedMembers++;
				}else if(reciprocal==1){
					state="UNIQUE_RECIPROCAL80_NOT_IDENTICAL";
					compatibleMembers++;
				}else if(reciprocal>1){
					state="AMBIGUOUS_RECIPROCAL80";
					ambiguousMembers++;
				}else{
					state=other>0 ? "OVERLAP_ONLY_UNRESOLVED" : "NO_QUALIFIED_COORDINATE_MATCH";
					unmatchedMembers++;
				}
				output.writer.print(b.clear().append(utf8(member.id)).tab().append(member.assembly.group).tab()
					.append(utf8(member.contig==null ? "." : member.contig)).tab().append(member.start).tab().append(member.stop).tab()
					.append(member.strand).tab().append(exact).tab().append(reciprocal).tab().append(other).tab().append(unmapped).tab().append(state).nl());
			}
		}
		assert(exactMembers+compatibleMembers+ambiguousMembers+unmatchedMembers+unmappedMembers==members.size())
			: "Coordinate match dispositions must conserve every retained member";
	}

	public static boolean reciprocal80(final int shared, final int a, final int b){
		assert(a>0 && b>0) : "Reciprocal coverage uses positive recorded genomic spans";
		return shared>0 && (long)shared*5>=4L*a && (long)shared*5>=4L*b;
	}

	private void writeLoci(){
		try(Output output=new Output(path("loci.tsv.gz"))){
			output.writer.println("locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tqualified\tprior_annotation_overlap\traw_events\texact_members\treciprocal80_members\tdisposition");
			for(final Locus locus:loci){
				if(locus.qualified){
					qualifiedLoci++;
					if(locus.prior){qualifiedPriorLoci++;}
					positiveGroups.put(locus.group, 1);
					final SourceFrame source=frames.get(locus.frame);
					require(source!=null, "Every native locus must retain its saved scan-frame census");
					source.gaLoci++;
					if(locus.exactMembers>0){lociWithExactMembers++;}
					if(locus.compatibleMembers>0){lociWithCompatibleMembers++;}
				}
				output.writer.print(b.clear().append(locus.id).tab().append(locus.group).tab().append(utf8(locus.frame)).tab()
					.append(utf8(locus.contig)).tab().append(locus.start).tab().append(locus.stop).tab().append(locus.strand).tab()
					.append(locus.qualified).tab().append(locus.prior).tab().append(locus.events).tab()
					.append(locus.exactMembers).tab().append(locus.compatibleMembers).tab()
					.append(locus.qualified ? "NATIVE_GA_COORDINATE_GROUP_OTHER_ALIASES_AND_OVERLAPS_UNRESOLVED" : "NO_NATIVE_GA_INCLUSION_RETAINED_AS_CONTEXT").nl());
			}
		}
		positiveAssemblies=positiveGroups.size();
	}

	/** Connected components of qualified, overlapping inclusive intervals in one
	 * proven coordinate frame. A chain joins transitively; adjacent bounds do not. */
	private void writeComponents(){
		final Locus[] qualified=new Locus[qualifiedLoci];
		int used=0;
		for(final Locus locus:loci){if(locus.qualified){qualified[used++]=locus;}}
		assert(used==qualified.length) : "Only native-GA coordinate groups enter the complete overlap roster";
		Arrays.sort(qualified, NcrnaDataReconciliation::compareCoordinates);
		try(Output output=new Output(path("merged_loci.tsv.gz")); Output links=new Output(path("locus_components.tsv.gz"))){
			output.writer.println("locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\texact_coordinate_groups\tnative_GA_events\tprotected_prior_witness");
			links.writer.println("exact_coordinate_group_id\tmerged_locus_id");
			int begin=0;
			while(begin<qualified.length){
				final Locus first=qualified[begin];
				int end=begin, stop=first.stop, nativeEvents=0;
				boolean prior=false;
				while(end<qualified.length){
					final Locus next=qualified[end];
					if(next.group!=first.group || !next.frame.equals(first.frame) || !next.contig.equals(first.contig)
						|| next.strand!=first.strand || next.start>stop){break;}
					stop=Math.max(stop, next.stop);
					prior|=next.prior;
					nativeEvents+=next.nativeEvents;
					links.writer.print(b.clear().append(next.id).tab().append(mergedLoci).nl());
					end++;
				}
				assert(end>begin && stop>=first.start) : "An overlap component must conserve at least its first qualified interval";
				output.writer.print(b.clear().append(mergedLoci).tab().append(first.group).tab().append(utf8(first.frame)).tab()
					.append(utf8(first.contig)).tab().append(first.start).tab().append(stop).tab().append(first.strand).tab()
					.append(end-begin).tab().append(nativeEvents).tab().append(prior).nl());
				mergedLoci++;
				mergedSpanBases+=(long)stop-first.start+1;
				maxMergedSpan=Math.max(maxMergedSpan, stop-first.start+1);
				begin=end;
			}
		}
		assert(mergedLoci<=qualifiedLoci) : "Overlap merging cannot create more loci than its qualified coordinate groups";
	}

	private static int compareCoordinates(final Locus a, final Locus b){
		int order=Integer.compare(a.group, b.group);
		if(order==0){order=a.frame.compareTo(b.frame);}
		if(order==0){order=a.contig.compareTo(b.contig);}
		if(order==0){order=Character.compare(a.strand, b.strand);}
		if(order==0){order=Integer.compare(a.start, b.start);}
		if(order==0){order=Integer.compare(a.stop, b.stop);}
		return order==0 ? Integer.compare(a.id, b.id) : order;
	}

	/** Cost inputs from saved census metadata, not a sequence-equivalence claim.
	 * Each canonical source frame contributes its bytes/bases exactly once. */
	private void writeFrames(){
		final IntHashMap seenGroups=new IntHashMap();
		try(Output output=new Output(path("source_frames.tsv.gz"))){
			output.writer.println("assembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\tcohort_input_rows\tframes_in_group\tnative_GA_exact_groups\tassembly_has_native_GA\tcoordinate_comparison_needed\tsource_available\treason");
			for(final SourceFrame frame:frameList){
				final int group=frame.assembly.group, count=frameCounts.get(group);
				final boolean pending=count>1 || uncertainGroups.contains(group), positive=positiveGroups.contains(group);
				final boolean available=new File(frame.assembly.literal).isFile();
				if(!available){unavailableFrames++;}
				frameBytes+=frame.bytes;
				frameBases+=frame.bases;
				if(pending){
					comparisonFrames++;
					comparisonBytes+=frame.bytes;
					comparisonBases+=frame.bases;
					if(positive){
						rosterComparisonFrames++;
						rosterComparisonBytes+=frame.bytes;
						rosterComparisonBases+=frame.bases;
					}
					if(!seenGroups.contains(group)){
						seenGroups.put(group, 1);
						comparisonGroups++;
						comparisonPairs+=(long)count*(count-1)/2;
						if(positive){rosterComparisonPairs+=(long)count*(count-1)/2;}
					}
				}
				output.writer.print(b.clear().append(group).tab().append(utf8(frame.assembly.canonical)).tab()
					.append(utf8(frame.assembly.literal)).tab().append(frame.bytes).tab().append(frame.bases).tab()
					.append(frame.records).tab().append(frame.inputs).tab().append(count).tab().append(frame.gaLoci).tab()
					.append(positive).tab().append(pending).tab().append(available).tab().append(count>1 ? "MULTIPLE_SOURCE_COORDINATE_FRAMES"
						: pending ? "UNRESOLVED_ASSEMBLY_AUDIT" : "SINGLE_FRAME_OTHER_ASSEMBLY_ALIASES_UNPROVED").nl());
			}
		}
	}

	private void writeSummary(){
		try(Output output=new Output(path("summary.tsv"))){
			output.writer.println("metric\tcount\tqualification");
			metric(output, "accepted_alias_audit_assembly_groups", auditedGroups.size(), "UNCHANGED_HISTORICAL_AUDIT_CHECK");
			metric(output, "scanned_assembly_groups", scannedGroups.size(), "ALL_DECLARED_COHORT_CENSUSES_BOUND_BEFORE_LIVE_INVENTORY");
			metric(output, "current_live_source_paths", livePaths, "HISTORICAL_51079_ROSTER_NOT_CONSULTED");
			metric(output, "current_live_assembly_groups", liveAssemblyGroups, "ACCEPTED_ALIASES_PLUS_CANONICAL_PATHS_OTHER_ALIASES_UNRESOLVED");
			metric(output, "unscanned_live_paths", missingPaths, "NO_ACCEPTED_SCAN_BINDING_NO_SEARCH_PERFORMED");
			metric(output, "unscanned_live_assembly_groups", unscannedGroups, "OTHER_ASSEMBLY_ALIASES_UNRESOLVED");
			metric(output, "original_scored_RF00174_reference_loci", originalScored, "EXISTING_REFERENCE_NOT_CONFIDENT_EVENT_COUNT");
			metric(output, "import_RF00174_confident_events", importConfident, "EVENTS_NOT_DISTINCT_LOCI");
			metric(output, "raw_RF00174_events", rawEvents, "ALL_SAVED_TIERS_RETAINED");
			metric(output, "native_GA_RF00174_events", nativeGaEvents, "INCLUDES_PARTIAL_AND_PRIOR_COVERED_HITS");
			metric(output, "native_GA_events_with_prior_annotation", nativeGaPriorEvents, "WITNESS_LINKED_OVERLAP");
			metric(output, "native_GA_events_without_prior_annotation", nativeGaEvents-nativeGaPriorEvents, "NO_LINKED_PROTECTED_PRIOR");
			metric(output, "marginal_RF00174_events", marginalEvents, "SAVED_TIER_NOT_ADDITIVE_TO_GA_COUNTS");
			metric(output, "marginal_native_GA_events", marginalGaEvents, "INTERSECTION_OF_MARGINAL_TIER_AND_NATIVE_GA");
			metric(output, "uncertain_RF00174_events", uncertainEvents, "SAVED_UNCERTAIN_FLOOR_TIER_NOT_ADDITIVE_TO_GA_COUNTS");
			metric(output, "uncertain_native_GA_events", uncertainGaEvents, "INTERSECTION_OF_UNCERTAIN_TIER_AND_NATIVE_GA");
			metric(output, "RF00174_prior_annotations_without_native_GA_hit", priorOnlyAnnotations, "OUTSIDE_CM_ROSTER");
			metric(output, "saved_scored_RF00174_annotation_rows", qualifiedAnnotations, "INCLUDES_PROTECTED_PRIORS_NOT_A_CM_HIT_COUNT");
			metric(output, "qualified_exact_coordinate_groups", qualifiedLoci, "NONIDENTICAL_OVERLAPS_NOT_MERGED");
			metric(output, "qualified_same_frame_overlap_loci", mergedLoci, "INCLUSIVE_CONNECTED_COMPONENTS_SAME_CONTIG_STRAND_FRAME_ONLY");
			metric(output, "coordinate_groups_collapsed_by_overlap", qualifiedLoci-mergedLoci, "EXACT_GROUPS_MINUS_COMPONENTS_NOT_A_CROSS_FRAME_DEDUPLICATION");
			metric(output, "merged_locus_span_bases", mergedSpanBases, "INCLUSIVE_COORDINATE_SPANS_NOT_EXTRACTED_SEQUENCE");
			metric(output, "largest_merged_locus_span", maxMergedSpan, "INCLUSIVE_COORDINATE_SPAN_NOT_A_COMPLETE_GENE_CLAIM");
			metric(output, "qualified_coordinate_groups_with_prior", qualifiedPriorLoci, "SAVED_WITNESS_OVERLAP");
			metric(output, "assemblies_with_qualified_coordinates", positiveAssemblies, "ACCEPTED_ASSEMBLY_GROUPS");
			metric(output, "kept_library_members", members.size(), "EXACT_ASSIGNMENT_LEDGER_JOIN_NOT_SEQUENCE_REVALIDATION");
			metric(output, "exact_member_matches", exactMembers, "CORRECTED_COORDINATES");
			metric(output, "unique_reciprocal80_member_matches", compatibleMembers, "NOT_IDENTICAL_BOUNDARIES");
			metric(output, "ambiguous_member_matches", ambiguousMembers, "MULTIPLE_RECIPROCAL80_LOCI");
			metric(output, "members_without_qualified_match", unmatchedMembers, "OVERLAP_ONLY_SEPARATE_IN_MEMBER_TABLE");
			metric(output, "unmapped_members", unmappedMembers, "NOT_ASSUMED_ABSENT");
			metric(output, "qualified_loci_with_exact_member", lociWithExactMembers, "CORRECTED_COORDINATES");
			metric(output, "qualified_loci_without_exact_member", qualifiedLoci-lociWithExactMembers, "NO_PROVEN_EXACT_COORDINATE_MATCH_UNRESOLVED_MAPPINGS_REMAIN");
			metric(output, "qualified_loci_with_reciprocal80_member", lociWithCompatibleMembers, "NOT_IDENTICAL_BOUNDARIES");
			metric(output, "distinct_scanned_source_frames", frameList.size(), "CANONICAL_PATHS_NOT_SEQUENCE_DEDUPLICATION");
			metric(output, "unavailable_scanned_source_frames", unavailableFrames, "NAMED_SOURCE_PATH_NOT_CURRENTLY_A_FILE_NO_PAYLOAD_READ");
			metric(output, "saved_source_frame_file_bytes", frameBytes, "SAVED_CENSUS_COMPRESSED_OR_PHYSICAL_FILE_BYTES_ONCE_PER_FRAME");
			metric(output, "saved_source_frame_nucleotide_bases", frameBases, "SAVED_CENSUS_BASES_ONCE_PER_FRAME");
			metric(output, "coordinate_comparison_assembly_groups", comparisonGroups, "MULTIPLE_FRAMES_OR_UNRESOLVED_AUDIT");
			metric(output, "coordinate_comparison_source_frames", comparisonFrames, "ALL_FRAMES_IN_THOSE_GROUPS");
			metric(output, "coordinate_comparison_frame_pairs", comparisonPairs, "UNORDERED_WITHIN_GROUP_PAIRS_NOT_A_SEQUENCE_RUN");
			metric(output, "coordinate_comparison_file_bytes", comparisonBytes, "SAVED_CENSUS_ONCE_PER_FRAME");
			metric(output, "coordinate_comparison_nucleotide_bases", comparisonBases, "SAVED_CENSUS_ONCE_PER_FRAME");
			metric(output, "roster_coordinate_comparison_source_frames", rosterComparisonFrames, "ALL_GROUP_FRAMES_WHERE_AT_LEAST_ONE_HAS_NATIVE_GA");
			metric(output, "roster_coordinate_comparison_frame_pairs", rosterComparisonPairs, "UNORDERED_WITHIN_POSITIVE_GROUP_PAIRS");
			metric(output, "roster_coordinate_comparison_file_bytes", rosterComparisonBytes, "SAVED_CENSUS_ONCE_PER_FRAME");
			metric(output, "roster_coordinate_comparison_nucleotide_bases", rosterComparisonBases, "SAVED_CENSUS_ONCE_PER_FRAME");
		}
	}

	private void metric(final Output output, final String name, final long count, final String qualification){
		assert(count>=0) : "Summary counts cannot silently encode missing values as negatives";
		output.writer.print(b.clear().append(name).tab().append(count).tab().append(qualification).nl());
	}
	private String path(final String name){return out+File.separator+name;}
	private String text(final int column){
		p.setBounds(column);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 must not change a metadata identity");
		return value;
	}
	private String decoded(final int column){
		p.setBounds(column);
		return decoder.decode(p.line(), p.a(), p.b());
	}
	/** Cache retained output identities once; ByteBuilder.append(String) is ASCII. */
	private byte[] utf8(final String value){
		byte[] bytes=encoded.get(value);
		if(bytes==null){
			bytes=value.getBytes(StandardCharsets.UTF_8);
			require(new String(bytes, StandardCharsets.UTF_8).equals(value), "Invalid Unicode must not change an output identity");
			encoded.put(value, bytes);
		}
		return bytes;
	}
	private char strand(final int column){
		require(p.termEquals("+", column) || p.termEquals("-", column), "Coordinates require an explicit strand");
		return p.termEquals("+", column) ? '+' : '-';
	}
	private static String bucket(final int group, final String contig, final char strand){return group+"\t"+contig+"\t"+strand;}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Cohort{
		Cohort(final String label_, final String census_, final String reference_, final String events_,
				final String witnesses_, final int expectedInputs_, final int expectedScored_){
			label=label_;
			census=census_;
			reference=reference_;
			events=events_;
			witnesses=witnesses_;
			expectedInputs=expectedInputs_;
			expectedScored=expectedScored_;
		}
		final String label, census, reference, events, witnesses;
		final int expectedInputs, expectedScored;
		final ObjectMap<String, Assembly> inputs=new ObjectMap<String, Assembly>(String.class, Assembly.class);
		final ArrayList<String> inputIds=new ArrayList<String>();
	}

	private static final class SourceFrame{
		SourceFrame(final Assembly assembly_, final long bytes_, final long bases_, final long records_){
			assembly=assembly_;
			bytes=bytes_;
			bases=bases_;
			records=records_;
		}
		final Assembly assembly;
		final long bytes, bases, records;
		int inputs, gaLoci;
	}

	public static final class Assembly{
		public Assembly(final int group_, final String literal_, final String canonical_, final boolean uncertain_){
			group=group_;
			literal=literal_;
			canonical=canonical_;
			uncertain=uncertain_;
		}
		final int group;
		final String literal, canonical;
		final boolean uncertain;
	}
	private static final class Annotation{
		Annotation(final Row row_, final Assembly assembly_, final boolean savedSupport_, final boolean prior_){
			row=row_;
			assembly=assembly_;
			savedSupport=savedSupport_;
			prior=prior_;
		}
		final Row row;
		final Assembly assembly;
		final boolean savedSupport, prior;
		int support, confident, nativeGa;
	}
	private static final class Event{
		Event(final byte[] row_, final Assembly assembly_, final boolean confident_, final boolean qualified_, final String input_, final String nativeTable_, final long nativeLine_){
			row=row_;
			assembly=assembly_;
			confident=confident_;
			qualified=qualified_;
			input=input_;
			nativeTable=nativeTable_;
			nativeLine=nativeLine_;
		}
		final byte[] row;
		final Assembly assembly;
		final boolean confident, qualified;
		boolean prior;
		final String input, nativeTable;
		final long nativeLine;
		int annotationLinks;
	}
	private static final class Locus{
		Locus(final int id_, final int group_, final String frame_, final String contig_, final int start_, final int stop_, final char strand_){
			id=id_;
			group=group_;
			frame=frame_;
			contig=contig_;
			start=start_;
			stop=stop_;
			strand=strand_;
		}
		final int id, group, start, stop;
		final String frame, contig;
		final char strand;
		int events, nativeEvents, exactMembers, compatibleMembers;
		boolean qualified, prior;
		Locus next, nextGroup;
	}
	private static final class Member{
		Member(final String id_, final String header_, final Assembly assembly_){
			id=id_;
			header=header_;
			assembly=assembly_;
		}
		final String id, header;
		final Assembly assembly;
		String contig;
		int start, stop, length;
		char strand='.';
		boolean bound, raw, corrected, known;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, out, groups, mapping, genomes, assignments, ledger, proxy, recovered;
	private String[] roots, overlays;
	private final PrintStream outstream;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ObjectMap<String, byte[]> encoded=new ObjectMap<String, byte[]>(String.class, byte[].class);
	private final ObjectMap<String, Assembly> assemblies=new ObjectMap<String, Assembly>(String.class, Assembly.class);
	private final ObjectMap<String, String> canonicalPaths=new ObjectMap<String, String>(String.class, String.class);
	private final ObjectMap<String, Locus> coordinates=new ObjectMap<String, Locus>(String.class, Locus.class);
	private final ObjectMap<String, Locus> byContig=new ObjectMap<String, Locus>(String.class, Locus.class);
	private final IntObjectMap<Locus> byGroup=new IntObjectMap<Locus>(Locus.class);
	private final ObjectMap<String, Member> memberIds=new ObjectMap<String, Member>(String.class, Member.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	private final ArrayList<Member> members=new ArrayList<Member>();
	private final ArrayList<Cohort> cohorts=new ArrayList<Cohort>();
	private final ObjectMap<String, SourceFrame> frames=new ObjectMap<String, SourceFrame>(String.class, SourceFrame.class);
	private final ArrayList<SourceFrame> frameList=new ArrayList<SourceFrame>();
	private final IntHashMap auditedGroups=new IntHashMap(), scannedGroups=new IntHashMap(), uncertainGroups=new IntHashMap();
	private final IntHashMap frameCounts=new IntHashMap(), positiveGroups=new IntHashMap();
	private int nextGroup, livePaths, liveAssemblyGroups, unscannedGroups, missingPaths;
	private int originalScored, importConfident, qualifiedAnnotations, qualifiedLoci, positiveAssemblies;
	private int exactMembers, compatibleMembers, ambiguousMembers, unmatchedMembers, unmappedMembers;
	private int lociWithExactMembers, lociWithCompatibleMembers;
	private int rawEvents, nativeGaEvents, nativeGaPriorEvents, priorOnlyAnnotations, qualifiedPriorLoci;
	private int marginalEvents, marginalGaEvents, uncertainEvents, uncertainGaEvents;
	private int mergedLoci, maxMergedSpan, comparisonGroups, comparisonFrames, rosterComparisonFrames, unavailableFrames;
	private long mergedSpanBases;
	private long frameBytes, frameBases, comparisonPairs, comparisonBytes, comparisonBases;
	private long rosterComparisonPairs, rosterComparisonBytes, rosterComparisonBases;
	private int expectedCohorts=2;
	private int expectedGroups=40723, expectedReference=25916, expectedConfident=57723, expectedKept=85726, expectedAssigned=85741;
	private int expectedGenomes=39979, expectedModels=508;
	private static final String FAMILY="RF00174";
	private static final String CENSUS_HEADER="input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records";
	private static final String COHORT_HEADER="cohort\tcensus\treference\tevents\twitnesses\texpected_inputs\texpected_scored_RF00174";
	private static final String ASSIGNMENT_HEADER="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus";
	private static final String LEDGER_HEADER="locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header";
	private static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation";
	private static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
}
