package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaHistoricalEvidence.Feature;
import prok.NcrnaHistoricalEvidence.Source;
import prok.NcrnaReferenceEvidence.Evidence;
import prok.NcrnaReferenceEvidence.Hit;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.Outputs;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Joins completed whole-genome searches and historical GFFs for reference review.
 * Evidence/report modes are descriptive. The explicit tiers mode applies the
 * declared construction policy; tierreport describes its saved results without
 * changing decisions. No mode infers biological absence from a missing annotation.
 * Derived from template.A_SampleByteFileMT: native parsing/output, static workers
 * and once-per-worker accumulation. Each worker owns one input's observations,
 * GFF graph and hit index at a time; input ownership is striped without a queue.
 * Original source files remain authoritative. Text fields in output TSVs escape
 * percent, tab, CR and LF as %HH and retain UTF-8; native raw lines are cited by
 * source path and one-based line, not rewritten into another raw-data copy.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaReference implements Accumulator<BuildNcrnaReference.ProcessThread>{

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BuildNcrnaReference.class, false);
		try{new BuildNcrnaReference(pp.args, pp.outstream).process(timer);}
		finally{Shared.closeStream(pp.outstream);}
	}

	private BuildNcrnaReference(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("census")){
				census=b;
			}else if(a.equals("memberships")){
				memberships=b;
			}else if(a.equals("scope")){
				scope=b;
			}else if(a.equals("assets")){
				assets=b;
			}else if(a.equals("queries")){
				queries=b;
			}else if(a.equals("attributes")){
				attributePrefix=b;
			}else if(a.equals("attributesevidence")){
				attributeEvidence=b;
			}else if(a.equals("marginalpercent")){
				marginalPercent=b;
			}else if(a.equals("collapsepercent")){
				collapsePercent=b;
			}else if(a.equals("samplemodel")){
				sampleModel=b;
			}else if(a.equals("examplespergroup")){
				examplesPerGroup=Parse.parseIntKMG(b);
				examplesSpecified=true;
			}else if(a.equals("priors")){
				priorManifest=b;
			}else if(a.equals("anchors")){
				anchorManifest=b;
			}else if(a.equals("sources")){
				sourcePrefix=b;
			}else if(a.equals("compare")){
				comparison=b;
			}else if(a.equals("searchspacemb")){
				batchSearchSpace=Double.parseDouble(b);
			}else if(a.equals("nativepaths")){
				nativePaths=b;
			}else if(a.equals("reference")){
				panelReference=b;
			}else if(a.equals("observed")){
				panelContigs=b;
			}else if(a.equals("sourcepanel")){
				panelSource=b;
			}else if(a.equals("selectioninventory")){
				selectionInventory=b;
			}else if(a.equals("selectionpolicy")){
				selectionPolicy=b;
			}else if(a.equals("selectiontaxonomy")){
				selectionTaxonomy=b;
			}else if(a.equals("selectionfamily")){
				selectionFamily=b;
			}else if(a.equals("inventorymembers")){
				inventoryMembers=b;
			}else if(a.equals("witnesses")){
				inventoryWitnesses=b;
			}else if(a.equals("events")){
				inventoryEvents=b;
			}else if(a.equals("maxinputs")){
				maxInputs=Parse.parseIntKMG(b);
				samplingSpecified=true;
			}else if(a.equals("mode")){
				require("evidence".equals(b) || "batchevidence".equals(b) || "profile".equals(b) || "scores".equals(b) || "retention".equals(b) || "tiers".equals(b) || "tierreport".equals(b)
					|| "geometry".equals(b) || "attributes".equals(b) || "strata".equals(b) || "prior".equals(b)
					|| "cost".equals(b) || "panels".equals(b) || "panel".equals(b) || "inventory".equals(b) || "selectpanel".equals(b) || "anchors".equals(b) || "sourceanchors".equals(b),
					"Unknown reference mode: "+b);
				batch=b.equals("batchevidence");
				profile=b.equals("profile");
				scores=b.equals("scores");
				retention=b.equals("retention");
				tiers=b.equals("tiers");
				tierReport=b.equals("tierreport");
				geometry=b.equals("geometry");
				attributes=b.equals("attributes");
				strata=b.equals("strata");
				prior=b.equals("prior");
				cost=b.equals("cost");
				panels=b.equals("panels");
				panel=b.equals("panel");
				inventory=b.equals("inventory");
				selectPanel=b.equals("selectpanel");
				anchors=b.equals("anchors");
				sourceAnchors=b.equals("sourceanchors");
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		input=parser.in1;
		prefix=parser.out1;
		require(selectPanel || selectionInventory==null && selectionPolicy==null && selectionTaxonomy==null && selectionFamily==null,
			"Selection-only options require mode=selectpanel; they must not silently affect another mode");
		if(selectPanel){
			require(census!=null && memberships!=null && panelReference!=null && inventoryMembers!=null
				&& selectionInventory!=null && selectionPolicy!=null && selectionTaxonomy!=null && selectionFamily!=null && prefix!=null,
				"selectpanel requires census, memberships, reference, inventorymembers, selectioninventory, selectionpolicy, selectiontaxonomy, selectionfamily and out");
			require(input==null && scope==null && assets==null && queries==null && comparison==null
				&& attributePrefix==null && attributeEvidence==null && priorManifest==null && panelSource==null && panelContigs==null
				&& marginalPercent==null && collapsePercent==null && anchorManifest==null && sourcePrefix==null && nativePaths==null
				&& inventoryWitnesses==null && inventoryEvents==null && !samplingSpecified && !examplesSpecified && sampleModel==null
				&& batchSearchSpace==-1 && !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
				&& parser.maxReads<0 && parser.samplerate==1,
				"Panel selection reads completed metadata only; no sequence reads, evidence reconstruction, sampling or overwrite");
			final ArrayList<String> files=new ArrayList<String>(Arrays.asList(census, memberships, panelReference, inventoryMembers,
				selectionPolicy, selectionTaxonomy, selectionInventory+".inputs.tsv", selectionInventory+".members.tsv", selectionInventory+".provenance.tsv"));
			require(Tools.testInputFiles(false, true, files.toArray(new String[files.size()])), "Selection metadata must be readable");
			for(final String suffix:NcrnaPanelSelector.SUFFIXES){
				require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Selection outputs must be fresh");
				files.add(prefix+suffix);
			}
			require(Tools.testForDuplicateFiles(true, files.toArray(new String[files.size()])), "Selection inputs and outputs must be distinct");
			return;
		}
		require(inventory || inventoryMembers==null && inventoryWitnesses==null && inventoryEvents==null,
			"Member inventory inputs require mode=inventory; they cannot silently affect another reference operation");
		if(inventory){
			require(census!=null && memberships!=null && panelReference!=null && inventoryMembers!=null
				&& inventoryWitnesses!=null && inventoryEvents!=null && prefix!=null,
				"Inventory requires census, memberships, reference, inventorymembers, witnesses, events and fresh out prefix");
			require(input==null && scope==null && assets==null && queries==null && comparison==null
				&& attributePrefix==null && attributeEvidence==null && priorManifest==null
				&& panelSource==null && panelContigs==null && marginalPercent==null && collapsePercent==null
				&& anchorManifest==null && sourcePrefix==null && nativePaths==null && !samplingSpecified
				&& !examplesSpecified && sampleModel==null && batchSearchSpace==-1 && !parser.overwrite && !parser.append
				&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
				"Inventory reads retained metadata only; no construction, sampling, reinterpretation or overwrite");
			final String[] files={census, memberships, panelReference, inventoryMembers, inventoryWitnesses, inventoryEvents};
			require(Tools.testInputFiles(false, true, files), "Inventory metadata inputs must be readable");
			final ArrayList<String> allPaths=new ArrayList<String>(java.util.Arrays.asList(files));
			for(final String suffix:NcrnaReferenceInventory.SUFFIXES){
				require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Inventory outputs must be fresh");
				allPaths.add(prefix+suffix);
			}
			require(Tools.testForDuplicateFiles(true, allPaths.toArray(new String[allPaths.size()])), "Inventory inputs and outputs must be distinct");
			return;
		}
		require(batch ? Double.isFinite(batchSearchSpace) && batchSearchSpace>0 && !samplingSpecified : batchSearchSpace==-1,
			"batchevidence requires explicit searchspacemb= matching native -Z and complete batch membership; other modes do not accept it");
		require(nativePaths==null || batch, "nativepaths= applies only to batchevidence; it changes read locations, not native provenance");
		require(anchors==(anchorManifest!=null), "anchors=manifest.tsv is required only for mode=anchors");
		require(sourceAnchors==(sourcePrefix!=null), "sources=saved_strata_prefix is required only for mode=sourceanchors");
		require(marginalPercent==null || tiers, "Only tiers mode selects marginalpercent");
		require(collapsePercent==null || tiers, "Only tiers mode selects collapsepercent");
		require(attributeEvidence==null || strata || tiers, "Only strata or tiers mode reuses an earlier attribute capture");
		require(tierReport || sampleModel==null && !examplesSpecified, "Example sampling options apply only to tierreport mode");
		require(panel || panelReference==null && panelContigs==null && panelSource==null,
			"reference, observed and sourcepanel apply only to selected-panel preparation");
		if(panel){
			require(input!=null && census!=null && memberships!=null && prefix!=null
				&& panelReference!=null && panelContigs!=null && panelSource!=null && !panelSource.isEmpty(),
				"Panel preparation requires in=flat_caller_config census= memberships= reference= observed= sourcepanel= out=prefix");
			require(scope==null && assets==null && queries==null && comparison==null && attributePrefix==null
				&& priorManifest==null && !samplingSpecified && !parser.overwrite && !parser.append
				&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
				"Panel preparation selects saved metadata only; no evidence construction or overwrite");
			require(Tools.testInputFiles(false, true, input, census, memberships, panelReference, panelContigs),
				"Selected-panel metadata inputs must be readable");
			for(final String suffix:NcrnaEvaluationPanel.SUFFIXES){
				require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Panel metadata outputs must be fresh");
			}
			return;
		}
		if(tierReport){
			require(input!=null && prefix!=null && new File(input).isDirectory(), "tierreport requires in=completed_reference_directory out=fresh_report_prefix");
			require(!parser.overwrite && !parser.append && !samplingSpecified && census==null && memberships==null
				&& scope==null && assets==null && queries==null && comparison==null && attributePrefix==null && priorManifest==null,
				"tierreport reads the completed reference outputs only and never changes decisions");
			final ArrayList<String> paths=new ArrayList<String>();
			for(final String name:NcrnaReferenceTierReport.INPUTS){paths.add(new File(input, name).getPath());}
			require(Tools.testInputFiles(false, true, paths.toArray(new String[paths.size()])), "A saved tier report input cannot be read");
			for(final String suffix:NcrnaReferenceTierReport.SUFFIXES){
				require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Tier report output must be fresh");
				paths.add(prefix+suffix);
			}
			require(Tools.testForDuplicateFiles(true, paths.toArray(new String[paths.size()])), "Tier report inputs and outputs must be distinct");
			return;
		}
		if(panels){
			require(census!=null && memberships!=null && prefix!=null,
				"Panel counts require census= memberships= out=new_prefix");
			require(input==null && scope==null && assets==null && queries==null && comparison==null
				&& attributePrefix==null && priorManifest==null && !samplingSpecified && !parser.overwrite && !parser.append,
				"Panel counts use the complete saved catalog, no other reference inputs, and fresh output");
			require(Tools.testOutputFiles(false, false, false, prefix+NcrnaReferencePanels.SUFFIX), "Panel output must be fresh");
			return;
		}
		if(cost){
			require(input!=null && prefix!=null && census!=null && new File(input).isDirectory(),
				"Cost mode requires in=completed_pilot_directory census=full_population_summary out=new_prefix");
			require(!parser.overwrite && !parser.append && !samplingSpecified && memberships==null
				&& scope==null && assets==null && queries==null && attributePrefix==null && priorManifest==null,
				"Cost mode accepts only pilot metadata and a full population summary; use fresh outputs");
			for(final String suffix:NcrnaReferenceCost.SUFFIXES){
				require(Tools.testOutputFiles(false, false, false, prefix+suffix), "Cost output must be fresh: "+prefix+suffix);
			}
			return;
		}
		require(retention==(comparison!=null), "compare=previous_saved_prefix is required for retention; otherwise compare= is supported only in cost mode");
		require(input!=null && prefix!=null && census!=null && memberships!=null && scope!=null && assets!=null && queries!=null,
			"Require in= census= memberships= scope= assets= queries= out=new_prefix; in is a collector table or saved profile prefix");
		require(!parser.overwrite && !parser.append && maxInputs>0, "Use a fresh output prefix, ow=f, append=f and positive maxinputs");
		require(!(profile || scores || retention || tiers || geometry || attributes || strata || prior || anchors || sourceAnchors) || !samplingSpecified, "Saved-table modes summarize the saved population and do not accept maxinputs");
		require((strata || tiers)==(attributePrefix!=null), "attributes=capture_prefix is required for strata or tiers");
		require((prior || tiers)==(priorManifest!=null), "priors= is required for prior import or tier construction");
		require(tiers==(marginalPercent!=null), "Tier construction requires an explicit marginalpercent=75..100; other modes do not select a floor");
		require(attributeEvidence==null || strata || tiers, "attributesevidence= applies only to strata or tier construction");
		final ArrayList<String> all=new ArrayList<String>();
		all.addAll(Arrays.asList(census, memberships, scope, assets, queries));
		if(nativePaths!=null){all.add(nativePaths);}
		if(profile || scores || retention || tiers || geometry || attributes || strata || prior || anchors || sourceAnchors){
			for(final String suffix:NcrnaReferenceProfiles.INPUT_SUFFIXES){all.add(input+suffix);}
			if(attributeEvidence!=null){
				all.add(attributeEvidence+".summary.tsv");
				all.add(attributeEvidence+".sources.tsv.gz");
			}
			if(retention){
				for(final String suffix:NcrnaReferenceProfiles.INPUT_SUFFIXES){all.add(comparison+suffix);}
				all.add(comparison+".hits.tsv.gz");
				all.add(comparison+".overlaps.tsv.gz");
			}
			if(prior || tiers){
				all.add(input+".contigs.tsv.gz");
				all.add(priorManifest);
			}
			if(tiers){all.add(input+".batch_sources.tsv");}
			if(anchors){
				all.add(anchorManifest);
				all.add(input+".contigs.tsv.gz");
			}
			if(sourceAnchors){
				all.add(sourcePrefix+".strata_summary.tsv");
				all.add(sourcePrefix+".source_links.tsv.gz");
			}
			if(scores || retention || tiers || geometry || strata || anchors || sourceAnchors){
				all.add(input+".hits.tsv.gz");
				all.add(input+".overlaps.tsv.gz");
			}
			if(attributes || strata || tiers){
				all.add(input+".historical.tsv.gz");
				all.add(input+".relations.tsv.gz");
			}
			if(strata || tiers){
				all.add(attributePrefix+".attribute_summary.tsv");
				all.add(attributePrefix+".attributes.tsv.gz");
			}
		}else{all.add(input);}
		require(Tools.testInputFiles(false, true, all.toArray(new String[all.size()])), "An evidence input cannot be read");
		if(tiers){
			require(!new File(prefix).exists(), "Tier construction requires a fresh output directory");
			require(Tools.testForDuplicateFiles(true, all.toArray(new String[all.size()])), "Tier input paths must be distinct");
			return;
		}
		final String[] suffixes=retention ? NcrnaReferenceRetention.SUFFIXES : scores ? NcrnaReferenceScores.SUFFIXES : sourceAnchors ? NcrnaReferenceSourceAnchors.SUFFIXES : anchors ? NcrnaReferenceAnchors.SUFFIXES : prior ? NcrnaReferencePrior.OUTPUT_SUFFIXES : strata ? NcrnaReferenceStrata.OUTPUT_SUFFIXES : attributes ? NcrnaReferenceAttributes.OUTPUT_SUFFIXES
			: geometry ? NcrnaReferenceGeometry.OUTPUT_SUFFIXES : profile ? NcrnaReferenceProfiles.OUTPUT_SUFFIXES : batch ? BATCH_SUFFIXES : SUFFIXES;
		for(final String suffix:suffixes){
			final String path=prefix+suffix;
			require(Tools.testOutputFiles(false, false, false, path), "Output is not fresh/writable: "+path);
			all.add(path);
		}
		require(Tools.testForDuplicateFiles(true, all.toArray(new String[all.size()])), "Input and output paths must be distinct");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input binding        ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates the entire collector/census join even for a bounded prefix diagnostic. */
	private void load(){
		catalog=NcrnaEvidenceCatalog.load(census, memberships);
		roster=NcrnaReferenceEvidence.loadModels(scope, assets, queries);
		for(int i=0; i<roster.queried(); i++){modelIndex.put(roster.query(i).accession, i);}
		if(batch){
			loadBatches();
			totals=new Counts(roster);
			return;
		}
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		final ObjectSet<String> directories=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, COLLECTOR_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==8, "Expected the eight-field native collector schema");
				final String id=p.parseString(0), fasta=p.parseString(1), directory=p.parseString(2);
				final Genome genome=catalog.get(id);
				require(genome!=null && genome.fasta.equals(fasta) && seen.add(id), "Collector input must bind once to its exact census FASTA: "+id);
				require(new File(directory).isAbsolute() && directories.add(directory), "Search result directories must be absolute and unique");
				final long hits=p.parseLong(4);
				require(p.parseInt(3)==roster.queried() && hits>=0, "Collector model/hit counts must agree with the searched roster");
				tasks.add(new Task(genome, directory, hits));
			}
		}
		require(seen.size()==catalog.size(), "Collector must cover every census entry, including zero-hit inputs");
		limit=Math.min(maxInputs, tasks.size());
		totals=new Counts(roster);
	}

	/** Binds the complete census once across prepared batches, including zero-hit inputs. */
	private void loadBatches(){
		final ObjectMap<String, NativeLocation> locations=loadNativePaths();
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), inputs=new ObjectSet<String>(String.class);
		final ObjectSet<String> prepared=new ObjectSet<String>(String.class), nativeDirs=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, BATCH_COLLECTOR_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==5 && p.parseInt(3)==roster.queried(), "Batch collector must declare every queried model");
				final String id=p.parseString(0), source=utf8(p, 1), results=utf8(p, 2);
				require(ids.add(id) && prepared.add(source) && nativeDirs.add(results), "Repeated batch ID, preparation or native result directory");
				String physical=results;
				if(locations!=null){
					final NativeLocation location=locations.get(id);
					require(location!=null, "Missing physical-native mapping for batch "+id);
					require(location.original.equals(results), "Original native directory differs from the collector for batch "+id);
					physical=location.physical;
				}
				final NcrnaReferenceBatch.Task task=NcrnaReferenceBatch.task(id, source, results, physical, p.parseLong(4), catalog);
				for(final Genome genome:task.genomes){require(inputs.add(genome.id), "A census input belongs to multiple batches: "+genome.id);}
				batchTasks.add(task);
			}
		}
		require(inputs.size()==catalog.size() && !batchTasks.isEmpty(), "Batch membership must cover the entire supplied census exactly once");
		require(locations==null || locations.size()==batchTasks.size(), "Physical-native mapping must cover exactly the imported batches, with no extras");
		limit=batchTasks.size();
	}

	/** Explicit complete relocation map; original directories may no longer exist.
	 * This declares read locations, not proof of archive identity or recovery. */
	private ObjectMap<String, NativeLocation> loadNativePaths(){
		if(nativePaths==null){return null;}
		final ObjectMap<String, NativeLocation> locations=new ObjectMap<String, NativeLocation>(String.class, NativeLocation.class);
		final ObjectSet<String> physicalDirs=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(nativePaths)){
			RrnaResourceIO.header(in.reader, NATIVE_PATH_HEADER, nativePaths);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==3, "Native relocation needs batch_id, original_native_directory and physical_native_directory");
				final String id=p.parseString(0), original=utf8(p, 1), physical=utf8(p, 2);
				require(id.matches("0|[1-9][0-9]*") && new File(original).isAbsolute()
					&& new File(physical).isAbsolute() && new File(physical).isDirectory(),
					"Native relocation requires a canonical batch ID, absolute paths and an existing physical directory");
				require(!original.equals(physical), "An explicit native relocation must name a different physical directory");
				require(!locations.contains(id) && physicalDirs.add(physical), "Duplicate native relocation batch or physical directory");
				locations.put(id, new NativeLocation(original, physical));
			}
		}
		return locations;
	}

	/** Original strings remain authoritative in native footers and emitted witnesses. */
	private static final class NativeLocation{
		NativeLocation(final String original_, final String physical_){
			original=original_;
			physical=physical_;
		}
		final String original, physical;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Thread lifecycle      ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		if(selectPanel){
			NcrnaPanelSelector.process(census, memberships, panelReference, selectionInventory,
				inventoryMembers, selectionPolicy, selectionTaxonomy, selectionFamily, prefix);
			timer.stop();
			outstream.println("Original-input panel selected from metadata; no sequence collapse, calling or adjudication. Time: "+timer);
			return;
		}
		if(inventory){
			NcrnaReferenceInventory.process(census, memberships, panelReference,
				inventoryMembers, inventoryWitnesses, inventoryEvents, prefix);
			timer.stop();
			outstream.println("Original-input reference inventory completed; no panel selection, FASTA reads or adjudication. Time: "+timer);
			return;
		}
		if(panel){
			NcrnaEvaluationPanel.process(NcrnaEvidenceCatalog.load(census, memberships), input,
				panelSource, panelReference, panelContigs, prefix);
			timer.stop();
			outstream.println("Selected panel metadata completed; reference dispositions unchanged. Registry authorization is external. Time: "+timer);
			return;
		}
		if(tierReport){
			NcrnaReferenceTierReport.process(input, prefix, sampleModel, examplesPerGroup);
			timer.stop();
			outstream.println("Saved reference reasons and prior measurement status reported; no decisions changed. Time: "+timer);
			return;
		}
		if(panels){
			NcrnaReferencePanels.process(NcrnaEvidenceCatalog.load(census, memberships), prefix);
			timer.stop();
			outstream.println("Saved panel input counts completed; no biological-reference decisions. Time: "+timer);
			return;
		}
		if(cost){
			NcrnaReferenceCost.process(input, census, comparison, prefix);
			timer.stop();
			outstream.println("Completed native timing report; projections are not full-scan measurements. Time: "+timer);
			return;
		}
		if(profile || scores || retention || tiers || geometry || attributes || strata || prior || anchors || sourceAnchors){
			final NcrnaEvidenceCatalog savedCatalog=NcrnaEvidenceCatalog.load(census, memberships);
			final Roster savedRoster=NcrnaReferenceEvidence.loadModels(scope, assets, queries);
			if(tiers){
				NcrnaReferenceTiers.process(input, attributePrefix, attributeEvidence, priorManifest,
					prefix, marginalPercent, collapsePercent==null ? "80" : collapsePercent,
					census, memberships, scope, assets, queries, savedCatalog, savedRoster);
			}else if(retention){
				NcrnaReferenceRetention.process(input, comparison, prefix, savedCatalog, savedRoster);
			}else if(scores){
				NcrnaReferenceScores.process(input, prefix, savedCatalog, savedRoster);
			}else if(sourceAnchors){
				NcrnaReferenceSourceAnchors.process(input, sourcePrefix, prefix, savedCatalog, savedRoster);
			}else if(anchors){
				NcrnaReferenceAnchors.process(input, anchorManifest, prefix, savedCatalog, savedRoster);
			}else if(prior){
				NcrnaReferencePrior.process(input, priorManifest, prefix, savedCatalog, savedRoster);
			}else if(strata){
				NcrnaReferenceStrata.process(input, attributePrefix, attributeEvidence, prefix, savedCatalog, savedRoster);
			}else if(attributes){
				NcrnaReferenceAttributes.process(input, prefix, savedCatalog, savedRoster);
			}else if(geometry){
				NcrnaReferenceGeometry.process(input, prefix, savedCatalog, savedRoster);
			}else{NcrnaReferenceProfiles.process(input, prefix, savedCatalog, savedRoster);}
			timer.stop();
			outstream.println("Reference "+(tiers ? "tier construction" : retention ? "native retention" : scores ? "score distributions" : sourceAnchors ? "annotation-concordant evidence" : anchors ? "prior-anchor evidence" : prior ? "prior-decision binding" : strata ? "evidence strata" : attributes ? "attribute join" : geometry ? "geometry" : "profiles")+" completed. Time: "+timer);
			return;
		}
		load();
		final String[] headers=batch ? BATCH_HEADERS : HEADERS;
		try(Outputs outputs=new Outputs(headers.length)){
			writers=new ByteStreamWriter[headers.length];
			for(int i=0; i<writers.length; i++){
				writers[i]=outputs.open(FileFormat.testOutput(prefix+SUFFIXES[i], FileFormat.TXT, null, true, false, false, false));
				writers[i].addJob(new ByteBuilder().append(headers[i]).nl());
			}
			spawnThreads();
		}
		// The summary is deliberately last: partial worker output is never a completed survey.
		if(batch){writeBatchSources();}
		writeCoverage();
		writeSummary();
		timer.stop();
		outstream.println("Inputs observed: "+totals.inputs+"; raw model-hit events: "+totals.hits+"; retained historical rows: "+totals.retained);
		outstream.println("Evidence only; no biological reference loci admitted. Time: "+timer);
	}

	private void spawnThreads(){
		assert(limit>0) : "The nonempty census/collector join supplies work before workers start";
		final int threads=Math.min(limit, Shared.threads());
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){workers.add(new ProcessThread(this, i, threads));}
		try{
			ThreadWaiter.startAndWait(workers, this);
		}catch(RuntimeException | Error failure){
			stop=true;
			final ArrayList<ProcessThread> started=new ArrayList<ProcessThread>();
			for(final ProcessThread worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
			ThreadWaiter.waitForThreadsToFinish(started);
			throw failure;
		}
		for(final ProcessThread worker:workers){
			if(worker.failure!=null){throw new IllegalStateException("Reference evidence worker failed", worker.failure);}
		}
		require(!stop && totals.inputs==(batch ? catalog.size() : limit), "Every selected input must finish before reporting success");
	}

	@Override
	public void accumulate(final ProcessThread worker){
		assert(worker.getState()==Thread.State.TERMINATED) : "Counters are merged only after joining their sole owner";
		totals.add(worker.counts);
	}
	@Override
	public boolean success(){return !stop;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	/** One independent input at a time; failures stop further input acquisition. */
	static final class ProcessThread extends Thread{
		ProcessThread(final BuildNcrnaReference owner_, final int id_, final int stride_){
			owner=owner_;
			id=id_;
			stride=stride_;
			counts=new Counts(owner.roster);
			for(int i=0; i<buffers.length; i++){buffers[i]=new ByteBuilder(BLOCK_SIZE);}
		}
		@Override
		public void run(){
			try{
				for(int i=id; i<owner.limit && !owner.stop; i+=stride){
					if(owner.batch){processBatch(owner.batchTasks.get(i));}else{process(owner.tasks.get(i));}
				}
				for(int i=0; i<buffers.length; i++){flush(i);}
			}catch(Throwable problem){
				failure=problem;
				owner.stop=true;
			}
		}

		private void process(final Task task){
			checkSearchInput(task);
			final Observation observed=owner.catalog.observe(task.genome);
			final Evidence evidence=NcrnaReferenceEvidence.read(task.directory+"/hits.tbl", task.directory+"/completion.tsv", observed, owner.roster);
			require(evidence.size()==task.hits, "Collector and imported raw-hit counts disagree for input "+task.genome.id);
			processEvidence(evidence);
		}

		/** One worker retains at most one complete batch; native reports are read once. */
		private void processBatch(final NcrnaReferenceBatch.Task task){
			final Evidence[] inputs=NcrnaReferenceBatch.read(task, owner.catalog, owner.roster, owner.batchSearchSpace);
			for(final Evidence evidence:inputs){
				if(owner.stop){return;}
				processEvidence(evidence);
			}
		}

		/** Shared observation/history/geometry path for legacy and captured batch inputs. */
		private void processEvidence(final Evidence evidence){
			final Observation observed=evidence.observed;
			final Genome genome=observed.genome;
			require(observed.captured==owner.batch, "Observation basis must match the declared evidence generation");
			for(int i=0; i<owner.roster.queried(); i++){
				if(evidence.hitsForQuery(i)==0){counts.zeroPairs++;}
			}
			final HitIndex index=new HitIndex(evidence);
			writeContigs(observed);
			writeHits(evidence, index);
			final ObjectIntMap<String> sources=new ObjectIntMap<String>(String.class);
			final ArrayList<String> paths=new ArrayList<String>();
			for(int i=0; i<genome.memberships(); i++){
				final String path=genome.membership(i).gff;
				if(!sources.contains(path)){
					sources.put(path, 1);
					paths.add(path);
				}else{sources.put(path, sources.get(path)+1);}
			}
			for(final String path:paths){writeSource(NcrnaHistoricalEvidence.read(path, observed), sources.get(path), index);}
			counts.inputs++;
			counts.bases+=genome.bases;
			counts.contigs+=observed.size();
		}

		/** Binds the native search directory itself, not just the later collector row. */
		private void checkSearchInput(final Task task){
			final String path=task.directory+"/input.tsv";
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, "input_id\tfasta\tnucleotide_bases\tFASTA_records", path);
				final byte[] row=input.reader.nextLine();
				require(row!=null, "Native search input binding is missing: "+path);
				p.set(row);
				require(p.terms()==4 && p.termEquals(task.genome.id, 0) && p.termEquals(task.genome.fasta, 1)
					&& p.parseLong(2)==task.genome.bases && p.parseLong(3)==task.genome.records,
					"Native search input binding disagrees with the exact census source: "+path);
				require(input.reader.nextLine()==null, "Native search input binding must contain exactly one input: "+path);
			}
		}

		private void writeContigs(final Observation observed){
			for(int i=0; i<observed.size(); i++){
				final ByteBuilder b=buffers[CONTIGS];
				field(b, observed.genome.id).tab();
				field(b, observed.genome.fasta).tab();
				field(b, observed.name(i)).tab().append(observed.length(i)).nl();
				flushIfFull(CONTIGS);
			}
		}

		private void writeHits(final Evidence evidence, final HitIndex index){
			for(int i=0; i<evidence.size(); i++){
				final Hit hit=evidence.get(i);
				require(hit.boundInclusion==owner.batch, "Saved schema must distinguish bound native inclusion from legacy printed scores");
				final int model=owner.modelIndex.get(hit.model.accession);
				counts.modelHits[model]++;
				counts.spans[model][hit.modelSpan*4+hit.truncation]++;
				final ByteBuilder b=buffers[HITS];
				key(b, hit.genome.id, hit.table, hit.line);
				field(b, hit.target).tab().append(hit.start).tab().append(hit.stop).tab().append(hit.strand).tab();
				field(b, hit.model.accession).tab();
				field(b, hit.model.release).tab();
				field(b, hit.model.scope.family).tab();
				field(b, hit.model.scope.role).tab();
				b.append(hit.modelFrom).tab().append(hit.modelTo).tab().append(hit.model.length).tab().append(hit.truncation).tab();
				b.append(hit.bits, 6).tab().append(hit.evalueText).tab().append(hit.gaPass).tab();
				field(b, hit.range).tab().append(hit.atLeftEdge).tab().append(hit.atRightEdge);
				if(owner.batch){field(b.tab(), hit.nativeTarget).tab().append(hit.included).tab().append(BOUND_GA_BASIS);}
				b.nl();
				flushIfFull(HITS);
			}
			counts.hits+=evidence.size();
			for(int i=0; i<index.hits.length; i++){
				final Hit a=index.hits[i];
				index.overlaps(a.target, a.start, a.stop, matches);
				for(int j=0; j<matches.size; j++){
					final int other=matches.get(j);
					if(other<=i){continue;}
					final Hit b=index.hits[other];
					link(a.genome.id, "MODEL_MODEL", a.table, a.line, b, overlap(a.start, a.stop, b.start, b.stop),
						a.strand==b.strand ? "SAME_STRAND_GEOMETRY" : "OPPOSITE_STRAND_GEOMETRY", false);
					counts.modelPairs++;
				}
			}
		}

		private void writeSource(final Source source, final int memberships, final HitIndex index){
			final String inputId=source.observed.genome.id;
			ByteBuilder b=buffers[SOURCES];
			field(b, inputId).tab();
			field(b, source.file).tab();
			field(b, source.state()).tab().append(memberships).tab().append(source.featureRows()).tab();
			b.append(source.size()).tab().append(source.excludedRows()).tab();
			b.append(source.excluded("UNSELECTED_CODING")).tab().append(source.excluded("UNSELECTED_GENERIC_STRUCTURE")).tab();
			b.append(source.excluded("UNSELECTED_REGION")).tab().append(source.missingParents()).tab().append(source.cycles()).nl();
			flushIfFull(SOURCES);
			counts.sources++;
			counts.features+=source.featureRows();
			counts.retained+=source.size();
			counts.excluded+=source.excludedRows();
			if(!source.state().equals("PARSED_GFF") && !source.state().equals("EMPTY_GFF")){counts.unavailableSources++;}
			final ObjectSet<String> missing=new ObjectSet<String>(String.class);
			for(int i=0; i<source.missingParents(); i++){missing.add(source.missingParent(i));}
			for(int i=0; i<source.cycles(); i++){
				final NcrnaHistoricalEvidence.Relation cycle=source.cycle(i);
				relationship(inputId, source.file, cycle.line, "CYCLE", cycle.child, cycle.parent, "BACK_EDGE_WITNESS");
			}
			for(int i=0; i<source.size(); i++){
				final Feature feature=source.get(i);
				for(final String parent:feature.parents){
					relationship(inputId, source.file, feature.line, "PARENT", feature.id, parent,
						missing.contains(parent) ? "DANGLING" : "PRESENT_ID");
				}
				for(final String rfam:feature.rfam){relationship(inputId, source.file, feature.line, "RFAM", feature.id, rfam, "SOURCE_CLAIM");}
				final String state=assessability(feature);
				Hit best=null;
				int same=0, opposite=0, invalidHit=0;
				matches.clear();
				if(state.equals("ASSESSABLE")){
					index.overlaps(feature.target, feature.start, feature.stop, matches);
					for(int j=0; j<matches.size; j++){
						final Hit hit=index.hits[matches.get(j)];
						if(!hit.range.equals("VALID_INTERVAL")){
							invalidHit++;
						}else if(hit.strand!=feature.strand){
							opposite++;
						}else{
							same++;
							if(best==null || NcrnaReferenceEvidence.compareSupport(hit, best, feature.start, feature.stop)<0){best=hit;}
						}
					}
					for(int j=0; j<matches.size; j++){
						final Hit hit=index.hits[matches.get(j)];
						final String relation=!hit.range.equals("VALID_INTERVAL") ? "UNASSESSABLE_HIT_RANGE"
							: hit.strand==feature.strand ? "SAME_STRAND_CANDIDATE" : "OPPOSITE_STRAND_GEOMETRY";
						link(inputId, "HISTORICAL_MODEL", source.file, feature.line, hit,
							overlap(feature.start, feature.stop, hit.start, hit.stop), relation, hit==best);
						counts.historicalPairs++;
					}
				}
				b=buffers[HISTORY];
				key(b, inputId, source.file, feature.line);
				field(b, feature.target).tab();
				field(b, feature.type).tab().append(feature.start).tab().append(feature.stop).tab().append(feature.strand).tab();
				field(b, feature.id==null ? "NA" : feature.id).tab().append(feature.selection()).tab();
				field(b, feature.range).tab().append(feature.incompatibleId()).tab();
				field(b, state).tab().append(feature.parents.length).tab().append(feature.rfam.length).tab();
				b.append(same).tab().append(opposite).tab().append(invalidHit).tab();
				b.append(best==null ? -1 : best.line);
				if(owner.batch){field(b.tab(), best==null ? "NA" : best.table);}
				b.nl();
				flushIfFull(HISTORY);
			}
		}

		/** Identifiers supplement row-level source pointers; Parent IDs are never unique row keys. */
		private void relationship(final String input, final String gff, final long line, final String kind,
			final String subject, final String object, final String state){
			final ByteBuilder b=buffers[RELATIONS];
			key(b, input, gff, line);
			field(b, kind).tab();
			field(b, subject==null ? "NA" : subject).tab();
			field(b, object).tab();
			field(b, state).nl();
			flushIfFull(RELATIONS);
		}

		private void link(final String inputId, final String kind, final String leftFile, final long leftLine,
			final Hit right, final long bases, final String relation, final boolean best){
			assert(bases>0) : "The interval index emits actual overlaps, never shared-window membership";
			final ByteBuilder b=buffers[LINKS];
			field(b, inputId).tab();
			field(b, kind).tab();
			field(b, leftFile).tab().append(leftLine).tab();
			field(b, right.table).tab().append(right.line).tab().append(bases).tab();
			field(b, relation).tab().append(best).nl();
			flushIfFull(LINKS);
		}

		private void flushIfFull(final int output){if(buffers[output].length()>=BLOCK_SIZE){flush(output);}}
		/** Ownership crosses the native writer queue once per block, never once per field. */
		private void flush(final int output){
			if(buffers[output].length()==0){return;}
			owner.writers[output].addJob(buffers[output]);
			buffers[output]=new ByteBuilder(BLOCK_SIZE);
		}
		final BuildNcrnaReference owner;
		final int id, stride;
		final Counts counts;
		final ByteBuilder[] buffers=new ByteBuilder[HEADERS.length];
		final IntList matches=new IntList();
		Throwable failure;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Interval index       ----------------*/
	/*--------------------------------------------------------------*/

	/** Sorted per-target intervals with prefix maximum ends; no locus merging. */
	public static final class HitIndex{
		public HitIndex(final Evidence evidence){
			hits=new Hit[evidence.size()];
			for(int i=0; i<hits.length; i++){hits[i]=evidence.get(i);}
			Arrays.sort(hits, HIT_ORDER);
			maxEnd=new int[hits.length];
			for(int first=0; first<hits.length; ){
				int end=first, maximum=0;
				while(end<hits.length && hits[first].target.equals(hits[end].target)){
					maximum=Math.max(maximum, hits[end].stop);
					maxEnd[end++]=maximum;
				}
				ranges.put(hits[first].target, new Slice(first, end));
				first=end;
			}
		}
		/** Returns sorted-array indices with inclusive overlap on either strand. */
		public void overlaps(final String target, final long start, final long stop, final IntList result){
			require(start>0 && stop>=start, "Interval lookup requires positive ordered coordinates");
			result.clear();
			final Slice slice=ranges.get(target);
			if(slice==null){return;}
			int low=slice.first, high=slice.end;
			while(low<high){
				final int mid=(low+high)>>>1;
				if(hits[mid].start<=stop){
					low=mid+1;
				}else{high=mid;}
			}
			for(int i=low-1; i>=slice.first && maxEnd[i]>=start; i--){
				if(hits[i].stop>=start){result.add(i);}
			}
		}
		public final Hit[] hits;
		private final int[] maxEnd;
		private final ObjectMap<String,Slice> ranges=new ObjectMap<String,Slice>(String.class, Slice.class);
	}
	private static final class Slice{
		Slice(final int first_, final int end_){
			first=first_;
			end=end_;
		}
		final int first, end;
	}

	/** Context ancestors are not genes; uncertain strand/ranges are not negative evidence. */
	public static String assessability(final Feature feature){
		if((feature.selection()&~NcrnaHistoricalEvidence.ANCESTOR)==0){return "STRUCTURAL_CONTEXT";}
		if(!feature.range.equals("VALID_INTERVAL")){return feature.range;}
		if(feature.strand!='+' && feature.strand!='-'){return "UNASSESSABLE_STRAND";}
		return "ASSESSABLE";
	}
	private static long overlap(final long a, final long b, final long c, final long d){return EukPanelCmEvidence.overlap(a, b, c, d);}
	private static final Comparator<Hit> HIT_ORDER=(a, b)->{
		int c=a.target.compareTo(b.target);
		if(c==0){c=Integer.compare(a.start, b.start);}
		if(c==0){c=Integer.compare(a.stop, b.stop);}
		if(c==0){c=a.table.compareTo(b.table);}
		return c==0 ? Long.compare(a.line, b.line) : c;
	};

	/*--------------------------------------------------------------*/
	/*----------------         Result accounting    ----------------*/
	/*--------------------------------------------------------------*/

	/** Declares captured map/FASTA provenance once per batch, without rereading payloads. */
	private void writeBatchSources(){
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(prefix+".batch_sources.tsv")){
			out.writer.println("batch_id\tprepared_directory\tnative_directory\tfasta\tfasta_sha80\tmapping\tmapping_sha80\tmembers\tmembers_sha80\tsearch_space_Mb\tobservation_basis");
			for(final NcrnaReferenceBatch.Task task:batchTasks){
				b.clear().append(task.id).tab();
				field(b, task.prepared).tab();
				field(b, task.results).tab();
				field(b, task.fasta).tab().append(task.fastaSha80).tab();
				field(b, task.mapping).tab().append(task.mapSha80).tab();
				field(b, task.members).tab().append(task.membersSha80).tab().append(batchSearchSpace, 6);
				out.writer.print(b.append("\tCAPTURED_PREPARATION_MAP_NOT_CURRENT_FASTA_READ\n"));
			}
		}
	}

	private static final class Task{
		Task(final Genome genome_, final String directory_, final long hits_){
			genome=genome_;
			directory=directory_;
			hits=hits_;
		}
		final Genome genome;
		final String directory;
		final long hits;
	}
	private static final class Counts{
		Counts(final Roster roster){
			modelHits=new long[roster.queried()];
			spans=new long[modelHits.length][];
			for(int i=0; i<spans.length; i++){spans[i]=new long[(roster.query(i).length+1)*4];}
		}
		void add(final Counts other){
			assert(modelHits.length==other.modelHits.length) : "All workers use the same immutable queried-model roster";
			inputs+=other.inputs;
			bases+=other.bases;
			contigs+=other.contigs;
			hits+=other.hits;
			zeroPairs+=other.zeroPairs;
			sources+=other.sources;
			features+=other.features;
			retained+=other.retained;
			excluded+=other.excluded;
			unavailableSources+=other.unavailableSources;
			historicalPairs+=other.historicalPairs;
			modelPairs+=other.modelPairs;
			for(int i=0; i<spans.length; i++){
				modelHits[i]+=other.modelHits[i];
				for(int j=0; j<spans[i].length; j++){spans[i][j]+=other.spans[i][j];}
			}
		}
		long inputs, bases, contigs, hits, zeroPairs, sources, features, retained, excluded, unavailableSources, historicalPairs, modelPairs;
		final long[] modelHits;
		final long[][] spans;
	}

	/** Exact model-position counts, not a universal biological completeness threshold. */
	private void writeCoverage(){
		try(Output output=new Output(prefix+".coverage.tsv")){
			final ByteBuilder b=new ByteBuilder();
			output.writer.println(COVERAGE_HEADER);
			for(int i=0; i<roster.queried(); i++){
				final Model model=roster.query(i);
				long sum=0;
				if(totals.modelHits[i]==0){
					field(b.clear(), model.accession).tab();
					field(b, model.release).tab().append(model.length).append("\tNA\tNA\t0\t0\t").append(totals.inputs).nl();
					output.writer.print(b);
				}
				for(int span=1; span<=model.length; span++){
					for(int truncation=0; truncation<4; truncation++){
						final long count=totals.spans[i][span*4+truncation];
						if(count==0){continue;}
						field(b.clear(), model.accession).tab();
						field(b, model.release).tab().append(model.length).tab().append(span).tab().append(truncation).tab();
						output.writer.print(b.append(count).tab().append(totals.modelHits[i]).tab().append(totals.inputs).nl());
						sum+=count;
					}
				}
				require(sum==totals.modelHits[i], "Exact span/truncation histogram must conserve every raw model hit");
			}
		}
	}
	private void writeSummary(){
		require(totals.features==totals.retained+totals.excluded, "Historical retained/excluded rows must conserve every parsed source feature");
		try(Output output=new Output(prefix+".summary.tsv")){
			output.writer.println(SUMMARY_HEADER);
			final ByteBuilder b=new ByteBuilder();
			b.append(totals.inputs).tab().append(catalog.size()).tab().append(totals.bases).tab().append(totals.contigs).tab();
			b.append(totals.hits).tab().append(totals.sources).tab().append(totals.unavailableSources).tab().append(totals.features).tab();
			b.append(totals.retained).tab().append(totals.excluded).tab().append(totals.historicalPairs).tab().append(totals.modelPairs).tab();
			b.append(roster.queried()).tab().append(totals.inputs*roster.queried()).tab().append(totals.zeroPairs).tab();
			output.writer.print(b.append(batch ? BATCH_STATUS : limit==tasks.size() ? "EVIDENCE_COMPLETE_NOT_ADJUDICATED" : "SAMPLE_EVIDENCE_NOT_ADJUDICATED").nl());
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Output framing       ----------------*/
	/*--------------------------------------------------------------*/

	private static void key(final ByteBuilder b, final String input, final String source, final long line){
		field(b, input).tab();
		field(b, source).tab().append(line).tab();
	}
	/** Ordinary ASCII fields allocate nothing; unusual Unicode uses a lossless UTF-8 fallback. */
	public static ByteBuilder field(final ByteBuilder b, final String text){
		assert(text!=null) : "Missing metadata is explicit NA, never an invented empty field";
		for(int i=0; i<text.length(); i++){
			if(text.charAt(i)>127){
				for(final byte value:text.getBytes(StandardCharsets.UTF_8)){fieldByte(b, value);}
				return b;
			}
		}
		for(int i=0; i<text.length(); i++){fieldByte(b, (byte)text.charAt(i));}
		return b;
	}
	private static void fieldByte(final ByteBuilder b, final byte value){
		if(value=='%' || value=='\t' || value=='\r' || value=='\n'){
			b.append('%').append(HEX.charAt((value&255)>>>4)).append(HEX.charAt(value&15));
		}else{b.append(value);}
	}
	/** Frames an existing UTF-8 field without allocating an intermediate String. */
	static ByteBuilder field(final ByteBuilder b, final byte[] text, final int from, final int to){
		assert(from>=0 && from<=to && to<=text.length) : "TSV escaping must stay within the source field";
		for(int i=from; i<to; i++){fieldByte(b, text[i]);}
		return b;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}
	private static String utf8(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private String input, prefix, census, memberships, scope, assets, queries, attributePrefix, attributeEvidence, marginalPercent, collapsePercent, priorManifest, comparison, anchorManifest, sourcePrefix, sampleModel;
	private String nativePaths;
	private String panelReference, panelContigs, panelSource;
	private boolean panel;
	private boolean inventory;
	private boolean selectPanel;
	private String selectionInventory, selectionPolicy, selectionTaxonomy, selectionFamily;
	private String inventoryMembers, inventoryWitnesses, inventoryEvents;
	private boolean profile, scores, retention, tiers, tierReport, geometry, attributes, strata, prior, cost, panels, anchors, sourceAnchors, samplingSpecified, batch, examplesSpecified;
	private int examplesPerGroup=3;
	private double batchSearchSpace=-1;
	private int maxInputs=Integer.MAX_VALUE, limit;
	private final PrintStream outstream;
	private NcrnaEvidenceCatalog catalog;
	private Roster roster;
	private Counts totals;
	private final ArrayList<Task> tasks=new ArrayList<Task>();
	private final ArrayList<NcrnaReferenceBatch.Task> batchTasks=new ArrayList<NcrnaReferenceBatch.Task>();
	private final ObjectIntMap<String> modelIndex=new ObjectIntMap<String>(String.class);
	private ByteStreamWriter[] writers;
	private volatile boolean stop;
	private static final int CONTIGS=0, HITS=1, SOURCES=2, HISTORY=3, LINKS=4, RELATIONS=5, BLOCK_SIZE=65536;
	private static final String HEX="0123456789ABCDEF";
	public static final String COLLECTOR_HEADER="input_id\tfasta\tresult_directory\tmodels\traw_model_hits\tsearch_CPU_seconds\tdecode_CPU_seconds\tmax_native_RSS_KiB";
	public static final String COVERAGE_HEADER="accession\trelease\tmodel_length\tmodel_span\tnative_truncation\traw_hit_events\taccession_hit_denominator\tinputs_observed";
	public static final String SUMMARY_HEADER="inputs_observed\tcatalog_inputs\tinput_bases\tcontigs\traw_hit_events\tunique_input_gff_claims\tunavailable_gff_claims\tsource_feature_rows\tretained_source_rows\texcluded_source_rows\thistorical_model_pairs\tmodel_model_pairs\tqueried_models\tcompleted_input_model_pairs\tzero_hit_pairs\tstatus";
	public static final String SOURCE_HEADER="input_id\tgff\tstate\tmemberships\tfeature_rows\tretained_rows\texcluded_rows\texcluded_coding\texcluded_generic_structure\texcluded_region\tmissing_parents\tcycle_witnesses";
	public static final String HIT_HEADER="input_id\tnative_table\tnative_line\ttarget\tstart\tstop\tstrand\taccession\trelease\treference_family\tevidence_role\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tbits\tevalue\tGA_pass\trange\tleft_contig_edge\tright_contig_edge";
	public static final String BATCH_HIT_HEADER=HIT_HEADER+"\tnative_target\tnative_included\tGA_decision_basis";
	public static final String BATCH_STATUS="BATCH_EVIDENCE_COMPLETE_NOT_ADJUDICATED";
	public static final String BOUND_GA_BASIS="ORIGINAL_GA_NATIVE_INCLUSION";
	public static final String BATCH_COLLECTOR_HEADER="batch_id\tprepared_directory\tresult_directory\tmodels\traw_model_hits";
	public static final String NATIVE_PATH_HEADER="batch_id\toriginal_native_directory\tphysical_native_directory";
	public static final String OVERLAP_HEADER="input_id\tpair_kind\tleft_source\tleft_line\tnative_table\tnative_line\toverlap_bases\trelation\tdeterministic_best_candidate";
	public static final String HISTORICAL_HEADER="input_id\tgff\tline\ttarget\ttype\tstart\tstop\tstrand\tID\tselection_bits\trange\tincompatible_ID\toverlap_state\tparent_count\tRFAM_count\tsame_strand_candidates\topposite_strand_overlaps\tinvalid_hit_overlaps\tbest_native_line";
	public static final String BATCH_HISTORICAL_HEADER=HISTORICAL_HEADER+"\tbest_native_table";
	public static final String RELATION_HEADER="input_id\tgff\tline\trelationship\tsubject_ID\tobject\tstate";
	private static final String[] SUFFIXES={".contigs.tsv.gz", ".hits.tsv.gz", ".sources.tsv.gz", ".historical.tsv.gz", ".overlaps.tsv.gz", ".relations.tsv.gz", ".coverage.tsv", ".summary.tsv"};
	private static final String[] BATCH_SUFFIXES={".contigs.tsv.gz", ".hits.tsv.gz", ".sources.tsv.gz", ".historical.tsv.gz", ".overlaps.tsv.gz", ".relations.tsv.gz", ".coverage.tsv", ".summary.tsv", ".batch_sources.tsv"};
	private static final String[] HEADERS={
		"input_id\tfasta\tliteral_target\tobserved_length",
		HIT_HEADER,
		SOURCE_HEADER,
		HISTORICAL_HEADER,
		OVERLAP_HEADER,
		RELATION_HEADER
	};
	private static final String[] BATCH_HEADERS={HEADERS[0], BATCH_HIT_HEADER, SOURCE_HEADER, BATCH_HISTORICAL_HEADER, OVERLAP_HEADER, RELATION_HEADER};
}
