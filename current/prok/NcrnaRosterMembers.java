package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.IntObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.LongList;
import structures.DoubleList;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** Extracts a complete reconciled RF00174 roster with construction-only exact
 * sequence representatives within accepted assembly groups. Original calls and
 * overlap components retain separate ledgers; neither disappears from recall.
 * Derived from template.A_SampleByteFileMT's prepass/worker/accumulator structure;
 * workers own synchronous native FASTA streams rather than shared text lines.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaRosterMembers implements Accumulator<NcrnaRosterMembers.ProcessThread>{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, NcrnaRosterMembers.class, false);
		try{
			final NcrnaRosterMembers tool=new NcrnaRosterMembers(pp.args, pp.outstream);
			tool.process();
			timer.stop();
			pp.outstream.println("NCRNA_ROSTER_MEMBERS_COMPLETE scope="+(tool.metadataOnly ? "METADATA_ONLY" : tool.canary ? "OPERATIONAL_CANARY" : "FULL_ROSTER")+" time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaRosterMembers(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		String mode="full";
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("expectedloci")){expectedLoci=Parse.parseIntKMG(b);}
			else if(a.equals("cohorts")){cohortManifest=b;}
			else if(a.equals("mode")){mode=b;}
			else if(a.equals("maxsources")){maxSources=Parse.parseIntKMG(b);}
			else if(a.equals("modeltolerance")){modelTolerance=Parse.parseIntKMG(b);}
			else if(a.equals("coveragefilter")){coverageFilter=Parse.parseBoolean(b);}
			else if(a.equals("maxcoveragefraction")){maxCoverageFraction=Double.parseDouble(b);}
			else{require(parser.parse(arg, a, b), "Unknown roster extraction argument: "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		canary="canary".equals(mode);
		metadataOnly="metadata".equals(mode);
		require(in!=null && new File(in).isDirectory() && cohortManifest!=null && out!=null && expectedLoci>0
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && (canary || metadataOnly || "full".equals(mode))
			&& modelTolerance>=0 && maxCoverageFraction>=0 && maxCoverageFraction<=1
			&& (canary ? maxSources>0 : maxSources==0),
			"Supply in=reconciliation_result cohorts=accepted_sources.tsv out=fresh_directory expectedloci=N; only explicit mode=canary permits maxsources=N");
		require(!Shared.AMINO_IN && !Read.IUPAC_TO_N && !Read.U_TO_T && !Read.DOT_DASH_X_TO_N
			&& !Read.LOWER_CASE_TO_N && Read.JUNK_MODE==Read.CRASH_JUNK,
			"Exact construction sequence equality forbids lossy nucleotide conversion");
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_PIGZ=false;
		ReadWrite.setZipThreads(1);
	}

	private void process(){
		loadSummary();
		loadFrames();
		loadLoci();
		loadLinks();
		loadOriginal();
		loadNativeEligibility();
		require(!new File(out).exists() && new File(out).mkdir(), "Extraction outputs require a fresh directory: "+out);
		writeEligibility();
		if(metadataOnly){return;}
		require(!coverageFilter || coverageFailures<=maxCoverageFraction*expectedOriginal,
			"Incomplete model coverage exceeded the allowed roster fraction; metadata tables were saved, no sequence was read");
		final ArrayList<Frame> selected=new ArrayList<Frame>();
		for(final Frame frame:frameList){
			if(frame.genes.isEmpty()){continue;}
			if(canary && selected.size()>=maxSources){continue;}
			frame.selected=true;
			selected.add(frame);
		}
		require(!selected.isEmpty(), "The complete roster must bind at least one source frame");
		final Frame[] jobs=selected.toArray(new Frame[0]);
		outstream.println("Roster extraction: "+expectedOriginal+" original called groups, "+genes.size()+" merged loci, "+jobs.length+" selected source frames.");
		final AtomicInteger next=new AtomicInteger();
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
		for(int i=0; i<Math.min(Shared.threads(), jobs.length); i++){workers.add(new ProcessThread(jobs, next));}
		final boolean complete=ThreadWaiter.startAndWait(workers, this);
		for(final ProcessThread worker:workers){
			if(worker.failure!=null){throw new IllegalStateException("Native roster extraction failed", worker.failure);}
		}
		require(complete && !errorState, "Every synchronous source worker must terminate before publishing gene sequences");
		selectRepresentatives();
		write();
		writeGenomicContext();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Metadata Prepass      ----------------*/
	/*--------------------------------------------------------------*/

	private void loadSummary(){
		try(Input input=new Input(input("summary.tsv"))){
			RrnaResourceIO.header(input.reader, "metric\tcount\tqualification", in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==3, "Reconciliation summary requires count and qualification");
				if(p.termEquals("qualified_same_frame_overlap_loci", 0)){
					require(p.parseInt(1)==expectedLoci, "Declare the complete overlap-locus population");
					summaryBound=true;
				}
				else if(p.termEquals("qualified_exact_coordinate_groups", 0)){expectedOriginal=p.parseInt(1);}
				else if(p.termEquals("native_GA_RF00174_events", 0)){nativeEvents=p.parseLong(1);}
				else if(p.termEquals("raw_RF00174_events", 0)){rawEvents=p.parseLong(1);}
			}
		}
		require(summaryBound && expectedOriginal>=expectedLoci && nativeEvents>=expectedOriginal && rawEvents>=nativeEvents,
			"Raw events, native exact groups and overlap loci must retain their separate complete populations");
	}

	private void loadFrames(){
		try(Input input=new Input(input("source_frames.tsv.gz"))){
			RrnaResourceIO.header(input.reader, FRAME_HEADER, in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==13, "Use the complete reconciler source-frame schema");
				final Frame frame=new Frame(frameList.size(), p.parseInt(0), text(1), text(2), p.parseLong(4), p.parseLong(5), text(12));
				require(frame.group>=0 && frame.canonical.startsWith("/") && frame.fasta.startsWith("/")
					&& frame.bases>0 && frame.records>0 && !frames.contains(frame.canonical), "Saved frames must be unique absolute nucleotide sources");
				frames.put(frame.canonical, frame);
				bindSource(frame.canonical, frame);
				bindSource(frame.fasta, frame);
				frameList.add(frame);
			}
		}
		require(!frameList.isEmpty(), "A complete native roster needs its source-frame census");
	}

	private void bindSource(final String path, final Frame frame){
		final Frame old=sourceNames.get(path);
		require(old==null || old==frame, "One saved source name cannot bind conflicting coordinate frames: "+path);
		if(old==null){sourceNames.put(path, frame);}
	}

	private void loadLoci(){
		try(Input input=new Input(input("merged_loci.tsv.gz"))){
			RrnaResourceIO.header(input.reader, MERGED_HEADER, in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==10 && p.parseInt(0)==genes.size(), "Merged locus IDs must retain their complete original ordinal");
				final Frame frame=frames.get(text(2));
				require(frame!=null && frame.group==p.parseInt(1), "Each merged locus must bind its exact saved frame and accepted group");
				final String contig=text(3);
				Target target=frame.targets.get(contig);
				if(target==null){
					target=new Target(contig);
					frame.targets.put(contig, target);
				}
				require(p.termEquals("+", 6) || p.termEquals("-", 6), "Sense extraction requires the original explicit strand");
				final Gene gene=new Gene(genes.size(), frame, target, p.parseInt(4), p.parseInt(5), p.termEquals("-", 6));
				require(gene.start>0 && gene.stop>=gene.start, "Merged loci require positive inclusive bounds");
				genes.add(gene);
				frame.genes.add(gene);
				target.genes.add(gene);
			}
		}
		require(genes.size()==expectedLoci, "No overlap locus may disappear before sequence extraction");
	}

	private void loadLinks(){
		try(Input input=new Input(input("locus_components.tsv.gz"))){
			RrnaResourceIO.header(input.reader, "exact_coordinate_group_id\tmerged_locus_id", in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				final int original=p.parseInt(0), merged=p.parseInt(1);
				require(p.terms()==2 && original>=0 && merged>=0 && merged<genes.size() && !links.contains(original),
					"Each original called coordinate group needs exactly one valid overlap component");
				links.put(original, merged);
			}
		}
		require(links.size()==expectedOriginal, "Component links must conserve all original called coordinate groups");
	}

	private void loadOriginal(){
		final IntHashMap seen=new IntHashMap();
		try(Input source=new Input(input("loci.tsv.gz"))){
			RrnaResourceIO.header(source.reader, ORIGINAL_HEADER, in);
			for(byte[] row=source.reader.nextLine(); row!=null; row=source.reader.nextLine()){
				p.set(row);
				require(p.terms()==13, "Retain all original called-coordinate columns");
				if(!p.termEquals("true", 7)){continue;}
				final int id=p.parseInt(0);
				require(links.contains(id) && !seen.contains(id), "Every original call needs exactly one declared component link");
				seen.put(id, 1);
				final Gene component=genes.get(links.get(id));
				final int start=p.parseInt(4), stop=p.parseInt(5);
				require(p.parseInt(1)==component.frame.group && text(2).equals(component.frame.canonical)
					&& text(3).equals(component.target.name) && start>=component.start && stop<=component.stop
					&& p.termEquals(component.minus ? "-" : "+", 6), "Original calls must stay inside their same-frame overlap component");
				final OriginalCall call=new OriginalCall(row, component, start, stop);
				for(final OriginalCall old:component.target.originals){
					require(old.start!=start || old.stop!=stop || old.component.minus!=component.minus,
						"Exact-coordinate source groups cannot contain repeated bounds");
				}
				component.target.originals.add(call);
				originalCalls.add(call);
			}
		}
		require(originalCalls.size()==expectedOriginal, "Original native called genes cannot disappear before eligibility binding");
	}

	/** Reads only native event fields. Saved annotations/library membership never
	 * choose a construction interval, strand or partial-hit disposition. */
	private void loadNativeEligibility(){
		final ArrayList<Cohort> cohorts=new ArrayList<Cohort>();
		try(Input manifest=new Input(cohortManifest)){
			RrnaResourceIO.header(manifest.reader, COHORT_HEADER, cohortManifest);
			for(byte[] row=manifest.reader.nextLine(); row!=null; row=manifest.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Use the accepted complete scan-cohort source manifest");
				final Cohort cohort=new Cohort(text(0));
				final String census=text(1);
				final int expected=p.parseInt(5);
				for(final Cohort old:cohorts){require(!old.label.equals(cohort.label), "Native cohort labels must be unique");}
				try(Input input=new Input(census)){
					RrnaResourceIO.header(input.reader, CENSUS_HEADER, census);
					for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
						p.set(line);
						require(p.terms()==5, "Native census bindings require complete source rows");
						final int id=p.parseInt(0);
						final Frame frame=sourceNames.get(text(1));
						require(frame!=null && id>=0 && !cohort.inputs.contains(id), "Native input IDs must bind the saved frame, without resolving a mutable alias");
						cohort.inputs.put(id, frame);
					}
				}
				require(cohort.inputs.size()==expected, "Every accepted cohort input must retain its source-frame binding");
				cohorts.add(cohort);
			}
		}
		long observed=0;
		try(Input input=new Input(input("events.tsv.gz"))){
			RrnaResourceIO.header(input.reader, EVENT_HEADER, in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==35 && p.termEquals("RF00174", 11), "Eligibility must use retained native RF00174 events only");
				if(!p.termEquals("true", 2)){continue;}
				require(p.termEquals("true", 22), "Qualified events must retain the original native GA inclusion");
				Cohort cohort=null;
				for(final Cohort candidate:cohorts){
					if(p.termEquals(candidate.label, 0)){
						cohort=candidate;
						break;
					}
				}
				require(cohort!=null, "Every native event requires its explicit accepted cohort");
				final Frame frame=cohort.inputs.get(p.parseInt(6));
				require(frame!=null && frame.group==p.parseInt(1), "A native event cannot change its accepted assembly/source binding");
				p.setBounds(7);
				final Target target=frame.targets.get(decoder.decode(p.line(), p.a(), p.b()));
				require(target!=null, "Every qualified native event must retain its original contig");
				final int start=p.parseInt(8), stop=p.parseInt(9);
				final boolean minus=p.termEquals("-", 10);
				require(minus || p.termEquals("+", 10), "Native strand must be explicit");
				OriginalCall call=null;
				for(final OriginalCall candidate:target.originals){
					if(candidate.start==start && candidate.stop==stop && candidate.component.minus==minus){
						call=candidate;
						break;
					}
				}
				require(call!=null, "A native hit's own bounds/strand must match one original called coordinate group");
				final int from=p.parseInt(13), to=p.parseInt(14), length=p.parseInt(15), truncation=p.parseInt(16);
				require(from>0 && to>=from && length>=to && truncation>=0 && truncation<=3, "Retain valid native model coordinates and documented0–3truncation bits");
				if(call.nativeHits==0){
					call.modelFrom=from;
					call.modelTo=to;
					call.modelLength=length;
				}else{
					require(call.modelFrom==from && call.modelTo==to && call.modelLength==length,
						"Repeated native events at one exact interval must retain consistent covariance-model coordinates");
				}
				call.nativeHits++;
				final double bits=p.parseDouble(19);
				require(Double.isFinite(bits), "Native printed bit scores must remain finite for endpoint-context diagnostics");
				call.bits=Double.isNaN(call.bits) ? bits : Math.max(call.bits, bits);
				call.truncation|=truncation;
				call.incomplete|=from>modelTolerance || to<length-modelTolerance;
				observed++;
			}
		}
		require(observed==nativeEvents, "Every native-GA hit must survive eligibility binding");
		for(final OriginalCall call:originalCalls){
			require(call.nativeHits>0, "Every original called gene needs native evidence, never annotation-only support");
			missingLeft.increment(call.modelFrom-1);
			missingRight.increment(call.modelLength-call.modelTo);
			groupCalls.increment(call.component.frame.group);
			if(call.incomplete){groupIncomplete.increment(call.component.frame.group);}
			if(call.incomplete){coverageFailures++;}
			if(call.truncation!=0 && call.incomplete){bothPartial++;}
			else if(call.truncation!=0){truncatedOnly++;}
			else if(call.incomplete){coverageOnly++;}
			else{neitherPartial++;}
			call.exclusion=partialReason(call);
			if(call.exclusion!=null){partialOriginal++;}
			if(call.start==call.component.start && call.stop==call.component.stop){
				call.component.nativeBoundsSeen=true;
				if(call.exclusion==null){
					call.component.nativeEligible=true;
					call.component.nativeCall=call;
				}
			}
		}
		require(truncatedOnly+coverageOnly+bothPartial+neitherPartial==expectedOriginal
			&& coverageFailures==coverageOnly+bothPartial,
			"Disjoint native eligibility classes must conserve every original call, including overlapping exclusion evidence");
	}

	private String partialReason(final OriginalCall call){
		return call.truncation!=0 ? (coverageFilter && call.incomplete ? "NATIVE_TRUNCATION_AND_INCOMPLETE_MODEL_SPAN" : "NATIVE_TRUNCATION")
			: coverageFilter && call.incomplete ? "INCOMPLETE_NATIVE_MODEL_SPAN" : null;
	}

	/** Metadata mode stops here: neither the FASTA reader nor clustering is started. */
	private void writeEligibility(){
		final ByteBuilder row=new ByteBuilder();
		try(Output table=new Output(output("eligibility.tsv")); Output distribution=new Output(output("missing_model_positions.tsv"))){
			table.writer.println("metric\tcount");
			metric(table, row, "original_called_coordinate_groups", expectedOriginal);
			metric(table, row, "infernal_truncated_only", truncatedOnly);
			metric(table, row, "incomplete_model_coverage_only", coverageOnly);
			metric(table, row, "infernal_truncated_and_incomplete_model_coverage", bothPartial);
			metric(table, row, "neither_partial_diagnostic", neitherPartial);
			metric(table, row, "construction_candidates", expectedOriginal-partialOriginal);
			int incompleteOnlyGroups=0;
			for(final int group:groupCalls.keys()){
				if(group>=0 && groupIncomplete.get(group)==groupCalls.get(group)){incompleteOnlyGroups++;}
			}
			metric(table, row, "assembly_groups_with_called_genes", groupCalls.size());
			metric(table, row, "assembly_groups_only_incomplete_model_coverage", incompleteOnlyGroups);
			metric(table, row, "infernal_truncated_total", truncatedOnly+bothPartial);
			metric(table, row, "incomplete_model_coverage_total", coverageFailures);
			metric(table, row, "excluded_total", partialOriginal);
			metric(table, row, "model_tolerance", modelTolerance);
			table.writer.println("maximum_coverage_failure_fraction\t"+maxCoverageFraction);
			table.writer.println("coverage_filter\t"+coverageFilter);
			table.writer.println("coverage_fraction_gate\t"+(!coverageFilter ? "INACTIVE_COVERAGE_FILTER_OFF"
				: coverageFailures<=maxCoverageFraction*expectedOriginal ? "PASS" : "STOP_BEFORE_SEQUENCE_EXTRACTION"));
			distribution.writer.println("end\tmissing_model_positions\toriginal_calls");
			writeDistribution(distribution, row, "left", missingLeft);
			writeDistribution(distribution, row, "right", missingRight);
		}
		try(Output calls=new Output(output("eligibility_calls.tsv.gz"))){
			calls.writer.println("original_locus_id\tmerged_locus_id\tnative_hit_count\tnative_truncation_bits\tmodel_from\tmodel_to\tmodel_length\tincomplete_model_coverage\tconstruction_eligibility");
			for(final OriginalCall call:originalCalls){
				p.set(call.row);
				calls.writer.print(row.clear().append(p.parseInt(0)).tab().append(call.component.id).tab().append(call.nativeHits)
					.tab().append(call.truncation).tab().append(call.modelFrom).tab().append(call.modelTo).tab().append(call.modelLength)
					.tab().append(call.incomplete).tab().append(call.exclusion==null ? "ELIGIBLE" : call.exclusion).nl());
			}
		}
		try(Output parameters=new Output(output("eligibility_parameters.tsv"))){
			parameters.writer.println("source\tNATIVE_CMSEARCH_RF00174_ONLY_NO_ANNOTATION_OR_OLD_LIBRARY_MEMBERS");
			parameters.writer.println("rule_authority\tG11_20261115_094204Z_WITHDREW_COVERAGE_EXCLUSION");
			parameters.writer.println("rule\tNATIVE_TRUNCATION_NONZERO_ONLY_UNLESS_COVERAGE_FILTER_EXPLICITLY_ENABLED");
			parameters.writer.println("coverage_diagnostic\tMODEL_FROM_GREATER_THAN_TOLERANCE_OR_MODEL_TO_LESS_THAN_MODEL_LENGTH_MINUS_TOLERANCE");
			parameters.writer.println("contig_edge_alone\tNOT_AN_EXCLUSION");
			parameters.writer.println("recall\tALL_ORIGINAL_CALLS_RETAINED");
		}
	}

	private void writeDistribution(final Output output, final ByteBuilder row, final String end, final LongList counts){
		long total=0;
		for(int i=0; i<counts.size; i++){
			final long count=counts.get(i);
			if(count==0){continue;}
			output.writer.print(row.clear().append(end).tab().append(i).tab().append(count).nl());
			total+=count;
		}
		require(total==expectedOriginal, "Each native model-end histogram must conserve every original call");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Source Extraction      ----------------*/
	/*--------------------------------------------------------------*/

	private static void extract(final Frame frame){
		assert(frame.selected && !frame.genes.isEmpty()) : "A worker owns one selected, target-bearing source frame";
		//Coordinates belong to the saved canonical frame. A supplied alias may
		//have been retargeted since the metadata pass, even to equal dimensions.
		final File file=new File(frame.canonical);
		if(!file.isFile()){
			frame.unresolved("SOURCE_UNAVAILABLE");
			return;
		}
		final long bytes=file.length(), modified=file.lastModified();
		final FileFormat format=FileFormat.testInput(frame.canonical, FileFormat.FASTA, null, true, true);
		require(format.fasta() && !format.amino() && !format.interleaved(), "Source frames must be unpaired nucleotide FASTA");
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 0))){
			input.stream.start();
			for(ListNum<Read> batch=input.stream.nextList(); batch!=null; batch=input.stream.nextList()){
				for(final Read read:batch.list){
					require(read.bases!=null && read.mate==null && read.length()>0, "Source records must retain complete unpaired bases");
					frame.readBases+=read.length();
					frame.readContigs++;
					final Target target=frame.targets.get(RrnaResourceIO.first(read.id));
					if(target==null){continue;}
					target.seen++;
					if(target.seen>1){
						for(final Gene gene:target.genes){gene.unresolved("AMBIGUOUS_CONTIG_ID");}
						continue;
					}
					for(final OriginalCall call:target.originals){
						if(call.stop>read.length()){continue;}
						final boolean minus=call.component.minus;
						call.leftAvailable=minus ? read.length()-call.stop : call.start-1;
						call.rightAvailable=minus ? call.start-1 : read.length()-call.stop;
						call.leftNs=adjacentNs(read.bases, minus ? call.stop : call.start-2, minus ? 1 : -1);
						call.rightNs=adjacentNs(read.bases, minus ? call.start-2 : call.stop, minus ? -1 : 1);
					}
					for(final Gene gene:target.genes){
						gene.contigLength=read.length();
						if(gene.stop>read.length()){
							gene.unresolved("OUT_OF_RANGE_COORDINATES");
							continue;
						}
						gene.sequence=Arrays.copyOfRange(read.bases, gene.start-1, gene.stop);
						Tools.toUpperCase(gene.sequence);
						for(final byte base:gene.sequence){
							if("ACGTRYSWKMBDHVN".indexOf(base)<0){
								throw new IllegalArgumentException("Invalid DNA IUPAC byte="+(base&255)+" in "+frame.canonical+" locus="+gene.id);
							}
						}
						if(gene.minus){AminoAcid.reverseComplementBasesInPlace(gene.sequence);}
					}
				}
				input.stream.returnList(batch);
			}
		}
		if(file.length()!=bytes || file.lastModified()!=modified){frame.unresolved("FRAME_CHANGED_DURING_READ");}
		else if(frame.readBases!=frame.bases || frame.readContigs!=frame.records){frame.unresolved("FRAME_DIMENSION_MISMATCH");}
		else{for(final Gene gene:frame.genes){if(gene.target.seen==0){gene.unresolved("CONTIG_NOT_FOUND");}}}
	}

	private static int adjacentNs(final byte[] bases, int position, final int step){
		int count=0;
		while(position>=0 && position<bases.length && (bases[position]=='N' || bases[position]=='n')){
			count++;
			position+=step;
		}
		return count;
	}

	private static int context(final OriginalCall call){
		if(call.component.sequence==null || call.leftAvailable<0 || call.rightAvailable<0){return 4;}
		final int left=call.modelFrom-1, right=call.modelLength-call.modelTo;
		if(left>call.leftAvailable || right>call.rightAvailable){return 0;}
		if(left>0 && call.leftNs>0 || right>0 && call.rightNs>0){return 1;}
		return call.leftAvailable>0 && call.rightAvailable>0 ? 2 : 3;
	}

	private void selectRepresentatives(){
		final ArrayList<Gene> extracted=new ArrayList<Gene>();
		for(final Gene gene:genes){
			if(gene.sequence!=null){
				extractedLoci++;
				if(gene.nativeEligible){extracted.add(gene);}
			}
		}
		final Gene[] sorted=extracted.toArray(new Gene[0]);
		Arrays.sort(sorted, (a, b)->{
			int order=Integer.compare(a.frame.group, b.frame.group);
			if(order==0){order=Integer.compare(a.sequence.length, b.sequence.length);}
			for(int i=0; order==0 && i<a.sequence.length; i++){order=Byte.compare(a.sequence[i], b.sequence[i]);}
			return order==0 ? Integer.compare(a.id, b.id) : order;
		});
		Gene representative=null;
		for(final Gene gene:sorted){
			if(representative==null || representative.frame.group!=gene.frame.group || !Arrays.equals(representative.sequence, gene.sequence)){
				representative=gene;
				representatives++;
			}
			gene.representative=representative.id;
		}
		eligibleExtracted=sorted.length;
		countAcrossGroups(sorted);
		assert(representatives<=eligibleExtracted && eligibleExtracted<=extractedLoci)
			: "Only complete native-eligible extracted genes may vote as construction representatives";
	}

	private void countAcrossGroups(final Gene[] sorted){
		Arrays.sort(sorted, (a, b)->{
			int order=Integer.compare(a.sequence.length, b.sequence.length);
			for(int i=0; order==0 && i<a.sequence.length; i++){order=Byte.compare(a.sequence[i], b.sequence[i]);}
			if(order==0){order=Integer.compare(a.frame.group, b.frame.group);}
			return order==0 ? Integer.compare(a.id, b.id) : order;
		});
		for(int from=0; from<sorted.length;){
			int to=from, groups=0, previous=-1;
			while(to<sorted.length && Arrays.equals(sorted[from].sequence, sorted[to].sequence)){
				final int group=sorted[to].frame.group;
				if(group!=previous){
					groups++;
					previous=group;
				}
				to++;
			}
			if(groups>1){
				crossGroupClasses++;
				crossGroupEligibleLoci+=to-from;
				crossGroupRepresentatives+=groups;
			}
			from=to;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		final ByteBuilder row=new ByteBuilder();
		try(Output fasta=new Output(output("members.fa.gz")); Output ledger=new Output(output("merged_members.tsv.gz"));
			Output original=new Output(output("original_calls.tsv.gz")); Output summary=new Output(output("summary.tsv"));
			Output parameters=new Output(output("parameters.tsv"))){
			parameters.writer.println("schema\tNCRNA_ROSTER_MEMBERS_V1");
			parameters.writer.println("scope\t"+(canary ? "OPERATIONAL_CANARY_NOT_COMPLETE_ROSTER" : "FULL_ROSTER"));
			parameters.writer.println("construction_equivalence\tSAME_ACCEPTED_GROUP_AND_EQUAL_FULL_CASE_NORMALIZED_STRAND_ORIENTED_DNA_IUPAC_ARRAY");
			parameters.writer.println("source_binding\tREAD_RECORDED_CANONICAL_COORDINATE_FRAME_NOT_RETARGETABLE_SUPPLIED_ALIAS");
			parameters.writer.println("construction_source\tNATIVE_CMSEARCH_RF00174_HIT_BOUNDS_STRAND_ONLY_NO_ANNOTATION_OR_OLD_LIBRARY_MEMBERS");
			parameters.writer.println("partial_rule\tSEE_ELIGIBILITY_PARAMETERS_AND_CONFIGURABLE_MODEL_TOLERANCE");
			parameters.writer.println("cross_group_identity\tCOUNT_ELIGIBLE_EXACT_SEQUENCE_CLASSES_AND_RETAIN_EVERY_ACCEPTED_GROUP_REPRESENTATIVE");
			parameters.writer.println("old_annotation_columns\tDESCRIPTIVE_ONLY_NEVER_CONSTRUCTION_ELIGIBILITY_OR_ENDPOINTS");
			parameters.writer.println("recall\tORIGINAL_CALLED_COORDINATES_RETAINED_THROUGH_OVERLAP_COMPONENTS_AND_REPRESENTATIVES_OR_EXPLICIT_UNRESOLVED_REASON");
			ledger.writer.println("merged_locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tcontig_length\tconstruction_representative\tconstruction_state\tassembly_alias_scope\tnative_construction_eligibility");
			for(final Gene gene:genes){
				final String state=state(gene);
				if(!gene.nativeEligible){excludedComponents++;}
				if(gene.sequence==null){
					if(gene.frame.selected){unresolvedLoci++;}
					else{skippedLoci++;}
				}
				if(gene.representative==gene.id){
					final OriginalCall nativeCall=gene.nativeCall;
					assert(nativeCall!=null) : "Every construction representative must retain its own eligible native model coordinates";
					fasta.writer.print(row.clear().append(">RF00174_X_").append(gene.id).append(" assembly_group=").append(gene.frame.group)
						.append(" source_frame=").append(gene.frame.id).append(" cmfrom=").append(nativeCall.modelFrom)
						.append(" cmto=").append(nativeCall.modelTo).append(" cmlen=").append(nativeCall.modelLength)
						.append(" incomplete_model_coverage=").append(nativeCall.incomplete).nl().append(gene.sequence).nl());
				}
				ledger.writer.print(row.clear().append(gene.id).tab().append(gene.frame.group).tab().append(gene.frame.canonicalUtf8).tab()
					.append(gene.target.utf8).tab().append(gene.start).tab().append(gene.stop).tab().append(gene.minus ? '-' : '+').tab()
					.append(gene.contigLength).tab().append(gene.representative).tab().append(state).tab().append(gene.frame.aliasScope)
					.tab().append(eligibility(gene)).nl());
			}
			writeOriginal(original, row);
			require(extractedLoci+unresolvedLoci+skippedLoci==expectedLoci && (canary || skippedLoci==0), "Every merged locus needs one explicit extraction disposition");
			summary.writer.println("metric\tcount");
			metric(summary, row, "raw_RF00174_events", rawEvents);
			metric(summary, row, "native_GA_events", nativeEvents);
			metric(summary, row, "original_called_coordinate_groups", expectedOriginal);
			metric(summary, row, "merged_loci", expectedLoci);
			metric(summary, row, "extracted_loci", extractedLoci);
			metric(summary, row, "native_eligible_extracted_loci", eligibleExtracted);
			metric(summary, row, "partial_or_truncated_original_calls", partialOriginal);
			metric(summary, row, "construction_ineligible_components", excludedComponents);
			metric(summary, row, "construction_representatives", representatives);
			metric(summary, row, "removed_construction_sequence_copies", eligibleExtracted-representatives);
			metric(summary, row, "eligible_identical_sequence_classes_across_groups", crossGroupClasses);
			metric(summary, row, "eligible_loci_in_cross_group_identical_classes", crossGroupEligibleLoci);
			metric(summary, row, "retained_representatives_in_cross_group_identical_classes", crossGroupRepresentatives);
			metric(summary, row, "unresolved_merged_loci", unresolvedLoci);
			metric(summary, row, "unresolved_original_calls", unresolvedOriginal);
			metric(summary, row, "skipped_operational_canary_loci", skippedLoci);
			metric(summary, row, "source_frames_processed", sourcesProcessed);
			metric(summary, row, "source_bases_read", basesRead);
		}
	}

	private void writeOriginal(final Output output, final ByteBuilder row){
		output.writer.println(ORIGINAL_HEADER+"\tmerged_locus_id\tconstruction_representative\tconstruction_state\tnative_hit_count\tnative_truncation_bits\tnative_model_from\tnative_model_to\tnative_model_length\tincomplete_model_span\tnative_construction_eligibility\tavailable_left_bases\tavailable_right_bases\tadjacent_left_N_run\tadjacent_right_N_run\tgenomic_context\tnative_printed_bits");
		for(final OriginalCall call:originalCalls){
			final Gene gene=call.component;
			if(gene.sequence==null && gene.frame.selected){unresolvedOriginal++;}
			output.writer.print(row.clear().append(call.row).tab().append(gene.id).tab().append(gene.representative).tab().append(state(gene))
				.tab().append(call.nativeHits).tab().append(call.truncation).tab().append(call.modelFrom).tab().append(call.modelTo)
				.tab().append(call.modelLength).tab().append(call.incomplete).tab().append(call.exclusion==null ? "ELIGIBLE_NATIVE_HIT" : call.exclusion)
				.tab().append(call.leftAvailable).tab().append(call.rightAvailable).tab().append(call.leftNs).tab().append(call.rightNs)
				.tab().append(CONTEXT_NAMES[context(call)]).tab().append(call.bits, 3).nl());
		}
		assert(originalCalls.size()==expectedOriginal) : "The immutable original-call ledger must conserve recall through output";
	}

	/** Sequence context is collected in the existing FASTA pass. Assembly counts
	 * use original calls, never construction-deduplicated representatives. */
	private void writeGenomicContext(){
		final long[][] counts=new long[6][CONTEXT_NAMES.length];
		final IntHashMap internalIncomplete=new IntHashMap(), unresolvedGroups=new IntHashMap();
		final DoubleList internalBits=new DoubleList(), fullBits=new DoubleList(), incompleteBits=new DoubleList();
		for(final OriginalCall call:originalCalls){
			final int state=context(call), group=call.component.frame.group;
			if(state==4){unresolvedGroups.increment(group);}
			if(call.modelFrom==1 && call.modelTo==call.modelLength){fullBits.add(call.bits);}
			if(!call.incomplete){counts[5][state]++;continue;}
			incompleteBits.add(call.bits);
			counts[0][state]++;
			final boolean left=call.modelFrom-1==44, right=call.modelLength-call.modelTo==20;
			if(left){counts[1][state]++;}
			if(right){counts[2][state]++;}
			if(left && right){counts[3][state]++;}
			if(!left && !right){counts[4][state]++;}
			if(state==2){
				internalIncomplete.increment(group);
				internalBits.add(call.bits);
			}
		}
		int onlyIncomplete=0, onlyInternal=0;
		for(final int group:groupCalls.keys()){
			if(group<0){continue;}
			if(groupIncomplete.get(group)==groupCalls.get(group)){onlyIncomplete++;}
			if(internalIncomplete.get(group)==groupCalls.get(group)){onlyInternal++;}
		}
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(output("genomic_context.tsv")); Output assembly=new Output(output("assembly_loss.tsv"));
			Output scores=new Output(output("bit_scores.tsv")); Output rules=new Output(output("context_parameters.tsv"))){
			output.writer.println("native_call_class\tgenomic_context\toriginal_calls\tbasis");
			for(int category=0; category<counts.length; category++){
				for(int state=0; state<CONTEXT_NAMES.length; state++){
					output.writer.print(row.clear().append(CONTEXT_CLASSES[category]).tab().append(CONTEXT_NAMES[state]).tab()
						.append(counts[category][state]).append("\tSEQUENCE_PLUS_NATIVE_MODEL_METADATA\n"));
				}
			}
			assembly.writer.println("metric\tcount\tbasis");
			assembly.writer.print(row.clear().append("assembly_groups_with_calls\t").append(groupCalls.size()).append("\tMETADATA\n"));
			assembly.writer.print(row.clear().append("groups_left_zero_removing_all_incomplete\t").append(onlyIncomplete).append("\tMETADATA\n"));
			assembly.writer.print(row.clear().append("groups_left_zero_removing_internal_incomplete\t").append(onlyInternal).append("\tSEQUENCE_PLUS_METADATA_CONFIRMED\n"));
			assembly.writer.print(row.clear().append("groups_with_unresolved_sequence_context\t").append(unresolvedGroups.size()).append("\tSEQUENCE_PLUS_METADATA\n"));
			scores.writer.println("native_call_class\tcount\tminimum\tq25\tmedian\tq75\tmaximum\tbasis");
			writeScores(scores, row, "INTERNAL_INCOMPLETE", internalBits, "SEQUENCE_CONTEXT_PLUS_NATIVE_PRINTED_BITS");
			writeScores(scores, row, "EXACT_FULL_MODEL_SPAN", fullBits, "NATIVE_METADATA_PRINTED_BITS");
			writeScores(scores, row, "ALL_INCOMPLETE", incompleteBits, "NATIVE_METADATA_PRINTED_BITS");
			rules.writer.println("coordinate_units\tGENOMIC_BASES_VERSUS_MISSING_MODEL_POSITIONS_APPROXIMATE_REACH_NOT_REALIGNMENT");
			rules.writer.println("precedence\tUNRESOLVED_THEN_MISSING_MODEL_POSITIONS_EXCEED_AVAILABLE_FLANK_THEN_IMMEDIATELY_ADJACENT_N_RUN_THEN_BOTH_FLANKS_PRESENT_ELSE_OTHER");
			rules.writer.println("strand\tLEFT_RIGHT_ARE_NATIVE_MODEL_ORIENTED");
			rules.writer.println("marginals\tLEFT44_AND_RIGHT20_ROWS_MAY_OVERLAP_JOINT_ROW_IS_EXPLICIT");
			rules.writer.println("scores\tMAXIMUM_RETAINED_NATIVE_PRINTED_BIT_SCORE_PER_ORIGINAL_EXACT_COORDINATE_CALL_QUANTILES_LINEAR_ON_SORTED_VALUES");
		}
	}

	private static void writeScores(final Output output, final ByteBuilder row, final String label, final DoubleList values, final String basis){
		values.sort();
		row.clear().append(label).tab().append(values.size);
		for(final double fraction:new double[]{0, 0.25, 0.5, 0.75, 1}){
			row.tab();
			if(values.size==0){row.append("NA");}
			else{
				final double index=fraction*(values.size-1);
				final int low=(int)index, high=Math.min(values.size-1, low+1);
				row.append(values.get(low)+(index-low)*(values.get(high)-values.get(low)), 3);
			}
		}
		output.writer.print(row.tab().append(basis).nl());
	}
	private static String eligibility(final Gene gene){
		return gene.nativeEligible ? "ELIGIBLE_NATIVE_HIT_WITH_OWN_BOUNDS" : gene.nativeBoundsSeen
			? "EXCLUDED_NATIVE_PARTIAL_OR_TRUNCATED" : "NO_ELIGIBLE_NATIVE_HIT_HAS_COMPONENT_BOUNDS";
	}

	private static String state(final Gene gene){
		if(!gene.frame.selected){return "SKIPPED_OPERATIONAL_CANARY";}
		if(gene.sequence==null){
			require(gene.reason!=null, "An unresolved gene needs a concrete extraction reason");
			return gene.reason;
		}
		if(!gene.nativeEligible){return eligibility(gene);}
		return gene.representative==gene.id ? "CONSTRUCTION_REPRESENTATIVE" : "EXACT_SEQUENCE_DUPLICATE_FOR_CONSTRUCTION_ONLY";
	}
	private static void metric(final Output out, final ByteBuilder row, final String name, final long count){out.writer.print(row.clear().append(name).tab().append(count).nl());}
	private String input(final String name){return in+File.separator+name;}
	private String output(final String name){return out+File.separator+name;}
	private String text(final int column){
		p.setBounds(column);
		final String value=new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 cannot silently change a source identity");
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------       Worker Accumulation    ----------------*/
	/*--------------------------------------------------------------*/

	@Override public void accumulate(final ProcessThread worker){
		sourcesProcessed+=worker.sources;
		basesRead+=worker.bases;
		errorState|=!worker.success;
	}
	@Override public boolean success(){return !errorState;}
	@Override public ReadWriteLock rwlock(){return rwlock;}

	static final class ProcessThread extends Thread{
		ProcessThread(final Frame[] jobs_, final AtomicInteger next_){
			jobs=jobs_;
			next=next_;
		}
		@Override public void run(){
			try{
				processInner();
				success=true;
			}catch(Throwable caught){
				failure=caught;
				next.set(jobs.length);
			}
		}
		private void processInner(){
			for(int index=next.getAndIncrement(); index<jobs.length; index=next.getAndIncrement()){
				final Frame frame=jobs[index];
				extract(frame);
				sources++;
				bases+=frame.readBases;
			}
		}
		final Frame[] jobs;
		final AtomicInteger next;
		long sources, bases;
		boolean success;
		Throwable failure;
	}

	private static final class Frame{
		Frame(final int id_, final int group_, final String canonical_, final String fasta_, final long bases_, final long records_, final String scope){
			id=id_;
			group=group_;
			canonical=canonical_;
			fasta=fasta_;
			bases=bases_;
			records=records_;
			canonicalUtf8=canonical.getBytes(StandardCharsets.UTF_8);
			aliasScope=scope.getBytes(StandardCharsets.UTF_8);
		}
		void unresolved(final String reason){for(final Gene gene:genes){gene.unresolved(reason);}}
		final int id, group;
		final String canonical, fasta;
		final byte[] canonicalUtf8, aliasScope;
		final long bases, records;
		long readBases, readContigs;
		boolean selected;
		final ArrayList<Gene> genes=new ArrayList<Gene>(4);
		final ObjectMap<String, Target> targets=new ObjectMap<String, Target>(4, String.class, Target.class);
	}
	private static final class Target{
		Target(final String name_){
			name=name_;
			utf8=name.getBytes(StandardCharsets.UTF_8);
		}
		final String name;
		final byte[] utf8;
		int seen;
		final ArrayList<Gene> genes=new ArrayList<Gene>(2);
		final ArrayList<OriginalCall> originals=new ArrayList<OriginalCall>(2);
	}
	private static final class Cohort{
		Cohort(final String label_){label=label_;}
		final String label;
		final IntObjectMap<Frame> inputs=new IntObjectMap<Frame>(Frame.class);
	}
	private static final class OriginalCall{
		OriginalCall(final byte[] row_, final Gene component_, final int start_, final int stop_){
			row=row_;
			component=component_;
			start=start_;
			stop=stop_;
		}
		final byte[] row;
		final Gene component;
		final int start, stop;
		int nativeHits, truncation, modelFrom, modelTo, modelLength;
		int leftAvailable=-1, rightAvailable=-1, leftNs, rightNs;
		double bits=Double.NaN;
		boolean incomplete;
		String exclusion;
	}
	private static final class Gene{
		Gene(final int id_, final Frame frame_, final Target target_, final int start_, final int stop_, final boolean minus_){
			id=id_;
			frame=frame_;
			target=target_;
			start=start_;
			stop=stop_;
			minus=minus_;
		}
		void unresolved(final String why){
			sequence=null;
			reason=why;
		}
		final int id, start, stop;
		final Frame frame;
		final Target target;
		final boolean minus;
		byte[] sequence;
		String reason;
		int contigLength, representative=-1;
		boolean nativeBoundsSeen, nativeEligible;
		OriginalCall nativeCall;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, out, cohortManifest;
	private int expectedLoci, expectedOriginal, maxSources, modelTolerance=5;
	private double maxCoverageFraction=0.25;
	private boolean canary, metadataOnly, coverageFilter, summaryBound, errorState;
	private long truncatedOnly, coverageOnly, bothPartial, neitherPartial, coverageFailures;
	private final LongList missingLeft=new LongList(), missingRight=new LongList();
	private long nativeEvents, rawEvents, sourcesProcessed, basesRead;
	private int extractedLoci, representatives, unresolvedLoci, unresolvedOriginal, skippedLoci;
	private int eligibleExtracted, partialOriginal, excludedComponents, crossGroupClasses, crossGroupEligibleLoci, crossGroupRepresentatives;
	private final PrintStream outstream;
	private final LineParser1 p=new LineParser1('\t');
	private final ObjectMap<String, Frame> frames=new ObjectMap<String, Frame>(String.class, Frame.class);
	private final ObjectMap<String, Frame> sourceNames=new ObjectMap<String, Frame>(String.class, Frame.class);
	private final ArrayList<Frame> frameList=new ArrayList<Frame>();
	private final ArrayList<Gene> genes=new ArrayList<Gene>();
	private final ArrayList<OriginalCall> originalCalls=new ArrayList<OriginalCall>();
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final IntHashMap links=new IntHashMap();
	private final IntHashMap groupCalls=new IntHashMap(), groupIncomplete=new IntHashMap();
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	private static final String FRAME_HEADER="assembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\tcohort_input_rows\tframes_in_group\tnative_GA_exact_groups\tassembly_has_native_GA\tcoordinate_comparison_needed\tsource_available\treason";
	private static final String MERGED_HEADER="locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\texact_coordinate_groups\tnative_GA_events\tprotected_prior_witness";
	private static final String ORIGINAL_HEADER="locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tqualified\tprior_annotation_overlap\traw_events\texact_members\treciprocal80_members\tdisposition";
	private static final String CENSUS_HEADER="input_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records";
	private static final String[] CONTEXT_NAMES={"CONTIG_END_REACH", "ADJACENT_N_RUN", "INTERNAL_BOTH_FLANKS", "OTHER_END_CONTEXT", "UNRESOLVED_OR_SKIPPED"};
	private static final String[] CONTEXT_CLASSES={"ALL_INCOMPLETE", "LEFT_MISSING44_MARGINAL", "RIGHT_MISSING20_MARGINAL", "LEFT44_AND_RIGHT20", "NEITHER_TARGET_END_CLASS", "WITHIN_MODEL_COVERAGE_TOLERANCE"};
	private static final String COHORT_HEADER="cohort\tcensus\treference\tevents\twitnesses\texpected_inputs\texpected_scored_RF00174";
	private static final String EVENT_HEADER="cohort\tassembly_group\tnative_GA_qualified\tprior_annotation_overlap\tqualifying_annotation_links\tevent_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
}
