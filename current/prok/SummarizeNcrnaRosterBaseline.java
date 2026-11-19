package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

import map.IntObjectMap;
import map.LongObjectMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceMatcher.Result;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/** Joins one complete saved native baseline to every original RF00174 gene.
 * Adapted from template.A_SampleBasic2's CLI/output lifecycle and existing byte
 * prepasses. No DNA, alignments, caller, shortlist replay or truth editing.
 * Source-frame identity is exact; contig titles use their original first token.
 * Final detection uses the existing maximum-cardinality one-to-one matcher.
 * Pre-path membership instead selects the best accepted proposal whose OWN span
 * meets reciprocal overlap >=1/2. A best proposal eligible for several original
 * genes is explicitly unresolved, never credited automatically to every gene.
 * Native event cost is summed once, separately from gene/window associations.
 * Window rows allocate no per-field Strings on their common nonaccepting path.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaRosterBaseline{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, SummarizeNcrnaRosterBaseline.class, false);
		try{
			final SummarizeNcrnaRosterBaseline tool=new SummarizeNcrnaRosterBaseline(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println((tool.sampleManifest!=null ? "NCRNA_ROSTER_BASELINE_SAMPLE_SUMMARY_COMPLETE"
				: tool.acceptanceOnly ? "NCRNA_ACCEPTANCE_ONLY_SUMMARY_COMPLETE" : "NCRNA_COMPLETE_ROSTER_BASELINE_SUMMARY_COMPLETE")+" time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private SummarizeNcrnaRosterBaseline(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("population")){population=b;}
			else if(a.equals("captures")){captureManifest=b;}
			else if(a.equals("assignments")){assignments=b;}
			else if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedframes")){expectedFrames=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmembers")){expectedMembers=Parse.parseIntKMG(b);}
			else if(a.equals("acceptanceonly")){acceptanceOnly=Parse.parseBoolean(b);}
			else if(a.equals("identity")){replayIdentity=Float.parseFloat(b);}
			else if(a.equals("parity")){parity=b;}
			else if(a.equals("samplemanifest")){sampleManifest=b;}
			else{require(parser.parse(arg, a, b), "Unknown baseline summary argument: "+arg);}
		}
		capture=parser.in1;
		prefix=parser.out1;
		require((capture!=null ^ captureManifest!=null) && population!=null && assignments!=null && prefix!=null
			&& expectedGenes>0 && expectedFrames>0 && expectedMembers>0 && !parser.overwrite
			&& !parser.append && parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply exactly one of in=capture_prefix or captures=manifest.tsv, population=prepared_prefix, assignments=builder.tsv, out=fresh_prefix and full expectedgenes/frames/members");
		require(acceptanceOnly ? Float.isFinite(replayIdentity) && replayIdentity>=0 && replayIdentity<=1
			&& (replayIdentity!=.64f || parity!=null) && (parity==null || replayIdentity==.64f)
			: Float.isNaN(replayIdentity) && parity==null,
			"Acceptance-only mode needs the actual replay identity;64% requires the original accepted-gene parity table, while ordinary grading accepts neither option");
		require(sampleManifest==null || !acceptanceOnly, "An explicit frame sample cannot relax the full-population acceptance-parity control");
		for(final String suffix:new String[]{".genes.tsv.gz", ".gene_windows.tsv.gz", ".calls.tsv.gz", ".model_gc.tsv", ".summary.tsv"}){
			require(!new File(prefix+suffix).exists(), "Baseline reports require fresh outputs: "+prefix+suffix);
		}
		if(acceptanceOnly){RrnaResourceIO.fresh(prefix, ".cohorts.tsv");}
		Shared.setThreads(1);
	}

	private void process(){
		if(sampleManifest!=null){loadSample();}
		loadFrames();
		loadGenes();
		loadAssignments();
		loadCaptures();
		loadInputs();
		readWindows();
		if(acceptanceOnly){checkAcceptanceParity();}
		else{
			readCalls();
			matchFinalCalls();
		}
		writeGenesAndFit();
		if(acceptanceOnly){writeAcceptanceSummary();}
		else{
			writeModels();
			writeSummary();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Complete Population   ----------------*/
	/*--------------------------------------------------------------*/

	/** Explicit frame/path selection changes the declared denominator, never truth. */
	private void loadSample(){
		try(Input in=new Input(sampleManifest)){
			RrnaResourceIO.header(in.reader, GradeNcrnaReference.MANIFEST_HEADER, sampleManifest);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==3, "Sample manifest requires exact input_id, fasta and gff fields");
				final String name=text(0), path=text(1);
				require(!name.isEmpty() && new File(path).isAbsolute() && sampleFrames.put(name, path)==null,
					"Sample frames need unique nonempty identities and explicit absolute source paths");
			}
		}
		require(sampleFrames.size()==expectedFrames, "The explicit sample manifest must contain every expected sample frame exactly once");
	}

	private void loadFrames(){
		try(Input in=new Input(population+".frames.tsv")){
			RrnaResourceIO.header(in.reader, FRAME_HEADER, population);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==9, "Use the complete prepared source-frame ledger");
				if(sampleManifest!=null){
					final String path=sampleFrames.get(text(0));
					if(path==null){continue;}
					require(p.termEquals(path, 2), "Sample frame path must equal the original coordinate-frame path");
				}
				final Frame frame=new Frame(frames.size(), text(0), p.parseInt(1), text(2), p.parseLong(5), p.parseLong(6), p.parseInt(7));
				frame.available=p.termEquals("AVAILABLE_RECORDED_FRAME", 8);
				require(frame.available || p.termEquals("UNRESOLVED_SOURCE_UNAVAILABLE", 8), "Unknown source availability cannot become zero calls");
				require(frame.group>=0 && frame.bases>0 && frame.contigs>0 && frame.expectedGenes>=0
					&& new File(frame.path).isAbsolute() && !frameNames.contains(frame.name), "Frames require unique identities, exact paths and saved dimensions");
				frameNames.put(frame.name, frame);
				frames.add(frame);
			}
		}
		require(frames.size()==expectedFrames, "Every recorded source frame, including unavailable and zero-gene frames, must remain present");
	}

	private void loadGenes(){
		try(Input in=new Input(population+".original_calls.tsv.gz")){
			RrnaResourceIO.header(in.reader, PrepareNcrnaRosterBaseline.ORIGINAL_HEADER+"\tbaseline_input_id\tbaseline_input_state", population);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==31 && p.termEquals("true", 7), "Every original native-GA gene must remain in the recall population");
				final Frame frame=frameNames.get(text(29));
				if(frame==null && sampleManifest!=null){continue;}
				require(frame!=null && frame.group==p.parseInt(1) && frame.path.equals(text(2))
					&& p.termEquals(frame.available ? "AVAILABLE_RECORDED_FRAME" : "UNRESOLVED_SOURCE_UNAVAILABLE", 30),
					"Original genes must retain their exact frame/group and unresolved-source state");
				final String name=text(3);
				Target target=frame.targets.get(name);
				if(target==null){
					target=new Target(name);
					frame.targets.put(name, target);
				}
				final int strand=p.termEquals("-", 6) ? 1 : 0;
				require(strand==1 || p.termEquals("+", 6), "Each original gene needs an explicit strand");
				final Gene gene=new Gene(genes.size(), p.parseInt(0), frame, target,
					p.parseInt(4)-1, p.parseInt(5)-1, strand, p.parseInt(14), p.parseInt(17), extractionState());
				if(acceptanceOnly){
					final boolean incomplete=p.termEquals("true", 21);
					require(incomplete || p.termEquals("false", 21), "Original incomplete-coverage flags must be explicitly true or false before population partitioning");
					gene.cohorts=1|(p.parseInt(18)==1 && p.parseInt(19)==p.parseInt(20) ? 2 : 0)
						|(incomplete ? 4 : 32)|(p.parseInt(18)-1==44 && p.parseInt(20)-p.parseInt(19)==20 ? 8 : 0)
						|(gene.truncation!=0 ? 16 : 0);
				}
				require(gene.id>=0 && !geneIds.contains(gene.id) && gene.start>=0 && gene.stop>=gene.start,
					"Original source coordinates and gene IDs must remain valid and unique");
				geneIds.put(gene.id, gene);
				genes.add(gene);
				target.genes.add(gene);
				frame.genes++;
			}
		}
		require(genes.size()==expectedGenes, "No construction exclusion or missing window may reduce original-gene recall");
		for(final Frame frame:frames){require(frame.genes==frame.expectedGenes, "Prepared per-frame gene counts must conserve the original population");}
	}

	private void loadAssignments(){
		try(Input in=new Input(assignments)){
			RrnaResourceIO.header(in.reader, "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus", assignments);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==6 && p.termStartsWith("RF00174_X_", 0), "Use the actual new roster builder assignments");
				p.setBounds(0);
				int end=p.a()+10;
				while(end<p.b() && p.line()[end]!=' '){end++;}
				final int representative=Parse.parseInt(p.line(), p.a()+10, end);
				require(representative>=0 && !construction.contains(representative), "One construction representative needs one builder disposition");
				final Model model;
				if(p.termEquals("kept", 5)){model=model(RrnaResourceIO.first(text(4)), true);}
				else{
					require(p.termEquals("", 4) && (p.termEquals("small", 5) || p.termEquals("orphan", 5)), "Nonvoting members cannot invent a construction model");
					model=UNMODELED;
				}
				construction.put(representative, model);
			}
		}
		require(construction.size()==expectedMembers && !models.isEmpty(), "Use all construction assignments and their actual retained model names");
		for(final Gene gene:genes){
			if(gene.representative<0){continue;}
			gene.constructed=construction.get(gene.representative);
			require(gene.constructed!=null, "Each original construction representative must bind its actual builder assignment");
		}
	}

	/** The completed native input table is the exact slot namespace for both
	 * event and final-call files; availability is never inferred from absence. */
	private void loadCaptures(){
		if(captureManifest==null){
			captures.add(new Capture("", capture, null));
			return;
		}
		final ObjectMap<String, Capture> names=new ObjectMap<String, Capture>(String.class, Capture.class);
		try(Input in=new Input(captureManifest)){
			RrnaResourceIO.header(in.reader, "capture_id\tcapture_prefix\treference_manifest", captureManifest);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==3, "A capture manifest binds each distinct namespace to its prefix and exact input list");
				final String name=text(0), path=text(1), reference=text(2);
				require(name.matches("[A-Za-z0-9_]+") && new File(path).isAbsolute()
					&& new File(reference).isAbsolute() && !names.contains(name),
					"Capture namespaces must be unique simple names and their saved paths absolute");
				final Capture part=new Capture(name, path, reference);
				names.put(name, part);
				captures.add(part);
			}
		}
		require(!captures.isEmpty(), "A sharded baseline needs at least one completed capture");
	}

	/** Check each shard against its actual ordered reference list as well as the
	 * global frame ledger. Empty inputs remain necessary observations. */
	private void loadInputs(){
		for(final Capture part:captures){loadInputs(part);}
		for(final Frame frame:frames){require(frame.seen==frame.available, "Every available source must complete, including zero-window inputs; unavailable sources stay explicit");}
	}

	private void loadInputs(final Capture part){
		final LineParser1 referenceParser=new LineParser1('\t');
		try(Input in=new Input(part.path+".inputs.tsv"); Input reference=part.reference==null ? null : new Input(part.reference)){
			RrnaResourceIO.header(in.reader, acceptanceOnly ? ReplayNcrnaWindows.INPUT_HEADER : INPUT_HEADER, part.path);
			if(reference!=null){RrnaResourceIO.header(reference.reader, GradeNcrnaReference.MANIFEST_HEADER, part.reference);}
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==(acceptanceOnly ? 9 : 10) && p.parseInt(0)==part.slots.size(), "Completed native input slots must be contiguous within each capture");
				final Frame frame=frameNames.get(text(1));
				require(frame!=null && frame.available && !frame.seen && field(2, frame.pathField)
					&& p.parseLong(3)==frame.contigs && p.parseLong(4)==frame.bases,
					"Every completed input must bind the declared exact source and saved contig/base census");
				if(reference!=null){
					final byte[] expected=reference.reader.nextLine();
					require(expected!=null, "A capture cannot add an input outside its declared shard");
					referenceParser.set(expected);
					require(referenceParser.terms()==3 && referenceParser.termEquals(frame.name, 0)
						&& referenceParser.termEquals(frame.path, 1), "A capture must retain its shard's exact ordered input bindings");
				}
				frame.seen=true;
				frame.capture=part;
				for(int i=0; i<(acceptanceOnly ? 4 : 5); i++){
					frame.expected[i]=p.parseLong(i+5);
					require(frame.expected[i]>=0, "Native counts cannot be negative");
				}
				require(frame.expected[0]<=frame.expected[1] && frame.expected[3]<=frame.expected[1], "Captured and accepted events cannot exceed scheduled windows");
				part.slots.add(frame);
				completedInputs++;
			}
			if(reference!=null){require(reference.reader.nextLine()==null, "A completed capture cannot omit even a zero-call shard input");}
		}
		if(acceptanceOnly){
			boolean identity=false;
			try(Input in=new Input(part.path+".replay.tsv")){
				RrnaResourceIO.header(in.reader, "metric\tvalue", part.path);
				for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
					p.set(line);
					require(p.terms()==2, "Each replay retains resolved scalar controls and counted work");
					if(p.termEquals("identity", 0)){
						require(!identity && p.parseFloat(1)==replayIdentity, "The grading identity must equal each replay's actual configured scalar gate");
						identity=true;
					}
				}
			}
			require(identity, "Every replay prefix must declare its actual identity, not inherit a filename label");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Native Observations  ----------------*/
	/*--------------------------------------------------------------*/

	private void readWindows(){
		try(Output associations=new Output(prefix+".gene_windows.tsv.gz")){
			associations.writer.println(ASSOCIATION_HEADER);
			for(final Capture part:captures){readWindows(part, associations);}
		}
	}

	private void readWindows(final Capture part, final Output associations){
		slots=part.slots;
		try(Input in=new Input(part.path+".results.tsv")){
			RrnaResourceIO.header(in.reader, acceptanceOnly ? ReplayNcrnaWindows.RESULT_HEADER : CaptureNcrnaWindows.RESULT_HEADER, part.path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==(acceptanceOnly ? 23 : 26), "Use the complete declared ordinary-capture or acceptance-replay schema");
				final Frame frame=slot(0, 1, 2);
				final long ordinal=p.parseLong(3);
				final int strand=p.parseInt(5), length=p.parseInt(6), start=p.parseInt(10), stop=p.parseInt(11);
				require(ordinal>=0 && ordinal<frame.contigs && (strand==0 || strand==1)
					&& length>0 && start>=0 && stop>=start && stop<length, "Scheduled windows require actual forward-contig coordinates");
				final int shift=acceptanceOnly ? 3 : 0;
				final long work=p.parseLong(17-shift);
				require(work>=0 && (acceptanceOnly || p.parseLong(14)>=0 && p.parseLong(15)>=0
					&& p.parseLong(16)>=0 && work==p.parseLong(14)+p.parseLong(15)+p.parseLong(16)),
					"Every event retains its complete actual alignment work exactly once; replay does not invent phase counters");
				final boolean captured=!p.termEquals("NA", 12), accepted=p.termEquals("ACCEPTED_PROPOSAL", 13);
				final int gc=p.parseInt(24-shift), acgt=p.parseInt(25-shift);
				if(captured){
					event.clear().append('i').append(p.parseInt(0)).append('e').append(frame.observed[0]);
					require(field(12, event.array, event.length()), "Captured event ordinals must occur once in each native input's execution order");
					frame.observed[0]++;
					if(acceptanceOnly){savedWindowBases+=stop-(long)start+1;}
					require(gc>=0 && acgt>=gc && acgt<=stop-start+1, "GC observations use known ACGT bases within the actual window");
				}else{require(!accepted && work==0 && gc==-1 && acgt==-1 && p.termEquals("UNCAPTURED_GATE_REJECTION", 13), "Only unaccepted zero-work gate rejections lack captured bases");}
				frame.observed[1]++;
				frame.observed[2]+=work;
				if(accepted){frame.observed[3]++;}
				else{require(p.termEquals(captured ? "REJECTED_VERIFICATION" : "UNCAPTURED_GATE_REJECTION", 13), "Unknown window status cannot become a verifier miss");}
				final int a=accepted ? p.parseInt(22-shift) : -1, b=accepted ? p.parseInt(23-shift) : -1;
				final float identity=accepted ? p.parseFloat(19-shift) : -1;
				require(!accepted || a>=0 && b>=a && b<length && Float.isFinite(identity) && identity>=0 && identity<=1,
					"Accepted proposals require their own measured identity and genomic span");
				require(ordinal>=frame.lastWindowOrdinal, "Window rows for one input must follow its actual contig execution order");
				frame.lastWindowOrdinal=ordinal;
				final Target target=binding(frame, ordinal, 4).target;
				if(target==null){
					if(acceptanceOnly && accepted){noRosterOverlap++;noSameStrandRosterOverlap++;}
					continue;
				}
				int matches=0, contained=0;
				boolean anyRosterOverlap=false, sameStrandRosterOverlap=false;
				Gene unique=null;
				for(final Gene gene:target.genes){
					if(acceptanceOnly && accepted && overlap(gene.start, gene.stop, a, b)>0){
						anyRosterOverlap=true;
						if(gene.strand==strand){sameStrandRosterOverlap=true;}
					}
					if(gene.strand!=strand){continue;}
					final boolean overlaps=overlap(gene.start, gene.stop, start, stop)>0;
					final boolean contains=start<=gene.start && stop>=gene.stop;
					final boolean qualifies=gene.resolved && accepted && qualifies(gene.start, gene.stop, a, b);
					if(overlaps){
						gene.windows++;
						gene.overlapWork+=work;
					}
					if(contains){
						gene.contained++;
						if(gene.resolved){contained++;}
					}
					if(qualifies){
						matches++;
						unique=gene;
					}
					if(!overlaps && !qualifies){continue;}
					// This cost is attached to an association, never added to the
					// global total again when another gene shares the same event.
					row.clear().append(gene.id).tab().append(frame.nameField).tab()
						.append(part.qualifier).append('i').append(frame.slot).append('s').append(frame.observed[1]-1).tab();
					if(captured){row.append(part.qualifier).append('i').append(frame.slot).append('e').append(frame.observed[0]-1);}
					else{row.append("NA");}
					associations.writer.print(row.tab().append(overlaps).tab().append(contains).tab().append(qualifies)
						.tab().append(work).tab().append(gene.resolved).nl());
					geneWindowAssociations++;
				}
				if(contained>0 && captured){trueWindowEvents++;}
				if(!accepted){continue;}
				if(acceptanceOnly){
					if(!anyRosterOverlap){noRosterOverlap++;}
					if(!sameStrandRosterOverlap){noSameStrandRosterOverlap++;}
				}
				if(matches>1){
					ambiguousEvents++;
					for(final Gene gene:target.genes){if(gene.resolved && gene.strand==strand && qualifies(gene.start, gene.stop, a, b)){gene.ambiguous=Math.max(gene.ambiguous, identity);}}
				}else if(matches==1){
					final Model model=model(RrnaResourceIO.first(decoded(18-shift)), false);
					final long id=frame.observed[0]-1;
					if(better(unique, identity, model, start, stop, id)){
						unique.best=model;
						unique.identity=identity;
						unique.event=id;
						unique.windowStart=start;
						unique.windowStop=stop;
						unique.gc=gc;
						unique.acgt=acgt;
						unique.selectedWork=work;
					}
				}
			}
		}
		for(final Frame frame:slots){for(int i=0; i<4; i++){require(frame.expected[i]==frame.observed[i], "Native input totals must equal all captured/scheduled/work/acceptance event rows");}}
	}

	private void readCalls(){
		for(final Capture part:captures){readCalls(part);}
	}

	private void readCalls(final Capture part){
		slots=part.slots;
		try(Input in=new Input(part.path+".final_calls.tsv")){
			RrnaResourceIO.header(in.reader, CaptureNcrnaWindows.FINAL_HEADER, part.path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==12, "Final path calls require the original input identity, not bare combined-GFF contig names");
				final Frame frame=slot(0, 1, 2);
				final long ordinal=p.parseLong(3);
				final int start=p.parseInt(5), stop=p.parseInt(6), strand=p.parseInt(7);
				final float identity=p.parseFloat(9);
				require(ordinal>=0 && ordinal<frame.contigs && start>=0 && stop>=start && stop<Integer.MAX_VALUE
					&& (strand==0 || strand==1) && Float.isFinite(identity) && identity>=0 && identity<=1,
					"Final path calls require explicit original-forward coordinates and measured identity");
				final String contig=binding(frame, ordinal, 4).name;
				final Model model=model(RrnaResourceIO.first(decoded(8)), false);
				calls.add(new Call(frame, new Locus(contig, FAMILY, start+1, stop+1, strand==0 ? '+' : '-'), model, identity));
				frame.observed[4]++;
			}
		}
		for(final Frame frame:slots){require(frame.expected[4]==frame.observed[4], "Every final target call must occur once in the input-bound final-call table");}
	}

	/** Sorting retained calls makes equivalent maximum-cardinality choices stable
	 * across input-worker emission order; it does not require replaying a caller. */
	private void matchFinalCalls(){
		Collections.sort(calls, CALL_ORDER);
		final ArrayList<Locus> query=new ArrayList<Locus>(calls.size()), reference=new ArrayList<Locus>(genes.size());
		final ArrayList<Gene> matchable=new ArrayList<Gene>();
		final int[] qi=new int[calls.size()];
		for(int i=0; i<calls.size(); i++){
			query.add(calls.get(i).locus);
			qi[i]=calls.get(i).frame.ordinal+1;
		}
		for(final Gene gene:genes){
			// Unresolved genes stay in the final denominator with no possible
			// matching credit. They are not masks or invented coordinate aliases.
			if(!gene.resolved){continue;}
			matchable.add(gene);
			reference.add(new Locus(gene.target.name, FAMILY, gene.start+1, gene.stop+1, gene.strand==0 ? '+' : '-'));
		}
		final int[] ri=new int[matchable.size()];
		for(int i=0; i<ri.length; i++){ri[i]=matchable.get(i).frame.ordinal+1;}
		unresolvedReferences=genes.size()-matchable.size();
		finalCallByGene=new int[genes.size()];
		Arrays.fill(finalCallByGene, -1);
		matched=NcrnaReferenceMatcher.match(query, qi, reference, ri, new boolean[ri.length], .5f);
		try(Output out=new Output(prefix+".calls.tsv.gz")){
			out.writer.println(CALLS_HEADER);
			for(int i=0; i<calls.size(); i++){
				final Call call=calls.get(i);
				final int referenceIndex=matched.callToReference[i];
				final Gene gene=referenceIndex<0 ? null : matchable.get(referenceIndex);
				if(gene!=null){finalCallByGene[gene.ordinal]=i;}
				row.clear().append(i).tab().append(call.frame.nameField).tab();
				BuildNcrnaReference.field(row, call.locus.seqid).tab().append(call.locus.start).tab().append(call.locus.stop).tab()
					.append(call.locus.strand).tab().append(call.model.field).tab().append(call.identity, 9).tab();
				if(gene!=null){row.append(gene.id);}else{row.append("NA");}
				out.writer.print(row.tab().append(matched.duplicateCalls[i]).tab().append(gene!=null ? "MATCHED" : "UNMATCHED_NOT_ADJUDICATED_FP").nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Reports          ----------------*/
	/*--------------------------------------------------------------*/

	/** Compare the actual accepted original-gene set, not aggregate counts or
	 * unrelated final calls. Every baseline row must bind one unchanged original. */
	private void checkAcceptanceParity(){
		assert(acceptanceOnly) : "Accepted-gene parity is separate from final path matching";
		if(parity==null){return;}
		try(Input in=new Input(parity)){
			RrnaResourceIO.header(in.reader, GENES_HEADER, parity);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==26, "The64% control must use the complete original ordinary-baseline gene table");
				final Gene gene=geneIds.get(p.parseInt(0));
				require(gene!=null && !gene.paritySeen && field(1, gene.frame.nameField)
					&& p.parseInt(3)==gene.start+1 && p.parseInt(4)==gene.stop+1
					&& p.termEquals(gene.strand==0 ? "+" : "-", 5),
					"Every original parity row must occur once with unchanged input, strand and source coordinates");
				BuildNcrnaReference.field(row.clear(), gene.target.name);
				require(field(2, row.array, row.length()), "The original parity table cannot substitute a different contig");
				boolean stateKnown=false;
				for(final String state:PROPOSAL_STATES){if(p.termEquals(state, 11)){stateKnown=true;break;}}
				require(stateKnown, "Unknown proposal states cannot silently become a parity miss");
				final boolean original=p.termEquals("UNIQUE_BEST_ACCEPTANCE", 11);
				final boolean replay=gene.resolved && gene.best!=null && gene.ambiguous<gene.identity;
				if(original!=replay){throw new IllegalStateException("The64% acceptance replay changed original gene "+gene.id
					+": baseline accepted="+original+", replay accepted="+replay);}
				gene.paritySeen=true;
				parityChecked++;
			}
		}
		require(parityChecked==genes.size(), "The64% parity check must retain every original gene, including baseline misses");
	}

	private void writeAcceptanceSummary(){
		assert(acceptanceOnly && cohortTotals[0]==genes.size() && cohortAccepted[0]==uniqueMembers)
			: "Acceptance cohorts must conserve the complete original denominator and attributed accepted genes";
		require(cohortTotals[2]+cohortTotals[5]==genes.size()
			&& cohortAccepted[2]+cohortAccepted[5]==uniqueMembers,
			"Within-coverage-tolerance and incomplete coverage must partition all original genes and their attributed accepted subset");
		long captured=0, scheduled=0, alignments=0, accepts=0, inputBases=0;
		for(final Frame frame:frames){
			captured+=frame.observed[0];
			scheduled+=frame.observed[1];
			alignments+=frame.observed[2];
			accepts+=frame.observed[3];
			if(frame.seen){inputBases+=frame.bases;}
		}
		try(Output out=new Output(prefix+".cohorts.tsv")){
			out.writer.println("cohort\toriginal_genes\taccepted_original_genes");
			for(int i=0; i<COHORTS.length; i++){
				require(cohortAccepted[i]<=cohortTotals[i], "Accepted subsets cannot exceed their original population");
				out.writer.print(row.clear().append(COHORTS[i]).tab().append(cohortTotals[i]).tab().append(cohortAccepted[i]).nl());
			}
		}
		try(Output out=new Output(prefix+".summary.tsv")){
			row.clear().append("metric\tvalue\nmode\tACCEPTANCE_ONLY_NO_FINAL_PATH\nidentity\t").append(Float.toString(replayIdentity))
				.append("\noriginal_genes\t").append(genes.size()).append("\nrecorded_frames\t").append(frames.size())
				.append("\ncompleted_inputs\t").append(completedInputs).append("\nunique_best_members\t").append(uniqueMembers)
				.append("\nscheduled_window_events\t").append(scheduled).append("\ncaptured_window_events\t").append(captured)
				.append("\nsaved_window_bases\t").append(savedWindowBases).append("\nrecorded_input_bases\t").append(inputBases)
				.append("\nactual_alignment_events\t").append(alignments).append("\naccepted_proposal_events\t").append(accepts)
				.append("\naccepted_events_without_any_strand_roster_overlap_unadjudicated\t").append(noRosterOverlap)
				.append("\naccepted_events_without_same_strand_roster_overlap_unadjudicated\t").append(noSameStrandRosterOverlap)
				.append("\nno_roster_overlap_events_per_saved_window\t");
			if(captured>0){row.append(Double.toString(noRosterOverlap/(double)captured));}else{row.append("NA");}
			row.append("\nno_roster_overlap_events_per_million_saved_window_bases\t");
			if(savedWindowBases>0){row.append(Double.toString(noRosterOverlap*1e6/savedWindowBases));}else{row.append("NA");}
			row.append("\nno_roster_overlap_events_per_million_recorded_input_bases\t");
			if(inputBases>0){row.append(Double.toString(noRosterOverlap*1e6/inputBases));}else{row.append("NA");}
			out.writer.print(row.append("\nambiguous_proposal_events\t").append(ambiguousEvents)
				.append("\ngenes_with_ambiguous_best_acceptance\t").append(ambiguousGenes)
				.append("\nparity_original_genes_checked\t").append(parityChecked).nl());
		}
	}

	private void writeGenesAndFit(){
		try(Output out=new Output(prefix+".genes.tsv.gz")){
			out.writer.println(acceptanceOnly ? ACCEPTANCE_GENES_HEADER : GENES_HEADER);
			for(final Gene gene:genes){
				final boolean usable=gene.resolved && gene.best!=null && gene.ambiguous<gene.identity;
				final String state=!gene.frame.available ? "UNRESOLVED_INPUT" : !gene.resolved ? "UNRESOLVED_EXTRACTION" : usable ? "UNIQUE_BEST_ACCEPTANCE"
					: gene.ambiguous>=0 ? "AMBIGUOUS_BEST_ACCEPTANCE" : gene.windows==0 ? "NO_SCHEDULED_WINDOW"
					: gene.contained==0 ? "PARTIAL_WINDOW_NO_ACCEPTANCE" : "WINDOW_WITHOUT_ACCEPTED_MODEL";
				if(usable){
					uniqueMembers++;
					if(gene.constructed!=null && gene.constructed!=UNMODELED){
						comparable++;
						if(gene.best!=gene.constructed){changed++;}
					}
					if(!acceptanceOnly){if(gene.acgt>0){gene.best.add(gene.gc/(double)gene.acgt);}else{gene.best.missing++;}}
				}else if(gene.ambiguous>=0){ambiguousGenes++;}
				if(!gene.frame.available){missingGenes++;}
				if(gene.frame.available && gene.windows==0){noWindowGenes++;}
				if(acceptanceOnly){
					for(int i=0; i<COHORTS.length; i++){
						if((gene.cohorts&(1<<i))==0){continue;}
						cohortTotals[i]++;
						if(usable){cohortAccepted[i]++;}
					}
				}
				row.clear().append(gene.id).tab().append(gene.frame.nameField).tab();
				BuildNcrnaReference.field(row, gene.target.name).tab().append(gene.start+1).tab().append(gene.stop+1).tab()
					.append(gene.strand==0 ? '+' : '-').tab().append(gene.truncation).tab().append(gene.representative).tab()
					.append(gene.constructed==null ? UNMODELED.field : gene.constructed.field).tab().append(gene.windows).tab()
					.append(gene.contained).tab().append(state).tab().append(gene.identity, 9).tab().append(gene.ambiguous, 9).tab();
				if(usable){row.append(gene.best.field).tab().append(gene.frame.capture.qualifier).append('i').append(gene.frame.slot).append('e').append(gene.event)
					.tab().append(gene.windowStart).tab().append(gene.windowStop).tab().append(gene.gc).tab().append(gene.acgt);}
				else{row.append("NA\tNA\tNA\tNA\tNA\tNA");}
				row.tab();
				if(usable && gene.constructed!=null && gene.constructed!=UNMODELED){row.append(gene.best!=gene.constructed);}else{row.append("NA");}
				if(!acceptanceOnly){
					final int ci=finalCallByGene[gene.ordinal];
					row.tab().append(ci>=0).tab();
					if(ci>=0){row.append(ci);}else{row.append("NA");}
				}
				row.tab().append(gene.overlapWork).tab();
				if(usable){row.append(gene.selectedWork);}else{row.append("NA");}
				out.writer.print(row.tab().append(EXTRACTION_STATES[gene.extraction]).nl());
			}
		}
	}

	private void writeModels(){
		long fitted=0;
		try(Output out=new Output(prefix+".model_gc.tsv")){
			out.writer.println("model\tknown_gc_members\tunknown_gc_members\tmean\tpopulation_sd\tminimum\tmaximum\tstate");
			for(final Model model:modelOrder){
				fitted+=model.n+model.missing;
				row.clear().append(model.field).tab().append(model.n).tab().append(model.missing).tab();
				if(model.n==0){row.append("NA\tNA\tNA\tNA\tNO_KNOWN_GC_MEMBERS");}
				else{
					// Rounded bounds can exclude their own input, particularly a ZERO_SD model.
					row.append(Double.toString(model.mean)).tab()
						.append(Double.toString(Math.sqrt(Math.max(0, model.m2/model.n)))).tab()
						.append(Double.toString(model.min)).tab().append(Double.toString(model.max)).tab()
						.append(model.m2==0 ? "ZERO_SD" : "FIT");
				}
				out.writer.print(row.nl());
			}
		}
		require(fitted==uniqueMembers, "Known and missing GC members must conserve every uniquely assigned original gene");
	}

	private void writeSummary(){
		assert(matched.matched+matched.falseNegatives+unresolvedReferences==genes.size())
			: "Resolved matches/misses plus explicitly unresolved genes must conserve the whole original denominator";
		long windows=0, captured=0, alignments=0, accepts=0;
		for(final Frame frame:frames){
			captured+=frame.observed[0];
			windows+=frame.observed[1];
			alignments+=frame.observed[2];
			accepts+=frame.observed[3];
		}
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.print(row.clear().append("metric\tcount\noriginal_genes\t").append(genes.size())
				.append("\nrecorded_frames\t").append(frames.size()).append("\ncompleted_inputs\t").append(completedInputs)
				.append("\nmatchable_original_genes\t").append(genes.size()-unresolvedReferences)
				.append("\nunresolved_original_genes\t").append(unresolvedReferences)
				.append("\nunresolved_input_genes\t").append(missingGenes).append("\navailable_input_genes_without_scheduled_overlap\t").append(noWindowGenes)
				.append("\nscheduled_window_events\t").append(windows).append("\ncaptured_window_events\t").append(captured)
				.append("\nactual_alignment_events\t").append(alignments).append("\naccepted_proposal_events\t").append(accepts)
				.append("\ngene_window_associations_nonadditive_cost\t").append(geneWindowAssociations)
				.append("\ncaptured_windows_containing_original_gene\t").append(trueWindowEvents)
				.append("\nambiguous_proposal_events\t").append(ambiguousEvents).append("\ngenes_with_ambiguous_best_acceptance\t").append(ambiguousGenes)
				.append("\nunique_best_members\t").append(uniqueMembers).append("\nconstruction_comparable_genes\t").append(comparable)
				.append("\nconstruction_model_changes\t").append(changed).append("\nfinal_call_rows\t").append(calls.size())
				.append("\nfinal_one_to_one_detected_genes\t").append(matched.matched).append("\nfinal_missed_genes\t").append(genes.size()-matched.matched)
				.append("\nunmatched_final_calls_not_adjudicated_fp\t").append(matched.falsePositives).append("\nexact_duplicate_final_calls\t").append(matched.duplicates).nl());
			if(sampleManifest!=null){
				out.writer.print(row.clear().append("scope\tEXPLICIT_INPUT_SAMPLE_NOT_FULL_POPULATION\nsample_manifest\t")
					.append(sampleManifest).nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	private Frame slot(final int number, final int name, final int path){
		final int index=p.parseInt(number);
		require(index>=0 && index<slots.size(), "Native rows must use their completed input-slot namespace");
		final Frame frame=slots.get(index);
		require(field(name, frame.nameField) && field(path, frame.pathField), "A native row cannot cross original input identity");
		frame.slot=index;
		return frame;
	}

	/** Names and ordinals form a bijection over observed contigs in one source.
	 * The same binding is checked during both the event and final-call passes. */
	private Contig binding(final Frame frame, final long ordinal, final int column){
		Contig bound=frame.lastContig;
		if(bound==null || bound.ordinal!=ordinal){bound=frame.ordinals.get(ordinal);}
		if(bound!=null){
			// Native ASCII separators take the allocation-free path. Use the
			// existing Unicode first-token rule for unusual whitespace/encodings.
			require(firstField(column, bound.field) || RrnaResourceIO.first(decoded(column)).equals(bound.name),
				"One source-frame contig ordinal cannot name another contig token");
		}else{
			final String name=RrnaResourceIO.first(decoded(column));
			require(!name.isEmpty(), "Observed contigs require a nonempty original first-token identity");
			final Contig old=frame.names.get(name);
			require(old==null || old.ordinal==ordinal, "One source-frame contig token cannot occur at two different ordinals");
			bound=new Contig(name, ordinal, frame.targets.get(name));
			frame.ordinals.put(ordinal, bound);
			frame.names.put(name, bound);
		}
		frame.lastContig=bound;
		return bound;
	}

	private int extractionState(){
		for(int i=0; i<EXTRACTION_STATES.length; i++){if(p.termEquals(EXTRACTION_STATES[i], 15)){return i;}}
		throw new IllegalArgumentException("Use the accepted full-extraction ledger; unknown or operationally skipped extraction states cannot become resolved truth");
	}

	private Model model(final String name, final boolean create){
		Model model=models.get(name);
		if(model==null && create){
			model=new Model(name);
			models.put(name, model);
			modelOrder.add(model);
		}
		require(model!=null && name.startsWith("cobalamin_RF00174_"), "Accepted and construction models must belong to this actual new RF00174 library");
		return model;
	}

	private static boolean better(final Gene gene, final float id, final Model model, final int start, final int stop, final long event){
		if(gene.best==null || id!=gene.identity){return gene.best==null || id>gene.identity;}
		final int order=model.name.compareTo(gene.best.name);
		return order<0 || order==0 && (start<gene.windowStart || start==gene.windowStart && (stop<gene.windowStop || stop==gene.windowStop && event<gene.event));
	}

	private static int overlap(final int a, final int b, final int c, final int d){return Math.max(0, Math.min(b, d)-Math.max(a, c)+1);}
	private static boolean qualifies(final int a, final int b, final int c, final int d){
		final long shared=overlap(a, b, c, d);
		return 2*shared>=b-(long)a+1 && 2*shared>=d-(long)c+1;
	}
	private boolean field(final int column, final byte[] expected){return field(column, expected, expected.length);}
	private boolean field(final int column, final byte[] expected, final int length){
		p.setBounds(column);
		if(p.b()-p.a()!=length){return false;}
		for(int i=0; i<length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		return true;
	}
	private boolean firstField(final int column, final byte[] expected){
		p.setBounds(column);
		if(p.b()-p.a()<expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(p.line()[p.a()+i]!=expected[i]){return false;}}
		final int end=p.a()+expected.length;
		return end==p.b() || p.line()[end]==' ' || end+2<p.b() && p.line()[end]=='%' && p.line()[end+1]=='0' && p.line()[end+2]=='9';
	}
	private String text(final int column){
		final String value=p.parseString(column);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 cannot change a source or model identity");
		return value;
	}
	private String decoded(final int column){
		p.setBounds(column);
		return decoder.decode(p.line(), p.a(), p.b());
	}
	private static byte[] escaped(final String text){return BuildNcrnaReference.field(new ByteBuilder(), text).toBytes();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------              Types           ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Frame{
		Frame(int ordinal_, String name_, int group_, String path_, long bases_, long contigs_, int genes_){
			ordinal=ordinal_;
			name=name_;
			group=group_;
			path=path_;
			bases=bases_;
			contigs=contigs_;
			expectedGenes=genes_;
			nameField=escaped(name);
			pathField=escaped(path);
		}
		final int ordinal, group, expectedGenes;
		final String name, path;
		final long bases, contigs;
		final byte[] nameField, pathField;
		final long[] expected=new long[5], observed=new long[5];
		final ObjectMap<String, Target> targets=new ObjectMap<String, Target>(4, String.class, Target.class);
		final ObjectMap<String, Contig> names=new ObjectMap<String, Contig>(4, String.class, Contig.class);
		final LongObjectMap<Contig> ordinals=new LongObjectMap<Contig>(4, Contig.class);
		int genes, slot=-1;
		boolean available, seen;
		long lastWindowOrdinal=-1;
		Contig lastContig;
		Capture capture;
	}

	/** Local native slots never escape their capture namespace. The qualifier is
	 * retained on event references so later saved-window consumers can locate DNA. */
	private static final class Capture{
		Capture(String name, String path_, String reference_){
			qualifier=name.isEmpty() ? "" : name+":";
			path=path_;
			reference=reference_;
		}
		final String qualifier, path, reference;
		final ArrayList<Frame> slots=new ArrayList<Frame>();
	}
	private static final class Contig{
		Contig(String name_, long ordinal_, Target target_){
			name=name_;
			ordinal=ordinal_;
			target=target_;
			field=escaped(name);
		}
		final String name;
		final long ordinal;
		final Target target;
		final byte[] field;
	}
	private static final class Target{
		Target(String name_){
			name=name_;
		}
		final String name;
		final ArrayList<Gene> genes=new ArrayList<Gene>();
	}
	private static final class Gene{
		Gene(int ordinal_, int id_, Frame frame_, Target target_, int start_, int stop_, int strand_, int representative_, int truncation_, int extraction_){
			ordinal=ordinal_;
			id=id_;
			frame=frame_;
			target=target_;
			start=start_;
			stop=stop_;
			strand=strand_;
			representative=representative_;
			truncation=truncation_;
			extraction=extraction_;
			resolved=frame.available && extraction<4;
		}
		final int ordinal, id, start, stop, strand, representative, truncation, extraction;
		final boolean resolved;
		final Frame frame;
		final Target target;
		Model constructed, best;
		float identity=-1, ambiguous=-1;
		long windows, contained, overlapWork, selectedWork, event=-1;
		int windowStart, windowStop, gc, acgt, cohorts;
		boolean paritySeen;
	}
	private static final class Model{
		Model(String name_){
			name=name_;
			field=escaped(name);
		}
		void add(final double gc){
			assert(gc>=0 && gc<=1) : "Known-base GC fractions must remain within their physical range";
			n++;
			final double delta=gc-mean;
			mean+=delta/n;
			m2+=delta*(gc-mean);
			min=Math.min(min, gc);
			max=Math.max(max, gc);
		}
		final String name;
		final byte[] field;
		long n, missing;
		double mean, m2, min=Double.POSITIVE_INFINITY, max=Double.NEGATIVE_INFINITY;
	}
	private static final class Call{
		Call(Frame frame_, Locus locus_, Model model_, float identity_){
			frame=frame_;
			locus=locus_;
			model=model_;
			identity=identity_;
		}
		final Frame frame;
		final Locus locus;
		final Model model;
		final float identity;
	}
	private static final Comparator<Call> CALL_ORDER=(a, b)->{
		int order=Integer.compare(a.frame.ordinal, b.frame.ordinal);
		if(order==0){order=a.locus.seqid.compareTo(b.locus.seqid);}
		if(order==0){order=Character.compare(a.locus.strand, b.locus.strand);}
		if(order==0){order=Integer.compare(a.locus.start, b.locus.start);}
		if(order==0){order=Integer.compare(a.locus.stop, b.locus.stop);}
		if(order==0){order=a.model.name.compareTo(b.model.name);}
		return order==0 ? Float.compare(b.identity, a.identity) : order;
	};
	private String capture, captureManifest, population, assignments, prefix, parity, sampleManifest;
	private final ObjectMap<String, String> sampleFrames=new ObjectMap<String, String>(String.class, String.class);
	private boolean acceptanceOnly;
	private float replayIdentity=Float.NaN;
	private long savedWindowBases, noRosterOverlap, noSameStrandRosterOverlap, parityChecked;
	private final long[] cohortTotals=new long[6], cohortAccepted=new long[6];
	private final IntObjectMap<Gene> geneIds=new IntObjectMap<Gene>(16, Gene.class);
	private int expectedGenes, expectedFrames, expectedMembers, completedInputs;
	private long ambiguousEvents, ambiguousGenes, uniqueMembers, comparable, changed, missingGenes, noWindowGenes;
	private long geneWindowAssociations, trueWindowEvents;
	private Result matched;
	private int unresolvedReferences;
	private int[] finalCallByGene;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder(), event=new ByteBuilder();
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ObjectMap<String, Frame> frameNames=new ObjectMap<String, Frame>(String.class, Frame.class);
	private final ArrayList<Frame> frames=new ArrayList<Frame>();
	private ArrayList<Frame> slots;
	private final ArrayList<Capture> captures=new ArrayList<Capture>();
	private final ArrayList<Gene> genes=new ArrayList<Gene>();
	private final ArrayList<Call> calls=new ArrayList<Call>();
	private final ObjectMap<String, Model> models=new ObjectMap<String, Model>(String.class, Model.class);
	private final ArrayList<Model> modelOrder=new ArrayList<Model>();
	private final IntObjectMap<Model> construction=new IntObjectMap<Model>(16, Model.class);
	private static final Model UNMODELED=new Model("NA");
	private static final String FAMILY="RF00174";
	public static final String ASSOCIATION_HEADER="original_gene\tinput_id\tcompleted_event\tcaptured_event\twindow_overlaps\tcontains_gene\taccepted_span_qualifies\tevent_alignments\treference_resolved";
	public static final String CALLS_HEADER="call_id\tinput_id\tcontig\tstart1\tstop1\tstrand\tmodel\tidentity\tmatched_original_gene\texact_duplicate\tunmatched_status";
	public static final String GENES_HEADER="original_gene\tinput_id\tcontig\tstart1\tstop1\tstrand\tnative_truncation_bits\tconstruction_representative\tconstruction_model\tscheduled_overlaps\tcontained_windows\tproposal_state\tbest_unique_identity\tbest_ambiguous_identity\tselected_model\tcaptured_event\twindow_start0\twindow_stop0\tgc_bases\tacgt_bases\tconstruction_changed\tfinal_detected\tfinal_call\toverlapping_event_alignments\tselected_event_alignments\textraction_state";
	public static final String ACCEPTANCE_GENES_HEADER=GENES_HEADER.replace("\tfinal_detected\tfinal_call", "");
	private static final String[] COHORTS={"ALL_ORIGINALS", "EXACT_FULL_MODEL_SPAN", "INCOMPLETE_COVERAGE", "MISSING44LEFT20RIGHT", "NATIVE_TRUNCATED", "WITHIN_COVERAGE_TOLERANCE"};
	private static final String[] PROPOSAL_STATES={"UNRESOLVED_INPUT", "UNRESOLVED_EXTRACTION", "UNIQUE_BEST_ACCEPTANCE",
		"AMBIGUOUS_BEST_ACCEPTANCE", "NO_SCHEDULED_WINDOW", "PARTIAL_WINDOW_NO_ACCEPTANCE", "WINDOW_WITHOUT_ACCEPTED_MODEL"};
	private static final String FRAME_HEADER="input_id\tassembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\toriginal_genes\tinput_state";
	private static final String INPUT_HEADER="input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tscheduled_windows\talignments\taccepted_proposals\tfinal_target_call_rows";
	// Only these four extractor states have proved the source coordinates. A
	// construction exclusion is deliberately different from an unresolved source.
	private static final String[] EXTRACTION_STATES={"CONSTRUCTION_REPRESENTATIVE", "EXACT_SEQUENCE_DUPLICATE_FOR_CONSTRUCTION_ONLY",
		"EXCLUDED_NATIVE_PARTIAL_OR_TRUNCATED", "NO_ELIGIBLE_NATIVE_HIT_HAS_COMPONENT_BOUNDS", "SOURCE_UNAVAILABLE",
		"AMBIGUOUS_CONTIG_ID", "OUT_OF_RANGE_COORDINATES", "FRAME_CHANGED_DURING_READ", "FRAME_DIMENSION_MISMATCH", "CONTIG_NOT_FOUND"};
}
