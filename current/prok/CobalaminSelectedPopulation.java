package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntObjectMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.KillSwitch;
import shared.Timer;
import structures.ByteBuilder;
import template.ThreadWaiter;
import template.Accumulator;

/** Recounts three cobalamin populations on an explicit assembly manifest.
 * Uses template.A_SampleBasic2's parse/process/output lifecycle and native ncRNA
 * text owners. No DNA, alignment, caller or biological reclassification occurs.
 * Native GA admission comes from retained event provenance, not reference tiers.
 * Exact same-member coordinates collapse repeated events within one input only;
 * original event identities remain in the output. Reference masks and protected
 * priors are counted separately and never select construction membership.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSelectedPopulation implements Accumulator<CobalaminSelectedPopulation.ProcessThread>{

	public static void main(final String[] args) throws IOException{
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, CobalaminSelectedPopulation.class, false);
		try{
			final CobalaminSelectedPopulation tool=new CobalaminSelectedPopulation(pp.args, pp.outstream);
			tool.process();
			timer.stop();
			pp.outstream.println("COBALAMIN_SELECTED_POPULATION_COMPLETE time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private CobalaminSelectedPopulation(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("cohorts")){
				cohortPath=b;
			}else if(a.equals("provenance")){
				provenance=b;
			}else if(a.equals("expectedgenomes")){
				expectedGenomes=Parse.parseIntKMG(b);
			}else{
				require(parser.parse(arg, a, b), "Unknown population argument: "+arg);
			}
		}
		manifest=parser.in1;
		directory=parser.out1;
		require(manifest!=null && cohortPath!=null && provenance!=null && directory!=null && expectedGenomes>0
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1,
			"Supply in=selected_MANIFEST.tsv provenance=original_MANIFEST.tsv cohorts=bindings.tsv out=fresh_directory expectedgenomes=N without sampling");
		require(!new File(directory).exists(), "Recount cannot overwrite a previous result: "+directory);
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_PIGZ=false;
		ReadWrite.setZipThreads(1);
	}

	private void process() throws IOException{
		loadCohorts();
		loadManifest();
		loadProvenance();
		for(final Cohort cohort:cohorts){
			bindSources(cohort);
			loadEvents(cohort);
			loadWitnesses(cohort);
		}
		require(new File(directory).mkdir(), "Create one fresh recount directory: "+directory);
		final String linkPath=new File(directory, "raw_reference_links.tsv.gz").getPath();
		try(Output links=new Output(FileFormat.testOutput(linkPath, FileFormat.TEXT, null, false, false, false, true))){
			links.writer.add(new ByteBuilder().append(LINK_HEADER).nl(), 0);
			final AtomicInteger next=new AtomicInteger();
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>();
			for(int i=0; i<Math.min(Shared.threads(), assemblies.size()); i++){
				workers.add(new ProcessThread(this, next, links.writer));
			}
			final boolean complete=ThreadWaiter.startAndWait(workers, this);
			writeProducerStatus();
			require(complete, "All reference readers must finish before population publication");
			try(Output examples=output("mismatch_examples.tsv")){
				examples.writer.println(LINK_HEADER);
				for(final ProcessThread worker:workers){
					for(final byte[] example:worker.examples){examples.writer.print(example);}
				}
			}
		}
		writePopulation();
		writeSummary();
	}

	/** ThreadWaiter.waitAndAccumulate joins every worker before merging counters.
	 * Ordered ByteStreamWriter jobs own all data merging; no ledger is retained here. */
	@Override
	public void accumulate(final ProcessThread worker){
		assert(Thread.currentThread()==ownerThread) : "ThreadWaiter.waitAndAccumulate merges counters on its joining caller without worker contention";
		synchronized(worker){
			if(worker.failure!=null){
				if(firstFailure==null){firstFailure=worker.failure;}
				return;
			}
			mismatches+=worker.mismatched;
			correspondences+=worker.linked;
		}
	}
	@Override
	public ReadWriteLock rwlock(){return null;}
	@Override
	public boolean success(){return firstFailure==null;}

	/** Statuses describe correspondence production, never biological negatives.
	 * Failed/skipped jobs carry no measured byte count. This ledger survives a
	 * failed run, but no population table can be emitted unless every input passed. */
	private void writeProducerStatus(){
		assert(Thread.currentThread()==ownerThread) : "Joined worker states must be inspected before any population counts are emitted";
		final long[] counts=new long[PRODUCER_STATES.length];
		for(int index=0; index<assemblies.size(); index++){
			final Assembly assembly=assemblies.get(index);
			require(assembly.state>=0 && assembly.state<PRODUCER_STATES.length,
				"Unknown producer state for input ordinal "+index+": "+assembly.state);
			counts[assembly.state]++;
		}
		validateProducerCounts(counts[COMPLETE_EMPTY], counts[COMPLETE_NONEMPTY], counts[FAILED], counts[SKIPPED], expectedGenomes);
		final ByteBuilder failures=new ByteBuilder();
		try(Output status=output("producer_status.tsv"); Output summary=output("producer_status_summary.tsv")){
			status.writer.println("input_ordinal0\twriter_job_id\tassembly_id\tinput_id\tstate\toutput_bytes\treason");
			for(int index=0; index<assemblies.size(); index++){
				final Assembly assembly=assemblies.get(index);
				buffer.clear().append(index).tab().append(index+1).tab().append(assembly.id).tab().append(assembly.input)
					.tab().append(PRODUCER_STATES[assembly.state]).tab().append(assembly.outputBytes).tab();
				if(assembly.state==FAILED){
					buffer.append(assembly.failureReason);
					failures.append("; input ordinal0=").append(index).append(" assembly=").append(assembly.id)
						.append(" reason=").append(assembly.failureReason);
				}else if(assembly.state==SKIPPED){
					buffer.append("after_failed_input_ordinal0=").append(assembly.failedAfter);
				}else{buffer.append('.');}
				status.writer.print(buffer.nl());
			}
			summary.writer.println("state\tinputs\tunit");
			for(int state=1; state<counts.length; state++){
				summary.writer.print(buffer.clear().append(PRODUCER_STATES[state]).tab().append(counts[state])
					.append("\tCORRESPONDENCE_PRODUCER_STATUS_NOT_BIOLOGICAL_POPULATION\n"));
			}
		}
		if(counts[FAILED]>0 || counts[SKIPPED]>0 || counts[UNPROCESSED]>0 || firstFailure!=null){
			throw new IllegalStateException("Reference producers failed: failed="+counts[FAILED]+", skipped_after_failure="+counts[SKIPPED]
				+", unprocessed="+counts[UNPROCESSED]+failures.toString()+"; details="+new File(directory, "producer_status.tsv"), firstFailure);
		}
	}

	/** Four exclusive terminal states must conserve the exact selected manifest.
	 * This always-on guard prevents an input outside all states silently vanishing
	 * from the downstream population denominator, even when assertions are disabled. */
	public static void validateProducerCounts(final long empty, final long nonempty, final long failed, final long skipped, final long manifestTotal){
		require(empty>=0 && nonempty>=0 && failed>=0 && skipped>=0 && manifestTotal>0,
			"Producer-state counts must be nonnegative against a positive selected manifest total");
		final long accounted=Math.addExact(Math.addExact(empty, nonempty), Math.addExact(failed, skipped));
		if(accounted!=manifestTotal){
			throw new IllegalStateException("Producer state conservation mismatch: accounted="+accounted+", selected_manifest_total="+manifestTotal
				+"; completed_empty="+empty+", completed_nonempty="+nonempty+", failed="+failed+", skipped_after_failure="+skipped
				+". This guards against an input escaping every state and silently vanishing from the population denominator.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Metadata Input       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadCohorts(){
		try(Input in=new Input(cohortPath)){
			RrnaResourceIO.header(in.reader, "cohort\treference\tevents\twitnesses", cohortPath);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==4, "Each cohort must name its exact reference, event and witness files");
				final Cohort cohort=new Cohort(p.parseString(0), p.parseString(1), p.parseString(2), p.parseString(3));
				require(!byCohort.contains(cohort.name), "Cohort names must be unique: "+cohort.name);
				byCohort.put(cohort.name, cohort);
				cohorts.add(cohort);
			}
		}
		require(!cohorts.isEmpty(), "A selected population requires its raw-event cohort bindings");
	}

	private void loadManifest(){
		final ObjectMap<String, Assembly> ids=new ObjectMap<String, Assembly>(String.class, Assembly.class);
		try(Input in=new Input(manifest)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, manifest);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==12, "Use the selected twelve-column assembly/new-GFF manifest");
				final Assembly assembly=new Assembly(p.parseString(0), p.parseInt(1), p.parseString(2),
					p.parseString(3), p.parseString(4), p.parseString(5), p.parseString(7), p.parseString(8));
				require(byCohort.contains(assembly.selectionCohort) || assembly.selectionCohort.equals("older_ANI95_addition"),
					"Unrecognized selection cohort '"+assembly.selectionCohort+"' for "+assembly.id);
				require(assembly.input>=0 && !ids.contains(assembly.id) && !byInput.contains(assembly.input)
					&& new File(assembly.fasta).isAbsolute() && new File(assembly.gff).isAbsolute(),
					"Manifest assembly IDs and cohort/input identities must be unique, with absolute paired paths");
				ids.put(assembly.id, assembly);
				byInput.put(assembly.input, assembly);
				assemblies.add(assembly);
			}
		}
		require(assemblies.size()==expectedGenomes, "Recount must include every selected genome, including zeros: "+assemblies.size());
	}

	/** The broader directory is an identity map only; it cannot add a selected input.
	 * ANI addition is a selection role, not the cohort that produced native calls. */
	private void loadProvenance() throws IOException{
		assert(assemblies.size()==expectedGenomes) : "Only the selected manifest creates population members";
		try(Input in=new Input(provenance)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER.substring(0, MANIFEST_HEADER.lastIndexOf('\t')), provenance);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==11, "Original corpus provenance retains eleven columns");
				final Assembly assembly=byInput.get(p.parseInt(1));
				if(assembly==null){continue;}
				final String origin=p.parseString(2);
				final Cohort cohort=byCohort.get(origin);
				require(cohort!=null, "Unrecognized provenance cohort '"+origin+"' for "+assembly.id);
				require(assembly.cohort==null && cohort!=null
					&& (assembly.selectionCohort.equals(cohort.name) || assembly.selectionCohort.equals("older_ANI95_addition")),
					"Each selected input must map once to its original producing cohort: "+assembly.id);
				require(samePath(assembly.source, p.parseString(7)) && samePath(assembly.sourceGff, p.parseString(8)),
					"Provenance binding requires exact physical source FASTA and new-GFF identities: "+assembly.id);
				assembly.cohort=cohort;
				cohort.inputs.put(assembly.input, assembly);
			}
		}
		for(final Assembly assembly:assemblies){require(assembly.cohort!=null, "Selected input lacks original provenance: "+assembly.id);}
	}

	/** Authenticate each selected input against its saved scan FASTA, not a taxID or basename. */
	private void bindSources(final Cohort cohort) throws IOException{
		assert(!assemblies.isEmpty()) : "Source binding follows complete manifest loading";
		try(Input in=new Input(cohort.reference)){
			RrnaResourceIO.header(in.reader, "input_id\tfasta\tgff", cohort.reference);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==3, "Scan reference manifest retains input/FASTA/GFF columns");
				final Assembly assembly=cohort.inputs.get(p.parseInt(0));
				if(assembly==null){continue;}
				require(!assembly.bound, "Repeated selected input in scan reference: "+assembly.id);
				final String source=p.parseString(1);
				require(source.equals(assembly.source) || new File(source).getCanonicalPath().equals(new File(assembly.source).getCanonicalPath()),
					"Selected assembly must bind its exact scanned physical FASTA: "+assembly.id);
				assembly.bound=true;
			}
		}
		for(final Assembly assembly:assemblies){
			if(assembly.cohort!=cohort){continue;}
			require(assembly.bound && new File(assembly.fasta).isFile()
				&& new File(assembly.fasta).getCanonicalPath().equals(new File(assembly.source).getCanonicalPath()),
				"Selected linked assembly must still resolve to its authenticated scanned source: "+assembly.id);
		}
	}

	private void loadEvents(final Cohort cohort){
		try(Input in=new Input(cohort.events)){
			RrnaResourceIO.header(in.reader, NcrnaSubtypeEvidence.EVENT_HEADER, cohort.events);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==30, "Native provenance requires all thirty fields");
				final int member=member(6);
				if(member<0 || !p.termEquals("true", 17)){continue;}
				final Assembly assembly=cohort.inputs.get(p.parseInt(1));
				if(assembly==null){continue;}
				final int eventId=p.parseInt(0), start=p.parseInt(3), stop=p.parseInt(4);
				final char strand=p.termEquals("+", 5) ? '+' : p.termEquals("-", 5) ? '-' : '.';
				require(eventId>0 && !cohort.eventsById.contains(eventId) && start>0 && stop>=start && strand!='.',
					"Native-GA events must retain unique IDs and valid source-bound inclusive coordinates");
				final String contig=decoded(2);
				final String key=member+"\t"+contig+"\t"+start+"\t"+stop+"\t"+strand;
				Gene gene=assembly.byLocus.get(key);
				if(gene==null){
					gene=new Gene(assembly, member, contig, start, stop, strand);
					assembly.byLocus.put(key, gene);
					assembly.genes.add(gene);
				}
				final boolean truncated=p.parseInt(11)!=0;
				final Event event=new Event(eventId, gene, truncated, p.parseString(26), p.parseLong(27));
				require(event.nativeLine>0, "Every admitted event retains its producing native table line");
				cohort.eventsById.put(eventId, event);
				assembly.events.add(event);
				gene.events++;
				if(!truncated){gene.untruncated++;}
			}
		}
	}

	/** Witness links establish correspondence; overlap alone never chooses a parent. */
	private void loadWitnesses(final Cohort cohort){
		try(Input in=new Input(cohort.witnesses)){
			RrnaResourceIO.header(in.reader, NcrnaSubtypeEvidence.WITNESS_HEADER, cohort.witnesses);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==11, "Annotation provenance retains all eleven witness columns");
				if(member(6)<0){continue;}
				final Assembly assembly=cohort.inputs.get(p.parseInt(1));
				if(assembly==null){continue;}
				final String id=decoded(0);
				Parent parent=assembly.parents.get(id);
				if(parent==null){
					parent=new Parent();
					assembly.parents.put(id, parent);
				}
				if(p.termEquals("PROTECTED_PRIOR", 2)){
					parent.prior=true;
				}else{
					require(p.termEquals("CM_HIT", 2), "Unknown annotation evidence kind");
					final Event event=cohort.eventsById.get(p.parseInt(3));
					if(event==null){continue;}// Non-native-GA observations never enter construction.
					require(event.gene.assembly==assembly && event.gene.member==member(6)
						&& event.gene.start==p.parseInt(7) && event.gene.stop==p.parseInt(8)
						&& event.nativeLine==p.parseLong(5) && p.termEquals(event.nativeTable, 4),
						"Reference witness must bind the exact native event, model, source and table line");
					parent.links.add(new Link(event, p.parseString(10)));
				}
			}
		}
	}

	private void readReference(final Assembly assembly, final ProcessThread worker, final ByteBuilder job){
		assert(assembly.bound) : "Reference comparison requires an authenticated source FASTA binding";
		final ByteBuilder rowBuffer=worker.rowBuffer;
		for(final Row row:worker.codec.read(assembly.gff, Integer.toString(assembly.input))){
			if(!row.locus.family.equals("cobalamin")){continue;}
			final Parent parent=assembly.parents.get(row.id);
			require(parent!=null, "Every cobalamin reference parent requires its saved witness provenance: "+row.id);
			assembly.categories[row.mask ? 2 : parent.prior ? 1 : 0]++;
			if(row.mask){continue;}
			for(final Link link:parent.links){
				final Event event=link.event;
				final Gene gene=event.gene;
				final boolean exact=gene.contig.equals(row.locus.seqid) && gene.start==row.locus.start
					&& gene.stop==row.locus.stop && gene.strand==row.locus.strand;
				if(exact){gene.exactScored++;}
				else{worker.mismatched++;}
				worker.linked++;
				rowBuffer.clear().append(assembly.id).tab().append(assembly.input).tab().append(MEMBERS[gene.member])
					.tab().append(event.id).tab().append(row.id).tab().append(link.association).tab().append(parent.prior)
					.tab().append(gene.contig).tab().append(gene.start).tab().append(gene.stop).tab().append(gene.strand)
					.tab().append(row.locus.seqid).tab().append(row.locus.start).tab().append(row.locus.stop).tab()
					.append(row.locus.strand).tab().append(exact).nl();
				job.append(rowBuffer);
				if(!exact && worker.examples.size()<3){worker.examples.add(rowBuffer.toBytes());}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writePopulation(){
		try(Output loci=output("native_loci.tsv.gz"); Output events=output("native_events.tsv.gz"); Output counts=output("inputs.tsv")){
			loci.writer.println("assembly_id\tinput_id\tcohort\tsource_fasta\tmember\tcontig\tstart\tstop\tstrand\tnative_events\tuntruncated_events\tconstruction_eligible\texact_scored_links");
			events.writer.println("assembly_id\tinput_id\tmember\tevent_id\tcontig\tstart\tstop\tstrand\tnative_truncated\tnative_table\tnative_line");
			counts.writer.println("assembly_id\tinput_id\tcohort\tmember\tnative_events\tnative_loci\tconstruction_loci\ttruncation_only_loci");
			for(final Assembly assembly:assemblies){
				final long[][] values=new long[3][4];
				for(final Gene gene:assembly.genes){
					final long[] v=values[gene.member];
					v[0]+=gene.events;
					v[1]++;
					v[gene.untruncated>0 ? 2 : 3]++;
					loci.writer.print(buffer.clear().append(assembly.id).tab().append(assembly.input).tab().append(assembly.cohort.name)
						.tab().append(assembly.source).tab().append(MEMBERS[gene.member]).tab().append(gene.contig).tab().append(gene.start)
						.tab().append(gene.stop).tab().append(gene.strand).tab().append(gene.events).tab().append(gene.untruncated)
						.tab().append(gene.untruncated>0).tab().append(gene.exactScored).nl());
				}
				for(final Event event:assembly.events){
					final Gene gene=event.gene;
					events.writer.print(buffer.clear().append(assembly.id).tab().append(assembly.input).tab().append(MEMBERS[gene.member])
						.tab().append(event.id).tab().append(gene.contig).tab().append(gene.start).tab().append(gene.stop).tab().append(gene.strand)
						.tab().append(event.truncated).tab().append(event.nativeTable).tab().append(event.nativeLine).nl());
				}
				for(int member=0; member<3; member++){
					final long[] v=values[member];
					require(v[1]==v[2]+v[3], "Truncation-only exclusions partition native coordinate loci");
					buffer.clear().append(assembly.id).tab().append(assembly.input).tab().append(assembly.cohort.name).tab().append(MEMBERS[member]);
					for(int column=0; column<4; column++){
						buffer.tab().append(v[column]);
						totals[member][column]+=v[column];
					}
					if(v[1]>0){totals[member][4]++;}
					counts.writer.print(buffer.nl());
				}
			}
		}
	}

	private void writeSummary(){
		try(Output summary=output("summary.tsv"); Output reference=output("reference_categories.tsv"); Output result=output("result.tsv")){
			summary.writer.println("member\tgenomes\tnative_events\tnative_loci\tconstruction_loci\ttruncation_only_loci\tpositive_genomes");
			for(int member=0; member<3; member++){
				buffer.clear().append(MEMBERS[member]).tab().append(assemblies.size());
				for(final long value:totals[member]){buffer.tab().append(value);}
				summary.writer.print(buffer.nl());
			}
			final long[] categories=new long[3];
			for(final Assembly assembly:assemblies){
				for(int i=0; i<3; i++){categories[i]+=assembly.categories[i];}
			}
			reference.writer.println("scope\tfresh_scored_CM_parents\tcarried_prior_scored_parents\tdont_care_parents\tunit");
			reference.writer.print(buffer.clear().append("cobalamin_union"));
			for(final long value:categories){reference.writer.print(buffer.clear().tab().append(value));}
			reference.writer.println("\tREFERENCE_PARENTS_NOT_RAW_NATIVE_CALLS");
			result.writer.println("metric\tvalue");
			result.writer.print(buffer.clear().append("raw_scored_witness_links\t").append(correspondences).nl()
				.append("nonidentical_raw_scored_links\t").append(mismatches).nl()
				.append("construction_status\t").append(mismatches==0 ? "AWAITING_SOURCE_REVIEW" : "BLOCKED_RAW_REFERENCE_COORDINATE_DIFFERENCES").nl());
		}
		outstream.println("Selected genomes="+assemblies.size()+"; nonidentical raw/scored witness links="+mismatches);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private int member(final int column){
		for(int i=0; i<MEMBERS.length; i++){if(p.termEquals(MEMBERS[i], column)){return i;}}
		return -1;
	}
	private String decoded(final int column){
		p.setBounds(column);
		return decoder.decode(p.line(), p.a(), p.b());
	}
	private Output output(final String name){return new Output(new File(directory, name).getPath());}
	private static boolean samePath(final String a, final String b) throws IOException{
		return a.equals(b) || new File(a).getCanonicalPath().equals(new File(b).getCanonicalPath());
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Retained Records      ----------------*/
	/*--------------------------------------------------------------*/

	static final class ProcessThread extends Thread{
		ProcessThread(final CobalaminSelectedPopulation t, final AtomicInteger n, final ByteStreamWriter w){
			tool=t;
			next=n;
			writer=w;
		}
		@Override
		public void run(){
			try{
				for(int index=next.getAndIncrement(); index<tool.assemblies.size(); index=next.getAndIncrement()){
					final Assembly assembly=tool.assemblies.get(index);
					final ByteBuilder job=new ByteBuilder(4096);
					try{
						if(failure==null){
							tool.readReference(assembly, this, job);
							assembly.state=job.length()==0 ? COMPLETE_EMPTY : COMPLETE_NONEMPTY;
							assembly.outputBytes=job.length();
						}else{
							assembly.state=SKIPPED;
							assembly.failedAfter=failedInput;
						}
					}catch(Throwable problem){
						failure=problem;
						failedInput=index;
						assembly.state=FAILED;
						assembly.failureReason=(problem.getClass().getSimpleName()+": "+problem.getMessage())
							.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
						job.clear();
					}finally{
						// A failed input still fills its ordered ID, so other producers
						// cannot block behind it. After failure, drain further claimed
						// ordinals with empty jobs; the accumulator rejects the run.
						writer.add(job, index+1L);
						// Ownership transferred: never clear or reuse this builder/array.
					}
				}
			}catch(Throwable problem){
				// Submission itself failed: a missing ordinal could block all other
				// producers before ThreadWaiter can join them. Terminate loudly.
				problem.printStackTrace(tool.outstream);
				KillSwitch.kill("Ordered correspondence submission failed before every input ordinal was supplied: "+problem);
			}
		}
		final AtomicInteger next;
		final CobalaminSelectedPopulation tool;
		final ByteStreamWriter writer;
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ByteBuilder rowBuffer=new ByteBuilder();
		final ArrayList<byte[]> examples=new ArrayList<byte[]>();
		long linked, mismatched;
		int failedInput=-1;
		Throwable failure;
	}

	private static final class Cohort{
		Cohort(final String n, final String r, final String e, final String w){
			name=n;
			reference=r;
			events=e;
			witnesses=w;
		}
		final String name, reference, events, witnesses;
		final IntObjectMap<Assembly> inputs=new IntObjectMap<Assembly>(Assembly.class);
		final IntObjectMap<Event> eventsById=new IntObjectMap<Event>(Event.class);
	}
	private static final class Assembly{
		Assembly(final String i, final int n, final String c, final String f, final String g, final String t, final String s, final String sg){
			id=i;
			input=n;
			selectionCohort=c;
			fasta=f;
			gff=g;
			taxids=t;
			source=s;
			sourceGff=sg;
		}
		final String id, fasta, gff, taxids, source, sourceGff, selectionCohort;
		final int input;
		Cohort cohort;
		boolean bound;
		int state=UNPROCESSED, failedAfter=-1, outputBytes=-1;
		String failureReason;
		final ObjectMap<String, Gene> byLocus=new ObjectMap<String, Gene>(String.class, Gene.class);
		final ObjectMap<String, Parent> parents=new ObjectMap<String, Parent>(String.class, Parent.class);
		final ArrayList<Gene> genes=new ArrayList<Gene>();
		final ArrayList<Event> events=new ArrayList<Event>();
		final long[] categories=new long[3];
	}
	private static final class Gene{
		Gene(final Assembly a, final int m, final String c, final int s, final int z, final char d){
			assembly=a;
			member=m;
			contig=c;
			start=s;
			stop=z;
			strand=d;
		}
		final Assembly assembly;
		final int member, start, stop;
		final String contig;
		final char strand;
		long events, untruncated, exactScored;
	}
	private static final class Event{
		Event(final int i, final Gene g, final boolean t, final String n, final long l){
			id=i;
			gene=g;
			truncated=t;
			nativeTable=n;
			nativeLine=l;
		}
		final int id;
		final Gene gene;
		final boolean truncated;
		final String nativeTable;
		final long nativeLine;
	}
	private static final class Link{
		Link(final Event e, final String a){
			event=e;
			association=a;
		}
		final Event event;
		final String association;
	}
	private static final class Parent{
		boolean prior;
		final ArrayList<Link> links=new ArrayList<Link>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String manifest, cohortPath, provenance, directory;
	private int expectedGenomes=-1;
	private long correspondences, mismatches;
	private Throwable firstFailure;
	private final Thread ownerThread=Thread.currentThread();
	private final PrintStream outstream;
	private final LineParser1 p=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectMap<String, Cohort> byCohort=new ObjectMap<String, Cohort>(String.class, Cohort.class);
	private final ArrayList<Cohort> cohorts=new ArrayList<Cohort>();
	private final ArrayList<Assembly> assemblies=new ArrayList<Assembly>();
	private final IntObjectMap<Assembly> byInput=new IntObjectMap<Assembly>(Assembly.class);
	private final long[][] totals=new long[3][5];
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final int UNPROCESSED=0, COMPLETE_EMPTY=1, COMPLETE_NONEMPTY=2, FAILED=3, SKIPPED=4;
	private static final String[] PRODUCER_STATES={"UNPROCESSED", "COMPLETE_EMPTY_OUTPUT", "COMPLETE_NONEMPTY_OUTPUT", "FAILED", "SKIPPED_AFTER_FAILURE"};
	private static final String LINK_HEADER="assembly_id\tinput_id\tmember\tevent_id\tparent_id\tassociation\tprior\traw_contig\traw_start\traw_stop\traw_strand\treference_contig\treference_start\treference_stop\treference_strand\tcoordinate_identity";
	private static final String MANIFEST_HEADER="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
}
