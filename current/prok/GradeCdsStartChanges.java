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

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashMap;
import map.IntHashSet;
import map.IntObjectMap;
import map.LongIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import prok.NcrnaObservationData.Symbols;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PanelAnnotationGrade.Feature;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Paired, source-relative CDS grading from saved per-input GFFs, without DNA.
 * Adapted from template.A_SampleByteFileMT: native byte I/O, coarse workers,
 * checked writer closure and one final counter merge. Retains one genome per
 * worker and delegates overlap/endpoint/frame metrics to PanelAnnotationGrade.
 * Same-stop start attribution is separate and excludes uncertain annotations.
 * @author Brian, Raiden
 */
public final class GradeCdsStartChanges implements Accumulator<GradeCdsStartChanges.ProcessThread>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final GradeCdsStartChanges tool=new GradeCdsStartChanges(args);
		try{tool.process();}
		finally{Shared.closeStream(tool.outstream);}
	}

	private GradeCdsStartChanges(final String[] args){
		final PreParser pp=new PreParser(args, GradeCdsStartChanges.class, false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("corpus")){corpus=b;}
			else if(a.equals("expectedinputs")){expectedInputs=Parse.parseIntKMG(b);}
			else if(a.equals("rnacallchanges")){rnaCallChanges=Parse.parseBoolean(b);}
			else if(a.equals("familymap")){familymap=b;}
			else{require(parser.parse(arg, a, b), "Unknown CDS grading option: "+arg);}
		}
		before=parser.in1;
		after=parser.in2;
		out=parser.out1;
		require(before!=null && after!=null && corpus!=null && out!=null && expectedInputs>0
			&& !parser.overwrite && !parser.append && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply complete in=before_manifest in2=after_manifest corpus=formal_manifest expectedinputs=N out=fresh_prefix");
		require(Tools.testInputFiles(false, true, before, after, corpus), "Saved-call manifests and formal assembly provenance must be readable");
		ReadWrite.setZipThreads(1);
		if(rnaCallChanges){
			require(familymap!=null && Tools.testInputFiles(false, true, familymap),
				"rnacallchanges=t requires the same explicit familymap= projection as RNA grading");
			rnaRules=GradeNcrnaReference.loadRules(familymap);
			suffixes=Arrays.copyOf(SUFFIXES, SUFFIXES.length+1);
			suffixes[SUFFIXES.length]=".rna_call_changes.tsv";
		}else{require(familymap==null, "familymap= belongs to explicit rnacallchanges=t, never a silent no-op");}
		RrnaResourceIO.fresh(out, suffixes);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Processing         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final Timer timer=new Timer();
		loadJobs();
		try(Output inputs=ordered(".inputs.tsv", INPUT_HEADER); Output changes=ordered(".changes.tsv", CHANGE_HEADER);
			Output others=ordered(".rna_counts.tsv", RNA_HEADER); Output sources=ordered(".sources.tsv", SOURCE_HEADER);
			Output rna=rnaCallChanges ? ordered(".rna_call_changes.tsv", RNA_CHANGE_HEADER) : null){
			final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>();
			for(int i=0, n=Tools.min(Shared.threads(), jobs.size()); i<n; i++){
				threads.add(new ProcessThread(this, inputs.writer, changes.writer, others.writer, sources.writer,
					rna==null ? null : rna.writer));
			}
			errorState|=!ThreadWaiter.startAndWait(threads, this);
		}
		require(!errorState && completed==jobs.size(), "Every declared input must finish before a CDS summary is valid");
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("arm\tmetric\tcount\tqualification");
			for(int a=0; a<2; a++){
				for(int i=0; i<STAT_NAMES.length; i++){
					summary.writer.print(new ByteBuilder(ARM[a]).tab().append(STAT_NAMES[i]).tab().append(totals[a][i]).tab().append(QUALIFICATION).nl());
				}
			}
			for(int i=0; i<CHANGE_NAMES.length; i++){
				summary.writer.print(new ByteBuilder("paired").tab().append(CHANGE_NAMES[i]).tab().append(changeTotals[i]).tab().append(QUALIFICATION).nl());
			}
		}
		timer.stop();
		outstream.println("CDS_START_CHANGES_GRADED inputs="+completed+" time="+timer+"; "+QUALIFICATION);
	}

	private void loadJobs(){
		final ArrayList<Binding> old=GradeNcrnaReference.loadManifest(before), next=GradeNcrnaReference.loadManifest(after);
		require(old.size()==expectedInputs && next.size()==old.size(), "Both call manifests must contain the full declared panel");
		final Symbols ids=new Symbols();
		final IntObjectMap<Job> byId=new IntObjectMap<Job>(16, Job.class);
		for(int i=0; i<old.size(); i++){
			final Binding b=old.get(i), a=next.get(i);
			require(b.id.equals(a.id) && b.fasta.equals(a.fasta), "Paired saved calls must bind the same exact input and FASTA in the same order");
			final Job job=new Job(b, a);
			final byte[] name=b.id.getBytes(StandardCharsets.UTF_8);
			byId.put(ids.intern(name, 0, name.length), job);
			jobs.add(job);
		}
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(corpus)){
			RrnaResourceIO.header(input.reader, CORPUS_HEADER, corpus);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==12, "Formal assembly provenance must retain its twelve-field schema");
				p.setBounds(1);
				final int id=ids.find(line, p.a(), p.b());
				final Job job=id<0 ? null : byId.get(id);
				if(job==null){continue;}
				require(job.original==null && p.termEquals(job.old.fasta, 3), "Each panel input needs exactly one matching formal assembly provenance row");
				final String fasta=p.parseString(7);
				require(fasta.endsWith(".fna.gz"), "Adjacent original GFF derivation requires the declared .fna.gz source convention: "+fasta);
				job.original=fasta.substring(0, fasta.length()-7)+".gff.gz";
			}
		}
		for(final Job job:jobs){
			require(job.original!=null && new File(job.original).isFile(), "Missing original CDS-containing annotation for "+job.old.id+": "+job.original);
			require(Tools.testInputFiles(false, true, job.old.gff, job.next.gff), "Both saved per-input GFFs must exist for "+job.old.id);
			for(final String suffix:suffixes){
				final String output=new File(out+suffix).getAbsolutePath();
				require(!output.equals(job.old.gff) && !output.equals(job.next.gff) && !output.equals(job.original)
					&& !output.equals(job.old.fasta), "An output cannot overwrite a bound genomic input");
			}
		}
	}

	private Output ordered(final String suffix, final String header){
		final Output output=new Output(new ByteStreamWriter(FileFormat.testOutput(out+suffix, FileFormat.TXT, null, false, false, false, true)));
		output.writer.add(new ByteBuilder(header).nl(), 0);
		return output;
	}

	@Override
	public void accumulate(final ProcessThread thread){
		errorState|=!thread.success;
		completed+=thread.completed;
		for(int a=0; a<2; a++){PanelAnnotationGrade.add(totals[a], thread.totals[a]);}
		PanelAnnotationGrade.add(changeTotals, thread.changes);
	}
	@Override
	public boolean success(){return !errorState;}
	@Override
	public ReadWriteLock rwlock(){return lock;}

	/*--------------------------------------------------------------*/
	/*----------------         Input Worker         ----------------*/
	/*--------------------------------------------------------------*/

	static final class ProcessThread extends Thread{
		ProcessThread(final GradeCdsStartChanges owner_, final ByteStreamWriter inputs_, final ByteStreamWriter changes_,
				final ByteStreamWriter others_, final ByteStreamWriter sources_, final ByteStreamWriter rna_){
			owner=owner_; inputs=inputs_; changeWriter=changes_; others=others_; sources=sources_;
			rna=rna_;
			projection=owner.rnaCallChanges ? new NcrnaReferenceCalls(owner.rnaRules) : null;
		}
		@Override
		public void run(){
			try{
				for(int index; (index=owner.nextInput.getAndIncrement())<owner.jobs.size();){
					final Job job=owner.jobs.get(index);
					final Genome genome=new Genome(job, projection);
					genome.read(job.original, -1);
					genome.validateReference();
					genome.read(job.old.gff, 0);
					genome.read(job.next.gff, 1);
					final ByteBuilder inputRows=new ByteBuilder(), changeRows=new ByteBuilder(), otherRows=new ByteBuilder();
					final long[][] inputStats=new long[2][PanelAnnotationGrade.NSTAT];
					for(final Contig contig:genome.order){grade(job.old.id, genome, contig, inputStats, changeRows);}
					for(int a=0; a<2; a++){
						require(inputStats[a][0]+inputStats[a][2]==genome.referenceRows, "Every source CDS row must be matched or missed within its own input");
						PanelAnnotationGrade.add(totals[a], inputStats[a]);
						inputRows.append(job.old.id).tab().append(ARM[a]);
						for(final long value:inputStats[a]){inputRows.tab().append(value);}
						inputRows.tab().append(QUALIFICATION).nl();
					}
					genome.rnaCounts(otherRows, changes);
					if(rna!=null){rna.add(genome.rnaChanges(new ByteBuilder()), index+1);}
					final ByteBuilder sourceRow=new ByteBuilder(job.old.id).tab().append(job.old.fasta).tab().append(job.original)
						.tab().append(genome.referenceRows).tab().append(genome.order.size()).append("\tDECLARED_GFF_TOPOLOGY_NO_DNA_REREAD\n");
					// Coarse ordered transfers, including EMPTY jobs, keep memory bounded and empty inputs reachable.
					inputs.add(inputRows, index+1);
					changeWriter.add(changeRows, index+1);
					others.add(otherRows, index+1);
					sources.add(sourceRow, index+1);
					completed++;
				}
				success=true;
			}catch(OutOfMemoryError failure){
				KillSwitch.memKill(failure);
			}catch(Throwable failure){
				failure.printStackTrace(owner.outstream);
				KillSwitch.kill("CDS grading worker failed; ordered output cannot remain valid: "+failure);
			}
		}

		private void grade(final String input, final Genome genome, final Contig c, final long[][] inputStats, final ByteBuilder changeRows){
			final int circle=c.circular ? c.length : 0;
			final MatchResult[] matched=new MatchResult[2];
			for(int a=0; a<2; a++){
				matched[a]=PanelAnnotationGrade.match(c.calls[a], c.truth, circle);
				final long[] stats=PanelAnnotationGrade.statistics(c.calls[a], c.truth, matched[a], true, circle);
				if(a==1){PanelAnnotationGrade.changes(c.calls[1], c.calls[0], stats);}
				PanelAnnotationGrade.add(inputStats[a], stats);
			}
			for(int i=0; i<c.truth.size(); i++){
				final boolean old=matched[0].truthMatchedCallIdx[i]>=0, next=matched[1].truthMatchedCallIdx[i]>=0;
				if(old && !next){changes[REFERENCE_LOST]++;}
				if(!old && next){changes[REFERENCE_GAINED]++;}
			}
			final LongIntMap old=index(c.calls[0], circle), next=index(c.calls[1], circle);
			for(final Feature b:c.calls[0]){
				final long key=stopKey(b, circle);
				final int n=next.get(key), previous=old.get(key);
				if(previous==-2 || n==-2){
					changes[AMBIGUOUS_CALL_STOP]++;
					change(input, genome, c, b, null, "AMBIGUOUS_BEFORE_STOP_MEMBER", changeRows);
					continue;
				}
				final Feature a=n<0 ? null : c.calls[1].get(n);
				if(previous>=0 && n>=0 && b.locus.start==a.locus.start && b.locus.stop==a.locus.stop){continue;}
				if(a==null){changes[STOP_LOST]++;}
				else{changes[START_CHANGED]++;}
				change(input, genome, c, b, a, a==null ? "STOP_NOT_CALLED" : "SAME_STOP_START_CHANGED", changeRows);
			}
			for(final Feature a:c.calls[1]){
				final long key=stopKey(a, circle);
				if(old.get(key)==-2 || next.get(key)==-2){
					changes[AMBIGUOUS_CALL_STOP]++;
					change(input, genome, c, null, a, "AMBIGUOUS_AFTER_STOP_MEMBER", changeRows);
					continue;
				}
				if(old.get(key)!=-1){continue;}
				changes[STOP_ADDED]++;
				change(input, genome, c, null, a, "ADDED_STOP", changeRows);
			}
		}

		private void change(final String input, final Genome genome, final Contig c, final Feature before,
				final Feature after, final String event, final ByteBuilder out){
			final Feature call=before==null ? after : before;
			final int circle=c.circular ? c.length : 0;
			Feature reference=null;
			int rows=0, firstStart=-1;
			boolean ambiguous=false, uncertain=false;
			for(int i=0; i<c.truth.size(); i++){
				final Feature r=c.truth.get(i);
				if(stopKey(r, circle)!=stopKey(call, circle)){continue;}
				rows++;
				final int start=physical(GradeCdsStartChanges.start(r), circle);
				if(firstStart>=0 && firstStart!=start){ambiguous=true;}
				firstStart=start;
				final int id=c.referenceIds.get(i);
				uncertain|=c.referenceFlags.get(i)!=0 || id<0 || genome.idCounts.get(id)!=1;
				reference=r;
			}
			final String status=rows==0 ? "NO_SAME_STOP_ANNOTATION" : ambiguous ? "AMBIGUOUS_ANNOTATED_START"
				: uncertain ? "UNCERTAIN_PARTIAL_OR_REPEATED_ID" : "UNIQUE_COMPLETE_ANNOTATED_START";
			final boolean reliable=rows>0 && !ambiguous && !uncertain && !event.startsWith("AMBIGUOUS_");
			final int oldDelta=!reliable || before==null ? 0 : startDelta(before, reference, circle);
			final int newDelta=!reliable || after==null ? 0 : startDelta(after, reference, circle);
			if(reliable && before!=null && after!=null){
				changes[RELIABLE_START_PAIR]++;
				if(oldDelta==0){changes[START_EXACT_BEFORE]++;}
				if(newDelta==0){changes[START_EXACT_AFTER]++;}
				if(Math.abs((long)newDelta)>Math.abs((long)oldDelta)){changes[START_ERROR_WORSE]++;}
				if(after.locus.stop-after.locus.start<before.locus.stop-before.locus.start && newDelta>0){changes[SHORTENED_PAST_START]++;}
			}
			out.append(input).tab().append(c.name).tab().append(call.locus.strand).tab().append(event);
			coordinates(out, before);
			coordinates(out, after);
			out.tab().append(rows).tab().append(status).tab();
			if(reliable && before!=null){out.append(oldDelta);}else{out.append("NA");}
			out.tab();
			if(reliable && after!=null){out.append(newDelta);}else{out.append("NA");}
			out.tab().append(circle).append("\tSTOP_NOT_CALLED_IS_NOT_PROVEN_GENE_LOSS;NO_STOP_MATCH_NOT_PROVEN_SPURIOUS\n");
		}

		private final GradeCdsStartChanges owner;
		private final ByteStreamWriter inputs, changeWriter, others, sources, rna;
		private final NcrnaReferenceCalls projection;
		private final long[][] totals=new long[2][PanelAnnotationGrade.NSTAT];
		private final long[] changes=new long[CHANGE_NAMES.length];
		private int completed;
		private boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Native GFF Reader     ----------------*/
	/*--------------------------------------------------------------*/

	/** One original CDS-containing annotation, parsed once without DNA or queries.
	 * Uses the same reader, physical/circular validation and per-input taxID aliases
	 * as the saved-call tool. The driver binds input/path to formal provenance;
	 * passing an ncRNA-only GFF here cannot establish complete CDS truth.
	 * No file remains open. One instance belongs to one worker; returned views
	 * expose immutable features and row metadata, never mutable genome collections.
	 */
	public static final class Reference{
		public Reference(final String input_, final String originalGff_){
			require(input_!=null && !input_.isEmpty() && originalGff_!=null && !originalGff_.isEmpty(),
				"An original CDS reference requires explicit input identity and annotation path");
			input=input_;
			originalGff=originalGff_;
			genome=new Genome(null);// Reference-only use never invokes saved-call output.
			genome.read(originalGff, -1);
			genome.validateReference();
		}
		/** Resolves only the original reader's literal or explicitly declared alias.
		 * Native headers use their emitted first token. Captured physical length must
		 * equal the original annotation's declaration; no clipping or inferred length.
		 */
		public ReferenceContig contig(final String nativeName, final int physicalLength){
			final ReferenceContig c=contig(nativeName);
			require(physicalLength>0 && physicalLength==c.length(),
				"Captured CDS contig length must equal its original annotation declaration");
			return c;
		}
		/** Resolves a saved GFF name through the original declared alias rules.
		 * The annotation already owns its physical length; no DNA is required.
		 */
		public ReferenceContig contig(final String nativeName){
			final byte[] name=NcrnaReferenceCalls.reportContig(nativeName).getBytes(StandardCharsets.UTF_8);
			final Contig c=genome.declaredQueryContig(name, 0, name.length);
			return new ReferenceContig(this, c);
		}
		/** Original contig order, including contigs with no CDS or selected query.
		 * Consumers must conserve reference-only contigs instead of dropping their FN.
		 */
		public ReferenceContig contig(final int index){return new ReferenceContig(this, genome.order.get(index));}
		public int contigs(){return genome.order.size();}
		public long cdsRows(){return genome.referenceRows;}
		/** Reads one saved pair through the grader's original contig, topology and
		 * GFF validation. A failed partial load cannot be retried on mutated state.
		 * No sequence input or annotation rematching is performed here.
		 */
		public void loadCalls(final String before, final String after){
			require(callState==0 && before!=null && after!=null,
				"A reference owns one complete saved-call pair; construct a fresh reference after a failed load");
			callState=1;
			genome.read(before, 0);
			genome.read(after, 1);
			callState=2;
		}
		public final String input, originalGff;
		private final Genome genome;
		private int callState;
	}

	/** Read-only original-order CDS rows for one declared physical contig. */
	public static final class ReferenceContig{
		private ReferenceContig(final Reference owner_, final Contig contig_){
			owner=owner_;
			contig=contig_;
			assert(contig.truth.size()==contig.referenceIds.size && contig.truth.size()==contig.referenceFlags.size)
				: "Original CDS row eligibility must remain aligned with the exact feature order retained by Genome.read";
		}
		public String input(){return owner.input;}
		public String name(){return contig.name;}
		public int length(){return contig.length;}
		public int circularLength(){return contig.circular ? contig.length : 0;}
		public int rows(){return contig.truth.size();}
		public Feature feature(final int index){return contig.truth.get(index);}
		/** Original selected-call order, as used by the one-to-one CDS grader. */
		public int calls(final int arm){
			checkCalls(arm);
			return contig.calls[arm].size();
		}
		public Feature call(final int arm, final int index){
			checkCalls(arm);
			return contig.calls[arm].get(index);
		}
		/** Returns exactly the grader's matching, including circular geometry,
		 * without exposing mutable feature lists to the diagnostic reader.
		 */
		public MatchResult matchCalls(final int arm){
			checkCalls(arm);
			return PanelAnnotationGrade.match(contig.calls[arm], contig.truth, circularLength());
		}
		private void checkCalls(final int arm){
			require(owner.callState==2 && arm>=0 && arm<2,
				"Saved CDS views require a completely loaded before/after pair and explicit arm0or1");
		}
		/** Row-level start eligibility only: complete flags, present ID and one
		 * occurrence in THIS input. Different eligible rows can still disagree on
		 * their same-stop start; retain the saved grader's grouped ambiguity test.
		 */
		public boolean completeUniqueRow(final int index){
			final int id=contig.referenceIds.get(index);
			return id>=0 && contig.referenceFlags.get(index)==0 && owner.genome.idCounts.get(id)==1;
		}
		private final Reference owner;
		private final Contig contig;
	}

	static final class Genome{
		Genome(final Job job_){this(job_, null);}
		Genome(final Job job_, final NcrnaReferenceCalls projection_){
			job=job_;
			projection=projection_;
		}
		void read(final String path, final int arm){
			assert(arm>=-1 && arm<2) : "One original annotation or one of the two declared call arms is required";
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(path)){
				long lineNumber=0;
				for(byte[] row; (row=input.reader.nextLine())!=null;){
					lineNumber++;
					if(row.length==0){continue;}
					if(row[0]=='#'){
						if(Tools.startsWith(row, "##FASTA")){require(arm<0, "Saved calls cannot embed DNA"); break;}
						if(arm<0 && Tools.startsWith(row, "##sequence-region")){region(row);}
						continue;
					}
					p.set(row);
					require(p.terms()==9, "Expected GFF9 in "+path+" at line "+lineNumber);
					final Contig c=arm<0 ? contig(p, 0) : declaredQueryContig(p);
					if(arm<0 && p.termEquals("region", 2)){
						final String circular=value(p, "Is_circular=", ';');
						if(circular!=null){
							require((circular.equals("true") || circular.equals("false")) && p.parseInt(3)==1,
								"Circular topology needs an explicit full-origin region and true/false value");
							final boolean yes=circular.equals("true");
							final int stop=p.parseInt(4);
							require(!c.topologyKnown || c.circular==yes && c.regionStop==stop,
								"Reference full-region declarations conflict for "+c.name);
							c.topologyKnown=true;
							c.circular=yes;
							c.regionStop=stop;
						}
					}
					if(!p.termEquals("CDS", 2)){
						if(arm>=0){countRna(p, arm, c);}
						continue;
					}
					final int start=p.parseInt(3), stop=p.parseInt(4);
					final char strand=p.termEquals("+", 6) ? '+' : p.termEquals("-", 6) ? '-' : '?';
					final int phase=p.termEquals(".", 7) ? -1 : p.parseInt(7);
					require(start>0 && stop>=start && strand!='?' && phase>=-1 && phase<=2, "CDS coordinates, strand and phase must be explicit and valid");
					if(arm>=0){
						require(p.termEquals("BBTools", 1) && stop<=c.length, "Saved native CDS must remain inside its declared physical contig");
						c.calls[arm].add(new Feature(new Locus(c.name, "CDS", start, stop, strand), phase, lineNumber));
					}else{
						c.truth.add(new Feature(new Locus(c.name, "CDS", start, stop, strand), phase, lineNumber));
						referenceMetadata(p);
						if(referenceId>=0){idCounts.increment(referenceId, 1);}
						c.referenceIds.add(referenceId);
						c.referenceFlags.add(referenceUncertain ? 1 : 0);
						referenceRows++;
					}
				}
			}
		}

		void validateReference(){
			for(final Contig c:order){
				require(c.length>0, "Reference contig lacks a declared physical length: "+c.name);
				require(!c.topologyKnown || c.regionStop==c.length,
					"Reference full-region length disagrees with sequence-region for "+c.name);
				for(final Feature f:c.truth){
					// GFF circular parts can lie wholly beyond the origin; retain their raw row coordinates.
					require((long)f.locus.stop-f.locus.start+1<=c.length && (f.locus.stop<=c.length || c.circular),
						"CDS must span at most one circumference; coordinates beyond the physical contig need declared circular geometry: "
						+c.name+" start="+f.locus.start+" stop="+f.locus.stop+" length="+c.length+" circular="+c.circular);
				}
			}
		}

		private void region(final byte[] row){
			final int[] bounds=regionBounds;
			int terms=0;
			for(int i=0; i<row.length;){
				while(i<row.length && row[i]<=32){i++;}
				if(i==row.length){break;}
				require(terms<4, "sequence-region requires exactly a directive, seqid, first and last base");
				bounds[2*terms]=i;
				while(i<row.length && row[i]>32){i++;}
				bounds[2*terms+1]=i;
				terms++;
			}
			require(terms==4 && bounds[1]==17 && Parse.parseInt(row, bounds[4], bounds[5])==1,
				"The exact sequence-region directive must declare one physical contig from base one");
			final Contig c=contig(row, bounds[2], bounds[3]);
			final int length=Parse.parseInt(row, bounds[6], bounds[7]);
			require(length>0 && (c.length==0 || c.length==length), "Conflicting declared contig lengths");
			c.length=length;
		}

		private void countRna(final LineParser1 p, final int arm, final Contig contig){
			final int type=p.termEquals("tRNA", 2) ? 0 : p.termEquals("rRNA", 2) ? 1 : p.termEquals("RNA", 2) ? 2 : -1;
			require(type>=0 && p.termEquals("BBTools", 1), "Saved calls must retain a known native feature type/source");
			long range=attributeRange(p, "family:", ',');
			if(range<0){
				require(type!=2, "Generic native RNA needs an explicit family tag for the other-family count tripwire");
				p.setBounds(8);
				int end=p.a();
				while(end<p.b() && p.line()[end]!=',' && p.line()[end]!=';'){end++;}
				range=((long)p.a()<<32)|end;
			}
			final int from=(int)(range>>>32), to=(int)range;
			require(to>from, "A native RNA count requires a nonempty family/type label");
			final int id=symbols.intern(p.line(), from, to);
			final int key=Math.addExact(Math.multiplyExact(3, id), type);
			if(rnaSeen.add(key)){rnaOrder.add(key);}
			rnaCounts[arm].increment(key, 1);
			if(projection==null){return;}
			int projected=id;
			if(type==2){
				final long model=attributeRange(p, "model:", ',');
				require(model>=0, "Exact generic-RNA reports need a model for the same explicit projection as training");
				final Rule rule=projection.modelRule(p.line(), (int)(model>>>32), (int)model);
				require(rule.matchesLabel(p.line(), from, to), "Saved RNA family label conflicts with its explicit model projection");
				if(rule.include){return;}
				projected=rule.internFamily(symbols);
			}
			final int start=p.parseInt(3), stop=p.parseInt(4);
			final int strand=p.termEquals("+", 6) ? 0 : p.termEquals("-", 6) ? 1 : -1;
			require(start>0 && stop>=start && stop<=contig.length && strand>=0,
				"Exact saved RNA reports must retain a strand and inclusive coordinates inside the declared physical contig");
			if(contig.rna==null){contig.rna=new NcrnaRnaCallChanges.Reported();}
			final int projectedKey=Math.addExact(Math.multiplyExact(3, projected), type);
			contig.rna.add(Math.addExact(Math.multiplyExact(2, projectedKey), strand), start, stop, arm);
		}
		/** Completes one explicit row per input, including inputs with no RNA.
		 * Equal family totals cannot hide losses at different positions or strands.
		 * Coordinates are exact reports; an endpoint shift is loss+gain, not proof
		 * of biological gene loss. Each contig's primitive groups are traversed once.
		 */
		ByteBuilder rnaChanges(final ByteBuilder out){
			require(projection!=null && job!=null, "Exact saved RNA changes require an opted-in, fully bound input pair");
			long before=0, after=0, lost=0, gained=0;
			for(final Contig contig:order){
				if(contig.rna==null){continue;}
				final NcrnaRnaCallChanges.Reported rna=contig.rna;
				rna.finish();
				before=Math.addExact(before, rna.before());
				after=Math.addExact(after, rna.after());
				lost=Math.addExact(lost, rna.lost());
				gained=Math.addExact(gained, rna.gained());
			}
			return out.append(job.old.id).tab().append(before).tab().append(after).tab().append(lost).tab().append(gained)
				.tab().append(RNA_CHANGE_QUALIFICATION).nl();
		}

		void rnaCounts(final ByteBuilder out, final long[] changes){
			for(int i=0; i<rnaOrder.size; i++){
				final int key=rnaOrder.get(i), id=key/3;
				final byte[] name=symbols.bytes(id);
				final String label=new String(name, StandardCharsets.UTF_8);
				boolean clan=false;
				for(final String member:CLAN_LABELS){clan|=member.equals(label);}
				final int old=Tools.max(0, rnaCounts[0].get(key)), next=Tools.max(0, rnaCounts[1].get(key));
				if(!clan && old!=next){changes[OTHER_RNA_COUNT_CHANGE]++;}
				out.append(job.old.id).tab().append(RNA_TYPES[key%3]).tab().append(name).tab().append(old).tab().append(next)
					.tab().append(clan).tab().append(old==next).nl();
			}
		}

		private Contig contig(final LineParser1 p, final int field){p.setBounds(field); return contig(p.line(), p.a(), p.b());}
		private Contig declaredQueryContig(final LineParser1 p){
			p.setBounds(0);
			return declaredQueryContig(p.line(), p.a(), p.b());
		}
		private Contig declaredQueryContig(final byte[] row, final int begin, final int end){
			final int from=accessionStart(row, begin, end);
			final int id=symbols.find(row, from, end);
			final Contig c=id<0 ? null : contigs.get(id);
			require(c!=null && c.length>0, "Saved query contig lacks a declared reference identity/length");
			if(from>begin){
				final int qualified=symbols.find(row, begin, end);
				require(qualified>=0 && qualified==c.qualifiedId, "A query taxID-qualified contig must equal the original annotation's qualified identity");
			}
			return c;
		}
		private Contig contig(final byte[] row, final int from, final int to){
			require(to>from && !(to-from==1 && row[from]=='.'), "A GFF contig must have a nonempty literal identity");
			final int accession=accessionStart(row, from, to);
			final int id=symbols.intern(row, accession, to);
			Contig c=contigs.get(id);
			if(c==null){
				final String name=new String(symbols.bytes(id), StandardCharsets.UTF_8);
				require(name.indexOf('\uFFFD')<0, "Malformed UTF-8 must not alter contig identity during CDS grading");
				c=new Contig(name);
				contigs.put(id, c);
				order.add(c);
			}
			if(accession>from){
				final int qualified=symbols.intern(row, from, to);
				require(c.qualifiedId<0 || c.qualifiedId==qualified, "Two taxID-qualified identities cannot share one declared physical contig: "+c.name);
				c.qualifiedId=qualified;
			}
			return c;
		}
		/** Older bacteria4 GFFs retain bare accessions in sequence-region but add
		 * tid|decimal_taxid| to feature seqids. Remove ONLY that recognized prefix
		 * within this one source; preserve a single exact qualified alias above.
		 * No accession/version trimming, cross-source join or inferred length. */
		private static int accessionStart(final byte[] row, final int from, final int to){
			if(to-from<7 || !prefix(row, from, to, "tid|")){return from;}
			int pos=from+4;
			final int first=pos;
			while(pos<to && row[pos]>='0' && row[pos]<='9'){pos++;}
			return pos>first && pos<to-1 && row[pos]=='|' ? pos+1 : from;
		}
		/** Scan source attributes once; retain only needed IDs, with no String per CDS field. */
		private void referenceMetadata(final LineParser1 p){
			p.setBounds(8);
			final byte[] row=p.line();
			referenceId=-1;
			referenceUncertain=false;
			for(int from=p.a(), end; from<p.b(); from=end+1){
				end=from;
				while(end<p.b() && row[end]!=';'){end++;}
				if(prefix(row, from, end, "ID=")){
					require(referenceId<0, "A CDS row cannot have two competing ID attributes");
					if(end>from+3){referenceId=symbols.intern(row, from+3, end);}
				}else if(end-from==12 && prefix(row, from, end, "partial=true")
					|| prefix(row, from, end, "start_range=") || prefix(row, from, end, "end_range=")
					|| prefix(row, from, end, "exception=") || prefix(row, from, end, "transl_except=")
					|| end-from==11 && prefix(row, from, end, "pseudo=true")){
					referenceUncertain=true;
				}
			}
		}
		private final Job job;
		private final NcrnaReferenceCalls projection;
		private final Symbols symbols=new Symbols();
		private final IntObjectMap<Contig> contigs=new IntObjectMap<Contig>(16, Contig.class);
		private final ArrayList<Contig> order=new ArrayList<Contig>();
		private final IntHashMap idCounts=new IntHashMap(16);
		private final IntHashMap[] rnaCounts={new IntHashMap(16), new IntHashMap(16)};
		private final IntHashSet rnaSeen=new IntHashSet(16);
		private final IntList rnaOrder=new IntList(16);
		private final int[] regionBounds=new int[8];
		private int referenceId;
		private boolean referenceUncertain;
		private long referenceRows;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private static LongIntMap index(final ArrayList<Feature> calls, final int circle){
		final LongIntMap map=new LongIntMap(Tools.max(16, calls.size()));
		for(int i=0; i<calls.size(); i++){
			final long key=stopKey(calls.get(i), circle);
			map.put(key, map.get(key)==-1 ? i : -2);
		}
		return map;
	}
	private static long stopKey(final Feature f, final int circle){return ((long)physical(f.locus.strand=='+' ? f.locus.stop : f.locus.start, circle)<<1)|(f.locus.strand=='-' ? 1 : 0);}
	private static int start(final Feature f){return f.locus.strand=='+' ? f.locus.start : f.locus.stop;}
	private static int physical(final int coordinate, final int circle){return circle==0 ? coordinate : (coordinate-1)%circle+1;}
	/** Positive means downstream/too short. A shared biological stop anchors the
	 * full sense interval, including a wrapped reference; nearest-half-circle
	 * offsets would give the wrong shortening sign for a very long CDS. */
	private static int startDelta(final Feature call, final Feature reference, final int circle){
		assert(stopKey(call, circle)==stopKey(reference, circle)) : "Length difference measures a start shift only with the same strand and biological stop";
		return (reference.locus.stop-reference.locus.start)-(call.locus.stop-call.locus.start);
	}
	private static void coordinates(final ByteBuilder out, final Feature f){
		out.tab();
		if(f==null){out.append("NA\tNA");}else{out.append(f.locus.start).tab().append(f.locus.stop);}
	}
	private static String value(final LineParser1 p, final String prefix, final char delimiter){
		final long range=attributeRange(p, prefix, delimiter);
		return range<0 ? null : new String(p.line(), (int)(range>>>32), (int)range-(int)(range>>>32), StandardCharsets.UTF_8);
	}
	private static long attributeRange(final LineParser1 p, final String prefix, final char delimiter){
		p.setBounds(8);
		final byte[] row=p.line();
		for(int from=p.a(), end; from<p.b(); from=end+1){
			end=from;
			while(end<p.b() && row[end]!=delimiter && row[end]!=';'){end++;}
			if(prefix(row, from, end, prefix)){return ((long)(from+prefix.length())<<32)|end;}
		}
		return -1;
	}
	private static boolean prefix(final byte[] row, final int from, final int to, final String prefix){
		if(to-from<prefix.length()){return false;}
		for(int i=0; i<prefix.length(); i++){if(row[from+i]!=prefix.charAt(i)){return false;}}
		return true;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private static final class Job{
		Job(final Binding old_, final Binding next_){old=old_; next=next_;}
		final Binding old, next;
		String original;
	}
	static final class Contig{
		@SuppressWarnings("unchecked")
		Contig(final String name_){name=name_; calls=new ArrayList[]{new ArrayList<Feature>(), new ArrayList<Feature>()};}
		final String name;
		int length, regionStop, qualifiedId=-1;
		boolean topologyKnown, circular;
		final ArrayList<Feature> truth=new ArrayList<Feature>();
		final ArrayList<Feature>[] calls;
		final IntList referenceIds=new IntList(16), referenceFlags=new IntList(16);
		NcrnaRnaCallChanges.Reported rna;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String before, after, corpus, out, familymap;
	private Rule[] rnaRules;
	private int expectedInputs=-1, completed;
	private boolean errorState, rnaCallChanges;
	private String[] suffixes=SUFFIXES;
	private final PrintStream outstream;
	private final ArrayList<Job> jobs=new ArrayList<Job>();
	private final AtomicInteger nextInput=new AtomicInteger();// One coarse claim per genome, never per feature.
	private final ReadWriteLock lock=new ReentrantReadWriteLock();
	private final long[][] totals=new long[2][PanelAnnotationGrade.NSTAT];
	private final long[] changeTotals=new long[CHANGE_NAMES.length];
	private static final int START_CHANGED=0, STOP_LOST=1, STOP_ADDED=2, RELIABLE_START_PAIR=3,
		START_EXACT_BEFORE=4, START_EXACT_AFTER=5, START_ERROR_WORSE=6, SHORTENED_PAST_START=7,
		REFERENCE_LOST=8, REFERENCE_GAINED=9, AMBIGUOUS_CALL_STOP=10, OTHER_RNA_COUNT_CHANGE=11;
	private static final String[] CHANGE_NAMES={"same_stop_start_changed", "stop_not_called_rows", "added_stop_rows", "reliable_changed_start_pairs",
		"same_stop_changed_start_exact_before", "same_stop_changed_start_exact_after", "changed_start_error_worsened", "shortened_past_annotated_start",
		"annotated_CDS_detection_lost", "annotated_CDS_detection_gained", "ambiguous_call_stop_rows", "other_RNA_family_count_changes"};
	private static final String[] STAT_NAMES=PanelAnnotationGrade.STAT_HEADER.split("\t");
	private static final String[] ARM={"before", "after"}, RNA_TYPES={"tRNA", "rRNA", "RNA"};
	private static final String[] CLAN_LABELS={NcrnaFamilyRegistry.displayName("cobalamin_rf00174"),
		NcrnaFamilyRegistry.displayName("cobalamin_rf01482"), NcrnaFamilyRegistry.displayName("cobalamin_rf01689")};
	private static final String[] SUFFIXES={".summary.tsv", ".inputs.tsv", ".changes.tsv", ".rna_counts.tsv", ".sources.tsv"};
	private static final String QUALIFICATION="SOURCE_GFF_FEATURE_ROWS_NOT_INDEPENDENT_BIOLOGICAL_TRUTH";
	private static final String INPUT_HEADER="input_id\tarm\t"+PanelAnnotationGrade.STAT_HEADER+"\tqualification";
	private static final String CHANGE_HEADER="input_id\tcontig\tstrand\tevent\tbefore_start1\tbefore_stop1\tafter_start1\tafter_stop1\tsame_stop_reference_rows\tannotation_status\tbefore_start_error_downstream_bp\tafter_start_error_downstream_bp\tcircular_length\tqualification";
	private static final String RNA_HEADER="input_id\ttype\tfamily\tbefore\tafter\tcobalamin_clan\tcount_unchanged";
	public static final String RNA_CHANGE_HEADER="input_id\tbefore_calls\tafter_calls\tlost_calls\tgained_calls\tqualification";
	public static final String RNA_CHANGE_QUALIFICATION="EXACT_REPORTED_OTHER_RNA_MULTISET;INCLUDED_PROJECTION_EXCLUDED;ENDPOINT_SHIFT_IS_LOSS_PLUS_GAIN";
	private static final String SOURCE_HEADER="input_id\tfasta\toriginal_gff\treference_CDS_rows\tcontigs_seen\tqualification";
	private static final String CORPUS_HEADER="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
}
