package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.nio.file.DirectoryStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.GeneticCode;
import dna.Data;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prot.ProteinFamilyAssigner;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;

/** Annotates existing contiguous CDS calls without changing caller coordinates.
 * FASTA/GFF extraction and ordered writing are serial; workers align proteins.
 * Uses the accepted expanded-family search rule, not a biochemical-function test.
 * @author Keqing
 */
public final class AnnotateGenes implements Accumulator<AnnotateGenes.ProcessThread> {

	public static void main(String[] args) throws IOException{
		final Timer timer=new Timer(); final AnnotateGenes tool=new AnnotateGenes(args);
		try{tool.process(timer);}finally{Shared.closeStream(tool.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public AnnotateGenes(String[] args) throws IOException{
		final PreParser pp=new PreParser(args, getClass(), false); outstream=pp.outstream;
		final Parser parser=new Parser(); parser.overwrite=true;
		boolean modelsSpecified=false, labelsSpecified=false;
		for(String arg : pp.args){
			final int eq=arg.indexOf('=');
			final String key=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase(Locale.ROOT);
			String value=eq<0 ? null : arg.substring(eq+1);
			if(value!=null && (value.isEmpty() || "null".equalsIgnoreCase(value))){value=null;}
			if(eq<0 && FileFormat.hasFastaExtension(arg)){
				require(parser.in1==null, "Multiple FASTA inputs: "+parser.in1+" and "+arg);
				parser.in1=arg;
			}else if(eq<0 && FileFormat.hasGffExtension(arg)){
				require(gff==null, "Multiple GFF inputs: "+gff+" and "+arg);
				gff=arg;
			}else if(key.equals("in") || key.equals("input") || key.equals("in1") || key.equals("input1")){
				require(value!=null, "FASTA input requires a filename: "+arg);
				require(parser.in1==null, "Multiple FASTA inputs: "+parser.in1+" and "+value);
				parser.parse(arg, key, value);
			}else if(key.equals("gff")){
				require(value!=null, "GFF input requires a filename: "+arg);
				require(gff==null, "Multiple GFF inputs: "+gff+" and "+value);
				gff=value;
			}else if(key.equals("models")){models=value; modelsSpecified=true;}
			else if(key.equals("labels")){labelFile=value; labelsSpecified=true;}
			else if(key.equals("evidence")){evidence=value;}
			else if(key.equals("transl_table") || key.equals("table")){table=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown AnnotateGenes parameter: "+arg);}
		}
		in=parser.in1; out=parser.out1; overwrite=parser.overwrite;
		require(in!=null && gff!=null && out!=null, "Required: genome.fna calls.gff out=annotated.gff (or in=genome.fna gff=calls.gff)");
		require(parser.in2==null && parser.out2==null && !parser.append && parser.maxReads<0,
			"AnnotateGenes requires whole unpaired FASTA input and does not support append, in2, out2 or maxreads");
		final String requestedModels=models;
		models=models==null ? null : Data.findPath(models, false);
		require(models!=null && new File(models).isDirectory(), "Model package directory not found: "+requestedModels
			+". Extract the BBGAP reference tar into BBTools/resources/ to create bbgap/models/ and bbgap/labels/, or specify models=DIR.");
		if(labelsSpecified){
			if(labelFile!=null){
				final String requestedLabels=labelFile;
				labelFile=Data.findPath(labelFile, false);
				require(labelFile!=null, "Family labels not found: "+requestedLabels);
			}
		}else{
			// Resolve labels beside the selected models, never from another resource installation.
			final File bundle=new File(models).getCanonicalFile().getParentFile();
			final File compressed=new File(bundle, "labels/family_labels.tsv.gz");
			final File raw=new File(bundle, "labels/family_labels.tsv");
			labelFile=compressed.isFile() ? compressed.getPath() : raw.isFile() ? raw.getPath() : null;
			require(modelsSpecified || labelFile!=null, "The installed BBGAP package is missing labels/family_labels.tsv.gz. "
				+"Extract the complete reference tar into BBTools/resources/, or use labels=null for family assignments without names.");
		}
		GeneticCode.forTable(table);
		require(Shared.threads()>0, "At least one alignment thread is required");
		require(Tools.testInputFiles(false, true, in, gff, labelFile), "An AnnotateGenes input is unreadable");
		require(Tools.testOutputFiles(overwrite, false, false, out, evidence), "Cannot write AnnotateGenes output");
		validatePaths();
	}

	/** Canonical paths catch aliases before any writer can truncate an input. */
	private void validatePaths() throws IOException{
		final Path modelRoot=new File(models).getCanonicalFile().toPath();
		final String[] sources={in, gff, labelFile}, targets={out, evidence};
		int standardInputs=0, standardOutputs=0;
		for(String source : sources){
			if(source!=null && FileFormat.testInput(source, FileFormat.TEXT, null, false, false).stdin()){standardInputs++;}
		}
		for(String target : targets){
			if(target!=null && FileFormat.testOutput(target, FileFormat.TEXT, null, false, overwrite, false, false).stdout()){standardOutputs++;}
		}
		require(standardInputs<=1 && standardOutputs<=1, "Multiple inputs or outputs alias a standard stream");
		for(int i=0; i<targets.length; i++){
			if(targets[i]==null){continue;}
			final File target=new File(targets[i]).getCanonicalFile();
			require(!target.toPath().startsWith(modelRoot), "Outputs must not overwrite the immutable model package: "+target);
			if(Files.isDirectory(modelRoot)){
				try(DirectoryStream<Path> members=Files.newDirectoryStream(modelRoot)){
					for(Path member : members){
						require(!sameFile(target, member.toFile()), "Output aliases a model package member: "+target);
					}
				}
			}
			for(String source : sources){
				if(source!=null){require(!sameFile(target, new File(source).getCanonicalFile()), "Output aliases an input: "+target);}
			}
			for(int j=0; j<i; j++){
				if(targets[j]!=null){require(!sameFile(target, new File(targets[j]).getCanonicalFile()), "Output files alias one another: "+target);}
			}
		}
		require(!sameFile(new File(in).getCanonicalFile(), new File(gff).getCanonicalFile()), "FASTA and GFF must be separate inputs");
	}
	private static boolean sameFile(final File a, final File b) throws IOException{
		return a.equals(b) || (a.exists() && b.exists() && Files.isSameFile(a.toPath(), b.toPath()));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer) throws IOException{
		final HashMap<String,byte[]> scaffolds=readScaffolds(in);
		final AnnotationGff document=AnnotationGff.load(gff, scaffolds, table);
		outstream.println("Loaded "+scaffolds.size()+" contigs and "+document.cds.size()+" CDS features.");
		final ProteinFamilyAssigner library=ProteinFamilyAssigner.load(Paths.get(models));
		final FamilyLabels labels=labelFile==null ? null : FamilyLabels.load(labelFile, library);
		final AtomicInteger next=new AtomicInteger(); final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final int count=Math.min(Shared.threads(), Math.max(1, document.cds.size()));
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(count);
		for(int i=0; i<count; i++){workers.add(new ProcessThread(document.cds, library, next, failure));}
		final boolean completed=ThreadWaiter.startAndWait(workers, this);
		if(!completed || failure.get()!=null){throw new RuntimeException("Protein alignment failed; annotation output was not opened", failure.get());}
		require(assigned+noHit+unsupported==document.cds.size(), "Every CDS must be assigned, no_hit or unsupported");
		require(scored+lengthRejected==50*(assigned+noHit), "Every searched CDS must account for all 50 shortlisted candidates");
		write(document, library, labels);
		timer.stop();
		outstream.println("CDS: "+document.cds.size()+"; assigned: "+assigned+"; no_hit: "+noHit+"; unsupported: "+unsupported);
		outstream.println("Scored candidates: "+scored+"; length rejected: "+lengthRejected+"; time: "+timer);
	}

	/** Existing FASTA reader, with synchronous conversion and no alignment workers. */
	static HashMap<String,byte[]> readScaffolds(final String path){
		final FileFormat format=FileFormat.testInput(path, FileFormat.FASTA, null, true, false);
		require(format.fasta() && !format.amino() && !Shared.AMINO_IN, "Expected nucleotide FASTA: "+path);
		final Streamer input=StreamerFactory.makeStreamer(format, null, 0, true, -1, false, true, 0);
		final HashMap<String,byte[]> sequences=new HashMap<String,byte[]>();
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null && batch.size()>0; batch=input.nextList()){
				for(Read read : batch){
					require(read.mate==null && read.id!=null && read.bases!=null && read.bases.length>0, "Expected a named, nonempty, unpaired FASTA contig");
					int end=0; while(end<read.id.length() && !Character.isWhitespace(read.id.charAt(end))){end++;}
					require(end>0, "Empty FASTA identifier"); final String id=read.id.substring(0, end);
					require(sequences.putIfAbsent(id, read.bases)==null, "Duplicate FASTA contig ID: "+id);
				}
			}
		}finally{input.close();}
		require(!input.errorState() && !sequences.isEmpty(), "FASTA read failed or contained no contigs: "+path);
		return sequences;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final AnnotationGff document, final ProteinFamilyAssigner library, final FamilyLabels labels){
		ByteStreamWriter calls=null, details=null;
		boolean error=false;
		try{
			calls=ByteStreamWriter.makeBSW(out, overwrite, false, true);
			if(evidence!=null){
				details=ByteStreamWriter.makeBSW(evidence, overwrite, false, true);
				details.println("#rule\t"+ProteinFamilyAssigner.RULE);
				details.println("#model_sha80\t"+library.modelPin());
				details.println("#family_identity_sha80\t"+library.identityPin());
				details.println("#consensus_sha80\t"+library.consensusPin());
				if(labels!=null){details.println("#label_version\t"+labels.version);}
				details.println("line\tseqid\tstart\tend\tstrand\tphase\ttranslation_table\tpartial_5\tpartial_3\tstatus\treason\tfamily_id\trep_id\tscore64\tidentity_pct\tcoverage_query\tcoverage_core\tmutual_coverage\tquery_length\tref_start\tref_end\tscored_candidates\tlength_rejected");
			}
			final ByteBuilder row=new ByteBuilder(1024); int next=0;
			for(int i=0; i<document.lines.size(); i++){
				if(next>=document.cds.size() || document.cds.get(next).line!=i){calls.println(document.lines.get(i)); continue;}
				final AnnotationGff.Cds cds=document.cds.get(next++);
				annotatedRow(row, document.lines.get(i), cds, library.modelPin(), labels); calls.println(row);
				if(details!=null){evidenceRow(row, cds); details.println(row);}
			}
			assert(next==document.cds.size()) : "Every ordered CDS result must be written to its original GFF row";
		}finally{
			if(calls!=null){error|=calls.poisonAndWait();}
			if(details!=null){error|=details.poisonAndWait();}
		}
		require(!error, "Annotation output failed; output files may be incomplete");
	}

	static void annotatedRow(final ByteBuilder row, final byte[] raw, final AnnotationGff.Cds cds,
			final String modelPin, final FamilyLabels labels){
		AnnotationGff.beginAnnotation(row, raw).append("bbgap_status=").append(cds.status());
		row.append(";bbgap_translation_table=").append(cds.translation.table);
		row.append(";bbgap_rule=").append(ProteinFamilyAssigner.RULE);
		row.append(";bbgap_package=").append(modelPin);
		if(labels!=null){row.append(";bbgap_label_version="); AnnotationGff.escaped(row, labels.version);}
		if(cds.reason!=null){row.append(";bbgap_reason=").append(cds.reason); return;}
		final ProteinFamilyAssigner.Hit hit=cds.hit;
		if(!hit.assigned()){return;}
		row.append(";bbgap_family_id=").append(hit.familyId).append(";bbgap_rep_id="); AnnotationGff.escaped(row, hit.representative);
		row.append(";bbgap_score64=").append(hit.score64).append(";bbgap_identity=").append(hit.identity, 8);
		row.append(";bbgap_query_cov=").append(hit.queryCoverage, 8).append(";bbgap_core_cov=").append(hit.coreCoverage, 8);
		appendFamilyLabel(row, labels, hit.rank);
	}

	/** Appends selected label fields without clearing existing GFF attributes. */
	static void appendFamilyLabel(final ByteBuilder row, final FamilyLabels labels, final int rank){
		assert(row!=null) : "Selected label fields must append to the existing annotated GFF row";
		final FamilyLabels.Label label=labels==null ? null : labels.get(rank);
		row.append(";bbgap_label_class=").append(label==null ? "unnamed" : label.decision);
		if(label!=null && label.product!=null){row.append(";bbgap_product="); AnnotationGff.escaped(row, label.product);}
		if(label!=null && label.display!=null){
			appendAttribute(row, "bbgap_label_display", label.display);
			appendAttribute(row, "bbgap_label_source", label.source);
			appendAttribute(row, "bbgap_label_support", label.support);
			appendAttribute(row, "bbgap_label_coverage", label.coverage);
			appendAttribute(row, "bbgap_label_score", label.score);
			appendAttribute(row, "bbgap_label_product_text", label.productText);
			appendAttribute(row, "bbgap_label_review", label.review);
			appendAttribute(row, "bbgap_label_hit", label.hit);
		}
	}

	private static void appendAttribute(final ByteBuilder row, final String key,
			final String value){
		row.append(';').append(key).append('=');
		AnnotationGff.escaped(row, value==null ? "" : value);
	}

	private static void evidenceRow(final ByteBuilder row, final AnnotationGff.Cds cds){
		row.clear().append(cds.line+1).tab().append(cds.feature.seqid).tab().append(cds.feature.start).tab().append(cds.feature.stop)
			.tab().append(cds.feature.strand==0 ? "+" : cds.feature.strand==1 ? "-" : ".").tab().append(cds.feature.phase)
			.tab().append(cds.translation.table).tab().append(cds.translation.partialFive).tab().append(cds.translation.partialThree)
			.tab().append(cds.status()).tab().append(cds.reason==null ? "." : cds.reason);
		final ProteinFamilyAssigner.Hit hit=cds.hit;
		if(hit==null || !hit.assigned()){
			row.append("\t.\t.\t.\t.\t.\t.\t.\t").append(cds.translation.peptide==null ? 0 : cds.translation.peptide.length)
				.append("\t.\t.\t").append(hit==null ? 0 : hit.scored).tab().append(hit==null ? 0 : hit.lengthRejected);
		}else{
			row.tab().append(hit.familyId).tab().append(hit.representative).tab().append(hit.score64).tab().append(hit.identity, 8)
				.tab().append(hit.queryCoverage, 8).tab().append(hit.coreCoverage, 8).tab().append(hit.mutualCoverage, 8)
				.tab().append(hit.queryLength).tab().append(hit.referenceStart).tab().append(hit.referenceEnd)
				.tab().append(hit.scored).tab().append(hit.lengthRejected);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Worker Threads        ----------------*/
	/*--------------------------------------------------------------*/

	static final class ProcessThread extends Thread {
		ProcessThread(ArrayList<AnnotationGff.Cds> records_, ProteinFamilyAssigner library,
				AtomicInteger next_, AtomicReference<Throwable> failure_){
			records=records_; worker=library.newWorker(); next=next_; failure=failure_;
		}
		@Override
		public void run(){
			try{
				for(int from=next.getAndAdd(8); from<records.size() && failure.get()==null; from=next.getAndAdd(8)){
					for(int i=from, end=Math.min(from+8, records.size()); i<end; i++){
						final AnnotationGff.Cds cds=records.get(i);
						if(cds.reason!=null){unsupported++; continue;}
						cds.hit=worker.assign("gff_line_"+(cds.line+1), cds.translation.peptide);
						if(cds.hit.assigned()){assigned++;}else{noHit++;}
						scored+=cds.hit.scored; lengthRejected+=cds.hit.lengthRejected;
					}
				}
				success=failure.get()==null;
			}catch(Throwable problem){failure.compareAndSet(null, problem);}
		}
		final ArrayList<AnnotationGff.Cds> records;
		final ProteinFamilyAssigner.Worker worker;
		final AtomicInteger next;
		final AtomicReference<Throwable> failure;
		long assigned, noHit, unsupported, scored, lengthRejected;
		boolean success;
	}

	@Override
	public void accumulate(ProcessThread worker){
		assigned+=worker.assigned; noHit+=worker.noHit; unsupported+=worker.unsupported;
		scored+=worker.scored; lengthRejected+=worker.lengthRejected; errorState|=!worker.success;
	}
	@Override
	public boolean success(){return !errorState;}
	@Override
	public ReadWriteLock rwlock(){return lock;}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private final ReadWriteLock lock=new ReentrantReadWriteLock();
	private final PrintStream outstream;
	private String in, gff, out, labelFile, evidence;
	private String models="?bbgap/models";
	private int table=11;
	private boolean overwrite=true, errorState=false;
	private long assigned, noHit, unsupported, scored, lengthRejected;
}
