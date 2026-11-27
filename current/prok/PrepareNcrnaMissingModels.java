package prok;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import tax.TaxTree;

/** Resolves saved inclusive FN references to admitted original gene sequences.
 * An explicit protected-prior supplement independently binds historical authority
 * and original assembly bytes without changing the native construction population.
 * No new grading, clustering, taxonomy inference or caller runs. Unresolved
 * originals remain in the audit and block FASTA publication rather than shrink
 * the declared primary population. Optional gradeprefixes joins disjoint input
 * rosters and keeps reference identity as the (input_id, reference_id) pair.
 * @author Brian Bushnell, Diona
 */
public final class PrepareNcrnaMissingModels{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PrepareNcrnaMissingModels.class, false);
		try{
			boolean overlay=false;
			for(final String arg:pp.args){if(Parse.splitOnFirst(arg, '=')[0].equalsIgnoreCase("overlayassignments")){overlay=true;}}
			if(overlay){new Overlay(pp.args).process();}
			else{new PrepareNcrnaMissingModels(pp.args).process();}
		}
		finally{Shared.closeStream(pp.outstream);}
	}

	private PrepareNcrnaMissingModels(final String[] args){
		// This extractor publishes original symbols; establish and enforce raw input
		// even when another native tool previously changed process-wide defaults.
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=false;
		Read.IUPAC_TO_N=false;
		Read.U_TO_T=false;
		Read.DOT_DASH_X_TO_N=false;
		Read.LOWER_CASE_TO_N=false;
		Read.FIX_HEADER=false;
		Read.JUNK_MODE=Read.CRASH_JUNK;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("sourceids")){sourceIds=b;}
			else if(a.equals("protectedpriors")){priorBindings=b;}
			else if(a.equals("priorfasta")){priorFasta=b;}
			else if(a.equals("priorloci")){priorLoci=b;}
			else if(a.equals("priorprovenance")){priorProvenance=b;}
			else if(a.equals("priorwitnesses")){priorWitnesses=b;}
			else if(a.equals("priorpanel")){priorPanel=b;}
			else if(a.equals("expectedpriors")){expectedPriors=Parse.parseIntKMG(b);}
			else if(a.equals("expectedpriororiginals")){expectedPriorOriginals=Parse.parseIntKMG(b);}
			else if(a.equals("gradeprefixes")){gradePrefixes=b.split(",", -1);}
			else if(a.equals("manifest")){manifest=b;}
			else if(a.equals("fasta")){fastas=b.split(",", -1);}
			else if(a.equals("expected")){expected=Parse.parseIntKMG(b);}
			else if(a.equals("expectedoriginals")){
				final String[] counts=b.split(",", -1);
				require(counts.length==3, "Declare all three original member counts");
				originals=new int[3];
				for(int i=0; i<3; i++){originals[i]=Parse.parseIntKMG(counts[i]);}
			}else{require(parser.parse(arg, a, b), "Unknown missing-model preparation option: "+arg);}
		}
		grade=parser.in1;
		out=parser.out1;
		require(!Shared.TRIM_READ_DESCRIPTION && !Read.TO_UPPER_CASE && !Read.IUPAC_TO_N && !Read.U_TO_T
			&& !Read.DOT_DASH_X_TO_N && !Read.LOWER_CASE_TO_N && !Read.FIX_HEADER && Read.JUNK_MODE==Read.CRASH_JUNK,
			"Original sequence/header preservation requires trd=f tuc=f itn=f utot=f dotdashxton=f lctn=f and no header/junk repair");
		require((grade!=null ^ gradePrefixes!=null) && sourceIds!=null && manifest!=null && fastas!=null && fastas.length==3
			&& originals!=null && expected>0 && out!=null && parser.in2==null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Bind either in=single_raw_grade or gradeprefixes=disjoint_prefixes, sourceids, formal manifest, three FASTAs, exact counts and a fresh output prefix; no sampling");
		if(gradePrefixes!=null){
			for(final String prefix:gradePrefixes){require(!prefix.isEmpty(), "Every grade partition needs its explicit prefix");}
		}
		for(int i=0; i<3; i++){require(originals[i]>0 && !fastas[i].isEmpty(), "Each declared source member needs its complete FASTA/count");}
		final boolean anyPrior=priorBindings!=null || priorFasta!=null || priorLoci!=null || priorProvenance!=null
			|| priorWitnesses!=null || priorPanel!=null || expectedPriors!=-1 || expectedPriorOriginals!=-1;
		require(!anyPrior || priorBindings!=null && priorFasta!=null && priorLoci!=null && priorProvenance!=null
			&& priorWitnesses!=null && priorPanel!=null && expectedPriors>0 && expectedPriorOriginals>=expectedPriors,
			"Protected-prior recovery requires all six prior inputs and both exact positive counts together");
		RrnaResourceIO.fresh(out, ".references.tsv", ".summary.tsv", ".success.txt",
			".RF00174.fa", ".RF01482.fa", ".RF01689.fa");
		if(anyPrior){RrnaResourceIO.fresh(out, ".protected_priors.tsv");}
		Shared.setThreads(1);
		fileIO.ReadWrite.setZipThreads(1);
		Shared.TRIM_READ_DESCRIPTION=false;
	}

	private void process(){
		readGrade();
		if(priorBindings!=null){
			priors=readProtectedPriors(priorBindings, expectedPriors);
			for(final ProtectedPrior prior:priors){
				prior.ref=byCoordinates.get(key(prior.input, prior.contig, prior.coreStart, prior.coreStop, prior.strand, prior.member));
				require(prior.ref!=null && prior.ref.id.equals(prior.reference), "Protected prior must name an exact selected FN reference: "+prior.id);
				priorById.put(prior.id, prior);
				priorInputs.add(prior.input);
				priorCoordinates.put(key(prior.input, prior.contig, prior.start, prior.stop, prior.strand, prior.member), prior);
			}
		}
		readManifest();
		readSources();
		for(int member=0; member<3; member++){readFasta(member);}
		if(priors!=null){recoverProtectedPriors();}
		int unresolved=0;
		for(final Reference ref:references){if(!ref.state().equals("RESOLVED")){unresolved++;}}
		writeAudit(unresolved);
		require(unresolved==0, "Unresolved primary references="+unresolved+"; complete audit retained, no cluster FASTAs released");
		if(priors!=null){
			try(Output output=new Output(out+".protected_priors.tsv")){
				output.writer.println(PRIOR_HEADER);
				for(final ProtectedPrior prior:priors){output.writer.println(prior.row);}
			}
		}
		assert(sequenceOrder.size()==references.size()) : "FASTAs must publish exactly one original vote per audited reference; references="+references.size()+" sequences="+sequenceOrder.size();
		for(int member=0; member<3; member++){
			try(Output output=new Output(out+"."+MEMBERS[member]+".fa")){
				for(final Reference ref:sequenceOrder){
					if(ref.member==member){output.writer.print(b.clear().append('>').append(ref.header).nl().append(ref.bases).nl());}
				}
			}
		}
		try(Output output=new Output(out+".success.txt")){
			output.writer.println("NCRNA_MISSING_MODELS_PREPARED references="+references.size()+" original_votes="+sequenceOrder.size());
		}
		System.err.println("NCRNA_MISSING_MODELS_PREPARED references="+references.size()+" original_votes="+sequenceOrder.size());
	}

	/*--------------------------------------------------------------*/
	/*----------------      Exact Saved Bindings     ----------------*/
	/*--------------------------------------------------------------*/

	private void readGrade(){
		assert(references.isEmpty() && byCoordinates.size()==0 && byInput.size()==0) : "The saved grade is consumed once; stale bindings could alias two primary populations";
		final ObjectSet<String> identities=new ObjectSet<String>(String.class);
		if(gradePrefixes==null){readGrade(grade, identities, null);}
		else{
			final ObjectSet<Target> allInputs=new ObjectSet<Target>(Target.class);
			for(final String prefix:gradePrefixes){
				final ObjectSet<Target> partitionInputs=readPartitionInputs(prefix+".inputs.tsv", allInputs);
				readGrade(prefix+".overlap.raw.tsv", identities, partitionInputs);
			}
		}
		require(references.size()==expected, "FN denominator changed: observed="+references.size()+" expected="+expected);
	}

	/** Include inputs with no raw rows or no FN, so an overlapping empty partition
	 * cannot pass by contributing no selected references. Keys own their bytes. */
	private ObjectSet<Target> readPartitionInputs(final String path, final ObjectSet<Target> allInputs){
		final ObjectSet<Target> partition=new ObjectSet<Target>(Target.class);
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, GradeNcrnaReference.INPUT_HEADER, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==10 && !p.termEquals("", 0) && !p.termEquals(".", 0), "Each grade roster row needs its original nonempty input ID");
				rosterProbe.set(p, 0);
				final Target id=rosterProbe.copy();
				require(partition.add(id) && allInputs.add(id), "Grade partitions repeat an input ID; empty and zero-FN inputs also require disjoint ownership");
			}
		}
		require(partition.size()>0, "Each grade partition requires its complete nonempty input roster");
		return partition;
	}

	/** Single-raw mode keeps its original identity rule. Explicit partition mode
	 * qualifies reference IDs by input, retaining both original fields in audit. */
	private void readGrade(final String path, final ObjectSet<String> identities, final ObjectSet<Target> partitionInputs){
		assert((partitionInputs==null)==(gradePrefixes==null)) : "Only explicit grade partitions supply roster-qualified reference identities";
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, GRADE, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==18, "Raw grade requires exactly eighteen fields");
				if(partitionInputs!=null){
					rosterProbe.set(p, 0);
					require(partitionInputs.contains(rosterProbe), "Raw grade row names an input outside its declared partition");
				}
				if(!p.termEquals("REFERENCE", 1) || !p.termEquals("FN", 7)){continue;}
				final String id=p.parseString(8);
				final String inputId=p.parseString(0);
				final String contig=p.parseString(3);
				final String strand=p.parseString(6);
				final int start=p.parseInt(4);
				final int stop=p.parseInt(5);
				final int member=member(p.parseString(15));
				require(p.termEquals("cobalamin", 2) && !id.isEmpty() && !id.equals(".")
					&& identities.add(partitionInputs==null ? id : inputId+'\t'+id)
					&& !inputId.isEmpty() && !contig.isEmpty() && start>0 && stop>=start && (strand.equals("+") || strand.equals("-")),
					"FN reference must retain unique identity, family, input, contig and inclusive strand coordinates: "+id);
				Binding binding=byInput.get(inputId);
				if(binding==null){
					binding=new Binding(inputId);
					byInput.put(inputId, binding);
				}
				final Reference ref=new Reference(id, binding, contig, start, stop, strand, member);
				require(byCoordinates.put(ref.key(), ref)==null, "Two FN identities cannot silently claim one original member locus: "+id);
				references.add(ref);
			}
		}
	}

	private void readManifest(){
		assert(references.size()==expected && byInput.size()>0) : "Taxonomy joins are restricted to inputs in the already validated complete FN population";
		try(Input input=new Input(manifest)){
			RrnaResourceIO.header(input.reader, MANIFEST, manifest);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==12, "Formal manifest requires twelve fields");
				final Binding binding=byInput.get(p.parseString(1));
				if(binding==null){continue;}
				require(binding.assembly==null && !p.termEquals("", 0), "Selected input needs exactly one formal assembly binding: "+binding.input);
				binding.assembly=p.parseString(0);
				binding.taxids=p.parseString(5);
				if(priors!=null){
					binding.fasta=p.parseString(3);
					binding.sourceFasta=p.parseString(7);
					binding.accession=p.parseString(6);
				}
			}
		}
		for(final Reference ref:references){require(ref.binding.assembly!=null, "Missing formal manifest input: "+ref.binding.input);}
	}

	private void readSources(){
		assert(byOriginal.size()==0 && sequenceOrder.isEmpty()) : "Source resolution must finish before any sequence acquires an original-gene vote";
		final ObjectSet<String> identities=new ObjectSet<String>(String.class);
		final int[] counts=new int[3];
		try(Input input=new Input(sourceIds)){
			RrnaResourceIO.header(input.reader, SOURCES, sourceIds);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==11, "Original source ledger requires eleven fields");
				final String id=p.parseString(0);
				final String assembly=p.parseString(2);
				final String inputId=p.parseString(3);
				final String contig=p.parseString(4);
				final String strand=p.parseString(7);
				final int member=member(p.parseString(1));
				final int start=p.parseInt(5);
				final int stop=p.parseInt(6);
				final int length=p.parseInt(8);
				require(!id.isEmpty() && identities.add(id) && id.startsWith(MEMBERS[member]+"_input_")
					&& !inputId.isEmpty() && !contig.isEmpty() && start>0 && stop>=start && length==(long)stop-start+1
					&& (strand.equals("+") || strand.equals("-")), "Original IDs/coordinates/member lengths must be unique and intact: "+id);
				if(priors!=null){require(!priorById.contains(id) && (!priorInputs.contains(inputId)
					|| !priorCoordinates.contains(key(inputId, contig, start, stop, strand, member))),
					"Protected prior collides with a frozen native source ID/locus: "+id);}
				counts[member]++;
				if(!byInput.contains(inputId)){continue;}
				final Reference ref=byCoordinates.get(key(inputId, contig, start, stop, strand, member));
				if(ref==null){continue;}
				require(assembly.equals(ref.binding.assembly), "Source/formal assembly disagreement: "+id);
				ref.candidates.add(id);
				if(ref.candidates.size()==1){
					ref.source=id;
					ref.sourceTaxid=p.parseString(9);
				}
			}
		}
		require(Arrays.equals(counts, originals), "Complete source ledger member counts changed");
		for(final Reference ref:references){
			if(ref.candidates.size()==1){require(byOriginal.put(ref.source, ref)==null, "Two primary references cannot create two votes for one original: "+ref.source);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------     Explicit Protected Priors   --------------*/
	/*--------------------------------------------------------------*/

	private void recoverProtectedPriors(){
		assert(priors!=null && priors.size()==expectedPriors) : "Only the declared, FN-bound prior population may append original votes";
		for(final ProtectedPrior prior:priors){
			final Binding binding=prior.ref.binding;
			require(prior.ref.state().equals("NO_ADMITTED_ORIGINAL"), "Protected prior cannot replace an admitted, ambiguous or missing-sequence original: "+prior.id);
			require(binding.assembly.equals(prior.assembly) && canonical(binding.fasta).equals(canonical(prior.fasta))
				&& canonical(binding.sourceFasta).equals(canonical(prior.fasta)), "Protected prior/formal original assembly binding disagrees: "+prior.id);
			prior.canonicalFasta=canonical(prior.fasta);
		}
		readPriorEvidence(priorPanel, 0);
		readPriorEvidence(priorProvenance, 1);
		readPriorEvidence(priorWitnesses, 2);
		readPriorEvidence(priorLoci, 3);
		final ObjectSet<String> historicalIds=new ObjectSet<String>(String.class);
		final int count=readOriginalSequences(priorFasta, read->{
			final String id=RrnaResourceIO.first(read.id);
			require(historicalIds.add(id), "Historical original FASTA repeats an ID: "+id);
			final ProtectedPrior prior=priorById.get(id);
			if(prior==null){return;}
			require(prior.bases==null && read.id.equals(prior.header) && read.length()==prior.length()
				&& sequenceIds.add(id), "Protected prior original header/length/identity changed: "+id);
			prior.bases=Arrays.copyOf(read.bases, read.bases.length);
		});
		require(count==expectedPriorOriginals, "Complete historical original FASTA count changed");
		for(final ProtectedPrior prior:priors){require(prior.bases!=null, "Protected prior historical sequence absent: "+prior.id);}
		final ObjectSet<String> assemblies=new ObjectSet<String>(String.class);
		for(final ProtectedPrior first:priors){
			if(!assemblies.add(first.canonicalFasta)){continue;}
			readOriginalSequences(first.canonicalFasta, read->{
				final String contig=RrnaResourceIO.first(read.id);
				for(final ProtectedPrior prior:priors){
					if(!prior.canonicalFasta.equals(first.canonicalFasta) || !prior.contig.equals(contig)){continue;}
					require(!prior.verified && read.length()==prior.contigLength, "Protected prior target must occur once at its declared contig length: "+prior.id);
					byte[] extracted=Arrays.copyOfRange(read.bases, prior.start-1, prior.stop);
					if(prior.strand.equals("-")){extracted=AminoAcid.reverseComplementBases(extracted);}
					require(Arrays.equals(extracted, prior.bases), "Protected prior original assembly sequence disagrees with historical gene: "+prior.id);
					prior.verified=true;
				}
			});
		}
		for(final ProtectedPrior prior:priors){
			require(prior.verified, "Protected prior literal assembly target absent: "+prior.id);
			final Reference ref=prior.ref;
			ref.candidates.add(prior.id);
			ref.source=prior.id;
			ref.sourceTaxid=prior.taxid;
			ref.header=prior.header;
			ref.bases=prior.bases;
			require(byOriginal.put(prior.id, ref)==null, "Protected prior cannot add a second vote for one original: "+prior.id);
			sequenceOrder.add(ref);
		}
	}

	/** Every authority is checked independently; no similarity or row-order join.
	 * Panel decision_line is the physical one-based line, including the header. */
	private void readPriorEvidence(final String path, final int kind){
		assert(kind>=0 && kind<PRIOR_EVIDENCE_HEADERS.length && priors!=null) : "The declared prior schema determines all field offsets and provenance checks";
		int lineNumber=1;
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, PRIOR_EVIDENCE_HEADERS[kind], path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				lineNumber++;
				p.set(line);
				require(p.terms()==PRIOR_EVIDENCE_FIELDS[kind], "Protected prior evidence schema changed: "+path);
				for(final ProtectedPrior prior:priors){
					final boolean matches=kind==1 ? p.parseInt(0)==prior.priorLine
						: kind==2 ? p.termEquals(prior.reference, 0) && p.termEquals(prior.input, 1) : p.termEquals(prior.id, 0);
					if(!matches){continue;}
					require(!prior.evidence[kind], "Protected prior evidence row is duplicated: "+prior.id+" file="+path);
					final String member=MEMBERS[prior.member];
					boolean valid;
					if(kind==0){
						valid=lineNumber==prior.decisionLine && p.termEquals(prior.genome, 1) && p.termEquals(prior.contig, 2)
							&& p.termEquals(member, 3) && p.termEquals(member, 4) && p.parseInt(5)==prior.start && p.parseInt(6)==prior.stop
							&& p.parseInt(7)==prior.coreStart && p.parseInt(8)==prior.coreStop && p.termEquals(prior.strand, 9);
					}else if(kind==1){
						valid=p.termEquals("PROTECTED_SCORED", 1) && p.termEquals("cobalamin", 2)
							&& canonical(p.parseString(3)).equals(canonical(priorPanel)) && p.parseInt(4)==prior.decisionLine
							&& p.termEquals(prior.id, 5) && p.termEquals(prior.genome, 6) && canonical(p.parseString(7)).equals(prior.canonicalFasta)
							&& p.termEquals(prior.input, 8) && canonical(p.parseString(9)).equals(prior.canonicalFasta)
							&& p.termEquals(prior.contig, 10) && p.parseInt(11)==prior.contigLength && p.termEquals(member, 12) && p.termEquals(member, 13)
							&& p.parseInt(14)==prior.start && p.parseInt(15)==prior.stop && p.parseInt(16)==prior.coreStart
							&& p.parseInt(17)==prior.coreStop && p.termEquals(prior.strand, 18) && p.termEquals("SOURCE_AND_REVIEWED_CORE", 19)
							&& p.termEquals(prior.authority, 20) && p.termEquals(prior.rowAuthority, 21);
					}else if(kind==2){
						valid=p.termEquals("PROTECTED_PRIOR", 2) && p.termEquals(".", 3) && p.termEquals(".", 4) && p.termEquals(".", 5)
							&& p.termEquals(member, 6) && p.parseInt(7)==prior.coreStart && p.parseInt(8)==prior.coreStop
							&& p.parseInt(9)==prior.priorLine && p.termEquals("ORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY", 10);
					}else{
						valid=p.termEquals(prior.genome, 1) && p.termEquals(prior.ref.binding.accession, 2) && p.termEquals(prior.taxid, 3)
							&& p.termEquals(prior.contig, 6) && p.termEquals(prior.strand, 7) && p.termEquals(prior.start+"-"+prior.stop, 8)
							&& p.parseInt(9)==prior.length() && p.termEquals("false", 12) && p.termEquals("false", 13)
							&& p.termEquals("ANNOTATION_FULL_COPY", 14) && p.parseInt(15)==prior.contigLength && p.termEquals("false", 16);
					}
					require(valid, "Protected prior authority/coordinates/source evidence disagrees: "+prior.id+" file="+path);
					prior.evidence[kind]=true;
				}
			}
		}
		if(kind==3){require(lineNumber-1==expectedPriorOriginals, "Complete historical original locus count changed");}
		for(final ProtectedPrior prior:priors){require(prior.evidence[kind], "Protected prior authority evidence absent: "+prior.id+" file="+path);}
	}

	/** Drain the bounded native streamer even after a consumer rejects evidence;
	 * its close operation is not cancellation of a live producer. */
	private static int readOriginalSequences(final String path, final java.util.function.Consumer<Read> consumer){
		assert(path!=null && consumer!=null) : "The explicit prior route binds one immutable sequence file and its validating consumer";
		int count=0;
		Throwable invalid=null;
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(
			FileFormat.testInput(path, FileFormat.FASTA, null, false, true), null, true, -1, false, true, 1))){
			input.stream.start();
			for(ListNum<Read> batch; (batch=input.stream.nextList())!=null;){
				try{
					for(final Read read:batch.list){
						count++;
						if(invalid!=null){continue;}
						try{
							require(read!=null && read.id!=null && !read.id.isEmpty() && read.bases!=null && read.length()>0 && read.mate==null,
								"Protected prior sequence input requires intact unpaired FASTA records");
							consumer.accept(read);
						}catch(RuntimeException | Error failure){invalid=failure;}
					}
				}finally{input.stream.returnList(batch);}
			}
			if(invalid instanceof RuntimeException){throw (RuntimeException)invalid;}
			if(invalid instanceof Error){throw (Error)invalid;}
			require(!input.stream.errorState(), "Protected prior native sequence reader failed: "+path);
		}
		return count;
	}

	/** Shared strict sidecar reader. Reporting validates the row against its audit;
	 * preparation additionally proves all authorities and assembly sequence. */
	static ArrayList<ProtectedPrior> readProtectedPriors(final String path, final int expected){
		require(path!=null && expected>0, "Protected-prior companion requires its exact positive count");
		final ArrayList<ProtectedPrior> rows=new ArrayList<ProtectedPrior>();
		final ObjectSet<String> ids=new ObjectSet<String>(String.class), refs=new ObjectSet<String>(String.class), loci=new ObjectSet<String>(String.class);
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, PRIOR_HEADER, path);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				parser.set(line);
				final ProtectedPrior prior=new ProtectedPrior(parser);
				require(ids.add(prior.id) && refs.add(prior.input+'\t'+prior.reference)
					&& loci.add(key(prior.input, prior.contig, prior.start, prior.stop, prior.strand, prior.member)),
					"Protected prior duplicate reference/original/source locus");
				rows.add(prior);
			}
		}
		require(rows.size()==expected, "Protected-prior companion count changed");
		return rows;
	}

	static final class ProtectedPrior{
		ProtectedPrior(final LineParser1 p){
			require(p.terms()==21, "Protected-prior binding requires twenty-one fields");
			reference=p.parseString(0);
			input=p.parseString(1);
			assembly=p.parseString(2);
			member=member(p.parseString(3));
			contig=p.parseString(4);
			strand=p.parseString(5);
			coreStart=p.parseInt(6);
			coreStop=p.parseInt(7);
			start=p.parseInt(8);
			stop=p.parseInt(9);
			contigLength=p.parseInt(10);
			id=p.parseString(11);
			header=p.parseString(12);
			taxid=p.parseString(13);
			priorLine=p.parseInt(14);
			decisionLine=p.parseInt(15);
			authority=p.parseString(17);
			rowAuthority=p.parseString(18);
			genome=p.parseString(19);
			fasta=p.parseString(20);
			row=p.line().clone();
			require(!reference.isEmpty() && !reference.equals(".") && !input.isEmpty() && !input.equals(".")
				&& !assembly.isEmpty() && !assembly.equals(".") && !contig.isEmpty() && !contig.equals(".")
				&& (strand.equals("+") || strand.equals("-")) && start>0 && coreStart>=start && coreStop>=coreStart
				&& stop>=coreStop && stop<=contigLength && id.startsWith(MEMBERS[member]+"_")
				&& !id.startsWith(MEMBERS[member]+"_input_") && id.equals(RrnaResourceIO.first(header))
				&& priorLine>0 && decisionLine>1 && p.termEquals("PROTECTED_PRIOR", 16)
				&& !authority.isEmpty() && !authority.equals(".") && !rowAuthority.isEmpty() && !genome.isEmpty() && !fasta.isEmpty(),
				"Protected prior requires intact historical identity, containing source/core bounds and explicit authority");
			final int parsed=TaxTree.parseTaxID(header);
			require(taxid.equals(parsed>0 ? Integer.toString(parsed) : "."), "Protected prior source TaxID must preserve the original header");
		}
		int length(){return stop-start+1;}
		boolean matchesAudit(final LineParser1 p){
			assert(p.terms()==14) : "The caller validates the reference audit schema before applying the source-length exception";
			return p.termEquals(reference, 0) && p.termEquals(input, 1) && p.termEquals(assembly, 2) && p.termEquals(contig, 3)
				&& p.parseInt(4)==coreStart && p.parseInt(5)==coreStop && p.termEquals(strand, 6) && p.termEquals(MEMBERS[member], 7)
				&& p.termEquals("RESOLVED", 9) && p.termEquals(id, 10) && p.termEquals(id, 11)
				&& p.termEquals(taxid, 12) && p.termEquals(header, 13);
		}
		final String reference, input, assembly, contig, strand, id, header, taxid, authority, rowAuthority, genome, fasta;
		final int member, coreStart, coreStop, start, stop, contigLength, priorLine, decisionLine;
		final byte[] row;
		// Slots match PRIOR_EVIDENCE_HEADERS: 0=reviewed panel, 1=prior provenance,
		// 2=annotation witnesses, 3=historical original-locus ledger.
		final boolean[] evidence=new boolean[4];
		String canonicalFasta;
		byte[] bases;
		boolean verified;
		Reference ref;
	}

	private static String canonical(final String path){
		try{return new File(path).getCanonicalPath();}
		catch(IOException error){throw new IllegalArgumentException("Cannot resolve original source path: "+path, error);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Original Sequences     ----------------*/
	/*--------------------------------------------------------------*/

	private void readFasta(final int member){
		assert(member>=0 && member<originals.length && fastas.length==originals.length) : "The fixed member order binds each complete FASTA to its own expected denominator; index="+member;
		int count=0;
		Throwable invalid=null;
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(
			FileFormat.testInput(fastas[member], FileFormat.FASTA, null, false, true), null, true, -1, false, true, 1))){
			input.stream.start();
			for(ListNum<Read> batch; (batch=input.stream.nextList())!=null;){
				try{
					for(final Read read:batch.list){
						count++;
						if(invalid!=null){continue;}
						try{retain(read, member);}
						catch(RuntimeException | Error failure){invalid=failure;}
					}
				}finally{input.stream.returnList(batch);}
			}
			// Drain bounded inputs even after a malformed record: native close is not cancellation.
			if(invalid instanceof RuntimeException){throw (RuntimeException)invalid;}
			if(invalid instanceof Error){throw (Error)invalid;}
			require(!input.stream.errorState(), "Native FASTA reader failed: "+fastas[member]);
		}
		require(count==originals[member], "Complete FASTA member count changed: "+MEMBERS[member]+" observed="+count+" expected="+originals[member]);
	}

	private void retain(final Read read, final int member){
		assert(references.size()==expected && byOriginal.size()<=expected) : "Retained sequences may supply only the prevalidated primary references, never extra original votes";
		require(read!=null && read.id!=null && !read.id.isEmpty() && read.bases!=null && read.length()>0 && read.mate==null,
			"Original genes require an intact FASTA header, nonempty sequence and unpaired gene-strand orientation");
		final String id=RrnaResourceIO.first(read.id);
		require(id.startsWith(MEMBERS[member]+"_input_") && sequenceIds.add(id), "Original FASTA identity/member duplicated or changed: "+id);
		final Reference ref=byOriginal.get(id);
		if(ref==null){return;}
		require(ref.member==member && ref.bases==null && read.length()==(long)ref.stop-ref.start+1,
			"Selected original sequence length/member changed: "+id);
		ref.header=read.id;
		ref.bases=Arrays.copyOf(read.bases, read.bases.length);
		sequenceOrder.add(ref);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Complete Audit        ----------------*/
	/*--------------------------------------------------------------*/

	private void writeAudit(final int unresolved){
		assert(references.size()==expected && unresolved>=0 && unresolved<=expected && sequenceOrder.size()+unresolved==expected) : "Every primary reference must be exactly resolved or explicitly unresolved before the complete audit is published";
		final int[] counts=new int[3];
		final int[] resolved=new int[3];
		try(Output audit=new Output(out+".references.tsv")){
			audit.writer.println("reference_id\tinput_id\tassembly_id\tcontig\tstart1\tstop1\tstrand\treference_member\tformal_taxids\tstate\toriginal_candidates\toriginal_source_id\toriginal_taxid\toriginal_header");
			for(final Reference ref:references){
				counts[ref.member]++;
				if(ref.state().equals("RESOLVED")){resolved[ref.member]++;}
				b.clear().append(ref.id).tab().append(ref.binding.input).tab().append(ref.binding.assembly).tab().append(ref.contig)
					.tab().append(ref.start).tab().append(ref.stop).tab().append(ref.strand).tab().append(MEMBERS[ref.member])
					.tab().append(ref.binding.taxids).tab().append(ref.state()).tab();
				for(int i=0; i<ref.candidates.size(); i++){
					if(i>0){b.append(';');}
					b.append(ref.candidates.get(i));
				}
				audit.writer.print(b.tab().append(ref.candidates.size()==1 ? ref.source : ".").tab()
					.append(ref.candidates.size()==1 ? ref.sourceTaxid : ".").tab().append(ref.header==null ? "." : ref.header).nl());
			}
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("member\tprimary_FN_references\tresolved_originals\tunresolved_references\tcohort_state\tinterpretation");
			for(int i=0; i<3; i++){
				summary.writer.print(b.clear().append(MEMBERS[i]).tab().append(counts[i]).tab().append(resolved[i]).tab()
					.append(counts[i]-resolved[i]).tab().append(unresolved==0 ? "COMPLETE" : "BLOCKED_UNRESOLVED")
					.append("\tINCLUSIVE_UNION_FN_NOT_PROOF_OF_NO_PHYSICAL_OVERLAP\n"));
			}
		}
	}

	private static int member(final String value){
		for(int i=0; i<MEMBERS.length; i++){if(MEMBERS[i].equals(value)){return i;}}
		throw new IllegalArgumentException("Unknown raw reference/source member: "+value);
	}
	private static String key(final String input, final String contig, final int start, final int stop, final String strand, final int member){
		return input+'\t'+contig+'\t'+start+'\t'+stop+'\t'+strand+'\t'+MEMBERS[member];
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/** Validated reassignment of a named subset while retaining the full original
	 * observation population. No sequence, profile, vote or caller is rebuilt. */
	private static final class Overlay{
		Overlay(final String[] args){
			Shared.TRIM_READ_DESCRIPTION=false;
			Read.FIX_HEADER=false;
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
				if(a.equals("overlayassignments")){require(Parse.parseBoolean(b), "Use overlayassignments=t for explicit assignment overlay");}
				else if(a.equals("assignments")){assignments=b;}
				else if(a.equals("replacements")){replacements=b;}
				else if(a.equals("audit")){audit=b;}
				else if(a.equals("panel")){panel=b;}
				else if(a.equals("member")){family=MEMBERS[member(b)];}
				else if(a.equals("expectedinput")){expectedInput=Parse.parseIntKMG(b);}
				else if(a.equals("expectedreplace")){expectedReplace=Parse.parseIntKMG(b);}
				else if(a.equals("protectedpriors")){protectedPriors=b;}
				else if(a.equals("expectedpriors")){expectedPriors=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown assignment-overlay option: "+arg);}
			}
			fasta=parser.in1;
			output=parser.out1;
			require(fasta!=null && output!=null && assignments!=null && replacements!=null && audit!=null && panel!=null && family!=null
				&& expectedInput>0 && expectedReplace>0 && expectedReplace<=expectedInput && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1
				&& !Shared.TRIM_READ_DESCRIPTION && !Read.FIX_HEADER,
				"Assignment overlay requires full in=FASTA assignments=original.tsv replacements=new.tsv audit=resolved.tsv panel=panel.sources.tsv member=RFxxxxx expectedinput=N expectedreplace=K out=fresh_prefix without sampling/header repair");
			require(protectedPriors==null ? expectedPriors==0 : expectedPriors>0 && expectedPriors<expectedInput
				&& expectedPriors<=expectedReplace, "Protected overlay requires companion/count and the full augmented population");
			RrnaResourceIO.fresh(output, ".assignments.tsv", ".overlay.tsv", ".clusters.tsv", ".success.txt");
			if(protectedPriors!=null){RrnaResourceIO.fresh(output, ".protected_priors.tsv");}
			Shared.setThreads(1);
			fileIO.ReadWrite.setZipThreads(1);
		}

		void process(){
			readAssignments(assignments, false);
			if(protectedPriors!=null){
				priors=readProtectedPriors(protectedPriors, expectedPriors);
				for(final ProtectedPrior prior:priors){
					require(MEMBERS[prior.member].equals(family) && byId.get(prior.id)==null,
						"Protected overlay source must belong to this member and remain distinct from native originals");
					final OverlayRow row=new OverlayRow(prior.id, prior.header, null);
					row.prior=prior;
					rows.add(row);
					byId.put(row.id, row);
				}
			}
			readFasta();
			readAudit();
			readPanel();
			readAssignments(replacements, true);
			int replaced=0;
			for(final OverlayRow row:rows){
				require(row.prior==null || row.replacement!=null, "Every protected original needs its audited replacement before output: "+row.id);
				require((row.replacement!=null)==(row.reference!=null), "Every audited original needs exactly one replacement: "+row.id);
				if(row.replacement==null){continue;}
				replaced++;
				OverlayGroup group=byGroup.get(row.cluster);
				if(group==null){
					group=new OverlayGroup(row.cluster, row.model, row.status);
					byGroup.put(row.cluster, group);
					groups.add(group);
				}
				require(group.model.equals(row.model) && group.status.equals(row.status), "A replacement cluster cannot change model or disposition");
				row.panelInput=panelInputs.contains(row.input);
				row.panelAssembly=panelAssemblies.contains(row.assembly);
				if(row.panelInput){require(panelInputs.get(row.input).equals(row.assembly), "Panel input/assembly disagrees with the accepted audit");}
				group.count++;
				group.inputs.add(row.input);
				group.assemblies.add(row.assembly);
				if(row.panelInput){group.panelInputs.add(row.input);group.panelInputRows++;}
				if(row.panelAssembly){group.panelAssemblies.add(row.assembly);group.panelAssemblyRows++;}
			}
			require(rows.size()==expectedInput && replaced==expectedReplace, "Full population and replacement counts must be conserved before output");
			write(replaced);
		}

		private void readAssignments(final String path, final boolean replacing){
			int count=0;
			try(Input in=new Input(path)){
				RrnaResourceIO.header(in.reader, ASSIGNMENTS, path);
				for(byte[] line; (line=in.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==6, "Assignment overlay requires the native six-field schema");
					final String header=p.parseString(0), id=RrnaResourceIO.first(header);
					final int taxid=TaxTree.parseTaxID(header);
					require((id.startsWith(family+"_input_") || replacing && byId.get(id)!=null && byId.get(id).prior!=null)
						&& p.parseInt(1)==(taxid>0 ? taxid : -1),
						"Assignment member/sourceTid must preserve its original header: "+id);
					final String cluster=p.parseString(3), model=p.parseString(4), status=p.parseString(5);
					require(!p.termEquals("", 2) && cluster.startsWith(family+"_")
						&& (status.equals("kept") ? RrnaResourceIO.first(model).startsWith("cobalamin_"+family+"_")
							: (status.equals("small") || status.equals("orphan")) && model.isEmpty()),
						"Assignment requires a member-bound cluster/model and known disposition");
					final OverlayRow old=byId.get(id);
					if(replacing){
						require(old!=null && old.reference!=null && old.replacement==null && old.header.equals(header),
							"Replacement must uniquely match an audited original ID/full header: "+id);
						old.replacement=line.clone();
						old.cluster=cluster;
						old.model=model;
						old.status=status;
					}else{
						require(old==null, "Original assignment identities must be unique: "+id);
						final OverlayRow row=new OverlayRow(id, header, line.clone());
						byId.put(id, row);
						rows.add(row);
					}
					count++;
				}
			}
			require(count==(replacing ? expectedReplace : expectedInput-expectedPriors), "Assignment row count differs from declared "+(replacing ? "replacement" : "original")+" population");
		}

		private void readFasta(){
			int count=0;
			Throwable invalid=null;
			try(SequenceInput in=new SequenceInput(StreamerFactory.makeStreamer(
				FileFormat.testInput(fasta, FileFormat.FASTA, null, false, true), null, true, -1, false, true, 1))){
				in.stream.start();
				for(ListNum<Read> batch; (batch=in.stream.nextList())!=null;){
					try{
						for(final Read read:batch.list){
							count++;
							if(invalid!=null){continue;}
							try{
								require(read!=null && read.id!=null && read.bases!=null && read.length()>0 && read.mate==null, "Overlay original FASTA must contain intact unpaired genes");
								final OverlayRow row=byId.get(RrnaResourceIO.first(read.id));
								require(row!=null && row.length==0 && row.header.equals(read.id), "Full FASTA must uniquely match every original assignment ID/header");
								require(row.prior==null || read.length()==row.prior.length(), "Protected overlay FASTA must retain its full source length");
								row.length=read.length();
							}catch(RuntimeException | Error failure){invalid=failure;}
						}
					}finally{in.stream.returnList(batch);}
				}
				if(invalid instanceof RuntimeException){throw (RuntimeException)invalid;}
				if(invalid instanceof Error){throw (Error)invalid;}
				require(!in.stream.errorState(), "Full original FASTA reader failed");
			}
			require(count==expectedInput, "Full original FASTA count changed");
			for(final OverlayRow row:rows){require(row.length>0, "Original assignment missing from full FASTA: "+row.id);}
		}

		private void readAudit(){
			int count=0;
			final ObjectSet<String> references=new ObjectSet<String>(String.class);
			try(Input in=new Input(audit)){
				RrnaResourceIO.header(in.reader, AUDIT, audit);
				for(byte[] line; (line=in.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==14, "Accepted original audit requires fourteen fields");
					if(!p.termEquals(family, 7)){continue;}
					final String id=p.parseString(11), header=p.parseString(13), sourceTaxid=p.parseString(12);
					final OverlayRow row=byId.get(id);
					// Saved grade partitions may reuse reference IDs in distinct inputs.
					require(row!=null && row.reference==null && p.termEquals("RESOLVED", 9) && p.termEquals(id, 10)
						&& row.header.equals(header) && !p.termEquals("", 0) && !p.termEquals(".", 0)
						&& references.add(p.parseString(1)+'\t'+p.parseString(0)),
						"Audit must uniquely resolve each selected original ID/full header");
					final int tid=TaxTree.parseTaxID(header), auditTid=sourceTaxid.equals(".") ? -1 : Parse.parseInt(sourceTaxid, 0);
					require(auditTid==(tid>0 ? tid : -1), "Audit original taxid differs from the original header; formal taxonomy is separate");
					row.reference=p.parseString(0);
					row.input=p.parseString(1);
					row.assembly=p.parseString(2);
					row.formalTaxids=p.parseString(8);
					require(!row.input.isEmpty() && !row.assembly.isEmpty() && !row.assembly.equals(".")
						&& (row.prior!=null ? row.prior.matchesAudit(p) && row.length==row.prior.length()
							: id.startsWith(family+"_input_"+row.input+"_locus_") && p.parseInt(4)>0
							&& (long)p.parseInt(5)-p.parseInt(4)+1==row.length && (p.termEquals("+", 6) || p.termEquals("-", 6))),
						"Audit input/assembly and original gene coordinates must remain intact");
					count++;
				}
			}
			require(count==expectedReplace, "Accepted audit selected-member count differs from declared replacement population");
		}

		private void readPanel(){
			int count=0;
			try(Input in=new Input(panel)){
				RrnaResourceIO.header(in.reader, "input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff", panel);
				for(byte[] line; (line=in.reader.nextLine())!=null;){
					p.set(line);
					require(p.terms()==5 && p.parseInt(0)==count && !p.termEquals("", 1) && !p.termEquals("", 2), "Panel requires ordered unique input/assembly identities");
					require(panelInputs.put(p.parseString(2), p.parseString(1))==null && panelAssemblies.add(p.parseString(1)), "Panel input/assembly identities must be unique");
					count++;
				}
			}
			require(count>0, "Panel intersection requires a nonempty saved roster");
		}

		private void write(final int replaced){
			final ByteBuilder b=new ByteBuilder();
			try(Output out=new Output(output+".assignments.tsv"); Output auditOut=new Output(output+".overlay.tsv")){
				out.writer.println(ASSIGNMENTS);
				auditOut.writer.println("ordinal0\toriginal_id\taction\treference_id\tinput_id\tassembly_id\tformal_taxids\tpanel_input_intersection\tpanel_assembly_intersection");
				for(int i=0; i<rows.size(); i++){
					final OverlayRow row=rows.get(i);
					out.writer.print(b.clear().append(row.replacement==null ? row.original : row.replacement).nl());
					b.clear().append(i).tab().append(row.id).tab();
					if(row.replacement==null){b.append("UNCHANGED\t.\t.\t.\t.\t.\t.");}
					else{
						b.append(row.prior==null ? "REPLACED" : "ADDED_PROTECTED_PRIOR").tab().append(row.reference).tab().append(row.input).tab().append(row.assembly)
							.tab().append(row.formalTaxids).tab().append(row.panelInput).tab().append(row.panelAssembly);
					}
					auditOut.writer.print(b.nl());
				}
			}
			try(Output out=new Output(output+".clusters.tsv")){
				out.writer.println("cluster_key\tmodel_label\tstatus\toriginals\tinputs\tassemblies\tpanel_input_originals\tpanel_assembly_originals\tpanel_inputs\tpanel_assemblies\tinterpretation");
				for(final OverlayGroup group:groups){
					out.writer.print(b.clear().append(group.cluster).tab().append(group.model).tab().append(group.status).tab().append(group.count)
						.tab().append(group.inputs.size()).tab().append(group.assemblies.size()).tab().append(group.panelInputRows).tab().append(group.panelAssemblyRows)
						.tab().append(group.panelInputs.size()).tab().append(group.panelAssemblies.size())
						.append("\tEXACT_INPUT_ASSEMBLY_INTERSECTIONS_NOT_SEQUENCE_ALIAS_INDEPENDENCE\n"));
				}
			}
			// Independent readback against original source bytes before a success marker.
			try(Input in=new Input(output+".assignments.tsv")){
				RrnaResourceIO.header(in.reader, ASSIGNMENTS, output);
				for(final OverlayRow row:rows){require(Arrays.equals(in.reader.nextLine(), row.replacement==null ? row.original : row.replacement), "Published overlay changed a source row or original order");}
				require(in.reader.nextLine()==null, "Published overlay added an original row");
			}
			if(priors!=null){
				try(Output out=new Output(output+".protected_priors.tsv")){
					out.writer.println(PRIOR_HEADER);
					for(final ProtectedPrior prior:priors){out.writer.println(prior.row);}
				}
			}
			final String success="NCRNA_ASSIGNMENT_OVERLAY_PASS originals="+rows.size()+" replaced="+replaced+" unchanged="+(rows.size()-replaced)
				+(priors==null ? "" : " protectedpriors="+expectedPriors);
			try(Output out=new Output(output+".success.txt")){out.writer.println(success);}
			System.err.println(success);
		}

		private String fasta, assignments, replacements, audit, panel, family, output, protectedPriors;
		private int expectedInput=-1, expectedReplace=-1, expectedPriors;
		private ArrayList<ProtectedPrior> priors;
		private final LineParser1 p=new LineParser1('\t');
		private final ArrayList<OverlayRow> rows=new ArrayList<OverlayRow>();
		private final ObjectMap<String, OverlayRow> byId=new ObjectMap<String, OverlayRow>(String.class, OverlayRow.class);
		private final ArrayList<OverlayGroup> groups=new ArrayList<OverlayGroup>();
		private final ObjectMap<String, OverlayGroup> byGroup=new ObjectMap<String, OverlayGroup>(String.class, OverlayGroup.class);
		private final ObjectMap<String, String> panelInputs=new ObjectMap<String, String>(String.class, String.class);
		private final ObjectSet<String> panelAssemblies=new ObjectSet<String>(String.class);
	}
	private static final class OverlayRow{
		OverlayRow(final String id_, final String header_, final byte[] original_){
			id=id_;
			header=header_;
			original=original_;
		}
		final String id, header;
		final byte[] original;
		ProtectedPrior prior;
		byte[] replacement;
		String reference, input, assembly, formalTaxids, cluster, model, status;
		int length;
		boolean panelInput, panelAssembly;
	}
	private static final class OverlayGroup{
		OverlayGroup(final String cluster_, final String model_, final String status_){
			cluster=cluster_;
			model=model_;
			status=status_;
		}
		final String cluster, model, status;
		int count, panelInputRows, panelAssemblyRows;
		final ObjectSet<String> inputs=new ObjectSet<String>(String.class), assemblies=new ObjectSet<String>(String.class),
			panelInputs=new ObjectSet<String>(String.class), panelAssemblies=new ObjectSet<String>(String.class);
	}

	private static final class Binding{
		Binding(final String id){input=id;}
		final String input;
		String assembly, taxids, fasta, sourceFasta, accession;
	}
	private static final class Reference{
		Reference(final String id_, final Binding binding_, final String contig_, final int start_, final int stop_, final String strand_, final int member_){
			id=id_;
			binding=binding_;
			contig=contig_;
			start=start_;
			stop=stop_;
			strand=strand_;
			member=member_;
		}
		String key(){return PrepareNcrnaMissingModels.key(binding.input, contig, start, stop, strand, member);}
		String state(){
			assert((bases==null)==(header==null) && (bases==null || candidates.size()==1)) : "Resolved sequence/header evidence belongs to exactly one admitted original, never an ambiguous locus";
			return candidates.isEmpty() ? "NO_ADMITTED_ORIGINAL" : candidates.size()>1 ? "AMBIGUOUS_ORIGINALS"
				: bases==null ? "ORIGINAL_SEQUENCE_ABSENT" : "RESOLVED";
		}
		final String id, contig, strand;
		final Binding binding;
		final int start, stop, member;
		final ArrayList<String> candidates=new ArrayList<String>();
		String source, sourceTaxid, header;
		byte[] bases;
	}

	private String grade, sourceIds, manifest, out;
	private String priorBindings, priorFasta, priorLoci, priorProvenance, priorWitnesses, priorPanel;
	private int expectedPriors=-1, expectedPriorOriginals=-1;
	private ArrayList<ProtectedPrior> priors;
	private final ObjectSet<String> priorInputs=new ObjectSet<String>(String.class);
	private final ObjectMap<String,ProtectedPrior> priorById=new ObjectMap<String,ProtectedPrior>(String.class, ProtectedPrior.class),
		priorCoordinates=new ObjectMap<String,ProtectedPrior>(String.class, ProtectedPrior.class);
	private String[] fastas, gradePrefixes;
	private int expected=-1;
	private int[] originals;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final Target rosterProbe=new Target();
	private final ArrayList<Reference> references=new ArrayList<Reference>(), sequenceOrder=new ArrayList<Reference>();
	private final ObjectMap<String,Reference> byCoordinates=new ObjectMap<String,Reference>(String.class, Reference.class),
		byOriginal=new ObjectMap<String,Reference>(String.class, Reference.class);
	private final ObjectMap<String,Binding> byInput=new ObjectMap<String,Binding>(String.class, Binding.class);
	private final ObjectSet<String> sequenceIds=new ObjectSet<String>(String.class);
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String GRADE="input_id\trole\tfamily\tseqid\tstart\tstop\tstrand\toutcome\treference_id\tquery_line\tmatched_reference_id\tmatched_query_line\tmask_reference_id\tduplicate\tscored_eligible\tmodel\td5\td3";
	private static final String SOURCES="source_id\tmember\tassembly_id\tinput_id\tcontig\tstart1\tstop1\tstrand\tlength\ttaxid\tchunk";
	private static final String MANIFEST="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
	private static final String ASSIGNMENTS="#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus";
	private static final String AUDIT="reference_id\tinput_id\tassembly_id\tcontig\tstart1\tstop1\tstrand\treference_member\tformal_taxids\tstate\toriginal_candidates\toriginal_source_id\toriginal_taxid\toriginal_header";
	private static final String PRIOR_HEADER="reference_id\tinput_id\tassembly_id\treference_member\tcontig\tstrand\tcore_start1\tcore_stop1\tsource_start1\tsource_stop1\tcontig_length\toriginal_id\toriginal_header\toriginal_taxid\tprior_line\tdecision_line\tsource_kind\tsource_authority\trow_authority\toriginal_genome\tsource_fasta";
	private static final int[] PRIOR_EVIDENCE_FIELDS={10, 23, 11, 21};
	private static final String[] PRIOR_EVIDENCE_HEADERS={
		"locus_id\tgenome\tseqid\tsource_family\tmember_family\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand",
		"prior_line\tretained_disposition\tsource_id\tdecision_file\tdecision_line\toriginal_locus\toriginal_genome\toriginal_fasta\tinput_id\tresolved_fasta\tliteral_target\tcontig_length\tsource_family\tassignment\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand\tinterval_semantics\tsource_authority\trow_authority\toriginal_row",
		"annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation",
		"locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines"
	};
}
