package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;
import java.util.function.Consumer;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntObjectMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Random;
import shared.Shared;
import shared.Timer;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Projects accepted raw native-CM loci for native CutGff, then independently
 * verifies its output against those loci. Every admitted original keeps one ID
 * and one vote, including equal DNA in different records. This class does not
 * extract sequences, choose clan assignments or interpret reference-GFF tiers.
 * Metadata lifecycle follows template.A_SampleByteFile; sequence validation uses
 * template.A_SampleStreamerMT's native Streamer/Writer batch interface with one
 * consumer, because validation has no alignment or other expensive per-read work.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminConstructionInputs{

	public static void main(final String[] args) throws IOException{
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, CobalaminConstructionInputs.class, false);
		try{
			final CobalaminConstructionInputs tool=new CobalaminConstructionInputs(pp.args, pp.outstream);
			tool.process();
			timer.stop();
			pp.outstream.println("COBALAMIN_CONSTRUCTION_INPUTS_COMPLETE time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private CobalaminConstructionInputs(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){
				mode=b;
			}else if(a.equals("loci")){
				lociPath=b;
			}else if(a.equals("sequences")){
				sequences=b;
			}else if(a.equals("flank")){
				flank=Parse.parseIntKMG(b);
			}else if(a.equals("member")){
				require(b!=null && !b.isEmpty(), "member= needs one exact producing-CM member");
				selectedName=b;
			}else if(a.equals("cores")){
				cores=b.split(",", -1);
			}else if(a.equals("expectedgenomes")){
				expectedGenomes=Parse.parseIntKMG(b);
			}else if(a.equals("expected")){
				final LineParser1 values=new LineParser1(',');
				values.set(b.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
				require(values.terms()==3, "expected= binds the three raw member counts");
				for(int i=0; i<3; i++){expected[i]=values.parseInt(i);}
			}else if(a.equals("chunkgenomes")){
				chunkGenomes=Parse.parseIntKMG(b);
			}else if(a.equals("chunk")){
				chunk=Parse.parseIntKMG(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else{require(parser.parse(arg, a, b), "Unknown construction-input option: "+arg);}
		}
		manifest=parser.in1;
		directory=parser.out1;
		require(("prepare".equals(mode) || "validate".equals(mode)) && manifest!=null && lociPath!=null
			&& directory!=null && expectedGenomes>0 && chunkGenomes>0 && chunk>=0 && seed>=0 && flank>=0
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1,
			"Use mode=prepare|validate in=selected_MANIFEST.tsv loci=accepted_native_loci.tsv.gz out=fresh_directory expectedgenomes=N expected=A,B,C chunkgenomes=N seed=N; validate also needs sequences= and chunk=");
		for(final int value:expected){require(value>=0, "Each member needs an explicit expected count");}
		if(selectedName!=null){
			for(int member=0; member<MEMBERS.length; member++){
				if(MEMBERS[member].equals(selectedName)){selectedMember=member;}
			}
			require(selectedMember>=0, "member= must exactly name RF00174, RF01482 or RF01689");
		}
		require(mode.equals("prepare") ? sequences==null : sequences!=null, "Only validate consumes extracted sequences");
		require(cores==null ? flank==0 || mode.equals("prepare") : flank>0 && mode.equals("validate") && cores.length==(selectedMember<0 ? 3 : 1),
			"Flanked validation requires all three cores= files, or exactly the selected member's one core file; zero-flank construction is unchanged");
		require(!new File(directory).exists(), "Outputs must be fresh: "+directory);
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_PIGZ=false;
		ReadWrite.setZipThreads(1);
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
	}

	private void process() throws IOException{
		loadManifest();
		loadLoci();
		assignChunks();
		require(new File(directory).mkdir(), "Could not create output directory: "+directory);
		if(mode.equals("prepare")){prepare();}
		else{
			if(flank>0){loadCores();}
			validate();
		}
	}

	/** The selected manifest alone creates inputs; taxIDs do not select genomes. */
	private void loadManifest() throws IOException{
		final ObjectMap<String, Assembly> ids=new ObjectMap<String, Assembly>(String.class, Assembly.class);
		try(Input in=new Input(manifest)){
			RrnaResourceIO.header(in.reader, MANIFEST_HEADER, manifest);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==12, "Selected assembly manifest requires twelve columns");
				final Assembly assembly=new Assembly(p.parseString(0), p.parseInt(1), p.parseString(3), p.parseString(4), p.parseString(7), taxid(5));
				require(assembly.input>=0 && !byInput.contains(assembly.input) && !ids.contains(assembly.id),
					"Duplicate selected assembly/input: "+assembly.id+" input="+assembly.input);
				validatePath(assembly.fasta);
				require(samePath(assembly.fasta, assembly.source), "Selected link must resolve to its scanned FASTA: "+assembly.id);
				byInput.put(assembly.input, assembly);
				ids.put(assembly.id, assembly);
				assemblies.add(assembly);
			}
		}
		require(assemblies.size()==expectedGenomes, "Selected manifest count="+assemblies.size()+" expected="+expectedGenomes);
	}

	/** Consumes only the accepted raw-locus ledger. Truncation is the sole admission
	 * exclusion; equal DNA is never examined here. One native event per coordinate
	 * locus was measured in this corpus and is required before assigning source IDs. */
	private void loadLoci() throws IOException{
		final int[] counts=new int[3];
		long ordinal=0;
		try(Input in=new Input(lociPath)){
			RrnaResourceIO.header(in.reader, LOCUS_HEADER, lociPath);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				require(p.terms()==13, "Accepted raw locus requires thirteen columns");
				final Assembly assembly=byInput.get(p.parseInt(1));
				require(assembly!=null && p.termEquals(assembly.id, 0), "Locus outside exact selected assembly set");
				require(samePath(assembly.source, p.parseString(3)),
					"Locus source must resolve to the selected scanned assembly: "+assembly.id);
				final int member=member(4), start=p.parseInt(6), stop=p.parseInt(7), events=p.parseInt(9), full=p.parseInt(10);
				final char strand=p.termEquals("+", 8) ? '+' : p.termEquals("-", 8) ? '-' : '.';
				require(member>=0 && start>0 && stop>=start && strand!='.' && events==1 && full>=0 && full<=events
					&& (p.termEquals("true", 11) ? full>0 : p.termEquals("false", 11) && full==0),
					"Original coordinate/event/admission invariant failed: assembly="+assembly.id+" raw_locus_ordinal0="+ordinal);
				final String contig=p.parseString(5);
				require(contig.length()>0 && contig.indexOf(' ')<0, "Raw contig must be an intact GFF seqid: "+contig);
				final String key=member+"\t"+contig+"\t"+start+"\t"+stop+"\t"+strand;
				require(!assembly.coordinates.contains(key), "Duplicate source locus: "+assembly.id+" "+key);
				assembly.coordinates.put(key, key);
				final String id=MEMBERS[member]+"_input_"+assembly.input+"_locus_"+ordinal;
				ordinal++;
				if(full==0){continue;}
				final Locus locus=new Locus(id, assembly, member, contig, start, stop, strand, flank);
				require(!byId.contains(id), "Duplicate admitted source ID: "+id);
				byId.put(id, locus);
				assembly.loci.add(locus);
				loci.add(locus);
				counts[member]++;
			}
		}
		for(int i=0; i<3; i++){
			require(counts[i]==expected[i], "Raw admitted "+MEMBERS[i]+" records="+counts[i]+" expected="+expected[i]);
		}
	}

	/** Fixed-seed Fisher-Yates shuffle over ALL selected genomes, including zero-locus
	 * inputs. Chunks are disjoint; the first chunk is retained after its cost assay. */
	private void assignChunks(){
		assert(assemblies.size()==expectedGenomes) : "Uniform chunks must include the complete selected roster";
		final Random random=Shared.random(seed);
		for(int i=assemblies.size()-1; i>0; i--){
			final int j=random.nextInt(i+1);
			final Assembly value=assemblies.get(i);
			assemblies.set(i, assemblies.get(j));
			assemblies.set(j, value);
		}
		for(int i=0; i<assemblies.size(); i++){assemblies.get(i).chunk=i/chunkGenomes;}
		require(chunk<chunkCount(), "Requested chunk outside selected roster: "+chunk);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Projection Output       ----------------*/
	/*--------------------------------------------------------------*/

	private void prepare(){
		final File projections=new File(directory, "projections"), chunks=new File(directory, "chunks");
		require(projections.mkdir() && chunks.mkdir(), "Could not create fresh extraction projection/chunk directories");
		try(Output ledger=output("source_ids.tsv"); Output inputs=output("inputs.tsv")){
			ledger.writer.println("source_id\tmember\tassembly_id\tinput_id\tcontig\tstart1\tstop1\tstrand\tlength\ttaxid\tchunk");
			inputs.writer.println("assembly_id\tinput_id\tchunk\tadmitted_loci\tlinked_fasta\tsource_fasta");
			for(final Assembly assembly:assemblies){
				final int selected=selectedCount(assembly);
				inputs.writer.print(buffer.clear().append(assembly.id).tab().append(assembly.input).tab().append(assembly.chunk)
					.tab().append(selected).tab().append(assembly.fasta).tab().append(assembly.source).nl());
				if(selected==0){continue;}
				try(Output gff=new Output(projection(assembly))){
					gff.writer.println("##gff-version 3");
					for(final Locus locus:assembly.loci){
						if(!wanted(locus.member)){continue;}
						gff.writer.print(buffer.clear().append(locus.contig).append("\traw_native_CM\tncRNA\t")
							.append(locus.start).tab().append(locus.stop).append("\t.\t").append(locus.strand)
							.append("\t.\t").append(locus.header()).nl());
						ledger.writer.print(buffer.clear().append(locus.id).tab().append(MEMBERS[locus.member]).tab().append(assembly.id)
							.tab().append(assembly.input).tab().append(locus.contig).tab().append(locus.start).tab().append(locus.stop)
							.tab().append(locus.strand).tab().append(locus.length()).tab().append(assembly.taxid).tab().append(assembly.chunk).nl());
					}
				}
			}
		}
		try(Output table=output("chunks.tsv")){
			table.writer.println("chunk\tselected_genomes\tgenomes_with_loci\tadmitted_loci\tconfig\textracted_fasta");
			for(int index=0; index<chunkCount(); index++){
				int genomes=0, positive=0, records=0;
				final String prefix=new File(chunks, Integer.toString(index)).getPath();
				try(Output config=new Output(prefix+".config")){
					for(final Assembly assembly:assemblies){
						if(assembly.chunk!=index){continue;}
						genomes++;
						final int selected=selectedCount(assembly);
						if(selected==0){continue;}
						positive++;
						records+=selected;
						config.writer.print(buffer.clear().append("in=").append(assembly.fasta).nl()
							.append("gff=").append(projection(assembly)).nl());
					}
					config.writer.print(buffer.clear().append("out=").append(prefix).append(".fa\n"));
				}
				table.writer.print(buffer.clear().append(index).tab().append(genomes).tab().append(positive)
					.tab().append(records).tab().append(prefix).append(".config\t").append(prefix).append(".fa\n"));
			}
		}
		writeComplete("PROJECTIONS_PREPARED_NOT_EXTRACTED", selectedMember<0 ? loci.size() : expected[selectedMember]);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Independent Validation   ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads only the admitted original cores, not genomes or new biological labels.
	 * Each locus gets exactly one original byte array for the subsequent context join.
	 */
	private void loadCores(){
		final int[] observed=new int[3];
		for(int member=0; member<3; member++){
			if(!wanted(member)){continue;}
			final int expectedMember=member;
			streamBatches(cores[selectedMember<0 ? member : 0], batch->{
				for(final Read read:batch.list){
					final Locus locus=byId.get(RrnaResourceIO.first(read.id));
					require(locus!=null && locus.member==expectedMember && locus.core==null
						&& locus.coreHeader().equals(read.id) && read.bases!=null && read.length()==locus.length(),
						"Admitted original core ID/member/header/length differs: "+read.id);
					locus.core=read.bases;
					observed[expectedMember]++;
				}
			});
		}
		for(int member=0; member<3; member++){
			if(!wanted(member)){continue;}
			require(observed[member]==expected[member], "Original core count differs for "+MEMBERS[member]);
		}
	}

	/** Native FASTA readers/writers retain record boundaries. Errors accumulate while
	 * the reader drains, so a consumer error cannot strand the reader's producer.
	 * Outputs are accepted only by the final completion file after exact-set checks. */
	private void validate(){
		final Writer[] writers=new Writer[3];
		try{
			for(int i=0; i<3; i++){
				if(!wanted(i)){continue;}
				final FileFormat format=FileFormat.testOutput(new File(directory, MEMBERS[i]+".fa").getPath(), FileFormat.FASTA, null, false, false, false, false);
				writers[i]=WriterFactory.makeWriter(format, null, 1, null, false);
				writers[i].start();
			}
			if(new File(sequences).isFile()){streamBatches(sequences, batch->{
				@SuppressWarnings("unchecked")
				final ArrayList<Read>[] groups=new ArrayList[3];
				for(int i=0; i<3; i++){groups[i]=new ArrayList<Read>();}
				for(final Read read:batch.list){
					final String id=RrnaResourceIO.first(read.id);
					final Locus locus=byId.get(id);
					if(locus==null || !wanted(locus.member) || locus.assembly.chunk!=chunk){
						unexpected.add(id);
						continue;
					}
					locus.seen++;
					locus.observedLength=read.length();
					if(flank==0){locus.headerMatches&=read.id.equals(locus.header());}
					else{validateContext(read, locus);}
					if(!locus.valid()){continue;}
					groups[locus.member].add(read);
				}
				for(int i=0; i<3; i++){if(!groups[i].isEmpty()){writers[i].add(groups[i], batch.id);}}
			});}
		}finally{
			closeWriters(writers);
		}
		int failures=unexpected.size(), records=0;
		for(final Assembly assembly:assemblies){
			if(assembly.chunk!=chunk){continue;}
			boolean bad=false;
			for(final Locus locus:assembly.loci){
				if(!wanted(locus.member)){continue;}
				if(!locus.valid()){
					bad=true;
					failures++;
				}
				records++;
			}
			if(bad){diagnose(assembly);}
		}
		try(Output differences=output("id_differences.tsv")){
			differences.writer.println("state\tsource_id\tassembly_id\tinput_id\tsource_fasta\tcontig\tstart1\tstop1\texpected_length\tobserved_length\toccurrences\tactual_contig_length");
			for(final String id:unexpected){differences.writer.print(buffer.clear().append("UNEXPECTED\t").append(id).append("\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"));}
			for(final Locus locus:loci){
				if(!wanted(locus.member) || locus.assembly.chunk!=chunk || locus.valid()){continue;}
				differences.writer.print(buffer.clear().append(locus.seen==0 ? "MISSING" : locus.seen>1 ? "DUPLICATE"
					: !locus.headerMatches ? "HEADER_MISMATCH" : locus.observedLength!=locus.length() ? "LENGTH_MISMATCH" : "SEQUENCE_MISMATCH")
					.tab().append(locus.id).tab().append(locus.assembly.id).tab().append(locus.assembly.input)
					.tab().append(locus.assembly.fasta).tab().append(locus.contig)
					.tab().append(locus.start).tab().append(locus.stop).tab().append(locus.length()).tab().append(locus.observedLength)
					.tab().append(locus.seen).tab().append(locus.contigLength<0 ? "ABSENT" : Integer.toString(locus.contigLength)).nl());
				outstream.println("Rejected extraction: source_id="+locus.id+" input="+locus.assembly.input+" source_fasta="+locus.assembly.fasta+" contig="+locus.contig
					+" start1="+locus.start+" stop1="+locus.stop+" expected_length="+locus.length()+" observed_length="+locus.observedLength
					+" occurrences="+locus.seen+" actual_contig_length="+(locus.contigLength<0 ? "ABSENT" : locus.contigLength));
			}
		}
		require(failures==0, "Extraction source-ID set/length validation failed: "+failures+" differences; see "+new File(directory, "id_differences.tsv"));
		if(flank>0){
			writeLearning();
			if(selectedMember>=0){writeClipping();}
		}
		writeComplete(flank>0 ? "NATIVE_CONTEXT_IDS_FLANKS_GC_AND_CORE_BYTES_PASS" : "EXACT_RAW_MEMBER_ID_SETS_AND_LENGTHS_PASS", records);
	}

	/** CutGff appends actual sense flanks and caller-GC in this exact order. The
	 * genomic-left flank is independently fixed by start1; the other is edge-clipped.
	 * Core bytes must agree with the separately accepted exact extraction.
	 */
	private void validateContext(final Read read, final Locus locus){
		assert(locus.core!=null) : "Flanked validation loads the admitted original core before accepting context bytes";
		final String header=read.id, prefix=locus.header()+" lflank=";
		final int right=header.indexOf(" rflank=", prefix.length());
		final int gc=right<0 ? -1 : header.indexOf(" contig_gc=", right+8);
		boolean valid=header.startsWith(prefix) && right>=prefix.length() && gc>right+8;
		if(valid){
			try{
				locus.left=NcrnaNumbers.unsignedInt(header, prefix.length(), right);
				locus.right=NcrnaNumbers.unsignedInt(header, right+8, gc);
				final float fraction=Float.parseFloat(header.substring(gc+11));
				valid=locus.left<=flank && locus.right<=flank && fraction>=0 && fraction<=1
					&& (locus.strand=='-' ? locus.right : locus.left)==Math.min(flank, locus.start-1);
			}catch(IllegalArgumentException failure){valid=false;}
		}
		locus.headerMatches&=valid;
		if(!valid){return;}
		locus.observedLength=read.length()-locus.left-locus.right;
		if(read.bases==null || locus.observedLength!=locus.length()){
			locus.sequenceMatches=false;
			return;
		}
		for(int i=0; i<locus.core.length; i++){
			if(locus.core[i]!=read.bases[locus.left+i]){
				locus.sequenceMatches=false;
				return;
			}
		}
	}

	/** Native-context metadata preserves original graph IDs and observed flanks.
	 * No contig census or phylum lookup was performed: those optional roster fields
	 * remain explicitly unknown, while source identity, TaxID and counts are exact.
	 */
	private void writeLearning(){
		assert(flank>0 && cores!=null) : "Only independently checked native contexts may create learning metadata";
		for(int member=0; member<3; member++){
			if(!wanted(member)){continue;}
			try(Output roster=output(MEMBERS[member]+".genomes.tsv"); Output ledger=output(MEMBERS[member]+".loci.tsv")){
				roster.writer.println(NcrnaLearningPopulation.ROSTER_HEADER);
				ledger.writer.println(NcrnaLearningPopulation.NATIVE_HEADER);
				int index=0;
				for(final Assembly assembly:assemblies){
					if(assembly.chunk!=chunk){continue;}
					int count=0;
					for(final Locus locus:assembly.loci){if(locus.member==member){count++;}}
					roster.writer.print(buffer.clear().append(index++).tab().append(assembly.id).tab().append(assembly.id)
						.tab().append(assembly.taxid>0 ? Integer.toString(assembly.taxid) : ".")
						.append("\tUNKNOWN\tUNKNOWN\t").append(assembly.gff).tab().append(assembly.fasta).tab().append(count).append("\t.\t.\n"));
					for(final Locus locus:assembly.loci){
						if(locus.member!=member){continue;}
						ledger.writer.print(buffer.clear().append(locus.id).tab().append(locus.id).tab().append(assembly.id).tab().append(assembly.id)
							.tab().append(assembly.taxid>0 ? Integer.toString(assembly.taxid) : ".")
							.append("\tUNKNOWN\tUNKNOWN\t").append(locus.contig).tab().append(locus.strand)
							.tab().append(locus.start).tab().append(locus.stop).tab().append(locus.length())
							.tab().append(locus.left).tab().append(locus.right).nl());
					}
				}
			}
		}
	}

	/** Observed physical context clipping, separate from later jittered-span failures. */
	private void writeClipping(){
		int records=0, left=0, right=0, both=0, full=0;
		for(final Locus locus:loci){
			if(!wanted(locus.member) || locus.assembly.chunk!=chunk){continue;}
			assert(locus.valid()) : "Clipping accounting follows exact source/core validation";
			records++;
			final boolean leftClipped=locus.left<flank, rightClipped=locus.right<flank;
			if(leftClipped){left++;}
			if(rightClipped){right++;}
			if(leftClipped && rightClipped){both++;}
			if(!leftClipped && !rightClipped){full++;}
		}
		assert(records==left+right-both+full) : "Every selected original has full flanks or observed clipping at one or both sense ends";
		try(Output out=output("context_geometry.tsv")){
			out.writer.println("member\trequested_flank_bp\toriginals\tsense_left_clipped\tsense_right_clipped\tboth_clipped\tfull_flanks");
			out.writer.print(buffer.clear().append(MEMBERS[selectedMember]).tab().append(flank).tab().append(records)
				.tab().append(left).tab().append(right).tab().append(both).tab().append(full).nl());
		}
	}

	/** Failure-only DNA read of affected assemblies. This observes contig lengths;
	 * it does not guess why a source record was absent or malformed. */
	private void diagnose(final Assembly assembly){
		assert(assembly.chunk==chunk) : "Only the failed extraction chunk can require targeted diagnostics";
		streamBatches(assembly.fasta, batch->{
			for(final Read read:batch.list){
				final String contig=RrnaResourceIO.first(read.id);
				for(final Locus locus:assembly.loci){
					if(locus.contig.equals(contig)){locus.contigLength=read.length();}
				}
			}
		});
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Drains native batches before propagating a consumer failure. No custom FASTA
	 * parsing, queue, dispatcher or sequence extraction is implemented here. */
	private static void streamBatches(final String path, final Consumer<ListNum<Read>> action){
		final Streamer reader=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, false, true), null, false, -1, false, true, 1);
		Throwable failure=null;
		reader.start();
		try{
			for(ListNum<Read> batch; (batch=reader.nextList())!=null;){
				if(failure==null){
					try{action.accept(batch);}
					catch(RuntimeException | Error problem){failure=problem;}
				}
				//As in A_SampleStreamerMT, submitted Reads remain borrowed until writer
				//completion; neither the batch nor its payload is recycled here.
			}
		}finally{require(!ReadWrite.closeStream(reader), "Sequence read failed: "+path);}
		if(failure instanceof Error){throw (Error)failure;}
		if(failure!=null){throw (RuntimeException)failure;}
	}

	/** Attempts every acquired writer's shutdown even if an earlier one fails. */
	private static void closeWriters(final Writer[] writers){
		Throwable failure=null;
		for(final Writer writer:writers){
			if(writer==null){continue;}
			try{require(!writer.poisonAndWait(), "Extracted member FASTA output did not close successfully");}
			catch(RuntimeException | Error problem){
				if(failure==null){failure=problem;}
				else if(problem!=failure){failure.addSuppressed(problem);}
			}
		}
		if(failure instanceof Error){throw (Error)failure;}
		if(failure!=null){throw (RuntimeException)failure;}
	}

	/** Multiple or absent positive taxIDs remain unknown; none is synthesized. */
	private int taxid(final int column){
		p.setBounds(column);
		int value=0;
		for(int i=p.a(); i<p.b(); i++){
			final byte b=p.line()[i];
			if(b<'0' || b>'9'){return -1;}
			if(value>(Integer.MAX_VALUE-(b-'0'))/10){return -1;}
			value=value*10+b-'0';
		}
		return value>0 ? value : -1;
	}
	private int member(final int column){
		for(int i=0; i<3; i++){if(p.termEquals(MEMBERS[i], column)){return i;}}
		return -1;
	}
	private String projection(final Assembly assembly){return new File(new File(directory, "projections"), assembly.input+".gff").getPath();}
	private int chunkCount(){return (assemblies.size()-1)/chunkGenomes+1;}
	/** Selection occurs only AFTER full-ledger IDs and all three counts are verified. */
	private boolean wanted(final int member){return selectedMember<0 || member==selectedMember;}
	private int selectedCount(final Assembly assembly){
		if(selectedMember<0){return assembly.loci.size();}
		int count=0;
		for(final Locus locus:assembly.loci){if(wanted(locus.member)){count++;}}
		return count;
	}
	private Output output(final String name){return new Output(new File(directory, name).getPath());}
	private void writeComplete(final String state, final int records){
		try(Output complete=output("complete.tsv")){
			complete.writer.println("state\tselected_genomes\tadmitted_records\tchunk_genomes\tshuffle_seed\tchunks"
				+(selectedMember<0 ? "" : "\tselected_member"));
			buffer.clear().append(state).tab().append(assemblies.size()).tab().append(records)
				.tab().append(chunkGenomes).tab().append(seed).tab().append(chunkCount());
			if(selectedMember>=0){buffer.tab().append(MEMBERS[selectedMember]);}
			complete.writer.print(buffer.nl());
		}
	}
	private static void validatePath(final String path){
		require(new File(path).isAbsolute() && path.indexOf('=')<0 && path.indexOf(',')<0 && path.indexOf('\n')<0,
			"CutGff configuration requires an absolute path without equals, comma or newline: "+path);
	}
	private static boolean samePath(final String a, final String b) throws IOException{return a.equals(b) || new File(a).getCanonicalPath().equals(new File(b).getCanonicalPath());}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Retained Records      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Assembly{
		Assembly(final String id_, final int input_, final String fasta_, final String gff_, final String source_, final int taxid_){
			id=id_;
			input=input_;
			fasta=fasta_;
			gff=gff_;
			source=source_;
			taxid=taxid_;
		}
		final String id, fasta, gff, source;
		final int input, taxid;
		int chunk=-1;
		final ArrayList<Locus> loci=new ArrayList<Locus>();
		final ObjectMap<String, String> coordinates=new ObjectMap<String, String>(String.class, String.class);
	}
	private static final class Locus{
		Locus(final String id_, final Assembly assembly_, final int member_, final String contig_, final int start_, final int stop_, final char strand_, final int flank_){
			id=id_;
			assembly=assembly_;
			member=member_;
			contig=contig_;
			start=start_;
			stop=stop_;
			strand=strand_;
			flank=flank_;
		}
		int length(){return stop-start+1;}
		boolean valid(){return seen==1 && observedLength==length() && headerMatches && sequenceMatches;}
		String coreHeader(){return assembly.taxid>0 ? id+" tid_"+assembly.taxid : id;}
		String header(){return flank>0 ? coreHeader()+" source="+assembly.id : coreHeader();}
		final String id, contig;
		final Assembly assembly;
		final int member, start, stop, flank;
		final char strand;
		int seen=0, observedLength=-1, contigLength=-1, left=-1, right=-1;
		boolean headerMatches=true, sequenceMatches=true;
		byte[] core;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String mode="prepare", manifest, lociPath, directory, sequences;
	private String[] cores;
	private String selectedName;
	private int selectedMember=-1;
	private int expectedGenomes=-1, chunkGenomes=1000, chunk=0, flank=0;
	private long seed=-1;
	private final int[] expected={-1, -1, -1};
	private final ArrayList<Assembly> assemblies=new ArrayList<Assembly>();
	private final IntObjectMap<Assembly> byInput=new IntObjectMap<Assembly>(Assembly.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	private final ObjectMap<String, Locus> byId=new ObjectMap<String, Locus>(String.class, Locus.class);
	private final ArrayList<String> unexpected=new ArrayList<String>();
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private final PrintStream outstream;
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String MANIFEST_HEADER="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
	private static final String LOCUS_HEADER="assembly_id\tinput_id\tcohort\tsource_fasta\tmember\tcontig\tstart\tstop\tstrand\tnative_events\tuntruncated_events\tconstruction_eligible\texact_scored_links";
}
