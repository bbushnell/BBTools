package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Plans the original gene sequences still needed after construction extraction.
 * Adapted from template.A_SampleBasic2 and PrepareNcrnaRosterBaseline.
 * Reads only the accepted prepared frame/original ledgers. Existing eligible
 * originals reuse their recorded same-group exact construction representatives;
 * a one-to-one original/component map is required before making that inference.
 * Excluded truncated originals become evaluation-only retrieval requests.
 * Source bytes, dimensions and availability are recorded metadata, not a new
 * file check or DNA identity claim. A retriever must validate the selected
 * canonical source bindings and preserve every requested original and strand.
 * No sequence, taxonomy, seed generation, caller or coverage is processed.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaSeedInputs{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, PrepareNcrnaSeedInputs.class, false);
		try{
			final PrepareNcrnaSeedInputs tool=new PrepareNcrnaSeedInputs(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println("NCRNA_SEED_SEQUENCE_INPUT_PLAN_READY_METADATA_ONLY time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private PrepareNcrnaSeedInputs(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedframes")){expectedFrames=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmembers")){expectedMembers=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmissing")){expectedMissing=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown sequence-plan argument: "+arg);}
		}
		input=parser.in1;
		prefix=parser.out1;
		require(input!=null && prefix!=null && expectedGenes>0 && expectedFrames>0
			&& expectedMembers>0 && expectedMembers<=expectedGenes && expectedMissing>=0
			&& expectedMissing<=expectedGenes && !parser.overwrite && !parser.append
			&& parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Supply in=prepared_population_prefix out=fresh_prefix expectedgenes=N expectedframes=N expectedmembers=N expectedmissing=N without sampling");
		require(Tools.testInputFiles(false, true, input+".frames.tsv", input+".original_calls.tsv.gz"),
			"The accepted prepared population is required; do not substitute construction-only gene metadata");
		prefix=new File(prefix).getAbsolutePath();
		for(final String suffix:SUFFIXES){require(!new File(prefix+suffix).exists(), "Sequence-plan output must be fresh: "+prefix+suffix);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Metadata Input       ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		loadFrames();
		loadOriginals();
		writeOutputs();
	}

	private void loadFrames(){
		try(Input in=new Input(input+".frames.tsv")){
			RrnaResourceIO.header(in.reader, FRAME_HEADER, input+".frames.tsv");
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==9 && p.termStartsWith("frame_", 0) && p.parseInt(0, 6)==frames.size(),
					"Prepared frame IDs must retain PrepareNcrnaRosterBaseline's complete canonical order");
				final boolean available=p.termEquals("AVAILABLE_RECORDED_FRAME", 8);
				require(available || p.termEquals("UNRESOLVED_SOURCE_UNAVAILABLE", 8),
					"Missing sources must retain the explicit prepared availability state");
				final Frame frame=new Frame(frames.size(), p.parseInt(1), p.parseString(2),
					p.parseLong(4), p.parseLong(5), p.parseInt(6), p.parseInt(7), available);
				require(frame.group>=0 && frame.path.startsWith("/") && frame.bytes>=0
					&& frame.bases>0 && frame.contigs>0 && frame.expectedGenes>=0,
					"Each recorded source frame needs an absolute canonical path and saved dimensions");
				frames.add(frame);
			}
		}
		require(frames.size()==expectedFrames,
			"Sequence planning must retain every prepared frame, including zero-gene and unavailable frames: observed="+frames.size()+" expected="+expectedFrames);
	}

	private void loadOriginals(){
		final boolean[] seen=new boolean[expectedGenes], representatives=new boolean[expectedGenes];
		final int[] representativeByMerged=new int[expectedGenes], groupByMerged=new int[expectedGenes];
		Arrays.fill(representativeByMerged, -2);
		try(Input in=new Input(input+".original_calls.tsv.gz")){
			RrnaResourceIO.header(in.reader, ORIGINAL_HEADER, input+".original_calls.tsv.gz");
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==31 && p.termEquals("true", 7) && p.termStartsWith("frame_", 29),
					"Every original native-GA gene must retain its complete prepared identity and eligibility fields");
				final int id=p.parseInt(0), merged=p.parseInt(13), representative=p.parseInt(14), frameId=p.parseInt(29, 6);
				require(id>=0 && id<expectedGenes && !seen[id] && merged>=0 && merged<expectedGenes
					&& representativeByMerged[merged]==-2 && frameId>=0 && frameId<frames.size(),
					"Exact representative reuse requires every original and merged ordinal once; component merges or duplicate IDs need a different sequence mapping: original="+id+" merged="+merged);
				final Frame frame=frames.get(frameId);
				require(p.parseInt(1)==frame.group && p.termEquals(frame.path, 2) && p.termEquals(frame.state(), 30),
					"Original genes must bind their exact canonical frame, group and recorded availability: original="+id+" frame="+frameId);
				final int start=p.parseInt(4), stop=p.parseInt(5), truncation=p.parseInt(17);
				require(!p.termEquals("", 3) && start>0 && stop>=start && (p.termEquals("+", 6) || p.termEquals("-", 6)),
					"Evaluation retrieval must preserve a nonempty contig, positive inclusive bounds and explicit strand: original="+id);
				seen[id]=true;
				representativeByMerged[merged]=representative;
				groupByMerged[merged]=frame.group;
				frame.observedGenes++;
				originalsRead++;
				if(representative>=0){
					require(representative<expectedGenes && truncation==0 && p.termEquals("ELIGIBLE_NATIVE_HIT", 22)
						&& p.termEquals(representative==merged ? "CONSTRUCTION_REPRESENTATIVE" : "EXACT_SEQUENCE_DUPLICATE_FOR_CONSTRUCTION_ONLY", 15),
						"A reusable saved sequence must be an eligible representative or its recorded exact construction copy: original="+id);
					if(!representatives[representative]){
						representatives[representative]=true;
						distinctRepresentatives++;
					}
					withSavedSequence++;
				}else{
					require(representative==-1 && truncation>0 && p.termEquals("EXCLUDED_NATIVE_PARTIAL_OR_TRUNCATED", 15)
						&& !p.termEquals("ELIGIBLE_NATIVE_HIT", 22),
						"Only resolved excluded truncations may become missing evaluation-sequence requests; do not reinterpret an extraction failure: original="+id);
					missing.add(new Missing(id, merged, frame, p.parseByteArray(3), start, stop, p.termEquals("-", 6), truncation));
					frame.missingGenes++;
					requestedBases=Math.addExact(requestedBases, (long)stop-start+1);
				}
			}
		}
		require(originalsRead==expectedGenes && distinctRepresentatives==expectedMembers && missing.size()==expectedMissing,
			"All-original and saved-member counts must be conserved before publishing requests: originals="+originalsRead+"/"+expectedGenes
			+" members="+distinctRepresentatives+"/"+expectedMembers+" missing="+missing.size()+"/"+expectedMissing);
		for(int merged=0; merged<expectedGenes; merged++){
			final int representative=representativeByMerged[merged];
			require(representative!=-2 && (representative<0 || representativeByMerged[representative]==representative
				&& groupByMerged[merged]==groupByMerged[representative]),
				"A recorded exact copy must resolve directly to an actual representative in the SAME accepted group: merged="+merged+" representative="+representative);
		}
		for(final Frame frame:frames){
			require(frame.observedGenes==frame.expectedGenes,
				"Original-gene counts must conserve each recorded canonical frame: frame="+frame.id+" observed="+frame.observedGenes+" expected="+frame.expectedGenes);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeOutputs(){
		assert(originalsRead==withSavedSequence+missing.size())
			: "Every original is either linked to saved construction DNA or retained as an evaluation-only sequence request";
		final ByteBuilder row=new ByteBuilder();
		try(Output requests=new Output(prefix+".requests.tsv"); Output sources=new Output(prefix+".sources.tsv");
			Output summary=new Output(prefix+".summary.tsv")){
			requests.writer.println(REQUEST_HEADER);
			for(final Missing gene:missing){
				row.clear().append(gene.id).tab().append(gene.merged).append("\tframe_").append(gene.frame.id)
					.tab().append(gene.frame.group).tab().append(gene.frame.path).tab().append(gene.contig)
					.tab().append(gene.start).tab().append(gene.stop).tab().append(gene.minus ? '-' : '+')
					.tab().append(gene.truncation).append("\ttrue\t").append(gene.frame.state()).nl();
				requests.writer.print(row);
			}
			sources.writer.println(SOURCE_HEADER);
			for(final Frame frame:frames){
				if(frame.missingGenes==0){continue;}
				missingSources++;
				if(!frame.available){unavailableMissingSources++;}
				sourceBytes=Math.addExact(sourceBytes, frame.bytes);
				sourceBases=Math.addExact(sourceBases, frame.bases);
				sources.writer.print(row.clear().append("frame_").append(frame.id).tab().append(frame.group)
					.tab().append(frame.path).tab().append(frame.bytes).tab().append(frame.bases).tab().append(frame.contigs)
					.tab().append(frame.missingGenes).tab().append(frame.state()).nl());
			}
			summary.writer.println("metric\tcount");
			metric(summary, row, "original_genes", originalsRead);
			metric(summary, row, "recorded_source_frames", frames.size());
			metric(summary, row, "originals_with_saved_sequences", withSavedSequence);
			metric(summary, row, "distinct_saved_representatives", distinctRepresentatives);
			metric(summary, row, "missing_original_sequences", missing.size());
			metric(summary, row, "missing_source_frames", missingSources);
			metric(summary, row, "recorded_unavailable_missing_sources", unavailableMissingSources);
			metric(summary, row, "requested_evaluation_bases", requestedBases);
			metric(summary, row, "planned_source_file_bytes", sourceBytes);
			metric(summary, row, "planned_source_nucleotide_bases", sourceBases);
		}
	}

	private static void metric(final Output out, final ByteBuilder row, final String name, final long value){
		out.writer.print(row.clear().append(name).tab().append(value).nl());
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Frame{
		Frame(final int id_, final int group_, final String path_, final long bytes_, final long bases_,
			final int contigs_, final int genes_, final boolean available_){
			id=id_;
			group=group_;
			path=path_;
			bytes=bytes_;
			bases=bases_;
			contigs=contigs_;
			expectedGenes=genes_;
			available=available_;
		}
		String state(){return available ? "AVAILABLE_RECORDED_FRAME" : "UNRESOLVED_SOURCE_UNAVAILABLE";}
		final int id, group, contigs, expectedGenes;
		final String path;
		final long bytes, bases;
		final boolean available;
		int observedGenes, missingGenes;
	}

	private static final class Missing{
		Missing(final int id_, final int merged_, final Frame frame_, final byte[] contig_, final int start_,
			final int stop_, final boolean minus_, final int truncation_){
			assert(frame_!=null && contig_!=null && contig_.length>0 && start_>0 && stop_>=start_ && truncation_>0)
				: "A retained evaluation request needs the original exact source, strand-oriented bounds and native truncation evidence";
			id=id_;
			merged=merged_;
			frame=frame_;
			contig=contig_;
			start=start_;
			stop=stop_;
			minus=minus_;
			truncation=truncation_;
		}
		final int id, merged, start, stop, truncation;
		final Frame frame;
		final byte[] contig;
		final boolean minus;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String input, prefix;
	private int expectedGenes=-1, expectedFrames=-1, expectedMembers=-1, expectedMissing=-1;
	private int originalsRead, withSavedSequence, distinctRepresentatives, missingSources, unavailableMissingSources;
	private long requestedBases, sourceBytes, sourceBases;
	private final ArrayList<Frame> frames=new ArrayList<Frame>();
	private final ArrayList<Missing> missing=new ArrayList<Missing>();
	private final LineParser1 p=new LineParser1('\t');
	private static final String[] SUFFIXES={".requests.tsv", ".sources.tsv", ".summary.tsv"};
	public static final String FRAME_HEADER="input_id\tassembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\toriginal_genes\tinput_state";
	public static final String ORIGINAL_HEADER=PrepareNcrnaRosterBaseline.ORIGINAL_HEADER+"\tbaseline_input_id\tbaseline_input_state";
	public static final String REQUEST_HEADER="locus_id\tmerged_locus_id\tinput_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tnative_truncation_bits\tevaluation_only\trecorded_input_state";
	public static final String SOURCE_HEADER="input_id\tassembly_group\tcoordinate_frame\tfile_bytes\tnucleotide_bases\tfasta_records\trequested_genes\trecorded_input_state";
}
