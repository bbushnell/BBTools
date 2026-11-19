package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.IntHashMap;

/** Prepares the complete recorded source-frame population for one-pass baselines.
 * Adapted from template.A_SampleBasic2's parser/output lifecycle and the native
 * NcrnaDataReconciliation byte prepass. No DNA, library or caller is processed.
 * Every recorded frame and original gene survives in explicit output ledgers;
 * missing source files are unresolved inputs, never zero-call observations.
 * The capture manifest's unused GFF column names placeholders, not reference files.
 * Source availability and size checks establish current metadata only, not DNA
 * identity or an equivalence between different assembly coordinate frames.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareNcrnaRosterBaseline{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, PrepareNcrnaRosterBaseline.class, false);
		try{
			final PrepareNcrnaRosterBaseline tool=new PrepareNcrnaRosterBaseline(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println("NCRNA_COMPLETE_ROSTER_BASELINE_INPUTS_READY time="+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private PrepareNcrnaRosterBaseline(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("originals")){originals=b;}
			else if(a.equals("expectedgenes")){expectedGenes=Parse.parseIntKMG(b);}
			else if(a.equals("expectedframes")){expectedFrames=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown baseline preparation argument: "+arg);}
		}
		input=parser.in1;
		prefix=parser.out1;
		require(input!=null && originals!=null && prefix!=null && expectedGenes>0 && expectedFrames>0
			&& !parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1,
			"Supply in=source_frames.tsv.gz originals=original_calls.tsv.gz out=fresh_prefix expectedgenes=N expectedframes=N without sampling");
		prefix=new File(prefix).getAbsolutePath();
		for(final String suffix:SUFFIXES){require(!new File(prefix+suffix).exists(), "Baseline preparation requires fresh output: "+prefix+suffix);}
	}

	private void process(){
		loadFrames();
		writeOriginals();
		writeInputs();
		writeSummary();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Metadata Input       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadFrames(){
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, FRAME_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==13, "Source-frame population must retain the complete reconciler schema");
				final Frame frame=new Frame(order.size(), p.parseInt(0), p.parseString(1), p.parseString(2),
					p.parseLong(3), p.parseLong(4), p.parseLong(5), p.parseInt(8));
				require(frame.group>=0 && new File(frame.path).isAbsolute() && new File(frame.alias).isAbsolute()
					&& frame.bytes>=0 && frame.bases>0 && frame.contigs>0 && frame.expectedGenes>=0
					&& !frames.contains(frame.path), "Recorded canonical frames require unique absolute paths and saved dimensions");
				require(frame.path.indexOf(',')<0 && frame.path.indexOf('\n')<0 && frame.path.indexOf('\r')<0
					&& frame.path.indexOf('\t')<0, "A complete CallGenes input list cannot silently split a source path: "+frame.path);
				final File source=new File(frame.path);
				frame.available=source.isFile() && source.canRead();
				if(frame.available){
					require(source.length()==frame.bytes, "Source size changed since the recorded coordinate-frame census: "+frame.path);
					available++;
					availableBases+=frame.bases;
				}
				frames.put(frame.path, frame);
				order.add(frame);
			}
		}
		require(order.size()==expectedFrames && available>0, "Declare every recorded canonical frame, including unavailable inputs");
	}

	/** Original records are read once; row IDs and per-frame counts conserve the
	 * complete native-call population independently of source availability. */
	private void writeOriginals(){
		final IntHashMap seen=new IntHashMap();
		try(Input in=new Input(originals); Output out=new Output(prefix+".original_calls.tsv.gz")){
			final byte[] header=in.reader.nextLine();
			require(header!=null && new String(header, StandardCharsets.UTF_8).equals(ORIGINAL_HEADER),
				"Use the complete29-column extraction ledger, not a construction-only member table");
			out.writer.print(row.clear().append(header).append("\tbaseline_input_id\tbaseline_input_state\n"));
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==29 && p.termEquals("true", 7), "Every original native-GA call, including truncations, must be retained");
				final int id=p.parseInt(0);
				final Frame frame=frames.get(p.parseString(2));
				require(id>=0 && !seen.contains(id) && frame!=null && frame.group==p.parseInt(1)
					&& p.parseInt(4)>0 && p.parseInt(5)>=p.parseInt(4)
					&& (p.termEquals("+", 6) || p.termEquals("-", 6)),
					"Each original gene needs a unique ID and its exact recorded source frame, strand and inclusive bounds");
				seen.put(id, 1);
				frame.genes++;
				genes++;
				if(!frame.available){unavailableGenes++;}
				out.writer.print(row.clear().append(line).append("\tframe_").append(frame.id).tab()
					.append(frame.available ? "AVAILABLE_RECORDED_FRAME" : "UNRESOLVED_SOURCE_UNAVAILABLE").nl());
			}
		}
		require(genes==expectedGenes, "Baseline recall must conserve every original call, regardless of construction eligibility");
		for(final Frame frame:order){
			require(frame.genes==frame.expectedGenes, "Original-gene counts must agree with each canonical frame's accepted census: "+frame.path);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void writeInputs(){
		try(Output ledger=new Output(prefix+".frames.tsv"); Output manifest=new Output(prefix+".reference.tsv");
			Output config=new Output(prefix+".inputs.config")){
			ledger.writer.println("input_id\tassembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\toriginal_genes\tinput_state");
			manifest.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			config.writer.print("in=");
			int written=0;
			for(final Frame frame:order){
				row.clear().append("frame_").append(frame.id).tab().append(frame.group).tab()
					.append(frame.path).tab().append(frame.alias).tab().append(frame.bytes).tab().append(frame.bases).tab()
					.append(frame.contigs).tab().append(frame.genes).tab()
					.append(frame.available ? "AVAILABLE_RECORDED_FRAME" : "UNRESOLVED_SOURCE_UNAVAILABLE").nl();
				ledger.writer.print(row);
				if(!frame.available){continue;}
				manifest.writer.print(row.clear().append("frame_").append(frame.id).tab().append(frame.path).tab()
					.append(prefix).append(".unused_reference_frame_").append(frame.id).append(".gff\n"));
				row.clear();
				if(written>0){row.append(',');}
				config.writer.print(row.append(frame.path));
				written++;
			}
			config.writer.print("\n");
			require(written==available, "Capture inputs include all and only currently readable recorded frames, including zero-gene frames");
		}
	}

	private void writeSummary(){
		assert(genes==expectedGenes && order.size()==expectedFrames) : "Completion records are published only after whole-population conservation";
		try(Output out=new Output(prefix+".summary.tsv")){
			out.writer.print(row.clear().append("metric\tcount\nrecorded_source_frames\t").append(order.size())
				.append("\navailable_source_frames\t").append(available)
				.append("\nunavailable_source_frames\t").append(order.size()-available)
				.append("\noriginal_genes\t").append(genes)
				.append("\ngenes_on_unavailable_sources\t").append(unavailableGenes)
				.append("\navailable_nucleotide_bases\t").append(availableBases).nl());
		}
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Frame{
		Frame(int id_, int group_, String path_, String alias_, long bytes_, long bases_, long contigs_, int genes_){
			id=id_;
			group=group_;
			path=path_;
			alias=alias_;
			bytes=bytes_;
			bases=bases_;
			contigs=contigs_;
			expectedGenes=genes_;
		}
		final int id, group, expectedGenes;
		final String path, alias;
		final long bytes, bases, contigs;
		boolean available;
		int genes;
	}
	private String input, originals, prefix;
	private int expectedGenes, expectedFrames, available, genes, unavailableGenes;
	private long availableBases;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder();
	private final ObjectMap<String, Frame> frames=new ObjectMap<String, Frame>(String.class, Frame.class);
	private final ArrayList<Frame> order=new ArrayList<Frame>();
	private static final String[] SUFFIXES={".original_calls.tsv.gz", ".frames.tsv", ".reference.tsv", ".inputs.config", ".summary.tsv"};
	public static final String FRAME_HEADER="assembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\tcohort_input_rows\tframes_in_group\tnative_GA_exact_groups\tassembly_has_native_GA\tcoordinate_comparison_needed\tsource_available\treason";
	public static final String ORIGINAL_HEADER="locus_id\tassembly_group\tcoordinate_frame\tcontig\tstart1\tstop1\tstrand\tqualified\tprior_annotation_overlap\traw_events\texact_members\treciprocal80_members\tdisposition\tmerged_locus_id\tconstruction_representative\tconstruction_state\tnative_hit_count\tnative_truncation_bits\tnative_model_from\tnative_model_to\tnative_model_length\tincomplete_model_span\tnative_construction_eligibility\tavailable_left_bases\tavailable_right_bases\tadjacent_left_N_run\tadjacent_right_N_run\tgenomic_context\tnative_printed_bits";
}
