package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CaptureNcrnaWindows;
import prok.ConservedRnaSeedIndex;
import prok.GradeNcrnaReference;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaJointPolicy;
import prok.NcrnaCallerProcess;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Exercises the native caller and capture writer on an asymmetric minus-strand
 * seed, with an independent literal coordinate expectation and a zero-window input.
 * The live capture assertion checks that the callback borrows its current Read's
 * array even while GeneCaller has reverse-complemented that array in place.
 * Adapted from A_SampleBasic2's output-only parser and native sequence I/O idioms.
 * @author Brian Bushnell, Raiden
 */
public final class CaptureNcrnaWindowsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CaptureNcrnaWindowsTest.class, false);
		try{
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "This fixture requires -ea to exercise borrowed-array identity on the caller worker");
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Supply out=fresh_fixture_directory without input or overwrite");
			new CaptureNcrnaWindowsTest(new File(parser.out1).getAbsoluteFile()).run();
		}finally{Shared.closeStream(pp.outstream);}
	}

	private CaptureNcrnaWindowsTest(final File root_){
		root=root_;
		check(root.mkdir(), "Fixture directory must be fresh: "+root);
	}

	/** Only a real forward seed is used; its reverse complement is absent from the seed set. */
	private void run(){
		final NcrnaFamily family=NcrnaFamilyLoader.load("cobalamin_rf00174");
		final byte[] forward=new byte[4096];
		Arrays.fill(forward, (byte)'N');
		final byte[] seed=minusSeed(family.kmerSet);
		System.arraycopy(seed, 0, forward, 500, seed.length);
		final byte[] oriented=forward.clone();
		AminoAcid.reverseComplementBasesInPlace(oriented);
		final ConservedRnaSeedIndex index=new ConservedRnaSeedIndex(new LongHashSet[]{family.kmerSet});
		check(index.scan(forward).hits(0).length==0, "The forward strand must contain no target seed");
		check(Arrays.equals(index.scan(oriented).hits(0), new int[]{3587}),
			"A17-mer at forward500..516 has reverse-complement center3587 on a4096-base contig");
		write("minus.fa", new ByteBuilder(">minus_contig\n").append(forward).nl());
		write("no_window.fa", new ByteBuilder(">no_seed\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"));
		final ByteBuilder manifest=new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl();
		manifest.append("minus\t").append(path("minus.fa")).tab().append(path("unused_minus.gff")).nl();
		manifest.append("empty\t").append(path("no_window.fa")).tab().append(path("unused_empty.gff")).nl();
		write("reference.tsv", manifest);
		checkDuplicateManifest();
		checkTaxonomyOverride("taxonomy=f\ntax=t\n", "short_alias");
		checkTaxonomyOverride("tax=f\ntaxonomy=t\n", "long_alias");
		final ByteBuilder config=new ByteBuilder();
		config.append("cdsonly=t\ncallcds=f\nncrna=f\ncobalamin=t\ntmrna=f\nsixs=f\nfmn=f\n")
			.append("passes=1\ntax=f\nncrnafamily=cobalamin_rf00174\nncrnapad=100\n")
			.append("ncrnaseedminhits=1\ncobalaminminlen=20\n")
			.append("ncrnaadaptiveminhits=f\nncrnaindexminhits=1000\nt=2\n");
		write("caller.config", withInputs(config, "minus.fa"));
		CaptureNcrnaWindows.main(new String[]{"caller="+path("caller.config"),
			"reference="+path("reference.tsv"), "out="+path("capture"),
			"expectedmodels="+family.library.length, "t=2"});
		checkMetadata();
		checkSequences(oriented);
		checkInputs();
		checkSummary();
		checkSelectedResults(family, config);
		checkReplayContext(family, config, oriented);
		System.err.println("NCRNA_WINDOW_CAPTURE_DRIVER_TEST_PASS");
	}

	/** The same real caller must retain completed work and both declared inputs
	 * without using the zero-window input as an implicit population filter. */
	private void checkSelectedResults(final NcrnaFamily family, final ByteBuilder config){
		final int[] voters=new int[family.library.length];
		Arrays.fill(voters, 1);
		NcrnaDdlProfileResource.writeSelected(path("selected.tsv"), new char[voters.length][512], voters,
			family.modelNames, family.library, 512, new LongHashSet(4));
		write("selected_policy.tsv", new ByteBuilder(NcrnaJointPolicy.HEADER).nl()
			.append("capture\tcobalamin_rf00174\tadaptive\t384\t2\t0\t0\t0\t8\n"));
		final ByteBuilder selected=new ByteBuilder().append(config.toBytes())
			.append("cobalaminddl=").append(path("selected.tsv")).nl()
			.append("ncrnapolicy=").append(path("selected_policy.tsv")).nl();
		write("selected.config", withInputs(selected, "minus.fa"));
		runCapture("verified", "selected.config", "reference.tsv", family.library.length);
		long alignments;
		try(Input input=new Input(path("verified.results.tsv"))){
			check(Arrays.equals(input.reader.nextLine(), CaptureNcrnaWindows.RESULT_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"Completed-window results must declare their pre-resolution status and counted work");
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.terms()==26 && p.parseInt(0)==0 && p.termEquals("minus", 1)
				&& p.parseInt(5)==1 && p.parseInt(6)==4096 && p.parseInt(7)==1
				&& p.parseInt(8)==3487 && p.parseInt(9)==3704 && p.parseInt(10)==391
				&& p.parseInt(11)==608 && p.termEquals("i0e0", 12),
				"The completed minus-strand event must retain its exact original input, capture ID and both coordinate systems");
			alignments=p.parseLong(17);
			check(alignments==p.parseLong(14)+p.parseLong(15)+p.parseLong(16)
				&& alignments>0 && p.parseInt(25)==17 && p.parseInt(24)>=0 && p.parseInt(24)<=17,
				"Completed work must conserve its primary/edge, rescue and other phases; GC must count the17known bases, not surrounding Ns");
			check(input.reader.nextLine()==null, "The zero-seed input must produce no invented scheduled-window result");
		}
		long finalRows=0;
		try(Input input=new Input(path("verified.attempts.tsv"))){
			check(Arrays.equals(input.reader.nextLine(), CaptureNcrnaWindows.ATTEMPT_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"Future captures must preserve raw attempt identity/span instead of discarding rejected evidence");
			final LineParser1 p=new LineParser1('\t');
			long observed=0;
			int rows=0;
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				check(p.terms()==13 && p.termEquals("minus", 1) && p.termEquals("i0e0", 2), "Raw stage summaries must retain their exact captured input/event");
				observed+=p.parseLong(4);
				if(p.termEquals("PRIMARY_RAW", 3)){
					check(p.parseLong(4)>0 && !p.termEquals("NA", 7) && p.parseInt(8)>0
						&& p.parseInt(9)>=391 && p.parseInt(10)<=608, "Rejected primary evidence must retain actual identity, span and original-forward minus coordinates");
				}
				rows++;
			}
			check(rows==3 && observed>0 && observed<=alignments, "Raw stage counts cannot invent alignment work or omit the three explicit stage states");
		}
		try(Input input=new Input(path("verified.calls.gff"))){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length>0 && line[0]!='#'){finalRows++;}
			}
		}
		try(Input input=new Input(path("verified.inputs.tsv"))){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.terms()==10 && p.termEquals("minus", 1) && p.parseInt(3)==1 && p.parseInt(4)==4096
				&& p.parseInt(5)==1 && p.parseInt(6)==1 && p.parseLong(7)==alignments && p.parseLong(9)==finalRows,
				"Input totals must conserve the actual completed window, all its alignment work and final native GFF rows separately");
			p.set(input.reader.nextLine());
			check(p.terms()==10 && p.termEquals("empty", 1) && p.parseInt(3)==1 && p.parseInt(4)==32,
				"The zero-window input must remain explicitly present after verification");
			for(int field=5; field<10; field++){check(p.parseLong(field)==0, "A zero-window input must retain real zero work and zero calls");}
			check(input.reader.nextLine()==null, "Both declared inputs must appear exactly once");
		}
		checkAcceptedResults(family, selected);
	}

	/** Positive and zero-window inputs deliberately reuse one contig name. The
	 * final TSV must retain the owning input while matching native GFF geometry. */
	private void checkAcceptedResults(final NcrnaFamily family, final ByteBuilder config){
		final ConservedRnaSeedIndex index=new ConservedRnaSeedIndex(new LongHashSet[]{family.kmerSet});
		byte[] model=null;
		for(final byte[] sequence:family.library){
			if(index.scan(sequence).hits(0).length>0){model=sequence;break;}
		}
		check(model!=null, "The real library must supply a consensus containing a declared discovery seed");
		final byte[] genome=new byte[model.length+400];
		Arrays.fill(genome, (byte)'N');
		System.arraycopy(model, 0, genome, 200, model.length);
		write("positive.fa", new ByteBuilder(">no_seed\n").append(genome).nl());
		write("positive_reference.tsv", new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl()
			.append("positive\t").append(path("positive.fa")).tab().append(path("unused_positive.gff")).nl()
			.append("empty\t").append(path("no_window.fa")).tab().append(path("unused_empty.gff")).nl());
		write("all_policy.tsv", new ByteBuilder(NcrnaJointPolicy.HEADER).nl()
			.append("all\tcobalamin_rf00174\tadaptive\t2147483647\t2\t0\t0\t0\t2147483647\n"));
		write("positive.config", withInputs(config, "positive.fa")
			.append("ncrnapad=512\ncobalaminddlcap=all\nncrnapolicy=").append(path("all_policy.tsv")).nl());
		runCapture("positive", "positive.config", "positive_reference.tsv", family.library.length);
		final LongHashSet coordinates=new LongHashSet(16);
		long gffRows=0, boundRows=0;
		try(Input input=new Input(path("positive.calls.gff"))){
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){continue;}
				p.set(line);
				check(p.terms()==9 && p.termEquals("no_seed", 0) && p.termEquals("RNA", 2),
					"The positive fixture must use its real native contig and target RNA output");
				check(coordinates.add(coordinate(p.parseInt(3)-1, p.parseInt(4)-1, p.termEquals("+", 6) ? 0 : 1)),
					"The finite positive fixture must retain unique native coordinate rows");
				gffRows++;
			}
		}
		try(Input input=new Input(path("positive.final_calls.tsv"))){
			check(Arrays.equals(input.reader.nextLine(), CaptureNcrnaWindows.FINAL_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"Retained calls must explicitly declare their original input binding");
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==12 && p.parseInt(0)==0 && p.termEquals("positive", 1)
					&& p.termEquals(path("positive.fa"), 2) && p.parseInt(3)==0 && p.termEquals("no_seed", 4)
					&& coordinates.remove(coordinate(p.parseInt(5), p.parseInt(6), p.parseInt(7))),
					"Each retained call must bind the positive file and match native GFF geometry; the same contig name in the empty file must not steal it");
				boundRows++;
			}
		}
		check(gffRows>0 && boundRows==gffRows && coordinates.size()==0,
			"At least one real seeded accepted call must survive, and native GFF and input-bound retained rows must conserve each other");
	}

	private static long coordinate(final int start, final int stop, final int strand){
		check(start>=0 && stop>=start && (strand==0 || strand==1), "Native coordinate keys require zero-based inclusive spans and explicit strand");
		return ((long)start<<33)|((long)stop<<1)|strand;
	}

	/** Shared policy contains no inputs; each scenario adds its pair exactly once. */
	private ByteBuilder withInputs(final ByteBuilder policy, final String first){
		check(first!=null && !first.equals("no_window.fa"), "Each scenario needs one distinct first input and its shared zero-window input");
		return new ByteBuilder().append(policy.toBytes()).append("in=").append(path(first)).append(',')
			.append(path("no_window.fa")).nl();
	}

	/** Policy changes run in fresh native JVMs through the installed shell launcher. */
	private void runCapture(final String name, final String caller, final String reference, final int models){
		runCapture(name, caller, reference, models, false);
	}

	private void runCapture(final String name, final String caller, final String reference,
			final int models, final boolean context){
		write(name+".driver.config", new ByteBuilder("caller=").append(path(caller)).nl()
			.append("reference=").append(path(reference)).nl().append("out=").append(path(name)).nl()
			.append("expectedmodels=").append(models).append("\nresults=t\nt=2\nreplaycontext=").append(context).nl());
		final String launcher=Data.findPath("?dev/capture_ncrna_windows.sh", false);
		check(launcher!=null, "The finite driver must use its installed native capture launcher");
		NcrnaCallerProcess.run(new File(launcher), new File(path(name+".driver.config")), 1024);
	}

	/** Read the actual new writer's products on both strands, independently of
	 * replay's manual context slicing. Legacy candidate-only capture stays valid. */
	private void checkReplayContext(final NcrnaFamily family, final ByteBuilder config, final byte[] oriented){
		check(!new File(path("capture.context.tsv")).exists() && !new File(path("capture.context.fa.gz")).exists(),
			"Default candidate-only capture must neither require supported replay routes nor create context sidecars");
		write("context_plus.fa", new ByteBuilder(">minus_contig\n").append(oriented).nl());
		write("context_plus_reference.tsv", new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl()
			.append("plus\t").append(path("context_plus.fa")).tab().append(path("unused_plus.gff")).nl()
			.append("empty\t").append(path("no_window.fa")).tab().append(path("unused_empty.gff")).nl());
		final ByteBuilder selected=new ByteBuilder().append(config.toBytes())
			.append("cobalaminddl=").append(path("selected.tsv")).nl()
			.append("ncrnapolicy=").append(path("selected_policy.tsv")).nl();
		write("context_plus.config", withInputs(selected, "context_plus.fa"));
		for(int strand=0; strand<2; strand++){
			final String prefix=strand==0 ? "context_plus" : "context_minus";
			runCapture(prefix, strand==0 ? "context_plus.config" : "selected.config",
				strand==0 ? "context_plus_reference.tsv" : "reference.tsv", family.library.length, true);
			try(Input input=new Input(path(prefix+".context.tsv"))){
				check(Arrays.equals(input.reader.nextLine(), CaptureNcrnaWindows.CONTEXT_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
					"Replay context metadata must explicitly separate its bounds from candidate bounds");
				final LineParser1 p=new LineParser1('\t');
				p.set(input.reader.nextLine());
				check(p.terms()==4 && p.termEquals("i0e0", 0) && p.parseInt(1)==3477
					&& p.parseInt(2)==3714 && p.parseInt(3)==10 && input.reader.nextLine()==null,
					"The actual context writer must retain exactly ten real bases around the independent3487..3704candidate");
			}
			final Read context=readOne(path(prefix+".context.fa.gz"));
			check(context.id.equals("i0e0") && context.length()==238,
				"The real context FASTA must bind its one exact event and extended length");
			for(int base=0; base<context.length(); base++){
				check(context.bases[base]==oriented[3477+base], "Real context must preserve the original caller-strand bytes on both strands");
			}
			try(Input input=new Input(path(prefix+".windows.tsv"))){
				input.reader.nextLine();
				final LineParser1 p=new LineParser1('\t');
				p.set(input.reader.nextLine());
				check(p.parseInt(6)==strand && p.parseInt(9)==3487 && p.parseInt(10)==3704
					&& p.parseInt(14)==218 && input.reader.nextLine()==null,
					"Adding a context sidecar must not enlarge candidate windows or change their strand/interval");
			}
		}
	}

	/** The last native alias must disable taxonomy before outputs or caller startup. */
	private void checkTaxonomyOverride(final String flags, final String name){
		write(name+".config", new ByteBuilder(flags));
		boolean rejected=false;
		try{
			CaptureNcrnaWindows.main(new String[]{"caller="+path(name+".config"),
				"reference="+path("reference.tsv"), "out="+path(name)});
		}catch(IllegalArgumentException failure){
			rejected="Capture requires final taxonomy=f or tax=f in the recorded caller configuration".equals(failure.getMessage());
		}
		check(rejected, "A final true taxonomy alias must fail at the capture configuration guard: "+name);
		for(final String suffix:new String[]{".windows.tsv", ".windows.fa.gz", ".inputs.tsv", ".summary.tsv"}){
			check(!new File(path(name+suffix)).exists(), "Taxonomy rejection must precede native capture output creation");
		}
	}

	/** Chooses the smallest eligible native DNA key, then writes its reverse complement. */
	private static byte[] minusSeed(final LongHashSet seeds){
		long selected=Long.MAX_VALUE;
		for(final long key:seeds.toArray()){
			if(key<selected && !seeds.contains(AminoAcid.reverseComplementBinaryFast(key, 17))){selected=key;}
		}
		check(selected!=Long.MAX_VALUE, "The real target resource needs one strand-specific17-mer for this fixture");
		final byte[] result=new byte[17];
		long key=selected;
		for(int i=result.length-1; i>=0; i--){
			result[i]=AminoAcid.numberToBase[(int)(key&3)];
			key>>>=2;
		}
		AminoAcid.reverseComplementBasesInPlace(result);
		return result;
	}

	/** Different IDs and GFFs must not hide an exact duplicate FASTA binding. */
	private void checkDuplicateManifest(){
		final ByteBuilder rows=new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl();
		rows.append("one\t").append(path("minus.fa")).tab().append(path("one.gff")).nl();
		rows.append("two\t").append(path("minus.fa")).tab().append(path("two.gff")).nl();
		write("duplicate.tsv", rows);
		boolean rejected=false;
		try{GradeNcrnaReference.loadManifest(path("duplicate.tsv"));}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().startsWith("Duplicate exact FASTA path");}
		check(rejected, "The existing manifest loader must reject duplicate FASTA paths before capture");
	}

	/** Native seed center3587 and pad100 imply inclusive3487..3704, forward391..608. */
	private void checkMetadata(){
		try(Input input=new Input(path("capture.windows.tsv"))){
			check(Arrays.equals(input.reader.nextLine(), CaptureNcrnaWindows.WINDOW_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"Window metadata must declare both coordinate systems");
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.terms()==15 && p.termEquals("i0e0", 0) && p.termEquals("minus", 2)
				&& p.termEquals("minus_contig", 5), "The event must retain its original input and contig binding");
			check(p.parseInt(1)==0 && p.parseInt(4)==0 && p.parseInt(6)==1 && p.parseInt(7)==4096
				&& p.parseInt(8)==1, "One minus-strand primary event must come from the first original contig");
			check(p.parseInt(9)==3487 && p.parseInt(10)==3704 && p.parseInt(11)==391
				&& p.parseInt(12)==608 && p.parseInt(14)==218,
				"Oriented and original-forward intervals must match the independent asymmetric coordinate oracle");
			check(input.reader.nextLine()==null, "The zero-seed input and rejecting index must create no extra captured events");
		}
	}

	/** A native FASTA reader verifies the exact strand-oriented captured bytes. */
	private void checkSequences(final byte[] oriented){
		final Read captured=readOne(path("capture.windows.fa.gz"));
		check(captured.id.equals("i0e0") && captured.length()==218, "Captured FASTA must bind exactly one metadata event");
		for(int i=0; i<captured.length(); i++){
			check(captured.bases[i]==oriented[3487+i], "Captured bytes must use the declared strand-oriented interval at offset"+i);
		}
	}

	private void checkInputs(){
		try(Input input=new Input(path("capture.inputs.tsv"))){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.terms()==6 && p.termEquals("minus", 1) && p.parseInt(3)==1
				&& p.parseInt(4)==4096 && p.parseInt(5)==1, "First input must conserve its contig, bases and captured event");
			p.set(input.reader.nextLine());
			check(p.terms()==6 && p.termEquals("empty", 1) && p.parseInt(3)==1
				&& p.parseInt(4)==32 && p.parseInt(5)==0, "An input with zero candidate windows must remain in the input denominator");
			check(input.reader.nextLine()==null, "Input accounting must contain exactly the two declared inputs");
		}
	}

	private void checkSummary(){
		boolean coordinates=false, retries=false, inputs=false, windows=false;
		try(Input input=new Input(path("capture.summary.tsv"))){
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				check(p.terms()==2, "Capture summary uses item/value pairs");
				if(p.termEquals("coordinates", 0)){coordinates=p.termEquals("ZERO_BASED_INCLUSIVE", 1);}
				if(p.termEquals("edge_retries_separately_observed", 0)){retries=p.termEquals("false", 1);}
				if(p.termEquals("inputs", 0)){inputs=p.parseInt(1)==2;}
				if(p.termEquals("captured_primary_windows", 0)){windows=p.parseInt(1)==1;}
			}
		}
		check(coordinates && retries && inputs && windows, "Final summary must state the coordinate, unobserved-retry and population contracts");
	}

	/** Drains every native batch so sequence-owner closure can report read errors. */
	private static Read readOne(final String path){
		final FileFormat format=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		Read result=null;
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1))){
			input.stream.start();
			for(ListNum<Read> list=input.stream.nextList(); list!=null; list=input.stream.nextList()){
				for(final Read read:list.list){
					check(result==null && read.mate==null, "Fixture FASTA must contain exactly one unpaired contig");
					result=read;
				}
			}
		}
		check(result!=null, "Fixture FASTA must contain its declared contig: "+path);
		return result;
	}

	private void write(final String name, final ByteBuilder bytes){
		try(Output output=new Output(path(name))){output.writer.print(bytes);}
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
	private final File root;
}
