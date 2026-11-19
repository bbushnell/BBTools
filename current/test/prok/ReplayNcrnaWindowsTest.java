package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import dna.AminoAcid;
import dna.Data;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CaptureNcrnaWindows;
import prok.GradeNcrnaReference;
import prok.NcrnaCallerProcess;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaJointPolicy;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec;
import prok.ReplayNcrnaWindows;
import shared.Shared;
import structures.ByteBuilder;

/** Captures the real caller once, then checks exact acceptance replay with new
 * saved flanks and legacy candidate-only data. Fresh native JVMs isolate caller
 * statics; both strands, a physical edge and a zero-window input are retained.
 * Adapted from the output-only parser in A_SampleBasic2 and the existing native
 * CaptureNcrnaWindowsTest process fixture. This is correctness, not a benchmark.
 * @author Brian Bushnell, Raiden
 */
public final class ReplayNcrnaWindowsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReplayNcrnaWindowsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown saved-window integration fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory without input or overwrite");
			new ReplayNcrnaWindowsTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_WINDOW_ACCEPTANCE_REPLAY_DRIVER_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private ReplayNcrnaWindowsTest(final File root_){
		root=root_;
		check(root.mkdir(), "The integration fixture needs a fresh output directory: "+root);
	}

	private void run(){
		resources();
		write("capture.driver.config", new ByteBuilder("caller=").append(path("caller.config"))
			.append("\nreference=").append(path("reference.tsv")).append("\nout=").append(path("capture"))
			.append("\nexpectedmodels=3\nresults=t\nreplaycontext=t\nt=2\n"));
		launch("capture_ncrna_windows.sh", "capture.driver.config");
		captureCounts();
		replay("capture", "with_context");
		checkReplay("with_context", 0);
		for(final String suffix:CAPTURE_SUFFIXES){link("capture"+suffix, "legacy"+suffix);}
		replay("legacy", "from_source");
		checkReplay("from_source", 1);
		for(final String cut:CUTS){
			sameRows("with_context."+cut+".results.tsv", "from_source."+cut+".results.tsv");
			sameRows("with_context."+cut+".inputs.tsv", "from_source."+cut+".inputs.tsv");
		}
		checkControlFailure();
		checkInputFailures();
		checkShuffle();
	}

	/** Same selected windows/permutation must retain all twelve actual acceptance
	 * outcomes and raw maxima across cached and uncached execution.
	 */
	private void checkShuffle(){
		for(final String cache:new String[]{"t", "f"}){
			final String name="shuffle_"+cache;
			write(name+".driver.config", new ByteBuilder("in=").append(path("capture"))
				.append("\ncaller=").append(path("caller.config")).append("\nout=").append(path(name))
				.append("\nexpectedinputs=2\nexpectedmodels=3\nexpectedwindows=").append(windows)
				.append("\nexpectedscheduled=").append(scheduled).append("\nshufflewindows=").append(windows)
				.append("\nshuffleseed=1701\nshufflecache=").append(cache).append("\nt=2\n"));
			launch("replay_ncrna_windows.sh", name+".driver.config");
		}
		try(Input cached=new Input(path("shuffle_t.shuffle.windows.tsv")); Input uncached=new Input(path("shuffle_f.shuffle.windows.tsv"))){
			cached.reader.nextLine();
			uncached.reader.nextLine();
			final LineParser1 a=new LineParser1('\t'), b=new LineParser1('\t');
			int count=0;
			for(byte[] line; (line=cached.reader.nextLine())!=null;){
				a.set(line);
				final byte[] other=uncached.reader.nextLine();
				check(other!=null, "Uncached shuffle cannot omit a sampled event");
				b.set(other);
				for(int field=0; field<5; field++){
					check(a.parseString(field).equals(b.parseString(field)), "Cached shuffle changed a sampled identity, raw maximum or actual acceptance mask");
				}
				check(a.parseLong(5)==b.parseLong(5) && a.parseLong(7)<=b.parseLong(7),
					"Exact reuse must preserve logical work and cannot increase physically executed kernels");
				count++;
			}
			check(count==windows && uncached.reader.nextLine()==null, "Every selected window must retain exactly one null-control row");
		}
	}

	/** Three small explicit models keep this fixture independent of library size.
 * Four spaced distinct17-mers make a clipped candidate with real outside flanks. */
	private void resources(){
		final byte[] model=Prok16sPartialRescueTest.random(320, 8173);
		final byte[] variant=model.clone();
		variant[170]=variant[170]=='A' ? (byte)'C' : (byte)'A';
		final byte[][] bases={model, variant, Prok16sPartialRescueTest.random(320, 1927)};
		final String[] names={"cobalamin_RF00174_fixture_0", "cobalamin_RF00174_fixture_1", "cobalamin_RF00174_fixture_2"};
		final ArrayList<BaseGraph> graphs=new ArrayList<BaseGraph>();
		final ByteBuilder fasta=new ByteBuilder();
		for(int i=0; i<bases.length; i++){
			fasta.append('>').append(names[i]).nl().append(bases[i]).nl();
			graphs.add(new BaseGraph(names[i], bases[i], null, i, 0));
		}
		write("consensus.fa", fasta);
		NucleotideHbmCodec.writeModels(graphs, path("models.hbmc.gz"), false);
		NcrnaDdlProfileResource.writeSelected(path("profiles.tsv"), new char[3][1024], new int[]{1, 1, 1},
			names, bases, 1024, new LongHashSet(4));
		final LongHashSet keys=new LongHashSet(8);
		final ByteBuilder seeds=new ByteBuilder();
		for(final int start:new int[]{150, 151, 230, 231}){
			long key=0;
			for(int i=start; i<start+17; i++){key=(key<<2)|AminoAcid.baseToNumber[model[i]];}
			check(keys.add(key), "Fixture discovery keys must be distinct before a two-key gate can be tested");
			seeds.append('>').append(start).nl().append(model, start, 17).nl();
		}
		write("seeds.fa", seeds);
		write("policy.tsv", new ByteBuilder(NcrnaJointPolicy.HEADER).nl()
			.append("all\tcobalamin_rf00174\tadaptive\t2147483647\t2\t0\t0\t0\t2147483647\n"));
		final byte[] full=Prok16sPartialRescueTest.embed(model, 400, 400), minus=full.clone();
		AminoAcid.reverseComplementBasesInPlace(minus);
		final byte[] edge=Prok16sPartialRescueTest.embed(Arrays.copyOfRange(model, 160, model.length), 0, 80);
		write("genome.fa", new ByteBuilder(">plus description\n").append(full).nl()
			.append(">minus description\n").append(minus).nl().append(">physical_left\n").append(edge).nl());
		write("zero.fa", new ByteBuilder(">plus unrelated input with reused name\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"));
		write("reference.tsv", new ByteBuilder(GradeNcrnaReference.MANIFEST_HEADER).nl()
			.append("real\t").append(path("genome.fa")).tab().append(path("unused.gff")).nl()
			.append("zero\t").append(path("zero.fa")).tab().append(path("unused_zero.gff")).nl());
		write("caller.config", new ByteBuilder("cdsonly=t\ncallcds=f\nncrna=f\ncobalamin=t\ntmrna=f\nsixs=f\nfmn=f\n")
			.append("sam1=f\ntpp=f\nglycine=f\npasses=1\ntaxonomy=f\nncrnafamily=cobalamin_rf00174\n")
			.append("ncrnaseedminhits=2\ncobalaminminlen=109\nncrnapad=150\n")
			.append("ncrnaconsensus=").append(path("consensus.fa")).append("\nncrnamodels=").append(path("models.hbmc.gz"))
			.append("\nncrnakmers=").append(path("seeds.fa")).append("\ncobalaminddl=").append(path("profiles.tsv"))
			.append("\ncobalaminddlcap=all\nncrnapolicy=").append(path("policy.tsv"))
			.append("\nin=").append(path("genome.fa")).append(',').append(path("zero.fa")).append("\nt=2\n"));
	}

	private void captureCounts(){
		int inputs=0;
		try(Input in=new Input(path("capture.inputs.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==10 && p.parseInt(0)==inputs, "The actual capture must retain original input order and dimensions");
				if(inputs==1){
					check(p.parseLong(3)==1 && p.parseLong(4)==32 && p.parseInt(5)==0 && p.parseLong(6)==0,
						"The reused-name zero-window original must survive with real zero work");
				}
				windows+=p.parseInt(5);
				scheduled+=p.parseLong(6);
				accepted+=p.parseLong(8);
				inputs++;
			}
		}
		check(inputs==2 && windows>=3 && scheduled>=windows && accepted>=2,
			"The actual fixture must exercise both full-model strands, a physical edge and accepted proposals");
	}

	private void replay(final String capture, final String output){
		replay(capture, output, "caller.config");
	}
	private void replay(final String capture, final String output, final String caller){
		write(output+".driver.config", new ByteBuilder("in=").append(path(capture))
			.append("\ncaller=").append(path(caller)).append("\nout=").append(path(output))
			.append("\nexpectedinputs=2\nexpectedmodels=3\nexpectedwindows=").append(windows)
			.append("\nexpectedscheduled=").append(scheduled).append("\nt=2\n"));
		launch("replay_ncrna_windows.sh", output+".driver.config");
	}

	/** The exact resource settings do not authorize a different original-input
	 * universe. Reject swapped and extra paths before any window is evaluated. */
	private void checkInputFailures(){
		write("extra.fa", new ByteBuilder(">extra\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"));
		for(int mode=0; mode<2; mode++){
			final String caller="wrong_inputs_"+mode+".config", output="wrong_inputs_"+mode;
			try(Input in=new Input(path("caller.config")); Output out=new Output(path(caller))){
				for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
					if(line.length>=3 && line[0]=='i' && line[1]=='n' && line[2]=='='){
						row.clear().append("in=");
						if(mode==0){row.append(path("zero.fa")).append(',').append(path("genome.fa"));}
						else{row.append(path("genome.fa")).append(',').append(path("zero.fa")).append(',').append(path("extra.fa"));}
					}else{row.clear().append(line);}
					out.writer.print(row.nl());
				}
			}
			boolean rejected=false;
			try{replay("capture", output, caller);}
			catch(RuntimeException expected){rejected=true;}
			check(rejected && !new File(path(output+".064.results.tsv")).exists(),
				"Swapped or extra configured original inputs must fail before window observations are written");
		}
	}

	private void checkReplay(final String prefix, final int sourceLoads){
		for(int cut=0; cut<CUTS.length; cut++){
			final String name=prefix+"."+CUTS[cut];
			metric(name, "identity", new float[]{.60f, .62f, .64f}[cut]);
			metric(name, "inputs", 2);
			metric(name, "captured_windows", windows);
			metric(name, "scheduled_windows", scheduled);
			metric(name, "shared_source_file_loads_NONADDITIVE_ACROSS_CUTS", sourceLoads);
			metric(name, "control_windows_checked", cut==2 ? windows : 0);
			check(!new File(path(name+".calls.gff")).exists(), "Acceptance replay must not fabricate final path output");
			try(Input in=new Input(path(name+".results.tsv"))){
				check(Arrays.equals(in.reader.nextLine(), ReplayNcrnaWindows.RESULT_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
					"Replay results must declare actual total work without invented alignment phases");
				long rows=0, accepts=0;
				for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
					p.set(line);
					check(p.terms()==23, "Replay rows must retain the declared23-field acceptance schema");
					if(p.termEquals("ACCEPTED_PROPOSAL", 13)){accepts++;}
					rows++;
				}
				check(rows==scheduled && (cut!=2 || accepts==accepted), "The unchanged control must conserve all actual accepted proposals and original scheduled events");
			}
		}
	}

	/** Alter one saved work count while retaining the exact captured bytes. A
 * false control must fail, even if the altered aggregate still looks plausible. */
	private void checkControlFailure(){
		for(final String suffix:new String[]{".windows.tsv", ".windows.fa.gz", ".inputs.tsv"}){link("capture"+suffix, "wrong"+suffix);}
		try(Input in=new Input(path("capture.results.tsv")); Output out=new Output(path("wrong.results.tsv"))){
			out.writer.print(row.clear().append(in.reader.nextLine()).nl());
			boolean changed=false;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				row.clear();
				for(int column=0; column<26; column++){
					if(column>0){row.tab();}
					if(!changed && column==17){row.append(p.parseLong(17)+1);}
					else{p.setBounds(column);row.append(line, p.a(), p.b()-p.a());}
				}
				changed=true;
				out.writer.print(row.nl());
			}
			check(changed, "The negative control must alter one actual observed event");
		}
		boolean rejected=false;
		try{replay("wrong", "must_fail");}
		catch(RuntimeException expected){rejected=true;}
		check(rejected && !new File(path("must_fail.064.replay.tsv")).exists(),
			"A changed64% work observation must fail before any completion summary");
	}

	private void launch(final String script, final String config){
		final String launcher=Data.findPath("?dev/"+script, false);
		check(launcher!=null, "The finite process fixture must use its installed native launcher: "+script);
		NcrnaCallerProcess.run(new File(launcher), new File(path(config)), 1024);
	}
	private void metric(final String name, final String key, final double expected){
		try(Input in=new Input(path(name+".replay.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals(key, 0)){
					check(Math.abs(p.parseDouble(1)-expected)<.000001, "Actual replay metric changed: "+key);
					return;
				}
			}
		}
		throw new AssertionError("Missing actual replay metric: "+key);
	}
	private void sameRows(final String first, final String second){
		try(Input a=new Input(path(first)); Input b=new Input(path(second))){
			for(byte[] line=a.reader.nextLine(); line!=null; line=a.reader.nextLine()){
				check(Arrays.equals(line, b.reader.nextLine()), "Saved flanks and exact source recovery must produce identical acceptance rows");
			}
			check(b.reader.nextLine()==null, "The alternate context route must not append events");
		}
	}
	private void link(final String source, final String target){
		try{Files.createSymbolicLink(new File(path(target)).toPath(), new File(path(source)).toPath());}
		catch(IOException failure){throw new IllegalStateException("Cannot bind the same finite capture data without copying a program tree", failure);}
	}
	private void write(final String name, final ByteBuilder bytes){try(Output out=new Output(path(name))){out.writer.print(bytes);}}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder row=new ByteBuilder();
	private int windows;
	private long scheduled, accepted;
	private static final String[] CUTS={"060", "062", "064"};
	private static final String[] CAPTURE_SUFFIXES={".windows.tsv", ".windows.fa.gz", ".results.tsv", ".inputs.tsv"};
}
