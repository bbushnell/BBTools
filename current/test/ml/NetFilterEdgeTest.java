package test.ml;

import ml.CellNet;
import ml.NetFilter;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Stream;

import fileIO.FileFormat;
import shared.Timer;

/**
 * Literal edge fixtures for {@link NetFilter} constructor and process edges:
 * frame score aggregation, paired pass/unmatched output expansion, file-routing
 * collision checks, overwrite guards and optional histogram output.
 * Run with {@code -ea}; most fixture oracles are Java assertions.
 * @author Brian, Jean
 */
public class NetFilterEdgeTest{

	private static final int AVERAGE=1;
	private static final int MAX=2;
	private static final int MIN=3;
	private static final byte[] SCORE_BASES=("ACGTACGTACGTACGTACGT"+
		"ACGTACGTACGTACGTACGT").getBytes(StandardCharsets.US_ASCII);
	private static final int SCORE_WIDTH=20;
	private static final int SCORE_STEP=8;
	private static final byte[] FIRST_SPAN_BASES="AAAAAAAAAAAAAAAAAAAAC"
		.getBytes(StandardCharsets.US_ASCII);
	private static final byte[] EXACT_END_BASES=("AAAAAAAAAAAAAAAAAAAA"+
		"CCCCCCCCCCCCCCCC").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] TAIL_BASES=("AAAAAAAAAAAAAAAAAAAA"+
		"CCCCCCCCCCCCCCCCGGGGGGGGG").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] SHORT_BASES="ACGT".getBytes(StandardCharsets.US_ASCII);
	private static final String FIRST_EXACT_WINDOW="AAAAAAAAAAAAAAAAAAAA";
	private static final String SECOND_EXACT_WINDOW="AAAACCCCCCCCCCCCCCCC";
	private static final String HIST_SEQUENCE="ACGTACGTACGTACGT";
	private static final String HIST_QUALITY="IIIIIIIIIIIIIIII";
	private static final String HIST_FASTQ="@neg result=0\n"+HIST_SEQUENCE+"\n+\n"+
		HIST_QUALITY+"\n@pos result=1\n"+HIST_SEQUENCE+"\n+\n"+
		HIST_QUALITY+"\n";
	
	/**
	 * Runs every fixture case or, when one argument is supplied, a single named
	 * case so baseline and repaired behavior can be compared row-by-row.
	 */
	public static void main(String[] args) throws Exception{
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts only one optional case name";
		run("score_average");
		run("score_min");
		run("score_max");
		run("score_exact_end_average");
		run("score_exact_end_min");
		run("score_exact_end_max");
		run("score_tail_omitted");
		run("score_short_clamped");
		run("score_first_span");
		run("score_explicit_span_control");
		run("unmatched_hash");
		run("main_hash");
		run("explicit_unmatched");
		run("unmatched_input_collision");
		run("unmatched_output_collision");
		run("unmatched_pair_collision");
		run("unmatched_overwrite_reject");
		run("unmatched_overwrite_control");
		run("unmatched_collision_control");
		run("overlap_zero");
		run("explicit_step");
		run("step_zero");
		run("step_negative");
		run("overlap_width");
		run("hist_write");
		run("hist_absent_control");
		run("hist_input_collision");
		run("hist_output_collision");
		run("hist_overwrite_reject");
	}
	
	/** Runs one selectable case in an isolated temporary directory. */
	private static void run(String name) throws Exception{
		final Path dir=Files.createTempDirectory("netfilter-edge-");
		try{
			if(name.equals("score_average")){
				final StubNet net=new StubNet(3, 6, 9);
				final float score=NetFilter.scoreFrames(SCORE_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, AVERAGE, SCORE_STEP, SCORE_WIDTH);
				assert(score==6f) : "Average must divide by the 3 scored frames, not 4: "+score;
				net.checkConsumed();
			}else if(name.equals("score_min")){
				final StubNet net=new StubNet(3, 6, 9);
				final float score=NetFilter.scoreFrames(SCORE_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, MIN, SCORE_STEP, SCORE_WIDTH);
				assert(score==3f) : "Minimum frame score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_max")){
				final StubNet net=new StubNet(3, 6, 9);
				final float score=NetFilter.scoreFrames(SCORE_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, MAX, SCORE_STEP, SCORE_WIDTH);
				assert(score==9f) : "Maximum frame score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_exact_end_average")){
				final WindowNet net=new WindowNet(
					new String[]{FIRST_EXACT_WINDOW, SECOND_EXACT_WINDOW}, 2, 8);
				final float score=NetFilter.scoreFrames(EXACT_END_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, AVERAGE, 16, SCORE_WIDTH);
				assert(score==5f) : "Exact-end average score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_exact_end_min")){
				final WindowNet net=new WindowNet(
					new String[]{FIRST_EXACT_WINDOW, SECOND_EXACT_WINDOW}, 2, 8);
				final float score=NetFilter.scoreFrames(EXACT_END_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, MIN, 16, SCORE_WIDTH);
				assert(score==2f) : "Exact-end minimum score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_exact_end_max")){
				final WindowNet net=new WindowNet(
					new String[]{FIRST_EXACT_WINDOW, SECOND_EXACT_WINDOW}, 2, 8);
				final float score=NetFilter.scoreFrames(EXACT_END_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, MAX, 16, SCORE_WIDTH);
				assert(score==8f) : "Exact-end maximum score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_tail_omitted")){
				final WindowNet net=new WindowNet(
					new String[]{FIRST_EXACT_WINDOW, SECOND_EXACT_WINDOW}, 4, 10);
				final float score=NetFilter.scoreFrames(TAIL_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, AVERAGE, 16, SCORE_WIDTH);
				assert(score==7f) : "Incomplete tail changed the average: "+score;
				net.checkConsumed();
			}else if(name.equals("score_short_clamped")){
				final WindowNet net=new WindowNet(new String[]{"ACGT"}, 11);
				final float score=NetFilter.scoreFrames(SHORT_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, AVERAGE, 16, SCORE_WIDTH);
				assert(score==11f) : "Short-read score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_first_span")){
				final RecordingNet net=new RecordingNet();
				final float score=NetFilter.scoreFrames(FIRST_SPAN_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, AVERAGE, 99, SCORE_WIDTH);
				assert(score==1f) : "First frame score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("score_explicit_span_control")){
				final RecordingNet net=new RecordingNet();
				final float score=NetFilter.score(FIRST_SPAN_BASES,
					new float[SCORE_WIDTH*4+4], 0, net, 0, SCORE_WIDTH-1);
				assert(score==1f) : "Explicit frame score drifted: "+score;
				net.checkConsumed();
			}else if(name.equals("unmatched_hash")){
				final NetFilter nf=filter(dir, "pass#.fq", "fail#.fq", null);
				check(names(nf, "ffoutu1", "ffoutu2"), "fail1.fq", "fail2.fq");
			}else if(name.equals("main_hash")){
				final NetFilter nf=filter(dir, "pass#.fq", "fail1.fq", "fail2.fq");
				check(names(nf, "ffout1", "ffout2"), "pass1.fq", "pass2.fq");
				check(names(nf, "ffoutu1", "ffoutu2"), "fail1.fq", "fail2.fq");
			}else if(name.equals("explicit_unmatched")){
				final NetFilter nf=filter(dir, "pass.fq", "failA.fq", "failB.fq");
				check(names(nf, "ffoutu1", "ffoutu2"), "failA.fq", "failB.fq");
			}else if(name.equals("unmatched_input_collision")){
				checkRejects(dir, "pass.fq", "in1.fq", "failB.fq",
					"was specified multiple times.");
			}else if(name.equals("unmatched_output_collision")){
				checkRejects(dir, "pass.fq", "pass.fq", "failB.fq",
					"was specified for multiple output streams.");
			}else if(name.equals("unmatched_pair_collision")){
				checkRejects(dir, "pass.fq", "failA.fq", "failA.fq",
					"was specified for multiple output streams.");
			}else if(name.equals("unmatched_overwrite_reject")){
				writeFastq(dir, "failA.fq");
				checkRejects(dir, "pass.fq", "failA.fq", "failB.fq",
					"exists and overwrite=false", "ow=f");
			}else if(name.equals("unmatched_overwrite_control")){
				writeFastq(dir, "failA.fq");
				filter(dir, "pass.fq", "failA.fq", "failB.fq");
			}else if(name.equals("unmatched_collision_control")){
				filter(dir, "pass.fq", "failA.fq", "failB.fq");
			}else if(name.equals("overlap_zero")){
				filter(dir, "pass.fq", "failA.fq", "failB.fq", "overlap=0");
			}else if(name.equals("explicit_step")){
				filter(dir, "pass.fq", "failA.fq", "failB.fq", "step=1");
			}else if(name.equals("step_zero")){
				checkRejects(dir, "step=0");
			}else if(name.equals("step_negative")){
				checkRejects(dir, "step=-1");
			}else if(name.equals("overlap_width")){
				checkRejects(dir, "overlap=1");
			}else if(name.equals("hist_write")){
				checkHistogram(dir, true);
			}else if(name.equals("hist_absent_control")){
				checkHistogram(dir, false);
			}else if(name.equals("hist_input_collision")){
				checkRejectsHistogram(dir, "was specified multiple times.",
					"hist="+dir.resolve("hist.fq"));
			}else if(name.equals("hist_output_collision")){
				checkRejectsHistogram(dir, "was specified for multiple output streams.",
					"hist="+dir.resolve("pass.fq"));
			}else if(name.equals("hist_overwrite_reject")){
				Files.write(dir.resolve("histogram.txt"), new byte[0]);
				checkRejectsHistogram(dir, "exists and overwrite=false",
					"hist="+dir.resolve("histogram.txt"), "ow=f");
			}else{
				throw new IllegalArgumentException("Unknown fixture mode: "+name);
			}
			System.err.println("NetFilterEdgeTest PASS "+name);
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}
	
	/** Constructs a paired-input filter with pass/unmatched outputs and optional extra flags. */
	private static NetFilter filter(Path dir, String out, String outu1, String outu2,
			String... extra) throws Exception{
		final Path net=writeNet(dir);
		final Path in1=writeFastq(dir, "in1.fq");
		final Path in2=writeFastq(dir, "in2.fq");
		final String[] args=new String[8+extra.length];
		args[0]="in="+in1;
		args[1]="in2="+in2;
		args[2]="net="+net;
		args[3]="out="+dir.resolve(out);
		args[4]="outu="+dir.resolve(outu1);
		args[5]=outu2==null ? "ordered=t" : "outu2="+dir.resolve(outu2);
		args[6]="ow=t";
		args[7]="width=1";
		System.arraycopy(extra, 0, args, 8, extra.length);
		return new NetFilter(args);
	}

	/** Runs a parsed two-read filter and checks ordinary output plus optional hist= emission. */
	private static void checkHistogram(Path dir, boolean writeHist) throws Exception{
		final Path net=writeNet(dir);
		final Path in=writeFastq(dir, "hist.fq", "@neg result=0", "@pos result=1");
		final Path hist=dir.resolve("histogram.txt");
		final Path out=dir.resolve("pass.fq");
		final String[] args=writeHist ?
			new String[]{"in="+in, "net="+net, "out="+out, "hist="+hist,
				"width=1", "parse=t", "filter=f", "t=1", "ow=t"} :
			new String[]{"in="+in, "net="+net, "out="+out,
				"width=1", "parse=t", "filter=f", "t=1", "ow=t"};
		new NetFilter(args).process(new Timer());
		checkHistogramOutput(out);
		if(writeHist){
			assert(Files.exists(hist)) : "Expected hist output at "+hist;
			final String actual=new String(Files.readAllBytes(hist), StandardCharsets.US_ASCII);
			final String expected="#count\tneg\tpos\n0\t1\t1\n";
			assert(actual.equals(expected)) : "Expected "+expected+", got "+actual;
		}else{
			assert(!Files.exists(hist)) : "hist= was absent but "+hist+" was created";
		}
	}

	/** Confirms invalid hist= paths still stop before processing. */
	private static void checkRejectsHistogram(Path dir, String message,
			String... extra) throws Exception{
		final Path net=writeNet(dir);
		final Path in=writeFastq(dir, "hist.fq", "@neg result=0", "@pos result=1");
		final String[] args=new String[5+extra.length];
		args[0]="in="+in;
		args[1]="net="+net;
		args[2]="out="+dir.resolve("pass.fq");
		args[3]="width=1";
		args[4]="parse=t";
		System.arraycopy(extra, 0, args, 5, extra.length);
		try{
			new NetFilter(args);
		}catch(RuntimeException|AssertionError expected){
			final String actual=String.valueOf(expected.getMessage());
			assert(actual.contains(message)) : "Expected rejection containing "+message+
				", got "+actual;
			return;
		}
		throw new AssertionError("Expected constructor rejection for "+
			Arrays.toString(extra));
	}

	/** Confirms invalid step geometry still stops in the constructor. */
	private static void checkRejects(Path dir, String extra) throws Exception{
		try{
			filter(dir, "pass.fq", "failA.fq", "failB.fq", extra);
		}catch(AssertionError expected){
			return;
		}
		throw new AssertionError("Expected constructor rejection for "+extra);
	}
	
	/** Confirms invalid file routing still stops during construction. */
	private static void checkRejects(Path dir, String out, String outu1, String outu2,
			String message, String... extra) throws Exception{
		try{
			filter(dir, out, outu1, outu2, extra);
		}catch(RuntimeException|AssertionError expected){
			final String actual=String.valueOf(expected.getMessage());
			assert(actual.contains(message)) : "Expected rejection containing "+message+
				", got "+actual;
			return;
		}
		throw new AssertionError("Expected constructor rejection for out="+out+
			", outu="+outu1+", outu2="+outu2);
	}

	/** Writes a minimal dense fixture network returning deterministic zero scores in NetFilter. */
	private static Path writeNet(Path dir) throws Exception{
		final Path path=dir.resolve("fixture.bbnet");
		Files.write(path, Arrays.asList(
			"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
			"#seed 1", "#layers 2", "#dims 8 1", "#edges 0", "#coding decimal",
			"##layer 1", "C9 LINEAR 0 0 0 0 0 0 0 0 0"
		), StandardCharsets.US_ASCII);
		return path;
	}
	
	/** Writes one valid four-base read so paired FileFormat probes can open FASTQ input. */
	private static Path writeFastq(Path dir, String name) throws Exception{
		final Path path=dir.resolve(name);
		Files.write(path, Arrays.asList("@r0", "ACGT", "+", "IIII"), StandardCharsets.US_ASCII);
		return path;
	}
	
	/** Writes named reads long enough to keep the histogram fixture away from short-span scoring. */
	private static Path writeFastq(Path dir, String name, String... ids) throws Exception{
		final Path path=dir.resolve(name);
		final String[] lines=new String[ids.length*4];
		for(int i=0, j=0; i<ids.length; i++){
			lines[j++]=ids[i];
			lines[j++]=HIST_SEQUENCE;
			lines[j++]="+";
			lines[j++]=HIST_QUALITY;
		}
		Files.write(path, Arrays.asList(lines), StandardCharsets.US_ASCII);
		return path;
	}

	/** Confirms the histogram fixture writes both parsed input reads unchanged. */
	private static void checkHistogramOutput(Path out) throws Exception{
		assert(Files.exists(out)) : "Expected FASTQ output at "+out;
		final String actual=new String(Files.readAllBytes(out), StandardCharsets.US_ASCII);
		assert(actual.equals(HIST_FASTQ)) : "Expected "+HIST_FASTQ+", got "+actual;
	}

	/** Extracts the normalized output names created by NetFilter. */
	private static String[] names(NetFilter nf, String a, String b) throws Exception{
		return new String[]{name(nf, a), name(nf, b)};
	}
	
	/** Extracts one FileFormat name from a private final field. */
	private static String name(NetFilter nf, String fieldName) throws Exception{
		final Field field=NetFilter.class.getDeclaredField(fieldName);
		field.setAccessible(true);
		final FileFormat ff=(FileFormat)field.get(nf);
		return Paths.get(ff.name()).getFileName().toString();
	}
	
	/** Asserts an exact two-name expansion result. */
	private static void check(String[] actual, String a, String b){
		assert(actual[0].equals(a) && actual[1].equals(b)) :
			"Expected ["+a+", "+b+"], got "+Arrays.toString(actual);
	}
	
	/** Asserts exact raw frame windows and feeds fixed scores. */
	private static final class WindowNet extends CellNet{
		WindowNet(String[] expected, float... scores){
			assert(expected.length==scores.length) :
				"Expected windows "+expected.length+" != scores "+scores.length;
			this.expected=expected;
			this.scores=scores;
		}

		@Override
		public CellNet applyInput(float[] valuesIn){
			assert(index<expected.length) : "Too many frames scored";
			final String window=expected[index];
			assert(valuesIn[0]==window.length()/(float)(SCORE_WIDTH+5)) :
				"Frame "+index+" length should be "+window.length()+", got "+valuesIn[0];
			for(int i=0; i<SCORE_WIDTH; i++){
				final int offset=4+i*4;
				final char base=i<window.length() ? window.charAt(i) : 0;
				checkBase(valuesIn, offset, base, index, i);
			}
			return this;
		}

		@Override
		public float feedForward(){
			assert(index<scores.length) : "Too many scores consumed";
			return scores[index++];
		}

		void checkConsumed(){
			assert(index==scores.length) : "Only scored "+index+" of "+scores.length+" frames";
		}

		/** Checks one raw one-hot base quartet at vector offset. */
		private static void checkBase(float[] valuesIn, int offset, char base,
				int frame, int pos){
			final int lane=(base=='A' ? 0 : base=='C' ? 1 : base=='G' ? 2 :
				base=='T' ? 3 : -1);
			assert(lane>=0 || base==0) : "Unexpected expected base "+base;
			for(int i=0; i<4; i++){
				final float expectedValue=(i==lane ? 1f : 0f);
				assert(valuesIn[offset+i]==expectedValue) :
					"Frame "+frame+" base "+pos+" lane "+i+" expected "+
					expectedValue+", got "+valuesIn[offset+i];
			}
		}

		private final String[] expected;
		private final float[] scores;
		private int index;
	}

	/** Asserts one 20-A vector and returns one fixed score. */
	private static final class RecordingNet extends CellNet{
		@Override
		public CellNet applyInput(float[] valuesIn){
			assert(!seen) : "Too many frames scored";
			seen=true;
			assert(valuesIn[0]==SCORE_WIDTH/(float)(SCORE_WIDTH+5)) :
				"First frame length should be exactly 20 bases: "+valuesIn[0];
			assert(valuesIn[1]==0f) : "First frame should have no GC bases: "+valuesIn[1];
			for(int i=0; i<SCORE_WIDTH; i++){
				final int offset=4+i*4;
				assert(valuesIn[offset]==1f) : "Expected A lane at base "+i;
				assert(valuesIn[offset+1]==0f && valuesIn[offset+2]==0f && valuesIn[offset+3]==0f) :
					"Unexpected non-A lane at base "+i;
			}
			return this;
		}

		@Override
		public float feedForward(){return 1f;}

		void checkConsumed(){
			assert(seen) : "No frame scored";
		}

		private boolean seen;
	}

	/** Feeds a fixed score sequence to scoreFrames without depending on trained weights. */
	private static final class StubNet extends CellNet{
		StubNet(float... scores){
			this.scores=scores;
		}
		
		@Override
		public CellNet applyInput(float[] valuesIn){return this;}
		
		@Override
		public float feedForward(){
			assert(index<scores.length) : "Too many frames scored";
			return scores[index++];
		}
		
		void checkConsumed(){
			assert(index==scores.length) : "Only scored "+index+" of "+scores.length+" frames";
		}
		
		private final float[] scores;
		private int index;
	}
}
