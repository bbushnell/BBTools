package test.prok;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.PreParser;
import prok.NcrnaBoundaryVectorGen;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Native optional-locus output ownership, including the two outer vector streams.
 * Uses forty saved synthetic contexts; literal row counts are engineering oracles,
 * not biological coverage or fitted production settings. Observations precede cleanup.
 * @author Raiden */
public final class NcrnaVectorOutputTest{
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaVectorOutputTest.class, false);
		String root=null, data=null, tables=null, expect=null;
		for(String arg:pp.args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("data=")){data=arg.substring(5);}
			else if(arg.startsWith("tables=")){tables=arg.substring(7);}
			else if(arg.startsWith("expect=")){expect=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown vector output fixture argument: "+arg);}
		}
		check(root!=null && new File(root).isAbsolute() && data!=null && tables!=null
			&& ("before".equals(expect) || "after".equals(expect)), "Declare absolute fresh output, fixed inputs and before/after expectation");
		check(new File(root).mkdir(), "Never replace fixture evidence");
		new NcrnaVectorOutputTest(root, data, tables, expect.equals("after")).run();
	}
	private NcrnaVectorOutputTest(String r, String d, String t, boolean a){
		root=r;
		data=d;
		tables=t;
		after=a;
	}
	private void run() throws Exception{
		shared.Shared.setThreads(1);fileIO.ReadWrite.setZipThreads(1);
		final ByteFile input=ByteFile.makeByteFile(data+"/contexts.fa", true);int changed=0;
		try{
			try(Output bad=new Output(FileFormat.testOutput(root+"/bad.fa", FileFormat.TEXT, null, true, false, false, false))){
				final ByteBuilder buffer=new ByteBuilder();
				for(byte[] row; (row=input.nextLine())!=null;){
					if(changed==0 && row.length>0 && row[0]=='>'){
						final String header=new String(row, StandardCharsets.UTF_8);
						check(header.contains("contig_gc=0.5"), "Saved synthetic header supplies the exact intended GC field");
						row=header.replace("contig_gc=0.5", "contig_gc=NaN").getBytes(StandardCharsets.UTF_8);changed++;
					}
					bad.writer.print(buffer.clear().append(row).nl());
				}
			}
		}finally{check(!input.close(), "Read the complete fixed corpus before fault controls");}
		check(changed==1, "Exactly one source header carries the deliberate body failure");
		final ByteBuilder observations=new ByteBuilder("case\tlocus_alive\tlocus_descriptors\tvector_descriptors\ttrainer_descriptors\tprimary\tsuppressed\n");
		for(String scenario:CASES){one(scenario, observations);}
		try(Output out=new Output(FileFormat.testOutput(root+"/observations.tsv", FileFormat.TEXT, null, true, false, false, false))){
			out.writer.print(observations);
		}
		// These finite outputs contain 481 lines each; compare their bytes without a hash subprocess.
		final byte[] expected=ReadWrite.readRaw(root+"/ordinary.vectors");
		check(Arrays.equals(expected, ReadWrite.readRaw(root+"/ordinary.trainer")), "Both normal sinks receive the identical vector rows");
		for(String suffix:new String[]{"vectors", "trainer"}){
			check(Arrays.equals(expected, ReadWrite.readRaw(root+"/omitted."+suffix)), "Omitting optional locus output cannot change vectors");
			check(lines(root+"/ordinary."+suffix)==481, "Forty contexts times twelve candidates plus one header");
		}
		final ByteFile locus=ByteFile.makeByteFile(root+"/ordinary.locus", true);final LineParser1 p=new LineParser1('\t');int n=0;
		try{
			check(new String(locus.nextLine(), StandardCharsets.UTF_8).equals("locus\tassembly\tmodel\tstart_rows\tstop_rows"), "Locus header contract");
			for(byte[] row; (row=locus.nextLine())!=null; n++){
				p.set(row);check(p.terms()==5 && p.termEquals("g"+n, 0) && p.termEquals("s"+n, 1) && p.termEquals("SAM1_fixture1", 2)
					&& p.parseInt(3)==5 && p.parseInt(4)==7, "Literal locus IDs, winning model and per-end candidate counts must survive");
			}
		}finally{check(!locus.close(), "Locus oracle read failed");}
		check(n==40 && !new File(root+"/omitted.locus").exists(), "All forty loci or explicitly no locus output");
		System.out.println("NCRNA_VECTOR_OUTPUT_PASS cases="+CASES.length+" expectation="+(after ? "after" : "before"));
	}
	private void one(String scenario, ByteBuilder observations) throws Exception{
		final String prefix=root+"/"+scenario;final Factory factory=new Factory(scenario);Throwable caught=null;
		boolean alive=false;int locusDescriptors=0, vectorDescriptors=0, trainerDescriptors=0;
		final ArrayList<String> args=new ArrayList<String>();
		java.util.Collections.addAll(args, "family=external", "externalname=fixture", "externalconsensus="+data+"/models.fa",
			"externalmodels="+data+"/models.hbmc", "externalseeds="+data+"/seeds.fa", "externalk=17", "externalseedminhits=1",
			"externalminlen=1", "externalwindowpad=20", "externalindexk=7", "externalfixedminhits=1", "externaltopn=3",
			"externalidpass=0.6", "externalidborderline=0.55", "externaladaptive=f", "externalrankedfallback=f", "externalstrictindexcutoff=f",
			"externalstartoffsets=-2,-1,0,1,2", "externalstopoffsets=-3,-2,-1,0,1,2,3", "meanlen=120", "leaveoneout=f", "boundaryfeatures=v1",
			"tablestart="+tables+"/tuning.start.tsv.gz", "tablestop="+tables+"/tuning.stop.tsv.gz",
			"fasta="+(scenario.startsWith("bad_input") ? root+"/bad.fa" : data+"/contexts.fa"), "out="+prefix+".vectors", "trainerout="+prefix+".trainer");
		if(!scenario.equals("omitted")){args.add("locusout="+prefix+".locus");}
		try{
			try{invoke(args.toArray(new String[0]), factory);}catch(Exception | Error failure){caught=failure;}
			alive=factory.writer!=null && factory.writer.isAlive();locusDescriptors=descriptors(prefix+".locus");
			vectorDescriptors=descriptors(prefix+".vectors");trainerDescriptors=descriptors(prefix+".trainer");
			check((caught==null)==(scenario.equals("ordinary") || scenario.equals("omitted")), "Unexpected vector invocation outcome: "+scenario+" / "+caught);
			check(!alive && vectorDescriptors==0 && trainerDescriptors==0 && locusDescriptors==(!after && scenario.equals("start") ? 1 : 0), "Observe native owners before fixture cleanup: "+scenario);
			if(scenario.equals("open") || scenario.equals("native_open") || scenario.equals("start")){check(caught==factory.primary, "Retain acquisition/start failure identity");}
			final boolean injected=scenario.startsWith("header") || scenario.startsWith("write");
			final boolean combined=scenario.contains("close") && (injected || scenario.startsWith("bad_input"));
			if(injected){
				if(after || !combined){check(caught==factory.body, "Keep original header/write failure");}
				else{check(caught!=factory.body, "Baseline finally masks the original error");}
			}
			if(scenario.startsWith("bad_input") && (after || !combined)){
				check(caught instanceof IllegalArgumentException && caught.getMessage().startsWith("Boundary GC must be finite"), "Exercise native malformed GC failure inside the owned body");
			}
			if(scenario.equals("close_exception") || (!after && (scenario.equals("header_close_exception") || scenario.equals("write_close_exception") || scenario.equals("bad_input_close")))){
				check(caught==factory.cleanup, "Retain actual cleanup exception identity");
			}
			check(caught==null || caught.getSuppressed().length==(after && combined ? 1 : 0), "Only repaired combined failures suppress cleanup");
		}finally{
			if(factory.writer!=null && !factory.writer.closed){try{if(factory.writer.getState()==Thread.State.NEW){factory.writer.closeBeforeStart();}else{factory.writer.poisonAndWait();}}catch(RuntimeException | Error cleanup){check(factory.writer.closed, "Injected cleanup follows real native closure");}}
		}
		check(descriptors(prefix+".locus")==0 && (factory.writer==null || !factory.writer.isAlive()), "Fixture closes deliberate baseline startup leak");
		observations.append(scenario).tab().append(alive).tab().append(locusDescriptors).tab().append(vectorDescriptors).tab().append(trainerDescriptors)
			.tab().append(caught==null ? "none" : caught.getClass().getName()).tab().append(caught==null ? 0 : caught.getSuppressed().length).nl();
	}
	private static void invoke(String[] args, Factory factory) throws Exception{
		final Method method=NcrnaBoundaryVectorGen.class.getDeclaredMethod("run", String[].class, Function.class);method.setAccessible(true);
		try{method.invoke(null, args, factory);}catch(InvocationTargetException wrapper){final Throwable cause=wrapper.getCause();if(cause instanceof Error){throw (Error)cause;}if(cause instanceof Exception){throw (Exception)cause;}throw new AssertionError(cause);}
	}
	private static int lines(String path){
		final ByteFile input=ByteFile.makeByteFile(path, true);int n=0;try{while(input.nextLine()!=null){n++;}}finally{check(!input.close(), "Complete vector row-count read");}return n;
	}
	private static int descriptors(String path) throws Exception{
		int n=0;try(DirectoryStream<Path> entries=Files.newDirectoryStream(new File("/proc/self/fd").toPath())){
			for(Path entry:entries){try{if(Files.readSymbolicLink(entry).toString().equals(path)){n++;}}catch(NoSuchFileException closedConcurrently){}}
		}return n;
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final class Factory implements Function<String, ByteStreamWriter>{
		Factory(String scenario_){scenario=scenario_;}
		@Override public ByteStreamWriter apply(String path){
			if(scenario.equals("open")){primary=new IllegalArgumentException("declared open failure");throw (RuntimeException)primary;}
			if(scenario.equals("native_open")){check(new File(path).mkdir(), "Native rejection after outer vector acquisition");}
			try{writer=new TrackedWriter(path, this);if(scenario.equals("start")){primary=writer.startup;}return writer;}
			catch(RuntimeException | Error failure){primary=failure;throw failure;}
		}
		final String scenario;TrackedWriter writer;Throwable primary;
		final IllegalArgumentException body=new IllegalArgumentException("declared locus write failure");
		final IllegalStateException cleanup=new IllegalStateException("declared locus close failure");
	}
	private static final class TrackedWriter extends ByteStreamWriter{
		TrackedWriter(String path, Factory factory_){super(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));factory=factory_;}
		@Override public synchronized void start(){if(factory.scenario.equals("start")){throw startup;}super.start();}
		@Override public ByteStreamWriter print(String row){if(factory.scenario.startsWith("header")){throw factory.body;}return super.print(row);}
		@Override public ByteStreamWriter print(ByteBuilder row){if(factory.scenario.startsWith("write")){throw factory.body;}return super.print(row);}
		@Override public synchronized boolean closeBeforeStart(){return finish(super.closeBeforeStart());}
		@Override public boolean poisonAndWait(){return finish(super.poisonAndWait());}
		private boolean finish(boolean error){closed=true;if(factory.scenario.endsWith("close_exception") || factory.scenario.equals("bad_input_close")){throw factory.cleanup;}return error || factory.scenario.endsWith("close_flag");}
		final Factory factory;boolean closed;final AssertionError startup=new AssertionError("declared locus startup failure");
	}
	private final String root, data, tables;
	private final boolean after;
	private static final String[] CASES={"open", "native_open", "start", "header", "header_close_exception", "write", "write_close_flag", "write_close_exception", "close_flag", "close_exception", "ordinary", "omitted", "bad_input", "bad_input_close"};
}
