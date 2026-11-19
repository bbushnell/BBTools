package test.prok;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import map.LongHashSet;
import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaBoundaryTrace;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaLocusComparison;
import prok.NcrnaPositionalBoundaryRefiner;
import prok.NcrnaScavenger;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import prok.PackagedBoundaryConfig;
import prok.ProkObject;
import prok.RrnaPositionalKmerTable;
import prok.TrnaBoundaryFeatures;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Literal trace lifecycle/frames/score-ownership oracles plus ordinary CLI plumbing.
 * Two native assemblies share a contig name, with opposite sequence orientation
 * and an input without target calls between them to expose source-ordinal drift.
 * Mock sweeps separately cover actual trim success/failure, both scalar orders,
 * physical ends, candidate/path rejection, unknown provenance and exception cleanup.
 * No biological training, genome benchmark or stochastic-output identity is tested.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryTraceTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaBoundaryTraceTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String installation=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("installation")){installation=b;}
				else{check(parser.parse(arg, a, b), "Unknown boundary-trace fixture argument: "+arg);}
			}
			check(installation!=null && parser.out1!=null && !parser.overwrite && !parser.append
				&& parser.in1==null && parser.in2==null && parser.maxReads<0 && parser.samplerate==1
				&& new File(parser.out1).mkdir(), "Supply installation= and out=fresh_fixture_directory");
			Shared.setThreads(1);
			final File root=new File(parser.out1);
			direct(root);
			cli(root, new File(installation));
			pp.outstream.println("NCRNA_BOUNDARY_TRACE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void direct(final File root){
		cutoffState();
		final String path=new File(root, "literal.tsv").getPath();
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		try(Output output=new Output(path)){
			output.writer.print(NcrnaBoundaryTrace.HEADER);
			final NcrnaBoundaryTrace.Writer writer=new NcrnaBoundaryTrace.Writer(output.writer, "sam1");
			for(int ordinal=0; ordinal<10; ordinal++){writer.input(ordinal, "literal_assembly_"+ordinal+".fa");}
			final NcrnaFamily forward=scalar(false), reverse=scalar(true);
			final NcrnaBoundaryTrace.Worker w0=writer.worker(0, 0);
			final NcrnaBoundaryTrace.Family f0=w0.family(forward);
			final TrimSpy c0=caller(forward, f0, 80, 139, true);
			w0.begin("shared", 240, 0);
			try{
				final Orf kept=run(c0, forward, bases, 0, 75, 145, 0);
				f0.returned(kept);
				final Orf pathRejected=run(c0, forward, bases, 0, 75, 145, 0);
				f0.returned(pathRejected);
				run(c0, forward, bases, 0, 75, 145, 0);// Scavenger rejects this accepted candidate later.
				c0.succeeds=false;
				final Orf failedTrim=run(c0, forward, bases, 0, 75, 145, 0);
				f0.returned(failedTrim);
				check(kept.start==79 && kept.stop==139 && failedTrim.start==75 && failedTrim.stop==145,
					"Ungated scalar inference moves only the real eligible trimmed candidate");
				final ArrayList<Orf> selected=new ArrayList<Orf>();
				selected.add(failedTrim);
				selected.add(kept);
				w0.finish(selected);
			}finally{w0.clear();}
			check(w0.empty(), "Successful trace completion must release every candidate and score buffer");

			final NcrnaBoundaryTrace.Worker w1=writer.worker(1, 0);
			final NcrnaBoundaryTrace.Family f1=w1.family(reverse);
			final TrimSpy c1=caller(reverse, f1, 80, 139, true);
			w1.begin("shared", 240, 0);
			try{
				final Orf kept=run(c1, reverse, bases, 1, 75, 145, 0);
				f1.returned(kept);
				check(kept.start==80 && kept.stop==140, "Stop-first scalar inference uses the expected original sense frame");
				kept.flip();
				check(kept.start==99 && kept.stop==159, "In-place reverse flip uses literal239-stop/239-start coordinates");
				w1.finish(one(kept));
			}finally{w1.clear();}

			physical(writer, reverse, bases, 2, 1);
			physical(writer, reverse, bases, 3, 3);
			physical(writer, reverse, bases, 9, 2);
			final NcrnaBoundaryTrace.Worker unknown=writer.worker(4, 0);
			final NcrnaBoundaryTrace.Family unknownFamily=unknown.family(forward);
			unknown.begin("shared", 240, 0);
			try{
				final Orf call=candidate(forward, bases, 0, 80, 139);
				unknownFamily.returned(call);// No eligibility event: retain known=0, never synthesize trim success.
				unknown.finish(one(call));
			}finally{unknown.clear();}
			cleanup(writer, forward, bases);
			final NcrnaBoundaryTrace.Worker empty=writer.worker(8, 0);
			empty.begin("empty", 0, 0);
			try{empty.finish(new ArrayList<Orf>());}finally{empty.clear();}
		}
		positional(root, bases);
		// Positional records have a separate stream so each stream has one declared target.
		final Counts scalar=read(path, true);
		check(scalar.inputs==10 && scalar.candidates==10 && scalar.scores==14 && scalar.summaries==8
			&& scalar.selected==8 && scalar.returned==9 && scalar.unknown==1,
			"Scalar literal records must conserve all ten candidates, eligibility failures and both rejection stages");
		check(scalar.first==6 && scalar.held==4 && scalar.moved==4 && scalar.independent==0 && scalar.reverse==1
			&& scalar.pathRejected==1 && scalar.scavengerRejected==1 && scalar.fixed==3 && scalar.afterAbort==1 && scalar.unavailableScores==1,
			"Both scalar orders, held/moved contexts, physical masks, reverse coordinates and exception reuse are explicit");
		final Counts positional=read(new File(root, "positional.tsv").getPath(), false);
		check(positional.candidates==3 && positional.scores==5 && positional.independent==5
			&& positional.selected==2 && positional.pathRejected==1 && positional.singleEnd==1,
			"Copied positional outputs must retain both candidate histories and explicitly missing five-prime scores");
	}

	/** A disabled finite cutoff and an enabled zero cutoff have different trace contracts. */
	private static void cutoffState(){
		final NcrnaFamily family=scalar(false);
		final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>(1);
		families.add(family);
		check(!family.nnCutoffEnabled && family.nnCutoff==0f, "The public default cutoff is explicitly disabled with a finite backing value");
		NcrnaBoundaryTrace.validateTarget("sam1", families);
		family.nnCutoffEnabled=true;
		boolean rejected=false;
		try{NcrnaBoundaryTrace.validateTarget("sam1", families);}
		catch(IllegalStateException expected){rejected=expected.getMessage().startsWith("Trace only ordinary ungated");}
		check(rejected, "Trace collection must reject an enabled zero cutoff rather than confuse it with the disabled default");
	}

	private static void physical(final NcrnaBoundaryTrace.Writer writer, final NcrnaFamily family,
			final byte[] bases, final int ordinal, final int fixed){
		final NcrnaBoundaryTrace.Worker owner=writer.worker(ordinal, 0);
		final NcrnaBoundaryTrace.Family trace=owner.family(family);
		// The short contig fixes the stop at its real edge and clips the start sweep's -1 slot.
		final byte[] input=fixed==2 ? Arrays.copyOf(bases, 60) : bases;
		final int stop=fixed==3 ? 239 : 59;
		final TrimSpy caller=caller(family, trace, 0, stop, true);
		owner.begin("shared", input.length, 0);
		try{
			final Orf call=run(caller, family, input, 0, 0, stop, fixed);
			trace.returned(call);
			check(call.start==0 && call.stop==(fixed==3 ? 239 : fixed==1 ? 60 : 59), "Physically fixed ends must never be scored or moved");
			owner.finish(one(call));
		}finally{owner.clear();}
	}

	private static void cleanup(final NcrnaBoundaryTrace.Writer writer, final NcrnaFamily family, final byte[] bases){
		final NcrnaBoundaryTrace.Worker owner=writer.worker(7, 0);
		final NcrnaBoundaryTrace.Family trace=owner.family(family);
		try{
			owner.begin("interrupted", 240, 0);
			trace.begin(candidate("interrupted", family, bases, 0, 80, 139), 0, true, true, 0, 0, .8f, .8f, true, 16, 1000, 0, 239, 60);
			throw new IllegalStateException("literal interrupted candidate");
		}catch(IllegalStateException expected){check(expected.getMessage().equals("literal interrupted candidate"), "Do not mask another cleanup-fixture failure");}
		finally{owner.clear();}
		check(owner.empty(), "Exception cleanup must release active-family state too");
		owner.begin("after_abort", 240, 1);
		try{
			final Orf call=new Orf("after_abort", 80, 139, 0, 0, bases, false, ProkObject.RNA);
			call.ncrnaFamily=family.name;
			call.trnaModel=family.modelNames[0];
			trace.begin(call, 0, false, false, 0, 0, .8f, .8f, true, 16, 1000, 0, 239, 60);
			trace.complete(call);
			trace.returned(call);
			owner.finish(one(call));
		}finally{owner.clear();}
	}

	private static void positional(final File root, final byte[] bases){
		final String path=new File(root, "positional.tsv").getPath();
		final byte[] reference=new byte[60];
		Arrays.fill(reference, (byte)'A');
		final NcrnaFamily family=new NcrnaFamily("cobalamin_rf00174", new byte[][]{reference}, null,
			new String[]{"selected_model"}, seeds(), 17, 16, 64, 6, 1, false, 1, 0, 0, 1, 0, 1, .6f, .6f);
		final RrnaPositionalKmerTable.Table five=table(bases, "5prime", 0, 80), three=table(bases, "3prime", -8, 139);
		family.boundaryRefinerFactory=new NcrnaPositionalBoundaryRefiner.Resources(family, five, three,
			new PositionNet(true), new PositionNet(false), 0);
		try(Output output=new Output(path)){
			output.writer.print(NcrnaBoundaryTrace.HEADER);
			final NcrnaBoundaryTrace.Writer writer=new NcrnaBoundaryTrace.Writer(output.writer, "cobalamin");
			for(int ordinal=0; ordinal<=6; ordinal++){writer.input(ordinal, "positional_"+ordinal+".fa");}
			for(int ordinal=5; ordinal<=6; ordinal++){
				final NcrnaBoundaryTrace.Worker owner=writer.worker(ordinal, 0);
				final NcrnaBoundaryTrace.Family trace=owner.family(family);
				final TrimSpy caller=caller(family, trace, 80, 139, true);
				final byte[] input=bases.clone();
				if(ordinal==6){input[80]='N';}
				owner.begin("shared", 240, 0);
				try{
					final Orf first=run(caller, family, input, 0, 75, 145, 0);
					trace.returned(first);
					if(ordinal==5){
						final Orf second=run(caller, family, input, 0, 75, 145, 0);
						trace.returned(second);
						check(first.start==79 && first.stop==140 && second.start==81 && second.stop==138,
							"Reused mock network buffers deliberately change their output between candidates");
					}else{check(first.start==80 && first.stop==140, "A missing five-prime word holds that end while the usable stop moves");}
					owner.finish(one(first));// Zero-margin validation also checks the earlier copied network buffer.
				}finally{owner.clear();}
			}
		}
	}

	private static Orf run(final TrimSpy caller, final NcrnaFamily family, final byte[] bases,
			final int strand, final int start, final int stop, final int fixed){
		final Orf result=candidate(family, bases, strand, start, stop);
		result.ncrnaPartial=fixed;
		caller.finishAcceptedBoundary(result, bases, 0, 0, bases.length-1, 0, .8f, stop-start+1,
			.8f, true, null, null, fixed, null);
		check(result.orfScore==42, "Diagnostic capture must not change the caller's detection score");
		return result;
	}
	private static Orf candidate(final NcrnaFamily family, final byte[] bases, final int strand, final int start, final int stop){
		return candidate("shared", family, bases, strand, start, stop);
	}
	private static Orf candidate(final String name, final NcrnaFamily family, final byte[] bases, final int strand, final int start, final int stop){
		final Orf result=new Orf(name, start, stop, strand, 0, bases, false, ProkObject.RNA);
		result.ncrnaFamily=family.name;
		result.trnaModel=family.modelNames[0];
		result.orfScore=42;
		return result;
	}
	private static ArrayList<Orf> one(final Orf orf){
		final ArrayList<Orf> result=new ArrayList<Orf>();
		result.add(orf);
		return result;
	}
	private static LongHashSet seeds(){
		final LongHashSet result=new LongHashSet();
		result.add(0);
		return result;
	}
	private static NcrnaFamily scalar(final boolean stopFirst){
		final byte[] reference=new byte[60];
		Arrays.fill(reference, (byte)'A');
		final TrnaBoundaryFeatures.NinemerTable table=(seq, pos)->0;
		final CellNet moving=new ScalarNet(false), constant=new ScalarNet(true);
		final int[] offsets={-1, 0, 1};
		return new NcrnaFamily("sam1", new byte[][]{reference}, null, new String[]{"sam1_fixture"}, seeds(),
			17, 16, 64, 6, 1, false, 1, 0, 0, 1, 0, 1, .6f, .6f, 1.01f, .9f,
			stopFirst ? constant : moving, stopFirst ? moving : constant, table, table,
			9, 0, 9, 0, 1000, offsets, offsets, 0, 0);
	}
	private static TrimSpy caller(final NcrnaFamily family, final NcrnaBoundaryTrace.Family trace,
			final int start, final int stop, final boolean trim){
		final TrimSpy caller=new TrimSpy(family, start, stop, trim);
		GeneCaller.applyNcrnaFamilyControls(caller, family);
		caller.setBoundaryTrace(trace);
		return caller;
	}
	private static RrnaPositionalKmerTable.Table table(final byte[] bases, final String end, final int offset, final int truth){
		final RrnaPositionalKmerTable.Table table=new RrnaPositionalKmerTable.Table("RF00174", end, 9, offset, 1.0/262144, 1);
		table.add(bases, truth);
		table.finish();
		return table;
	}

	private static final class TrimSpy extends NcrnaScavenger{
		TrimSpy(final NcrnaFamily f, final int start_, final int stop_, final boolean succeeds_){
			super(f.library, f.models, f.modelNames, f.kmerSet, f.kLong, f.minLen, f.windowPad, f.indexK, f.indexTopN,
				f.adaptive, f.adaptFloor, f.adaptTopFrac, f.adaptQFrac, f.fixedMinHits, f.scoreA, f.scoreB, f.idPass, f.idBorderline,
				f.boundary5NetTemplate, f.boundary3NetTemplate, f.boundaryStartTable, f.boundaryStopTable,
				f.boundaryStartInside, f.boundaryStartOutside, f.boundaryStopInside, f.boundaryStopOutside,
				f.boundaryMeanLen, f.boundaryStartOffsets, f.boundaryStopOffsets, f.boundaryMarginStart, f.boundaryMarginStop);
			start=start_;
			stop=stop_;
			succeeds=succeeds_;
		}
		@Override
		public boolean trimToAlignmentExtent(final Orf orf, final byte[] bases, final int model,
				final int first, final int last, final float identity){
			if(succeeds){
				orf.start=start;
				orf.stop=stop;
			}
			return succeeds;
		}
		final int start, stop;
		boolean succeeds;
	}
	private static final class ScalarNet extends CellNet{
		ScalarNet(final boolean constant_){constant=constant_;}
		@Override
		public CellNet copy(final boolean deltas){return new ScalarNet(constant);}
		@Override
		public CellNet applyInput(final float[] values){
			NcrnaBoundaryTraceTest.check(values.length==10, "Scalar fixture receives the native V1 layout");
			length=Math.round(values[8]*1000);
			return this;
		}
		@Override
		public float feedForward(){return constant ? 1 : length>60 ? .9f : length<60 ? .8995f : .2f;}
		private final boolean constant;
		private int length;
	}
	private static final class PositionNet extends CellNet{
		PositionNet(final boolean start_){start=start_;}
		@Override
		public CellNet copy(final boolean deltas){return new PositionNet(start);}
		@Override
		public int numInputs(){return 6;}
		@Override
		public int numOutputs(){return 3;}
		@Override
		public CellNet applyInput(final float[] values){NcrnaBoundaryTraceTest.check(values.length==6, "Positional fixture uses radius1");return this;}
		@Override
		public float feedForward(){
			Arrays.fill(output, .25f);
			final boolean first=(calls++&1)==0;
			output[start==first ? 0 : 2]=.875f;
			return .875f;
		}
		@Override
		public float[] getOutput(){return output;}
		private final boolean start;
		private final float[] output=new float[3];
		private int calls;
	}

	private static Counts read(final String path, final boolean literal){
		final Counts result=new Counts();
		final LineParser1 lp=new LineParser1('\t');
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.termEquals("I", 0)){
					check(lp.terms()==3, "Input record schema");
					result.inputs++;
					check(lp.parseInt(1)>=0 && lp.parseInt(1)<31, "Literal fixture input ordinals fit its small mask");
					result.inputMask|=1<<lp.parseInt(1);
				}
				else if(lp.termEquals("C", 0)){
					check(lp.terms()==37, "Candidate trace schema");
					check(lp.termEquals("after_abort", 5) ? lp.parseInt(22)==1000 : lp.length(22)==0,
						"The explicit after-abort maximum stays1000; unbounded fixture families omit their maximum");
					result.candidates++;
					result.returned+=lp.parseInt(31);
					result.selected+=lp.parseInt(32);
					if(lp.parseInt(32)!=0){result.selectedInputMask|=1<<lp.parseInt(1);}
					if(lp.parseInt(14)==0){
						result.unknown++;
						for(final int field:new int[]{12, 17, 18, 19, 20, 34, 35, 36}){
							check(lp.length(field)==0, "Unknown provenance must omit unmeasured metadata");
						}
					}
					if(lp.parseInt(32)==0){
						check(lp.length(25)==0 && lp.length(26)==0 && lp.length(27)==0,
							"Rejected candidates have no final selected coordinates or index");
					}
					if(lp.parseInt(16)!=0){result.eligible++;}
					if(lp.parseInt(14)!=0 && lp.parseInt(12)>0){result.fixed++;}
					if(lp.termEquals("PATH_REJECTED", 33)){result.pathRejected++;}
					if(lp.termEquals("SCAVENGER_REJECTED", 33)){result.scavengerRejected++;}
					if(lp.termEquals("after_abort", 5)){result.afterAbort++;}
					if(lp.parseInt(30)==2){result.singleEnd++;}
					if(literal && lp.parseInt(1)==1){
						check(lp.parseInt(6)==1 && lp.parseInt(10)==80 && lp.parseInt(11)==139
							&& lp.parseInt(23)==80 && lp.parseInt(24)==140 && lp.parseInt(25)==99 && lp.parseInt(26)==159,
							"Original sense and final genomic coordinates must remain distinct on the reverse strand");
						result.reverse++;
					}
				}else if(lp.termEquals("S", 0)){
					check(lp.terms()==9, "Score-window schema");
					result.scores++;
					final LineParser1 values=new LineParser1(',').set(lp.parseByteArray(8));
					check(values.terms()==3, "The radius-one fixture retains every positional slot, including missing positions");
					for(int i=0; i<values.terms(); i++){
						if(values.length(i)==0){result.unavailableScores++;}
						else{check(Float.isFinite(Float.intBitsToFloat(values.parseInt(i))), "Every emitted score is a finite measurement");}
					}
					final int context=lp.parseInt(5);
					if(context==NcrnaBoundaryTrace.FIRST){result.first++;}
					else if(context==NcrnaBoundaryTrace.HELD){result.held++;}
					else if(context==NcrnaBoundaryTrace.MOVED){result.moved++;}
					else if(context==NcrnaBoundaryTrace.INDEPENDENT){result.independent++;}
					else{throw new AssertionError("Unknown trace context");}
				}else if(lp.termEquals("T", 0)){
					check(lp.terms()==11, "Contig conservation schema");
					result.summaries++;
				}
				else{throw new AssertionError("Unknown trace record");}
			}
		}
		return result;
	}
	private static final class Counts{
		int inputs, inputMask, selectedInputMask, candidates, scores, summaries, returned, selected, unknown, eligible, first, held, moved,
			independent, reverse, pathRejected, scavengerRejected, fixed, afterAbort, singleEnd, unavailableScores;
	}

	private static void cli(final File root, final File installation) throws Exception{
		final File directory=new File(root, "cli");
		check(directory.mkdir(), "Native CLI fixture requires fresh outputs");
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final NcrnaFamilyLoader.Options off=new NcrnaFamilyLoader.Options();
		off.boundaryOverride=false;
		final NcrnaFamily baseline=NcrnaFamilyLoader.load("sam1", off, resolver);
		final ByteBuilder manifest=new ByteBuilder("#ncrna_family\t1\tsam1\nrole\tmodel\tend\tpath\tformat\tbytes\tsha80\n");
		for(final String role:new String[]{"consensus", "models"}){
			final File source=new File(resolver.resolve("sam1", role, null).path).getAbsoluteFile();
			final String name=role+"_"+source.getName();
			Files.createSymbolicLink(new File(directory, name).toPath(), source.toPath());
			manifest.append(role).append("\t.\t.\t").append(name).tab().append(role.equals("models") ? "hbm" : "fasta").append("\t-\t-\n");
		}
		manifest.append("seeds\t.\t.\tseeds.fa\tfasta\t-\t-\nboundary_config\t.\t.\tboundary.config\tconfig\t-\t-\n");
		write(new File(directory, "manifest.tsv"), manifest.toString());
		final byte[] model=baseline.library[0];
		final ByteBuilder seedText=new ByteBuilder();
		for(int start=0; start+17<=model.length; start++){
			boolean valid=true;
			for(int j=start; j<start+17; j++){
				if(AminoAcid.baseToNumber[model[j]]<0){
					valid=false;
					break;
				}
			}
			if(valid){seedText.append('>').append(start).nl().append(model, start, 17).nl();}
		}
		write(new File(directory, "seeds.fa"), seedText.toString());
		write(new File(directory, "boundary.config"), "family=sam1\nenabled=t\nfeatures=v1\nfeatureorder="+PackagedBoundaryConfig.FEATURE_ORDER
			+"\nmeanlen="+PackagedBoundaryConfig.upperMedianLength(baseline.library)+"\nstartradius=1\nstopradius=1\nmarginstart=0\nmarginstop=0\ncutoff=none\nsharedconfidence=t\nnet="
			+resolver.resolve("sam1", "boundary", null).path+"\ntablestart="+resolver.resolve("sam1", "boundary_start", null).path
			+"\ntablestop="+resolver.resolve("sam1", "boundary_stop", null).path+"\n");
		final byte[] genome=new byte[model.length+200];
		Arrays.fill(genome, (byte)'A');
		System.arraycopy(model, 0, genome, 100, model.length);
		write(new File(directory, "forward.fa"), new ByteBuilder(">shared\n").append(genome).nl().toString());
		final byte[] noTarget=new byte[240];
		Arrays.fill(noTarget, (byte)'N');
		write(new File(directory, "no_target.fa"), new ByteBuilder(">no_target\n").append(noTarget).nl().toString());
		AminoAcid.reverseComplementBasesInPlace(genome);
		write(new File(directory, "reverse.fa"), new ByteBuilder(">shared\n").append(genome).nl().toString());
		for(int mode=0; mode<5; mode++){
			final String name=mode==0 ? "plain" : mode==1 ? "trace_seq" : mode==2 ? "trace_files" : mode==3 ? "bad_target" : "bad_policy";
			final ArrayList<String> command=new ArrayList<String>();
			command.add("bash");
			command.add(new File(installation, "callgenes.sh").getPath());
			for(final String arg:new String[]{"in="+new File(directory, "forward.fa")+","+new File(directory, "no_target.fa")+","+new File(directory, "reverse.fa"),
				"out="+new File(directory, name+".gff"), "sam1root="+directory, "sam1=t", "ncrna=t", "fmn=f", "cobalamin=f",
				"taxonomy=f", "callcds=f", "euk=f", "gffinputordinal=t", "ow=f", "t="+(mode==2 ? 2 : 1), "-Xmx1g"}){command.add(arg);}
			if(mode>0){
				command.add("ncrnaboundarytrace="+new File(directory, name+".tsv"));
				command.add("ncrnaboundarytracefamily="+(mode==3 ? "fmn" : "sam1"));
			}
			if(mode==4){command.add("sam1boundarynn=f");}
			final File log=new File(directory, name+".log");
			final int exit=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log).start().waitFor();
			if(mode>=3){
				check(exit!=0 && !new File(directory, name+".gff").exists() && !new File(directory, name+".tsv").exists()
					&& contains(log, mode==3 ? "Boundary trace target is not enabled" : "Trace only ordinary ungated"),
					"Invalid trace selection must fail before opening writers; inspect "+log);
				continue;
			}
			check(exit==0, "Native trace CLI failed; inspect "+log);
			if(mode>0){
				final Counts counts=read(new File(directory, name+".tsv").getPath(), false);
				check(counts.inputs==3 && counts.inputMask==7 && counts.selectedInputMask==5
					&& counts.selected>=2 && counts.eligible>0 && counts.scores>0,
					"The middle input with no target calls must preserve both flanking ordinals and their real selected/eligible traces");
				join(new File(directory, name+".tsv"), new File(directory, name+".gff"), false);
				boolean rejected=false;
				try{join(new File(directory, name+".tsv"), new File(directory, name+".gff"), true);}
				catch(AssertionError expected){rejected=expected.getMessage().contains("source-qualified GFF");}
				check(rejected, "A deliberately swapped source ordinal must break the actual selected-call trace join");
			}
		}
		// The saved finite worker also uses raw cmp, preserving exact GFF bytes.
		equal(new File(directory, "plain.gff"), new File(directory, "trace_seq.gff"));
		equal(new File(directory, "plain.gff"), new File(directory, "trace_files.gff"));
	}

	/** Exact one-to-one association to real GFF rows; the negative oracle swaps only sources. */
	private static void join(final File trace, final File gff, final boolean swapSources){
		final ArrayList<NcrnaLocusComparison.Locus> calls=NcrnaLocusComparison.readCalls(gff.getPath(), NcrnaFamilyRegistry.displayName("sam1"), 3);
		final boolean[] used=new boolean[calls.size()];
		final LineParser1 lp=new LineParser1('\t');
		int matched=0;
		try(Input input=new Input(trace.getPath())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(!lp.termEquals("C", 0) || lp.parseInt(32)==0){continue;}
				final int recorded=lp.parseInt(1), ordinal=swapSources ? 2-recorded : recorded;
				final int strand=lp.parseInt(6), start=lp.parseInt(25)+1, stop=lp.parseInt(26)+1;
				int found=-1;
				for(int i=0; i<calls.size(); i++){
					final NcrnaLocusComparison.Locus call=calls.get(i);
					if(!used[i] && call.inputOrdinal==ordinal && call.strand==strand && call.start==start && call.stop==stop
						&& lp.termEquals(call.seqid, 5) && lp.termEquals(call.model, 8)){
						check(found<0, "Selected trace is ambiguous against source-qualified GFF rows");
						found=i;
					}
				}
				check(found>=0, "Selected trace has no matching source-qualified GFF row");
				used[found]=true;
				matched++;
			}
		}
		check(matched==calls.size() && matched>=2, "Every source-qualified GFF call needs exactly one selected trace");
	}
	private static boolean contains(final File file, final String text){
		try(Input input=new Input(file.getPath())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(Tools.indexOf(line, text, 0)>=0){return true;}
			}
		}
		return false;
	}
	private static void equal(final File left, final File right){
		try(Input a=new Input(left.getPath()); Input b=new Input(right.getPath())){
			for(byte[] x=a.reader.nextLine(), y=b.reader.nextLine(); x!=null || y!=null; x=a.reader.nextLine(), y=b.reader.nextLine()){
				check(Arrays.equals(x, y), "Trace-disabled and trace-enabled GFF records must be identical: "+left+" "+right);
			}
		}
	}
	private static void write(final File file, final String text){try(Output out=new Output(file.getPath())){out.writer.print(text);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
