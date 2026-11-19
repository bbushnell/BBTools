package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.function.Consumer;

import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Projects saved 69-input endpoint trials into the retained inside-K9 layout.
 * The sidecar, trial order, selected model and raw/truth spans never change.
 * Only the 21 table scores are recomputed; the three common extras and labels
 * are copied as original text. No alignment, table teaching or jitter drawing.
 * CLI/byte I/O adapted from template.A_SampleByteFile; sequence consumption uses
 * the owned A_SampleStreamerMT pattern through NcrnaWorkers.
 * @author Brian Bushnell, Raiden
 */
public final class ProjectNcrnaPositionalInputs{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Projects one declared end, then closes the diagnostic stream. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final ProjectNcrnaPositionalInputs tool=new ProjectNcrnaPositionalInputs(args);
		try{
			final RrnaPositionalKmerTable.Table table=RrnaPositionalKmerTable.loadScores(tool.table, tool.member);
			final long trials=project(tool.in, tool.rows, tool.contexts, tool.roster, tool.ref, table,
				tool.end, tool.expected, tool.out);
			timer.stop();
			tool.outstream.println(Tools.thingsProcessed(timer.elapsed, trials, 8, "Trials"));
			tool.outstream.println("SAVED_TRIAL_PROJECTION_COMPLETE");
		}finally{Shared.closeStream(tool.outstream);}
	}

	/** Accepts fixed input bindings; output is always fresh and unfiltered. */
	public ProjectNcrnaPositionalInputs(String[] args){
		final PreParser pre=new PreParser(args, null, false);
		args=pre.args;
		outstream=pre.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(a.equals("rows")){
				rows=b;
			}else if(a.equals("contexts")){
				contexts=b;
			}else if(a.equals("roster")){
				roster=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("table")){
				table=b;
			}else if(a.equals("member")){
				member=b;
			}else if(a.equals("expected")){
				expected=Parse.parseKMG(b);
			}else if(a.equals("end")){
				if("5prime".equals(b)){end=0;}
				else if("3prime".equals(b)){end=1;}
				else{throw new IllegalArgumentException("end= must be 5prime or 3prime");}
			}else if(parser.parse(arg, a, b)){
				// Native input/output and shared options.
			}else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		if(in==null || out==null || rows==null || contexts==null || roster==null || ref==null || table==null
				|| member==null || end<0 || expected<1 || parser.append || parser.overwrite || parser.in2!=null
				|| parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
			throw new IllegalArgumentException("Require in=69matrix rows= contexts= roster= ref= table= member= end= expected= out=fresh24matrix");
		}
		if(!Tools.testInputFiles(false, true, in, rows, contexts, roster, ref, table)
				|| !Tools.testOutputFiles(false, false, false, out)
				|| !Tools.testForDuplicateFiles(true, in, rows, contexts, roster, ref, table, out)){
			throw new IllegalArgumentException("Projection needs distinct readable inputs and a fresh output");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Projection         ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the complete trial count. The borrowed table must be finalized.
	 * The original sidecar is reused directly; this method never rewrites it.
	 * A failed invocation may leave a partial matrix, which must not be scored.
	 */
	public static long project(final String matrix, final String sidecar, final String contexts,
			final String roster, final String ref, final RrnaPositionalKmerTable.Table table,
			final int end, final long expected, final String out){
		if(table==null || table.k!=9 || table.radius!=10 || table.sites!=21 || !table.finished
				|| expected<1 || end<0 || end>1 || !table.end.equals(end==0 ? "5prime" : "3prime")
				|| table.offset!=(end==0 ? 0 : -8)){
			throw new IllegalArgumentException("Saved-trial projection requires the finalized radius10 fully-inside K9 table for this end");
		}
		final ProjectNcrnaPositionalInputs tool=new ProjectNcrnaPositionalInputs();
		tool.readTrials(sidecar, expected, end);
		tool.checkRoster(roster);
		tool.readSequences(contexts, read->tool.keepContext(read));
		for(final Locus locus:tool.loci){
			if(locus.bases==null){throw new IllegalArgumentException("Missing source context for "+locus.id);}
		}
		final ArrayList<String> models=new ArrayList<String>();
		final ObjectIntMap<String> modelIndex=new ObjectIntMap<String>(String.class);
		tool.readSequences(ref, read->{
			final String name=NcrnaSupport.recordName(read.id);
			if(modelIndex.contains(name)){throw new IllegalArgumentException("Duplicate selected-library model "+name);}
			modelIndex.put(name, models.size());
			models.add(name);
		});
		final LineParser1 values=new LineParser1('\t'), metadata=new LineParser1('\t');
		final float[] features=new float[69], scores=new float[21];
		final ByteBuilder buffer=new ByteBuilder();
		long trials=0;
		Locus locus=null;
		try(final NcrnaTextIO.Input data=new NcrnaTextIO.Input(matrix);
				final NcrnaTextIO.Input identities=new NcrnaTextIO.Input(sidecar);
				final NcrnaTextIO.Output output=new NcrnaTextIO.Output(out)){
			EvaluateNcrnaPositionalNet.header(data.reader.nextLine(), EvaluateNcrnaPositionalNet.matrixHeader(69), matrix);
			EvaluateNcrnaPositionalNet.header(identities.reader.nextLine(), EvaluateNcrnaPositionalNet.ROW_HEADER, sidecar);
			output.writer.print(buffer.append(EvaluateNcrnaPositionalNet.matrixHeader(24)).nl());
			for(byte[] line=data.reader.nextLine(); line!=null; line=data.reader.nextLine()){
				final byte[] row=identities.reader.nextLine();
				if(row==null){throw new IllegalArgumentException("Sidecar ended before projected trial "+(trials+1));}
				values.set(line);
				final int truth=EvaluateNcrnaPositionalNet.parseVector(values, features, trials+1);
				metadata.set(row);
				EvaluateNcrnaPositionalNet.validateIdentity(metadata, truth, end, trials+1);
				locus=tool.locus(metadata, locus);
				locus.check(metadata);
				final int model=metadata.parseInt(4), rawStart=metadata.parseInt(7), rawStop=metadata.parseInt(8);
				if(model>=models.size() || !metadata.termEquals(models.get(model), 5) || rawStop>=locus.bases.length){
					throw new IllegalArgumentException("Saved selected model or raw span disagrees with fixed library/context at trial "+(trials+1));
				}
				table.scoreWindow(locus.bases, (end==0 ? rawStart : rawStop)+table.offset, scores);
				buffer.clear();
				for(int i=0; i<scores.length; i++){
					if(!Float.isFinite(scores[i])){throw new IllegalArgumentException("A saved trial has unavailable inside-K9 geometry: "+(trials+1));}
					if(i>0){buffer.tab();}
					buffer.appendSlow(scores[i]);
				}
				// Preserve the richer producer's exact common-feature and target tokens.
				for(int i=63; i<66; i++){appendField(buffer.tab(), values, i);}
				for(int i=69; i<90; i++){appendField(buffer.tab(), values, i);}
				output.writer.print(buffer.nl());
				trials++;
			}
			if(identities.reader.nextLine()!=null || trials!=expected){
				throw new IllegalArgumentException("Projection must conserve all "+expected+" paired trials; observed "+trials);
			}
		}
		assert(trials==expected) : "No draw, filtering or model selection is permitted when projecting saved trials";
		return trials;
	}

	/** Retains only the distinct original-locus identities needed by these trials. */
	private void readTrials(final String path, final long expected, final int end){
		final LineParser1 parser=new LineParser1('\t');
		long count=0;
		Locus previous=null;
		try(final NcrnaTextIO.Input input=new NcrnaTextIO.Input(path)){
			EvaluateNcrnaPositionalNet.header(input.reader.nextLine(), EvaluateNcrnaPositionalNet.ROW_HEADER, path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.terms()!=11){throw new IllegalArgumentException("Expected eleven saved-trial metadata fields");}
				final long label=(end==0 ? (long)parser.parseInt(9)-parser.parseInt(7)
					: (long)parser.parseInt(10)-parser.parseInt(8))+10;
				if(label<0 || label>20){throw new IllegalArgumentException("Saved evaluation label is outside its original radius");}
				EvaluateNcrnaPositionalNet.validateIdentity(parser, (int)label, end, ++count);
				if(previous==null || !parser.termEquals(previous.id, 0)){
					final String id=parser.parseString(0);
					int index=byId.get(id);
					if(index<0){
						index=loci.size();
						byId.put(id, index);
						loci.add(new Locus(id, parser));
					}
					previous=loci.get(index);
				}
				previous.check(parser);
			}
		}
		if(count!=expected){throw new IllegalArgumentException("Expected "+expected+" saved trials, found "+count);}
	}

	/** Verifies source-to-assembly joins without requiring or inventing taxonomy. */
	private void checkRoster(final String path){
		final ObjectIntMap<String> sources=new ObjectIntMap<String>(String.class);
		final ArrayList<String> assemblies=new ArrayList<String>();
		final LineParser1 parser=new LineParser1('\t');
		try(final NcrnaTextIO.Input input=new NcrnaTextIO.Input(path)){
			EvaluateNcrnaPositionalNet.header(input.reader.nextLine(),
				NcrnaLearningPopulation.ROSTER_HEADER.getBytes(java.nio.charset.StandardCharsets.US_ASCII), path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.terms()!=11 || parser.parseInt(0)!=assemblies.size()){
					throw new IllegalArgumentException("Original source roster must retain eleven fields and dense indices");
				}
				final String source=parser.parseString(1);
				if(sources.contains(source)){throw new IllegalArgumentException("Duplicate original source "+source);}
				sources.put(source, assemblies.size());
				assemblies.add(parser.parseString(2));
			}
		}
		for(final Locus locus:loci){
			final int index=sources.get(locus.source);
			if(index<0 || !locus.assembly.equals(assemblies.get(index))){
				throw new IllegalArgumentException("Trial source/assembly disagrees with original roster: "+locus.id);
			}
		}
	}

	/** Keeps one native base array per needed locus; all other reads are discarded. */
	private void keepContext(final Read read){
		final int index=byId.get(NcrnaSupport.recordName(read.id));
		if(index<0){return;}
		final Locus locus=loci.get(index);
		final int left=Integer.parseInt(RrnaResourceIO.field(read.id, "lflank"));
		final int right=Integer.parseInt(RrnaResourceIO.field(read.id, "rflank"));
		if(locus.bases!=null || read.bases==null || read.mate!=null || left!=locus.start
				|| right<0 || read.length()-right-1!=locus.stop
				|| !RrnaResourceIO.field(read.id, "source").equals(locus.source)){
			throw new IllegalArgumentException("Duplicate or inconsistent saved source context: "+locus.id);
		}
		locus.bases=read.bases;
	}

	/** Fast path reuses the retained identity while the saved locus is unchanged. */
	private Locus locus(final LineParser1 metadata, final Locus previous){
		if(previous!=null && metadata.termEquals(previous.id, 0)){return previous;}
		final int index=byId.get(metadata.parseString(0));
		if(index<0){throw new IllegalArgumentException("Saved sidecar changed during projection");}
		return loci.get(index);
	}

	/** Copies one retained numeric token without parsing/reformatting its spelling. */
	private static void appendField(final ByteBuilder buffer, final LineParser1 parser, final int term){
		parser.setBounds(term);
		buffer.append(parser.line(), parser.a(), parser.b()-parser.a());
	}

	/** One joined consumer owns retained arrays; failure terminates loudly instead
	 * of abandoning a bounded producer queue. No worker survives stream closure.
	 */
	private void readSequences(final String path, final Consumer<Read> action){
		final boolean trim=Shared.TRIM_READ_DESCRIPTION;
		Shared.TRIM_READ_DESCRIPTION=false;
		try(NcrnaWorkers workers=new NcrnaWorkers();
				NcrnaTextIO.SequenceInput owner=new NcrnaTextIO.SequenceInput(StreamerFactory.makeStreamer(
					FileFormat.testInput(path, FileFormat.FASTA, null, true, true), null, false, -1, false, true, 1))){
			final ArrayList<SequenceThread> threads=new ArrayList<SequenceThread>(1);
			final SequenceThread consumer=new SequenceThread(owner.stream, action);
			threads.add(consumer);
			owner.stream.start();
			workers.startAndJoin(threads, owner.stream);
			if(!consumer.success){throw new IllegalStateException("Sequence input did not complete: "+path);}
		}finally{Shared.TRIM_READ_DESCRIPTION=trim;}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One retained source-qualified identity, shared by all its saved trial rows. */
	private static final class Locus{
		Locus(final String id_, final LineParser1 row){
			id=id_;
			source=row.parseString(1);
			assembly=row.parseString(2);
			start=row.parseInt(9);
			stop=row.parseInt(10);
		}
		void check(final LineParser1 row){
			if(!row.termEquals(id, 0) || !row.termEquals(source, 1) || !row.termEquals(assembly, 2)
					|| row.parseInt(9)!=start || row.parseInt(10)!=stop){
				throw new IllegalArgumentException("Saved trials disagree on original source/truth: "+id);
			}
		}
		final String id, source, assembly;
		final int start, stop;
		byte[] bases;
	}

	private static final class SequenceThread extends Thread{
		SequenceThread(final Streamer input_, final Consumer<Read> action_){
			input=input_;
			action=action_;
		}
		@Override
		public void run(){
			try{
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(final Read read:batch){action.accept(read);}
				}
				success=true;
			}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		}
		final Streamer input;
		final Consumer<Read> action;
		boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ProjectNcrnaPositionalInputs(){}
	private final ObjectIntMap<String> byId=new ObjectIntMap<String>(String.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();
	private String in, rows, contexts, roster, ref, table, member, out;
	private int end=-1;
	private long expected=-1;
	private PrintStream outstream=System.err;
}
