package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaCallerEvaluation;
import prok.NcrnaCallerProcess;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Ordinary-caller/config integration and finite direct-child lifecycle controls.
 * probe=t runs only the public caller; the external fixture records forwarded JVM arguments.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerEvaluationTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, null, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}

	/** Main closes the selected log after both the probe and full fixture paths. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		boolean probe=false;
		String callerConfig=null;
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("probe")){probe=Parse.parseBoolean(b);}
			else if(a.equals("callerconfig")){callerConfig=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown caller fixture argument: "+arg);}
		}
		if(probe){
			NcrnaCallerProcess.runInstalled(new File(callerConfig), 1024);
			return;
		}
		if(parser.out1==null){throw new IllegalArgumentException("Specify a fresh out= directory");}
		final Path root=Files.createDirectory(new File(parser.out1).toPath());
		final String input=new FamilyResourceResolver().resolve("sixs_rf01685", "consensus", null).path;
		final File installed=new File(input).getParentFile().getParentFile();
		final Path genomes=root.resolve("genomes.txt");
		final Path truth=root.resolve("truth.tsv");
		final Path policy=root.resolve("policy.config");
		write(genomes, new ByteBuilder(input).nl());
		write(policy, new ByteBuilder("taxonomy=f\ncallcds=f\n"));
		truth(input, truth);
		final NcrnaCallerEvaluation evaluation=new NcrnaCallerEvaluation(policy.toString(), genomes.toString(),
			truth.toString(), "sixs_rf01685", "fixture");
		final File calls=evaluation.run(installed, root.resolve("ordinary trial").toFile(), 1, 1024);
		final Result result=evaluation.grade(calls, null);
		assert(result.matched==19 && result.unmatched==0 && result.missed==0)
			: "The nineteen installed consensus contigs are nineteen distinct complete reference loci";
		policyCases();
		final Path duplicate=root.resolve("duplicate.txt");
		write(duplicate, new ByteBuilder(input).nl().append(new File(input).getParent()+"/./"+new File(input).getName()).nl());
		fails(()->new NcrnaCallerEvaluation(policy.toString(), duplicate.toString(), truth.toString(), "sixs_rf01685", "fixture"), "Duplicate evaluation genome");
		final Path commaGenome=root.resolve("genome,comma.fa");
		write(commaGenome, new ByteBuilder(">fixture\nACGT\n"));
		final Path commaList=root.resolve("comma.txt");
		write(commaList, new ByteBuilder(commaGenome.toAbsolutePath().toString()).nl());
		fails(()->NcrnaCallerEvaluation.readGenomes(commaList.toString()), "without commas");
		lifecycle(root, policy.toFile());
		System.out.println("Caller evaluation passed: nineteen complete detections, native config routes, policy ownership, duplicate genomes and child failure/interruption cleanup");
	}

	private static void policyCases(){
		NcrnaCallerEvaluation.validatePolicy(new String[]{"taxonomy=f", "callcds=f", "trnaboundarynet=f"}, "sixs_rf01685");
		NcrnaCallerEvaluation.validatePolicy(new String[]{"fmnroot=/unrelated"}, "rnasep");
		for(String arg:new String[]{"in=/alternate", "out=/alternate", "reads=1", "config=/nested",
			"ncrnafamilyroot=/alternate", "ncrnaboundarynet=f", "sixsrf01685boundarynn=t", "sixsrf01685kmers=/alternate"}){
			fails(()->NcrnaCallerEvaluation.validatePolicy(new String[]{arg}, "sixs_rf01685"), null);
		}
	}

	/** Truth is intrinsic to these synthetic contigs: one complete consensus sequence per contig. */
	private static void truth(final String fasta, final Path target){
		int count=0;
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(
				FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true), null, false, -1, false, true, 1));
				Output output=new Output(target.toString())){
			final Streamer reader=input.stream;
			final ByteStreamWriter writer=output.writer;
			reader.start();
			writer.print("seqid\tfamily\tstart\tstop\tstrand\n");
			final ByteBuilder row=new ByteBuilder();
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				for(Read read:batch){
					final int space=read.id.indexOf(' ');
					final String id=space<0 ? read.id : read.id.substring(0, space);
					row.clear();
					writer.print(row.append(id).append("\tfixture\t1\t").append(read.length()).append("\t+\n"));
					count++;
				}
			}
		}
		assert(count==19) : "Fixture population must be the known nineteen RF01685 consensuses";
	}

	private static void lifecycle(final Path root, final File config) throws Exception{
		final Path failed=mock(root, "failed", "exit 7\n");
		fails(()->NcrnaCallerProcess.run(failed.toFile(), config, 256), "exited 7");
		gone(failed.getParent());
		final Path pre=mock(root, "preinterrupted", "exec sleep 60\n");
		Thread.currentThread().interrupt();
		try{
			fails(()->NcrnaCallerProcess.run(pre.toFile(), config, 256), "before launch");
			assert(Thread.currentThread().isInterrupted() && !Files.exists(pre.getParent().resolve("pid")))
				: "An existing interrupt must be retained and must prevent launch";
		}finally{Thread.interrupted();}
		final Path blocked=mock(root, "interrupted", "exec sleep 60\n");
		final AtomicBoolean retained=new AtomicBoolean();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Thread worker=new Thread(()->{
			try{NcrnaCallerProcess.run(blocked.toFile(), config, 256);}
			catch(Throwable caught){failure.set(caught);}
			finally{retained.set(Thread.currentThread().isInterrupted());}
		});
		worker.start();
		final Path pid=blocked.getParent().resolve("pid");
		final long deadline=System.nanoTime()+5_000_000_000L;
		while(!Files.exists(pid) && worker.isAlive() && System.nanoTime()<deadline){Thread.sleep(10);}
		try(AutoCloseable cleanup=()->{
			worker.interrupt();
			worker.join(5000);
		}){
			assert(Files.exists(pid)) : "Owned child must announce startup before interruption";
		}
		assert(!worker.isAlive() && retained.get() && failure.get()!=null)
			: "Interruption must finish child cleanup, preserve interruption and report failure";
		gone(blocked.getParent());
	}

	private static Path mock(Path root, String name, String body) throws Exception{
		final Path directory=Files.createDirectory(root.resolve(name));
		final Path script=directory.resolve("caller.sh");
		write(script, new ByteBuilder("#!/bin/bash\nprintf '%s\\n' \"$$\" > ")
			.append(quote(directory.resolve("pid").toString())).nl().append(body));
		return script;
	}

	private static String quote(String value){return "'"+value.replace("'", "'\\''")+"'";}

	private static void gone(final Path directory) throws Exception{
		final ByteFile input=ByteFile.makeByteFile(directory.resolve("pid").toString(), false);
		final byte[] row=input.nextLine();
		if(input.close() || row==null){throw new IllegalStateException("Could not read owned mock PID");}
		final String pid=new String(row, StandardCharsets.US_ASCII);
		assert(pid.matches("[0-9]+")) : "The mock writes only its own numeric PID";
		final int result=new ProcessBuilder("kill", "-0", pid).redirectError(directory.resolve("probe.err").toFile()).start().waitFor();
		assert(result!=0) : "The owned child must already be gone when its caller returns";
	}

	private static void fails(Runnable action, String message){
		boolean rejected=false;
		try{action.run();}
		catch(IllegalArgumentException | IllegalStateException expected){
			rejected=message==null || expected.getMessage().contains(message);
		}
		assert(rejected) : "The intended invalid input or child failure must be reported";
	}

	private static void write(Path path, ByteBuilder text){
		try(Output owner=new Output(path.toString())){
			owner.writer.print(text);
		}
	}
}
