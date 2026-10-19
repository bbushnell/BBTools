package consensus;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Runs conservative contig-tip substitution polishing through isolated tools.
 * All evidence is generated from the supplied assembly and unmodified reads;
 * no reference genome or truth labels are used. Contexts remain separate mapping
 * universes, matching the validated experimental workflow.
 * @author Brian Bushnell, Fischl
 */
public final class TipPolish {

	/** Expands config files with the standard preparser before orchestration. */
	public static void main(final String[] args) throws IOException, InterruptedException{
		final PreParser pp=new PreParser(args, TipPolish.class, false);
		try{new TipPolish(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns I/O and the frozen recipe; arbitrary child overrides are not accepted. */
	private TipPolish(final String[] args, final PrintStream log){
		assert(args!=null && log!=null) : "PreParser must supply expanded arguments and a diagnostic stream.";
		outstream=log;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("assembly") || a.equals("ref")){assembly=path(b);}
			else if(a.equals("report")){report=path(b);}
			else if(a.equals("temp")){temp=path(b);}
			else if(a.equals("bbtools")){bbtools=path(b);}
			else if(a.equals("childheap")){heap=b;}
			else if(a.equals("interleaved") || a.equals("int")){interleaved=Parse.parseBoolean(b);}
			else if(a.equals("minoverlap")){minOverlap=Parse.parseIntKMG(b);}
			else if(a.equals("maxoverlap")){maxOverlap=Parse.parseIntKMG(b);}
			else if(a.equals("indexed")){indexed=Parse.parseBoolean(b);}
			else if(standardOption(a) && parser.parse(arg, a, b)){
				// Use native aliases and t=auto semantics. Other process-wide Parser
				// options would affect only this coordinator, not its child JVMs.
			}else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		reads=parser.in1==null ? null : path(parser.in1);
		reads2=parser.in2==null ? null : path(parser.in2);
		output=parser.out1==null ? null : path(parser.out1);
		threads=Shared.threads();
		overwrite=parser.overwrite;
		require(!parser.append, "Appending would combine independent assemblies and edit ledgers.");
		require(assembly!=null && reads!=null && output!=null && bbtools!=null,
				"Required: assembly=contigs.fa in=reads.fq out=polished.fa; use tippolish.sh.");
		if(report==null){report=path(output+".edits.tsv");}
		require(threads>0 && minOverlap>=33 && maxOverlap>=minOverlap, "Require t>0 and33<=minoverlap<=maxoverlap.");
		require(heap!=null && heap.matches("[1-9][0-9]*[kKmMgG]?"), "Invalid childheap: "+heap);
		require(reads2==null || !interleaved, "in2 and interleaved=t are mutually exclusive.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Execution          ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds all artifacts privately and publishes only after every stage succeeds. */
	private void process() throws IOException, InterruptedException{
		validate();
		work=Files.createTempDirectory(temp, "tippolish-").toRealPath();
		outstream.println("TipPolish work directory: "+work);
		try{
			final HashSet<String> ambiguous=selectTargets();
			run(work, "prepare", "dev/tipcontextpanel.sh", "mode=prepare", "in="+assembly,
					"targets="+work.resolve("targets.header"), "dir="+work.resolve("panel"),
					"indexed="+indexed, "minoverlap="+minOverlap, "maxoverlap="+maxOverlap);
			boolean mapped=false;
			int cases=0, skipped=0;
			final LineParser1 parser=new LineParser1('\t');
			for(byte[] line : ByteFile.toLines(work.resolve("panel/targets.tsv").toString())){
				parser.set(line);
				if(parser.termEquals("id", 0)){continue;}
				require(parser.terms()==6, "Unexpected prepared case schema.");
				final String id=parser.parseString(0), name=parser.parseString(1);
				require(id.matches("case[0-9]+"), "Invalid prepared case identifier: "+id);
				final int length=parser.parseInt(3), left=parser.parseInt(4), right=parser.parseInt(5);
				if(left<0 && right<0){continue;}
				final Path dir=work.resolve("panel").resolve(id);
				if(ambiguous.contains(name)){
					// Full original-kmer evidence cannot cover an ambiguous base. Abstain,
					// retaining this record as a possible donor rather than deleting it.
					if(left>=0){Files.copy(dir.resolve("original.fa"), dir.resolve("left.guarded.fa"));}
					if(right>=0){Files.copy(dir.resolve("original.fa"), dir.resolve("right.guarded.fa"));}
					Files.write(dir.resolve("abstained.txt"), Arrays.asList("Ambiguous original sequence; no edits permitted."), StandardCharsets.UTF_8);
					skipped++;
					continue;
				}
				if(!mapped){map(work, "whole", assembly, work.resolve("whole.sam")); mapped=true;}
				run(dir, "extract", "samstreamer.sh", "in="+work.resolve("whole.sam"), "out=original.sam",
						"bed=target.bed", "mappedonly=t", "secondary=f", "qfail=f");
				consensus(dir, "original", "original.sam", "original.fa");
				final ArrayList<String> counts=readArgs();
				counts.addAll(Arrays.asList("ref=original.fa", "scaffold="+name, "regions=32-"+length,
						"k=32", "out=original.kmers.tsv"));
				run(dir, "kmers", "dev/kmerwindowdepth.sh", counts.toArray(new String[0]));
				if(left>=0){polishSide(dir, "left", left);}
				if(right>=0){polishSide(dir, "right", right);}
				cases++;
			}
			final Path result=work.resolve(output.toString().endsWith(".gz") ? "polished.fa.gz" : "polished.fa");
			run(work, "merge", "dev/tipcontextpanel.sh", "mode=merge", "in="+assembly,
					"dir="+work.resolve("panel"), "out="+result, "report="+work.resolve("edits.tsv"));
			// Publication follows all successful children. Complete private results
			// remain available if either final copy fails; inputs are always protected.
			validate();
			publish(work.resolve("edits.tsv"), report);
			publish(result, output);
			outstream.println("TIP_POLISH_PASS polished_cases="+cases+" ambiguous_cases="+skipped+" out="+output+" report="+report);
			outstream.println("Inputs unchanged; commands, logs and all evidence retained in "+work);
		}catch(IOException | InterruptedException | RuntimeException e){
			outstream.println("TipPolish failed; evidence retained in "+work);
			throw e;
		}
	}

	/** Records every complete contig and detects unsafe names before mapping. */
	private HashSet<String> selectTargets(){
		assert(work!=null) : "Target names must be written only inside the new run directory.";
		FastaReadInputStream.SPLIT_READS=false;
		final FileFormat format=FileFormat.testInput(assembly.toString(), FileFormat.FASTA, null, true, true);
		require(format.fasta(), "assembly must be FASTA.");
		final Streamer input=StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1);
		final ByteStreamWriter writer=new ByteStreamWriter(work.resolve("targets.header").toString(), false, false, false);
		final HashSet<String> seen=new HashSet<String>(), ambiguous=new HashSet<String>();
		final ByteBuilder row=new ByteBuilder();
		input.start(); writer.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read : batch){
					require(read.id!=null && !read.id.isEmpty() && seen.add(read.id), "Empty or duplicate contig name: "+read.id);
					for(int i=0; i<read.id.length(); i++){
						require(!Character.isWhitespace(read.id.charAt(i)), "Contig names must be whitespace-free: "+read.id);
					}
					require(read.bases!=null && read.length()>0 && read.mate==null, "Empty or paired assembly record: "+read.id);
					for(byte base : read.bases){if(AminoAcid.baseToNumber[base]<0){ambiguous.add(read.id); break;}}
					writer.print(row.clear().append(read.id).nl());
				}
			}
		}finally{
			input.close();
			require(!writer.poisonAndWait() && !input.errorState(), "Assembly input or target-list output failed.");
		}
		require(!seen.isEmpty(), "The assembly contains no contigs.");
		return ambiguous;
	}

	/** Only real local files are accepted, because the evidence is read repeatedly. */
	private void validate() throws IOException{
		assert(output!=null && report!=null) : "Constructor supplies both publication paths.";
		for(Path input : new Path[]{assembly, reads, reads2}){
			if(input==null){continue;}
			require(Files.isRegularFile(input) && Files.isReadable(input), "Unreadable input: "+input);
			require(input.toString().indexOf(',')<0, "Comma-list input paths are not supported: "+input);
		}
		for(Path destination : new Path[]{output, report}){
			require(overwrite || !Files.exists(destination, LinkOption.NOFOLLOW_LINKS), "Output already exists: "+destination);
			require(!Files.isSymbolicLink(destination) && !Files.isDirectory(destination), "Output must not be a symlink or directory: "+destination);
			require(Files.isDirectory(destination.getParent()) && Files.isWritable(destination.getParent()), "Missing/unwritable output parent: "+destination);
			if(Files.exists(destination)){
				for(Path input : new Path[]{assembly, reads, reads2}){
					require(input==null || !Files.isSameFile(input, destination), "Output aliases an input: "+destination);
				}
				for(String tool : TOOLS){require(!Files.isSameFile(destination, bbtools.resolve(tool)), "Output aliases a required launcher.");}
			}
		}
		require(!output.getParent().toRealPath().resolve(output.getFileName()).equals(
				report.getParent().toRealPath().resolve(report.getFileName())), "Assembly output and report must differ.");
		if(Files.exists(output) && Files.exists(report)){require(!Files.isSameFile(output, report), "Output and report are hard-linked.");}
		require(Files.isDirectory(temp) && Files.isWritable(temp), "temp must be an existing writable directory: "+temp);
		for(String tool : TOOLS){require(Files.isRegularFile(bbtools.resolve(tool)), "Missing BBTools launcher: "+tool);}
	}

	/** Uses the identical read pairing and alignment policy for original and context. */
	private void map(final Path dir, final String label, final Path ref, final Path sam) throws IOException, InterruptedException{
		assert(ref!=null && sam!=null) : "Mapping evidence needs explicit reference and output paths.";
		final ArrayList<String> args=readArgs();
		args.addAll(Arrays.asList("fast=t", "nodisk=t", "ambiguous=best", "secondary=f", "ref="+ref, "out="+sam));
		run(dir, label+"_map", "bbmap.sh", args.toArray(new String[0]));
	}

	/** Builds trusted local, unpadded alignment models with the frozen support policy. */
	private void consensus(final Path dir, final String stem, final String sam, final String ref) throws IOException, InterruptedException{
		assert(stem!=null) : "Each model must be named for its original or separately mapped context.";
		run(dir, stem+"_consensus", "consensus.sh", "in="+sam, "ref="+ref, "out="+stem+".polished.fa",
				"outm="+stem+".alm", "mindepth=2", "mafsub=0.5", "mafins=0.5", "mafdel=0.5",
				"usemapq=f", "trimdepth=0", "trimns=f", "noindels=t", "-Xss64m");
	}

	/** The context adds mapping support, never sequence to the final assembly. */
	private void polishSide(final Path dir, final String side, final int offset) throws IOException, InterruptedException{
		assert(offset>=0) : "Prepared context offsets locate the unchanged original target.";
		map(dir, side, dir.resolve(side+".context.fa"), dir.resolve(side+".sam"));
		consensus(dir, side, side+".sam", side+".context.fa");
		run(dir, side+"_guard", "dev/tipcontextguard.sh", "originalmodel=original.alm", "contextmodel="+side+".alm",
				"assembly="+assembly, "originalkmers=original.kmers.tsv", "offset="+offset, "side="+side,
				"minoverlap="+minOverlap, "maxoverlap="+maxOverlap, "out="+side+".guarded.fa",
				"report="+side+".decisions.tsv", "-Xss64m");
	}

	/** Separate children prevent tool statics or kmer tables from leaking between stages. */
	private void run(final Path dir, final String label, final String tool, final String... args) throws IOException, InterruptedException{
		assert(work!=null && dir.startsWith(work)) : "All child artifacts must remain in the owned run directory.";
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", bbtools.resolve(tool).toString(),
				"-Xmx"+heap, "-Xms32m", "-ea", "-eoom", "t="+threads, "overwrite=f"));
		command.addAll(Arrays.asList(args));
		final StringBuilder replay=new StringBuilder("(cd ").append(quote(dir.toString())).append(" &&");
		for(String part : command){replay.append(' ').append(quote(part));}
		replay.append(")");
		Files.write(work.resolve("commands.sh"), Arrays.asList(replay.toString()), StandardCharsets.UTF_8,
				StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		outstream.println("Starting "+dir.getFileName()+"/"+label);
		final long start=System.nanoTime();
		final Process child=new ProcessBuilder(command).directory(dir.toFile()).redirectErrorStream(true)
				.redirectOutput(dir.resolve(label+".log").toFile()).start();
		child.getOutputStream().close();
		final int exit;
		try{exit=child.waitFor();}
		catch(InterruptedException e){child.destroy(); Thread.currentThread().interrupt(); throw e;}
		final String timing=String.format(Locale.ROOT, "%s\t%s\t%.3f\t%d", work.relativize(dir), label, (System.nanoTime()-start)*1e-9, exit);
		Files.write(work.resolve("timing.tsv"), Arrays.asList(timing), StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		if(exit!=0){throw new IOException("Stage "+label+" exited "+exit+"; see "+dir.resolve(label+".log"));}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Only standard settings actually consumed by this coordinator are accepted. */
	private static boolean standardOption(final String key){
		return key.equals("in") || key.equals("in1") || key.equals("input") || key.equals("input1")
				|| key.equals("in2") || key.equals("input2") || key.equals("out") || key.equals("out1")
				|| key.equals("output") || key.equals("output1") || key.equals("t") || key.equals("threads")
				|| key.equals("ow") || key.equals("overwrite") || key.equals("append") || key.equals("app");
	}

	/** Honors explicit overwrite only after validation excludes all evidence aliases. */
	private void publish(final Path source, final Path destination) throws IOException{
		assert(source.startsWith(work)) : "Only completed private results may be published.";
		if(overwrite){Files.copy(source, destination, StandardCopyOption.REPLACE_EXISTING);}
		else{Files.copy(source, destination);}
	}

	/** Reuses the exact same input evidence for every mapping and count pass. */
	private ArrayList<String> readArgs(){
		assert(reads!=null) : "Every evidence stage requires the user's unchanged reads.";
		final ArrayList<String> result=new ArrayList<String>();
		result.add("in="+reads);
		if(reads2!=null){result.add("in2="+reads2);}
		result.add("interleaved="+interleaved);
		return result;
	}

	/** Paths are absolute before children change directory; values are never shell code. */
	private static Path path(final String value){
		require(value!=null && !value.isEmpty(), "A file or directory argument needs a value.");
		return Paths.get(value).toAbsolutePath().normalize();
	}

	/** Quoting is for the reproduction log only; ProcessBuilder receives unmodified arguments. */
	private static String quote(final String value){return "'"+value.replace("'", "'\\''")+"'";}

	/** User input errors must remain fatal even under disabled assertions. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private Path assembly, reads, reads2, output, report, bbtools, work;
	private Path temp=path(System.getProperty("java.io.tmpdir"));
	private String heap="8g";
	private int threads=Shared.threads(), minOverlap=33, maxOverlap=200;
	private boolean interleaved=false, indexed=true, overwrite=false;
	private final PrintStream outstream;
	private static final String[] TOOLS={"dev/tipcontextpanel.sh", "dev/tipcontextguard.sh",
			"dev/kmerwindowdepth.sh", "bbmap.sh", "samstreamer.sh", "consensus.sh"};
}
