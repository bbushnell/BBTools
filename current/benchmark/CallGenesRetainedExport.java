package benchmark;

import java.io.IOException;
import java.io.OutputStream;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.IntFunction;

import gff.GffLine;
import prok.Orf;
import prok.ProkObject;
import structures.ByteBuilder;

/** Bounded retained-feature export with explicit per-process compression options.
 * No ReadWrite compression setting is read or changed by an export worker.
 * Payload buffers belong to one active task through compressor completion;
 * failures join active processes and retain any partial files for diagnosis.
 * @author Sayu, Fischl
 */
public final class CallGenesRetainedExport{

	/*--------------------------------------------------------------*/
	/*----------------         Task records         ----------------*/
	/*--------------------------------------------------------------*/

	/** Explicit compressor paths, with no executable-content stamps. */
	public static final class Tools{
		/** Bind canonical executable paths; neither tool is started here. */
		public Tools(Path bgzip_, Path pigz_){
			bgzip=boundTool(bgzip_);
			pigz=boundTool(pigz_);
		}
		/** Compatibility only; historical pin arguments are ignored, never hashed. */
		public Tools(Path bgzip_, String bgzipPin, Path pigz_, String pigzPin){
			this(bgzip_, pigz_);
		}
		/** Canonical executable paths shared read-only by export workers. */
		public final Path bgzip, pigz;
	}

	/** Task-owned output binding borrowing its mutable text buffer without a copy.
	 * Keep the buffer unchanged until its compressor has joined. The constructor
	 * checks directory shape; writeOne separately rejects an existing target file. */
	public static final class Payload{
		/** Validate an existing owned output directory, genome/arm key and borrowed payload. */
		public Payload(Path directory, String genome_, String arm_, ByteBuilder text_, long rows_){
			require(directory!=null && directory.isAbsolute() && Files.isDirectory(directory)
				&& !Files.isSymbolicLink(directory), "Retained output needs a fresh owned directory");
			require(genome_!=null && genome_.matches("g[0-9]{4}")
				&& ("native".equals(arm_) || "formula".equals(arm_)), "Unsafe retained genome/policy key");
			require(text_!=null && text_.length()>0 && rows_>=0, "Retained payload and feature count are required");
			file=directory.resolve(genome_+"."+arm_+".gff.gz");
			genome=genome_;
			arm=arm_;
			text=text_;
			rows=rows_;
		}
		/** Absolute output path for this genome/arm. */
		public final Path file;
		/** Validated genome and native/formula arm labels. */
		public final String genome, arm;
		/** Borrowed buffer held unchanged until export completion. */
		public final ByteBuilder text;
		/** Declared feature count, excluding headers; not recounted by the compressor. */
		public final long rows;
	}

	/** One receipt per file, merged in original task order after all workers join. */
	public static final class Receipt{
		/** Capture completed task measurements and defensively copy the actual process argv. */
		private Receipt(Payload payload, long stored_, long nanos_, boolean bgzf_, int exit_, String[] argv_){
			file=payload.file;
			genome=payload.genome;
			arm=payload.arm;
			rows=payload.rows;
			raw=payload.text.length();
			stored=stored_;
			nanos=nanos_;
			bgzf=bgzf_;
			exit=exit_;
			argv=argv_.clone();
		}
		/** Preserve the existing CDS per-file cost schema. Timing is the measured file task. */
		public void appendCost(ByteBuilder out){
			require(out!=null && raw>0 && stored>0 && nanos>=0, "Incomplete retained-output receipt");
			out.append(genome).tab().append(arm).tab().append(raw).tab().append(stored)
				.tab().appendSlow(nanos*1e-9).tab().append(bgzf ? "BGZF_level9" : "pigz_11").tab();
			appendArgv(out);
			out.nl();
		}
		/** Record the joined process argv and exit beside its exact output population. */
		public void appendProcess(ByteBuilder out, String population){
			require(out!=null && exit==0 && ("CDS".equals(population) || "non_CDS".equals(population)),
				"Only joined successful compressors have complete receipts");
			out.append(genome).tab().append(arm).tab().append(population).tab().append(file.toString()).tab().append(exit).tab();
			appendArgv(out);
			out.nl();
		}
		/** Quote each actual argv token for a readable TSV receipt; no shell executes it. */
		private void appendArgv(ByteBuilder out){
			for(int i=0; i<argv.length; i++){
				if(i>0){out.append(' ');}
				out.append('\'').append(argv[i].replace("'", "'\\''")).append('\'');
			}
		}
		/** Absolute output path for this genome/arm. */
		public final Path file;
		/** Validated genome and native/formula arm labels. */
		public final String genome, arm;
		/** Feature count, uncompressed bytes, stored bytes and whole-task elapsed nanoseconds. */
		public final long rows, raw, stored, nanos;
		/** True for the large-payload BGZF path; false for pigz. */
		public final boolean bgzf;
		/** Joined compressor exit code, zero for a successfully returned receipt. */
		public final int exit;
		/** Owned copy of actual process arguments, used only for reporting. */
		private final String[] argv;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Feature serialization     ----------------*/
	/*--------------------------------------------------------------*/

	/** Capture the original emitted RNA bytes before LeanCalls discards non-CDS rows.
	 * Orf.appendGff preserves scores, partials, models, anticodons and family names. */
	public static final class Rna{
		/** Append non-CDS features in input order using Orf.appendGff.
		 * Null batches add nothing; invalid input can leave a partial capture. */
		public void append(String contig, List<Orf> calls){
			require(contig!=null && !contig.isEmpty() && contig.indexOf(' ')<0 && contig.indexOf('\t')<0,
				"RNA capture requires the emitted GFF contig key");
			if(calls==null){return;}
			for(Orf call:calls){
				require(call!=null, "An emitted feature cannot disappear during RNA capture");
				if(call.type==ProkObject.CDS){continue;}
				require(call.scafName!=null && call.scafName.startsWith(contig)
					&& (call.scafName.length()==contig.length() || call.scafName.charAt(contig.length())==' '
					|| call.scafName.charAt(contig.length())=='\t'), "RNA emitted contig differs from its input binding");
				call.appendGff(text);
				text.nl();
				rows=Math.addExact(rows, 1);
			}
		}
		/** Number of captured non-CDS feature rows, excluding the GFF header. */
		public long rows(){return rows;}
		/** Return the borrowed mutable buffer, not a snapshot. After calling joins,
		 * the owner must keep it unchanged while an export task uses it. */
		public ByteBuilder text(){return text;}
		/** Owned mutable GFF buffer, including its version header. */
		private final ByteBuilder text=new ByteBuilder("##gff-version 3\n");
		/** Captured non-CDS feature count. */
		private long rows;
	}

	/** Return a new headered buffer using the legacy serial retained-CDS formatter.
	 * Borrow the input list and records without mutating them; reject non-CDS rows. */
	public static ByteBuilder cdsText(List<GffLine> lines){
		require(lines!=null, "Retained CDS records are required");
		final ByteBuilder text=new ByteBuilder("##gff-version 3\n");
		for(GffLine line:lines){
			require(line!=null && "CDS".equals(line.type), "The CDS exporter cannot consume RNA features");
			line.appendTo(text);
			text.nl();
		}
		return text;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Bounded export        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reuse the accepted bounded dispatcher and its join-on-failure semantics.
	 * A payload is constructed only by its active task; queued tasks hold no text.
	 * Return receipts in task order only after all workers join. Per-file timing
	 * includes source.apply, compression, writing and process joining; it is not
	 * a compressor-only benchmark. Payload paths must be unique for this invocation. */
	public static ArrayList<Receipt> write(int count, int workers, IntFunction<Payload> source, Tools tools){
		require(count>0 && workers>0 && workers<=count && source!=null && tools!=null,
			"Retained export needs finite tasks, workers, payload source and explicit compressors");
		final ConcurrentHashMap<Path, Boolean> claimed=new ConcurrentHashMap<Path, Boolean>();
		return CallGenesLeanWorkers.run(count, workers, index->{
			final long start=System.nanoTime();
			final Payload payload=source.apply(index);
			require(payload!=null && claimed.putIfAbsent(payload.file, Boolean.TRUE)==null,
				"Missing or duplicate retained output path at task "+index);
			return writeOne(payload, tools, start);
		});
	}

	/** Own and join the external process on every exit, including interrupted writes.
	 * Use BGZF level 9 above 200,000 uncompressed bytes and pigz level 11 otherwise,
	 * each with one worker. The target must be absent; failures can leave partial
	 * output. Restore interruption after process cleanup. start precedes payload creation. */
	private static Receipt writeOne(Payload payload, Tools tools, long start){
		final boolean bgzf=payload.text.length()>200000;
		final String[] argv=bgzf ? new String[]{tools.bgzip.toString(), "-c", "-f", "-l", "9", "-@", "1"}
			: new String[]{tools.pigz.toString(), "-c", "-p", "1", "-11", "-b", "128"};
		Process process=null;
		Throwable failure=null;
		boolean interrupted=false;
		try{
			require(!Files.exists(payload.file) && !Files.isSymbolicLink(payload.file),
				"Retained output collision: "+payload.file);
			Files.createFile(payload.file);
			// Explicit argv avoids both shell interpretation and mutable ReadWrite defaults.
			final ProcessBuilder builder=new ProcessBuilder(argv);
			builder.redirectOutput(payload.file.toFile());
			builder.redirectError(ProcessBuilder.Redirect.INHERIT);
			process=builder.start();
			try(OutputStream stdin=process.getOutputStream()){
				stdin.write(payload.text.array, 0, payload.text.length());
			}
			final int exit=process.waitFor();
			require(exit==0, "Retained compressor exit="+exit+" output="+payload.file);
			final long stored=Files.size(payload.file);
			require(stored>0, "Successful compressor produced no bytes: "+payload.file);
			return new Receipt(payload, stored, System.nanoTime()-start, bgzf, exit, argv);
		}catch(IOException problem){
			failure=problem;
			throw new UncheckedIOException("Retained compression failed: "+payload.file, problem);
		}catch(InterruptedException problem){
			failure=problem;
			interrupted=true;
			throw new IllegalStateException("Interrupted retained export: "+payload.file, problem);
		}catch(RuntimeException | Error problem){
			failure=problem;
			throw problem;
		}finally{
			if(process!=null){
				if(failure!=null){process.destroyForcibly();}
				for(boolean joined=false; !joined; ){
					try{
						process.waitFor();
						joined=true;
					}catch(InterruptedException again){
						interrupted=true;
						process.destroyForcibly();
						if(failure!=null){failure.addSuppressed(again);}
					}
				}
				if(failure!=null){
					failure.addSuppressed(new IllegalStateException("Compressor exit="+process.exitValue()+
						" argv="+java.util.Arrays.toString(argv)+" output="+payload.file));
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Bind an executable's canonical path without reading its payload. */
	private static Path boundTool(Path path){
		try{
			require(path!=null && path.isAbsolute() && Files.isRegularFile(path) && Files.isExecutable(path)
				&& path.equals(path.toRealPath()), "Explicit canonical compressor executable required: "+path);
			return path;
		}catch(IOException failure){throw new UncheckedIOException("Cannot bind retained compressor: "+path, failure);}
	}
	/** Reject invalid export state even when assertions are disabled. */
	private static void require(boolean value, String message){
		if(!value){throw new IllegalArgumentException(message);}
	}
	/** Static export entry points; state belongs to individual captures and tasks. */
	private CallGenesRetainedExport(){}
}
