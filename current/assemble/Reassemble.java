package assemble;

import java.io.File;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import bin.BinObject;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.GiToTaxid;

/** Runs Tadpole separately for each input genome, publishing contigs to one output.
 * Sequential per-file assembly avoids mixing those files in a single assembly;
 * it does not establish taxonomic correctness or eliminate within-file misassembly.
 * Positive IDs parsed from filenames or the first header are passed to Tadpole.
 * Missing IDs are permitted and reported; labels are not resolved against a tree.
 *
 * <p>Input accepts comma-separated paths and one level of FASTA directory children.
 * Existing bare paths are also accepted. Directory order is not sorted, duplicate
 * inputs are retained, and this class does not expand wildcard patterns.</p>
 *
 * <p>Default append mode clears a replaceable existing output once, then appends
 * each genome with a cumulative contig-ID offset. Specifying {@code tempdir=DIR}
 * instead retains per-genome temporary assemblies until final concatenation succeeds.
 * Reported assembly failures produce a final error even with {@code failfast=f};
 * successfully written output may be partial. The object is intended for one serial
 * run and changes shared parser/Tadpole settings; it is not reentrant.</p>
 *
 * <pre>
 * reassemble.sh -Xmx8g in=genomes/ out=assembled.fa k=155
 * reassemble.sh in=genomes/ out=assembled.fa tempdir=assembly-temp/ failfast=t
 * </pre>
 *
 * <p>{@code out=} and at least one input are required; {@code k=} defaults to 31.
 * {@code delete=t} permits temporary-file deletion only after every assembly and
 * final output completion succeed. {@code verbose=t} enables additional diagnostics.
 * Unhandled options pass to Tadpole, but per-genome input/output, append mode,
 * ID offset and statistics settings are subsequently imposed by this wrapper.
 * Sparse-genome defaults are {@code mcs=1}, {@code mce=1}, {@code mincontig=1};
 * explicit user options follow them and retain Tadpole's last-value-wins behavior.</p>
 *
 * @author Brian Bushnell
 * @contributor Noire
 * @date February 9, 2026
 */
public class Reassemble{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs one wrapper and closes its diagnostic stream after successful construction.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		Reassemble x=new Reassemble(args);
		try{x.process(t);}
		finally{Shared.closeStream(x.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Constructor           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Expands input specifications, validates named output aliases and records options.
	 * Construction may create a requested temporary directory but does not assemble.
	 * Unhandled options retain their input order for later Tadpole parsing.
	 * @param args Command-line arguments, expanded through PreParser for config/help
	 * @throws IllegalArgumentException If input/output/K requirements or output protection fail
	 */
	public Reassemble(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		parser.overwrite=true;
		tadpoleArgs=new ArrayList<String>();

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 && !split[1].isEmpty() ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("in1")){
				// Handle multiple input specifications
				if(b!=null){
					String[] files=b.split(",");
					for(String s:files){
						// Expand one directory level; this helper does not expand wildcards.
						ArrayList<String> expanded=Tools.getFileOrFiles(s, true, false, false, false);
						in.addAll(expanded);
					}
				}
			}else if(b==null && new File(arg).exists()){
				// Bare filename as positional argument
				ArrayList<String> expanded=Tools.getFileOrFiles(arg, true, false, false, false);
				in.addAll(expanded);
			}else if(a.equals("out") || a.equals("out1")){
				out=b;
			}else if(a.equals("failfast") || a.equals("strict")){
				failfast=Parse.parseBoolean(b);
			}else if(a.equals("delete") || a.equals("deletetemp")){
				deleteTemp=Parse.parseBoolean(b);
			}else if(a.equals("tempdir") || a.equals("temp")){
				tempdir=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k") || a.equals("kmer")){
				k=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				// Parser handles standard args, but we still need to pass to Tadpole
				tadpoleArgs.add(arg);
			}else{
				// Unknown to us, might be Tadpole-specific
				tadpoleArgs.add(arg);
			}
		}

		// Validate inputs
		if(in.isEmpty() || out==null || k<1){
			throw new IllegalArgumentException("Specify input files, an output, and a positive k-mer length.");
		}
		overwrite=parser.overwrite;
		validateOutput();

		// Determine mode: append (default) or temp files (if tempdir specified)
		useAppendMode=(tempdir==null);

		// Setup temp directory if using temp file mode
		if(!useAppendMode){
			if(!tempdir.endsWith("/")){tempdir=tempdir+"/";}
			// Ensure temp directory exists
			File td=new File(tempdir);
			if(!td.exists()){td.mkdirs();}
		}

		// Setup output file format
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false);
	}

	/** Checks named output replacement and aliases before any output initialization.
	 * Canonical paths and existing hard links are checked against all named inputs;
	 * stream names are skipped. This is a preflight check, not filesystem locking. */
	private void validateOutput(){
		if(FileFormat.isStdio(out)){return;}
		try{
			final File target=new File(out).getCanonicalFile();
			if(target.exists() && (!overwrite || !target.isFile())){
				throw new IllegalArgumentException("Output exists and cannot be replaced: "+out+"; overwrite="+overwrite);
			}
			for(String input:in){
				if(FileFormat.isStdio(input)){continue;}
				final File source=new File(input).getCanonicalFile();
				if(source.equals(target) || (source.exists() && target.exists() && Files.isSameFile(source.toPath(), target.toPath()))){
					throw new IllegalArgumentException("Assembly output aliases an input: "+input+" -> "+out);
				}
			}
		}catch(IOException e){throw new RuntimeException("Could not verify assembly output aliases: "+out, e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Main Process          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Processes genomes serially, retaining successful counters and publication state.
	 * With failfast=false, caught assembly Exceptions permit later genomes to run;
	 * any failed genome still makes the final result an error. Errors and exceptions
	 * outside the per-genome assembly try block propagate immediately.
	 * Temporary outputs are deleted only after every genome and concatenation succeed;
	 * individual deletion return values are not checked. Partial final output may remain.
	 * @param t Timer for tracking execution time
	 * @throws RuntimeException If assembly, output initialization or final output completion fails
	 */
	void process(Timer t){
		outstream.println("Reassemble start");
		outstream.println("Input genomes: "+in.size());
		outstream.println("Output: "+out);
		outstream.println("Mode: "+(useAppendMode ? "append (direct output)" : "temp files (concatenate after)"));

		//[assemble/Reassemble#001 FIXED 2026-07-02] In append mode every genome (incl. the first) writes with app=t; nothing
		//used to clear a pre-existing out file, so re-running with the same out= silently appended the new assembly onto stale
		//contigs. Delete it once up front so the run starts clean (temp-file mode already overwrites out in concatenateOutputs).
		//Matches BBTools' default-overwrite convention.
		validateOutput();
		if(useAppendMode && !FileFormat.isStdio(out)){
			File of=new File(out);
			if(of.exists() && !of.delete()){throw new RuntimeException("Could not clear assembly output: "+out);}
		}

		ArrayList<GenomeRecord> records=new ArrayList<GenomeRecord>();
		long cumulativeContigOffset=0;
		boolean failed=false;

		// Process each genome individually
		for(int i=0; i<in.size(); i++){
			String inputFile=in.get(i);
			outstream.println("\n=== Processing genome "+(i+1)+"/"+in.size()+": "+inputFile+" ===");

			GenomeRecord record=processGenome(inputFile, cumulativeContigOffset);
			records.add(record);

			if(record.success){
				// ASM056: final filtering can leave gaps, so emitted counts do not reserve the ID range.
				cumulativeContigOffset=record.nextContigOffset;
			}

			if(!record.success){
				failed=true;
				if(failfast){throw new RuntimeException("Assembly failed for "+inputFile+": "+record.errorMessage);}
			}
		}

		// If using temp files, concatenate all successful outputs
		if(!useAppendMode){
			concatenateOutputs(records);

			// Clean up temp files
			if(deleteTemp && !failed){
				for(GenomeRecord r:records){
					if(r.tempOutput!=null){
						File f=new File(r.tempOutput);
						if(f.exists()){f.delete();}
					}
				}
			}
		}

		// Print summary statistics
		printSummary(records);

		t.stop();
		outstream.println("\nTime: "+t);
		if(failed){throw new RuntimeException("One or more genome assemblies failed; output may be incomplete. Temporary assemblies were retained.");}
		outstream.println("Reassemble finish");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Output Concatenation  ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies successful temporary assemblies to a newly opened final output, in order.
	 * Input compression is decoded and final-output compression is handled by ReadWrite.
	 * No successful genomes means no final output is opened. Otherwise a missing
	 * successful temporary file is an error, while a zero-byte file contributes no bytes.
	 * Final stream/subprocess completion precedes the caller's temporary-file deletion;
	 * a failure can leave partial final output but preserves those temporary assemblies.
	 * @param records Completed per-genome results in input order
	 */
	private void concatenateOutputs(ArrayList<GenomeRecord> records){
		outstream.println("\n=== Concatenating outputs to "+out+" ===");

		// Count successful assemblies
		int successCount=0;
		for(GenomeRecord r:records){
			if(r.success){successCount++;}
		}
		outstream.println("Successful assemblies: "+successCount+"/"+records.size());

		if(successCount==0){
			outstream.println("ERROR: No successful assemblies to concatenate.");
			return;
		}

		// Concatenate using binary stream copying
		try(AssemblyOutput output=new AssemblyOutput(out)){
			final byte[] buf=new byte[32768];  // 32KB buffer

			for(GenomeRecord r:records){
				if(r.success && r.tempOutput!=null){
					File f=new File(r.tempOutput);
					if(!f.isFile()){throw new IOException("Completed temporary assembly is missing: "+r.tempOutput);}
					if(f.length()>0){
						if(verbose){
							outstream.println("Concatenating: "+r.tempOutput);
						}

						try(InputStream is=ReadWrite.getInputStream(r.tempOutput, false, true, true)){
							for(int lim=is.read(buf); lim>0; lim=is.read(buf)){
								output.stream.write(buf, 0, lim);
							}
						}
					}
				}
			}

		}catch(IOException e){throw new RuntimeException("Failed to concatenate assemblies to "+out, e);}
		outstream.println("Concatenation complete: "+out);
	}

	/**
	 * Prints per-genome success/error rows without changing their state.
	 * Append mode reports contig/base counters; temporary mode reports physical file
	 * bytes. These measures are not interchangeable.
	 * @param records Completed results in processing order
	 */
	private void printSummary(ArrayList<GenomeRecord> records){
		outstream.println("\n=== Assembly Summary ===");
		if(useAppendMode){
			outstream.println("Genome\tTaxID\tStatus\tContigs\tBases");
			for(GenomeRecord r:records){
				outstream.print(new File(r.inputFile).getName()+"\t");
				outstream.print((r.taxID>0 ? r.taxID : "N/A")+"\t");
				outstream.print((r.success ? "SUCCESS" : "FAILED")+"\t");
				if(r.success){
					outstream.println(r.contigsWritten+"\t"+r.basesWritten);
				}else{
					outstream.println("0\t"+r.errorMessage);
				}
			}
		}else{
			outstream.println("Genome\tTaxID\tStatus\tOutputSize");
			for(GenomeRecord r:records){
				outstream.print(new File(r.inputFile).getName()+"\t");
				outstream.print((r.taxID>0 ? r.taxID : "N/A")+"\t");
				outstream.print((r.success ? "SUCCESS" : "FAILED")+"\t");
				if(r.success){
					outstream.println(r.outputSize+" bytes");
				}else{
					outstream.println(r.errorMessage);
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Per-Genome Process    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Selects the label/output path, invokes Tadpole and captures assembly Exceptions.
	 * Filename parsing and argument construction occur outside that exception handler.
	 * A successful record copies the completed Tadpole contig/base counts and exclusive
	 * contig-ID bound, including retained contigs removed by final filtering; temporary
	 * mode additionally observes the file's physical length. No temporary file is deleted.
	 * @param inputFile Path to genome file
	 * @param contigOffset Starting contig ID offset for this genome
	 * @return New success/failure record; a caught Exception leaves an error message
	 */
	private GenomeRecord processGenome(String inputFile, long contigOffset){
		GenomeRecord record=new GenomeRecord();
		record.inputFile=inputFile;

		// Parse taxID from filename or header
		record.taxID=parseTaxIDFromFilename(inputFile);
		if(record.taxID>0){
			outstream.println("Detected taxID: "+record.taxID);
		}else{
			outstream.println("Warning: Could not parse taxID from filename or header: "+inputFile);
		}

		// Determine output file (append mode or temp file)
		String outputFile;
		if(useAppendMode){
			outputFile=out;
		}else{
			// Generate unique temp output filename
			String basename=new File(inputFile).getName();
			// Remove extension
			if(basename.contains(".")){
				basename=basename.substring(0, basename.lastIndexOf('.'));
			}
			record.tempOutput=tempdir+"reassemble_temp_"+basename+"_"+System.nanoTime()+".fa";
			outputFile=record.tempOutput;
		}

		// Sparse-genome defaults precede explicit options, including aliases and combined mincount.
		ArrayList<String> args=new ArrayList<String>();
		args.add("mcs=1");
		args.add("mce=1");
		args.add("mincontig=1");
		args.addAll(tadpoleArgs);
		args.add("in="+inputFile);
		args.add("out="+outputFile);
		//[assemble/Reassemble#001] app=t is correct here now that process() deletes any pre-existing out up front (see there).
		if(useAppendMode){
			args.add("app=t");  // Append mode
		}
		args.add("k="+k);  // K-mer size for assembly
		args.add("idoffset="+contigOffset);  // Unique contig IDs
		args.add("showstats=f");
		if(record.taxID>0){
			args.add("tid="+record.taxID);  // Pass taxID to Tadpole
		}

		String[] argsArray=args.toArray(new String[0]);

		// Log command for reproducibility
		if(verbose){
			outstream.println("Tadpole command: "+String.join(" ", args));
		}

		// Run Tadpole with exception handling
		Tadpole tadpole=null;
		try{
			System.gc();  // Memory management between runs (from TadpoleWrapper pattern)

			// Create Tadpole instance and process
			tadpole=Tadpole.makeTadpole(argsArray, true);
			Timer tadpoleTimer=new Timer();
			tadpole.process(tadpoleTimer);

			record.success=true;
			record.contigsWritten=tadpole.contigsWritten;
			record.basesWritten=tadpole.basesWritten;
			record.nextContigOffset=tadpole.nextContigIDOffset();
			outstream.println("Assembly succeeded: "+record.contigsWritten+" contigs, "+record.basesWritten+" bases");
		}catch(Exception e){
			record.success=false;
			record.errorMessage=e.getMessage();
			outstream.println("ERROR: Assembly failed for "+inputFile);
			outstream.println("Exception: "+e.getMessage());
			if(verbose){
				e.printStackTrace(outstream);
			}
		}

		// Gather statistics from temp output (if using temp files)
		if(!useAppendMode && record.success && record.tempOutput!=null){
			File f=new File(record.tempOutput);
			if(f.exists()){
				record.outputSize=f.length();
			}
		}

		return record;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tries the legacy path parser after removing recognized format/compression
	 * suffixes, then tid_/tid| plus digits in the original basename,
	 * then BinObject parsing of the first logical input line. Only positive IDs count.
	 * Header-read Exceptions are tolerated (logged only when verbose); Errors and
	 * failures from the first legacy parser propagate. No taxonomy-tree resolution occurs.
	 * @param filename Literal input name or path, also used to open the header
	 * @return Positive parsed ID, or -1 when all attempted fallbacks find none
	 */
	private int parseTaxIDFromFilename(String filename){
		// The header-token parser requires digits or field separators, not a filename extension (ASM038).
		int tid=GiToTaxid.parseTaxidNumber(ReadWrite.stripExtension(filename), '_');
		if(tid>0){return tid;}

		// Fallback: simple manual parsing for tid_<digits> in filename
		String name=new File(filename).getName();
		int pos=name.indexOf("tid_");
		if(pos<0){pos=name.indexOf("tid|");}
		if(pos>=0){
			int start=pos+4;
			int end=start;
			while(end<name.length() && Character.isDigit(name.charAt(end))){
				end++;
			}
			if(end>start){
				try{
					tid=Integer.parseInt(name.substring(start, end));
					if(tid>0){return tid;}
				}catch(NumberFormatException e){
					// Continue to header parsing
				}
			}
		}

		// If not found in filename, try parsing from FASTA header (first line)
		try{
			String[] lines=TextFile.toStringLines(filename, 1);
			if(lines!=null && lines.length>0){
				tid=BinObject.parseTaxID(lines[0]);
				if(tid>0){return tid;}
			}
		}catch(Exception e){
			// If file can't be read, just return -1
			if(verbose){
				outstream.println("Warning: Could not read file to parse header: "+e.getMessage());
			}
		}

		return -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Mutable result owned by one serial genome invocation; retained for final reporting. */
	private static class GenomeRecord{
		/** Original input name, used for report labeling. */
		String inputFile;
		/** Positive parsed label, or -1 when no label was found. */
		int taxID=-1;
		/** Per-genome output in temporary mode; null in direct-append mode. */
		String tempOutput;
		/** True only after Tadpole.process returns normally. */
		boolean success=false;
		/** Physical temporary-file size observed after success, rather than a base count. */
		long outputSize=0;
		/** Tadpole contig/base output counters captured after successful processing. */
		long contigsWritten=0;
		long basesWritten=0;
		/** Exclusive absolute ID bound from successful output, including filtered retained contigs. */
		long nextContigOffset=0;
		/** Caught Exception message; may be null, empty on unfailed records. */
		String errorMessage="";
	}

	/** Owns final-output completion so close/subprocess errors prevent temporary-file deletion. */
	private static final class AssemblyOutput implements Closeable{
		/** Opens the final native backend; its owner must close it before deleting any inputs. */
		AssemblyOutput(final String name_){
			name=name_;
			stream=ReadWrite.getOutputStream(name, false, true, true);
		}
		/** Completes the backend and any registered output process, reporting retained errors. */
		@Override
		public void close() throws IOException{
			if(ReadWrite.finishWriting(null, stream, name, true)){throw new IOException("Failed to finish assembly output: "+name);}
		}
		/** Native output backend; bulk writes go directly to it without a per-byte wrapper. */
		final OutputStream stream;
		/** Filename used to wait for any registered compression process. */
		final String name;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Fields                ----------------*/
	/*--------------------------------------------------------------*/

	/** Expanded inputs in supplied/directory-listing order; duplicates are retained. */
	private ArrayList<String> in=new ArrayList<String>();

	/** Required final output name, possibly a stream name. */
	private String out=null;

	/** Optional temporary directory with trailing slash; null selects direct append. */
	private String tempdir=null;

	/** K-mer size for assembly (default 31) */
	private int k=31;

	/** Derived from tempdir during construction; true publishes each genome directly. */
	private boolean useAppendMode=true;

	/** Abort on the first caught assembly failure; false still reports a final error. */
	private boolean failfast=false;

	/** Delete temporary files only after all genomes and final output succeed (default true). */
	private boolean deleteTemp=true;
	/** Whether an existing final output may be replaced; never permits an input alias. */
	private boolean overwrite=true;

	/** Ordered pass-through options, copied before per-genome overrides are appended. */
	private ArrayList<String> tadpoleArgs;

	/** Construction-time output descriptor; actual writes use per-genome/final names. */
	private FileFormat ffout;

	/** PreParser diagnostic stream; the CLI entry point closes it after processing. */
	private PrintStream outstream=System.err;

	/** Shared diagnostic verbosity; persists across instances in the same JVM. */
	public static boolean verbose=false;
}
