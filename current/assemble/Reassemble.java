package assemble;

import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

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

/**
 * Processes multiple genome files individually through Tadpole assembler while
 * preserving taxonomic ID labels. Reads genome files with taxID in filename
 * (pattern: tid_<number>_...), runs Tadpole on each, and concatenates results
 * to a single output file.
 *
 * <p>This tool eliminates the need for coassembly, preventing chimeric contigs and
 * simplifying the workflow for metagenomic binning evaluation datasets.</p>
 *
 * <h2>Problem Solved:</h2>
 * Traditional metagenomic assembly pipelines coassemble all genomes simultaneously,
 * which leads to:
 * <ul>
 * <li>High memory usage (1000+ genomes at once)</li>
 * <li>Chimeric contigs that confound binning evaluation</li>
 * <li>Complex post-processing: mapping reads → SAM parsing → contig renaming</li>
 * <li>Ambiguous read mappings causing incorrect contig labels</li>
 * </ul>
 *
 * <p>Reassemble processes each genome independently, maintaining ground truth labels
 * directly from filenames without requiring read mapping.</p>
 *
 * <h2>Input Requirements:</h2>
 * Genome files must have taxID in filename following pattern:
 * <ul>
 * <li>tid_<number>_... (e.g., tid_910964_1282131.Streptococcus.asm.fa)</li>
 * <li>tid|<number>|... (alternative delimiter)</li>
 * </ul>
 *
 * <h2>Usage Examples:</h2>
 * <pre>
 * # Basic usage with directory input
 * java -ea -Xmx8g driver.Reassemble in=genomes/ out=assembled.fa k=155
 *
 * # Comma-delimited file list
 * java -ea -Xmx8g driver.Reassemble in=tid_123.fa,tid_456.fa out=output.fa k=31
 *
 * # With custom parameters and fail-fast mode
 * java -ea -Xmx8g driver.Reassemble in=genomes/*.fa out=output.fa k=155 \\
 *     mcs=2 mce=2 mincontig=200 failfast=t tempdir=/tmp/
 * </pre>
 *
 * <h2>Key Features:</h2>
 * <ul>
 * <li>Processes genomes sequentially to avoid memory accumulation</li>
 * <li>Extracts taxID from filename automatically</li>
 * <li>Prevents chimeric contigs through isolated assembly</li>
 * <li>Configurable failure handling (continue vs fail-fast)</li>
 * <li>Automatic temp file management and cleanup</li>
 * <li>Comprehensive statistics reporting per genome</li>
 * </ul>
 *
 * <h2>Parameters:</h2>
 * <h3>Reassemble-specific:</h3>
 * <ul>
 * <li>in=<file|dir> - Input files (comma-delimited, directories, wildcards supported)</li>
 * <li>out=<file> - Output file path (required)</li>
 * <li>failfast=<t/f> - Abort on first failure (default: false)</li>
 * <li>delete=<t/f> - Delete temp files after success (default: true)</li>
 * <li>tempdir=<path> - Temporary file directory (default: output file's directory)</li>
 * <li>verbose=<t/f> - Verbose logging (default: false)</li>
 * </ul>
 *
 * <h3>Tadpole pass-through:</h3>
 * All other parameters are passed to Tadpole. Common ones:
 * <ul>
 * <li>k=<int> - K-mer size for assembly (required)</li>
 * <li>mcs=<int> - minCountSeed (default: 1 in code for sparse genomes)</li>
 * <li>mce=<int> - minCountExtend (default: 1 in code for sparse genomes)</li>
 * <li>mincontig=<int> - Minimum contig length (default: 1 in code)</li>
 * <li>prefilter=<int> - Prefilter level</li>
 * <li>mode=<contig|extend> - Assembly mode</li>
 * </ul>
 *
 * @author Brian Bushnell
 * @contributor Noire
 * @date February 9, 2026
 */
public class Reassemble {

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point. Creates Reassemble instance and executes processing.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		Reassemble x=new Reassemble(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Constructor           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a Reassemble instance by parsing command-line arguments.
	 * Separates Reassemble-specific arguments from Tadpole pass-through arguments.
	 * @param args Command-line arguments
	 */
	public Reassemble(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		tadpoleArgs=new ArrayList<String>();

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("in1")){
				// Handle multiple input specifications
				if(b!=null){
					String[] files=b.split(",");
					for(String s : files){
						// Expand directories and wildcards using Tools.getFileOrFiles
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
		assert(in.size()>0) : "No input files specified. Use in=<file> or in=<directory>";
		assert(out!=null) : "No output file specified. Use out=<file>";
		assert(k>0) : "K-mer size must be specified and positive. Use k=<value>";

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
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, null, true, true, false, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Main Process          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main processing orchestration method.
	 * Processes each input genome through Tadpole, using append mode or temp files.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		outstream.println("Reassemble start");
		outstream.println("Input genomes: "+in.size());
		outstream.println("Output: "+out);
		outstream.println("Mode: "+(useAppendMode ? "append (direct output)" : "temp files (concatenate after)"));

		ArrayList<GenomeRecord> records=new ArrayList<GenomeRecord>();
		long cumulativeContigOffset=0;

		// Process each genome individually
		for(int i=0; i<in.size(); i++){
			String inputFile=in.get(i);
			outstream.println("\n=== Processing genome "+(i+1)+"/"+in.size()+": "+inputFile+" ===");

			GenomeRecord record=processGenome(inputFile, cumulativeContigOffset);
			records.add(record);

			if(record.success){
				// Update cumulative offset for next genome
				cumulativeContigOffset+=record.contigsWritten;
			}

			if(!record.success && failfast){
				outstream.println("ERROR: Assembly failed for "+inputFile+" and failfast=true. Aborting.");
				System.exit(1);
			}
		}

		// If using temp files, concatenate all successful outputs
		if(!useAppendMode){
			concatenateOutputs(records);

			// Clean up temp files
			if(deleteTemp){
				for(GenomeRecord r : records){
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
		outstream.println("Reassemble finish");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Output Concatenation  ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Concatenates all successful assembly outputs into single output file.
	 * Uses binary stream copying pattern from ConcatenateFiles.
	 * @param records List of genome processing records
	 */
	private void concatenateOutputs(ArrayList<GenomeRecord> records){
		outstream.println("\n=== Concatenating outputs to "+out+" ===");

		// Count successful assemblies
		int successCount=0;
		for(GenomeRecord r : records){
			if(r.success){successCount++;}
		}
		outstream.println("Successful assemblies: "+successCount+"/"+records.size());

		if(successCount==0){
			outstream.println("ERROR: No successful assemblies to concatenate.");
			return;
		}

		// Concatenate using binary stream copying
		try{
			final byte[] buf=new byte[32768];  // 32KB buffer
			final OutputStream os=ReadWrite.getOutputStream(out, false, true, true);

			for(GenomeRecord r : records){
				if(r.success && r.tempOutput!=null){
					File f=new File(r.tempOutput);
					if(f.exists() && f.length()>0){
						if(verbose){
							outstream.println("Concatenating: "+r.tempOutput);
						}

						InputStream is=ReadWrite.getInputStream(r.tempOutput, false, true, true);

						for(int lim=is.read(buf); lim>0; lim=is.read(buf)){
							os.write(buf, 0, lim);
						}

						is.close();
					}
				}
			}

			ReadWrite.close(os);
			outstream.println("Concatenation complete: "+out);
		}catch(Exception e){
			outstream.println("ERROR: Failed to concatenate outputs");
			e.printStackTrace(outstream);
		}
	}

	/**
	 * Prints summary table of all genome processing results.
	 * @param records List of genome processing records
	 */
	private void printSummary(ArrayList<GenomeRecord> records){
		outstream.println("\n=== Assembly Summary ===");
		if(useAppendMode){
			outstream.println("Genome\tTaxID\tStatus\tContigs\tBases");
			for(GenomeRecord r : records){
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
			for(GenomeRecord r : records){
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
	 * Processes a single genome file through Tadpole assembler.
	 * @param inputFile Path to genome file
	 * @param contigOffset Starting contig ID offset for this genome
	 * @return GenomeRecord containing processing results
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

		// Build Tadpole argument array
		ArrayList<String> args=new ArrayList<String>(tadpoleArgs);
		args.add("in="+inputFile);
		args.add("out="+outputFile);
		if(useAppendMode){
			args.add("app=t");  // Append mode
		}
		args.add("k="+k);  // K-mer size for assembly
		args.add("mcs=1");  // minCountSeed=1 for sparse genomes
		args.add("mce=1");  // minCountExtend=1 for sparse genomes
		args.add("mincontig=1");  // Keep all contigs initially
		args.add("idoffset="+contigOffset);  // Unique contig IDs
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
	 * Parses taxID from filename, falling back to FASTA header if not found in filename.
	 * First tries filename patterns (tid_<number> or tid|<number>).
	 * If not found, reads first line of file to parse from FASTA header.
	 * @param filename File name or path containing taxID
	 * @return taxID number, or -1 if not found
	 */
	private int parseTaxIDFromFilename(String filename){
		// Try GiToTaxid.parseTaxidNumber (handles tid_ and tid| with various delimiters)
		int tid=GiToTaxid.parseTaxidNumber(filename, '_');
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

	/**
	 * Records processing results for a single genome.
	 */
	private static class GenomeRecord{
		String inputFile;
		int taxID=-1;
		String tempOutput;
		boolean success=false;
		long outputSize=0;
		long contigsWritten=0;
		long basesWritten=0;
		String errorMessage="";
	}

	/*--------------------------------------------------------------*/
	/*----------------        Fields                ----------------*/
	/*--------------------------------------------------------------*/

	/** Input genome file paths */
	private ArrayList<String> in=new ArrayList<String>();

	/** Output file path */
	private String out=null;

	/** Temporary file directory */
	private String tempdir=null;

	/** K-mer size for assembly (default 31) */
	private int k=31;

	/** Use append mode (true) or temp files (false) */
	private boolean useAppendMode=true;

	/** Abort on first failure (default false) */
	private boolean failfast=false;

	/** Delete temporary files after concatenation (default true) */
	private boolean deleteTemp=true;

	/** Arguments to pass through to Tadpole */
	private ArrayList<String> tadpoleArgs;

	/** Output file format */
	private FileFormat ffout;

	/** Output stream for logging */
	private PrintStream outstream=System.err;

	/** Verbose logging flag */
	public static boolean verbose=false;
}
