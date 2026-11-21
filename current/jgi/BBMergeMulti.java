package jgi;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

import shared.Shared;
import shared.Tools;

/**
 * Runs a nonnegative kmer schedule as sequential BBMerge child JVMs.
 * Zero entries request ordinary overlap merging; positive entries request
 * bridge-only mode with that kmer length. Reciprocal matching follows the
 * setting forwarded to each child.
 * <p>
 * Each pass writes its merged reads to an owned temporary file and supplies
 * its residual reads to the next pass. The last child writes the final residual
 * outputs directly. After all children succeed, merged files are copied in pass
 * order to an aggregate beside the named merged output, then moved into place.
 * This stages the merged output; it does not stage the final residual outputs.
 * Cleanup attempts all registered pass paths on success or failure.
 *
 * @author Brian, Noelle, Sucrose
 * @date August 29, 2026
 */
public class BBMergeMulti{

	/*--------------------------------------------------------------*/
	/*----------------        Entry Points          ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses the wrapper arguments and runs each child before starting the next.
	 * @param args Nonnull command-line array containing nonnull tokens */
	public static void main(final String[] args){
		new BBMergeMulti(args).process();
	}

	/** Returns true when a comma-delimited kmer series requests iterative mode.
	 * Detects the routing syntax; parsing and validation occur in the constructor.
	 * @param args Argument array; null arrays and null entries are ignored
	 * @return Whether a k or kmer option contains a comma after its first equals sign */
	public static boolean isMultiK(final String[] args){
		if(args==null){return false;}
		for(final String arg : args){
			final int equals=(arg==null ? -1 : arg.indexOf('='));
			if(equals<1){continue;}
			final String key=arg.substring(0, equals).toLowerCase();
			if((key.equals("k") || key.equals("kmer")) && arg.indexOf(',', equals+1)>=0){return true;}
		}
		return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains the argument array for later child forwarding, without modifying it.
	 * @param args_ Command-line arguments borrowed for the lifetime of this wrapper */
	private BBMergeMulti(final String[] args_){
		args=args_;
		parse();
	}

	/** Collects wrapper settings, records positional inputs and checks output policy.
	 * Expands a missing paired input template and the primary residual output template
	 * when applicable. Other options remain available for child argument forwarding.
	 * @throws RuntimeException If required settings or supported mode/output rules fail */
	private void parse(){
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final int equals=arg.indexOf('=');
			if(equals<1){
				if(looksLikeInput(arg)){
					if(in1==null){in1=arg; positionalIndices.add(i);}
					else if(in2==null){in2=arg; positionalIndices.add(i);}
				}else if(arg.equalsIgnoreCase("ecco") || arg.equalsIgnoreCase("ecc") ||
						arg.equalsIgnoreCase("errorcorrect") || arg.equalsIgnoreCase("mix")){
					throw new RuntimeException("Iterative BBMerge does not support "+arg+".");
				}
				continue;
			}
			final String key=arg.substring(0, equals).toLowerCase();
			final String value=arg.substring(equals+1);
			if(key.equals("k") || key.equals("kmer")){
				if(value.indexOf(',')>=0){parseKmers(value);}
			}else if(key.equals("in") || key.equals("in1")){
				in1=value;
			}else if(key.equals("in2")){
				in2=value;
			}else if(isMergedOutput(key)){
				outMerged=value;
			}else if(isUnmergedOutput1(key)){
				outUnmerged1=value;
			}else if(isUnmergedOutput2(key)){
				outUnmerged2=value;
			}else if(key.equals("extra")){
				for(final String s : value.split(",")){if(s.length()>0){extra.add(s);}}
			}else if(key.equals("tmpdir") || key.equals("tempdir")){
				tempDir=value;
			}else if(key.equals("overwrite") || key.equals("ow")){
				overwrite=!isFalse(value);
			}else if(key.equals("append") || key.equals("app")){
				if(!isFalse(value)){throw new RuntimeException("Iterative BBMerge does not support append=t.");}
			}else if(isUnsupportedOutput(key)){
				throw new RuntimeException("Iterative BBMerge does not yet support "+key+"; use primary out/outu files only.");
			}else if(key.equals("join") || key.equals("merge")){
				if(isFalse(value)){throw new RuntimeException("Iterative BBMerge requires merge=t.");}
			}else if(key.equals("ecco") || key.equals("ecc") || key.equals("errorcorrect") || key.equals("mix")){
				if(!isFalse(value)){throw new RuntimeException("Iterative BBMerge does not support "+key+".");}
			}
		}
		if(kmers==null || kmers.length<2){throw new RuntimeException("Iterative BBMerge requires at least two comma-delimited kmer lengths.");}
		if(in1==null){throw new RuntimeException("Iterative BBMerge requires an input file.");}
		if(in2==null && in1.indexOf('#')>=0 && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		if(outUnmerged2==null && outUnmerged1!=null && outUnmerged1.indexOf('#')>=0){
			outUnmerged2=outUnmerged1.replaceFirst("#", "2");
			outUnmerged1=outUnmerged1.replaceFirst("#", "1");
		}
		if(outMerged==null || isStdout(outMerged)){
			throw new RuntimeException("Iterative BBMerge requires a named cumulative merged output file.");
		}
		if(outUnmerged1==null || isStdout(outUnmerged1)){
			throw new RuntimeException("Iterative BBMerge requires a named final unmerged output file.");
		}
		if(outMerged.indexOf('#')>=0){throw new RuntimeException("Merged output must be a single file: "+outMerged);}
		if(unsupportedConcatenation(outMerged)){
			throw new RuntimeException("Iterative BBMerge does not support cumulative output compression for "+outMerged);
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, outMerged, outUnmerged1, outUnmerged2)){
			throw new RuntimeException("Some iterative BBMerge input and output file names were specified multiple times.");
		}
		if(!Tools.testOutputFiles(overwrite, false, false, outMerged, outUnmerged1, outUnmerged2)){
			throw new RuntimeException("Cannot write an iterative BBMerge output file with overwrite="+overwrite+".");
		}
	}

	/** Parses the comma-separated integer schedule and rejects negative entries.
	 * The caller separately requires at least two parsed entries.
	 * @param value Kmer lengths in bases, with zero denoting an overlap pass */
	private void parseKmers(final String value){
		final String[] split=value.split(",");
		kmers=new int[split.length];
		for(int i=0; i<split.length; i++){
			kmers[i]=Integer.parseInt(split[i]);
			if(kmers[i]<0){throw new RuntimeException("Kmer lengths must be nonnegative: "+value);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the schedule, stages cumulative merged bytes and attempts pass-temp cleanup.
	 * Pass files are registered before later allocations. A successful next pass
	 * retires the preceding residual file; finally attempts cleanup of the remaining
	 * registered paths. Child interruption is reported only after waiting for that child;
	 * owner interrupt status is restored after cleanup and primary-failure handling.
	 * @throws RuntimeException If a child exits nonzero, I/O fails or this failure type is retained
	 * @throws Error If an Error is retained as the primary failure */
	private void process(){
		final ArrayList<File> residualTemps=new ArrayList<File>();
		final ArrayList<File> mergedTemps=new ArrayList<File>();
		String current1=in1, current2=in2;
		Throwable failure=null;
		boolean interrupted=false;
		try{
			try{
				for(int pass=0; pass<kmers.length; pass++){
					final boolean last=(pass==kmers.length-1);
					final File residualTemp=(last ? null : makeTemp("residual", pass, outUnmerged1));
					// Register each owned temp before allocating its sibling (JGI-085).
					if(residualTemp!=null){residualTemps.add(residualTemp);}
					final File mergedTemp=makeTemp("merged", pass, outMerged);
					mergedTemps.add(mergedTemp);
					final ArrayList<String> command=makeCommand(pass, current1, current2,
							mergedTemp.getAbsolutePath(), last ? outUnmerged1 : residualTemp.getAbsolutePath(),
							last ? outUnmerged2 : null, mergedTemps);
					System.err.println("\nBBMerge iterative pass "+(pass+1)+"/"+kmers.length+", k="+kmers[pass]);
					final Process process=new ProcessBuilder(command).inheritIO().start();
					final int exit=waitFor(process);
					if(exit!=0){throw new RuntimeException("BBMerge iterative pass "+(pass+1)+" failed with exit status "+exit+".");}
					if(pass>0){deleteTemp(residualTemps.get(pass-1));}
					current1=(last ? null : residualTemp.getAbsolutePath());
					current2=null;
				}
				concatenate(mergedTemps);
			}catch(IOException e){
				failure=new RuntimeException("Could not launch an iterative BBMerge pass.", e);
			}catch(InterruptedException e){
				interrupted=true;
				failure=new RuntimeException("Interrupted while waiting for an iterative BBMerge pass.", e);
			}catch(RuntimeException e){
				failure=e;
			}catch(Error e){
				failure=e;
			}finally{
				failure=addSuppressed(failure, cleanupTemps(residualTemps, mergedTemps));
			}
			rethrow(failure);
		}finally{
			if(interrupted){Thread.currentThread().interrupt();}
		}
	}

	/** Builds a child argument vector using the current Java executable and classpath.
	 * Inherits recorded JVM arguments except wrapper heap-control properties and heap
	 * flags replaced by them. The bbmerge.child.xmx/xms properties supply complete JVM
	 * arguments. Wrapper-owned input/output/k/extra options replace their original
	 * aliases; other arguments pass through. Child output uses append=f/overwrite=t
	 * after the wrapper's output checks. Extension and correction modes are disabled.
	 * Extra inputs are added after the first pass or whenever its kmer length is positive.
	 * @param pass Zero-based index in the kmer schedule
	 * @param input1 Current primary input path
	 * @param input2 Optional mate input path; null after the first pass
	 * @param merged Owned merged temporary output for this pass
	 * @param residual1 Primary residual output, temporary except on the final pass
	 * @param residual2 Optional final mate residual output
	 * @param mergedTemps Merged paths registered through this pass, including its output
	 * @return New ProcessBuilder argument list, with no shell interpretation */
	private ArrayList<String> makeCommand(final int pass, final String input1, final String input2,
			final String merged, final String residual1, final String residual2, final ArrayList<File> mergedTemps){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(System.getProperty("java.home")+File.separator+"bin"+File.separator+"java");
		final String childXmx=System.getProperty("bbmerge.child.xmx");
		final String childXms=System.getProperty("bbmerge.child.xms");
		for(final String s : Shared.JVM_ARGS()){
			if(s.startsWith("-Dbbmerge.child.xm")){continue;}
			if(childXmx!=null && s.startsWith("-Xmx")){continue;}
			if(childXms!=null && s.startsWith("-Xms")){continue;}
			command.add(s);
		}
		if(childXmx!=null){command.add(childXmx);}
		if(childXms!=null){command.add(childXms);}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(BBMerge.class.getName());
		for(int i=0; i<args.length; i++){
			if(positionalIndices.contains(i)){continue;}
			final String arg=args[i];
			final int equals=arg.indexOf('=');
			if(equals>0 && isReplacedArgument(arg.substring(0, equals).toLowerCase())){continue;}
			command.add(arg);
		}
		command.add("in="+input1);
		if(input2!=null){command.add("in2="+input2);}
		command.add("out="+merged);
		command.add("outu="+residual1);
		if(residual2!=null){command.add("outu2="+residual2);}
		command.add("append=f");
		command.add("overwrite=t");
		if(pass>0 || kmers[pass]>0){
			final String extras=makeExtra(pass, mergedTemps);
			if(extras.length()>0){command.add("extra="+extras);}
		}
		if(kmers[pass]==0){
			command.add("k=31");
			command.add("bridgeonly=f");
			command.add("kmerbridge=f");
			command.add("useoverlap=t");
		}else{
			command.add("k="+kmers[pass]);
			command.add("bridgeonly=t");
		}
		command.add("extend=0");
		command.add("extend2=0");
		command.add("rem=f");
		command.add("rsem=f");
		command.add("ecct=f");
		command.add("eccb=f");
		command.add("testmerge=f");
		command.add("kfilter=0");
		return command;
	}

	/** Joins original extra inputs and merged outputs from strictly earlier passes.
	 * @param pass Number of preceding merged pass files to include
	 * @param mergedTemps Registered merged files, whose current-pass entry is excluded
	 * @return Comma-separated paths, or an empty string when no extras are available */
	private String makeExtra(final int pass, final ArrayList<File> mergedTemps){
		final StringBuilder sb=new StringBuilder();
		for(final String s : extra){
			if(sb.length()>0){sb.append(',');}
			sb.append(s);
		}
		for(int i=0; i<pass; i++){
			if(sb.length()>0){sb.append(',');}
			sb.append(mergedTemps.get(i).getAbsolutePath());
		}
		return sb.toString();
	}

	/** Creates a pass file in tmpdir when set, otherwise beside the associated output.
	 * The caller registers the returned file before allocating another pass file.
	 * @param type Role used in the temporary filename
	 * @param pass Zero-based pass index, displayed one-based in the filename
	 * @param output Output path supplying the default directory and suffix
	 * @return Newly created, owned temporary file
	 * @throws IOException If the directory or temporary file cannot be created */
	private File makeTemp(final String type, final int pass, final String output) throws IOException{
		final File dir=(tempDir==null ? outputDirectory(output) : new File(tempDir));
		if(!dir.exists() && !dir.mkdirs()){
			throw new IOException("Could not create temporary directory "+dir);
		}
		return File.createTempFile("bbmerge_"+type+"_pass"+(pass+1)+"_", tempSuffix(output), dir);
	}

	/** Retains the basename from its first dot, preserving compound suffixes.
	 * @param output Output path whose basename supplies the suffix
	 * @return Suffix beginning at the first dot, or .fq when no dot occurs */
	private static String tempSuffix(final String output){
		final String name=new File(output).getName();
		final int dot=name.indexOf('.');
		return dot<0 ? ".fq" : name.substring(dot);
	}

	/** Copies inputs through an owned aggregate, preserving primary failures during cleanup.
	 * Copies raw bytes in list order; no decoding or re-encoding occurs. The aggregate
	 * is created beside outMerged, independently of tmpdir. After the overwrite check,
	 * Files.move replaces the destination; no atomic-move option is requested. An
	 * unsuccessful aggregate has deletion attempted, with failures suppressed onto a
	 * prior failure when one exists.
	 * @param inputs Owned merged pass files, in cumulative output order
	 * @throws IOException If an aggregate file operation fails with an IOException */
	private void concatenate(final ArrayList<File> inputs) throws IOException{
		final File output=new File(outMerged);
		if(output.exists() && !overwrite){throw new IOException("Merged output exists and overwrite=f: "+output);}
		final File aggregateDir=outputDirectory(outMerged);
		if(!aggregateDir.exists() && !aggregateDir.mkdirs()){
			throw new IOException("Could not create output directory "+aggregateDir);
		}
		final File aggregate=File.createTempFile("bbmerge_aggregate_", tempSuffix(outMerged), aggregateDir);
		Throwable failure=null;
		boolean moved=false;
		try{
			final byte[] buffer=new byte[1<<20];
			try(FileOutputStream out=new FileOutputStream(aggregate)){
				for(final File input : inputs){
					try(FileInputStream in=new FileInputStream(input)){
						for(int len=in.read(buffer); len>0; len=in.read(buffer)){out.write(buffer, 0, len);}
					}
				}
			}
			Files.move(aggregate.toPath(), output.toPath(), StandardCopyOption.REPLACE_EXISTING);
			moved=true;
		}catch(IOException|RuntimeException|Error e){
			failure=e;
			throw e;
		}finally{
			if(!moved){
				try{Files.deleteIfExists(aggregate.toPath());}
				catch(IOException|RuntimeException|Error e){
					if(failure==null){throw e;}
					addSuppressed(failure, e);
				}
			}
		}
	}

	/** Returns the absolute parent directory, falling back to the current directory. */
	private static File outputDirectory(final String output){
		final File parent=new File(output).getAbsoluteFile().getParentFile();
		return parent==null ? new File(".") : parent;
	}

	/** Tests the selected suffix denylist for cumulative raw-byte concatenation. */
	private static boolean unsupportedConcatenation(final String output){
		final String lower=output.toLowerCase();
		return lower.endsWith(".zip") || lower.endsWith(".bz2") || lower.endsWith(".xz") ||
				lower.endsWith(".dsrc") || lower.endsWith(".fqz") || lower.endsWith(".ac");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Waits for one child, canceling and quiescing it before reporting owner interruption.
	 * @param process Child handle returned by this wrapper's launch
	 * @return Child exit status after normal completion
	 * @throws InterruptedException First owner interruption, with later failures suppressed */
	private static int waitFor(final Process process) throws InterruptedException{
		try{return process.waitFor();}
		catch(InterruptedException e){
			try{
				process.destroyForcibly();
			}catch(RuntimeException e2){
				e.addSuppressed(e2);
			}catch(Error e2){
				e.addSuppressed(e2);
			}
			awaitTermination(process, e);
			throw e;
		}
	}

	/** Waits for the same child after cancellation, preserving later interruptions.
	 * @param process The same owned child, even when cancellation threw
	 * @param failure First interruption, which receives subsequent interruptions */
	private static void awaitTermination(final Process process, final InterruptedException failure){
		boolean waiting=true;
		while(waiting){
			try{
				process.waitFor();
				waiting=false;
			}catch(InterruptedException e){
				failure.addSuppressed(e);
			}
		}
	}

	/** Cleans registered pass temps and preserves the first cleanup RuntimeException or Error.
	 * @param residualTemps Registered residual paths, some already retired between passes
	 * @param mergedTemps Registered merged paths
	 * @return First cleanup failure with later failures suppressed, or null */
	private static Throwable cleanupTemps(final ArrayList<File> residualTemps,
			final ArrayList<File> mergedTemps){
		Throwable failure=null;
		failure=cleanupTemps(residualTemps, failure);
		return cleanupTemps(mergedTemps, failure);
	}

	/** Attempts every extant temp deletion and suppresses later failures onto the first.
	 * @param files Owned temporary paths to inspect in list order
	 * @param failure Existing cleanup failure, or null
	 * @return Existing failure or the first new failure, with subsequent failures suppressed */
	private static Throwable cleanupTemps(final ArrayList<File> files, Throwable failure){
		for(final File f : files){
			try{
				if(f.exists()){deleteTemp(f);}
			}catch(RuntimeException e){
				failure=addSuppressed(failure, e);
			}catch(Error e){
				failure=addSuppressed(failure, e);
			}
		}
		return failure;
	}

	/** Preserves an existing primary failure and suppresses a secondary failure onto it.
	 * @param failure Primary failure, or null
	 * @param secondary Later failure, or null; self-suppression is ignored
	 * @return Primary failure when present, otherwise the secondary failure */
	private static Throwable addSuppressed(final Throwable failure, final Throwable secondary){
		if(secondary==null || failure==secondary){return failure;}
		if(failure==null){return secondary;}
		failure.addSuppressed(secondary);
		return failure;
	}

	/** Rethrows a preserved unchecked primary with its original RuntimeException/Error type.
	 * @param failure RuntimeException or Error to rethrow; null returns normally */
	private static void rethrow(final Throwable failure){
		if(failure==null){return;}
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		throw (Error)failure;
	}

	/** Deletes an extant owned temp, throwing when the deletion reports failure. */
	private static void deleteTemp(final File f){
		if(f.exists() && !f.delete()){throw new RuntimeException("Could not delete temporary file "+f);}
	}

	/** Recognizes lowercase options replaced with wrapper-selected child arguments. */
	private static boolean isReplacedArgument(final String key){
		return key.equals("k") || key.equals("kmer") || key.equals("in") || key.equals("in1") || key.equals("in2") ||
				key.equals("extra") || key.equals("append") || key.equals("app") || key.equals("overwrite") || key.equals("ow") ||
				isMergedOutput(key) || isUnmergedOutput1(key) || isUnmergedOutput2(key);
	}

	/** Recognizes aliases for the one cumulative merged output. */
	private static boolean isMergedOutput(final String key){
		return key.equals("out") || key.equals("out1") || key.equals("outm") || key.equals("outm1") ||
				key.equals("outgood") || key.equals("outgood1") || key.equals("outmerged") || key.equals("outmerged1");
	}

	/** Recognizes aliases for the primary final residual output. */
	private static boolean isUnmergedOutput1(final String key){
		return key.equals("outu") || key.equals("outu1") || key.equals("outb") || key.equals("outb1") ||
				key.equals("outbad") || key.equals("outbad1") || key.equals("outunmerged") || key.equals("outunmerged1");
	}

	/** Recognizes aliases for the optional mate final residual output. */
	private static boolean isUnmergedOutput2(final String key){
		return key.equals("outu2") || key.equals("outb2") || key.equals("outbad2") || key.equals("outunmerged2");
	}

	/** Recognizes unsupported secondary merged outputs and report options. */
	private static boolean isUnsupportedOutput(final String key){
		return key.equals("out2") || key.equals("outm2") || key.equals("outgood2") || key.equals("outmerged2") ||
				key.startsWith("outinsert") || key.startsWith("outlength") || key.equals("outi") || key.equals("ihist") ||
				key.equals("hist") || key.equals("histogram") || key.equals("outhist") || key.equals("outa") ||
				key.equals("outadapter") || key.equals("outc") || key.equals("outcardinality");
	}

	/** Recognizes f, false and 0; other values, including null, are not false tokens. */
	private static boolean isFalse(final String value){
		return value!=null && (value.equalsIgnoreCase("f") || value.equalsIgnoreCase("false") || value.equals("0"));
	}

	/** Recognizes positional stdin aliases, hash templates or existing regular files. */
	private static boolean looksLikeInput(final String arg){
		if(arg.startsWith("-") || arg.indexOf('=')>=0){return false;}
		final String lower=arg.toLowerCase();
		return lower.equals("stdin") || lower.startsWith("stdin.") || arg.indexOf('#')>=0 || new File(arg).isFile();
	}

	/** Recognizes stdout and stdout-suffix aliases, ignoring case. */
	private static boolean isStdout(final String s){
		final String lower=s.toLowerCase();
		return lower.equals("stdout") || lower.startsWith("stdout.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed command-line array, retained without modification for child forwarding. */
	private final String[] args;
	/** Original positional input indices, omitted when forwarding child arguments. */
	private final ArrayList<Integer> positionalIndices=new ArrayList<Integer>();
	/** External extra-read paths, preceding earlier merged-pass files in child extras. */
	private final ArrayList<String> extra=new ArrayList<String>();
	/** Nonnegative kmer lengths in bases, one child per entry; zero requests overlap. */
	private int[] kmers;
	/** Initial primary and optional mate input paths. */
	private String in1, in2;
	/** Cumulative merged output and primary/optional mate final residual paths. */
	private String outMerged, outUnmerged1, outUnmerged2;
	/** Optional pass-file directory override; the aggregate uses the merged-output parent. */
	private String tempDir;
	/** Wrapper-level output replacement policy, checked before child work and assembly. */
	private boolean overwrite=true;
}
