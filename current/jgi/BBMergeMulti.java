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
 * Runs progressive BBMerge overlap and target-bridge passes in isolated JVMs.
 * A zero kmer length performs ordinary overlap merging; positive lengths
 * perform reciprocal target-kmer bridging on the preceding residual reads.
 *
 * @author Brian, Noelle
 * @date August 29, 2026
 */
public class BBMergeMulti {

	public static void main(String[] args){
		new BBMergeMulti(args).process();
	}

	/** Returns true when a comma-delimited kmer series requests iterative mode. */
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

	private BBMergeMulti(final String[] args_){
		args=args_;
		parse();
	}

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
			}else if(key.equals("in") || key.equals("in1")){in1=value;}
			else if(key.equals("in2")){in2=value;}
			else if(isMergedOutput(key)){outMerged=value;}
			else if(isUnmergedOutput1(key)){outUnmerged1=value;}
			else if(isUnmergedOutput2(key)){outUnmerged2=value;}
			else if(key.equals("extra")){
				for(final String s : value.split(",")){if(s.length()>0){extra.add(s);}}
			}else if(key.equals("tmpdir") || key.equals("tempdir")){tempDir=value;}
			else if(key.equals("overwrite") || key.equals("ow")){overwrite=!isFalse(value);}
			else if(key.equals("append") || key.equals("app")){
				if(!isFalse(value)){throw new RuntimeException("Iterative BBMerge does not support append=t.");}
			}
			else if(isUnsupportedOutput(key)){
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

	private void parseKmers(final String value){
		final String[] split=value.split(",");
		kmers=new int[split.length];
		for(int i=0; i<split.length; i++){
			kmers[i]=Integer.parseInt(split[i]);
			if(kmers[i]<0){throw new RuntimeException("Kmer lengths must be nonnegative: "+value);}
		}
	}

	private void process(){
		final ArrayList<File> residualTemps=new ArrayList<File>();
		final ArrayList<File> mergedTemps=new ArrayList<File>();
		String current1=in1, current2=in2;
		try{
			for(int pass=0; pass<kmers.length; pass++){
				final boolean last=(pass==kmers.length-1);
				final File residualTemp=(last ? null : makeTemp("residual", pass, outUnmerged1));
				final File mergedTemp=makeTemp("merged", pass, outMerged);
				if(residualTemp!=null){residualTemps.add(residualTemp);}
				mergedTemps.add(mergedTemp);
				final ArrayList<String> command=makeCommand(pass, current1, current2,
						mergedTemp.getAbsolutePath(), last ? outUnmerged1 : residualTemp.getAbsolutePath(),
						last ? outUnmerged2 : null, mergedTemps);
				System.err.println("\nBBMerge iterative pass "+(pass+1)+"/"+kmers.length+", k="+kmers[pass]);
				final Process process=new ProcessBuilder(command).inheritIO().start();
				final int exit=process.waitFor();
				if(exit!=0){throw new RuntimeException("BBMerge iterative pass "+(pass+1)+" failed with exit status "+exit+".");}
				if(pass>0){deleteTemp(residualTemps.get(pass-1));}
				current1=(last ? null : residualTemp.getAbsolutePath());
				current2=null;
			}
			concatenate(mergedTemps);
		}catch(IOException e){
			throw new RuntimeException("Could not launch an iterative BBMerge pass.", e);
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new RuntimeException("Interrupted while waiting for an iterative BBMerge pass.", e);
		}finally{
			for(final File f : residualTemps){if(f.exists()){deleteTemp(f);}}
			for(final File f : mergedTemps){if(f.exists()){deleteTemp(f);}}
		}
	}

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

	private File makeTemp(final String type, final int pass, final String output) throws IOException{
		final File dir=(tempDir==null ? outputDirectory(output) : new File(tempDir));
		if(!dir.exists() && !dir.mkdirs()){
			throw new IOException("Could not create temporary directory "+dir);
		}
		return File.createTempFile("bbmerge_"+type+"_pass"+(pass+1)+"_", tempSuffix(output), dir);
	}

	private static String tempSuffix(final String output){
		final String name=new File(output).getName();
		final int dot=name.indexOf('.');
		return dot<0 ? ".fq" : name.substring(dot);
	}

	private void concatenate(final ArrayList<File> inputs) throws IOException{
		final File output=new File(outMerged);
		if(output.exists() && !overwrite){throw new IOException("Merged output exists and overwrite=f: "+output);}
		final File aggregateDir=outputDirectory(outMerged);
		if(!aggregateDir.exists() && !aggregateDir.mkdirs()){
			throw new IOException("Could not create output directory "+aggregateDir);
		}
		final File aggregate=File.createTempFile("bbmerge_aggregate_", tempSuffix(outMerged), aggregateDir);
		final byte[] buffer=new byte[1<<20];
		final FileOutputStream out=new FileOutputStream(aggregate);
		try{
			for(final File input : inputs){
				final FileInputStream in=new FileInputStream(input);
				try{
					for(int len=in.read(buffer); len>0; len=in.read(buffer)){out.write(buffer, 0, len);}
				}finally{in.close();}
			}
		}finally{out.close();}
		Files.move(aggregate.toPath(), output.toPath(), StandardCopyOption.REPLACE_EXISTING);
	}

	private static File outputDirectory(final String output){
		final File parent=new File(output).getAbsoluteFile().getParentFile();
		return parent==null ? new File(".") : parent;
	}

	private static boolean unsupportedConcatenation(final String output){
		final String lower=output.toLowerCase();
		return lower.endsWith(".zip") || lower.endsWith(".bz2") || lower.endsWith(".xz") ||
				lower.endsWith(".dsrc") || lower.endsWith(".fqz") || lower.endsWith(".ac");
	}

	private static void deleteTemp(final File f){
		if(f.exists() && !f.delete()){throw new RuntimeException("Could not delete temporary file "+f);}
	}

	private static boolean isReplacedArgument(final String key){
		return key.equals("k") || key.equals("kmer") || key.equals("in") || key.equals("in1") || key.equals("in2") ||
				key.equals("extra") || key.equals("append") || key.equals("app") || key.equals("overwrite") || key.equals("ow") ||
				isMergedOutput(key) || isUnmergedOutput1(key) || isUnmergedOutput2(key);
	}

	private static boolean isMergedOutput(final String key){
		return key.equals("out") || key.equals("out1") || key.equals("outm") || key.equals("outm1") ||
				key.equals("outgood") || key.equals("outgood1") || key.equals("outmerged") || key.equals("outmerged1");
	}

	private static boolean isUnmergedOutput1(final String key){
		return key.equals("outu") || key.equals("outu1") || key.equals("outb") || key.equals("outb1") ||
				key.equals("outbad") || key.equals("outbad1") || key.equals("outunmerged") || key.equals("outunmerged1");
	}

	private static boolean isUnmergedOutput2(final String key){
		return key.equals("outu2") || key.equals("outb2") || key.equals("outbad2") || key.equals("outunmerged2");
	}

	private static boolean isUnsupportedOutput(final String key){
		return key.equals("out2") || key.equals("outm2") || key.equals("outgood2") || key.equals("outmerged2") ||
				key.startsWith("outinsert") || key.startsWith("outlength") || key.equals("outi") || key.equals("ihist") ||
				key.equals("hist") || key.equals("histogram") || key.equals("outhist") || key.equals("outa") ||
				key.equals("outadapter") || key.equals("outc") || key.equals("outcardinality");
	}

	private static boolean isFalse(final String value){
		return value!=null && (value.equalsIgnoreCase("f") || value.equalsIgnoreCase("false") || value.equals("0"));
	}

	private static boolean looksLikeInput(final String arg){
		if(arg.startsWith("-") || arg.indexOf('=')>=0){return false;}
		final String lower=arg.toLowerCase();
		return lower.equals("stdin") || lower.startsWith("stdin.") || arg.indexOf('#')>=0 || new File(arg).isFile();
	}

	private static boolean isStdout(final String s){
		final String lower=s.toLowerCase();
		return lower.equals("stdout") || lower.startsWith("stdout.");
	}

	private final String[] args;
	private final ArrayList<Integer> positionalIndices=new ArrayList<Integer>();
	private final ArrayList<String> extra=new ArrayList<String>();
	private int[] kmers;
	private String in1, in2;
	private String outMerged, outUnmerged1, outUnmerged2;
	private String tempDir;
	private boolean overwrite=true;
}
