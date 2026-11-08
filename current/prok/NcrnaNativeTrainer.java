package prok;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;

/** Runs one native boundary trainer through its public shell entry point.
 * Recipes contain training choices; this driver owns data routes, output and
 * thread budget. Regression training shards stream through native decompression
 * into stdin, with no merged bulk copy. The caller supplies explicit monitoring
 * data; for final refit it is the same single all-eligible input file.
 * Failure/interruption terminates only the trainer's owned process group,
 * including descendants of an unchanged public launcher that waits for Java.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaNativeTrainer{

	/** The two supported native ten-feature endpoint learners. */
	public enum Method{CLASSIFIER, REGRESSION}

	/** Binds the public installation and the maximum worker budget for one trainer. */
	public NcrnaNativeTrainer(final File installation_, final int threads_){
		if(installation_==null || !installation_.isDirectory() || threads_<1){
			throw new IllegalArgumentException("Native training requires the BBTools installation directory and positive thread budget");
		}
		installation=installation_.getAbsoluteFile();
		threads=threads_;
	}

	/** Trains a fresh network, requiring an explicit ten-input/one-output finite recipe.
	 * @param method Native classifier or regression
	 * @param recipe Native config file with dims, seed, epochs/batches and vfraction=0
	 * @param training Complete training shard paths, nonempty
	 * @param validation One monitoring file; same as training[0] for single-file refit
	 * @param output Fresh network path under the caller's output directory
	 */
	public void train(final Method method, final String recipe, final String[] training,
			final String validation, final File output){
		if(method==null || training==null || training.length==0 || output==null || output.exists()){
			throw new IllegalArgumentException("Native training needs a method, nonempty inputs and a fresh network path");
		}
		final String[] arguments=recipe(method, recipe);
		final File script=new File(installation, method==Method.CLASSIFIER ? "train.sh" : "regressiontrainer.sh");
		readable(script.getPath());
		final String[] paths=new String[training.length];
		for(int i=0; i<paths.length; i++){paths[i]=readable(training[i]);}
		final String monitor=readable(validation);
		// The native regression constructor rejects identical in/valin paths.
		// Streaming even a single refit file preserves all-data monitoring without a duplicate file.
		final boolean pipe=method==Method.REGRESSION;
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(script.getPath());
		for(String arg:arguments){command.add(arg);}
		command.add("in="+(pipe ? "stdin.txt" : fileList(paths)));
		command.add((method==Method.CLASSIFIER ? "validate=" : "valin=")+monitor);
		command.add("out="+output.getAbsolutePath());
		command.add("t="+threadCount(method, arguments, threads));
		command.add("vfraction=0");
		// Keep child heap under the same explicit -Xmx setting as the invoking tool.
		// Initial heap must also override train.sh's larger standalone default.
		final long heapMiB=Math.max(1, Runtime.getRuntime().maxMemory()/(1024L*1024));
		NcrnaProcess.addJvmOptions(command, heapMiB);
		final NcrnaProcess.InputWriter input=pipe ? sink->feed(paths, sink) : null;
		NcrnaProcess.run(command, input, "Native trainer");
		if(!output.isFile()){throw new IllegalStateException("Native trainer completed without producing "+output);}
	}

	/** Loads native config syntax once; validates only the builder's scientific contract. */
	private static String[] recipe(final Method method, final String path){
		readable(path);
		final String[] arguments=Parser.parseConfig(new String[]{"config="+path});
		int[] dims=null;
		long seed=-1;
		int epochs=0, nets=1, cycles=1;
		double fraction=Double.NaN;
		for(String arg:arguments){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(owned(a)){throw new IllegalArgumentException("Trainer recipe must not override builder-owned routing/population: "+a);}
			if(a.equals("dims") || a.equals("dimensions")){
				try{
					dims=Parse.parseIntArray(NcrnaNumbers.value(a, b), ",", "x");
				}catch(NumberFormatException failure){
					throw new IllegalArgumentException("Invalid numeric option "+a+"="+b, failure);
				}
			}else if(a.equals("seed") || (method==Method.CLASSIFIER && a.equals("netseed"))){
				seed=NcrnaNumbers.longValue(a, b);
			}else if(a.equals("epochs") || (method==Method.CLASSIFIER && (a.equals("maxepochs") || a.equals("batches")))){
				epochs=method==Method.CLASSIFIER ? NcrnaNumbers.intKMG(a, b) : NcrnaNumbers.intValue(a, b);
			}else if(a.equals("vfraction") || a.equals("validatefraction")){
				fraction=NcrnaNumbers.doubleValue(a, b);
			}else if(a.equals("nets") || a.equals("numnets") || a.equals("networks") || a.equals("numnetworks")){
				nets=NcrnaNumbers.intKMG(a, b);
			}else if(a.equals("cycles")){
				cycles=NcrnaNumbers.intKMG(a, b);
			}
		}
		if(dims==null || dims.length<2 || dims[0]!=10 || dims[dims.length-1]!=1
				|| seed<0 || epochs<1 || fraction!=0 || nets!=1 || cycles!=1){
			throw new IllegalArgumentException("Boundary recipe requires dims=10,...,1, nonnegative seed, positive epochs/batches, vfraction=0 and one net/cycle");
		}
		for(int width:dims){
			if(width<1){throw new IllegalArgumentException("Recipe layer widths must be fixed positive integers");}
		}
		return arguments;
	}

	/** Retains a recipe's worker count within the tool budget; scalar regression defaults to one. */
	private static int threadCount(final Method method, final String[] arguments, final int budget){
		int requested=method==Method.REGRESSION ? 1 : budget;
		for(String arg:arguments){
			final String[] split=Parse.splitOnFirst(arg, '=');
			if(split[0].equalsIgnoreCase("t") || split[0].equalsIgnoreCase("threads")){
				requested=NcrnaNumbers.intKMG(split[0], split.length>1 ? split[1] : null);
			}
		}
		if(requested<1){throw new IllegalArgumentException("Recipe thread count must be positive");}
		return Math.min(requested, budget);
	}

	/** These flags can replace data, outputs or the complete declared population. */
	private static boolean owned(final String a){
		if(a.startsWith("out") || a.startsWith("-xm")){return true;}
		switch(a){
			case "config": case "in": case "input": case "in1": case "in2":
			case "data": case "din": case "dataout":
			case "netout": case "netoutbest": case "networkout":
			case "validate": case "validateset": case "vin": case "vset":
			case "valin": case "validationin":
			case "net": case "network": case "netin": case "networkin": case "innet":
			case "lines": case "maxlines": case "samples": case "maxsamples":
			case "vlines": case "vsamples": case "reads": case "maxreads":
			case "samplerate": case "skipreads": case "train": case "training":
			case "test": case "evaluate": return true;
			default: return false;
		}
	}

	/** Streams each compressed shard once; a separator also handles a missing final newline. */
	private static void feed(final String[] paths, final OutputStream sink) throws IOException{
		assert(paths.length>0) : "Regression stdin must receive every declared training shard, including a one-file refit";
		final byte[] buffer=new byte[65536];
		for(String path:paths){
			try(InputStream source=ReadWrite.getInputStream(path, true, false)){
				for(int n=source.read(buffer); n>=0; n=source.read(buffer)){
					if(n>0){sink.write(buffer, 0, n);}
				}
			}
			sink.write('\n');
		}
	}

	private static String readable(final String path){
		if(path==null){throw new IllegalArgumentException("Native trainer requires explicit input and monitoring paths");}
		final File file=new File(path).getAbsoluteFile();
		if(!file.isFile() || !file.canRead()){throw new IllegalArgumentException("Cannot read trainer input: "+file);}
		return file.getPath();
	}

	/** DataLoader uses commas only when the full joined string is not an existing file. */
	private static String fileList(final String[] paths){
		if(paths.length==1){return paths[0];}
		final StringBuilder list=new StringBuilder();
		for(String path:paths){
			if(path.indexOf(',')>=0){throw new IllegalArgumentException("Multi-file classifier input paths cannot contain commas: "+path);}
			if(list.length()>0){list.append(',');}
			list.append(path);
		}
		return list.toString();
	}

	/** Installation containing the public classifier and regression launchers. */
	private final File installation;
	/** Maximum workers for one child; a smaller explicit recipe count is retained. */
	private final int threads;
}
