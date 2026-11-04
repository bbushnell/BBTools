package prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectSet;
import parse.Parse;
import parse.Parser;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaTextIO.Input;

/** Immutable evaluation inputs for ordinary combined-genome caller trials.
 * Each trial gets one closed native configuration and one GFF in a fresh directory.
 * Repeated in= records avoid OS argument-size limits and preserve the full genome list.
 * Independent run() calls may execute concurrently with explicit thread/heap shares;
 * grade after OFF completes to protect its matched reference IDs during comparison.
 * The genome/truth population is never filtered by learning exclusions.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerEvaluation{

	/** Loads immutable policy, complete genome paths and the selected truth-family population. */
	public NcrnaCallerEvaluation(final String callerConfig, final String genomes, final String truth,
			final String family_, final String truthFamily){
		if(callerConfig==null || genomes==null || truth==null || truthFamily==null){
			throw new IllegalArgumentException("Caller evaluation requires policy, genome list, truth and an explicit truth-family label");
		}
		family=CallGenes.parseNcrnaFamily(family_);
		NcrnaFamilyRegistry.sharedBoundaryMember(family);
		displayFamily=NcrnaFamilyRegistry.displayName(family);
		final String config=NcrnaCallerProcess.configurationPath(new File(callerConfig));
		arguments=Parser.parseConfig(new String[]{"config="+config});
		validatePolicy(arguments, family);
		genomePaths=readGenomes(genomes);
		references=NcrnaLocusComparison.readTruth(truth, truthFamily);
	}

	/** The resource directory already declares OFF or a candidate's complete learned settings. */
	public File run(final File familyRoot, final File directory, final int threads, final long heapMiB){
		if(familyRoot==null || !familyRoot.isDirectory() || directory==null || threads<1 || heapMiB<1){
			throw new IllegalArgumentException("A caller trial needs a family directory and positive thread/heap shares");
		}
		final String root=recordValue(familyRoot.getAbsolutePath());
		final File config=new File(directory, "caller.config").getAbsoluteFile();
		final File gff=new File(directory, "calls.gff").getAbsoluteFile();
		recordValue(gff.getPath());
		if(!directory.mkdir()){throw new IllegalArgumentException("Caller output directory must be fresh: "+directory);}
		NcrnaCallerProcess.writeConfiguration(config, arguments, genomePaths, gff, family, threads,
			"ncrnafamilyroot="+root);
		NcrnaCallerProcess.runInstalled(config, heapMiB);
		if(!gff.isFile()){throw new IllegalStateException("Ordinary caller exited without its GFF: "+gff);}
		return gff;
	}

	/** Uses the same complete truth object for every trial and the public family-aware parser. */
	public Result grade(final File gff, final Result off){
		if(gff==null || !gff.isFile()){throw new IllegalArgumentException("A completed caller GFF is required for grading");}
		return NcrnaBoundaryComparison.grade(references, NcrnaLocusComparison.readCalls(gff.getPath(), displayFamily), off);
	}

	/** Selected resources and refinement belong to the trial; unrelated family policies remain valid. */
	public static void validatePolicy(final String[] args, final String family){
		NcrnaCallerProtection.validatePolicy(args);
		final String prefix=family.replace("_", "");
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			if(a.equals("config") || a.equals("ncrnafamily") || a.equals("ncrnafamilyroot")
					|| a.equals("ncrnakmers") || a.startsWith("ncrnaboundary")
					|| a.startsWith(prefix+"boundary") || a.equals(prefix+"root")
					|| a.equals(prefix+"consensus") || a.equals(prefix+"models") || a.equals(prefix+"kmers")
					|| NcrnaFamilyRegistry.isCobalamin(family) && a.equals("cobalaminroot")){
				throw new IllegalArgumentException("Boundary trials own selected-family resource and refinement options: "+a);
			}
		}
	}

	/** Reads the complete evaluation population shared by reduction and boundary trials.
	 * Paths must be readable and absolute, with no commas or surrounding whitespace:
	 * ordinary CallGenes splits each in= value on commas. Canonical-path uniqueness
	 * prevents evaluating one canonical path twice through alternate spellings.
	 * @param path One genome FASTA path per line; empty lines and # comments are skipped
	 * @return Owned list in the supplied order, preserving each accepted path spelling
	 */
	public static ArrayList<String> readGenomes(final String path){
		final ArrayList<String> genomes=new ArrayList<String>();
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				final String value=recordValue(new String(line, StandardCharsets.UTF_8));
				final File file=new File(value);
				if(!file.isAbsolute() || !file.isFile() || !file.canRead() || value.indexOf(',')>=0){
					throw new IllegalArgumentException("genomes= requires readable absolute FASTA paths without commas: "+value);
				}
				final String canonical;
				try{
					canonical=file.getCanonicalPath();
				}catch(IOException failure){
					throw new IllegalArgumentException("Could not resolve genome path: "+value, failure);
				}
				if(!seen.add(canonical)){throw new IllegalArgumentException("Duplicate evaluation genome: "+value);}
				genomes.add(value);
			}
		}
		if(genomes.isEmpty()){throw new IllegalArgumentException("Evaluation genome list is empty");}
		return genomes;
	}

	/** A native configuration has one trimmed flag=value record per line. */
	private static String recordValue(final String value){
		if(value==null || value.isEmpty() || !value.equals(value.trim()) || value.indexOf('\n')>=0 || value.indexOf('\r')>=0){
			throw new IllegalArgumentException("A caller path cannot contain line breaks or leading/trailing whitespace: "+value);
		}
		return value;
	}

	private final String family, displayFamily;
	private final String[] arguments;
	private final ArrayList<String> genomePaths;
	private final ArrayList<Locus> references;
}
