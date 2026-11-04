package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Proposes a smaller ncRNA consensus library and optionally protects caller detections.
 * Forked from template.A_SampleStreamerMT. The shared alignment worker stage is
 * NcrnaSupport; this driver owns CLI parsing, proposal selection and output.
 * Original resources are never changed. A proposal without full-genome calling
 * is explicitly labeled as such and is not a validated caller replacement.
 * @author Brian Bushnell, Raiden
 */
public final class ReduceNcrnaModels{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs reduction and optional caller protection, then closes the diagnostic stream. */
	public static void main(final String[] args){
		final Timer t=new Timer();
		final ReduceNcrnaModels tool=new ReduceNcrnaModels(args);
		try{
			tool.process(t);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Parses native key=value options and validates full-corpus unpaired mode and input/output routes.
	 * Failure after preprocessing closes the selected diagnostic log; successful
	 * construction leaves it open for processing and main's final cleanup.
	 * @param args Command-line arguments, including any config= files
	 */
	public ReduceNcrnaModels(String[] args){
		{
			final PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		// Main cannot own cleanup until this constructor returns successfully.
		try{
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			final Parser parser=parse(args);
			in=parser.in1;
			out=parser.out1;
			threads=Tools.max(1, parser.workers()>0 ? parser.workers() : Shared.threads());
			threadsIn=parser.threadsIn;
			if(parser.append || parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Reduction requires the complete unpaired corpus and a fresh output directory");
			}
			ReadWrite.setZipThreads(threads);
			validateParams();
			in=Tools.fixExtension(in);
			ref=Tools.fixExtension(ref);
			hbm=Tools.fixExtension(hbm);
			if(!Tools.testInputFiles(false, true, in, ref, hbm, genomes, truth, callerConfig)){
				throw new IllegalArgumentException("Cannot read a required reducer input");
			}
			ffin=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
			if(!ffin.fasta()){throw new IllegalArgumentException("in= must contain sense-oriented gene FASTA");}
		}catch(RuntimeException | Error failure){
			try{Shared.closeStream(outstream);}catch(RuntimeException | Error cleanup){
				if(cleanup!=failure){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
	}

	/** Parses reducer-specific options while retaining the standard BBTools option handling. */
	private Parser parse(final String[] args){
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("ref")){
				ref=b;
			}else if(a.equals("hbm")){
				hbm=b;
			}else if(a.equals("cutoff")){
				cutoff=NcrnaNumbers.floatValue(a, b);
				requestedCutoff=NcrnaNumbers.doubleValue(a, b);
			}else if(a.equals("target")){
				target=NcrnaNumbers.doubleValue(a, b);
			}else if(a.equals("genomes")){
				genomes=b;
			}else if(a.equals("truth")){
				truth=b;
			}else if(a.equals("family")){
				family=b;
			}else if(a.equals("truthfamily")){
				truthFamily=b;
			}else if(a.equals("callerconfig")){
				callerConfig=b;
			}else if(a.equals("callerheap")){
				callerHeap=b;
			}else if(a.equals("maxcalls")){
				maxCalls=NcrnaNumbers.intKMG(a, b);
			}else if(a.equals("regenerate")){
				regenerate=Parse.parseBoolean(b);
			}else if(a.equals("consensusconfig")){
				consensusConfig=b;
			}else if(a.equals("consensusheap")){
				consensusHeap=NcrnaNumbers.binaryValue(a, b)/(1024L*1024);
			}else if(a.equals("maxregencalls")){
				maxRegenCalls=NcrnaNumbers.intKMG(a, b);
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools options.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		return parser;
	}

	/** Requires a fresh destination and the complete option set when caller protection is requested. */
	private void validateParams(){
		if(in==null || ref==null || out==null || !(cutoff>=0 && cutoff<=1) || !(target>0 && target<=1)){
			throw new IllegalArgumentException("Require in=genes.fa ref=consensus.fa out=freshdir cutoff=[0,1] target=(0,1]");
		}
		if((genomes==null)!=(truth==null)){
			throw new IllegalArgumentException("genomes= and truth= must be supplied together for caller protection");
		}
		if(genomes!=null && (hbm==null || family==null || truthFamily==null || callerConfig==null || maxCalls<2)){
			throw new IllegalArgumentException("Caller protection also requires hbm=, family=, truthfamily=, callerconfig= and maxcalls>=2");
		}
		if(regenerate && (genomes==null || consensusConfig==null || consensusHeap<1 || maxRegenCalls<1)){
			throw new IllegalArgumentException("regenerate=t requires caller protection, consensusconfig, positive consensusheap and maxregencalls");
		}
		if(new File(out).exists()){throw new IllegalArgumentException("out= must name a fresh directory: "+out);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Scans the full corpus once, proposes a subset, optionally restores detections and writes results. */
	private void process(final Timer timer){
		final NcrnaModelRegeneration regeneration=regenerate ? new NcrnaModelRegeneration(consensusConfig,
			consensusHeap, maxRegenCalls, threads) : null;
		final NcrnaCallerProtection protection=genomes==null ? null : new NcrnaCallerProtection(callerConfig,
			genomes, truth, family, truthFamily, threads, callerHeap, maxCalls);
		final Library library=NcrnaSupport.loadLibrary(ref);
		final NcrnaSupport.ScanResult scan=NcrnaSupport.scanWithStats(ffin, library, cutoff, threads, threadsIn);
		final Graph original=scan.graph;
		final NcrnaModelSubset.Result proposal=NcrnaModelSubset.reduce(library.names.length, original.starts, original.models, target);
		final File directory=new File(out).getAbsoluteFile();
		if(!directory.mkdir()){throw new IllegalArgumentException("Could not create fresh output directory: "+directory);}
		NcrnaSupportIO.write(original, new File(directory, "original.support.tsv").getPath());
		boolean[] kept=proposal.kept;
		String status="Proposal only; caller protection not run";
		if(protection!=null){
			kept=protection.protect(library, ref, hbm, kept, directory);
			status=protection.status();
		}
		Graph reduced=original.subset(kept);
		if(regeneration==null){
			writeLibrary(library, hbm, kept, directory);
		}else{
			final NcrnaModelPairs baseline=new NcrnaModelPairs(library, hbm, kept);
			regeneration.process(in, baseline, reduced, protection, new File(directory, "regeneration"), outstream);
			regeneration.pairs.write(directory, modelFilename(hbm));
			reduced=regeneration.graph;
			status=status+"; "+regeneration.status;
		}
		NcrnaSupportIO.write(reduced, new File(directory, "members.tsv").getPath());
		writeSummary(directory, original, reduced, status, scan, regeneration==null ? 0 : regeneration.alignments());
		timer.stop();
		outstream.println("Models retained: "+reduced.modelNames.length+" of "+original.modelNames.length);
		outstream.println("Result: "+status);
		outstream.println("Time: "+timer);
	}

	/** Exports a paired subset with original headers and HBM count data, and upper-case sequences.
	 * Each HBM is checked against the same-position consensus before selection.
	 * This is a library subset, not an installable family resource directory.
	 */
	static void writeLibrary(final Library library, final String hbm,
			final boolean[] kept, final File directory){
		new NcrnaModelPairs(library, hbm, kept).write(directory, hbm==null ? null : modelFilename(hbm));
	}

	/** Preserves the source HBM encoding and gzip suffix in the output filename. */
	static String modelFilename(final String hbm){
		final String suffix=NucleotideHbmCodec.isCompact(hbm) ? ".hbmc" : ".hbm";
		return "models"+suffix+(hbm.endsWith(".gz") ? ".gz" : "");
	}

	/** Reports full-corpus coverage and measured work without implying genomic false-positive validation. */
	private void writeSummary(final File directory, final Graph original, final Graph reduced,
			final String status, final NcrnaSupport.ScanResult scan, final long regenerationAlignments){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("status\t").append(status).nl();
		bb.append("models_original\t").append(original.modelNames.length).nl();
		bb.append("models_retained\t").append(reduced.modelNames.length).nl();
		bb.append("genes\t").append(original.queryNames.length).nl();
		bb.append("coverage_original\t").append(coverage(original), 8, true).nl();
		bb.append("coverage_retained\t").append(coverage(reduced), 8, true).nl();
		bb.append("identity_cutoff\t").append(requestedCutoff, 8, true).nl();
		bb.append("coverage_target\t").append(target, 8, true).nl();
		bb.append("pair_alignments\t").append(scan.alignments+regenerationAlignments).nl();
		bb.append("initial_pair_alignments\t").append(scan.alignments).nl();
		bb.append("regeneration_pair_alignments\t").append(regenerationAlignments).nl();
		bb.append("input_batches\t").append(scan.batches).nl();
		bb.append("active_workers\t").append(scan.activeWorkers).nl();
		bb.append("false_positive_rate\tUNMEASURED\n");
		bb.append("method\tGreedy subset proposal; no global minimum claim.\n");
		if(genomes==null){bb.append("limitation\tNo full-genome caller protection was performed.\n");}
		try(Output output=new Output(new File(directory, "summary.txt").getPath())){
			output.writer.print(bb);
		}
	}

	/** Fraction of all input genes having at least one retained model meeting the alignment cutoff. */
	private static double coverage(final Graph graph){
		assert(graph.queryNames.length>0) : "Coverage denominator includes the complete nonempty query population; genes="+graph.queryNames.length;
		int covered=0;
		for(int gene=0; gene<graph.queryNames.length; gene++){
			if(graph.starts[gene+1]>graph.starts[gene]){covered++;}
		}
		return covered/(double)graph.queryNames.length;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Gene FASTA, consensus FASTA, optional paired HBM, and fresh output directory. */
	private String in, ref, hbm, out;
	/** Optional protection inputs: genome list, truth table, caller/truth labels and scientific config. */
	private String genomes, truth, family, truthFamily, callerConfig;
	/** Native JVM heap value for each sequential caller child, such as 2g. */
	private String callerHeap="2g";
	/** Minimum native float alignment identity for a gene-to-model support edge. */
	private float cutoff=Float.NaN;
	/** Original double value for human reports; alignments retain native float semantics. */
	private double requestedCutoff=Double.NaN;
	/** Required supported fraction of ALL input genes, including initially unsupported genes. */
	private double target=0.98;
	/** Alignment-worker budget shared by support and caller stages. */
	private int threads=1;
	/** Native input-reader worker setting, separate from alignment workers. */
	private int threadsIn=1;
	/** Initial caller-protection budget, including its baseline call. */
	private int maxCalls=4;
	/** Whether to rebuild retained models after the initial caller protection passes. */
	private boolean regenerate=false;
	/** Scientific consensus-builder recipe used only for regeneration. */
	private String consensusConfig;
	/** Consensus-builder child heap budget in MiB. */
	private long consensusHeap=1024;
	/** Additional regeneration caller budget; the protected baseline is reused. */
	private int maxRegenCalls=5;
	/** Validated complete unpaired gene FASTA input descriptor. */
	private FileFormat ffin;
	/** Selected diagnostic stream; constructor failure or main closes it. */
	private PrintStream outstream=System.err;
}
