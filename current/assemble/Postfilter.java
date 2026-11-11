package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import align2.BBMap;
import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.FilterByCoverage;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;

/**
 * Maps reads to an assembly with BBMap, then filters contigs with FilterByCoverage.
 * Coverage statistics are generated before optional end trimming. Passing the
 * configured coverage thresholds is not a guarantee of assembly correctness.
 * Instances are intended for one sequential run; child tools and the diagnostic
 * stream use shared static state, so this is not an isolated concurrent pipeline.
 *
 * @author Brian Bushnell
 * @date Jul 8, 2015
 */
public class Postfilter{

	/**
	 * Program entry point that creates a Postfilter instance and executes the
	 * two-phase filtering pipeline.
	 * @param args Command line arguments containing input/output paths and filtering parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		
		//Create and run the mapping/filtering pipeline.
		Postfilter x=new Postfilter(args, true);
		
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructs a Postfilter with command-line arguments and optional default configuration.
	 * Parses all filtering parameters, validates input/output files, and configures
	 * BBMap and FilterByCoverage settings.
	 * Retains equals signs inside values. Checks rejected-output collisions and
	 * overwrite refusal before mapping can create coverage output.
	 * Rejects effective NaN coverage thresholds after all arguments have been parsed.
	 * Normalizes option names independently of locale without changing their values.
	 *
	 * @param args Command line arguments for filtering configuration
	 * @param setDefaults Whether to apply default compression and threading settings
	 */
	public Postfilter(String[] args, boolean setDefaults){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		if(setDefaults){
			/* Set global defaults */
			ReadWrite.ZIPLEVEL=8;
			ReadWrite.USE_UNPIGZ=true;
			ReadWrite.USE_PIGZ=true;
			if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
				ByteFile.FORCE_MODE_BF2=true;
			}
		}
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("ref") || a.equals("contigs") || a.equals("assembly")){
				ref=b;
			}else if(a.equals("out") || a.equals("out1")){
				out=b;
			}else if(a.equals("outdirty") || a.equals("outd") || a.equals("outbad")){
				outdirty=b;
			}else if(a.equals("showstats")){
				showStats=Parse.parseBoolean(b);
			}else if(a.equals("covstats") || a.equals("cov")){
				covstats=b;
			}else if(a.equals("maxindel")){
				maxIndel=Integer.parseInt(b);
			}else if(a.equals("minhits")){
				minHits=Integer.parseInt(b);
			}else if(a.equals("minc") || a.equals("mincov") || a.equals("mincoverage")){
				minCoverage=Double.parseDouble(b);
			}else if(a.equals("minp") || a.equals("minpercent")){
				minCoveredPercent=Double.parseDouble(b);
			}else if(a.equals("minr") || a.equals("minreads")){
				minReads=Parse.parseKMG(b);
			}else if(a.equals("minl") || a.equals("minlen") || a.equals("minlength")){
				minLength=Integer.parseInt(b);
			}else if(a.equals("rescue")){
				rescue=Parse.parseBoolean(b);
			}else if(a.equals("trim") || a.equals("trimends")){
				if(b==null || Character.isLetter(b.charAt(0))){
					trimEnds=Parse.parseBoolean(b) ? 100 : 0;
				}else{
					trimEnds=Integer.parseInt(b);
				}
				trimEnds=Tools.max(trimEnds, 0);
			}else{
				mapArgs.add(arg);
			}
		}
		
		if(Double.isNaN(minCoverage) || Double.isNaN(minCoveredPercent)){
			throw new IllegalArgumentException("mincov and minpercent must not be NaN.");
		}
		if(in2==null && in1!=null && in1.contains("#") && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		
		// Check all path collisions before output validation can assert on duplicates.
		if(!Tools.testForDuplicateFiles(true, in1, in2, covstats, out, ref, outdirty)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		// The optional rejected-output stream opens only after mapping; fail before that work.
		if(!overwrite && !append && Tools.isOutputFileName(outdirty) && new File(outdirty).exists()){
			throw new RuntimeException("File "+outdirty+" exists and overwrite=false.");
		}
		if(!Tools.testOutputFiles(overwrite, append, false, covstats, out, outdirty)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					covstats+", "+out+", "+outdirty+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		assert(in1!=null);
		assert(out!=null);
		assert(ref!=null);
		assert(covstats!=null);
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Runs BBMap, unloads its reference, and invokes FilterByCoverage on the original
	 * contigs with the generated coverage file. Unrecognized constructor arguments
	 * are forwarded only to BBMap; wrapper mapping options are appended afterward.
	 * Reuses the argument list for the filtering call, so this method is intended
	 * to run once per instance. Optionally rereads non-stdio output for assembly stats.
	 *
	 * @param t Legacy timer parameter; not read or stopped by this method
	 */
	public void process(Timer t){

//		bbmap.sh in=../reads.fq.gz ref=contigs.fasta nodisk ambig=all maxindel=100 covstats=covstats.txt minhits=2;
//		filterbycoverage.sh in=contigs.fasta out=filtered.fasta cov=covstats.txt mincov=2 minr=6 minp=95 minl=400

		mapArgs.add("in="+in1);
		if(in2!=null){mapArgs.add("in2="+in2);}
		mapArgs.add("ref="+ref);
		mapArgs.add("covstats="+covstats);
		mapArgs.add("ambig=all");
		mapArgs.add("minhits="+minHits);
		mapArgs.add("maxindel="+maxIndel);
		mapArgs.add("nodisk");
		mapArgs.add("append="+append);
		mapArgs.add("ow="+overwrite);
		mapArgs.add("bw="+bw);
		mapArgs.add("tipsearch="+tipsearch);
		mapArgs.add("rescue="+rescue);
		BBMap.main(mapArgs.toArray(new String[0]));
		Data.unloadAll();
		
		mapArgs.clear();
		mapArgs.add("in="+ref);
		mapArgs.add("out="+out);
		if(outdirty!=null){mapArgs.add("outdirty="+outdirty);}
		mapArgs.add("covstats="+covstats);
		mapArgs.add("mincov="+minCoverage);
		mapArgs.add("minr="+minReads);
		mapArgs.add("minp="+minCoveredPercent);
		mapArgs.add("minl="+minLength);
		mapArgs.add("trim="+trimEnds);
		mapArgs.add("append="+append);
		mapArgs.add("ow="+overwrite);
		FilterByCoverage.main(mapArgs.toArray(new String[0]));
		
		if(showStats && out!=null && !FileFormat.isStdio(out)){
			outstream.println();
			jgi.AssemblyStats2.main(new String[] {"in="+out});
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Forwarded mapping arguments, then reused for the filtering invocation. */
	private ArrayList<String> mapArgs=new ArrayList<String>();
	
	/** Primary reads to map, with an optional mate-placeholder expansion. */
	private String in1=null;
	/** Optional separate mate-read file. */
	private String in2=null;
	/** Assembly used as mapping reference and subsequent filtering input. */
	private String ref=null;
	/** Destination for retained contigs. */
	private String out="filtered.fa";
	/** Optional destination for rejected contigs. */
	private String outdirty=null;
	/** Coverage table produced by mapping and consumed by filtering. */
	private String covstats="covstats.txt";

	/** Maximum indel setting passed to BBMap. */
	private int maxIndel=0;
	/** Minimum seed-hit setting passed to BBMap. */
	private int minHits=2;
	/** Fixed mapping bandwidth appended after forwarded arguments. */
	private int bw=20;
	/** Fixed mapping tip-search setting appended after forwarded arguments. */
	private int tipsearch=0;
	/** Whether BBMap may attempt mate rescue. */
	private boolean rescue=false;
	
	/** Bases trimmed from each contig end during filtering. */
	private int trimEnds=0;
	
	/** Minimum average mapped depth required by the coverage filter. */
	private double minCoverage=2;
	/** Minimum covered-base percentage required by the coverage filter. */
	private double minCoveredPercent=95;
	/** Minimum mapped-read count required by the coverage filter. */
	private long minReads=6;
	/** Minimum contig length after optional end trimming. */
	private int minLength=400;
	
	/** Whether to reread a non-stdio output for assembly statistics. */
	boolean showStats=true;
	
	/** Append policy forwarded to mapping and filtering outputs. */
	boolean append=false;
	/** Overwrite policy checked before mapping and forwarded to both child tools. */
	boolean overwrite=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Shared diagnostic stream selected by the most recent constructor. */
	private static PrintStream outstream=System.err;
	
}
