package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import align2.BBMap;
import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.FilterByCoverage;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Jul 8, 2015
 *
 */
public class Postfilter {

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		
		//Create a new CountKmersExact instance
		Postfilter x=new Postfilter(args, true);
		
		///And run it
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
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
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
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
		
		if(in2==null && in1!=null && in1.contains("#") && !new File(in1).exists()){
			in2=in1.replaceFirst("#", "2");
			in1=in1.replaceFirst("#", "1");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, covstats, out)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					covstats+", "+out+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2, ref)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, covstats, out, ref)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		assert(in1!=null);
		assert(out!=null);
		assert(ref!=null);
		assert(covstats!=null);
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
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
	
	private ArrayList<String> mapArgs=new ArrayList<String>();
	
	private String in1=null;
	private String in2=null;
	private String ref=null;
	private String out="filtered.fa";
	private String outdirty=null;
	private String covstats="covstats.txt";

	private int maxIndel=0;
	private int minHits=2;
	private int bw=20;
	private int tipsearch=0;
	private boolean rescue=false;
	
	private int trimEnds=0;
	
	private double minCoverage=2;
	private double minCoveredPercent=95;
	private long minReads=6;
	private int minLength=400;
	
	boolean showStats=true;
	
	boolean append=false;
	boolean overwrite=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	
}
