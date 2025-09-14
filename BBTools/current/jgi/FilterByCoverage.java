package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date October 8, 2014
 *
 */
public class FilterByCoverage {

	public static void main(String[] args){
		Timer t=new Timer();
		FilterByCoverage x=new FilterByCoverage(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FilterByCoverage(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		Shared.capBufferLen(20);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("coverage") || a.equals("cov") || a.equals("covstats") || a.equals("coverage1") || a.equals("cov1") || a.equals("covstats1")){
				covStatsAfter=b;
			}else if(a.equals("coverage0") || a.equals("cov0") || a.equals("covstats0")){
				covStatsBefore=b;
			}else if(a.equals("minc") || a.equals("mincov") || a.equals("mincoverage")){
				minCoverage=Double.parseDouble(b);
			}else if(a.equals("minp") || a.equals("minpercent")){
				minCoveredPercent=Double.parseDouble(b);
			}else if(a.equals("minr") || a.equals("minreads")){
				minReads=Parse.parseKMG(b);
			}else if(a.equals("minratio") || a.equals("ratio")){
				minRatio=Double.parseDouble(b);
			}else if(a.equals("basesundermin")){
				basesUnderMin=Integer.parseInt(b);
			}else if(a.equals("minl") || a.equals("minlen") || a.equals("minlength")){
				minLength=Integer.parseInt(b);
			}else if(a.equals("trim") || a.equals("trimends")){
				if(b==null || Character.isLetter(b.charAt(0))){
					trimEnds=Parse.parseBoolean(b) ? 100 : 0;
				}else{
					trimEnds=Integer.parseInt(b);
				}
				trimEnds=Tools.max(trimEnds, 0);
			}else if(a.equals("appendresults") || a.equals("logappend") || a.equals("appendlog") || a.equals("appendtolog")){
				logappend=Parse.parseBoolean(b);
			}else if(a.equals("log") || a.equals("results")){
				logfile=b;
			}else if(a.equals("logheader")){
				logheader=Parse.parseBoolean(b);
			}else if(a.equals("outd") || a.equals("outdirty")){
				outdirty=b;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					parser.in1=arg.replace("#", "1");
					parser.in2=arg.replace("#", "2");
				}
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			if(parser.minReadLength>0){minLength=parser.minReadLength;}
			
			in1=parser.in1;
			qfin1=parser.qfin1;

			outclean=parser.out1;
			qfoutclean=parser.qfout1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		minLength=Tools.max(1, minLength);
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		name=ReadWrite.stripToCore(in1);
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(outclean!=null && outclean.equalsIgnoreCase("null")){outclean=null;}
		if(outdirty!=null && outdirty.equalsIgnoreCase("null")){outdirty=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outclean, outdirty)){
			outstream.println((outclean==null)+", "+outclean+", "+(outdirty==null)+", "+outdirty);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outclean+", "+outdirty+"\n");
		}

		ffoutclean=FileFormat.testOutput(outclean, FileFormat.FASTA, extout, true, overwrite, append, false);
		ffoutdirty=FileFormat.testOutput(outdirty, FileFormat.FASTA, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffCov0=FileFormat.testInput(covStatsBefore, FileFormat.TEXT, ".txt", true, false);
		ffCov1=FileFormat.testInput(covStatsAfter, FileFormat.TEXT, ".txt", true, false);
		
		assert(covStatsAfter!=null) : "No coverage file specified.";
	}
	
	void process(Timer t){

		final HashMap<String, CovStatsLine> cslMap0=new HashMap<String, CovStatsLine>(1024);
		final HashMap<String, CovStatsLine> cslMap1=new HashMap<String, CovStatsLine>(1024);
		if(ffCov0!=null){
			TextFile tf=new TextFile(ffCov0);
			int i=0;
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(i==0){
					assert(s.startsWith("#")) : "Expected a header line starting with #";
					CovStatsLine.initializeHeader(s);
				}else{
					CovStatsLine csl=new CovStatsLine(s);
					CovStatsLine old=cslMap0.put(csl.id, csl);
					assert(old==null);
				}
				i++;
			}
			tf.close();
		}
		if(ffCov1!=null){
			TextFile tf=new TextFile(ffCov1);
			int i=0;
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(i==0){
					assert(s.startsWith("#")) : "Expected a header line starting with #";
					CovStatsLine.initializeHeader(s);
				}else{
					CovStatsLine csl=new CovStatsLine(s);
					CovStatsLine old=cslMap1.put(csl.id, csl);
					assert(old==null);
				}
				i++;
			}
			tf.close();
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		assert(!cris.paired());
		
		final ConcurrentReadOutputStream rosClean;
		if(outclean!=null){
			final int buff=4;
			
			assert(!outclean.equalsIgnoreCase(in1) && !outclean.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosClean=ConcurrentReadOutputStream.getStream(ffoutclean, null, qfoutclean, null, buff, null, false);
			rosClean.start();
		}else{rosClean=null;}
		
		final ConcurrentReadOutputStream rosDirty;
		if(outdirty!=null){
			final int buff=4;
			
			assert(!outdirty.equalsIgnoreCase(in1) && !outdirty.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosDirty=ConcurrentReadOutputStream.getStream(ffoutdirty, null, qfoutdirty, null, buff, null, false);
			rosDirty.start();
		}else{rosDirty=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;

		long basesTrimmed=0;
		
		long readsOut=0;
		long basesOut=0;
		
		long readsFiltered=0;
		long basesFiltered=0;
		
		final TextStreamWriter tsw=(logfile==null ? null : new TextStreamWriter(logfile, (overwrite && !logappend), logappend, true));
//		System.err.println("***** overwrite="+overwrite+", logappend="+logappend+", combined="+(overwrite && !logappend));
		if(tsw!=null){
			tsw.start();
			if(logheader){tsw.print("#assembly\tcontig\tcontam\tlength\tavgFold\treads\tpercentCovered"+(ffCov0==null ? "" : "\tavgFold0\treads0\tnormRatio")+"\n");}
		}
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				final ArrayList<Read> cleanList=new ArrayList<Read>(reads.size());
				final ArrayList<Read> dirtyList=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					assert(r1.mate==null);
					
					final int initialLength1=r1.length();
					
					readsProcessed++;
					basesProcessed+=initialLength1;
					
					if(trimEnds>0){
						if(initialLength1-trimEnds*2<minLength){
							r1.bases=r1.quality=null;
						}else{
							TrimRead.trimByAmount(r1, trimEnds, trimEnds, 0);
						}
					}
					final int length=r1.length();
					basesTrimmed+=(initialLength1-length);
					
					final double covRatio;
					final boolean contam;
					
					final CovStatsLine csl0=cslMap0.get(r1.id);
					final CovStatsLine csl1=cslMap1.get(r1.id);
					if(csl1!=null){
						
						if(csl0!=null){
							covRatio=csl0.avgFold/Tools.max(0.01, csl1.avgFold);
							int underMin=csl0.underMin-csl1.underMin;
							
							if(csl1.reads()<minReads || length<minLength || csl1.coveredPercent()<minCoveredPercent){
								contam=true;
							}else if((csl1.avgFold<minCoverage && covRatio>minRatio) || csl1.avgFold<0.5){
								contam=true;
							}else if(basesUnderMin>0 && underMin>basesUnderMin){
								contam=true;
							}else{
								contam=false;
							}
						}else{
							covRatio=0;
							int underMin=csl1.underMin;
							
							if(csl1.reads()<minReads || length<minLength || csl1.coveredPercent()<minCoveredPercent || csl1.avgFold<minCoverage){
								contam=true;
							}else if(basesUnderMin>0 && underMin>basesUnderMin){
								contam=true;
							}else{
								contam=false;
							}
						}
						
					}else{
						contam=true;
						covRatio=0;
					}
					
					if(!contam){
						cleanList.add(r1);
						readsOut++;
						basesOut+=length;
					}else{
						dirtyList.add(r1);
						readsFiltered++;
						basesFiltered+=length;
					}
					if(tsw!=null && (length>=minLength || PRINT_SHORT_CONTIG_RESULTS)){
						if(csl1==null){
							if(ffCov0==null){
								tsw.print(Tools.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\n", name, r1.id, contam ? "1" : "0", length, 0, 0, 0));
							}else{
								tsw.print(Tools.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\n",
										name, r1.id, contam ? "1" : "0", length, 0, 0, 0, 0, 0, 0));
							}
							
						}else if(csl0==null){
							tsw.print(Tools.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\n", name, csl1.id, contam ? "1" : "0", length,
									csl1.avgFold, csl1.plusReads+csl1.minusReads, csl1.coveredPercent()));
						}else{
							tsw.print(Tools.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\n", name, csl1.id, contam ? "1" : "0", length,
									csl1.avgFold, csl1.plusReads+csl1.minusReads, csl1.coveredPercent(), csl0.avgFold, csl0.plusReads+csl0.minusReads, covRatio));
						}
					}
				}
				
				if(rosClean!=null){rosClean.add(cleanList, ln.id);}
				if(rosDirty!=null){rosDirty.add(dirtyList, ln.id);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, rosClean, rosDirty);
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:               "+t);
		outstream.println("Reads In:           "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases In:           "+basesProcessed+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		outstream.println("Reads Filtered:     "+readsFiltered);
		outstream.println("Bases Filtered:     "+basesFiltered);
		if(trimEnds>0){
			outstream.println("Bases Trimmed:      "+basesTrimmed);
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String covStatsBefore=null;
	private String covStatsAfter=null;
	private String name=null;
	
	private String qfin1=null;

	private String outclean=null;
	private String outdirty=null;

	private String qfoutclean=null;
	private String qfoutdirty=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	private long maxReads=-1;

	/** Scaffolds shorter than this will be discarded. */
	private int minLength=0;
	/** Scaffolds with fewer mapped reads will be discarded. */
	private long minReads=0;
	/** Scaffolds with lower average coverage will be discarded. */
	private double minCoverage=5;
	/** Scaffolds with a lower percent of covered bases will be discarded. */
	private double minCoveredPercent=40;
	/** Scaffolds will NOT be discarded based on low coverage unless the coverage dropped by at least this factor. */
	private double minRatio=0;
	/** Scaffolds will be discarded if there are at least this many bases in windows below a coverage cutoff. */
	private int basesUnderMin=-1;
	
	/** Trim this much from sequence ends */
	private int trimEnds=0;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffCov0;
	private final FileFormat ffCov1;

	private final FileFormat ffoutclean;
	private final FileFormat ffoutdirty;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean logappend=false;
	private String logfile=null;
	private boolean logheader=true;
	private static boolean PRINT_SHORT_CONTIG_RESULTS=false;
	
}
