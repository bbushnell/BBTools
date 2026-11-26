package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.GiToTaxid;
import tracker.ReadStats;

/**
 * Keeps the best copy of an SSU per organism based on quality metrics.
 * Processes sequence reads and maintains only the highest quality copy for each
 * taxonomic ID, comparing reads by defined sequence length and ambiguous base count.
 * Selection, accounting, and {@code maxlen} checks use each list root; mates
 * remain attached to retained roots but are not ranked independently.
 *
 * @author Brian Bushnell
 * @date Oct 4, 2019
 */
public class KeepBestCopy{
	
	/**
	 * Parses arguments, runs the filter, reports throughput, and closes the
	 * diagnostic stream if it was redirected.
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		KeepBestCopy x=new KeepBestCopy(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Applies command-line configuration, checks the required input name and
	 * output availability, and records file formats for read and write streams.
	 * @param args Command-line arguments
	 */
	public KeepBestCopy(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		Shared.capBuffers(4); //Only for singlethreaded programs
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("maxlen")){
				maxLen=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}
			
			else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
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
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			qfin1=parser.qfin1;
			
			out1=parser.out1;
			qfout1=parser.qfout1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}

		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
	}
	
	/**
	 * Opens and starts the configured primary input stream.
	 * The stream may report paired input, but downstream selection uses
	 * only the root reads in each returned list.
	 */
	ConcurrentReadInputStream makeCris(){
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
		cris.start();
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/**
	 * Consumes roots returned by the configured stream, up to the optional
	 * maximum read limit, and retains at most one root per taxonomic ID.
	 * Read counters count roots; base counters and {@code maxlen} comparisons
	 * use root lengths only. Mates remain attached to retained roots but are
	 * not ranked independently. If output is configured, retained roots are
	 * written in first-seen taxonomic insertion order.
	 * @param t Timer to stop and report after streams are closed
	 */
	void process(Timer t){
		final ConcurrentReadInputStream cris=makeCris();
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		long readsProcessed=0, readsOut=0;
		long basesProcessed=0, basesOut=0;
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				final ArrayList<Read> listOut=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					final boolean keep=process(r1);
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		errorState|=ReadWrite.closeStreams(cris);

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, qfout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		{
			ArrayList<Read> list=new ArrayList<Read>(200);
			long ln=0;
			for(Entry<Integer, Read> e : map.entrySet()){
				Read r=e.getValue();
				list.add(r);
				readsOut++;
				basesOut+=r.length();
				if(list.size()>=200){
					if(ros!=null){ros.add(list, ln);}
					ln++;
					list=new ArrayList<Read>(200);
				}
			}
			if(!list.isEmpty()){
				if(ros!=null){ros.add(list, ln);}
			}
		}
		
		errorState|=ReadStats.writeAll();
		
		errorState|=ReadWrite.closeStream(ros);
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses the root read ID for a taxonomic key, skips negative or
	 * unrecognized TaxIDs, and retains the root if no earlier root used that
	 * key or this root outranks the incumbent.
	 * Replacing a value does not move its key in the insertion-ordered map.
	 * @param r Root read to classify and possibly retain
	 * @return {@code true} if this root became the retained root for its key
	 */
	private boolean process(Read r){
		int tid=GiToTaxid.parseTaxidNumber(r.id, '|');
		if(tid<0){return false;}
		Integer key=tid;
		Read old=map.get(key);
		if(old==null || isBetterThan(r, old)){
			map.put(key, r);
			return true;
		}
		return false;
	}
	
	/**
	 * Ranks candidate roots by length and ambiguity.
	 * A null incumbent is always replaced. If lengths differ and either
	 * exceeds {@code maxLen}, the shorter root wins; otherwise, the candidate
	 * must have more fully defined bases or an equal number of fully defined
	 * bases with fewer undefined bases. Exact metric ties keep the incumbent.
	 * @param r Candidate root
	 * @param old Incumbent root, or {@code null}
	 * @return {@code true} if the candidate should replace the incumbent
	 */
	private boolean isBetterThan(Read r, Read old){
		if(old==null){return true;}
		int oldNs=old.countUndefined();
		int Ns=r.countUndefined();
		int oldDef=old.length()-oldNs;
		int def=r.length()-Ns;
		if(old.length()>maxLen && r.length()<old.length()){return true;}
		if(r.length()>maxLen && old.length()<r.length()){return false;}
		return def>oldDef || (def==oldDef && Ns<oldNs);
	}
	
	/*--------------------------------------------------------------*/
	
	/** Input sequence filename. */
	private String in1=null;
	
	/** Optional input quality filename. */
	private String qfin1=null;

	/** Optional output sequence filename. */
	private String out1=null;

	/** Optional output quality filename. */
	private String qfout1=null;
	
	/** Optional forced input file extension. */
	private String extin=null;
	
	/** Optional forced output file extension. */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Length above which shorter unequal roots are preferred. */
	int maxLen=1600;
	
	/** Maximum root reads to consume from the input stream. */
	private long maxReads=-1;
	
	/** Retained root read per taxonomic ID, in first-seen TaxID order. */
	private LinkedHashMap<Integer, Read> map=new LinkedHashMap<Integer, Read>();
	
	/*--------------------------------------------------------------*/
	
	/** Parsed primary input format. */
	private final FileFormat ffin1;
	
	/** Parsed primary output format, or null when output is disabled. */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	/** Diagnostic stream, possibly redirected by the preparser. */
	private PrintStream outstream=System.err;
	
	/** Enables diagnostic output in this class and selected readers. */
	public static boolean verbose=false;
	
	/** Set when a stream or statistics writer reports a failure. */
	public boolean errorState=false;
	
	/** Allows existing output files to be replaced. */
	private boolean overwrite=true;
	
	/** Appends to existing output files instead of truncating them. */
	private boolean append=false;
	
}
