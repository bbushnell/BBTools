package jgi;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;

/**
 * Filters root reads by clipping bounds, total substitutions plus at least one
 * substitution in an inclusive quality-score interval unless minSubs is below
 * one, optional counted indel evidence, or optional zero substitutions and zero
 * counted indels. Used for manually examining specific reads that may have
 * incorrectly calibrated quality scores.
 *
 * @author Brian Bushnell
 * @date May 5, 2015
 */
public class FilterReadsWithSubs{

	/** Program entry point: constructs FilterReadsWithSubs, runs processing, reports timing.
	 * @param args Command-line arguments for input/output and filter settings */
	public static void main(String[] args){
		Timer t=new Timer();
		
		FilterReadsWithSubs x=new FilterReadsWithSubs(args);
		
		x.process(t);
		
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command-line arguments, configures substitution/quality/clipping thresholds,
	 * and initializes SAM input/output file formats for filtering.
	 * @param args Command-line arguments for filter configuration and I/O
	 */
	public FilterReadsWithSubs(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		SamLine.SET_FROM_OK=true;
		minq=0;
		maxq=99;
		minSubs=1;
		countIndels=true;
		keepPerfect=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//do nothing
			}else if(a.equals("minq")){
				minq=Parse.parseIntKMG(b);
			}else if(a.equals("maxq")){
				maxq=Parse.parseIntKMG(b);
			}else if(a.equals("keepperfect")){
				keepPerfect=Parse.parseBoolean(b);
			}else if(a.equals("countindels")){
				countIndels=Parse.parseBoolean(b);
			}else if(a.equals("minsubs")){
				minSubs=Integer.parseInt(b);
			}else if(a.equals("clips") || a.equalsIgnoreCase("maxclips") || a.equalsIgnoreCase("maxclip")){
				maxClips=Integer.parseInt(b);
				if(maxClips<0){maxClips=Integer.MAX_VALUE;}
			}else if(a.equalsIgnoreCase("minclips") || a.equalsIgnoreCase("minclip")){
				minClips=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
	}
	
	/**
	 * Streams root reads, applies filtering to each list root, writes retained
	 * roots, and reports throughput statistics counted with pairCount() and
	 * pairLength(). Mates contribute to counters and output as part of their root
	 * read; they are not evaluated separately.
	 * @param t Timer for tracking execution
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}
		boolean paired=cris.paired();

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, ffout1.samOrBam());
			ros.start();
		}else{ros=null;}
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					
					boolean keep=processRead(r1);
					if(keep){
						readsOut+=r1.pairCount();
						basesOut+=r1.pairLength();
					}else{reads.set(idx, null);}
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		boolean errorState=ReadWrite.closeStreams(cris, ros);
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		if(verbose){outstream.println("Finished.");}
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
	}
	
	/**
	 * Evaluates a root read against clipping and three alternate retention paths:
	 * total substitutions at least minSubs with at least one substitution in the
	 * inclusive minq/maxq interval unless minSubs is below one, counted indel
	 * evidence when countIndels is enabled, or zero substitutions and zero counted
	 * indels when keepPerfect is enabled. The indel and keepPerfect alternatives
	 * are not gated by minSubs or by the minq/maxq interval.
	 * @param r1 Read to evaluate
	 * @return true if mandatory checks pass and any retention route accepts the read
	 */
	private boolean processRead(Read r1){
		final byte[] quals=r1.quality, bases=r1.bases;
		final byte[] match=(r1.match==null ? null : !r1.shortmatch() ? r1.match : Read.toLongMatchString(r1.match));
		if(match==null || quals==null || bases==null){return false;}
		
		int subs=0, passingSubs=0;
		int indels=0;
		int clips=0;
		for(int qpos=0, mpos=0; mpos<match.length; mpos++){
			
			final byte m=match[mpos];
			final byte mprev=match[Tools.max(mpos-1, 0)];
			final byte mnext=match[Tools.min(mpos+1, match.length-1)];
			if(m=='D' && qpos==bases.length && qpos==quals.length){
				// Terminal deletions consume no query base after the read is exhausted.
				continue;
			}
			
			final byte q=quals[qpos];
			final byte b=bases[qpos];
			
			if(m=='S'){
				subs++;
				if(q>=minq && q<=maxq){passingSubs++;}
			}else if(m=='I'){
				indels++;
			}else if(m=='m'){
				if(mprev=='D' || mnext=='D'){
					indels++;
				}
			}else if(m=='D'){
				//do nothing
			}else if(m=='N'){
				//do nothing
			}else if(m=='C'){
				clips++;
			}else{
				throw new RuntimeException("Bad symbol m='"+((char)m)+"'\n"+new String(match)+"\n"+new String(bases)+"\n");
			}
			if(clips>maxClips){return false;}
			
			if(m!='D'){qpos++;}
		}
		if(clips>maxClips || clips<minClips){return false;}
//		System.err.println("match="+new String(match));
//		System.err.println("clips="+clips);
//		assert(false) : clips+", "+subs+", "+indels+", "+maxClips+", "+minClips;
		if(subs>=minSubs && (passingSubs>0 || minSubs<1)){return true;}
		if(countIndels && indels>0){return true;}
		return keepPerfect && subs==0 && indels==0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           I/O Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------       Tracking Fields        ----------------*/
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long readsProcessed=0;
	private long readsOut=0;
	private long basesProcessed=0;
	private long basesOut=0;

	/**
	 * Substitution-route minimum and maximum quality scores and minimum total
	 * substitutions; ignored by indel and keepPerfect retention.
	 */
	public int minq, maxq, minSubs;
	/** Maximum number of clipped bases allowed */
	public int maxClips=Integer.MAX_VALUE;
	/** Minimum number of clipped bases required */
	public int minClips=0;
	/** Whether reads with insertion/deletion evidence can be retained */
	public boolean countIndels;
	/** Whether reads with no substitutions or counted indels can be retained */
	public boolean keepPerfect;
	
	/*--------------------------------------------------------------*/
	/*----------------       Output Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
