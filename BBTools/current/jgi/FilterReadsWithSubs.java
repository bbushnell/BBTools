package jgi;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;

/**
 * Filters to select only reads with substitution errors
 * for bases with quality scores in a certain interval.
 * Used for manually examining specific reads that may have
 * incorrectly calibrated quality scores, for example.
 * @author Brian Bushnell
 * @date May 5, 2015
 *
 */
public class FilterReadsWithSubs {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		FilterReadsWithSubs x=new FilterReadsWithSubs(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
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
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here

//				System.err.println("Calling parseArgument("+arg+","+a+","+b+")");
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
				if(maxClips<0) {maxClips=Integer.MAX_VALUE;}
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
					if(keep) {
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
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
	}
	
	private boolean processRead(Read r1) {
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
			if(clips>maxClips) {return false;}
			
			if(m!='D'){qpos++;}
		}
		if(clips>maxClips || clips<minClips) {return false;}
//		System.err.println("match="+new String(match));
//		System.err.println("clips="+clips);
//		assert(false) : clips+", "+subs+", "+indels+", "+maxClips+", "+minClips;
		if(subs>=minSubs && (passingSubs>0 || minSubs<1)) {return true;}
		if(countIndels && indels>0) {return true;}
		return keepPerfect && subs==0 && indels==0;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long readsProcessed=0;
	private long readsOut=0;
	private long basesProcessed=0;
	private long basesOut=0;

	public int minq, maxq, minSubs;
	public int maxClips=Integer.MAX_VALUE;
	public int minClips=0;
	public boolean countIndels;
	public boolean keepPerfect;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
