package jgi;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.KmerTracker;

/**
 * @author Brian Bushnell
 * @date Oct 2, 2025
 *
 */
public class Scalars {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Scalars x=new Scalars(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public Scalars(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("header") || a.equals("colheader") || a.equals("columnheader")){
				header=Parse.parseBoolean(b);
			}else if(a.equals("rowheader")){
				rowheader=Parse.parseBoolean(b);
			}else if(a.equals("window")){
				window=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		dimers=new KmerTracker(2, window);
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		long readsProcessed=0, basesProcessed=0;
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
					
					if(window<1) {
						dimers.add(r1.bases);
						if(r1.mate!=null) {dimers.add(r1.mate.bases);}
					}else {
						addWindowed(r1.bases);
						if(r1.mate!=null) {addWindowed(r1.mate.bases);}
					}
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}
		
		outputResults();
		
		t.stop();
		if(printTime) {
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+readsProcessed+" \t"+
				Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		}
		assert(!errorState) : "An error was encountered.";
	}
	
	void addWindowed(byte[] bases) {
		for(byte b : bases) {
			boolean newValid=dimers.addWindowed(b);
			if(newValid) {
				hist[0][(int)(dimers.GC()*1024)]++;
				hist[1][(int)(dimers.strandedness()*1024)]++;
				hist[2][(int)(dimers.HH()*1024)]++;
				hist[3][(int)(dimers.PP()*1024)]++;
				hist[4][(int)(dimers.AAAT()*1024)]++;
				hist[5][(int)(dimers.CCCG()*1024)]++;
				hist[6][(int)(dimers.HMH()*1024)]++;
				hist[7][(int)(dimers.HHPP()*1024)]++;
			}
		}
	}
	
	private void outputResults(){
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();
		if(header){
			if(rowheader) {bb.append("Header\t");}
			bb.append("GC\tSTR\tHH\tPP\tAAAT\tCCCG\tHMH\tHHPP\n");
		}
		if(window<1) {
			if(rowheader) {bb.append("Mean\t");}
			bb.appendt(dimers.GC(), 5);
			bb.appendt(dimers.strandedness(), 5);
			bb.appendt(dimers.HH(), 5);
			bb.appendt(dimers.PP(), 5);
			bb.appendt(dimers.AAAT(), 5);
			bb.appendt(dimers.CCCG(), 5);
			bb.appendt(dimers.HMH(), 5);
			bb.append(dimers.HHPP(), 5);
			bb.nl();
		}else {
			if(rowheader) {bb.append("Mean\t");}
			for(int i=0; i<hist.length; i++) {
				bb.appendt(Tools.averageHistogram(hist[i])/1024, 5);
			}
			bb.set(bb.length()-1, '\n');
			if(rowheader) {bb.append("STDev\t");}
			for(int i=0; i<hist.length; i++) {
				bb.appendt(Tools.standardDeviationHistogram(hist[i])/1024, 5);
			}
			bb.set(bb.length()-1, '\n');
		}
		bsw.print(bb);
		
		errorState=bsw.poisonAndWait() | errorState;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	private int window=0;
	private final KmerTracker dimers;
	private boolean header=false;
	private boolean rowheader=false;
	private boolean printTime=false;
	private long[][] hist=new long[8][1025];
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private boolean errorState=false;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
