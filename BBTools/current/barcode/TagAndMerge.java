package barcode;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.LineParserS4;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * This is intended to run simultaneously on all files from a demultiplexed lane,
 * labeling each read with the barcode of the input file, then printing them
 * to a single output file.
 * For example, all reads from input file 52866.4.475040.GAGGCCGCCA-TTATCTAGCT.filter-DNA.fastq.gz
 * would have "(tab)GAGGCCGCCA+TTATCTAGCT" appended to their headers.
 * @author Brian Bushnell
 * @date May 10, 2024
 *
 */
public class TagAndMerge {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		TagAndMerge x=new TagAndMerge(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public TagAndMerge(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("trim")){
				trimLen=(b==null ? 1 : Integer.parseInt(b));
			}else if(a.equals("r1only") || a.equals("read1only") || a.equals("dropr2")){
				dropR2=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("shrinkHeader") || a.equals("shrink")){
				shrinkHeader=Parse.parseBoolean(b);
			}else if(a.equals("barcodes") || a.equals("barcodesout")){
				barcodesOutFile=b;
			}else if(a.equals("remap")){
				remap=(b==null ? null : b.toCharArray());
			}else if(a.equals("in")){
				in.clear();
				String[] b2=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
				for(String b3 : b2){in.add(b3);}
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, false);
	}
	
	void process(Timer t){

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;

			for(String s : in) {
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		for(String s : in) {
			processInner(s, ros);
		}
		
		ReadWrite.closeStream(ros);
		outstream.println("Finished processing "+in.size()+" files.");
		
		if(barcodesOutFile!=null) {
			TextStreamWriter tsw=new TextStreamWriter(barcodesOutFile, overwrite, append, false);
			tsw.start();
			for(String tag : barcodes) {tsw.println(tag);}
			tsw.poisonAndWait();
		}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
	}
	
	void processInner(String fname, ConcurrentReadOutputStream ros) {
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		String tag=Barcode.parseBarcodeFromFname(fname);
		assert(tag!=null) : "Can't find barcode in filename "+fname;
		if(remap!=null) {
			for(int i=0; i<remap.length; i+=2) {
				tag=tag.replace(remap[i], remap[i+1]);
			}
		}
		if(tag!=null && Barcode.isBarcode(tag)) {
			barcodes.add(tag);
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			cris.start();
		}
		boolean paired=cris.paired();
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin==null || ffin.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					processPair(r1, tag);
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
		ReadWrite.closeStream(cris);
	}
	
	private void processPair(Read r, String tag) {
		readsProcessed+=r.pairCount();
		basesProcessed+=r.pairLength();
		if(dropR2) {r.mate=null;}
		processRead(r, tag);
		processRead(r.mate, tag);
		readsOut+=r.pairCount();
		basesOut+=r.pairLength();
	}
	
	private void processRead(Read r, String tag) {
		if(r==null) {return;}
		if(trimLen>=0 && trimLen<r.length()) {
			final int x=TrimRead.trimToPosition(r, 0, trimLen-1, 1);
			assert(x>0) : x+", "+trimLen+", "+r.length();
		}
		bb.clear();
		if(!shrinkHeader) {
			bb.append(r.id).tab().append(tag);
		}else {
			//LH00223:36:22JWNJLT3:8:1101:31031:1064 1:N:0:AAAGACCAGG+GAGTTTGAGA
//			ihp.parse(r.id);
//			final LineParser lp=ihp.lp();
			lp.set(r.id);
			bb.colon().colon().colon().appendTerm(lp, 3).colon().appendTerm(lp, 4).colon();
			bb.space().appendTerm(lp, 7);
//			bb.colon().appendTerm(lp,  8).colon().appendTerm(lp, 9).colon().appendTerm(lp,  10);
			bb.tab().append(tag);
		}
		r.id=bb.toString();
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	private String out1=null;
	private String barcodesOutFile=null;
	private LinkedHashSet<String> barcodes=new LinkedHashSet<String>();
	
	private int trimLen=-1;
	private boolean dropR2=false;
	private boolean shrinkHeader=false;
	private char[] remap=new char[] {'-', '+'};
	
	private final FileFormat ffout1;
	private boolean overwrite=true;
	private boolean append=false;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long readsProcessed=0, basesProcessed=0;
	private long readsOut=0, basesOut=0;
	private final ByteBuilder bb=new ByteBuilder();
//	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	private final LineParserS4 lp=new LineParserS4(":::::: ");
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
