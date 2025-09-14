package hiseq;

import java.util.ArrayList;

import barcode.Barcode;
import barcode.BarcodeCounter;
import barcode.PCRMatrixHDist;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class PlotReadPosition {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		PlotReadPosition x=new PlotReadPosition(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public PlotReadPosition(String[] args){
		
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

			if(a.equals("expected") || a.equals("names") || a.equals("barcodes")){
				expectedPath=b;
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
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

		IlluminaHeaderParser1.PARSE_COMMENT=true;
		ffout1=FileFormat.testOutput(out1, FileFormat.HEADER, ".header", true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		final int delimiter=ffin1.barcodeDelimiter();
		final int len1=ffin1.barcodeLength(1);
		final int len2=ffin1.barcodeLength(2);
		matrix=new PCRMatrixHDist(len1, len2, delimiter, true);
		ArrayList<String> expected=BarcodeCounter.loadBarcodes(expectedPath, delimiter);
		matrix.populateExpected(expected);
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

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, "x\ty\thdist\n", false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0, basesProcessed=0;
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			ByteBuilder bb=new ByteBuilder(128);
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					r1.mate=null;
					
					//  *********  Process reads here  *********
					ihp.parse(r1);
//					System.err.println(ihp);
					String barcode=ihp.barcode();
					Barcode bc=matrix.findClosest(barcode, 99, 0);
					int dist=barcode.length();
					if(bc!=null) {
						dist=bc.hdist(barcode);
					}
					bb.clear();
					bb.append(ihp.xPos()).tab().append(ihp.yPos()).tab().append(dist);
					r1.id=bb.toString();
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
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String expectedPath=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private final PCRMatrixHDist matrix;
	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
