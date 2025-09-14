package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

import dna.Data;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Grab reads with specified numbers from a file.
 * TODO Note that much of this is ripped directly from ReformatReads, but is incorrect, because this class does not support dual output files.
 * @author Brian Bushnell
 * @date Jul 10, 2013
 *
 */
public class GetReads {
	
	public static void main(String[] args){
		GetReads x=new GetReads(args);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public GetReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Timer t=new Timer();

		Parser parser=new Parser();
		String in1=null;
		String in2=null;
		
		String qfin1=null;
		String qfin2=null;

		String out1=null;
		String out2=null;

		String qfout1=null;
		String qfout2=null;
		
		boolean errorState=false;
		long maxReads=-1;
		int passes=1;
		boolean testsize=false;
		boolean overwrite=false, append=false;
		float samplerate=1f;
		long sampleseed=1;

		
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		HashSet<Long> table=new HashSet<Long>();
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("id") || a.equals("number")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] b2=b.split(",");
				for(String c : b2){
					final long x, y;
					if(c.indexOf('-')>=0){
						String[] c2=c.split("-");
						assert(c2.length==2) : c;
						x=Long.parseLong(c2[0]);
						y=Long.parseLong(c2[1]);
					}else{
						x=y=Long.parseLong(c);
					}
					for(long z=x; z<=y; z++){
						table.add(z);
					}
				}
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				assert(b!=null) : "Bad parameter: "+arg;
				in1=b;
				if(b.indexOf('#')>-1 && !new File(b).exists()){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				assert(b!=null) : "Bad parameter: "+arg;
				out1=b;
				if(b.indexOf('#')>-1){
					out1=b.replace("#", "1");
					out2=b.replace("#", "2");
				}
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfin2")){
				qfin2=b;
			}else if(a.equals("qfout2")){
				qfout2=b;
			}else if(a.equals("testsize")){
				testsize=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					in1=arg.replace("#", "1");
					in2=arg.replace("#", "2");
				}
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(out1==null){
			if(out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
			out1="stdout";
		}
		
		if(!parser.setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		

		FileFormat ffin=FileFormat.testInput(in1, 0, null, true, true);
		FileFormat ffout=FileFormat.testOutput(out1, 0, null, true, overwrite, append, false);
		
		
		final boolean useSharedHeader=(ffin!=null && ffout!=null && ffin.samOrBam() && ffout.samOrBam());
		
		if(ffin!=null && ffout!=null && ffin.samOrBam() && (ffout.samOrBam() || ffout.bread())){
			throw new RuntimeException("\nDirect conversion of sam to sam or bread are not currently supported.\nAll other conversions are possible.");
		}
		
		
		ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ff1, ff2);
		}
		
		cris.setSampleRate(samplerate, sampleseed);
		outstream.println("Input is "+(cris.paired() ? "paired" : "unpaired"));
		cris.start(); //4567

		TextStreamWriter tsw=new TextStreamWriter(out1, overwrite, false, false);
		tsw.start();
		
		
		long readsProcessed=0;
		long basesProcessed=0;

		for(int pass=1; pass<=passes; pass++){
//			outstream.println("pass="+pass);
			if(pass>1){
				FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
				FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ff1, ff2);
				cris.setSampleRate(samplerate, sampleseed);
				cris.start();
			}
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin==null || ffin.samOrBam() || (r.mate!=null)==cris.paired());//ffin cannot be null
			}

			while(reads!=null && reads.size()>0 && !table.isEmpty() && ln!=null){//ln!=null is implied

				for(Read r1 : reads){
					{
						readsProcessed++;
						basesProcessed+=r1.length();
					}
					Read r2=r1.mate;
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=r2.length();
					}
					
					if(table.remove(r1.numericID)){
						tsw.println(r1);
						if(r2!=null){tsw.println(r2);}
						if(table.isEmpty()){break;}
					}
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
			errorState|=ReadWrite.closeStream(cris);
		}

		if(tsw!=null){
			tsw.poisonAndWait();
		}
		
		errorState|=(cris.errorState());
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length()))*passes;
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+Tools.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(errorState){
			throw new RuntimeException("GetReads terminated in an error state; the output may be corrupt.");
		}
		
	}

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
