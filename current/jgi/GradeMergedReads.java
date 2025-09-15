package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.CustomHeader;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date May 20, 2014
 *
 */
public class GradeMergedReads {


	public static void main(String[] args){
		Timer t=new Timer();
		GradeMergedReads x=new GradeMergedReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public GradeMergedReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		FASTQ.DETECT_QUALITY=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("raw") || a.equals("raw1")){
				assert(b!=null) : "Bad parameter: "+arg;
				raw1=b;
				if(b.indexOf("#")>=0){
					raw1=b.replace('#', '1');
					raw2=b.replace('#', '2');
				}
			}else if(a.equals("raw2")){
				raw2=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				System.err.println("Unknown parameter "+i+": "+args[i]);
				assert(false) : "Unknown parameter "+i+": "+args[i];
//					+"\n"+arg+", "+parser.in1+", "+arg.contains("=")+", "+(arg.toLowerCase().startsWith("stdin")+", "+new File(arg).exists()+", "+new File(arg).getAbsolutePath());
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in=parser.in1;
			extin=parser.extin;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		FASTQ.PARSE_CUSTOM=false;

		ffin=FileFormat.testInput(in, FileFormat.FASTQ, extin, true, true);
	}
	
	void process(Timer t){
		
		long mergeable=0, total=0;
		if(raw1!=null){
			FileFormat ffraw1=FileFormat.testInput(raw1, FileFormat.FASTQ, extin, true, true);
			FileFormat ffraw2=FileFormat.testInput(raw2, FileFormat.FASTQ, extin, true, true);
			ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffraw1, ffraw2, null, null);
			cris.start();
			
			{
				
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);

				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					for(int idx=0; idx<reads.size(); idx++){
						final Read r1=reads.get(idx);
						String s=r1.id;
						total++;
						final int insert=parseInsert(r1.id);
						if(insert>0 && insert<r1.pairLength()){
							mergeable++;
						}
					}

					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				errorState|=ReadWrite.closeStream(cris);
			}
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		long readsProcessed=0;
		long basesProcessed=0;

		long correct=0;
		long tooLong=0;
		long tooShort=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin==null || ffin.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					final int insert=parseInsert(r1.id);
					
					int delta=insert-initialLength1;
					if(delta==0){
						correct++;
					}else if(delta>0){
						tooLong++;
					}else{
						tooShort++;
					}
					
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
		
		errorState|=ReadWrite.closeStream(cris);
		
		t.stop();
		
		long incorrect=tooShort+tooLong;
		double snr=10*Math.log10((correct+incorrect+0.0001)/(incorrect+0.0001));
		
		if(total>0){
			outstream.println("Input Total:            \t"+total+" pairs");
			outstream.println("Input Overlapping:      \t"+Tools.format("%.5f",mergeable*100.0/total)+"%\t"+mergeable+" reads");
		}
		outstream.println("Correct:                \t"+Tools.format("%.5f",correct*100.0/readsProcessed)+"%\t"+correct+" reads");
		outstream.println("Incorrect:              \t"+Tools.format("%.5f",incorrect*100.0/readsProcessed)+"%\t"+incorrect+" reads");
		outstream.println("Too Short:              \t"+Tools.format("%.5f",tooShort*100.0/readsProcessed)+"%\t"+tooShort+" reads");
		outstream.println("Too Long:               \t"+Tools.format("%.5f",tooLong*100.0/readsProcessed)+"%\t"+tooLong+" reads");
		outstream.println("SNR:                    \t"+Tools.format("%.3f",snr));

		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			throw new RuntimeException("GradeMergedReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	public static int parseInsert(String s){
		if(s.startsWith("SYN")){
			CustomHeader h=new CustomHeader(s, 0);
			return h.insert;
		}
		assert(s.startsWith("insert=")) : "Can't parse insert size for header "+s;
		
//		int space=s.indexOf(' ');
//		if(space<0){space=s.length();}
		int space=s.length();
		int equals=s.indexOf('=');
		for(int i=equals+1;i<s.length();i++){//For programs that rename my reads!
			if(!Tools.isDigit(s.charAt(i))){
				space=i;
				break;
			}
		}
		s=s.substring(equals+1, space);
		int insert=Integer.parseInt(s);
		return insert;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in=null;
	
	private String extin=null;
	
	private String raw1=null;
	private String raw2=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	

}
