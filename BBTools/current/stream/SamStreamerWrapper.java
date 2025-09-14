package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;
import var2.SamFilter;
import var2.ScafMap;

public class SamStreamerWrapper {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		SamStreamerWrapper x=new SamStreamerWrapper(args);

		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	SamStreamerWrapper(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=true;
		
		filter=new SamFilter();
		boolean doFilter=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
				Parser.parseSam(arg, a, b);
				fixCigar=true;
			}
			
			//Filter
			else if(a.equals("filter")){
				doFilter=Parse.parseBoolean(b);
			}else if(filter.parse(arg, a, b)){
				//do nothing
			}
				
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if(!doFilter){filter=null;}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, true, false, true);
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		
		if(!forceParse && !fixCigar && (ffout1==null || !ffout1.samOrBam())){
			SamLine.PARSE_2=false;
			SamLine.PARSE_5=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=false;
		}
	}
	
	void process(Timer t){
		
		if(ref!=null){
			ScafMap.loadReference(ref, true);
			SamLine.RNAME_AS_BYTES=false;
		}
		
		boolean useSharedHeader=(ffout1!=null && ffout1.samOrBam());
		final SamReadStreamer ss=new SamReadStreamer(ffin1, ordered ? 1 : SamStreamer.DEFAULT_THREADS, useSharedHeader, maxReads);
		ss.start();

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, useSharedHeader);
			ros.start();
		}else{ros=null;}

		long readsProcessed=0, readsOut=0;
		long basesProcessed=0, basesOut=0;
		for(ListNum<Read> ln=ss.nextReads(); ln!=null && ln.size()>0; ln=ss.nextReads()){
			ArrayList<Read> list=ln.list;
			if(verbose){outstream.println("Got list of size "+ln.size());}
			ArrayList<Read> out=new ArrayList<Read>(list.size());
			for(Read r : list){
				final int len=r.length();
				readsProcessed++;
				basesProcessed+=len;
				final SamLine sl=r.samline;
				boolean keep=filter==null || filter.passesFilter(sl);
				if(keep){
					if(sl.cigar!=null){
						if(fixCigar){
							if(SamLine.VERSION==1.3f){
								sl.cigar=SamLine.toCigar13(sl.cigar);
							}else{
								r.toLongMatchString(false);
								int start=sl.pos-1;
								int stop=start+Read.calcMatchLength(r.match)-1;
								sl.cigar=SamLine.toCigar14(r.match, start, stop, Integer.MAX_VALUE, r.bases);
							}
						}
//						if(SamLine.MAKE_MD_TAG){
//							//No code for this currently except from ChromArray
//							if(sl.optional==null){sl.optional=new ArrayList<String>(1);}
//							for(int i=0; i<sl.optional.size(); i++){
//								if(sl.optional.get(i).startsWith("MD:")){sl.optional.remove(i);}
//							}
//							String md=makeMdTag(r.chrom, r.start, r.match, r.bases, scafloc, scaflen);
//							if(md!=null){sl.optional.add(md);}
//						}
					}
					
					out.add(r);
					readsOut++;
					basesOut+=len;
				}
			}
			if(ros!=null){ros.add(out, ln.id);}
		}
		errorState|=ss.errorState;
		errorState|=ReadWrite.closeStream(ros);
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+Tools.format("%.2f Mbp/sec", (basesProcessed/(double)(t.elapsed))*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getSimpleName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	SamFilter filter;
	
	private String in1=null;
	private String out1=null;
	private String ref=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	public boolean errorState=false;
	public boolean ordered=true;
	private long maxReads=-1;
	private boolean forceParse;
	private boolean fixCigar;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
