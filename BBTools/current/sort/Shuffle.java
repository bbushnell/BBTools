package sort;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
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
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SamReadInputStream;
import structures.ListNum;
import tracker.ReadStats;


/**
 * Randomizes the order of reads.
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */
public class Shuffle {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		Shuffle x=new Shuffle(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public Shuffle(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, printClass ? getClass() : null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLine.SET_FROM_OK=true;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=true;
		
		int mode_=SHUFFLE;
		
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
			}else if(a.equals("shuffle")){
				mode_=SHUFFLE;
			}else if(a.equals("name")){
				mode_=SORT_NAME;
			}else if(a.equals("coordinate")){
				mode_=SORT_COORD;
			}else if(a.equals("sequence")){
				mode_=SORT_SEQ;
			}else if(a.equals("id")){
				mode_=SORT_ID;
			}else if(a.equals("mode")){
				if(b==null){
					throw new RuntimeException("mode must be shuffle, name, coordinate, sequence, or id.");
				}else if(b.equals("shuffle")){
					mode_=SHUFFLE;
				}else if(b.equals("name")){
					mode_=SORT_NAME;
				}else if(b.equals("coordinate")){
					mode_=SORT_COORD;
				}else if(b.equals("sequence")){
					mode_=SORT_SEQ;
				}else if(b.equals("id")){
					mode_=SORT_ID;
				}else{
					throw new RuntimeException("mode must be shuffle, name, coordinate, sequence, or id.");
				}
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		mode=mode_;
		assert(mode>=1 && mode<=5) : "mode must be shuffle, name, coordinate, sequence, or id.";
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		if(!setInterleaved){
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
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		useSharedHeader=(ffout1.samOrBam() && ffout1!=null && ffout1.samOrBam());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		ArrayList<Read> bigList=new ArrayList<Read>(65530);
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					bigList.add(r1);
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
		errorState|=ReadStats.writeAll();
		
		if(mode==SHUFFLE){
			Collections.shuffle(bigList);
		}else if(mode==SORT_NAME){
			Shared.sort(bigList, ReadComparatorName.comparator);
		}else if(mode==SORT_SEQ){
			Shared.sort(bigList, ReadComparatorTopological.comparator);
		}else if(mode==SORT_COORD){
			Shared.sort(bigList, new ReadComparatorMapping());
		}else if(mode==SORT_ID){
			Shared.sort(bigList, ReadComparatorID.comparator);
		}else{
			assert(false) : "No mode set.";
		}
		
		if(ffout1!=null){
			final ByteStreamWriter bsw1, bsw2;
			if(ffout1!=null){
				bsw1=new ByteStreamWriter(ffout1);
				bsw1.start();
				if(useSharedHeader){writeHeader(bsw1);}
			}else{bsw1=null;}
			if(ffout2!=null){
				bsw2=new ByteStreamWriter(ffout2);
				bsw2.start();
				if(useSharedHeader){writeHeader(bsw1);}
			}else{bsw2=null;}
			final boolean b=(bsw2==null);
			for(int i=0, lim=bigList.size(); i<lim; i++){
				final Read r1=bigList.set(i, null);
				final Read r2=r1.mate;
				bsw1.println(r1, b);
				if(r2!=null && !b){bsw2.println(r2);}
			}
			if(bsw1!=null){errorState|=bsw1.poisonAndWait();}
			if(bsw2!=null){errorState|=bsw2.poisonAndWait();}
		}
		
		t.stop();
		
		if(showSpeed){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeHeader(ByteStreamWriter bsw){
		ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
		if(list==null){
			System.err.println("Header was null.");
		}else{
			for(byte[] line : list){
				bsw.print(line).print('\n');
			}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static class ShuffleThread extends Thread{
		
		public ShuffleThread(String in1_, String in2_, String out1_, String out2_, int mode_, boolean ow_){
			in1=in1_;
			in2=in2_;
			out1=out1_;
			out2=out2_;
			mode=mode_;
			ow=ow_;
		}
		
		@Override
		public synchronized void start(){
			addThread(1);
			super.start();
		}
		
		@Override
		public void run(){
			ArrayList<String> list=new ArrayList<String>();
			if(in1!=null){list.add("in1="+in1);}
			if(in2!=null){list.add("in1="+in2);}
			if(out1!=null){list.add("out1="+out1);}
			if(out2!=null){list.add("out2="+out2);}
			list.add("mode="+MODES[mode]);
			list.add("ow="+ow);
			try{
				Shuffle.main(list.toArray(new String[0]));
			}catch(Throwable e){
				System.err.println("Failed to shuffle "+in1+"\nException:"+e+"\n");
			}
			addThread(-1);
		}
		
		final String in1, in2;
		final String out1, out2;
		final int mode;
		final boolean ow;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	private final int mode;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	private final boolean useSharedHeader;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int maxShuffleThreads=1;
	private static int currentShuffleThreads=0;
	
	public static void setMaxThreads(final int x){
		assert(x>0);
		synchronized(SHUFFLE_LOCK){
			maxShuffleThreads=x;
		}
	}
	
	public static int addThread(final int x){
		synchronized(SHUFFLE_LOCK){
			while(x>0 && currentShuffleThreads>=maxShuffleThreads){
				try {
					SHUFFLE_LOCK.wait(2000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			currentShuffleThreads+=x;
			if(currentShuffleThreads<maxShuffleThreads){SHUFFLE_LOCK.notify();}
			return currentShuffleThreads;
		}
	}
	
	public static void waitForFinish(){
		synchronized(SHUFFLE_LOCK){
			while(currentShuffleThreads>=maxShuffleThreads){
				try {
					SHUFFLE_LOCK.wait(2000);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private static String SHUFFLE_LOCK=new String("SHUFFLE_LOCK");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	public static boolean showSpeed=true;
	public static boolean printClass=true;

	public static final int SHUFFLE=1, SORT_NAME=2, SORT_SEQ=3, SORT_COORD=4, SORT_ID=5;
	public static final String[] MODES={"shuffle", "name", "sequence", "coordinate", "id"};
	
	
}
