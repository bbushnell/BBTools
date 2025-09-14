package synth;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
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
import structures.ListNum;


/**
 * Fuses reads together randomly to make chimeric reads.
 * @author Brian Bushnell
 * @date Oct 7, 2014
 *
 */
public class MakeChimeras {

	public static void main(String[] args){
		Timer t=new Timer();
		MakeChimeras x=new MakeChimeras(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public MakeChimeras(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
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
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(a.equals("forcelength")){
				forceLength=Integer.parseInt(b);
			}else if(a.equals("readsout") || a.equals("chimeras")){
				readsOut=Parse.parseKMG(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			readsIn=parser.maxReads;
			
			in1=parser.in1;
			qfin1=parser.qfin1;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}
	
	void process(Timer t){
		assert(readsOut>0) : "Please set the 'readsout' flag to a positive integer.";
		
		ArrayList<Read> source=new ArrayList<Read>();
		{
			final ConcurrentReadInputStream cris;
			{
				cris=ConcurrentReadInputStream.getReadInputStream(readsIn, false, ffin1, null, qfin1, null);
				if(verbose){outstream.println("Started cris");}
				cris.start(); //4567
			}
			assert(!cris.paired());
			
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
						assert(r1.mate==null);

						final int initialLength1=r1.length();
						
						if(initialLength1>0){
							source.add(r1);
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
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		
		
		if(readsOut>=0){
			t.start();
			
			final TextStreamWriter tsw;
			if(ffout1==null){
				tsw=null;
			}else{
				tsw=new TextStreamWriter(ffout1);
				tsw.start();
			}
			
			final Random randy=Shared.threadLocalRandom();
			
			long readsProcessed=0;
			long basesProcessed=0;
			final int mod=source.size();
			for(long i=0; i<readsOut; i++){
				Read a=source.get(randy.nextInt(mod));
				Read b=source.get(randy.nextInt(mod));
				Read c=makeChimera(a, b, randy, i);
				if(c==null){
					i--;
				}else{
					if(tsw!=null && c!=null){
						tsw.println(c);
						readsProcessed++;
						basesProcessed+=c.length();
					}
				}
			}
			
			if(tsw!=null){errorState|=tsw.poisonAndWait();}
			
			t.stop();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * @param a
	 * @param b
	 * @param randy
	 * @return
	 */
	private Read makeChimera(Read a, Read b, Random randy, long numericID) {
		final String id=a.id+" ~ "+b.id;

		final Read a2, b2;
		if(forceLength>0){
//			if(a.length()>b.length()){
//				Read c=a;
//				a=b;
//				b=c;
//			}
			a2=getPiece(a, randy, forceLength);
			b2=getPiece(b, randy, b.length()-forceLength);
			if(a2==null || b2==null){return null;}
		}else{
			a2=getPiece(a, randy);
			b2=getPiece(b, randy);
		}
		
		a=b=null;
		
		final byte[] abases=a2.bases, bbases=b2.bases, aquals=a2.quality, bquals=b2.quality;
		final int alen=a2.length(), blen=b2.length();
		final int len=a2.length()+b2.length();
		byte[] bases=new byte[len];
		byte[] quals=(aquals==null || bquals==null) ? null : new byte[len];

		for(int i=0; i<alen; i++){
			bases[i]=abases[i];
			if(quals!=null){quals[i]=aquals[i];}
		}
		for(int i=alen, j=0; j<blen; i++, j++){
			bases[i]=bbases[j];
			if(quals!=null){quals[i]=bquals[j];}
		}
		
		Read r=new Read(bases, quals, id, numericID);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}

	/**
	 * @param value
	 * @param randy
	 * @return
	 */
	private static Read getPiece(Read a, Random randy) {
		int len=randy.nextInt(a.length())+1;
		
		final int start;
		if(Tools.nextBoolean(randy)){
			if(Tools.nextBoolean(randy)){
				start=0;
			}else{
				start=a.length()-len;
			}
		}else{
			int range=a.length()-len;
			start=randy.nextInt(range+1);
		}
		
		byte[] bases=KillSwitch.copyOfRange(a.bases, start, start+len);
		byte[] quals=a.quality==null ? null : KillSwitch.copyOfRange(a.quality, start, start+len);
		
		Read r=new Read(bases, quals, a.id, a.numericID);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}

	/**
	 * @param value
	 * @param randy
	 * @return
	 */
	private Read getPiece(Read a, Random randy, int len) {
		len=Tools.min(len, a.length());
		if(len<1){return null;}
		
		final int start;
		if(Tools.nextBoolean(randy)){
			if(Tools.nextBoolean(randy)){
				start=0;
			}else{
				start=a.length()-len;
			}
		}else{
			int range=a.length()-len;
			start=randy.nextInt(range+1);
		}
		
		byte[] bases=KillSwitch.copyOfRange(a.bases, start, start+len);
		byte[] quals=a.quality==null ? null : KillSwitch.copyOfRange(a.quality, start, start+len);
		
		Read r=new Read(bases, quals, a.id, a.numericID);
		if(Tools.nextBoolean(randy)){r.reverseComplement();}
		return r;
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	
	private String qfin1=null;

	private String out1=null;
	
	private String extin=null;
	private String extout=null;

	private int forceLength=0;
	
	/*--------------------------------------------------------------*/

	private long readsIn=-1;
	private long readsOut=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;

	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
