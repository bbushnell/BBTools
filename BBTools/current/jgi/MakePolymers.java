package jgi;

import java.io.PrintStream;

import dna.AminoAcid;
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
import stream.FASTQ;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class MakePolymers {
	
	public static void main(String[] args){
		Timer t=new Timer();
		MakePolymers x=new MakePolymers(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public MakePolymers(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("k")){
				mink=maxk=Integer.parseInt(b);
			}else if(a.equals("mink")){
				mink=Integer.parseInt(b);
			}else if(a.equals("maxk")){
				maxk=Integer.parseInt(b);
			}else if(a.equals("len") || a.equals("minlen")){
				minLen=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);
	}
	
	void process(Timer t){

		final ByteStreamWriter bsw;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		for(int i=mink; i<=maxk; i++){
			writeSequence(i, bsw);
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private void writeSequence(int k, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		
		final int minLen2=((minLen+k-1)/k)*k;
		final int minCount;
		if(minLen2-minLen>=k-1){
			minCount=minLen2/k;
		}else{
			minCount=minLen2/k+1;
		}
		
		final long max=(1<<(2*k))-1;
		for(long kmer=0; kmer<=max; kmer++){
			bb.append('>').append(k).append('_').append(kmer).append('\n');
			for(int i=0; i<minCount; i++){
				basesProcessed+=k;
				toBytes(kmer, k, bb);
			}
			readsProcessed++;
			bb.nl();
			if(bb.length>=16384){
				bsw.print(bb);
				bb.clear();
			}
		}
		if(bb.length>0){
			bsw.print(bb);
			bb.clear();
		}
	}
	
	public static final void toBytes(long kmer, int k, ByteBuilder bb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			bb.append(AminoAcid.numberToBase[x]);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private long readsProcessed=0;
	private long basesProcessed=0;
	
	private int mink=1, maxk=1;
	
	private int minLen=31;
	
	private String out1=null;

	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
