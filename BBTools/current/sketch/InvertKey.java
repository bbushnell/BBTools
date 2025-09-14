package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.AminoAcid;
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
import structures.ListNum;
import structures.LongHashSet;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class InvertKey extends SketchObject {
	
	public static void main(String[] args){
		Timer t=new Timer();
		InvertKey x=new InvertKey(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public InvertKey(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		int k_=32, k2_=0;
		
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
			}else if(a.equals("key")){
				keyString=b;
			}else if(a.equals("out")){
				out1=b;
			}else if(a.equalsIgnoreCase("k")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.indexOf(',')>=0){
					String[] bsplit=b.split(",");
					assert(bsplit.length==2) : "Bad argument "+arg;
					int x=Integer.parseInt(bsplit[0]);
					int y=Integer.parseInt(bsplit[1]);
					k_=Tools.max(x, y);
					k2_=Tools.min(x, y);
					if(k_==k2_){k2_=0;}
				}else{
					k_=Integer.parseInt(b);
					k2_=0;
				}
			}else if(a.equalsIgnoreCase("printonce")){
				printOnce=Parse.parseBoolean(b);
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				out1=arg;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		k=k_;
		k2=k2_;
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift)); //Conditional allows K=32
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
		
		SketchObject.postParse();
		
		if(keyString.indexOf(',')>0){
			String[] split=keyString.split(",");
			set=new LongHashSet(split.length*2);
			for(String s : split){
				long x=Long.MAX_VALUE-Sketch.parseA48(s);
				set.add(x);
//				assert(set.contains(x)) : x+", "+set.size()+", "+set.toStringListView();
			}
			key0=-1;
//			System.err.println(set.toStringListView()+", "+set.size());
			assert(!set.isEmpty());
		}else if(keyString.endsWith(".sketch")){
			SketchTool tool=new SketchTool(10000, 0, false, false, rcomp);
			Sketch sk=tool.loadSketchesFromFile(keyString, null, 0, 1000000, SketchObject.ONE_SKETCH, 1f, 0f, 0f, (byte)0, false).get(0);
			set=new LongHashSet(sk.length()*2);
			for(long x : sk.keys){set.add(Long.MAX_VALUE-x);}
			key0=-1;
//			System.err.println(set.toStringListView()+", "+set.size());
			assert(!set.isEmpty());
		}else{
			key0=Long.MAX_VALUE-Sketch.parseA48(keyString);
			set=null;
//			System.err.println(key0);
		}
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
//		if(verbose){
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
//		}

		final ByteStreamWriter bsw;
		if(out1!=null){
			fasta=ffout1.fasta() && !out1.endsWith(".txt");
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}else{bsw=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		boolean finished=false;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(reads!=null && reads.size()>0 && !finished){
				
				for(int idx=0; idx<reads.size() && !finished; idx++){
					final Read r1=reads.get(idx);

					finished=invert(key0, r1, bsw);
					
					final int initialLength1=r1.length();
					
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
		
		errorState|=(ReadWrite.closeStream(cris));
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState && !finished && maxReads<1){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private boolean invert(long key2, Read r, ByteStreamWriter bsw) {
		final byte[] bases=r.bases;
		
		long kmer=0;
		long rkmer=0;
		int len=0;
		

//		System.err.println("Looking for "+key+"\t"+Sketch.toA48(key)+"\t"+Sketch.toA48(Long.MAX_VALUE-key));
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){len=0; rkmer=0;}else{len++;}
			if(len>=k){
				kmersProcessed++;
				final long hashcode=hash(kmer, rkmer);
				boolean found=(key0>=0 ? hashcode==key0 : set.contains(hashcode));
				if(found){
					if(fasta){bsw.println(">"+Sketch.toA48(Long.MAX_VALUE-hashcode)+" "+(i-k+1)+" "+r.id);}
					bsw.println(AminoAcid.kmerToString(Tools.min(kmer, rkmer), k));
					if(printOnce){
						if(key0>=0){return true;}
						else{
							set.remove(hashcode);
							return set.isEmpty();
						}
					}
				}
			}
		}
		return false;
	}
	
	/*--------------------------------------------------------------*/
	
	final long key0;
	final LongHashSet set;
	
	final int shift;
	final int shift2;
	final long mask;
	
	boolean printOnce=true;
	long kmersProcessed=0;
	
	private String in1=null;
	boolean fasta;
	boolean sketch;
	private String keyString=null;

	private String out1="stdout.fa";
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
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
