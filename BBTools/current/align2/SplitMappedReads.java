package align2;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.zip.ZipOutputStream;

import dna.Data;
import fileIO.ReadWrite;
import shared.Timer;
import stream.ConcurrentLegacyReadInputStream;
import stream.RTextInputStream;
import stream.Read;
import stream.SiteScore;
import structures.ListNum;

public class SplitMappedReads {
	
	
	public static void main(String[] args){

		String reads1=args[0];
		String reads2=args[1].equalsIgnoreCase("null") ?  null : args[1];
		String outname=args[2].equalsIgnoreCase("null") ?  "" : args[2];
		
		int minChrom=1;
		int maxChrom=25;
		if(args.length>3){
			minChrom=maxChrom=Byte.parseByte(args[3]);
			if(args.length>4){
				maxChrom=Byte.parseByte(args[4]);
			}
		}
		assert(minChrom<=maxChrom && minChrom>=0);
		
		SplitMappedReads smr=new SplitMappedReads(reads1, reads2, outname, minChrom, maxChrom);
		smr.process();
		
	}
	
	public SplitMappedReads(String fname1, String fname2, String outname_, int minChrom, int maxChrom){
		this(new RTextInputStream(fname1, fname2, -1), outname_, minChrom, maxChrom);
		assert(fname2==null || !fname1.equals(fname2)) : "Error - input files have same name.";
	}
	
	public SplitMappedReads(RTextInputStream stream_, String outname_, int minChrom, int maxChrom){
		stream=stream_;
		outname=outname_;
		paired=stream.paired();
//		assert(outname.contains("#")) : "Output file name must contain the character '#' to be used for chromosome number.";
		
		MIN_CHROM=minChrom;
		MAX_CHROM=maxChrom;
		assert(MIN_CHROM>=0);
		assert(MAX_CHROM>=MIN_CHROM);

		outArraySingle1=new OutputStream[maxChrom+1];
		printArraySingle1=new PrintWriter[maxChrom+1];
		bufferArraySingle1=new ArrayList[maxChrom+1];
		for(int i=minChrom; i<outArraySingle1.length; i++){
			bufferArraySingle1[i]=new ArrayList<Read>(WRITE_BUFFER);
			outArraySingle1[i]=ReadWrite.getOutputStream(outname.replace("#", "single_1_chr"+i), false, true, false);
			printArraySingle1[i]=new PrintWriter(outArraySingle1[i]);
			printArraySingle1[i].println("#Chromosome "+i+" Read 1 Singletons");
			printArraySingle1[i].println("#"+Read.header());
		}
		
		if(!paired){
			outArraySingle2=null;
			printArraySingle2=null;
			bufferArraySingle2=null;
			outArrayPaired1=null;
			printArrayPaired1=null;
			bufferArrayPaired1=null;
			outArrayPaired2=null;
			printArrayPaired2=null;
			bufferArrayPaired2=null;
		}else{
			
			outArraySingle2=new OutputStream[maxChrom+1];
			printArraySingle2=new PrintWriter[maxChrom+1];
			bufferArraySingle2=new ArrayList[maxChrom+1];
			for(int i=minChrom; i<outArraySingle2.length; i++){
				bufferArraySingle2[i]=new ArrayList<Read>(WRITE_BUFFER);
				outArraySingle2[i]=ReadWrite.getOutputStream(outname.replace("#", "single_2_chr"+i), false, true, false);
				printArraySingle2[i]=new PrintWriter(outArraySingle2[i]);
				printArraySingle2[i].println("#Chromosome "+i+" Read 2 Singletons");
				printArraySingle2[i].println("#"+Read.header());
			}
			
			outArrayPaired1=new OutputStream[maxChrom+1];
			printArrayPaired1=new PrintWriter[maxChrom+1];
			bufferArrayPaired1=new ArrayList[maxChrom+1];
			for(int i=minChrom; i<outArrayPaired1.length; i++){
				bufferArrayPaired1[i]=new ArrayList<Read>(WRITE_BUFFER);
				outArrayPaired1[i]=ReadWrite.getOutputStream(outname.replace("#", "paired_1_chr"+i), false, true, false);
				printArrayPaired1[i]=new PrintWriter(outArrayPaired1[i]);
				printArrayPaired1[i].println("#Chromosome "+i+" Read 1 Paired");
				printArrayPaired1[i].println("#"+Read.header());
			}
			
			outArrayPaired2=new OutputStream[maxChrom+1];
			printArrayPaired2=new PrintWriter[maxChrom+1];
			bufferArrayPaired2=new ArrayList[maxChrom+1];
			for(int i=minChrom; i<outArrayPaired2.length; i++){
				bufferArrayPaired2[i]=new ArrayList<Read>(WRITE_BUFFER);
				outArrayPaired2[i]=ReadWrite.getOutputStream(outname.replace("#", "paired_2_chr"+i), false, true, false);
				printArrayPaired2[i]=new PrintWriter(outArrayPaired2[i]);
				printArrayPaired2[i].println("#Chromosome "+i+" Read 2 Paired");
				printArrayPaired2[i].println("#"+Read.header());
			}
			
		}
		
		cris=(USE_CRIS ? new ConcurrentLegacyReadInputStream(stream, -1) : null);
	}
	
	public void process(){
		
		Timer t=new Timer();
		
		if(cris!=null){
			cris.start();
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				processReads(reads);
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}else{
			ArrayList<Read> reads=stream.nextList();
			while(reads!=null && reads.size()>0){
				processReads(reads);
				reads=stream.nextList();
			}
		}
		
		synchronized(this){this.notifyAll();}
		
		finish();
		
		t.stop();
		Data.sysout.println("Time:\t"+t);
	}
	

	
	private void processReads(ArrayList<Read> reads){
		for(Read r : reads){
			addRead(r, 1);
			if(r.mate!=null){
				addRead(r.mate, 2);
			}
		}
	}

	
	private void addRead(Read r, int side){
		
		if(r.chrom<1 && r.numSites()>0){
			SiteScore ss=r.topSite(); //Should not be necessary
			r.start=ss.start;
			r.stop=ss.stop;
			r.chrom=ss.chrom;
			r.setStrand(ss.strand);
		}
		
		//Ensure no superfluous data is written
		r.sites=null;
		r.originalSite=null;
		r.samline=null;
		
//		System.err.println("Adding to chrom "+r.chrom+", side "+side+", paired="+r.paired+", "+(r.list==null ? "null" : r.list.size()));
		if(r.chrom<MIN_CHROM || r.chrom>MAX_CHROM){return;}
		
		final PrintWriter writer;
		final ArrayList<Read> list;
		
		if(side==1){
			if(r.paired()){
				writer=printArrayPaired1[r.chrom];
				list=bufferArrayPaired1[r.chrom];
			}else{
				writer=printArraySingle1[r.chrom];
				list=bufferArraySingle1[r.chrom];
			}
		}else{
			assert(side==2);
			if(r.paired()){
				writer=printArrayPaired2[r.chrom];
				list=bufferArrayPaired2[r.chrom];
			}else{
				writer=printArraySingle2[r.chrom];
				list=bufferArraySingle2[r.chrom];
			}
		}
		
		assert(list.size()<WRITE_BUFFER);
		list.add(r);
		
		if(list.size()>=WRITE_BUFFER){
			writeList((ArrayList<Read>)list.clone(), writer);
			list.clear();
		}
	}
	
	
	private static void writeList(ArrayList<Read> list, PrintWriter writer){
		synchronized(writer){
			for(Read r : list){
				writer.println(r.toText(true));
			}
		}
	}
	
	
	public void finish(){

		final PrintWriter[][] writers=new PrintWriter[][] {printArraySingle1, printArraySingle2, printArrayPaired1, printArrayPaired2};
		final OutputStream[][] streams=new OutputStream[][] {outArraySingle1, outArraySingle2, outArrayPaired1, outArrayPaired2};
		final ArrayList<Read>[][] buffers=new ArrayList[][] {bufferArraySingle1, bufferArraySingle2, bufferArrayPaired1, bufferArrayPaired2};
		

		for(int x=0; x<buffers.length; x++){


			PrintWriter[] printArray=writers[x];
			ArrayList<Read>[] bufferArray=buffers[x];

			for(int i=0; printArray!=null && i<printArray.length; i++){
				PrintWriter writer=printArray[i];
				ArrayList<Read> list=bufferArray[i];

				if(list!=null && !list.isEmpty()){
					writeList(list, writer);
					list=null;
				}
			}
		}
		
		//TODO: Wait for writing to finish, if it is done in threads.
		
		
		for(int x=0; x<writers.length; x++){


			PrintWriter[] printArray=writers[x];
			OutputStream[] outArray=streams[x];
			
			for(int i=0; printArray!=null && i<printArray.length; i++){
				if(printArray[i]!=null){
					synchronized(printArray[i]){
						printArray[i].flush();
						if(outArray[i].getClass()==ZipOutputStream.class){
							ZipOutputStream zos=(ZipOutputStream)outArray[i];
							try {
								zos.closeEntry();
								zos.finish();
							} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
						printArray[i].close();
						try {
							outArray[i].close();
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
		}
		
//		if(cris!=null){cris.shutdown();}
//		stream.shutdown();
		
		if(cris!=null){ReadWrite.closeStream(cris);}
		else{stream.close();}
	}
	
	
	public final String outname;
	private final RTextInputStream stream;
	private final ConcurrentLegacyReadInputStream cris;
	
	private final OutputStream[] outArraySingle1;
	private final PrintWriter[] printArraySingle1;
	private final ArrayList<Read>[] bufferArraySingle1;
	
	private final OutputStream[] outArraySingle2;
	private final PrintWriter[] printArraySingle2;
	private final ArrayList<Read>[] bufferArraySingle2;
	
	private final OutputStream[] outArrayPaired1;
	private final PrintWriter[] printArrayPaired1;
	private final ArrayList<Read>[] bufferArrayPaired1;
	
	private final OutputStream[] outArrayPaired2;
	private final PrintWriter[] printArrayPaired2;
	private final ArrayList<Read>[] bufferArrayPaired2;

	private final int MIN_CHROM;
	private final int MAX_CHROM;
	
	public final boolean paired;
	
	public static boolean USE_CRIS=true; //Similar speed either way.  "true" may be better with many threads.
	
	public static final int WRITE_BUFFER=400; //Bigger number uses more memory, for less frequent writes.
	
	
}
