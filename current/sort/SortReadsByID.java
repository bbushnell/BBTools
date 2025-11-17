package sort;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Nov 1, 2012
 *
 */
public class SortReadsByID {
	
	/**
	 * Program entry point for command-line execution.
	 * Parses arguments and initiates the read sorting process.
	 * @param args Command-line arguments including input/output files and options
	 */
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		String in1=null;
		String in2=null;
		String out="raw_idsorted#.txt.gz";
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("i") || a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				assert(b!=null) : "Bad parameter: "+arg;
				in1=b;
				if(b.indexOf('#')>=0){
					in1=b.replaceFirst("#", "1");
					in2=b.replaceFirst("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("o") || a.equals("out") || a.equals("output")){
				out=b;
			}else if(a.endsWith("renumber")){
				RENUMBER=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
				Data.sysout.println("Set overwrite to "+overwrite);
			}else if(a.endsWith("blocksize")){
				BLOCKSIZE=Integer.parseInt(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		if(in1==null){throw new RuntimeException("Please specify input file.");}
		if(out==null){throw new RuntimeException("Please specify output file.");}
		if(in1.equalsIgnoreCase(in2) || in1.equalsIgnoreCase(out) || (in2!=null && in2.equalsIgnoreCase(out))){
			throw new RuntimeException("Duplicate filenames.");
		}

		if(out!=null && !out.contains("#")){
			throw new RuntimeException("Output filename must contain '#' symbol.");
		}
		
		SortReadsByID srid=new SortReadsByID(in1, in2, out);
		srid.process();
	}
	
	
	/**
	 * Executes the external sorting algorithm for reads.
	 * Phase 1: Reads input and distributes reads into bins based on ID ranges.
	 * Phase 2: Sorts each bin individually and merges into final sorted output.
	 * Uses temporary files to handle datasets larger than available memory.
	 */
	public void process(){

		Timer tRead=new Timer();
		Timer tSort=new Timer();
		Timer tAll=new Timer();
		
		tRead.start();
		tAll.start();
		
		final long maxReads=-1;
		ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
			cris.start(); //4567
		}
		
		HashMap<Integer, Block> map=new HashMap<Integer, Block>();
		
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				//System.err.println("reads.size()="+reads.size());
				for(Read r : reads){

					int bin=(int)(r.numericID/BLOCKSIZE);
					Block b=map.get(bin);
					if(b==null){
						String o1=out.replaceFirst("#", "_bin"+bin+"_1");
						String o2=(cris.paired() && !OUT_INTERLEAVED) ? out.replaceFirst("#", "_bin"+bin+"_2") : null;
						b=new Block(o1, o2);
						map.put(bin, b);
					}
					b.add(r);
				}
				//System.err.println("returning list");
				cris.returnList(ln);
				//System.err.println("fetching list");
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			
			cris.returnList(ln);
			ReadWrite.closeStream(cris);
		}
		
		for(Block b : map.values()){b.close();}
		
		tRead.stop();
		Data.sysout.println("Read time:   \t"+tRead);
		tSort.start();
		
		String o1=out.replaceFirst("#", "1");
		String o2=(cris.paired() && !OUT_INTERLEAVED) ? out.replaceFirst("#", "2") : null;
		Block sorted=new Block(o1, o2);
		
		long count=0;
		
		ArrayList<Integer> keys=new ArrayList<Integer>();
		keys.addAll(map.keySet());
		Shared.sort(keys);
		for(Integer key : keys){
			Block b=map.get(key);
			b.join();
			map.remove(key);
			{
				FileFormat ff1=FileFormat.testInput(b.out1, FileFormat.FASTQ, null, true, true);
				FileFormat ff2=FileFormat.testInput(b.out2, FileFormat.FASTQ, null, true, true);
				cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2);
				cris.start(); //4567
			}
			ArrayList<Read> reads2=new ArrayList<Read>((int)b.count);
			count+=b.count;
			
			{
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);

				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					reads2.addAll(reads);
					//System.err.println("returning list");
					cris.returnList(ln);
					//System.err.println("fetching list");
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				
				cris.returnList(ln);
				ReadWrite.closeStream(cris);
			}
			
			Shared.sort(reads2, ReadComparatorID.comparator);
			for(Read r : reads2){sorted.add(r);}
			new File(b.out1).delete();
			if(b.out2!=null){new File(b.out2).delete();}
		}
		
		sorted.close();
		sorted.join();
		
		tSort.stop();
		tAll.stop();
		
		Data.sysout.println("Total reads: \t"+count);
		Data.sysout.println("Sort time:   \t"+tSort);
		Data.sysout.println("Total time:  \t"+tAll);
		
	}
	
	/**
	 * Constructs a SortReadsByID instance with input/output file paths.
	 * Determines output format (FASTQ, FASTA, or text) from output filename.
	 *
	 * @param in1_ Primary input file path
	 * @param in2_ Secondary input file path (may be null for single-end)
	 * @param out_ Output file path template containing '#' placeholder
	 */
	public SortReadsByID(String in1_, String in2_, String out_) {
		in1=in1_;
		in2=in2_;
		out=out_;
		
		FileFormat ff=FileFormat.testOutput(out, FileFormat.BREAD, null, true, false, append, false);
		outFastq=ff.fastq();
		outFasta=ff.fasta();
		outText=ff.bread();
	}

	/** Primary input file path */
	public String in1;
	/** Secondary input file path for paired reads */
	public String in2;
	/** Output file path template containing '#' placeholder */
	public String out;

	/** Whether output should be in text format */
	private final boolean outText;
	/** Whether output should be in FASTA format */
	private final boolean outFasta;
	/** Whether output should be in FASTQ format */
	private final boolean outFastq;
	
	/** Number of reads per bin for external sorting algorithm */
	public static int BLOCKSIZE=8000000;
	/** Whether to overwrite existing output files */
	public static boolean overwrite=true;
	/** Whether to append to existing output files */
	public static boolean append=false;
	/** Whether to renumber reads sequentially after sorting */
	public static boolean RENUMBER=false;
	/** Whether to output paired reads in interleaved format */
	public static boolean OUT_INTERLEAVED=false;
	
	/** Represents a temporary file block containing reads from a specific ID range.
	 * Manages writing reads to temporary files during the binning phase of external sort. */
	private class Block{
		
		/**
		 * Creates a new block with specified output file paths.
		 * Initializes text stream writers for the temporary files.
		 * @param out1_ Output path for first reads in pairs
		 * @param out2_ Output path for second reads in pairs (may be null)
		 */
		public Block(String out1_, String out2_){
			out1=out1_;
			out2=out2_;
			
			tsw1=new TextStreamWriter(out1, overwrite, false, false);
			tsw2=(out2==null ? null : new TextStreamWriter(out2, overwrite, false, false));
			
			tsw1.start();
			if(tsw2!=null){tsw2.start();}
		}
		
		/**
		 * Adds a read and its mate (if present) to this block's temporary files.
		 * Converts reads to the appropriate output format before writing.
		 * @param r The read to add to this block
		 */
		public void add(Read r){
			count++;
			Read r2=r.mate;
			
			ByteBuilder sb1=outText ? r.toText(true) : outFastq ? r.toFastq() : outFasta ? r.toFasta() : null;
			ByteBuilder sb2=r2==null ? null : outText ? r2.toText(true) : outFastq ? r2.toFastq() : outFasta ? r2.toFasta() : null;
			
			tsw1.print(sb1.append('\n'));
			if(sb2!=null){
				if(tsw2!=null){
					tsw2.print(sb2.append('\n'));
				}else{
					tsw1.print(sb2.append('\n')); //Interleaved
				}
			}
			
		}
		
		/** Closes the text stream writers by sending poison signals.
		 * Must be called before reading back the temporary files. */
		public void close(){
			tsw1.poison();
			if(tsw2!=null){tsw2.poison();}
		}
		
		/** Waits for the text stream writers to finish writing and close.
		 * Ensures all data is flushed to disk before proceeding. */
		public void join(){
			tsw1.waitForFinish();
			if(tsw2!=null){tsw2.waitForFinish();}
		}
		
		/** Output file path for first reads in pairs */
		String out1;
		/** Output file path for second reads in pairs */
		String out2;
		
		/** Text stream writer for first reads in pairs */
		TextStreamWriter tsw1;
		/** Text stream writer for second reads in pairs */
		TextStreamWriter tsw2;
		
		/** Number of reads processed by this block */
		long count=0;
		
	}
	
	
}
