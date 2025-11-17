package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.regex.Pattern;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import hiseq.IlluminaHeaderParser2;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.Quantizer;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Aug 23, 2013
 *
 */
public class RenameReads {
	
	/**
	 * Program entry point.
	 * Creates RenameReads instance, processes reads, and closes output streams.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameReads x=new RenameReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs RenameReads with command-line argument parsing.
	 * Parses input/output files, renaming options, and processing parameters.
	 * Sets up file formats and validates output file accessibility.
	 * @param args Command-line arguments specifying input files and renaming options
	 */
	public RenameReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
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
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equalsIgnoreCase("fixSRA")){
				fixSRA=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("prefix") || a.equals("p")){
				prefix=b;
			}else if(a.equals("suffix")){
				suffix=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
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
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("trimright")){
				trimRight=Integer.parseInt(b);
			}else if(a.equals("trimleft")){
				trimLeft=Integer.parseInt(b);
			}else if(a.equals("trimbeforesymbol")){
				trimBeforeSymbol=Integer.parseInt(b);
			}else if(a.equals("symbol")){
				String s=Parse.parseSymbol(b);
				assert(s.length()==1) : "'"+s+"'";
				symbol=s.charAt(0);
			}
			
			else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("renamebyinsert")){
				renameByInsert=Parse.parseBoolean(b);
			}else if(a.equals("renamebytrim")){
				renameByTrim=Parse.parseBoolean(b);
			}else if(a.equals("renamebycoordinates") || a.equals("coordinates") || a.equals("coords")){
				renameByCoords=Parse.parseBoolean(b);
			}else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}else if(a.equals("addprefix")){
				addPrefix=Parse.parseBoolean(b);
			}else if(a.equals("addpairnum")){
				addPairnum=Parse.parseBoolean(b);
			}else if(a.equals("prefixonly")){
				prefixOnly=Parse.parseBoolean(b);
			}else if(a.equals("underscore") || a.equals("addunderscore")){
				addUnderscore=Parse.parseBoolean(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
			
			renameByMapping=FASTQ.TAG_CUSTOM;
		}
		
		if(addUnderscore && prefix!=null && !prefix.endsWith("_") && !prefixOnly){prefix+="_";}
		if(!addPairnum){pairnums=new String[] {"",""};}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
//		assert(false) : prefix;
		if(prefix==null || prefix.length()<1){prefix="";}
//		else if(!prefix.endsWith("_") && !prefixOnly){
//			prefix=prefix+"_";
//		}

		if(renameByInsert){
			prefix="insert=";
			FASTQ.PARSE_CUSTOM=true;
		}else if(renameByTrim){
			prefix="";
			FASTQ.PARSE_CUSTOM=true;
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
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		if(renameByMapping){
			assert(ffout1==null || ffout1.fastq()) : "Currently renameByMapping requires fastq output.";
		}
	}
	
	/**
	 * Main processing method that reads input sequences and applies renaming.
	 * Creates input/output streams, processes reads in batches, and applies
	 * the configured renaming strategy to each read header.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //4567
		}
		
//		TextStreamWriter tsw=new TextStreamWriter(args[2], false, false, true);
//		tsw.start();
		
		ConcurrentReadOutputStream ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
		ByteBuilder bb=new ByteBuilder();
		
		long x=0;
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				final Read r2=r1.mate;
//				assert(false) : trimBeforeSymbol+", '"+symbol+"'";
				
				if(quantizeQuality){
					final byte[] quals1=r1.quality, quals2=(r2==null ? null : r2.quality);
					Quantizer.quantize(quals1);
					Quantizer.quantize(quals2);
				}
				
				if(trimRight>0 || trimLeft>0 || trimBeforeSymbol>0) {
					if(trimRight>0 || trimLeft>0){
						r1.id=trim(r1.id, trimLeft, trimRight);
						if(r2!=null) {r2.id=trim(r2.id, trimLeft, trimRight);}
					}
					if(trimBeforeSymbol>0){
						r1.id=trimBeforeSymbol(r1.id, trimBeforeSymbol, symbol);
						if(r2!=null) {r2.id=trimBeforeSymbol(r2.id, trimBeforeSymbol, symbol);}
					}
				}else if(fixSRA){
					fixSRA(r1);
					fixSRA(r2);
				}else if(renameByCoords){
					bb.clear().colon().colon().colon();
					ihp.parse(r1);
					ihp.appendCoordinates(bb).space().append(1).colon();
					r1.id=bb.toString();
					if(r2!=null) {
						bb.set(bb.length-2, (byte)'2');
						r2.id=bb.toString();
					}
				}else if(renameByMapping){
					//Should be handled automatically, if output is fastq.
				}else if(r2!=null && (renameByInsert || renameByTrim)){
					
					r1.setMapped(true);
					r2.setMapped(true);
					x=Read.insertSizeMapped(r1, r2, false);
					if(verbose){System.err.println("True Insert: "+x);}
					if(renameByTrim){
						r1.id=r1.numericID+"_"+r1.length()+"_"+Tools.min(x, r1.length())+pairnums[0];
						r2.id=r2.numericID+"_"+r2.length()+"_"+Tools.min(x, r2.length())+pairnums[1];
					}else{
						String s=prefix+x;
						r1.id=s+(addPairnum ? " 1:"+r1.numericID : "");
						if(r2!=null){
							r2.id=s+(addPairnum ? " 2:"+r1.numericID : "");
						}
					}
					
				}else if(prefixOnly){
					r1.id=prefix;
					if(r2!=null){
						r2.id=prefix;
					}
					x++;
				}else if(addPrefix){
					r1.id=prefix+r1.id;
					if(r2!=null){
						r2.id=prefix+r2.id;
					}
					x++;
				}else if(suffix!=null) {
					r1.id+="\t"+suffix;
					if(r2!=null) {r2.id+="\t"+suffix;}
				}else{
					r1.id=prefix+x;
					if(r2!=null){
						r1.id=r1.id+pairnums[0];
						r2.id=prefix+x+pairnums[1];
					}
					x++;
				}
			}
			if(ros!=null){ros.add(reads, ln.id);}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		System.err.println("Time: "+t);
	}
	
	/**
	 * Fixes SRA format read headers by extracting the Illumina portion.
	 * Converts format "SRR1726611.11001 HWI-ST797:117:D091UACXX:4:1101:21093:8249 length=101"
	 * to "HWI-ST797:117:D091UACXX:4:1101:21093:8249" with pair number suffix.
	 * @param r The read to fix, may be null
	 */
	private void fixSRA(Read r){
		//SRR1726611.11001 HWI-ST797:117:D091UACXX:4:1101:21093:8249 length=101
		if(r==null){return;}
		String id=r.id;
		String[] split=spacePattern.split(id);
		assert(split.length==3) : "Unrecognized format: "+id;
		assert(split.length>0 && split[0].indexOf(':')<0) : "Unrecognized format: "+id;
		if(split.length>1){
			r.id=split[1]+pairnums[r.pairnum()];
		}
	}
	
	/**
	 * Trims characters from both ends of a string.
	 * Removes specified number of characters from left and right sides.
	 *
	 * @param s The string to trim
	 * @param left Number of characters to remove from the left
	 * @param right Number of characters to remove from the right
	 * @return Trimmed string, empty string if trimming exceeds string length
	 */
	private static String trim(String s, int left, int right) {
		assert(left>=0 && right>=0) : left+", "+right;
		assert(left>0 || right>0) : left+", "+right;
		int len=s.length()-left-right;
		if(len<0) {return "";}
		else if(len>=s.length()) {return s;}
		return s.substring(left, left+len);
	}
	
	/**
	 * Trims characters before the last occurrence of a symbol.
	 * Removes specified number of characters before the rightmost symbol position.
	 *
	 * @param s The string to trim
	 * @param right Number of characters to remove before the symbol
	 * @param symbol The symbol to search for from the right
	 * @return Modified string with characters removed before the symbol
	 */
	private static String trimBeforeSymbol(String s, int right, char symbol) {
		int pos=s.lastIndexOf(symbol);
//		assert(false) : pos+", '"+symbol+"', "+right;
		if(pos<0) {return s;}
		int len=Tools.max(0, pos-right);
//		assert(false) : pos+", "+len+", "+right;
		String ret=s.substring(0, len)+s.substring(pos);
//		assert(false) : pos+", "+len+", "+right+"\n"+s+"\n"+ret+"\n";
		return ret;
	}
	
	/** Output stream for status messages */
	private PrintStream outstream=System.err;
	
	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path for paired reads */
	private String in2=null;
	
	/** Quality file for primary input */
	private String qfin1=null;
	/** Quality file for secondary input */
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path for paired reads */
	private String out2=null;

	/** Quality output file for primary reads */
	private String qfout1=null;
	/** Quality output file for secondary reads */
	private String qfout2=null;
	
	/** Input file extension override */
	private String extin=null;
	/** Output file extension override */
	private String extout=null;

	/** Prefix to add to read names */
	private String prefix=null;
	/** Suffix to add to read names */
	private String suffix=null;
	
	/** File format for primary input */
	private final FileFormat ffin1;
	/** File format for secondary input */
	private final FileFormat ffin2;

	/** File format for primary output */
	private final FileFormat ffout1;
	/** File format for secondary output */
	private final FileFormat ffout2;

	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	/** Whether to enable verbose output */
	private boolean verbose=false;
	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Flag indicating whether an error occurred during processing */
	public boolean errorState=false;

	/** Number of characters to trim before a specified symbol */
	public int trimBeforeSymbol=0;
	/** Symbol used for trimBeforeSymbol operation */
	public char symbol;
	/** Number of characters to trim from right side of read names */
	public int trimRight=0;
	/** Number of characters to trim from left side of read names */
	public int trimLeft=0;

	/** Whether to add underscore after prefix */
	public boolean addUnderscore=true;
	/** Whether to rename reads based on mapping information */
	public boolean renameByMapping=false;
	/** Whether to rename paired reads using insert size */
	public boolean renameByInsert=false;
	/** Whether to rename reads using numeric ID, length, and insert size */
	public boolean renameByTrim=false;
	/**
	 * Whether to rename reads using coordinate information from Illumina headers
	 */
	public boolean renameByCoords=false;
	/** Whether to add prefix to existing read names instead of replacing them */
	public boolean addPrefix=false;
	/** Whether to replace read names entirely with just the prefix */
	public boolean prefixOnly=false;
	/** Whether to fix SRA-format headers by extracting Illumina portion */
	public boolean fixSRA=false;
	/** Whether to add pair number suffixes to paired reads */
	public boolean addPairnum=true;
	/** Whether to quantize quality scores */
	public boolean quantizeQuality=false;
	/** Pair number suffixes for read 1 and read 2 */
	private String[] pairnums={" 1:", " 2:"};

	/** Pattern for splitting on whitespace characters */
	private static final Pattern spacePattern=Pattern.compile("\\s+");
	/** Pattern for splitting on single space characters */
	private static final Pattern whitespacePattern=Pattern.compile(" ");
	
}
