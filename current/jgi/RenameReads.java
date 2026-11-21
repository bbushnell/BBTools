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
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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
 * Renames sequence read headers while preserving the configured input/output stream lifecycle.
 * <p>
 * The constructor accepts several ID transformations. Quality quantization is independent; ID
 * modes then take this precedence: fixed or symbol trimming, SRA repair, coordinate renaming,
 * FASTQ custom mapping, paired insert-size or trim-length renaming, fixed prefix-only IDs,
 * prefix addition, suffix addition, then default numeric IDs with an optional prefix.
 *
 * @author Brian Bushnell
 * @date Aug 23, 2013
 */
public class RenameReads {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RenameReads x=new RenameReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		if(x.errorState){
			throw new RuntimeException("RenameReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Parses arguments, resolves paired file templates and interleaving, validates input/output
	 * path collisions and output permissions, and constructs the input and output descriptors.
	 *
	 * @param args command-line arguments
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
				// parsed by shared static options
			}else if(Parser.parseZip(arg, a, b)){
				// parsed by compression options
			}else if(Parser.parseQuality(arg, a, b)){
				// parsed by quality options
			}else if(Parser.parseFasta(arg, a, b)){
				// parsed by FASTA options
			}else if(parser.parseInterleaved(arg, a, b)){
				// parsed by interleaving options
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
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
			}else if(a.equals("renamebycoords") || a.equals("renamebycoordinates") || a.equals("coordinates") || a.equals("coords")){
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
			}
			
			renameByMapping=FASTQ.TAG_CUSTOM;
		}
		
		if(addUnderscore && prefix!=null && !prefix.endsWith("_") && !prefixOnly){prefix+="_";}
		if(!addPairnum){pairnums=new String[] {"",""};}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		if(prefix==null || prefix.length()<1){prefix="";}

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
		
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("Input and output file names must be distinct.");
		}

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
	 * Streams reads, applies the selected ID transformation to each list, writes requested output,
	 * and, on normal completion, returns all lists and records close errors in {@link #errorState}.
	 *
	 * @param t timer used for final elapsed-time reporting
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start();
		}
		
		
		ConcurrentReadOutputStream ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

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
		//Main reports retained stream errors after closing the diagnostic stream.
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		System.err.println("Time: "+t);
	}
	
	/**
	 * Converts an SRA-style three-token ID to the embedded instrument ID plus pair suffix.
	 *
	 * @param r read to rename, or {@code null} for an absent mate
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
	 * Trims fixed numbers of characters from both ends of an ID.
	 *
	 * @param s source ID
	 * @param left characters to trim from the left
	 * @param right characters to trim from the right
	 * @return trimmed ID, or an empty ID when the trim consumes the whole ID
	 */
	private static String trim(String s, int left, int right){
		assert(left>=0 && right>=0) : left+", "+right;
		assert(left>0 || right>0) : left+", "+right;
		//[JGI-128] Subtract right alone; comparing left avoids combined-width overflow.
		final int end=s.length()-right;
		return left>=end ? "" : s.substring(left, end);
	}
	
	/**
	 * Removes up to {@code right} characters immediately before the last requested symbol.
	 *
	 * @param s source ID
	 * @param right maximum characters to remove before {@code symbol}
	 * @param symbol delimiter whose last occurrence marks the preserved suffix
	 * @return ID with the shortened prefix, or the original ID when {@code symbol} is absent
	 */
	private static String trimBeforeSymbol(String s, int right, char symbol) {
		int pos=s.lastIndexOf(symbol);
		if(pos<0) {return s;}
		int len=Tools.max(0, pos-right);
		String ret=s.substring(0, len)+s.substring(pos);
		return ret;
	}
	
	// Diagnostic output after PreParser handles redirection.
	private PrintStream outstream=System.err;
	
	// Primary and optional mate input paths.
	private String in1=null;
	private String in2=null;
	
	// Optional input quality-sidecar paths.
	private String qfin1=null;
	private String qfin2=null;

	// Primary and optional mate output paths.
	private String out1=null;
	private String out2=null;

	// Optional output quality-sidecar paths.
	private String qfout1=null;
	private String qfout2=null;
	
	// Format extension overrides.
	private String extin=null;
	private String extout=null;

	// ID text to add or substitute.
	private String prefix=null;
	private String suffix=null;
	
	// Resolved input formats.
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	// Resolved output formats.
	private final FileFormat ffout1;
	private final FileFormat ffout2;

	// Output policy, progress and terminal state.
	private boolean overwrite=true;
	private boolean append=false;
	private boolean verbose=false;
	private long maxReads=-1;
	public boolean errorState=false;

	// ID trimming controls.
	public int trimBeforeSymbol=0;
	public char symbol;
	public int trimRight=0;
	public int trimLeft=0;

	// Renaming modes and mate suffixes.
	public boolean addUnderscore=true;
	public boolean renameByMapping=false;
	public boolean renameByInsert=false;
	public boolean renameByTrim=false;
	public boolean renameByCoords=false;
	public boolean addPrefix=false;
	public boolean prefixOnly=false;
	public boolean fixSRA=false;
	public boolean addPairnum=true;
	public boolean quantizeQuality=false;
	private String[] pairnums={" 1:", " 2:"};

	// Header tokenization patterns.
	private static final Pattern spacePattern=Pattern.compile("\\s+");
	private static final Pattern whitespacePattern=Pattern.compile(" ");
	
}
