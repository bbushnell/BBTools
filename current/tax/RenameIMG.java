package tax;

import java.io.File;
import java.io.PrintStream;

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
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.IntHashSet;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class RenameIMG {
	
	/** Program entry point for IMG sequence renaming.
	 * @param args Command-line arguments specifying input files and options */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameIMG x=new RenameIMG(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs RenameIMG instance and parses command-line arguments.
	 * Configures input/output files, IMG database paths, and processing options.
	 * Sets up file format handlers and validates output directory permissions.
	 * @param args Command-line arguments containing file paths and configuration
	 */
	public RenameIMG(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("img")){
				imgFile=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}//TODO: why are these set to the same default?
		if("auto".equalsIgnoreCase(in1)){in1=TaxTree.defaultImgFile();}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FA, null, true, overwrite, append, false);
	}
	
	/**
	 * Processes array of IMG records and copies sequences with renamed headers.
	 * For each IMG record, extracts taxonomy ID and processes the associated FASTA file.
	 * Tracks known and unknown taxonomy IDs during processing.
	 * @param array Array of IMG records containing file paths and metadata
	 */
	void copyFiles(ImgRecord[] array){
		if(useSet){set=new IntHashSet(10000);}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();
		for(ImgRecord ir : array){
			if(ir.taxID>0){set.add(ir.taxID);}
			else{unknownTaxid++;}
			FileFormat ffin=FileFormat.testInput(ir.path(), FileFormat.FA, null, true, true);
			process_inner(ffin, bsw, ir.imgID);
		}
		knownTaxid=set.size();
		set=null;
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
	}
	
	/**
	 * Main processing pipeline that loads IMG database and renames sequences.
	 * Loads IMG records from input file, initializes taxonomy tree, processes files,
	 * and prints comprehensive statistics about the operation.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){
		ImgRecord[] array=ImgRecord.toArray(in1, TaxTree.IMG_HQ);
		if(imgFile==null){
			TaxTree.loadIMG(array);
		}else{
			ImgRecord[] array2=ImgRecord.toArray(imgFile, TaxTree.IMG_HQ);
			TaxTree.loadIMG(array2);
		}
		
		copyFiles(array);
		
		t.stop();

		final int spaces=8;
		String fpstring=""+filesProcessed;
		String cpstring=Tools.padKMB(sequencesProcessed, spaces);
		String bapstring=Tools.padKMB(basesProcessed, spaces);
		String tpstring=""+knownTaxid;
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Files Processed:    "+fpstring);
		outstream.println("Contigs Processed:  "+cpstring);
		outstream.println("Bases Processed:    "+bapstring);
		if(useSet){outstream.println("TaxIDs Processed:   "+tpstring+" \t"+"("+unknownTaxid+" unknown)");}
		outstream.println(Tools.linesBytesProcessed(t.elapsed, linesProcessed, bytesProcessed, spaces));
		
		outstream.println();
		outstream.println("Valid Files:       \t"+filesValid);
		outstream.println("Invalid Files:     \t"+(filesProcessed-filesValid));
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Processes individual FASTA file and renames sequence headers with taxonomy information.
	 * Reads each line from input file, identifies sequence headers (starting with '>'),
	 * and prepends "tid|{taxid}|img|{imgid}" prefix to create taxonomically annotated headers.
	 *
	 * @param ffin Input file format specification
	 * @param bsw Output stream writer for processed sequences
	 * @param img IMG identifier for the current genome
	 */
	void process_inner(final FileFormat ffin, final ByteStreamWriter bsw, final long img){
		
		filesProcessed++;
		{
			File f=new File(ffin.name());
			if(!f.exists() || !f.canRead()){
				System.err.println("Can't find "+f);
				errorState=true;
				return;
			}
		}
		final int tid=TaxTree.imgToTaxid(img);
		ByteFile bf=ByteFile.makeByteFile(ffin);
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;

				linesValid++;
				if(line[0]=='>'){
					sequencesProcessed++;
					bb.append('>');
					if(tid>=0){
						bb.append("tid|");
						bb.append(tid);
						bb.append('|');
					}
					bb.append("img|");
					bb.append(img);
					bb.append(' ');
					for(int i=1; i<line.length; i++){
						bb.append(line[i]);
					}
				}else{
					basesProcessed+=line.length;
					bb.append(line);
				}
				bb.nl();
				bsw.print(bb.toBytes());
				bb.clear();
			}
			line=bf.nextLine();
		}
		
		filesValid++;
		errorState|=bf.close();
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file path containing IMG record list or sequence data */
	private String in1=null;
	/** Output file path for sequences with renamed headers */
	private String out1=null;
	/** Optional IMG database file path for taxonomy information lookup */
	private String imgFile=null;
	
	/*--------------------------------------------------------------*/
	
	/** Hash set for tracking unique taxonomy IDs encountered during processing */
	private IntHashSet set=null;
	/** Count of sequences with known taxonomy IDs */
	private int knownTaxid=0;
	/** Count of sequences with unknown or invalid taxonomy IDs */
	private int unknownTaxid=0;
	/** Whether to use hash set for tracking unique taxonomy IDs */
	private boolean useSet=true;
	
	/** Total number of lines read from all input files */
	private long linesProcessed=0;
	/** Number of valid lines successfully processed */
	private long linesValid=0;
	/** Total bytes read from all input files including headers and sequences */
	private long bytesProcessed=0;

	/** Total number of sequence bases processed (excludes header lines) */
	private long basesProcessed=0;
	/** Total number of sequences (FASTA records) processed */
	private long sequencesProcessed=0;
	/** Total number of files attempted for processing */
	private long filesProcessed=0;
	/** Number of files successfully processed without errors */
	private long filesValid=0;
	
	/** Maximum number of lines to process before stopping */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	
	/** Output file format specification for renamed sequences */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and statistics */
	private PrintStream outstream=System.err;
	/** Global flag controlling verbosity of diagnostic output */
	public static boolean verbose=false;
	/** Flag indicating whether any errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
