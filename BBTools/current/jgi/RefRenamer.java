package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Converts reference sequence names in SAM files.
 * Updates both @SQ header lines and RNAME/RNEXT fields in alignment records.
 * @author Brian Bushnell
 * @author Isla
 * @date July 8, 2025
 *
 */
public class RefRenamer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		RefRenamer x=new RefRenamer(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public RefRenamer(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			in1=parser.in1;
			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.SAM, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		ffMapping=FileFormat.testInput(mappingFile, FileFormat.TXT, null, true, true);

		loadMapping();
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse arguments from the command line */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Set any necessary Parser defaults here
		parser.out1="stdout";

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("mapping") || a.equals("map") || a.equals("ref")){
				mappingFile=b;
			}else if(a.equals("strict")){
				strict=Parse.parseBoolean(b);
			}else if(a.equals("invert") || a.equals("swap")){
				invert=Parse.parseBoolean(b);
			}else if(a.equals("lines")){
				maxLines=Parse.parseKMG(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(i<2 && new File(arg).isFile()){
				if(i==0 && parser.in1==null) {parser.in1=arg;}
				else if(i==1 && mappingFile==null) {mappingFile=arg;}
				else {throw new RuntimeException("Unknown parameter "+arg);}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		return parser;
	}

	/** Load the reference name mapping from TSV file */
	private void loadMapping(){
		if(mappingFile==null){
			outstream.println("Warning: No mapping file specified. No reference names will be changed.");
			return;
		}

		if(verbose){outstream.println("Loading reference mapping from "+mappingFile);}

		ByteFile bf=ByteFile.makeByteFile(ffMapping);
		LineParser1 lp=new LineParser1('\t');
		byte[] line=bf.nextLine();

		while(line!=null){
			if(line.length>0 && line[0]!='#'){
				lp.set(line);
				if(lp.terms()>=2){
					String oldName=lp.parseString(0);
					String newName=lp.parseString(1);
					if(!invert) {refMap.put(oldName, newName);}
					else {refMap.put(newName, oldName);}
					if(verbose){
						outstream.println("Mapping: "+oldName+" -> "+newName);
					}
				}
			}
			line=bf.nextLine();
		}

		bf.close();
		if(verbose){outstream.println("Loaded "+refMap.size()+" reference mappings.");}
	}

	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		mappingFile=Tools.fixExtension(mappingFile);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, mappingFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, mappingFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}

	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
		if(mappingFile==null && strict){
			throw new RuntimeException("Error: strict mode requires a mapping file.");
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create streams and process all data */
	void process(Timer t){

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);

		processInner(bf, bsw);

		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		t.stop();

		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));

		outstream.println();
		outstream.println("Lines Processed:   \t"+linesProcessed);
		outstream.println("Headers Processed: \t"+headersProcessed);
		outstream.println("Headers Converted: \t"+headersConverted);
		outstream.println("Records Converted: \t"+recordsConverted);
		if(unknownRefs.size()>0){
			outstream.println("Unknown References:\t"+unknownsProcessed);
			outstream.println("Unique Unknowns:   \t"+unknownRefs.size());
			for(String ref : unknownRefs){outstream.println(ref);}
		}

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Process the input file line by line, converting reference names based on file type.
	 * Automatically detects FASTA vs SAM format and processes accordingly.
	 * @param bf Input ByteFile to read from
	 * @param bsw Output ByteStreamWriter to write to
	 */
	private void processInner(ByteFile bf, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		LineParser1 lp=new LineParser1('\t');

		final boolean fasta=FileFormat.hasFastaExtension(in1);
		final boolean vcf=FileFormat.hasVcfExtension(in1);
		final boolean gff=FileFormat.hasGffExtension(in1);
		final boolean sam0=FileFormat.hasSamOrBamExtension(in1);
		final boolean sam=!(fasta || vcf || gff);
		if(sam) {
			if(sam0) {outstream.println("Processing sam file.");}
			else {outstream.println("Processing unknown file as sam.");}
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);

				bb.clear();
				if(line[0]=='@'){processSamHeaderLine(line, bb, lp);}
				else{processAlignmentLine(line, bb, lp);}
				bsw.print(bb);
			}
		}else if(fasta){
			outstream.println("Processing fasta file.");
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				bsw.print(processFastaLine(line, bb.clear()));
			}
		}else if(vcf){
			outstream.println("Processing vcf file.");
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				bsw.print(processVcfLine(line, bb.clear(), lp));
			}
		}else if(gff){
			outstream.println("Processing gff file.");
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				bsw.print(processGffLine(line, bb.clear(), lp));
			}
		}else {throw new RuntimeException("Unsupported file type: "+in1);}
	}

	/**
	 * Process a SAM header line, converting @SQ reference names if mappings exist.
	 * Other header lines pass through unchanged.
	 * @param line Raw header line bytes
	 * @param bb ByteBuilder for output construction
	 * @param lp LineParser for field parsing
	 */
	private void processSamHeaderLine(byte[] line, ByteBuilder bb, LineParser1 lp){
		if(Tools.startsWith(line, "@SQ")){
			headersProcessed++;
			//Process @SQ header line - need to convert SN: field
			lp.set(line);
			bb.append("@SQ");

			for(int i=1; i<lp.terms(); i++){
				bb.tab();
				String field=lp.parseString(i);
				if(field.startsWith("SN:")){
					String oldRef=field.substring(3);
					String newRef=refMap.get(oldRef);
					if(newRef!=null){
						bb.append("SN:").append(newRef);
						headersConverted++;
					}else{
						bb.append(field); //Keep original if no mapping
						handleUnknownRef(oldRef);
					}
				}else{
					bb.append(field);
				}
			}
			bb.nl();
		}else{
			//Other header lines pass through unchanged
			bb.append(line).nl();
		}
	}

	/**
	 * Process a SAM alignment line, converting RNAME and RNEXT fields if mappings exist.
	 * Maintains all other fields unchanged while tracking conversion statistics.
	 * @param line Raw alignment line bytes  
	 * @param bb ByteBuilder for output construction
	 * @param lp LineParser for field parsing
	 */
	private void processAlignmentLine(byte[] line, ByteBuilder bb, LineParser1 lp){
		lp.set(line);
		if(lp.terms()<11){
			//Invalid SAM line, pass through unchanged
			bb.append(line).nl();
			return;
		}

		//Get original RNAME and RNEXT
		String rname=lp.parseString(2);
		String rnext=lp.parseString(6);

		//Convert if mappings exist
		String newRname=refMap.getOrDefault(rname, rname);
		String newRnext=refMap.getOrDefault(rnext, rnext);

		//Track unknown references
		if(!rname.equals("*") && !refMap.containsKey(rname)){
			handleUnknownRef(rname);
		}
		if(!rnext.equals("*") && !rnext.equals("=") && !refMap.containsKey(rnext)){
			handleUnknownRef(rnext);
		}

		//Rebuild the line with converted references
		lp.appendTerm(bb, 0).tab();      // QNAME
		lp.appendTerm(bb, 1).tab();      // FLAG
		bb.append(newRname).tab();       // RNAME - converted
		lp.appendTerm(bb, 3).tab();      // POS
		lp.appendTerm(bb, 4).tab();      // MAPQ
		lp.appendTerm(bb, 5).tab();      // CIGAR
		bb.append(newRnext).tab();       // RNEXT - converted
		lp.appendTerm(bb, 7);            // PNEXT

		//Add remaining fields (TLEN, SEQ, QUAL, and optional fields)
		for(int i=8; i<lp.terms(); i++){
			bb.tab();
			lp.appendTerm(bb, i);
		}
		bb.nl();

		//Track conversion
		if(!newRname.equals(rname) || !newRnext.equals(rnext)){
			recordsConverted++;
		}
	}
	
	/**
	 * Process a FASTA header line, converting reference names while preserving descriptions.
	 * Handles both full header replacement and partial replacement with whitespace.
	 * Sequence lines pass through unchanged.
	 * @param line Raw FASTA line bytes
	 * @param bb ByteBuilder for output construction
	 * @return Updated ByteBuilder with processed line
	 */
	private ByteBuilder processFastaLine(byte[] line, ByteBuilder bb){
		final boolean header=Tools.startsWith(line, ">");
		if(!header || line.length<2){
			//Could assert that it starts with a valid symbol but lots of fastas are misformatted
			return bb.append(line).nl();
		}
		headersProcessed+=(header ? 1 : 0);
		
		int limit=Tools.indexOfWhitespace(line);
		if(limit<0) {limit=line.length;}
		String oldRef=new String(line, 1, line.length-1);
		String newRef=(refMap.get(oldRef));
		String name=(newRef!=null ? newRef : limit==line.length ? oldRef : null);
		if(name!=null) {
			if(name==newRef) {headersConverted++;}
			else {handleUnknownRef(oldRef);}
			return bb.append('>').append(name).nl();
		}
		
		oldRef=new String(line, 1, limit);
		newRef=(refMap.get(oldRef));
		name=(newRef!=null ? newRef : oldRef);
		if(name==newRef) {headersConverted++;}
		else {handleUnknownRef(oldRef);}
		return bb.append('>').append(name).append(line, limit, line.length-limit).nl();
	}
	
	/**
	 * Process a VCF line, converting chromosome names in headers and data records.
	 * Handles ##contig header lines and CHROM field in variant records.
	 * @param line Raw VCF line bytes
	 * @param bb ByteBuilder for output construction
	 * @param lp LineParser for field parsing
	 * @return Updated ByteBuilder with processed line
	 */
	private ByteBuilder processVcfLine(byte[] line, ByteBuilder bb, LineParser1 lp){
	    if(Tools.startsWith(line, "##contig=")){
	        //Handle ##contig=<ID=chr1,length=249250621> lines
			headersProcessed++;
	        String lineStr=new String(line);
	        String oldStr=lineStr;
	        for(String oldRef : refMap.keySet()){
	            String newRef=refMap.get(oldRef);
	            lineStr=lineStr.replace("ID="+oldRef+",", "ID="+newRef+",");
	        }
	        bb.append(lineStr).nl();
	        if(!lineStr.equals(oldStr)){headersConverted++;}
	    }else if(line.length>0 && line[0]!='#'){
	        //Handle variant data lines - convert CHROM field (field 0)
	        lp.set(line);
	        if(lp.terms()>=8){
	            String chrom=lp.parseString(0);
	            String newChrom=refMap.getOrDefault(chrom, chrom);
	            
	            bb.append(newChrom); //CHROM - converted
	            for(int i=1; i<lp.terms(); i++){
	                bb.tab();
	                lp.appendTerm(bb, i);
	            }
	            bb.nl();
	            
	            if(!newChrom.equals(chrom)){recordsConverted++;}
	            if(!refMap.containsKey(chrom)){handleUnknownRef(chrom);}
	        }else{
	            bb.append(line).nl(); //Invalid line, pass through
	        }
	    }else{
	        bb.append(line).nl(); //Other headers pass through
	    }
	    return bb;
	}
	
	/**
	 * Process a GFF line, converting sequence names in annotation records.
	 * Comment lines and headers pass through unchanged.
	 * Converts the seqname field (field 0) in feature records.
	 * @param line Raw GFF line bytes
	 * @param bb ByteBuilder for output construction  
	 * @param lp LineParser for field parsing
	 * @return Updated ByteBuilder with processed line
	 */
	private ByteBuilder processGffLine(byte[] line, ByteBuilder bb, LineParser1 lp){
	    if(line.length==0 || line[0]=='#'){
	        //Comments and headers pass through unchanged
	        return bb.append(line).nl();
	    }
	    
	    lp.set(line);
	    if(lp.terms()>=9){
	        //Valid GFF line with all required fields
	        String seqname=lp.parseString(0);
	        String newSeqname=refMap.getOrDefault(seqname, seqname);
	        
	        bb.append(newSeqname);  //seqname - converted
	        for(int i=1; i<lp.terms(); i++){
	            bb.tab();
	            lp.appendTerm(bb, i);
	        }
	        bb.nl();
	        
	        //Track conversions and unknowns
	        if(!newSeqname.equals(seqname)){recordsConverted++;}
	        if(!seqname.equals(".") && !refMap.containsKey(seqname)){
	            handleUnknownRef(seqname);
	        }
	    }else{
	        //Invalid GFF line, pass through unchanged
	        bb.append(line).nl();
	    }
	    
	    return bb;
	}

	/**
	 * Handle unknown reference names by tracking them and optionally warning or crashing.
	 * Each unknown reference is only reported once to avoid spam.
	 * @param ref Unknown reference name encountered
	 */
	private void handleUnknownRef(String ref){
		unknownsProcessed++;
		if(!unknownRefs.contains(ref)){
			unknownRefs.add(ref);
			if(strict){
				throw new RuntimeException("Unknown reference: "+ref+". Use strict=false to continue.");
			}else if(verbose){
				outstream.println("Warning: Unknown reference: "+ref);
			}
		}
	}

	/**
	 * Create a ByteStreamWriter for the specified FileFormat.
	 * @param ff FileFormat to create writer for
	 * @return ByteStreamWriter or null if ff is null
	 */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Reference mapping file path */
	private String mappingFile=null;

	/** Crash on unknown references */
	private boolean strict=false;

	/** Flip the order of the map file names */
	private boolean invert=false;

	/*--------------------------------------------------------------*/

	private final HashMap<String,String> refMap=new HashMap<String,String>();
	private final HashSet<String> unknownRefs=new HashSet<String>();

	private long linesProcessed=0;
	private long bytesProcessed=0;
	private long headersProcessed=0;
	private long unknownsProcessed=0;
	private long headersConverted=0;
	private long recordsConverted=0;

	private long maxLines=Long.MAX_VALUE;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input File */
	private final FileFormat ffin1;
	/** Output File */
	private final FileFormat ffout1;
	/** Mapping File */
	private final FileFormat ffMapping;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;

}