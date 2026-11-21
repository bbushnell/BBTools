package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Translates reference names in FASTA, SAM, VCF and GFF text using a tab-separated map.
 * Format-specific helpers select the header or record fields to translate; they do not
 * rewrite every occurrence of a reference name. Unknown-name checks are also format-specific.
 * Output may be disabled while processing, counters and applicable strict checks still run.
 *
 * @author Brian Bushnell
 * @author Isla
 * @date July 8, 2025
 */
public class RefRenamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads the mapping, processes the input and, on normal completion, closes redirected diagnostics.
	 * @param args command-line input/output paths, mapping and processing options
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
	 * Parses options, validates and resolves paths, prepares format descriptors and loads the map.
	 * The main input and optional output streams are created later by {@link #process(Timer)}.
	 * @param args command-line input/output paths, mapping and processing options
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
		fixExtensions(); //Resolve existing input and mapping path variants
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

	/**
	 * Parses tool-specific options and delegates shared options to Parser.
	 * Only the first equals sign separates an option from its value. Output defaults to stdout;
	 * an explicit null value disables it. A negative line limit means no configured limit.
	 * @param args arguments after preprocessing
	 * @return shared parser containing configured input/output options
	 */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Write by default; an explicit out=null still disables output.
		parser.out1="stdout";

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Only the first '=' separates the option; filenames may contain more.
			final int equals=arg.indexOf('=');
			String a=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			String b=(equals<0 || equals==arg.length()-1 ? null : arg.substring(equals+1));
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
				if(i==0 && parser.in1==null){parser.in1=arg;}
				else if(i==1 && mappingFile==null){mappingFile=arg;}
				else{throw new RuntimeException("Unknown parameter "+arg);}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		return parser;
	}

	/**
	 * Loads the first two tab-separated fields of each eligible mapping line.
	 * Skips empty lines, comments beginning with {@code #}, and rows with fewer than two fields.
	 * Inversion swaps key and value before insertion; later duplicate keys replace earlier values.
	 * With no mapping file, leaves the map empty and reports a warning. On normal completion,
	 * closes its owned mapping reader; its close status is not currently propagated to errorState.
	 */
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
					if(!invert){refMap.put(oldName, newName);}
					else{refMap.put(newName, oldName);}
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

	/** Resolves input/mapping paths through Tools.fixExtension and requires a primary input. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		mappingFile=Tools.fixExtension(mappingFile);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/** Validates that input files can be read and output files can be written.
	 * Checks for file permissions, duplicate file specifications, and path validity. */
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

	/** Selects BF2 when neither BF1 nor BF2 is forced and the configured thread count exceeds two. */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}

	/**
	 * Requires a mapping file when strict unknown-name handling is enabled.
	 * @return true when this check passes; otherwise throws an exception
	 */
	private boolean validateParams(){
		if(mappingFile==null && strict){
			throw new RuntimeException("Error: strict mode requires a mapping file.");
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates the input and optional output, processes lines and reports counters and elapsed time.
	 * On normal completion, closes the input and waits for an output writer when present.
	 * This method does not provide exceptional-path cleanup for every processing failure.
	 * @param t timer started by the caller, stopped after normal stream completion
	 */
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

		//Report retained stream errors after the normal completion statistics.
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Processes nonempty lines up to maxLines using extension-selected FASTA, VCF or GFF handling.
	 * Other extensions use the SAM branch. Each helper runs even when output is disabled, so
	 * conversion counters and applicable unknown-name checks do not depend on writing.
	 * The caller owns stream creation and normal completion; this method borrows both streams.
	 * @param bf input reader
	 * @param bsw output writer, or null to process and count without writing
	 */
	private void processInner(ByteFile bf, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		LineParser1 lp=new LineParser1('\t');

		final boolean fasta=FileFormat.hasFastaExtension(in1);
		final boolean vcf=FileFormat.hasVcfExtension(in1);
		final boolean gff=FileFormat.hasGffExtension(in1);
		final boolean sam0=FileFormat.hasSamOrBamExtension(in1);
		final boolean sam=!(fasta || vcf || gff);
		if(sam){
			if(sam0){outstream.println("Processing sam file.");}
			else{outstream.println("Processing unknown file as sam.");}
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);

				bb.clear();
				if(line[0]=='@'){processSamHeaderLine(line, bb, lp);}
				else{processAlignmentLine(line, bb, lp);}
				if(bsw!=null){bsw.print(bb);}
			}
		}else if(fasta){
			outstream.println("Processing fasta file.");
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				processFastaLine(line, bb.clear());
				if(bsw!=null){bsw.print(bb);}
			}
		}else if(vcf){
			outstream.println("Processing vcf file.");
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				processVcfLine(line, bb.clear(), lp);
				if(bsw!=null){bsw.print(bb);}
			}
		}else if(gff){
			outstream.println("Processing gff file.");
			for(byte[] line=bf.nextLine(); line!=null && linesProcessed<maxLines; line=bf.nextLine()){
				if(line.length==0){continue;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				processGffLine(line, bb.clear(), lp);
				if(bsw!=null){bsw.print(bb);}
			}
		}else {throw new RuntimeException("Unsupported file type: "+in1);}
	}

	/**
	 * Rebuilds a SAM header, translating SN fields in recognized {@code @SQ} lines.
	 * Other header lines pass through unchanged. Each mapped SN field increments the header
	 * conversion counter; unmapped SN fields use the shared unknown-name handler.
	 *
	 * @param line nonempty input line
	 * @param bb caller-owned output builder, cleared before this call
	 * @param lp reusable tab-field parser
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
	 * Translates SAM RNAME and RNEXT fields and appends all other fields unchanged.
	 * Lines with fewer than eleven fields pass through. Unknown-name tracking excludes
	 * {@code *} in either reference field and {@code =} in RNEXT; a record counts as converted
	 * when either translated field differs from its original value.
	 *
	 * @param line nonempty input line
	 * @param bb caller-owned output builder, cleared before this call
	 * @param lp reusable tab-field parser
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
	 * Appends a FASTA line, looking up a complete header before its first whitespace-delimited name.
	 * A complete-header mapping replaces its description too; a name-only mapping preserves the
	 * remaining description bytes. Nonheaders and single-character header markers pass through.
	 * An applied mapping increments the header conversion counter even for an identical value.
	 *
	 * @param line input line
	 * @param bb caller-owned output builder, cleared before this call
	 * @return the supplied builder with the processed line and newline appended
	 */
	private ByteBuilder processFastaLine(byte[] line, ByteBuilder bb){
		final boolean header=Tools.startsWith(line, ">");
		if(!header || line.length<2){
			//Could assert that it starts with a valid symbol but lots of fastas are misformatted
			return bb.append(line).nl();
		}
		headersProcessed+=(header ? 1 : 0);

		int limit=Tools.indexOfWhitespace(line);
		if(limit<0){limit=line.length;}
		String oldRef=new String(line, 1, line.length-1);
		String newRef=(refMap.get(oldRef));
		String name=(newRef!=null ? newRef : limit==line.length ? oldRef : null);
		if(name!=null){
			if(name==newRef){headersConverted++;}
			else{handleUnknownRef(oldRef);}
			return bb.append('>').append(name).nl();
		}

		oldRef=new String(line, 1, limit-1);
		newRef=(refMap.get(oldRef));
		name=(newRef!=null ? newRef : oldRef);
		if(name==newRef){headersConverted++;}
		else{handleUnknownRef(oldRef);}
		return bb.append('>').append(name).append(line, limit, line.length-limit).nl();
	}

	/**
	 * Appends a VCF line, translating a recognized contig ID or a record's CHROM field once.
	 * Contig handling looks for the first literal {@code ID=} marker followed by a comma;
	 * it is not a general attribute parser. Headers count only recognized contig lines, and only
	 * changed contig text counts as converted. Unknown-name handling applies to records with at
	 * least eight fields, not headers.
	 *
	 * @param line input line
	 * @param bb caller-owned output builder, cleared before this call
	 * @param lp reusable tab-field parser
	 * @return the supplied builder with the processed line and newline appended
	 */
	private ByteBuilder processVcfLine(byte[] line, ByteBuilder bb, LineParser1 lp){
		if(Tools.startsWith(line, "##contig=")){
			//Handle ##contig=<ID=chr1,length=249250621> lines
			headersProcessed++;
			String lineStr=new String(line);
			String oldStr=lineStr;
			//Look up the original comma-terminated ID once; never remap its replacement.
			final int keyStart=lineStr.indexOf("ID=");
			final int nameStart=keyStart+3;
			final int nameEnd=(keyStart<0 ? -1 : lineStr.indexOf(',', nameStart));
			if(nameEnd>=nameStart){
				final String oldRef=lineStr.substring(nameStart, nameEnd);
				final String newRef=refMap.get(oldRef);
				if(newRef!=null && !newRef.equals(oldRef)){
					lineStr=lineStr.substring(0, nameStart)+newRef+lineStr.substring(nameEnd);
				}
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
	 * Appends a GFF line, translating the sequence field in records with at least nine fields.
	 * Empty lines, comments and shorter records pass through without header/conversion counting.
	 * An unmapped sequence name other than {@code .} uses the shared unknown-name handler.
	 *
	 * @param line input line
	 * @param bb caller-owned output builder, cleared before this call
	 * @param lp reusable tab-field parser
	 * @return the supplied builder with the processed line and newline appended
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
	 * Counts every call for an unknown name and records distinct names for the final report.
	 * On the first occurrence of a name, strict mode throws; otherwise verbose mode warns.
	 * Call sites determine which format fields participate in this check.
	 * @param ref unknown reference name
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
	 * Creates and starts an output writer when an output descriptor is present.
	 * @param ff output descriptor, or null to disable writing
	 * @return started writer owned by the caller, or null when output is disabled
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

	/** Primary text input path. */
	private String in1=null;

	/** Output path; null disables writing without disabling processing. */
	private String out1=null;

	/** Optional tab-separated source/replacement name mapping. */
	private String mappingFile=null;

	/** Throw on the first unknown name encountered by a participating format helper. */
	private boolean strict=false;

	/** Swap mapping keys and values while loading. */
	private boolean invert=false;

	/*--------------------------------------------------------------*/

	/** Owned name map; later insertions replace earlier entries for the same key. */
	private final HashMap<String,String> refMap=new HashMap<String,String>();
	/** Owned set of distinct names passed to handleUnknownRef. */
	private final HashSet<String> unknownRefs=new HashSet<String>();

	/** Nonempty input lines processed, including recognized headers and comments. */
	private long linesProcessed=0;
	/** Sum of processed line lengths plus one byte per processed line. */
	private long bytesProcessed=0;
	/** Recognized FASTA headers, SAM SQ lines and VCF contig lines; excludes GFF comments. */
	private long headersProcessed=0;
	/** Calls to handleUnknownRef, including repeated occurrences of a name. */
	private long unknownsProcessed=0;
	/** Mapped FASTA/SAM header fields or changed VCF contig lines; see the format helpers. */
	private long headersConverted=0;
	/** Records with at least one changed selected reference field. */
	private long recordsConverted=0;

	/** Maximum number of nonempty input lines to process. */
	private long maxLines=Long.MAX_VALUE;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input descriptor. */
	private final FileFormat ffin1;
	/** Output descriptor, or null when writing is disabled. */
	private final FileFormat ffout1;
	/** Mapping input descriptor used when a mapping file is configured. */
	private final FileFormat ffMapping;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostics selected by PreParser; main closes redirected output on normal completion. */
	private PrintStream outstream=System.err;
	/** Shared mapping-detail and unknown-name warning verbosity. */
	public static boolean verbose=false;
	/** Retained main-input/output stream error flag checked after normal completion. */
	public boolean errorState=false;
	/** Whether an existing output may be overwritten. */
	private boolean overwrite=true;
	/** Whether output is appended to an existing destination. */
	private boolean append=false;

}
