package tax;

import java.io.PrintStream;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
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
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tracker.ReadStats;

/**
 * Compresses and optimizes accession number storage and assembly summaries.
 * NCBI four-column accession-to-taxID rows drop accession.version and optionally
 * GI numbers; two-column accession.version rows retain the full first column and
 * insert an empty second output column. Assembly mode keeps accessions, taxIDs
 * and selected assembly metrics in the same three-leading-column layout.
 *
 * @author Brian Bushnell
 * @date April 4, 2017
 */
public class ShrinkAccession{
	
	/** Runs accession-table shrinking from command-line arguments.
	 * @param args BBTools-style arguments selecting input, output and mode */
	public static void main(String[] args){
		Timer t=new Timer();
		ShrinkAccession x=new ShrinkAccession(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/** Parses command-line arguments and validates input/output files.
	 * Positional arguments are treated as input then output when the parser has
	 * not already consumed those paths.
	 * @param args Command-line arguments in BBTools key=value form */
	public ShrinkAccession(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		if(Data.PIGZ()){
			ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
		}
		
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
				ReadWrite.verbose=verbose;
			}else if(a.equals("seq") || a.equals("sequence")){
				mode=(Parse.parseBoolean(b) ? SEQUENCE : ASSEMBLY);
			}else if(a.equals("asm") || a.equals("assembly")){
				mode=(!Parse.parseBoolean(b) ? SEQUENCE : ASSEMBLY);
			}else if(a.equals("gi")){
				KEEP_GI_NUMBERS=Parse.parseBoolean(b);
			/* tax/ShrinkAccession cleanup: removed the dead
			 * outgi/giout/gi->giOut->ffoutGi cluster; its gi alias was
			 * shadowed by the branch above, and ffoutGi was never written. */
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in=parser.in1;

			out=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		
	}
	
	/** Converts the configured accession input into the smaller output form.
	 * Attempts reader closure and output cleanup on normal and exceptional
	 * row-processing exits. Writers that fail before worker startup are closed
	 * directly; writers with a started worker are poisoned and joined.
	 * @param t Timer for throughput reporting */
	void process(Timer t){
		
		ByteFile bf=null;
		ByteStreamWriter bsw=null;
		boolean bswStarted=false;
		Throwable failure=null;
		try{
			bf=ByteFile.makeByteFile(ffin);
			bsw=makeWriter(ffout);
			bsw.start();
			bswStarted=true;

			if(mode==SEQUENCE){
				processSeq(bf, bsw);
			}else if(mode==ASSEMBLY){
				processAsm(bf, bsw);
			}else{
				throw new RuntimeException("Unknown mode: "+mode);
			}
		}catch(RuntimeException|Error e){
			failure=e;
		}finally{
			if(bf!=null){
				try{
					errorState|=bf.close();
				}catch(RuntimeException|Error e){
					failure=recordCleanupFailure(failure, e);
				}
			}
			if(bsw!=null){
				try{
					errorState|=closeWriter(bsw, bswStarted);
				}catch(RuntimeException|Error e){
					failure=recordCleanupFailure(failure, e);
				}
			}
		}
		
		if(failure!=null){rethrow(failure);}
		
		t.stop();
		outstream.println("Discarded "+badLines+" lines.\n");
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Shrinks NCBI accession2taxid rows into accession, taxID and optional GI columns.
	 * @param bf Input accession table
	 * @param bsw Output writer */
	void processSeq(ByteFile bf, ByteStreamWriter bsw){
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder(10000);
		int columns=4;
		while(line!=null){
			if(Tools.startsWith(line, "accession\t")){
				bb.append(line);
				bb.nl();
			}else if(Tools.startsWith(line, "accession.version\ttaxid")){
				columns=2;
				bb.append("accession\t\ttaxid\t");//dummy header
				bb.nl();
			}else{
				charsProcessed+=line.length+1;
				linesProcessed++;
				
				final int tid=(columns==4 ? AccessionToTaxid.parseLineToTaxid(line, (byte)'\t') : 
					AccessionToTaxid.parseLineToTaxid_2col(line, (byte)'\t'));
				if(tid<1){
					badLines++;
				}else{
					int i=0;
					
					while(i<line.length){//Accession
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					
					if(columns==4){
						while(i<line.length){//Accession with decimal
							byte b=line[i];
							i++;
							if(b=='\t'){break;}
						}
					}
					bb.append('\t');
					
					while(i<line.length){//Taxid
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					
					if(KEEP_GI_NUMBERS){
						if(line.length>i && Tools.isDigit(line[i])){//GI number or "na"
							while(i<line.length){
								byte b=line[i];
								bb.append(b);
								i++;
							}
						}
					}
					bb.nl();
				}
			}
			if(bb.length()>8000){
				bsw.print(bb);
				bb.clear();
			}
			line=bf.nextLine();
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
	}
	
	/** Shrinks assembly_summary rows into accession, taxID and assembly metrics.
	 * Copies columns 0 #assembly_accession, 5 taxid, 23 assembly_type,
	 * 25 genome_size, 26 genome_size_ungapped, 28 replicon_count,
	 * 29 scaffold_count, 30 contig_count, 34 total_gene_count,
	 * 35 protein_coding_gene_count and 36 non_coding_gene_count.
	 * @param bf Input assembly_summary table
	 * @param bsw Output writer */
	void processAsm(ByteFile bf, ByteStreamWriter bsw){
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder(20000);
		LineParser1 lp=new LineParser1('\t');
		while(line!=null && Tools.startsWith(line, "#")){
			if(Tools.startsWith(line, "#assembly_accession\t")){
				lp.set(line);
				bb.append(lp.parseByteArray(0)).tab().tab();//accession and blank compatibility column
				bb.append(lp.parseByteArray(5)).tab();//tid
				bb.append(lp.parseByteArray(23)).tab();//ploidy?
				bb.append(lp.parseByteArray(25)).tab();//size
				bb.append(lp.parseByteArray(26)).tab();//size_ungapped
				bb.append(lp.parseByteArray(28)).tab();//chromosomes
				bb.append(lp.parseByteArray(29)).tab();//scafs
				bb.append(lp.parseByteArray(30)).tab();//contigs
				bb.append(lp.parseByteArray(34)).tab();//genes
				bb.append(lp.parseByteArray(35)).tab();//coding
				bb.append(lp.parseByteArray(36));//noncoding
				bb.nl();
			}
			charsProcessed+=line.length+1;
			linesProcessed++;
			line=bf.nextLine();
		}
		while(line!=null){
			charsProcessed+=line.length+1;
			linesProcessed++;

			/* tax/ShrinkAccession#001: comment rows must advance before
			 * continuing; assembly_summary inputs can concatenate mid-stream
			 * headers, and this branch exists to skip them. */
			if(Tools.startsWith(line, "#")){line=bf.nextLine(); continue;}
			lp.set(line);

			bb.append(lp.parseByteArray(0)).tab().tab();//accession
			//Blank column to match seq files
			bb.append(parseNum(lp, 5)).tab();//tid
			bb.append(lp.parseByteArray(23)).tab();//ploidy?
			bb.append(parseNum(lp, 25)).tab();//size
			bb.append(parseNum(lp, 26)).tab();//size_ungapped
			bb.append(parseNum(lp, 28)).tab();//chromosomes
			bb.append(parseNum(lp, 29)).tab();//scafs
			bb.append(parseNum(lp, 30)).tab();//contigs
			bb.append(parseNum(lp, 34)).tab();//genes
			bb.append(parseNum(lp, 35)).tab();//coding
			bb.append(parseNum(lp, 36));//noncoding
			bb.nl();
			if(bb.length()>16000){
				bsw.print(bb);
				bb.clear();
			}
			line=bf.nextLine();
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
	}
	
	/** Parses a numeric assembly-summary column.
	 * @param lp Parser positioned on the current assembly_summary row
	 * @param field Zero-based field index
	 * @return Field value, or -1 when the field starts with a letter */
	private long parseNum(LineParser1 lp, int field){
		long x=-1;
		if(lp.termStartsWithLetter(field)){return -1;}
		x=lp.parseLong(field);
		return x;
	}
	
	/** Creates the output writer; fixtures can inject start failures.
	 * @param ff Output file descriptor
	 * @return Writer for {@code ff} */
	protected ByteStreamWriter makeWriter(final FileFormat ff){return new ByteStreamWriter(ff);}

	/** Closes an output writer whether or not the worker started.
	 * @param bsw Writer whose constructor already opened the output stream
	 * @param started Whether {@link ByteStreamWriter#start()} returned
	 * @return True if writer cleanup reported an error */
	private static boolean closeWriter(final ByteStreamWriter bsw, final boolean started){
		return started || bsw.getState()!=Thread.State.NEW ?
				bsw.poisonAndWait() : bsw.closeBeforeStart();
	}

	/** Retains the first failure and suppresses distinct cleanup failures on it.
	 * @param primary First failure, or null if cleanup is the first failure
	 * @param cleanupFailure Failure raised while closing an owned resource
	 * @return The failure that should be rethrown */
	private static Throwable recordCleanupFailure(Throwable primary, Throwable cleanupFailure){
		if(primary==null){return cleanupFailure;}
		if(cleanupFailure!=primary){primary.addSuppressed(cleanupFailure);}
		return primary;
	}

	/** Rethrows a recorded unchecked failure.
	 * @param failure RuntimeException or Error to rethrow */
	private static void rethrow(Throwable failure){
		if(failure instanceof Error){throw (Error)failure;}
		throw (RuntimeException)failure;
	}

	/*--------------------------------------------------------------*/
	
	/** Input accession or assembly-summary file. */
	private String in=null;
	/** Output shrunk accession table. */
	private String out=null;
	/** Input interpretation mode; either SEQUENCE or ASSEMBLY. */
	private int mode=SEQUENCE;

	/** Nonheader sequence rows parsed, or assembly rows and skipped comments scanned. */
	long linesProcessed=0;
	/** Bytes in rows counted by {@link #linesProcessed}, including one newline per row. */
	long charsProcessed=0;
	/** Sequence-mode rows rejected because their parsed taxID was below one. */
	long badLines=0;
	
	/*--------------------------------------------------------------*/
	
	/** Validated input descriptor. */
	private final FileFormat ffin;
	/** Validated output descriptor. */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	
	/** Destination for status and error messages. */
	private PrintStream outstream=System.err;
	/** Enables verbose file IO logging. */
	public static boolean verbose=false;
	/** Whether sequence-mode output keeps the source GI column. */
	public static boolean KEEP_GI_NUMBERS=true;
	/** NCBI accession2taxid input mode. */
	public static final int SEQUENCE=1, ASSEMBLY=2;
	/** True if input close or output-writer completion reported an error. */
	public boolean errorState=false;
	/** Whether existing output files may be overwritten. */
	private boolean overwrite=true;
	/** Whether existing output files may be appended. */
	private boolean append=false;
	
}
