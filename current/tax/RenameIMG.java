package tax;

import java.io.File;
import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.IntHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;

/**
 * Concatenates legacy IMG FASTA files while prefixing every sequence header
 * with an IMG genome ID and, when available, the IMG record's NCBI taxID.
 * The primary {@code in=} file selects which IMG records to read and the
 * optional {@code img=} file can provide a broader taxonomy map for those IDs.
 *
 * @author Brian Bushnell
 * @date May 9, 2016
 */
public class RenameIMG{

	/** Parses arguments, rewrites selected IMG FASTA records, and prints summary statistics.
	 * @param args BBTools-style {@code key=value} arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		RenameIMG x=new RenameIMG(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses command-line flags for an IMG dump input, optional alternate IMG
	 * taxID map, and combined FASTA output.
	 * @param args BBTools-style {@code key=value} arguments
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

		//TODO: Probable bug - missing output passes here, then ByteStreamWriter
		//dereferences the null FileFormat created by FileFormat.testOutput.
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FA, null, true, overwrite, append, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Streams every selected IMG FASTA through one output writer.
	 * Each input path is derived from {@link ImgRecord#path()}.
	 * @param array IMG records to concatenate
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
	 * Loads the selection and taxID-mapping dumps, copies all selected FASTA
	 * files, and reports file, line, contig, base and taxID counts.
	 * @param t Timer started by {@link #main(String[])}
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
	 * Copies one IMG FASTA file into the shared output stream.
	 * Header lines receive {@code >tid|N|img|IMGID } when the IMG ID resolved
	 * to a nonnegative taxID, or {@code >img|IMGID } when it did not.
	 * Non-header lines are copied unchanged.
	 *
	 * @param ffin FASTA file to read
	 * @param bsw Shared output writer
	 * @param img IMG ID assigned to the emitted headers
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
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** IMG dump path selecting records whose FASTA files will be concatenated. */
	private String in1=null;
	/** Required combined FASTA output path; use {@code stdout.fa} for standard output. */
	private String out1=null;
	/** Optional alternate IMG dump path used only to resolve IMG IDs to taxIDs. */
	private String imgFile=null;

	/** Distinct known taxIDs observed in the selected {@link #in1} dump. */
	private IntHashSet set=null;
	/** Number of distinct positive taxIDs in the selected {@link #in1} dump. */
	private int knownTaxid=0;
	/** Number of selected records with no positive taxID in {@link #in1}. */
	private int unknownTaxid=0;
	/** Whether to collect distinct taxID and unknown-taxID summary counts. */
	private boolean useSet=true;

	/** Nonempty FASTA lines copied before the global {@link #maxLines} cap. */
	private long linesProcessed=0;
	/** Nonempty FASTA lines queued for output. */
	private long linesValid=0;
	/** Input bytes from nonempty lines, excluding line terminators. */
	private long bytesProcessed=0;

	/** FASTA base lines' byte count, excluding line terminators. */
	private long basesProcessed=0;
	/** FASTA header lines queued for output. */
	private long sequencesProcessed=0;
	/** IMG records attempted. */
	private long filesProcessed=0;
	/** IMG FASTA files processed through the read loop. */
	private long filesValid=0;

	/** Maximum number of nonempty FASTA lines to copy; zero or negative means unlimited. */
	private long maxLines=Long.MAX_VALUE;

	/** Output FASTA descriptor. */
	private final FileFormat ffout1;

	/** Destination for status messages and summary statistics. */
	private PrintStream outstream=System.err;
	/** Whether verbose diagnostics are enabled in lower-level I/O classes. */
	public static boolean verbose=false;
	/** Set when an input is missing or unreadable, or a stream close reports failure. */
	public boolean errorState=false;
	/** Whether output files may be overwritten. */
	private boolean overwrite=true;
	/** Whether output should append to existing files. */
	private boolean append=false;

}
