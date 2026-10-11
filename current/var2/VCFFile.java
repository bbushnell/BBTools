package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;

/**
 * Loads and parses VCF (Variant Call Format) files into memory.
 * Retains headers and sample names, builds an insertion-ordered variant map,
 * and can construct scaffold mappings from contig declarations. Variant parsing
 * inherits VCFLine's shared trimming/cache settings; this is not a structural validator.
 * Duplicate VCFLine identities replace the stored value without moving the first
 * key's insertion position. Identity is not full record-byte equality.
 *
 * Used primarily for VCF file comparison and processing operations
 * rather than direct variant calling.
 * Instances expose mutable collections and are not synchronized. Loading can install
 * the process-wide default scaffold map when absent, then asserts that it is nonempty.
 * Reader errors are retained in errorState; callers must inspect it (VAR2-046).
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date January 14, 2017
 */
public class VCFFile{
	
	/**
	 * Main method for standalone VCF file processing.
	 * Enables shared compression and VCFLine caching, forces BF1 only when BF2 was not forced,
	 * and prints load statistics. Positive lines= limits nonempty rows including headers;
	 * zero or a negative value requests unlimited loading. Recorded input errors fail main.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		VCFLine.AUTOCACHE=true;
		
		Timer t=new Timer();
		String in=null;
		long maxLines=Long.MAX_VALUE;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("lines")){
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
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			in=parser.in1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}
		
		VCFFile vf=new VCFFile(in, maxLines);//VAR2-044: install the limit before loading.
		t.stop();
		vf.printTime(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
		if(vf.errorState){throw new RuntimeException("Error reading or closing VCF input: "+in);}//VAR2-046
	}
	
	/** Constructs a VCFFile from an input filename and immediately loads it.
	 * Uses TEXT as the descriptor fallback and allows subprocess decompression.
	 * @param s Input VCF filename */
	public VCFFile(String s){
		this(s, Long.MAX_VALUE);
	}

	/** Installs main's row limit before loading; public filename construction remains unlimited.
	 * @param s Input VCF filename
	 * @param maxLines_ Positive nonempty-row limit including headers, or zero for unlimited */
	private VCFFile(String s, final long maxLines_){
		assert(maxLines_>=0) : "main normalizes negative limits to Long.MAX_VALUE before constructing VCFFile: "+maxLines_;
		in1=s;
		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, false);
		maxLines=maxLines_;
		load();
	}
	
	/** Borrows a FileFormat and immediately loads its input; the descriptor is not copied.
	 * @param ff FileFormat specifying the input file */
	public VCFFile(FileFormat ff){
		in1=ff.name();
		ffin1=ff;
		load();
	}
	
	/**
	 * Loads all nonempty rows, retaining headers even after variant data. Repeated
	 * column headers append sample names without deduplication. An absent global
	 * scaffold map is initialized from the leading headers only; an existing map is reused.
	 * Retains variants under VCFLine identity, replacing duplicate values in place.
	 * Repeated calls append to existing collections/counters rather than clearing them.
	 * Closes owned input in finally. Normal read/close errors accumulate in errorState;
	 * an earlier parsing failure remains primary if closing also throws.
	 */
	void load(){
		
		final ByteFile bf=ByteFile.makeByteFile(ffin1);
		Throwable failure=null;
		try{
			byte[] line=bf.nextLine();

			// Process header lines first
			while(line!=null && (line.length==0 || line[0]=='#')){
				if(line.length>0){
					if(maxLines>0 && linesProcessed>=maxLines){break;}
					linesProcessed++;
					bytesProcessed+=line.length;
					header.add(line);
					// Extract sample names from column header line
					if(Tools.startsWith(line, CHROM_POS)){
						String[] split=new String(line).split("\t");
						for(int i=9; i<split.length; i++){
							sampleNames.add(split[i]);
						}
					}
				}
				line=bf.nextLine();
			}

			// VAR2-044: A positive header-only limit can return before any contig exists.
			if(maxLines>0 && linesProcessed>=maxLines){return;}

			// Initialize default scaffold map from VCF header if needed
			if(ScafMap.defaultScafMap()==null){ScafMap.setDefaultScafMap(toScafMap(null), in1);}
			assert(ScafMap.defaultScafMap().size()>0) : ScafMap.defaultScafMap()+"\n"+headerToString();

			// Process remaining lines (should be variant data)
			while(line!=null){
				if(line.length>0){
					if(maxLines>0 && linesProcessed>=maxLines){break;}
					linesProcessed++;
					bytesProcessed+=line.length;
					final boolean isHeader=(line[0]=='#');
					if(isHeader){
						header.add(line);
						// Handle additional header lines that might appear after variants
						if(Tools.startsWith(line, CHROM_POS)){
							String[] split=new String(line).split("\t");
							for(int i=9; i<split.length; i++){
								sampleNames.add(split[i]);
							}
						}
					}else{
						// Parse variant line and store in map
						VCFLine vline=new VCFLine(line);
						map.put(vline, vline);
					}
				}
				line=bf.nextLine();
			}
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			//VAR2-045: close on every exit while retaining the original parsing failure.
			try{errorState|=bf.close();}
			catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}
		//TODO: Probable bug [VAR2-046] - GradeVCF.loadTruthSet and ApplyVariants.process ignore errorState;
		//main now rejects it and CompareVCF.getSet accumulates it. Their caller coverage remains separate.
	}
	
	/** Prints timing, header count, and variant count statistics for the loaded VCF.
	 * @param t Timer containing elapsed time */
	void printTime(Timer t){
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Header Lines:      \t"+header.size());
		outstream.println("Variant Lines:     \t"+map.size());
	}
	
	/**
	 * Returns map values in first-key insertion order, borrowing unsplit VCFLine objects.
	 * Ordinary multi-allelic entries request allele splitting alone. Complex entries
	 * still request the unsupported complex mode (VAR2-047); success is not guaranteed.
	 * @param simplify Request splitting for multi-allelic or complex entries
	 * @return New list containing retained values or the dependency's split results
	 */
	public ArrayList<VCFLine> lines(boolean simplify){
		ArrayList<VCFLine> lines=new ArrayList<VCFLine>(map.size());
		for(Entry<VCFLine, VCFLine> e : map.entrySet()){
			VCFLine line=e.getValue();
			if(simplify && (line.isMulti() || line.isComplex())){
				//VAR2-047: ordinary MULTI entries do not require the unsupported complex mode.
				//TODO: Probable bug [VAR2-047] - complex mode still asserts; null split results are not handled here.
				ArrayList<VCFLine> list=line.split(true, line.isComplex(), false);
				lines.addAll(list);
			}else{
				lines.add(line);
			}
		}
		return lines;
	}
	
	/**
	 * Creates a ScafMap from VCF contig header lines.
	 * Delegates exact {@code ##contig=<ID=} prefixes to ScafMap.addFromVcf; adds to the supplied
	 * map without clearing it. Does not itself install the returned map as the default.
	 * @param sm Existing ScafMap to add to (null to create new one)
	 * @return ScafMap containing scaffold information from VCF header
	 */
	public ScafMap toScafMap(ScafMap sm){
		if(sm==null){sm=new ScafMap();}
		for(byte[] line : header){
			if(Tools.startsWith(line, "##contig=<ID=")){
				sm.addFromVcf(line);
			}
		}
		return sm;
	}
	
	/** Renders stored headers with one LF each by casting bytes directly to chars; not UTF-8 decoding. */
	public String headerToString(){
		StringBuilder sb=new StringBuilder();
		for(byte[] line : header){
			for(byte b : line){
				sb.append((char)b);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns nonempty rows consumed, including headers and duplicate variant records.
	 * @return Nonempty input row count */
	public long linesProcessed(){
		return linesProcessed;
	}
	
	/** Returns consumed nonempty row bytes, excluding line terminators removed by the reader.
	 * @return Retained line-byte lengths summed before duplicate removal */
	public long bytesProcessed(){
		return bytesProcessed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Nonempty consumed rows, including headers and duplicates. */
	private long linesProcessed=0;
	/** Sum of consumed row lengths, excluding line terminators. */
	private long bytesProcessed=0;
	/** Positive row limit including headers; zero disables the limit. Public constructors use Long.MAX_VALUE. */
	private long maxLines=Long.MAX_VALUE;
	
	/** Mutable retained header arrays in encounter order, without line terminators. */
	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	/** Mutable sample-name list; every matching column header appends names. */
	public ArrayList<String> sampleNames=new ArrayList<String>();
	/**
	 * Mutable map keyed by VCFLine equality: first equal key keeps its position; latest equal record becomes
	 * the value. Identity fields of retained keys must not change while used for lookup.
	 */
	public LinkedHashMap<VCFLine, VCFLine> map=new LinkedHashMap<VCFLine, VCFLine>();
	
	/** Input name used when installing the default scaffold map. */
	private String in1=null;
	/** Borrowed input descriptor, opened by load. */
	private final FileFormat ffin1;
	
	/*--------------------------------------------------------------*/
	/*----------------       Static Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Exact column-header prefix recognized when extracting sample names. */
	public static final String CHROM_POS="#CHROM\tPOS\t";
	
	/** Shared diagnostics destination selected by main's PreParser. */
	private static PrintStream outstream=System.err;
	/** Main's parsed verbosity, also propagated to shared reader classes. */
	public static boolean verbose=false;
	/** Accumulated normal input-close error flag; constructors currently return despite it. */
	public boolean errorState=false;
}
