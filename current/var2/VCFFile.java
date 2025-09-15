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
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;

/**
 * Loads and parses VCF (Variant Call Format) files into memory.
 * Handles VCF headers, sample names, and variant lines with support for
 * complex variant splitting and scaffold mapping generation.
 * 
 * Used primarily for VCF file comparison and processing operations
 * rather than direct variant calling.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @date January 14, 2017
 */
public class VCFFile {
	
	/**
	 * Main method for standalone VCF file processing.
	 * Supports command-line parsing and file loading with timing statistics.
	 * 
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
		
		VCFFile vf=new VCFFile(in);
		t.stop();
		vf.printTime(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructs a VCFFile from a filename string.
	 * 
	 * @param s Input VCF filename
	 */
	public VCFFile(String s){
		in1=s;
		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, false);
		load();
	}
	
	/**
	 * Constructs a VCFFile from a FileFormat object.
	 * 
	 * @param ff FileFormat specifying the input file
	 */
	public VCFFile(FileFormat ff){
		in1=ff.name();
		ffin1=ff;
		load();
	}
	
	/**
	 * Loads VCF data from the input file. Parses header lines to extract
	 * scaffold information and sample names, then loads variant lines into
	 * a LinkedHashMap for preservation of order and fast lookup.
	 */
	void load(){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		
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
		
		errorState|=bf.close();
	}
	
	/**
	 * Prints timing and statistics information about the loaded VCF file.
	 * 
	 * @param t Timer object containing elapsed time information
	 */
	void printTime(Timer t){
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Header Lines:      \t"+header.size());
		outstream.println("Variant Lines:     \t"+map.size());
	}
	
	/**
	 * Returns all VCF lines as a list, with optional splitting of complex variants.
	 * 
	 * @param simplify Whether to split multi-allelic and complex variants into simple variants
	 * @return List of VCFLine objects
	 */
	public ArrayList<VCFLine> lines(boolean simplify){
		ArrayList<VCFLine> lines=new ArrayList<VCFLine>(map.size());
		for(Entry<VCFLine, VCFLine> e : map.entrySet()){
			VCFLine line=e.getValue();
			if(simplify && (line.isMulti() || line.isComplex())){
				ArrayList<VCFLine> list=line.split(true, true, false);
				lines.addAll(list);
			}else{
				lines.add(line);
			}
		}
		return lines;
	}
	
	/**
	 * Creates a ScafMap from VCF contig header lines.
	 * Parses ##contig=<ID=...> lines to build scaffold mapping.
	 * 
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
	
	/**
	 * Converts VCF header to a formatted string.
	 * 
	 * @return Complete VCF header as string with newlines
	 */
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
	
	/** Returns number of lines processed during loading
	 * @return Lines processed count */
	public long linesProcessed() {
		return linesProcessed;
	}
	
	/** Returns number of bytes processed during loading  
	 * @return Bytes processed count */
	public long bytesProcessed() {
		return bytesProcessed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of lines processed */
	private long linesProcessed=0;
	/** Number of bytes processed */
	private long bytesProcessed=0;
	/** Maximum lines to process (for testing/debugging) */
	private long maxLines=Long.MAX_VALUE;
	
	/** VCF header lines */
	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	/** Sample names extracted from VCF header */
	public ArrayList<String> sampleNames=new ArrayList<String>();
	/** Map of VCF variant lines (preserves order, allows fast lookup) */
	public LinkedHashMap<VCFLine, VCFLine> map=new LinkedHashMap<VCFLine, VCFLine>();
	
	/** Input filename */
	private String in1=null;
	/** Input file format */
	private final FileFormat ffin1;
	
	/*--------------------------------------------------------------*/
	/*----------------       Static Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** VCF column header prefix for identification */
	public static final String CHROM_POS="#CHROM\tPOS\t";
	
	/** Output stream for messages */
	private static PrintStream outstream=System.err;
	/** Verbose output flag */
	public static boolean verbose=false;
	/** Error state flag */
	public boolean errorState=false;
}