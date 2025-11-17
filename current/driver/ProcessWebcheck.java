package driver;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

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
import structures.IntList;

/**
 * @author Brian Bushnell
 * @date May 9, 2016
 *
 */
public class ProcessWebcheck {
	
	/** Program entry point for web log analysis.
	 * @param args Command-line arguments specifying input files and options */
	public static void main(String[] args){
		Timer t=new Timer();
		ProcessWebcheck x=new ProcessWebcheck(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs ProcessWebcheck instance with parsed command-line arguments.
	 * Configures input/output files, processing parameters, and validation settings.
	 * Sets up file formats and validates output file accessibility.
	 * @param args Command-line arguments array containing file paths and options
	 */
	public ProcessWebcheck(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, null, false);
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

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("extendedstats") || a.equals("es")){
				extendedStats=Parse.parseBoolean(b);
			}else if(a.equals("invalid") || a.equals("outinvalid")){
				outInvalid=b;
			}else if(a.equals("fail") || a.equals("outfail")){
				outFail=b;
			}else if(a.equals("ms") || a.equals("millis")){
				boolean x=Parse.parseBoolean(b);
				ms=(x ? "ms" : "");
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
			}else if(!arg.contains("=")){
				for(String s : arg.split(",")){
					in1.add(s);
				}
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			if(parser.in1!=null){
				for(String s : parser.in1.split(",")){
					in1.add(s);
				}
			}
			
			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutFail=FileFormat.testOutput(outFail, FileFormat.TXT, null, true, overwrite, append, false);
		
		ffin1=new ArrayList<FileFormat>(in1.size());
		for(String s : in1){
			FileFormat ff=FileFormat.testInput(s, FileFormat.TXT, null, true, true);
			assert(ff!=null) : "Cannot read file "+s;
			ffin1.add(ff);
		}
		assert(ffin1.size()>0) : "No input files.";
	}
	
	/**
	 * Main processing method that analyzes all input log files.
	 * Creates output streams for results, failed requests, and invalid lines.
	 * Processes each input file sequentially and generates comprehensive statistics.
	 * Outputs aggregated results including latency metrics and failure codes.
	 *
	 * @param t Timer for tracking total execution time
	 */
	void process(Timer t){

		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}

		ByteStreamWriter bswFail=null;
		if(ffoutFail!=null){
			bswFail=new ByteStreamWriter(ffoutFail);
			bswFail.start();
		}
		
		for(FileFormat ff : ffin1){
			ByteFile bf=ByteFile.makeByteFile(ff);
			process2(bf, bswFail, bswInvalid);
		}

		passLatency.shrink();
		failLatency.shrink();
		failCode.sort();
		failCode.shrinkToUnique();

		StringBuilder sb=new StringBuilder();

		ArrayList<String> list=new ArrayList<String>();
		list.addAll(map.keySet());
		Shared.sort(list);
		for(String s : list){
			sb.append(s+"\t"+map.get(s)[0]+"\n");
		}

		if(extendedStats){
			sb.append('\n');
			sb.append("Lines_Processed:\t"
					+ linesProcessed).append('\n');
			sb.append("Invalid_Lines:\t"
					+ (linesProcessed-linesValid)).append('\n');
			sb.append("Passing:\t"
					+ passLatency.size).append('\n');
			sb.append("Failing:\t"
					+ failLatency.size).append('\n');
			sb.append("Avg_Pass_Latency:\t"
					+ (passLatency.size>0 ? Tools.averageInt(passLatency.array) : 0)+ms).append('\n');
			sb.append("Max_Pass_Latency:\t"
					+ (passLatency.size>0 ? Tools.max(passLatency.array) : 0)+ms).append('\n');
			sb.append("Avg_Fail_Latency:\t"
					+ (failLatency.size>0 ? Tools.averageInt(failLatency.array) : 0)+ms).append('\n');
			sb.append("Max_Fail_Latency:\t"
					+ (failLatency.size>0 ? Tools.max(failLatency.array) : 0)+ms).append('\n');
			sb.append("Observed_Fail_Codes:");
			for(int i=0; i<failCode.size; i++){
				sb.append('\t').append(failCode.get(i));
			}
			sb.append('\n');
		}

		outstream.print(sb);

		if(bsw!=null){
			bsw.print(sb);
			errorState|=bsw.poisonAndWait();
		}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		if(bswFail!=null){errorState|=bswFail.poisonAndWait();}

		t.stop();

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
		
	/**
	 * Processes individual log file parsing line-by-line entries.
	 * Validates log format (pipe-separated fields), extracts HTTP codes and latency.
	 * Categorizes entries as passing (200 codes) or failing, tracks statistics.
	 * Writes failed and invalid entries to separate output streams if configured.
	 *
	 * @param bf Input ByteFile containing web check log entries
	 * @param bswFail Output stream for failed requests (non-200 response codes)
	 * @param bswInvalid Output stream for malformed/unparseable log lines
	 */
	private void process2(ByteFile bf, ByteStreamWriter bswFail, ByteStreamWriter bswInvalid){
		
		byte[] line=bf.nextLine();
		
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;
				
				boolean valid=(line[0]!='#' && Tools.isDigit(line[line.length-1]));
				String[] split=null;
				if(valid){
					split=new String(line).split("\\|");
					valid=(split.length==4);
				}
//				assert(false) : (char)line[0]+", "+(char)line[line.length-1]+", "+new String(line).split("\\|").length;
				
				int code=-1;
				float latency=-1;
				if(valid){
					try {
						code=Integer.parseInt(split[2].substring(0, split[2].indexOf(' ')));
						latency=Float.parseFloat(split[3]);
					} catch (Exception e) {
						valid=false;
					}
				}
				
				if(valid){
					linesValid++;
					
					long[] cnt=map.get(split[2]);
					if(cnt==null){
						cnt=new long[1];
						map.put(split[2], cnt);
					}
					cnt[0]++;
					
					int latency2=(int)(latency*1000);
					
					if(code==200){
						passLatency.add(latency2);
					}else{
						failLatency.add(latency2);
						failCode.add(code);
						if(bswFail!=null){
							bswFail.println(line);
						}
					}
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}
		
		errorState|=bf.close();
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	/** List of input log file paths to process */
	private ArrayList<String> in1=new ArrayList<String>();
	/** Primary output file path for aggregated statistics */
	private String out1=null;
	/** Output file path for invalid/unparseable log entries */
	private String outInvalid=null;
	/** Output file path for failed requests (non-200 response codes) */
	private String outFail=null;
	
	/** Units suffix for latency display ("ms" for milliseconds or empty) */
	private String ms="ms";
	
	/*--------------------------------------------------------------*/
	
	/** Total number of log lines processed across all input files */
	private long linesProcessed=0;
	/** Number of successfully parsed and valid log entries */
	private long linesValid=0;
	/** Total bytes processed from all input log files */
	private long bytesProcessed=0;

	/** Collection of latency values for successful requests (200 response codes) */
	private IntList passLatency=new IntList();
	/** Collection of latency values for failed requests (non-200 response codes) */
	private IntList failLatency=new IntList();
	/** Collection of HTTP response codes for failed requests */
	private IntList failCode=new IntList();
	
	/** Maps response code strings to occurrence counts for statistics generation */
	private HashMap<String, long[]> map=new HashMap<String, long[]>();
	
	/** Maximum number of lines to process from input files (default unlimited) */
	private long maxLines=Long.MAX_VALUE;
	/**
	 * Whether to include extended statistics in output (latency metrics, fail codes)
	 */
	private boolean extendedStats=false;
	
	/*--------------------------------------------------------------*/
	
	/** File format objects for all input log files */
	private final ArrayList<FileFormat> ffin1;
	/** File format object for primary statistics output file */
	private final FileFormat ffout1;
	/** File format object for invalid entries output file */
	private final FileFormat ffoutInvalid;
	/** File format object for failed requests output file */
	private final FileFormat ffoutFail;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for console/redirected output (default System.out) */
	private PrintStream outstream=System.out;
	/** Global verbose logging flag for detailed processing information */
	public static boolean verbose=false;
	/** Tracks if any errors occurred during processing operations */
	public boolean errorState=false;
	/** Whether to overwrite existing output files (default true) */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
