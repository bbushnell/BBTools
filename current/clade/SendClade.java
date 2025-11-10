package clade;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import bin.AdjustEntropy;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import server.ServerTools;
import shared.KillSwitch;
import shared.LineParserS1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.SendSketch;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringNum;
import tracker.ReadStats;

/**
 * Client-side tool for sending clade queries to a remote taxonomic classification server.
 * Processes FASTA sequences by generating clade objects containing genomic statistics
 * and signature data, then transmits them to a clade server for taxonomic identification.
 * Supports both single-sequence aggregation and per-contig analysis modes.
 * Mirrors the SendSketch architecture for consistent client-server communication patterns.
 * Results are returned in either human-readable or machine-parseable oneline format.
 *
 * @author Chloe
 * @date September 14, 2025
 */
public class SendClade extends CladeObject {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Entry point for command line execution. Creates a SendClade instance,
	 * processes all input files, and ensures clean shutdown with error checking.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		assert(args != null) : "Command line arguments cannot be null";

		Timer t=new Timer();
		assert(t != null) : "Timer initialization failed";

		//Create an instance of this class
		SendClade x=new SendClade(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);

		assert(!x.errorState) : "This program ended in an error state.";
	}

	/**
	 * Constructor that parses command line arguments and initializes all configuration.
	 * Sets up input/output handling, validates file accessibility, configures server
	 * communication parameters, and initializes entropy models for clade calculations.
	 * @param args Command line arguments
	 */
	public SendClade(String[] args){
		assert(args != null) : "Constructor arguments cannot be null";

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		assert(outstream != null) : "Output stream initialization failed";

		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		//Initialize entropy model for clade calculations
		if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150){
			AdjustEntropy.load(4, 150);
		}
		assert(AdjustEntropy.kLoaded == 4 && AdjustEntropy.wLoaded == 150) :
			"Entropy model initialization failed: k=" + AdjustEntropy.kLoaded + " w=" + AdjustEntropy.wLoaded;

		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";

		//Set defaults
		banSelf=false;
		printQTID=false;
		heapSize=1;
		perContig=false;
		oneline=false;
		hits=1;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			assert(arg != null) : "Null argument at position " + i;

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("local")){
				local=Parse.parseBoolean(b);
				if(local){address=localAddress;}
			}else if(a.equals("address") || a.equals("server")){
				address=b;
			}else if(a.equals("percontig") || a.equals("persequence")){
				perContig=Parse.parseBoolean(b);
			}else if(a.equals("mode")){
				if("perseq".equals(b) || "persequence".equals(b) || "percontig".equals(b) || "sequence".equals(b)){
					perContig=true;
				}else{
					perContig=false;
				}
			}else if(a.equals("oneline") || a.equals("machine")){
				oneline=Parse.parseBoolean(b);
			}else if(a.equals("format")){
				oneline="oneline".equals(b) || "machine".equals(b);
			}else if(a.equals("hits")){
				hits=Integer.parseInt(b);
				assert(hits > 0) : "Hits must be positive: " + hits;
			}else if(a.equals("printqtid") || a.equals("qtid")){
				printQTID=Parse.parseBoolean(b);
			}else if(a.equals("banself")){
				banSelf=Parse.parseBoolean(b);
			}else if(a.equals("heap")){
				heapSize=Integer.parseInt(b);
				assert(heapSize > 0) : "Heap size must be positive: " + heapSize;
			}else if(a.equals("minlen") || a.equals("mincontig")){
				minlen=Integer.parseInt(b);
				assert(minlen >= 0) : "Minimum length cannot be negative: " + minlen;
			}else if(a.equals("concise")){
				Clade.CONCISE=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("a48")){
				Clade.outputCoding=Parse.parseBoolean(b) ? Clade.A48 : Clade.DECIMAL;
			}else if(a.equals("in")){
				Tools.getFileOrFiles(b, in, true, false, false, false);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(new File(arg).canRead()){
				Tools.getFileOrFiles(arg, in, true, false, false, false);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			if(parser.in1!=null){in.add(parser.in1);}
			out1=parser.out1;
		}

		//Ensure there is an input file
		if(in==null || in.isEmpty()){
			throw new RuntimeException("Error - at least one input file is required.");
		}
		assert(!in.isEmpty()) : "Input file list should not be empty after validation";
		assert(in.size() > 0) : "Input file list size should be positive";

		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		ArrayList<String> allFiles=new ArrayList<String>(in);
		if(out1!=null){allFiles.add(out1);}
		if(!Tools.testForDuplicateFiles(true, allFiles.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		//Create output file - use stdout as default
		ffout=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);

		//Ensure address is valid
		assert(address!=null) : "No server address specified. Use address=<address> or local=t";
		assert(address.length() > 0) : "Server address cannot be empty";
		if(!address.startsWith("http://") && !address.startsWith("https://")){address="http://"+address;}
		assert(address.startsWith("http://") || address.startsWith("https://")) : "Address should start with http:// or https:// after normalization";
		assert(address.length() > 7) : "Address too short after normalization: " + address;

		if(verbose){
			System.err.println("[" + new java.util.Date() + "] SendClade configured - server: " + address + ", hits: " + hits + ", format: " + (oneline ? "oneline" : "human"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Processing           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Primary processing method that coordinates the entire workflow.
	 * Sets up unified output handling, processes all input files sequentially,
	 * and reports comprehensive timing and statistics upon completion.
	 * @param t Timer for tracking total processing time
	 */
	private void process(Timer t){
		assert(t != null) : "Timer cannot be null";
		if(verbose){System.err.println("[" + new java.util.Date() + "] SendClade.process() starting");}
		Clade.MAKE_FREQUENCIES=false;

		//Set up unified output stream - always use TextStreamWriter for consistency
		final TextStreamWriter tsw;
		if(ffout!=null) {
			tsw=new TextStreamWriter(ffout);
			tsw.start();
		}else {tsw=null;}
		
		if(oneline && tsw!=null){
			//Print header for oneline format
			tsw.println("#QueryName\tQ_GC\tQ_Bases\tQ_Contigs\tRefName\tR_TaxID\tR_GC\tR_Bases\tR_Contigs\tR_Level\tGCdif\tSTRdif\tk3dif\tk4dif\tk5dif\tlineage");
		}

		//Process all input files
		long sequencesProcessed=0;
		assert(in.size() > 0) : "No input files to process";
		if(verbose){System.err.println("[" + new java.util.Date() + "] Processing " + in.size() + " input files");}

		for(String fname : in){
			assert(fname != null) : "Null filename in input list";
			if(verbose){System.err.println("[" + new java.util.Date() + "] Processing file: " + fname);}
			FileFormat ff=FileFormat.testInput(fname, FileFormat.UNKNOWN, null, true, true);
			assert(ff != null) : "FileFormat creation failed for " + fname;
			long seqs = process_inner(ff, tsw);
			assert(seqs >= 0) : "Invalid sequence count: " + seqs;
			if(verbose){System.err.println("[" + new java.util.Date() + "] Processed " + seqs + " sequences from " + fname);}
			sequencesProcessed+=seqs;
		}

		//Clean up
		if(tsw!=null) {tsw.poisonAndWait();}
		assert(sequencesProcessed >= 0) : "Invalid total sequence count: " + sequencesProcessed;

		//Report timing and results
		if(verbose){
			t.stop();
			outstream.println("\nTime: \t"+t);
			outstream.println("Sequences Processed: "+sequencesProcessed);
		}

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Processes a single input file by detecting format and routing appropriately.
	 * Routes .clade/.spectra files to processCladeFile() and FASTA/FASTQ to processFastaFile().
	 * @param ff Input file format descriptor
	 * @param tsw Unified output writer for results
	 * @return Number of sequences successfully processed
	 */
	private long process_inner(FileFormat ff, TextStreamWriter tsw){
		assert(ff != null) : "FileFormat cannot be null";
		if(verbose){System.err.println("[" + new java.util.Date() + "] process_inner() starting for " + ff.name());}

		// Check file format and route to appropriate processing method
		if (ff.clade()) {
			return processCladeFile(ff, tsw);
		} else {
			return processFastaFile(ff, tsw);
		}
	}

	/**
	 * Processes .clade/.spectra files by loading pre-computed clades directly.
	 * Uses CladeLoader to read clade objects and sends them in batches to the server.
	 * @param ff Input file format descriptor for .clade/.spectra file
	 * @param tsw Unified output writer for results
	 * @return Number of clades successfully processed
	 */
	private long processCladeFile(FileFormat ff, TextStreamWriter tsw){
		assert(ff != null) : "FileFormat cannot be null";
		assert(ff.clade()) : "File format must be .clade or .spectra: " + ff.name();
		if(verbose){System.err.println("[" + new java.util.Date() + "] processCladeFile() starting for " + ff.name());}

		//Load clades directly from clade file
		ArrayList<Clade> allClades = CladeLoader.loadCladesFromClade(ff);
		assert(allClades != null) : "CladeLoader returned null for " + ff.name();
		if(verbose){System.err.println("[" + new java.util.Date() + "] Loaded " + allClades.size() + " clades from " + ff.name());}

		//Send clades in batches
		long cladesProcessed = 0;
		ArrayList<Clade> batch = new ArrayList<Clade>();

		for(Clade clade : allClades){
			batch.add(clade);
			cladesProcessed++;

			//Send batch if buffer is full
			if(batch.size() >= MAX_CLADES_PER_BATCH){
				if(verbose){System.err.println("[" + new java.util.Date() + "] Batch full with " + batch.size() + " clades, sending to server");}
				sendAndPrint(batch, tsw);
				batch.clear();
			}
		}

		//Send any remaining clades
		if(!batch.isEmpty()){
			if(verbose){System.err.println("[" + new java.util.Date() + "] Sending final batch of " + batch.size() + " clades");}
			sendAndPrint(batch, tsw);
		}

		if(verbose){outstream.println("Processed " + cladesProcessed + " clades from " + ff.name());}
		return cladesProcessed;
	}

	/**
	 * Processes FASTA/FASTQ files by reading sequences and organizing them into clades.
	 * Supports both per-contig mode (each sequence becomes separate clade) and aggregation
	 * mode (all sequences combined into single clade). Handles batching to prevent memory
	 * overflow and maintains detailed progress tracking.
	 * @param ff Input file format descriptor for FASTA/FASTQ file
	 * @param tsw Unified output writer for results
	 * @return Number of sequences successfully processed
	 */
	private long processFastaFile(FileFormat ff, TextStreamWriter tsw){
		assert(ff != null) : "FileFormat cannot be null";
		if(verbose){System.err.println("[" + new java.util.Date() + "] processFastaFile() starting for " + ff.name());}

		//Load sequences using standard BBTools ConcurrentReadInputStream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
			cris.start();
		}
		assert(cris != null) : "ConcurrentReadInputStream creation failed for " + ff.name();
		ArrayList<Clade> clades=new ArrayList<Clade>();

		//Track progress
		long sequencesProcessed=0;

		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r != null) : "Null read encountered";
					if(perContig){
						//Process each contig separately
						if(r.bases==null || r.bases.length<minlen){continue;}
						assert(r.bases.length >= minlen) : "Sequence shorter than minimum length: " + r.bases.length + " < " + minlen;

						Clade clade=new Clade(0, 0, r.id.replace('\t', ' '));
						assert(clade != null) : "Clade creation failed for " + r.id;
						clade.add(r.bases, null);
						clade.finish();
						clades.add(clade);
						sequencesProcessed++;
						assert(sequencesProcessed > 0) : "Sequence count should be positive";

						//Send batch if buffer is full
						if(clades.size()>=MAX_CLADES_PER_BATCH){
							assert(clades.size() <= MAX_CLADES_PER_BATCH) : "Clade batch size exceeded maximum: " + clades.size() + " > " + MAX_CLADES_PER_BATCH;
							if(verbose){System.err.println("[" + new java.util.Date() + "] Batch full with " + clades.size() + " clades, sending to server");}
							sendAndPrint(clades, tsw);
							clades.clear();
							assert(clades.isEmpty()) : "Clade list should be empty after clearing";
						}
					}else{
						//Accumulate all sequences into one clade
						if(clades.isEmpty()){
							clades.add(new Clade(0, 0, ff.simpleName()));
							assert(!clades.isEmpty()) : "Clade list should not be empty after adding";
						}
						assert(clades.size() == 1) : "Single-sequence mode should have exactly one clade: " + clades.size();
						Clade clade=clades.get(0);
						assert(clade != null) : "Null clade in single-sequence mode";
						clade.add(r.bases, null);
						sequencesProcessed++;
					}
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		//Finish any accumulated clade
		if(!perContig && !clades.isEmpty()){
			Clade clade=clades.get(0);
			clade.finish();
		}

		//Send any remaining clades
		if(!clades.isEmpty()){
			if(verbose){System.err.println("[" + new java.util.Date() + "] Sending final batch of " + clades.size() + " clades");}
			sendAndPrint(clades, tsw);
		}

		ReadWrite.closeStream(cris);

		if(verbose){outstream.println("Processed "+sequencesProcessed+" sequences from "+ff.name());}

		return sequencesProcessed;
	}
	
	/**
	 * Transmits a batch of clades to the remote server and processes the response.
	 * Constructs the request message with all configuration parameters, handles server
	 * communication with comprehensive error checking, and routes the response through
	 * the unified output system. Includes detailed timing and debugging information.
	 * @param clades List of clades to transmit
	 * @param tsw Unified output writer for server response
	 */
	private void sendAndPrint(ArrayList<Clade> clades, TextStreamWriter tsw) {
		String s=sendClades(clades);
		//Write response to unified output
		if(tsw!=null) {tsw.print(s);}
	}
	
	public static boolean sendAndLabel(List<Clade> clades) {
		String response=sendClades(clades, null, true, 1, false, false, 1, false);
		ArrayList<Comparison> comps=null;
		try{
			comps=SendClade.responseToComparisons(response);
			assert(comps.size()==clades.size());
			for(int i=0; i<clades.size(); i++) {
				Clade clade=clades.get(i);
				Comparison comp=comps.get(i);
				if(clade!=null && comp!=null) {
					clade.name=comp.ref.name;
					clade.taxID=comp.ref.taxID;
					clade.lineage=comp.ref.lineage;
				}
			}
		}catch(Throwable e){
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	public static ArrayList<Comparison> responseToComparisons(String s) {

		
//		#Query1
//		f_0_c_0_s_0_p_6294_i_40_tid_1002367 1:	0.475	40	1	Cortinarius geophilus var. subauroreus	2764306	0.416	609	1	1	0.059	0.167	0.579	1.000	1.000	sk__Eukaryota;k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Cortinariaceae;g__Cortinarius;s__Cortinarius geophilus
//		#Query2
//		f_0_c_0_s_0_p_11130_i_40_tid_1002367 1:	0.475	40	1	Passion fruit yellow mosaic virus	185692	0.540	1115	1	1	0.065	0.013	0.423	1.000	1.000	sk__Viruses;k__Orthornavirae;p__Kitrinoviricota;c__Alsuviricetes;o__Tymovirales;f__Tymoviridae;g__Tymovirus;s__Tymovirus passiflorae;ss__Passion fruit yellow mosaic virus
//		#Query3
		String[] lines=s.split("\n");
		ArrayList<Comparison> list=new ArrayList<Comparison>(lines.length/2+1);
		LineParserS1 lp=new LineParserS1('\t');
		for(int i=0; i<lines.length; i++) {
			String a=lines[i], b=null;
			Comparison c=null;
			assert(a.startsWith("#Query")) : "\n"+a+"\n"+b+"\n";
			if(i+1<lines.length) {b=lines[i+1];}
			if(b!=null && !b.startsWith("#Query")) {
				c=new Comparison(b, lp);
				i++;
			}
			list.add(c);
		}
		return list;
	}
	
	private String sendClades(Collection<Clade> clades) {
		return sendClades(clades, null, oneline, hits, 
			printQTID, banSelf, heapSize, verbose);
	}

	/**
	 * Transmits a batch of clades to the remote server and processes the response.
	 * Constructs the request message with all configuration parameters, handles server
	 * communication with comprehensive error checking, and routes the response through
	 * the unified output system. Includes detailed timing and debugging information.
	 * @param clades List of clades to transmit
	 * @param tsw Unified output writer for server response
	 */
	public static String sendClades(Collection<Clade> clades, String address, boolean oneline, int hits, 
			boolean printQTID, boolean banSelf, int heapSize, boolean verbose){
			if(clades.size()<=MAX_CLADES_PER_BATCH) {
				return sendBatch(clades, address, oneline, hits, 
					printQTID, banSelf, heapSize, verbose);
			}
			StringBuilder sb=new StringBuilder();
			ArrayList<Clade> temp=new ArrayList<Clade>(MAX_CLADES_PER_BATCH);
			for(Clade c : clades) {
				temp.add(c);
				if(temp.size()>=MAX_CLADES_PER_BATCH) {
					String s=sendBatch(temp, address, oneline, hits, 
						printQTID, banSelf, heapSize, verbose);
					sb.append(s);
					temp.clear();
				}
			}
			if(temp.size()>0) {
				String s=sendBatch(temp, address, oneline, hits, 
					printQTID, banSelf, heapSize, verbose);
				sb.append(s);
				temp.clear();
			}
			return sb.toString();
	}
	
	private String sendBatch(Collection<Clade> clades) {
		return sendBatch(clades, null, oneline, hits, 
			printQTID, banSelf, heapSize, verbose);
	}

	/**
	 * Transmits a batch of clades to the remote server and processes the response.
	 * Constructs the request message with all configuration parameters, handles server
	 * communication with comprehensive error checking, and routes the response through
	 * the unified output system. Includes detailed timing and debugging information.
	 * @param clades List of clades to transmit
	 * @param tsw Unified output writer for server response
	 */
	private static String sendBatch(Collection<Clade> clades, String address, boolean oneline, int hits, 
			boolean printQTID, boolean banSelf, int heapSize, boolean verbose){
		if(clades==null || clades.isEmpty()){return null;}
		assert(clades.size()<=MAX_CLADES_PER_BATCH) : clades.size()+">"+MAX_CLADES_PER_BATCH;
		if(verbose){System.err.println("[" + new java.util.Date() + "] sendClades() called with " + clades.size() + " clades");}
		if(address==null) {address=defaultAddress;}
		
		//Send to server
		Timer t=new Timer();
		byte[] message=toMessage(clades, oneline, hits, banSelf, printQTID, heapSize);
		if(verbose){
			t.stopAndStart("toMessage bytes: "+message.length+", time:");
			System.err.println("[" + new java.util.Date() + "] Sending " + clades.size() + " clades (" + message.length + " bytes) to " + address);
			if(message.length < 500) {
				System.err.println("[" + new java.util.Date() + "] Message content: " + new String(message));
			}
			System.err.println("Sending "+clades.size()+" clades ("+message.length+" bytes) to "+address);
		}
		String response=sendMessage(message, address, verbose);
		if(verbose) {t.stopAndStart("sendClades time:");}
		return response;
	}

	/**
	 * Transmits a batch of clades to the remote server and processes the response.
	 * Constructs the request message with all configuration parameters, handles server
	 * communication with comprehensive error checking, and routes the response through
	 * the unified output system. Includes detailed timing and debugging information.
	 * @param clades List of clades to transmit
	 * @param tsw Unified output writer for server response
	 */
	public static byte[] toMessage(Collection<Clade> clades, boolean oneline, int hits, 
			boolean printQTID, boolean banSelf, int heapSize){
		if(clades==null || clades.isEmpty()){return null;}
		//Build message
		ByteBuilder bb=new ByteBuilder();

		//Add parameters
		bb.append("format=").append(oneline ? "oneline" : "human").append('/');
		bb.append("hits=").append(hits).append('/');
		if(printQTID){bb.append("printqtid=t/");}
		if(banSelf){bb.append("banself=t/");}
		bb.append("heap=").append(heapSize).append('/');
		bb.append('\n');

		//Add clades
		for(Clade clade : clades){clade.toBytes(bb);}

		//Send to server
		byte[] message=bb.toBytes();
		assert(message != null) : "Message creation failed";
		assert(message.length > 0) : "Empty message created";
		assert(message.length < 100000000) : "Message too large: " + message.length + " bytes";
		return message;
	}

	/**
	 * Transmits a batch of clades to the remote server and processes the response.
	 * Constructs the request message with all configuration parameters, handles server
	 * communication with comprehensive error checking, and routes the response through
	 * the unified output system. Includes detailed timing and debugging information.
	 * @param clades List of clades to transmit
	 * @param tsw Unified output writer for server response
	 */
	public static String sendMessage(byte[] message, String address, boolean verbose){
		if(message==null || message.length==0){return null;}
		if(verbose){System.err.println("[" + new java.util.Date() + "] sendClades() called with " + message.length + " bytes");}
		if(address==null) {address=defaultAddress;}

		try{
			if(verbose){System.err.println("[" + new java.util.Date() + "] Calling ServerTools.sendAndReceive()");}
			Timer sendTimer=new Timer();
			StringNum result=sendAndReceive(message, address);
			assert(result != null) : "Server returned null result";
			assert(result.n >= 100 && result.n < 600) : "Invalid HTTP status code: " + result.n;
			sendTimer.stop();
			long sendTime=sendTimer.elapsed;
			assert(sendTime >= 0) : "Invalid timing measurement: " + sendTime;
			if(verbose){
				System.err.println("[" + new java.util.Date() + "] Server responded in " + (sendTime/1000000) + "ms with code " + result.n);
				System.err.println("[" + new java.util.Date() + "] Response length: " + (result.s != null ? result.s.length() : 0) + " chars");
				if(result.s != null && result.s.length() < 500) {
					System.err.println("[" + new java.util.Date() + "] Response: " + result.s);
				}
			}
			if(!ServerTools.suppressErrors && (result.n<200 || result.n>299)){
				System.err.println("ERROR: Server returned code "+result.n+" and this message:\n"+result.s);
				KillSwitch.kill();
			}
			return result.s;
		}catch(Exception e){
			if(verbose){
				System.err.println("[" + new java.util.Date() + "] ERROR in sendClades: " + e.getMessage());
				System.err.println("[" + new java.util.Date() + "] Stack trace:");
			}
			e.printStackTrace();
		}
		return null;
	}
	
	private static StringNum sendAndReceive(byte[] message, String address) {
		StringNum sn=null;
		if(sync) {
			synchronized(SendSketch.class) {
				sn=ServerTools.sendAndReceive(message, address);
			}
		}else {
			while(concurrency.addAndGet(1)>maxConcurrency) {
				concurrency.addAndGet(-1);
				try{Thread.sleep(20);}
				catch(InterruptedException e){}
			}
			sn=ServerTools.sendAndReceive(message, address);
			concurrency.addAndGet(-1);
		}
		return sn;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input files */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file name */
	private String out1=null;
	/** Output file format */
	private FileFormat ffout=null;

	/** Server address */
	private String address=defaultAddress;
	/** Use local server */
	private boolean local=false;

	/** Process each contig separately */
	private boolean perContig=false;
	/** Print in one-line format */
	private boolean oneline=false;
	/** Number of hits to return */
	private int hits=1;
	/** Heap size for comparisons */
	private int heapSize=1;
	/** Minimum contig length */
	private int minlen=0;
	/** Maximum reads to process */
	private long maxReads=-1;

	/** Print query TaxID */
	private boolean printQTID=false;
	/** Ban self-matches */
	private boolean banSelf=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output stream */
	private final PrintStream outstream;

	/** True if an error was encountered */
	private boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;

	/** Display verbose output */
	private boolean verbose=false;
	
	private static AtomicInteger concurrency=new AtomicInteger(0);
	public static boolean sync=false;
	public static int maxConcurrency=8;

	/** Default server address */
	static final String defaultAddress="https://bbmapservers.jgi.doe.gov/quickclade";
	/** Local server address */
	private static final String localAddress="http://localhost:5002";
	/** Maximum clades to send in one batch */
	private static final int MAX_CLADES_PER_BATCH=4000;

}