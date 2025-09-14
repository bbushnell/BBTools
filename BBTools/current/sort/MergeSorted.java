package sort;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.SamLine;
import structures.ListNum;
import tax.AccessionToTaxid;
import tax.GiToTaxid;
import tax.TaxTree;
import tracker.ReadStats;
import var2.ScafMap;

/**
 * Sorts reads by name, potentially from multiple input files.
 * 
 * @author Brian Bushnell
 * @date September 21, 2016
 *
 */
public class MergeSorted {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		final boolean oldFI=FASTQ.FORCE_INTERLEAVED, oldTI=FASTQ.TEST_INTERLEAVED;
		MergeSorted x=new MergeSorted(args);
		x.process(t);
		FASTQ.FORCE_INTERLEAVED=oldFI;
		FASTQ.TEST_INTERLEAVED=oldTI;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public MergeSorted(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		boolean ascending=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbose2")){
				assert(false) : "Verbose2 is disabled.";
//				verbose2=Parse.parseBoolean(b);
			}else if(a.equals("allowsubprocess")){
				allowInputSubprocess=Parse.parseBoolean(b);
			}else if(a.equals("delete")){
				delete=Parse.parseBoolean(b);
			}else if(a.equals("ascending")){
				ascending=Parse.parseBoolean(b);
			}else if(a.equals("descending")){
				ascending=!Parse.parseBoolean(b);
			}else if(a.equals("length")){
				if(Parse.parseBoolean(b)){
					comparator=ReadLengthComparator.comparator;
				}
			}else if(a.equals("name")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorName.comparator;
				}
			}else if(a.equals("quality")){
				if(Parse.parseBoolean(b)){
					comparator=ReadQualityComparator.comparator;
				}
			}else if(a.equals("position")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorPosition.comparator;
				}
			}else if(a.equals("list") || a.equals("names")){
				comparator=new ReadComparatorList(b);
			}else if(a.equals("random") || a.equals("shuffle")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorRandom.comparator;
				}
			}else if(a.equals("taxa")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorTaxa.comparator;
				}
			}else if(a.equals("topo") || a.equals("topological") || a.equals("alphabetic") || a.equals("sequence") || a.equals("lexographic")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorTopological.comparator;
				}
			}else if(a.equals("flowcell")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorFlowcell.comparator;
				}
			}else if(a.equalsIgnoreCase("clump") || a.equalsIgnoreCase("clumpify")){
				if(Parse.parseBoolean(b)){
					comparator=ReadComparatorClump.comparator;
				}
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				if(b==null || "ignore".equalsIgnoreCase(b) || "skip".equalsIgnoreCase(b)){
					giTableFile=null;
					TaxTree.CRASH_IF_NO_GI_TABLE=false;
				}else{giTableFile=b;}
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("in") || a.equals("in1")){
				for(String s : b.split(",")){in1.add(s);}
			}else if(a.equals("maxfiles") || a.equals("files")){
				maxFiles=Integer.parseInt(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(b==null && new File(arg).exists()){
				in1.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		if(comparator==ReadComparatorTopological.comparator && genKmer){
			comparator=ReadComparatorTopological5Bit.comparator;
		}else{
			genKmer=false;
		}
		
		comparator.setAscending(ascending);
		SamLine.SET_FROM_OK=true;
		
		if(comparator==ReadComparatorRandom.comparator){
			ListNum.setDeterministicRandomSeed(-1);
			ListNum.setDeterministicRandom(true);
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;

			out1=parser.out1;
			out2=parser.out2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		if(comparator==ReadComparatorPosition.comparator){
			if(ReadComparatorPosition.scafMap==null){
				ReadComparatorPosition.scafMap=ScafMap.loadSamHeader(in1.get(0));
			}
		}
		
		if((comparator==ReadComparatorTaxa.comparator)){
			if(giTableFile!=null){
				outstream.println("Loading gi table.");
				GiToTaxid.initialize(giTableFile);
			}
			if(accessionFile!=null){
				outstream.println("Loading accession table.");
				AccessionToTaxid.load(accessionFile);
			}
			if(taxTreeFile!=null){
				ReadComparatorTaxa.tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, false);
				assert(ReadComparatorTaxa.tree.nameMap()!=null);
			}else{
				throw new RuntimeException("No tree specified.");
			}
		}
		
		{
			FileFormat ffin=null;
			if(in1!=null && in1.size()>0){
				ffin=FileFormat.testInput(in1.get(0), null, false);
			}
			tempExt=Tools.getTempExt(ffin, ffout1, extout);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		ArrayList<String> currentList=mergeRecursive(in1);
		merge(currentList, ffout1, ffout2);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Report timing and results
		t.stop("Time: \t");
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ArrayList<String> mergeRecursive(final ArrayList<String> inList){
		assert(maxFiles>1);
		ArrayList<String> currentList=inList;
		final int oldZL=ReadWrite.ZIPLEVEL;
		while(currentList.size()>maxFiles){
			ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, 4);
			final int size=currentList.size();
			final int groups=(size+maxFiles-1)/maxFiles;
			assert(groups>0 && groups<size);
			ArrayList<ArrayList<String>> listList=new ArrayList<ArrayList<String>>();
			ArrayList<String> outList=new ArrayList<String>();
			for(int i=0; i<groups; i++){
				listList.add(new ArrayList<String>());
			}
			for(int i=0; i<size; i++){
				listList.get(i%groups).add(currentList.get(i));
			}
			for(ArrayList<String> subList : listList){
				String temp=getTempFile();
				FileFormat ff=FileFormat.testOutput(temp, FileFormat.FASTQ, null, true, false, false, false);
				merge(subList, ff, null);
				outList.add(temp);
			}
			currentList=outList;
		}
		ReadWrite.ZIPLEVEL=oldZL;
		return currentList;
	}
	
	public String getTempFile(){
		String temp;
		File dir=new File(".");//(Shared.tmpdir()==null ? null : new File(Shared.tmpdir()));
		if(dir!=null && !dir.exists()){dir.mkdirs();}
		try {
			temp=File.createTempFile("sort_temp_", tempExt, dir).toString();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			KillSwitch.kill(e.getMessage());
			return null;
		}
		return temp;
	}
	
	public void merge(ArrayList<String> inList, FileFormat ff1, FileFormat ff2){
		errorState|=SortByName.mergeAndDump(inList, /*null, */ff1, ff2, delete, useSharedHeader, allowInputSubprocess, outstream, SortByName.maxSizeObservedStatic());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in1=new ArrayList<String>();

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	private String giTableFile=null;;
	private String taxTreeFile=null;
	private String accessionFile=null;
	
	/*--------------------------------------------------------------*/

//	long maxLengthObserved=0;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	
	private int maxFiles=16;
	
	private boolean delete=true;
	
	private boolean useSharedHeader=false;
	
	private boolean allowInputSubprocess=true;
	
	private String tempExt=null;
	
	private boolean genKmer=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	private static ReadComparator comparator=ReadComparatorName.comparator;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** Print verbose messages */
	public static final boolean verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
