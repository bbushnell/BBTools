package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import tax.TaxNode;
import tax.TaxTree;

/**
 * @author Brian Bushnell
 * @date Oct 17, 2014
 *
 */
public class SummarizeContamReport {
	
	/** Program entry point for contamination report summarization.
	 * @param args Command-line arguments including input files and parameters */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SummarizeContamReport x=new SummarizeContamReport(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs a SummarizeContamReport instance and parses command-line arguments.
	 * Initializes file formats, taxonomy tree, and processing parameters.
	 * @param args Command-line arguments containing input files, output options, and filters
	 */
	public SummarizeContamReport(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("minreads")){
				minReads=Long.parseLong(b);
			}else if(a.equals("minsequnits") || a.equals("minunits") || a.equals("minseqs")){
				minSeqUnits=Long.parseLong(b);
			}else if(a.equals("in")){
				for(String term : b.split(",")){
					in.add(term);
				}
			}else if(a.equals("tree")){
				treeFile=b;
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if("auto".equalsIgnoreCase(sizeFile)){sizeFile=TaxTree.defaultSizeFile();}
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=parser.overwrite;
			append=parser.append;

			out1=parser.out1;
		}
		
		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffinArray=new FileFormat[in.size()];
		for(int i=0; i<in.size(); i++){
			ffinArray[i]=FileFormat.testInput(in.get(i), FileFormat.TEXT, null, false, false);
		}
		
		tree=TaxTree.loadTaxTree(treeFile, System.err, true, false);
		if(tree!=null){tree.loadSizeFile(sizeFile);}
	}
	
	/**
	 * Main processing method that reads input files and generates summary output.
	 * Processes each input file, aggregates contamination data, and prints results.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){
		
		for(FileFormat ff : ffinArray){
			if(ff.canRead()){
				processOneFile(ff);
			}else{
				System.err.println("Skipping unreadable file "+ff.name());
			}
		}
		
		printOutput();
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Processes a single contamination report file.
	 * Parses the file line by line, extracting taxonomy and count information.
	 * @param ff FileFormat object representing the input file to process
	 */
	void processOneFile(FileFormat ff){

		
		final TextFile tf;
		{
			tf=new TextFile(ff);
			if(verbose){outstream.println("Started tf");}
		}
		
		{
			String line;
			
			line=tf.nextLine();
			if(line.startsWith("CONTAM SUMMARY")){
				line=tf.nextLine();
			}
			assert(line.startsWith("Examined ")) : line;
			line=tf.nextLine();
			assert(line.startsWith("|Taxonomy")) : line;
			
			while(/*(maxReads<0 || linesProcessed<maxReads) && */(line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				if(line.startsWith("|")){
					if(line.startsWith("|TOTAL")){break;}
					processLine(line);
				}else{
					break;
				}
			}
		}
		errorState|=tf.close();
	}
	
	/**
	 * Generates and writes the summarized contamination report output.
	 * Sorts results by sequence units, applies filters, and enriches with taxonomy data.
	 * Output includes taxonomy names, counts, tax IDs, and genome size information.
	 */
	private void printOutput(){
		final TextStreamWriter tsw;
		{
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
			tsw.println("#Name\tSeqUnits\tReads\tTaxID\tClade\tsize\tcSize\tseqs\tcSeqs\tcNodes");
		}
		
		ArrayList<StringLongLong> list=new ArrayList<StringLongLong>(map.size());
		list.addAll(map.values());
		Collections.sort(list, new ComparatorA());
		boolean filterA=minSeqUnits>1, filterB=minReads>1;
		boolean filter=filterA || filterB;
		for(StringLongLong sll : list){
			if(sll.a>=minSeqUnits && sll.b>=minReads){
//			if(!filter || (filterA && sll.a>=minSeqUnits) || (filterB && sll.b>=minReads)){
				int tid=-1;
				TaxNode tn;
				//			TaxNode tn=null;
				TaxNode ancestor=null;

				long size=0;
				long cumulative_size=0;
				long seqs=0;
				long cumulative_seqs=0;
				long cumulative_nodes=0;
				if(tree!=null){
					tid=tree.parseNameToTaxid(sll.s);
					if(tid>=0){
						tn=tree.getNode(tid);
						ancestor=tree.getNodeAtLevelExtended(tid, TaxTree.SUPERKINGDOM_E);
						size=tree.toSize(tn);
						cumulative_size=tree.toSizeC(tn);
						seqs=tree.toSeqs(tn);
						cumulative_seqs=tree.toSeqsC(tn);
						cumulative_nodes=tree.toNodes(tn);
					}
				}
				
				tsw.println(sll.s+"\t"+sll.a+"\t"+sll.b+"\t"+tid+"\t"+(ancestor==null ? "null" : ancestor.name)+
						"\t"+size+"\t"+cumulative_size+"\t"+seqs+"\t"+cumulative_seqs+"\t"+cumulative_nodes);
			}
		}
		errorState|=tsw.poisonAndWait();
	}
	
	/**
	 * Parses a single line from a contamination report file.
	 * Extracts taxonomy name and numeric values, aggregating counts in the global map.
	 * @param line Input line containing pipe-separated contamination data
	 */
	private void processLine(String line){
		String[] split=line.split("\\|");
		String[] split2=split[1].split(";");
		String name=split2[split2.length-1];
		try {
			long a=Long.parseLong(split[2]);
			long b=Long.parseLong(split[3]);
			StringLongLong p=map.get(name);
			if(p==null){
				p=new StringLongLong(name, a, b);
				map.put(name, p);
			}else{
				p.a+=a;
				p.b+=b;
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.err.println(line);
			System.err.println(Arrays.toString(split));
			System.err.println(Arrays.toString(split2));
			shared.KillSwitch.kill();
		}
	}
	
	
	/*--------------------------------------------------------------*/
	
	/** Comparator for sorting StringLongLong objects primarily by field 'a' (sequence units).
	 * Sorts in descending order by sequence units, then by reads, then alphabetically by name. */
	class ComparatorA implements Comparator<StringLongLong> {

		@Override
		public int compare(StringLongLong x, StringLongLong y) {
			if(x.a!=y.a){return  x.a<y.a ? 1 : -1;}
			if(x.b!=y.b){return  x.b<y.b ? 1 : -1;}
			return x.s.compareTo(y.s);
		}
		
	}
	
	/** Comparator for sorting StringLongLong objects primarily by field 'b' (read counts).
	 * Sorts in descending order by reads, then by sequence units, then alphabetically by name. */
	class ComparatorB implements Comparator<StringLongLong> {

		@Override
		public int compare(StringLongLong x, StringLongLong y) {
			if(x.b!=y.b){return  x.b<y.b ? 1 : -1;}
			if(x.a!=y.a){return  x.a<y.a ? 1 : -1;}
			return x.s.compareTo(y.s);
		}
		
	}
	
	/** Container class holding a string name and two long values.
	 * Used to store taxonomy names with associated sequence unit and read counts. */
	class StringLongLong {
		
		/** Constructs a StringLongLong with only a string value.
		 * @param s_ The string value to store */
		StringLongLong(String s_){
			s=s_;
		}
		
		/**
		 * Constructs a StringLongLong with string and two long values.
		 * @param s_ The string value (typically taxonomy name)
		 * @param a_ First long value (typically sequence units)
		 * @param b_ Second long value (typically read count)
		 */
		StringLongLong(String s_, long a_, long b_){
			s=s_;
			a=a_;
			b=b_;
		}
		
		/** String value, typically containing taxonomy name */
		final String s;
		/** First numeric value, typically sequence units count */
		long a;
		/** Second numeric value, typically read count */
		long b;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** List of input file paths to process */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file path for the summarized report */
	private String out1=null;
	/** Path to taxonomy tree file, defaults to "auto" */
	private String treeFile="auto";
	/** Path to genome size file, defaults to "auto" */
	private String sizeFile="auto";
	
	/** Taxonomy tree for resolving taxonomic information */
	TaxTree tree=null;
	/** Map storing aggregated counts by taxonomy name */
	private HashMap<String, StringLongLong> map=new HashMap<String, StringLongLong>();
	
	/*--------------------------------------------------------------*/

	/** Minimum read count threshold for output filtering */
	long minReads=0;
	/** Minimum sequence units threshold for output filtering */
	long minSeqUnits=0;
	
	/** Counter for total lines processed from input files */
	long linesProcessed=0;
	/** Counter for total characters processed from input files */
	long charsProcessed=0;
	
	/*--------------------------------------------------------------*/
	
	/** Array of FileFormat objects for input files */
	private final FileFormat ffinArray[];
	/** FileFormat object for the primary output file */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Print stream for status messages and logging output */
	private PrintStream outstream=System.err;
	/** Flag controlling verbosity of status output */
	public static boolean verbose=false;
	/** Flag indicating if an error occurred during processing */
	public boolean errorState=false;
	/** Flag controlling whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Flag controlling whether to append to existing output files */
	private boolean append=false;
	
}
