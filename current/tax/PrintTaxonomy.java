package tax;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Filters sequences and prints taxonomy information based on sequence identifiers.
 * Accepts FASTA/FASTQ/SAM/text inputs labeled with gi numbers or NCBI taxIDs and can emit either entire hierarchies or a specific taxonomic level.
 * Supports column-based translation, direct name lists, or streamed reads while tracking per-node counts.
 * @author Brian Bushnell
 * @date November 23, 2015
 */
public class PrintTaxonomy {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Program entry point that constructs PrintTaxonomy, runs process(Timer), and closes shared streams.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		PrintTaxonomy x=new PrintTaxonomy(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Parses command-line arguments, loads taxonomy tables/trees, configures Parser overrides, and validates file formats.
	 * Initializes optional gi/accession tables and builds the TaxTree used for lookups.
	 * @param args Command-line arguments for configuration
	 */
	public PrintTaxonomy(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		//Create a parser object
		Parser parser=new Parser();
		
		int taxLevel=0, minLevel=0, maxLevel=TaxTree.LIFE;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("out")){
				out1=b;
			}else if(a.equals("counts")){
				countFile=b;
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("level") || a.equals("lv") || a.equals("taxlevel") || a.equals("tl")){
				taxLevel=TaxTree.parseLevel(b);
			}else if(a.equals("minlevel")){
				minLevel=TaxTree.parseLevel(b);
			}else if(a.equals("maxlevel")){
				maxLevel=TaxTree.parseLevel(b);
			}else if(a.equals("printname")){
				printName=Parse.parseBoolean(b);
			}else if(a.equals("reverse")){
				reverseOrder=Parse.parseBoolean(b);
			}else if(a.equals("silva")){
				TaxTree.SILVA_MODE=Parse.parseBoolean(b);
			}else if(a.equals("unite")){
				TaxTree.UNITE_MODE=Parse.parseBoolean(b);
			}else if(a.equals("simple")){
				skipNonCanonical=Parse.parseBoolean(b);
			}else if(a.equals("column")){
				keyColumn=Integer.parseInt(b);
			}else if(b!=null && (a.equals("name") || a.equals("names") || a.equals("id") || a.equals("ids"))){
				for(String s : b.split(",")){
					names.add(s);
				}
			}else{
				names.add(arg);
			}
		}
		
		if(taxTreeFile==null || "auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		taxLevelExtended=TaxTree.levelToExtended(taxLevel);
		minLevelExtended=TaxTree.levelToExtended(minLevel);
		maxLevelExtended=TaxTree.levelToExtended(maxLevel);
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			maxReads=parser.maxReads;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);
		
		ffcount=FileFormat.testOutput(countFile, FileFormat.TEXT, null, true, overwrite, append, false);
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, false);
		
		if(giTableFile!=null){
			outstream.println("Loading gi table.");
			GiToTaxid.initialize(giTableFile);
		}
		if(accessionFile!=null){
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
		}
		if(taxTreeFile!=null){
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, true, true);
			assert(tree.nameMap!=null);
		}else{
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates writers, dispatches to the appropriate processing path (reads, text input, or command-line names), and writes optional count summaries.
	 * Closes writers, propagates error flags, and throws if any worker reported failure.
	 * @param t Timer tracking the total execution time
	 */
	void process(Timer t){
		
		TextStreamWriter tsw=null;
		if(ffout1!=null){
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
		}
		
		if(ffin1!=null){
			if(ffin1.fasta() || ffin1.fastq() || ffin1.samOrBam() || ffin1.scarf()){
				processReads(tsw);
			}else{
				processFile(new TextFile(ffin1), tsw);
			}
		}else{
			processNames(tsw);
		}
		
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		if(ffcount!=null){
			TextStreamWriter tswc=new TextStreamWriter(ffcount);
			tswc.start();
			for(TaxNode tn : tree.nodes){
				if(tn!=null && tn.countRaw>0){
					tswc.println(tn.countRaw+"\t"+tn.name);
				}
			}
			errorState|=tswc.poisonAndWait();
		}
		
		t.stop();
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterates over the collected name/id list from the CLI, printing either the full taxonomy or a single level per entry.
	 * @param tsw Output writer for taxonomy results (may be null when writing to stdout) */
	void processNames(final TextStreamWriter tsw){
		for(String name : names){
			if(taxLevelExtended>0){
				printTaxLevel(name, tsw);
			}else{
				printTaxonomy(name, tsw);
			}
		}
	}
	
	/**
	 * Streams a tab-delimited text file, translating the configured keyColumn or printing taxonomy/level information per line.
	 * @param tf Input text file containing names or identifiers
	 * @param tsw Output writer for taxonomy results
	 */
	void processFile(final TextFile tf, final TextStreamWriter tsw){
		for(String name=tf.nextLine(); name!=null; name=tf.nextLine()){
			
			if(keyColumn>=0){
				String result=translateLine(name, keyColumn);
				tsw.print(result);
			}else if(taxLevelExtended>0){
				printTaxLevel(name, tsw);
			}else{
				printTaxonomy(name, tsw);
			}
		}
	}
	
	/** Consumes sequence files through ConcurrentReadInputStream, resolves taxonomy from read headers, and prints hierarchy/level data per read.
	 * @param tsw Output writer for taxonomy results */
	void processReads(final TextStreamWriter tsw){
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){System.err.println("Started cris");}
			cris.start();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				if(keyColumn>=0){
					String result=translateLine(r1.id, keyColumn);
					tsw.println(result);
				}else if(taxLevelExtended>0){
					printTaxLevel(r1.id, tsw);
				}else{
					printTaxonomy(r1.id, tsw);
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Translates a tab-delimited line by replacing the specified column with formatted taxonomy text or appending NOT_FOUND markers.
	 * @param line Input line containing tab-delimited tokens
	 * @param col Column index that holds the sequence identifier
	 * @return Line rebuilt with taxonomy strings in the requested column
	 */
	String translateLine(String line, int col){
		StringBuilder sb=new StringBuilder();
		String[] split=line.split("\t");
		assert(split.length>col) : "Too few columns in line:\n"+line+"\n->\n"+Arrays.toString(split);
		
		if(col<split.length){
			String name=split[col];
			while(name.startsWith(">") || name.startsWith("@")){name=name.substring(1);}
			
			TaxNode tn=parseNodeFromHeader(name);
			if(tn!=null){
				String tl=makeTaxLine(tree, tn, minLevelExtended, maxLevelExtended, skipNonCanonical, reverseOrder).toString();
				split[col]=tl;
			}else{
				List<TaxNode> list=tree.getNodesByNameExtended(name);
				if(list!=null){
					String tab="";
					for(TaxNode tn2 : list){
						sb.append(tab);
						sb.append(makeTaxLine(tree, tn2, minLevelExtended, maxLevelExtended, skipNonCanonical, reverseOrder).toString());
						tab="\t";
					}
				}else{
					split[col]=split[col]+"_***NOT_FOUND***";
				}
			}
		}
		
		for(int i=0; i<split.length; i++){
			if(i>0){sb.append('\t');}
			sb.append(split[i]);
		}
		sb.append('\n');
		return sb.toString();
	}
	
	/**
	 * Resolves a sequence name to one or more TaxNodes and prints the complete hierarchy for each hit, including optional name headers.
	 * @param name Sequence name or identifier stripped of FASTA/FASTQ markers
	 * @param tsw Output writer that receives the formatted taxonomy
	 */
	void printTaxonomy(String name, final TextStreamWriter tsw){
		while(name.startsWith(">") || name.startsWith("@")){name=name.substring(1);}
		tsw.print("\n");
		if(printName){tsw.print(name+":\n");}
		TaxNode tn=parseNodeFromHeader(name);
		if(tn!=null){
			printTaxonomy(tn, tsw);
			return;
		}else{
			List<TaxNode> list=tree.getNodesByNameExtended(name);
			if(list!=null){
				String nl="";
				for(TaxNode tn2 : list){
					tsw.print(nl);
					printTaxonomy(tn2, tsw);
					nl="\n";
				}
				return;
			}
		}
		tsw.println("Could not find node" + (printName ? "." : " for '"+name+"'"));
		return;
	}
	
	/**
	 * Emits only the target taxonomy level for a sequence name, expanding multiple matches when necessary.
	 * @param name Sequence name or identifier stripped of FASTA/FASTQ markers
	 * @param tsw Output writer for taxonomy results
	 */
	void printTaxLevel(String name, final TextStreamWriter tsw){
		while(name.startsWith(">") || name.startsWith("@")){name=name.substring(1);}
		tsw.print("\n");
		if(printName){tsw.print(name+":\n");}
		TaxNode tn=parseNodeFromHeader(name);
		if(tn!=null){
			printTaxLevel(tn, tsw);
			return;
		}else{
			List<TaxNode> list=tree.getNodesByNameExtended(name);
			if(list!=null){
				for(TaxNode tn2 : list){
					printTaxLevel(tn2, tsw);
				}
				return;
			}
		}
		tsw.println("Could not find node" + (printName ? "." : " for '"+name+"'"));
		return;
	}
	
//	void printTaxCounts(String name, final TextStreamWriter tsw){
//		TaxNode tn=null;
//		tn=tree.getNode(name);
//		if(tn==null){tn=tree.getNodeByName(name);}
//		if(tn==null){tn=unknown;}
//		while(tn!=null && tn.id!=tn.pid && tn.level<taxLevel){tn=tree.getNode(tn.pid);}
//		if(tsw!=null)tsw.println(tn.name);
//		tn.incrementRaw(1);
//	}
	
	/**
	 * Traverses upward from a resolved TaxNode, printing each node between min/max level constraints and incrementing raw counts for nodes at or under the target level.
	 * @param tn Taxonomy node to print
	 * @param tsw Output writer for taxonomy rows
	 */
	void printTaxonomy(TaxNode tn, final TextStreamWriter tsw){
//		assert(false) : tn.levelExtended+", "+taxLevelExtended+", "+minLevelExtended+", "+maxLevelExtended;
		assert(tn!=null);
//		tsw.print("\n");
		do{
			if(tn.levelExtended<=taxLevelExtended){tn.incrementRaw(1);}
			if(tn.levelExtended>=minLevelExtended && tn.levelExtended<=maxLevelExtended){
				if(!tn.cellularOrganisms() && (!skipNonCanonical || tn.isSimple())){
					tsw.println(tn.levelStringExtended(false)+"\t"+tn.id+"\t"+tn.name);
				}
			}
			tn=tree.getNode(tn.pid);
		}while(tn!=null && tn.id!=tn.pid);
	}
	
	/**
	 * Builds a semicolon-delimited taxonomy string with level prefixes (e.g., k__, p__) between the requested min/max levels, optionally reversing the order.
	 * @param tree Taxonomy tree used for parent traversal
	 * @param tn Starting taxonomy node
	 * @param minLevelE Minimum extended taxonomy level to include
	 * @param maxLevelE Maximum extended taxonomy level to include
	 * @param skipNonCanonical If true, omit non-simple nodes
	 * @param reverseOrder If true, emit levels from leaf to root
	 * @return Mutable StringBuilder containing the taxonomy line
	 */
	public static StringBuilder makeTaxLine(TaxTree tree, TaxNode tn, int minLevelE, int maxLevelE, boolean skipNonCanonical, boolean reverseOrder){
//		assert(false) : tn+", "+minLevelE+", "+maxLevelE;
		assert(tn!=null);
		StringBuilder sb=new StringBuilder();
		
		if(reverseOrder){
			ArrayList<TaxNode> list=new ArrayList<TaxNode>();
			while(tn.levelExtended<=maxLevelE){
				if(tn.levelExtended>=minLevelE){
					if(!tn.cellularOrganisms() && (!skipNonCanonical || tn.isSimple())){
						list.add(tn);
					}
				}
				if(tn.id==tn.pid){break;}
				tn=tree.getNode(tn.pid);
			}
			
			String semi="";
			Collections.reverse(list);
			for(TaxNode tn2 : list){
				sb.append(semi);
				sb.append(tn2.levelToStringShort());
				sb.append("__");
				sb.append(tn2.name);
				semi=";";
			}
		}else{
			String semi="";
			while(tn.levelExtended<=maxLevelE){
				if(tn.levelExtended>=minLevelE && !tn.cellularOrganisms() && (!skipNonCanonical || tn.isSimple())){
					sb.append(semi);
					sb.append(tn.levelToStringShort());
					sb.append("__");
					sb.append(tn.name);
					semi=";";
				}
				if(tn.id==tn.pid){break;}
				tn=tree.getNode(tn.pid);
			}
		}
		
		return sb;
	}
	
//	public static void printTaxonomy(TaxNode tn, final StringBuilder sb, final TaxTree tree, final int maxLevel, boolean skipNonCanonical){
//		final int maxLevelE=maxLevel<0 ? maxLevel : TaxTree.levelToExtended(maxLevel);
//		assert(tn!=null);
////		tsw.print("\n");
//		do{
//			if(!tn.cellularOrganisms() && (!skipNonCanonical || tn.isSimple())){
//				sb.append(tn.levelStringExtended(false)+"\t"+tn.id+"\t"+tn.name+"\n");
//			}
//			tn=tree.getNode(tn.pid);
//		}while(tn!=null && tn.id!=tn.pid && tn.levelExtended<=maxLevelE);
//	}
	
	/**
	 * Appends a full taxonomy hierarchy to a ByteBuilder with tab-separated columns, honoring the maxLevel constraint and optional canonical filtering.
	 * @param tn Starting taxonomy node
	 * @param sb Buffer to append level/id/name triples
	 * @param tree TaxTree instance used for parent lookups
	 * @param maxLevel Highest taxonomy level to include (negative uses extended levels)
	 * @param skipNonCanonical Whether to skip non-simple nodes
	 */
	public static void printTaxonomy(TaxNode tn, final ByteBuilder sb, final TaxTree tree, final int maxLevel, boolean skipNonCanonical){
		final int maxLevelE=maxLevel<0 ? maxLevel : TaxTree.levelToExtended(maxLevel);
		assert(tn!=null);
//		tsw.print("\n");
		do{
			if(!tn.cellularOrganisms() && (!skipNonCanonical || tn.isSimple())){
				sb.append(tn.levelStringExtended(false)).append('\t').append(tn.id).append('\t').append(tn.name).append('\n');
			}
			tn=tree.getNode(tn.pid);
		}while(tn!=null && tn.id!=tn.pid && tn.levelExtended<=maxLevelE);
	}
	
//	public static void printTaxonomy(TaxNode tn, final TextStreamWriter tsw, final TaxTree tree, final int maxLevel){
//		final int maxLevelE=maxLevel<0 ? maxLevel : TaxTree.levelToExtended(maxLevel);
//		assert(tn!=null);
////		tsw.print("\n");
//		do{
//			if(!skipNonCanonical || tn.isSimple()){
//				tsw.println(tn.levelStringExtended(false)+"\t"+tn.id+"\t"+tn.name);
//			}
//			tn=tree.getNode(tn.pid);
//		}while(tn!=null && tn.id!=tn.pid && tn.levelExtended<=maxLevelE);
//	}
	
	/**
	 * Walks a taxonomy node up to the configured target level, prints the node name, and increments its raw count for downstream summaries.
	 * @param tn Taxonomy node to resolve to the configured level (falls back to UNKNOWN)
	 * @param tsw Destination writer for the selected level
	 */
	void printTaxLevel(TaxNode tn, final TextStreamWriter tsw){
		if(tn==null){tn=unknown;}
		while(tn.id!=tn.pid && tn.levelExtended<taxLevelExtended){tn=tree.getNode(tn.pid);}
		if(tsw!=null){tsw.println(tn.name);}
		tn.incrementRaw(1);
	}
	
//	void printTaxCounts(TaxNode tn, final TextStreamWriter tsw){
//		if(tn==null){tn=unknown;}
//		while(tn!=null && tn.id!=tn.pid && tn.level<taxLevel){tn=tree.getNode(tn.pid);}
//		if(tsw!=null)tsw.println(tn.name);
//		tn.incrementRaw(1);
//	}
	
	/**
	 * Delegates to the TaxTree to parse gi/taxid/accession data from a header string, returning the matching TaxNode when available.
	 * @param header Sequence header containing taxonomy identifiers
	 * @return Matching TaxNode or null if no match is found
	 */
	public TaxNode parseNodeFromHeader(String header){
		if(tree==null){return null;}
		return tree.parseNodeFromHeader(header, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Optional input file path for staged taxonomy translation */
	private String in1=null;

	/** Primary output file path (defaults to stdout.txt) */
	private String out1="stdout.txt";
	
	/** Optional output path for aggregated taxonomy counts */
	private String countFile=null;

	/** Source file for gi-to-taxid mappings (auto uses TaxTree defaults) */
	private String giTableFile=null;
	/** Raw taxonomy tree file to load when constructing the TaxTree */
	private String taxTreeFile=null;
	/** Source file for accession-to-taxid mappings (optional) */
	private String accessionFile=null;
	
	/** Taxonomy tree instance used for all lookups and traversals */
	private final TaxTree tree;
	
//	/** Level to print */
//	private int taxLevel=-1;//TaxTree.stringToLevel("phylum");
//
//	/** Min level to print */
//	private int minLevel=-1;
//
//	/** Max level to print */
//	private int maxLevel=TaxTree.stringToLevel("life");
	
	private final int taxLevelExtended, minLevelExtended, maxLevelExtended;
	
	/** Whether taxonomy strings are emitted in reverse (leaf-to-root) order */
	private boolean reverseOrder=true;
	
	/**
	 * Names or identifiers supplied via CLI to process when no input file is given
	 */
	private ArrayList<String> names=new ArrayList<String>();
	
	/**
	 * Maximum number of reads to process when streaming sequence files (-1 for unlimited)
	 */
	private long maxReads=-1;
	
	/** Whether to echo the original sequence name ahead of taxonomy output */
	boolean printName=true;
	/** Skip non-canonical (non-simple) taxonomy nodes when building hierarchies */
	boolean skipNonCanonical=false;
	
	/**
	 * Zero-based column index for extracting taxonomy identifiers from tab-delimited inputs (-1 disables column mode)
	 */
	int keyColumn=-1;
//	Deprecated.  Description from shellscript:
//	column=-1       If set to a non-negative integer, parse the taxonomy
//            information from this column in a tab-delimited file.
//            Example if column=1:
//            read1 TAB gi|944259871|gb|KQL24128.1| TAB score:42
//            becomes
//            read1 TAB  k__Viridiplantae;p__Streptophyta;... TAB score:42
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input file descriptor (null when reading from CLI names) */
	private final FileFormat ffin1;
	
	/** Primary output file descriptor for taxonomy text */
	private final FileFormat ffout1;
	
	/** Writer descriptor for taxonomy count summaries (optional) */
	private final FileFormat ffcount;
	
	/**
	 * Fallback TaxNode used when a sequence cannot be resolved to any taxonomy entry
	 */
	private final TaxNode unknown=new TaxNode(-99, -99, TaxTree.LIFE, TaxTree.LIFE_E, "UNKNOWN");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status and logging messages */
	private PrintStream outstream=System.err;
	/** Enable verbose logging during processing */
	public static boolean verbose=false;
	/** Set to true if any worker thread or writer failed */
	public boolean errorState=false;
	/** Overwrite behavior for output writers */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of truncating */
	private boolean append=false;
	
}
