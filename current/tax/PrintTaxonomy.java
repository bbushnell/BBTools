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
 * Resolves command-line names, text columns, or read headers through local GI,
 * accession, and taxonomy-tree data, then prints either full lineages or a
 * selected taxonomy level.
 * <p>
 * Static helpers also format taxonomy lineages for sketch and clade tools.
 *
 * @author Brian Bushnell
 * @date November 23, 2015
 */
public class PrintTaxonomy{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs taxonomy printing with the supplied command-line options.
	 * @param args BBTools-style {@code flag=value} arguments plus optional bare names or IDs */
	public static void main(String[] args){
		Timer t=new Timer();
		PrintTaxonomy x=new PrintTaxonomy(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/** Parses arguments, initializes optional GI/accession tables, and loads the required {@link TaxTree}.
	 * @param args BBTools-style {@code flag=value} arguments plus optional bare names or IDs
	 * @throws IllegalArgumentException If the primary output destination is disabled */
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
		boolean useTree=true;

		int taxLevel=0, minLevel=0, maxLevel=TaxTree.LIFE;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("tree") || a.equals("usetree") || a.equals("taxtree")){
				useTree=TaxTree.parseTreeFlag(b);
			}else if(a.equals("out")){
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
		if(ffout1==null){
			throw new IllegalArgumentException("PrintTaxonomy requires out=<file> or out=stdout.txt; out=null is not supported.");
		}

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
		if(!useTree){
			throw new RuntimeException("PrintTaxonomy requires a taxonomy tree; tree=f is not supported.");
		}
		tree=TaxTree.loadTaxTree(outstream, true, true);
		if(tree==null){
			throw new RuntimeException("PrintTaxonomy requires a taxonomy tree; none could be loaded.");
		}
		assert(tree.nameMap!=null);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Dispatches to read, text-column, or command-line-name processing and writes optional count summaries.
	 * Finishes the started primary writer even when processing fails, preserving
	 * the processing exception if writer teardown also throws.
	 * @param t active timer for the processing phase */
	void process(Timer t){

		TextStreamWriter tsw=null;
		if(ffout1!=null){
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
		}

		Throwable failure=null;
		try{
			if(ffin1!=null){
				if(ffin1.fasta() || ffin1.fastq() || ffin1.samOrBam() || ffin1.scarf()){
					processReads(tsw);
				}else{
					processFile(new TextFile(ffin1), tsw);
				}
			}else{
				processNames(tsw);
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			if(tsw!=null){
				try{
					errorState|=tsw.poisonAndWait();
				}catch(RuntimeException | Error e){
					errorState=true;
					if(failure==null){throw e;}
					if(failure!=e){failure.addSuppressed(e);}
				}
			}
		}

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

		//Throw an exception if there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Prints taxonomy output for every name or ID collected from the command line.
	 * @param tsw non-null output writer for taxonomy results */
	void processNames(final TextStreamWriter tsw){
		for(String name : names){
			if(taxLevelExtended>0){
				printTaxLevel(name, tsw);
			}else{
				printTaxonomy(name, tsw);
			}
		}
	}

	/** Streams text input, translating one tab-delimited column or resolving each full line as an identifier.
	 * Closes the owned input on normal and exceptional exits, accumulating reported
	 * close errors. A thrown close failure is suppressed behind a processing failure.
	 * @param tf owned input text file containing names or identifiers
	 * @param tsw output writer for taxonomy results */
	void processFile(final TextFile tf, final TextStreamWriter tsw){
		Throwable failure=null;
		try{
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
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			try{
				errorState|=tf.close();
			}catch(RuntimeException | Error e){
				errorState=true;
				if(failure==null){throw e;}
				if(failure!=e){failure.addSuppressed(e);}
			}
		}
	}

	/** Resolves taxonomy from read headers and prints hierarchy or level data per read.
	 * Processes only first-mate headers from the input stream.
	 * Returns the current borrowed batch before attempting input closure on every
	 * exit after startup. Reported close errors accumulate in errorState; unchecked
	 * cleanup failures are suppressed behind an earlier processing failure.
	 * @param tsw output writer for taxonomy results */
	void processReads(final TextStreamWriter tsw){
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){System.err.println("Started cris");}
			cris.start();
		}

		ListNum<Read> ln=null;
		Throwable failure=null;
		try{
			ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){
				for(Read r1 : reads){
					if(keyColumn>=0){
						String result=translateLine(r1.id, keyColumn);
						tsw.print(result);
					}else if(taxLevelExtended>0){
						printTaxLevel(r1.id, tsw);
					}else{
						printTaxonomy(r1.id, tsw);
					}
				}
				final ListNum<Read> returned=ln;
				ln=null;//Do not retry a batch return if returnList itself throws.
				cris.returnList(returned);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			try{
				cris.returnList(ln);
			}catch(RuntimeException | Error e){
				errorState=true;
				if(failure==null){failure=e; throw e;}
				if(failure!=e){failure.addSuppressed(e);}
			}finally{
				try{
					errorState|=ReadWrite.closeStreams(cris);
				}catch(RuntimeException | Error e){
					errorState=true;
					if(failure==null){throw e;}
					if(failure!=e){failure.addSuppressed(e);}
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces the selected tab-delimited column with a semicolon lineage or a {@code _***NOT_FOUND***} suffix.
	 * @param line input line containing tab-delimited tokens
	 * @param col zero-based column index holding the sequence identifier
	 * @return line rebuilt with the selected column translated */
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
					//FIXED [tax/PrintTaxonomy#001]: multi-match now REPLACES split[col] with the tab-joined taxlines, matching
					//the single-match branch above (split[col]=tl). It previously appended to sb while leaving split[col] as the
					//original name, so the final column-rebuild emitted the taxlines PLUS the full unreplaced original line.
					StringBuilder sb2=new StringBuilder();
					String tab="";
					for(TaxNode tn2 : list){
						sb2.append(tab);
						sb2.append(makeTaxLine(tree, tn2, minLevelExtended, maxLevelExtended, skipNonCanonical, reverseOrder).toString());
						tab="\t";
					}
					split[col]=sb2.toString();
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

	/** Resolves a sequence name to one or more taxonomy nodes and prints a full hierarchy for each hit.
	 * Unresolved names print a not-found message rather than using {@link #unknown}.
	 * @param name sequence name or identifier, with or without FASTA/FASTQ prefix markers
	 * @param tsw output writer that receives the formatted taxonomy */
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

	/** Resolves a sequence name to one or more taxonomy nodes and prints the configured target level for each hit.
	 * Unresolved names print a not-found message rather than using {@link #unknown}.
	 * @param name sequence name or identifier, with or without FASTA/FASTQ prefix markers
	 * @param tsw output writer for taxonomy results */
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


	/** Traverses a node toward the root, printing nodes within the configured min/max levels.
	 * Increments raw counts for every visited node at or below the target level,
	 * independently of min/max and canonical-output filtering.
	 * @param tn taxonomy node to print
	 * @param tsw output writer for taxonomy rows */
	void printTaxonomy(TaxNode tn, final TextStreamWriter tsw){
		assert(tn!=null);
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

	/** Builds a semicolon-delimited lineage with short level prefixes such as {@code p__}.
	 * @param tree taxonomy tree used for parent traversal
	 * @param tn starting taxonomy node
	 * @param minLevelE minimum extended taxonomy level to include
	 * @param maxLevelE maximum extended taxonomy level to include
	 * @param skipNonCanonical true to omit non-simple nodes
	 * @param reverseOrder true to emit from root to leaf instead of leaf to root
	 * @return mutable builder containing the lineage */
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


	/**
	 * Appends tab-delimited level, taxID, and name rows while traversing a node
	 * toward the root. Always considers the starting node; {@code maxLevel}
	 * bounds only parent traversal. This intentional best-effort exception keeps
	 * the known starting classification even when it is already above the limit;
	 * traversal cannot descend to infer a more specific classification.
	 * @param tn starting taxonomy node
	 * @param sb buffer to receive level/id/name triples
	 * @param tree taxonomy tree used for parent lookups
	 * @param maxLevel taxonomy level whose extended value bounds parent traversal; negative values stop after the starting node
	 * @param skipNonCanonical true to omit non-simple nodes */
	public static void printTaxonomy(TaxNode tn, final ByteBuilder sb, final TaxTree tree, final int maxLevel, boolean skipNonCanonical){
		final int maxLevelE=maxLevel<0 ? maxLevel : TaxTree.levelToExtended(maxLevel);
		assert(tn!=null);
		do{
			if(!tn.cellularOrganisms() && (!skipNonCanonical || tn.isSimple())){
				sb.append(tn.levelStringExtended(false)).append('\t').append(tn.id).append('\t').append(tn.name).append('\n');
			}
			tn=tree.getNode(tn.pid);
		}while(tn!=null && tn.id!=tn.pid && tn.levelExtended<=maxLevelE);
	}


	/** Prints and counts a best-effort match for the configured target level.
	 * Walks toward the root until the first node at or above the requested rank,
	 * or the root if no such node is found. Missing ranks in irregular lineages
	 * therefore resolve to a broader ancestor. A starting node already at or above
	 * the requested rank is retained; this method never descends the tree.
	 * @param tn taxonomy node to resolve, or null from direct helper calls to count {@link #unknown}
	 * @param tsw destination writer for the selected level */
	void printTaxLevel(TaxNode tn, final TextStreamWriter tsw){
		if(tn==null){tn=unknown;}
		while(tn.id!=tn.pid && tn.levelExtended<taxLevelExtended){tn=tree.getNode(tn.pid);}
		if(tsw!=null){tsw.println(tn.name);}
		tn.incrementRaw(1);
	}


	/** Delegates GI, taxID and accession extraction from a header string to the loaded taxonomy tree.
	 * @param header sequence header containing a supported taxonomy identifier
	 * @return matching taxonomy node, or null when no match is found */
	public TaxNode parseNodeFromHeader(String header){
		if(tree==null){return null;}
		return tree.parseNodeFromHeader(header, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Optional input file path for text or sequence-file taxonomy translation. */
	private String in1=null;

	/** Primary output file path; the default is an explicit stdout descriptor. */
	private String out1="stdout.txt";

	/** Optional output path for per-node raw count summaries. */
	private String countFile=null;

	/** Source file for GI-to-taxID mappings; {@code auto} uses TaxTree defaults. */
	private String giTableFile=null;
	/** Optional accession-to-taxID source path. */
	private String accessionFile=null;

	/** Taxonomy tree used for all lookups and traversals. */
	private final TaxTree tree;


	/** Extended target level requested for output filtering. */
	private final int taxLevelExtended;
	/** Minimum extended level requested for output filtering. */
	private final int minLevelExtended;
	/** Maximum extended level requested for output filtering. */
	private final int maxLevelExtended;

	/** Emit semicolon lineages from root to leaf instead of leaf to root. */
	private boolean reverseOrder=true;

	/** Names or identifiers supplied on the command line when no input file is given. */
	private ArrayList<String> names=new ArrayList<String>();

	/** Maximum number of root fragments to request from CRIS, or -1 for unlimited. */
	private long maxReads=-1;

	/** Echo the original sequence name ahead of taxonomy output. */
	boolean printName=true;
	/** Skip non-canonical, non-simple taxonomy nodes when building hierarchies. */
	boolean skipNonCanonical=false;

	/** Zero-based column index for tab-delimited translation, or -1 to resolve each full line. */
	int keyColumn=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Input file descriptor; null when reading command-line names. */
	private final FileFormat ffin1;

	/** Primary taxonomy-text output descriptor. */
	private final FileFormat ffout1;

	/** Optional taxonomy-count output descriptor. */
	private final FileFormat ffcount;

	/** Fallback node counted by direct null-node {@link #printTaxLevel(TaxNode, TextStreamWriter)} calls. */
	private final TaxNode unknown=new TaxNode(-99, -99, TaxTree.LIFE, TaxTree.LIFE_E, "UNKNOWN");

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for status and logging messages. */
	private PrintStream outstream=System.err;
	/** Enable verbose logging during processing. */
	public static boolean verbose=false;
	/** True if a checked input or output close reported an error. */
	public boolean errorState=false;
	/** Overwrite existing output files. */
	private boolean overwrite=true;
	/** Append to existing output files instead of truncating them. */
	private boolean append=false;

}
