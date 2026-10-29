package tax;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
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
import structures.ByteBuilder;
import structures.IntList;
import tracker.ReadStats;

/**
 * Finds taxonomic ancestors and majority classifications for sequences based on GI numbers.
 * Processes input files containing GI numbers, converts them to NCBI taxonomy IDs, and determines common ancestors and majority classifications using a taxonomy tree.
 * @author Brian Bushnell
 * @date May 9, 2016
 */
public class FindAncestor{

	/** Launches GI-list ancestor classification and closes the status stream.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		FindAncestor x=new FindAncestor(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/** Parses parameters, prepares IO formats and loads GI/taxonomy lookups.
	 * @param args Command-line arguments */
	public FindAncestor(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		boolean useTree=true;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("tree") || a.equals("usetree") || a.equals("taxtree")){
				useTree=TaxTree.parseTreeFlag(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("invalid")){
				outInvalid=b;
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
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}

		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;

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
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);

		if(giTableFile!=null){
			outstream.println("Loading gi table.");
			GiToTaxid.initialize(giTableFile);
		}
		if(!useTree){
			throw new RuntimeException("FindAncestor requires a taxonomy tree; tree=f is not supported.");
		}
		tree=TaxTree.loadTaxTree(outstream, true, true);
		if(tree==null){
			throw new RuntimeException("FindAncestor requires a taxonomy tree; none could be loaded.");
		}
		assert(tree.nameMap!=null);
		lifeNode=tree.getNodeByName("life");
	}

	/** Processes all input lines, writing per-line LCA/majority classifications.
	 * @param t Timer tracking elapsed processing time */
	void process(Timer t){

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		bsw.print("#Name\tAncestor\tMajority\tTaxonomy...\n".getBytes());

		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}

		final IntList giList=new IntList();
		final IntList tidList=new IntList();
		final IntList traversal=new IntList();

		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();

		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=line.length;

				giList.clear();
				tidList.clear();
				traversal.clear();

				final int giCount=getGiNumbers(line, giList, ',');
				final int ncbiCount=getTaxidNumbers(giList, tidList);

				taxaCounted+=giCount;
				taxaValid+=ncbiCount;
				final boolean valid=(ncbiCount>0);

				if(valid){
					linesValid++;
					int ancestor=findAncestor(tidList);
					int majority=findMajority(tidList);

					for(int i=0; i<line.length && line[i]!='\t'; i++){
						bb.append(line[i]);
					}
					bb.tab();
					bb.append(ancestor);
					bb.tab();
					bb.append(majority);
					bb.tab();

					fillTraversal(majority, traversal, true);
					writeTraversal(traversal, bb);
					bb.nl();

					for(int i=0; i<tidList.size; i++){
						final int id=tidList.get(i);
						fillTraversal(id, traversal, true);
						writeTraversal(traversal, bb);
						bb.nl();
					}

					bsw.print(bb.toBytes());
					bb.clear();
				}else{
					if(bswInvalid!=null){
						bswInvalid.println(line);
					}
				}
			}
			line=bf.nextLine();
		}

		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));

		outstream.println();
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Fills a leaf-to-root path for one taxonomic ID.
	 * @param id Taxonomic ID to trace
	 * @param traversal Destination list, cleared before use
	 * @param addLife When true, include the root life node even when the lineage is
	 * nonempty; life is always emitted for an otherwise empty traversal */
	private void fillTraversal(int id, IntList traversal, boolean addLife){
		traversal.clear();
		for(TaxNode node=tree.getNode(id); node!=null && node!=lifeNode; node=tree.getNode(node.pid)){
			traversal.add(node.id);
		}
		if(addLife || traversal.size==0){traversal.add(lifeNode.id);}
	}

	/** Appends a root-to-leaf lineage path to a byte builder.
	 * @param traversal Leaf-to-root list of taxonomic IDs
	 * @param bb Destination builder */
	private void writeTraversal(IntList traversal, ByteBuilder bb){
		for(int i=traversal.size-1; i>=0; i--){
			final int id=traversal.get(i);
			if(id>=0){
				TaxNode tn=tree.getNode(id);
				bb.append(tn.name);
				if(i>0){bb.tab();}
			}
		}
	}

	/** Parses comma-delimited GI numbers from the tab-delimited line body.
	 * @param line Input line whose first tab-delimited column is skipped
	 * @param list Destination GI list, appended in input order
	 * @param delimiter GI-list delimiter
	 * @return Number of parsed GI IDs now present in {@code list} */
	private int getGiNumbers(final byte[] line, final IntList list, final char delimiter){
		int i=0;

		//Skip name
		while(i<line.length && line[i]!='\t'){i++;}

		//Skip whitespaces
		while(i<line.length && Character.isWhitespace(line[i])){i++;}

		while(i<line.length){
			while(i<line.length && line[i]==delimiter){i++;}
			int start=i;
			while(i<line.length && line[i]!=delimiter){i++;}
			final int stop=i;
			if(Tools.startsWith(line, prefix, start)){start+=3;}
			// RELAXED [tax/FindAncestor]: the old assert(start<stop)
			// rejected benign trailing/empty delimiter tokens such as
			// "gi1,gi2," or bare "gi|". GI numbers are legacy; best-effort
			// parsing is preferable, and the guard below skips empty tokens.
			if(start<stop){
				final int number=Parse.parseInt(line, start, stop);
				list.add(number);
			}
		}
		return list.size;
	}

	/** Converts GI IDs to taxonomic IDs, discarding missing mappings.
	 * @param giList GI numbers to look up
	 * @param ncbiList Destination taxonomic-ID list
	 * @return Total number of taxonomic IDs in {@code ncbiList} after appending */
	private static int getTaxidNumbers(final IntList giList, final IntList ncbiList){
		final int size=giList.size;
		for(int i=0; i<size; i++){
			final int gi=giList.get(i);
			final int ncbi=GiToTaxid.getID(gi);
			if(ncbi>=0){ncbiList.add(ncbi);}
		}
		return ncbiList.size;
	}

	/** Finds the common ancestor for this instance's taxonomy tree.
	 * @param list Nonempty valid taxonomic IDs from this instance's connected tree
	 * @return Common ancestor taxonomic ID */
	private int findAncestor(IntList list){
		return findAncestor(tree, list);
	}

	/** Finds the common ancestor of all taxonomic IDs in a list.
	 * Empty lists assert under {@code -ea} and return -1 otherwise. With
	 * assertions disabled, intermediate disconnected-pair results are ignored so
	 * later IDs can still be compared to the last connected ancestor.
	 * @param tree Connected taxonomy tree used for LCA checks
	 * @param list Nonempty valid taxonomic IDs to compare
	 * @return Last connected common ancestor taxonomic ID, or -1 for an empty list */
	public static int findAncestor(TaxTree tree, IntList list){
		if(list.size<1){
			assert(false);
			return -1;
		}
		int ancestor=list.get(0);
		for(int i=1; i<list.size && ancestor>-1; i++){
			final int id=list.get(i);
			int x=tree.commonAncestor(ancestor, id);
			if(x>-1){
				ancestor=x;
			}
		}
		return ancestor;
	}

	/** Finds the lowest-level taxon containing a strict majority of input IDs.
	 * Temporarily increments each traversed node's {@code countSum}; relevant nodes
	 * must start with zero {@code countSum} values and be exclusively owned during
	 * the call.
	 * @param list Taxonomic IDs to count
	 * @return Majority taxonomic ID by {@code levelExtended}, or the life node when
	 * none is more specific */
	private int findMajority(IntList list){
		if(list.size<3){return findAncestor(list);}
		final int majority=list.size/2+1;

		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			assert(tn!=null) : "No node for id "+id;
			if(tn!=null){
				tree.percolateUp(tn, 1);
			}
		}

		TaxNode best=lifeNode;
		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			while(tn!=null && tn!=lifeNode){
				if(tn.countSum>=majority && tn.levelExtended<best.levelExtended){
					best=tn;
					break;
				}
				tn=tree.getNode(tn.pid);
			}
		}


		for(int i=0; i<list.size; i++){
			final int id=list.get(i);
			TaxNode tn=tree.getNode(id);
			if(tn!=null){
				tree.percolateUp(tn, -1);
			}
		}

		return best.id;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Configuration        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input table containing record names and comma-delimited GI numbers. */
	private String in1=null;

	/** Primary output path for per-record ancestor and majority calls. */
	private String out1=null;

	/** Optional output path for records without valid mapped taxonomic IDs. */
	private String outInvalid=null;

	/** GI-to-taxID table path, or null to use the initialized default table. */
	private String giTableFile=null;

	/** Taxonomy tree used for all lineage and LCA operations. */
	private final TaxTree tree;

	/** Root node used to terminate and default lineage traversals. */
	private final TaxNode lifeNode;

	/*--------------------------------------------------------------*/

	/** GI IDs parsed from input. */
	private long taxaCounted=0;

	/** GI IDs that resolved to NCBI taxonomic IDs. */
	private long taxaValid=0;

	/** Nonempty input lines processed. */
	private long linesProcessed=0;

	/** Input lines containing at least one mapped taxonomic ID. */
	private long linesValid=0;

	/** Input bytes processed, excluding line terminators. */
	private long bytesProcessed=0;

	/** Maximum number of nonempty input lines to process; {@code maxLines<=0} is unlimited. */
	private long maxLines=Long.MAX_VALUE;

	/** Legacy placeholder for optional count-table output. */
	private boolean countTable=true;

	/** Prefix stripped from GI tokens when present. */
	private final String prefix="gi|";

	/*--------------------------------------------------------------*/
	/*----------------       Format Handles        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input file format descriptor. */
	private final FileFormat ffin1;

	/** Primary output file format descriptor. */
	private final FileFormat ffout1;

	/** Invalid-record output format descriptor, or null when disabled. */
	private final FileFormat ffoutInvalid;

	/*--------------------------------------------------------------*/
	/*----------------        Output State         ----------------*/
	/*--------------------------------------------------------------*/

	/** Status stream, potentially redirected by the pre-parser. */
	private PrintStream outstream=System.err;

	/** Emit verbose library diagnostics. */
	public static boolean verbose=false;

	/** True when the input close or any output close reports a failure. */
	public boolean errorState=false;

	/** Allow overwriting existing outputs. */
	private boolean overwrite=true;

	/** Append to existing outputs rather than replacing them. */
	private boolean append=false;

}
