package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Greedily selects query names from a tab-delimited ANI comparison table.
 * Edges are considered in descending ANI/KID/name order. For each still-valid pair,
 * the endpoint with the larger initial outgoing score is removed; names break ties.
 * Removal stops when the requested node count or ANI threshold is reached,
 * or the supplied edge list is exhausted.
 * This heuristic does not guarantee a globally optimal subset or an exact output count.
 * Unlike RepresentativeSet, it does not use taxonomy or centroids.
 * <p>The first three columns are query, reference and ANI. Query rows define nodes;
 * self comparisons define a node but do not contribute an edge. Sparse comparisons are
 * accepted, but a reference needed during selection or final reporting must have a query
 * row among the records read. The lines option can truncate that query coverage.
 * Recognized #Query/Ref/ANI or #Query/Ref/AAI headers locate optional KID and SSU columns;
 * headerless input uses zero-based positions 9 and 11. Omitted scores use zero.
 * ANI and KID are used in their input units, without percentage conversion.
 * <p>Each instance is intended for one processing run and is not thread-safe.
 * Construction changes shared I/O settings. Cleanup attempts to close every acquired
 * reader, writer and adopted diagnostic stream; borrowed System streams
 * remain open. A failure can leave partial output, including names written before the
 * final maximum-ANI scan fails.
 * @author Brian Bushnell
 * @date July 30, 2024
 */
public class PickSubset{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs selection and closes adopted diagnostics even when processing fails.
	 * @param args Command-line flag=value arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		PickSubset x=new PickSubset(args);
		Throwable failure=null;
		try{
			x.process(t);
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			x.closeDiagnostics(failure);
		}
	}

	/** Parses options, validates file paths and prepares format descriptors.
	 * Does not create the processing ByteFile or selected/discarded output writers.
	 * Diagnostics adopted during preprocessing are closed if construction fails.
	 * @param args Command-line flag=value arguments */
	public PickSubset(String[] args){
		boolean success=false;
		Throwable failure=null;
		try{
			{//Preprocess config files and diagnostic-stream selection.
				PreParser pp=new PreParser(args, /*getClass()*/null, false);
				args=pp.args;
				outstream=pp.outstream;
			}
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			ReadWrite.setZipThreads(Shared.threads());
			{
				final Parser parser=parse(args);
				overwrite=parser.overwrite;
				append=parser.append;
				in1=parser.in1;
				out1=parser.out1;
			}
			validateParams();
			fixExtensions();
			checkFileExistence();
			checkStatics();
			ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
			ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
			ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
			success=true;
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			//The caller cannot close adopted diagnostics when construction fails.
			if(!success){closeDiagnostics(failure);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses tool options and delegates shared options to Parser.
	 * Negative lines values become unlimited; zero is also unlimited during reading.
	 * @param args Preprocessed arguments
	 * @return Parser carrying shared input/output and overwrite/append settings */
	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			//Split once so filenames retain later equals signs.
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("invalid") || a.equals("outu")){
				outInvalid=b;
			}else if(a.equals("nodes") || a.equals("files") || a.equals("genomes")){
				desiredNodes=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("ani") || a.equalsIgnoreCase("maxani")){
				desiredANI=Float.parseFloat(b);
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Retained template option parses its value but has no operational effect.
			}else if(parser.parse(arg, a, b)){
				//Shared option handled by Parser.
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		return parser;
	}

	/** Resolves a supported input compression suffix and requires an input path. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/** Checks input readability, output writability and collisions before opening data streams. */
	private void checkFileExistence(){
		//Check both outputs before opening either writer or the input stream.
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outInvalid)){
			outstream.println(out1+", "+outInvalid);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+outInvalid+"\n");
		}

		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		if(!Tools.testForDuplicateFiles(true, in1, out1, outInvalid)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Selects BF2 for larger thread counts unless either ByteFile mode was explicitly forced. */
	private static void checkStatics(){
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}

	/** Asserts that at least one positive stopping criterion was supplied.
	 * @return Always true when the assertion passes or assertions are disabled */
	private boolean validateParams(){
		assert(desiredANI>0 || desiredNodes>0) : "Desired ANI or genomes must be greater than 0.";
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads comparisons, selects query nodes and writes kept/discarded names.
	 * Reports the largest supplied ANI between retained nodes, with zero as its floor.
	 * Acquired data streams are closed before normal statistics are printed. A failure
	 * during selection or the final ANI scan can leave partial output files.
	 * @param t Running timer started before construction */
	void process(Timer t){
		ByteFile bf=null;
		ByteStreamWriter bsw=null, bswInvalid=null;
		ArrayList<Node> chosen;
		float maxANI=0;
		Throwable failure=null;
		try{
			bf=ByteFile.makeByteFile(ffin1);
			bsw=makeBSW(ffout1);
			bswInvalid=makeBSW(ffoutInvalid);
			ArrayList<Edge> edges=readFile(bf);
			chosen=pickNodes(edges, desiredNodes);
			if(bsw!=null){
				for(Node n : chosen){
					bsw.println(n.name);
				}
			}
			if(bswInvalid!=null){
				for(Entry<String, Node> e : nodeMap.entrySet()){
					Node n=e.getValue();
					if(!n.valid){
						bswInvalid.println(n.name);
					}
				}
			}
			for(Node n : chosen){
				for(Edge e : n.edges){
					Node r=referenceNode(e);
					if(r.valid){
						maxANI=Tools.max(maxANI, e.ani);
					}
				}
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			finishStreams(bf, bsw, bswInvalid, failure);
		}

		t.stop();

		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));

		outstream.println();
		outstream.println("Nodes In:          \t"+nodeMap.size());
		outstream.println("Nodes Out:         \t"+chosen.size());
		outstream.println("Maximum ANI:       \t"+maxANI);

		//Close methods can report asynchronous I/O errors without throwing.
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/** Closes each acquired resource even if an earlier cleanup throws.
	 * Preserves the processing failure, or throws the first unchecked cleanup failure.
	 * Factory-internal acquisition and underlying blocking I/O retain their own limits.
	 * @param bf Acquired input, or null
	 * @param bsw Started selected-output writer, or null
	 * @param bswInvalid Started discarded-output writer, or null
	 * @param primary Existing processing failure, or null after normal processing */
	private void finishStreams(final ByteFile bf, final ByteStreamWriter bsw,
			final ByteStreamWriter bswInvalid, final Throwable primary){
		Throwable failure=primary;
		if(bf!=null){
			try{errorState|=bf.close();}
			catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		}
		if(bsw!=null){
			try{errorState|=bsw.poisonAndWait();}
			catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		}
		if(bswInvalid!=null){
			try{errorState|=bswInvalid.poisonAndWait();}
			catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		}
		if(primary==null && failure!=null){
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			throw (Error)failure;
		}
	}

	/** Closes adopted diagnostics while preserving System streams and a primary failure.
	 * @param primary Earlier failure, or null on normal completion */
	private void closeDiagnostics(final Throwable primary){
		try{Shared.closeStream(outstream);}
		catch(RuntimeException | Error e){
			if(primary==null){throw e;}
			retainFailure(primary, e);
		}
	}

	/** Keeps the first failure and suppresses a distinct later cleanup failure.
	 * @param primary Earlier failure, or null
	 * @param secondary Nonnull cleanup failure
	 * @return The earlier failure when present, otherwise the cleanup failure */
	private static Throwable retainFailure(final Throwable primary, final Throwable secondary){
		if(primary==null){return secondary;}
		if(primary!=secondary){primary.addSuppressed(secondary);}
		return primary;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads query nodes and directed non-self comparison edges.
	 * A query's first row supplies its name and optional size/base metadata. Self rows
	 * still create nodes, but their ANI/KID values do not enter the retained-edge check.
	 * Retained non-self ANI and KID scores must be finite. Recognized headers update the
	 * optional-column layout for subsequent rows; unrelated comments leave it unchanged.
	 * <p>The positive lines limit counts nonempty records, including comments and headers.
	 * Blank records do not count. Reported bytes sum each counted record's length plus one,
	 * rather than measuring physical file bytes or compression size.
	 * @param bf Acquired input, borrowed here and closed by process
	 * @return Mutable list of retained non-self edges, also referenced by their query nodes */
	private ArrayList<Edge> readFile(ByteFile bf){
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();//Retained template allocation; currently unused.

		ArrayList<Edge> edges=new ArrayList<Edge>();
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);

				if(line[0]!='#'){
					Edge e=new Edge(line);
					Node query=nodeMap.get(e.query);
					//Query rows define nodes; referenceNode diagnoses unavailable references when needed.
					if(query==null){
						query=new Node(line);
						nodeMap.put(e.query, query);
					}
					if(!e.query.equals(e.ref)){//Self rows define nodes without adding edges.
						if(!Float.isFinite(e.ani) || !Float.isFinite(e.kid)){
							throw new IllegalArgumentException("Nonfinite ANI or KID on input line "+linesProcessed+
									": "+e.query+" -> "+e.ref+", ANI="+e.ani+", KID="+e.kid);
						}
						query.edges.add(e);
						edges.add(e);
					}
				}else{
					parseHeader(line);
				}
			}
			line=bf.nextLine();
		}
		return edges;
	}

	/** Reads optional score positions from a producer query/reference header.
	 * Unrelated comments leave the current layout unchanged; omitted optional columns use zero.
	 * @param line Comment or header bytes supplied by readFile
	 * @throws IllegalArgumentException If a recognized header repeats a score column */
	private void parseHeader(final byte[] line){
		assert(line!=null && line.length>0 && line[0]=='#') : "readFile routes only comment/header lines to parseHeader";
		lp.set(line);
		if(lp.terms()<3 || !lp.termEquals("#Query", 0) || !lp.termEquals("Ref", 1) ||
				(!lp.termEquals("ANI", 2) && !lp.termEquals("AAI", 2))){return;}
		int kid=-1, ssu=-1;
		for(int column=3; column<lp.terms(); column++){
			if(lp.termEquals("KID", column)){
				if(kid>=0){throw new IllegalArgumentException("Duplicate KID column in input header on line "+linesProcessed);}
				kid=column;
			}else if(lp.termEquals("SSU", column)){
				if(ssu>=0){throw new IllegalArgumentException("Duplicate SSU column in input header on line "+linesProcessed);}
				ssu=column;
			}
		}
		kidColumn=kid;
		ssuColumn=ssu;
	}

	/** Greedily invalidates endpoints of the highest-ranked remaining edges.
	 * Computes each node's outgoing score once, then sorts and reverses the supplied list.
	 * For a valid pair, removes the endpoint with the larger score/name comparison.
	 * Scores are not recomputed as nodes disappear. Either stopping criterion ends removal:
	 * remaining nodes at most desiredNodes, or the next edge's ANI at most desiredANI.
	 * Exhausting the supplied edges also ends selection. The result may contain more
	 * than the requested number of nodes, or fewer if the input defines fewer nodes.
	 * @param list Mutable edge list; reordered and its edge-valid flags consumed
	 * @param desiredNodes Count-based stopping threshold
	 * @return Retained nodes in nodeMap iteration order, not sorted or input order */
	ArrayList<Node> pickNodes(ArrayList<Edge> list, int desiredNodes){
		int remaining=nodeMap.size();
		for(Node n : nodeMap.values()){n.calcSum();}
		Collections.sort(list);
		Collections.reverse(list);
		for(Edge e : list){
			if(remaining<=desiredNodes || e.ani<=desiredANI){break;}
			if(e.valid){
				Node q=nodeMap.get(e.query);
				Node r=referenceNode(e);
				if(q.valid && r.valid){
					if(q.compareTo(r)>=0){
						q.setInvalid();
					}else{
						r.setInvalid();
					}
					remaining--;
				}
				e.valid=false;
			}
		}
		ArrayList<Node> chosen=new ArrayList<Node>(remaining);
		for(Node n : nodeMap.values()){
			if(n.valid){chosen.add(n);}
		}
		return chosen;
	}

	/** Resolves a comparison's reference among query nodes already read.
	 * @param edge Stored comparison supplied by process or pickNodes
	 * @return Existing reference node, which may already be invalidated
	 * @throws IllegalArgumentException If no query row was read for the reference */
	private Node referenceNode(final Edge edge){
		assert(edge!=null) : "process and pickNodes must pass stored non-null comparison edges";
		final Node node=nodeMap.get(edge.ref);
		if(node==null){
			throw new IllegalArgumentException("No query row was read for reference '"+edge.ref+
					"' (query '"+edge.query+"'). Check input coverage and lines= limits.");
		}
		return node;
	}

	/** Creates and starts an optional name-output writer.
	 * @param ff Output descriptor, or null to disable this output
	 * @return Started writer owned by process, or null when disabled */
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}

	/** Query-name lookup and mutable node state; reference-only names are not inserted. */
	private HashMap<String, Node> nodeMap=new HashMap<String, Node>();

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Query-defined graph node with its original outgoing edges and cached ranking score. */
	private class Node implements Comparable<Node>{

		/** Takes metadata from the first row encountered for this query.
		 * Columns 3 and 5, when present, are parsed as longs even though selection does not
		 * use them. Missing metadata defaults to zero.
		 * @param line Tab-delimited comparison containing the query name in column 0 */
		Node(byte[] line){
			lp.set(line);
			name=lp.parseString(0);
			//The compact input format contains only query, reference and ANI.
			size=(lp.terms()>3 ? lp.parseLong(3) : 0);
			bases=(lp.terms()>5 ? lp.parseLong(5) : 0);
		}

		/** Caches the initial sum of outgoing KID plus ANI weighted by 0.0001.
		 * Includes all stored edges, without consulting their validity. The score is
		 * computed once before selection and is not updated after removals.
		 * @return Cached float score used to choose which endpoint to remove
		 * @throws IllegalArgumentException If the accumulated ranking score is nonfinite */
		float calcSum(){
			sum=0;
			for(Edge e : edges){
				sum+=e.kid;
				sum+=e.ani*0.0001;//In case kid is absent
			}
			if(!Float.isFinite(sum)){
				throw new IllegalArgumentException("Nonfinite accumulated score for node "+name+": "+sum);
			}
			return sum;
		}

		/** Invalidates this node and its outgoing edges when removed from selection.
		 * Incoming edges are left stored; selection checks the reference node's validity. */
		void setInvalid(){
			assert(valid);
			valid=false;
			for(Edge e : edges){e.valid=false;}
		}

		/** Orders finite initial scores ascending, then breaks ties by query name.
		 * calcSum rejects nonfinite totals before the selection loop compares nodes.
		 * @param o Other node with its initial score already calculated
		 * @return Negative, zero or positive according to score then name ordering */
		@Override
		public int compareTo(Node o){
			if(sum!=o.sum){return sum>o.sum ? 1 : -1;}
			return name.compareTo(o.name);
		}


		/** Node name/identifier from input file */
		final String name;
		/** Node size metadata from input, or zero when absent; not used for ranking. */
		long size;
		/** Node base-count metadata from input, or zero when absent; not used for ranking. */
		long bases;
		/** Initial outgoing score, calculated before selection and not updated by removals. */
		float sum;
		/** Whether this node is still valid (not removed from subset) */
		boolean valid=true;
		/** Directed non-self edges for which this node is the query. */
		ArrayList<Edge> edges=new ArrayList<Edge>();
	}

//	#Query	Ref	ANI	QSize	RefSize	QBases	RBases	QTaxID	RTaxID	KID	WKID	SSU
//	a.fa	b.fa	80.980	11178684	8473676	12055201	8514625	-1	-1	0.237	0.313	94.233

	/** One directed comparison, sharing the outer parser and optional-column layout. */
	private class Edge implements Comparable<Edge>{

		/** Parses fixed query/reference/ANI columns and optional KID/SSU scores.
		 * Missing optional scores use zero; a lone dot is also accepted for SSU.
		 * readFile performs the retained non-self edge's finite ANI/KID check afterward.
		 * @param line Tab-delimited comparison with at least the first three columns */
		Edge(byte[] line){
			lp.set(line);
			query=lp.parseString(0);
			ref=lp.parseString(1);
			ani=lp.parseFloat(2);
			//Recognized headers map optional scores; headerless input uses the legacy positions.
			kid=(kidColumn>=0 && lp.terms()>kidColumn) ? lp.parseFloat(kidColumn) : 0;
			//An absent SSU column or dot uses zero.
			if(ssuColumn<0 || lp.terms()<=ssuColumn || lp.termEquals('.', ssuColumn)){
				ssu=0;
			}else{
				ssu=lp.parseFloat(ssuColumn);
			}
		}

		/** Orders edges ascending by ANI, KID, query name and reference name.
		 * readFile rejects nonfinite scores before non-self edges enter the sortable list.
		 * pickNodes reverses this ordering before selecting pairs to remove.
		 * @param o Other retained edge with finite ANI and KID
		 * @return Negative, zero or positive according to the first differing key */
		@Override
		public int compareTo(Edge o){
			if(ani!=o.ani){return ani>o.ani ? 1 : -1;}
			if(kid!=o.kid){return kid>o.kid ? 1 : -1;}
			int x=query.compareTo(o.query);
			if(x!=0){return x;}
			return ref.compareTo(o.ref);
		}

		/** Query sequence name in the ANI comparison */
		final String query;
		/** Reference sequence name in the ANI comparison */
		final String ref;
		/** Average Nucleotide Identity value between query and reference */
		float ani;
		/** K-mer identity in input units; used for edge ordering and node scores. */
		float kid;
		/** Parsed SSU similarity; currently unused by selection and reporting. */
		float ssu;
		/** Whether selection may still consider this edge; consumed or query-invalidated edges are false. */
		boolean valid=true;
	}

	/** Tab-delimited line parser for processing input file lines */
	private LineParser1 lp=new LineParser1('\t');
	/** KID position in the current layout, or -1 when a recognized header omits it. */
	private int kidColumn=9;
	/** SSU position in the current layout, or -1 when a recognized header omits it. */
	private int ssuColumn=11;

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Optional output path for names of nodes removed by selection. */
	private String outInvalid=null;

	/** Count threshold: removal stops once the remaining count is at most this value. */
	private int desiredNodes=0;
	/** Edge-ANI stopping threshold in input units; may stop removal before the count target. */
	private float desiredANI=0;

	/*--------------------------------------------------------------*/

	/** Number of edges processed (currently unused) */
	private long edgesProcessed=0;
	/** Number of nodes processed (currently unused) */
	private long nodesProcessed=0;
	/** Number of edges output (currently unused) */
	private long edgesOut=0;
	/** Number of nodes output (currently unused) */
	private long nodesOut=0;

	/** Nonempty input records counted, including comments and headers. */
	private long linesProcessed=0;
	/** Number of output lines written (currently unused) */
	private long linesOut=0;
	/** Sum of counted record lengths plus one per record; not physical/compressed file size. */
	private long bytesProcessed=0;
	/** Number of output bytes written (currently unused) */
	private long bytesOut=0;

	/** Positive limit on counted nonempty records; zero disables the limit. */
	private long maxLines=Long.MAX_VALUE;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Comparison-input descriptor; its ByteFile is acquired by process. */
	private final FileFormat ffin1;
	/** Optional selected-name output descriptor. */
	private final FileFormat ffout1;
	/** Optional discarded-name output descriptor. */
	private final FileFormat ffoutInvalid;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostics, borrowed from System.err by default or adopted during preprocessing. */
	private PrintStream outstream=System.err;
	/** Shared verbosity flag, also copied into ByteFile1, ByteFile2 and ReadWrite when parsed. */
	public static boolean verbose=false;
	/** Aggregates errors reported by acquired reader/writer close methods. */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;

}
