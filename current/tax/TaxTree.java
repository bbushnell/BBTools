package tax;

import java.io.File;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.TextFile;
import structures.IntHashMap;
import structures.IntLongHashMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import template.ThreadWaiter;

/**
 * NCBI taxonomy topology, rank indexes, name lookup and mutable abundance counters.
 * Durable TaxIDs normally index nodes directly; auxiliary IDs instead use a packed
 * mapping, so clients must use getNode rather than indexing nodes by an arbitrary ID.
 * Parent ranks need not increase monotonically. Counts and optional indexes are
 * mutable: complete initialization before sharing read-only use, and externally
 * coordinate mutations. The process-wide shared-tree cache reuses the first tree.
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxTree implements Serializable{
	
	/** Compatibility identifier for legacy serialized taxonomy resources. */
	private static final long serialVersionUID=5894416521711540017L;
	
	/*--------------------------------------------------------------*/
	/*----------------              Main            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * This is not called normally, only when converting NCBI text files
	 * into a portable TSV or legacy serialized representation.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		assert(args.length>=4) : "TaxTree syntax:\ntaxtree.sh names.dmp nodes.dmp merged.dmp taxtree.tsv.gz\n";
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=(Shared.threads()>2 ? 11 : 9);
		ReadWrite.PIGZ_BLOCKSIZE=256;
		ReadWrite.PIGZ_ITERATIONS=60;
		
		Timer t=new Timer();
		TaxTree tree=new TaxTree(args);
		t.stop();
		
		outstream.println("Retained "+tree.nodeCount+" nodes:");
		
		for(int i=tree.treeLevelsExtended.length-1; i>=0; i--){
			outstream.print(tree.nodesPerLevelExtended[i]+"\t"+taxLevelNamesExtended[i]);
			if(verbose){
				int lim=10;
				for(int j=0; j<lim && j<tree.treeLevelsExtended[i].length; j++){
					TaxNode n=tree.treeLevelsExtended[i][j];
					outstream.print("\n"+n+" -> "+tree.getNode(n.pid));
				}
				for(int j=tree.treeLevelsExtended[i].length-lim; j<tree.treeLevelsExtended[i].length; j++){
					if(j>=lim){
						TaxNode n=tree.treeLevelsExtended[i][j];
						outstream.print("\n"+n+" -> "+tree.getNode(n.pid));
					}
				}
			}
			outstream.println();
		}
		
		
		outstream.println();
		outstream.println("Time: \t"+t);
		
		if(args.length>2){//Write a tree
			writeTaxTree(tree, args[3], true);
		}
	}
	
	/** Parses shared flags while allowing the first four positional filenames.
	 * @param args Original command-line arguments; optional flags follow the filenames
	 * @return Parsed shared settings; tree-specific state is also updated */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(i>3){
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor using filenames from command line arguments, in the format of:
	 * {names, nodes, merged}
	 * @param args Command line arguments
	 */
	private TaxTree(String[] args){
		this(args[0], args[1], args[2], args);
	}
	
	/**
	 * Builds a tree from NCBI dump files, then infers ranks and builds rank indexes.
	 * Existing rank inversions can be reported without rejecting the entire tree.
	 * @param namesFile NCBI names.dmp scientific-name records
	 * @param nodesFile NCBI nodes.dmp parent and rank records
	 * @param mergedFile NCBI merged.dmp old-to-current IDs, or null
	 * @param args Optional shared command-line settings, or null
	 */
	private TaxTree(String namesFile, String nodesFile, String mergedFile, String[] args){
		
		if(args!=null){
			Parser parser=parse(args);
		}
		
		nodes=getNames(namesFile);
		getNodes(nodesFile, nodes);
		
		mergedMap=getMerged(mergedFile);
		externalToInternalMap=null;
		firstAuxIndex=nodes.length;
		
		countChildren();
		outstream.println("Counted children.");
		int rounds=percolate();
		outstream.println("Percolated "+rounds+" rounds to fixpoint.");
		
		if(assignStrains){
			assignStrains();
			rounds=percolate();
			outstream.println("Percolated "+rounds+" rounds to fixpoint.");
		}
		
		if(simplify){
			if(verbose){outstream.println("Simplifying.");}
			int removed=simplify(nodes);
			if(verbose){outstream.println("Removed "+removed+" nodes.");}
			rounds=percolate();
			outstream.println("Percolated "+rounds+" rounds to fixpoint.");
		}
		int errors=test(nodes);
//		assert(errors==0); //Not possible since the tree is wrong.
		if(errors>0){
			System.err.println("Found "+errors+" errors in tree.");
		}
		
		nodeCount=populateLevelData();

	}

	/** Adopts prepared native-ID nodes and builds rank indexes with one thread.
	 * Parents, child counts and inferred bounds must already be populated. Non-null
	 * nodes must occupy their TaxID slot and contain valid normal/extended ranks.
	 * Arrays, nodes and the merge map are retained, not copied. Configuration flags
	 * are recorded without rerunning simplification or inference.
	 * @param nodes_ Prepared sparse node array
	 * @param mergedMap_ Old-to-current TaxID mapping, or null
	 * @param minValidTaxa_ Legacy minimum-count metadata
	 * @param simplify_ Recorded simplification setting
	 * @param reassign_ Recorded no-rank reassignment setting
	 * @param skipNorank_ Recorded no-rank removal setting
	 * @param inferRankLimit_ Recorded rank-inference limit */
	public TaxTree(TaxNode[] nodes_, IntHashMap mergedMap_, int minValidTaxa_, boolean simplify_,
			boolean reassign_, boolean skipNorank_, int inferRankLimit_){
		this(nodes_, mergedMap_, minValidTaxa_, simplify_, reassign_, skipNorank_, inferRankLimit_, 1, null, nodes_.length);
	}

	/** Adopts prepared native-ID nodes with optional parallel rank indexing.
	 * @param nodes_ Prepared sparse array, indexed by native TaxID
	 * @param mergedMap_ Old-to-current TaxID mapping, or null
	 * @param minValidTaxa_ Legacy minimum-count metadata
	 * @param simplify_ Recorded simplification setting, not an instruction to simplify
	 * @param reassign_ Recorded no-rank reassignment setting
	 * @param skipNorank_ Recorded no-rank removal setting
	 * @param inferRankLimit_ Recorded rank-inference limit
	 * @param threads Requested indexing workers; small arrays use the serial path */
	TaxTree(TaxNode[] nodes_, IntHashMap mergedMap_, int minValidTaxa_, boolean simplify_,
			boolean reassign_, boolean skipNorank_, int inferRankLimit_, int threads){
		this(nodes_, mergedMap_, minValidTaxa_, simplify_, reassign_, skipNorank_, inferRankLimit_, threads, null, nodes_.length);
	}
	
	/** Adopts prepared topology, including optional packed auxiliary nodes.
	 * Every non-null node must resolve back to its array position through nodeIndex.
	 * Derived child/rank bounds are supplied by the caller; only rank lists/counts
	 * are rebuilt here. Mutating retained nodes afterward can invalidate those indexes.
	 * @param nodes_ Sparse native nodes followed by any packed auxiliary nodes
	 * @param mergedMap_ Old-to-current native TaxID map, or null
	 * @param minValidTaxa_ Legacy minimum-count metadata
	 * @param simplify_ Recorded simplification setting
	 * @param reassign_ Recorded no-rank reassignment setting
	 * @param skipNorank_ Recorded no-rank removal setting
	 * @param inferRankLimit_ Recorded rank-inference limit
	 * @param threads Requested indexing workers
	 * @param externalToInternalMap_ Auxiliary durable IDs to packed positions, or null
	 * @param firstAuxIndex_ First private packed position; native-only trees use nodes_.length */
	TaxTree(TaxNode[] nodes_, IntHashMap mergedMap_, int minValidTaxa_, boolean simplify_,
			boolean reassign_, boolean skipNorank_, int inferRankLimit_, int threads,
			IntHashMap externalToInternalMap_, int firstAuxIndex_){
		nodes=nodes_;
		mergedMap=mergedMap_;
		externalToInternalMap=externalToInternalMap_;
		firstAuxIndex=firstAuxIndex_;
		minValidTaxa=minValidTaxa_;
		simplify=simplify_;
		reassign=reassign_;
		skipNorank=skipNorank_;
		inferRankLimit=inferRankLimit_;
		nodeCount=populateLevelData(threads);
	}

	/** Counts ranks and fills stable node lists, serially below 100,000 array slots.
	 * Called once on zeroed index counts. Worker failures are rethrown after joining;
	 * disjoint output ranges preserve increasing array-slot order within each rank.
	 * @param threads Requested worker count, clamped to at least one
	 * @return Number of non-null nodes */
	private int populateLevelData(int threads){
		final int threadCount=Tools.max(1, threads);
		if(threadCount<2 || nodes.length<100000){return populateLevelData();}

		ArrayList<LevelCountThread> counters=new ArrayList<LevelCountThread>(threadCount);
		for(int tid=0; tid<threadCount; tid++){
			int from=(int)(((long)nodes.length*tid)/threadCount);
			int to=(int)(((long)nodes.length*(tid+1))/threadCount);
			counters.add(new LevelCountThread(nodes, from, to, this));
		}
		ThreadWaiter.startAndWait(counters);

		int count=0;
		for(LevelCountThread ct : counters){
			rethrow(ct.failure);
			if(!ct.success){throw new RuntimeException("TaxTree level-count thread terminated abnormally.");}
			count+=ct.count;
			for(int level=0; level<nodesPerLevel.length; level++){
				nodesPerLevel[level]+=ct.counts[level];
			}
			for(int level=0; level<nodesPerLevelExtended.length; level++){
				nodesPerLevelExtended[level]+=ct.countsExtended[level];
			}
		}
		for(int level=0; level<nodesPerLevelExtended.length; level++){
			treeLevelsExtended[level]=new TaxNode[nodesPerLevelExtended[level]];
		}

		int[] offsets=new int[nodesPerLevelExtended.length];
		ArrayList<LevelFillThread> fillers=new ArrayList<LevelFillThread>(threadCount);
		for(LevelCountThread ct : counters){
			fillers.add(new LevelFillThread(nodes, ct.from, ct.to, treeLevelsExtended, offsets.clone()));
			for(int level=0; level<offsets.length; level++){
				offsets[level]+=ct.countsExtended[level];
			}
		}
		ThreadWaiter.startAndWait(fillers);
		for(LevelFillThread ft : fillers){
			rethrow(ft.failure);
			if(!ft.success){throw new RuntimeException("TaxTree level-fill thread terminated abnormally.");}
		}
		return count;
	}

	/** Rethrows a worker failure, preserving errors and runtime exceptions.
	 * @param t Failure to propagate; null does nothing */
	private static void rethrow(Throwable t){
		if(t instanceof Error){throw (Error)t;}
		if(t instanceof RuntimeException){throw (RuntimeException)t;}
		if(t!=null){throw new RuntimeException(t);}
	}

	/** Validates node positions/ranks and counts one disjoint array range. */
	private static final class LevelCountThread extends Thread{
		/** Assigns a read-only node range and allocates per-worker rank counters.
		 * @param nodes_ Shared prepared nodes
		 * @param from_ Inclusive range start
		 * @param to_ Exclusive range end
		 * @param tree_ Tree resolving durable IDs to array positions */
		LevelCountThread(TaxNode[] nodes_, int from_, int to_, TaxTree tree_){
			nodes=nodes_;
			from=from_;
			to=to_;
			tree=tree_;
			counts=new int[taxLevelNames.length];
			countsExtended=new int[taxLevelNamesExtended.length];
		}

		@Override
		public void run(){
			try{
				for(int i=from; i<to; i++){
					TaxNode n=nodes[i];
					if(n!=null){
						if(tree.nodeIndex(n.id)!=i){throw new IllegalArgumentException("TaxNode "+n.id+" is stored at index "+i);}
						if(n.level<0 || n.level>=counts.length){
							throw new IllegalArgumentException("Invalid taxonomic level "+n.level+" for TaxID "+n.id);
						}
						if(n.levelExtended<0 || n.levelExtended>=countsExtended.length){
							throw new IllegalArgumentException("Invalid extended taxonomic level "+n.levelExtended+" for TaxID "+n.id);
						}
						counts[n.level]++;
						countsExtended[n.levelExtended]++;
						count++;
					}
				}
				success=true;
			}catch(Throwable t){failure=t;}
		}

		/** Shared read-only node array. */
		final TaxNode[] nodes;
		/** Owning tree used for packed-ID resolution. */
		final TaxTree tree;
		/** Inclusive start and exclusive end of this worker's range. */
		final int from, to;
		/** Private counts indexed by normal and extended rank. */
		final int[] counts, countsExtended;
		/** Number of non-null nodes observed. */
		int count=0;
		/** Captured failure, read after joining the worker. */
		Throwable failure=null;
		/** True only after the entire assigned range was counted. */
		boolean success=false;
	}

	/** Fills preassigned non-overlapping rank-list ranges from one input range. */
	private static final class LevelFillThread extends Thread{
		/** Adopts private offsets and shared output lists for this worker.
		 * @param nodes_ Validated node array
		 * @param from_ Inclusive input start
		 * @param to_ Exclusive input end
		 * @param levels_ Shared per-rank output arrays
		 * @param offsets_ Private mutable write cursors, one per extended rank */
		LevelFillThread(TaxNode[] nodes_, int from_, int to_, TaxNode[][] levels_, int[] offsets_){
			nodes=nodes_;
			from=from_;
			to=to_;
			levels=levels_;
			offsets=offsets_;
		}

		@Override
		public void run(){
			try{
				for(int i=from; i<to; i++){
					TaxNode n=nodes[i];
					if(n!=null){levels[n.levelExtended][offsets[n.levelExtended]++]=n;}
				}
				success=true;
			}catch(Throwable t){failure=t;}
		}

		/** Shared read-only node array. */
		final TaxNode[] nodes;
		/** Inclusive start and exclusive end of the input range. */
		final int from, to;
		/** Shared output arrays, with this worker's ranges disjoint from all others. */
		final TaxNode[][] levels;
		/** Private mutable write position for each extended rank. */
		final int[] offsets;
		/** Captured failure, read after joining. */
		Throwable failure=null;
		/** True only after the assigned range was filled. */
		boolean success=false;
	}

	/** Serial rank-index construction from prepared nodes and zeroed rank counters.
	 * Validates positions and rank bounds, but does not validate parent topology.
	 * @return Number of non-null nodes */
	private int populateLevelData(){
		int count=0;
		for(int i=0; i<nodes.length; i++){
			TaxNode n=nodes[i];
			if(n!=null){
				if(nodeIndex(n.id)!=i){throw new IllegalArgumentException("TaxNode "+n.id+" is stored at index "+i);}
				if(n.level<0 || n.level>=nodesPerLevel.length){
					throw new IllegalArgumentException("Invalid taxonomic level "+n.level+" for TaxID "+n.id);
				}
				if(n.levelExtended<0 || n.levelExtended>=nodesPerLevelExtended.length){
					throw new IllegalArgumentException("Invalid extended taxonomic level "+n.levelExtended+" for TaxID "+n.id);
				}
				nodesPerLevel[n.level]++;
				nodesPerLevelExtended[n.levelExtended]++;
				count++;
			}
		}
		for(int i=0; i<nodesPerLevelExtended.length; i++){
			treeLevelsExtended[i]=new TaxNode[nodesPerLevelExtended[i]];
		}
		int[] temp=new int[nodesPerLevelExtended.length];
		for(TaxNode n : nodes){
			if(n!=null){
				int level=n.levelExtended;
				treeLevelsExtended[level][temp[level]++]=n;
			}
		}
		return count;
	}

	/** Detects auxiliary mappings or optional runtime indexes omitted by TaxTree2.
	 * Node abundance counters are checked separately by the text writer.
	 * @return Whether any unsupported map is present, even if empty */
	final boolean hasTextExcludedState(){
		return externalToInternalMap!=null || nameMap!=null || nameMapLower!=null || childMap!=null || refseqSizeMap!=null
				|| refseqSizeMapC!=null || refseqSeqMap!=null || refseqSeqMapC!=null || nodeMapC!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*---------------            Loaders            ----------------*/
	/*--------------------------------------------------------------*/


	/**
	 * Load a tax tree from disk.
	 * Uses the process-wide treeFile setting; null disables loading.
	 * @param outstream Optional status destination
	 * @param hashNames Hash nodes using names as keys
	 * @param hashDotFormat Hash nodes using abbreviations, e.g. H.sapiens
	 * @return Loaded/shared tree, or null when treeFile is null
	 */
	public static final TaxTree loadTaxTree(PrintStream outstream, boolean hashNames, 
			boolean hashDotFormat){
		if(treeFile==null){return null;}
		return loadTaxTree(treeFile, null, null, null, outstream, hashNames, hashDotFormat);
	}

	/**
	 * Load a tax tree from disk.
	 * @param treeFile Portable TSV, legacy serialized tree, or auto; null disables loading
	 * @param outstream Optional status destination
	 * @param hashNames Hash nodes using names as keys
	 * @param hashDotFormat Hash nodes using abbreviations, e.g. H.sapiens
	 * @return Loaded/shared tree, or null for a null treeFile
	 */
	public static final TaxTree loadTaxTree(String treeFile, PrintStream outstream, boolean hashNames, 
			boolean hashDotFormat){
		if(treeFile==null){return null;}
		return loadTaxTree(treeFile, null, null, null, outstream, hashNames, hashDotFormat);
	}
	
	/**
	 * Load a tax tree from disk, either from a binary tree file,
	 * or from NCBI text files.
	 * When sharing is enabled, an existing cached tree takes precedence over a new
	 * filename. Direct dump construction uses names and nodes only when no tree file
	 * is selected. Name hashing mutates the loaded tree before it is returned.
	 * @param taxTreeFile Portable TSV or serialized representation, taking precedence over dumps
	 * @param taxNameFile NCBI names.txt
	 * @param taxNodeFile NCBI nodes.txt
	 * @param taxMergedFile NCBI merged.txt
	 * @param outstream Optional status destination
	 * @param hashNames Hash nodes using names as keys
	 * @param hashDotFormat Hash nodes using abbreviations, e.g. H.sapiens
	 * @return The loaded tree
	 */
	public static final TaxTree loadTaxTree(String taxTreeFile, String taxNameFile, String taxNodeFile, 
			String taxMergedFile, PrintStream outstream, boolean hashNames, boolean hashDotFormat){
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=defaultTreeFile();}
		if(taxTreeFile!=null || taxNodeFile==null){
			TaxTree tree=sharedTree(taxTreeFile, hashNames, hashDotFormat, outstream);
			if(tree!=null){return tree;}
		}
		assert(taxTreeFile!=null || (taxNameFile!=null && taxNodeFile!=null)) : 
			"Must specify both taxname and taxnode files.";
		Timer t=new Timer();
		if(outstream!=null){outstream.print("\nLoading tax tree; ");}
		final TaxTree tree;
		if(taxTreeFile!=null){
			tree=loadTaxTreeFile(taxTreeFile);
		}else{
			tree=new TaxTree(taxNameFile, taxNodeFile, taxMergedFile, null);
		}
		t.stop();
		if(hashNames){
			if(outstream!=null){outstream.println("Hashing names.");}
			tree.hashNames(hashDotFormat);
		}
		if(outstream!=null){
			outstream.println("Time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		if(ALLOW_SHARED_TREE){sharedTree=tree;}
		return tree;
	}

	/** Loads portable TSV or legacy serialization without using the shared-tree cache.
	 * Applies the configured auxiliary overlay when LOAD_FAKE_TREE is enabled.
	 * @param fname Tree path, optionally followed by a comma and an auxiliary path
	 * @return Loaded tree, with any requested overlay */
	public static TaxTree loadTaxTreeFile(String fname){
		return loadTaxTreeFile(fname, LOAD_FAKE_TREE ? fakeTreeFile : null);
	}
	
	/** Loads one base and applies an inline overlay before the optional global overlay.
	 * @param fname Base path, optionally followed by a comma and an auxiliary path
	 * @param auxFile Additional auxiliary path, or null
	 * @return Loaded tree after overlays */
	private static TaxTree loadTaxTreeFile(String fname, String auxFile){
		int comma=fname==null ? -1 : fname.indexOf(',');
		if(comma>=0){
			TaxTree tree=loadTaxTreeFile(fname.substring(0, comma), null);
			tree=loadAuxTree(tree, fname.substring(comma+1));
			return auxFile==null ? tree : loadAuxTree(tree, auxFile);
		}
		final TaxTree tree;
		if(isTextTreeFile(fname)){
			tree=TaxTreeText.load(fname);
		}else{
			tree=ReadWrite.read(TaxTree.class, fname, true);
		}
		return auxFile==null ? tree : loadAuxTree(tree, auxFile);
	}
	
	/** Overlays an auxiliary file using BBTools resource-path semantics.
	 * @param tree Base taxonomy
	 * @param fname Auxiliary resource or filesystem path
	 * @return Overlay tree, or the base tree when the auxiliary table has no rows */
	private static TaxTree loadAuxTree(TaxTree tree, String fname){
		String path=Data.findPath(fname);
		if(path==null){throw new RuntimeException("Could not find auxiliary TaxTree file "+fname);}
		return TaxTreeAux.apply(tree, path);
	}

	/** Writes portable TSV or legacy serialization according to the raw extension.
	 * TSV serialization rejects runtime state that its format cannot preserve.
	 * @param tree Tree to write
	 * @param fname Destination path; a raw tsv extension selects portable text
	 * @param overwrite Permit replacement of an existing destination */
	public static void writeTaxTree(TaxTree tree, String fname, boolean overwrite){
		if(isTextTreeFile(fname)){TaxTreeText.write(tree, fname, overwrite);}
		else{ReadWrite.write(tree, fname, overwrite);}
	}

	/** Tests whether the extension after removing compression suffixes is tsv.
	 * @param fname Path to classify, or null
	 * @return Whether the path selects portable taxonomy text */
	public static boolean isTextTreeFile(String fname){
		return fname!=null && "tsv".equals(ReadWrite.rawExtension(fname));
	}

	/**
	 * Parse a component-local tree flag while preserving legacy tree=path syntax.
	 * Boolean values only control the caller; a path selects the shared source and enables the caller.
	 * @param value Boolean token, path, null or empty
	 * @return Whether the caller should enable its tree-dependent operation
	 */
	public static boolean parseTreeFlag(String value){
		if(value==null || value.isEmpty()){return true;}
		if(Parse.isBoolean(value) || "1".equals(value) || "0".equals(value)
				|| "null".equalsIgnoreCase(value) || "none".equalsIgnoreCase(value)){
			return Parse.parseBoolean(value);
		}
		treeFile=value;
		return true;
	}
	
	/** Parses the global auxiliary fake-TaxID tree flag, updating its enable flag/path.
	 * @param value Boolean token, explicit path, null or empty for the default enabled path
	 * @return Updated LOAD_FAKE_TREE value */
	public static boolean parseFakeTreeFlag(String value){
		if(value==null || value.isEmpty()){
			LOAD_FAKE_TREE=true;
			fakeTreeFile=defaultFakeTreeFile;
		}else if(Parse.isBoolean(value) || "1".equals(value) || "0".equals(value)
				|| "null".equalsIgnoreCase(value) || "none".equalsIgnoreCase(value)){
			LOAD_FAKE_TREE=Parse.parseBoolean(value);
			if(LOAD_FAKE_TREE){fakeTreeFile=defaultFakeTreeFile;}
		}else{
			LOAD_FAKE_TREE=true;
			fakeTreeFile=value;
		}
		return LOAD_FAKE_TREE;
	}
	
	/*--------------------------------------------------------------*/
	/*---------------      Constructor Helpers      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Finds unranked nodes in the archaeal and bacterial kingdoms.
	 * If these are below species level, have a ranked parent,
	 * and have no ranked children, they are assigned strain or substrain.
	 */
	private void assignStrains(){

		outstream.println("Assigning strains.");
		int strains=0, substrains=0;
		TaxNode bacteria=getNode(BACTERIA_ID); //Can't do a name lookup since the names are not hashed yet
		TaxNode archaea=getNode(ARCHAEA_ID);
		assert(bacteria.name.equalsIgnoreCase("Bacteria"));
		assert(archaea.name.equalsIgnoreCase("Archaea"));

		ArrayList<TaxNode> bactList=new ArrayList<TaxNode>();
		ArrayList<TaxNode> archList=new ArrayList<TaxNode>();
		for(TaxNode tn : nodes){
			if(tn!=null && tn.originalLevel()==NO_RANK && tn.minParentLevelExtended<=SPECIES_E){
				if(descendsFrom(tn, bacteria)){
					bactList.add(tn);
				}else if(descendsFrom(tn, archaea)){
					archList.add(tn);
				}
			}
		}
		
		ArrayList<TaxNode> prokList=new ArrayList<TaxNode>(bactList.size()+archList.size());
		prokList.addAll(bactList);
		prokList.addAll(archList);
		
		for(TaxNode tn : prokList){
			if(tn.maxDescendantLevelIncludingSelf()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==SPECIES_E || parent.levelExtended==SUBSPECIES_E){
					tn.levelExtended=STRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(STRAIN_E);
					strains++;
				}
			}
		}
		
//		outstream.println("Assigned "+strains+" strains.");
		for(TaxNode tn : prokList){
			if(tn.maxDescendantLevelIncludingSelf()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==STRAIN_E){
					tn.levelExtended=SUBSTRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(SUBSTRAIN_E);
					substrains++;
				}
			}
		}
//		outstream.println("Assigned "+substrains+" substrains.");
	}
	
	/** @deprecated Legacy strain assignment method */
	@Deprecated
	private void assignStrainsOld(){

		outstream.println("Assigning strains.");
		int strains=0, substrains=0;
		TaxNode bacteria=getNode(BACTERIA_ID); //Can't do a name lookup since the names are not hashed
		assert(bacteria.name.equalsIgnoreCase("Bacteria"));
		for(TaxNode tn : nodes){
			if(tn!=null && tn.originalLevel()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==SPECIES_E && commonAncestor(parent, bacteria)==bacteria){
//					nodesPerLevelExtended[STRAIN_E]++;
//					nodesPerLevelExtended[tn.levelExtended]--;
					tn.levelExtended=STRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(STRAIN_E);
					strains++;
				}
			}
		}
//		outstream.println("Assigned "+strains+" strains.");
		for(TaxNode tn : nodes){
			if(tn!=null && tn.originalLevel()==NO_RANK){
				TaxNode parent=nodes[tn.pid];
				if(parent.levelExtended==STRAIN_E && commonAncestor(parent, bacteria)==bacteria){
//					nodesPerLevelExtended[SUBSTRAIN_E]++;
//					nodesPerLevelExtended[tn.levelExtended]--;
					tn.levelExtended=SUBSTRAIN_E;
					tn.level=SUBSPECIES;
					tn.setOriginalLevel(SUBSTRAIN_E);
					substrains++;
				}
			}
		}
//		outstream.println("Assigned "+substrains+" substrains.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Create tax nodes using names in the designated file.
	 * @param fname NCBI names.txt
	 * @return Array of created nodes, where array[x] contains the node with TaxID x.
	 */
	private static TaxNode[] getNames(String fname){
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(200000);
		int max=0;
		
		//TODO: Probable bug - this owned TextFile is never closed, omitting cleanup and close-status validation.
		TextFile tf=new TextFile(fname, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.contains("scientific name")){
				String[] split=delimiter.split(s, 3);
				assert(split.length==3) : s;
				int id=Integer.parseInt(split[0]);
				String name=split[1];
				if(id==1 && name.equalsIgnoreCase("root")){name="Life";}
				max=Tools.max(max, id);
				list.add(new TaxNode(id, name));
			}
		}
		
		TaxNode[] nodes=new TaxNode[max+1];
		for(TaxNode n : list){
			assert(nodes[n.id]==null || nodes[n.id].equals(n)) : nodes[n.id]+" -> "+n;
			nodes[n.id]=n;
		}
		
		return nodes;
	}
	
	/**
	 * Adds parent IDs, rank values and canonical flags from NCBI node records.
	 * TODO: Consider combining the two dump-loading passes where their inputs permit it.
	 * @param fname NCBI nodes.dmp file, not the names file
	 * @param nodes Mutable array previously populated from scientific names
	 * @return The supplied array after node fields are filled
	 */
	private static TaxNode[] getNodes(String fname, TaxNode[] nodes){
		
		int max=0;
		
		LinkedHashMap<String, int[]> oddNames=new LinkedHashMap<String, int[]>();
		
		//TODO: Probable bug - this owned TextFile is never closed, omitting cleanup and close-status validation.
		TextFile tf=new TextFile(fname, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			String[] split=delimiter.split(s, 4);
			assert(split.length==4) : s;
			int id=-1, pid=-1, level=-1, levelExtended=-1;
			
			id=Integer.parseInt(split[0]);
			try{
				pid=Integer.parseInt(split[1]);
			}catch(NumberFormatException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println("Bad line: "+s+"\n"+Arrays.toString(split));
			}
			boolean alt=false;
			{
				String key=split[2];
				Integer obj0=levelMap.get(key);
				Integer obj=levelMapExtended.get(key);
				if(obj==null){
					if(key.equalsIgnoreCase("cellular root") || key.equalsIgnoreCase("acellular root")){
						obj=DOMAIN_E;
						obj0=DOMAIN;
					}
					assert(obj!=null) : "No level found for "+key+"; line="+Arrays.toString(split);
				}
				
				if(obj0==null){
					obj0=altLevelMap.get(key);
					alt=true;
				}
				if(obj0!=null){
					level=obj0;
					levelExtended=obj;
					if(id==pid){
						level=LIFE;
						levelExtended=LIFE_E;
						alt=false;
					}
				}else{
					if(id==pid){
						level=LIFE;
						levelExtended=LIFE_E;
						alt=false;
					}else{
						int[] count=oddNames.get(key);
						if(count==null){
							count=new int[1];
							oddNames.put(key, count);
						}
						count[0]++;
					}
				}
			}
			max=Tools.max(max, id);
			TaxNode n=nodes[id];
			assert(n!=null && n.pid<0) : n+" -> "+s;
			n.pid=pid;
			n.level=level;
			n.levelExtended=levelExtended;
			n.setOriginalLevel(levelExtended);
			n.setCanonical(!alt);
			assert(n.canonical()==n.isSimple() || n.levelExtended==NO_RANK_E) : n.canonical()+", "+n.isSimple()+", "+n.level+", "+n.levelExtended+"\n"+n.toString()+"\n";
		}
		
		if(oddNames.size()>0){
			outstream.println("Found "+oddNames.size()+" unknown taxonomic levels:");
			if(verbose){
				for(String s : oddNames.keySet()){
					outstream.println(oddNames.get(s)[0]+"\t"+s);
				}
			}
		}
		
		return nodes;
	}
	
	/**
	 * Count child nodes of each node.
	 * This can be used to size arrays or determine which nodes are leaves.
	 */
	private void countChildren(){
		for(TaxNode child : nodes){
			if(child!=null && child.pid!=child.id){
				TaxNode parent=getNode(child.pid);
				if(parent!=child){parent.numChildren++;}
			}
		}
	}
	
	//TODO - This could be finished in 2 passes using the childTable.
	/**
	 * Fill derived fields minParentLevelExtended and maxChildLevelExtended
	 * by percolating information through the tree until a fixpoint is reached.
	 * TODO: This could be finished in two passes using a child table.
	 * @return Number of rounds required to reach fixpoint.
	 */
	private int percolate(){
		boolean changed=true;
		int rounds=0;
		while(changed){
			changed=false;
			rounds++;
			for(TaxNode child : nodes){
				if(child!=null && child.pid!=child.id){
					TaxNode parent=getNode(child.pid);
					changed=(child.discussWithParent(parent) | changed);
				}
			}
			
			if(!changed){break;}
			changed=false;
			rounds++;
			for(int i=nodes.length-1; i>=0; i--){
				TaxNode child=nodes[i];
				if(child!=null && child.pid!=child.id){
					TaxNode parent=getNode(child.pid);
					changed=(child.discussWithParent(parent) | changed);
				}
			}
		}
		return rounds;
	}
	
	/**
	 * Load nodes into the nameMap and nameMapLower, mapped to their names.
	 * @param genusDotSpecies Also hash abbreviations such as E.coli.
	 */
	public synchronized void hashNames(boolean genusDotSpecies){
		if(nameMap!=null){return;}
		assert(nameMap==null);
		assert(nameMapLower==null);
		final int size=((int)Tools.mid(2, (nodes.length+(genusDotSpecies ? nodesPerLevelExtended[SPECIES_E] : 0))*1.5, Shared.MAX_ARRAY_LEN));
		nameMap=new HashMap<String, ArrayList<TaxNode>>(size);
		nameMapLower=new HashMap<String, ArrayList<TaxNode>>(size);
		
		//Hash the names, both lowercase and uppercase
		for(TaxNode n : nodes){
			if(n!=null){
				String name=n.name;
				if(name.indexOf('_')>=0){
					name=name.replace('_', ' ').trim();
				}
				if(name!=null && !name.equals("environmental samples")){
					{
						ArrayList<TaxNode> list=nameMap.get(name);
						if(list==null){
							list=new ArrayList<TaxNode>();
							nameMap.put(name, list);
						}
						list.add(n);
					}
					{
						String lc=name.toLowerCase();
						ArrayList<TaxNode> list=nameMapLower.get(lc);
						if(list==null){
							list=new ArrayList<TaxNode>();
							nameMapLower.put(lc, list);
						}
						list.add(n);
					}
				}
			}
		}
		
		//Hash G.species versions of the names, both lowercase and uppercase
		if(genusDotSpecies){
			ByteBuilder bb=new ByteBuilder(64);
			for(TaxNode n : nodes){
				if(n!=null && n.levelExtended==SPECIES_E){
					String name=n.name;
					if(name.indexOf('_')>=0){
						name=name.replace('_', ' ').trim();
					}
					if(name!=null && !name.equals("environmental samples")){
						final String dotFormat=dotFormat(name, bb);
						if(dotFormat!=null){
							{
								ArrayList<TaxNode> list=nameMap.get(dotFormat);
								if(list==null){
									list=new ArrayList<TaxNode>();
									nameMap.put(dotFormat, list);
								}
								list.add(n);
							}
							{
								String lc=dotFormat.toLowerCase();
								ArrayList<TaxNode> list=nameMapLower.get(lc);
								if(list==null){
									list=new ArrayList<TaxNode>();
									nameMapLower.put(lc, list);
								}
								list.add(n);
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * Generate the "dot format" name of a node.
	 * For example, transform "Homo sapiens" to "H.sapiens"
	 * @param name Node name
	 * @param buffer A ByteBuilder that may be modified
	 * @return Dot format
	 */
	private static String dotFormat(String name, ByteBuilder buffer){
		if(name==null || name.indexOf('.')>=0){return null;}
		final int firstSpace=name.indexOf(' ');
		if(firstSpace<0 || firstSpace>=name.length()-1){return null;}
		final int lastSpace=name.lastIndexOf(' ');
		if(firstSpace!=lastSpace){return null;}
		final String a=name.substring(0, firstSpace);
		final String b=name.substring(lastSpace+1);
		final char ca=a.charAt(0);
		final char cb=b.charAt(0);
		if(!Tools.isUpperCase(ca) || !Tools.isLowerCase(cb)){return null;}
		if(buffer==null){buffer=new ByteBuilder(2+b.length());}
		else{buffer.clear();}
		buffer.append(ca).append('.').append(b);
		return buffer.toString();
	}
	
	/**
	 * Fill childMap, which maps nodes to their children.
	 */
	public synchronized void hashChildren(){
		assert(childMap==null);
		int nodesWithChildren=0;
		for(TaxNode tn : nodes){
			if(tn!=null && tn.numChildren>0){nodesWithChildren++;}
		}
		childMap=new HashMap<TaxNode, ArrayList<TaxNode>>((int)Tools.mid(2, nodesWithChildren*1.5, Shared.MAX_ARRAY_LEN));
		for(TaxNode tn : nodes){
			if(tn!=null){
				if(tn.numChildren>0){
					childMap.put(tn, new ArrayList<TaxNode>(tn.numChildren));
				}
			}
		}
		for(TaxNode tn : nodes){
			if(tn!=null){
				if(tn.id!=tn.pid){
					ArrayList<TaxNode> list=childMap.get(getNode(tn.pid));
					if(list!=null){list.add(tn);}
				}
			}
		}
	}
	
	/**
	 * Fetches this node's children, borrowing the cached list when available.
	 * @param parent Non-null node whose children are requested
	 * @return Child list, or null when the recorded child count is nonpositive or the cache
	 * has no entry; without a cache, a new list is collected from the node array
	 */
	public ArrayList<TaxNode> getChildren(TaxNode parent){
		if(parent.numChildren<1){return null;}
		if(childMap!=null){return childMap.get(parent);}
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(parent.numChildren);
		for(TaxNode tn : nodes){
			if(tn!=null && tn.id!=tn.pid && tn.pid==parent.id){
				list.add(tn);
			}
		}
		return list;
	}
	
	/**
	 * Load a map of old to new TaxIDs.
	 * @param mergedFile NCBI merged.txt.
	 * @return Map of old to new TaxIDs
	 */
	private static IntHashMap getMerged(String mergedFile){
		if(mergedFile==null){return null;}
		String[] lines=TextFile.toStringLines(mergedFile);
		if(lines.length<1){return null;}
		IntHashMap map=new IntHashMap((int)(lines.length*1.3));
		for(String line : lines){
			String[] split=delimiterTab.split(line);
			int a=Integer.parseInt(split[0]);
			int b=Integer.parseInt(split[2]);
			map.put(a, b);
		}
		return map;
	}
	
	/**
	 * Simplify the tree by assigning ranks to unranked nodes,
	 * where possible, through inference.
	 * Optionally removes unranked nodes based on the skipNorank field.
	 * @param nodes Array of all TaxNodes.
	 * @return Number of nodes removed.
	 */
	private int simplify(TaxNode nodes[]){
		
		int failed=test(nodes);
		
		int removed=0;
		int reassigned=0;
		
		if(reassign){
			boolean changed=true;
			int changedCount=0;
			while(changed){
				changed=false;
				for(int i=0; i<nodes.length; i++){
					TaxNode n=nodes[i];
					if(n!=null && n.levelExtended<1){
						int pid=n.pid;
						TaxNode parent=nodes[pid];
						assert(parent!=null) : n;
						if(n.maxDescendantLevelIncludingSelf()<SUBSPECIES_E &&
								parent.minAncestorLevelIncludingSelf()<=SPECIES_E && parent.minAncestorLevelIncludingSelf()>=SUBSPECIES_E){
//						if(parent.levelExtended==SPECIES_E || parent.levelExtended==SUBSPECIES_E){
							changed=true;
							n.levelExtended=SUBSPECIES_E;
							n.level=SUBSPECIES;
							changedCount++;
						}
					}
				}
			}
			System.err.println("Assigned levels to "+changedCount+" unranked nodes.");
		}
		
		
		if(skipNorank){//Skip nodes with unknown taxa
			if(verbose){outstream.println("A0");}
			
			for(int i=0; i<nodes.length; i++){
				TaxNode n=nodes[i];
				if(n!=null){
					int pid=n.pid;
					TaxNode parent=nodes[pid];
					assert(parent!=null) : n;
					assert(parent!=n || pid==1) : n+", "+pid;
					while(parent.levelExtended<1 && n.levelExtended>parent.levelExtended){
						//System.err.println("Reassigned from "+parent);
						assert(parent.id!=parent.pid);
						parent=nodes[parent.pid];
						n.pid=parent.id;
						reassigned++;
					}
				}
			}
			
			for(int i=0; i<nodes.length; i++){
				if(nodes[i]!=null && nodes[i].levelExtended<0){
					System.err.println("Removed "+nodes[i]);
					nodes[i]=null;
					removed++;
				}
			}
			if(verbose){outstream.println("Skipped "+reassigned+" unranked parents, removed "+removed+" invalid nodes.");}
		}
		
		if(inferRankLimit>0){//Infer level for unset nodes (from "no rank")
			if(verbose){outstream.println("A");}
			int changed=1;
			while(changed>0){
				changed=0;
				for(final TaxNode n : nodes){
					if(n!=null){
						if(n.levelExtended==0){
							TaxNode parent=nodes[n.pid];
							if(n!=parent && parent.levelExtended>0 && parent.levelExtended<=inferRankLimit+1){
								n.levelExtended=Tools.max(1, parent.levelExtended-1);
								assert(n.levelExtended>0 && n.levelExtended<=parent.levelExtended && n.levelExtended<=inferRankLimit);
								changed++;
							}
						}
					}
				}
				if(verbose){outstream.println("changed: "+changed);}
			}
			
//			outstream.println("B");
//			for(TaxNode n : nodes){
//				if(n!=null && n.level==0){
//					n.level=-1;
//				}
//			}
		}
		
		failed=test(nodes);
		
//		if(reassign){//Skip nodes with duplicate taxa
//			if(verbose){outstream.println("D");}
//			int changed=1;
//			while(changed>0){
//				changed=0;
//				for(final TaxNode n : nodes){
//					if(n!=null){
//						TaxNode parent=nodes[n.pid];
//						TaxNode grandparent=nodes[parent.pid];
//						assert(n.level<=parent.level || parent.level<1 || !parent.canonical()) : n+" -> "+parent+" -> "+grandparent;
//						assert(parent.level<=grandparent.level || grandparent.level<1 || !grandparent.canonical()) : n+" -> "+parent+" -> "+grandparent;
//
//						while(parent!=grandparent && (parent.level<0 || (parent.level==grandparent.level && !parent.canonical()) ||
//								n.level>parent.level || (n.level==parent.level))){
//							parent=grandparent;
//							grandparent=nodes[parent.pid];
//							n.pid=parent.id;
//							reassigned++;
//							changed++;
//						}
//					}
//				}
//				if(verbose){outstream.println("changed: "+changed);}
//			}
//			if(verbose){outstream.println("E");}
//			for(int i=0; i<nodes.length; i++){
//				if(nodes[i]!=null && nodes[i].level<0){
//					nodes[i]=null;
//					removed++;
//				}
//			}
//		}
		
		failed=test(nodes);

		if(verbose){outstream.println("F");}
		{//Ensure the tree is now clean
			for(int i=0; i<nodes.length; i++){
				TaxNode n=nodes[i];
				if(n!=null){
					TaxNode parent=nodes[n.pid];
					TaxNode grandparent=nodes[parent.pid];
					assert(n==parent || n.levelExtended<=parent.levelExtended || !n.canonical() || n.levelExtended<1 || parent.levelExtended<1) : n+" -> "+parent+" -> "+grandparent;
					assert(parent==grandparent || parent.levelExtended<=grandparent.levelExtended || !parent.canonical() || parent.levelExtended<1 || grandparent.levelExtended<1) : n+" -> "+parent+" -> "+grandparent;
				}
			}
		}
		
//		if(verbose){System.err.println("Reassignments: "+reassigned);}
		
		return removed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Ensure tree has monotonically increasing (or nondescending) ranks.
	 * @param nodes All TaxNodes.
	 * @return Number of violations.
	 */
	private static int test(TaxNode[] nodes){
		int failed=0;
		for(final TaxNode n : nodes){
			if(n!=null){
				TaxNode parent=nodes[n.pid];
				try{
					assert(n==parent || n.level<=parent.level || parent.level<1 || !parent.canonical()) : 
						"\n"+n+" -> "+parent+", level="+n.level+", plevel="+parent.level+", pcanon="+parent.canonical()+"\n"
								+ "levelE="+n.levelExtended+", plevelE="+parent.levelExtended;
					assert(n==parent || n.levelExtended<=parent.levelExtended || parent.levelExtended<1) : n+" -> "+parent;
//				assert(n==parent || n.level<parent.level || parent.level<1 || !n.canonical() || !parent.canonical()) : n+" -> "+parent;
					if(n!=parent && n.level>parent.level && parent.level>=1 && n.canonical() && parent.canonical()){
						if(verbose){outstream.println("Error: "+n+" -> "+parent);}
						failed++;
					}else if(n!=parent && parent.levelExtended>=1 && n.levelExtended>=parent.levelExtended){
//					if(verbose){outstream.println("Error: "+n+" -> "+parent);}
//					failed++;
					}
					assert(n!=parent || n.id<=1) : n;
				}catch(Throwable e){
					// TODO Auto-generated catch block
					e.printStackTrace();
					failed++;
				}
			}
		}
		if(verbose || failed>0){outstream.println(failed+" nodes failed.");}
		return failed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Print Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Format full name in semicolon format, e.g. 
	 * "SK:Bacteria;P:Protobacteria;..."
	 * @param tn0 Base node
	 * @param skipNonCanonical Ignore noncanonical (aka "nonsimple") levels like Tribe.
	 * @param mononomial Remove a matching parent-name prefix when rendering a node name
	 * @return Resultant String
	 */
	public String toSemicolon(final TaxNode tn0, boolean skipNonCanonical, boolean mononomial){
		StringBuilder sb=new StringBuilder();
		if(tn0==null){return "Not found";}
		String semi="";
		ArrayList<TaxNode> list=toAncestors(tn0, skipNonCanonical);
		boolean addTaxLevel=true;
		for(int i=list.size()-1; i>=0; i--){
			sb.append(semi);
			TaxNode tn=list.get(i);
			if(tn.id!=LIFE_ID || list.size()==1){
				if(addTaxLevel && tn.canonical() && !tn.levelChanged() && tn.isSimple()){
					sb.append(tn.levelToStringShort()).append(':');
				}
				sb.append(mononomial ? mononomial(tn) : tn.name);
				semi=";";
			}
		}
		return sb.toString();
	}
	
	/**
	 * Collects selected TaxIDs from the starting node toward the available root.
	 * Omits the cellular-organisms node when it is an ancestor, and optionally
	 * omits noncanonical ranks. If filtering removes every candidate, retains tn0.
	 * @param tn0 Starting node, or null
	 * @param skipNonCanonical Omit noncanonical ranks such as tribe
	 * @return New list in start-to-root order, or null for a null starting node
	 */
	public IntList toAncestorIds(final TaxNode tn0, boolean skipNonCanonical){
		if(tn0==null){return null;}
		IntList list=new IntList(8);
		
		TaxNode tn=tn0;
		while(tn!=null){
			if(!skipNonCanonical || tn.isSimple()){
				if(tn.id!=CELLULAR_ORGANISMS_ID || tn==tn0){list.add(tn.id);}
			}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		if(list.isEmpty()){list.add(tn0.id);}
		return list;
	}
	
	/**
	 * Collects selected nodes from the starting node toward the available root.
	 * Omits the cellular-organisms node when it is an ancestor, and optionally
	 * omits noncanonical ranks. If filtering removes every candidate, retains tn0.
	 * @param tn0 Starting node, or null
	 * @param skipNonCanonical Omit noncanonical ranks such as tribe
	 * @return New list of borrowed nodes in start-to-root order, or null for null input
	 */
	public ArrayList<TaxNode> toAncestors(final TaxNode tn0, boolean skipNonCanonical){
		if(tn0==null){return null;}
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(8);
		
		TaxNode tn=tn0;
		while(tn!=null){
			if(!skipNonCanonical || tn.isSimple()){
				if(tn.id!=CELLULAR_ORGANISMS_ID || tn==tn0){list.add(tn);}
			}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		if(list.isEmpty()){list.add(tn0);}
		return list;
	}
	
	/**
	 * Generate a path to the genome of an organism on the filesystem;
	 * used by ExplodeTree.  Intended for internal JGI use.
	 * @param node Starting node; its ancestry must resolve to TaxID1
	 * @param root Location of the exploded tree.
	 * @return Path to a genome.
	 */
	public String toDir(TaxNode node, String root){
		StringBuilder sb=new StringBuilder();
		if(root==null){root="";}
		sb.append(root);
		if(root.length()>0 && !root.endsWith("/")){sb.append('/');}
		IntList list=toAncestorIds(node, false);
		list.reverse();
		assert(list.get(0)==1) : list + "," +getNode(list.get(0));
		for(int i=0; i<list.size(); i++){
			sb.append(list.get(i));
			sb.append('/');
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Use various techniques to get a TaxID from an unknown String, such as parsing,
	 * name lookups, and accession translation.
	 * @param s String to process.
	 * @return Decoded TaxID.
	 */
	public static int getID(String s){return GiToTaxid.getID(s);}
	
	/**
	 * Use various techniques to get a TaxID from an unknown byte[], such as parsing,
	 * name lookups, and accession translation.
	 * @param s Byte array to process.
	 * @return Decoded TaxID.
	 */
	public static int getID(byte[] s){return GiToTaxid.getID(s);}
	
	/** Resolves a header and ascends until the current extended rank meets the bound.
	 * A root or missing ancestor also stops traversal.
	 * @param s Header/name accepted by parseNodeFromHeader with best-effort lookup
	 * @param minLevelExtended Minimum desired extended rank
	 * @return First qualifying ancestor, stopping root, or null if unresolved */
	public TaxNode getNode(String s, int minLevelExtended){
		TaxNode tn=parseNodeFromHeader(s, true);
		while(tn!=null && tn.levelExtended<minLevelExtended && tn.pid!=tn.id){
			tn=getNode(tn.pid);
		}
		return tn;
	}
	
	/**
	 * Determine whether a node is a descendant of another.
	 * @param child Possible child TaxID.
	 * @param parent Possible parent TaxID.
	 * @return true iff child descends from parent.
	 */
	public boolean descendsFrom(final int child, final int parent){
		TaxNode cn=getNode(child), pn=getNode(parent);
		assert(cn!=null) : "Invalid taxID: "+child;
		assert(pn!=null) : "Invalid taxID: "+parent;
		return descendsFrom(cn, pn);
	}
	
	/**
	 * Determine whether a node is a descendant of another.
	 * @param child Possible child node.
	 * @param parent Possible parent node.
	 * @return true iff child descends from parent.
	 */
	public boolean descendsFrom(TaxNode child, TaxNode parent){
		assert(child!=null && parent!=null) : "Null parameters.";
		if(child==null || parent==null){return false;}
		
		//FIXED [tax/TaxTree#001]: DROPPED the level early-exit (child.levelExtended<=parent.levelExtended). It assumed
		//levels are non-decreasing up the path, but the real NCBI tree is heavily non-monotonic: (a) no-rank nodes have
		//levelExtended==0 (a no-rank PARENT made any ranked child return a wrong FALSE, e.g. descendsFrom(Bacteria,
		//cellularOrganisms)), AND (b) genuine rank inversions are COMMON. Empirically, the minimal `parent.levelExtended<1`
		//guard still left 76,633 wrong-FALSE pairs in a 1.58M-pair sample (predict-then-run vs descendsFrom2), so no
		//level-based early-exit is safe here. Now a plain parent-chain walk (identical logic to descendsFrom2); depth is
		//~30 so the lost optimization is negligible. Post-fix: 0 mismatches vs descendsFrom2 across the whole tree.
		while(child!=parent && child.id!=child.pid){
			child=getNode(child.pid);
		}
		return child==parent;
	}
	
	/**
	 * Tests lineage membership using the legacy walk that excludes a self-parent root.
	 * Unlike descendsFrom(TaxNode, TaxNode), this never matches the root, even when
	 * the starting ID is that root. The starting node and its parent chain must exist.
	 * @param taxID Valid starting organism ID
	 * @param ancestorID Non-root ancestor ID to seek
	 * @return True if a non-root node on the lineage has ancestorID
	 */
	public boolean descendsFrom2(int taxID, final int ancestorID){
		TaxNode tn=getNode(taxID);
		//NOTE [tax/TaxTree]: the loop's root guard runs before its equality check,
		//so a self-parent root is never tested. descendsFrom2(x, LIFE_ID) is false.
		//The four isPlant/isAnimal/isFungus/isEukaryote wrappers seek non-root IDs.
		//For a root ancestor use descendsFrom(TaxNode,TaxNode), which includes it.
		while(tn.id!=tn.pid){
			if(tn.id==ancestorID){return true;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	/** Tests membership in the green-plant lineage.
	 * @param taxID Resolving ID with an intact parent chain
	 * @return True when the legacy lineage walk reaches VIRIDIPLANTAE_ID */
	public boolean isPlant(int taxID){return descendsFrom2(taxID, VIRIDIPLANTAE_ID);}
	
	/** Tests membership in the animal lineage.
	 * @param taxID Resolving ID with an intact parent chain
	 * @return True when the legacy lineage walk reaches METAZOA_ID */
	public boolean isAnimal(int taxID){return descendsFrom2(taxID, METAZOA_ID);}
	
	/** Tests membership in the fungal lineage.
	 * @param taxID Resolving ID with an intact parent chain
	 * @return True when the legacy lineage walk reaches FUNGI_ID */
	public boolean isFungus(int taxID){return descendsFrom2(taxID, FUNGI_ID);}
	
	/** Tests membership in the eukaryotic lineage.
	 * @param taxID Resolving ID with an intact parent chain
	 * @return True when the legacy lineage walk reaches EUKARYOTA_ID */
	public boolean isEukaryote(int taxID){return descendsFrom2(taxID, EUKARYOTA_ID);}
	
	/** Tests whether a lineage reaches Bacteria or Archaea before a self-parent root.
	 * An unresolved in-range ID prints a warning and returns false. getNode's range
	 * assertion still applies; resolved nodes must have intact parent chains.
	 * @param taxID Starting TaxID
	 * @return True when either prokaryotic domain is encountered */
	public boolean isProkaryote(int taxID){
		TaxNode tn=getNode(taxID);
		if(tn==null){
			System.err.println("*** Warning: Can't find node "+taxID+" ***");
			return false;
		}
		while(tn.id!=tn.pid){
			if(tn.id==BACTERIA_ID || tn.id==ARCHAEA_ID){return true;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	/**
	 * Returns the normal rank of the nearest shared ancestor.
	 * Positive IDs must resolve to nodes in this tree; assertions diagnose missing
	 * IDs or disconnected components. Parent chains must terminate at roots.
	 * @param a First taxID, including an auxiliary/merged ID supported by getNode
	 * @param b Second taxID
	 * @return Normal rank, or -1 for a nonpositive ID or no ancestor with assertions disabled
	 */
	public int commonAncestorLevel(final int a, final int b){
		if(a<1 || b<1){return -1;}
		TaxNode an=getNode(a), bn=getNode(b);
		assert(an!=null) : "Invalid taxID: "+a;
		assert(bn!=null) : "Invalid taxID: "+b;
		TaxNode cn=commonAncestor(an, bn);
		assert(cn!=null) : "No common ancestor: "+an+", "+bn;
		if(cn==null){return -1;}
		return cn.level;
	}
	
	/**
	 * Returns the taxID of the nearest shared ancestor of two valid IDs.
	 * ID resolution follows getNode, including its auxiliary and merged-ID maps.
	 * Assertions diagnose missing IDs and disconnected components.
	 * @param a First valid taxID
	 * @param b Second valid taxID
	 * @return Shared ancestor ID, or -1 for disconnected components with assertions disabled
	 */
	public int commonAncestor(final int a, final int b){
		TaxNode an=getNode(a), bn=getNode(b);
		assert(an!=null) : "Invalid taxID: "+a;
		assert(bn!=null) : "Invalid taxID: "+b;
		TaxNode cn=commonAncestor(an, bn);
		assert(cn!=null) : "No common ancestor: "+an+", "+bn;
		if(cn==null){return -1;}
		return cn.id;
	}
	
	/**
	 * Finds the nearest shared node on two parent chains.
	 * Strictly increasing ranks use the fast rank walk. Equal or inverted parent
	 * ranks trigger an exact depth-aligned fallback, so unranked ancestors are not
	 * skipped. This method allocates no per-call collection and does not mutate nodes.
	 * Both inputs must be canonical node objects from this tree, with finite parent
	 * chains ending at self-parent roots; topology must not change during a call.
	 * @param a First non-null node
	 * @param b Second non-null node
	 * @return Nearest shared node, or null for separate rooted components; with
	 * assertions disabled a null input returns the other input
	 */
	public TaxNode commonAncestor(TaxNode a, TaxNode b){
		assert(a!=null && b!=null) : "Null parameters.";
		if(a==null){return b;}
		if(b==null){return a;}
		final TaxNode originalA=a, originalB=b;
		while(a!=b){
			if(a.levelExtended<b.levelExtended){
				TaxNode parent=getNode(a.pid);
				if(parent.levelExtended<=a.levelExtended){return commonAncestorByDepth(originalA, originalB);}
				a=parent;
			}else{
				TaxNode parent=getNode(b.pid);
				if(parent.levelExtended<=b.levelExtended){return commonAncestorByDepth(originalA, originalB);}
				b=parent;
			}
		}
		return a;
	}

	/** Finds the exact parent-chain intersection without assuming ordered ranks.
	 * TAX-001: fallback must restart from the original inputs, since the rank walk
	 * may already have passed the shared ancestor before seeing a non-increasing edge.
	 * No per-call collections are allocated. Parent chains must terminate at roots.
	 * Stops early when either input is an ancestor of the other.
	 * @param a Original node from this tree
	 * @param b Other original node from this tree
	 * @return Nearest shared ancestor, or null for separate rooted components */
	private TaxNode commonAncestorByDepth(TaxNode a, TaxNode b){
		TaxNode rootA=a, rootB=b;
		int depthA=0, depthB=0;
		while(rootA!=b && rootA.pid!=rootA.id){rootA=getNode(rootA.pid); depthA++;}
		if(rootA==b){return b;}
		while(rootB!=a && rootB.pid!=rootB.id){rootB=getNode(rootB.pid); depthB++;}
		if(rootB==a){return a;}
		if(rootA!=rootB){return null;}
		while(depthA>depthB){a=getNode(a.pid); depthA--;}
		while(depthB>depthA){b=getNode(b.pid); depthB--;}
		while(a!=b){a=getNode(a.pid); b=getNode(b.pid);}
		return a;
	}
	
	/**
	 * Walks a finite parent chain to its self-parent root, normally Life.
	 * @param a Non-null node whose complete parent chain exists in this tree
	 * @return Borrowed root node
	 */
	public TaxNode highestAncestor(TaxNode a){
		assert(a!=null);
		while(a.id!=a.pid){a=getNode(a.pid);}
		return a;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Header Parsing        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Determine the TaxID of a String,
	 * without a loaded TaxTree.
	 * This only works if the literal TaxID is embedded in the String.
	 * @param header Typically a sequence header
	 * @return Decoded TaxID, or -1 if unsuccessful
	 */
	public static int parseHeaderStatic(String header){
		if(header.length()<3){return -1;}
		if(header.charAt(0)=='>'){header=header.substring(1);}
		if(!header.startsWith("tid|")){return -1;}
		int idx=3;
		int idx2=header.indexOf('|', 4);
		if(idx2<5){return -1;}
		int id=-1;
		try{
			id=Parse.parseInt(header, idx+1, idx2);
//			System.err.println("d"+", "+header.substring(idx+1, idx2));
		}catch(Throwable e){
//			System.err.println("e"+", "+header.substring(idx+1, idx2));
			//ignore
		}
		return id;
	}
	
	/** Extracts a literal taxID first, then optionally resolves other header forms.
	 * A literal ID is returned without requiring that it exists in the tree.
	 * @param header Sequence header or supported filename, or null
	 * @param tree Optional tree for accession/name/header lookup
	 * @return Parsed/resolved ID, or -1 when neither path resolves */
	public static int parseHeaderStatic2(String header, TaxTree tree){
		if(header==null){return -1;}
		int tid=parseTaxID(header);
		if(tid>=0 || tree==null){return tid;}
		TaxNode tn=tree.parseNodeFromHeader(header, true);
		return tn==null ? tid : tn.id;
	}
	
	/**
	 * Parse a taxonomic ID from a header or filename string.
	 * Recognizes "tid|N", "tid_N", "ncbi|N", and "ncbi_N" formats.
	 * Returns -1 if no valid (positive) taxid is found.
	 * Retries on false substring matches (e.g. "plastid_" contains "tid_").
	 * Prefix types are tried in the listed order, not their order within the input.
	 * @param s Header or filename text, or null
	 * @return First accepted ID for the preferred prefix type, or -1
	 */
	public static int parseTaxID(String s){
		if(s==null){return -1;}
		int tid=parseTaxIDPrefix(s, "tid|");
		if(tid<1){tid=parseTaxIDPrefix(s, "tid_");}
		if(tid<1){tid=parseTaxIDPrefix(s, "ncbi|");}
		if(tid<1){tid=parseTaxIDPrefix(s, "ncbi_");}
		return tid;
	}

	/** Searches successive occurrences of one prefix for a decimal ID.
	 * @param s Non-null input text
	 * @param prefix Literal marker preceding decimal digits
	 * @return Accepted positive ID below Integer.MAX_VALUE, or -1 */
	private static int parseTaxIDPrefix(String s, String prefix){
		int pos=s.indexOf(prefix);
		while(pos>=0){
			int start=pos+prefix.length();
			long id=0;
			int end=start;
			for(; end<s.length(); end++){
				char c=s.charAt(end);
				if(c<'0' || c>'9'){break;}
				id=id*10+(c-'0');
			}
			if(id>0 && id<Integer.MAX_VALUE){
				if(end-start<=10){return (int)id;}
				//Only long tokens need the leading-zero scan. Ten significant digits
				//fit in long; larger spans reject even if accumulation wrapped.
				while(start<end && s.charAt(start)=='0'){start++;}
				if(end-start<=10){return (int)id;}
			}
			pos=s.indexOf(prefix, pos+1);
		}
		return -1;
	}

	/** Parses the first legacy marker followed by one delimiter and a signed integer.
	 * Delegates to Parse's unchecked int accumulator: digit values after an optional
	 * leading minus are assertion-checked, and overflow is not checked.
	 * @param header Header containing the marker
	 * @param term Marker without its delimiter, such as tid
	 * @return Parsed integer, or -1 when the marker/delimiter/digit pattern is absent */
	private static int parseTaxIDLegacy(String header, String term){
		int loc=header.indexOf(term);
		if(loc<0){return -1;}
		int start=loc+term.length();
		if(header.length()<start+2){return -1;}
		final char c=header.charAt(start);
		if(Tools.isDigitOrSign(c)){return -1;}
		start++;
		final char d=header.charAt(start);
		if(!Tools.isDigitOrSign(d)){return -1;}
		int stop=start+1;
		while(stop<header.length() && Tools.isDigitOrSign(header.charAt(stop))){stop++;}
		return Parse.parseInt(header, start, stop);
	}
	
	/**
	 * Determine the TaxID of a String.
	 * @param header Typically a sequence header or file name
	 * @param bestEffort In some cases, try certain substrings if the name is not found.
	 * @return Borrowed resolved node, or null if no supported lookup resolves
	 */
	public TaxNode parseNodeFromHeader(String header, boolean bestEffort){
		if(header==null || header.length()<2){return null;}
		if(header.charAt(0)=='>'){header=header.substring(1);}
		if(header.endsWith(".gz")){header=header.substring(0, header.lastIndexOf('.'));}
		if(header.endsWith(".fa") || header.endsWith(".fna") || header.endsWith(".fasta") || header.endsWith(".faa")){
			header=header.substring(0, header.lastIndexOf('.'));
		}
		TaxNode tn;
		if(SILVA_MODE){
			tn=getNodeSilva(header, bestEffort);
		}else if(UNITE_MODE){
			tn=getNodeUnite(header, bestEffort);
		}else{
			final char delimiter=ncbiHeaderDelimiter(header);
			if(delimiter==' '){
				tn=getNodeNewStyle(header);
			}else{
				tn=getNodeOldStyle(header, delimiter);
				if(tn==null && delimiter=='|'){
//					System.err.println("A: "+header);
					int id=-1;
					String[] split=header.split("\\|");
					if(AccessionToTaxid.LOADED()){
						for(int i=0; i<split.length && id<0; i++){//Try accessions first
							if(AccessionToTaxid.isValidAccession(split[i])){
								id=AccessionToTaxid.get(split[i]);
							}
						}
					}
					for(int i=0; i<split.length && id<0; i++){//Then names
						id=parseNameToTaxid(split[i]);
					}
//					System.err.println("E: "+id);
					if(id>=0){tn=getNode(id);}
//					System.err.println("F: "+tn);
				}
			}
		}
		if(tn==null){
			int tid=parseTaxID(header);
			if(tid>0){tn=getNode(tid);}
		}
		return tn;
	}

	/**
	 * Guess the delimiter character in a String;
	 * typically assumed to be '|', '~', or ' '.
	 * The first pipe, tilde or whitespace wins; absence defaults to a space.
	 * @param header Non-null header
	 * @return Pipe, tilde or space; a leading pipe/tilde triggers an assertion
	 */
	public static char ncbiHeaderDelimiter(String header){
		for(int i=0; i<header.length(); i++){
			final char c=header.charAt(i);
			if(c=='|' || c=='~'){
				assert(i>0) : "i="+i+"; malformatted header '"+header+"'";
				return c;
			}else if(Character.isWhitespace(c)){
				return ' ';
			}
		}
		return ' ';
	}
	
	/**
	 * Parse a Silva header to a Node.
	 * @param s Silva header.
	 * @param bestEffort Try certain substrings if the name is not found.
	 * @return Node
	 */
	TaxNode getNodeSilva(String s, boolean bestEffort){
		if(s==null){return null;}
		if(s.length()>=5 && s.startsWith("tid") && (s.charAt(3)=='|' || s.charAt(3)=='~') && Tools.isDigit(s.charAt(4))){
			return getNodeOldStyle(s, s.charAt(3));
		}
		String[] split=Tools.semiPattern.split(s);
		
		int number=-1;
//		final boolean chloroplast=(split.length>1 && split[split.length-1].equals("Chloroplast"));
//		if(chloroplast){return null;}
		for(int i=split.length-1; number<0 && i>=0; i--){
			String last=split[i];
			int paren=last.indexOf('(');
			if(paren>=0){last=last.substring(0, paren);}
			last=last.trim();
			
			if(!last.startsWith("uncultured") && !last.startsWith("unidentified")){
				number=parseNameToTaxid(last);
			}
			
			if(number>=0){return getNode(number);}
			else if(!bestEffort){break;}
		}
		return null;
	}

	/**
	 * Parse a Unite header to a Node.
	 * @param s Unite header.
	 * @param bestEffort Try certain substrings if the name is not found.
	 * @return Node
	 */
	TaxNode getNodeUnite(String s, boolean bestEffort){
		if(s==null){return null;}
		if(s.length()>=5 && s.startsWith("tid") && (s.charAt(3)=='|' || s.charAt(3)=='~') && Tools.isDigit(s.charAt(4))){
			return getNodeOldStyle(s, s.charAt(3));
		}
		String[] split=Tools.pipePattern.split(s);

		int number=-1;
		String name=split[0];
		String acc=split[1];
		if(AccessionToTaxid.LOADED() && acc.length()>0){
			number=AccessionToTaxid.get(acc);
		}
		if(number<1){
			TaxNode tn=getNodeByName(name);
			if(tn!=null){number=tn.id;}
		}
		
		if(number>=0){return getNode(number);}
		return null;
	}
	
	/** Resolves a legacy delimited header by GI, literal taxID, IMG ID or fallback fields.
	 * Failed GI parsing is rethrown only if accession/name fallback also fails.
	 * IMG lookup requires the IMG map; name fallback requires initialized name indexes.
	 * @param s Non-null header without a leading FASTA marker
	 * @param delimiter Preferred separator; missing separators fall back to tilde then underscore
	 * @return Resolved tree node, or null when no lookup resolves */
	private TaxNode getNodeOldStyle(final String s, char delimiter){
		{
			int index=s.indexOf(delimiter);
			if(index<0){
				delimiter='~';
				index=s.indexOf(delimiter);
				if(index<0){
					delimiter='_';
					index=s.indexOf(delimiter);
				}
			}
			int number=-1;
			
			Throwable e=null;
			
			if(index==2 && s.length()>3 && s.startsWith("gi") && Tools.isDigit(s.charAt(3))){
//				System.err.println("Parsing gi number.");
				
				if(GiToTaxid.isInitialized()){
					try{
						number=GiToTaxid.parseGiToTaxid(s, delimiter);
					}catch(Throwable e2){
						e=e2;
					}
				}else{
					assert(!CRASH_IF_NO_GI_TABLE) : "To use gi numbers, you must load a gi table.\n"+s;
				}
//				if(number!=-1){System.err.println("number="+number);}
			}else if(index==3 && s.length()>4 && s.startsWith("tid") && Tools.isDigit(s.charAt(4))){
//				System.err.println("Parsing ncbi number.");
				number=GiToTaxid.parseTaxidNumber(s, delimiter);
			}else if(index==3 && s.length()>4 && s.startsWith("img") && Tools.isDigit(s.charAt(4))){
//				System.err.println("Parsing ncbi number.");
				long img=parseDelimitedNumber(s, delimiter);
				ImgRecord record=imgMap.get(img);
				number=(record==null ? -1 : record.taxID);
			}else if(index==4 && s.length()>5 && s.startsWith("ncbi") && Tools.isDigit(s.charAt(5))){//obsolete
//				System.err.println("Parsing ncbi number.");
				number=GiToTaxid.parseTaxidNumber(s, delimiter);
			}
			
			if(number<0 && index>=0 && (delimiter=='|' || delimiter=='~')){
				String[] split=(delimiter=='|' ? delimiterPipe.split(s) : delimiterTilde.split(s));
				if(AccessionToTaxid.LOADED()){
					number=parseAccessionToTaxid(split);
				}
				if(number<0){
					number=parseHeaderNameToTaxid(split);
				}
			}
			
			if(number<0 && e!=null){
				assert(false) : e;
				throw new RuntimeException(e);
			}
			
			//TaxServer code could go here...
			
			if(number>=0){return getNode(number);}
		}
		if(verbose){System.err.println("Can't process name "+s);}
		if(Tools.isNumeric(s)){return getNode(Integer.parseInt(s), true);}
		return null;
	}
	
	/** Parses an unsigned decimal after the first separator, stopping at separator, space or tab.
	 * Subsequent characters are assertion-checked digits; arithmetic overflow is unchecked.
	 * @param s Header text, or null
	 * @param delimiter Separator preceding the number and optionally ending it
	 * @return Parsed value, or -1 for null text, a missing number or a nondigit first character */
	static long parseDelimitedNumber(String s, char delimiter){
		if(s==null){return -1;}
		int i=0;
		while(i<s.length() && s.charAt(i)!=delimiter){i++;}
		i++;
		if(i>=s.length() || !Tools.isDigit(s.charAt(i))){return -1;}
		
		long number=0;
		while(i<s.length()){
			char c=s.charAt(i);
			if(c==delimiter || c==' ' || c=='\t'){break;}
			assert(Tools.isDigit(c)) : c+"\n"+s;
			number=(number*10)+(c-'0');
			i++;
		}
		return number;
	}
	
	/** Resolves an accession-first header, then tries its organism-name suffix.
	 * An unresolved bare numeric header uses merged-ID lookup; underscore-only text
	 * is retried with underscores changed to spaces.
	 * @param s Non-null header text without a leading FASTA marker
	 * @return Resolved tree node, or null when neither accession nor name resolves */
	private TaxNode getNodeNewStyle(final String s){
		
		int space=s.indexOf(' ');
		int number=-1;
		
		if(AccessionToTaxid.LOADED()){
			if(space>0){
				number=AccessionToTaxid.get(s.substring(0, space));
			}else{
				number=AccessionToTaxid.get(s);
			}
		}
		
		if(number<0 && space<0 && Tools.isNumeric(s)){return getNode(Integer.parseInt(s), true);}
		
		if(number<0 && space>0){
			number=parseNameToTaxid(s.substring(space+1));
		}
		
		if(number>-1){return getNode(number);}
		if(space<0 && s.indexOf('_')>0){
			return getNodeNewStyle(s.replace('_', ' '));
		}
		return null;
	}
	
	/** Resolves the fourth field of a split legacy header as an accession.
	 * @param split Non-null legacy header fields; lookup state must already be loaded
	 * @return Accession lookup result, or -1 when fewer than four fields are present */
	public int parseAccessionToTaxid(String[] split){
		if(split.length<4){
			return -1;
		}
		int ncbi=AccessionToTaxid.get(split[3]);
		return ncbi;
	}
	
	/** Resolves the fifth field of a split legacy header as an organism name.
	 * @param split Non-null legacy header fields; name indexes must already be initialized
	 * @return Unique name's TaxID, or -1 for a missing field, unknown name or ambiguity */
	public int parseHeaderNameToTaxid(String[] split){
		if(split.length<5){
			return -1;
		}
		return parseNameToTaxid(split[4]);
	}
	
	/**
	 * Returns the TaxID from the organism's scientific name (e.g. "Homo sapiens").
	 * Ambiguous results are rejected; use getNodesByNameExtended to inspect them.
	 * @param name Organism name.
	 * @return Unique organism TaxID, or -1 when missing or ambiguous
	 */
	public int parseNameToTaxid(String name){
//		assert(false) : name+", "+(nameMap==null)+", "+(nameMap==null ? 0 : nameMap.size());
		List<TaxNode> list=null;
		
		list=getNodesByNameExtended(name);
		
		if(list==null || list.size()>1){return -1;}
		return list.get(0).id;
	}
	
	/**
	 * Fetch nodes indicated by this name.
	 * @param name A taxonomic name delimited by space or underscore.
	 * @return Nodes corresponding to the name.
	 */
	public List<TaxNode> getNodesByNameExtended(String name){
		List<TaxNode> list=null;
		
		list=getNodesByName(name);
		if(list!=null){return list;}
		
		name=name.replaceAll("_", " ").trim();
		list=getNodesByName(name);
		if(list!=null){return list;}
		
		String[] split2=name.split(" ");
		
		if(split2.length>7){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6]+" "+split2[7];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		
		if(split2.length>6){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		
		if(split2.length>5){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		if(split2.length>4){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		if(split2.length>3){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		if(split2.length>2){
			String term=split2[0]+" "+split2[1]+" "+split2[2];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		if(split2.length>1){
			String term=split2[0]+" "+split2[1];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		if(split2.length>0){
			String term=split2[0];
			list=getNodesByName(term);
			if(list!=null){return list;}
		}
		
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Assorted Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Ascends toward a normal rank without stepping beyond it or onto Life.
	 * If the requested rank is absent, this may return a lower-ranked node.
	 * @param taxID Starting TaxID; unresolved/nonpositive values are returned unchanged
	 * @param taxLevel Desired normal rank, not an extended-rank index
	 * @return Selected TaxID, or the input when no node resolves
	 */
	public int promote(final int taxID, int taxLevel){
		TaxNode tn=null;
		tn=(taxID<1 ? null : getNode(taxID));
		tn=promote(tn, taxLevel);
		return (tn==null ? taxID : tn.id);
	}
	
	/**
	 * Ascends while below the desired normal rank, stopping before an overshoot,
	 * Life, a missing parent, or a self-parent. The result may remain below the target.
	 * @see #getNodeAtLevelExtended(int, int)
	 * @param tn Node in question
	 * @param taxLevel Desired minimum level
	 * @return Selected node, or null for a null input
	 */
	public TaxNode promote(TaxNode tn, int taxLevel){
		while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
			TaxNode temp=getNode(tn.pid);
			if(temp==null || temp==tn || temp.level>=TaxTree.LIFE || temp.level>taxLevel){break;}
			tn=temp;
		}
		return tn;
	}
	
	/**
	 * Resolves a native, auxiliary or merged ID and returns its stored parent ID.
	 * The out-of-range assertion follows getNode; unresolved IDs return -1.
	 * @param id TaxID of child node
	 * @return Parent's durable TaxID, or -1 when no node resolves
	 */
	public int getParentID(int id){
		final int idx=nodeIndex(id);
		assert(idx<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(idx<0 || idx>=nodes.length){return -1;}
		TaxNode tn=nodes[idx];
		if(tn==null && mergedMap!=null){tn=getNode(mergedMap.get(id), true);}
		return tn==null ? -1 : tn.pid;
	}
	
	/**
	 * Resolves native IDs, durable auxiliary IDs and chains of merged aliases.
	 * Private packed auxiliary slots are not accepted as public IDs. An ID beyond
	 * the addressable native range triggers an assertion, or returns null under -da.
	 * @param id External TaxID
	 * @return Borrowed node, or null for a missing/unaddressable ID
	 */
	public TaxNode getNode(int id){
		final int idx=nodeIndex(id);
		assert(idx<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(idx<0 || idx>=nodes.length){return null;}
		TaxNode tn=nodes[idx];
		if(tn!=null || mergedMap==null){return tn;}
		return getNode(mergedMap.get(id), true);
	}
	
	/**
	 * Gets the normal taxonomic rank of a resolved ID.
	 * @param id TaxID
	 * @return Taxonomic level or -1 if not found
	 */
	public int toLevel(int id){
		TaxNode tn=getNode(id);
		return tn==null ? -1 : tn.level;
	}
	
	/**
	 * Applies one merged-ID mapping without validating that the target node exists.
	 * Known durable auxiliary IDs are preserved. This does not recursively flatten
	 * alias chains; getNode performs recursive resolution when fetching a node.
	 * @param id TaxID to resolve
	 * @return Direct mapped ID, or the unchanged input when no mapping applies
	 */
	public int resolveID(int id){
		if(externalToInternalMap!=null && externalToInternalMap.contains(id)){return id;}
		if(mergedMap==null){return id;}
		int x=mergedMap.get(id);
		return x<0 ? id : x;
	}
	
	/**
	 * Resolves a node with optional suppression of the out-of-range assertion.
	 * Missing slots follow merged-ID chains, suppressing that assertion on later hops.
	 * @param id External TaxID
	 * @param skipAssertion Permit a too-new/out-of-range ID to return null under -ea
	 * @return Borrowed resolved node, or null when absent
	 */	
	public TaxNode getNode(int id, boolean skipAssertion){
		final int idx=nodeIndex(id);
		assert(skipAssertion || idx<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(idx<0 || idx>=nodes.length){return null;}
		TaxNode tn=nodes[idx];
		if(tn!=null || mergedMap==null){return tn;}
		return getNode(mergedMap.get(id), true);
	}

	/** Maps durable auxiliary IDs to private slots while hiding those slots as IDs.
	 * Native IDs otherwise pass through; this method does not resolve merged aliases
	 * or check the upper array bound.
	 * @param id External TaxID
	 * @return Candidate array index, or -1 for a negative ID or exposed private slot */
	int nodeIndex(int id){
		if(id<0){return -1;}
		if(externalToInternalMap!=null){
			int idx=externalToInternalMap.get(id);
			if(idx>=0){return idx;}
			if(id>=firstAuxIndex && id<nodes.length){return -1;}
		}
		return id;
	}
	
	/**
	 * Promotes toward a normal minimum rank, with DOMAIN as the ascent ceiling.
	 * @param id TaxID
	 * @param minLevel Minimum level
	 * @return Selected node, possibly below the requested minimum, or null
	 */
	public TaxNode getNodeAtLevel(int id, int minLevel){
		return getNodeAtLevel(id, minLevel, DOMAIN);
	}
	
	/**
	 * Promotes toward an extended minimum rank, with DOMAIN_E as the ascent ceiling.
	 * @param id TaxID
	 * @param minLevelE Minimum extended level
	 * @return Selected node, possibly below the requested minimum, or null
	 */
	public TaxNode getNodeAtLevelExtended(int id, int minLevelE){
		return getNodeAtLevelExtended(id, minLevelE, DOMAIN_E);
	}
	
	/**
	 * Converts normal rank indices to extended ranks and applies bounded promotion.
	 *
	 * @param id TaxID
	 * @param minLevel Desired minimum normal rank, in the normal-rank index range
	 * @param maxLevel Normal rank ceiling on parent steps
	 * @return Selected node, possibly outside the requested interval, or null
	 */
	public TaxNode getNodeAtLevel(int id, int minLevel, int maxLevel){
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
		return getNodeAtLevelExtended(id, minLevelExtended, maxLevelExtended);
	}
	
	/**
	 * Ascends while below the extended minimum, stopping before a parent exceeds
	 * the maximum. Missing parents and self-parent roots also stop the walk. The
	 * result can remain below the minimum; an initial node above the maximum is
	 * returned unchanged. Parent chains must be finite.
	 *
	 * @param id TaxID
	 * @param minLevelE Minimum extended level
	 * @param maxLevelE Maximum extended level
	 * @return Last accepted node, or null when the starting ID does not resolve
	 */
	public TaxNode getNodeAtLevelExtended(int id, int minLevelE, int maxLevelE){
		TaxNode tn=getNode(id);
		while(tn!=null && tn.pid!=tn.id && tn.levelExtended<minLevelE){
			TaxNode temp=getNode(tn.pid);
			if(temp==null || temp.levelExtended>maxLevelE){break;}
			tn=temp;
		}
		return tn;
	}
	
	/**
	 * Ascends toward an extended rank without accepting an overshooting parent.
	 * Negative rank limits return the input without lookup. Existing nodes must have
	 * intact parent chains; unlike promote, this method does not guard missing parents.
	 * @param taxID Starting TaxID
	 * @param taxLevelExtended Desired extended level
	 * @return Last accepted ID, or the input when no step is accepted
	 */
	public int getIdAtLevelExtended(int taxID, int taxLevelExtended){
		if(taxLevelExtended<0){return taxID;}
		TaxNode tn=getNode(taxID);
		while(tn!=null && tn.id!=tn.pid && tn.levelExtended<taxLevelExtended){
			tn=getNode(tn.pid);
			if(tn.levelExtended>taxLevelExtended){break;}
			taxID=tn.id;
		}
		return taxID;
	}

	/**
	 * Looks up a name, trying original case before lowercase matching.
	 * Multiple matches trigger an assertion; under -da the smallest
	 * minAncestorLevelIncludingSelf value wins, retaining the first match on ties.
	 * @param s Organism name.
	 * @return Borrowed matching node, or null when no match exists
	 */
	public TaxNode getNodeByName(String s){
		List<TaxNode> list=getNodesByName(s, false);
		if(list==null){list=getNodesByName(s, true);}
		if(list==null || list.isEmpty()){return null;}
		if(list.size()==1){return list.get(0);}
		assert(false) : "Found multiple nodes for '"+s+"':\n"+list+"\n";
		TaxNode a=list.get(0);
		for(int i=1; i<list.size(); i++){
			TaxNode b=list.get(i);
			//Keep the most specific node
//			if(a==null || (b!=null && b.minAncestorLevelIncludingSelf()<a.minAncestorLevelIncludingSelf())){//not necessary
			if(b.minAncestorLevelIncludingSelf()<a.minAncestorLevelIncludingSelf()){
				a=b;
			}
		}
		return a;
	}
	
	/**
	 * Looks up a name with original-case precedence and a lowercase fallback.
	 * Underscores are normalized to spaces; multiword inputs can also match ancestry.
	 * Returned lists may be shared with the name index and must not be modified.
	 * @param s Organism name.
	 * @return Matching borrowed nodes in an existing or filtered list, or null
	 */
	public List<TaxNode> getNodesByName(String s){
		List<TaxNode> list=getNodesByName(s, false);
		if(list==null){list=getNodesByName(s, true);}
		return list;
	}
	
	/**
	 * Fetch a map of names to nodes.  If absent, create it first.
	 * @param lowercase If true, return the map with lowercase keys.
	 * @return Map of names to nodes.
	 */
	private HashMap<String, ArrayList<TaxNode>> getMap(boolean lowercase){
		HashMap<String, ArrayList<TaxNode>> map=(lowercase ? nameMapLower : nameMap);
		if(map==null){
			synchronized(this){hashNames(true);}
			map=(lowercase ? nameMapLower : nameMap);
		}
		assert(map!=null) : "Tax names were not hashed.";
		return map;
	}
	
	/**
	 * Internal method for name-based node lookup.
	 * @param s Name to search for
	 * @param lowercase Use lowercase matching
	 * @return List of matching nodes
	 */
	private List<TaxNode> getNodesByName(String s, boolean lowercase){
		if(s==null){return null;}
		if(s.indexOf('_')>=0){s=s.replace('_', ' ');}
		if(lowercase){s=s.toLowerCase();}
		final HashMap<String, ArrayList<TaxNode>> map=getMap(lowercase);
		ArrayList<TaxNode> list=map.get(s);
		if(list!=null){return list;}
		
		
		if(s.indexOf('_')<0 && s.indexOf(' ')<0){return null;}
		String[] split=delimiter2.split(lowercase ? s.toLowerCase() : s, 8);
		list=map.get(split[split.length-1]);
		if(list==null){return list;}
		
		int matchCount=0;
		for(TaxNode tn : list){
			if(tn.matchesName(split, split.length-1, this)){matchCount++;}
		}
		if(matchCount==list.size()){return list;}
		if(matchCount<1){return null;}
		ArrayList<TaxNode> hits=new ArrayList<TaxNode>(matchCount);
		for(TaxNode tn : list){
			if(tn.matchesName(split, split.length-1, this)){hits.add(tn);}
		}
		return hits;
	}
	/**
	 * Collects the starting node and its ancestors in leaf-to-root order, excluding
	 * the self-parent root. Missing IDs and the root itself produce an empty list.
	 * @param id TaxID
	 * @return Newly allocated list containing borrowed node references
	 */
	public ArrayList<TaxNode> getAncestors(int id){
		TaxNode current=getNode(id);
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		while(current!=null && current.pid!=current.id){//ignores root
			list.add(current);
			current=getNode(current.pid);
		}
		//optionally add root here
		return list;
	}
	
	/**
	 * Increments counts for multiple TaxIDs after sorting and deduplicating ids in place.
	 * The two lists must be distinct and exclusively owned during this operation.
	 * @param ids TaxIDs, replaced with their sorted unique values
	 * @param counts Scratch list replaced with matching multiplicities
	 * @param sync Synchronize raw counter updates on this tree
	 */
	public void increment(IntList ids, IntList counts, boolean sync){
		
		ids.sort();
		ids.getUniqueCounts(counts);
		
		if(!sync){
			for(int i=0; i<ids.size; i++){
				int id=ids.get(i);
				int count=counts.get(i);
				incrementRaw(id, count);
			}
		}else{
			synchronized(this){
				for(int i=0; i<ids.size; i++){
					int id=ids.get(i);
					int count=counts.get(i);
					incrementRaw(id, count);
				}
			}
		}
	}
	
	/**
	 * Increment raw count for a TaxID.
	 * @param id TaxID
	 * @param amt Amount to increment
	 */
	public void incrementRaw(int id, long amt){
		final int idx=nodeIndex(id);
		assert(idx>=0 && idx<nodes.length) : "TaxID "+id+" is out of range."+(idx<0 ? "" : "  Possibly the taxonomy data needs to be updated.");
		assert(nodes[idx]!=null) : "No node for TaxID "+id+"; possibly the taxonomy data needs to be updated.";
		nodes[idx].incrementRaw(amt);
	}
	
	/** Percolate counts upward through the entire tree */
	public void percolateUp(){
		for(int i=0; i<treeLevelsExtended.length; i++){
			percolateUp(i);
		}
	}
	
	/** Percolate counts upward from a specific level.
	 * @param fromLevel Starting level */
	public void percolateUp(final int fromLevel){
		//FIXED [tax/TaxTree#002]: the old stratified push (n.incrementSum(countRaw); parent.incrementSum(n.countSum))
		//assumed a parent always sits in a strictly-higher stratum than its children, so that a child's countSum was
		//complete before its parent forwarded it. NO-RANK internal nodes (levelExtended==0, the lowest stratum,
		//processed FIRST) break that: a no-rank node forwarded its incomplete sum, then its ranked children (later
		//strata) pushed into it too late to reach the no-rank node's ancestors. Empirically this dropped 92% of counts
		//on the real NCBI tree (root got 216,768 of 2,835,394). Reroute each node's own countRaw through the correct
		//order-independent parent-chain walk percolateUp(node,amt); since every node is in exactly one stratum,
		//percolateUp() (which calls this for every level) still visits each node exactly once -> no double count.
		final TaxNode[] stratum=treeLevelsExtended[fromLevel];
		for(final TaxNode n : stratum){
			if(n!=null){percolateUp(n, n.countRaw);}
		}
	}
	
	/** Adds a signed amount to countSum for this node and every ancestor, including root.
	 * Nonzero updates require a complete, acyclic parent chain ending at a self-parent
	 * root. Zero is a no-op without inspecting the node.
	 * @param node Starting node, non-null for a nonzero update
	 * @param amt Signed count adjustment, including negative subtraction */
	public void percolateUp(TaxNode node, long amt){
		if(amt==0){return;}
		if(verbose){System.err.println("percolateUp("+amt+") node: "+node);}
		while(node.id!=node.pid){
			node.incrementSum(amt);
			node=getNode(node.pid);
		}
		node.incrementSum(amt);
	}
	
	/**
	 * Gather nodes with count at least the limit.
	 * Subtracts each selected node's cumulative count from its ancestors.
	 * @param limit Minimum count threshold
	 * @return List of qualifying nodes
	 */
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit){
		return gatherNodesAtLeastLimit(limit, NO_RANK, LIFE);
	}
	
	/**
	 * Gather nodes with count at least the limit within level range.
	 * Visits extended levels in ascending order and subtracts each selected node's
	 * cumulative count from its ancestors; this mutates the tree's count state.
	 *
	 * @param limit Minimum count threshold
	 * @param minLevel Minimum taxonomic level
	 * @param maxLevel Maximum taxonomic level
	 * @return List of qualifying nodes
	 */
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit, final int minLevel, final int maxLevel){
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
//		assert(false) : limit+", "+minLevel+", "+maxLevel+", "+minLevelExtended+", "+maxLevelExtended;
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		for(int i=minLevelExtended; i<nodesPerLevelExtended.length && i<=maxLevelExtended; i++){
			list.addAll(gatherNodesAtLeastLimitExtended(i, limit));
		}
		Shared.sort(list, TaxNode.countComparator);
		return list;
	}
	
	/**
	 * Gather nodes at extended level with count at least the limit.
	 * Subtracts each selected node's cumulative count from its ancestors.
	 * Mutates the tree's count state.
	 * @param fromLevelExtended Extended level to gather from
	 * @param limit Minimum count threshold
	 * @return List of qualifying nodes
	 */
	public ArrayList<TaxNode> gatherNodesAtLeastLimitExtended(final int fromLevelExtended, final long limit){
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		final TaxNode[] stratum=treeLevelsExtended[fromLevelExtended];
		for(final TaxNode n : stratum){
			if(n.countSum>=limit){
				list.add(n);
				TaxNode parent=getNode(n.pid);
				if(n!=parent){
					percolateUp(parent, -n.countSum);//123 This was negative for some reason
				}
			}
		}
		Shared.sort(list, TaxNode.countComparator);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Static Initializers      ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds normal-rank aliases in lowercase and uppercase, including clade as no-rank.
	 * @return Newly allocated rank-name map */
	private static HashMap<String, Integer> makeLevelMap(){
		HashMap<String, Integer> map=new HashMap<String, Integer>(31);
		for(int i=0; i<taxLevelNames.length; i++){
			map.put(taxLevelNames[i], i);
			map.put(taxLevelNames[i].toUpperCase(), i);
		}
		map.put("clade", NO_RANK);
		map.put("clade".toUpperCase(), NO_RANK);
		return map;
	}

	/** Builds lowercase and uppercase abbreviations for normal ranks.
	 * @return Newly allocated short-name map */
	private static HashMap<String, Integer> makeShortLevelMap(){
		HashMap<String, Integer> map=new HashMap<String, Integer>(31);
		for(int i=0; i<taxLevelNamesShort.length; i++){
			map.put(taxLevelNamesShort[i], i);
			map.put(taxLevelNamesShort[i].toUpperCase(), i);
		}
		return map;
	}

	/** Builds lowercase and uppercase extended-rank names, including clade as no-rank.
	 * @return Newly allocated extended-rank map */
	private static HashMap<String, Integer> makeLevelMapExtended(){
		HashMap<String, Integer> map=new HashMap<String, Integer>(129);
		for(int i=0; i<taxLevelNamesExtended.length; i++){
			map.put(taxLevelNamesExtended[i], i);
			map.put(taxLevelNamesExtended[i].toUpperCase(), i);
		}
		map.put("clade", NO_RANK_E);
		map.put("clade".toUpperCase(), NO_RANK_E);
		return map;
	}

	/** Maps normal and intermediate rank names to their containing normal rank.
	 * Includes lowercase and uppercase spellings and clade as no-rank.
	 * @return Newly allocated full-name-to-normal-rank map */
	private static HashMap<String, Integer> makeAltLevelMap(){
		HashMap<String, Integer> map=new HashMap<String, Integer>(129);
		for(int i=0; i<taxLevelNames.length; i++){
			map.put(taxLevelNames[i], i);
			map.put(taxLevelNames[i].toUpperCase(), i);
		}
		map.put("clade", NO_RANK);
		map.put("clade".toUpperCase(), NO_RANK);
		
		//Add synonyms
//		map.put("subfamily", map.get("family"));
//		map.put("tribe", map.get("family"));
//		map.put("varietas", map.get("subspecies"));
//		map.put("subgenus", map.get("genus"));
//		map.put("forma", map.get("subspecies"));
//		map.put("species group", map.get("genus"));
//		map.put("species subgroup", map.get("genus"));
//		map.put("cohort", map.get("class"));
//		map.put("subclass", map.get("class"));
//		map.put("infraorder", map.get("order"));
//		map.put("superorder", map.get("class"));
//		map.put("subphylum", map.get("phylum"));
//		map.put("infraclass", map.get("class"));
//		map.put("superkingdom", map.get("division"));
//		map.put("parvorder", map.get("order"));
//		map.put("superclass", map.get("phylum"));
//		map.put("superphylum", map.get("kingdom"));
//		map.put("subkingdom", map.get("kingdom"));
//		map.put("superfamily", map.get("order"));
//		map.put("superkingdom", map.get("domain"));
//		map.put("suborder", map.get("order"));
//		map.put("subtribe", map.get("family"));
		
		for(String[] array : taxLevelNamesExtendedMatrix){
			String head=array[array.length-1];
			Integer value=map.get(head);
			assert(value!=null) : head+", "+value+", "+Arrays.toString(array);
			for(String key : array){
				if(key!=head){
					assert(!map.containsKey(key)) : "Map already contains key "+key+": "+Arrays.toString(array);
					map.put(key, value);
					map.put(key.toUpperCase(), value);
				}
			}
		}
		
		return map;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Size             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the loaded base count assigned directly to a node.
	 * @param tn Node whose ID indexes the loaded size table, or null
	 * @return Zero for null; -1 if the table is unloaded; otherwise the stored count,
	 * with absent or negative values reported as zero */
	public long toSize(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSizeMap==null){return -1L;}
		final long x=refseqSizeMap.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Returns the loaded cumulative base count for a node and its descendants.
	 * This reads the summary; it does not traverse or recompute ancestry.
	 * @param tn Node whose ID indexes the loaded summary, or null
	 * @return Zero for null; -1 if the table is unloaded; otherwise the stored count,
	 * with absent or negative values reported as zero */
	public long toSizeC(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSizeMapC==null){return -1L;}
		final long x=refseqSizeMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Returns the loaded sequence count assigned directly to a node.
	 * @param tn Node whose ID indexes the loaded summary, or null
	 * @return Zero for null; -1 if the table is unloaded; otherwise the stored count,
	 * with absent or negative values reported as zero */
	public int toSeqs(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSeqMap==null){return -1;}
		final int x=refseqSeqMap.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Returns the loaded cumulative sequence count for a node and its descendants.
	 * @param tn Node whose ID indexes the loaded summary, or null
	 * @return Zero for null; -1 if the table is unloaded; otherwise the stored count,
	 * with absent or negative values reported as zero */
	public long toSeqsC(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSeqMapC==null){return -1L;}
		final long x=refseqSeqMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Returns the loaded cumulative node count for a TaxID.
	 * @param tn Node whose ID indexes the loaded summary, or null
	 * @return Zero for null; -1 if the table is unloaded; otherwise the stored count,
	 * with absent or negative values reported as zero */
	public int toNodes(TaxNode tn){
		if(tn==null){return 0;}
		if(nodeMapC==null){return -1;}
		final int x=nodeMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	/**
	 * Fills refseqSizeMap, refseqSizeMapC, etc. from a file containing the summary.
	 * @param fname Size file name
	 */
	public void loadSizeFile(String fname){
		if(fname==null){return;}
		assert(refseqSizeMap==null);
		refseqSizeMap=new IntLongHashMap();
		refseqSizeMapC=new IntLongHashMap();
		refseqSeqMap=new IntHashMap();
		refseqSeqMapC=new IntLongHashMap();
		nodeMapC=new IntHashMap();
		
		final ByteFile bf=ByteFile.makeByteFile(fname, true);
		final byte delimiter='\t';
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0 && line[0]!='#'){
				int a=0, b=0;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 0: "+new String(line);
				int tid=Parse.parseInt(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 1: "+new String(line);
				long size=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 2: "+new String(line);
				long csize=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 3: "+new String(line);
				int seqs=Parse.parseInt(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 4: "+new String(line);
				long cseqs=Parse.parseLong(line, a, b);
				b++;
				a=b;

				while(b<line.length && line[b]!=delimiter){b++;}
				assert(b>a) : "Missing field 5: "+new String(line);
				int cnodes=Parse.parseInt(line, a, b);
				b++;
				a=b;

				if(refseqSizeMap!=null && size>0){refseqSizeMap.put(tid, size);}
				if(refseqSizeMapC!=null && csize>0){refseqSizeMapC.put(tid, csize);}
				if(refseqSeqMap!=null && seqs>0){refseqSeqMap.put(tid, seqs);}
				if(refseqSeqMapC!=null && cseqs>0){refseqSeqMapC.put(tid, cseqs);}
				if(nodeMapC!=null && cnodes>0){nodeMapC.put(tid, cnodes);}
			}
		}
		bf.close();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             IMG              ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Convert IMG ID to TaxID.
	 * @param img IMG identifier
	 * @return TaxID or -1 if not found
	 */
	public static int imgToTaxid(long img){
		ImgRecord ir=imgMap.get(img);
		return ir==null ? -1 : ir.taxID;
	}
	
	/**
	 * Convert IMG ID to TaxNode.
	 * @param img IMG identifier
	 * @return TaxNode or null if not found
	 */
	public TaxNode imgToTaxNode(long img){
		int tid=imgToTaxid(img);
		return tid<1 ? null : getNode(tid);
	}
	
	
	/**
	 * Load IMG mapping data from file.
	 *
	 * @param fname IMG file name
	 * @param storeName Store organism names
	 * @param outstream Print stream for output
	 * @return Number of records loaded
	 */
	public static int loadIMG(String fname, boolean storeName, PrintStream outstream){
		assert(imgMap==null);
		if(fname==null){return 0;}
		ImgRecord.storeName=storeName;
		if(outstream!=null){System.err.println("Loading IMG.");}
		Timer t=new Timer(outstream, false);
		ImgRecord[] array=ImgRecord.toArray(fname, IMG_HQ);
		int x=loadIMG(array);
		t.stopAndPrint();
		return x;
	}
	
	/**
	 * Load IMG mapping data from array.
	 * @param array Array of IMG records
	 * @return Number of records loaded
	 */
	public static int loadIMG(ImgRecord[] array){
		assert(imgMap==null);
		imgMap=new HashMap<Long, ImgRecord>((int)(array.length*1.5));
		for(ImgRecord record : array){
			imgMap.put(record.imgID, record);
		}
		return imgMap.size();
	}
	
	/**
	 * @deprecated Parse taxonomic level string
	 * @param b Level string
	 * @return Level number
	 */
	@Deprecated
	public static int parseLevel(String b){
		final int level;
		if(b==null){level=-1;}
		else if(Tools.isNumeric(b.charAt(0))){
			level=Integer.parseInt(b);
		}else{
			level=stringToLevel(b.toLowerCase());
		}
		return level;
	}
	
	/**
	 * Parse extended taxonomic level string.
	 * @param b Level string
	 * @return Extended level number
	 */
	public static int parseLevelExtended(String b){
		final int level;
		if(b==null){level=-1;}
		else if(Tools.isNumeric(b.charAt(0))){
			level=levelToExtended(Integer.parseInt(b));
		}else{
			level=stringToLevelExtended(b.toLowerCase());
		}
		return level;
	}
	
	/**
	 * Check if TaxID represents unclassified organism.
	 * @param tid TaxID to check
	 * @return true if unclassified
	 */
	public boolean isUnclassified(int tid){
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			if(tn.isUnclassified()){return true;}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	/**
	 * Check if TaxID represents environmental sample.
	 * @param tid TaxID to check
	 * @return true if environmental sample
	 */
	public boolean isEnvironmentalSample(int tid){
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			if(tn.isEnvironmentalSample()){return true;}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	/**
	 * Check if TaxID represents a virus.
	 * @param tid TaxID to check
	 * @return true if virus
	 */
	public boolean isVirus(int tid){
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			if(tn.id==VIRUSES_ID){return true;}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	/**
	 * Get bit mask of defined taxonomic levels for TaxID.
	 * @param tid TaxID
	 * @return Bitmask of defined levels
	 */
	public long definedLevels(int tid){
		long levels=0;
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			levels=levels|(1L<<tn.level);
			tn=getNode(tn.pid);
		}
		return levels;
	}
	
	/**
	 * Get bit mask of defined extended taxonomic levels for TaxID.
	 * @param tid TaxID
	 * @return Bitmask of defined extended levels
	 */
	public long definedLevelsExtended(int tid){
		long levels=0;
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			levels=levels|(1L<<tn.levelExtended);
			tn=getNode(tn.pid);
		}
		return levels;
	}
	
	/** Removes a resolved node's parent-name prefix when followed by one space.
	 * For example, "Homo sapiens" with parent "Homo" becomes "sapiens".
	 * Case is preserved; a name without the exact parent prefix is unchanged.
	 * @param tid TaxID resolved through getNode
	 * @return Remaining name, unchanged name, or null if no node resolves
	 */
	public String mononomial(int tid){return mononomial(getNode(tid));}
	/**
	 * Removes an exact parent-name prefix and its following space, preserving case.
	 * Missing parents and names without that prefix leave the name unchanged.
	 * @param tn Node whose name is shortened, or null
	 * @return Remaining name, unchanged name, or null for a null node
	 */
	public String mononomial(TaxNode tn){
		if(tn==null){return null;}
		String name=tn.name;
		if(name.indexOf(' ')<0){return name;}
		TaxNode parent=getNode(tn.pid);
		if(parent==null){return name;}
		String pname=parent.name;
		if(name.length()>pname.length() && name.charAt(pname.length())==' ' && name.startsWith(pname)){
			name=name.substring(pname.length()+1);
		}
		return name;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Nodes indexed by native TaxID or private packed auxiliary slot; use getNode for IDs. */
	public final TaxNode[] nodes;
	
	/** Number of nodes per normal level */
	public final int[] nodesPerLevel=new int[taxLevelNames.length];
	
	/** Number of nodes per extended level */
	public final int[] nodesPerLevelExtended=new int[taxLevelNamesExtended.length];
	
	/** Number of nodes in the tree */
	public final int nodeCount;
	
	/** Maps old TaxIDs to new TaxIDs */
	public final IntHashMap mergedMap;

	/** Maps durable auxiliary TaxIDs to their packed runtime node indices. */
	public final IntHashMap externalToInternalMap;

	/** First packed auxiliary index; packed IDs are private and not accepted as TaxIDs. */
	public final int firstAuxIndex;

	/** Arrays of all nodes at a given taxonomic level (extended) */
	public final TaxNode[][] treeLevelsExtended=new TaxNode[taxLevelNamesExtended.length][];
	
	/** Map of names to nodes */
	HashMap<String, ArrayList<TaxNode>> nameMap;
	/** Map of lowercase names to nodes */
	HashMap<String, ArrayList<TaxNode>> nameMapLower;
	/** Map of nodes to child nodes */
	HashMap<TaxNode, ArrayList<TaxNode>> childMap;
	/** Returns the current name index without constructing it.
	 * The returned map and its lists are borrowed mutable state.
	 * @return Existing exact-name index, or null before names are hashed */
	public HashMap<String, ArrayList<TaxNode>> nameMap(){return nameMap;}
	
	/** @deprecated Minimum valid taxa count */
	@Deprecated
	public int minValidTaxa=0; //TODO: Remove (will break serialization)
	
	/** Infer ranks for no-rank nodes, when possible */
	public boolean simplify=true;
	/** See simplify() for details, works in conjunction with simplify */
	public boolean reassign=true;
	/** Discard no-rank nodes */
	public boolean skipNorank=false;
	/** Maximum level for rank inference */
	public int inferRankLimit=0;//levelMap.get("species");
	
	//Node Statistics
	/** Number of bases assigned to this TaxID in RefSeq */
	private IntLongHashMap refseqSizeMap;
	/** Number of bases assigned to this TaxID and descendants in RefSeq */
	private IntLongHashMap refseqSizeMapC;
	/** Number of sequences assigned to this TaxID in RefSeq */
	private IntHashMap refseqSeqMap;
	/** Number of sequences assigned to this TaxID and descendants in RefSeq */
	private IntLongHashMap refseqSeqMapC;
	/** Number of descendant nodes, inclusive, for each TaxID */
	private IntHashMap nodeMapC;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Assign levels to unranked nodes below species level, when possible */
	public static boolean assignStrains=true;
	/** Assume headers are in Silva format */
	public static boolean SILVA_MODE=false;
	/** Assume headers are in Unite format */
	public static boolean UNITE_MODE=false;
	/** Probably unnecessary at this point...  present for legacy reasons */
	public static boolean CRASH_IF_NO_GI_TABLE=true;

	/** Enable verbose output */
	public static boolean verbose=false;
	/** Show warning messages */
	public static boolean SHOW_WARNINGS=false;
	
	/** Maps IMG IDs to records from the dump file */
	private static HashMap<Long, ImgRecord> imgMap;

	/** Set to false if the tree is expected to be mutated.
	 * TODO: Separate mutable fields such as counters from the shared topology.
	 */
	public static boolean ALLOW_SHARED_TREE=true;
	
	/** Load the bundled auxiliary fake-TaxID tree when loading a main TaxTree. */
	public static boolean LOAD_FAKE_TREE=false;
	
	/** Bundled auxiliary fake-TaxID tree file. */
	public static final String defaultFakeTreeFile="?TaxTreeExtra.tsv.gz";
	
	/** Current auxiliary fake-TaxID tree path. */
	public static String fakeTreeFile=defaultFakeTreeFile;
	
	/** Universal location of the shared TaxTree used by various classes */
	private static TaxTree sharedTree;
	
	/** Returns the currently published shared tree without loading or hashing it.
	 * @return Shared tree, or null before publication */
	public static TaxTree getTree(){return sharedTree;}
	
	/** Loads or reuses the shared tree, requesting name and dot-format indexes.
	 * @return Shared tree, or null when shared-tree caching is disabled */
	public static TaxTree sharedTree(){return sharedTree(defaultTreeFile(), true, true, System.err);}
	
	/**
	 * Reuses the published shared tree or loads one when no tree is cached.
	 * An existing tree takes precedence over fname. Name hashing is added only
	 * when requested and the exact-name index is absent.
	 * @param fname Initial tree path, auto for the default, or null to avoid loading
	 * @param hashNames Request name indexes
	 * @param hashDotFormat Include dot-format names when constructing those indexes
	 * @param outstream Optional destination for loading/hash timing messages
	 * @return Shared tree, or null when caching is disabled or none is loaded
	 * TODO: Review safe publication of the lazily initialized shared tree and indexes.
	 */
	private static TaxTree sharedTree(String fname, boolean hashNames, boolean hashDotFormat, PrintStream outstream){
		if(!ALLOW_SHARED_TREE){return null;}
		if(sharedTree==null && fname!=null){
			if("auto".equalsIgnoreCase(fname)){fname=defaultTreeFile();}
			synchronized(TaxTree.class){
				if(sharedTree==null){
//					if(outstream!=null){outstream.println("Loading tax tree.");}
					if(outstream!=null){outstream.println("Loading tax tree from "+fname);} //Useful for debugging
					Timer t=new Timer(outstream, false);
					setSharedTree(loadTaxTreeFile(fname), hashNames, hashDotFormat);
					t.stopAndPrint();
				}
			}
		}
		if(sharedTree!=null && hashNames && sharedTree.nameMap==null){
			synchronized(sharedTree){
				if(sharedTree.nameMap==null){
					if(outstream!=null){outstream.println("Hashing names.");}
					Timer t=new Timer(outstream, false);
					sharedTree.hashNames(hashDotFormat);
					t.stopAndPrint();
				}
			}
		}
		return sharedTree;
	}
	
	/** Publishes the initial shared tree and optionally builds missing name indexes.
	 * Requires shared caching enabled and no prior shared tree; assertions enforce
	 * those preconditions. With assertions disabled, a prior tree is replaced.
	 * @param tree Non-null tree to publish
	 * @param hashNames Build name indexes if absent
	 * @param hashDotFormat Include dot-format names when constructing the indexes
	 */
	private static synchronized void setSharedTree(TaxTree tree, boolean hashNames, boolean hashDotFormat){
		assert(ALLOW_SHARED_TREE);
		assert(sharedTree==null);
		sharedTree=tree;
		if(hashNames && sharedTree.nameMap==null){
			synchronized(sharedTree){
				if(sharedTree.nameMap==null){
					sharedTree.hashNames(hashDotFormat);
				}
			}
		}
	}
	
	/**
	 * Tests whether an extended rank is the standard representative of its normal rank.
	 * For example, {@code isSimple(PHYLUM_E)} is true, while
	 * {@code isSimple(stringToLevelExtended("subphylum"))} and
	 * {@code isSimple(NO_RANK_E)} are false.
	 * @param levelExtended The extended level to test.
	 * @return True if this level is not no-rank, and the names of the normal and extended levels match.
	 */
	public static boolean isSimple(int levelExtended){
		int level=extendedToLevel(levelExtended);
		return levelExtended!=NO_RANK_E && (levelExtended==levelToExtended(level));
	}
	
	/**
	 * Tests standard ranks plus the explicit subspecies, strain and substrain ranks.
	 * For example, {@code isSimple2(PHYLUM_E)} and {@code isSimple2(SUBSTRAIN_E)}
	 * are true; {@code isSimple2(NO_RANK_E)} is false.
	 * @param levelExtended The extended level to test.
	 * @return True for a standard rank or one of those explicit additional ranks
	 */
	public static boolean isSimple2(int levelExtended){
		int level=extendedToLevel(levelExtended);
		return levelExtended!=NO_RANK_E && (levelExtended==levelToExtended(level) 
				|| levelExtended==STRAIN_E || levelExtended==SUBSPECIES_E || levelExtended==SUBSTRAIN_E);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Maps an exact accepted full rank name to its normal rank group.
	 * @param s Full rank name, including supported intermediate ranks
	 * @return Normal rank number
	 * @throws NullPointerException If the rank name is absent from the map */
	public static final int stringToLevel(String s){return altLevelMap.get(s);}
	/**
	 * Get the number for the short level name.
	 * @param s Short level name
	 * @return Level number
	 * @throws NullPointerException If the abbreviation is absent from the map
	 */
	public static final int shortStringToLevel(String s){return shortLevelMap.get(s);}
	/**
	 * Check if extended level map contains this name.
	 * @param s Level name
	 * @return true if contained
	 */
	public static final boolean levelMapExtendedContains(String s){return levelMapExtended.containsKey(s);}
	/** Maps an exact extended rank name to its extended rank number.
	 * @param s Extended rank name
	 * @return Extended rank number
	 * @throws NullPointerException If the rank name is absent from the map */
	public static final int stringToLevelExtended(String s){return levelMapExtended.get(s);}
	/** Returns the name of a normal rank.
	 * @param x Normal rank index in the supported range
	 * @return Corresponding normal rank name */
	public static final String levelToString(int x){return taxLevelNames[x];}
	/** Returns the name of an extended rank.
	 * @param x Extended rank index in the supported range
	 * @return Corresponding extended rank name */
	public static final String levelToStringExtended(int x){return taxLevelNamesExtended[x];}
	/** Returns the short name of a normal rank, such as s for species.
	 * @param x Normal rank index in the supported range
	 * @return Corresponding rank abbreviation */
	public static final String levelToStringShort(int x){return taxLevelNamesShort[x];}
	
	/** Normal, aka canonical, aka simple tax level names */
	private static final String[] taxLevelNames=new String[]{
		"no rank", "subspecies", "species", "genus",
		"family", "order", "class", "phylum",
		"kingdom", "superkingdom", "domain", "life"
	};
	/** Number of supported normal rank slots, including no-rank and Life. */
	public static final int numTaxLevelNames=taxLevelNames.length;
	
	/** 
	 * Definitive representation of all NCBI taxonomic level names.
	 * All levels used by NCBI must be present here, or parsing a new NCBI tax tree will crash.
	 * The first dimension maps normal ranks to extended ranks.
	 * Both dimensions are ordered ascending.
	 * TODO: More than 63 names would exceed the original-rank flag and level-bitmask capacity.
	 */
	private static final String[][] taxLevelNamesExtendedMatrix=new String[][]{
		{"no rank"},
		{"subgenotype", "genotype", "substrain", "isolate", "strain", "pathotype", "pathogroup", 
			"biotype", "serotype", "serogroup", "morph", "forma specialis", "forma", "subvariety", "varietas", 
			"subspecies"},
		{"species"},
		{"species subgroup", "species group", "series", "subsection", "section", "subgenus", "genus"},
		{"subtribe", "tribe", "subfamily", "family"},
		{"superfamily", "parvorder", "infraorder", "suborder", "order"},
		{"superorder", "subcohort", "cohort", "infraclass", "subclass", "class"},
		{"superclass", "subdivision", "division", "subphylum", "phylum"},
		{"superphylum", "subkingdom", "kingdom"},
		{"superkingdom"},
		{"realm", "domain"},
		{"acellular root", "cellular root", "life"}
	};
	
	/** Extended tax level names as a 1D array */
	private static final String[] taxLevelNamesExtended=makeNamesExtended();
	/** Number of extended tax levels */
	public static final int numTaxLevelNamesExtended=taxLevelNamesExtended.length;
	
	/** Flattens rank groups in their declared order.
	 * @return Newly allocated array of all extended rank names */
	private static final String[] makeNamesExtended(){
		ArrayList<String> list=new ArrayList<String>();
		for(String[] s : taxLevelNamesExtendedMatrix){
			for(String ss : s){
				list.add(ss);
			}
		}
		return list.toArray(new String[0]);
	}
	
	/** Abbreviations of tax level names, mainly for semicolon form */
	private static final String[] taxLevelNamesShort=new String[]{
			"nr", "ss", "s", "g",
			"f", "o", "c", "p",
			"k", "sk", "d", "l"
	};
	
	/** Normal tax level numbers as constants */
	public static final int NO_RANK=0, SUBSPECIES=1, SPECIES=2, GENUS=3,
			FAMILY=4, ORDER=5, CLASS=6, PHYLUM=7, KINGDOM=8, SUPERKINGDOM=9, DOMAIN=10, LIFE=11;
	
	/** TaxID of Life node */
	public static final int LIFE_ID=1;
	/** TaxID of Cellular Organisms node */
	public static final int CELLULAR_ORGANISMS_ID=131567;
	/** TaxID of Bacteria node */
	public static final int BACTERIA_ID=2; //Is this safe?  Who knows...
	/** TaxID of Archaea node */
	public static final int ARCHAEA_ID=2157;
	/** TaxID of Euk node */
	public static final int EUKARYOTA_ID=2759;
	/** TaxID of Animal node */
	public static final int METAZOA_ID=33208, ANIMALIA_ID=33208;
	/** TaxID of Plant node */
	public static final int VIRIDIPLANTAE_ID=33090, PLANTAE_ID=33090;
	/** TaxID of Fungi node */
	public static final int FUNGI_ID=4751;
	/** TaxID of Virus node */
	public static final int VIRUSES_ID=10239;
	/** TaxID of Viroids node (now defunct) */
	public static final int VIROIDS_ID=12884;

	/** Maps normal level names to normal level numbers */
	private static final HashMap<String, Integer> levelMap=makeLevelMap();
	/** Maps short level names to normal level numbers */
	private static final HashMap<String, Integer> shortLevelMap=makeShortLevelMap();
	/** Maps extended level names to extended level numbers */
	private static final HashMap<String, Integer> levelMapExtended=makeLevelMapExtended();
	/** Maps extended level names to normal level numbers */
	private static final HashMap<String, Integer> altLevelMap=makeAltLevelMap();
	
	/** Common extended level numbers as constants */
	public static final int NO_RANK_E=NO_RANK,
			SUBSTRAIN_E=stringToLevelExtended("substrain"), STRAIN_E=stringToLevelExtended("strain"),
			SUBSPECIES_E=stringToLevelExtended("subspecies"),
			SPECIES_E=stringToLevelExtended("species"), GENUS_E=stringToLevelExtended("genus"),
			FAMILY_E=stringToLevelExtended("family"), ORDER_E=stringToLevelExtended("order"),
			CLASS_E=stringToLevelExtended("class"), PHYLUM_E=stringToLevelExtended("phylum"),
			KINGDOM_E=stringToLevelExtended("kingdom"), SUPERKINGDOM_E=stringToLevelExtended("superkingdom"),
			DOMAIN_E=stringToLevelExtended("domain"), 
//			CELLULAR_ROOT_E=stringToLevelExtended("cellular root"),
//			ACELLULAR_ROOT_E=stringToLevelExtended("acellular root"),
			LIFE_E=stringToLevelExtended("life");
	
	/** Map of normal to extended level numbers */
	private static final int[] levelToExtended=new int[]{
			NO_RANK_E, SUBSPECIES_E, SPECIES_E, GENUS_E, FAMILY_E,
			ORDER_E, CLASS_E, PHYLUM_E, KINGDOM_E, SUPERKINGDOM_E, DOMAIN_E, LIFE_E
		};
	
	/** Map of extended to normal level numbers */
	private static final int[] extendedToLevel=makeExtendedToLevel();
	
	/** Creates one normal-group index per extended rank during class initialization.
	 * @return Newly allocated extended-to-normal mapping */
	private static int[] makeExtendedToLevel(){
		int len=0;
		for(String[] array : taxLevelNamesExtendedMatrix){
			len+=array.length;
		}
		int[] ret=new int[len];
		
		int pos=0;
		for(int level=0; level<taxLevelNamesExtendedMatrix.length; level++){
			String[] array=taxLevelNamesExtendedMatrix[level];
			for(int i=0; i<array.length; i++){
				ret[pos]=level;
				pos++;
			}
		}
		return ret;
	}
	
	/** Converts a normal rank to its corresponding extended rank.
	 * @param level Normal rank index; negative values are returned unchanged
	 * @return Extended rank, or the supplied negative sentinel */
	public static final int levelToExtended(int level){
		return level<0 ? level : levelToExtended[level];
	}

	/** Converts an extended rank to its normal rank group.
	 * Intermediate ranks map to their containing normal group, for example
	 * subphylum to phylum.
	 * @param extended Extended rank index; any negative value is accepted
	 * @return Normal rank, or -1 for a negative input */
	public static final int extendedToLevel(int extended){
		return extended<0 ? -1 : extendedToLevel[extended];
	}
	
	/** Tab delimiter for legacy text records. */
	private static final Pattern delimiterTab = Pattern.compile("\t");
	/** Tab-pipe-tab separator used in NCBI dump records. */
	private static final Pattern delimiter = Pattern.compile("\t\\|\t");
	/** Literal pipe separator for legacy sequence headers. */
	private static final Pattern delimiterPipe = Pattern.compile("\\|");
	/** Literal tilde separator for alternate legacy headers. */
	private static final Pattern delimiterTilde = Pattern.compile("\\~");
	/** One or more whitespace or underscore characters used in name matching. */
	private static final Pattern delimiter2 = Pattern.compile("[\\s_]+");
	
	/** Restrict IMG dump loading to records marked high quality when requested by callers. */
	public static boolean IMG_HQ=false;
	/** Configured main-tree path; auto requests the default resource lookup. */
	public static String treeFile="auto";
	
	/** Default taxonomy directory on NERSC systems. */
	private static final String defaultTaxPathNersc="/global/cfs/cdirs/bbtools/tax/latest";
	/** Default taxonomy directory on the AWS installation. */
	private static final String defaultTaxPathAws="/test1/tax/latest";
	/** Default taxonomy directory on the IGB virtual-machine installation. */
	private static final String defaultTaxPathIGBVM="/data/tax/latest";
	/** Default taxonomy directory on Dori. */
	private static final String defaultTaxPathDori="/clusterfs/jgi/groups/gentech/homes/bbushnell/tax/latest";
	/** NERSC prokaryotic 16S consensus FASTA default. */
	private static final String default16SFileNersc="/global/cfs/cdirs/bbtools/silva/16S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	/** Dori prokaryotic 16S consensus FASTA default. */
	private static final String default16SFileDori="/clusterfs/jgi/groups/gentech/homes/bbushnell/ribo/all_prok_16S_consensus_species_maxns10_taxsorted.fa.gz";
	/** AWS prokaryotic 16S consensus FASTA default. */
	private static final String default16SFileAws="/test1/16S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	/** IGB virtual-machine prokaryotic 16S consensus FASTA default. */
	private static final String default16SFileIGBVM="/data/sketch/silva/16S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	/** NERSC eukaryotic 18S consensus FASTA default. */
	private static final String default18SFileNersc="/global/cfs/cdirs/bbtools/silva/18S_consensus_silva_maxns10_taxsorted.fa.gz";
	/** Dori eukaryotic 18S consensus FASTA default. */
	private static final String default18SFileDori="/clusterfs/jgi/groups/gentech/homes/bbushnell/ribo/all_euk_18S_consensus_species_maxns10_taxsorted.fa.gz";
	/** AWS eukaryotic 18S consensus FASTA default. */
	private static final String default18SFileAws="/test1/18S_consensus_silva_maxns10_taxsorted.fa.gz";
	/** IGB virtual-machine eukaryotic 18S consensus FASTA default. */
	private static final String default18SFileIGBVM="/data/sketch/silva/18S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	
	/** IMG dump path template; TAX_PATH is replaced by the configured taxonomy directory. */
	private static final String defaultImgFile="TAX_PATH/imgDump.txt";
	/** Legacy one-dimensional GI table path template. */
	private static final String defaultTableFileInt="TAX_PATH/gitable.int1d.gz";
	/** Two-dimensional GI table path template. */
	private static final String defaultTableFile="TAX_PATH/gitable.int2d.gz";
	/** Serialized taxonomy tree fallback path template. */
	private static final String defaultTreeFile="TAX_PATH/tree.taxtree.gz";
	/** Accession digitization-pattern path template. */
	private static final String defaultPatternFile="TAX_PATH/patterns.txt";
	/** Per-taxID sequence-size table path template. */
	private static final String defaultSizeFile="TAX_PATH/taxsize.tsv.gz";

	/** Comma-separated accession mapping path templates, each rooted at TAX_PATH. */
	private static final String defaultAccessionFile=
		"TAX_PATH/shrunk.assembly_summary_genbank.txt.gz"
		+ ",TAX_PATH/shrunk.assembly_summary_refseq.txt.gz"
		+ ",TAX_PATH/shrunk.dead_nucl.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.dead_prot.accession2taxid.FULL.gz"
		+ ",TAX_PATH/shrunk.dead_wgs.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.nuclE_wgs.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.nucl_gb.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.nucl_wgs.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.pdb.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.prot.accession2taxid.gz"
		+ ",TAX_PATH/shrunk.wgs.accession2taxid.gz";


	/** Selects a taxonomy root from host flags, falling back to bundled resources.
	 * Priority is AWS outside NERSC, IGB VM, Dori, then Perlmutter/NERSC.
	 * @return Default directory path; existence is not checked */
	public static final String defaultTaxPath(){
		String s=(Shared.AWS && !Shared.NERSC) ? 
				defaultTaxPathAws : 
					Shared.IGBVM ? defaultTaxPathIGBVM : 
						Shared.DORI ? defaultTaxPathDori : 
							(Shared.PERLMUTTER || Shared.NERSC) ? defaultTaxPathNersc :
								Data.RESOURCES();
		return s;
	}

	/** Selects the host-specific prokaryotic 16S consensus path.
	 * @return Dori, AWS, IGB VM or fallback NERSC path; existence is not checked */
	public static final String default16SFile(){
		return Shared.DORI ? default16SFileDori : (Shared.AWS && !Shared.NERSC)
			? default16SFileAws : Shared.IGBVM ? default16SFileIGBVM : default16SFileNersc;
	}

	/** Selects the host-specific eukaryotic 18S consensus path.
	 * @return Dori, AWS, IGB VM or fallback NERSC path; existence is not checked */
	public static final String default18SFile(){
		return Shared.DORI ? default18SFileDori : (Shared.AWS && !Shared.NERSC)
			? default18SFileAws : Shared.IGBVM ? default18SFileIGBVM : default18SFileNersc;
	}

	/** Path to all taxonomy files, substituted in to make specific file paths */
	public static String TAX_PATH=defaultTaxPath();
	
	/** Resolves the two-dimensional GI lookup table template.
	 * @return GI table path under the configured TAX_PATH; existence is not checked */
	public static final String defaultTableFile(){return defaultTableFile.replaceAll("TAX_PATH", TAX_PATH);}
	
	/** Location of taxtree.tsv.gz.
	 * The tree now ships with BBTools, so resources/ wins; TAX_PATH is only a legacy fallback for
	 * installations that predate the bundled tree.  The other taxonomy files (gitable, accession2taxid,
	 * imgDump, ...) are far too large to bundle and still live under TAX_PATH -- this priority applies
	 * to the TREE ALONE.
	 * @return First existing bundled text-tree path, otherwise the configured serialized
	 * tree fallback path even when that fallback does not exist */
	public static final String defaultTreeFile(){
		String s2=Data.RESOURCES()+"taxtree.tsv.gz";
		if(new File(s2).exists()){return s2;}
		s2=Data.ROOT()+"taxtree.tsv.gz";
		if(new File(s2).exists()){return s2;}
		String s=defaultTreeFile.replaceAll("TAX_PATH", TAX_PATH);
		if(new File(s).exists()){return s;}
		return s;
	}
	
	/** Whether to prepend the full-protein accession table to the default input list. */
	public static boolean protFull=false;
	
	/** Resolves the configured accession table list, optionally prepending full proteins.
	 * @return Comma-separated paths under TAX_PATH; existence is not checked */
	public static final String defaultAccessionFile(){
		String s=(protFull ? "TAX_PATH/shrunk.protF.accession2taxid.gz," : "")+defaultAccessionFile;
		return s.replaceAll("TAX_PATH", TAX_PATH);
	}
	/** Resolves the observed accession-pattern table path.
	 * @return patterns.txt under TAX_PATH; existence is not checked */
	public static final String defaultPatternFile(){return defaultPatternFile.replaceAll("TAX_PATH", TAX_PATH);}
	/** Resolves the IMG-to-NCBI taxonomy dump path.
	 * @return imgDump.txt under TAX_PATH; existence is not checked */
	public static final String defaultImgFile(){return defaultImgFile.replaceAll("TAX_PATH", TAX_PATH);}
	/** Resolves the table of sequence amounts associated with each TaxID.
	 * @return taxsize.tsv.gz under TAX_PATH; existence is not checked */
	public static final String defaultSizeFile(){return defaultSizeFile.replaceAll("TAX_PATH", TAX_PATH);}
	
	/** Screen output gets printed here */
	private static PrintStream outstream=System.out;
	
}
