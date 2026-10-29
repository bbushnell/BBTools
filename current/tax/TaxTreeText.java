package tax;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.ListNum;
import template.ThreadWaiter;

/** Reads and writes versioned TaxTree2 text with A48-encoded numeric fields.
 * Text persists native node topology, rank metadata, flags and merged IDs; it
 * excludes auxiliary packed IDs, optional runtime indexes and nonzero counters.
 * Loading reconstructs rank indexes and preserves node order within each rank.
 * This class validates stored fields/counts, not every parent-chain invariant.
 * @author Brian Bushnell, Barbara
 * @date August 13, 2026
 */
public final class TaxTreeText{

	/** Creates a stateless holder; all conversion methods are static. */
	public TaxTreeText(){}

	/** Converts a tree, reloads its output and compares persisted/reconstructed fields.
	 * Defaults to the configured tree resource and an output name derived from it.
	 * @param args Input/output paths or BBTools key=value options */
	public static void main(String[] args){
		PreParser pp=new PreParser(args, TaxTreeText.class, false);
		args=pp.args;
		outstream=pp.outstream;

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		String oldFile=null, textFile=null;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("in") || a.equals("old") || a.equals("tree")){oldFile=b;}
			else if(a.equals("out") || a.equals("new") || a.equals("text")){textFile=b;}
			else if(parser.parse(arg, a, b)){}
			else if(split.length==1 && oldFile==null){oldFile=arg;}
			else if(split.length==1 && textFile==null){textFile=arg;}
			else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		if(oldFile==null){oldFile=TaxTree.defaultTreeFile();}
		if(textFile==null){textFile=defaultOutputName(oldFile);}
		if(!Tools.testInputFiles(false, true, oldFile)){
			throw new RuntimeException("Cannot read input file "+oldFile);
		}
		if(!Tools.testOutputFiles(parser.overwrite, false, false, textFile)){
			throw new RuntimeException("Cannot write output file "+textFile+"; overwrite="+parser.overwrite);
		}
		if(!Tools.testForDuplicateFiles(true, oldFile, textFile)){
			throw new RuntimeException("Input and output files must be different.");
		}

		Timer total=new Timer();
		Timer t=new Timer();
		TaxTree oldTree=TaxTree.loadTaxTreeFile(oldFile);
		t.stop();
		outstream.println("Input load:\t"+t);

		t.start();
		TaxTree.writeTaxTree(oldTree, textFile, parser.overwrite);
		t.stop();
		outstream.println("Output write:\t"+t);

		t.start();
		TaxTree textTree=TaxTree.loadTaxTreeFile(textFile);
		t.stop();
		outstream.println("Output load:\t"+t);

		t.start();
		validate(oldTree, textTree);
		t.stop();
		outstream.println("Validation:\t"+t);
		total.stop();
		int mergedCount=textTree.mergedMap==null ? 0 : textTree.mergedMap.size();
		outstream.println("Validated "+textTree.nodeCount+" nodes and "+mergedCount+" merged TaxIDs.");
		outstream.println("Total time:\t"+total);
		Shared.closeStream(outstream);
	}

	/** Loads TaxTree2 text with one to four parser workers based on available threads.
	 * @param fname Input text path, optionally compressed
	 * @return Reconstructed tree with fresh runtime counters/index state */
	public static TaxTree load(String fname){
		return load(fname, Tools.mid(1, 4, Shared.threads()/2));
	}

	/** Loads text using one to eight parser workers and joins before returning.
	 * Header parsing starts serially; later batches are claimed by workers. Count,
	 * duplicate-ID and rank-position checks follow merging. Reader-close errors and
	 * worker failures are propagated. No global shared-tree cache is changed.
	 * @param fname Input text path, optionally compressed
	 * @param threads Requested parser workers, clamped to the supported range
	 * @return Reconstructed tree */
	public static TaxTree load(String fname, int threads){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, true, true);
		ByteFile bf=ByteFile.makeByteFile(ff);
		BatchSource source=new BatchSource(bf);
		Header h=new Header();
		final int threadCount=Tools.mid(1, threads, MAX_THREADS);
		ParseResult initial=new ParseResult(threadCount);
		ArrayList<ParseThread> workers=null;
		boolean error=false;
		try{
			boolean dataStarted=false;
			while(!dataStarted){
				TaxBatch batch=source.nextBatch();
				if(batch==null){break;}
				dataStarted=parseInitialBatch(batch, h, initial);
			}
			if(!dataStarted){h.validate(source.recordsRead());}
			if(dataStarted){
				workers=new ArrayList<ParseThread>(threadCount);
				for(int i=0; i<threadCount; i++){
					workers.add(new ParseThread(source, h, threadCount, i));
				}
				ThreadWaiter.startAndWait(workers);
			}
		}finally{
			error|=bf.close();
		}
		if(workers!=null){
			for(ParseThread pt : workers){
				if(pt.failure instanceof Error){throw (Error)pt.failure;}
				if(pt.failure instanceof RuntimeException){throw (RuntimeException)pt.failure;}
				if(pt.failure!=null){throw new RuntimeException(pt.failure);}
				if(!pt.success){throw new RuntimeException("TaxTree parser thread "+pt.tid+" terminated abnormally.");}
			}
		}
		if(error){throw new RuntimeException("Error while reading "+fname);}

		TaxNode[] nodes=new TaxNode[h.nodesLength];
		IntHashMap mergedMap=h.mergedPresent==1 ? makeMergedMap(h.mergedCount) : null;
		int nodesRead=0, mergedRead=0;
		ArrayList<ParseResult> results=new ArrayList<ParseResult>(workers==null ? 1 : workers.size()+1);
		results.add(initial);
		if(workers!=null){for(ParseThread pt : workers){results.add(pt.result);}}
		for(ParseResult result : results){
			nodesRead+=result.nodesRead;
			mergedRead+=result.mergedRead;
		}
		if(threadCount<2){
			for(ParseResult result : results){
				for(TaxNode n : result.nodes){
					if(nodes[n.id]!=null){throw malformed("Duplicate TaxID "+n.id+" detected during merge");}
					nodes[n.id]=n;
				}
			}
		}else{
			final int mergeThreads=Tools.min(threadCount, results.size());
			ArrayList<NodeMergeThread> mergers=new ArrayList<NodeMergeThread>(mergeThreads);
			for(int tid=0; tid<mergeThreads; tid++){
				mergers.add(new NodeMergeThread(nodes, results, tid, mergeThreads));
			}
			ThreadWaiter.startAndWait(mergers);
			for(NodeMergeThread mt : mergers){
				if(mt.failure instanceof Error){throw (Error)mt.failure;}
				if(mt.failure instanceof RuntimeException){throw (RuntimeException)mt.failure;}
				if(mt.failure!=null){throw new RuntimeException(mt.failure);}
				if(!mt.success){throw new RuntimeException("TaxTree node-merge thread terminated abnormally.");}
			}
		}
		for(ParseResult result : results){
			if(result.mergedMap!=null){
				int[] keys=result.mergedMap.keys(), values=result.mergedMap.values();
				int invalid=result.mergedMap.invalid();
				for(int i=0; i<keys.length; i++){
					int oldID=keys[i];
					if(oldID!=invalid){
						if(mergedMap.contains(oldID)){
							throw malformed("Duplicate merged TaxID "+oldID+" detected during parallel merge");
						}
						mergedMap.put(oldID, values[i]);
					}
				}
			}
		}
		if(nodesRead!=h.nodeCount){
			throw malformed(source.recordsRead(), "Expected "+h.nodeCount+" nodes but read "+nodesRead);
		}
		if(mergedRead!=h.mergedCount || (mergedMap!=null && mergedMap.size()!=h.mergedCount)){
			throw malformed(source.recordsRead(), "Expected "+h.mergedCount+" merged IDs but read "+mergedRead);
		}
		TaxTree tree=new TaxTree(nodes, mergedMap, h.minValidTaxa, h.simplify==1,
				h.reassign==1, h.skipNorank==1, h.inferRankLimit, threadCount);
		if(tree.nodeCount!=h.nodeCount){
			throw malformed(source.recordsRead(), "Reconstructed node count is "+tree.nodeCount+", expected "+h.nodeCount);
		}
		return tree;
	}

	/** Reads header records and any first data records in the same batch.
	 * @param batch Numbered input batch
	 * @param h Mutable header accumulated before the first data row
	 * @param result Initial serial parse output
	 * @return Whether this batch contained a data row */
	private static boolean parseInitialBatch(TaxBatch batch, Header h, ParseResult result){
		LineParser1 lp=new LineParser1('\t');
		boolean dataStarted=false;
		for(int i=0; i<batch.lines.size(); i++){
			byte[] line=batch.lines.get(i);
			long lineNumber=batch.lines.firstRecordNum+i+1;
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(line.length>1 && line[1]=='#'){
					if(dataStarted){throw malformed(lineNumber, "Header appears after data");}
					parseHeader(line, h, lp, lineNumber);
				}
				continue;
			}
			if(!dataStarted){h.validate(lineNumber); dataStarted=true;}
			parseDataLine(line, h, result, lp, lineNumber);
		}
		return dataStarted;
	}

	/** Parses a post-header batch; blank lines and single-# comments are skipped.
	 * @param batch Numbered input batch
	 * @param h Validated read-only header
	 * @param result Worker-local output
	 * @param lp Reused worker-local tab parser */
	private static void parseBatch(TaxBatch batch, Header h, ParseResult result, LineParser1 lp){
		for(int i=0; i<batch.lines.size(); i++){
			byte[] line=batch.lines.get(i);
			long lineNumber=batch.lines.firstRecordNum+i+1;
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(line.length>1 && line[1]=='#'){throw malformed(lineNumber, "Header appears after data");}
				continue;
			}
			parseDataLine(line, h, result, lp, lineNumber);
		}
	}

	/** Decodes one recognized double-# metadata row into the header.
	 * @param line Raw header bytes
	 * @param h Mutable header state
	 * @param lp Reused tab parser
	 * @param lineNumber One-based input record number for diagnostics */
	private static void parseHeader(byte[] line, Header h, LineParser1 lp, long lineNumber){
		lp.set(line);
		if(lp.termEquals("##TaxTree2", 0)){
			requireTerms(lp, 2, lineNumber);
			h.version=parseInt(lp, 1, lineNumber, "version");
		}else if(lp.termEquals("##Encoding", 0)){
			requireTerms(lp, 2, lineNumber);
			if(!lp.termEquals("A48", 1)){throw malformed(lineNumber, "Unsupported numeric encoding");}
			h.encoding=true;
		}else if(lp.termEquals("##Tree", 0)){
			requireTerms(lp, 3, lineNumber);
			h.nodesLength=parseInt(lp, 1, lineNumber, "nodes length");
			h.nodeCount=parseInt(lp, 2, lineNumber, "node count");
		}else if(lp.termEquals("##Config", 0)){
			requireTerms(lp, 6, lineNumber);
			h.minValidTaxa=parseInt(lp, 1, lineNumber, "minimum valid taxa");
			h.simplify=parseBoolean(lp, 2, lineNumber, "simplify");
			h.reassign=parseBoolean(lp, 3, lineNumber, "reassign");
			h.skipNorank=parseBoolean(lp, 4, lineNumber, "skip no-rank");
			h.inferRankLimit=parseInt(lp, 5, lineNumber, "rank inference limit");
		}else if(lp.termEquals("##MergedMap", 0)){
			requireTerms(lp, 3, lineNumber);
			h.mergedPresent=parseBoolean(lp, 1, lineNumber, "merged-map presence");
			h.mergedCount=parseInt(lp, 2, lineNumber, "merged-map size");
		}else if(lp.termEquals("##LevelExtended", 0)){
			requireTerms(lp, NODE_TERMS, lineNumber);
			h.columns=true;
		}else{throw malformed(lineNumber, "Unknown header");}
	}

	/** Parses either a two-column merged-ID row or an eight-column node row.
	 * @param line Raw data bytes
	 * @param h Validated header supplying array/count metadata
	 * @param result Worker-local accumulated nodes and mappings
	 * @param lp Reused tab parser
	 * @param lineNumber One-based input record number for diagnostics */
	private static void parseDataLine(byte[] line, Header h,
			ParseResult result, LineParser1 lp, long lineNumber){
		lp.set(line);
		final int terms=lp.terms();
		if(terms==2){
			if(h.mergedPresent<0){throw malformed(lineNumber, "Merged-map data precedes its header");}
			if(h.mergedPresent==0){throw malformed(lineNumber, "Merged-map data present for a null map");}
			int oldID=parseInt(lp, 0, lineNumber, "old TaxID");
			int newID=parseInt(lp, 1, lineNumber, "new TaxID");
			IntHashMap map=result.mergedMap(h.mergedCount);
			if(map.contains(oldID)){throw malformed(lineNumber, "Duplicate merged TaxID "+oldID);}
			map.put(oldID, newID);
			result.mergedRead++;
			return;
		}

		if(h.nodesLength<0){throw malformed(lineNumber, "Node data precedes the tree header");}
		if(terms!=NODE_TERMS){throw malformed(lineNumber, "Expected "+NODE_TERMS+" columns but found "+terms);}
		int levelExtended=parseInt(lp, 0, lineNumber, "extended level");
		if(levelExtended>=TaxTree.numTaxLevelNamesExtended){
			throw malformed(lineNumber, "Invalid extended level: "+levelExtended);
		}
		int level=TaxTree.extendedToLevel(levelExtended);
		int id=parseInt(lp, 1, lineNumber, "TaxID");
		if(id<0 || id>=h.nodesLength){throw malformed(lineNumber, "TaxID outside node array: "+id);}
		int pid=parseInt(lp, 2, lineNumber, "parent TaxID");
		int numChildren=parseInt(lp, 3, lineNumber, "number of children");
		int minParentLevelExtended=parseInt(lp, 4, lineNumber, "minimum parent level");
		int maxChildLevelExtended=parseInt(lp, 5, lineNumber, "maximum child level");
		long flag=lp.parseLongA48(6);
		String name=parseName(lp, 7, lineNumber);
		TaxNode n=new TaxNode(id, pid, level, levelExtended, name);
		n.numChildren=numChildren;
		n.minParentLevelExtended=minParentLevelExtended;
		n.maxChildLevelExtended=maxChildLevelExtended;
		n.setRawFlag(flag);
		result.nodes.add(n);
		result.nodesRead++;
	}

	/** Writes text and waits for the output writer to finish.
	 * Rejects unsupported runtime state rather than silently discarding it. Numeric
	 * fields use A48; names must be single-line ASCII without tabs.
	 * @param tree Native-ID tree to serialize
	 * @param fname Output path, optionally compressed
	 * @param overwrite Permit replacing existing output */
	public static void write(TaxTree tree, String fname, boolean overwrite){
		if(tree==null){throw new NullPointerException("tree");}
		if(tree.hasTextExcludedState()){
			throw new IllegalArgumentException("TaxTree2 stores tree structure, not populated runtime name or statistics maps.");
		}
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ff);
		boolean error;
		try{
			writeInner(tree, bsw);
		}finally{
			error=bsw.poisonAndWait();
		}
		if(error){throw new RuntimeException("Error while writing "+fname);}
	}

	/** Writes metadata, nodes by descending rank, and mappings by ascending old ID.
	 * @param tree Tree whose runtime state passed the outer guard
	 * @param bsw Started writer owned and closed by the caller */
	private static void writeInner(TaxTree tree, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder(BUFFER_SIZE);
		bb.append("##TaxTree2").tab().appendA48(VERSION).nl();
		bb.append("##Encoding\tA48\n");
		bb.append("##Tree").tab().appendA48(tree.nodes.length).tab().appendA48(tree.nodeCount).nl();
		bb.append("##Config").tab().appendA48(tree.minValidTaxa).tab().appendA48(tree.simplify ? 1 : 0)
				.tab().appendA48(tree.reassign ? 1 : 0).tab().appendA48(tree.skipNorank ? 1 : 0)
				.tab().appendA48(tree.inferRankLimit).nl();
		int mergedCount=tree.mergedMap==null ? 0 : tree.mergedMap.size();
		bb.append("##MergedMap").tab().appendA48(tree.mergedMap==null ? 0 : 1).tab().appendA48(mergedCount).nl();
		bb.append("##LevelExtended\tTID\tPID\tNumChildren\tMinParentLevelExtended")
				.append("\tMaxChildLevelExtended\tFlag\tName\n");

		int written=0;
		for(int level=tree.treeLevelsExtended.length-1; level>=0; level--){
			bb.append('#');
			appendTitle(bb, TaxTree.levelToStringExtended(level));
			bb.nl();
			for(TaxNode n : tree.treeLevelsExtended[level]){
				if(n.level!=TaxTree.extendedToLevel(n.levelExtended)){
					throw new IllegalArgumentException("Normal level is not derivable for TaxID "+n.id);
				}
				if(n.countRaw!=0 || n.countSum!=0){
					throw new IllegalArgumentException("TaxTree2 does not store runtime abundance counts; TaxID "+n.id+" is nonzero.");
				}
				bb.appendA48(n.levelExtended).tab().appendA48(n.id).tab().appendA48(n.pid).tab().appendA48(n.numChildren)
						.tab().appendA48(n.minParentLevelExtended).tab().appendA48(n.maxChildLevelExtended)
						.tab().appendA48(n.rawFlag()).tab();
				appendName(bb, n.name);
				bb.nl();
				written++;
				if(bb.length()>=BUFFER_SIZE){bb=submit(bb, bsw);}
			}
		}
		if(written!=tree.nodeCount){
			throw new IllegalStateException("Expected to write "+tree.nodeCount+" nodes but wrote "+written);
		}
		bb.append("#MergedMap\n");
		if(tree.mergedMap!=null){
			int[] oldIDs=tree.mergedMap.toArray();
			Arrays.sort(oldIDs);
			for(int oldID : oldIDs){
				bb.appendA48(oldID).tab().appendA48(tree.mergedMap.get(oldID)).nl();
				if(bb.length()>=BUFFER_SIZE){bb=submit(bb, bsw);}
			}
		}
		if(bb.length()>0){bsw.print(bb, true);}
	}

	/** Compares supported persisted state, rank indexes and counters without mutation.
	 * @param a Expected tree
	 * @param b Reconstructed tree
	 * @throws IllegalStateException If a compared field differs
	 * @throws IllegalArgumentException If unsupported runtime maps are present */
	public static void validate(TaxTree a, TaxTree b){
		if(a.hasTextExcludedState() || b.hasTextExcludedState()){
			throw new IllegalArgumentException("Cannot validate TaxTree2 with populated runtime name or statistics maps.");
		}
		check(a.nodes.length==b.nodes.length, "nodes.length", a.nodes.length, b.nodes.length);
		check(a.nodeCount==b.nodeCount, "nodeCount", a.nodeCount, b.nodeCount);
		check(a.minValidTaxa==b.minValidTaxa, "minValidTaxa", a.minValidTaxa, b.minValidTaxa);
		check(a.simplify==b.simplify, "simplify", a.simplify, b.simplify);
		check(a.reassign==b.reassign, "reassign", a.reassign, b.reassign);
		check(a.skipNorank==b.skipNorank, "skipNorank", a.skipNorank, b.skipNorank);
		check(a.inferRankLimit==b.inferRankLimit, "inferRankLimit", a.inferRankLimit, b.inferRankLimit);
		check(Arrays.equals(a.nodesPerLevel, b.nodesPerLevel), "nodesPerLevel", Arrays.toString(a.nodesPerLevel), Arrays.toString(b.nodesPerLevel));
		check(Arrays.equals(a.nodesPerLevelExtended, b.nodesPerLevelExtended), "nodesPerLevelExtended",
				Arrays.toString(a.nodesPerLevelExtended), Arrays.toString(b.nodesPerLevelExtended));
		validateMergedMap(a.mergedMap, b.mergedMap);
		for(int i=0; i<a.nodes.length; i++){
			TaxNode x=a.nodes[i], y=b.nodes[i];
			if(x==null || y==null){check(x==y, "nodes["+i+"]", x, y);}
			else{validateNode(x, y);}
		}
		check(a.treeLevelsExtended.length==b.treeLevelsExtended.length, "treeLevelsExtended.length",
				a.treeLevelsExtended.length, b.treeLevelsExtended.length);
		for(int level=0; level<a.treeLevelsExtended.length; level++){
			TaxNode[] x=a.treeLevelsExtended[level], y=b.treeLevelsExtended[level];
			check(x.length==y.length, "treeLevelsExtended["+level+"].length", x.length, y.length);
			for(int i=0; i<x.length; i++){
				check(x[i].id==y[i].id, "treeLevelsExtended["+level+"]["+i+"]", x[i].id, y[i].id);
			}
		}
	}

	/** Compares one pair of non-null nodes, reporting the first differing field.
	 * @param a Expected node
	 * @param b Reconstructed node */
	private static void validateNode(TaxNode a, TaxNode b){
		if(a.id!=b.id){throw nodeMismatch(a.id, "id", a.id, b.id);}
		if(!Objects.equals(a.name, b.name)){throw nodeMismatch(a.id, "name", a.name, b.name);}
		if(a.pid!=b.pid){throw nodeMismatch(a.id, "pid", a.pid, b.pid);}
		if(a.level!=b.level){throw nodeMismatch(a.id, "level", a.level, b.level);}
		if(a.levelExtended!=b.levelExtended){throw nodeMismatch(a.id, "levelExtended", a.levelExtended, b.levelExtended);}
		if(a.numChildren!=b.numChildren){throw nodeMismatch(a.id, "numChildren", a.numChildren, b.numChildren);}
		if(a.minParentLevelExtended!=b.minParentLevelExtended){
			throw nodeMismatch(a.id, "minParentLevelExtended", a.minParentLevelExtended, b.minParentLevelExtended);
		}
		if(a.maxChildLevelExtended!=b.maxChildLevelExtended){
			throw nodeMismatch(a.id, "maxChildLevelExtended", a.maxChildLevelExtended, b.maxChildLevelExtended);
		}
		if(a.rawFlag()!=b.rawFlag()){throw nodeMismatch(a.id, "flag", a.rawFlag(), b.rawFlag());}
		if(a.countRaw!=b.countRaw){throw nodeMismatch(a.id, "countRaw", a.countRaw, b.countRaw);}
		if(a.countSum!=b.countSum){throw nodeMismatch(a.id, "countSum", a.countSum, b.countSum);}
	}

	/** Compares map nullness, size, keys and values without relying on slot order.
	 * @param a Expected old-to-new ID map, possibly null
	 * @param b Reconstructed map, possibly null */
	private static void validateMergedMap(IntHashMap a, IntHashMap b){
		if(a==null || b==null){check(a==b, "mergedMap", a, b); return;}
		check(a.size()==b.size(), "mergedMap.size", a.size(), b.size());
		int[] keys=a.keys();
		int invalid=a.invalid();
		for(int key : keys){
			if(key!=invalid){
				check(b.contains(key), "mergedMap key "+key, true, false);
				check(a.get(key)==b.get(key), "mergedMap["+key+"]", a.get(key), b.get(key));
			}
		}
	}

	/** Transfers the filled buffer to the writer and allocates the next owned buffer.
	 * @param bb Buffer that must not be reused after submission
	 * @param bsw Started asynchronous writer
	 * @return Fresh empty output buffer */
	private static ByteBuilder submit(ByteBuilder bb, ByteStreamWriter bsw){
		bsw.print(bb, true);
		return new ByteBuilder(BUFFER_SIZE);
	}

	/** Appends a rank title with its first character capitalized.
	 * @param bb Destination buffer
	 * @param s Non-null title; empty appends nothing */
	private static void appendTitle(ByteBuilder bb, String s){
		if(s.length()==0){return;}
		bb.append(Character.toUpperCase(s.charAt(0)));
		bb.append(s, 1, s.length());
	}

	/** Validates and appends a non-null single-line ASCII name without tab delimiters.
	 * @param bb Destination buffer
	 * @param s Name to append; empty names are permitted */
	private static void appendName(ByteBuilder bb, String s){
		if(s==null){throw new IllegalArgumentException("Taxonomic names may not be null.");}
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c>127 || c=='\t' || c=='\n' || c=='\r'){
				throw new IllegalArgumentException("TaxTree2 names must be single-line ASCII without tabs: "+s);
			}
		}
		bb.append(s);
	}

	/** Parses an A48 value and requires its decoded result to fit a nonnegative int.
	 * @param lp Parser positioned on the current row
	 * @param term Zero-based field index
	 * @param lineNumber Input record number for diagnostics
	 * @param field Human-readable field description
	 * @return Decoded nonnegative value */
	private static int parseInt(LineParser1 lp, int term, long lineNumber, String field){
		long x=lp.parseLongA48(term);
		if(x<0 || x>Integer.MAX_VALUE){throw malformed(lineNumber, "Invalid "+field+": "+x);}
		return (int)x;
	}

	/** Extracts an ASCII name field after rejecting bytes with their high bit set.
	 * @param lp Parser positioned on the current row
	 * @param term Zero-based name field
	 * @param lineNumber Input record number for diagnostics
	 * @return Parsed name, possibly empty */
	private static String parseName(LineParser1 lp, int term, long lineNumber){
		lp.setBounds(term);
		byte[] line=lp.line();
		for(int i=lp.a(); i<lp.b(); i++){
			if(line[i]<0){throw malformed(lineNumber, "Non-ASCII taxonomic name");}
		}
		return lp.parseStringFromCurrentField();
	}

	/** Parses an A48 boolean restricted to zero or one.
	 * @param lp Parser positioned on the current row
	 * @param term Zero-based field index
	 * @param lineNumber Input record number for diagnostics
	 * @param field Human-readable field description
	 * @return Zero or one */
	private static int parseBoolean(LineParser1 lp, int term, long lineNumber, String field){
		int x=parseInt(lp, term, lineNumber, field);
		if(x>1){throw malformed(lineNumber, "Invalid "+field+" boolean: "+x);}
		return x;
	}

	/** Allocates an old-to-new mapping with approximately 50 percent spare capacity.
	 * @param size Expected number of mappings
	 * @return Empty map */
	private static IntHashMap makeMergedMap(int size){
		long capacity=Math.max(2L, size*3L/2L+2L);
		if(capacity>Integer.MAX_VALUE){throw new IllegalArgumentException("Merged map is too large: "+size);}
		return new IntHashMap((int)capacity);
	}

	/** Requires an exact column count for the current row.
	 * @param lp Current row parser
	 * @param expected Required number of fields
	 * @param lineNumber Input record number for diagnostics */
	private static void requireTerms(LineParser1 lp, int expected, long lineNumber){
		if(lp.terms()!=expected){throw malformed(lineNumber, "Expected "+expected+" columns but found "+lp.terms());}
	}

	/** Throws a diagnostic mismatch when a compared field is unequal.
	 * @param condition Whether the values match
	 * @param field Field description
	 * @param a Expected value
	 * @param b Observed value */
	private static void check(boolean condition, String field, Object a, Object b){
		if(!condition){throw new IllegalStateException(field+" differs: "+a+" != "+b);}
	}

	/** Builds a mismatch exception without throwing it.
	 * @param id Expected TaxID
	 * @param field Field description
	 * @param a Expected value
	 * @param b Observed value
	 * @return Diagnostic exception */
	private static IllegalStateException nodeMismatch(int id, String field, Object a, Object b){
		return new IllegalStateException("TaxID "+id+" "+field+" differs: "+a+" != "+b);
	}

	/** Builds a row-specific input error.
	 * @param lineNumber One-based input record number
	 * @param message Reason the input is rejected
	 * @return Diagnostic exception */
	private static IllegalArgumentException malformed(long lineNumber, String message){
		return new IllegalArgumentException("Malformed TaxTree2 at line "+lineNumber+": "+message);
	}

	/** Builds an input error not tied to one record.
	 * @param message Reason the input is rejected
	 * @return Diagnostic exception */
	private static IllegalArgumentException malformed(String message){
		return new IllegalArgumentException("Malformed TaxTree2: "+message);
	}

	/** Derives a conversion output name, preserving compression suffix handling.
	 * @param oldFile Input path
	 * @return Sibling TSV name or a .copy.tsv.gz name for existing text input */
	private static String defaultOutputName(String oldFile){
		String suffix="tree.taxtree.gz";
		if(oldFile.endsWith(suffix)){return oldFile.substring(0, oldFile.length()-suffix.length())+"taxtree.tsv.gz";}
		if(TaxTree.isTextTreeFile(oldFile)){return ReadWrite.stripExtension(oldFile)+".copy.tsv.gz";}
		return oldFile+"taxtree.tsv.gz";
	}

	/** Serializes only batch acquisition and record numbering. */
	private static final class BatchSource{
		/** Borrows a reader whose close is owned by the outer load operation.
		 * @param bf_ Shared input reader */
		BatchSource(ByteFile bf_){bf=bf_;}

		/** Acquires the next batch and validates contiguous record numbering.
		 * @return Numbered batch, or null at EOF */
		synchronized TaxBatch nextBatch(){
			ListNum<byte[]> lines=bf.nextList();
			if(lines==null){return null;}
			if(lines.firstRecordNum<0){lines.firstRecordNum=recordsRead;}
			else if(lines.firstRecordNum!=recordsRead){
				throw new IllegalStateException("Noncontiguous ByteFile list "+lines.id+": starts at "+
						lines.firstRecordNum+", expected "+recordsRead);
			}
			recordsRead+=lines.size();
			return new TaxBatch(lines);
		}

		/** Returns records acquired so far, including blank/comment records.
		 * @return Record count, read serially or after workers are joined */
		long recordsRead(){return recordsRead;}

		/** Shared input, accessed only while holding this object's monitor. */
		private final ByteFile bf;
		/** Number of records assigned to batches. */
		private long recordsRead=0;
	}

	/** Numbered batch handed to one parser worker. */
	private static final class TaxBatch{
		/** Retains the reader's list and line arrays without copying.
		 * @param lines_ Acquired input list */
		TaxBatch(ListNum<byte[]> lines_){
			lines=lines_;
		}

		/** Owned batch payload; its record numbering includes skipped comments. */
		final ListNum<byte[]> lines;
	}

	/** Parses claimed batches into private lists and maps, capturing any failure. */
	private static final class ParseThread extends Thread{
		/** Configures one parser, with no input acquired until run.
		 * @param source_ Shared synchronized batch source
		 * @param h_ Validated read-only header
		 * @param threads Total parser workers, used to size private maps
		 * @param tid_ Worker identifier used in diagnostics */
		ParseThread(BatchSource source_, Header h_, int threads, int tid_){
			source=source_;
			h=h_;
			tid=tid_;
			result=new ParseResult(threads);
			setName("TaxTreeText-"+tid);
		}

		@Override
		public void run(){
			LineParser1 lp=new LineParser1('\t');
			try{
				for(TaxBatch batch=source.nextBatch(); batch!=null; batch=source.nextBatch()){
					parseBatch(batch, h, result, lp);
				}
				success=true;
			}catch(Throwable t){failure=t;}
		}

		/** Shared source coordinating batch acquisition. */
		final BatchSource source;
		/** Immutable-during-run parsed header. */
		final Header h;
		/** Output owned by this worker until it is joined. */
		final ParseResult result;
		/** Diagnostic worker number. */
		final int tid;
		/** Captured failure, read after joining. */
		Throwable failure=null;
		/** True only after normal EOF. */
		boolean success=false;
	}

	/** Accumulated nodes and mappings from one serial/worker parse stream. */
	private static final class ParseResult{
		/** Stores the worker count used for later map sizing.
		 * @param threads_ Positive total parser count */
		ParseResult(int threads_){threads=threads_;}

		/** Lazily allocates this parser's merged-ID map.
		 * @param totalSize Total mappings expected across all parsers
		 * @return Private mutable map */
		IntHashMap mergedMap(int totalSize){
			if(mergedMap==null){mergedMap=makeMergedMap(Tools.max(2, totalSize/threads+2));}
			return mergedMap;
		}

		/** Parsed nodes in this parser's encounter order. */
		final ArrayList<TaxNode> nodes=new ArrayList<TaxNode>(4096);
		/** Total parser count used for initial map capacity. */
		final int threads;
		/** Private old-to-new ID map, absent until a mapping row is seen. */
		IntHashMap mergedMap=null;
		/** Counts of successfully parsed node and mapping records. */
		int nodesRead=0, mergedRead=0;
	}

	/** Places parsed nodes in their native-ID slots; final counts detect lost duplicates. */
	private static final class NodeMergeThread extends Thread{
		/** Assigns a strided set of parser result lists.
		 * @param nodes_ Shared sparse destination array
		 * @param results_ Joined parser results, read-only
		 * @param tid_ First result-list index for this worker
		 * @param threads_ Stride between result-list indices */
		NodeMergeThread(TaxNode[] nodes_, ArrayList<ParseResult> results_, int tid_, int threads_){
			nodes=nodes_;
			results=results_;
			tid=tid_;
			threads=threads_;
		}

		@Override
		public void run(){
			try{
				for(int i=tid; i<results.size(); i+=threads){
					for(TaxNode n : results.get(i).nodes){nodes[n.id]=n;}
				}
				success=true;
			}catch(Throwable t){failure=t;}
		}

		/** Shared destination; valid input gives every TaxID exactly one writer. */
		final TaxNode[] nodes;
		/** Read-only parser outputs. */
		final ArrayList<ParseResult> results;
		/** First assigned list and stride between assigned lists. */
		final int tid, threads;
		/** Captured failure, read after joining. */
		Throwable failure=null;
		/** True only after every assigned node was stored. */
		boolean success=false;
	}

	/** Mutable metadata accumulated before data parsing; negative fields mean missing. */
	private static final class Header{
		/** Creates an empty header whose required fields are initially absent. */
		Header(){}

		/** Requires a supported version and all structural/configuration headers.
		 * @param lineNumber Current record number for missing/invalid header diagnostics */
		void validate(long lineNumber){
			if(version!=VERSION){throw malformed(lineNumber, "Unsupported version "+version);}
			if(!encoding){throw malformed(lineNumber, "Missing A48 encoding header");}
			if(nodesLength<0 || nodeCount<0){throw malformed(lineNumber, "Missing tree header");}
			if(minValidTaxa<0 || simplify<0 || reassign<0 || skipNorank<0 || inferRankLimit<0){
				throw malformed(lineNumber, "Missing configuration header");
			}
			if(mergedPresent<0 || mergedCount<0){throw malformed(lineNumber, "Missing merged-map header");}
			if(mergedPresent==0 && mergedCount!=0){throw malformed(lineNumber, "Null merged map has nonzero size");}
			if(!columns){throw malformed(lineNumber, "Missing column header");}
		}

		/** Format version, native array length and expected non-null node count. */
		int version=-1, nodesLength=-1, nodeCount=-1;
		/** Persisted configuration values; booleans are encoded as zero or one. */
		int minValidTaxa=-1, simplify=-1, reassign=-1, skipNorank=-1, inferRankLimit=-1;
		/** Whether a merged map exists and its expected record count. */
		int mergedPresent=-1, mergedCount=-1;
		/** Presence of the encoding declaration and eight-column node declaration. */
		boolean encoding=false, columns=false;
	}

	/** Supported TaxTree2 schema version. */
	private static final int VERSION=1;
	/** Number of columns in a node record. */
	private static final int NODE_TERMS=8;
	/* ln004, 2026-08-21: raw 138,747,082-byte tree, forced BF1, 5 warmups plus
	 * 20 loads: 1/2/4/8-thread medians were 0.558/0.313/0.196/0.185 seconds. */
	/** Upper bound on parser workers. */
	private static final int MAX_THREADS=8;
	/** Initial output-buffer capacity and submission threshold in bytes. */
	private static final int BUFFER_SIZE=65536;
	/** CLI status destination; changed by command-line preparsing. */
	private static PrintStream outstream=System.err;
}
