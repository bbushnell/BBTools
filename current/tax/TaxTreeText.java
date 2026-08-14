package tax;

import java.io.PrintStream;
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

/** 
 * Converts serialized TaxTree files to a compact, versioned A48 text format.
 * Loads TaxTree objects from text format.
 * @author Brian Bushnell, Barbara
 * @date August 13, 2026
 */
public final class TaxTreeText {

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
		TaxTree oldTree=ReadWrite.read(TaxTree.class, oldFile, true);
		t.stop();
		outstream.println("Serialized load:\t"+t);

		t.start();
		write(oldTree, textFile, parser.overwrite);
		t.stop();
		outstream.println("Text write:\t"+t);

		t.start();
		TaxTree textTree=load(textFile);
		t.stop();
		outstream.println("Text load:\t"+t);

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

	/** Load a TaxTree from a TaxTree2 text file. */
	public static TaxTree load(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, true, true);
		ByteFile bf=ByteFile.makeByteFile(ff);
		LineParser1 lp=new LineParser1('\t');
		Header h=new Header();
		TaxNode[] nodes=null;
		IntHashMap mergedMap=null;
		boolean readingMerged=false;
		long lineNumber=0;
		int nodesRead=0, mergedRead=0;
		boolean error=false;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				if(line[0]=='#'){
					if(line.length>1 && line[1]=='#'){
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
						}else{
							throw malformed(lineNumber, "Unknown header");
						}
						continue;
					}
					readingMerged=equals(line, "#MergedMap");
					continue;
				}

				lp.set(line);
				if(readingMerged){
					requireTerms(lp, 2, lineNumber);
					if(h.mergedPresent<0){throw malformed(lineNumber, "Merged-map data precedes its header");}
					if(h.mergedPresent==0){throw malformed(lineNumber, "Merged-map data present for a null map");}
					if(mergedMap==null){mergedMap=makeMergedMap(h.mergedCount);}
					int oldID=parseInt(lp, 0, lineNumber, "old TaxID");
					int newID=parseInt(lp, 1, lineNumber, "new TaxID");
					if(mergedMap.contains(oldID)){throw malformed(lineNumber, "Duplicate merged TaxID "+oldID);}
					mergedMap.put(oldID, newID);
					mergedRead++;
				}else{
					if(h.nodesLength<0){throw malformed(lineNumber, "Node data precedes the tree header");}
					if(nodes==null){nodes=new TaxNode[h.nodesLength];}
					requireTerms(lp, NODE_TERMS, lineNumber);
					int levelExtended=parseInt(lp, 0, lineNumber, "extended level");
					if(levelExtended>=TaxTree.numTaxLevelNamesExtended){
						throw malformed(lineNumber, "Invalid extended level: "+levelExtended);
					}
					int level=TaxTree.extendedToLevel(levelExtended);
					int id=parseInt(lp, 1, lineNumber, "TaxID");
					if(id<0 || id>=nodes.length){throw malformed(lineNumber, "TaxID outside node array: "+id);}
					if(nodes[id]!=null){throw malformed(lineNumber, "Duplicate TaxID "+id);}
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
					nodes[id]=n;
					nodesRead++;
				}
			}
		}finally{
			error|=bf.close();
		}
		if(error){throw new RuntimeException("Error while reading "+fname);}
		h.validate(lineNumber);
		if(nodes==null){nodes=new TaxNode[h.nodesLength];}
		if(nodesRead!=h.nodeCount){
			throw malformed(lineNumber, "Expected "+h.nodeCount+" nodes but read "+nodesRead);
		}
		if(h.mergedPresent==1 && mergedMap==null){mergedMap=makeMergedMap(h.mergedCount);}
		if(mergedRead!=h.mergedCount){
			throw malformed(lineNumber, "Expected "+h.mergedCount+" merged IDs but read "+mergedRead);
		}
		TaxTree tree=new TaxTree(nodes, mergedMap, h.minValidTaxa, h.simplify==1,
				h.reassign==1, h.skipNorank==1, h.inferRankLimit);
		if(tree.nodeCount!=h.nodeCount){
			throw malformed(lineNumber, "Reconstructed node count is "+tree.nodeCount+", expected "+h.nodeCount);
		}
		return tree;
	}

	/** Write a TaxTree in TaxTree2 text format. */
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

	/** Validate every persisted and reconstructed field. */
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

	private static ByteBuilder submit(ByteBuilder bb, ByteStreamWriter bsw){
		bsw.print(bb, true);
		return new ByteBuilder(BUFFER_SIZE);
	}

	private static void appendTitle(ByteBuilder bb, String s){
		if(s.length()==0){return;}
		bb.append(Character.toUpperCase(s.charAt(0)));
		bb.append(s, 1, s.length());
	}

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

	private static int parseInt(LineParser1 lp, int term, long lineNumber, String field){
		long x=lp.parseLongA48(term);
		if(x<0 || x>Integer.MAX_VALUE){throw malformed(lineNumber, "Invalid "+field+": "+x);}
		return (int)x;
	}

	private static String parseName(LineParser1 lp, int term, long lineNumber){
		lp.setBounds(term);
		byte[] line=lp.line();
		for(int i=lp.a(); i<lp.b(); i++){
			if(line[i]<0){throw malformed(lineNumber, "Non-ASCII taxonomic name");}
		}
		return lp.parseStringFromCurrentField();
	}

	private static int parseBoolean(LineParser1 lp, int term, long lineNumber, String field){
		int x=parseInt(lp, term, lineNumber, field);
		if(x>1){throw malformed(lineNumber, "Invalid "+field+" boolean: "+x);}
		return x;
	}

	private static IntHashMap makeMergedMap(int size){
		long capacity=Math.max(2L, size*3L/2L+2L);
		if(capacity>Integer.MAX_VALUE){throw new IllegalArgumentException("Merged map is too large: "+size);}
		return new IntHashMap((int)capacity);
	}

	private static void requireTerms(LineParser1 lp, int expected, long lineNumber){
		if(lp.terms()!=expected){throw malformed(lineNumber, "Expected "+expected+" columns but found "+lp.terms());}
	}

	private static boolean equals(byte[] line, String s){
		if(line.length!=s.length()){return false;}
		for(int i=0; i<line.length; i++){if(line[i]!=(byte)s.charAt(i)){return false;}}
		return true;
	}

	private static void check(boolean condition, String field, Object a, Object b){
		if(!condition){throw new IllegalStateException(field+" differs: "+a+" != "+b);}
	}

	private static IllegalStateException nodeMismatch(int id, String field, Object a, Object b){
		return new IllegalStateException("TaxID "+id+" "+field+" differs: "+a+" != "+b);
	}

	private static IllegalArgumentException malformed(long lineNumber, String message){
		return new IllegalArgumentException("Malformed TaxTree2 at line "+lineNumber+": "+message);
	}

	private static String defaultOutputName(String oldFile){
		String suffix="tree.taxtree.gz";
		if(oldFile.endsWith(suffix)){return oldFile.substring(0, oldFile.length()-suffix.length())+"taxtree.tsv.gz";}
		return oldFile+"taxtree.tsv.gz";
	}

	private static final class Header {
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

		int version=-1, nodesLength=-1, nodeCount=-1;
		int minValidTaxa=-1, simplify=-1, reassign=-1, skipNorank=-1, inferRankLimit=-1;
		int mergedPresent=-1, mergedCount=-1;
		boolean encoding=false, columns=false;
	}

	private static final int VERSION=1;
	private static final int NODE_TERMS=8;
	private static final int BUFFER_SIZE=65536;
	private static PrintStream outstream=System.err;
}
