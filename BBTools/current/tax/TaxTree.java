package tax;

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
import shared.LineParser1;
import shared.LineParserS1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.IntLongHashMap;

/**
 * Represents a taxonomic tree.
 * Usually just one of these needs to be created for a process.
 * Designed for NCBI's taxdmp.zip file contents.
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxTree implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 5894416521711540017L;
	
	/*--------------------------------------------------------------*/
	/*----------------              Main            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * This is not called normally, only when converting NCBI text files
	 * into a binary representation and writing it to disk.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, outstream, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		assert(args.length>=4) : "TaxTree syntax:\ntaxtree.sh names.dmp nodes.dmp merged.dmp tree.taxtree.gz\n";
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
					outstream.print("\n"+n+" -> "+tree.nodes[n.pid]);
				}
				for(int j=tree.treeLevelsExtended[i].length-lim; j<tree.treeLevelsExtended[i].length; j++){
					if(j>=lim){
						TaxNode n=tree.treeLevelsExtended[i][j];
						outstream.print("\n"+n+" -> "+tree.nodes[n.pid]);
					}
				}
			}
			outstream.println();
		}
		
		
		outstream.println();
		outstream.println("Time: \t"+t);
		
		if(args.length>2){//Write a tree
			ReadWrite.write(tree, args[3], true);
		}
	}
	
	/** Parse arguments from the command line */
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
	 * @param namesFile NCBI names.txt
	 * @param nodesFile NCBI nodes.txt
	 * @param mergedFile NCBI merged.txt
	 * @param args
	 */
	private TaxTree(String namesFile, String nodesFile, String mergedFile, String[] args){
		
		if(args!=null) {
			Parser parser=parse(args);
		}
		
		nodes=getNames(namesFile);
		getNodes(nodesFile, nodes);
		
		mergedMap=getMerged(mergedFile);
		
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
		if(errors>0) {
			System.err.println("Found "+errors+" errors in tree.");
		}
		
		for(TaxNode n : nodes){
			if(n!=null){
				nodesPerLevel[n.level]++;
				nodesPerLevelExtended[n.levelExtended]++;
			}
		}
		
//		for(int i=0; i<nodesPerLevel.length; i++){
//			treeLevels[i]=new TaxNode[nodesPerLevel[i]];
//		}
		for(int i=0; i<nodesPerLevelExtended.length; i++){
			treeLevelsExtended[i]=new TaxNode[nodesPerLevelExtended[i]];
		}
		
//		{
//			int[] temp=new int[nodesPerLevel.length];
//			for(TaxNode n : nodes){
//				if(n!=null){
//					int level=n.level;
//					treeLevels[level][temp[level]]=n;
//					temp[level]++;
//				}
//			}
//		}
		
		{
			int[] temp=new int[nodesPerLevelExtended.length];
			for(TaxNode n : nodes){
				if(n!=null){
					int level=n.levelExtended;
					treeLevelsExtended[level][temp[level]]=n;
					temp[level]++;
				}
			}
		}
		nodeCount=(int)shared.Vector.sum(nodesPerLevelExtended);
		
	}
	
	/*--------------------------------------------------------------*/
	/*---------------            Loaders            ----------------*/
	/*--------------------------------------------------------------*/


	/**
	 * Load a tax tree from disk.
	 * @param taxTreeFile Serialized TaxTree.
	 * @param hashNames Hash nodes using names as keys
	 * @param hashDotFormat Hash nodes using abbreviations, e.g. H.sapiens
	 * @return
	 */
	public static final TaxTree loadTaxTree(String taxTreeFile, PrintStream outstream, boolean hashNames, 
			boolean hashDotFormat){
		if(taxTreeFile==null){return null;}
		return loadTaxTree(taxTreeFile, null, null, null, outstream, hashNames, hashDotFormat);
	}
	
	/**
	 * Load a tax tree from disk, either from a binary tree file,
	 * or from NCBI text files.
	 * @param taxTreeFile Binary representation; mutually exclusive with other files.
	 * @param taxNameFile NCBI names.txt
	 * @param taxNodeFile NCBI nodes.txt
	 * @param taxMergedFile NCBI merged.txt
	 * @param hashNames Hash nodes using names as keys
	 * @param hashDotFormat Hash nodes using abbreviations, e.g. H.sapiens
	 * @return The loaded tree
	 */
	public static final TaxTree loadTaxTree(String taxTreeFile, String taxNameFile, String taxNodeFile, 
			String taxMergedFile, PrintStream outstream, boolean hashNames, boolean hashDotFormat){
		if(taxTreeFile!=null || taxNodeFile==null){
			TaxTree tree=sharedTree(taxTreeFile, hashNames, hashDotFormat, outstream);
			if(tree!=null){return tree;}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=defaultTreeFile();}
		assert(taxTreeFile!=null || (taxNameFile!=null && taxNodeFile!=null)) : "Must specify both taxname and taxnode files.";
		Timer t=new Timer();
		if(outstream!=null){outstream.print("\nLoading tax tree; ");}
		final TaxTree tree;
		if(taxTreeFile!=null){
			tree=ReadWrite.read(TaxTree.class, taxTreeFile, true);
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
	 * Parses names file a second time to fill in additional information.
	 * Should really be merged into getNames.
	 * @TODO Merge into getNames
	 */
	private static TaxNode[] getNodes(String fname, TaxNode[] nodes){
		
		int max=0;
		
		LinkedHashMap<String, int[]> oddNames=new LinkedHashMap<String, int[]>();
		
		TextFile tf=new TextFile(fname, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			String[] split=delimiter.split(s, 4);
			assert(split.length==4) : s;
			int id=-1, pid=-1, level=-1, levelExtended=-1;
			
			id=Integer.parseInt(split[0]);
			try {
				pid=Integer.parseInt(split[1]);
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println("Bad line: "+s+"\n"+Arrays.toString(split));
			}
			boolean alt=false;
			{
				String key=split[2];
				Integer obj0=levelMap.get(key);
				Integer obj=levelMapExtended.get(key);
				if(obj==null) {
					if(key.equalsIgnoreCase("cellular root") || key.equalsIgnoreCase("acellular root")) {
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
	 * @TODO This could be finished in 2 passes using the childTable.
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
	 * Fetch this node's children.
	 * @param parent Node in question
	 * @return List of child nodes
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
	private static IntHashMap getMerged(String mergedFile) {
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
				try {
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
				} catch (Throwable e) {
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
	 * Return a list of TaxIDs of all ancestors.
	 * @param tn0 Base node
	 * @param skipNonCanonical Ignore noncanonical (aka "nonsimple") levels like Tribe.
	 * @return List of TaxIDs.
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
	 * Return a list of all ancestors.
	 * @param tn0 Base node
	 * @param skipNonCanonical Ignore noncanonical (aka "nonsimple") levels like Tribe.
	 * @return List of ancestor nodes.
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
	 * @param s String to process.
	 * @return Decoded TaxID.
	 */
	public static int getID(byte[] s){return GiToTaxid.getID(s);}
	
	/** Return the lowest ancestor of the named node with taxonomic level at least minLevel */
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
		
		while(child!=parent && child.levelExtended<=parent.levelExtended && child.id!=child.pid){
			child=getNode(child.pid);
		}
		return child==parent;
	}
	
	/**
	 * Determine whether an organism is classified as X.
	 * @param taxID taxID of organism.
	 * @param ancestorID taxID of possible ancestor.
	 * @return true if the organism is an X.
	 */
	public boolean descendsFrom2(int taxID, final int ancestorID) {
		TaxNode tn=getNode(taxID);
		while(tn.id!=tn.pid){
			if(tn.id==ancestorID){return true;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	/** Determine whether an organism is classified as a plant. */
	public boolean isPlant(int taxID) {return descendsFrom2(taxID, VIRIDIPLANTAE_ID);}
	
	/** Determine whether an organism is classified as an animal. */
	public boolean isAnimal(int taxID) {return descendsFrom2(taxID, METAZOA_ID);}
	
	/** Determine whether an organism is classified as a fungus. */
	public boolean isFungus(int taxID) {return descendsFrom2(taxID, FUNGI_ID);}
	
	/** Determine whether an organism is classified as a eukaryote. */
	public boolean isEukaryote(int taxID) {return descendsFrom2(taxID, EUKARYOTA_ID);}
	
	/** Determine whether an organism is classified as a prokaryote. */
	public boolean isProkaryote(int taxID) {
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
	 * Calculate the taxLevel of the common ancestor of two nodes.
	 * @param a TaxID of a node.
	 * @param b TaxID of a node.
	 * @return Common ancestor ID of a and b.
	 */
	public int commonAncestorLevel(final int a, final int b){
		if(a<1 || b<1) {return -1;} 
		TaxNode an=getNode(a), bn=getNode(b);
		assert(an!=null) : "Invalid taxID: "+a;
		assert(bn!=null) : "Invalid taxID: "+b;
		TaxNode cn=commonAncestor(an, bn);
		assert(cn!=null) : "No common ancestor: "+an+", "+bn;
		if(cn==null){return -1;}
		return cn.level;
	}
	
	/**
	 * Calculate the common ancestor of two nodes.
	 * @param a TaxID of a node.
	 * @param b TaxID of a node.
	 * @return Common ancestor ID of a and b.
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
	 * Calculate the common ancestor of two nodes.
	 * @param a A node.
	 * @param b A node.
	 * @return Common ancestor of a and b.
	 */
	public TaxNode commonAncestor(TaxNode a, TaxNode b){
		assert(a!=null && b!=null) : "Null parameters.";
		if(a==null){return b;}
		if(b==null){return a;}
		
		while(a!=b){
			if(a.levelExtended<b.levelExtended){
				a=getNode(a.pid);
			}else{
				b=getNode(b.pid);
			}
		}
		return a;
	}
	
	/**
	 * Identify the highest ancestor of a node;
	 * this will presumably be "Life".
	 * @param a Node
	 * @return Highest ancestor
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
		try {
			id=Parse.parseInt(header, idx+1, idx2);
//			System.err.println("d"+", "+header.substring(idx+1, idx2));
		} catch (Throwable e) {
//			System.err.println("e"+", "+header.substring(idx+1, idx2));
			//ignore
		}
		return id;
	}
	
	/** New and improved.  Tree is optional. */
	public static int parseHeaderStatic2(String header, TaxTree tree) {
		int tid=-1;
		if(header==null) {return -1;}
		if(tid<0) {tid=parseTaxID(header, "tid");}
		if(tid<0) {tid=parseTaxID(header, "ncbi");}
		if(tid>=0 || tree==null) {return tid;}
		TaxNode tn=tree.parseNodeFromHeader(header, true);
		return tn==null ? tid : tn.id;
	}
	
	/** Header should contain something like "tid|1234" somewhere; 
	 * delimiter is autodetected. */
	private static int parseTaxID(String header, String term) {
		int loc=header.indexOf(term);
//		System.err.println("a) loc="+loc+", term="+term+", header="+header);
		if(loc<0) {return -1;}
		int start=loc+term.length();
//		System.err.println("b) start="+start);
		if(header.length()<start+2) {return -1;}
		final char c=header.charAt(start);
//		System.err.println("c) c="+c);
		if(Tools.isDigitOrSign(c)) {return -1;}
		start++;
		final char d=header.charAt(start);
//		System.err.println("d) d="+d);
		if(!Tools.isDigitOrSign(d)) {return -1;}
		int stop=start+1;
//		System.err.println("e) stop="+stop);
		while(stop<header.length() && Tools.isDigitOrSign(header.charAt(stop))) {stop++;}
//		System.err.println("f) stop="+stop);
//		System.err.println("g) start="+start+", stop="+stop+", int="+Parse.parseInt(header, start, stop));
		return Parse.parseInt(header, start, stop);
	}
	
	/**
	 * Determine the TaxID of a String.
	 * @param header Typically a sequence header or file name
	 * @param bestEffort In some cases, try certain substrings if the name is not found.
	 * @return
	 */
	public TaxNode parseNodeFromHeader(String header, boolean bestEffort){
		if(header==null || header.length()<2){return null;}
		if(header.charAt(0)=='>'){header=header.substring(1);}
		if(header.endsWith(".gz")) {header=header.substring(0, header.lastIndexOf('.'));}
		if(header.endsWith(".fa") || header.endsWith(".fna") || header.endsWith(".fasta") || header.endsWith(".faa")) {
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
		return tn;
	}
	
	/** 
	 * Guess the delimiter character in a String;
	 * typically assumed to be '|', '~', or ' '.
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
	
	/** Parses sequence headers using NCBI's old-style header system, prior to Accessions. */
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
					try {
						number=GiToTaxid.parseGiToTaxid(s, delimiter);
					} catch (Throwable e2) {
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
		if(Tools.isDigit(s.charAt(0)) && s.length()<=9){
			try {
				return getNode(Integer.parseInt(s));
			} catch (NumberFormatException e) {
				//ignore
			}
		}
		return null;
	}
	
	/** Parse a delimited number from a header, or return -1 if formatted incorrectly. */
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
	
	/** Parses sequence headers using NCBI's current header system, with Accessions. */
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
		
		if(number<0 && Tools.isDigit(s.charAt(0)) && s.length()<=9 && space<0){
			try {
				return getNode(Integer.parseInt(s));
			} catch (NumberFormatException e) {
				//ignore
			}
		}
		
		if(number<0 && space>0){
			number=parseNameToTaxid(s.substring(space+1));
		}
		
		if(number>-1){return getNode(number);}
		if(space<0 && s.indexOf('_')>0){
			return getNodeNewStyle(s.replace('_', ' '));
		}
		return null;
	}
	
	/**
	 * For parsing old-style NCBI headers.
	 */
	public int parseAccessionToTaxid(String[] split){
		if(split.length<4){
			return -1;
		}
		int ncbi=AccessionToTaxid.get(split[3]);
		return ncbi;
	}
	
	/**
	 * For parsing old-style NCBI headers.
	 */
	public int parseHeaderNameToTaxid(String[] split){
		if(split.length<5){
			return -1;
		}
		return parseNameToTaxid(split[4]);
	}
	
	/**
	 * Returns the TaxID from the organism's scientific name (e.g. "Homo sapiens").
	 * If multiple nodes share the same name, returns the first; to get the full list,
	 * use getNodesByNameExtended.
	 * @param name Organism name.
	 * @return Organism TaxID, or -1 if not found.
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
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		if(split2.length>6){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		if(split2.length>5){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>4){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4];
			list=getNodesByName(term);
//			System.err.println("5:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>3){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3];
			list=getNodesByName(term);
//			System.err.println("4:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>2){
			String term=split2[0]+" "+split2[1]+" "+split2[2];
			list=getNodesByName(term);
//			System.err.println("3:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>1){
			String term=split2[0]+" "+split2[1];
			list=getNodesByName(term);
//			System.err.println("2:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>0){
			String term=split2[0];
			list=getNodesByName(term);
//			System.err.println("1:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		return null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Assorted Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Return the TaxID of the lowest ancestor node at least the specified level,
	 * including this node itself.  Level is the normal (non-extended) level.
	 * @param taxID
	 * @param taxLevel
	 * @return
	 */
	public int promote(final int taxID, int taxLevel){
		TaxNode tn=null;
		tn=(taxID<1 ? null : getNode(taxID));
		tn=promote(tn, taxLevel);
		return (tn==null ? taxID : tn.id);
	}
	
	/**
	 * Fetch the first node in this node's lineage of at least the indicated level.
	 * This can be the node itself or an ancestor.
	 * @see getNodeAtLevelExtended
	 * @param tn Node in question
	 * @param taxLevel Desired minimum level
	 * @return A node at the desired level
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
	 * Determine the TaxID of the node's parent.
	 * @param id TaxID of child node
	 * @return Parent TaxID
	 */
	public int getParentID(int id){
		assert(id<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(id<0 || id>=nodes.length){return -1;}
		TaxNode tn=nodes[id];
		if(tn==null && mergedMap!=null){tn=getNode(mergedMap.get(id), true);}
		return tn==null ? -1 : tn.pid;
	}
	
	/**
	 * Fetch the node with this TaxID.
	 * @param id TaxID
	 * @return Node
	 */
	public TaxNode getNode(int id){
		assert(id<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(id<0 || id>=nodes.length){return null;}
		TaxNode tn=nodes[id];
		if(tn!=null || mergedMap==null){return tn;}
		return getNode(mergedMap.get(id), true);
	}
	
	public int toLevel(int id) {
		TaxNode tn=getNode(id);
		return tn==null ? -1 : tn.level;
	}
	
	public int resolveID(int id) {
		if(mergedMap==null) {return id;}
		int x=mergedMap.get(id);
		return x<0 ? id : x;
	}
	
	/**
	 * Fetch the node with this TaxID, but don't throw assertions upon failure.
	 * @param id TaxID
	 * @return Node
	 */	
	public TaxNode getNode(int id, boolean skipAssertion){
		assert(skipAssertion || id<nodes.length) : id+", "+nodes.length+"\nYou have encountered a TaxID more recent than your NCBI dump."
				+ "\nPlease redownload it and regenerate the taxtree.";
		if(id<0 || id>=nodes.length){return null;}
		TaxNode tn=nodes[id];
		if(tn!=null || mergedMap==null){return tn;}
		return getNode(mergedMap.get(id), true);
	}
	
	public TaxNode getNodeAtLevel(int id, int minLevel){
		return getNodeAtLevel(id, minLevel, DOMAIN);
	}
	
	public TaxNode getNodeAtLevelExtended(int id, int minLevelE){
		return getNodeAtLevelExtended(id, minLevelE, DOMAIN_E);
	}
	
	public TaxNode getNodeAtLevel(int id, int minLevel, int maxLevel){
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
		return getNodeAtLevelExtended(id, minLevelExtended, maxLevelExtended);
	}
	
	public TaxNode getNodeAtLevelExtended(int id, int minLevelE, int maxLevelE){
		TaxNode tn=getNode(id);
		while(tn!=null && tn.pid!=tn.id && tn.levelExtended<minLevelE){
			TaxNode temp=getNode(tn.pid);
			if(temp==null || temp.levelExtended>maxLevelE){break;}
			tn=temp;
		}
		return tn;
	}
	
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
	 * Fetch the node with this name.
	 * Throw an assertion if there are multiple such nodes.
	 * @param s Organism name.
	 * @return Node with given name.
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
	 * Fetch a list of all nodes with this name.
	 * @param s Organism name.
	 * @return Nodes with given name.
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
	
	private List<TaxNode> getNodesByName(String s, boolean lowercase){
		if(s==null){return null;}
		if(s.indexOf('_')>=0){s=s.replace('_', ' ');}
		if(lowercase){s=s.toLowerCase();}
//		System.err.println("Searching for "+s);
		final HashMap<String, ArrayList<TaxNode>> map=getMap(lowercase);
		ArrayList<TaxNode> list=map.get(s);
		if(list!=null){return list;}
//		System.err.println("No matches for '"+s+"'");
		
//		assert(false) : nameMap.containsKey(s)+", "+nameMapLower.containsKey(s);
		
		if(s.indexOf('_')<0 && s.indexOf(' ')<0){return null;}
		String[] split=delimiter2.split(lowercase ? s.toLowerCase() : s, 8);
//		System.err.println("Array: "+Arrays.toString(split));
		list=map.get(split[split.length-1]);
		if(list==null){return list;}
//		System.err.println(list==null ? "No matches for "+split[split.length-1] : "Found list( "+list.size()+")");
		
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
	
	public void incrementRaw(int id, long amt){
		assert(id>=0 && id<nodes.length) : "TaxID "+id+" is out of range."+(id<0 ? "" : "  Possibly the taxonomy data needs to be updated.");
		assert(nodes[id]!=null) : "No node for TaxID "+id+"; possibly the taxonomy data needs to be updated.";
		nodes[id].incrementRaw(amt);
	}
	
	public void percolateUp(){
		for(int i=0; i<treeLevelsExtended.length; i++){
			percolateUp(i);
		}
	}
	
	public void percolateUp(final int fromLevel){
		final TaxNode[] stratum=treeLevelsExtended[fromLevel];
		for(final TaxNode n : stratum){
			n.incrementSum(n.countRaw);
			TaxNode parent=nodes[n.pid];
			if(n!=parent){
				parent.incrementSum(n.countSum);
			}
		}
	}
	
	/** Add this amount to the node and all its ancestors. */
	public void percolateUp(TaxNode node, long amt){
		if(amt==0){return;}
		if(verbose){System.err.println("percolateUp("+amt+") node: "+node);}
		while(node.id!=node.pid){
			node.incrementSum(amt);
			node=nodes[node.pid];
		}
		node.incrementSum(amt);
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit){
		return gatherNodesAtLeastLimit(limit, 0, nodesPerLevelExtended.length-1);
	}
	
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
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimitExtended(final int fromLevelExtended, final long limit){
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		final TaxNode[] stratum=treeLevelsExtended[fromLevelExtended];
		for(final TaxNode n : stratum){
			if(n.countSum>=limit){
				list.add(n);
				TaxNode parent=nodes[n.pid];
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

	/**
	 * Generate the name to level number map.
	 */
	private static HashMap<String, Integer> makeLevelMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(31);
		for(int i=0; i<taxLevelNames.length; i++){
			map.put(taxLevelNames[i], i);
			map.put(taxLevelNames[i].toUpperCase(), i);
		}
		map.put("clade", NO_RANK);
		map.put("clade".toUpperCase(), NO_RANK);
		return map;
	}

	/**
	 * Generate the name to level number map.
	 */
	private static HashMap<String, Integer> makeShortLevelMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(31);
		for(int i=0; i<taxLevelNamesShort.length; i++){
			map.put(taxLevelNamesShort[i], i);
			map.put(taxLevelNamesShort[i].toUpperCase(), i);
		}
		return map;
	}

	/**
	 * Generate the name to extended level number map.
	 */
	private static HashMap<String, Integer> makeLevelMapExtended() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(129);
		for(int i=0; i<taxLevelNamesExtended.length; i++){
			map.put(taxLevelNamesExtended[i], i);
			map.put(taxLevelNamesExtended[i].toUpperCase(), i);
		}
		map.put("clade", NO_RANK_E);
		map.put("clade".toUpperCase(), NO_RANK_E);
		return map;
	}

	/**
	 * I think this maps normal and extend names to normal level numbers.
	 */
	private static HashMap<String, Integer> makeAltLevelMap() {
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
	
	/** Number of bp associated with this node in RefSeq */
	public long toSize(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSizeMap==null){return -1L;}
		final long x=refseqSizeMap.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Number of bp associated with this node and descendants in RefSeq */
	public long toSizeC(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSizeMapC==null){return -1L;}
		final long x=refseqSizeMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Number of sequences associated with this node in RefSeq */
	public int toSeqs(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSeqMap==null){return -1;}
		final int x=refseqSeqMap.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Number of sequences associated with this node and descandants in RefSeq */
	public long toSeqsC(TaxNode tn){
		if(tn==null){return 0;}
		if(refseqSeqMapC==null){return -1L;}
		final long x=refseqSeqMapC.get(tn.id);
		return Tools.max(0, x);
	}
	
	/** Number of descendants of this node */
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
	
	public static int imgToTaxid(long img){
		ImgRecord ir=imgMap.get(img);
//		assert(false) : "\n"+img+"\n"+imgMap.get(img)+"\n"+562+"\n"+imgMap.get(562)+"\n"+imgMap.size()+"\n"+IMGHQ+"\n"+defaultImgFile()+"\n";
		return ir==null ? -1 : ir.taxID;
	}
	
	public TaxNode imgToTaxNode(long img){
		int tid=imgToTaxid(img);
		return tid<1 ? null : getNode(tid);
	}
	
//	public static int loadIMGOld(String fname, boolean storeName, PrintStream outstream){
//		assert(imgMap==null);
//		if(fname==null){return 0;}
//		ImgRecord2.storeName=storeName;
//		if(outstream!=null){System.err.println("Loading IMG.");}
//		Timer t=new Timer(outstream, false);
//		ImgRecord2[] array=ImgRecord2.toArray(fname);
//		int x=loadIMG(array);
//		t.stopAndPrint();
//		return x;
//	}
	
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
	
	public static int loadIMG(ImgRecord[] array){
		assert(imgMap==null);
		imgMap=new HashMap<Long, ImgRecord>((int)(array.length*1.5));
		for(ImgRecord record : array){
			imgMap.put(record.imgID, record);
		}
		return imgMap.size();
	}
	
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
	
	public boolean isUnclassified(int tid){
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			if(tn.isUnclassified()){return true;}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	public boolean isEnvironmentalSample(int tid){
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			if(tn.isEnvironmentalSample()){return true;}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	public boolean isVirus(int tid){
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			if(tn.id==VIRUSES_ID){return true;}
			if(tn.pid==tn.id){break;}
			tn=getNode(tn.pid);
		}
		return false;
	}
	
	public long definedLevels(int tid){
		long levels=0;
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			levels=levels|(1L<<tn.level);
		}
		return levels;
	}
	
	public long definedLevelsExtended(int tid){
		long levels=0;
		TaxNode tn=getNode(tid);
		while(tn!=null && tn.id!=tn.pid){
			levels=levels|(1L<<tn.levelExtended);
		}
		return levels;
	}
	
	/** 
	 * Generates the mononomial name for this taxonomic level based on the scientific name.
	 * For example, "Homo sapiens" -> "Sapiens"
	 * @param tid TaxID
	 * @return Correct name for this node.
	 */
	public String mononomial(int tid){return mononomial(getNode(tid));}
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
	
	/** All nodes in the tree in a flat array, indexed by TaxiD */
	public final TaxNode[] nodes;
	
	/** Number of nodes per normal level */
	public final int[] nodesPerLevel=new int[taxLevelNames.length];
	
	/** Number of nodes per extended level */
	public final int[] nodesPerLevelExtended=new int[taxLevelNamesExtended.length];
	
	/** Number of nodes in the tree */
	public final int nodeCount;
	
	/** Maps old TaxIDs to new TaxIDs */
	public final IntHashMap mergedMap;

	/** Arrays of all nodes at a given taxonomic level (extended) */
	public final TaxNode[][] treeLevelsExtended=new TaxNode[taxLevelNamesExtended.length][];
	
	/** Map of names to nodes */
	HashMap<String, ArrayList<TaxNode>> nameMap;
	/** Map of lowercase names to nodes */
	HashMap<String, ArrayList<TaxNode>> nameMapLower;
	/** Map of nodes to child nodes */
	HashMap<TaxNode, ArrayList<TaxNode>> childMap;
	public HashMap<String, ArrayList<TaxNode>> nameMap(){return nameMap;}
	
	@Deprecated
	public int minValidTaxa=0; //TODO: Remove (will break serialization)
	
	/** Infer ranks for no-rank nodes, when possible */
	public boolean simplify=true;
	/** See simplify() for details, works in conjunction with simplify */
	public boolean reassign=true;
	/** Discard no-rank nodes */
	public boolean skipNorank=false;
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

	public static boolean verbose=false;
	public static boolean SHOW_WARNINGS=false;
	
	/** Maps IMG IDs to records from the dump file */
	private static HashMap<Long, ImgRecord> imgMap;

	/** Set to false if the tree is expected to be mutated.
	 * @TODO Remove mutable fields from the tree (like counters).
	 */
	public static boolean ALLOW_SHARED_TREE=true;
	
	/** Universal location of the shared TaxTree used by various classes */
	private static TaxTree sharedTree;
	
	/** A simpler and probably less safe version of sharedTree(...) */
	public static TaxTree getTree(){return sharedTree;}
	
	/**
	 * Fetch the shared tree, loading it from file if not present.
	 * @return A tree.
	 * @TODO: Check proper-construction of double-checked synchronize
	 */
	private static TaxTree sharedTree(String fname, boolean hashNames, boolean hashDotFormat, PrintStream outstream) {
		if(!ALLOW_SHARED_TREE){return null;}
		if(sharedTree==null && fname!=null){
			if("auto".equalsIgnoreCase(fname)){fname=defaultTreeFile();}
			synchronized(TaxTree.class){
				if(sharedTree==null){
//					if(outstream!=null){outstream.println("Loading tax tree.");}
					if(outstream!=null){outstream.println("Loading tax tree from "+fname);} //Useful for debugging
					Timer t=new Timer(outstream, false);
					setSharedTree(ReadWrite.read(TaxTree.class, fname, true), hashNames, hashDotFormat);
					t.stopAndPrint();
				}
			}
		}
		if(hashNames && sharedTree.nameMap==null){
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
	
	/** 
	 * For initialization.  Normally only one tree is needed by a process so it is set here.
	 * If the tree is already set nothing will happen, unless additional hashing is needed.
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
	 * Determine whether a taxonomic level is standard. e.g.:<br>
	 * isSimple("phylum")=true<br>
	 * isSimple("subphylum")=false<br>
	 * isSimple("no-rank")=false
	 * @param levelExtended The extended level to test.
	 * @return True if this level is not no-rank, and the names of the normal and extended levels match.
	 */
	public static boolean isSimple(int levelExtended){
		int level=extendedToLevel(levelExtended);
		return levelExtended!=NO_RANK_E && (levelExtended==levelToExtended(level));
	}
	
	/**
	 * Determine whether a taxonomic level is standard, but allows substrain and lower. e.g.:<br>
	 * isSimple("phylum")=true<br>
	 * isSimple("substrain")=true<br>
	 * isSimple("subphylum")=false<br>
	 * isSimple("no-rank")=false
	 * @param levelExtended The extended level to test.
	 * @return True if this level is not no-rank, and the names of the normal and extended levels match.
	 */
	public static boolean isSimple2(int levelExtended){
		int level=extendedToLevel(levelExtended);
		return levelExtended!=NO_RANK_E && (levelExtended==levelToExtended(level) 
				|| levelExtended==STRAIN_E || levelExtended==SUBSPECIES_E || levelExtended==SUBSTRAIN_E);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Get the number for the normal level of this name */
	public static final int stringToLevel(String s){return altLevelMap.get(s);}
	public static final int shortStringToLevel(String s){return shortLevelMap.get(s);}
	public static final boolean levelMapExtendedContains(String s){return levelMapExtended.containsKey(s);}
	/** Get the number for the extended level of this name */
	public static final int stringToLevelExtended(String s){return levelMapExtended.get(s);}
	/** Get the normal name for this normal level */
	public static final String levelToString(int x){return taxLevelNames[x];}
	/** Get the extended name for this extended level */
	public static final String levelToStringExtended(int x){return taxLevelNamesExtended[x];}
	/** Get the abbreviated name for this normal level */
	public static final String levelToStringShort(int x){return taxLevelNamesShort[x];}
	
	/** Normal, aka canonical, aka simple tax level names */
	private static final String[] taxLevelNames=new String[] {
		"no rank", "subspecies", "species", "genus",
		"family", "order", "class", "phylum",
		"kingdom", "superkingdom", "domain", "life"
	};
	public static final int numTaxLevelNames=taxLevelNames.length;
	
	/** 
	 * Definitive representation of all NCBI taxonomic level names.
	 * All levels used by NCBI must be present here, or parsing a new NCBI tax tree will crash.
	 * The first dimension maps normal ranks to extended ranks.
	 * Both dimensions are ordered ascending.
	 * @TODO Note! If this goes over 63 names it will cause a problem with getDefinedLevels().
	 */
	//TODO See @TODO
	private static final String[][] taxLevelNamesExtendedMatrix=new String[][] {
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
	
	/** Flatten the extended tax level names matrix to a 1D array */
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
	private static final String[] taxLevelNamesShort=new String[] {
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
	private static final int[] levelToExtended=new int[] {
			NO_RANK_E, SUBSPECIES_E, SPECIES_E, GENUS_E, FAMILY_E,
			ORDER_E, CLASS_E, PHYLUM_E, KINGDOM_E, SUPERKINGDOM_E, DOMAIN_E, LIFE_E
		};
	
	/** Map of extended to normal level numbers */
	private static final int[] extendedToLevel=makeExtendedToLevel();
	
	/** Creates extendedToLevel from taxaNamesExtendedMatrix during initialization. */
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
	
	/** Convert a standard level number (like KINGDOM) to extended (like KINGDOM_E). */
	public static final int levelToExtended(int level){
		return level<0 ? level : levelToExtended[level];
	}

	/** Convert an extended level number (like PHYLUM_E) to extended (like PHYLUM).
	 * Non-standard levels will be converted to the next higher standard level;
	 * e.g., subphylum -> phylum */
	public static final int extendedToLevel(int extended){
		return extended<0 ? -1 : extendedToLevel[extended];
	}
	
	/* Pre-compiled delimiters to save time when splitting lines */
	private static final Pattern delimiterTab = Pattern.compile("\t");
	private static final Pattern delimiter = Pattern.compile("\t\\|\t");
	private static final Pattern delimiterPipe = Pattern.compile("\\|");
	private static final Pattern delimiterTilde = Pattern.compile("\\~");
	private static final Pattern delimiter2 = Pattern.compile("[\\s_]+");
	
	public static boolean IMG_HQ=false;
	
	/* For these fields, see the corresponding functions, below.
	 * They define the default paths to various data on NERSC. */

	private static final String defaultTaxPathNersc="/global/cfs/cdirs/bbtools/tax/latest";
	private static final String defaultTaxPathAws="/test1/tax/latest";
	private static final String defaultTaxPathIGBVM="/data/tax/latest";
	private static final String defaultTaxPathDori="/clusterfs/jgi/groups/gentech/homes/bbushnell/tax/latest";
	private static final String default16SFileNersc="/global/cfs/cdirs/bbtools/silva/16S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	private static final String default16SFileAws="/test1/16S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	private static final String default16SFileIGBVM="/data/sketch/silva/16S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	private static final String default18SFileNersc="/global/cfs/cdirs/bbtools/silva/18S_consensus_silva_maxns10_taxsorted.fa.gz";
	private static final String default18SFileAws="/test1/18S_consensus_silva_maxns10_taxsorted.fa.gz";
	private static final String default18SFileIGBVM="/data/sketch/silva/18S_consensus_with_silva_maxns10_taxsorted.fa.gz";
	
	private static final String defaultImgFile="TAX_PATH/imgDump.txt";
	private static final String defaultTableFileInt="TAX_PATH/gitable.int1d.gz";
	private static final String defaultTableFile="TAX_PATH/gitable.int2d.gz";
	private static final String defaultTreeFile="TAX_PATH/tree.taxtree.gz";
	private static final String defaultPatternFile="TAX_PATH/patterns.txt";
	private static final String defaultSizeFile="TAX_PATH/taxsize.tsv.gz";

	private static final String defaultAccessionFile=
			//"TAX_PATH/shrunk.protF.accession2taxid.gz," +
			"TAX_PATH/shrunk.prot.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.nucl_wgs.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.nucl_gb.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.dead_prot.accession2taxid.gz,"
//			+ "TAX_PATH/shrunk.nucl_est.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.dead_wgs.accession2taxid.gz,"
//			+ "TAX_PATH/shrunk.nucl_gss.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.dead_nucl.accession2taxid.gz,"
			+ "TAX_PATH/shrunk.pdb.accession2taxid.gz";

	/** For setting TAX_PATH, the root to taxonomy files */
	public static final String defaultTaxPath(){
		return (Shared.AWS && !Shared.NERSC) ? 
				defaultTaxPathAws : 
					Shared.IGBVM ? defaultTaxPathIGBVM : 
						Shared.DORI ? defaultTaxPathDori : 
							(Shared.PERLMUTTER || Shared.NERSC) ? defaultTaxPathNersc :
								Data.ROOT()+"/resources/";
	}

	/** 16S consensus sequences per TaxID */
	public static final String default16SFile(){
		return (Shared.AWS && !Shared.NERSC) ? default16SFileAws : Shared.IGBVM ? default16SFileIGBVM : default16SFileNersc;
	}

	/** 18S consensus sequences per TaxID */
	public static final String default18SFile(){
		return (Shared.AWS && !Shared.NERSC) ? default18SFileAws : Shared.IGBVM ? default18SFileIGBVM : default18SFileNersc;
	}

	/** Path to all taxonomy files, substituted in to make specific file paths */
	public static String TAX_PATH=defaultTaxPath();
	
	/** Location of gitable.int2d.gz for gi lookups */
	public static final String defaultTableFile(){return defaultTableFile.replaceAll("TAX_PATH", TAX_PATH);}
	/** Location of tree.taxtree.gz */
	public static final String defaultTreeFile(){return defaultTreeFile.replaceAll("TAX_PATH", TAX_PATH);}
	
	//Use the prot.FULL.gz ncbi file.
	public static boolean protFull=false;
	
	/** Location of shrunk.*.accession2taxid.gz (all accession files, comma-delimited) */
	public static final String defaultAccessionFile(){
		String s=(protFull ? "TAX_PATH/shrunk.protF.accession2taxid.gz," : "")+defaultAccessionFile;
		return s.replaceAll("TAX_PATH", TAX_PATH);
	}
	/** Location of patterns.txt, which holds information about observed accession string formats */
	public static final String defaultPatternFile(){return defaultPatternFile.replaceAll("TAX_PATH", TAX_PATH);}
	/** Location of imgDump.txt, which translates IMG to NCBI IDs for internal JGI use */
	public static final String defaultImgFile(){return defaultImgFile.replaceAll("TAX_PATH", TAX_PATH);}
	/** Location of taxsize.tsv, which indicates the amount of sequence associated with a TaxID */
	public static final String defaultSizeFile(){return defaultSizeFile.replaceAll("TAX_PATH", TAX_PATH);}
	
	/** Screen output gets printed here */
	private static PrintStream outstream=System.out;
	
}
