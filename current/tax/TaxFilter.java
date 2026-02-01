package tax;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Pattern;

import fileIO.TextFile;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 30, 2015
 *
 */
public class TaxFilter {
	
	/** Test method for regex pattern matching.
	 * @param args Command line arguments: regex pattern and test string */
	public static void main(String[] args){
		String regex=args[0];
		String s=args[1];
		Pattern regexPattern=Pattern.compile(regex);
		boolean b=regexPattern.matcher(s).matches();
		System.err.println(regex);
		System.err.println(s);
		System.err.println(b);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public static TaxFilter makeFilter(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		String names=null;
		String ids=null;

		String giTableFile=null;
		String taxTreeFile=null;
		String accessionFile=null;

		int taxLevelE=-1;
		int reqLevels=0;
		boolean include=false;
		boolean promote=true;
		String regex=null;
		String contains=null;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("table") || a.equals("gi")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("level") || a.equals("lv") || a.equals("taxlevel") || a.equals("tl")){
				taxLevelE=TaxTree.parseLevelExtended(b);
//				System.err.println("Set taxLevelE = "+TaxTree.levelToStringExtended(taxLevelE)); //123
			}else if(a.equals("reqlevel") || a.equals("requiredlevel") || a.equals("reqlevels") || a.equals("requiredlevels")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.toLowerCase().split(",");
				reqLevels=0;
				for(String s : split2){
					int level=TaxTree.parseLevel(s);
					reqLevels|=(1<<level);
				}
			}else if(a.equals("name") || a.equals("names")){
				names=b;
			}else if(a.equals("regex")){
				regex=b;
			}else if(a.equals("contains")){
				contains=b;
			}else if(a.equals("printnodesadded") || a.equals("printnodes")){
				printNodesAdded=Parse.parseBoolean(b);
			}else if(a.equals("promote")){
				promote=Parse.parseBoolean(b);
			}else if(a.equals("include")){
				include=Parse.parseBoolean(b);
			}else if(a.equals("exclude")){
				include=!Parse.parseBoolean(b);
			}else if(a.equals("requirepresent")){
				REQUIRE_PRESENT=Parse.parseBoolean(b);
				TaxTree.SHOW_WARNINGS=REQUIRE_PRESENT;
			}else if(a.equals("id") || a.equals("ids") || a.equals("taxid") || a.equals("tid") || a.equals("taxids")){
				ids=b;
			}
		}
		
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		
		TaxTree tree=loadTree(taxTreeFile);
		loadGiTable(giTableFile);
		loadAccession(accessionFile, tree);
		
		TaxFilter filter=new TaxFilter(tree, taxLevelE, reqLevels, include, promote, null, regex, contains);
		filter.addNames(names);
		filter.addNumbers(ids, true);
		return filter;
	}
	
	/**
	 * Constructor.
	 * @param tree_
	 */
	public TaxFilter(TaxTree tree_, boolean include_){
		this(tree_, -1, 0, include_, true, null, null, null);
	}
	
	/**
	 * Constructor.
	 */
	public TaxFilter(TaxTree tree_, int taxLevelE_, int reqLevels_, boolean include_, boolean promote_,
			HashSet<Integer> taxSet_, String regex_, String contains_){
		tree=tree_;
		taxLevelE=taxLevelE_;
		reqLevels=reqLevels_;
		include=include_;
		promote=promote_;
		taxSet=(taxSet_==null ? new HashSet<Integer>() : taxSet_);
		regex=regex_;
		regexPattern=(regex==null ? null : Pattern.compile(regex));
		containsString=contains_;
	}
	
	/** Alter taxSet and taxLevel to ensure they intersect with the file contents. */
	public void reviseByBestEffort(String fname){
		HashSet<Integer> desired=new HashSet<Integer>();
		int currentLevelE=taxLevelE;
		for(int i : taxSet){
			int x=tree.getIdAtLevelExtended(i, currentLevelE);
			desired.add(x);
		}
		HashSet<Integer> present=new HashSet<Integer>();
		TextFile tf=new TextFile(fname);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith(">")){
				TaxNode tn=tree.parseNodeFromHeader(line.substring(1), true);
				if(tn!=null){
					present.add(tree.getIdAtLevelExtended(tn.id, currentLevelE));
				}
			}
		}
		
		while(currentLevelE<TaxTree.LIFE_E){
			//Intersect
			for(Integer i : desired){
				if(present.contains(i)){
					if(currentLevelE!=taxLevelE){System.err.println("Widened filter from "+
							TaxTree.levelToStringExtended(taxLevelE)+" to "+TaxTree.levelToStringExtended(currentLevelE));}
					taxLevelE=currentLevelE;
					taxSet=desired;
					return;
				}
			}
			currentLevelE++;
			
			HashSet<Integer> desired2=new HashSet<Integer>();
			for(int i : desired){
				desired2.add(tree.getIdAtLevelExtended(i, currentLevelE));
			}
			desired=desired2;
			
			HashSet<Integer> present2=new HashSet<Integer>();
			for(int i : present){
				present2.add(tree.getIdAtLevelExtended(i, currentLevelE));
			}
			present=present2;
		}
	}
	
	/**
	 * Checks if a parameter name is valid for TaxFilter configuration.
	 * @param a Parameter name to validate
	 * @return true if the parameter is recognized, false otherwise
	 */
	public static boolean validArgument(String a){
		if(a.equals("table") || a.equals("gi")){
		}else if(a.equals("tree") || a.equals("taxtree")){
		}else if(a.equals("accession")){
		}else if(a.equals("taxpath")){
		}else if(a.equals("level") || a.equals("lv") || a.equals("taxlevel") || a.equals("tl")){
		}else if(a.equals("name") || a.equals("names")){
		}else if(a.equals("regex")){
		}else if(a.equals("contains")){
		}else if(a.equals("besteffort")){
		}else if(a.equals("include")){
		}else if(a.equals("promote")){
		}else if(a.equals("exclude")){
		}else if(a.equals("printnodesadded") || a.equals("printnodes")){
		}else if(a.equals("id") || a.equals("ids") || a.equals("taxid") || a.equals("tid") || a.equals("taxids")){
		}else if(a.equals("requirepresent")){
		}else if(a.equals("reqlevel") || a.equals("requiredlevel") || a.equals("reqlevels") || a.equals("requiredlevels")){
		}else{
			return false;
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads NCBI accession-to-taxid mapping table.
	 * @param accessionFile Path to accession table file
	 * @param tree Taxonomy tree to associate with accession data
	 */
	private static void loadAccession(String accessionFile, TaxTree tree){
		if(accessionFile!=null){
			AccessionToTaxid.tree=tree;
			assert(tree!=null);
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			System.gc();
		}
	}
	
	/** Loads NCBI gi-to-taxid mapping table.
	 * @param fname Path to gi table file */
	public static void loadGiTable(String fname){
		if(fname==null){return;}
		if(PRINT_STUFF){outstream.println("Loading gi table.");}
		GiToTaxid.initialize(fname);
	}
	
	/**
	 * Loads NCBI taxonomy tree from file.
	 * @param fname Path to taxonomy tree file
	 * @return Loaded TaxTree instance or null if no file specified
	 */
	public static TaxTree loadTree(String fname){
		if(fname==null){return null;}
		TaxTree tt=TaxTree.loadTaxTree(fname, PRINT_STUFF ? outstream : null, true, false);
		assert(tt.nameMap!=null);
		return tt;
	}
	
	/**
	 * Adds comma-separated taxonomy names or IDs to the filter set.
	 * Automatically detects whether each entry is a name or numeric ID.
	 * @param names Comma-separated list of taxonomy names or IDs
	 * @param promote Whether to promote matches up the taxonomic tree
	 */
	public void addNamesOrNumbers(String names, boolean promote){
		if(names==null){return;}
		String[] array=names.split(",");
		for(String name : array){
			addNameOrNumber(name, promote);
		}
	}
	
	/**
	 * Adds a single taxonomy name or ID to the filter set.
	 * Determines if the string is numeric and routes to appropriate handler.
	 * @param s Taxonomy name or numeric ID
	 * @param promote Whether to promote matches up the taxonomic tree
	 */
	public void addNameOrNumber(String s, boolean promote){
		if(s==null || s.length()<1){return;}
		if(Tools.isDigit(s.charAt(0))){addNumber(Integer.parseInt(s), promote);}
		else{addName(s);}
	}
	
	/** Adds comma-separated taxonomy names to the filter set.
	 * @param names Comma-separated list of taxonomy names */
	public void addNames(String names){
		if(names==null){return;}
		String[] array=names.split(",");
		for(String name : array){
			addName(name);
		}
	}
	
	/**
	 * Adds a taxonomy name to the filter set.
	 * First tries header parsing, then searches by extended name matching.
	 * @param name Taxonomy name to add
	 * @return true if at least one node was added successfully
	 */
	public boolean addName(String name){
		{
			TaxNode tn=tree.parseNodeFromHeader(name, true);
			if(tn!=null){return addNode(tn);}
		}
		List<TaxNode> list=tree.getNodesByNameExtended(name);
		boolean success=false;
		assert(list!=null) : "Could not find a node for '"+name+"'";
		if(list==null){throw new RuntimeException("Could not find a node for '"+name+"'");}
		for(TaxNode tn : list){
			success=addNode(tn)|success;
		}
		return success;
	}
	
	/**
	 * Adds taxonomy IDs from string or file to the filter set.
	 * Supports comma-separated IDs, file paths, or organism names.
	 * @param numbers Comma-separated taxonomy IDs, file path, or organism names
	 * @param promote Whether to promote matches up the taxonomic tree
	 */
	public void addNumbers(String numbers, boolean promote){
		if(numbers==null){return;}
		String[] array=numbers.split(",");
		for(String s : array){
			if(Tools.isDigit(s.charAt(0))){
				final int x=Integer.parseInt(s);
				addNumber(x, promote);
			}else if(new File(s).exists()){
				TextFile tf=new TextFile(s);
				for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
					addNumbers(line, promote);
				}
			}else{
				assert(tree!=null) : "Using organism names requires a taxonomic tree loaded; please use a numeric NCBI taxID.";
				final int x=tree.parseNameToTaxid(s);
				assert(x>0) : "Can't find a tax node for "+s;
				addNumber(x, promote);
			}
		}
	}
	
	/**
	 * Adds a numeric taxonomy ID to the filter set.
	 * @param taxID NCBI taxonomy ID to add
	 * @param promote Whether to promote the node and ancestors to the filter level
	 * @return true if the ID was added successfully
	 */
	public boolean addNumber(int taxID, boolean promote){
		if(promote){
			TaxNode tn=tree.getNode(taxID);
			assert(tn!=null) : "Could not find a node for '"+taxID+"'";
			return addNode(tn);
		}else{
			maxChildLevelExtended=TaxTree.LIFE_E;
			return taxSet.add(taxID);
		}
	}
	
	/**
	 * Adds a taxonomy node and its ancestors to the filter set.
	 * Promotes the node up to the specified taxonomic level if needed.
	 * @param tn Taxonomy node to add
	 * @return true if the node was added successfully
	 */
	public boolean addNode(TaxNode tn){
		if(tn==null){return false;}
		taxSet.add(tn.id);
		maxChildLevelExtended=Tools.max(maxChildLevelExtended, tn.maxChildLevelExtended, tn.levelExtended);
		if(printNodesAdded){System.err.println("Added node "+tn);}//123
		while(tn.id!=tn.pid && tn.levelExtended<taxLevelE){
			tn=tree.getNode(tn.pid);
			if(tn.levelExtended<=taxLevelE){
				if(printNodesAdded){System.err.println("Added node "+tn);}//123
				taxSet.add(tn.id);
			}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Tests if a sequence read passes the taxonomic filter.
	 * @param r Sequence read with taxonomic header information
	 * @return true if the read passes the filter criteria
	 */
	public boolean passesFilter(final Read r){
		return passesFilter(r.id);
	}
	
	/**
	 * Tests if a sequence name passes regex and substring filters only.
	 * Does not perform taxonomic ID lookups.
	 * @param name Sequence name to test
	 * @return true if the name matches regex/substring criteria
	 */
	public boolean passesFilterByNameOnly(final String name){
		if(regexPattern!=null){
			boolean b=matchesRegex(name);
			if(b!=include){return false;}
		}
		if(containsString!=null){
			boolean b=containsString(name);
			if(b!=include){return false;}
		}
		return true;
	}
	
	/**
	 * Tests if a sequence name passes all filter criteria.
	 * Combines name-based filtering with taxonomic ID matching.
	 * @param name Sequence name with taxonomic information
	 * @return true if the sequence passes all filter criteria
	 */
	public boolean passesFilter(final String name){
		if(!passesFilterByNameOnly(name)){return false;}
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		TaxNode tn=tree.parseNodeFromHeader(name, true);
		if(tn==null){tn=tree.getNodeByName(name);}
//		assert(tn!=null || !REQUIRE_PRESENT) : "Could not find node for '"+name+"'";
		
		if(REQUIRE_PRESENT && tn==null){
			KillSwitch.kill("ERROR: Could not find node for '"+name+"'"
					+ "\nTo bypass this error, add the flag 'requirepresent=f'");
		}
		
//		assert(false) : passesFilter(tn);
		return passesFilter(tn);
	}
	
	/**
	 * Tests if a taxonomy ID passes the filter criteria.
	 * @param id NCBI taxonomy ID to test
	 * @return true if the ID passes the filter criteria
	 */
	public boolean passesFilter(final int id){
//		if((taxSet==null || taxSet.isEmpty()) && reqLevels==0){return !include;}
		if((taxSet==null || taxSet.isEmpty()) && reqLevels==0){return true;}
		TaxNode tn=tree.getNode(id);
//		assert(tn!=null || !REQUIRE_PRESENT) : "Could not find node number "+id;
		
		if(tn==null){
			if(REQUIRE_PRESENT){ 
				KillSwitch.kill("ERROR: Could not find node for "+id);
			}else if(WARN_ABSENT_NODE){
				System.err.println("Warning: Could not find node for "+id);
				WARN_ABSENT_NODE=false;
			}
		}
		
		return passesFilter(tn);
	}
	
	/**
	 * Core filtering logic for taxonomy nodes.
	 * Checks taxonomic set membership and required ancestor levels.
	 * Supports promotion up the taxonomic tree.
	 *
	 * @param tn0 Taxonomy node to test
	 * @return true if the node passes filter criteria
	 */
	boolean passesFilter(final TaxNode tn0){
		TaxNode tn=tn0;
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		if(tn==null){
			assert(!REQUIRE_PRESENT) : "Null TaxNode.";
			return !include && reqLevels==0;
		}
		boolean found=taxSet.contains(tn.id);
//		System.err.println("found="+found+", node="+tn);
		int levels=1<<tn.level;
		while((!found || (levels&reqLevels)!=reqLevels) && tn.id!=tn.pid){
			tn=tree.getNode(tn.pid);
			levels|=(1<<tn.level);
			if(promote){found=found||taxSet.contains(tn.id);}
//			System.err.println("found="+found+", node="+tn);
		}
//		assert(false) : levels+", "+reqLevels+", "+tn0+", "+tree.getAncestors(tn0.pid);
		return include==found && (levels&reqLevels)==reqLevels;
	}
	
	/**
	 * Fast filtering method that skips required level checking.
	 * Optimized for performance when ancestor level requirements are not needed.
	 * @param id NCBI taxonomy ID to test
	 * @return true if the ID passes the filter criteria
	 */
	public boolean passesFilterFast(final int id){
		assert(reqLevels==0);
		if(taxSet==null || taxSet.isEmpty()){return true;}
		TaxNode tn=tree.getNode(id);
		
		if(tn==null){
			if(REQUIRE_PRESENT){ 
				KillSwitch.kill("ERROR: Could not find node for "+id);
			}else if(WARN_ABSENT_NODE){
				System.err.println("Warning: Could not find node for "+id);
				WARN_ABSENT_NODE=false;
			}
		}
		
		return passesFilterFast(tn);
	}
	
	/**
	 * Fast filtering logic that uses maximum child level optimization.
	 * Stops early when nodes cannot contain target taxa.
	 * @param tn0 Taxonomy node to test
	 * @return true if the node passes filter criteria
	 */
	boolean passesFilterFast(final TaxNode tn0){
		TaxNode tn=tn0;
		if(taxSet==null || taxSet.isEmpty()){return true;}
		if(tn==null){
			assert(!REQUIRE_PRESENT) : "Null TaxNode.";
			return !include;
		}
		boolean found=taxSet.contains(tn.id);
		if(!promote){return include==found;}
		
		while(!found && tn.maxChildLevelExtended<=maxChildLevelExtended && tn.id!=tn.pid){
			tn=tree.getNode(tn.pid);
			found=found||taxSet.contains(tn.id);
		}
		return include==found;
	}
	
	/**
	 * Tests if a string matches the configured regular expression pattern.
	 * @param s String to test against regex
	 * @return true if the string matches the regex pattern
	 */
	boolean matchesRegex(String s){
		return regexPattern.matcher(s).matches();
	}
	
	/**
	 * Tests if a string contains the configured substring (case-insensitive).
	 * @param s String to test for substring presence
	 * @return true if the string contains the target substring
	 */
	boolean containsString(String s){
		return s.toLowerCase().contains(containsString);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Clone             ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a deep copy of this TaxFilter with independent taxonomy set.
	 * @return Deep copy of this TaxFilter instance */
	public TaxFilter deepCopy() {
		TaxFilter copy=null;
		try {
			copy=(TaxFilter) this.clone();
			if(taxSet!=null){
				copy.taxSet=new HashSet<Integer>();
				copy.taxSet.addAll(taxSet);
			}
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return copy;
	}
	
	/** Clears the taxonomy ID set by setting it to null */
	public void clearSet(){
		taxSet=null;
	}
	
	/** Initializes the taxonomy ID set to an empty HashSet */
	public void makeSet(){
		taxSet=new HashSet<Integer>();
	}
	
//	public void setInclude(boolean b){
//		include=b;
//	}
	
	/**
	 * Changes the taxonomic filtering level.
	 * Level can only be increased when the taxonomy set is non-empty.
	 * @param newLevel New taxonomy level for filtering
	 * @param promote Whether to promote existing nodes to the new level
	 */
	public void setLevel(final int newLevel, boolean promote){
		final int newLevelE=TaxTree.levelToExtended(newLevel);
		assert(newLevelE>=taxLevelE || newLevelE<1 || taxSet==null && taxSet.isEmpty()) : "taxLevel may only be increased when the set is non-empty.";
		taxLevelE=newLevelE;
		if(promote){promote();}
	}
	
	/** Promotes all nodes in the taxonomy set to the current filtering level.
	 * Replaces specific nodes with their ancestors at the target level. */
	public void promote(){
		if(taxSet!=null && !taxSet.isEmpty() && taxLevelE>0){
			ArrayList<Integer> list=new ArrayList<Integer>(taxSet.size());
			list.addAll(taxSet);
			taxSet.clear();
			for(Integer i : list){addNumber(i, true);}
		}
	}
	
	/** Gets the number of taxonomy IDs in the filter set */
	public int size(){return (taxSet==null ? 0 : taxSet.size());}

	/** Gets the current taxonomic filtering level */
	public int taxLevel(){return TaxTree.extendedToLevel(taxLevelE);}
	/** Gets array of taxonomy IDs in the filter set, or null if empty */
	public Integer[] taxSet(){
		return (taxSet==null || taxSet.isEmpty()) ? null : taxSet.toArray(new Integer[0]);
	}
	/** Gets whether the filter is in include mode (true) or exclude mode (false) */
	public boolean include(){return include;}
	/** Sets the taxonomy tree for this filter.
	 * @param tree_ New taxonomy tree to use */
	public void setTree(TaxTree tree_){tree=tree_;}
	/** Gets the taxonomy tree used by this filter */
	public TaxTree tree(){return tree;}
	/** Sets the substring matching pattern (converted to lowercase).
	 * @param s Substring pattern to match, or null to disable */
	public void setContainsString(String s){containsString=(s==null ? null : s.toLowerCase());}
	/** Gets the current substring matching pattern */
	public String containsString(){return containsString;}
	@Override
	public String toString(){return ""+taxSet;}
	
	/** Returns true if the taxonomy set is empty */
	public boolean isEmpty() {return taxSet.isEmpty();}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Taxonomy tree for taxonomic lookups and navigation */
	private TaxTree tree;
	
	/** Level at which to filter */
	private int taxLevelE;
	
	/** Branch must contain ancestors at these levels (bitflag) */
	private final int reqLevels;
	
	/** Set of numeric NCBI TaxIDs */
	private HashSet<Integer> taxSet;
	
	/**
	 * Maximum extended level of child nodes in the taxonomy set for optimization
	 */
	private int maxChildLevelExtended=TaxTree.LIFE_E;
	/** If true, include matching sequences; if false, exclude them */
	private final boolean include;
	/** Whether to promote matches up the taxonomic tree */
	private boolean promote;
	
	/** Regular expression pattern for name matching */
	private String regex;
	/** Compiled regular expression pattern for efficient matching */
	private final Pattern regexPattern;
	/** Substring pattern for name matching (stored in lowercase) */
	private String containsString;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private static PrintStream outstream=System.err;
	
	/** Print loading messages */
	static boolean PRINT_STUFF=true;
	
	/** Whether to require all taxonomy nodes to be present in the tree */
	public static boolean REQUIRE_PRESENT=true;
	/** Whether to warn about absent taxonomy nodes (one-time warning) */
	public static boolean WARN_ABSENT_NODE=true;
	/** Whether to print debug information when nodes are added to the filter */
	public static boolean printNodesAdded=true;

}
