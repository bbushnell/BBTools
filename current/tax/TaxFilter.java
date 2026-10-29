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

/** Reusable taxonomy include/exclude filter used by read, assembly and sketch
 * tools. A filter can match numeric NCBI taxIDs, names, regexes, substrings and
 * required ancestor ranks against a loaded {@link TaxTree}.
 *
 * @author Brian Bushnell
 * @date Nov 30, 2015
 */
public class TaxFilter implements Cloneable{

	/** Tests a regex pattern against a string and prints the match result.
	 * @param args regex pattern and test string */
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

	/** Builds a filter from command-line arguments and loads requested taxonomy
	 * resources before adding configured names and IDs.
	 * @param args command-line arguments
	 * @return configured taxonomy filter */
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

	/** Constructs a default tree-backed include or exclude filter.
	 * @param tree_ taxonomy tree used for lookups and ancestor traversal
	 * @param include_ true to include matching taxa; false to exclude them */
	public TaxFilter(TaxTree tree_, boolean include_){
		this(tree_, -1, 0, include_, true, null, null, null);
	}

	/** Constructs a fully configured include or exclude filter.
	 * @param tree_ taxonomy tree used for lookups and ancestor traversal
	 * @param taxLevelE_ extended rank to promote selected IDs to
	 * @param reqLevels_ bit mask of ranks that a passing node's ancestry must include
	 * @param include_ true to include matching taxa; false to exclude them
	 * @param promote_ true to search ancestors when looking for selected IDs
	 * @param taxSet_ optional borrowed set of selected NCBI taxIDs; this filter may mutate it
	 * @param regex_ optional regex for name-only filtering
	 * @param contains_ optional lowercase substring for name-only filtering */
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
		setContainsString(contains_);
	}

	/** Widens {@link #taxLevelE} upward until promoted selected IDs intersect
	 * promoted taxIDs observed in a FASTA input file.
	 * @param fname FASTA file to scan for taxonomic headers */
	public void reviseByBestEffort(String fname){
		HashSet<Integer> desired=new HashSet<Integer>();
		int currentLevelE=taxLevelE;
		for(int i : taxSet){
			int x=tree.getIdAtLevelExtended(i, currentLevelE);
			desired.add(x);
		}
		HashSet<Integer> present=new HashSet<Integer>();
		TextFile tf=new TextFile(fname);
		Throwable thrown=null;
		try{
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(line.startsWith(">")){
					TaxNode tn=tree.parseNodeFromHeader(line.substring(1), true);
					if(tn!=null){
						present.add(tree.getIdAtLevelExtended(tn.id, currentLevelE));
					}
				}
			}
		}catch(RuntimeException|Error e){
			thrown=e;
			throw e;
		}finally{
			closeOrSuppress(tf, fname, thrown);
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

	/** Checks if a parameter name is valid for TaxFilter configuration.
	 * @param a parameter name to validate, without its value
	 * @return true if the parameter is recognized, false otherwise */
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

	/** Loads an NCBI accession-to-taxid mapping table and binds its tree.
	 * @param accessionFile path to accession table file
	 * @param tree taxonomy tree to associate with accession data */
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
	 * @param fname path to gi table file */
	public static void loadGiTable(String fname){
		if(fname==null){return;}
		if(PRINT_STUFF){outstream.println("Loading gi table.");}
		GiToTaxid.initialize(fname);
	}

	/** Loads an NCBI taxonomy tree from file.
	 * @param fname path to taxonomy tree file
	 * @return loaded TaxTree instance, or null if no file was specified */
	public static TaxTree loadTree(String fname){
		if(fname==null){return null;}
		TaxTree tt=TaxTree.loadTaxTree(fname, PRINT_STUFF ? outstream : null, true, false);
		assert(tt.nameMap!=null);
		return tt;
	}

	/** Adds comma-separated taxonomy names or IDs to the filter set.
	 * Automatically detects whether each entry is a name or numeric ID.
	 * @param names comma-separated list of taxonomy names or IDs
	 * @param promote whether to promote matches up the taxonomic tree */
	public void addNamesOrNumbers(String names, boolean promote){
		if(names==null){return;}
		String[] array=names.split(",");
		for(String name : array){
			addNameOrNumber(name, promote);
		}
	}

	/** Adds a single taxonomy name or ID to the filter set.
	 * Determines if the string is numeric and routes to the appropriate handler.
	 * @param s taxonomy name or numeric ID
	 * @param promote whether to promote matches up the taxonomic tree */
	public void addNameOrNumber(String s, boolean promote){
		if(s==null || s.length()<1){return;}
		if(Tools.isDigit(s.charAt(0))){addNumber(Integer.parseInt(s), promote);}
		else{addName(s);}
	}

	/** Adds comma-separated taxonomy names to the filter set.
	 * @param names comma-separated list of taxonomy names */
	public void addNames(String names){
		if(names==null){return;}
		String[] array=names.split(",");
		for(String name : array){
			addName(name);
		}
	}

	/** Adds a taxonomy name to the filter set.
	 * First tries header parsing, then searches by extended name matching.
	 * @param name taxonomy name to add
	 * @return true if at least one node was added successfully */
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

	/** Adds taxonomy IDs from a comma-separated string or file to the filter set.
	 * Tokens can be numeric taxIDs, existing filenames containing more tokens, or
	 * organism names when a taxonomy tree is loaded.
	 * @param numbers comma-separated taxonomy IDs, file paths or organism names
	 * @param promote whether to promote matches up the taxonomic tree */
	public void addNumbers(String numbers, boolean promote){
		if(numbers==null){return;}
		String[] array=numbers.split(",");
		for(String s : array){
			// FIXED [tax/TaxFilter#001]: skip leading/interior empty tokens,
			// such as taxids=,5 or taxids=3,,4. The sibling addNameOrNumber
			// already guards empty names, and TextFile.nextLine drops blank file
			// lines; String.split already drops trailing empty tokens, so only
			// this comma-token path could hit s.charAt(0).
			if(s.isEmpty()){continue;}
			if(Tools.isDigit(s.charAt(0))){
				final int x=Integer.parseInt(s);
				addNumber(x, promote);
			}else if(new File(s).exists()){
				TextFile tf=new TextFile(s);
				Throwable thrown=null;
				try{
					for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
						addNumbers(line, promote);
					}
				}catch(RuntimeException|Error e){
					thrown=e;
					throw e;
				}finally{
					closeOrSuppress(tf, s, thrown);
				}
			}else{
				assert(tree!=null) : "Using organism names requires a taxonomic tree loaded; please use a numeric NCBI taxID.";
				final int x=tree.parseNameToTaxid(s);
				assert(x>0) : "Can't find a tax node for "+s;
				addNumber(x, promote);
			}
		}
	}

	/** Closes a TextFile, preserving any primary parse failure.
	 * @param tf TextFile to close
	 * @param fname Source path for diagnostics
	 * @param thrown Primary exception or error, or null if the parse succeeded */
	private static void closeOrSuppress(TextFile tf, String fname, Throwable thrown){
		try{
			if(tf.close()){throw new RuntimeException("Error reading "+fname);}
		}catch(RuntimeException|Error closeFailure){
			if(thrown==null){throw closeFailure;}
			thrown.addSuppressed(closeFailure);
		}
	}

	/** Adds a numeric taxonomy ID to the filter set.
	 * @param taxID NCBI taxonomy ID to add
	 * @param promote whether to add the node and eligible ancestors to the filter level
	 * @return true if the ID was added successfully */
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

	/** Adds a taxonomy node and eligible ancestors to the filter set.
	 * Preserves the original node ID and any accepted ancestor IDs up to the configured
	 * extended taxonomy level.
	 * @param tn taxonomy node to add
	 * @return true if the node was added successfully */
	public boolean addNode(TaxNode tn){
		if(tn==null){return false;}
		taxSet.add(tn.id);
		maxChildLevelExtended=Tools.max(maxChildLevelExtended, tn.maxChildLevelExtended, tn.levelExtended);
		if(printNodesAdded){System.err.println("Added node "+tn);}
		while(tn.id!=tn.pid && tn.levelExtended<taxLevelE){
			tn=tree.getNode(tn.pid);
			if(tn.levelExtended<=taxLevelE){
				if(printNodesAdded){System.err.println("Added node "+tn);}
				taxSet.add(tn.id);
			}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests if a read ID passes the filter.
	 * @param r read with taxonomic header information
	 * @return true if the read passes the filter criteria */
	public boolean passesFilter(final Read r){
		return passesFilter(r.id);
	}

	/** Tests if a sequence name passes regex and substring filters only.
	 * Does not perform taxonomic ID lookups.
	 * @param name sequence name to test
	 * @return true if the name passes all name-only criteria */
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

	/** Tests if a sequence name passes all filter criteria.
	 * Combines name-based filtering with taxonomic ID matching.
	 * @param name sequence name with taxonomic information
	 * @return true if the sequence passes all filter criteria */
	public boolean passesFilter(final String name){
		if(!passesFilterByNameOnly(name)){return false;}
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		TaxNode tn=tree.parseNodeFromHeader(name, true);
		if(tn==null){tn=tree.getNodeByName(name);}

		if(REQUIRE_PRESENT && tn==null){
			KillSwitch.kill("ERROR: Could not find node for '"+name+"'"
					+ "\nTo bypass this error, add the flag 'requirepresent=f'");
		}

		return passesFilter(tn);
	}

	/** Tests if a taxonomy ID passes the filter criteria.
	 * @param id NCBI taxonomy ID to test
	 * @return true if the ID passes the filter criteria */
	public boolean passesFilter(final int id){
		// NOTE [tax/TaxFilter]: empty filters intentionally return TRUE here
		// so sketch DisplayParams.postParse can build an include-mode,
		// name-only whitelist: the ID check passes first, then
		// passesFilterByNameOnly does the real contains-string filtering.
		// Do not unify this with passesFilter(String)/(TaxNode), whose
		// empty-filter result is !include by design.
		if((taxSet==null || taxSet.isEmpty()) && reqLevels==0){return true;}
		TaxNode tn=tree.getNode(id);

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

	/** Core filtering logic for taxonomy nodes.
	 * Checks taxonomic set membership and required ancestor levels, optionally
	 * promoting the node to its ancestors while looking for a selected ID.
	 *
	 * @param tn0 taxonomy node to test
	 * @return true if the node passes filter criteria */
	boolean passesFilter(final TaxNode tn0){
		TaxNode tn=tn0;
		if(taxSet.isEmpty() && reqLevels==0){return !include;}
		if(tn==null){
			assert(!REQUIRE_PRESENT) : "Null TaxNode.";
			return !include && reqLevels==0;
		}
		boolean found=taxSet.contains(tn.id);
		int levels=1<<tn.level;
		while((!found || (levels&reqLevels)!=reqLevels) && tn.id!=tn.pid){
			tn=tree.getNode(tn.pid);
			levels|=(1<<tn.level);
			if(promote){found=found||taxSet.contains(tn.id);}
		}
		return include==found && (levels&reqLevels)==reqLevels;
	}

	/** Fast filtering method that skips required level checking.
	 * Optimized for performance when ancestor level requirements are not needed.
	 * @param id NCBI taxonomy ID to test
	 * @return true if the ID passes the filter criteria */
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

	/** Fast filtering logic that uses maximum child level optimization.
	 * Stops early when nodes cannot contain target taxa.
	 * @param tn0 taxonomy node to test
	 * @return true if the node passes filter criteria */
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

	/** Tests if a string matches the configured regular expression pattern.
	 * @param s string to test against regex
	 * @return true if the string matches the regex pattern */
	boolean matchesRegex(String s){
		return regexPattern.matcher(s).matches();
	}

	/** Tests if a string contains the configured lowercase substring.
	 * @param s string to test for substring presence
	 * @return true if the string contains the target substring */
	boolean containsString(String s){
		return s.toLowerCase().contains(containsString);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Clone             ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a copy with an independent taxonomy set.
	 * The copy shares the taxonomy tree, immutable regex pattern and substring data.
	 * @return Deep copy of this TaxFilter instance */
	public TaxFilter deepCopy(){
		TaxFilter copy=null;
		try{
			copy=(TaxFilter) this.clone();
			if(taxSet!=null){
				copy.taxSet=new HashSet<Integer>();
				copy.taxSet.addAll(taxSet);
			}
		}catch(CloneNotSupportedException e){
			e.printStackTrace();
		}
		return copy;
	}

	/** Clears the taxonomy ID set by setting it to null. */
	public void clearSet(){
		taxSet=null;
	}

	/** Initializes the taxonomy ID set to an empty HashSet. */
	public void makeSet(){
		taxSet=new HashSet<Integer>();
	}


	/** Changes the taxonomic filtering level.
	 * The level can only be lowered when the taxonomy set is empty.
	 * @param newLevel new taxonomy level for filtering
	 * @param promote whether to promote existing nodes to the new level */
	public void setLevel(final int newLevel, boolean promote){
		final int newLevelE=TaxTree.levelToExtended(newLevel);
		assert(newLevelE>=taxLevelE || newLevelE<1 || taxSet==null || taxSet.isEmpty()) : "taxLevel may only be increased when the set is non-empty.";
		taxLevelE=newLevelE;
		if(promote){promote();}
	}

	/** Adds eligible ancestors for all nodes in the taxonomy set.
	 * Original IDs are preserved while ancestors up to the current filter level are added. */
	public void promote(){
		if(taxSet!=null && !taxSet.isEmpty() && taxLevelE>0){
			ArrayList<Integer> list=new ArrayList<Integer>(taxSet.size());
			list.addAll(taxSet);
			taxSet.clear();
			for(Integer i : list){addNumber(i, true);}
		}
	}

	/** Gets the number of taxonomy IDs in the filter set.
	 * @return number of selected TaxIDs, or zero when the backing set is null */
	public int size(){return (taxSet==null ? 0 : taxSet.size());}

	/** Gets the current taxonomic filtering level.
	 * @return regular taxonomy rank for the current extended filter level */
	public int taxLevel(){return TaxTree.extendedToLevel(taxLevelE);}
	/** Gets a snapshot array of taxonomy IDs in the filter set.
	 * @return new array containing selected TaxIDs, or null when the backing set is null or empty */
	public Integer[] taxSet(){
		return (taxSet==null || taxSet.isEmpty()) ? null : taxSet.toArray(new Integer[0]);
	}
	/** Gets whether the filter is in include mode.
	 * @return true when matching taxa pass, false when matching taxa are excluded */
	public boolean include(){return include;}

	/** Sets the taxonomy tree for this filter.
	 * @param tree_ new taxonomy tree to use */
	public void setTree(TaxTree tree_){tree=tree_;}

	/** Gets the taxonomy tree used by this filter.
	 * @return currently bound taxonomy tree */
	public TaxTree tree(){return tree;}

	/** Sets the substring matching pattern (converted to lowercase).
	 * @param s substring pattern to match, or null to disable */
	public void setContainsString(String s){containsString=(s==null ? null : s.toLowerCase());}

	/** Gets the current lowercase substring matching pattern.
	 * @return lowercase substring pattern, or null when substring filtering is disabled */
	public String containsString(){return containsString;}

	/** Converts the backing taxonomy set to text.
	 * @return string representation of the backing taxonomy ID set */
	@Override
	public String toString(){return ""+taxSet;}

	/** Tests whether the initialized taxonomy set has no entries.
	 * @return true when the backing taxonomy ID set is empty */
	public boolean isEmpty(){return taxSet.isEmpty();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Taxonomy tree for taxonomic lookups and navigation. */
	private TaxTree tree;

	/** Extended taxonomy level to which selected IDs are promoted. */
	private int taxLevelE;

	/** Bit mask of regular ranks required in a passing node's ancestry. */
	private final int reqLevels;

	/** Set of selected NCBI TaxIDs. */
	private HashSet<Integer> taxSet;

	/** Maximum extended child level in the taxonomy set for fast ancestor search. */
	private int maxChildLevelExtended=TaxTree.LIFE_E;
	/** True to include matching taxa; false to exclude them. */
	private final boolean include;
	/** Whether to promote matches up the taxonomic tree. */
	private boolean promote;

	/** Regular expression pattern for name matching. */
	private String regex;
	/** Compiled regular expression pattern for efficient matching. */
	private final Pattern regexPattern;
	/** Substring pattern for name matching. */
	private String containsString;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream. */
	private static PrintStream outstream=System.err;

	/** Print loading messages. */
	static boolean PRINT_STUFF=true;

	/** Whether to require all taxonomy nodes to be present in the tree. */
	public static boolean REQUIRE_PRESENT=true;
	/** Whether to warn once about absent taxonomy nodes. */
	public static boolean WARN_ABSENT_NODE=true;
	/** Whether to print debug information when nodes are added to the filter. */
	public static boolean printNodesAdded=true;

}
