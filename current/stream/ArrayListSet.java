package stream;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;

import tax.TaxNode;
import tax.TaxTree;

/**
 * Associates reads with named lists.
 * Designed for dynamically demultiplexing reads into output streams with MultiCros.
 * This class is not thread-safe; one should be instantiated per thread.
 * @author Brian Bushnell
 * @date Apr 2, 2015
 *
 */
public class ArrayListSet {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an ArrayListSet with default taxonomy settings.
	 * Uses phylum level as the default taxonomic level.
	 * @param ordered_ Whether input order should be maintained (unimplemented)
	 */
	public ArrayListSet(boolean ordered_){
		this(ordered_, null, TaxTree.stringToLevelExtended("phylum"));
	}

	/**
	 * Create an ArrayListSet with an optional TaxTree and level.
	 * The tree is to assign reads to a list based on the taxonomy of the name,
	 * rather than the name itself.
	 * @param ordered_ Whether input order should be maintained.  Unimplemented.
	 * @param tree_ A taxonomic tree.
	 * @param taxLevelE_ The minimum level in the tree to stop.
	 */
	public ArrayListSet(boolean ordered_, TaxTree tree_, int taxLevelE_){
		ordered=ordered_;
		tree=tree_;
		taxLevelE=taxLevelE_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a read to lists corresponding to multiple names.
	 * @param r The read to add
	 * @param names Collection of names to associate with the read
	 */
	public void add(Read r, Iterable<String> names){
		for(String s : names){add(r, s);}
	}
	
	/**
	 * Adds a read to the list corresponding to the specified name.
	 * @param r The read to add
	 * @param name The name/key to associate with the read
	 */
	public void add(Read r, String name){
		final Pack p=getPack(name, true);
		p.add(r);
	}
	
	/**
	 * Adds a read to the list corresponding to the specified numeric ID.
	 * @param r The read to add
	 * @param id The numeric ID to associate with the read
	 */
	public void add(Read r, int id){
		final Pack p=getPack(id, true);
		p.add(r);
	}
	
	/**
	 * Retrieves and clears the read list for the specified name.
	 * @param name The name/key to retrieve reads for
	 * @return The list of reads associated with the name, or null if none exist
	 */
	public ArrayList<Read> getAndClear(String name){
		final Pack p=getPack(name, false);
		return p==null ? null : p.getAndClear();
	}
	
	/**
	 * Retrieves and clears the read list for the specified numeric ID.
	 * @param id The numeric ID to retrieve reads for
	 * @return The list of reads associated with the ID, or null if none exist
	 */
	public ArrayList<Read> getAndClear(int id){
		final Pack p=getPack(id, false);
		return p==null ? null : p.getAndClear();
	}
	
	/** Returns the collection of all names currently tracked.
	 * @return Collection of all registered names */
	public Collection<String> getNames(){
		return nameList;
	}
	
	/** Returns the number of named lists currently tracked */
	public int size(){return nameList.size();}
	
	/*--------------------------------------------------------------*/
	/*----------------        TaxId Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Look up the sequence name, which should start with a gi or ncbi number, and
	 * associate the read with the ancestor node at some taxonomic level.
	 * @param r
	 * @param name
	 */
	public void addByTaxid(Read r, String name){
		addByTaxid(r, nameToTaxid(name));
	}
	
	/**
	 * Associates a read with a taxonomic ID.
	 * @param r The read to add
	 * @param taxid The taxonomic ID to associate with the read
	 */
	public void addByTaxid(Read r, int taxid){
		String key=Integer.toString(taxid);
		final Pack p=getPack(key, true);
		p.add(r);
	}
	
	/**
	 * Associates a read with multiple taxonomic names.
	 * Handles empty, single, or multiple name cases efficiently.
	 * @param r The read to add
	 * @param names List of names containing taxonomic identifiers
	 */
	public void addByTaxid(Read r, ArrayList<String> names){
		if(names.size()==0){return;}
		else if(names.size()==1){addByTaxid(r, names.get(0));}
		else{addByTaxid(r, (Iterable<String>)names);}
	}
	
	/**
	 * Associates a read with multiple taxonomic names from an iterable.
	 * Uses thread-local storage for efficient taxonomic ID deduplication.
	 * @param r The read to add
	 * @param names Iterable of names containing taxonomic identifiers
	 */
	public void addByTaxid(Read r, Iterable<String> names){
		HashSet<Integer> idset=tls.get();
		if(idset==null){
			idset=new HashSet<Integer>();
			tls.set(idset);
		}
		assert(idset.isEmpty());
		for(String s : names){
			idset.add(nameToTaxid(s));
		}
		for(Integer i : idset){
			addByTaxid(r, i);
		}
		idset.clear();
	}
	
	/**
	 * Converts a sequence name to a taxonomic ID using the taxonomy tree.
	 * @param name Sequence name to look up
	 * @return Taxonomic ID, or -1 if not found
	 */
	private int nameToTaxid(String name){
		TaxNode tn=tree.getNode(name, taxLevelE);
		return (tn==null ? -1 :tn.id);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Retrieves or creates a Pack for the specified name.
	 * @param name The name to look up
	 * @param add Whether to create a new Pack if none exists
	 * @return The Pack associated with the name, or null if not found and add is false
	 */
	private Pack getPack(String name, boolean add){
		Pack p=stringMap.get(name);
		if(p==null && add){p=new Pack(name);}
		return p;
	}
	
	/**
	 * Retrieves or creates a Pack for the specified numeric ID.
	 * @param id The numeric ID to look up
	 * @param add Whether to create a new Pack if none exists
	 * @return The Pack associated with the ID, or null if not found and add is false
	 */
	private Pack getPack(int id, boolean add){
		Pack p=packList.size()>id ? packList.get(id) : null;
		if(p==null && add){p=new Pack(id);}
		return p;
	}
	
	@Override
	public String toString(){
		return nameList.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Internal container class that holds reads associated with a specific name or ID.
	 * Manages the mapping between names/IDs and their corresponding read lists. */
	private class Pack {
		
		/**
		 * Creates a Pack with a string name.
		 * Registers the Pack in both name and ID mappings.
		 * @param s The name to associate with this Pack
		 */
		Pack(String s){
			assert(s==null || !stringMap.containsKey(s));
			name=s;
			id=packList.size();
			nameList.add(s);
			packList.add(this);
			if(s!=null){stringMap.put(s, this);}
		}
		
		/**
		 * Creates a Pack with a numeric ID only.
		 * Expands the pack list as needed to accommodate the ID.
		 * @param x The numeric ID to associate with this Pack
		 */
		Pack(int x){
			name=null;
			id=x;
			while(packList.size()<=x){packList.add(null);}
			assert(packList.get(x)==null);
			packList.set(x, this);
		}
		
		/**
		 * Adds a read to this Pack's read list.
		 * Creates the list if it doesn't exist yet.
		 * @param r The read to add
		 */
		public void add(Read r){
			if(list==null){list=new ArrayList<Read>();}
			list.add(r);
		}
		
		/** Retrieves and clears the read list from this Pack.
		 * @return The list of reads, or null if no reads were added */
		public ArrayList<Read> getAndClear(){
			ArrayList<Read> temp=list;
			list=null;
			return temp;
		}
		
		@Override
		public String toString(){
			return "Pack "+name;
		}
		
		/** The name associated with this Pack */
		final String name;
		/** The numeric ID associated with this Pack */
		@SuppressWarnings("unused")
		final int id;
		/** The list of reads stored in this Pack */
		private ArrayList<Read> list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Whether input order should be maintained (unimplemented) */
	private final boolean ordered;
	/** List of all registered names in order of addition */
	private final ArrayList<String> nameList=new ArrayList<String>();
	/** List of all Packs indexed by numeric ID */
	private final ArrayList<Pack> packList=new ArrayList<Pack>();
	/** Map from names to their corresponding Packs */
	private final LinkedHashMap<String, Pack> stringMap=new LinkedHashMap<String, Pack>();
	/** The minimum taxonomic level to stop at when traversing the tree */
	private final int taxLevelE;
	/** Taxonomic tree used for name-to-taxid lookups */
	private final TaxTree tree;
	/**
	 * Thread-local storage for taxonomic ID deduplication during batch operations
	 */
	private final ThreadLocal<HashSet<Integer>> tls=new ThreadLocal<HashSet<Integer>>();
	
	
}
