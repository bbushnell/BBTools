package stream;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;

import tax.TaxNode;
import tax.TaxTree;

/**
 * Collects Read references into reusable buckets for dynamic demultiplexing.
 * String-created buckets also occupy integer indexes; integer lookup addresses that
 * shared index array, not a decimal string key. Integer-only buckets are not added
 * to the name list consumed by MultiCros. Reads are stored as references, not copied.
 * Draining detaches a bucket's current list while retaining its name/index metadata.
 * This class is not thread-safe; use one instance per thread or coordinate all access.
 * The stored ordered option does not change collection or output ordering here.
 *
 * @author Brian Bushnell
 * @date Apr 2, 2015
 */
public class ArrayListSet{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a set without a taxonomy tree and with phylum as the requested extended rank.
	 * Name-based taxonomic lookup requires a tree and is unavailable with this constructor;
	 * direct string, integer-index and numeric-taxid grouping do not need one.
	 * @param ordered_ Retained compatibility option; currently unused by this class
	 */
	public ArrayListSet(boolean ordered_){
		this(ordered_, null, TaxTree.stringToLevelExtended("phylum"));
	}

	/**
	 * Creates a bucket set with optional name-to-taxid lookup support.
	 * @param ordered_ Retained compatibility option; does not enable ordering logic
	 * @param tree_ Tree used by name-based taxonomy methods, or null when those are unused
	 * @param taxLevelE_ Requested minimum extended rank passed to TaxTree.getNode
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
	 * Adds the same Read reference once for every supplied name, without deduplication.
	 * Repeated names append the reference repeatedly to the same nonnull-name bucket.
	 * @param r Read reference to store without copying
	 * @param names Nonnull iterable of bucket names
	 */
	public void add(Read r, Iterable<String> names){
		for(String s : names){add(r, s);}
	}
	
	/**
	 * Appends a reference to a named bucket, creating the bucket on first use.
	 * A null name creates a new bucket without a string-map entry on every call; null is not
	 * stored in the string map and cannot retrieve that bucket through string lookup.
	 * @param r Read reference to store without copying
	 * @param name Bucket name; a null entry is still appended to the name list
	 */
	public void add(Read r, String name){
		final Pack p=getPack(name, true);
		p.add(r);
	}
	
	/**
	 * Appends a reference to an indexed bucket, extending the index array if needed.
	 * An existing string-created bucket at this index is shared. A newly created
	 * integer-only bucket does not contribute to getNames() or size().
	 * @param r Read reference to store without copying
	 * @param id Nonnegative bucket index, not a taxonomic ID or decimal string name
	 */
	public void add(Read r, int id){
		final Pack p=getPack(id, true);
		p.add(r);
	}
	
	/**
	 * Detaches the current list for a string key without removing the bucket metadata.
	 * Later additions allocate a new list; the returned list and Read references are not copied.
	 * @param name String lookup key; null has no string-map entry
	 * @return Detached list, or null if the key is absent or its bucket has no current list
	 */
	public ArrayList<Read> getAndClear(String name){
		final Pack p=getPack(name, false);
		return p==null ? null : p.getAndClear();
	}
	
	/**
	 * Detaches the current list at an index without removing the bucket or changing names.
	 * @param id Nonnegative index in the shared string/integer bucket array
	 * @return Detached list, or null if the slot is absent or has no current list
	 */
	public ArrayList<Read> getAndClear(int id){
		final Pack p=getPack(id, false);
		return p==null ? null : p.getAndClear();
	}
	
	/**
	 * Exposes the live name list in string-bucket creation order, including null names.
	 * Integer-only buckets are omitted. Callers must treat the returned collection as
	 * read-only because mutations do not update the bucket array or string map.
	 * @return The backing name list, not a snapshot of nonempty buckets
	 */
	public Collection<String> getNames(){
		return nameList;
	}
	
	/** @return Number of name-list entries, including drained buckets; not total buckets or reads */
	public int size(){return nameList.size();}
	
	/*--------------------------------------------------------------*/
	/*----------------        TaxId Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Resolves a name through the configured TaxTree and adds under the resulting decimal taxid key.
	 * Lookup requests taxLevelE; a missing node uses the shared "-1" bucket.
	 * @param r Read reference to store without copying
	 * @param name Name/header accepted by the configured tree's lookup; a nonnull tree is required
	 */
	public void addByTaxid(Read r, String name){
		addByTaxid(r, nameToTaxid(name));
	}
	
	/**
	 * Adds under a decimal string key without consulting the tree or changing taxonomic rank.
	 * This uses the same key as add(r, Integer.toString(taxid)). The integer overload of
	 * add interprets the number as a bucket index, which may happen to address that bucket.
	 * @param r Read reference to store without copying
	 * @param taxid ID used verbatim in the key, including negative values
	 */
	public void addByTaxid(Read r, int taxid){
		String key=Integer.toString(taxid);
		final Pack p=getPack(key, true);
		p.add(r);
	}
	
	/**
	 * Adds to distinct resolved taxid buckets, with empty and single-name fast paths.
	 * @param r Read reference to store without copying
	 * @param names Nonnull name list; a tree is required when a name is resolved
	 */
	public void addByTaxid(Read r, ArrayList<String> names){
		if(names.size()==0){return;}
		else if(names.size()==1){addByTaxid(r, names.get(0));}
		else{addByTaxid(r, (Iterable<String>)names);}
	}
	
	/**
	 * Resolves names and adds the reference once per distinct resulting taxid for this call.
	 * Multiple names mapping to the same node or to -1 are deduplicated. Scratch IDs are
	 * thread-local, but bucket mutation is not synchronized. Set traversal does not promise
	 * the input-name order for newly created buckets; deduplication does not span calls.
	 * @param r Read reference to store without copying
	 * @param names Nonnull iterable; a tree is required when a name is resolved
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
	 * Delegates name resolution and ancestor selection to the configured tree.
	 * @param name Name/header to resolve; tree must be nonnull
	 * @return Resolved node ID, or -1 if lookup returns no node
	 */
	private int nameToTaxid(String name){
		TaxNode tn=tree.getNode(name, taxLevelE);
		return (tn==null ? -1 :tn.id);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Finds a string-mapped bucket, optionally creating one; null names are never mapped.
	 * @param name String lookup key
	 * @param add Whether to create a bucket when lookup finds none
	 * @return Existing/new bucket, or null when absent and creation is disabled
	 */
	private Pack getPack(String name, boolean add){
		Pack p=stringMap.get(name);
		if(p==null && add){p=new Pack(name);}
		return p;
	}
	
	/**
	 * Finds an indexed bucket, optionally extending the index array and creating one.
	 * @param id Nonnegative index in the shared bucket array
	 * @param add Whether to create a bucket when the slot is absent
	 * @return Existing/new bucket, or null when absent and creation is disabled
	 */
	private Pack getPack(int id, boolean add){
		Pack p=packList.size()>id ? packList.get(id) : null;
		if(p==null && add){p=new Pack(id);}
		return p;
	}
	
	/** @return String representation of the live name list, without bucket contents */
	@Override
	public String toString(){
		return nameList.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** One persistent bucket identity whose detachable read list is allocated on demand. */
	private class Pack{
		
		/** Appends a string-created bucket to both lists and maps its name when nonnull. */
		Pack(String s){
			assert(s==null || !stringMap.containsKey(s));
			name=s;
			id=packList.size();
			nameList.add(s);
			packList.add(this);
			if(s!=null){stringMap.put(s, this);}
		}
		
		/** Occupies a nonnegative index, extending with null slots and leaving names unchanged. */
		Pack(int x){
			name=null;
			id=x;
			while(packList.size()<=x){packList.add(null);}
			assert(packList.get(x)==null);
			packList.set(x, this);
		}
		
		/** Appends the supplied reference, allocating a list for the first addition after a drain. */
		public void add(Read r){
			if(list==null){list=new ArrayList<Read>();}
			list.add(r);
		}
		
		/** @return Current list, detached without copying; null when no list is allocated */
		public ArrayList<Read> getAndClear(){
			ArrayList<Read> temp=list;
			list=null;
			return temp;
		}
		
		/** @return Diagnostic bucket name; integer-only buckets report null as their name */
		@Override
		public String toString(){
			return "Pack "+name;
		}
		
		/** String name, or null for unnamed/integer-only buckets. */
		final String name;
		/** Assigned index in the outer bucket array; retained for metadata only. */
		@SuppressWarnings("unused")
		final int id;
		/** Current accumulated references, or null before allocation/after detachment. */
		private ArrayList<Read> list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Compatibility flag stored but unused; output ordering is handled outside this class. */
	private final boolean ordered;//NOTE: set from the ctor but never read anywhere - vestigial (ordered demux was never wired in here; MultiCros handles ordering).
	/** Live string-created name entries in creation order; draining does not remove them. */
	private final ArrayList<String> nameList=new ArrayList<String>();
	/** Shared index space for named and integer-only buckets, possibly containing null slots. */
	private final ArrayList<Pack> packList=new ArrayList<Pack>();
	/** Lookup map for nonnull names, including decimal keys made by addByTaxid. */
	private final LinkedHashMap<String, Pack> stringMap=new LinkedHashMap<String, Pack>();
	/** Requested extended rank for name-based TaxTree lookup. */
	private final int taxLevelE;
	/** Optional name-resolution tree; numeric taxid grouping does not consult it. */
	private final TaxTree tree;
	/** Per-thread scratch deduplication set; does not protect the outer collections. */
	private final ThreadLocal<HashSet<Integer>> tls=new ThreadLocal<HashSet<Integer>>();
	
	
}
