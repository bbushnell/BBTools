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
	
	public void add(Read r, Iterable<String> names){
		for(String s : names){add(r, s);}
	}
	
	public void add(Read r, String name){
		final Pack p=getPack(name, true);
		p.add(r);
	}
	
	public void add(Read r, int id){
		final Pack p=getPack(id, true);
		p.add(r);
	}
	
	public ArrayList<Read> getAndClear(String name){
		final Pack p=getPack(name, false);
		return p==null ? null : p.getAndClear();
	}
	
	public ArrayList<Read> getAndClear(int id){
		final Pack p=getPack(id, false);
		return p==null ? null : p.getAndClear();
	}
	
	public Collection<String> getNames(){
		return nameList;
	}
	
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
	
	public void addByTaxid(Read r, int taxid){
		String key=Integer.toString(taxid);
		final Pack p=getPack(key, true);
		p.add(r);
	}
	
	public void addByTaxid(Read r, ArrayList<String> names){
		if(names.size()==0){return;}
		else if(names.size()==1){addByTaxid(r, names.get(0));}
		else{addByTaxid(r, (Iterable<String>)names);}
	}
	
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
	
	private int nameToTaxid(String name){
		TaxNode tn=tree.getNode(name, taxLevelE);
		return (tn==null ? -1 :tn.id);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private Pack getPack(String name, boolean add){
		Pack p=stringMap.get(name);
		if(p==null && add){p=new Pack(name);}
		return p;
	}
	
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
	
	private class Pack {
		
		Pack(String s){
			assert(s==null || !stringMap.containsKey(s));
			name=s;
			id=packList.size();
			nameList.add(s);
			packList.add(this);
			if(s!=null){stringMap.put(s, this);}
		}
		
		Pack(int x){
			name=null;
			id=x;
			while(packList.size()<=x){packList.add(null);}
			assert(packList.get(x)==null);
			packList.set(x, this);
		}
		
		public void add(Read r){
			if(list==null){list=new ArrayList<Read>();}
			list.add(r);
		}
		
		public ArrayList<Read> getAndClear(){
			ArrayList<Read> temp=list;
			list=null;
			return temp;
		}
		
		@Override
		public String toString(){
			return "Pack "+name;
		}
		
		final String name;
		@SuppressWarnings("unused")
		final int id;
		private ArrayList<Read> list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final boolean ordered;
	private final ArrayList<String> nameList=new ArrayList<String>();
	private final ArrayList<Pack> packList=new ArrayList<Pack>();
	private final LinkedHashMap<String, Pack> stringMap=new LinkedHashMap<String, Pack>();
	private final int taxLevelE;
	private final TaxTree tree;
	private final ThreadLocal<HashSet<Integer>> tls=new ThreadLocal<HashSet<Integer>>();
	
	
}
