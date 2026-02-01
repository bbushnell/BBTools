package sketch;

import java.io.File;
import java.util.ArrayList;

import dna.Data;
import kmer.AbstractKmerTable;
import parse.PreParser;

/**
 * Manages blacklist k-mer databases for filtering common or contaminating sequences.
 * Provides access to pre-built blacklists for various databases (NT, RefSeq, Silva, etc.)
 * and functionality to load custom blacklist files for sketch-based contamination filtering.
 *
 * @author Brian Bushnell
 * @date June 2, 2025
 */
public class Blacklist extends SketchObject {
	
	/**
	 * Parses blacklist arguments and sets the appropriate blacklist database.
	 * Recognizes predefined database names (nt, silva, refseq, img, nr, etc.)
	 * and sets the corresponding blacklist file path.
	 *
	 * @param arg The full argument string
	 * @param a The argument key
	 * @param b The argument value
	 * @return true if blacklist was successfully parsed and set
	 */
	public static boolean parseBlacklist(String arg, String a, String b) {
		String x;
		boolean force=false;
		if(a.equals("bl") || a.equals("blacklist")){
			force=true;
			x=b;
		}else{
			x=a;
		}
		
		if(force || blacklist==null){
			if(x.equalsIgnoreCase("nt")){
				blacklist=Blacklist.ntBlacklist();
			}else if(x.equalsIgnoreCase("silva") || x.equalsIgnoreCase("ribo")){
				blacklist=Blacklist.silvaBlacklist();
			}else if(x.equalsIgnoreCase("refseq") || x.equalsIgnoreCase("refseqbig")){
				blacklist=Blacklist.refseqBlacklist();
			}else if(x.equalsIgnoreCase("img")){
				blacklist=Blacklist.imgBlacklist();
			}else if(x.equalsIgnoreCase("nr")){
				blacklist=Blacklist.nrBlacklist();
			}else if(x.equalsIgnoreCase("refseqprot") || x.equalsIgnoreCase("prokprot") || x.equalsIgnoreCase("prokprotbig") 
					|| x.equalsIgnoreCase("protein") || x.equalsIgnoreCase("protien") || x.equalsIgnoreCase("prot")){
				blacklist=Blacklist.prokProtBlacklist();
			}else if(x.equalsIgnoreCase("mito")){
				blacklist=Blacklist.mitoBlacklist();
			}else if(x.equalsIgnoreCase("fungi")){
				blacklist=Blacklist.fungiBlacklist();
			}else if(force){
				blacklist=b;
			}else{
				return false;
			}
			return true;
		}
		return false;
	}
	
	/**
	 * Adds blacklist files to the k-mer table for filtering.
	 * Handles comma-separated lists of files and resolves predefined database names.
	 * @param fname Filename or comma-separated list of filenames to add
	 */
	public static void addFiles(String fname){
		if(fname==null){return;}
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		if(fname.indexOf(',')<0 || new File(fname).exists()){
			ArrayList<Sketch> temp=addFile(fname);
			sketches.addAll(temp);
		}else{
			String[] split=fname.split(",");
			for(String s : split){
				if(s!=null && !"null".equalsIgnoreCase(s)){
					ArrayList<Sketch> temp=addFile(s);
					sketches.addAll(temp);
				}
			}
		}
		addSketches(sketches);
	}
	
	/**
	 * Loads sketches from a single blacklist file.
	 * Resolves predefined database names to actual file paths and loads
	 * the sketches using SketchTool.
	 *
	 * @param fname Filename to load sketches from
	 * @return ArrayList of loaded sketches
	 */
	private static ArrayList<Sketch> addFile(String fname){
		if(fname==null){return null;}
		if(!new File(fname).exists()){
			if("nt".equalsIgnoreCase(fname)){
				fname=Blacklist.ntBlacklist();
			}else if(("silva".equalsIgnoreCase(fname) || "ribo".equalsIgnoreCase(fname))){
				fname=Blacklist.silvaBlacklist();
			}else if("refseq".equalsIgnoreCase(fname)){
				fname=Blacklist.refseqBlacklist();
			}else if("img".equalsIgnoreCase(fname)){
				fname=Blacklist.imgBlacklist();
			}else if("prokprot".equalsIgnoreCase(fname) || "protein".equalsIgnoreCase(fname) || "protien".equalsIgnoreCase(fname) || "prot".equalsIgnoreCase(fname)){
				fname=Blacklist.prokProtBlacklist();
			}else if("fungi".equalsIgnoreCase(fname)){
				fname=Blacklist.fungiBlacklist();
			}else if("mito".equalsIgnoreCase(fname)){
				fname=Blacklist.mitoBlacklist();
			}
		}
		if(!PreParser.silent){System.err.println("Adding "+fname+" to blacklist.");}
		assert(!added.contains(fname));
		added.add(fname);
		SketchTool tool=new SketchTool(1000000, 1, false, false, rcomp);
		ArrayList<Sketch> sketches=tool.loadSketchesFromFile(fname, null, 1, -1, ONE_SKETCH, 1f, 0f, 0f, (byte)0, false);
		return sketches;
	}
	
	/**
	 * Adds loaded sketches to the k-mer table for blacklist filtering.
	 * Initializes the k-mer table if necessary and stores all k-mers
	 * from the sketches with inverted keys (Long.MAX_VALUE - key).
	 * @param sketches List of sketches to add to the blacklist
	 */
	private static void addSketches(ArrayList<Sketch> sketches){
		if(sketches==null || sketches.isEmpty()){return;}
		long size=0;
		for(Sketch sk : sketches){
			size+=sk.length();
		}
		long size2=(size*4)/3;
		assert(size2>0 && size2+1000<Integer.MAX_VALUE) : size2;
		if(keySets==null){
			keySets=AbstractKmerTable.preallocate(ways, AbstractKmerTable.ARRAY1D, new int[] {(int)size2}, -1L);
			
//			keySets=AbstractKmerTable.preallocate(ways, AbstractKmerTable.ARRAY1D, (int)size2, -1L, true);
		}
		for(Sketch sk : sketches){
			for(long key : sk.keys){
				increment(Long.MAX_VALUE-key);
			}
		}
	}
	
	/**
	 * Increments the count for a blacklist k-mer.
	 * Uses way 0 for all keys in the current implementation.
	 * @param key The k-mer key to increment
	 * @return The new count for this k-mer
	 */
	public static int increment(long key){
		int way=0;//(int)(key%ways);
		return keySets[way].increment(key, 1);
	}
	
	/**
	 * Checks if a k-mer is present in the blacklist.
	 * Returns false if no blacklist is loaded.
	 * @param key The k-mer key to check
	 * @return true if the k-mer is blacklisted, false otherwise
	 */
	public static boolean contains(long key){
		if(keySets==null){return false;}
		int way=0;//(int)(key%ways);
		return keySets[way].getValue(key)>0;
	}
	
	/** Checks if a blacklist has been loaded.
	 * @return true if a blacklist k-mer table exists, false otherwise */
	public static boolean exists(){
		return keySets!=null;
	}
	
	/**
	 * Converts a blacklist name to its corresponding file path.
	 * Recognizes predefined database names and returns the appropriate
	 * file path using Data.findPath().
	 *
	 * @param b Blacklist name or file path
	 * @return Resolved file path for the blacklist
	 */
	public static final String toBlacklist(String b){
		String blacklist=b;
		if(b==null){
			//do nothing
		}else if(b.equalsIgnoreCase("nt")){
			blacklist=ntBlacklist();
		}else if(b.equalsIgnoreCase("refseq")){
			blacklist=refseqBlacklist();
		}else if(b.equalsIgnoreCase("silva") || b.equalsIgnoreCase("ribo")){
			blacklist=silvaBlacklist();
		}else if(b.equalsIgnoreCase("img")){
			blacklist=imgBlacklist();
		}else if(b.equalsIgnoreCase("refseqprot") || b.equalsIgnoreCase("prokprot") 
				|| b.equalsIgnoreCase("protein") || b.equalsIgnoreCase("protien") || b.equalsIgnoreCase("prot")){
			blacklist=prokProtBlacklist();
		}else if(b.equalsIgnoreCase("refseqmito") || b.equalsIgnoreCase("mito")){
			blacklist=mitoBlacklist();
		}else if(b.equalsIgnoreCase("refseqfungi") || b.equalsIgnoreCase("fungi")){
			blacklist=fungiBlacklist();
		}
		return blacklist;
	}
	
	/** Gets the NT (nucleotide) blacklist file path, loading it lazily */
	static synchronized String ntBlacklist(){return ntBlacklist!=null ? ntBlacklist : (ntBlacklist=Data.findPath("?blacklist_nt_merged.sketch"));}
	/** Gets the Silva ribosomal RNA blacklist file path, loading it lazily */
	static synchronized String silvaBlacklist(){return silvaBlacklist!=null ? silvaBlacklist : (silvaBlacklist=Data.findPath("?blacklist_silva_merged.sketch"));}
	/** Gets the RefSeq blacklist file path, loading it lazily */
	static synchronized String refseqBlacklist(){return refseqBlacklist!=null ? refseqBlacklist : (refseqBlacklist=Data.findPath("?blacklist_refseq_merged.sketch"));}
	/**
	 * Gets the IMG (Integrated Microbial Genomes) blacklist file path, loading it lazily
	 */
	static synchronized String imgBlacklist(){return imgBlacklist!=null ? imgBlacklist : (imgBlacklist=Data.findPath("?blacklist_img_species_300.sketch"));}
	/** Gets the NR (non-redundant) blacklist file path - currently returns null */
	static synchronized String nrBlacklist(){return null;}//Data.findPath("?blacklist_nr_species_1000.sketch");
	/** Gets the prokaryotic protein blacklist file path, loading it lazily */
	static synchronized String prokProtBlacklist(){return prokProtBlacklist!=null ? prokProtBlacklist : (prokProtBlacklist=Data.findPath("?blacklist_prokprot_merged.sketch"));}
	/** Gets the mitochondrial blacklist file path, loading it lazily */
	static synchronized String mitoBlacklist(){return mitoBlacklist!=null ? mitoBlacklist : (mitoBlacklist=Data.findPath("?blacklist_refseq_merged.sketch"));}
	/** Gets the fungi blacklist file path, loading it lazily */
	static synchronized String fungiBlacklist(){return fungiBlacklist!=null ? fungiBlacklist : (fungiBlacklist=Data.findPath("?blacklist_refseq_merged.sketch"));}
	
	/** Cached path to NT blacklist file */
	private static String ntBlacklist;
	/** Cached path to Silva ribosomal RNA blacklist file */
	private static String silvaBlacklist;
	/** Cached path to RefSeq blacklist file */
	private static String refseqBlacklist;
	/** Cached path to IMG blacklist file */
	private static String imgBlacklist;
	/** Cached path to prokaryotic protein blacklist file */
	private static String prokProtBlacklist;
	/** Cached path to NR blacklist file */
	private static String nrBlacklist;
	/** Cached path to mitochondrial blacklist file */
	private static String mitoBlacklist;
	/** Cached path to fungi blacklist file */
	private static String fungiBlacklist;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	static AbstractKmerTable[] keySets;
	/** Number of hash table ways for distributing k-mers (currently 1) */
	private static final int ways=1;
//	private static final int initialSize=16000;
	/** List of blacklist files that have been added to prevent duplicates */
	private static ArrayList<String> added=new ArrayList<String>();
	
}
