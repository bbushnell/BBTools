package sketch;

import java.io.File;
import java.util.ArrayList;

import dna.Data;
import kmer.AbstractKmerTable;
import shared.PreParser;

public class Blacklist extends SketchObject {
	
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
	
	public static int increment(long key){
		int way=0;//(int)(key%ways);
		return keySets[way].increment(key, 1);
	}
	
	public static boolean contains(long key){
		if(keySets==null){return false;}
		int way=0;//(int)(key%ways);
		return keySets[way].getValue(key)>0;
	}
	
	public static boolean exists(){
		return keySets!=null;
	}
	
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
	
	static synchronized String ntBlacklist(){return ntBlacklist!=null ? ntBlacklist : (ntBlacklist=Data.findPath("?blacklist_nt_merged.sketch"));}
	static synchronized String silvaBlacklist(){return silvaBlacklist!=null ? silvaBlacklist : (silvaBlacklist=Data.findPath("?blacklist_silva_merged.sketch"));}
	static synchronized String refseqBlacklist(){return refseqBlacklist!=null ? refseqBlacklist : (refseqBlacklist=Data.findPath("?blacklist_refseq_merged.sketch"));}
	static synchronized String imgBlacklist(){return imgBlacklist!=null ? imgBlacklist : (imgBlacklist=Data.findPath("?blacklist_img_species_300.sketch"));}
	static synchronized String nrBlacklist(){return null;}//Data.findPath("?blacklist_nr_species_1000.sketch");
	static synchronized String prokProtBlacklist(){return prokProtBlacklist!=null ? prokProtBlacklist : (prokProtBlacklist=Data.findPath("?blacklist_prokprot_merged.sketch"));}
	static synchronized String mitoBlacklist(){return mitoBlacklist!=null ? mitoBlacklist : (mitoBlacklist=Data.findPath("?blacklist_refseq_merged.sketch"));}
	static synchronized String fungiBlacklist(){return fungiBlacklist!=null ? fungiBlacklist : (fungiBlacklist=Data.findPath("?blacklist_refseq_merged.sketch"));}
	
	private static String ntBlacklist;
	private static String silvaBlacklist;
	private static String refseqBlacklist;
	private static String imgBlacklist;
	private static String prokProtBlacklist;
	private static String nrBlacklist;
	private static String mitoBlacklist;
	private static String fungiBlacklist;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	static AbstractKmerTable[] keySets;
	private static final int ways=1;
//	private static final int initialSize=16000;
	private static ArrayList<String> added=new ArrayList<String>();
	
}
