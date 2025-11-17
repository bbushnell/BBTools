package prok;

import java.io.File;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ListNum;
import structures.LongHashSet;

/** Contains a lot of statics and static methods for gene-calling */
public abstract class ProkObject {
	
	/**
	 * Parses command-line arguments for ribosomal RNA and strand processing parameters.
	 * Handles slop settings, identity thresholds, sequence loading flags, and k-mer parameters.
	 *
	 * @param arg The full argument string (unused)
	 * @param a The parameter name
	 * @param b The parameter value
	 * @return true if the parameter was recognized and parsed, false otherwise
	 */
	public static boolean parse(String arg, String a, String b){
		if(a.equalsIgnoreCase("16sstartslop") || a.equalsIgnoreCase("ssustartslop")){
			ssuStartSlop=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("23sstartslop") || a.equalsIgnoreCase("lsustartslop")){
			lsuStartSlop=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("5sstartslop")){
			r5SStartSlop=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("16sstopslop") || a.equalsIgnoreCase("ssustopslop")){
			ssuStopSlop=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("23sstopslop") || a.equalsIgnoreCase("lsustopslop")){
			lsuStopSlop=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("5sstopslop")){
			r5SStopSlop=Integer.parseInt(b);
		}else if(a.equals("plus")){
			PROCESS_PLUS_STRAND=Parse.parseBoolean(b);
		}else if(a.equals("minus")){
			PROCESS_MINUS_STRAND=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("min16SIdentity") || a.equalsIgnoreCase("min16SId")) {
			min16SIdentity=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("min18SIdentity") || a.equalsIgnoreCase("min18SId")) {
			min18SIdentity=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("min23SIdentity") || a.equalsIgnoreCase("min23SId")) {
			min23SIdentity=Float.parseFloat(b);
		}else if(a.equalsIgnoreCase("min5SIdentity") || a.equalsIgnoreCase("min5SId")) {
			min5SIdentity=Float.parseFloat(b);
		}			
		
		else if(a.equalsIgnoreCase("align16s") || a.equalsIgnoreCase("load16SSequence")){
			load16SSequence=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("align23s") || a.equalsIgnoreCase("load23SSequence")){
			load23SSequence=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("align18s") || a.equalsIgnoreCase("load18SSequence")){
			load18SSequence=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("align5s") || a.equalsIgnoreCase("load5SSequence")){
			load5SSequence=Parse.parseBoolean(b);
		}
		
		else if(a.equalsIgnoreCase("load16skmers") || a.equalsIgnoreCase("load18skmers") || a.equalsIgnoreCase("loadssukmers")){
			loadSSUkmers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("load23skmers") || a.equalsIgnoreCase("load28skmers") || a.equalsIgnoreCase("loadlsukmers")){
			loadLSUkmers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("load5skmers")){
			load5Skmers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("loadtrnakmers")){
			loadtRNAkmers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("klongtrna")){
			kLongTRna=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("longkmers")){
			loadSSUkmers=loadLSUkmers=load5Skmers=loadtRNAkmers=Parse.parseBoolean(b);
		}else if(a.equalsIgnoreCase("klong5s")){
			kLong5S=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("klong16s") || a.equalsIgnoreCase("klong18s") || a.equalsIgnoreCase("klongssu")){
			kLongSSU=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("klong23s") || a.equalsIgnoreCase("klong28s") || a.equalsIgnoreCase("klonglsu")){
			kLongLSU=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("klongtrna")){
			kLongTRna=Integer.parseInt(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Determines whether a given gene type should be processed based on current settings.
	 * @param type The gene type constant (CDS, r16S, r23S, r18S, r5S, tRNA)
	 * @return true if this type should be processed, false otherwise
	 */
	public static boolean processType(int type){
		return (type==CDS ? callCDS : type==r16S ? call16S : type==r23S ? call23S : type==r18S ? call18S : type==r5S ? call5S : type==tRNA ? calltRNA : true);
	}
	
	/**
	 * Gets the start position tolerance for a given ribosomal RNA type.
	 * @param type The ribosomal RNA type constant
	 * @return The start slop value in bases, or 9999 for non-rRNA types
	 */
	public static int startSlop(int type) {
		int slop=(type==r16S ? ssuStartSlop : type==r23S ? lsuStartSlop : type==r18S ? ssuStartSlop : type==r5S ? r5SStartSlop : 9999);
		return slop;
	}
	
	/**
	 * Gets the stop position tolerance for a given ribosomal RNA type.
	 * @param type The ribosomal RNA type constant
	 * @return The stop slop value in bases, or 9999 for non-rRNA types
	 */
	public static int stopSlop(int type) {
		int slop=(type==r16S ? ssuStopSlop : type==r23S ? lsuStopSlop : type==r18S ? ssuStopSlop : type==r5S ? r5SStopSlop : 9999);
		return slop;
	}
	
	/**
	 * Gets the minimum identity threshold for a given ribosomal RNA type.
	 * @param type The ribosomal RNA type constant
	 * @return The minimum identity threshold as a fraction, or 0 for non-rRNA types
	 */
	public static float minID(int type) {
		float minIdentity=(type==r16S ? min16SIdentity : type==r23S ? min23SIdentity : type==r18S ? min18SIdentity : type==r5S ? min5SIdentity : 0);
		return minIdentity;
	}
	
	/**
	 * Gets the consensus sequence reads for a given ribosomal RNA type.
	 * @param type The ribosomal RNA type constant
	 * @return Array of consensus reads, or null if not available for this type
	 */
	public static Read[] consensusReads(int type) {
		Read[] consensusReads=(type==r16S ? r16SSequence : type==r23S ? r23SSequence : type==r18S ? r18SSequence : type==r5S ? r5SSequence : null);
		return consensusReads;
	}
	
	/**
	 * Gets the k-mer set for a given gene type.
	 * @param type The gene type constant
	 * @return The LongHashSet containing k-mers for this type, or null if not available
	 */
	public static LongHashSet kmerSet(int type) {
		LongHashSet set=(type==tRNA ? trnaKmers : type==r16S ? ssuKmers : type==r23S ? lsuKmers : type==r5S ? r5SKmers : type==r18S ? ssuKmers : null);
		return set;
	}
	
	/**
	 * Gets the long k-mer length for a given gene type.
	 * @param type The gene type constant
	 * @return The k-mer length for this type, or -1 if not applicable
	 */
	public static int kLongLen(int type) {
		int kLongLen=(type==tRNA ? kLongTRna : type==r16S ? kLongSSU : type==r23S ? kLongLSU : type==r5S ? kLong5S : type==r18S ? kLongSSU : -1);
		return kLongLen;
	}
	
	/**
	 * Converts a bit flag to the corresponding gene type constant.
	 * Uses trailing zero count to decode the flag position.
	 * @param flag The bit flag
	 * @return The gene type constant
	 */
	public static int flagToType(int flag) {
		return Integer.numberOfTrailingZeros(flag)+1;
	}
	
	/**
	 * Converts a gene type constant to the corresponding bit flag.
	 * @param type The gene type constant (must be <= 6)
	 * @return The bit flag as a byte
	 */
	public static byte typeToFlag(int type) {
		assert(type<=6);
		return (byte)(1<<(type-1));
	}
	
	/**
	 * Determines whether a specific gene type should be called based on current settings.
	 * @param type The gene type constant
	 * @return true if this type should be called, false otherwise
	 */
	public static boolean callType(int type){//TODO: Turn these functions into array lookups
		if(type==CDS){return callCDS;}
		else if(type==tRNA){return calltRNA;}
		else if(type==r16S){return call16S;}
		else if(type==r23S){return call23S;}
		else if(type==r5S){return call5S;}
		else if(type==r18S){return call18S;}
		assert(false) : type;
		return false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Long Kmers          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads long k-mer sets for all enabled gene types from reference files.
	 * Thread-safe method that loads SSU, LSU, 5S, and tRNA k-mers as configured.
	 * Only loads once per execution to avoid redundant file I/O.
	 */
	public static synchronized void loadLongKmers(){
//		assert(ssuKmers==null);
//		assert(false) : load5Skmers+", "+kLong5s;
		if(loadedLongKmers){return;}
		if(loadSSUkmers){ssuKmers=loadLongKmersByType(kLongSSU, "ssu");}
		if(loadLSUkmers){lsuKmers=loadLongKmersByType(kLongLSU, "lsu");}
		if(load5Skmers){r5SKmers=loadLongKmersByType(kLong5S, "5S");}
		if(loadtRNAkmers){trnaKmers=loadLongKmersByType(kLongTRna, "tRNA");}
		loadedLongKmers=true;
	}
	
//	private static LongHashSet loadLongKmers(StatsContainer sc, int k, String prefix){
//		String fname=Data.findPath("?"+prefix+"_"+k+"mers.fa");
//		if(!new File(fname).exists()){
//			fname=fname+".gz";
//			if(!new File(fname).exists()){
//				System.err.println("Can't find "+fname);
//				return null;
//			}
//		}
//		LongHashSet set=loadLongKmers(fname, k);
//		sc.kmerSet=set;
//		sc.kLongLen=k;
//		return set;
//	}
	
	/**
	 * Loads long k-mers from a specific reference file type.
	 * Searches for files with pattern "prefix_kmers.fa" in the data directory.
	 *
	 * @param k The k-mer length
	 * @param prefix The file prefix (e.g., "ssu", "lsu", "5S", "tRNA")
	 * @return LongHashSet containing the k-mers, or null if file not found
	 */
	private static LongHashSet loadLongKmersByType(int k, String prefix){
		String fname=Data.findPath("?"+prefix+"_"+k+"mers.fa", true);
		if(!new File(fname).exists()){
			fname=fname+".gz";
			if(!new File(fname).exists()){
				System.err.println("Can't find "+fname);
				return null;
			}
		}
		LongHashSet set=loadLongKmers(fname, k);
		return set;
	}
	
	/**
	 * Loads k-mers from a FASTA file into a LongHashSet.
	 * Processes all sequences and extracts all k-mers of the specified length.
	 *
	 * @param fname The input FASTA filename
	 * @param k The k-mer length
	 * @return LongHashSet containing all k-mers from the file
	 */
	private static LongHashSet loadLongKmers(String fname, int k){//TODO: Consider making this a LongHashSet.  No reason not to...
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FA, null, false, false);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start(); //Start the stream
//		if(verbose){outstream.println("Started cris");}
		
		LongHashSet set=new LongHashSet(1000);
		ListNum<Read> ln=cris.nextList();
		while(ln!=null && ln.size()>0){
			processList(ln, set, k);
			cris.returnList(ln);
			ln=cris.nextList();
		}
		if(ln!=null){cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());}
		ReadWrite.closeStream(cris);
		return set;
	}
	
	/**
	 * Extracts all k-mers from a list of reads and adds them to the set.
	 * Uses rolling hash to efficiently generate k-mers from sequence data.
	 *
	 * @param ln List of reads to process
	 * @param set The LongHashSet to add k-mers to
	 * @param k The k-mer length
	 * @return The modified LongHashSet
	 */
	private static LongHashSet processList(ListNum<Read> ln, LongHashSet set, int k){
		final long mask=~((-1L)<<(2*k));
		for(Read r : ln){
			final byte[] bases=r.bases;
			long kmer=0;
			int len=0;
			for(byte b : bases){
				final int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					len++;
					kmer=((kmer<<2)|num)&mask;
					if(len>=k){
						set.add(kmer);
					}
				}else{
					len=0;
				}
			}
		}
		return set;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Consensus Sequence      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads consensus sequences for all enabled ribosomal RNA types from reference files.
	 * Thread-safe method that optionally filters out mitochondrial and chloroplast sequences.
	 * Only loads once per execution to avoid redundant file I/O.
	 *
	 * @param removeMito Whether to remove mitochondrial sequences
	 * @param removeChloro Whether to remove chloroplast sequences
	 */
	public static synchronized void loadConsensusSequenceFromFile(boolean removeMito, boolean removeChloro){
		if(loadedConsensusSequence){return;}
//		assert(r16SSequence==null);
		if(load16SSequence){r16SSequence=loadConsensusSequenceType("16S", removeMito, removeChloro);}
		if(load18SSequence){r18SSequence=loadConsensusSequenceType("18S", removeMito, removeChloro);}
		if(load23SSequence){r23SSequence=loadConsensusSequenceType("23S", removeMito, removeChloro);}
		if(load5SSequence){r5SSequence=loadConsensusSequenceType("5S", removeMito, removeChloro);}
		if(loadtRNASequence){trnaSequence=loadConsensusSequenceType("tRNA", removeMito, removeChloro);}
		loadedConsensusSequence=true;
	}
	
	/**
	 * Loads consensus sequences for a specific RNA type from reference files.
	 * Searches for FASTQ files first, then FASTA files as fallback.
	 *
	 * @param prefix The RNA type prefix (e.g., "16S", "23S", "5S", "18S", "tRNA")
	 * @param removeMito Whether to remove mitochondrial sequences
	 * @param removeChloro Whether to remove chloroplast sequences
	 * @return Array of consensus reads, or null if file not found
	 */
	public static Read[] loadConsensusSequenceType(String prefix, boolean removeMito, boolean removeChloro){
		String fname=null;
		fname=Data.findPath("?"+prefix+"_consensus_sequence.fq", false);
		if(fname!=null && (fname.endsWith(".jar") || new File(fname).exists())){
			fname=Tools.fixExtension(fname);
		}else{
			fname=Data.findPath("?"+prefix+"_consensus_sequence.fa", true);
			fname=Tools.fixExtension(fname);
			if(!fname.endsWith(".jar") && !new File(fname).exists()){
				System.err.println("Can't find "+fname);
				return null;
			}
		}
		Read[] array=loadConsensusSequence(fname);
		if(removeMito){array=stripOrganelle(array, "mito");}
		if(removeChloro){array=stripOrganelle(array, "plastid");}
		return array;
	}
	
	/**
	 * Loads consensus sequences from a FASTA or FASTQ file.
	 * @param fname The input filename
	 * @return Array of reads from the file
	 */
	private static Read[] loadConsensusSequence(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FA, null, false, false);
		Read[] array=ReadInputStream.toReadArray(ff, -1);
		return array;
	}
	
	/**
	 * Removes organellar sequences from a read array based on header prefix.
	 * Filters out reads whose IDs start with the specified organelle key.
	 *
	 * @param array The input read array
	 * @param key The organelle prefix to filter (e.g., "mito", "plastid")
	 * @return Filtered array with organellar sequences removed
	 */
	private static Read[] stripOrganelle(Read[] array, String key){
		int removed=0;
		for(int j=0; j<array.length; j++){
			if(array[j].id.toLowerCase().startsWith(key)) {
				array[j]=null;
				removed++;
			}
		}
		if(removed>0){array=Tools.condenseStrict(array);}
		return array;
	}
	
	/*--------------------------------------------------------------*/
	
	public static final int CDS=0, tRNA=1, r16S=2, r23S=3, r5S=4, r18S=5, r28S=6, RNA=7;
	public static String[] typeStrings=new String[] {"CDS", "tRNA", "16S", "23S", "5S", "18S", "28S", "RNA"};
	public static String[] typeStrings2=new String[] {"CDS", "tRNA", "rRNA", "rRNA", "rRNA", "rRNA", "rRNA", "RNA"};
	public static String[] specialTypeStrings=new String[] {null, "tRNA", "16S", "23S", "5S", "18S", "28S", null};
	/**
	 * Determines if a type string represents a special RNA type (not CDS).
	 * @param type The type string to check
	 * @return true if this is a special RNA type, false otherwise
	 */
	public static boolean isSpecialType(String type){
		if(type==null){return false;}
		for(String s : specialTypeStrings){
			if(type.equalsIgnoreCase(s)){return true;}
		}
		return false;
	}

	/** K-mer length for inner RNA region analysis */
	public static int kInnerRNA=6;
	/** K-mer length for RNA start region analysis */
	public static int kStartRNA=3;
	/** K-mer length for RNA stop region analysis */
	public static int kStopRNA=3;

	/** Long k-mer length for small subunit ribosomal RNA (16S/18S) analysis */
	public static int kLongSSU=15;
	/** Long k-mer length for large subunit ribosomal RNA (23S/28S) analysis */
	public static int kLongLSU=15;
	/** Long k-mer length for 5S ribosomal RNA analysis */
	public static int kLong5S=15;
	/** Long k-mer length for transfer RNA analysis */
	public static int kLongTRna=15;
	
	/** Minimum identity threshold for 16S ribosomal RNA classification */
	public static float min16SIdentity=0.62f;
	/** Minimum identity threshold for 23S ribosomal RNA classification */
	public static float min23SIdentity=0.60f;
	/** Minimum identity threshold for 5S ribosomal RNA classification */
	public static float min5SIdentity=0.60f;
	/** Minimum identity threshold for 18S ribosomal RNA classification */
	public static float min18SIdentity=0.60f;
	
	/** Start position tolerance for small subunit ribosomal RNA detection */
	static int ssuStartSlop=200;
	/** Stop position tolerance for small subunit ribosomal RNA detection */
	static int ssuStopSlop=0;
	/** Start position tolerance for large subunit ribosomal RNA detection */
	static int lsuStartSlop=220;
	/** Stop position tolerance for large subunit ribosomal RNA detection */
	static int lsuStopSlop=0;
	/** Start position tolerance for 5S ribosomal RNA detection */
	static int r5SStartSlop=50;
	/** Stop position tolerance for 5S ribosomal RNA detection */
	static int r5SStopSlop=50;

	/** Whether to call coding sequences during gene annotation */
	public static boolean callCDS=true;
	/** Whether to call transfer RNA genes during gene annotation */
	public static boolean calltRNA=true;
	/** Whether to call 16S ribosomal RNA genes during gene annotation */
	public static boolean call16S=true;
	/** Whether to call 23S ribosomal RNA genes during gene annotation */
	public static boolean call23S=true;
	/** Whether to call 5S ribosomal RNA genes during gene annotation */
	public static boolean call5S=true;
	/** Whether to call 18S ribosomal RNA genes during gene annotation */
	public static boolean call18S=false;

	/** K-mer set for small subunit ribosomal RNA detection */
	public static LongHashSet ssuKmers=null;
	/** K-mer set for large subunit ribosomal RNA detection */
	public static LongHashSet lsuKmers=null;
	/** K-mer set for 5S ribosomal RNA detection */
	public static LongHashSet r5SKmers=null;
	/** K-mer set for transfer RNA detection */
	public static LongHashSet trnaKmers=null;

	/** Consensus sequences for transfer RNA reference */
	public static Read[] trnaSequence=null;
	/** Consensus sequences for 16S ribosomal RNA reference */
	public static Read[] r16SSequence=null;
	/** Consensus sequences for 23S ribosomal RNA reference */
	public static Read[] r23SSequence=null;
	/** Consensus sequences for 5S ribosomal RNA reference */
	public static Read[] r5SSequence=null;
	/** Consensus sequences for 18S ribosomal RNA reference */
	public static Read[] r18SSequence=null;

	/** Whether to process the forward DNA strand during gene annotation */
	public static boolean PROCESS_PLUS_STRAND=true;
	/** Whether to process the reverse DNA strand during gene annotation */
	public static boolean PROCESS_MINUS_STRAND=true;

	/** Whether to load small subunit ribosomal RNA k-mers */
	public static boolean loadSSUkmers=true;
	/** Whether to load large subunit ribosomal RNA k-mers */
	public static boolean loadLSUkmers=true;
	/** Whether to load 5S ribosomal RNA k-mers */
	public static boolean load5Skmers=true;
	/** Whether to load transfer RNA k-mers */
	public static boolean loadtRNAkmers=true;
	/** Flag indicating whether long k-mer sets have been loaded */
	private static boolean loadedLongKmers=false;

	/** Whether to load transfer RNA consensus sequences */
	public static boolean loadtRNASequence=false;
	/** Whether to load 16S ribosomal RNA consensus sequences */
	public static boolean load16SSequence=true;
	/** Whether to load 23S ribosomal RNA consensus sequences */
	public static boolean load23SSequence=true;
	/** Whether to load 5S ribosomal RNA consensus sequences */
	public static boolean load5SSequence=true;
	/** Whether to load 18S ribosomal RNA consensus sequences */
	public static boolean load18SSequence=true;
	/** Flag indicating whether consensus sequences have been loaded */
	private static boolean loadedConsensusSequence=false;
	
}
