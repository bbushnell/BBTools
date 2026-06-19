package prok;

import java.io.File;

import dna.AminoAcid;
import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import parse.Parse;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ListNum;

/**
 * Abstract base class containing static utilities and configuration for prokaryotic gene calling.
 * Centralizes ribosomal RNA (16S, 18S, 23S, 5S), tRNA, and CDS calling parameters, including slop, identity, and k-mer settings.
 * Handles consensus sequence and k-mer loading, parameter parsing, and type-specific helper methods for gene annotation tools.
 * @author Brian Bushnell
 * @date 2013
 */
public abstract class ProkObject {
	
	/** Parses one prok-specific argument (a=key lowercased, b=value): rRNA slop/identity thresholds, plus/minus strand
	 * flags, per-type sequence/kmer load toggles, and klong* kmer lengths. @return true if {@code a} was recognized. */
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
		//TODO: Possible bug [prok/ProkObject#002] - the "klongtrna" branch just below is DEAD: it duplicates the identical
		//branch ~10 lines up (in the load*kmers group), which catches it first. Harmless (both set kLongTRna) but redundant;
		//recommend removing one - Brian's call which. Copy-paste-drift fingerprint.
		}else if(a.equalsIgnoreCase("klongtrna")){
			kLongTRna=Integer.parseInt(b);
		}
		
		else{
			return false;
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Whether this element type is enabled for calling (its call-flag): callCDS/call16S/.../calltRNA by type.
	 * Unlisted types (r28S, RNA) default to true. Compare callType(), which instead asserts on unlisted types. */
	public static boolean processType(int type){
		return (type==CDS ? callCDS : type==r16S ? call16S : type==r23S ? call23S : type==r18S ? call18S : type==r5S ? call5S : type==tRNA ? calltRNA : true);
	}
	
	/** Allowed start-coordinate wiggle for an RNA type after alignment (SSU types 16S/18S share ssuStartSlop; 23S→lsu;
	 * 5S→r5S). Non-RNA types get 9999 (effectively unbounded — they don't use alignment slop). */
	public static int startSlop(int type) {
		int slop=(type==r16S ? ssuStartSlop : type==r23S ? lsuStartSlop : type==r18S ? ssuStartSlop : type==r5S ? r5SStartSlop : 9999);
		return slop;
	}
	
	/** Allowed stop-coordinate wiggle for an RNA type after alignment (SSU 16S/18S share ssuStopSlop; 23S→lsu; 5S→r5S).
	 * Non-RNA types get 9999 (unbounded). */
	public static int stopSlop(int type) {
		int slop=(type==r16S ? ssuStopSlop : type==r23S ? lsuStopSlop : type==r18S ? ssuStopSlop : type==r5S ? r5SStopSlop : 9999);
		return slop;
	}
	
	/** Minimum alignment identity to accept an rRNA hit, per type. Note: NOT SSU-shared — 16S and 18S keep separate
	 * thresholds (min16SIdentity vs min18SIdentity), unlike the slop/kmer settings. Non-rRNA types get 0. */
	public static float minID(int type) {
		float minIdentity=(type==r16S ? min16SIdentity : type==r23S ? min23SIdentity : type==r18S ? min18SIdentity : type==r5S ? min5SIdentity : 0);
		return minIdentity;
	}
	
	/** Consensus reference sequences for an rRNA type (used for alignment), or null for types without one (CDS/tRNA). */
	public static Read[] consensusReads(int type) {
		Read[] consensusReads=(type==r16S ? r16SSequence : type==r23S ? r23SSequence : type==r18S ? r18SSequence : type==r5S ? r5SSequence : null);
		return consensusReads;
	}
	
	/** Long-kmer set for an rRNA/tRNA type, or null if none loaded. SSU sharing: 16S and 18S both map to ssuKmers
	 * (homologous small-subunit rRNA), so SSU kmers are loaded/tuned once; 23S→lsuKmers, 5S→r5SKmers, tRNA→trnaKmers. */
	public static LongHashSet kmerSet(int type) {
		LongHashSet set=(type==tRNA ? trnaKmers : type==r16S ? ssuKmers : type==r23S ? lsuKmers : type==r5S ? r5SKmers : type==r18S ? ssuKmers : null);
		return set;
	}
	
	/** Long-kmer length for a type (SSU 16S/18S share kLongSSU; 23S→kLongLSU; 5S→kLong5S; tRNA→kLongTRna), or -1 if none. */
	public static int kLongLen(int type) {
		int kLongLen=(type==tRNA ? kLongTRna : type==r16S ? kLongSSU : type==r23S ? kLongLSU : type==r5S ? kLong5S : type==r18S ? kLongSSU : -1);
		return kLongLen;
	}
	
	/** Inverse of typeToFlag(): recovers the type from its single-bit flag via numberOfTrailingZeros+1 (domain: one set
	 * bit, for types 1..6). Currently uncalled. flagToType(0) would give 33 (numberOfTrailingZeros(0)==32) — see #001. */
	public static int flagToType(int flag) {
		return Integer.numberOfTrailingZeros(flag)+1;
	}
	
	/** Packs an RNA element type (tRNA=1 .. r28S=6) into a single bit: 1<<(type-1). Inverse of flagToType().
	 * Domain is 1..6 — CDS(0) has no type-flag (it uses the 3-frame coding bits instead). */
	public static byte typeToFlag(int type) {
		//TODO: [prok/ProkObject#001] FIXED - assert previously guarded only type<=6, not type>=1; typeToFlag(0)=1<<-1=(byte)0
		//(silent-wrong). Sole caller (GeneModel:483) is guarded by type>0 so CDS never reaches it (latent), but now it crashes
		//loud if it ever does. flagToType is currently uncalled.
		assert(type>=1 && type<=6) : type;
		return (byte)(1<<(type-1));
	}
	
	/** Like processType() but STRICTER: returns the call-flag for CDS/tRNA/16S/23S/5S/18S and asserts(false) on any other
	 * type (e.g. r28S, RNA) — so an unexpected type crashes loud here rather than silently defaulting true. */
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
	
	/** Loads the enabled long-kmer reference sets (SSU/LSU/5S/tRNA) from resource FASTA once; idempotent via
	 * loadedLongKmers. synchronized so concurrent callers load exactly once. */
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
	
	/** Resolves the {@code prefix}_{k}mers.fa[.gz] resource for one type and loads its kmers; returns null (with a stderr
	 * note) if the file is absent. */
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
	
	/** Streams a FASTA file and collects all length-k forward kmers into a LongHashSet. */
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
	
	/** Adds every length-k forward kmer from each read in {@code ln} to {@code set} (no reverse-complement canonicalization).
	 * An invalid base resets len; the kmer self-corrects after k more valid bases via the mask + len>=k guard. */
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
	
	/** Loads the enabled rRNA/tRNA consensus reference sequences once (idempotent via loadedConsensusSequence).
	 * @param removeMito drop mitochondrial entries @param removeChloro drop plastid entries. */
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
	
	/** Loads one type's consensus sequences, preferring {@code prefix}_consensus_sequence.fq then .fa (also resolves inside
	 * a .jar); optionally strips mito/plastid entries. Returns null if the file can't be found. */
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
	
	/** Reads all sequences from a FASTA file into a Read array. */
	private static Read[] loadConsensusSequence(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FA, null, false, false);
		Read[] array=ReadInputStream.toReadArray(ff, -1);
		return array;
	}
	
	/** Nulls out entries whose id starts with {@code key} (case-insensitive), then compacts; used to drop organellar
	 * (mito/plastid) consensus sequences. */
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
	/** True if {@code type} names a non-CDS special element (tRNA/16S/23S/5S/18S/28S). Null-safe: the null slots in
	 * specialTypeStrings (CDS, RNA) are skipped because String.equalsIgnoreCase(null) is false. */
	public static boolean isSpecialType(String type){
		if(type==null){return false;}
		for(String s : specialTypeStrings){
			if(type.equalsIgnoreCase(s)){return true;}
		}
		return false;
	}

	public static int kInnerRNA=6;
	public static int kStartRNA=3;
	public static int kStopRNA=3;

	public static int kLongSSU=15;
	public static int kLongLSU=15;
	public static int kLong5S=15;
	public static int kLongTRna=15;
	
	public static float min16SIdentity=0.62f;
	public static float min23SIdentity=0.60f;
	public static float min5SIdentity=0.60f;
	public static float min18SIdentity=0.60f;
	
	static int ssuStartSlop=200;
	static int ssuStopSlop=0;
	static int lsuStartSlop=220;
	static int lsuStopSlop=0;
	static int r5SStartSlop=50;
	static int r5SStopSlop=50;

	public static boolean callCDS=true;
	public static boolean calltRNA=true;
	public static boolean call16S=true;
	public static boolean call23S=true;
	public static boolean call5S=true;
	public static boolean call18S=false;

	public static LongHashSet ssuKmers=null;
	public static LongHashSet lsuKmers=null;
	public static LongHashSet r5SKmers=null;
	public static LongHashSet trnaKmers=null;

	public static Read[] trnaSequence=null;
	public static Read[] r16SSequence=null;
	public static Read[] r23SSequence=null;
	public static Read[] r5SSequence=null;
	public static Read[] r18SSequence=null;

	public static boolean PROCESS_PLUS_STRAND=true;
	public static boolean PROCESS_MINUS_STRAND=true;

	public static boolean loadSSUkmers=true;
	public static boolean loadLSUkmers=true;
	public static boolean load5Skmers=true;
	public static boolean loadtRNAkmers=true;
	private static boolean loadedLongKmers=false;

	public static boolean loadtRNASequence=false;
	public static boolean load16SSequence=true;
	public static boolean load23SSequence=true;
	public static boolean load5SSequence=true;
	public static boolean load18SSequence=true;
	private static boolean loadedConsensusSequence=false;
	
}
