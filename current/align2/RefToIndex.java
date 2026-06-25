package align2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;

import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.SummaryFile;
import shared.Shared;
import shared.Tools;

/**
 * Utility class for creating and managing reference genome indices.
 * Handles conversion of FASTA reference files to indexed chromosome arrays
 * and manages associated metadata files for the BBTools suite.
 *
 * @author Brian Bushnell
 * @date Sep 25, 2013
 */
public class RefToIndex {
	
	/** Clears cached chromosome arrays from the last index build to free memory. */
	public static final void clear(){
		chromlist=null;
	}
	
	/**
	 * Constructs the file path for the genome summary file.
	 * Transforms index directory path to genome directory and appends summary.txt.
	 * @param build The genome build number
	 * @return Path to the summary.txt file for this build
	 */
	public static String summaryLoc(int build){
		String s=IndexMaker4.fname(1, 1, 13, 1, build);
		String dir=new File(s).getParent();
		dir=dir.replace('\\', '/');
		dir=dir.replace("ref/index/", "ref/genome/");
		String sf=dir+"/summary.txt";
		return sf;
	}
	
	/**
	 * Constructs the file path for the bloom.serial file for a given genome build.
	 * @param build Genome build number
	 * @return Path to the bloom filter file
	 */
	public static String bloomLoc(int build){
		return Data.ROOT_INDEX+build+"/bloom.serial";
	}
	
	/**
	 * Creates reference genome index from FASTA file.
	 * Validates input file format, manages existing index cleanup, and delegates
	 * to FastaToChromArrays2 for chromosome array generation. Handles directory
	 * structure creation and logging.
	 *
	 * @param reference Path to input FASTA reference file
	 * @param build Build number for this reference genome
	 * @param sysout Output stream for status messages
	 * @param keylen K-mer length for indexing
	 * @throws RuntimeException if reference file is invalid or directories cannot be created
	 */
	public static void makeIndex(String reference, int build, PrintStream sysout, int keylen){
		assert(reference!=null);
		{
			File f=new File(reference);
			if(!f.exists() || !f.isFile() || !f.canRead()){
				if(!reference.startsWith("stdin")){
					throw new RuntimeException("Cannot read file "+f.getAbsolutePath());
				}
			}else{
				FileFormat ff=FileFormat.testInput(reference, FileFormat.FA, null, false, true, true);
				if(!ff.fasta()){
					throw new RuntimeException("Reference file is not in fasta format: "+reference+"\n"+ff);
				}
			}
		}

		//Comprehension (Eru 2026-06-24): makeIndex = one-time reference setup, NOT the hot align path. Validates the
		//FASTA (L72-82), skips rebuild if an up-to-date summary already exists (L94), else deletes stale genome+index
		//dirs (overwrite-gated) and delegates the actual build to FastaToChromArrays2.main2 (L174). dir here =
		//".../ref/index/<GENOME_BUILD>" (IndexMaker4.fname parent; ROOT_INDEX ends in "index/").
		String s=IndexMaker4.fname(1, 1, keylen, 1);
		String dir=new File(s).getParent();
		dir=dir.replace('\\', '/');
		//TODO: Possible bug [align2/RefToIndex#001] (LOW, in deep-research baseline) - the fixed -7 trim assumes a
		//SINGLE-DIGIT build: dir ends "index/<build>", and "index/"=6 + 1 digit = 7. For build>=10 it over-trims into
		//the path (e.g. "ref/index/12" -> "ref/i"). `base` is used only for the indexlog name (L90-91) and base.mkdirs()
		//(L97/L153), BOTH gated behind LOG (default false) -> only a misnamed log + a stray empty dir, no data
		//corruption. Reachable via build>=10 (no <10 guard; BBSplitter parses build with plain parseInt). Brian's call.
		final String base=dir.substring(0, dir.length()-7);
		final String args=(Shared.COMMAND_LINE==null ? "null" : Arrays.toString(Shared.COMMAND_LINE));
		final String indexlog=base+"build"+build+"_"+
				(System.nanoTime()&Long.MAX_VALUE)+"."+((args==null ? (reference==null ? "null" : reference) : args).hashCode()&Integer.MAX_VALUE)+".log";
		dir=dir.replace("ref/index/", "ref/genome/");
		String sf=dir+"/summary.txt";
		if(FORCE_READ_ONLY || (!NODISK && new File(sf).exists() && SummaryFile.compare(sf, reference))){
			//do nothing
			if(LOG && !NODISK){
				if(!new File(base).exists()){new File(base).mkdirs();}
				ReadWrite.writeString(new Date()+"\nFound an already-written genome for build "+build+".\n"+args+"\n", indexlog, true);
			}
			sysout.println("NOTE:\tIgnoring reference file because it already appears to have been processed.");
			sysout.println("NOTE:\tIf you wish to regenerate the index, please manually delete "+dir+"/summary.txt");
		}else{
			if(NODISK){}
			else{//Delete old data if present
				File f=new File(dir);
				if(f.exists()){
					File[] f2=f.listFiles();
					if(f2!=null && f2.length>0){
						if(overwrite || f2[0].getAbsolutePath().equals(new File(reference).getAbsolutePath())){
							sysout.println("NOTE:\tDeleting contents of "+dir+" because reference is specified and overwrite="+overwrite);
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nDeleting genome for build "+build+".\n"+args+"\n", indexlog, true);}
							for(File f3 : f2){
								if(f3.isFile()){
									String f3n=f3.getName();
									if((f3n.contains(".chrom") || f3n.endsWith(".txt") || f3n.endsWith(".txt.gz")) && !f3n.endsWith("list.txt")){
										f3.delete();
									}
								}
							}
						}else{
							sysout.println(Arrays.toString(f2));
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nFailed to overwrite genome for build "+build+".\n"+args+"\n", indexlog, true);}
							throw new RuntimeException("\nThere is already a reference at location '"+f.getAbsolutePath()+"'.  " +
									"Please delete it (and the associated index), or use a different build ID, " +
									"or remove the 'reference=' parameter from the command line, or set overwrite=true.");
						}
					}
				}

				dir=dir.replace("ref/genome/", "ref/index/");
				f=new File(dir);
				if(f.exists()){
					File[] f2=f.listFiles();
					if(f2!=null && f2.length>0){
						if(overwrite){
							sysout.println("NOTE:\tDeleting contents of "+dir+" because reference is specified and overwrite="+overwrite);
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nDeleting index for build "+build+".\n"+args+"\n", indexlog, true);}
							for(File f3 : f2){
								if(f3.isFile()){f3.delete();}
							}
						}else{
							if(LOG && !NODISK){ReadWrite.writeString(new Date()+"\nFailed to overwrite index for build "+build+".\n"+args+"\n", indexlog, true);}
							throw new RuntimeException("\nThere is already an index at location '"+f.getAbsolutePath()+"'.  " +
									"Please delete it, or use a different build ID, or remove the 'reference=' parameter from the command line.");
						}
					}
				}
			}

			if(!NODISK){
				sysout.println("Writing reference.");
				if(LOG && !NODISK){
					if(!new File(base).exists()){new File(base).mkdirs();}
					ReadWrite.writeString(new Date()+"\nWriting genome for build "+build+".\n"+args+"\n", indexlog, true);
				}
			}

			int oldzl=ReadWrite.ZIPLEVEL;
			ReadWrite.ZIPLEVEL=Tools.max(4, ReadWrite.ZIPLEVEL);

			//assert(false) : "minScaf="+minScaf+", midPad="+midPad+", maxChromLen="+maxChromLen+
			//		", startPad="+startPad+", stopPad="+stopPad+", FastaToChromArrays2.END_PADDING="+FastaToChromArrays2.END_PADDING;
			
			//TODO: Possible bug [align2/RefToIndex#002] (LOW, in deep-research baseline) - the else branch (reached
			//only when maxChromLen unset AND AUTO_CHROMBITS=false, i.e. user set cbits=N) computes
			//(1L<<(31-chrombits))-200000 with NO Tools.max(1,..) clamp, so it goes NEGATIVE for chrombits>=14
			//(14 -> -68928). A chromosome length can't be negative. Reachable via cbits>=14 (legal range; auto caps at
			//16). With -ea it dies loud at FastaToChromArrays2's `assert(len>0..)` maxlen guard (crash-loud, acceptable);
			//with -ea off it throws "scaffold exceeds maximum" or yields 0 chroms. Advanced-param-gated -> LOW. Brian's call.
			maxChromLen=maxChromLen>0 ? maxChromLen : AUTO_CHROMBITS ? FastaToChromArrays2.MAX_LENGTH : ((1L<<(31-(chrombits<0 ? 2 : chrombits)))-200000);
			minScaf=minScaf>-1 ? minScaf : FastaToChromArrays2.MIN_SCAFFOLD;
			midPad=midPad>-1 ? midPad : FastaToChromArrays2.MID_PADDING;
			startPad=startPad>-1 ? startPad : FastaToChromArrays2.START_PADDING;
			stopPad=stopPad>-1 ? stopPad : FastaToChromArrays2.END_PADDING;
			
			String[] ftcaArgs=new String[] {reference, ""+build, "writeinthread=false", "genscaffoldinfo="+genScaffoldInfo, "retain", "waitforwriting=false",
					"gz="+(Data.CHROMGZ), "maxlen="+maxChromLen,
					"writechroms="+(!NODISK), "minscaf="+minScaf, "midpad="+midPad, "startpad="+startPad, "stoppad="+stopPad, "nodisk="+NODISK};
			
			chromlist=FastaToChromArrays2.main2(ftcaArgs);

			ReadWrite.ZIPLEVEL=oldzl;
		}

	}

	public static boolean AUTO_CHROMBITS=true;
	public static boolean LOG=false;
	public static boolean NODISK=false;
	public static boolean FORCE_READ_ONLY=false;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean genScaffoldInfo=true;
	
	public static long maxChromLen=-1;
	
	//align2/RefToIndex#003 DOC FIXED (Eru 2026-06-24): the three stacked javadocs here were orphaned (Java attaches
	//only the LAST javadoc to a declaration) AND misordered relative to the fields (minScaf was undocumented).
	//Corrected to the actual declaration order; all default -1 means "use the FastaToChromArrays2 default".
	/** minScaf: minimum scaffold length to keep; midPad: padding inserted between joined scaffolds;
	 * stopPad: padding bases added at chromosome ends; startPad: padding bases added at chromosome starts. (-1=default) */
	public static int minScaf=-1, midPad=-1, stopPad=-1, startPad=-1;
	public static int chrombits=-1;
//	public static int minScaf=FastaToChromArrays2.MIN_SCAFFOLD;
//	public static int midPad=FastaToChromArrays2.MID_PADDING;
//	public static int startPad=FastaToChromArrays2.START_PADDING;
//	public static int stopPad=FastaToChromArrays2.END_PADDING;
	
	public static ArrayList<ChromosomeArray> chromlist=null;
	
}
