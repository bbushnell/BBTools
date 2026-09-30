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
 * Prepares cached chromosome arrays and reference metadata for index construction.
 * Handles conversion of FASTA reference files to indexed chromosome arrays
 * and manages associated metadata files for the BBTools suite.
 * Actual k-mer blocks are built separately by the IndexMaker classes. Settings
 * and retained chromosome state are process-global; setup calls must be serialized.
 * Rebuilding may delete existing cache files when overwrite permits it.
 *
 * @author Brian Bushnell
 * @date Sep 25, 2013
 */
public class RefToIndex{

	/** Clears cached chromosome arrays from the last index build to free memory. */
	public static final void clear(){
		chromlist=null;
	}

	/**
	 * Constructs the file path for the genome summary file.
	 * Uses the configured genome root and the explicitly supplied build number.
	 * @param build The genome build number
	 * @return Path to the summary.txt file for this build
	 */
	public static String summaryLoc(int build){
		return Data.ROOT_GENOME+build+"/summary.txt";
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
	 * Prepares reference chromosome arrays from a FASTA file.
	 * Validates input file format, manages existing index cleanup, and delegates
	 * to FastaToChromArrays2 for chromosome array generation. Handles directory
	 * structure creation and logging.
	 *
	 * @param reference Path to input FASTA reference file
	 * @param build Build number for this reference genome
	 * @param sysout Output stream for status messages
	 * @param keylen Retained for caller compatibility; reference directory names do not depend on k
	 * @throws RuntimeException if reference file is invalid or directories cannot be created
	 */
	public static void makeIndex(String reference, int build, PrintStream sysout, int keylen){
		assert(reference!=null) : "Reference conversion requires a FASTA path or stdin name";
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

		// Use the same explicit build and configured roots as the downstream converter.
		// Deriving paths from Data.GENOME_BUILD can target a different cache at this API boundary.
		String dir=Data.ROOT_GENOME+build;
		final String base=Data.ROOT_REF;
		final String args=(Shared.COMMAND_LINE==null ? "null" : Arrays.toString(Shared.COMMAND_LINE));
		final String indexlog=base+"build"+build+"_"+
				(System.nanoTime()&Long.MAX_VALUE)+"."+((args==null ? (reference==null ? "null" : reference) : args).hashCode()&Integer.MAX_VALUE)+".log";
		final String sf=summaryLoc(build);
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
				//TODO: Probable bug - delete() failures are ignored below, so a rebuild can retain stale cache files.
				//Failure handling and partial-rebuild recovery are not repaired by this path/ZIPLEVEL change.
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

				dir=Data.ROOT_INDEX+build;
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

			final int oldzl=ReadWrite.ZIPLEVEL;
			ReadWrite.ZIPLEVEL=Tools.max(4, ReadWrite.ZIPLEVEL);
			try{
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

				String[] ftcaArgs=new String[]{reference, ""+build, "writeinthread=false", "genscaffoldinfo="+genScaffoldInfo, "retain", "waitforwriting=false",
						"gz="+(Data.CHROMGZ), "maxlen="+maxChromLen,
						"writechroms="+(!NODISK), "minscaf="+minScaf, "midpad="+midPad, "startpad="+startPad, "stoppad="+stopPad, "nodisk="+NODISK};

				chromlist=FastaToChromArrays2.main2(ftcaArgs);
			}finally{
				// A validation/conversion failure must not change later writers' compression level.
				ReadWrite.ZIPLEVEL=oldzl;
			}
		}

	}

	public static boolean AUTO_CHROMBITS=true;
	public static boolean LOG=false;
	public static boolean NODISK=false;
	/** Skips conversion without requiring a matching summary; optional logging may still write. */
	public static boolean FORCE_READ_ONLY=false;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean genScaffoldInfo=true;

	public static long maxChromLen=-1;

	/** minScaf: minimum scaffold length to keep; midPad: padding inserted between joined scaffolds;
	 * stopPad: padding bases added at chromosome ends; startPad: padding bases added at chromosome starts. (-1=default) */
	public static int minScaf=-1, midPad=-1, stopPad=-1, startPad=-1;
	public static int chrombits=-1;

	/** Retained arrays from the last successful conversion; clear() releases this reference. */
	public static ArrayList<ChromosomeArray> chromlist=null;

}
