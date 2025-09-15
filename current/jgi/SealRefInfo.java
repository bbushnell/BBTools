package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLongArray;

import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Tools;
import structures.IntList;


/**
 * Handles duties related to tracking Seal's reference set and scaffolds.
 * @author BBushnell
 *
 */
public class SealRefInfo {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * @param outstream Typically stderr
	 */
	public SealRefInfo(ArrayList<String> refs_, String[] literals_, PrintStream outstream){
		
		refs=refs_;
		literals=literals_;
		
		{//Add a fake first scaffold so that the true first scaffold has ID 1, not 0.
			scaffoldNames.add("");
			scaffoldLengths.add(0);
			scaffoldKmers.add(0);
			scaffolds.add(null);
		}
		refNames.add(null);//Fake first reference also
		
		if(!refs.isEmpty()){
			ArrayList<String> temp=new ArrayList<String>();
			for(String s : refs){
				if(s==null){
					assert(false) : "Null reference file.";
				}else if(new File(s).exists()){
					Tools.getFileOrFiles(s, temp, true, false, false, false);
				}else{
					String fname=null;
					if("phix".equalsIgnoreCase(s)){
						fname=Data.findPath("?phix2.fa.gz");
					}else if("polya".equalsIgnoreCase(s) || "polyt".equalsIgnoreCase(s)){
						fname=Data.findPath("?polyA.fa.gz");
					}else if("lambda".equalsIgnoreCase(s)){
						fname=Data.findPath("?lambda.fa.gz");
					}else if("kapa".equalsIgnoreCase(s)){
						fname=Data.findPath("?kapatags.L40.fa");
					}else if("pjet".equalsIgnoreCase(s)){
						fname=Data.findPath("?pJET1.2.fa");
					}else if("mtst".equalsIgnoreCase(s)){
						fname=Data.findPath("?mtst.fa");
					}else if("adapters".equalsIgnoreCase(s)){
						fname=Data.findPath("?adapters.fa");
					}else if("truseq".equalsIgnoreCase(s)){
						fname=Data.findPath("?truseq.fa.gz");
					}else if("nextera".equalsIgnoreCase(s)){
						fname=Data.findPath("?nextera.fa.gz");
					}else if("artifacts".equalsIgnoreCase(s)){
						fname=Data.findPath("?sequencing_artifacts.fa.gz");
					}else if("crisprs".equalsIgnoreCase(s)){
						s=Data.findPath("?crisprs.fa.gz");
					}else{
						assert(false) : "Can't find reference file "+s;
					}
					temp.add(fname);
				}
			}
			refs.clear();
			refs.addAll(temp);
			refNames.addAll(temp);
		}
		
		if(literals!=null){refNames.add("literal");}
		refScafCounts=new int[refNames.size()];
		
		if(!Tools.testInputFiles(true, true, refs)){
			throw new RuntimeException("\nCan't read some reference files.\n");
		}
		
		if(refs.isEmpty() && literals==null){
			outstream.println("ERROR: No reference sequences specified.  Use the -da flag to run anyway.");
			assert(false) : "Please specify a reference.";
		}
				
		if(refs!=null){
			for(String s0 : refs){
				assert(s0!=null) : "Specified a null reference.";
				String s=s0.toLowerCase();
				assert(s==null || s.startsWith("stdin") || s.startsWith("standardin") || new File(s0).exists()) : "Can't find "+s0;
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	boolean isEmpty(){
		return refs.isEmpty() && (literals==null || literals.length<1);
	}
	
	int numScaffolds(){
		return Tools.max(0, scaffoldNames.size()-1);
	}
	
	/**
	 * Clear stored sequence data.
	 */
	public void unloadScaffolds(){
		if(scaffoldNames!=null && !scaffoldNames.isEmpty()){
			scaffoldNames.clear();
			scaffoldNames.trimToSize();
		}
		scaffoldReadCounts=null;
		scaffoldFragCounts=null;
		scaffoldBaseCounts=null;
		scaffoldLengths=null;
		scaffoldKmers=null;
		scaffolds=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Printing           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Write statistics on a per-reference basis.
	 */
	void writeRefStats(String fname, boolean ow, boolean useRefNames, boolean printNonZeroOnly,
			String in1, String in2, long readsIn){
		if(fname==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(fname, ow, false, false);
		tsw.start();
		
		/* Count mapped reads */
		long mapped=0;
		for(int i=0; i<scaffoldReadCounts.length(); i++){
			mapped+=scaffoldReadCounts.get(i);
		}
		
		final int numRefs=refNames.size();
		long[] refReadCounts=new long[numRefs];
		long[] refFragCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		long[] refLengths=new long[numRefs];
		
		for(int r=1, s=1; r<numRefs; r++){
			final int lim=s+(useRefNames ? 1 : refScafCounts[r]);
			while(s<lim){
				refReadCounts[r]+=scaffoldReadCounts.get(s);
				refFragCounts[r]+=scaffoldFragCounts.get(s);
				refBaseCounts[r]+=scaffoldBaseCounts.get(s);
				refLengths[r]+=scaffoldLengths.get(s);
				s++;
			}
		}
		
		/* Print header */
		tsw.print("#File\t"+in1+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mapped));
		tsw.print(Tools.format("#References\t%d\n",refNames.size()-1));
		tsw.print("#Name\tLength\tScaffolds\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\n");
		
		final float mult=1000000000f/Tools.max(1, mapped);
		
		/* Print data */
		for(int i=1; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long frags=refFragCounts[i];
			final long bases=refBaseCounts[i];
			final long len=refLengths[i];
			final int scafs=refScafCounts[i];
			final String name=ReadWrite.stripToCore(refNames.get(i));
			final double invlen=1.0/Tools.max(1, len);
			final double mult2=mult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\n",name,len,scafs,bases,bases*invlen,reads,reads*mult2,frags,frags*mult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	/**
	 * Write statistics on a per-reference basis.
	 */
	void writeRefStats_BBSplitStyle(String fname, boolean ow, boolean useRefNames, 
			boolean printNonZeroOnly, long totalReads){
		if(fname==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(fname, ow, false, false);
		tsw.start();
		
		final int numRefs=refNames.size();
		long[] refReadCounts=new long[numRefs];
		long[] refBaseCounts=new long[numRefs];
		
		for(int r=1, s=1; r<numRefs; r++){
			final int lim=s+(useRefNames ? 1 : refScafCounts[r]);
			while(s<lim){
				refReadCounts[r]+=scaffoldReadCounts.get(s);
				refBaseCounts[r]+=scaffoldBaseCounts.get(s);
				s++;
			}
		}
		
		/* Print header */
		tsw.print("#name\t%unambiguousReads\tunambiguousMB\t%ambiguousReads\tambiguousMB\tunambiguousReads\tambiguousReads\n");

		final float rmult=100f/Tools.max(1, totalReads);
		
		/* Print data */
		for(int i=1; i<refNames.size(); i++){
			final long reads=refReadCounts[i];
			final long bases=refBaseCounts[i];
			final float unambigMB=bases*0.000001f;
			
			final long ambigReads=0; //TODO but not urgent
			final long ambigBases=0; //TODO but not urgent
			final float ambigMB=ambigBases*0.000001f;
			
			final String name=ReadWrite.stripToCore(refNames.get(i));

			final double unambigReadP=rmult*reads;
			final double ambigReadP=rmult*ambigReads;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\n",name,unambigReadP,unambigMB,ambigReadP,ambigMB,reads,ambigReads));
			}
		}
		tsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

//	/** Array of reference files from which to load kmers */
	final ArrayList<String> refs;
//	/** Array of literal strings from which to load kmers */
	final String[] literals;
	
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	final ArrayList<String> scaffoldNames=new ArrayList<String>();
	/** Names of reference files (refNames[0] is valid). */
	final ArrayList<String> refNames=new ArrayList<String>();
	/** Number of scaffolds per reference. */
	final int[] refScafCounts;
	
	/** scaffoldCounts[id] stores the number of reads with kmer matches to that scaffold */
	AtomicLongArray scaffoldReadCounts;
	/** scaffoldFragCounts[id] stores the number of fragments (reads or pairs) with kmer matches to that scaffold */
	AtomicLongArray scaffoldFragCounts;
	/** scaffoldBaseCounts[id] stores the number of bases with kmer matches to that scaffold */
	AtomicLongArray scaffoldBaseCounts;
	
	/** scaffoldLengths[id] stores the length of that scaffold */
	IntList scaffoldLengths=new IntList();
	/** scaffoldLengths[id] stores the number of kmers in that scaffold (excluding mutants) */
	IntList scaffoldKmers=new IntList();
	/** scaffolds[id] stores the number of kmers in that scaffold */
	ArrayList<byte[]> scaffolds=new ArrayList<byte[]>();
	
}
