package assemble;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import jgi.AssemblyStats2;
import jgi.ReformatReads;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import ukmer.Kmer;

/**
 * Wrapper for Tadpole assembler that tests multiple k-mer lengths to find optimal assembly.
 * Iterates through specified k-mer values, runs Tadpole for each, and selects the k-mer
 * that produces the best assembly metrics (L50, L90, max contig length). Optionally
 * performs bisection search or expansion to refine the optimal k-mer value.
 *
 * @author Brian Bushnell
 * @date Oct 15, 2015
 */
public class TadpoleWrapper {
	

	public static void main(String[] args){
		process(args);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
		
	
	/**
	 * Main processing method that runs Tadpole with multiple k-mer values and selects the best.
	 * Parses arguments to extract k-mer values, runs Tadpole for each k-mer, compares assembly
	 * statistics, and optionally performs bisection/expansion to refine the optimal k-mer.
	 * Handles output file naming and cleanup of intermediate assemblies if requested.
	 *
	 * @param args Command-line arguments containing k-mer values, output paths, and options
	 * @return The optimal k-mer value that produced the best assembly metrics
	 */
	public static int process(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		HashSet<Integer> set=new HashSet<Integer>();
		ArrayList<String> argList=new ArrayList<String>();
		String contigsName="contigs%.fa";
		String outFinal=null;
		bestAssembly=null;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("k")){
				assert(b!=null) : "Bad parameter: "+arg;
				for(String s2 : b.split(",")){
					int x=Integer.parseInt(s2);
					x=Kmer.getKbig(x);
					set.add(x);
				}
			}else if(a.equals("out")){
				assert(b!=null) : "Bad parameter: "+arg;
				contigsName=b;
				assert(b.contains("%")) : "Output name must contain % symbol.";
			}else if(a.equals("outfinal")){
				outFinal=b;
			}else if(a.equals("quitearly")){
				quitEarly=Parse.parseBoolean(b);
			}else if(a.equals("delete")){
				delete=Parse.parseBoolean(b);
			}else if(a.equals("bisect")){
				bisect=Parse.parseBoolean(b);
			}else if(a.equals("expand") /*|| a.equals("extend")*/){
				expand=Parse.parseBoolean(b);
			}else{
				argList.add(arg);
			}
		}
		
		if(outFinal==null && contigsName!=null && contigsName.indexOf('%')<0){
			outFinal=contigsName;
			contigsName="contigs_k%.fa";
		}
		
		if(set.isEmpty()){
			kmers=new int[] {31};
		}
		
		{
			kmers=new int[set.size()];
			int i=0;
			for(Integer x : set){
				kmers[i]=x;
				i++;
			}
			Arrays.sort(kmers);
		}
		
		ArrayList<Record> records=new ArrayList<Record>();
		
		argList.add("");
		argList.add("");
		
		int best=0;
		Record bestRecord=null;
		
		for(int i=0; i<kmers.length; i++){
			int k=kmers[i];
			argList.set(argList.size()-2, "k="+k);
			argList.set(argList.size()-1, "out="+contigsName.replace("%", ""+k));
			String[] args2=argList.toArray(new String[0]);
			System.gc();
			Tadpole.main(args2);
			
			Record r=new Record(k, AssemblyStats2.lastL50, AssemblyStats2.lastL90,
					AssemblyStats2.lastSize, AssemblyStats2.lastContigs, AssemblyStats2.lastMaxContig);
			records.add(r);
			
			if(bestRecord==null){
				bestRecord=r;
				best=i;
			}else{
				int x=bestRecord.compareTo(r);
				if(x<0){
					bestRecord=r;
					best=i;
				}else if(x>0 && quitEarly){
					System.err.println("Metrics stopped improving; quit early.");
					break;
				}
			}
		}
		
		if(expand || bisect){
			best=bisect(records, best, argList, contigsName, true);
			bestRecord=records.get(best);
		}
		
		assert(bestRecord!=null);
		int bestK=bestRecord.k;
		bestAssembly=(contigsName.replace("%", ""+bestK));
		System.err.println("Recommended K:\t"+bestK);
		
		if(outFinal!=null){
			File f=new File(bestAssembly);
			boolean success=(Shared.WINDOWS ? false : f.renameTo(new File(outFinal)));
			if(!success && f.exists()){
				try {
					ReformatReads.main(new String[] {"in="+bestAssembly, "out="+outFinal, "ow"});
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(f.exists()){f.delete();}
				
//				try {
//					Files.move(f.toPath(), new File(outFinal).toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
//				} catch (IOException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
			}
		}
		
		if(delete){
			for(Record r : records){
				if(r.k!=bestK){
					String s=contigsName.replace("%", ""+r.k);
					File f=new File(s);
					if(f.exists()){
						f.delete();
					}
				}
			}
		}
		
		return bestK;
	}
	
	/**
	 * Expands search to smaller k-mer values when the best k-mer is at the left boundary.
	 * Tests a k-mer that is 70% of the current best k-mer value. If the new k-mer produces
	 * better results, recursively continues expanding leftward.
	 *
	 * @param list List of assembly records for different k-mer values
	 * @param best Index of current best assembly in the list
	 * @param argList Command-line arguments template for Tadpole execution
	 * @param contigsName Output filename template with % placeholder for k-mer value
	 * @param recur Whether to recursively continue expansion if improvement found
	 * @return Updated index of the best assembly after expansion
	 */
	static int expandLeft(ArrayList<Record> list, int best, ArrayList<String> argList, String contigsName, boolean recur){
		assert(best==0);
		Record mid=list.get(best);
		int k1=Kmer.getKbig((int)(mid.k*0.7f));
		outstream.println("ExpandLeft: best="+mid.k+", next="+k1);
		if(k1>=mid.k || k1<1){return best;}
		
		{
			argList.set(argList.size()-2, "k="+k1);
			argList.set(argList.size()-1, "out="+contigsName.replace("%", ""+k1));
			String[] args2=argList.toArray(new String[0]);
			System.gc();
			Tadpole.main(args2);

			Record r=new Record(k1, AssemblyStats2.lastL50, AssemblyStats2.lastL90,
					AssemblyStats2.lastSize, AssemblyStats2.lastContigs, AssemblyStats2.lastMaxContig);
			list.add(best, r);
			best++;
			
			if(r.compareTo(mid)>0){
				return expandLeft(list, best-1, argList, contigsName, true);
			}
		}
		return best;
	}
	
	/**
	 * Expands search to larger k-mer values when the best k-mer is at the right boundary.
	 * Tests a k-mer that is min(current+40, current*1.25+10). If the new k-mer produces
	 * better results, recursively continues expanding rightward.
	 *
	 * @param list List of assembly records for different k-mer values
	 * @param best Index of current best assembly in the list
	 * @param argList Command-line arguments template for Tadpole execution
	 * @param contigsName Output filename template with % placeholder for k-mer value
	 * @param recur Whether to recursively continue expansion if improvement found
	 * @return Updated index of the best assembly after expansion
	 */
	static int expandRight(ArrayList<Record> list, int best, ArrayList<String> argList, String contigsName, boolean recur){
		assert(best==list.size()-1);
		Record mid=list.get(best);
		int k1=Kmer.getKbig(Tools.min(mid.k+40, (int)(10+mid.k*1.25f)));
		outstream.println("ExpandRight: best="+mid.k+", next="+k1);
		if(k1<=mid.k){return best;}
		
		{
			argList.set(argList.size()-2, "k="+k1);
			argList.set(argList.size()-1, "out="+contigsName.replace("%", ""+k1));
			String[] args2=argList.toArray(new String[0]);
			System.gc();
			Tadpole.main(args2);

			Record r=new Record(k1, AssemblyStats2.lastL50, AssemblyStats2.lastL90,
					AssemblyStats2.lastSize, AssemblyStats2.lastContigs, AssemblyStats2.lastMaxContig);
			list.add(r);
			
			if(r.compareTo(mid)>0){
				return expandRight(list, best+1, argList, contigsName, true);
			}
		}
		return best;
	}
	
	/**
	 * Performs bisection search to find optimal k-mer between tested values.
	 * If expand mode is enabled, first expands boundaries. Then tests intermediate k-mer
	 * values between the best k-mer and its neighbors. If better results are found,
	 * recursively continues bisection search to refine the optimal k-mer.
	 *
	 * @param list List of assembly records for different k-mer values
	 * @param best Index of current best assembly in the list
	 * @param argList Command-line arguments template for Tadpole execution
	 * @param contigsName Output filename template with % placeholder for k-mer value
	 * @param recur Whether to recursively continue bisection if improvement found
	 * @return Updated index of the best assembly after bisection
	 */
	static int bisect(ArrayList<Record> list, int best, ArrayList<String> argList, String contigsName, boolean recur){
		
		if(expand){
			if(best==0){best=expandLeft(list, best, argList, contigsName, recur);}
			if(best==list.size()-1){best=expandRight(list, best, argList, contigsName, recur);}
		}
		
		if(!bisect || best==0 || best==list.size()-1){return best;}
		
		Record left=list.get(best-1);
		Record mid=list.get(best);
		Record right=list.get(best+1);
		
		int k1=Kmer.getKbig((left.k+mid.k+1)/2);
		int k2=Kmer.getKbig((mid.k+right.k+1)/2);
		
		outstream.println("Bisect: left="+left.k+", k1="+k1+", best="+mid.k+", k2="+k2+", right="+right.k);
		
		if(k1==left.k || k1==mid.k){return best;}
		if(k2==mid.k || k2==right.k){return best;}
		
//		Record r1=null, r2=null;
		
		{
			argList.set(argList.size()-2, "k="+k1);
			argList.set(argList.size()-1, "out="+contigsName.replace("%", ""+k1));
			String[] args2=argList.toArray(new String[0]);
			System.gc();
			Tadpole.main(args2);

			Record r=new Record(k1, AssemblyStats2.lastL50, AssemblyStats2.lastL90,
					AssemblyStats2.lastSize, AssemblyStats2.lastContigs, AssemblyStats2.lastMaxContig);
			list.add(best, r);
			best++;
			
			if(r.compareTo(mid)>0){
				return bisect(list, best-1, argList, contigsName, true);
			}
		}

		{
			argList.set(argList.size()-2, "k="+k2);
			argList.set(argList.size()-1, "out="+contigsName.replace("%", ""+k2));
			String[] args2=argList.toArray(new String[0]);
			System.gc();
			Tadpole.main(args2);

			Record r=new Record(k2, AssemblyStats2.lastL50, AssemblyStats2.lastL90,
					AssemblyStats2.lastSize, AssemblyStats2.lastContigs, AssemblyStats2.lastMaxContig);
			list.add(best+1, r);
			
			if(r.compareTo(mid)>0){
				return bisect(list, best+1, argList, contigsName, true);
			}
		}
		
		if(recur){
			best=bisect(list, best, argList, contigsName, false);
		}
		
		return best;
	}
	
	/**
	 * Stores assembly statistics and quality metrics for a specific k-mer value.
	 * Used to compare assemblies from different k-mer lengths and determine the optimal k-mer.
	 * Implements Comparable to rank assemblies by quality using multiple metrics with tolerances.
	 */
	private static class Record implements Comparable<Record>{
		
		public Record(){}
		
		public Record(int k_, long L50_, long L90_, long contigLen_, long contigs_, long maxContig_){
			k=k_;
			L50=L50_;
			L90=L90_;
			contigLen=contigLen_;
			contigs=contigs_;
			maxContig=maxContig_;
		}
		
		/**
		 * Compares assembly quality between two Record objects using multiple metrics.
		 * Uses hierarchical comparison with tolerance thresholds: first compares L50 and L90
		 * with 1% tolerance, then maximum contig length, then number of contigs (fewer is better),
		 * then uses tighter 0.2% tolerances, and finally prefers smaller k-mer values if tied.
		 *
		 * @param b The Record to compare against
		 * @return Negative if this record is worse, positive if better, or k-mer difference if tied
		 */
		@Override
		public int compareTo(Record b){
			if(b==null){return 1;}
			
//			if(L50==b.L50 && L90==b.L90 && contigLen==b.contigLen && contigs==b.contigs && maxContig==b.maxContig){
//				return 0;
//			}
			
			if(L50<b.L50*0.99){return -1;}
			if(L50>b.L50*1.01){return 1;}
			
			if(L90<b.L90*0.99){return -1;}
			if(L90>b.L90*1.01){return 1;}
			
			if(maxContig<b.maxContig*0.99){return -1;}
			if(maxContig>b.maxContig*1.01){return 1;}
			
			if(contigs>b.contigs){return -1;}
			
			if(L50<b.L50*0.998){return -1;}
			if(L50>b.L50*1.002){return 1;}
			
			if(L90<b.L90*0.998){return -1;}
			if(L90>b.L90*1.002){return 1;}
			
			if(maxContig<b.maxContig*0.998){return -1;}
			if(maxContig>b.maxContig*1.002){return 1;}
			
			return b.k-k;
		}
		
		int k;
		long L50;
		long L90;
		long contigLen;
		long contigs;
		long maxContig;
		
	}
	
	private static int[] kmers;
	private static boolean quitEarly=false;
	private static boolean delete=false;
	private static boolean expand=false;
	private static boolean bisect=false;
	private static String bestAssembly=null;
	
	private static PrintStream outstream=System.err;
	
}
