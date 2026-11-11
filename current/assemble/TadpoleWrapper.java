package assemble;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.AssemblyStats2;
import jgi.ReformatReads;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.NullOutputStream;
import ukmer.Kmer;

/**
 * Runs serial Tadpole trials and selects an assembly using a length-based heuristic.
 * Requested K values are rounded to supported representations and deduplicated.
 * Optional boundary expansion and neighbor bisection explore additional values;
 * this search does not guarantee a globally optimal assembly.
 * Trials use Tadpole's normal policy/rescue dispatch and share its JVM state.
 * Each trial computes statistics from its final named output even when reporting
 * is disabled. Static search options also persist between
 * invocations; this class is not an independent or concurrent assembler instance.
 *
 * @author Brian Bushnell
 * @date Oct 15, 2015
 */
public class TadpoleWrapper{

	/** Runs the wrapper and closes its diagnostic stream, including on failure. */
	public static void main(String[] args){
		try{process(args);}
		finally{Shared.closeStream(outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Expands configuration, runs trials, and optionally publishes the selected file.
	 * Initial K trials run in ascending order; an empty K set defaults to 31.
	 * The output pattern must contain '%' with assertions enabled. Unrecognized wrapper
	 * options are forwarded to Tadpole. Final publication refuses input aliases and
	 * forbidden replacement, retaining the selected assembly when publication fails.
	 * With delete enabled, nonselected intermediate files other than the final destination
	 * are removed on a best-effort basis after selection and any requested publication.
	 * This method leaves stream closure to main.
	 *
	 * @param args Command-line arguments containing k-mer values, output paths, and options
	 * @return The supported K selected by the comparison heuristic
	 * @throws RuntimeException If the selected assembly cannot be published; a failed copy retains its source
	 */
	public static int process(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object(){}.getClass().getEnclosingClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}

		HashSet<Integer> set=new HashSet<Integer>();
		ArrayList<String> argList=new ArrayList<String>();
		String contigsName="contigs%.fa";
		String outFinal=null;
		bestAssembly=null;
		showStats=true;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 && !split[1].isEmpty() ? split[1] : null;

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
			}else if(a.equals("showstats") || a.equals("stats")){
				showStats=Parse.parseBoolean(b);
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
		}else{//Keep the default trial when no K was supplied; an empty array skipped every trial (#001 fixed).
			kmers=new int[set.size()];
			int i=0;
			for(Integer x : set){
				kmers[i]=x;
				i++;
			}
			Arrays.sort(kmers);
		}

		ArrayList<Record> records=new ArrayList<Record>();

		//The wrapper measures the final product once, after any scaffold/rescue stage.
		argList.add("showstats=f");
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

			Record r=readStatistics(k, contigsName.replace("%", ""+k));
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
			if(!f.isFile()){throw new RuntimeException("Selected assembly is missing: "+bestAssembly);}
			final boolean alreadyPublished=checkPublication(outFinal, f, argList);
			boolean success=alreadyPublished || (Shared.WINDOWS ? false : f.renameTo(new File(outFinal)));
			if(!success && f.exists()){
				try{
					ReformatReads.main(new String[] {"in="+bestAssembly, "out="+outFinal, "ow="+Tadpole.overwrite});
				}catch(Exception e){
					throw new RuntimeException("Could not publish assembly to "+outFinal+"; selected assembly retained at "+bestAssembly, e);
				}
				if(f.exists() && !f.delete()){throw new RuntimeException("Assembly copied, but source could not be removed: "+bestAssembly);}

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
					if(f.exists() && (outFinal==null || !sameFile(f, new File(outFinal)))){
						f.delete();
					}
				}
			}
		}

		return bestK;
	}

	/** Computes current final-output metrics independently of report visibility.
	 * Native reporting is disabled for every trial, so ordinary assembly still makes
	 * one statistics pass. The input must be a named, completed sequence file.
	 * Temporary report/compression settings are restored even if analysis fails;
	 * last* metrics deliberately remain available to existing serial callers.
	 */
	private static Record readStatistics(final int k, final String filename){
		if(FileFormat.isStdio(filename) || !new File(filename).isFile()){
			throw new IllegalArgumentException("Wrapper statistics require a completed named assembly: "+filename);
		}
		final int format=AssemblyStats2.FORMAT;
		final boolean speed=AssemblyStats2.showspeed, extended=AssemblyStats2.printExtended;
		final boolean l90=AssemblyStats2.printL90, pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final PrintStream diagnostics=AssemblyStats2.outstream;
		final PrintStream quiet=showStats ? null : new PrintStream(new NullOutputStream());
		try{
			if(!showStats){
				AssemblyStats2.FORMAT=0;
				AssemblyStats2.showspeed=false;
				AssemblyStats2.outstream=quiet;
			}else{outstream.println();}
			//Use the instance API; the diagnostic stream is borrowed, not ours to close.
			new AssemblyStats2(new String[]{"in="+filename, "printextended"}).process();
			return new Record(k, AssemblyStats2.lastL50, AssemblyStats2.lastL90,
					AssemblyStats2.lastSize, AssemblyStats2.lastContigs, AssemblyStats2.lastMaxContig);
		}finally{
			AssemblyStats2.FORMAT=format;
			AssemblyStats2.showspeed=speed;
			AssemblyStats2.printExtended=extended;
			AssemblyStats2.printL90=l90;
			AssemblyStats2.outstream=diagnostics;
			ReadWrite.USE_PIGZ=pigz;
			ReadWrite.USE_UNPIGZ=unpigz;
			if(quiet!=null){quiet.close();}
		}
	}

	/** Checks the final destination before rename/copy; refusal retains the selected assembly.
	 * Uses the native overwrite setting and protects both table-loading and read-processing
	 * input lists. This is preflight validation, not protection against concurrent path changes.
	 * @return True when the selected assembly already occupies the final destination
	 */
	private static boolean checkPublication(final String output, final File selected, final ArrayList<String> args){
		final File target=new File(output);
		final ArrayList<String> inputs=publicationInputs(args);
		for(String input:inputs){
			if(!FileFormat.isStdio(input) && sameFile(new File(input), target)){
				throw new IllegalArgumentException("Final assembly output aliases an input: "+input+" -> "+output+
						"; selected assembly retained at "+selected);
			}
		}
		if(sameFile(selected, target)){return true;}
		if(target.exists() && (!Tadpole.overwrite || !target.isFile())){
			throw new IllegalArgumentException("Final assembly output exists and cannot be replaced: "+output+
					"; overwrite="+Tadpole.overwrite+"; selected assembly retained at "+selected);
		}
		return false;
	}

	/** Resolves the final input options without replaying parsers that mutate global settings.
	 * Tadpole's ine aliases affect read processing but not the table loader; protect both.
	 * Primary/extra lists expand an absent final '#' to mates, while explicit mate lists do not.
	 */
	private static ArrayList<String> publicationInputs(final ArrayList<String> args){
		String primary=null, mate=null, extra=null, reads=null, readMate=null;
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase(Locale.ROOT), value=split.length>1 ? split[1] : null;
			if(key.equals("in") || key.equals("in1")){primary=reads=value;
			}else if(key.equals("in2")){mate=readMate=value;
			}else if(key.equals("ine") || key.equals("ine1")){reads=value;
			}else if(key.equals("ine2")){readMate=value;
			}else if(key.equals("extra")){extra=value;}
		}
		final ArrayList<String> result=new ArrayList<String>();
		addPublicationInputs(result, primary, true, false);
		addPublicationInputs(result, mate, false, false);
		addPublicationInputs(result, extra, true, false);
		addPublicationInputs(result, reads, true, true);
		addPublicationInputs(result, readMate, false, true);
		return result;
	}

	/** Appends one resolved list, matching loader expansion and optional Tadpole extension repair. */
	private static void addPublicationInputs(final ArrayList<String> inputs, final String value,
			final boolean expandPairs, final boolean fixExtension){
		if(value==null){return;}
		for(String input:value.split(",")){
			final int pound=input.lastIndexOf('#');
			if(expandPairs && pound>=0 && !new File(input).exists()){
				final String prefix=input.substring(0, pound), suffix=input.substring(pound+1);
				inputs.add(fixExtension ? Tools.fixExtension(prefix+1+suffix) : prefix+1+suffix);
				inputs.add(fixExtension ? Tools.fixExtension(prefix+2+suffix) : prefix+2+suffix);
			}else{inputs.add(fixExtension ? Tools.fixExtension(input) : input);}
		}
	}

	/** Compares canonical paths and existing hard links; inability to check is a loud failure. */
	private static boolean sameFile(final File a, final File b){
		try{
			return a.getCanonicalFile().equals(b.getCanonicalFile()) ||
					(a.exists() && b.exists() && Files.isSameFile(a.toPath(), b.toPath()));
		}catch(IOException e){throw new RuntimeException("Could not verify assembly paths: "+a+" and "+b, e);}
	}

	/**
	 * Tests a smaller supported K derived from 70% of the leftmost trial's K.
	 * Inserts the new trial in K order and follows strict improvements recursively.
	 * The list and the final two argument slots are mutated; each trial runs Tadpole.
	 *
	 * @param list List of assembly records for different k-mer values
	 * @param best Index of current best assembly in the list
	 * @param argList Command-line arguments template for Tadpole execution
	 * @param contigsName Output filename template with % placeholder for k-mer value
	 * @param recur Retained compatibility parameter; expansion does not consult it
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

			Record r=readStatistics(k1, contigsName.replace("%", ""+k1));
			list.add(best, r);
			best++;

			if(r.compareTo(mid)>0){
				return expandLeft(list, best-1, argList, contigsName, true);
			}
		}
		return best;
	}

	/**
	 * Tests a larger supported K derived from min(current+40, current*1.25+10).
	 * Appends the trial and follows strict improvements recursively. The list and the
	 * final two argument slots are mutated; each trial runs Tadpole.
	 *
	 * @param list List of assembly records for different k-mer values
	 * @param best Index of current best assembly in the list
	 * @param argList Command-line arguments template for Tadpole execution
	 * @param contigsName Output filename template with % placeholder for k-mer value
	 * @param recur Retained compatibility parameter; expansion does not consult it
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

			Record r=readStatistics(k1, contigsName.replace("%", ""+k1));
			list.add(r);

			if(r.compareTo(mid)>0){
				return expandRight(list, best+1, argList, contigsName, true);
			}
		}
		return best;
	}

	/**
	 * Optionally expands boundary winners, then tests midpoints around an interior winner.
	 * The neighboring trials are adjacent in K order, not necessarily the second-best
	 * assemblies. Rounded midpoints must differ from both interval ends. Improvements
	 * recurse; after neither midpoint improves, recur allows one additional refinement.
	 * Mutates the ordered trial list and the final two command-line argument slots.
	 *
	 * @param list List of assembly records for different k-mer values
	 * @param best Index of current best assembly in the list
	 * @param argList Command-line arguments template for Tadpole execution
	 * @param contigsName Output filename template with % placeholder for k-mer value
	 * @param recur Allows one additional refinement after neither midpoint improves; improvements always recurse
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

			Record r=readStatistics(k1, contigsName.replace("%", ""+k1));
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

			Record r=readStatistics(k2, contigsName.replace("%", ""+k2));
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

	/*--------------------------------------------------------------*/
	/*----------------          Trial Record        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Snapshot of reported assembly statistics for one supported K.
	 * Comparison prefers having contigs, then orders fixed 64-base length bins,
	 * count and K. It ranks recorded metrics without claiming biological quality.
	 */
	private static class Record implements Comparable<Record>{

		/** Creates an empty snapshot with Java's zero defaults. */
		public Record(){}

		/** Copies scalar metrics; L50/L90 are scaffold length thresholds, not counts. */
		public Record(int k_, long L50_, long L90_, long contigLen_, long contigs_, long maxContig_){
			k=k_;
			L50=L50_;
			L90=L90_;
			contigLen=contigLen_;
			contigs=contigs_;
			maxContig=maxContig_;
		}

		/**
		 * Prefers an assembly with contigs over one without, including in the zero-length bin.
		 * Then compares L50, L90 and maximum length in order, rounding each down to a
		 * multiple of 64 bases. Larger bins win, followed by fewer contigs and smaller K.
		 * Fixed per-record bins preserve transitivity; pairwise percentage tolerances
		 * did not (ASM012). Raw lengths within a bin and contigLen are not tie-breakers.
		 * Records must not be modified during a comparison sequence.
		 *
		 * @param b The Record to compare against
		 * @return Negative for worse, positive for better (including null b), or zero for equal comparison keys
		 */
		@Override
		public int compareTo(Record b){
			if(b==null){return 1;}
			int x=Boolean.compare(contigs>0, b.contigs>0);
			if(x!=0){return x;}
			x=Long.compare(L50&~63L, b.L50&~63L);
			if(x!=0){return x;}
			x=Long.compare(L90&~63L, b.L90&~63L);
			if(x!=0){return x;}
			x=Long.compare(maxContig&~63L, b.maxContig&~63L);
			if(x!=0){return x;}
			x=Long.compare(b.contigs, contigs);
			if(x!=0){return x;}
			return Integer.compare(b.k, k);
		}

		/** Supported K used by this trial. */
		int k;
		/** Scaffold length at the 50% base threshold, copied from AssemblyStats2.lastL50. */
		long L50;
		/** Scaffold length at the 90% base threshold, copied from AssemblyStats2.lastL90. */
		long L90;
		/** Total contig bases; retained in the snapshot but unused by compareTo. */
		long contigLen;
		/** Number of contigs in the reported assembly. */
		long contigs;
		/** Longest reported contig, in bases. */
		long maxContig;

	}

	/*--------------------------------------------------------------*/
	/*----------------        Shared State          ----------------*/
	/*--------------------------------------------------------------*/

	/** Initial supported K values in ascending order; rebuilt by process. */
	private static int[] kmers;
	/** Stops the initial ascending scan on a strictly worse trial, not a tie. */
	private static boolean quitEarly=false;
	/** Requests best-effort deletion of nonselected intermediates, excluding the final destination. */
	private static boolean delete=false;
	/** Enables recursive exploration beyond the tested K boundaries. */
	private static boolean expand=false;
	/** Enables refinement between the winner and adjacent tested K values. */
	private static boolean bisect=false;
	/** Report visibility for this invocation; metrics are always computed for selection. */
	private static boolean showStats=true;
	/** Selected intermediate pathname; publication may move or remove this path. */
	private static String bestAssembly=null;

	/** Configured diagnostics; several legacy messages still write directly to System.err. */
	private static PrintStream outstream=System.err;

}
