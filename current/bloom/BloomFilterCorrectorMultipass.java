package bloom;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import parse.Parse;
import shared.Shared;
import shared.Timer;

/**
 * Runs BBCMS once, or filters reads before a second correction invocation.
 * The default {@code passes=1} delegates to {@link BloomFilterCorrectorWrapper}.
 * With {@code passes=2}, first-pass settings default to k=31, two-bit counts,
 * junk removal enabled and correction disabled. Surviving reads go through one
 * temporary compressed FASTQ before the second invocation builds its own filter.
 * Paired survivors share that temporary file as interleaved records.
 *
 * Both wrapper invocations run sequentially in the same JVM; static parser and
 * I/O settings can carry over. This class does not isolate them or support
 * concurrent pipeline invocations. It copies argument references into its own
 * list and leaves the caller's argument array unchanged. The first pass receives
 * only its explicit settings and the flags in {@link #FORWARD_TO_PASS1}; the
 * second receives the remaining user settings with its input and read limit
 * adjusted. The first input stream's actual pairing state determines whether
 * the second invocation reads the temporary FASTQ as interleaved.
 *
 * Design: Amber.  Implementation: Nepgear.
 *
 * @author Brian Bushnell, Nepgear
 * @date August 21, 2026
 */
public class BloomFilterCorrectorMultipass{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses the pipeline settings and runs the requested invocation sequence.
	 * Exceptions thrown by either wrapper invocation propagate to the launcher.
	 * @param args Non-null array of BBCMS flag=value arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		BloomFilterCorrectorMultipass x=new BloomFilterCorrectorMultipass(args);
		x.process(t);
	}

	/**
	 * Separates pipeline flags, input paths and arguments for the wrapper.
	 * Recognizes passes, k0, bits0, hcf0, mincount0, tossjunk0 and ecc0 here;
	 * other flags remain in their original order for downstream parsing. The last
	 * reads/maxreads argument is also saved for the first pass. No files are opened.
	 * Supported pass counts are 1 and 2, checked by an assertion.
	 * @param args Non-null array of flag=value arguments; entries must be non-null
	 */
	public BloomFilterCorrectorMultipass(String[] args){
		int passes_=1;
		int k0_=31;
		int bits0_=2;
		float hcf0_=0.4f;
		int mincount0_=2;
		boolean tossjunk0_=true;
		boolean ecc0_=false;

		String in1_=null, in2_=null;
		String maxReadsArg_=null;
		ArrayList<String> rest=new ArrayList<String>();

		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("passes")){
				passes_=Integer.parseInt(b);
			}else if(a.equals("k0")){
				k0_=Integer.parseInt(b);
			}else if(a.equals("bits0")){
				bits0_=Integer.parseInt(b);
			}else if(a.equals("hcf0")){
				hcf0_=Float.parseFloat(b);
			}else if(a.equals("mincount0")){
				mincount0_=Integer.parseInt(b);
			}else if(a.equals("tossjunk0")){
				tossjunk0_=Parse.parseBoolean(b);
			}else if(a.equals("ecc0")){
				ecc0_=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1_=b;
			}else if(a.equals("in2") || a.equals("input2")){
				in2_=b;
			}else if(a.equals("reads") || a.equals("maxreads")){
				//Kept in pass2Args too so the passes=1 delegate path is unaffected;
				//buildPass2Args strips it back out for the real 2-pass case (the
				//design's edge case: forward to pass 1 only, input is already
				//bounded by the time pass 2 sees it).
				maxReadsArg_=arg;
				rest.add(arg);
			}else{
				rest.add(arg);
			}
		}

		passes=passes_;
		k0=k0_;
		bits0=bits0_;
		hcf0=hcf0_;
		mincount0=mincount0_;
		tossjunk0=tossjunk0_;
		ecc0=ecc0_;
		in1=in1_;
		in2=in2_;
		maxReadsArg=maxReadsArg_;
		pass2Args=rest;

		assert(passes==1 || passes==2) : "passes= must be 1 or 2, not "+passes;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs the wrapper once or runs both passes through a temporary FASTQ.
	 * Two-pass mode creates the temporary file under Shared.tmpdir(), when set,
	 * otherwise the JVM temporary directory. After file creation, cleanup is
	 * attempted even when a wrapper throws; deletion failure is not reported.
	 * A first-pass failure prevents the second invocation. The timer is stopped
	 * only after normal completion, and its elapsed time is not printed here.
	 * @param t Non-null timer started at program entrance
	 */
	void process(Timer t){
		if(passes<=1){
			//pass2Args deliberately excludes in=/in2= (needed separately for the
			//passes=2 swap) -- add them back before the single wrapper invocation.
			ArrayList<String> list=new ArrayList<String>(pass2Args.size()+2);
			if(in1!=null){list.add("in="+in1);}
			if(in2!=null){list.add("in2="+in2);}
			list.addAll(pass2Args);
			BloomFilterCorrectorWrapper.main(list.toArray(new String[0]));
			t.stop();
			return;
		}

		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Shared.tmpdir() is null unless TMPDIR is set (same convention as
		//TadPipe.java): a null File arg to createTempFile means "use the JVM
		//default (java.io.tmpdir)".
		final File tmpDir=(Shared.tmpdir()==null ? null : new File(Shared.tmpdir()));
		if(tmpDir!=null && !tmpDir.exists()){tmpDir.mkdirs();}

		//Always a SINGLE interleaved temp file, regardless of whether the
		//user's own input/output is twin or interleaved (Brian, 2026-08-21:
		//"intermediate files should generally be interleaved... it's just
		//cleaner"). BLM201: carry the input stream's actual pairing state into
		//pass 2; original input flags and header-based guessing are insufficient.
		final String temp1;
		try{
			temp1=File.createTempFile("bbcms_p1_", ".fq.gz", tmpDir).getAbsolutePath();
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		try{
			System.err.println("\n----- BBCMS Pass 1: junk removal (k="+k0+") -----");
			final boolean paired=BloomFilterCorrectorWrapper.run(buildPass1Args(temp1));

			//Pass 1's filter (can be gigabytes) is unreachable now but not yet
			//collected -- Runtime.freeMemory() still reports it as used, so pass
			//2's memFraction-based cell sizing (BloomFilter -> KCountArray7MTA)
			//undersizes and can crash on a negative cell count (verified live:
			//AssertionError in Primes.primeAtMost at -Xmx4.5g). Request collection
			//before pass 2 reads free memory; the JVM can ignore System.gc().
			System.gc();

			System.err.println("\n----- BBCMS Pass 2: correction -----");
			BloomFilterCorrectorWrapper.main(buildPass2Args(temp1, paired));
		}finally{
			delete(temp1);
		}
		t.stop();
	}

	/**
	 * Builds first-pass arguments from the original input and the *0 settings.
	 * Output is the single temporary FASTQ, with overwrite enabled and smooth=3.
	 * Appends the last read-limit argument and whitelisted user flags in input
	 * order. Options absent from this list use the wrapper's defaults or shared
	 * process state; for example, the user's merge flag is not forwarded here.
	 * @param temp1 Non-null path of the already-created temporary output file
	 * @return Fresh argument array for BloomFilterCorrectorWrapper.main
	 */
	private String[] buildPass1Args(String temp1){
		ArrayList<String> list=new ArrayList<String>();
		list.add("in="+in1);
		if(in2!=null){list.add("in2="+in2);}
		list.add("out="+temp1);
		list.add("k="+k0);
		list.add("bits="+bits0);
		list.add("tossjunk="+tossjunk0);
		list.add("ecc="+ecc0);
		list.add("mincount="+mincount0);
		list.add("hcf="+hcf0);
		list.add("smooth=3");
		list.add("ow=t");
		if(maxReadsArg!=null){list.add(maxReadsArg);}
		for(String arg : pass2Args){
			String a=arg.split("=")[0].toLowerCase();
			if(FORWARD_TO_PASS1.contains(a)){list.add(arg);}
		}
		return list.toArray(new String[0]);
	}

	/**
	 * Builds second-pass arguments using the first pass's single output as input.
	 * Drops every reads/maxreads argument because the first pass already limited
	 * input. Replaces original interleaved/int flags with the first input stream's
	 * actual pairing state. Preserves other stored arguments and their order;
	 * adds ow=t only if no ow/overwrite flag was given. Does not modify the stored list.
	 * @param temp1 Non-null path of the first pass's completed temporary FASTQ
	 * @param paired Whether the first input stream processed paired reads
	 * @return Fresh argument array for BloomFilterCorrectorWrapper.main
	 */
	private String[] buildPass2Args(String temp1, boolean paired){
		ArrayList<String> list=new ArrayList<String>();
		list.add("in="+temp1);
		boolean sawOw=false;
		for(String arg : pass2Args){
			//pass2Args never contains in2=/input2= -- captured separately in the
			//constructor (merged into the single interleaved temp1 above).
			String a=arg.split("=")[0].toLowerCase();
			if(a.equals("reads") || a.equals("maxreads")){continue;}
			if(a.equals("interleaved") || a.equals("int")){continue;}
			if(a.equals("ow") || a.equals("overwrite")){sawOw=true;}
			list.add(arg);
		}
		if(!sawOw){list.add("ow=t");}
		list.add("interleaved="+paired);
		return list.toArray(new String[0]);
	}

	/**
	 * Attempts to delete the pipeline's temporary file if it exists.
	 * Ignores File.delete()'s false result; it does not guarantee removal or
	 * suppress an exception thrown by the filesystem security checks.
	 * @param s Temporary path owned by this pipeline, or null for no operation
	 */
	private static void delete(String s){
		if(s!=null){
			File f=new File(s);
			if(f.exists()){f.delete();}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** User flag names forwarded to pass 1; this list does not translate format or pairing settings. */
	private static final HashSet<String> FORWARD_TO_PASS1=new HashSet<String>(Arrays.asList(
			"t", "threads", "rcomp", "memmult", "memfraction", "memratio", "cells", "seed",
			"ordered", "verbose", "interleaved", "int", "qfin", "qfin1", "qfin2"));

	/** Number of passes; only 1 (delegate) and 2 (junk-removal + correction) are supported. */
	final int passes;
	/** First-pass k-mer length in bases; default 31. */
	final int k0;
	/** First-pass count-cell width in bits; default 2. */
	final int bits0;
	/** First-pass high-count-fraction flag; default 0.4. */
	final float hcf0;
	/** First-pass minimum-count flag; default 2. */
	final int mincount0;
	/** Whether the first pass requests junk-read removal; default true. */
	final boolean tossjunk0;
	/** Whether the first pass requests error correction; default false. */
	final boolean ecc0;
	/** User's original primary input file. */
	final String in1;
	/** User's secondary input path, or null when no second path was specified. */
	final String in2;
	/** Last reads/maxreads argument, verbatim, or null; applied to pass 1 only in two-pass mode. */
	final String maxReadsArg;
	/** Non-input, non-multipass args; includes reads/maxreads until buildPass2Args removes them. */
	final ArrayList<String> pass2Args;

}
