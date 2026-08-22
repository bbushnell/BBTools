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
 * Thin orchestrator for BBCMS two-pass mode (Item 3a).  With passes=2, runs a
 * cheap K=31 junk-removal pass first (survivors written to a temp file), then
 * runs the user's normal correction pass on the survivors.  With passes=1
 * (default), delegates straight to BloomFilterCorrectorWrapper with zero
 * overhead.  Never modifies BloomFilterCorrectorWrapper; each pass is a
 * complete, independent invocation of it, each building its own filter.
 *
 * Design: Amber.  Implementation: Nepgear.
 *
 * @author Brian Bushnell, Nepgear
 * @date August 21, 2026
 */
public class BloomFilterCorrectorMultipass {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		BloomFilterCorrectorMultipass x=new BloomFilterCorrectorMultipass(args);
		x.process(t);
	}

	/**
	 * Constructor.  Parses the multipass-specific flags (passes=, k0=, bits0=,
	 * hcf0=, mincount0=, tossjunk0=, ecc0=) and partitions everything else into
	 * in1/in2 (recorded separately, since both passes need to control them) and
	 * pass2Args (forwarded to the correction pass verbatim).
	 * @param args Command line arguments
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
	 * Runs either a single delegated pass or the full two-pass pipeline.
	 * @param t Timer started at program entrance
	 */
	void process(Timer t){
		if(passes<=1){
			//pass2Args deliberately excludes in=/in2= (needed separately for the
			//passes=2 swap) -- add them back for a true zero-overhead delegate.
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
		//cleaner"). Verified live: twin in=/in2= with only out= (no out2=) set
		//makes BloomFilterCorrectorWrapper write interleaved output, and a
		//plain in= with no in2= on the reread correctly auto-detects as paired.
		final String temp1;
		try{
			temp1=File.createTempFile("bbcms_p1_", ".fq.gz", tmpDir).getAbsolutePath();
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		try{
			System.err.println("\n----- BBCMS Pass 1: junk removal (k="+k0+") -----");
			BloomFilterCorrectorWrapper.main(buildPass1Args(temp1));

			//Pass 1's filter (can be gigabytes) is unreachable now but not yet
			//collected -- Runtime.freeMemory() still reports it as used, so pass
			//2's memFraction-based cell sizing (BloomFilter -> KCountArray7MTA)
			//undersizes and can crash on a negative cell count (verified live:
			//AssertionError in Primes.primeAtMost at -Xmx4.5g). Force the
			//collection before pass 2 reads free memory.
			System.gc();

			System.err.println("\n----- BBCMS Pass 2: correction -----");
			BloomFilterCorrectorWrapper.main(buildPass2Args(temp1));
		}finally{
			delete(temp1);
		}
		t.stop();
	}

	/**
	 * Builds the arg array for the junk-removal pass: original input, a
	 * hardcoded (k0-overridable) filter/threshold set, a single interleaved
	 * temp output, and a small whitelist of user args safe to carry over
	 * (see FORWARD_TO_PASS1).
	 * @param temp1 Pass-1 interleaved output path
	 * @return Argument array for BloomFilterCorrectorWrapper.main
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
	 * Builds the arg array for the correction pass: the user's original args
	 * verbatim, with in=/in2= swapped to the single interleaved pass-1 output
	 * and maxreads= dropped (pass 1 already bounded the input).
	 * @param temp1 Pass-1 interleaved output path to read as pass-2 input
	 * @return Argument array for BloomFilterCorrectorWrapper.main
	 */
	private String[] buildPass2Args(String temp1){
		ArrayList<String> list=new ArrayList<String>();
		list.add("in="+temp1);
		boolean sawOw=false;
		for(String arg : pass2Args){
			//pass2Args never contains in2=/input2= -- captured separately in the
			//constructor (merged into the single interleaved temp1 above).
			String a=arg.split("=")[0].toLowerCase();
			if(a.equals("reads") || a.equals("maxreads")){continue;}
			if(a.equals("ow") || a.equals("overwrite")){sawOw=true;}
			list.add(arg);
		}
		if(!sawOw){list.add("ow=t");}
		return list.toArray(new String[0]);
	}

	/**
	 * Deletes a temp file if it exists.
	 * @param s Path to delete, or null (no-op)
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

	/** Flags safe to carry from the user's args into the pass-1 (junk-removal) invocation. */
	private static final HashSet<String> FORWARD_TO_PASS1=new HashSet<String>(Arrays.asList(
			"t", "threads", "rcomp", "memmult", "memfraction", "memratio", "cells", "seed",
			"ordered", "verbose", "interleaved", "int", "qfin", "qfin1", "qfin2"));

	/** Number of passes; only 1 (delegate) and 2 (junk-removal + correction) are supported. */
	final int passes;
	/** Pass-1 kmer length override (default 31). */
	final int k0;
	/** Pass-1 bits-per-cell override (default 2). */
	final int bits0;
	/** Pass-1 high-count-fraction override (default 0.4). */
	final float hcf0;
	/** Pass-1 mincount override (default 2). */
	final int mincount0;
	/** Pass-1 tossjunk override (default true). */
	final boolean tossjunk0;
	/** Pass-1 ecc override (default false). */
	final boolean ecc0;
	/** User's original primary input file. */
	final String in1;
	/** User's original secondary input file, or null if unpaired/interleaved. */
	final String in2;
	/** The user's original reads=/maxreads= arg, verbatim, or null if unset. */
	final String maxReadsArg;
	/** Every user arg except in=/in2=/maxreads=/the passes=&*0= multipass flags; forwarded to pass 2 as-is. */
	final ArrayList<String> pass2Args;

}
