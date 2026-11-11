package test.assemble;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Random;

import assemble.ShaveObject;
import assemble.Contig;
import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import kmer.KmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import stream.ConcurrentReadInputStream;
import structures.ByteBuilder;
import structures.IntList;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Finite status/sequence contracts for a real counted path reaching its length cap.
 * Each invocation uses a fresh JVM, one concrete representation and new output directory.
 * Includes an actual dead end and a competing owner at the cap as unaffected controls.
 * Optional caller=table|read also tests the real two-direction contig builders.
 * This does not test complete assembly or graph postprocessing.
 * @author Brian Bushnell, Shinobu
 */
public final class TadpoleLengthLimitTest{

	/** Accepts k, ownership, case and outdir; caller=direct|table|read defaults to direct. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TadpoleLengthLimitTest.class, false);
		try{
			int k=0;
			boolean owned=false, ownershipSet=false;
			String name=null, caller="direct";
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(pair[1]);}
				else if(pair[0].equalsIgnoreCase("ownership")){owned=Parse.parseBoolean(pair[1]); ownershipSet=true;}
				else if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else if(pair[0].equalsIgnoreCase("caller")){caller=pair[1];}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			check(root!=null && ownershipSet && (k==31 || k==63) &&
					("cap".equals(name) || "dead-end".equals(name) || "owner-at-cap".equals(name)),
					"Specify k=31|63 ownership=t|f case=cap|dead-end|owner-at-cap outdir=NEW.");
			check(!"owner-at-cap".equals(name) || owned, "The competing-owner control requires ownership=t.");
			check("direct".equals(caller) || "table".equals(caller) || "read".equals(caller), "Unknown caller: "+caller);
			check("direct".equals(caller) || !"owner-at-cap".equals(name), "Caller cases require cap or dead-end.");
			Files.createDirectory(root);
			run(root, k, owned, name, caller);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Counts a unique reference, checks both branch directions, and observes one extension. */
	private static void run(final Path root, final int k, final boolean owned, final String name, final String caller) throws Exception{
		final byte[] truth=new byte[2*k+17];
		final Random random=new Random(706119L);
		for(int i=0; i<truth.length; i++){truth[i]=BASES[random.nextInt(4)];}
		final String sequence=new String(truth, StandardCharsets.US_ASCII);
		final StringBuilder reference=new StringBuilder();
		for(int copy=0; copy<DEPTH; copy++){reference.append('>').append(copy).append('\n').append(sequence).append('\n');}
		final Path input=root.resolve("counts.fa");
		write(input, reference.toString());
		final boolean deadEnd=name.equals("dead-end"), conflict=name.equals("owner-at-cap");
		final int limit=deadEnd ? truth.length+8 : k+4;
		String[] options={"in="+input, "out=null", "mode=contig", "k="+k, "packed=f",
				"ownership="+owned, "maxlen="+limit, "extendright=-1", "t=1", "buildthreads=1",
				"ways=1", "initialsize=1024", "prealloc=f", "prefilter=0", "minprob=0",
				"minprobmain=f", "rcomp=t", "maskmiddle=f", "maskcore=f", "qin=33",
				"mcs=2", "mce=2", "pop=f", "shave=f", "rinse=f", "interleaved=f",
				"pigz=f", "unpigz=f", "showstats=f", "ecc=f", "ecco=f", "merge=f"};
		if(!caller.equals("direct")){
			final int n=options.length;
			options=Arrays.copyOf(options, n+4);
			options[n]="minextension=0"; options[n+1]="mincontig=0";
			options[n+2]="trimends=0"; options[n+3]="trimcircular=f";
		}
		Files.write(root.resolve("table.args"), Arrays.asList(options), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final Tadpole tad=Tadpole.makeTadpole(options, false);
		try{
			check(tad.tables().kbig()==k && tad.maxContigLen==limit && !Kmer.PACKED && !Kmer.MASK_CORE &&
					(k==31 ? tad instanceof Tadpole1 : tad instanceof Tadpole2), "Unexpected fixture geometry or cap.");
			tad.loadKmers(new Timer());
			verifyPath(tad, truth, k);
			if(owned){tad.tables().initializeOwnership();}
			if(!caller.equals("direct")){
				runCaller(root, tad, truth, k, owned, name, caller, limit);
				return;
			}
			if(conflict){
				final byte[] window=Arrays.copyOfRange(truth, limit-k, limit);
				final boolean claimed=k==31 ? ((KmerTableSet)tad.tables()).claim(window, k, 1, true) :
						((KmerTableSetU)tad.tables()).claim(window, k, 1, true, new Kmer(k));
				check(claimed, "The cap-ending word must belong to competing owner1 before extending as owner0.");
			}
			final ByteBuilder bases=new ByteBuilder(truth.length);
			bases.append(truth, 0, k);
			int status=Integer.MIN_VALUE;
			AssertionError failure=null;
			try{
				status=k==31 ? ((Tadpole1)tad).extendToRight(bases, new int[4], new int[4], 0) :
						((Tadpole2)tad).extendToRight(bases, new int[4], new int[4], 0, new Kmer(k));
			}catch(AssertionError e){failure=e;}
			if(failure!=null){
				final StringWriter text=new StringWriter();
				final PrintWriter writer=new PrintWriter(text);
				failure.printStackTrace(writer); writer.flush();
				write(root.resolve("extension-error.txt"), text.toString());
			}
			final int expected=deadEnd ? ShaveObject.DEAD_END : conflict ? ShaveObject.BAD_OWNER : ShaveObject.TOO_LONG;
			final byte[] wanted=Arrays.copyOf(truth, deadEnd ? truth.length : limit);
			final byte[] observed=bases.toBytes();
			final boolean sequenceMatches=Arrays.equals(wanted, observed);
			final boolean passed=failure==null && status==expected && sequenceMatches;
			write(root.resolve("observed.fa"), ">observed\n"+new String(observed, StandardCharsets.US_ASCII)+"\n");
			final String report="case\tk\townership\texpected_status\tactual_status\tassertion\tsequence\tpassed\n"+
					name+"\t"+k+"\t"+owned+"\t"+expected+"\t"+status+"\t"+(failure!=null)+"\t"+sequenceMatches+"\t"+passed+"\n";
			write(root.resolve("checks.tsv"), report);
			System.out.print(report);
			System.out.println("LENGTH_LIMIT assertions="+TadpoleLengthLimitTest.class.desiredAssertionStatus()+" setup=true passed="+passed);
			if(!passed){throw new AssertionError("Extension length-limit status or sequence contract failed.");}
		}finally{tad.clearData();}
	}

	/** Invokes the actual unstarted worker's private builder; no substituted extension method. */
	private static void runCaller(final Path root, final Tadpole tad, final byte[] truth, final int k,
			final boolean owned, final String name, final String caller, final int limit) throws Exception{
		final boolean table=caller.equals("table"), deadEnd=name.equals("dead-end");
		final byte[] seed=Arrays.copyOf(truth, k), savedSeed=seed.clone();
		final Method factory=tad.getClass().getDeclaredMethod("makeBuildThread", int.class, int.class, ConcurrentReadInputStream[].class);
		factory.setAccessible(true);
		final Object worker=factory.invoke(tad, 0, ShaveObject.contigMode, null);
		check(worker instanceof Thread && ((Thread)worker).getState()==Thread.State.NEW, "Caller fixture must not start a worker.");
		final Field builderField=worker.getClass().getSuperclass().getDeclaredField("builderT");
		builderField.setAccessible(true);
		final ByteBuilder builder=(ByteBuilder)builderField.get(worker);
		final Method method;
		final Object[] arguments;
		if(table && k==31){
			long forward=0, reverse=0;
			for(int i=0; i<k; i++){
				forward=(forward<<2)|baseNumber(seed[i]);
				reverse=(reverse<<2)|(3-baseNumber(seed[k-1-i]));
			}
			method=worker.getClass().getDeclaredMethod("makeContig", long.class, ByteBuilder.class, boolean.class);
			arguments=new Object[]{Math.max(forward, reverse), builder, false};
		}else if(table){
			final Kmer word=new Kmer(k);
			for(byte b : seed){word.addRight(b);}
			if(!Arrays.equals(word.array1(), word.key())){word.rcomp();}
			method=worker.getClass().getDeclaredMethod("makeContig", ByteBuilder.class, Kmer.class, boolean.class);
			arguments=new Object[]{builder, word, false};
		}else if(k==31){
			method=worker.getClass().getDeclaredMethod("makeContig", byte[].class, ByteBuilder.class, long.class);
			arguments=new Object[]{seed, builder, 0L};
		}else{
			method=worker.getClass().getDeclaredMethod("makeContig", byte[].class, ByteBuilder.class, long.class, Kmer.class);
			arguments=new Object[]{seed, builder, 0L, new Kmer(k)};
		}
		method.setAccessible(true);
		Object result=null;
		Throwable failure=null;
		try{result=method.invoke(worker, arguments);}
		catch(InvocationTargetException e){failure=e.getCause();}
		if(failure!=null){
			final StringWriter text=new StringWriter();
			final PrintWriter writer=new PrintWriter(text);
			failure.printStackTrace(writer); writer.flush();
			write(root.resolve("caller-error.txt"), text.toString());
		}
		byte[] wanted=Arrays.copyOf(truth, deadEnd ? truth.length : limit);
		int expectedLeft=ShaveObject.DEAD_END, expectedRight=deadEnd ? ShaveObject.DEAD_END : ShaveObject.TOO_LONG;
		if(table){
			final byte[] reverse=reverseComplement(wanted);
			if(compare(wanted, reverse)>0){
				wanted=reverse;
				final int tmp=expectedLeft; expectedLeft=expectedRight; expectedRight=tmp;
			}
		}
		final byte[] observed=result==null ? new byte[0] : table ? ((Contig)result).bases : (byte[])result;
		int left=-1, right=-1;
		if(table && result!=null){
			final Field l=Contig.class.getDeclaredField("leftCode"), r=Contig.class.getDeclaredField("rightCode");
			l.setAccessible(true); r.setAccessible(true);
			left=l.getInt(result); right=r.getInt(result);
		}
		final boolean sequence=Arrays.equals(wanted, observed), seedUnchanged=Arrays.equals(seed, savedSeed);
		final boolean endpoints=!table || (result!=null && left==expectedLeft && right==expectedRight);
		final boolean passed=failure==null && result!=null && sequence && seedUnchanged && endpoints;
		write(root.resolve("caller-expected.fa"), ">expected\n"+new String(wanted, StandardCharsets.US_ASCII)+"\n");
		write(root.resolve("caller-observed.fa"), ">observed\n"+new String(observed, StandardCharsets.US_ASCII)+"\n");
		write(root.resolve("caller-builder.fa"), ">builder\n"+builder.toString()+"\n");
		final String report="caller\tcase\tk\townership\treturned\texception\tsequence\tseed\tleft\tright\tendpoints\tpassed\n"+
				caller+"\t"+name+"\t"+k+"\t"+owned+"\t"+(result!=null)+"\t"+
				(failure==null ? "none" : failure.getClass().getName())+"\t"+sequence+"\t"+seedUnchanged+"\t"+
				left+"\t"+right+"\t"+endpoints+"\t"+passed+"\n";
		write(root.resolve("caller-checks.tsv"), report);
		System.out.print(report);
		System.out.println("LENGTH_CALLER assertions="+TadpoleLengthLimitTest.class.desiredAssertionStatus()+" setup=true passed="+passed);
		if(!passed){throw new AssertionError("Contig builder length-limit contract failed.");}
	}

	/** Independent literal DNA encoding for constructing the short seed and expected strand. */
	private static int baseNumber(final byte b){
		switch(b){case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3;
			default: throw new IllegalArgumentException("Non-DNA fixture base: "+b);}
	}

	/** Computes the fixture oracle without calling the contig canonicalization implementation. */
	private static byte[] reverseComplement(final byte[] bases){
		final byte[] answer=new byte[bases.length];
		for(int i=0; i<bases.length; i++){answer[i]=BASES[3-baseNumber(bases[bases.length-1-i])];}
		return answer;
	}

	/** Lexicographic comparison of equal-length fixture DNA sequences. */
	private static int compare(final byte[] a, final byte[] b){
		check(a.length==b.length, "Strand comparison requires equal lengths.");
		for(int i=0; i<a.length; i++){if(a[i]!=b[i]){return a[i]-b[i];}}
		return 0;
	}

	/** Rejects repeated words and every competing one-base predecessor or successor. */
	private static void verifyPath(final Tadpole tad, final byte[] truth, final int k){
		final IntList counts=new IntList(truth.length);
		final Kmer scratch=new Kmer(k);
		tad.tables().fillCounts(truth, counts, scratch);
		check(counts.size==truth.length-k+1, "Reference lacks one count per complete window.");
		for(int i=0; i<counts.size; i++){check(counts.get(i)==DEPTH, "Repeated or missing reference word at "+i);}
		for(int start=0; start+k<=truth.length; start++){
			final byte[] window=Arrays.copyOfRange(truth, start, start+k);
			for(int end : new int[]{0, k-1}){
				final byte original=window[end];
				for(byte replacement : BASES){if(replacement!=original){
					window[end]=replacement;
					tad.tables().fillCounts(window, counts, scratch);
					check(counts.size==1 && counts.get(0)==0, "Competing predecessor/successor invalidates fixture at "+start);
				}}
				window[end]=original;
			}
		}
	}

	/** Saves a fresh tiny UTF-8 artifact without replacing earlier evidence. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Setup checks remain active in the assertion-disabled comparison. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Independently chosen positive support depth, above both explicit thresholds. */
	private static final int DEPTH=8;
	/** Seeded DNA alphabet in fixed random-draw order. */
	private static final byte[] BASES={'A', 'C', 'G', 'T'};
}
