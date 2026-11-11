package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import assemble.AbstractExploreThread;
import assemble.AbstractShaveThread;
import assemble.Shaver;
import kmer.AbstractKmerTableSet;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Injected worker-failure contracts for Shaver's ordinary explore/remove phases.
 * Uses a tiny real count table but throws only from test workers.
 * @author Brian Bushnell, Jean
 */
public final class ShaverFailureTest{

	/** Checks normal completion and propagation of ordinary worker failures. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ShaverFailureTest.class, false);
		try{
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else{throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(root==null){throw new IllegalArgumentException("Specify new outdir=DIR.");}
			Files.createDirectory(root);
			final Path input=root.resolve("reads.fa");
			Files.write(input, FASTA.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			final String[] options={"in="+input, "k=63", "packed=t", "rcomp=t", "ways=1",
					"initialsize=128", "prealloc=f", "prefilter=f", "minprob=0", "minprobmain=f",
					"interleaved=f", "t=1"};
			final KmerTableSetU table=new KmerTableSetU(options, 0);
			table.process(new Timer());
			verifyLoaded(table);
			try{
				int failed=0;
				final StringBuilder checks=new StringBuilder("case\tpassed\tdetail\n");
				failed+=control(checks, table);
				failed+=injected(checks, table, "explore-runtime", new RuntimeException("explore-runtime"), null);
				failed+=injected(checks, table, "explore-assertion", new AssertionError("explore-assertion"), null);
				failed+=injected(checks, table, "remove-runtime", null, new RuntimeException("remove-runtime"));
				failed+=injected(checks, table, "remove-assertion", null, new AssertionError("remove-assertion"));
				Files.write(root.resolve("checks.tsv"), checks.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
				System.out.println("SHAVER_FAILURE_TEST checks=5 failed="+failed);
				if(failed>0){throw new AssertionError("Shaver worker-failure contract failures; inspect "+root);}
			}finally{table.clear();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No-op workers should finish both phases normally and report no removals. */
	private static int control(final StringBuilder checks, final AbstractKmerTableSet table){
		table.clearOwnership();
		final InjectedShaver shaver=new InjectedShaver(table, null, null);
		try{
			final long removed=shaver.shave();
			return check(checks, "control", removed==0 && shaver.removeFactoryCalls==1, "removed="+removed+", removeFactories="+shaver.removeFactoryCalls);
		}catch(Throwable t){return check(checks, "control", false, "unexpected "+t.getClass().getSimpleName());}
	}

	/** A fixed owner should propagate each worker failure and skip removal after explorer failure. */
	private static int injected(final StringBuilder checks, final AbstractKmerTableSet table,
			final String name, final Throwable explore, final Throwable remove){
		table.clearOwnership();
		final InjectedShaver shaver=new InjectedShaver(table, explore, remove);
		try{
			shaver.shave();
			if(!shaver.injectionRan()){
				throw new AssertionError("Injected worker hook did not run for "+name);
			}
			return check(checks, name, false, "shave returned; removeFactories="+shaver.removeFactoryCalls);
		}catch(Throwable t){
			final Throwable expected=explore==null ? remove : explore;
			final Throwable cause=t.getCause()==null ? t : t.getCause();
			if(cause!=expected){throw new AssertionError("Unexpected exception before injected failure for "+name, t);}
			final boolean propagated=cause==expected;
			final boolean skipped=explore==null || shaver.removeFactoryCalls==0;
			return check(checks, name, propagated && skipped,
					"cause="+cause.getClass().getSimpleName()+", removeFactories="+shaver.removeFactoryCalls);
		}
	}

	/** Keeps every scientific expectation active without depending on JVM assertions. */
	private static int check(final StringBuilder checks, final String name,
			final boolean passed, final String detail){
		checks.append(name).append('\t').append(passed ? 1 : 0).append('\t').append(detail).append('\n');
		return passed ? 0 : 1;
	}

	/** Verifies that the tiny real graph loaded at least one expected kmer. */
	private static void verifyLoaded(final KmerTableSetU table){
		final Kmer kmer=new Kmer(table.kbig(), table.hollow);
		for(int i=0; i<table.kbig(); i++){kmer.addRight(FASTA_BASES.charAt(i));}
		if(table.getCount(kmer)<=0){throw new IllegalStateException("Fixture table did not load the first expected kmer.");}
	}

	/** Test-only Shaver whose workers either finish immediately or fail on their first hook. */
	private static final class InjectedShaver extends Shaver{
		InjectedShaver(final AbstractKmerTableSet tables_, final Throwable exploreFailure_,
				final Throwable removeFailure_){
			super(tables_, 1, 1, 1, 1, 1, 3, 100, 100, true, true);
			tables=tables_;
			exploreFailure=exploreFailure_;
			removeFailure=removeFailure_;
		}

		@Override
		protected AbstractExploreThread makeExploreThread(final int id){
			return new AbstractExploreThread(id, tables.kbig()){
				@Override
				protected boolean processNextTable(final Kmer kmer, final Kmer temp){
					exploreHookRan=true;
					throwIf(exploreFailure);
					return false;
				}
				@Override
				protected boolean processNextVictims(final Kmer kmer, final Kmer temp){return false;}
			};
		}

		@Override
		protected AbstractShaveThread makeShaveThread(final int id){
			removeFactoryCalls++;
			return new AbstractShaveThread(id){
				@Override
				protected boolean processNextTable(){
					removeHookRan=true;
					throwIf(removeFailure);
					return false;
				}
			};
		}

		@Override
		protected AbstractKmerTableSet tables(){return tables;}

		boolean injectionRan(){return exploreFailure==null ? removeHookRan : exploreHookRan;}

		private final AbstractKmerTableSet tables;
		private final Throwable exploreFailure, removeFailure;
		private boolean exploreHookRan, removeHookRan;
		private int removeFactoryCalls;
	}

	/** Throws RuntimeException and Error injections without wrapping them in the worker. */
	private static void throwIf(final Throwable failure){
		if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
		if(failure instanceof Error){throw (Error)failure;}
		if(failure!=null){throw new RuntimeException(failure);}
	}

	private static final String FASTA_BASES="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
			+"ACGTACGTACGTACGTACGT";
	private static final String FASTA=">a\n"+FASTA_BASES+"\n";
}
