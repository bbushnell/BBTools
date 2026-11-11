package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

import assemble.Rollback;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.IntList;

/** Finite one-use restore contracts with literal expected read and count state.
 * Records failures before exiting, so the same oracle can capture a baseline and
 * verify a repair. Snapshot reuse and mutation of public counts0 are not tested.
 * @author Brian Bushnell, Shinobu, Sucrose
 */
public final class RollbackContractTest{

	/** Runs nine independent cases and optional single-use transaction timings.
	 * outdir must be new; iterations=0 disables timing and warmup defaults to100000.
	 * Functional checks remain active under -da and are saved before optional timing. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, RollbackContractTest.class, false);
		try{
			Path root=null;
			int iterations=0, warmup=100000;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
				if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else if(pair[0].equalsIgnoreCase("iterations")){iterations=Parse.parseIntKMG(pair[1]);}
				else if(pair[0].equalsIgnoreCase("warmup")){warmup=Parse.parseIntKMG(pair[1]);}
				else{
					throw new IllegalArgumentException("Expected outdir, iterations or warmup: "+arg);
				}
			}
			if(root==null){throw new IllegalArgumentException("Specify outdir=NEW_DIRECTORY.");}
			if(iterations<0 || warmup<0){throw new IllegalArgumentException("iterations and warmup must be nonnegative.");}
			Files.createDirectory(root);
			final StringBuilder rows=new StringBuilder("case\tpassed\tdetail\n");
			int failed=0;
			for(int mode=0; mode<NAMES.length; mode++){failed+=runCase(mode, rows);}
			Files.write(root.resolve("checks.tsv"), rows.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			if(iterations>0){timeTransactions(root, iterations, warmup);}
			System.out.println("ROLLBACK_CONTRACTS checks="+NAMES.length+" failed="+failed);
			if(failed!=0){throw new AssertionError("Rollback contract failures; inspect "+root);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Mutates saved fields and an unsaved numeric ID, then performs one restore.
	 * Modes 0,1,8 are unaffected controls; the other cases vary destination storage.
	 * @return Zero on the literal contract, or one after recording a failed restore
	 */
	private static int runCase(final int mode, final StringBuilder rows){
		assert(mode>=0 && mode<NAMES.length) : "Each restore scenario has one fixed result-table name.";
		final boolean noQualities=mode==6 || mode==8, noCounts=mode==8;
		final Read read=new Read(new byte[]{'A', 'C', 'G', 'T'},
				noQualities ? null : new byte[]{10, 20, 30, 40}, "saved", 41, false);
		read.flags=0;
		IntList counts=new IntList(8);
		counts.add(2); counts.add(3); counts.add(5);
		final Rollback point=noCounts ? new Rollback(read) : new Rollback(read, counts);
		final byte[] originalBases=read.bases, originalQualities=read.quality;
		final int[] originalCounts=counts.array;
		Arrays.fill(read.bases, (byte)'N');
		if(read.quality!=null){Arrays.fill(read.quality, (byte)1);}
		Arrays.fill(counts.array, -7);
		read.id="changed"; read.flags=1; read.numericID=99;
		if(mode==1){
			read.bases=new byte[]{'N'}; read.quality=new byte[]{1};
			counts.clear(); counts.add(99);
		}else if(mode==2){counts.setSize(1);}
		else if(mode==3){counts.add(99); counts.add(101);}
		else if(mode==4){counts=new IntList(1); counts.add(99);}
		else if(mode==5){read.quality=null;}
		else if(mode==6){read.quality=new byte[]{1, 1, 1, 1};}
		else if(mode==7){read.quality=new byte[]{1};}
		RuntimeException failure=null;
		try{
			if(noCounts){point.rollback(read);}else{point.rollback(read, counts);}
		}catch(RuntimeException caught){failure=caught;}
		final boolean bases=Arrays.equals(read.bases, new byte[]{'A', 'C', 'G', 'T'});
		final boolean qualities=Arrays.equals(read.quality, noQualities ? null : new byte[]{10, 20, 30, 40});
		final boolean values=noCounts || (counts.size==3 && counts.array.length>=3 &&
				counts.array[0]==2 && counts.array[1]==3 && counts.array[2]==5);
		final boolean metadata="saved".equals(read.id) && read.flags==0 && read.numericID==99;
		final boolean reuse=mode!=0 || (read.bases==originalBases &&
				read.quality==originalQualities && counts.array==originalCounts);
		final boolean passed=failure==null && bases && qualities && values && metadata && reuse;
		rows.append(NAMES[mode]).append('\t').append(passed ? 1 : 0).append('\t')
				.append("bases=").append(bases).append(" qualities=").append(qualities)
				.append(" counts=").append(values).append(" metadata=").append(metadata)
				.append(" reuse=").append(reuse).append(" exception=")
				.append(failure==null ? "none" : failure.getClass().getSimpleName()).append('\n');
		if(failure!=null){System.err.println(NAMES[mode]); failure.printStackTrace(System.err);}
		return passed ? 0 : 1;
	}

	/** Times warmup and measured transactions on reused live storage, with a fresh
	 * snapshot each iteration. Full state checks run outside each timed loop. */
	private static void timeTransactions(final Path root, final int iterations, final int warmup) throws Exception{
		final byte[] bases=new byte[150], qualities=new byte[150];
		Arrays.fill(bases, (byte)'A'); Arrays.fill(qualities, (byte)30);
		final Read read=new Read(bases, qualities, "timing", 41, false);
		read.flags=0;
		final IntList counts=new IntList(120);
		for(int i=0; i<120; i++){counts.add(7);}
		final int[] countStorage=counts.array;
		final StringBuilder rows=new StringBuilder("phase\titerations\telapsed_nanos\tchecksum\texpected_checksum\trestored\n");
		boolean valid=true;
		for(int phase=0; phase<2; phase++){
			final int n=phase==0 ? warmup : iterations;
			final long start=System.nanoTime();
			final long checksum=transactions(read, counts, n);
			final long nanos=System.nanoTime()-start;
			final long expected=102L*n;//ASCII A(65)+Q30+count7 after each restore.
			final boolean restored=timingStateRestored(read, counts, bases, qualities, countStorage);
			valid&=checksum==expected && restored;
			rows.append(phase==0 ? "warmup" : "measured").append('\t').append(n).append('\t')
					.append(nanos).append('\t').append(checksum).append('\t').append(expected)
					.append('\t').append(restored ? 1 : 0).append('\n');
		}
		Files.write(root.resolve("timing.tsv"), rows.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(!valid){throw new AssertionError("Timed transactions did not restore literal state/checksum; inspect "+root);}
	}

	/** Executes aligned one-use snapshots; no snapshot is reused after a restore.
	 * Mutations touch rotating live positions, never the snapshot's saved arrays. */
	private static long transactions(final Read read, final IntList counts, final int iterations){
		assert(read.length()==150 && read.quality.length==150 && counts.size==120) :
				"Matched timing uses a150-base read, aligned qualities and120counts for every new snapshot.";
		long checksum=0;
		int base=0, count=0;
		for(int i=0; i<iterations; i++){
			final Rollback point=new Rollback(read, counts);
			read.bases[base]='N'; read.quality[base]=0; counts.array[count]=-7;
			read.id="changed"; read.flags=1;
			point.rollback(read, counts);
			checksum+=read.bases[base]+read.quality[base]+counts.array[count]+read.flags+
					("timing".equals(read.id) ? 0 : 1);
			if(++base==150){base=0;}
			if(++count==120){count=0;}
		}
		return checksum;
	}

	/** Checks full live contents, metadata and original buffer identities outside timing. */
	private static boolean timingStateRestored(final Read read, final IntList counts,
			final byte[] bases, final byte[] qualities, final int[] countStorage){
		if(read.length()!=150 || read.quality==null || read.quality.length!=150 ||
				counts.size!=120 || read.flags!=0 || !"timing".equals(read.id) || read.numericID!=41 ||
				read.bases!=bases || read.quality!=qualities || counts.array!=countStorage){return false;}
		for(int i=0; i<150; i++){if(read.bases[i]!='A' || read.quality[i]!=30){return false;}}
		for(int i=0; i<120; i++){if(counts.array[i]!=7){return false;}}
		return true;
	}

	/** Ordered scenarios; baseline/fixed runner validates the complete set. */
	private static final String[] NAMES={"stable-control", "length-control", "counts-short",
			"counts-long", "counts-capacity", "qualities-null", "qualities-added",
			"qualities-short", "plain-control"};
}
