package test.clump;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import clump.KmerSort;
import stream.Read;

/** Finite oracle for the single-group KmerSort namesort final-ordering branch.
 * @author Jean
 */
public final class KmerSortNameOrderTest{

	/** Records one-group baseline and fixed name-order expectations. */
	public static void main(String[] args) throws Exception{
		check(args.length==1 || args.length==3,
				"Expected phase, with optional runtime and fresh fixture directory");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown phase: "+args[0]);

		System.out.println("case\torder\tpairs\tbaseline_ok\tfixed_ok");
		boolean baselineOk=true, fixedOk=true;
		final Result[] results=new Result[] {
				checkNamesort("unpaired", fixed, unpaired(), false, false),
				checkNamesort("attached-pairs", fixed, attachedPairs(), false, false),
				checkNamesort("attached-repair", fixed, attachedRepair(), true, false),
				checkNamesort("repair-lexical", fixed, repairLexical(), true, true),
				checkNamesort("repair-duplicate-names", fixed, repairDuplicateNames(), true,
						true),
				checkRepairOnly()
		};
		for(Result result : results){
			System.out.println(result);
			baselineOk&=result.baselineOk;
			fixedOk&=result.fixedOk;
		}
		if(args.length>1){
			final Result result=runCli(Paths.get(args[1]).toAbsolutePath(),
					Paths.get(args[2]).toAbsolutePath());
			System.out.println(result);
			baselineOk&=result.baselineOk;
			fixedOk&=result.fixedOk;
		}
		check(fixed ? fixedOk : baselineOk,
				"Unexpected "+args[0]+" single-group namesort behavior");
	}

	/** Invokes the final namesort helper selected for the requested phase. */
	private static Result checkNamesort(String label, boolean fixed, Fixture fixture,
			boolean repair, boolean unpair) throws Exception{
		ArrayList<Read> out=null;
		Throwable thrown=null;
		try{
			out=invoke(method(fixed ? "nameSortById" : "idSort"), fixture.reads,
					fixed ? repair && unpair : repair);
		}catch(Throwable t){thrown=unwrap(t);}
		return evaluate(label, out, fixture, thrown);
	}

	/** Verifies that repair without namesort still leaves clump order untouched. */
	private static Result checkRepairOnly() throws Exception{
		final Fixture fixture=repairOnly();
		return evaluate("repair-only", invoke(READ1_ONLY, fixture.reads), fixture, null);
	}

	/** Computes exact identity and mate expectations without throwing immediately. */
	private static Result evaluate(String label, ArrayList<Read> out, Fixture fixture,
			Throwable thrown){
		final String order=thrown==null ? order(out) : thrown.getClass().getSimpleName();
		final boolean pairs=thrown==null && linkedPairs(fixture.pairs);
		final boolean baselineOk=fixture.baselineThrows ? thrown instanceof AssertionError :
				thrown==null && matches(out, fixture.baseline) && pairs;
		final boolean fixedOk=thrown==null && matches(out, fixture.fixed) && pairs;
		return new Result(label, order, pairs, baselineOk, fixedOk);
	}

	/** Builds an unpaired namesort case whose lexical order opposes input IDs. */
	private static Fixture unpaired(){
		final Read alpha=read("alpha", 3, 0);
		final Read beta=read("beta", 2, 0);
		final Read zeta=read("zeta", 1, 0);
		return new Fixture(list(alpha, beta, zeta), reads(zeta, beta, alpha),
				reads(alpha, beta, zeta), new Read[0][], false);
	}

	/** Builds an attached-pair case that must sort by read-one name and keep mates. */
	private static Fixture attachedPairs(){
		final Read alpha1=read("alpha/1", 3, 0), alpha2=read("alpha/2", 3, 1);
		final Read beta1=read("beta/1", 2, 0), beta2=read("beta/2", 2, 1);
		final Read zeta1=read("zeta/1", 1, 0), zeta2=read("zeta/2", 1, 1);
		link(alpha1, alpha2); link(beta1, beta2); link(zeta1, zeta2);
		return new Fixture(list(alpha1, beta1, zeta1), reads(zeta1, beta1, alpha1),
				reads(alpha1, beta1, zeta1), pairs(alpha1, alpha2, beta1, beta2,
						zeta1, zeta2), false);
	}

	/** Builds an attached-pair repair case that the old idSort repair path asserts. */
	private static Fixture attachedRepair(){
		final Read alpha1=read("alpha/1", 3, 0), alpha2=read("alpha/2", 3, 1);
		final Read beta1=read("beta/1", 2, 0), beta2=read("beta/2", 2, 1);
		final Read zeta1=read("zeta/1", 1, 0), zeta2=read("zeta/2", 1, 1);
		link(alpha1, alpha2); link(beta1, beta2); link(zeta1, zeta2);
		return new Fixture(list(alpha1, beta1, zeta1), null,
				reads(alpha1, beta1, zeta1), pairs(alpha1, alpha2, beta1, beta2,
						zeta1, zeta2), true);
	}

	/** Builds a detached-pair repair case whose pair name order opposes input IDs. */
	private static Fixture repairLexical(){
		final Read alpha1=read("alpha/1", 3, 0), alpha2=read("alpha/2", 3, 1);
		final Read beta1=read("beta/1", 2, 0), beta2=read("beta/2", 2, 1);
		final Read zeta1=read("zeta/1", 1, 0), zeta2=read("zeta/2", 1, 1);
		return new Fixture(list(alpha1, beta1, zeta1, alpha2, beta2, zeta2),
				reads(zeta1, beta1, alpha1), reads(alpha1, beta1, zeta1),
				pairs(alpha1, alpha2, beta1, beta2, zeta1, zeta2), false);
	}

	/** Builds a repair case that catches name-based pairing of nonunique names. */
	private static Fixture repairDuplicateNames(){
		final Read first1=read("duplicate", 1, 0), first2=read("duplicate", 1, 1);
		final Read second1=read("duplicate", 2, 0), second2=read("duplicate", 2, 1);
		return new Fixture(list(first1, second1, first2, second2),
				reads(first1, second1), reads(first1, second1),
				pairs(first1, first2, second1, second2), false);
	}

	/** Builds a repair-only control where clump order differs from input order. */
	private static Fixture repairOnly(){
		final Read alpha1=read("alpha/1", 1, 0), alpha2=read("alpha/2", 1, 1);
		final Read zeta1=read("zeta/1", 2, 0), zeta2=read("zeta/2", 2, 1);
		link(zeta1, zeta2); link(alpha1, alpha2);
		return new Fixture(list(zeta1, alpha1), reads(zeta1, alpha1),
				reads(zeta1, alpha1), pairs(zeta1, zeta2, alpha1, alpha2), false);
	}

	/** Runs one groups=1 CLI dispatch case through clumpify.sh when requested. */
	private static Result runCli(Path runtime, Path root) throws Exception{
		Files.createDirectory(root);
		final Path dir=root.resolve("cli-attached-repair");
		Files.createDirectory(dir);
		final Path left=dir.resolve("left.fq"), right=dir.resolve("right.fq");
		write(left, ZETA_1+ALPHA_1);
		write(right, ZETA_2+ALPHA_2);
		final Path output=dir.resolve("out#.fq");
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(runtime.resolve("clumpify.sh").toString());
		command.add("-Xmx512m");
		command.add("-Xms64m");
		command.add("t=1");
		command.add("k=7");
		command.add("seed=1");
		command.add("hashes=4");
		command.add("border=0");
		command.add("qin=33");
		command.add("qout=33");
		command.add("rcomp=f");
		command.add("ecco=f");
		command.add("correct=f");
		command.add("dedupe=f");
		command.add("condense=f");
		command.add("quantize=f");
		command.add("reorder=f");
		command.add("overwrite=f");
		command.add("deleteinput=f");
		command.add("groups=1");
		command.add("namesort=t");
		command.add("repair=t");
		command.add("unpair=f");
		command.add("in="+left);
		command.add("in2="+right);
		command.add("out="+output);
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8,
				StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.directory(dir.toFile());
		builder.environment().put("PATH", System.getProperty("java.home")+"/bin:"+
				builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		final int status=builder.start().waitFor();
		final boolean exact=equalsFile(dir.resolve("out1.fq"), ALPHA_1+ZETA_1) &&
				equalsFile(dir.resolve("out2.fq"), ALPHA_2+ZETA_2);
		final boolean inputExact=equalsFile(left, ZETA_1+ALPHA_1) &&
				equalsFile(right, ZETA_2+ALPHA_2);
		final String stderr=new String(Files.readAllBytes(dir.resolve("stderr.txt")),
				StandardCharsets.UTF_8);
		final boolean baselineOk=status!=0 && inputExact &&
				stderr.contains("AssertionError") &&
				stderr.contains("clump.KmerSort.idSort");
		final boolean fixedOk=status==0 && inputExact && exact;
		return new Result("cli-attached-repair", "exit="+status, exact && inputExact,
				baselineOk, fixedOk);
	}

	/** True only when the top-level list has exact object identity and order. */
	private static boolean matches(ArrayList<Read> out, Read[] expected){
		if(out.size()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			if(out.get(i)!=expected[i]){return false;}
		}
		return true;
	}

	/** True only when every expected mate pair remains linked by identity. */
	private static boolean linkedPairs(Read[][] pairs){
		for(Read[] pair : pairs){
			if(pair[0].mate!=pair[1] || pair[1].mate!=pair[0]){return false;}
			if(pair[0].pairnum()!=0 || pair[1].pairnum()!=1){return false;}
		}
		return true;
	}

	/** Returns read identifiers in final top-level order for the result table. */
	private static String order(ArrayList<Read> out){
		final StringBuilder sb=new StringBuilder();
		for(Read read : out){
			if(sb.length()>0){sb.append(',');}
			sb.append(read.id).append('#').append(read.numericID);
		}
		return sb.toString();
	}

	/** Invokes one static KmerSort helper with a repair argument. */
	@SuppressWarnings("unchecked")
	private static ArrayList<Read> invoke(Method method, ArrayList<Read> reads,
			boolean repair) throws Exception{
		return (ArrayList<Read>)method.invoke(null, reads, Boolean.valueOf(repair));
	}

	/** Invokes the static read-one projection helper. */
	@SuppressWarnings("unchecked")
	private static ArrayList<Read> invoke(Method method, ArrayList<Read> reads)
			throws Exception{
		return (ArrayList<Read>)method.invoke(null, reads);
	}

	/** Makes a tiny read with deterministic pair metadata. */
	private static Read read(String id, long numericID, int pairnum){
		final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null,
				id, numericID);
		read.setPairnum(pairnum);
		return read;
	}

	/** Connects existing read objects as mates. */
	private static void link(Read read1, Read read2){
		read1.mate=read2;
		read2.mate=read1;
	}

	/** Packs reads into a mutable list. */
	private static ArrayList<Read> list(Read... reads){
		final ArrayList<Read> list=new ArrayList<Read>(reads.length);
		for(Read read : reads){list.add(read);}
		return list;
	}

	/** Packs reads into an array. */
	private static Read[] reads(Read... reads){
		return reads;
	}

	/** Packs alternating read1/read2 arguments into expected mate pairs. */
	private static Read[][] pairs(Read... reads){
		check(reads.length%2==0, "Expected an even number of pair reads");
		final Read[][] pairs=new Read[reads.length/2][2];
		for(int i=0, j=0; i<reads.length; i+=2, j++){
			pairs[j][0]=reads[i];
			pairs[j][1]=reads[i+1];
		}
		return pairs;
	}

	/** Lightweight assertion independent of JVM assertion flags. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Finds an otherwise package-private static helper. */
	private static Method method(String name){
		try{
			final Method method=KmerSort.class.getDeclaredMethod(name, ArrayList.class,
					boolean.class);
			method.setAccessible(true);
			return method;
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Finds the otherwise package-private static read-one projection helper. */
	private static Method read1OnlyMethod(){
		try{
			final Method method=KmerSort.class.getDeclaredMethod("read1Only",
					ArrayList.class);
			method.setAccessible(true);
			return method;
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Unwraps reflective invocation failures. */
	private static Throwable unwrap(Throwable thrown){
		if(thrown instanceof java.lang.reflect.InvocationTargetException){
			return ((java.lang.reflect.InvocationTargetException)thrown).getCause();
		}
		return thrown;
	}

	/** Writes one FASTQ fixture. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
	}

	/** Compares exact bytes when a file is present. */
	private static boolean equalsFile(Path path, String expected) throws Exception{
		return Files.exists(path) && java.util.Arrays.equals(Files.readAllBytes(path),
				expected.getBytes(StandardCharsets.US_ASCII));
	}

	/** Recorded fixture inputs and exact expectations. */
	private static final class Fixture{

		Fixture(ArrayList<Read> reads_, Read[] baseline_, Read[] fixed_, Read[][] pairs_,
				boolean baselineThrows_){
			reads=reads_;
			baseline=baseline_;
			fixed=fixed_;
			pairs=pairs_;
			baselineThrows=baselineThrows_;
		}

		final ArrayList<Read> reads;
		final Read[] baseline;
		final Read[] fixed;
		final Read[][] pairs;
		final boolean baselineThrows;
	}

	/** Recorded result for one single-group final-ordering case. */
	private static final class Result{

		Result(String label_, String order_, boolean pairs_, boolean baselineOk_,
				boolean fixedOk_){
			label=label_;
			order=order_;
			pairs=pairs_;
			baselineOk=baselineOk_;
			fixedOk=fixedOk_;
		}

		@Override
		public String toString(){
			return label+'\t'+order+'\t'+pairs+'\t'+baselineOk+'\t'+fixedOk;
		}

		final String label;
		final String order;
		final boolean pairs;
		final boolean baselineOk;
		final boolean fixedOk;
	}

	private static final Method READ1_ONLY=read1OnlyMethod();

	private static final String ALPHA_1="@alpha/1\nTTGCAACCGGTT\n+\nHHHHHHHHHHHH\n";
	private static final String ALPHA_2="@alpha/2\nAAAACCCCTTTT\n+\nFFFFFFFFFFFF\n";
	private static final String ZETA_1="@zeta/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n";
	private static final String ZETA_2="@zeta/2\nTTTTCCCCAAAA\n+\nGGGGGGGGGGGG\n";

	private KmerSortNameOrderTest(){}

}
