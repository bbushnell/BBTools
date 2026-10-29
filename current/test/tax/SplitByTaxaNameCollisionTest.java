package test.tax;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ReadWrite;
import shared.Timer;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tax.SplitByTaxa;
import tax.TaxTree;

/** Finite baseline for SplitByTaxa raw-name keys colliding after filename sanitization.
 * @author Sucrose */
public final class SplitByTaxaNameCollisionTest{

	/** Runs a healthy distinct-path control and the collision case.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If setup, reflective invocation or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final CaseResult healthy=runCase("gamma delta", "gamma-delta", 3, 4, false);
		check(healthy.failure==null, "Distinct sanitized paths must process: "+healthy.failure);
		check(healthy.rawWriterCount==2, "Distinct raw names must get separate writers");
		check(healthy.uniquePathCount==2, "Distinct sanitized paths must stay distinct");
		check(healthy.returnedLists==2, "Healthy input stream should see two returned lists");

		final CaseResult collision=runCase("alpha beta", "alpha_beta", 5, 6, true);
		if(baseline){
			check(collision.failure==null,
					"Baseline must expose the silent sanitized-path collision: "+collision.failure);
			check(collision.rawWriterCount==2,
					"Baseline must open one writer per raw colliding name");
			check(collision.uniquePathCount==1,
					"Baseline colliding writers must target the same sanitized path");
			check(collision.returnedLists==2,
					"Baseline collision stream should see two returned lists");
		}else{
			check(collision.failure instanceof IllegalArgumentException,
					"Fixed SplitByTaxa should reject with IllegalArgumentException: "+
					collision.failure);
			checkCollisionMessage(collision.failure.getMessage(), "alpha beta",
					"alpha_beta", "split-alpha_beta.fq");
			checkPublicProcessCleanup();
			checkPairedPathCollision();
		}
		System.out.println("SPLIT_BY_TAXA_NAME_COLLISION_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private SplitByTaxaNameCollisionTest(){}

	/** Runs processInner over two read names and captures the writers it constructs.
	 * @param name1 First taxon name
	 * @param name2 Second taxon name
	 * @param taxID1 TaxID matching name1 in the fixture tree
	 * @param taxID2 TaxID matching name2 in the fixture tree
	 * @param colliding Whether the case should collide after filename sanitization
	 * @return Writer counts and any processInner failure
	 * @throws Exception If setup or cleanup fails */
	private static CaseResult runCase(final String name1, final String name2,
			final int taxID1, final int taxID2, final boolean colliding)
			throws Exception{
		final Path root=Files.createTempDirectory("split-by-taxa-collision-");
		final HashMap<String, ConcurrentReadOutputStream> map=
				new HashMap<String, ConcurrentReadOutputStream>();
		final FakeCris cris=new FakeCris(taxID1, taxID2);
		Throwable failure=null;
		try{
			final Path tree=root.resolve("tree.tsv");
			final Path input=root.resolve("input.fa");
			Files.write(tree, treeText().toBytes());
			Files.write(input, new byte[0]);
			final SplitByTaxa split=new SplitByTaxa(new String[] {
					"in="+input,
					"out="+root.resolve("split-%.fq"),
					"tree="+tree,
					"level=species",
					"ow=t"
			});
			final TaxTree loaded=(TaxTree)TREE_FIELD.get(split);
			check(name1.equals(loaded.getNode(taxID1).name),
					"Loaded tree did not preserve taxID "+taxID1+
					" literal name");
			check(name2.equals(loaded.getNode(taxID2).name),
					"Loaded tree did not preserve taxID "+taxID2+
					" literal name");
			check(colliding==(sanitize(name1).equals(sanitize(name2))),
					"Fixture collision flag does not match literal names");
			try{
				PROCESS_INNER.invoke(split, cris, map);
			}catch(InvocationTargetException e){
				failure=e.getCause();
			}
			final HashSet<String> paths=new HashSet<String>();
			for(ConcurrentReadOutputStream ros : map.values()){paths.add(ros.ff1.name());}
			return new CaseResult(failure, map.size(), paths.size(), cris.returnedLists);
		}finally{
			boolean closeError=false;
			try{
				for(ConcurrentReadOutputStream ros : map.values()){
					closeError|=ReadWrite.closeStream(ros);
				}
			}finally{
				deleteTempTree(root);
			}
			check(!closeError, "A SplitByTaxa test writer failed to close cleanly");
		}
	}

	/**
	 * Verifies the public {@code process()} path closes earlier writers when a
	 * later taxon collides before its writer can open.
	 * @throws Exception If fixture setup or cleanup fails
	 */
	private static void checkPublicProcessCleanup() throws Exception{
		final Path root=Files.createTempDirectory("split-by-taxa-process-collision-");
		try{
			final Path tree=root.resolve("tree.tsv");
			final Path input=root.resolve("input.fa");
			Files.write(tree, treeText().toBytes());
			writeReads(input, 3, 5, 6);
			final SplitByTaxa split=new SplitByTaxa(new String[] {
					"in="+input,
					"out="+root.resolve("split-%.fa"),
					"tree="+tree,
					"level=species",
					"ow=t"
			});
			final Throwable failure=processFailure(split);
			check(failure instanceof IllegalArgumentException,
					"Public process should reject colliding output paths: "+failure);
			checkCollisionMessage(failure.getMessage(), "alpha beta", "alpha_beta",
					"split-alpha_beta.fa");
			check(readAscii(root.resolve("split-gamma_delta.fa")).equals(
					">tid|3|\nACGT\n"),
					"Noncolliding writer should contain the complete first record");
			check(readAscii(root.resolve("split-alpha_beta.fa")).equals(
					">tid|5|\nTGCA\n"),
					"First colliding writer should contain the complete second record");
		}finally{
			deleteTempTree(root);
		}
	}

	/**
	 * Verifies a single new taxon writer cannot point {@code out1} and
	 * {@code out2} at the same normalized path.
	 * @throws Exception If fixture setup or cleanup fails
	 */
	private static void checkPairedPathCollision() throws Exception{
		final Path root=Files.createTempDirectory("split-by-taxa-paired-collision-");
		try{
			final Path tree=root.resolve("tree.tsv");
			final Path input=root.resolve("input.fa");
			Files.write(tree, treeText().toBytes());
			writeReads(input, 3, 3);
			final SplitByTaxa split=new SplitByTaxa(new String[] {
					"in="+input,
					"out1="+root.resolve("%_gamma_delta.fq"),
					"out2="+root.resolve("gamma_delta_%.fq"),
					"tree="+tree,
					"level=species",
					"ow=t"
			});
			final Throwable failure=processFailure(split);
			check(failure instanceof IllegalArgumentException,
					"Paired output should reject converging mate paths: "+failure);
			final String message=failure.getMessage();
			check(message.contains("mate output paths resolve to the same file"),
					"Paired collision message should identify mate outputs: "+message);
			check(message.contains("gamma delta") &&
					message.contains("gamma_delta_gamma_delta.fq"),
					"Paired collision message should identify taxon/path: "+message);
			check(!Files.exists(root.resolve("gamma_delta_gamma_delta.fq")),
					"Paired collision should fail before opening the shared writer");
		}finally{
			deleteTempTree(root);
		}
	}

	/** Runs public process and returns its failure, or null on success.
	 * @param split Configured splitter to run
	 * @return Process failure, or null if it returned normally */
	private static Throwable processFailure(final SplitByTaxa split){
		try{
			split.process(new Timer());
			return null;
		}catch(Throwable failure){
			return failure;
		}
	}

	/** Writes finite FASTA input with headers parseable as explicit taxIDs.
	 * @param input Destination input path
	 * @param taxIDs TaxIDs to write as FASTA read headers
	 * @throws Exception If writing fails */
	private static void writeReads(final Path input, final int... taxIDs)
			throws Exception{
		final ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<taxIDs.length; i++){
			bb.append(">tid|").append(taxIDs[i]).append("|").nl();
			bb.append(i%2==0 ? "ACGT" : "TGCA").nl();
		}
		Files.write(input, bb.toBytes());
	}

	/** Reads an ASCII text file.
	 * @param path Path to read
	 * @return File text decoded as ASCII
	 * @throws Exception If reading fails */
	private static String readAscii(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Builds one tree containing both test cases under the same parent.
	 * @return Serialized TaxTreeText bytes */
	private static ByteBuilder treeText(){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t7\t6\n");
		bb.append("##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n");
		bb.append("##LevelExtended\tTID\tPID\tNumChildren\tMinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		row(bb, 1, 1, TaxTree.LIFE_E, 1, "Life");
		row(bb, 2, 1, TaxTree.DOMAIN_E, 4, "Bacteria");
		row(bb, 3, 2, TaxTree.SPECIES_E, 0, "gamma delta");
		row(bb, 4, 2, TaxTree.SPECIES_E, 0, "gamma-delta");
		row(bb, 5, 2, TaxTree.SPECIES_E, 0, "alpha beta");
		row(bb, 6, 2, TaxTree.SPECIES_E, 0, "alpha_beta");
		return bb;
	}

	/** Mirrors SplitByTaxa's output-name substitution.
	 * @param name Raw taxon name
	 * @return Sanitized path component */
	private static String sanitize(final String name){
		return name.replaceAll("\\s+", "_").replaceAll("[/\\\\]", "");
	}

	/** Appends one A48 TaxTreeText row for the fixed test tree.
	 * @param bb Destination tree text
	 * @param id Node TaxID
	 * @param parent Parent TaxID
	 * @param rank Extended taxonomic rank
	 * @param children Number of direct children
	 * @param name Literal taxon name, including any filename-sensitive characters */
	private static void row(final ByteBuilder bb, final int id, final int parent,
			final int rank, final int children, final String name){
		bb.appendA48(rank).tab().appendA48(id).tab().appendA48(parent)
				.tab().appendA48(children).tab().appendA48(TaxTree.LIFE_E)
				.tab().appendA48(0).tab().appendA48(rank).tab().append(name).nl();
	}

	/** Deletes an owned temporary tree after the tested writers have closed.
	 * @param root Temporary directory whose descendants are owned by this fixture
	 * @throws Exception If directory traversal or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps test checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Checks fixed collision diagnostics include both raw names and the path.
	 * @param message Exception message
	 * @param name1 First raw taxon name expected in the message
	 * @param name2 Second raw taxon name expected in the message
	 * @param path Output path expected in the message */
	private static void checkCollisionMessage(final String message,
			final String name1, final String name2, final String path){
		check(message.startsWith(
				"Distinct taxon names produce the same SplitByTaxa output path: "),
				"Collision message should identify SplitByTaxa output collision: "+
				message);
		check(message.contains(name1) && message.contains(name2) &&
				message.contains(path),
				"Collision message should identify both taxon names and path: "+
				message);
	}

	/** Reflected package-private SplitByTaxa.processInner method. */
	private static final Method PROCESS_INNER;
	/** Reflected private SplitByTaxa.tree field. */
	private static final java.lang.reflect.Field TREE_FIELD;

	static{
		try{
			PROCESS_INNER=SplitByTaxa.class.getDeclaredMethod("processInner",
					ConcurrentReadInputStream.class, HashMap.class);
			PROCESS_INNER.setAccessible(true);
			TREE_FIELD=SplitByTaxa.class.getDeclaredField("tree");
			TREE_FIELD.setAccessible(true);
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Captured behavior from one processInner invocation. */
	private static final class CaseResult{
		/** Captures the outcome without changing the tested streams.
		 * @param failure_ Cause from reflective invocation, or null on normal return
		 * @param rawWriterCount_ Number of raw-name entries in the output map
		 * @param uniquePathCount_ Number of distinct first-mate output path strings
		 * @param returnedLists_ Number of returnList calls observed by the fake input */
		CaseResult(final Throwable failure_, final int rawWriterCount_,
				final int uniquePathCount_, final int returnedLists_){
			failure=failure_;
			rawWriterCount=rawWriterCount_;
			uniquePathCount=uniquePathCount_;
			returnedLists=returnedLists_;
		}
		/** Invocation failure, or null if processInner returned normally. */
		final Throwable failure;
		/** Output-map entries, counted before stream cleanup. */
		final int rawWriterCount;
		/** Distinct first-mate path strings among the constructed output streams. */
		final int uniquePathCount;
		/** Input-list returns observed during processing. */
		final int returnedLists;
	}

	/** Finite read stream that yields exactly one nonempty list and one empty tail. */
	private static final class FakeCris extends ConcurrentReadInputStream{
		/** Creates two unpaired reads with explicit taxonomy headers.
		 * @param taxID1 TaxID embedded in the first read's header
		 * @param taxID2 TaxID embedded in the second read's header */
		FakeCris(final int taxID1, final int taxID2){
			super("SplitByTaxaNameCollisionTest");
			reads.add(new Read("ACGT".getBytes(StandardCharsets.US_ASCII),
					null, "tid|"+taxID1+"|", 1));
			reads.add(new Read("TGCA".getBytes(StandardCharsets.US_ASCII),
					null, "tid|"+taxID2+"|", 2));
		}
		@Override
		public ListNum<Read> nextList(){
			if(next==0){next++; return new ListNum<Read>(reads, 0);}
			if(next==1){next++; return new ListNum<Read>(new ArrayList<Read>(), 1);}
			return null;
		}
		@Override
		public void returnList(final long listNum, final boolean poison){returnedLists++;}
		@Override
		public void run(){}
		@Override
		public void shutdown(){}
		@Override
		public void restart(){next=0;}
		@Override
		public void close(){}
		@Override
		public boolean paired(){return false;}
		@Override
		public Object[] producers(){return new Object[0];}
		@Override
		public boolean errorState(){return false;}
		@Override
		public void setSampleRate(final float rate, final long seed){}
		@Override
		public long basesIn(){return 0;}
		@Override
		public long readsIn(){return 0;}
		@Override
		public boolean verbose(){return false;}

		/** Owned two-read payload returned as the first list. */
		final ArrayList<Read> reads=new ArrayList<Read>();
		/** Delivery cursor and cumulative returnList call count; restart resets only the cursor. */
		int next=0, returnedLists=0;
	}

}
