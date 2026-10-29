package test.tax;

import java.lang.reflect.Constructor;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Stream;

import structures.IntHashMap;
import tax.TaxNode;
import tax.TaxTree;

/** Captures invalid QuickClade parent handling for leaf-only FakeTree nodes.
 * @author Jean
 */
public final class FakeTreeLeafTaxonomyInvalidQcTest{

	/** Creates an instance with no retained state. */
	public FakeTreeLeafTaxonomyInvalidQcTest(){}

	/** Positive in-range taxID slot left absent from the tiny tree. */
	private static final int UNKNOWN_TID=TaxTree.BACTERIA_ID+1;
	/** Usable domain-level QuickClade evidence for the single fake leaf. */
	private static final String USABLE_QC_ROW=
			"9000001\t2\tBacteria\tProteobacteria\tdomain\tdomain:0.99\t99.0\n";

	/** Runs invalid-evidence checks and optional broad-SSU precedence checks.
	 * @param args Ignored BBTools root, expect=baseline/fixed, and optional unknown=baseline/fixed
	 * @throws Exception If temporary fixture creation, child execution or output parsing fails */
	public static void main(String[] args) throws Exception{
		check(args.length==2 || args.length==3,
				"Expected <BBTools dir> expect=baseline|fixed [unknown=baseline|fixed]");
		final boolean fixed=parseSelector(args[1], "expect");
		final boolean runUnknown=(args.length==3);
		final boolean unknownFixed=runUnknown && parseSelector(args[2], "unknown");
		final Path root=Files.createTempDirectory("faketree-invalid-qc-");

		final Result qcEa=run(root, true, Case.INVALID_QC);
		final Result ssuEa=run(root, true, Case.ZERO_SSU);
		if(fixed){
			assertFallback(qcEa, "invalid-QC -ea");
			assertValidQc(ssuEa, "invalid-SSU -ea");
		}else{
			assertInvalidParentFailure(qcEa, "invalid-QC -ea");
			assertInvalidParentFailure(ssuEa, "invalid-SSU -ea");
		}

		final Result qcDa=run(root, false, Case.INVALID_QC);
		final Result ssuDa=run(root, false, Case.ZERO_SSU);
		if(fixed){
			assertFallback(qcDa, "invalid-QC -da");
			assertValidQc(ssuDa, "invalid-SSU -da");
		}else{
			assertBlankParent(qcDa, "invalid-QC -da", "quickclade");
			assertBlankParent(ssuDa, "invalid-SSU -da", "findssu");
		}

		final Result positiveEa, positiveDa;
		if(runUnknown){
			positiveEa=run(root, true, Case.POSITIVE_UNKNOWN_SSU);
			positiveDa=run(root, false, Case.POSITIVE_UNKNOWN_SSU);
			if(unknownFixed){
				assertValidQc(positiveEa, "positive-unknown-SSU -ea");
				assertValidQc(positiveDa, "positive-unknown-SSU -da");
			}else{
				assertUnknownSsuLife(positiveEa, "positive-unknown-SSU -ea");
				assertUnknownSsuLife(positiveDa, "positive-unknown-SSU -da");
			}
		}else{
			positiveEa=positiveDa=null;
		}

		check(qcEa.inputsPreserved && ssuEa.inputsPreserved && qcDa.inputsPreserved &&
				ssuDa.inputsPreserved && (!runUnknown ||
				(positiveEa.inputsPreserved && positiveDa.inputsPreserved)),
				"Taxonomy assignment must not mutate input tables.");
		printReceipt(Case.INVALID_QC, true, qcEa);
		printReceipt(Case.ZERO_SSU, true, ssuEa);
		if(runUnknown){printReceipt(Case.POSITIVE_UNKNOWN_SSU, true, positiveEa);}
		printReceipt(Case.INVALID_QC, false, qcDa);
		printReceipt(Case.ZERO_SSU, false, ssuDa);
		if(runUnknown){printReceipt(Case.POSITIVE_UNKNOWN_SSU, false, positiveDa);}
		if(runUnknown){
			for(Case testCase : Case.values()){
				if(testCase.baselineParent==0){continue;}
				for(boolean assertions : new boolean[]{true, false}){
					final Result result=run(root, assertions, testCase);
					assertParentPolicy(result, testCase, unknownFixed);
					printReceipt(testCase, assertions, result);
				}
			}
		}
		deleteRecursively(root);
		System.out.println("FAKETREE_INVALID_QC_TEST_PASS expect="+(fixed ? "fixed" : "baseline")+
				(runUnknown ? " unknown="+(unknownFixed ? "fixed" : "baseline") : ""));
	}

	/** Verifies -ea still exposes the captured blank-parent baseline.
	 * @param result Captured child exit and logs
	 * @param name Assertion-mode label for diagnostics */
	private static void assertInvalidParentFailure(final Result result, final String name){
		check(result.exit!=0, "Baseline "+name+" should catch a blank invalid parent: "+result.log);
		check(result.logText.contains("No proposed parent for fake_taxid 9000001"),
				"Baseline "+name+" must explain the invalid parent assertion: "+result.log);
	}

	/** Verifies -da exposes the captured blank-parent baseline.
	 * @param result Captured child output tables and logs
	 * @param name Assertion-mode label for diagnostics
	 * @param source Expected evidence_source from the accepted invalid row */
	private static void assertBlankParent(final Result result, final String name, final String source){
		check(result.exit==0, "Baseline "+name+" should expose the emitted blank parent: "+result.log);
		check("".equals(result.proposal.get("proposed_parent_taxid")),
				"Baseline "+name+" must emit the missing parent before the repair: "+result.proposals);
		check(source.equals(result.proposal.get("evidence_source")),
				"Baseline "+name+" should reveal that "+source+" was selected: "+result.proposals);
	}

	/** Verifies invalid QuickClade taxids are treated like absent usable evidence.
	 * @param result Captured child output tables and logs
	 * @param name Assertion-mode label for diagnostics */
	private static void assertFallback(final Result result, final String name){
		check(result.exit==0, name+" invalid-QC fallback should succeed: "+result.log);
		check("2".equals(result.proposal.get("proposed_parent_taxid")),
				name+" invalid QC should fall back to Bacteria: "+result.proposals);
		check("fallback".equals(result.proposal.get("evidence_source")),
				name+" invalid QC should not be accepted as QuickClade evidence: "+result.proposals);
		check("unresolved".equals(result.proposal.get("decision_class")),
				name+" invalid QC should use the normal unresolved fallback decision: "+result.proposals);
		check("2".equals(result.node.get("parent_taxid")),
				name+" emitted node parent must match the fallback proposal: "+result.nodes);
	}

	/** Verifies invalid SSU taxids are skipped so usable QuickClade evidence can win.
	 * @param result Captured child output tables and logs
	 * @param name Assertion-mode label for diagnostics */
	private static void assertValidQc(final Result result, final String name){
		check(result.exit==0, name+" invalid-SSU fallback should succeed: "+result.log);
		check("2".equals(result.proposal.get("proposed_parent_taxid")),
				name+" invalid SSU should fall back to the valid QC row: "+result.proposals);
		check("quickclade".equals(result.proposal.get("evidence_source")),
				name+" invalid SSU should not block QuickClade evidence: "+result.proposals);
		check("quickclade_only".equals(result.proposal.get("decision_class")),
				name+" valid QC should use the QuickClade-only decision: "+result.proposals);
		check("2".equals(result.node.get("parent_taxid")),
				name+" emitted node parent must match the QC proposal: "+result.nodes);
	}

	/** Verifies a positive absent SSU taxid resolves to Life and beats usable QuickClade.
	 * @param result Captured child output tables and logs
	 * @param name Assertion-mode label for diagnostics */
	private static void assertUnknownSsuLife(final Result result, final String name){
		check(result.exit==0, name+" positive-unknown SSU baseline should succeed: "+result.log);
		check("1".equals(result.proposal.get("ssu_parent_taxid")),
				name+" baseline should map unknown positive SSU to Life: "+result.proposals);
		check("2".equals(result.proposal.get("qc_parent_taxid")),
				name+" baseline should keep the usable QC parent visible: "+result.proposals);
		check("1".equals(result.proposal.get("proposed_parent_taxid")),
				name+" baseline should select the Life SSU parent: "+result.proposals);
		check("findssu".equals(result.proposal.get("evidence_source")),
				name+" baseline should prefer the unknown positive SSU hit: "+result.proposals);
		check("ssu_with_qc_conflict".equals(result.proposal.get("decision_class")),
				name+" baseline should expose the SSU/QC conflict decision: "+result.proposals);
		check("1".equals(result.node.get("parent_taxid")),
				name+" emitted node parent must match the Life proposal: "+result.nodes);
	}

	/** Checks literal parent/source expectations against the emitted proposal and node.
	 * @param result Captured child output
	 * @param testCase Independent expected IDs for the tiny lineage
	 * @param fixed True for QuickClade precedence over broad SSU parents */
	private static void assertParentPolicy(final Result result, final Case testCase,
			final boolean fixed){
		final int expected=fixed ? testCase.fixedParent : testCase.baselineParent;
		final String source=fixed ? testCase.fixedEvidence : "findssu";
		check(result.exit==0, testCase.receiptName+" child failed: "+result.logText);
		check(Integer.toString(expected).equals(result.proposal.get("proposed_parent_taxid")),
				testCase.receiptName+" expected parent "+expected+": "+result.proposals);
		check(Integer.toString(expected).equals(result.node.get("parent_taxid")),
				testCase.receiptName+" node must use the selected proposal parent: "+result.nodes);
		check(source.equals(result.proposal.get("evidence_source")),
				testCase.receiptName+" expected evidence "+source+": "+result.proposals);
		check(result.inputsPreserved, testCase.receiptName+" mutated an input table");
	}

	/** Runs the tool in a child JVM on a one-row fixture with selected invalid evidence.
	 * @param root Parent temp directory for this fixture run
	 * @param assertions True for a child JVM with assertions enabled
	 * @param testCase Evidence case to write
	 * @return Captured output tables, logs, exit code and input preservation state
	 * @throws Exception If fixture writes, child execution or output parsing fails */
	private static Result run(final Path root, final boolean assertions, final Case testCase)
			throws Exception{
		final Path dir=Files.createDirectory(root.resolve(testCase.dirPrefix+
				(assertions ? "ea" : "da")));
		final Path leaf=dir.resolve("leaf.tsv");
		final Path ssu=dir.resolve("ssu.tsv");
		final Path qc=dir.resolve("qc.tsv");
		final Path tree=dir.resolve("tree.tsv");
		final Path proposals=dir.resolve("proposals.tsv");
		final Path nodes=dir.resolve("nodes.tsv");
		final Path log=dir.resolve("child.log");
		final byte[] leafBytes=("fake_taxid\tcollection\torigin\tsource_id\tgenome_id\tparent_taxid\tdomain\tphylum\n"+
				"9000001\tfixture\torigin\tgenomeA\tgenomeA\t2\tBacteria\tProteobacteria\n").getBytes(StandardCharsets.UTF_8);
		final byte[] ssuBytes=("fake_taxid\thit_tid\thit_name\thit_phylum\tani\twkid\tmatches\n"+
				testCase.ssuRow).getBytes(StandardCharsets.UTF_8);
		final byte[] qcBytes=("fake_taxid\tref_taxid\tref_name\tqc_hit_phylum\tsketch_lca\tconf_level\tANI\n"+
				testCase.qcRow).getBytes(StandardCharsets.UTF_8);
		Files.write(leaf, leafBytes);
		Files.write(ssu, ssuBytes);
		Files.write(qc, qcBytes);
		writeTinyTree(tree);

		final List<String> command=javaCommand(assertions, "in="+leaf, "ssu="+ssu,
				"quickclade="+qc, "out="+proposals, "outnodes="+nodes,
				"tree="+tree, "t=1", "usepigz=f", "useunpigz=f");
		final Process child=new ProcessBuilder(command).redirectErrorStream(true).redirectOutput(log.toFile())
				.start();
		final int exit;
		try{
			child.getOutputStream().close();
			exit=waitForChild(child, log);
		}finally{
			if(child.isAlive()){terminateChild(child, log);}
		}
		final Map<String, String> proposal=Files.exists(proposals) && Files.size(proposals)>0 ?
				readFirstRow(proposals) : new LinkedHashMap<String, String>();
		final Map<String, String> node=Files.exists(nodes) && Files.size(nodes)>0 ?
				readFirstRow(nodes) : new LinkedHashMap<String, String>();
		return new Result(exit, log, new String(Files.readAllBytes(log), StandardCharsets.UTF_8),
				Files.exists(proposals) ? new String(Files.readAllBytes(proposals), StandardCharsets.UTF_8) : "",
				Files.exists(nodes) ? new String(Files.readAllBytes(nodes), StandardCharsets.UTF_8) : "",
				proposal, node, Arrays.equals(leafBytes, Files.readAllBytes(leaf)) &&
				Arrays.equals(ssuBytes, Files.readAllBytes(ssu)) &&
				Arrays.equals(qcBytes, Files.readAllBytes(qc)) && Files.size(tree)>0);
	}

	/** Writes a sparse lineage: species10, genus11, family12, Bacteria2, Life1.
	 * TaxID3 remains an in-range hole for missing-node cases.
	 * @param tree Destination serialized TaxTree path
	 * @throws Exception If reflective TaxTree construction or serialization fails */
	private static void writeTinyTree(final Path tree) throws Exception{
		final TaxNode[] nodes=new TaxNode[13];
		nodes[TaxTree.LIFE_ID]=new TaxNode(TaxTree.LIFE_ID, TaxTree.LIFE_ID,
				TaxTree.LIFE, TaxTree.LIFE_E, "Life");
		nodes[TaxTree.BACTERIA_ID]=new TaxNode(TaxTree.BACTERIA_ID, TaxTree.LIFE_ID,
				TaxTree.DOMAIN, TaxTree.DOMAIN_E, "Bacteria");
		nodes[10]=new TaxNode(10, 11, TaxTree.SPECIES, TaxTree.SPECIES_E, "Fixture species");
		nodes[11]=new TaxNode(11, 12, TaxTree.GENUS, TaxTree.GENUS_E, "Fixture genus");
		nodes[12]=new TaxNode(12, TaxTree.BACTERIA_ID, TaxTree.FAMILY, TaxTree.FAMILY_E,
				"Fixture family");
		final Constructor<TaxTree> constructor=TaxTree.class.getDeclaredConstructor(TaxNode[].class,
				IntHashMap.class, int.class, boolean.class, boolean.class, boolean.class, int.class);
		constructor.setAccessible(true);
		final TaxTree taxTree=constructor.newInstance(nodes, null, 0, true, true, false, 0);
		TaxTree.writeTaxTree(taxTree, tree.toString(), true);
	}

	/** Reads the first data row by header name from a tiny TSV.
	 * @param table TSV containing a header and at least one data row
	 * @return First data row keyed by header
	 * @throws Exception If the file cannot be read */
	private static Map<String, String> readFirstRow(final Path table) throws Exception{
		final List<String> lines=Files.readAllLines(table, StandardCharsets.UTF_8);
		check(lines.size()>=2, "Expected a header and data row in "+table);
		final String[] header=lines.get(0).split("\t", -1);
		final String[] values=lines.get(1).split("\t", -1);
		final LinkedHashMap<String, String> row=new LinkedHashMap<String, String>();
		for(int i=0; i<header.length; i++){row.put(header[i], i<values.length ? values[i] : "");}
		return row;
	}

	/** Prints one successful child receipt before its temporary tables are deleted.
	 * @param testCase Evidence case that was run
	 * @param assertions Whether the child had Java assertions enabled
	 * @param result Child output snapshot */
	private static void printReceipt(final Case testCase, final boolean assertions, final Result result){
		System.out.println("child case="+testCase.receiptName+" assertions="+
				(assertions ? "ea" : "da")+" exit="+result.exit+" proposed_parent="+
				result.proposal.get("proposed_parent_taxid")+" evidence="+
				result.proposal.get("evidence_source")+" decision="+
				result.proposal.get("decision_class"));
	}

	/** Parses one baseline/fixed selector.
	 * @param value Argument text
	 * @param key Expected argument key
	 * @return True for fixed mode, false for baseline mode */
	private static boolean parseSelector(final String value, final String key){
		if((key+"=fixed").equals(value)){return true;}
		if((key+"=baseline").equals(value)){return false;}
		throw new AssertionError("Expected "+key+"=baseline or "+key+"=fixed, not "+value);
	}

	/** Test expectations stay fatal with assertions disabled in the parent JVM.
	 * @param condition Required test predicate
	 * @param message AssertionError message if the predicate fails */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Waits for one child JVM with a finite timeout and verifies termination.
	 * @param child Child process
	 * @param log Combined stdout/stderr destination for diagnostics
	 * @return Child exit code
	 * @throws Exception If interrupted while waiting */
	private static int waitForChild(final Process child, final Path log) throws Exception{
		try{
			if(child.waitFor(30, TimeUnit.SECONDS)){return child.exitValue();}
			terminateChild(child, log);
			throw new AssertionError("Timed out running child: "+log);
		}catch(InterruptedException e){
			try{
				terminateChild(child, log);
			}finally{
				Thread.currentThread().interrupt();
			}
			throw e;
		}
	}

	/** Forcibly terminates a child JVM and verifies it exits.
	 * @param child Child process
	 * @param log Combined stdout/stderr destination for diagnostics
	 * @throws InterruptedException If interrupted while joining the child */
	private static void terminateChild(final Process child, final Path log) throws InterruptedException{
		if(child.isAlive()){child.destroyForcibly();}
		check(child.waitFor(5, TimeUnit.SECONDS), "Timed out waiting for child termination: "+log);
	}

	/** Deletes a successful fixture tree.
	 * @param root Root directory to delete
	 * @throws Exception If the tree cannot be walked or deleted */
	private static void deleteRecursively(final Path root) throws Exception{
		try(Stream<Path> paths=Files.walk(root)){
			final Object[] sorted=paths.sorted(Comparator.reverseOrder()).toArray();
			for(Object path : sorted){Files.deleteIfExists((Path)path);}
		}
	}

	/** Builds a deterministic direct FakeTreeLeafTaxonomy child JVM command.
	 * @param assertions Whether to enable assertions
	 * @param args Tool arguments for FakeTreeLeafTaxonomy
	 * @return Fully expanded child JVM command */
	private static List<String> javaCommand(final boolean assertions, final String... args){
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(assertions ? "-ea" : "-da");
		command.add("-Xmx256m");
		command.add("-Xms256m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-XX:+ExitOnOutOfMemoryError");
		addVectorModule(command);
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add("tax.FakeTreeLeafTaxonomy");
		command.addAll(Arrays.asList(args));
		return command;
	}

	/** Creates one SSU evidence row with explicit taxID and ANI.
	 * @param tid Hit taxID in the fixture tree
	 * @param ani Literal alignment identity
	 * @return One tab-delimited evidence row */
	private static String ssuRow(final int tid, final String ani){
		return "9000001\t"+tid+"\tfixture-hit\tProteobacteria\t"+ani+"\t0.900\t100\n";
	}

	/** Adds the Vector incubator module when this JDK exposes it.
	 * @param command Child JVM command to extend */
	private static void addVectorModule(final ArrayList<String> command){
		try{
			final Class<?> moduleLayer=Class.forName("java.lang.ModuleLayer");
			final Object boot=moduleLayer.getMethod("boot").invoke(null);
			final Object optional=moduleLayer.getMethod("findModule", String.class)
					.invoke(boot, "jdk.incubator.vector");
			final boolean present=((Boolean)optional.getClass()
					.getMethod("isPresent").invoke(optional)).booleanValue();
			if(present){
				command.add("--add-modules");
				command.add("jdk.incubator.vector");
			}
		}catch(ClassNotFoundException e){
			return;
		}catch(Exception e){
			throw new RuntimeException("Could not probe jdk.incubator.vector", e);
		}
	}

	/** Finite invalid-evidence scenarios for the one-row FakeTree child fixture. */
	private enum Case{
		/** Invalid QuickClade evidence with no SSU row. */
		INVALID_QC("invalid-qc", "invalid-qc-", "",
				"9000001\t0\tbad\tProteobacteria\tphylum\tphylum:0.99\t99.0\n"),
		/** Zero SSU taxID with a usable QuickClade row. */
		ZERO_SSU("zero-ssu", "invalid-ssu-",
				"9000001\t0\tbad\tProteobacteria\t0.900\t0.900\t100\n",
				"9000001\t2\tBacteria\tProteobacteria\tdomain\tdomain:0.99\t99.0\n"),
		/** Positive SSU taxID absent from the tree with a usable QuickClade row. */
		POSITIVE_UNKNOWN_SSU("positive-unknown-ssu", "positive-unknown-ssu-",
				"9000001\t"+UNKNOWN_TID+"\tunknown-positive\tProteobacteria\t0.900\t0.900\t100\n",
				USABLE_QC_ROW),
		/** High-ANI species hit still selects its genus parent in both versions. */
		SPECIES_HIT("species-hit", ssuRow(10, "0.990"), USABLE_QC_ROW, 11, 11, "findssu"),
		/** Genus is the precise-parent boundary and retains SSU preference. */
		GENUS_SSU("genus-ssu", ssuRow(11, "0.990"), USABLE_QC_ROW, 11, 11, "findssu"),
		/** Family-level SSU yields even to broader but usable domain-level QC. */
		FAMILY_SSU("family-ssu", ssuRow(12, "0.990"), USABLE_QC_ROW, 12, 2, "quickclade"),
		/** A species hit promoted to family by its ANI is also broad evidence. */
		PROMOTED_SPECIES("promoted-species", ssuRow(10, "0.960"), USABLE_QC_ROW, 12, 2, "quickclade"),
		/** Domain SSU and QC agree on the ID, but the selected evidence changes. */
		DOMAIN_SSU("domain-ssu", ssuRow(2, "0.990"), USABLE_QC_ROW, 2, 2, "quickclade"),
		/** Life is broader than genus even though it has no tracked rank index. */
		LIFE_SSU("life-ssu", ssuRow(1, "0.990"), USABLE_QC_ROW, 1, 2, "quickclade"),
		/** Broad SSU remains usable when no QC row exists. */
		NO_QC("no-qc", ssuRow(12, "0.990"), "", 12, 12, "findssu"),
		/** A zero QC identifier cannot displace a broad SSU parent. */
		INVALID_QC_WITH_SSU("invalid-qc-with-ssu", ssuRow(12, "0.990"),
				"9000001\t0\tinvalid\tProteobacteria\tdomain\tdomain:0.99\t99.0\n", 12, 12, "findssu"),
		/** An unresolved positive QC identifier also cannot displace broad SSU. */
		MISSING_QC_WITH_SSU("missing-qc-with-ssu", ssuRow(12, "0.990"),
				"9000001\t3\tmissing\tProteobacteria\tdomain\tdomain:0.99\t99.0\n", 12, 12, "findssu");

		/** Records a precedence case with literal expected parent IDs.
		 * @param name Case and directory prefix
		 * @param ssu SSU evidence row
		 * @param qc QuickClade evidence row
		 * @param baseline Expected selected parent before the policy change
		 * @param fixed Expected selected parent after the policy change
		 * @param fixedEvidence_ Literal expected evidence source after the change */
		Case(final String name, final String ssu, final String qc,
				final int baseline, final int fixed, final String fixedEvidence_){
			receiptName=name;
			dirPrefix=name+"-";
			ssuRow=ssu;
			qcRow=qc;
			baselineParent=baseline;
			fixedParent=fixed;
			fixedEvidence=fixedEvidence_;
		}

		/** Records one invalid-evidence case.
		 * @param receiptName_ Case name for emitted receipts
		 * @param dirPrefix_ Prefix for the temporary child directory
		 * @param ssuRow_ FindSSU data row, or empty
		 * @param qcRow_ QuickClade data row */
		Case(final String receiptName_, final String dirPrefix_, final String ssuRow_,
				final String qcRow_){
			receiptName=receiptName_;
			dirPrefix=dirPrefix_;
			ssuRow=ssuRow_;
			qcRow=qcRow_;
			baselineParent=fixedParent=0;
			fixedEvidence=null;
		}

		/** Case name for emitted receipts. */
		final String receiptName;
		/** Prefix for the temporary child directory. */
		final String dirPrefix;
		/** FindSSU and QuickClade data rows. */
		final String ssuRow, qcRow;
		/** Literal parent expectations; zero denotes the older dedicated checks. */
		final int baselineParent, fixedParent;
		/** Literal evidence-source expectation for a precedence case. */
		final String fixedEvidence;
	}

	/** Captured output snapshot from one child JVM. */
	private static final class Result{
		/** Records one child outcome.
		 * @param exit_ Process exit code
		 * @param log_ Combined stdout/stderr log path
		 * @param logText_ Combined stdout/stderr log contents
		 * @param proposals_ Proposal TSV text, or empty if absent
		 * @param nodes_ Node TSV text, or empty if absent
		 * @param proposal_ First proposal row keyed by header, or empty if absent
		 * @param node_ First node row keyed by header, or empty if absent
		 * @param inputsPreserved_ True if the leaf, SSU and QuickClade inputs stayed byte-identical */
		Result(final int exit_, final Path log_, final String logText_, final String proposals_,
				final String nodes_, final Map<String, String> proposal_, final Map<String, String> node_,
				final boolean inputsPreserved_){
			exit=exit_;
			log=log_;
			logText=logText_;
			proposals=proposals_;
			nodes=nodes_;
			proposal=proposal_;
			node=node_;
			inputsPreserved=inputsPreserved_;
		}
		/** Child process exit code. */
		final int exit;
		/** Combined stdout/stderr log path. */
		final Path log;
		/** Combined child log, proposal TSV text, and node TSV text. */
		final String logText, proposals, nodes;
		/** First proposal row and node row keyed by their TSV headers. */
		final Map<String, String> proposal, node;
		/** True when the leaf, SSU and QuickClade inputs stayed byte-identical and the tree survived nonempty. */
		final boolean inputsPreserved;
	}
}
