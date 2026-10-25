package test.tax;

import java.io.File;
import java.lang.reflect.Constructor;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import structures.IntHashMap;
import tax.TaxNode;
import tax.TaxTree;

/** Captures invalid QuickClade parent handling for leaf-only FakeTree nodes.
 * @author Jean
 */
public final class FakeTreeLeafTaxonomyInvalidQcTest{

	/** Runs the captured baseline or the corrected invalid-QC fallback contract. */
	public static void main(String[] args) throws Exception{
		check(args.length==2, "Expected <BBTools dir> and expect=baseline or expect=fixed");
		final boolean fixed=args[1].equals("expect=fixed");
		check(fixed || args[1].equals("expect=baseline"), "Unknown expectation: "+args[1]);
		final Path code=Paths.get(args[0]).toAbsolutePath();
		final Path root=Files.createTempDirectory("faketree-invalid-qc-");

		final Result qcEa=run(code, root, true, false);
		final Result ssuEa=run(code, root, true, true);
		if(fixed){
			assertFallback(qcEa, "invalid-QC -ea");
			assertValidQc(ssuEa, "invalid-SSU -ea");
		}else{
			assertInvalidParentFailure(qcEa, "invalid-QC -ea");
			assertInvalidParentFailure(ssuEa, "invalid-SSU -ea");
		}

		final Result qcDa=run(code, root, false, false);
		final Result ssuDa=run(code, root, false, true);
		if(fixed){
			assertFallback(qcDa, "invalid-QC -da");
			assertValidQc(ssuDa, "invalid-SSU -da");
		}else{
			assertBlankParent(qcDa, "invalid-QC -da", "quickclade");
			assertBlankParent(ssuDa, "invalid-SSU -da", "findssu");
		}
		check(qcEa.inputsPreserved && ssuEa.inputsPreserved && qcDa.inputsPreserved && ssuDa.inputsPreserved,
				"Taxonomy assignment must not mutate input tables.");
		System.out.println("FAKETREE_INVALID_QC_TEST_PASS expect="+(fixed ? "fixed" : "baseline"));
	}

	/** Verifies -ea still exposes the captured blank-parent baseline. */
	private static void assertInvalidParentFailure(final Result result, final String name){
		check(result.exit!=0, "Baseline "+name+" should catch a blank invalid parent: "+result.log);
		check(result.logText.contains("No proposed parent for fake_taxid 9000001"),
				"Baseline "+name+" must explain the invalid parent assertion: "+result.log);
	}

	/** Verifies -da exposes the captured blank-parent baseline. */
	private static void assertBlankParent(final Result result, final String name, final String source){
		check(result.exit==0, "Baseline "+name+" should expose the emitted blank parent: "+result.log);
		check("".equals(result.proposal.get("proposed_parent_taxid")),
				"Baseline "+name+" must emit the missing parent before the repair: "+result.proposals);
		check(source.equals(result.proposal.get("evidence_source")),
				"Baseline "+name+" should reveal that "+source+" was selected: "+result.proposals);
	}

	/** Verifies invalid QuickClade taxids are treated like absent usable evidence. */
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

	/** Verifies invalid SSU taxids are skipped so usable QuickClade evidence can win. */
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

	/** Runs the real wrapper on a one-row fixture with invalid QC or SSU evidence. */
	private static Result run(final Path code, final Path root, final boolean assertions,
			final boolean invalidSsu) throws Exception{
		final Path dir=Files.createDirectory(root.resolve((invalidSsu ? "invalid-ssu-" : "invalid-qc-")+
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
				(invalidSsu ? "9000001\t0\tbad\tProteobacteria\t0.900\t0.900\t100\n" : ""))
				.getBytes(StandardCharsets.UTF_8);
		final byte[] qcBytes=("fake_taxid\tref_taxid\tref_name\tqc_hit_phylum\tsketch_lca\tconf_level\tANI\n"+
				(invalidSsu ? "9000001\t2\tBacteria\tProteobacteria\tdomain\tdomain:0.99\t99.0\n" :
				"9000001\t0\tbad\tProteobacteria\tphylum\tphylum:0.99\t99.0\n"))
				.getBytes(StandardCharsets.UTF_8);
		Files.write(leaf, leafBytes);
		Files.write(ssu, ssuBytes);
		Files.write(qc, qcBytes);
		writeTinyTree(tree);

		final List<String> command=new ArrayList<String>(Arrays.asList("bash",
				code.resolve("faketreeleaves.sh").toString(), "in="+leaf, "ssu="+ssu,
				"quickclade="+qc, "out="+proposals, "outnodes="+nodes,
				"tree="+tree, "t=1", "usepigz=f", "useunpigz=f",
				"-Xmx256m", "-Xms256m", assertions ? "-ea" : "-da"));
		final Process child=makeChild(command, log).start();
		child.getOutputStream().close();
		final int exit=child.waitFor();
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

	/** Writes the smallest real TaxTree this invalid-taxid fixture can use. */
	private static void writeTinyTree(final Path tree) throws Exception{
		final TaxNode[] nodes=new TaxNode[TaxTree.BACTERIA_ID+1];
		nodes[TaxTree.LIFE_ID]=new TaxNode(TaxTree.LIFE_ID, TaxTree.LIFE_ID,
				TaxTree.LIFE, TaxTree.LIFE_E, "Life");
		nodes[TaxTree.BACTERIA_ID]=new TaxNode(TaxTree.BACTERIA_ID, TaxTree.LIFE_ID,
				TaxTree.DOMAIN, TaxTree.DOMAIN_E, "Bacteria");
		final Constructor<TaxTree> constructor=TaxTree.class.getDeclaredConstructor(TaxNode[].class,
				IntHashMap.class, int.class, boolean.class, boolean.class, boolean.class, int.class);
		constructor.setAccessible(true);
		final TaxTree taxTree=constructor.newInstance(nodes, null, 0, true, true, false, 0);
		TaxTree.writeTaxTree(taxTree, tree.toString(), true);
	}

	/** Reads the first data row by header name from a tiny TSV. */
	private static Map<String, String> readFirstRow(final Path table) throws Exception{
		final List<String> lines=Files.readAllLines(table, StandardCharsets.UTF_8);
		check(lines.size()>=2, "Expected a header and data row in "+table);
		final String[] header=lines.get(0).split("\t", -1);
		final String[] values=lines.get(1).split("\t", -1);
		final LinkedHashMap<String, String> row=new LinkedHashMap<String, String>();
		for(int i=0; i<header.length; i++){row.put(header[i], i<values.length ? values[i] : "");}
		return row;
	}

	/** Test expectations stay fatal with assertions disabled in the parent JVM. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Pins wrapper subprocesses to this JVM's Java and a deterministic tiny runtime. */
	private static ProcessBuilder makeChild(final List<String> command, final Path log){
		final ProcessBuilder pb=new ProcessBuilder(command).redirectErrorStream(true)
				.redirectOutput(log.toFile());
		final String javaHome=System.getProperty("java.home");
		final String oldPath=pb.environment().get("PATH");
		pb.environment().put("JAVA_HOME", javaHome);
		pb.environment().put("PATH", Paths.get(javaHome, "bin")+File.pathSeparator+
				(oldPath==null ? "" : oldPath));
		return pb;
	}

	private static final class Result{
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
		final int exit;
		final Path log;
		final String logText, proposals, nodes;
		final Map<String, String> proposal, node;
		final boolean inputsPreserved;
	}
}
