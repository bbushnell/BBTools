package prot;

import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import stream.bam.BgzfInputStream;
import structures.ByteBuilder;

/** Native topology fixture for {@link HbmTextDump}. */
public final class HbmTextDumpTest {

	public static void main(final String[] args) throws Exception{
		final Path root=Files.createTempDirectory("hbm-text-dump-test-");
		final AAGraph graph=fixture();
		final byte[] pivot=graph.pivot.clone();
		final ArrayList<HbmBundleBuilder.FamilyInput> originalFamilies=new ArrayList<HbmBundleBuilder.FamilyInput>();
		originalFamilies.add(new HbmBundleBuilder.FamilyInput("fam%0;A", pivot, graph));
		final byte[][] provenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][];
		for(int i=0; i<provenance.length; i++){
			provenance[i]=HbmBundleFormat.sha256(new byte[]{(byte)(i+1)}, 0, 1);
		}
		final Path original=root.resolve("original.mqhb");
		HbmBundleBuilder.build(original, originalFamilies, provenance);
		HbmRareStateProbe.prune(graph, 200);
		final Path candidate=root.resolve("candidate.mqhb");
		HbmBundleBuilder.build(candidate, originalFamilies, provenance);
		final ByteBuilder fasta=new ByteBuilder();
		fasta.append(">fam%0;A\n");
		for(final byte b : pivot){fasta.append(AminoAcid.numberToAcid[b]);}
		fasta.nl();
		final Path ref=root.resolve("consensus.faa");
		Files.write(ref, fasta.toBytes());

		final Path sparse=root.resolve("sparse.txt");
		HbmTextDump.run(candidate, original, ref.toString(), sparse, HbmTextDump.Encoding.SPARSE,
			HbmTextDump.CoordMode.EXPLICIT, 1, false);
		final String text=new String(Files.readAllBytes(sparse), StandardCharsets.UTF_8);
		require(text.contains("#format\thbm_text_v1\n"), "format header missing");
		require(text.contains("#start_end_counts\tunavailable\n"), "start/end reservation missing");
		require(text.contains("f\t0\tfam%250%3BA\t3\t"), "encoded family id missing");
		require(text.contains("r\t0\tr\t?[\tA1\tR?V\tX3\n"), "REF original sum or sparse A48 counts wrong");
		require(text.contains("d\t0\t:\n"), "surviving DEL row missing original A48 sum");
		require(text.contains("i\t0\t0\tr\t:\tA:\n"), "REF insertion A48 row missing");
		require(text.contains("di\t0\t0\tr\t:\tA:\n"), "DEL insertion A48 row missing");
		require(!containsLinePrefix(text, "i\t0\t1\t"), "pruned insertion suffix survived");

		final Path sparseBgz=root.resolve("sparse.bgz");
		HbmTextDump.run(candidate, original, ref.toString(), sparseBgz, HbmTextDump.Encoding.SPARSE,
			HbmTextDump.CoordMode.EXPLICIT, 1, true);
		require(text.equals(readBgzf(sparseBgz)), "sparse BGZF differs from plain sparse text");

		final Path sparseMin2=root.resolve("sparse_min2.txt");
		HbmTextDump.run(candidate, original, ref.toString(), sparseMin2, HbmTextDump.Encoding.SPARSE,
			HbmTextDump.CoordMode.EXPLICIT, 2, false);
		final String min2Text=new String(Files.readAllBytes(sparseMin2), StandardCharsets.UTF_8);
		require(min2Text.contains("#min_count_emitted\t2\n"), "mincount header missing");
		require(min2Text.contains("r\t0\tr\t?[\tR?V\tX3\n"), "mincount=2 sparse row wrong");

		final Path dense=root.resolve("dense.txt");
		HbmTextDump.run(candidate, original, ref.toString(), dense, HbmTextDump.Encoding.DENSE,
			HbmTextDump.CoordMode.EXPLICIT, 1, false);
		final String denseText=new String(Files.readAllBytes(dense), StandardCharsets.UTF_8);
		require(denseText.contains(denseRef0()), "dense REF row missing native 22-slot order");

		final Path denseBgz=root.resolve("dense.bgz");
		HbmTextDump.run(candidate, original, ref.toString(), denseBgz, HbmTextDump.Encoding.DENSE,
			HbmTextDump.CoordMode.EXPLICIT, 1, true);
		require(denseText.equals(readBgzf(denseBgz)), "dense BGZF differs from plain dense text");

		final Path adaptiveText=root.resolve("adaptive.txt");
		HbmTextDump.run(candidate, candidate, ref.toString(), adaptiveText, HbmTextDump.Encoding.ADAPTIVE,
			HbmTextDump.CoordMode.EXPLICIT, 1, false);
		final Path adaptive=root.resolve("adaptive.bgz");
		HbmTextDump.run(candidate, candidate, ref.toString(), adaptive, HbmTextDump.Encoding.ADAPTIVE,
			HbmTextDump.CoordMode.EXPLICIT, 1, true);
		require(new String(Files.readAllBytes(adaptiveText), StandardCharsets.UTF_8).equals(readBgzf(adaptive)),
			"adaptive BGZF differs from plain adaptive text");

		final Path implicit=root.resolve("implicit.txt");
		HbmTextDump.run(candidate, original, ref.toString(), implicit, HbmTextDump.Encoding.SPARSE,
			HbmTextDump.CoordMode.IMPLICIT, 1, false);
		final String implicitText=new String(Files.readAllBytes(implicit), StandardCharsets.UTF_8);
		require(implicitText.contains("#coordinate_mode\timplicit\n"), "implicit coordinate header missing");
		require(implicitText.contains("r\tr\t?[\tA1\tR?V\tX3\n"), "implicit REF sparse row wrong");
		require(implicitText.contains("d\t:\n"), "implicit DEL row wrong");
		require(implicitText.contains("i\tr\t:\tA:\n"), "implicit REF insertion row wrong");
		require(implicitText.contains("di\tr\t:\tA:\n"), "implicit DEL insertion row wrong");
		require(!implicitText.contains("r\t0\tr\t"), "explicit REF anchor survived implicit mode");
		require(!implicitText.contains("i\t0\t0\t"), "explicit insertion coordinates survived implicit mode");
		System.out.println("HbmTextDumpTest PASS sparse/dense/adaptive BGZF, mincount, explicit/implicit coords, REF/DEL insertion topology, encoded IDs, original totals");
		System.out.println("fixture_directory\t"+root);
	}

	private static AAGraph fixture(){
		final byte[] pivot={0, 1, 2};
		final AAGraph graph=new AAGraph(pivot, 0);
		set(graph.ref[0], 0, 1); add(graph.ref[0], 1, 998); add(graph.ref[0], 2, 1); add(graph.ref[0], 21, 3);
		set(graph.ref[1], 1, 996); add(graph.ref[1], 2, 4);
		set(graph.ref[2], 2, 1000);
		graph.del[0].countSum=graph.del[0].weightSum=10;
		graph.del[1].countSum=graph.del[1].weightSum=1;
		graph.ref[0].insEdge=ins(1, 10, 0);
		graph.ref[0].insEdge.insEdge=ins(1, 5, 0);
		graph.del[0].insEdge=ins(1, 10, 0);
		return graph;
	}

	private static AAGraphNode ins(final int position, final int count, final int residue){
		final AAGraphNode node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, position);
		set(node, residue, count);
		return node;
	}

	private static void set(final AAGraphNode node, final int residue, final int count){
		Arrays.fill(node.count, 0);
		Arrays.fill(node.weight, 0);
		node.countSum=node.weightSum=0;
		add(node, residue, count);
	}

	private static void add(final AAGraphNode node, final int residue, final int count){
		node.count[residue]+=count;
		node.weight[residue]+=count;
		node.countSum+=count;
		node.weightSum+=count;
	}

	private static String denseRef0(){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("r\t0\td\t?[");
		bb.tab().appendA48(1).tab().appendA48(998);
		for(int i=2; i<Blosum62.X_CODE; i++){bb.tab().appendA48(0);}
		bb.tab().appendA48(3).nl();
		return new String(bb.toBytes(), StandardCharsets.UTF_8);
	}

	private static String readBgzf(final Path path) throws Exception{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		try(InputStream in=new BgzfInputStream(new BufferedInputStream(Files.newInputStream(path)))){
			final byte[] buffer=new byte[8192];
			for(int len=in.read(buffer); len>=0; len=in.read(buffer)){
				if(len>0){out.write(buffer, 0, len);}
			}
		}
		return new String(out.toByteArray(), StandardCharsets.UTF_8);
	}

	private static boolean containsLinePrefix(final String text, final String prefix){
		return text.startsWith(prefix) || text.contains("\n"+prefix);
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private HbmTextDumpTest(){}
}
