package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import parse.LineParser1;
import shared.Shared;
import stream.bam.BgzfOutputStream;

/** Exact topology and malformed-input fixtures for the experimental dense reader. @author Collei */
public final class HbmDenseTextLoaderTest {

	public static void main(final String[] args) throws Exception{
		final Path root=Files.createTempDirectory("hbm-dense-test-");
		final byte[] pivot={0, 1, 21};
		final AAGraph graph=new AAGraph(pivot.clone(), 0);
		graph.ref[0].count[0]=graph.ref[0].weight[0]=Integer.MAX_VALUE;
		graph.ref[0].countSum=graph.ref[0].weightSum=Integer.MAX_VALUE;
		graph.ref[0].insEdge=insertion(1, 3);
		graph.ref[0].insEdge.insEdge=insertion(1, 7);
		graph.ref[1].count[2]=graph.ref[1].weight[2]=3;
		graph.ref[1].countSum=graph.ref[1].weightSum=4;
		// Crucial: a zero-count DEL may still carry its own insertion chain.
		graph.del[0].insEdge=insertion(1, 9);
		graph.del[1].countSum=graph.del[1].weightSum=4;
		graph.del[1].insEdge=insertion(2, 7);
		graph.del[1].insEdge.insEdge=insertion(2, 1);
		final List<String> ids=Arrays.asList("fam%0;A", "family1");
		final ArrayList<HbmBundleBuilder.FamilyInput> families=new ArrayList<HbmBundleBuilder.FamilyInput>();
		for(final String id : ids){families.add(new HbmBundleBuilder.FamilyInput(id, pivot, graph));}
		final byte[][] provenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][];
		for(int i=0; i<provenance.length; i++){provenance[i]=HbmBundleFormat.sha256(new byte[]{(byte)(i+1)}, 0, 1);}
		final Path binary=root.resolve("fixture.mqhb"), fasta=root.resolve("ref.faa"), text=root.resolve("dense.txt"), compressed=root.resolve("dense.txt.gz");
		HbmBundleBuilder.build(binary, families, provenance);
		Files.write(fasta, ">fam%0;A\nARX\n>family1\nARX\n".getBytes(StandardCharsets.UTF_8));
		HbmTextDump.run(binary, null, fasta.toString(), text, HbmTextDump.Encoding.DENSE, HbmTextDump.CoordMode.EXPLICIT, 1, false);
		HbmTextDump.run(binary, null, fasta.toString(), compressed, HbmTextDump.Encoding.DENSE, HbmTextDump.CoordMode.EXPLICIT, 1, true);
		final String pin=DigestSuffix.file(binary.toString());
		final HbmBundleLoader.Loaded expected=HbmBundleLoader.load(binary, ids, id->pivot, provenance);
		for(final int threads : new int[]{1, 4}){
			Shared.setThreads(threads);
			HbmDenseTextLoader.load(text.toString(), ids, id->pivot, pin, threads).assertStructuralEquivalent(expected);
			HbmDenseTextLoader.load(compressed.toString(), ids, id->pivot, pin, threads).assertStructuralEquivalent(expected);
		}
		verified(root, binary, fasta, ids, pivot, provenance, expected);
		final String good=new String(Files.readAllBytes(text), StandardCharsets.UTF_8);
		final String[] bad={
			good.replace("#sum_source\tcandidate", "#sum_source\toriginal"),
			good.replace("#min_count_emitted\t1", "#min_count_emitted\t2"),
			good.replace("#coordinate_mode\texplicit", "#coordinate_mode\timplicit"),
			good.substring(0, good.length()-2), good+"trailing\n",
			good.replace("r\t0\td", "r\t1\td"),
			good.replace("i\t0\t1\td", "i\t0\t2\td"),
			good.replace("d\t1\t4\n", "d\t1\t4\nd\t1\t4\n"),
			good.replace("r\t2\td\t1", "r\t2\td\t2"),
			good.replace("f\t1\tfamily1", "f\t0\tfamily1"),
			good.replace("\tARX\n", "\tARY\n"),
			good.replace("fam%250%3BA", "fam%250%3BZ")
		};
		for(int i=0; i<bad.length; i++){
			if(bad[i].equals(good)){throw new AssertionError("Fixture mutation did not change input: "+i);}
			final Path broken=root.resolve("broken"+i+".txt");
			Files.write(broken, bad[i].getBytes(StandardCharsets.UTF_8));
			for(final int threads : new int[]{1, 4}){
				boolean rejected=false;
				try{HbmDenseTextLoader.load(broken.toString(), ids, id->pivot, pin, threads);}
				catch(IllegalArgumentException | java.io.IOException e){rejected=true;}
				if(!rejected){throw new AssertionError("Accepted malformed fixture "+i+" t="+threads);}
			}
		}
		for(final String invalid : new String[]{"", "00", "/", "p", "200000", "oooooo"}){
			final LineParser1 lp=new LineParser1('\t');
			lp.set(invalid.getBytes(StandardCharsets.US_ASCII));
			boolean rejected=false;
			try{HbmDenseTextLoader.integer(lp, 0);}catch(IllegalArgumentException e){rejected=true;}
			if(!rejected){throw new AssertionError("Accepted invalid integer "+invalid);}
		}
		System.out.println("HbmDenseTextLoaderTest PASS: native topology/counts, zero-DEL insertions, BGZF, t1/t4, malformed input, A48 bounds");
		System.out.println("fixture_directory\t"+root);
	}

	/** Uses the actual packer, production dispatch and externally supplied pins. */
	private static void verified(final Path root, final Path binary, final Path fasta,
			final List<String> ids, final byte[] pivot, final byte[][] provenance,
			final HbmBundleLoader.Loaded expected) throws Exception{
		final String[] keys={"format_spec", "loader_source", "loader_class", "aagraph",
			"aagraphnode", "aagraphscorer", "glocalaminolinear", "blosum62", "builder_source",
			"builder_class", "roster", "consensus_ref", "source_corpus", "cluster_membership", "member_policy", "member_manifest"};
		final StringBuilder metadata=new StringBuilder();
		for(int i=0; i<keys.length; i++){metadata.append(keys[i]).append("_sha80\t").append(DigestSuffix.fromDigest(provenance[i])).append('\n');}
		final Path manifest=root.resolve("provenance.tsv"), packed=root.resolve("verified.hbmt"), bgzf=root.resolve("verified.hbmt.gz");
		Files.write(manifest, metadata.toString().getBytes(StandardCharsets.US_ASCII));
		HbmDenseTextPacker.pack(binary, fasta.toString(), manifest.toString(), packed);
		final byte[] content=Files.readAllBytes(packed);
		try(BgzfOutputStream out=new BgzfOutputStream(Files.newOutputStream(bgzf), 9)){out.write(content); out.writeEOF();}
		for(final int threads : new int[]{1, 4}){
			Shared.setThreads(threads);
			HbmBundleLoader.load(packed, ids, id->pivot, provenance).assertStructuralEquivalent(expected);
			HbmBundleLoader.load(bgzf, ids, id->pivot, Arrays.copyOf(provenance, 8)).assertStructuralEquivalent(expected);
			for(final int minCount : new int[]{2, 4, 8, 10, Integer.MAX_VALUE}){
				final HbmBundleLoader.Loaded nativeFiltered=HbmBundleLoader.load(binary, ids, id->pivot, provenance, minCount);
				HbmBundleLoader.load(packed, ids, id->pivot, provenance, minCount).assertStructuralEquivalent(nativeFiltered);
				HbmBundleLoader.load(bgzf, ids, id->pivot, provenance, minCount).assertStructuralEquivalent(nativeFiltered);
			}
		}
		final String good=new String(content, StandardCharsets.UTF_8);
		final String[] corrupt={
			good.replace("i\t0\t0\td\t3\t3", "i\t0\t0\td\t4\t4"), // Sum stays valid; checksum must catch it.
			good.replace("#row_f\tf rank", "#row_f\tf RANK"), // Nonsemantic header corruption must also be detected.
			good.replace(HbmDenseTextBundle.CONTRACT, "unrecognized_contract"),
			good.substring(0, good.lastIndexOf("z\t")),
			good+"trailing\n"
		};
		for(int i=0; i<corrupt.length; i++){
			if(corrupt[i].equals(good)){throw new AssertionError("V2 fixture mutation missed: "+i);}
			final Path broken=root.resolve("corrupt"+i+".hbmt");
			Files.write(broken, corrupt[i].getBytes(StandardCharsets.UTF_8));
			for(final int threads : new int[]{1, 4}){
				Shared.setThreads(threads);
				rejectVerified(broken, ids, pivot, provenance, 1);
			}
		}
		final byte[][] wrong=provenance.clone(); wrong[0]=provenance[0].clone(); wrong[0][31]^=1;
		rejectVerified(packed, ids, pivot, wrong, 1);
		rejectVerified(packed, ids, pivot, provenance, 0);
		final Path crlf=root.resolve("crlf.hbmt");
		Files.write(crlf, good.replace("\n", "\r\n").getBytes(StandardCharsets.UTF_8));
		HbmBundleLoader.load(crlf, ids, id->pivot, provenance).assertStructuralEquivalent(expected);
		System.out.println("VERIFIED_DENSE_V2_FIXTURES_PASS native parity, t1/t4, BGZF, semantic pins, family/root checksums, CRLF, minCount2/4/8/10/MAX");
	}
	private static void rejectVerified(final Path file, final List<String> ids, final byte[] pivot,
			final byte[][] provenance, final int minCount) throws Exception{
		boolean rejected=false;
		try{HbmBundleLoader.load(file, ids, id->pivot, provenance, minCount);}
		catch(IllegalArgumentException | java.io.IOException e){rejected=true;}
		if(!rejected){throw new AssertionError("Accepted invalid verified text: "+file);}
	}
	private static AAGraphNode insertion(final int position, final int count){
		final AAGraphNode node=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, position);
		node.count[0]=node.weight[0]=node.countSum=node.weightSum=count;
		return node;
	}
	private HbmDenseTextLoaderTest(){}
}
