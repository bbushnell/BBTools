package ml;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ByteArrayOutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/** Analytic fixtures for ModifyNN growth, pruning, no-op copies, and RNG parity. */
public final class ModifyNNTest {

	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("modifynn_").toFile();
		try{
			checkNoop(dir);
			checkDuplicateNoop(dir);
			checkUnsafeNewWeight(dir);
			checkDenseGrowth(dir);
			checkZeroAbsentAccounting(dir);
			checkZeroToEpsilon(dir);
			checkCompressedGrowth(dir);
			checkCompressedNoop(dir);
			checkNoopEncodingMatrix(dir);
			checkNoopTransportRejection(dir);
			checkReadonlyInputDirectory(dir);
			checkMetadataArchive(dir);
			checkPartitionRejection(dir);
			checkMalformedA48(dir);
			checkHeterogeneousPrune(dir);
			checkSparsePruneAndGrowth(dir);
			checkActualRng(dir);
			checkInputNormalizationZeroPreservation(dir);
			checkOutputOnlyPreservation(dir);
			checkSigmoidHiddenWidening(dir);
		}finally{
			delete(dir);
		}
		System.out.println("ModifyNNTest PASS");
	}

	private static void checkNoop(File dir) throws Exception{
		final File in=new File(dir, "noop.bbnet"), out=new File(dir, "noop_out.bbnet");
		final File report=new File(dir, "noop_report.tsv");
		final CellNet net=denseLinear(new int[]{2, 2, 1});
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "report="+report});
		final byte[] a=Files.readAllBytes(in.toPath()), b=Files.readAllBytes(out.toPath());
		check(Arrays.equals(a, b), "no-op must copy original bytes");
		final String text=new String(Files.readAllBytes(report.toPath()));
		check(text.contains("retained_edges\t6\n"), "no-op report must count original active edges");
	}

	private static void checkDuplicateNoop(File dir) throws Exception{
		final File in=new File(dir, "duplicate.bbnet");
		final byte[] original=denseLinear(new int[]{2, 2, 1}).toBytes().toString().getBytes();
		Files.write(in.toPath(), original);
		boolean failed=false;
		try{
			ModifyNN.main(new String[]{"in="+in, "out="+in});
		}catch(RuntimeException e){
			failed=true;
		}
		check(failed, "duplicate in/out must be rejected");
		check(Arrays.equals(original, Files.readAllBytes(in.toPath())), "duplicate no-op must not delete input");
	}

	private static void checkUnsafeNewWeight(File dir) throws Exception{
		final File in=new File(dir, "unsafe.bbnet"), out=new File(dir, "unsafe_out.bbnet");
		Files.write(in.toPath(), denseLinear(new int[]{2, 1}).toBytes().toString().getBytes());
		boolean failed=false;
		try{
			ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,1", "newweight=1e-8"});
		}catch(IllegalArgumentException e){
			failed=true;
		}
		check(failed, "FP16-unsafe newweight must be rejected");
	}

	private static void checkDenseGrowth(File dir) throws Exception{
		final File in=new File(dir, "dense.bbnet"), out=new File(dir, "dense_out.bbnet");
		final File report=new File(dir, "dense_report.tsv");
		Files.write(in.toPath(), denseLinear(new int[]{2, 2, 1}, true).toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,3,2", "seed=7", "report="+report});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		check(grown.dims[0]==3 && grown.dims[1]==3 && grown.dims[2]==2, "dense dims not grown");
		check(grown.net[1][0].weights[0]==0 && grown.net[1][0].weights[1]==13f, "old dense weights changed");
		check(grown.net[1][0].weights[2]!=0 && grown.net[1][2].weights[0]!=0, "new dense edges not initialized");
		check(grown.net[1][2].bias==0 && grown.net[2][1].bias==0, "new biases must be zero");
		check(grown.net[2][0].weights[0]==0 && grown.net[2][0].weights[1]==19f, "old output row changed");
		final String text=new String(Files.readAllBytes(report.toPath()));
		check(text.contains("retained_edges\t3\n"), "dense zero slots must be absent from retained edge counts");
		check(text.matches("(?s).*parent_sha80\t[0-9a-f]{20}\n.*"), "parent sha80 missing from report");
		check("append_only".equals(grown.getTag("index_preservation"))==false, "report key leaked to net tags");
		check(grown.getTag("parent_sha80")!=null && grown.getTag("parent_sha80").matches("[0-9a-f]{20}"),
			"parent sha80 missing from child provenance");
	}

	private static void checkZeroAbsentAccounting(File dir) throws Exception{
		final File noopIn=new File(dir, "zero_absent_noop.bbnet");
		final File noopOut=new File(dir, "zero_absent_noop_out.bbnet");
		final File noopReport=new File(dir, "zero_absent_noop_report.tsv");
		Files.write(noopIn.toPath(), denseLinear(new int[]{2, 2, 1}, true).toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+noopIn, "out="+noopOut, "report="+noopReport});
		check(Arrays.equals(Files.readAllBytes(noopIn.toPath()), Files.readAllBytes(noopOut.toPath())),
			"zero-absent no-op must still copy original bytes");
		check(new String(Files.readAllBytes(noopReport.toPath())).contains("retained_edges\t3\n"),
			"zero-absent no-op must count only nonzero dense edges");

		final File legacyIn=new File(dir, "zero_legacy_noop.bbnet");
		final File legacyOut=new File(dir, "zero_legacy_noop_out.bbnet");
		final File legacyReport=new File(dir, "zero_legacy_noop_report.tsv");
		Files.write(legacyIn.toPath(), unflaggedDenseLegacy(new int[]{2, 2, 1}, true).getBytes());
		ModifyNN.main(new String[]{"in="+legacyIn, "out="+legacyOut, "report="+legacyReport});
		check(new String(Files.readAllBytes(legacyReport.toPath())).contains("retained_edges\t6\n"),
			"unflagged dense no-op must preserve legacy zero-is-active edge counts");

		final File denseIn=new File(dir, "zero_absent_dense.bbnet");
		final File denseOut=new File(dir, "zero_absent_dense_out.bbnet");
		final File denseReport=new File(dir, "zero_absent_dense_report.tsv");
		final CellNet dense=denseLinear(new int[]{2, 2, 1}, true);
		dense.net[1][0].weights[1]=0.005f;
		Files.write(denseIn.toPath(), dense.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+denseIn, "out="+denseOut, "pruneabs=0.01", "report="+denseReport});
		final CellNet densePruned=CellNetParser.load(denseOut.toString(), false);
		check(densePruned.net[1][0].inputs.length==0, "dense zero slots must be omitted from sparse rewrite");
		final String denseText=new String(Files.readAllBytes(denseReport.toPath()));
		check(denseText.contains("retained_edges\t2\n"), "dense zero slots must not count as retained");
		check(denseText.contains("deleted_edges\t1\n"), "only nonzero subthreshold edges should count as deleted");
		check(denseText.contains("layer_1_disconnected_neurons\t1\n"), "all-zero rows must count as disconnected");

		final File sparseIn=new File(dir, "zero_absent_sparse.bbnet");
		final File sparseOut=new File(dir, "zero_absent_sparse_out.bbnet");
		final File sparseReport=new File(dir, "zero_absent_sparse_report.tsv");
		CellNet.DENSE=false;
		final CellNet sparse=new CellNet(new int[]{3, 2}, 2, 1f, 0f, 1, new ArrayList<String>());
		sparse.net[1][0].function=Function.getFunction(Function.LINEAR);
		sparse.net[1][0].bias=1;
		sparse.net[1][0].inputs=new int[]{0, 2};
		sparse.net[1][0].weights=new float[]{0, 0.5f};
		sparse.net[1][0].deltas=new float[2];
		sparse.net[1][1].function=Function.getFunction(Function.LINEAR);
		sparse.net[1][1].bias=2;
		sparse.net[1][1].inputs=new int[]{1};
		sparse.net[1][1].weights=new float[]{0.125f};
		sparse.net[1][1].deltas=new float[1];
		CellNet.makeOutputSets(sparse.net);
		sparse.makeWeightMatrices();
		Files.write(sparseIn.toPath(), sparse.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+sparseIn, "out="+sparseOut, "dims=4,3", "report="+sparseReport});
		final CellNet sparseGrown=CellNetParser.load(sparseOut.toString(), false);
		check(Arrays.equals(sparseGrown.net[1][0].inputs, new int[]{0, 2, 3}),
			"explicit sparse zero must be retained on rewrite");
		final String sparseText=new String(Files.readAllBytes(sparseReport.toPath()));
		check(sparseText.contains("retained_edges\t3\n"), "sparse zero entry must count as retained");
		check(sparseText.contains("deleted_edges\t0\n"), "sparse zero entry must not count as deleted");
	}

	/** Checks deterministic zero revival, separate accounting and preservation of absent sparse edges. */
	private static void checkZeroToEpsilon(File dir) throws Exception{
		final File denseIn=new File(dir, "zero2eps_dense.bbnet");
		final File denseOut=new File(dir, "zero2eps_dense_out.bbnet");
		final File denseOut2=new File(dir, "zero2eps_dense_out2.bbnet");
		final File denseReport=new File(dir, "zero2eps_dense_report.tsv");
		Files.write(denseIn.toPath(), unflaggedDenseLegacy(new int[]{2, 2, 1}, true).getBytes());
		ModifyNN.main(new String[]{"in="+denseIn, "out="+denseOut, "seed=3", "zero2epsilon=t", "report="+denseReport});
		ModifyNN.main(new String[]{"in="+denseIn, "out="+denseOut2, "seed=3", "zero2epsilon=t"});
		check(Arrays.equals(Files.readAllBytes(denseOut.toPath()), Files.readAllBytes(denseOut2.toPath())),
			"zero2epsilon must be deterministic for a fixed seed");
		final CellNet denseGrown=CellNetParser.load(denseOut.toString(), false);
		checkEpsilon(denseGrown.net[1][0].weights[0], "first dense active zero was not changed to epsilon");
		checkEpsilon(denseGrown.net[1][1].weights[0], "second dense active zero was not changed to epsilon");
		checkEpsilon(denseGrown.net[2][0].weights[0], "output dense active zero was not changed to epsilon");
		check(denseGrown.net[1][0].weights[1]==13f && denseGrown.net[2][0].weights[1]==19f,
			"zero2epsilon changed retained dense nonzeros");
		final String denseText=new String(Files.readAllBytes(denseReport.toPath()));
		check(denseText.contains("retained_edges\t3\n"), "dense zero2epsilon retained count changed");
		check(denseText.contains("converted_edges\t3\n"), "dense active zero converted count missing");
		check(denseText.contains("deleted_edges\t0\n"), "dense active zeros must not be deleted");

		final File absentIn=new File(dir, "zero2eps_absent.bbnet");
		final File absentOut=new File(dir, "zero2eps_absent_out.bbnet");
		final File absentReport=new File(dir, "zero2eps_absent_report.tsv");
		Files.write(absentIn.toPath(), denseLinear(new int[]{2, 2, 1}, true).toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+absentIn, "out="+absentOut, "seed=3", "zero2epsilon=t", "report="+absentReport});
		final CellNet absentGrown=CellNetParser.load(absentOut.toString(), false);
		check(absentGrown.net[1][0].weights[0]==0 && absentGrown.net[2][0].weights[0]==0,
			"zero2epsilon must not convert absent dense zero-mask slots");
		check(new String(Files.readAllBytes(absentReport.toPath())).contains("converted_edges\t0\n"),
			"absent dense zero-mask slots must not count as converted");

		final File sparseIn=new File(dir, "zero2eps_sparse.bbnet");
		final File sparseOut=new File(dir, "zero2eps_sparse_out.bbnet");
		final File sparseReport=new File(dir, "zero2eps_sparse_report.tsv");
		CellNet.DENSE=false;
		final CellNet sparse=new CellNet(new int[]{3, 2}, 2, 1f, 0f, 1, new ArrayList<String>());
		sparse.net[1][0].function=Function.getFunction(Function.LINEAR);
		sparse.net[1][0].bias=1;
		sparse.net[1][0].inputs=new int[]{0, 2};
		sparse.net[1][0].weights=new float[]{0, 0.5f};
		sparse.net[1][0].deltas=new float[2];
		sparse.net[1][1].function=Function.getFunction(Function.LINEAR);
		sparse.net[1][1].bias=2;
		sparse.net[1][1].inputs=new int[]{1};
		sparse.net[1][1].weights=new float[]{0.125f};
		sparse.net[1][1].deltas=new float[1];
		CellNet.makeOutputSets(sparse.net);
		sparse.makeWeightMatrices();
		Files.write(sparseIn.toPath(), sparse.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+sparseIn, "out="+sparseOut, "dims=4,3", "seed=3",
			"zero2epsilon=t", "report="+sparseReport});
		final CellNet sparseGrown=CellNetParser.load(sparseOut.toString(), false);
		check(Arrays.equals(sparseGrown.net[1][0].inputs, new int[]{0, 2, 3}),
			"zero2epsilon must preserve explicit sparse zero while leaving absent sparse inputs absent");
		checkEpsilon(sparseGrown.net[1][0].weights[0], "explicit sparse zero was not changed to epsilon");
		check(sparseGrown.net[1][0].weights[1]==0.5f, "zero2epsilon changed retained sparse nonzero");
		final String sparseText=new String(Files.readAllBytes(sparseReport.toPath()));
		check(sparseText.contains("retained_edges\t2\n"), "sparse zero2epsilon retained count changed");
		check(sparseText.contains("converted_edges\t1\n"), "sparse explicit zero converted count missing");
		check(sparseText.contains("added_edges\t6\n"), "growth adds one column to each old sink and four edges to the new sink");
	}

	private static void checkCompressedGrowth(File dir) throws Exception{
		final File in=new File(dir, "compressed.bbnet"), out=new File(dir, "compressed_out.bbnet.gz");
		final File report=new File(dir, "compressed_report.tsv.gz");
		Files.write(in.toPath(), denseLinear(new int[]{2, 2, 1}).toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,2,1", "report="+report});
		final byte[] bytes=Files.readAllBytes(out.toPath());
		check((bytes[0]&0xFF)==0x1F && (bytes[1]&0xFF)==0x8B, "compressed output lost gzip encoding");
		final CellNet grown=CellNetParser.load(out.toString(), false);
		check(grown.dims[0]==3, "compressed grown output did not reload natively");
		check(readGzip(report).contains("added_edges\t2\n"), "compressed report lost gzip encoding");
	}

	private static void checkCompressedNoop(File dir) throws Exception{
		final File in=new File(dir, "noop.bbnet.gz"), out=new File(dir, "noop_out.bbnet.gz");
		final byte[] raw=denseLinear(new int[]{2, 2, 1}).toBytes().toString().getBytes();
		writeGzip(in, raw);
		ModifyNN.main(new String[]{"in="+in, "out="+out});
		check(Arrays.equals(Files.readAllBytes(in.toPath()), Files.readAllBytes(out.toPath())),
			"compressed no-op must copy original bytes");
		CellNetParser.load(out.toString(), false);
	}

	private static void checkNoopEncodingMatrix(File dir) throws Exception{
		final int[] bits={18, 24, 32, 32};
		final boolean[] a48={true, true, true, false};
		for(boolean sparse : new boolean[]{false, true}){
			for(boolean gzip : new boolean[]{false, true}){
				for(int i=0; i<bits.length; i++){
					final CellNet net=sparse ? sparseLinear() : denseLinear(new int[]{3, 2});
					final String name="noop_matrix_"+(sparse ? "sparse" : "dense")+"_"+
						(a48[i] ? "a48_" : "decimal_")+bits[i]+(gzip ? ".bbnet.gz" : ".bbnet");
					final File in=new File(dir, "in_"+name), out=new File(dir, "out_"+name);
					writeEncoded(in, net, a48[i], bits[i], gzip);
					ModifyNN.main(new String[]{"in="+in, "out="+out});
					check(Arrays.equals(Files.readAllBytes(in.toPath()), Files.readAllBytes(out.toPath())),
						"no-op matrix changed bytes for "+name);
					final CellNet parsed=CellNetParser.load(out.toString(), false);
					check(parsed.weightBits()==bits[i], "no-op matrix changed weightbits for "+name);
				}
			}
		}
	}

	private static void checkNoopTransportRejection(File dir) throws Exception{
		final byte[] raw=denseLinear(new int[]{2, 2, 1}).toBytes().toString().getBytes();
		final File gzipIn=new File(dir, "transport_in.bbnet.gz");
		final File plainOut=new File(dir, "transport_plain_out.bbnet");
		writeGzip(gzipIn, raw);
		boolean failed=false;
		try{
			ModifyNN.main(new String[]{"in="+gzipIn, "out="+plainOut});
		}catch(IllegalArgumentException e){
			failed=true;
		}
		check(failed, "gzip-to-plain no-op must be rejected");
		check(!plainOut.exists(), "gzip-to-plain no-op must not publish");

		final File plainIn=new File(dir, "transport_in.bbnet");
		final File gzipOut=new File(dir, "transport_gzip_out.bbnet.gz");
		Files.write(plainIn.toPath(), raw);
		failed=false;
		try{
			ModifyNN.main(new String[]{"in="+plainIn, "out="+gzipOut});
		}catch(IllegalArgumentException e){
			failed=true;
		}
		check(failed, "plain-to-gzip no-op must be rejected");
		check(!gzipOut.exists(), "plain-to-gzip no-op must not publish");
	}

	private static void checkReadonlyInputDirectory(File dir) throws Exception{
		final File inputDir=new File(dir, "readonly_input");
		final File in=new File(inputDir, "source.bbnet"), out=new File(dir, "readonly_out.bbnet");
		inputDir.mkdirs();
		Files.write(in.toPath(), denseLinear(new int[]{2, 2, 1}).toBytes().toString().getBytes());
		inputDir.setWritable(false, false);
		try{
			ModifyNN.main(new String[]{"in="+in, "out="+out});
		}finally{
			inputDir.setWritable(true, false);
		}
		check(out.exists(), "read-only input directory no-op did not publish");
	}

	private static void checkMetadataArchive(File dir) throws Exception{
		final File in=new File(dir, "metadata.bbnet"), out=new File(dir, "metadata_out.bbnet");
		final CellNet net=denseLinear(new int[]{2, 2, 1});
		net.setTag("output_contract", "{\"outputs\":[\"old\"]}");
		net.setTag("training_dims", "2 2 1");
		net.setTag("magqc_layout", "stale_layout");
		net.setTag("descriptive", "keep");
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,2,1"});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		check(grown.getTag("output_contract")==null, "stale output contract stayed active");
		check(grown.getTag("training_dims")==null, "stale training dims stayed active");
		check(grown.getTag("magqc_layout")==null, "stale MAG-QC resource binding stayed active");
		check("{\"outputs\":[\"old\"]}".equals(grown.getTag("parent_output_contract")), "parent contract not archived");
		check("2 2 1".equals(grown.getTag("parent_training_dims")), "parent training dims not archived");
		check("stale_layout".equals(grown.getTag("parent_magqc_layout")), "parent MAG-QC binding not archived");
		check("keep".equals(grown.getTag("descriptive")), "ordinary descriptive metadata not preserved");
	}

	private static void checkSparsePruneAndGrowth(File dir) throws Exception{
		CellNet.DENSE=false;
		final File in=new File(dir, "sparse.bbnet"), out=new File(dir, "sparse_out.bbnet");
		final File report=new File(dir, "sparse_report.tsv");
		final CellNet net=new CellNet(new int[]{3, 2}, 2, 1f, 0f, 1, new ArrayList<String>());
		net.net[1][0].function=Function.getFunction(Function.LINEAR);
		net.net[1][0].bias=1;
		net.net[1][0].inputs=new int[]{0, 2};
		net.net[1][0].weights=new float[]{0, 0.5f};
		net.net[1][0].deltas=new float[2];
		net.net[1][1].function=Function.getFunction(Function.LINEAR);
		net.net[1][1].bias=2;
		net.net[1][1].inputs=new int[]{1};
		net.net[1][1].weights=new float[]{0.02f};
		net.net[1][1].deltas=new float[1];
		CellNet.makeOutputSets(net.net);
		net.makeWeightMatrices();
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=4,3", "seed=9", "pruneabs=0.02", "report="+report});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		check(Arrays.equals(grown.net[1][0].inputs, new int[]{0, 2, 3}), "sparse pruning/growth for old sink failed");
		check(grown.net[1][0].weights[0]==0 && grown.net[1][0].weights[1]==0.5f &&
			grown.net[1][0].weights[2]!=0, "sparse old/new weights failed");
		check(Arrays.equals(grown.net[1][2].inputs, new int[]{0, 1, 2, 3}), "new sparse sink not fully connected");
		final String text=new String(Files.readAllBytes(report.toPath()));
		check(text.contains("retained_edges\t3\n"), "explicit sparse zero must count as retained");
		check(text.contains("deleted_edges\t0\n"), "explicit sparse zero must not count as deleted");
	}

	private static void checkPartitionRejection(File dir) throws Exception{
		final File in=new File(dir, "partition.bbnet"), out=new File(dir, "partition_out.bbnet");
		final CellNet net=denseLinear(new int[]{2, 2, 1});
		net.setTag("output_partition", "0:1");
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		boolean failed=false;
		try{
			ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,3,2"});
		}catch(IllegalArgumentException e){
			failed=true;
		}
		check(failed, "output_partition rewrites must fail loudly");
		check(!out.exists(), "partition failure must not publish output");
	}

	private static void checkHeterogeneousPrune(File dir) throws Exception{
		final File in=new File(dir, "heterogeneous.bbnet"), out=new File(dir, "heterogeneous_out.bbnet");
		final CellNet net=denseLinear(new int[]{2, 2});
		net.net[1][1].function=Function.getFunction(Function.TANH);
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "pruneabs=0.01"});
		final CellNet pruned=CellNetParser.load(out.toString(), false);
		check(pruned.net[1][0].function==net.net[1][0].function, "first old activation changed");
		check(pruned.net[1][1].function==net.net[1][1].function, "second old activation changed");
		check(pruned.net[1][0].inputs.length==2, "dense prune did not emit explicit sparse rows");
	}

	private static void checkMalformedA48(File dir) throws Exception{
		final File nan=new File(dir, "nan.bbnet"), nanOut=new File(dir, "nan_out.bbnet");
		final CellNet nanNet=denseLinear(new int[]{2, 1});
		nanNet.net[1][0].weights[0]=Float.NaN;
		writeA48(nan, nanNet);
		boolean failed=false;
		try{
			ModifyNN.main(new String[]{"in="+nan, "out="+nanOut});
		}catch(IllegalArgumentException e){
			failed=true;
		}
		check(failed, "32-bit A48 NaN must be rejected before no-op publication");
		check(!nanOut.exists(), "malformed no-op must not publish");

		final File inf=new File(dir, "inf.bbnet"), infOut=new File(dir, "inf_out.bbnet");
		final CellNet infNet=denseLinear(new int[]{2, 1});
		infNet.net[1][0].bias=Float.POSITIVE_INFINITY;
		writeA48(inf, infNet);
		failed=false;
		try{
			ModifyNN.main(new String[]{"in="+inf, "out="+infOut, "dims=3,1"});
		}catch(IllegalArgumentException e){
			failed=true;
		}
		check(failed, "32-bit A48 Inf bias must be rejected before growth");
		check(!infOut.exists(), "malformed growth must not publish");
	}

	private static void checkActualRng(File dir) throws Exception{
		final File in=new File(dir, "rng.bbnet"), out=new File(dir, "rng_out.bbnet");
		Files.write(in.toPath(), denseLinear(new int[]{1, 1, 1}).toBytes().toString().getBytes());
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,3,3", "seed=1"});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		final int[] bits=new int[]{
			Float.floatToRawIntBits(grown.net[1][0].weights[1]),
			Float.floatToRawIntBits(grown.net[1][0].weights[2]),
			Float.floatToRawIntBits(grown.net[1][1].weights[0]),
			Float.floatToRawIntBits(grown.net[1][1].weights[1]),
			Float.floatToRawIntBits(grown.net[1][1].weights[2]),
			Float.floatToRawIntBits(grown.net[1][2].weights[0]),
			Float.floatToRawIntBits(grown.net[1][2].weights[1]),
			Float.floatToRawIntBits(grown.net[1][2].weights[2])};
		final int[] expected=new int[]{
			-1174124527, -1169408077, -1172514882, 975520799,
			973337228, -1173498822, -1172383905, -1172877678};
		check(Arrays.equals(bits, expected), "actual added-edge RNG fixture changed: "+Arrays.toString(bits));
	}

	private static void checkInputNormalizationZeroPreservation(File dir) throws Exception{
		final File in=new File(dir, "norm.bbnet"), out=new File(dir, "norm_out.bbnet");
		final CellNet net=denseLinear(new int[]{2, 2, 1});
		net.setInputNormalization(new float[]{10f, 20f}, new float[]{0.5f, 2f});
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		final float oldPrediction=predict(net, new float[]{12f, 21f})[0];
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,2,1"});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		check(Arrays.equals(grown.inputMeanCopy(), new float[]{10f, 20f, 0f}), "input means were not extended");
		check(Arrays.equals(grown.inputInverseStdCopy(), new float[]{0.5f, 2f, 1f}), "inverse deviations were not extended");
		check(predict(grown, new float[]{12f, 21f, 0f})[0]==oldPrediction,
			"zero-normalized appended input changed the old prediction");
	}

	private static void checkOutputOnlyPreservation(File dir) throws Exception{
		final File in=new File(dir, "outonly.bbnet"), out=new File(dir, "outonly_out.bbnet");
		final CellNet net=denseLinear(new int[]{2, 2, 1});
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		final float oldPrediction=predict(net, new float[]{0.25f, -0.5f})[0];
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=2,2,2"});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		check(predict(grown, new float[]{0.25f, -0.5f})[0]==oldPrediction,
			"appended output row changed existing output");
	}

	private static void checkSigmoidHiddenWidening(File dir) throws Exception{
		final File in=new File(dir, "sigmoid.bbnet"), out=new File(dir, "sigmoid_out.bbnet");
		final CellNet net=denseLinear(new int[]{3, 4, 2});
		for(Cell cell : net.net[1]){cell.function=Function.getFunction(Function.SIG);}
		Files.write(in.toPath(), net.toBytes().toString().getBytes());
		final float[] oldZero=predict(net, new float[]{0, 0, 0});
		final float[] oldOut=predict(net, new float[]{0.1f, -0.2f, 0.3f});
		ModifyNN.main(new String[]{"in="+in, "out="+out, "dims=3,7,2", "seed=1", "newweight=0.001"});
		final CellNet grown=CellNetParser.load(out.toString(), false);
		final float[] newZero=predict(grown, new float[]{0, 0, 0});
		final float[] newOut=predict(grown, new float[]{0.1f, -0.2f, 0.3f});
		final float[] zeroDelta=maxDelta(oldZero, newZero), nonzeroDelta=maxDelta(oldOut, newOut);
		System.out.println("sigmoid_hidden_widening\tzero\tmaxabs="+zeroDelta[0]+
			"\tmaxrel_floor1e-6="+zeroDelta[1]);
		System.out.println("sigmoid_hidden_widening\tnonzero\tmaxabs="+nonzeroDelta[0]+
			"\tmaxrel_floor1e-6="+nonzeroDelta[1]);
		check(nonzeroDelta[0]>1e-4f && nonzeroDelta[0]<0.01f,
			"sigmoid hidden widening drift outside expected scale: "+nonzeroDelta[0]);
	}

	private static CellNet denseLinear(int[] dims){
		return denseLinear(dims, false);
	}

	private static String unflaggedDenseLegacy(int[] dims, boolean activeZero){
		return denseLinear(dims, activeZero).toBytes().toString().replaceAll("#zeroabsent true\n#liveedges [0-9]+\n", "");
	}

	private static CellNet denseLinear(int[] dims, boolean activeZero){
		CellNet.DENSE=true;
		final CellNet net=new CellNet(dims, 1, 1f, 0f, 1, new ArrayList<String>());
		float value=11;
		for(int layer=1; layer<dims.length; layer++){
			for(Cell cell:net.net[layer]){
				cell.function=Function.getFunction(Function.LINEAR);
				cell.bias=value++;
				cell.weights=new float[dims[layer-1]];
				cell.deltas=new float[cell.weights.length];
				for(int i=0; i<cell.weights.length; i++){cell.weights[i]=value++;}
				if(activeZero && cell.weights.length>0){cell.weights[0]=0;}
			}
		}
		net.makeWeightMatrices();
		return net;
	}

	private static CellNet sparseLinear(){
		CellNet.DENSE=false;
		final CellNet net=new CellNet(new int[]{3, 2}, 2, 1f, 0f, 1, new ArrayList<String>());
		net.net[1][0].function=Function.getFunction(Function.LINEAR);
		net.net[1][0].bias=1;
		net.net[1][0].inputs=new int[]{0, 2};
		net.net[1][0].weights=new float[]{0.5f, -0.25f};
		net.net[1][0].deltas=new float[2];
		net.net[1][1].function=Function.getFunction(Function.LINEAR);
		net.net[1][1].bias=2;
		net.net[1][1].inputs=new int[]{1};
		net.net[1][1].weights=new float[]{0.125f};
		net.net[1][1].deltas=new float[1];
		CellNet.makeOutputSets(net.net);
		net.makeWeightMatrices();
		return net;
	}

	private static float[] maxDelta(float[] a, float[] b){
		float maxAbs=0, maxRel=0;
		for(int i=0; i<a.length; i++){
			final float abs=Math.abs(b[i]-a[i]);
			maxAbs=Math.max(maxAbs, abs);
			maxRel=Math.max(maxRel, abs/Math.max(Math.abs(a[i]), 1e-6f));
		}
		return new float[]{maxAbs, maxRel};
	}

	private static float[] predict(CellNet net, float[] in){
		net.applyInput(in);
		net.feedForward();
		return net.getOutput();
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Requires a converted edge to be live and inside the default initialization range. */
	private static void checkEpsilon(final float weight, final String message){
		check(weight!=0 && Math.abs(weight)>=0.5e-3f && Math.abs(weight)<=1e-3f, message+": "+weight);
	}

	private static void writeA48(File file, CellNet net) throws Exception{
		writeEncoded(file, net, true, 32, false);
	}

	private static void writeEncoded(File file, CellNet net, boolean a48, int weightBits, boolean gzip) throws Exception{
		final boolean old=CellNet.codingA48Out;
		CellNet.codingA48Out=a48;
		net.setWeightBits(weightBits);
		try{
			final byte[] bytes=net.toBytes().toString().getBytes();
			if(gzip){writeGzip(file, bytes);}
			else{Files.write(file.toPath(), bytes);}
		}finally{
			CellNet.codingA48Out=old;
		}
	}

	private static void writeGzip(File file, byte[] bytes) throws Exception{
		final GZIPOutputStream out=new GZIPOutputStream(new FileOutputStream(file));
		out.write(bytes);
		out.close();
	}

	private static String readGzip(File file) throws Exception{
		final GZIPInputStream in=new GZIPInputStream(new FileInputStream(file));
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		final byte[] buffer=new byte[4096];
		for(int len; (len=in.read(buffer))>=0; ){
			out.write(buffer, 0, len);
		}
		in.close();
		return out.toString();
	}

	private static void delete(File file){
		if(file==null || !file.exists()){return;}
		if(file.isDirectory()){for(File child:file.listFiles()){delete(child);}}
		file.delete();
	}
}
