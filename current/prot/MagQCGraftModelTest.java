package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

/** Descriptor parsing/dispatch fixture; does not claim serialized-network prediction parity.
 * @author Brian Bushnell, Yoimiya */
public final class MagQCGraftModelTest {

	public static void main(String[] args) throws Exception{
		final Path root=Files.createTempDirectory("magqc-graft-descriptor-");
		final String text="#magqc_frozen_graft_v1\nbase\tbase.bbnet.gz\nbase_sha80\t00000000000000000000\n"+
			"private\tprivate.bbnet.gz\nprivate_sha80\t11111111111111111111\n"+
			"heads\theads.bbnet.gz\nheads_sha80\t22222222222222222222\n"+
			"upstream_layer\t4\nshared_nodes\t80\nprivate_per_head\t16\nhidden_activation\ttanh\n"+
			"output_contract\t{}\nprediction_contract\t{}\n";
		final Path descriptor=root.resolve("model.graft");
		Files.write(descriptor, text.getBytes(StandardCharsets.UTF_8));
		final MagQCGraftModel.Descriptor parsed=MagQCGraftModel.describe(descriptor.toString());
		check(parsed!=null && parsed.upstreamLayer==4 && parsed.shared==80 && parsed.privatePerHead==16,
			"Descriptor must preserve declared bypass geometry");
		check(parsed.baseFile.equals(root.resolve("base.bbnet.gz").toString()), "Components must resolve next to the descriptor");
		final Path ordinary=root.resolve("ordinary.bbnet"), compressed=root.resolve("ordinary.bbnet.gz");
		final byte[] model="##bbnet\n#version 1\n".getBytes(StandardCharsets.US_ASCII);
		Files.write(ordinary, model);
		try(GZIPOutputStream output=new GZIPOutputStream(Files.newOutputStream(compressed))){output.write(model);}
		check(MagQCGraftModel.describe(ordinary.toString())==null, "Plain native models must retain their ordinary loader");
		check(MagQCGraftModel.describe(compressed.toString())==null, "Gzip native models must retain their ordinary loader");
		for(String invalid:new String[]{text+"shared_nodes\t80\n", text.replace("base.bbnet.gz", "../base.bbnet.gz"),
			text.replace("private_per_head\t16", "private_per_head\t0"), text.replace("shared_nodes\t80\n", "")}){
			Files.write(descriptor, invalid.getBytes(StandardCharsets.UTF_8));
			boolean failed=false;
			try{MagQCGraftModel.describe(descriptor.toString());}catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Invalid descriptor geometry/path/fields must fail before component loading");
		}
		System.out.println("MagQCGraftModelTest PASS: descriptor schema, sibling paths, duplicate rejection and ordinary plain/gzip dispatch; no inference claim");
	}

	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
