package prot;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.zip.GZIPOutputStream;
import java.util.zip.GZIPInputStream;

import shared.Shared;
import structures.ByteBuilder;

/** Native graph and byte-image oracles for the lossless production transport.
 * @author Yoimiya
 */
public final class HbmA48CodecTest {

	/** Uses a tiny real bundle, not a codec-generated expected answer. */
	public static void main(String[] args) throws IOException{
		final Path root=Files.createTempDirectory("hbm-a48-test-");
		final byte[] pivot={0, 1, 2};
		final AAGraph graph=new AAGraph(pivot, 0);
		graph.ref[0].count[1]=graph.ref[0].weight[1]=65536;
		graph.ref[0].countSum+=65536; graph.ref[0].weightSum+=65536;
		graph.del[1].countSum=graph.del[1].weightSum=63;
		final AAGraphNode insertion=new AAGraphNode(Blosum62.X_CODE, AAGraphNode.INS, 1);
		insertion.count[2]=insertion.weight[2]=64;
		insertion.countSum=insertion.weightSum=64;
		graph.ref[0].insEdge=insertion;
		final byte[][] provenance=new byte[HbmBundleFormat.PROVENANCE_COUNT][];
		for(int i=0; i<provenance.length; i++){
			provenance[i]=HbmBundleFormat.sha256(new byte[]{(byte)i}, 0, 1);
		}
		final Path raw=root.resolve("native.mqhb"), encoded=root.resolve("encoded.mqhb.bgz"), decoded=root.resolve("decoded.mqhb");
		HbmBundleBuilder.build(raw, Arrays.asList(new HbmBundleBuilder.FamilyInput("family0", pivot, graph)), provenance);
		final byte[] expected=Files.readAllBytes(raw);
		HbmBundleBuilder.encodeA48(raw, encoded);
		HbmA48Codec.decode(encoded, decoded);
		check(Arrays.equals(expected, Files.readAllBytes(decoded)), "Decoded native bytes differ");
		checkReadModes(root, encoded, expected);
		checkNativeDigest(root, encoded);
		HbmBundleLoader.load(raw, Arrays.asList("family0"), id->pivot, provenance).assertStructuralMatch(0, graph);
		HbmBundleLoader.load(encoded, Arrays.asList("family0"), id->pivot, provenance).assertStructuralMatch(0, graph);
		check(!HbmA48Codec.isEncoded(raw) && HbmA48Codec.isEncoded(encoded), "Native/A48 dispatch differs");
		for(int value:new int[]{0, 1, 63, 64, 65535, 65536, Integer.MAX_VALUE}){
			final ByteBuilder token=new ByteBuilder(); token.appendA48(value).tab();
			check(HbmA48Codec.integer(new DataInputStream(new ByteArrayInputStream(token.toBytes())))==value,
				"Integer boundary changed: "+value);
		}
		for(String token:new String[]{"\t", "00\t", "p\t", "/\t", "oooooo\t", "1000000\t", "1"}){
			boolean rejected=false;
			try{HbmA48Codec.integer(new DataInputStream(new ByteArrayInputStream(token.getBytes(StandardCharsets.US_ASCII))));}
			catch(IOException expectedFailure){rejected=true;}
			check(rejected, "Malformed integer accepted: "+token);
		}
		boolean rejected=false;
		try{HbmA48Codec.decode(encoded, raw);}catch(IOException expectedFailure){rejected=true;}
		check(rejected && Arrays.equals(expected, Files.readAllBytes(raw)), "Existing output was changed");
		final Path malformed=root.resolve("wrong_encoding.bgz");
		try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(malformed))){
			out.write("MQHA\t1\n#encoding\tu16\n".getBytes(StandardCharsets.US_ASCII));
		}
		reject(malformed, root.resolve("invalid.mqhb"));
		final byte[] compressed=Files.readAllBytes(encoded);
		final Path truncated=root.resolve("truncated.bgz");
		Files.write(truncated, Arrays.copyOf(compressed, compressed.length/2));
		reject(truncated, root.resolve("truncated.mqhb"));
		System.out.println("HbmA48CodecTest PASS exact native bytes, graph parity, integer boundaries, legacy dispatch, malformed/truncated rejection, fresh-output preservation");
		System.out.println("fixture_directory\t"+root);
	}

	/** Both gzip transports and scalar/parallel input must preserve native bytes. */
	private static void checkReadModes(final Path root, final Path encoded, final byte[] expected) throws IOException{
		final Path plain=root.resolve("plain-gzip.mqhb.bgz");
		try(InputStream in=new GZIPInputStream(Files.newInputStream(encoded));
				GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(plain))){
			final byte[] buffer=new byte[8192];
			for(int count=in.read(buffer); count>=0; count=in.read(buffer)){
				if(count>0){out.write(buffer, 0, count);}
			}
		}
		final int oldThreads=Shared.threads();
		try{
			for(int threads:new int[]{1, 4}){
				Shared.setThreads(threads);
				for(Path input:new Path[]{encoded, plain}){
					final Path output=root.resolve(input.getFileName()+".t"+threads+".decoded");
					HbmA48Codec.decode(input, output);
					check(Arrays.equals(expected, Files.readAllBytes(output)),
						"Gzip transport changed native bytes: "+input+", threads="+threads);
				}
			}
		}finally{Shared.setThreads(oldThreads);}
	}

	/** Valid gzip with a corrupt native digest must still be rejected after buffered writes. */
	private static void checkNativeDigest(final Path root, final Path encoded) throws IOException{
		final ByteArrayOutputStream decoded=new ByteArrayOutputStream();
		try(InputStream in=new GZIPInputStream(Files.newInputStream(encoded))){
			final byte[] buffer=new byte[8192];
			for(int count=in.read(buffer); count>=0; count=in.read(buffer)){
				if(count>0){decoded.write(buffer, 0, count);}
			}
		}
		final byte[] bytes=decoded.toByteArray();
		check(bytes.length>=36,"Fixture must include the native digest and sentinel");
		bytes[bytes.length-36]^=1;
		final Path bad=root.resolve("bad-native-digest.bgz");
		try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(bad))){out.write(bytes);}
		reject(bad,root.resolve("bad-native-digest.mqhb"));
	}

	/** Failed decodes must not leave a publishable native output. */
	private static void reject(final Path input, final Path output) throws IOException{
		boolean rejected=false;
		try{HbmA48Codec.decode(input, output);}catch(IOException expected){rejected=true;}
		check(rejected && !Files.exists(output), "Malformed encoded model left an accepted output");
	}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private HbmA48CodecTest(){}
}
