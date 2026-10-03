package prot;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import ml.CellNet;
import ml.CellNetParser;

/**
 * Packs real CellNet fixtures with all four observation types and checks their
 * independent release binding, inference and malformed semantic definitions.
 * IDs deliberately differ from type names, so ID-based dispatch cannot pass.
 * @author Yoimiya
 */
public final class MagQCNetBundleSpecialTest {

	/** Runs through the product pack/load path; fixture files live under one temporary root. */
	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("Special bundle fixture accepts no arguments");}
		checkHeaderBounds();
		final Path root=Files.createTempDirectory("magqc-special-bundle-");
		final String[] types={"famset","ncrna","rrna","trna_anticodon"};
		final int[] counts={2,5,4,65};
		final StringBuilder release=new StringBuilder(),aggregate=new StringBuilder();
		write(root.resolve("families.tsv"),"#rank\trep_id\tocc_total\n0\tf0\t1\n1\tf1\t1\n");
		write(root.resolve("subset.txt"),"0\n1\n");
		write(root.resolve("taxonomy.tsv"),"1\tP\n");
		for(int i=0; i<types.length; i++){
			final String id="entry_"+i,subset=i==0 ? "subset.txt" : "-";
			final int inputs=counts[i]+19;
			MagQCVectorMakerSixFeatureTest.writeFourOutputNet(root.resolve(id+".bbnet").toFile(),inputs,100+i);
			release.append(id).append('\t').append(types[i]).append('\t').append(subset).append("\toutputs=4\n");
			aggregate.append(id).append('\t').append(counts[i]).append('\t').append(inputs).append('\t')
				.append(subset).append("\t-\t").append(id).append(".bbnet\t").append(types[i]).append('\n');
		}
		final Path releaseFile=root.resolve("release.tsv"),aggregateFile=root.resolve("aggregate.tsv");
		write(releaseFile,release.toString());
		write(aggregateFile,aggregate.toString());
		final Map<String,String> config=new LinkedHashMap<String,String>();
		config.put("subnetmanifest",releaseFile.toString());
		config.put("subnetmanifestsha80",sha80(Files.readAllBytes(releaseFile)));
		config.put("aggmanifest",aggregateFile.toString());
		config.put("familylist",root.resolve("families.tsv").toString());
		config.put("taxpgm",root.resolve("taxonomy.tsv").toString());
		config.put("netroot",root.toString());
		final Path output=root.resolve("mixed.bbnets");
		config.put("out",output.toString());
		MagQCNetBundle.pack(config);
		final MagQCNetBundle bundle=MagQCNetBundle.loadMultiOutput(output);
		bundle.requireReleaseManifest(releaseFile,config.get("subnetmanifestsha80"));
		check(bundle.size()==4,"Every declared type must survive the round trip");
		for(int i=0; i<types.length; i++){
			final MagQCNetBundle.Subnet subnet=bundle.subnet(i);
			check(types[i].equals(subnet.type) && counts[i]==subnet.numObs,"Type and observation count at "+i);
			check((i==0)==(subnet.familyRanks!=null),"Only famset carries ranks");
			final CellNet loose=CellNetParser.load(root.resolve("entry_"+i+".bbnet").toString(),false);
			final float[] input=new float[subnet.expectedInputs];
			for(int j=0; j<input.length; j++){input[j]=(j%7)*0.25f;}
			loose.applyInput(input);
			loose.feedForward();
			final float[] actual=subnet.scoreAll(input);
			for(int j=0; j<actual.length; j++){
				check(Float.floatToIntBits(actual[j])==Float.floatToIntBits(loose.getOutput(j)),
					"Loose/bundled inference differs at subnet="+i+" output="+j);
			}
		}
		// Removing explicit type from a new special row must not reinterpret its ID.
		config.put("out",root.resolve("rejected.bbnets").toString());
		write(aggregateFile,aggregate.toString().replace("entry_2.bbnet\trrna", "entry_2.bbnet"));
		expectPackFailure(config,"type mismatch");
		write(aggregateFile,aggregate.toString().replace("entry_2\t4\t23", "entry_2\t5\t23"));
		expectPackFailure(config,"observation count");
		write(aggregateFile,aggregate.toString().replace("entry_3.bbnet\ttrna_anticodon", "entry_3.bbnet\tunknown_kind"));
		expectPackFailure(config,"Unknown MAG-QC observation type");
		write(aggregateFile,aggregate.toString().replace("entry_2\t4\t23\t-", "entry_2\t4\t23\tsubset.txt"));
		expectPackFailure(config,"subset must be -");
		write(aggregateFile,aggregate.toString());
		checkDefinitionTamper(root,output);
		// Legacy six-column famset rows continue to work in the same new mixed roster.
		write(aggregateFile,aggregate.toString().replace("entry_0.bbnet\tfamset", "entry_0.bbnet"));
		config.put("out",root.resolve("legacy-famset.bbnets").toString());
		MagQCNetBundle.pack(config);
		check(MagQCNetBundle.loadMultiOutput(root.resolve("legacy-famset.bbnets")).size()==4,"Legacy famset row accepted");
		System.out.println("MAGQC_SPECIAL_BUNDLE PASS: four types, exact inference, explicit dispatch, malformed contracts");
	}

	/** Differential against the prior full-text header reader, including EOF/CRLF boundaries. */
	private static void checkHeaderBounds() throws Exception{
		final java.lang.reflect.Method method=MagQCNetBundle.class.getDeclaredMethod("headerLines", byte[].class);
		method.setAccessible(true);
		final String[] prefixes={"", "\n", "\r\n", "\r", "#dense\n", "#dense\r\n#sparse\r\n\r\n",
			"#dims 2 4\n\n#seed 17\n", "##name forêt\n#density1 0.5\n"};
		final String[] bodies={"", "\n", "\r\n", "C0\n#seed 999\n", "W1 0.5\r\n#sparse\r\n", "\r\r\n"};
		for(String prefix:prefixes){
			for(String body:bodies){
				final byte[] bytes=(prefix+body).getBytes(StandardCharsets.UTF_8);
				final String[] all=new String(bytes,StandardCharsets.UTF_8).split("\\r?\\n");
				int end=0;
				while(end<all.length && (all[end].isEmpty() || all[end].startsWith("#"))){end++;}
				final String[] actual=(String[])method.invoke(null,(Object)bytes);
				check(Arrays.equals(Arrays.copyOf(all,end),actual),"Header boundary differs for "+Arrays.toString(all));
			}
		}
		System.out.println("MAGQC_HEADER_BOUNDARIES_PASS cases="+(prefixes.length*bodies.length));
	}

	/** A width-compatible but false observation description must fail before hashing. */
	private static void checkDefinitionTamper(Path root,Path output) throws Exception{
		final ByteArrayOutputStream plain=new ByteArrayOutputStream();
		try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(output))){
			final byte[] buffer=new byte[8192];
			for(int n; (n=in.read(buffer))!=-1;){plain.write(buffer,0,n);}
		}
		final String text=new String(plain.toByteArray(),StandardCharsets.UTF_8);
		final String changed=text.replace("ncrna_obs_definition=r16,r23,r5,rother", "ncrna_obs_definition=wrong");
		check(!text.equals(changed),"rRNA definition must be serialized explicitly");
		final Path tampered=root.resolve("bad-definition.bbnets");
		try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(tampered))){
			out.write(changed.getBytes(StandardCharsets.UTF_8));
		}
		try{MagQCNetBundle.loadMultiOutput(tampered); throw new AssertionError("False observation description accepted");}
		catch(IOException expected){check(expected.getMessage().contains("observation definition"),"Wrong tamper rejection");}
	}

	/** A bad aggregate row cannot leave a published bundle. */
	private static void expectPackFailure(Map<String,String> config,String reason) throws Exception{
		try{MagQCNetBundle.pack(config); throw new AssertionError("Malformed aggregate accepted: "+reason);}
		catch(IOException expected){check(expected.getMessage().contains(reason),"Unexpected rejection: "+expected.getMessage());}
		check(!Files.exists(java.nio.file.Paths.get(config.get("out"))),"Rejected pack published output");
	}

	/** Writes only deterministic small fixture metadata. */
	private static void write(Path path,String text) throws IOException{
		assert(text!=null) : "Test fixtures must specify their complete intended metadata";
		Files.write(path,text.getBytes(StandardCharsets.UTF_8));
	}

	/** Computes the release pin without emitting a full digest. */
	private static String sha80(byte[] bytes) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(bytes);
		final StringBuilder out=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){
			out.append(Character.forDigit((digest[i]>>>4)&15,16)).append(Character.forDigit(digest[i]&15,16));
		}
		assert(out.length()==20) : "The independent release binding requires exactly 80 digest bits";
		return out.toString();
	}

	/** Enforces test invariants under both -ea and -da. */
	private static void check(boolean condition,String message){
		if(!condition){throw new AssertionError(message);}
	}
}
