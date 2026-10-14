package prot;

import java.io.IOException;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * Adds ProkCC native-harness resource bindings to a bare explicit-input final net.
 * The adapter preserves exported weights and normalization headers verbatim.
 * @author Yelan
 */
public final class MagQCFinalNetBinder {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args) throws Exception{
		final HashMap<String,String> o=new HashMap<String,String>();
		for(String arg:args){
			final int eq=arg.indexOf('=');
			require(eq>0, "Expected key=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT);
			require(KEYS.contains(key) && o.put(key, arg.substring(eq+1))==null, "Unknown or duplicate option: "+key);
		}
		if(parse.Parse.parseBoolean(o.getOrDefault("selftest", "f"))){selftest();return;}
		final String in=required(o, "in"), out=required(o, "out"), config=required(o, "config");
		final HashMap<String,String> options=MagQCAssemblyBatch.parseOptions(new String[]{"config="+config});
		bind(in, out, options);
		System.out.println("MAGQC_FINAL_NET_BINDER_PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Binding            ----------------*/
	/*--------------------------------------------------------------*/

	private static void bind(String in, String out, HashMap<String,String> options) throws IOException{
		final String[] contracts=scan(in);
		final String[] names=parseExplicitNames(contracts[0]);
		if(contracts[1]!=null){MagQCPredictionContract.parse(contracts[1], names);}
		final ByteBuilder tags=new ByteBuilder();
		tags.append("##magqc_harness_contract prepared_raw_v1\n");
		tags.append("##magqc_layout ").append(MagQCNetworkHarness.expectedLayout(options)).nl();
		tags.append("##magqc_encoding raw\n");
		for(String key:PIN_KEYS){
			tags.append("##magqc_").append(key).append(' ').append(required(options, key)).nl();
		}
		tags.append("##magqc_output_mode explicit\n");
		tags.append("##magqc_output_names ").append(String.join(",", names)).nl();
		tags.append("##magqc_input_normalization embedded\n");
		copyWithTags(in, out, tags);
	}

	/** Prediction metadata is exported before resource binding and must survive verbatim. */
	private static String[] scan(String in) throws IOException{
		final ByteFile bf=ByteFile.makeByteFile(in, true);
		String contract=null;
		String prediction=null;
		long lineNumber=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lineNumber++;
				if(lineNumber==1 && !equals(line, BBNET)){
					throw new IllegalArgumentException("Final net is missing ##bbnet header");
				}
				if(startsWith(line, PREDICTION_CONTRACT)){
					if(prediction!=null){throw new IllegalArgumentException("Duplicate magqc_prediction_contract");}
					prediction=new String(line, PREDICTION_CONTRACT.length, line.length-PREDICTION_CONTRACT.length, StandardCharsets.UTF_8);
				}else if(startsWith(line, MAGQC)){
					throw new IllegalArgumentException("Final net is already bound: "+new String(line, StandardCharsets.UTF_8));
				}
				if(startsWith(line, OUTPUT_CONTRACT)){
					if(contract!=null){throw new IllegalArgumentException("Duplicate output_contract");}
					contract=new String(line, OUTPUT_CONTRACT.length, line.length-OUTPUT_CONTRACT.length, StandardCharsets.UTF_8);
				}
			}
		}finally{
			if(bf.close()){throw new IllegalStateException("Error reading "+in);}
		}
		if(lineNumber<1){throw new IllegalArgumentException("Final net is missing ##bbnet header");}
		if(contract==null){throw new IllegalArgumentException("Missing line: ##output_contract");}
		return new String[]{contract, prediction};
	}

	private static void copyWithTags(String in, String out, ByteBuilder tags) throws IOException{
		final Path target=Paths.get(out).toAbsolutePath();
		final Path dir=target.getParent();
		Files.createDirectories(dir);
		final Path tmp=dir.resolve(tempName(target.getFileName().toString()));
		ByteStreamWriter writer=null;
		boolean success=false;
		try{
			writer=new ByteStreamWriter(tmp.toString(), true, false, true);
			writer.start();
			final ByteFile bf=ByteFile.makeByteFile(in, true);
			try{
				for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
					writer.println(line);
					if(equals(line, BBNET)){writer.print(tags);}
				}
			}finally{
				if(bf.close()){throw new IllegalStateException("Error reading "+in);}
			}
			if(writer.poisonAndWait()){throw new IllegalStateException("Error writing "+tmp);}
			writer=null;
			move(tmp, target);
			success=true;
		}finally{
			if(writer!=null){writer.poisonAndWait();}
			if(!success){Files.deleteIfExists(tmp);}
		}
	}

	private static void move(Path tmp, Path target) throws IOException{
		try{
			Files.move(tmp, target, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE);
		}catch(AtomicMoveNotSupportedException e){
			Files.move(tmp, target, StandardCopyOption.REPLACE_EXISTING);
		}
	}

	private static String tempName(String name){
		final String suffix=name.endsWith(".gz") ? ".gz" : "";
		final String base=suffix.isEmpty() ? name : name.substring(0, name.length()-suffix.length());
		return base+".tmp."+System.nanoTime()+suffix;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	private static String[] parseExplicitNames(String contract){
		require(contract.indexOf("\"mode\":\"explicit\"")>=0 &&
			contract.indexOf("\"schema\":\"explicit_targets_v1\"")>=0,
			"Final net output_contract must be explicit_targets_v1");
		final ArrayList<String> names=new ArrayList<String>();
		final HashSet<String> seen=new HashSet<String>();
		for(int start=contract.indexOf("\"name\":\""); start>=0; start=contract.indexOf("\"name\":\"", start)){
			start+=8;
			final int end=contract.indexOf('"', start);
			require(end>=0, "Malformed output_contract name");
			final String name=contract.substring(start, end);
			require(name.length()>0 && name.indexOf('\\')<0 && seen.add(name), "Malformed or duplicate output name: "+name);
			names.add(name);
		}
		require(!names.isEmpty(), "Empty explicit output_contract");
		return names.toArray(new String[names.size()]);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Self Test           ----------------*/
	/*--------------------------------------------------------------*/

	private static void selftest() throws IOException{
		final Path root=Files.createTempDirectory("magqc-final-net-binder-");
		final Path net=root.resolve("bare.bbnet"), config=root.resolve("release.config"), out=root.resolve("bound.bbnet");
		Files.write(net, ("##bbnet\n#version 1\n#dims 3 2\n##normalization explicit-input\n##mean_sd_folded false\n"+
			"##output_contract {\"mode\":\"explicit\",\"schema\":\"explicit_targets_v1\",\"targets\":[{\"name\":\"a\"},{\"name\":\"b\"}]}\n").getBytes(StandardCharsets.UTF_8));
		Files.write(config, ("bundle=b\nbundlesha80=11111111111111111111\nfamilylist=f\nfamilylistsha80=22222222222222222222\n"+
			"subnetmanifest=s\nsubnetmanifestsha80=33333333333333333333\nexpectedcopytable=e\nexpectedcopytablesha80=44444444444444444444\n"+
			"subnetpopulations=p\nsubnetpopulationssha80=55555555555555555555\nrequiredness=r\nrequirednesssha80=66666666666666666666\n"+
			"sizemarkers=z\nsizemarkerssha80=77777777777777777777\n").getBytes(StandardCharsets.UTF_8));
		bind(net.toString(), out.toString(), MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}));
		final String bound=new String(Files.readAllBytes(out), StandardCharsets.UTF_8);
		require(bound.contains("##magqc_output_names a,b\n"), "Output names were not bound");
		require(bound.contains("##magqc_bundlesha80 11111111111111111111\n"), "Resource pins were not bound");
		require(bound.contains("##normalization explicit-input\n"), "Normalization was not preserved");
		boolean rejected=false;
		try{bind(out.toString(), root.resolve("bad.bbnet").toString(), MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}));}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Already-bound net was not rejected");
		rejected=false;
		final Path duplicate=root.resolve("duplicate.bbnet");
		Files.write(duplicate, ("##bbnet\n##output_contract {\"mode\":\"explicit\",\"schema\":\"explicit_targets_v1\",\"targets\":[{\"name\":\"a\"}]}\n"+
			"##output_contract {\"mode\":\"explicit\",\"schema\":\"explicit_targets_v1\",\"targets\":[{\"name\":\"b\"}]}\n").getBytes(StandardCharsets.UTF_8));
		try{bind(duplicate.toString(), root.resolve("duplicate-bound.bbnet").toString(), MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}));}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Duplicate output_contract was not rejected");
		rejected=false;
		final Path missing=root.resolve("missing.bbnet");
		Files.write(missing, "##bbnet\n#version 1\n".getBytes(StandardCharsets.UTF_8));
		try{bind(missing.toString(), root.resolve("missing-bound.bbnet").toString(), MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}));}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Missing output_contract was not rejected");
		final String[] declared={"gene_completeness", "gene_contamination", "magnified_completeness", "magnified_contamination"};
		final String metadata=MagQCPredictionContractTest.metadata(declared, 6, 10);
		final StringBuilder declaration=new StringBuilder("{\"mode\":\"explicit\",\"schema\":\"explicit_targets_v1\",\"heads\":[");
		for(int i=0; i<declared.length; i++){
			if(i>0){declaration.append(',');}
			declaration.append("{\"name\":\"").append(declared[i]).append("\"}");
		}
		declaration.append("]}");
		final String bare="##bbnet\n#dims 3 4\n##output_contract "+declaration+
			"\n##magqc_prediction_contract "+metadata+"\n";
		final Path interpreted=root.resolve("interpreted.bbnet"), interpretedOut=root.resolve("interpreted-bound.bbnet");
		Files.write(interpreted, bare.getBytes(StandardCharsets.UTF_8));
		bind(interpreted.toString(), interpretedOut.toString(), MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}));
		require(new String(Files.readAllBytes(interpretedOut), StandardCharsets.UTF_8).contains("##magqc_prediction_contract "+metadata+"\n"),
			"Resource binding must preserve the exported prediction metadata byte-for-byte");
		for(String invalid:new String[]{bare.replace("magqc_prediction_contract_v1", "unknown_schema"),
			bare+"##magqc_prediction_contract "+metadata+"\n"}){
			Files.write(interpreted, invalid.getBytes(StandardCharsets.UTF_8));
			rejected=false;
			try{bind(interpreted.toString(), root.resolve("invalid-prediction.bbnet").toString(),
				MagQCAssemblyBatch.parseOptions(new String[]{"config="+config}));}
			catch(IllegalArgumentException expected){rejected=true;}
			require(rejected && !Files.exists(root.resolve("invalid-prediction.bbnet")),
				"Invalid prediction metadata must fail before publishing a bound net");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private static String required(HashMap<String,String> o, String key){
		final String value=o.get(key);
		require(value!=null && !value.isEmpty(), "Required: "+key);
		return value;
	}
	private static void require(boolean good, String why){if(!good){throw new IllegalArgumentException(why);}}

	private static boolean equals(byte[] line, byte[] prefix){
		if(line.length!=prefix.length){return false;}
		for(int i=0; i<prefix.length; i++){if(line[i]!=prefix[i]){return false;}}
		return true;
	}

	private static boolean startsWith(byte[] line, byte[] prefix){
		if(line.length<prefix.length){return false;}
		for(int i=0; i<prefix.length; i++){if(line[i]!=prefix[i]){return false;}}
		return true;
	}

	private static final byte[] BBNET="##bbnet".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] MAGQC="##magqc_".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] OUTPUT_CONTRACT="##output_contract ".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] PREDICTION_CONTRACT="##magqc_prediction_contract ".getBytes(StandardCharsets.US_ASCII);
	private static final String[] PIN_KEYS={"bundlesha80", "familylistsha80", "subnetmanifestsha80",
		"expectedcopytablesha80", "subnetpopulationssha80"};
	private static final HashSet<String> KEYS=new HashSet<String>(java.util.Arrays.asList("in", "out", "config", "selftest"));
}
