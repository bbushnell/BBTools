package prot;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import parse.Parse;
import parse.PreParser;
import structures.ByteBuilder;

/**
 * Streams prepared whole-bin statistics through the frozen subnet ensemble to input-only composite vectors.
 * Optional final-network scoring consumes the exact serialized input values, without a second transform.
 * @author Nilou
 */
public final class MagQCNetworkHarness {

	private MagQCNetworkHarness(){}

	/** Resource pins are SHA80 of file bytes; all outputs must be fresh files. */
	public static void main(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		//TODO: Probable bug - PreParser consumes outstream and can truncate a file
		//before this developer harness validates its output paths. The public client
		//uses side-effect-free Parser.parseConfig; apply the same separation here.
		for(String arg:new PreParser(args, null, false).args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected key=value argument");}
			final String key=arg.substring(0, equals).toLowerCase(java.util.Locale.ROOT);
			if(!ALLOWED.containsKey(key) || options.put(key, arg.substring(equals+1))!=null){
				throw new IllegalArgumentException("Unknown or duplicate harness argument: "+key);
			}
		}
		// Match MagQCVectorMaker deterministic=t before subnet or final-network inference.
		if(options.containsKey("deterministic") && Parse.parseBoolean(required(options, "deterministic"))){
			shared.Shared.SIMD_FMA=false; shared.Shared.SIMD_FEED_FORWARD=false;
		}
		final boolean assembly=options.containsKey("fasta");
		if(assembly==options.containsKey("in")){
			throw new IllegalArgumentException("Specify exactly one of in=prepared.tsv or fasta=assembly.fa");
		}
		if(!assembly){
			for(String key:MagQCAssemblyInput.OPTIONS){
				if(options.containsKey(key)){throw new IllegalArgumentException(key+" requires fasta= assembly input");}
			}
		}
		final String input=assembly ? null : required(options, "in");
		final String output=options.containsKey("out") ? required(options, "out") : null;
		final String scores=options.containsKey("scores") ? required(options, "scores") : null;
		if(output==null && scores==null){throw new IllegalArgumentException("Required: out= vectors or scores= destination");}
		if(output==null && options.containsKey("idsout")){throw new IllegalArgumentException("idsout= requires vector out=");}
		if(scores==null && (options.containsKey("net") || options.containsKey("netsha80") || options.containsKey("dummy"))){
			throw new IllegalArgumentException("Final-network arguments require scores=");
		}
		final String ledger=output==null ? null : options.containsKey("idsout") ? required(options, "idsout") : output+".ids.tsv";
		final String bundle=pinned(options, "bundle", "bundlesha80");
		final String family=pinned(options, "familylist", "familylistsha80");
		final String release=pinned(options, "subnetmanifest", "subnetmanifestsha80");
		final String table=pinned(options, "expectedcopytable", "expectedcopytablesha80");
		final String populations=pinned(options, "subnetpopulations", "subnetpopulationssha80");
		// Classification/calling precedes subnet initialization: no classifier can alter loaded NN state.
		final MagQCPreparedBin prepared=assembly ? MagQCAssemblyInput.build(options) : null;
		final MagQCVectorMaker vm=MagQCVectorMaker.initializePrepared(bundle, family, release,
			options.get("subnetmanifestsha80"), table, options.get("expectedcopytablesha80"),
			populations, options.get("subnetpopulationssha80"));
		final Scorer scorer=scores==null ? null : new Scorer(options, vm.preparedInputWidth());
		final long rows=write(input, prepared, output, ledger, scores, vm, scorer, options);
		System.err.println("MagQCNetworkHarness PASS: rows="+rows+" inputs="+vm.preparedInputWidth()+
			" targets=0 scores="+(scorer!=null)+(scorer!=null && scorer.dummy ? " PLUMBING_ONLY" : ""));
	}

	/** Returns a required nonempty value without accepting invisible defaults for resource bindings. */
	static String required(HashMap<String,String> options, String key){
		final String value=options.get(key);
		if(value==null || value.isEmpty()){throw new IllegalArgumentException("Required harness argument: "+key);}
		return value;
	}

	/** Verifies operator-facing hashes before any writer or prepared row is opened. */
	static String pinned(HashMap<String,String> options, String key, String hashKey){
		final String path=required(options, key), pin=required(options, hashKey);
		if(!pin.matches("[0-9a-f]{20}")){throw new IllegalArgumentException(hashKey+" requires 20 lowercase hex characters");}
		final boolean text=key.equals("familylist") || key.equals("subnetmanifest") ||
			key.equals("expectedcopytable") || key.equals("subnetpopulations");
		final String actual=text ? MagQCTextResource.sha80(path) : sha80(path);
		if(!actual.equals(pin)){throw new IllegalArgumentException("Harness resource sha80 mismatch: "+key);}
		return path;
	}

	/** Validates the whole input before publication; owned partials are cleaned even on parse or writer failure. */
	private static long write(String input, MagQCPreparedBin prepared, String output, String ledger, String scores, MagQCVectorMaker vm, Scorer scorer,
			HashMap<String,String> options) throws IOException{
		final boolean stdout="stdout".equals(scores) || "-".equals(scores);
		final Path vectors=path(output), ids=path(ledger), scoreFile=stdout ? null : path(scores);
		final HashSet<Path> destinations=new HashSet<Path>();
		for(Path destination:new Path[]{vectors, ids, scoreFile}){
			if(destination!=null && (!destinations.add(destination) || Files.exists(destination, LinkOption.NOFOLLOW_LINKS))){
				throw new IllegalArgumentException("Vector, identity and score outputs must be distinct fresh paths");
			}
		}
		long rows=0, identities=0;
		try(OwnedPartials owned=new OwnedPartials();
			PreparedReader reader=new PreparedReader(input, prepared, options.get("familylistsha80"), vm.preparedFamilyCount())){
			final Path vectorPart=vectors==null ? null : owned.create(vectors, ".magqc-vectors-");
			final Path idPart=ids==null ? null : owned.create(ids, ".magqc-identities-");
			final Path scorePart=scorer==null ? null : owned.create(scoreFile, ".magqc-scores-");
			try(Output vectorWriter=vectorPart==null ? null : new Output(vectorPart);
				Output idWriter=idPart==null ? null : new Output(idPart);
				Output scoreWriter=scorePart==null ? null : new Output(scorePart)){
				final ByteBuilder header=new ByteBuilder();
				if(vectorWriter!=null){
					header.append("#dims\t").append(vm.preparedInputWidth()).append("\t0\t0\n");
					appendBindings(header, options); vectorWriter.writer.print(header);
					header.clear(); header.append("#schema_version\tmagqc_composite_rows_v1\n");
					header.append("#input_width\t").append(vm.preparedInputWidth()).nl().append("#target_width\t0\n");
					appendBindings(header, options); header.append("#columns\trow_index\tbin_id\n");
					idWriter.writer.print(header);
				}
				if(scoreWriter!=null){header.clear(); scorer.header(header, options); scoreWriter.writer.print(header);}
				final ByteBuilder row=new ByteBuilder(vm.preparedInputWidth()*8+64), identity=new ByteBuilder();
				final ByteBuilder prediction=scorer==null ? null : new ByteBuilder(256);
				for(MagQCPreparedBin bin=reader.next(); bin!=null; bin=reader.next()){
					row.clear(); vm.formatPreparedBin(row, bin);
					if(scorer!=null){
						prediction.clear(); scorer.score(row, rows, bin.id, prediction); scoreWriter.writer.print(prediction);
					}
					if(vectorWriter!=null){
						vectorWriter.writer.print(row);
						identity.clear(); identity.append(identities).tab().append(bin.id).nl();
						idWriter.writer.print(identity); identities++;
					}
					rows++;
				}
				if(idWriter!=null && rows!=identities){throw new IllegalStateException("Composite vector/identity row conservation failed");}
			}
			// Close the input before publication too: a failed input close invalidates the run.
			reader.close();
			if(ids!=null){Files.move(idPart, ids);}
			if(scoreFile!=null){Files.move(scorePart, scoreFile);}
			// The vector file is the final publication point, after its identity ledger exists.
			if(vectors!=null){Files.move(vectorPart, vectors);}
			if(stdout){
				Files.copy(scorePart, System.out); System.out.flush();
				if(System.out.checkError()){throw new IOException("I/O error publishing scores to stdout");}
			}
		}
		return rows;
	}

	/** Null denotes an output that was not requested. */
	private static Path path(String value){return value==null ? null : Paths.get(value).toAbsolutePath().normalize();}

	/** Both paired outputs carry the same frozen resource identities and layout definition. */
	static void appendBindings(ByteBuilder out, HashMap<String,String> options){
		assert(options!=null) : "Output provenance must use the validated inference configuration";
		out.append("#layout\t").append(MagQCVectorMaker.JOINT_COMPOSITE_LAYOUT).nl();
		for(String key:PIN_KEYS){out.append('#').append(key).tab().append(options.get(key)).nl();}
		MagQCAssemblyInput.appendProvenance(out, options);
	}

	/** Uses the established native digest implementation without exposing the full digest. */
	static String sha80(String path){return ReferenceCdsSurvivalLabelReader.sha256File(path).substring(44);}

	/** Reads model text through native decompression; resource pins still cover file bytes. */
	private static byte[] readNetworkBytes(String path) throws IOException{
		final ByteFile input=ByteFile.makeByteFile(path, true);
		final ByteBuilder bytes=new ByteBuilder(65536);
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				bytes.append(line).nl();
			}
		}finally{
			if(input.close()){throw new IOException("I/O error reading final network: "+path);}
		}
		return bytes.toBytes();
	}

	/** Owns only temporary paths created by this invocation, never final or preexisting outputs. */
	private static final class OwnedPartials implements AutoCloseable {
		/** A null destination stages stdout in the configured temporary directory. */
		Path create(Path destination, String prefix) throws IOException{
			final Path temporary;
			if(destination==null){temporary=Files.createTempFile(prefix, ".partial");}
			else{
				Files.createDirectories(destination.getParent());
				temporary=Files.createTempFile(destination.getParent(), prefix, ".partial");
			}
			paths.add(temporary); return temporary;
		}
		/** Try-with-resources suppresses cleanup errors behind an earlier parse/write failure. */
		@Override public void close() throws IOException{
			IOException failure=null;
			for(Path temporary:paths){
				try{Files.deleteIfExists(temporary);}
				catch(IOException e){if(failure==null){failure=e;}else{failure.addSuppressed(e);}}
			}
			if(failure!=null){throw failure;}
		}
		final ArrayList<Path> paths=new ArrayList<Path>(3);
	}

	/** Validated final network, with one private worker and one reusable rounded input array. */
	static final class Scorer {
		Scorer(HashMap<String,String> options, int width) throws Exception{
			this(options, width, true);
		}

		/** Shared serving can discard the construction master after validating its one worker. */
		Scorer(HashMap<String,String> options, int width, boolean retainModel) throws Exception{
			this(options, width, retainModel, false);
		}

		/** Loads independently of the subnet bundle; the batch checks their widths after joining. */
		static Scorer loadConfigured(HashMap<String,String> options, boolean retainModel) throws Exception{
			return new Scorer(options, 0, retainModel, true);
		}

		/** Infers width only for parallel resource loading, never for ordinary supplied vectors. */
		private Scorer(HashMap<String,String> options, int width, boolean retainModel, boolean inferWidth) throws Exception{
			final String file=pinned(options, "net", "netsha80");
			final String flag=options.containsKey("dummy") ? required(options, "dummy") : "f";
			if(!Arrays.asList("t", "true", "f", "false").contains(flag)){
				throw new IllegalArgumentException("dummy= must be t or f");
			}
			dummy=flag.equals("t") || flag.equals("true");
			final MagQCNetBundle.InferenceModel parsed=MagQCNetBundle.parseInferenceModel(
				readNetworkBytes(file), "harness final net");
			model=retainModel ? parsed : null; sharedBinding=null;
			if(inferWidth){width=parsed.inputs();}
			if(parsed.inputs()!=width || parsed.outputs()!=(dummy ? 3 : 6)){
				throw new IllegalArgumentException("Final net requires input width "+width+" and "+(dummy ? 3 : 6)+" outputs");
			}
			worker=parsed.newWorker(); input=new float[width];
			output=new float[dummy ? 3 : 6];
			requireTag("magqc_harness_contract", "prepared_raw_v1");
			requireTag("magqc_layout", MagQCVectorMaker.JOINT_COMPOSITE_LAYOUT);
			requireTag("magqc_encoding", "raw");
			for(String key:PIN_KEYS){requireTag("magqc_"+key, options.get(key));}
			names=dummy ? DUMMY_NAMES : SIX_NAMES;
			requireTag("magqc_output_mode", dummy ? "dummy3" : "live_six");
			requireTag("magqc_output_names", String.join(",", names));
			requireTag("magqc_input_normalization", dummy ? "identity" : "embedded");
			if(dummy){
				if(worker.hasInputNormalization() || worker.getTag("mean_sd_folded")!=null || worker.getTag("normalization")!=null ||
					worker.getTag("affine_prefix_layers")!=null || worker.getTag("output_contract")!=null){
					throw new IllegalArgumentException("Dummy identity net cannot declare trained normalization/output semantics");
				}
			}else{
				// bbnet_export.py emits sorted compact JSON; pin its entire live_six semantics.
				requireTag("output_contract", SIX_CONTRACT);
				final String normalization=worker.getTag("normalization");
				final boolean headers=worker.hasInputNormalization();
				final boolean folded=!headers && "true".equals(worker.getTag("mean_sd_folded")) &&
					(normalization==null || "folded".equals(normalization)) && worker.getTag("affine_prefix_layers")==null;
				final boolean affine=!headers && "explicit-affine".equals(normalization) &&
					"false".equals(worker.getTag("mean_sd_folded")) && "2".equals(worker.getTag("affine_prefix_layers"));
				// CellNet.applyInput applies these validated float32 arrays once, before
				// trained edges. Reject mixed declarations that would normalize twice.
				final boolean explicitInput=headers && "explicit-input".equals(normalization) &&
					"false".equals(worker.getTag("mean_sd_folded")) && worker.getTag("affine_prefix_layers")==null;
				if(!folded && !affine && !explicitInput){
					throw new IllegalArgumentException("Final net must embed one consistent exported normalization representation");
				}
				if(affine){
					final int[] dims=parsed.dimensions();
					if(dims.length<4 || dims[1]!=width || dims[2]!=width){
						throw new IllegalArgumentException("Explicit-affine net lacks two input-width prefix layers");
					}
				}
			}
		}

		/** Copies inference scratch from one validated model without rereading the model or resource files. */
		private Scorer(Scorer source){
			model=source.model; sharedBinding=source.sharedBinding;
			worker=sharedBinding==null ? model.newWorker() : null;
			input=new float[source.input.length]; output=new float[source.output.length];
			names=source.names; dummy=source.dummy;
			assert((sharedBinding!=null || worker!=source.worker) && input!=source.input && output!=source.output) :
				"Batch final-network inference requires private activations and input scratch";
		}

		/** No model load: the first synchronized API call validates and loads this binding. */
		Scorer(MagQCCompositeInference.Binding binding){
			sharedBinding=binding; model=null; worker=null; dummy=false; names=SIX_NAMES;
			input=new float[binding.inputs()]; output=new float[6];
		}

		/** Private result storage keeps report formatting outside the composite lock. */
		float getOutput(int index){return output[index];}

		/** One caller retains this private scorer for its entire worker lifetime. */
		Scorer newWorker(){return new Scorer(this);}

		/** Resource identity and semantic mismatches fail before any output is opened. */
		private void requireTag(String key, String expected){
			if(!expected.equals(worker.getTag(key))){throw new IllegalArgumentException("Final net contract mismatch: "+key);}
		}

		/** Names, resource pins and dummy status travel with both file and stdout scores. */
		void header(ByteBuilder out, HashMap<String,String> options){
			out.append("#schema_version\tmagqc_prepared_scores_v1\n");
			appendBindings(out, options);
			out.append("#encoding\traw\n#netsha80\t").append(options.get("netsha80")).nl();
			out.append("#input_width\t").append(input.length).nl();
			out.append("#output_mode\t").append(dummy ? "dummy3" : "live_six").nl();
			out.append("#dummy\t").append(dummy ? "PLUMBING_ONLY" : "false").nl();
			if("synthetic_weights".equals(worker.getTag("magqc_fixture"))){out.append("#fixture\tSYNTHETIC_WEIGHTS_NOT_ACCURACY\n");}
			if(!dummy){out.append("#error_estimates\tpredicted absolute fraction-point errors; not confidence intervals\n");}
			out.append("#columns\trow_index\tbin_id");
			for(String name:names){out.tab().append(name);} out.nl();
		}

		/** Parses only used formatter bytes, retaining the exact rounded values seen by training. */
		void score(ByteBuilder row, long ordinal, String id, ByteBuilder out){
			evaluate(row, id);
			out.append(ordinal).tab().append(id);
			for(int i=0; i<names.length; i++){
				out.tab().appendSlow(getOutput(i));
			}
			out.nl();
		}

		/** Shared inference arithmetic for the developer harness and public assembly client. */
		void evaluate(ByteBuilder row, String id){
			assert(row.length()>0 && row.array[row.length()-1]=='\n') : "Prepared formatter must terminate one row";
			int first=0, term=0;
			for(int pos=0; pos<row.length(); pos++){
				final byte value=row.array[pos];
				if(value=='\t' || value=='\n'){
					if(term>=input.length || pos==first || (value=='\n' && pos!=row.length()-1)){
						throw new IllegalStateException("Malformed formatted composite row: "+id);
					}
					// Native decimal parsing accumulates digits in longs; unusually wide fields need its safe slow path.
					final float parsed=pos-first<=18 ? Parse.parseFloat(row.array, first, pos) :
						(float)Parse.parseDoubleSlow(row.array, first, pos);
					if(!Float.isFinite(parsed)){throw new IllegalArgumentException("Nonfinite final-net input: "+id);}
					input[term++]=parsed; first=pos+1;
				}
			}
			if(first!=row.length() || term!=input.length){throw new IllegalStateException("Formatted composite width mismatch: "+id);}
			evaluateInput(id);
		}

		/** Runs validated input scratch and captures outputs before another bin can overwrite them. */
		private void evaluateInput(String id){
			if(sharedBinding!=null){output=prok.ProkCC.calcCCSynced(input, sharedBinding); return;}
			worker.applyInput(input); worker.feedForward();
			for(int i=0; i<names.length; i++){
				final float value=worker.getOutput(i);
				output[i]=value;
				if(!Float.isFinite(value)){throw new IllegalStateException("Nonfinite final-net output "+i+": "+id);}
			}
		}
		/** Copies caller values into owned scratch and returns an independent snapshot. */
		float[] evaluateCopy(float[] values){
			if(values==null || values.length!=input.length){throw new IllegalArgumentException("Composite input width mismatch");}
			for(int i=0; i<values.length; i++){
				final float value=values[i];
				if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite composite input at "+i);}
				input[i]=value;
			}
			evaluateInput("synchronized API");
			return output.clone();
		}
		private final MagQCCompositeInference.Binding sharedBinding;
		private float[] output;

		final CellNet worker;
		private final MagQCNetBundle.InferenceModel model;
		final float[] input;
		final String[] names;
		final boolean dummy;
	}

	/** ByteFile-backed prepared schema reader; reuses one typed row across the entire stream. */
	private static final class PreparedReader implements AutoCloseable {
		PreparedReader(String path, MagQCPreparedBin prepared, String familyPin, int familyCount){
			if(prepared!=null){
				if(path!=null || prepared.families.length!=familyCount){
					throw new IllegalArgumentException("Assembly family width must match the bound subnet roster");
				}
				prepared.validate(); input=null; bin=prepared; return;
			}
			input=ByteFile.makeByteFile(path, true); bin=new MagQCPreparedBin(familyCount);
			try{
				requireLine("#schema_version\t"+MagQCPreparedBin.SCHEMA);
				requireLine("#family_list_sha80\t"+familyPin);
				requireLine("#columns\t"+MagQCPreparedBin.COLUMNS);
			}catch(RuntimeException|Error failure){
				try{close();}catch(RuntimeException cleanup){failure.addSuppressed(cleanup);}
				throw failure;
			}
		}
		/** Header order and exact vocabulary are part of the versioned input contract. */
		private void requireLine(String expected){
			final byte[] actual=input.nextLine();
			if(actual==null || !Arrays.equals(actual, expected.getBytes(java.nio.charset.StandardCharsets.UTF_8))){
				throw new IllegalArgumentException("Prepared schema/header binding mismatch");
			}
		}
		/** Exposes one checked reusable row; empty lines and extra headers are invalid data. */
		MagQCPreparedBin next(){
			if(closed){throw new IllegalStateException("Prepared input is closed");}
			if(input==null){if(consumed){return null;} consumed=true; return bin;}
			final byte[] line=input.nextLine();
			if(line==null){return null;}
			bin.parse(line); return bin;
		}
		/** Input close errors prevent final publication, including under try-with-resources cleanup. */
		@Override public void close(){
			if(closed){return;} closed=true;
			if(input!=null && input.close()){throw new RuntimeException("I/O error reading prepared bins");}
		}
		final ByteFile input;
		final MagQCPreparedBin bin;
		boolean closed;
		boolean consumed;
	}

	/** Makes the native asynchronous writer participate in suppressed-exception-safe cleanup. */
	private static final class Output implements AutoCloseable {
		Output(Path file){writer=new ByteStreamWriter(file.toString(), true, false, true); writer.start();}
		@Override public void close(){
			if(writer.poisonAndWait()){throw new RuntimeException("I/O error writing composite vectors/identities");}
		}
		final ByteStreamWriter writer;
	}

	private static final String[] PIN_KEYS={"bundlesha80", "familylistsha80", "subnetmanifestsha80",
		"expectedcopytablesha80", "subnetpopulationssha80"};
	private static final String[] DUMMY_NAMES={"dummy_0", "dummy_1", "dummy_2"};
	private static final String[] SIX_NAMES={"gene_completeness", "gene_contamination",
		"abs_residual_gene_completeness", "abs_residual_gene_contamination", "contamination_ge_0.01", "completeness_gt_0.99"};
	private static final String SIX_CONTRACT="{\"error_target_rule\":\"same_forward_absolute_residual_stop_gradient\",\"mode\":\"live_six\",\"outputs\":[\""+
		String.join("\",\"", SIX_NAMES)+"\"],\"selection\":\"value_only\",\"source_labels\":2}";
	private static final HashMap<String,Boolean> ALLOWED=new HashMap<String,Boolean>();
	static{
		for(String key:new String[]{"in", "out", "idsout", "bundle", "familylist", "subnetmanifest",
			"expectedcopytable", "subnetpopulations", "net", "netsha80", "scores", "dummy", "deterministic"}){ALLOWED.put(key, true);}
		for(String key:PIN_KEYS){ALLOWED.put(key, true);}
		for(String key:MagQCAssemblyInput.OPTIONS){ALLOWED.put(key, true);}
	}
}
