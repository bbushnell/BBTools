package prot;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.HashMap;

import fileIO.ByteFile;
import ml.CellNet;
import ml.CellNetGraft;

/** A small relocatable descriptor binds an unchanged base and two native graft components.
 * Component files are pinned independently; the release's netsha80 pins the descriptor.
 * Loading validates every component before any worker can publish predictions.
 * @author Brian Bushnell, Yoimiya */
final class MagQCGraftModel {

	/** Null means an ordinary native network; an explicit graft header requires the complete schema. */
	static Descriptor describe(final String file) throws IOException{
		if(!hasDescriptorHeader(file)){return null;}
		final ByteFile input=ByteFile.makeByteFile(file, true);
		try{
			final byte[] first=input.nextLine();
			if(first==null || !HEADER.equals(new String(first, StandardCharsets.UTF_8))){return null;}
			final HashMap<String,String> fields=new HashMap<String,String>();
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				final String text=new String(line, StandardCharsets.UTF_8);
				final int tab=text.indexOf('\t');
				if(tab<1 || tab==text.length()-1 || fields.put(text.substring(0, tab), text.substring(tab+1))!=null){
					throw new IllegalArgumentException("Malformed or duplicate graft descriptor field");
				}
			}
			return new Descriptor(Paths.get(file).toAbsolutePath().getParent(), fields);
		}finally{
			if(input.close()){throw new IOException("Error reading graft descriptor: "+file);}
		}
	}

	/** Ordinary gzip models must not start a decompressor merely to inspect an optional descriptor. */
	private static boolean hasDescriptorHeader(final String file) throws IOException{
		try(InputStream input=Files.newInputStream(Paths.get(file))){
			for(int i=0; i<HEADER.length(); i++){if(input.read()!=HEADER.charAt(i)){return false;}}
			final int end=input.read();
			return end=='\n' || end=='\r';
		}
	}

	/** The base is loaded by Scorer through its ordinary path; only the two small extra masters live here. */
	MagQCGraftModel(final Descriptor descriptor) throws Exception{
		description=descriptor;
		privateModel=load(descriptor.privateFile, descriptor.privatePin, "graft private affine");
		headModel=load(descriptor.headFile, descriptor.headPin, "graft output heads");
		if(privateModel.outputs()!=headModel.outputs()*descriptor.privatePerHead || headModel.outputs()!=2){
			throw new IllegalArgumentException("Graft component outputs disagree with the declared two-head partition");
		}
	}

	/** Preserves native per-instance dense/sparse dispatch and private activation arrays. */
	CellNetGraft newWorker(final CellNet base){
		return new CellNetGraft(base, privateModel.newWorker(), headModel.newWorker(),
			description.upstreamLayer, description.shared, description.privatePerHead, description.activation);
	}

	/** Existing native parsing checks every weight; the descriptor also binds the complete file bytes. */
	private static MagQCNetBundle.InferenceModel load(final String path, final String pin, final String label) throws Exception{
		verify(path, pin);
		final MagQCNetBundle.InferenceModel model=MagQCNetBundle.parseInferenceModel(MagQCNetworkHarness.readNetworkBytes(path), label);
		requireNative18(model.newWorker());
		return model;
	}

	/** Native18 is mandatory for all serving components, including the unmodified base. */
	static void requireNative18(final CellNet model){
		if(model.weightBits()!=18){throw new IllegalArgumentException("Graft components require native18 precision");}
	}

	/** Pins are validated before reading a component as a model. */
	static void verify(final String path, final String pin){
		if(!pin.matches("[0-9a-f]{20}") || !pin.equals(MagQCNetworkHarness.sha80(path))){
			throw new IllegalArgumentException("Graft component file pin mismatch: "+path);
		}
	}

	/** Immutable interpretation metadata with filenames confined to the descriptor's directory. */
	static final class Descriptor {
		Descriptor(final Path directory, final HashMap<String,String> fields){
			if(fields.size()!=12){throw new IllegalArgumentException("Graft descriptor requires exactly twelve fields");}
			baseFile=component(directory, fields, "base"); basePin=required(fields, "base_sha80");
			privateFile=component(directory, fields, "private"); privatePin=required(fields, "private_sha80");
			headFile=component(directory, fields, "heads"); headPin=required(fields, "heads_sha80");
			upstreamLayer=positive(fields, "upstream_layer"); shared=positive(fields, "shared_nodes");
			privatePerHead=positive(fields, "private_per_head"); activation=required(fields, "hidden_activation");
			outputContract=required(fields, "output_contract"); predictionContract=required(fields, "prediction_contract");
			if(baseFile.equals(privateFile) || baseFile.equals(headFile) || privateFile.equals(headFile)){
				throw new IllegalArgumentException("Graft components must be distinct files");
			}
		}
		private static String required(final HashMap<String,String> fields, final String key){
			final String value=fields.get(key);
			if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing graft field: "+key);}
			return value;
		}
		private static int positive(final HashMap<String,String> fields, final String key){
			final int value=Integer.parseInt(required(fields, key));
			if(value<1){throw new IllegalArgumentException("Graft geometry must be positive: "+key);}
			return value;
		}
		private static String component(final Path directory, final HashMap<String,String> fields, final String key){
			final String value=required(fields, key);
			if(value.equals(".") || value.equals("..") || value.indexOf('/')>=0 || value.indexOf('\\')>=0){
				throw new IllegalArgumentException("Graft components must be sibling filenames: "+key);
			}
			return directory.resolve(value).toString();
		}
		final String baseFile, privateFile, headFile, basePin, privatePin, headPin;
		final String activation, outputContract, predictionContract;
		final int upstreamLayer, shared, privatePerHead;
	}

	private final Descriptor description;
	private final MagQCNetBundle.InferenceModel privateModel, headModel;
	private static final String HEADER="#magqc_frozen_graft_v1";
}
