package assemble;

import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;

import ml.CellNet;

/** Default-off and fail-loud CLI/config contract for the frozen local-edit model. */
public final class TadpoleNeuralConfigTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private TadpoleNeuralConfigTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args) throws Exception{
		if(args.length!=2){throw new IllegalArgumentException("Expected BBTools root and selected.bbnet path.");}
		ml.Function.normalizeTypeRates();
		final Path dir=Files.createTempDirectory("tadpole-dual-neural-config.");
		final Path low39=dir.resolve("low39.bbnet"), high39=dir.resolve("high39.bbnet");
		writeNet(low39, new CellNet(new int[]{39, 8, 1}, 17, 1, 1, 1, new ArrayList<String>()));
		writeNet(high39, new CellNet(new int[]{39, 8, 1}, 19, 1, 1, 1, new ArrayList<String>()));
		final String low39Sha80=sha80(low39), high39Sha80=sha80(high39);
		final Tadpole valid=make(args[0], args[1], "fixindels=t", "k=62");
		if(!valid.localEdit || valid.localEditNeuralNet==null || valid.localEditPairs || valid.localEditBatch){
			throw new AssertionError("Valid frozen neural config was not retained exactly.");
		}
		if(Float.floatToRawIntBits(valid.localEditNeuralCutoff)!=Float.floatToRawIntBits(0.7466576f)){
			throw new AssertionError("Explicit deployment cutoff was not retained exactly.");
		}
		final Tadpole dual=makeWithModel(args[0], low39.toString(), low39Sha80, "0.25", "fixindels=t", "k=62",
				"fixindelshighnet="+high39, "fixindelshighsha80="+high39Sha80, "fixindelshighcutoff=0.75");
		if(dual.localEditNeuralNet==null || dual.localEditNeuralHighNet==null ||
				Float.floatToRawIntBits(dual.localEditNeuralHighCutoff)!=Float.floatToRawIntBits(0.75f)){
			throw new AssertionError("Valid dual frozen neural config was not retained exactly.");
		}
		reject(args[0], args[1], "k=62");
		reject(args[0], args[1], "fixindels=t", "k=31");
		reject(args[0], args[1], "fixindels=t", "k=62", "fixindelspairs=t");
		reject(args[0], args[1], "fixindels=t", "k=62", "fixindelsbatch=t", "fixindelsstride=1");
		reject(args[0], args[0]+"/README.md", "fixindels=t", "k=62");
		reject(args[0], args[1], "fixindels=t", "k=62", "fixindelscutoff=2");
		reject(args[0], args[1], "fixindels=t", "k=62", "fixindelscutoff=NaN");
		reject(args[0], args[1], "fixindels=t", "k=62", "fixindelssha80=BAD");
		rejectRaw(args[0], "fixindelshighcutoff=NaN");
		rejectRaw(args[0], "fixindels=t", "k=62", "fixindelshighnet="+high39,
				"fixindelshighsha80="+high39Sha80, "fixindelshighcutoff=0.75");
		rejectModel(args[0], low39.toString(), low39Sha80, "0.25", "fixindels=t", "k=62", "fixindelshighnet="+high39);
		rejectModel(args[0], low39.toString(), low39Sha80, "0.25", "fixindels=t", "k=62",
				"fixindelshighnet="+high39, "fixindelshighsha80=00000000000000000000", "fixindelshighcutoff=0.75");
		rejectModel(args[0], args[1], "05d0812fb6b4cbf1c5fe", "0.7466576", "fixindels=t", "k=62",
				"fixindelshighnet="+high39, "fixindelshighsha80="+high39Sha80, "fixindelshighcutoff=0.75");
		checkPacBio(args[0], dir);
		writeSmokeReads(java.nio.file.Paths.get(args[0], "pacbio_smoke.fa"));
		System.out.println("TADPOLE_NEURAL_CONFIG_TEST_OK valid=2 rejected=13 pacbio=PASS");
	}

	/** Check bundled defaults, override order, config expansion, and explicit failures. */
	private static void checkPacBio(final String root, final Path dir) throws Exception{
		checkPreset(makeRaw(root, "pacbio"), 0.93f, 0.96f);
		checkPreset(makeRaw(root, "pacbio=t"), 0.93f, 0.96f);
		checkPreset(makeRaw(root, "k=62", "ecc", "pacbio"), 0.93f, 0.96f);
		checkPreset(makeRaw(root, "pacbio", "ecc=t", "k=62"), 0.93f, 0.96f);
		checkPreset(makeRaw(root, "ecct", "pacbio"), 0.93f, 0.96f);
		final Tadpole legacy=makeRaw(root, "ecc", "pacbio=f");
		if(!legacy.ecc || legacy.localEdit || legacy.localEditNeuralNet!=null){
			throw new AssertionError("ecc without PacBio must keep the legacy correction path.");
		}
		checkPreset(makeRaw(root, "-pacbio"), 0.93f, 0.96f);
		checkPreset(makeRaw(root, "--pacbio=t"), 0.93f, 0.96f);
		checkPreset(makeRaw(root, "fixindelscutoff=0.91", "fixindelshighcutoff=0.95", "pacbio"), 0.91f, 0.95f);
		checkPreset(makeRaw(root, "pacbio", "fixindelscutoff=0.92", "fixindelshighcutoff=0.94"), 0.92f, 0.94f);
		final Path config=dir.resolve("pacbio.config");
		Files.write(config, java.util.Arrays.asList("ecc", "pacbio=t", "fixindelscutoff=0.92"));
		checkPreset(makeRaw(root, "config="+config, "fixindelshighcutoff=0.95"), 0.92f, 0.95f);
		final Tadpole off=makeRaw(root, "pacbio=f");
		if(off.localEdit || off.localEditNeuralNet!=null || off.localEditNeuralHighNet!=null || off.kbig!=31){
			throw new AssertionError("pacbio=f must preserve the legacy default-off correction and K=31.");
		}
		final Tadpole disabled=makeRaw(root, "pacbio", "pacbio=f");
		if(disabled.localEdit || disabled.localEditNeuralNet!=null){
			throw new AssertionError("The last pacbio boolean must select whether defaults are applied.");
		}
		rejectRaw(root, "pacbio", "k=31");
		rejectRaw(root, "k=31", "pacbio");
		rejectRaw(root, "pacbio", "fixindelspairs=t");
		rejectRaw(root, "pacbio", "fixindelsbatch=t", "fixindelsstride=1");
		rejectRaw(root, "pacbio", "fixindels=f");
		rejectRaw(root, "pacbio", "ecc", "mode=contig");
		rejectRaw(root, "mode=contig", "ecc", "pacbio");
		rejectRaw(root, "pacbio", "fixindelscutoff=1.01");
		rejectRaw(root, "pacbio", "fixindelshighcutoff=NaN");
		rejectRaw(root, "pacbio", "fixindelssha80=00000000000000000000");
		rejectRaw(root, "pacbio", "fixindelsnet=?missing_pacbio_test_model.bbnet");
	}

	/** Verify dispatch and the exact cutoffs actually consumed by worker gates. */
	private static void checkPreset(final Tadpole tadpole, final float low, final float high){
		assert(tadpole!=null) : "The preset test must inspect the constructed Tadpole, not just expanded arguments.";
		if(!(tadpole instanceof Tadpole2) || tadpole.kbig!=62 || !tadpole.localEdit || tadpole.ecc ||
				tadpole.processingMode!=Tadpole.correctMode ||
				tadpole.localEditNeuralNet==null || tadpole.localEditNeuralHighNet==null ||
				tadpole.localEditNeuralCutoff!=low || tadpole.localEditNeuralHighCutoff!=high){
			throw new AssertionError("PacBio preset lost K=62, correction mode, bundled models or explicit cutoffs.");
		}
		if(tadpole.localEditNeuralNet.numInputs()!=39 || tadpole.localEditNeuralHighNet.numInputs()!=39){
			throw new AssertionError("PacBio depth routing requires both bundled 39-input models.");
		}
	}

	/** Generate two depth families with known single-base errors for launcher parity. */
	private static void writeSmokeReads(final Path path) throws Exception{
		final java.util.Random random=new java.util.Random(20260930L);
		final structures.ByteBuilder reads=new structures.ByteBuilder();
		for(final int depth:new int[]{20, 40}){
			final StringBuilder sequence=new StringBuilder(400);
			for(int i=0; i<400; i++){sequence.append("ACGT".charAt(random.nextInt(4)));}
			final String truth=sequence.toString();
			for(int i=0; i<depth; i++){
				reads.append(">truth_").append(depth).append('_').append(i).nl().append(truth).nl();
			}
			final char changed=truth.charAt(200)=='A' ? 'C' : 'A';
			reads.append(">sub_").append(depth).nl().append(truth.substring(0, 200))
					.append(changed).append(truth.substring(201)).nl();
			reads.append(">ins_").append(depth).nl().append(truth.substring(0, 200))
					.append('A').append(truth.substring(200)).nl();
			reads.append(">del_").append(depth).nl().append(truth.substring(0, 200))
					.append(truth.substring(201)).nl();
		}
		Files.write(path, reads.toBytes());
	}

	/** Assert that a standard fixture plus explicit raw options is rejected. */
	private static void rejectRaw(final String root, final String... options){
		try{makeRaw(root, options);}
		catch(final IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid frozen neural config did not fail loudly.");
	}

	private static void reject(final String root, final String model, final String... options){
		try{make(root, model, options);}
		catch(final IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid frozen neural config did not fail loudly.");
	}

	/** Assert that an invalid low-model trio plus extra options is rejected. */
	private static void rejectModel(final String root, final String model, final String sha80, final String cutoff, final String... options){
		try{makeWithModel(root, model, sha80, cutoff, options);}
		catch(final IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid frozen neural config did not fail loudly.");
	}

	private static Tadpole make(final String root, final String model, final String... options){
		return makeWithModel(root, model, "05d0812fb6b4cbf1c5fe", "0.7466576", options);
	}

	/** Build a Tadpole from one explicit frozen low-model trio plus extra options. */
	private static Tadpole makeWithModel(final String root, final String model, final String sha80, final String cutoff, final String... options){
		final ArrayList<String> args=baseArgs(root);
		args.add("fixindelsnet="+model);
		args.add("fixindelssha80="+sha80);
		args.add("fixindelscutoff="+cutoff);
		for(final String option:options){args.add(option);}
		return Tadpole.makeTadpole(args.toArray(new String[args.size()]), true);
	}

	/** Build a Tadpole from options without implying a frozen low-model trio. */
	private static Tadpole makeRaw(final String root, final String... options){
		final ArrayList<String> args=baseArgs(root);
		for(final String option:options){args.add(option);}
		return Tadpole.makeTadpole(args.toArray(new String[args.size()]), true);
	}
	/** Supply the small single-threaded input fixture shared by all config cases. */
	private static ArrayList<String> baseArgs(final String root){
		final ArrayList<String> args=new ArrayList<String>();
		args.add("in="+root+"/testdata/crossk_left_bridge_reads.fa");
		args.add("out=null"); args.add("t=1"); args.add("prealloc=f"); args.add("prefilter=f"); args.add("pop=f");
		return args;
	}

	/** Write a deterministic temporary network through CellNet's native serializer. */
	private static void writeNet(final Path path, final CellNet net) throws Exception{
		net.randomize();
		Files.write(path, net.toBytes().toBytes());
	}

	/** Return the lowercase 20-hex SHA-256 suffix used by frozen model binding. */
	private static String sha80(final Path path) throws Exception{
		final byte[] hash=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(path));
		final StringBuilder sb=new StringBuilder(20);
		for(int i=hash.length-10; i<hash.length; i++){
			final int b=hash[i]&0xff;
			if(b<16){sb.append('0');}
			sb.append(Integer.toHexString(b));
		}
		return sb.toString();
	}
}
