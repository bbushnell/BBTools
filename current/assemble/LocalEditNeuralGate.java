package assemble;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import dna.Data;
import ml.CellNet;
import ml.CellNetParser;
import structures.IntList;
import ukmer.Kmer;

/** Worker-local model veto for an already unique, fully verified local edit.
 * This gate scores supplied proposals; it neither discovers nor applies edits.
 * Counts must stay stable and use the same K/hash/canonicalization as loading.
 * Networks and count provider are retained, not copied. Frozen weights do not
 * make mutable inference scratch thread-safe: use one gate/network set per worker.
 * Feature layouts and candidate order are defined by LocalEditNeuralFeatures.
 */
public final class LocalEditNeuralGate{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves a literal path or ?resource and checks its trailing SHA80 before loading.
	 * The file must remain stable between hashing and parsing. Validates a finite
	 * [0,1] cutoff and supported dimensions; a matching caller-supplied hash alone
	 * does not establish model calibration or suitability for a particular dataset.
	 * @return Newly loaded network with 28 or 39 inputs and one output
	 */
	public static CellNet loadFrozen(final String path, final String expectedSha80, final float cutoff){
		if(path==null || path.length()==0){throw new IllegalArgumentException("A frozen local-edit model path is required.");}
		if(!isSha80(expectedSha80)){throw new IllegalArgumentException("A lowercase 20-character model SHA80 is required.");}
		validateCutoff(cutoff);
		final String resolved=Data.findPath(path, false);
		if(resolved==null){throw new IllegalArgumentException("Cannot find local-edit model: "+path);}
		final String observed=sha80(resolved);
		if(!expectedSha80.equals(observed)){
			throw new IllegalArgumentException("Local-edit model SHA80 mismatch: expected "+expectedSha80+", observed "+observed+'.');
		}
		final CellNet net=CellNetParser.load(resolved, false);
		validateModel(net);
		return net;
	}

	/** Constructs a single-model gate; a 39-input model still needs original-read context. */
	public LocalEditNeuralGate(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final CellNet net_, final float cutoff_){
		this(k_, lookup_, net_, cutoff_, null, Float.NaN);
	}
	/** Constructs a K=62 gate and enables simdFF on the retained networks.
	 * A single model accepts either supported layout. A high model requires both
	 * networks to use 39 inputs; its cutoff must be finite in [0,1], or NaN when absent.
	 * This constructor checks dimensions/cutoffs, not file identity; loadFrozen owns
	 * the optional file-hash check. Callers supply separate mutable worker instances.
	 */
	public LocalEditNeuralGate(final int k_, final HomopolymerIndelProposal.CountLookup lookup_,
			final CellNet lowNet_, final float lowCutoff_, final CellNet highNet_, final float highCutoff_){
		if(k_!=MODEL_K || lookup_==null || lowNet_==null){
			throw new IllegalArgumentException("The frozen local-edit model requires K=62, immutable counts, and a loaded network.");
		}
		validateModel(lowNet_); validateCutoff(lowCutoff_);
		if((highNet_==null)!=Float.isNaN(highCutoff_)){
			throw new IllegalArgumentException("High-depth local-edit model and cutoff must be configured together.");
		}
		if(highNet_!=null){
			validateDepthModel(lowNet_, "low-depth");
			validateDepthModel(highNet_, "high-depth");
			validateCutoff(highCutoff_);
		}
		k=k_;
		lookup=lookup_;
		net=lowNet_;
		highNet=highNet_;
		dualNet=highNet_!=null;
		features=new double[dualNet ? LocalEditNeuralFeatures.DEPTH_FEATURE_COUNT : net.numInputs()];
		input=new float[features.length];
		cutoff=lowCutoff_;
		highCutoff=highCutoff_;
		net.simdFF=true;
		if(highNet!=null){highNet.simdFF=true;}
		key=new Kmer(k);
		probe=new LocalEditKmerProbe(k, lookup, 1);
		if(key.kbig!=k){throw new IllegalArgumentException("Local-edit neural K must equal the count-table K.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks only the supported feature/output dimensions of a loaded network. */
	private static void validateModel(final CellNet net){
		if(net==null || !LocalEditNeuralFeatures.supportedCount(net.numInputs()) || net.numOutputs()!=1){
			throw new IllegalArgumentException("Local-edit model dimensions must be 28 or 39 inputs and one output.");
		}
	}

	/** High-depth routing requires both loaded models to share the 39-input feature set. */
	private static void validateDepthModel(final CellNet net, final String label){
		if(net==null || net.numInputs()!=LocalEditNeuralFeatures.DEPTH_FEATURE_COUNT || net.numOutputs()!=1){
			throw new IllegalArgumentException(label+" local-edit model dimensions must be 39 inputs and one output.");
		}
	}

	/** Checks both routed layouts when a high model is present; otherwise does nothing. */
	static void validateRoutedPair(final CellNet lowNet, final CellNet highNet){
		if(highNet!=null){
			validateDepthModel(lowNet, "low-depth");
			validateDepthModel(highNet, "high-depth");
		}
	}
	/** Requires a finite inclusive score threshold in [0,1]. */
	private static void validateCutoff(final float cutoff){
		if(!Float.isFinite(cutoff) || cutoff<0 || cutoff>1){throw new IllegalArgumentException("fixindelscutoff must be finite in [0,1].");}
	}

	/** Accepts only a finite score at or above the selected model's inclusive cutoff. */
	public boolean accept(final byte[] bases, final int position, final LocalSingleBaseEdit.Operation operation, final byte base){
		final float score=score(bases, position, operation, base);
		return Float.isFinite(score) && score>=lastCutoff;
	}

	/** Borrows the original read for lazy depth percentiles; does not copy its bytes.
	 * Keep them unchanged throughout the transaction, even after edits install new arrays.
	 * Call before each new read transaction, including reused arrays. Resets depth/routing
	 * state and lastCutoff, but not cumulative counters or lastScore.
	 */
	public void beginRead(final byte[] bases){
		originalBases=bases;
		depthsReady=false;
		depthsAvailable=false;
		lastHighDepth=false;
		lastCutoff=Float.NaN;
	}

	/** Calculates original-read deciles once per transaction for a 39-input layout.
	 * A 28-input layout needs no whole-read context. Zero valid windows causes abstention;
	 * measured zero counts remain valid context. Lookup exceptions propagate.
	 */
	private boolean prepareReadDepths(){
		if(features.length==LocalEditNeuralFeatures.FEATURE_COUNT){return true;}
		if(originalBases==null){throw new IllegalStateException("A depth-aware neural gate requires beginRead before scoring.");}
		if(!depthsReady){
			final int queries=LocalEditNeuralFeatures.readDepthPercentiles(originalBases, key, lookup, depthScratch, readPercentiles);
			windowQueries+=queries;
			depthsAvailable=queries>0;
			depthsReady=true;
		}
		return depthsAvailable;
	}

	/** Scores a proposal on the supplied current sequence without editing it.
	 * Position addresses an existing base, or the gap before it for insertion;
	 * the terminal insertion gap is unsupported. Base is ignored for deletion.
	 * Requires all eight candidate depths and the feature encoder's context checks;
	 * unavailable evidence or nonfinite network output returns NaN, while exceptions propagate.
	 * Local evidence uses these current bytes; read-depth features/routing retain the
	 * original sequence supplied to beginRead. Resets lastScore/lastCutoff on entry.
	 */
	public float score(final byte[] bases, final int position, final LocalSingleBaseEdit.Operation operation, final byte base){
		lastScore=Float.NaN;
		lastCutoff=Float.NaN;
		if(bases==null || bases.length<k || position<0 || position>=bases.length || operation==null){return lastScore;}
		if(!prepareReadDepths()){return lastScore;}
		final int original=LocalEditNeuralFeatures.baseIndex(bases[position]);
		final int candidateBase=operation==LocalSingleBaseEdit.Operation.DELETION ? -1 :
				LocalEditNeuralFeatures.baseIndex(base);
		final int operationCode=operation==LocalSingleBaseEdit.Operation.SUBSTITUTION ? LocalEditNeuralFeatures.SUBSTITUTION :
				operation==LocalSingleBaseEdit.Operation.DELETION ? LocalEditNeuralFeatures.DELETION :
				LocalEditNeuralFeatures.INSERTION;
		final int candidate=LocalEditNeuralFeatures.candidateIndex(original, operationCode, candidateBase);
		if(candidate<0){return lastScore;}
		final int start=Math.max(0, Math.min(position-k/2, bases.length-k));
		probe.beginRead(bases);
		probe.substitutions(bases, start, position);
		probe.indels(); probeQueries+=probe.queries;
		if(!probe.substitutionsAvailable || !probe.deletionAvailable || !probe.insertionsAvailable ||
				!completeSubstitutions(probe.substitutionDepth, original) || probe.deletionDepth<0 || !complete(probe.insertionDepth)){
			return lastScore;
		}
		int next=0;
		for(int b=0; b<4; b++){if(b!=original){candidateDepths[next++]=probe.substitutionDepth[b];}}
		candidateDepths[next++]=probe.deletionDepth;
		for(int b=0; b<4; b++){candidateDepths[next++]=probe.insertionDepth[b];}
		assert(next==LocalEditNeuralFeatures.CANDIDATE_COUNT) : "Frozen candidate order has eight outcomes.";
		final int current=countWindow(bases, start);
		final int left=countWindow(bases, position-k);
		final int right=countWindow(bases, operation==LocalSingleBaseEdit.Operation.INSERTION ? position : position+1);
		final int run=operation==LocalSingleBaseEdit.Operation.SUBSTITUTION ?
				LocalEditNeuralFeatures.substitutionRun(bases, position, candidateBase) :
				operation==LocalSingleBaseEdit.Operation.DELETION ? LocalEditNeuralFeatures.deletionRun(bases, position) :
				LocalEditNeuralFeatures.insertionRun(bases, position, candidateBase);
		if(!LocalEditNeuralFeatures.fill(features, current, candidateDepths, candidate, original, left, right, run,
				LocalEditNeuralFeatures.lowComplexity(bases, start, k, position), bases.length, position)){
			return lastScore;
		}
		if(features.length==LocalEditNeuralFeatures.DEPTH_FEATURE_COUNT &&
				!LocalEditNeuralFeatures.appendReadDepths(features, readPercentiles)){return lastScore;}
		if(!LocalEditNeuralFeatures.toModelInput(features, input)){return lastScore;}
		lastScore=selectedNet().applyInput(input).feedForwardDense(); evaluations++;
		return Float.isFinite(lastScore) ? lastScore : Float.NaN;
	}

	/** Select low-depth for P50&lt;=30 and high-depth for P50>30 within one read transaction. */
	private CellNet selectedNet(){
		if(!dualNet){lastHighDepth=false; lastCutoff=cutoff; return net;}
		lastHighDepth=readPercentiles[MEDIAN_PERCENTILE_INDEX]>LOW_MAX_MEDIAN_DEPTH;
		lastCutoff=lastHighDepth ? highCutoff : cutoff;
		return lastHighDepth ? highNet : net;
	}

	/** Counts one complete defined window; returns -1 for bounds/undefined bases.
	 * A table's -1 absent sentinel becomes zero; depths below -1 are rejected.
	 */
	private int countWindow(final byte[] bases, final int start){
		if(start<0 || (long)start+k>bases.length){return -1;}
		key.clearFast();
		for(int i=start; i<start+k; i++){
			final int b=LocalEditNeuralFeatures.baseIndex(bases[i]);
			if(b<0){return -1;}
			key.addRightNumeric(b);
		}
		final int depth=lookup.count(key); windowQueries++;
		if(depth<-1){throw new IllegalStateException("Invalid count depth: "+depth);}
		return Math.max(0, depth);
	}

	/** Requires three measured substitutions and a negative original-base slot. */
	private static boolean completeSubstitutions(final int[] depths, final int original){
		if(depths==null || depths.length!=4 || original<0 || original>3 || depths[original]>=0){return false;}
		for(int b=0; b<4; b++){if(b!=original && depths[b]<0){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires one measured nonnegative insertion depth for each of A/C/G/T. */
	private static boolean complete(final int[] depths){
		if(depths==null || depths.length!=4){return false;}
		for(final int depth:depths){if(depth<0){return false;}}
		return true;
	}

	/** Hashes the file and returns the lowercase final ten SHA-256 bytes. */
	private static String sha80(final String path){
		try(InputStream in=Files.newInputStream(Paths.get(path))){
			final MessageDigest digest=MessageDigest.getInstance("SHA-256");
			final byte[] buffer=new byte[65536];
			for(int len=in.read(buffer); len>=0; len=in.read(buffer)){if(len>0){digest.update(buffer, 0, len);}}
			final byte[] bytes=digest.digest();
			final StringBuilder sb=new StringBuilder(20);
			for(int i=bytes.length-10; i<bytes.length; i++){
				final int value=bytes[i]&0xff;
				if(value<16){sb.append('0');}
				sb.append(Integer.toHexString(value));
			}
			return sb.toString();
		}catch(final IOException e){throw new IllegalArgumentException("Could not read local-edit model: "+path, e);
		}catch(final NoSuchAlgorithmException e){throw new AssertionError("SHA-256 is required by the Java runtime.", e);}
	}
	/** Checks the exact lowercase 20-hex-digit spelling, without accessing a file. */
	private static boolean isSha80(final String value){
		if(value==null || value.length()!=20){return false;}
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(!((c>='0' && c<='9') || (c>='a' && c<='f'))){return false;}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Cumulative completed inferences and reported probe work. beginRead does not
	 * reset these; exceptions may leave partial totals. */
	public long evaluations, probeQueries;
	/** Cumulative direct-window and original-read depth queries; beginRead does not
	 * reset it, and exceptions may leave partial totals. */
	public long windowQueries;
	/** Raw last network output, possibly nonfinite; NaN also marks early abstention.
	 * Reset on score entry, not by beginRead. Observe after the worker's score call;
	 * this gate owns the diagnostic value, which callers should not modify or read concurrently. */
	public float lastScore=Float.NaN;
	/** Selected cutoff, or NaN if the latest score call did not reach network selection. */
	public float lastCutoff=Float.NaN;
	/** Last selected route, reset by beginRead; an early abstention can leave it stale. */
	public boolean lastHighDepth=false;
	/** Serial scratch in the fixed substitution/deletion/insertion candidate order. */
	private final int[] candidateDepths=new int[LocalEditNeuralFeatures.CANDIDATE_COUNT];
	/** Raw feature vector and its rounded float representation for inference. */
	private final double[] features;
	private final float[] input;
	/** Reused original-read count scratch and nearest-rank deciles. */
	private final IntList depthScratch=new IntList();
	private final int[] readPercentiles=new int[LocalEditNeuralFeatures.PERCENTILE_COUNT];
	/** Borrowed transaction-start sequence; never replaced by an edited sequence. */
	private byte[] originalBases;
	/** Whether deciles were computed, and whether any defined windows were counted. */
	private boolean depthsReady, depthsAvailable;
	/** Effective table K; construction requires MODEL_K. */
	private final int k;
	/** Single/low and optional high inclusive thresholds; absent high cutoff is NaN. */
	final float cutoff, highCutoff;
	/** Whether original-read median depth selects between two models. */
	final boolean dualNet;
	/** Stable borrowed count source; queried keys are mutable scratch, not retained values. */
	private final HomopolymerIndelProposal.CountLookup lookup;
	/** Borrowed networks with worker-local inference state. */
	private final CellNet net, highNet;
	/** Reused direct-window key and eight-candidate probe. */
	private final Kmer key;
	private final LocalEditKmerProbe probe;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** K and routing boundary fixed by the trained feature contract. */
	static final int MODEL_K=62;
	private static final int MEDIAN_PERCENTILE_INDEX=5;
	private static final int LOW_MAX_MEDIAN_DEPTH=30;
}
