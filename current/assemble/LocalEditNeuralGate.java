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

/** Worker-local frozen-model veto for an already unique, fully verified local edit. */
final class LocalEditNeuralGate {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolve a literal path or ?resource, then verify its identity before loading. */
	static CellNet loadFrozen(final String path, final String expectedSha80, final float cutoff){
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

	LocalEditNeuralGate(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final CellNet net_, final float cutoff_){
		this(k_, lookup_, net_, cutoff_, null, Float.NaN);
	}
	/** Build a gate whose optional high-depth model is paired with the low-depth model. */
	LocalEditNeuralGate(final int k_, final HomopolymerIndelProposal.CountLookup lookup_,
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

	/** Validate the low/high pair as soon as Tadpole has loaded the frozen models. */
	static void validateRoutedPair(final CellNet lowNet, final CellNet highNet){
		if(highNet!=null){
			validateDepthModel(lowNet, "low-depth");
			validateDepthModel(highNet, "high-depth");
		}
	}
	private static void validateCutoff(final float cutoff){
		if(!Float.isFinite(cutoff) || cutoff<0 || cutoff>1){throw new IllegalArgumentException("fixindelscutoff must be finite in [0,1].");}
	}

	boolean accept(final byte[] bases, final int position, final LocalSingleBaseEdit.Operation operation, final byte base){
		final float score=score(bases, position, operation, base);
		return Float.isFinite(score) && score>=lastCutoff;
	}

	/** Freeze the source read for lazy percentile calculation before any edits.
	 * The caller starts a new transaction even when it reuses a bases array. */
	void beginRead(final byte[] bases){
		originalBases=bases;
		depthsReady=false;
		depthsAvailable=false;
		lastHighDepth=false;
		lastCutoff=Float.NaN;
	}

	/** Calculate whole-read context once, only if a 39-input model scores a locus. */
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

	/** Return NaN to abstain when any frozen feature is unavailable. */
	float score(final byte[] bases, final int position, final LocalSingleBaseEdit.Operation operation, final byte base){
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

	/** Select low-depth for P50<=30 and high-depth for P50>30 within one read transaction. */
	private CellNet selectedNet(){
		if(!dualNet){lastHighDepth=false; lastCutoff=cutoff; return net;}
		lastHighDepth=readPercentiles[MEDIAN_PERCENTILE_INDEX]>LOW_MAX_MEDIAN_DEPTH;
		lastCutoff=lastHighDepth ? highCutoff : cutoff;
		return lastHighDepth ? highNet : net;
	}

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

	private static boolean completeSubstitutions(final int[] depths, final int original){
		if(depths==null || depths.length!=4 || original<0 || original>3 || depths[original]>=0){return false;}
		for(int b=0; b<4; b++){if(b!=original && depths[b]<0){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean complete(final int[] depths){
		if(depths==null || depths.length!=4){return false;}
		for(final int depth:depths){if(depth<0){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

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
		}catch(final IOException e){throw new IllegalArgumentException("Could not read local-edit model: "+path, e);}
		catch(final NoSuchAlgorithmException e){throw new AssertionError("SHA-256 is required by the Java runtime.", e);}
	}
	private static boolean isSha80(final String value){
		if(value==null || value.length()!=20){return false;}
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(!((c>='0' && c<='9') || (c>='a' && c<='f'))){return false;}
		}
		return true;
	}

	long evaluations, probeQueries, windowQueries;
	float lastScore=Float.NaN;
	float lastCutoff=Float.NaN;
	boolean lastHighDepth=false;
	private final int[] candidateDepths=new int[LocalEditNeuralFeatures.CANDIDATE_COUNT];
	private final double[] features;
	private final float[] input;
	private final IntList depthScratch=new IntList();
	private final int[] readPercentiles=new int[LocalEditNeuralFeatures.PERCENTILE_COUNT];
	private byte[] originalBases;
	private boolean depthsReady, depthsAvailable;
	private final int k;
	final float cutoff, highCutoff;
	final boolean dualNet;
	private final HomopolymerIndelProposal.CountLookup lookup;
	private final CellNet net, highNet;
	private final Kmer key;
	private final LocalEditKmerProbe probe;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int MODEL_K=62;
	private static final int MEDIAN_PERCENTILE_INDEX=5;
	private static final int LOW_MAX_MEDIAN_DEPTH=30;
}
