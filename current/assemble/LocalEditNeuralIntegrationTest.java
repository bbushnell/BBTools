package assemble;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Random;

import dna.AminoAcid;
import ml.CellNet;
import stream.Read;
import ukmer.HashArrayU1D;
import ukmer.Kmer;

/** Original-orientation neural-veto wiring over one controlled verified substitution. */
public final class LocalEditNeuralIntegrationTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private LocalEditNeuralIntegrationTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		if(args.length!=3){throw new IllegalArgumentException("Expected model path, SHA80, and cutoff.");}
		final float cutoff=Float.parseFloat(args[2]);
		final int k=62, position=120;
		final String truth=sequence();
		final char replacement=truth.charAt(position)=='T' ? 'A' : 'T';
		final String query=truth.substring(0, position)+replacement+truth.substring(position+1);
		final Counts counts=new Counts(k); counts.add(truth, 40); counts.add(query, 1);
		final boolean generated=args[0].equals("generated39");
		final CellNet model;
		final CellNet highModel;
		if(generated){
			ml.Function.normalizeTypeRates();
			model=new CellNet(new int[]{39, 8, 1}, 17, 1, 1, 1, new ArrayList<String>());
			highModel=new CellNet(new int[]{39, 8, 1}, 19, 1, 1, 1, new ArrayList<String>());
			model.randomize(); highModel.randomize();
		}else{
			model=LocalEditNeuralGate.loadFrozen(args[0], args[1], cutoff);
			highModel=null;
		}
		int acceptedCount=0;
		for(final boolean rc:new boolean[]{false, true}){
			final byte[] raw=rc ? AminoAcid.reverseComplementBases(bytes(query)) : bytes(query);
			final byte[] expected=rc ? AminoAcid.reverseComplementBases(bytes(truth)) : bytes(truth);
			final int rawPosition=rc ? raw.length-1-position : position;
			final byte rawBase=rc ? complement((byte)truth.charAt(position)) : (byte)truth.charAt(position);
			final LocalEditNeuralGate oracle=new LocalEditNeuralGate(k, counts, model.copy(false), cutoff);
			oracle.beginRead(raw);
			final float expectedScore=oracle.score(raw, rawPosition, LocalSingleBaseEdit.Operation.SUBSTITUTION, rawBase);
			if(!Float.isFinite(expectedScore)){throw new AssertionError("Controlled complete evidence must score.");}
			final LocalEditNeuralGate integratedGate=new LocalEditNeuralGate(k, counts, model.copy(false), cutoff);
			integratedGate.beginRead(raw);
			final LocalEditCorrector corrector=new LocalEditCorrector(k, counts, 1, true, 1, integratedGate);
			final Read read=new Read(raw.clone(), null, "neural-orientation", 0, false);
			final int changed=corrector.correctOne(read);
			if(Float.floatToIntBits(integratedGate.lastScore)!=Float.floatToIntBits(expectedScore)){
				throw new AssertionError("Integrated gate did not score the original sequencing orientation.");
			}
			final boolean accepted=expectedScore>=cutoff;
			acceptedCount+=accepted ? 1 : 0;
			if(corrector.neuralEvaluations!=1 || corrector.neuralAcceptedLoci!=(accepted ? 1 : 0) ||
					corrector.neuralRejectedLoci!=(accepted ? 0 : 1) || changed!=(accepted ? 1 : 0)){
				throw new AssertionError("Neural veto accounting disagrees with the frozen score.");
			}
			if(accepted && !Arrays.equals(read.bases, expected)){
				throw new AssertionError("Accepted neural-veto path did not restore truth.");
			}
			if(!accepted && !Arrays.equals(read.bases, raw)){
				throw new AssertionError("Rejected neural-veto path changed the read.");
			}
			final LocalEditEngine transaction=new LocalEditEngine(k, new LocalEditEngine.CountLookup(){
				@Override public int count(final Kmer key){return counts.count(key);}
			}, 1, true, 1, model.copy(false), cutoff);
			final Read transactionalRead=new Read(raw.clone(), null, "neural-transaction", 2, false);
			final int edits=transaction.correct(transactionalRead, 8);
			if(edits!=(accepted ? 1 : 0) || !Arrays.equals(transactionalRead.bases, accepted ? expected : raw)){
				throw new AssertionError("Whole-read neural transaction disagrees with the direct gate.");
			}
			if(generated){
				final LocalEditEngine lowHighTransaction=new LocalEditEngine(k, new LocalEditEngine.CountLookup(){
					@Override public int count(final Kmer key){return counts.count(key);}
				}, 1, true, 1, model.copy(false), 0f, highModel.copy(false), 1f);
				final Read lowHighRead=new Read(raw.clone(), null, "dual-veto", 3, false);
				if(lowHighTransaction.correct(lowHighRead, 8)!=0 || !Arrays.equals(lowHighRead.bases, raw) ||
						lowHighTransaction.neuralEvaluations!=1 || lowHighTransaction.neuralRejectedLoci!=1){
					throw new AssertionError("Dual high-depth neural veto stopped rejecting after complete heuristic verification.");
				}
			}
		}
		if(acceptedCount<1){throw new AssertionError("Controlled high-contrast correction never exercised transactional apply.");}
		constructorDefenses(k, counts, model, cutoff);
		final LocalEditEngine engine=new LocalEditEngine(k, new LocalEditEngine.CountLookup(){
			@Override public int count(final Kmer key){return counts.count(key);}
		}, 1, true, 1, model.copy(false), cutoff);
		final byte[] blockedBases=bytes(query);
		final Read blocked=new Read(blockedBases, null, "pair-bypass", 1, false);
		try{engine.correct(blocked, 8, true);}
		catch(final IllegalArgumentException expected){
			if(blocked.bases!=blockedBases){throw new AssertionError("Rejected pair mode touched the read before failing.");}
			System.out.println("LOCAL_EDIT_NEURAL_INTEGRATION_TEST_OK orientations=2 accepted="+acceptedCount+" pair_bypass=rejected");
			return;
		}
		throw new AssertionError("Neural-enabled engine permitted pair lookahead to bypass the veto.");
	}

	private static void constructorDefenses(final int k, final Counts counts, final CellNet model, final float cutoff){
		for(int mode=0; mode<2; mode++){
			final LocalEditNeuralGate gate=new LocalEditNeuralGate(k, counts, model.copy(false), cutoff);
			try{new LocalEditCorrector(k, counts, mode==0 ? 3 : 1, mode==0, 1, gate);}
			catch(final IllegalArgumentException expected){continue;}
			throw new AssertionError("Package-level neural corrector accepted unsupported probe/competition configuration.");
		}
	}

	private static String sequence(){
		final Random random=new Random(202609211140L);
		final StringBuilder bb=new StringBuilder(244);
		for(int i=0; i<244; i++){bb.append("ACGT".charAt(random.nextInt(4)));}
		return bb.toString();
	}

	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static byte complement(final byte base){
		return base=='A' ? (byte)'T' : base=='C' ? (byte)'G' : base=='G' ? (byte)'C' : (byte)'A';
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact kmer-count fixture populated from controlled truth and query sequences. */
	private static final class Counts implements HomopolymerIndelProposal.CountLookup {
		Counts(final int k_){
			k=k_;
			final Kmer key=new Kmer(k);
			table=new HashArrayU1D(new int[]{2003}, key.k, k);
		}
		void add(final String text, final int copies){
			final Kmer key=new Kmer(k);
			for(int copy=0; copy<copies; copy++){
				key.clearFast();
				for(int i=0; i<text.length(); i++){key.addRight(text.charAt(i)); if(key.len()>=k){table.increment(key);}}
			}
		}
		@Override public int count(final Kmer key){return table.getValue(key);}
		final int k;
		final HashArrayU1D table;
	}
}
