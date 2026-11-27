package prok;

import java.util.Arrays;

import dna.GeneticCode;
import parse.LineParser1;
import structures.ByteBuilder;

/** Effective native path settings and the actual loaded input PGM.
 * One worker owns this immutable policy and its reusable comparison scratch.
 * No model is synthesized and no PGM statistics are reconstructed from GFFs.
 * @author Brian, Raiden
 */
public final class GeneCallerPathPolicy{

	GeneCallerPathPolicy(final GeneCaller caller, final GeneticCode code_){
		if(caller.pgm==null){throw new IllegalArgumentException("Path replay requires the actual loaded input PGM");}
		model=caller.pgm;
		code=code_;
		instance(caller, state);
		statics(state);
		families=GeneCaller.ncrnaFamilies.toArray(new NcrnaFamily[GeneCaller.ncrnaFamilies.size()]);
	}

	/** Fresh counters and mutable caller scratch, with the captured model/code. */
	public GeneCaller newCaller(){
		validateStatic();
		final GeneCaller caller=new GeneCaller(state[0], state[1], state[2],
			value(state[3]), value(state[4]), value(state[5]), value(state[6]), value(state[7]), model, code);
		validate(caller);
		return caller;
	}

	/** Fail before a trial if caller, global geometry or family registration changed. */
	public void validate(final GeneCaller caller){
		if(!sameModel(caller)){
			throw new IllegalArgumentException("Replay must use the same actual input PGM and genetic code");
		}
		caller.verifyPathReplayPolicy();
		instance(caller, scratch);
		statics(scratch);
		compare();
		registry();
	}
	public boolean sameModel(final GeneCaller caller){return caller!=null && caller.pgm==model && caller.pathGeneticCode()==code;}

	private void validateStatic(){
		System.arraycopy(state, 0, scratch, 0, 8);
		statics(scratch);
		compare();
		registry();
	}
	private void compare(){
		if(!Arrays.equals(state, scratch)){
			for(int i=0; i<state.length; i++){
				if(state[i]!=scratch[i]){
					throw new IllegalStateException("Native path setting changed after capture: "+FIELDS[i]);
				}
			}
		}
	}
	private void registry(){
		if(families.length!=GeneCaller.ncrnaFamilies.size()){
			throw new IllegalStateException("Replay must retain the native registered-family counter axis");
		}
		for(int i=0; i<families.length; i++){
			if(families[i]!=GeneCaller.ncrnaFamilies.get(i)){
				throw new IllegalStateException("Registered family binding changed at index "+i);
			}
		}
	}

	/** Versioned numeric field order; float values are raw IEEE bits. */
	public ByteBuilder append(final ByteBuilder row){
		assert(row!=null) : "The capture writer owns a reusable output buffer";
		row.append("POLICY\t1");
		for(final int number:state){row.tab().append(number);}
		return row;
	}
	/** A saved graph is admitted only under the matching loaded runtime settings. */
	public void validateRecord(final LineParser1 row){
		if(row.terms()!=state.length+2 || !row.termEquals("POLICY", 0) || row.parseInt(1)!=1){
			throw new IllegalArgumentException("Unsupported or incomplete native path policy record");
		}
		validateStatic();
		for(int i=0; i<state.length; i++){
			if(state[i]!=row.parseInt(i+2)){
				throw new IllegalArgumentException("Saved graph has a different effective path setting: "+FIELDS[i]);
			}
		}
	}

	private static void instance(final GeneCaller caller, final int[] out){
		assert(caller!=null && out.length==FIELDS.length) : "Policy columns bind the native constructor's complete path configuration";
		out[0]=caller.minLen;
		out[1]=caller.maxOverlapSameStrand;
		out[2]=caller.maxOverlapOppositeStrand;
		out[3]=bits(caller.minStartScore);
		out[4]=bits(caller.minStopScore);
		out[5]=bits(caller.minInnerScore);
		out[6]=bits(caller.minOrfScore);
		out[7]=bits(caller.minAvgScore);
	}
	private static void statics(final int[] out){
		assert(out.length==FIELDS.length) : "Saved policy field order must cover every input to native overlap/transition scoring";
		out[8]=GeneModel.kInnerCDS;
		out[9]=bits(Orf.e1);
		out[10]=bits(Orf.e2);
		out[11]=bits(Orf.e3);
		out[12]=bits(Orf.f1);
		out[13]=bits(Orf.f2);
		out[14]=bits(Orf.f3);
		out[15]=GeneCaller.lookbackPlus;
		out[16]=GeneCaller.lookbackMinus;
		out[17]=bits(GeneCaller.p0);
		out[18]=bits(GeneCaller.p1);
		out[19]=bits(GeneCaller.p2);
		out[20]=bits(GeneCaller.p3);
		out[21]=bits(GeneCaller.p4);
		out[22]=bits(GeneCaller.p5);
		out[23]=bits(GeneCaller.p6);
		out[24]=bits(GeneCaller.q1);
		out[25]=bits(GeneCaller.q2);
		out[26]=bits(GeneCaller.q3);
		out[27]=bits(GeneCaller.q4);
		out[28]=bits(GeneCaller.q5);
		out[29]=bits(GeneCaller.pathMinorStartMultiplier());
		out[30]=GeneCaller.FLAT_RNA_TRANSITIONS ? 1 : 0;
		out[31]=GeneCaller.EUK5S_ARRAY_BONUS ? 1 : 0;
		out[32]=GeneCaller.EUK5S_ARRAY_MAX_GAP;
	}
	private static int bits(final float value){return Float.floatToRawIntBits(value);}
	private static float value(final int bits){return Float.intBitsToFloat(bits);}
	private final GeneModel model;
	private final GeneticCode code;
	private final NcrnaFamily[] families;
	private final int[] state=new int[FIELDS.length], scratch=new int[FIELDS.length];
	private static final String[] FIELDS={"minimum_length", "same_strand_overlap", "opposite_strand_overlap",
		"minimum_start", "minimum_stop", "minimum_inner", "minimum_orf", "minimum_average", "inner_k",
		"e1", "e2", "e3", "f1", "f2", "f3", "lookback_plus", "lookback_minus",
		"p0", "p1", "p2", "p3", "p4", "p5", "p6", "q1", "q2", "q3", "q4", "q5",
		"minor_start_multiplier", "flat_rna", "euk5s_array_bonus", "euk5s_array_gap"};
}
