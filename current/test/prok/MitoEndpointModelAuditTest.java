package test.prok;


import java.util.Arrays;
import java.util.Random;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;

/** Terminal-edit immutability and actual native terminal-tie witness.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointModelAuditTest {
	public static void main(String[] args){
		final byte[] q=new byte[128],alphabet={'A','C','G','T'};final Random random=new Random(1729);
		for(int i=0;i<q.length;i++){q[i]=alphabet[random.nextInt(4)];}q[126]='G';q[127]='C';
		final byte[] original=q.clone();
		for(int end:new int[]{3,5}){final byte[][] v=MitoEndpointModelAudit.terminalVariants(q,end);
			check(Arrays.equals(q,original) && Arrays.equals(v[0],q) && v[0]!=q,"Diagnostic alternatives must not mutate/share baseline arrays");
			check(Arrays.equals(v[1],Arrays.copyOfRange(q,end==5?1:0,q.length-(end==3?1:0))),"Trim exactly one base from the declared end");
			for(int i=0;i<4;i++){check(v[i+2].length==129 && v[i+2][end==5?0:128]==alphabet[i],"All four appended/prepended bases must remain explicit alternatives");
				check(Arrays.equals(Arrays.copyOfRange(v[i+2],end==5?1:0,end==5?129:128),q),"Every original base is unchanged by a terminal extension");}
		}
		final byte[] ref=new byte[512];Arrays.fill(ref,(byte)'A');System.arraycopy(q,0,ref,80,q.length);
		final AlignmentStats exact=new AlignmentStats(true);QuantumAligner.alignAndTraceStatic(q,ref,exact);
		check(exact.rStart==80 && exact.rStop==207,"An exact control must retain both inclusive source ends");
		final byte[] extended=MitoEndpointModelAudit.terminalEdit(q,3,0,(byte)'A');final AlignmentStats extension=new AlignmentStats(true);
		QuantumAligner.alignAndTraceStatic(extended,ref,extension);
		check(extension.rStart==80 && extension.rStop==208 && extension.identity==1f,"A real matching extension changes the aligned end through sequence alignment, not post-hoc coordinate addition");
		final byte[] neutral=MitoEndpointModelAudit.terminalEdit(q,3,0,(byte)'N');
		for(byte base:alphabet){ref[208]=base;final AlignmentStats n=new AlignmentStats(true);final float returned=QuantumAligner.alignAndTraceStatic(neutral,ref,n);
			check(n.rStart==80 && n.rStop==208 && n.matchString[n.matchString.length-1]=='N' && n.ns==1 && !n.opposingIndelsWereRepaired(),
				"Quantum N_SCORE=0 beats INS=-1 on the exact shared prefix for every terminal reference base; it also consumes a flank beyond the old exact endpoint");
			check(n.score==128 && Math.abs(returned-128.5f/129)<0.000001f && Math.abs(n.identity-128f/129)<0.000001f,
				"Quantum return identity gives N half credit; trace stats identity gives zero match credit, after endpoint selection");
			check(MitoEndpointModelAudit.terminalNMargin(neutral.length,ref.length,207)==1,"Same-prefix final N consumption wins by one score unit where final-row positional penalties agree");
		}
		ref[208]='A';check(Arrays.equals(q,original),"N experiment must copy the model without changing any original base");
		final byte[] trimmed=MitoEndpointModelAudit.terminalEdit(q,5,2,(byte)0);final AlignmentStats trimming=new AlignmentStats(true);
		QuantumAligner.alignAndTraceStatic(trimmed,ref,trimming);
		check(trimming.rStart==82 && trimming.rStop==207 && Arrays.equals(q,original),"Trimming two leading model bases must move its true-match start by two without mutating the library");
		check(MitoEndpointModelAudit.supportStatus(false,0,(char)exact.matchString[exact.matchString.length-1],false,1,q[127],ref[207],ref[208]).equals("ALIGNED_MODEL_WITH_ANNOTATED_OVERHANG"),
			"A fully mapped model tip with an actual extra annotated A supplies extension evidence");
		ref[207]='T';final AlignmentStats mismatch=new AlignmentStats(true);QuantumAligner.alignAndTraceStatic(q,ref,mismatch);
		check(mismatch.rStart==80 && mismatch.rStop==206 && mismatch.matchString[mismatch.matchString.length-1]=='I',"Native equal-cost terminal substitution/insertion must expose the earlier-column tie in this fixture");
		check(MitoEndpointModelAudit.terminalTie(q,ref,mismatch.rStop),"Replacing final I by S consumes the truth base for the same terminal score");
		check(MitoEndpointModelAudit.supportStatus(false,1,'I',true,1,q[127],ref[206],ref[207]).equals("TERMINAL_I_EQUAL_SCORE_WITNESS"),"A terminal query insertion is not a completely aligned model tip with extension support");
		ref[207]='C';check(!MitoEndpointModelAudit.terminalTie(q,ref,206),"Matching next base is better, not an equal SUB/INS witness");
		ref[207]='N';check(!MitoEndpointModelAudit.terminalTie(q,ref,206),"N has a different score and cannot witness the substitution tie");
		q[127]='R';ref[207]='R';final AlignmentStats ambiguous=new AlignmentStats(true);QuantumAligner.alignAndTraceStatic(q,ref,ambiguous);
		check(ambiguous.rStop==207 && ambiguous.matchString[ambiguous.matchString.length-1]=='m',"Quantum treats equal non-N IUPAC bases as a match");
		check(MitoEndpointModelAudit.supportStatus(false,0,'m',false,1,q[127],ref[207],ref[208]).equals("AMBIGUOUS_TIP_OR_OVERHANG"),"An IUPAC terminal must not become unambiguous evidence merely because its native operation is m");
		System.out.println("MITO_MODEL_TERMINAL_FIXTURE_PASS immutable_edits_inclusive_exact_control_native_terminal_tie");
		System.out.println("MITO_MODEL_SOURCE_SUPPORT_FIXTURE_PASS mapped_overhang_terminal_insertion_iupac_tip");
		System.out.println("MITO_MODEL_NEUTRAL_N_FIXTURE_PASS all_four_bases_consumed_margin_one_identity_conventions_old_endpoint_overshoot");
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
