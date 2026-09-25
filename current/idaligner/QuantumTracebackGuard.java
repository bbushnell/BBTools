package idaligner;

import java.util.Arrays;
import java.util.Random;

import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Deterministic opposing-indel and trace/statistics regression cases.
 * record=t preserves unmodified output without suppressing structural failures.
 * @author Nilou
 */
public class QuantumTracebackGuard{

	public static void main(String[] args){
		args=new PreParser(args, System.err, null, false, true, false).args;
		boolean record=false;
		boolean wide=false;
		int trials=500;
		long seed=20260925L;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split[0].equals("record")){record=Boolean.parseBoolean(split[1]) || split[1].equals("t");}
			else if(split[0].equals("trials")){trials=Integer.parseInt(split[1]);}
			else if(split[0].equals("seed")){seed=Long.parseLong(split[1]);}
			else if(split[0].equals("wide")){wide=Boolean.parseBoolean(split[1]) || split[1].equals("t");}
			else if(split[0].equals("simd")){
				final boolean enable=Boolean.parseBoolean(split[1]) || split[1].equals("t");
				if(enable && !simd.Vector.simd256){throw new IllegalArgumentException("SIMD requested on unsupported hardware/JVM");}
				Shared.SIMD=enable;
			}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		assert(trials>0) : "At least one random trial is required to exercise the sparse traceback";
		final QuantumTracebackGuard guard=new QuantumTracebackGuard(record);
		final Random random=new Random(seed);
		if(!record){normalizationCases();}
		System.out.println("case\tquery\tref\tstart\tstop\treturnedIdentity\ttraceIdentity\tscore\tops\tplainIdentity\tplainStart\tplainStop\tplainScore\tplainDels");
		guard.run("single_match", "A".getBytes(), "A".getBytes());
		guard.run("single_sub", "A".getBytes(), "C".getBytes());
		guard.run("single_N", "N".getBytes(), "N".getBytes());
		guard.run("IUPAC", "ACGTRYKMSWBDHVNACGTNNRY".getBytes(), "ACGTRYKMSWBDHVNACGTNNRY".getBytes());
		guard.run("prok203_capture_26418567",
			"TCCTGGTGGCCATAGCGAGGTGGAAACACCCGTTCCCATCCCGAACACGGAAGTTAAGCCCCTCAGCGCCGATGGTACTGTGGGGTCACCCCGTGGGAGAGTAGGTCGCCGCCAGGC".getBytes(),
			"CAGAGCGGTAGGGAAACACCCGTACCCATTCCGAACACGGAAGTTAAGCCTACCAGCGTATCGTGAAGTACTGGAGTGAGCGATCCTCTGGGAACCACGAGTCGCCGCCTGCCCACTCATACGGGAACATTCCCACCACAGCCCACGGACAACGCGTCCGTGGGCTTTCTTCATATCCAACGAGCCGCAGGCCC".getBytes());
		for(int trial=0; trial<trials; trial++){
			final int len=wide ? 200+random.nextInt(1400) : 32+random.nextInt(225);
			final int pad=wide ? (trial%3)*16 : 32;
			final byte[] ref=randomSeq(random, len+(wide ? pad*2 : 128));
			final ByteBuilder query=new ByteBuilder();
			final double subRate=(trial%6)*0.1;
			for(int p=pad; p<pad+len; p++){
				final int event=random.nextInt(100);
				if(event<8){p+=random.nextInt(12);}
				else{
					if(event<16){
						final int ins=1+random.nextInt(12);
						for(int i=0; i<ins; i++){query.append(BASES[random.nextInt(BASES.length)]);}
					}
					query.append(random.nextDouble()<subRate ? BASES[random.nextInt(BASES.length)] : ref[p]);
				}
			}
			byte[] q=query.toBytes();
			if(q.length==0){q=new byte[]{'A'};}
			// Respect the sparse fill's documented qLen<=rLen operating assumption.
			final byte[] r=q.length<=ref.length ? ref : Arrays.copyOf(ref, q.length+16);
			if(r!=ref){for(int i=ref.length; i<r.length; i++){r[i]=BASES[random.nextInt(4)];}}
			guard.run("seed_"+seed+"_trial_"+trial, q, r);
		}
		System.err.println("cases="+guard.cases+" opposing="+guard.opposing+" identityMismatch="+guard.identityMismatch+
			" returnEstimateMismatch="+guard.returnEstimateMismatch+" rawRepairs="+guard.rawRepairs+
			" seed="+seed+" simd="+Shared.SIMD+" record="+record);
		if(!record && (guard.opposing>0 || guard.returnEstimateMismatch>0)){
			throw new AssertionError("Quantum trace must be canonical; return identity must preserve neutral-N scoring");
		}
	}

	QuantumTracebackGuard(final boolean record_){record=record_;}

	private void run(final String name, final byte[] query, final byte[] ref){
		final AlignmentStats stats=new AlignmentStats(true);
		final float identity=QuantumAligner.alignAndTraceStatic(query, ref, stats);
		final int[] pos=new int[4];
		final float plain=QuantumAligner.alignStatic(query, ref, pos);
		if(stats.opposingIndelsRepaired){
			rawRepairs++;
			assert(stats.score>pos[2] && stats.dels<pos[3]) : name+": each canceled opposing pair must improve score and remove a deletion";
		}else{
			assert(stats.score==pos[2] && stats.dels==pos[3] && Float.floatToIntBits(identity)==Float.floatToIntBits(plain)) :
				name+": unmodified trace must exactly preserve its non-trace score, deletions and return";
		}
		assert(!name.equals("prok203_capture_26418567") || stats.opposingIndelsRepaired) :
			"The saved production fixture must actually exercise opposing-gap repair";
		final String problem=TracerReconstructionGuard.checkReconstruction(stats, query, ref);
		assert(problem==null) : name+": "+problem;
		assert(stats.rStart==pos[0] && stats.rStop==pos[1]) : name+": trace must retain sparse-fill endpoints";
		int matches=0, subs=0, ins=0, dels=0, ns=0, qi=0, ri=stats.rStart;
		boolean mixed=false;
		byte previous=0;
		for(byte op : stats.matchString){
			mixed|=(previous=='I' && op=='D') || (previous=='D' && op=='I');
			previous=op;
			if(op=='I'){ins++; qi++;}
			else if(op=='D'){dels++; ri++;}
			else{
				final byte expected=(byte)(query[qi]=='N' || ref[ri]=='N' ? 'N' : query[qi]==ref[ri] ? 'm' : 'S');
				assert(op==expected) : name+": op "+(char)op+" != base-derived "+(char)expected+" at "+qi+","+ri;
				if(op=='m'){matches++;}
				else if(op=='N'){ns++;}
				else{subs++;}
				qi++; ri++;
			}
		}
		assert(stats.matches==matches && stats.subs==subs && stats.ins==ins && stats.dels==dels && stats.ns==ns) :
			name+": reported counts must agree with an independent trace walk";
		assert(stats.score==matches-subs-ins-dels) : name+": score must use Quantum's +1/-1/0 scoring";
		final float expectedIdentity=(matches+0.5f*ns)/stats.matchString.length;
		assert(Float.floatToIntBits(stats.identity)==Float.floatToIntBits(matches/(float)stats.matchString.length)) :
			name+": keep AlignmentStats.setFromMatchString's existing no-match-credit N convention";
		cases++;
		if(mixed){opposing++;}
		if(Float.floatToIntBits(identity)!=Float.floatToIntBits(stats.identity)){identityMismatch++;}
		if(Float.floatToIntBits(identity)!=Float.floatToIntBits(expectedIdentity)){returnEstimateMismatch++;}
		if(record || mixed){
			System.out.println(name+"\t"+new String(query)+"\t"+new String(ref)+"\t"+stats.rStart+"\t"+stats.rStop+
				"\t"+identity+"\t"+stats.identity+"\t"+stats.score+"\t"+new String(stats.matchString)+
				"\t"+plain+"\t"+pos[0]+"\t"+pos[1]+"\t"+pos[2]+"\t"+pos[3]);
		}
	}

	/** Hand-derived full-block fixtures include both orders, alternating runs and real bases. */
	private static void normalizationCases(){
		final String[][] cases={
			{"mDDIIm", "ACCA", "AGGA", "mSSm"},
			{"mDDDIm", "ACA", "AGGGA", "mSDDm"},
			{"mIIIDm", "ACCCA", "AGA", "mSIIm"},
			{"mIDDDIm", "ACCA", "AGGGA", "mSSDm"},
			{"mIDIDIDm", "ACCCA", "AGGGA", "mSSSm"},
			{"mDDDIIm", "ACCA", "AGGGA", "mSSDm"},
			{"mIIDDm", "ACCA", "AGGA", "mSSm"},
			{"mDIm", "ACA", "ACA", "mmm"},
			{"mIDm", "ANA", "AGA", "mNm"},
			{"mDIm", "ACA", "ANA", "mNm"},
			{"mIDm", "ANA", "ANA", "mNm"},
			{"mDIm", "AYA", "AYA", "mmm"},
			{"mDIm", "AYA", "AGA", "mSm"},
			{"mDIm", "AaA", "AaA", "mmm"},
			{"ID", "C", "G", "S"},
			{"DI", "C", "C", "m"},
			{"DII", "AC", "A", "mI"},
			{"IID", "NA", "A", "NI"},
			{"mIDmDIm", "ACACA", "AGAGA", "mSmSm"},
			{"mIIDDDm", "ACNA", "AGNGA", "mSNDm"},
			{"mmImmDDm", "AAAAAA", "AAAAAAA", "mmImmDDm"}
		};
		for(String[] c : cases){
			for(int start=0; start<=3; start+=3){
				final byte[] q=c[1].getBytes();
				final byte[] r=(start==0 ? c[2] : "TTT"+c[2]+"TT").getBytes();
				final byte[] original=c[0].getBytes();
				final byte[] normalized=Tracer.normalizeOpposingIndels(original.clone(), q, r, start);
				assert(Arrays.equals(normalized, c[3].getBytes())) : c[0]+" became "+new String(normalized)+", expected "+c[3];
				final AlignmentStats before=new AlignmentStats(), after=new AlignmentStats();
				before.setFromMatchString(original);
				after.setFromMatchString(normalized);
				assert(before.matches+before.subs+before.ns+before.ins==after.matches+after.subs+after.ns+after.ins) :
					"Normalization changed query consumption for "+c[0];
				assert(before.matches+before.subs+before.ns+before.dels==after.matches+after.subs+after.ns+after.dels) :
					"Normalization changed reference consumption for "+c[0];
				assert(after.score>=before.score) : "Replacing opposing gaps must not reduce the +1/-1/0 score";
				assert(!Tracer.hasOpposingIndels(normalized)) : "BaseGraph cannot consume opposing gaps: "+c[0];
				assert(Arrays.equals(normalized, Tracer.normalizeOpposingIndels(normalized.clone(), q, r, start))) :
					"Canonicalization must be idempotent for "+c[0];
			}
		}
		System.err.println("normalizationFixtures="+(cases.length*2));
	}

	private static byte[] randomSeq(final Random random, final int len){
		assert(len>0) : "Nonempty inputs are required by Quantum's first/last base access";
		final byte[] out=new byte[len];
		for(int i=0; i<len; i++){out[i]=BASES[random.nextInt(BASES.length)];}
		return out;
	}

	private final boolean record;
	private int cases, opposing, identityMismatch, returnEstimateMismatch, rawRepairs;
	private static final byte[] BASES={'A', 'C', 'G', 'T'};
}
