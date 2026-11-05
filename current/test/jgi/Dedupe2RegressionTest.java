package test.jgi;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

import align2.BandedAligner;
import align2.BandedAlignerConcrete;
import dna.AminoAcid;
import jgi.Dedupe2;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import stream.Read;

/** Saved-overlap and real HashThread failure regression controls.
 * Reflection supplies finite internal states without a production fault-injection flag.
 * Run workerfail in a separate JVM: successful failure propagation halts that JVM.
 * @author Brian Bushnell, Raiden
 */
public final class Dedupe2RegressionTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, Dedupe2RegressionTest.class, false);
		String mode="bounds", input=null, output=null;
		for(String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0], b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){mode=b;}
			else if(a.equals("in")){input=b;}
			else if(a.equals("out")){output=b;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		check(Shared.EA(), "Regression checks require assertions enabled");
		check(input!=null && output!=null, "Supply an input fixture and fresh output path");
		final Dedupe2 tool=new Dedupe2(new String[]{"in="+input, "out="+output,
			"k=13", "nam=3", "t=1", "e=3", "bandwidth=9", "minidentity=99",
			"am=f", "ac=f", "fo=t", "cluster=t", "mo=1", "mop=95", "pc=f", "pto=f", "ow=f"});
		if(mode.equals("workerfail")){
			workerFailure(tool);
		}else if(mode.equals("bounds")){
			check((Integer)field(Dedupe2.class, "bandwidth").get(tool)==7,
				"The ordinary edit-three band remains seven columns");
			pair(tool, FIRST_A, FIRST_B, 111, 108);
			pair(tool, SECOND_A, SECOND_B, 101, 98);
			boundedApi(FIRST_A, FIRST_B, 111, 108);
			boundedApi(SECOND_A, SECOND_B, 101, 98);
			fixedEndControls();
			System.err.println("DEDUPE2_BOUNDED_OVERLAP_TEST_PASS");
		}else{throw new IllegalArgumentException("Unknown mode: "+mode);}
		Shared.closeStream(pp.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Tests            ----------------*/
	/*--------------------------------------------------------------*/

	private static void pair(final Dedupe2 tool, final String a, final String b,
		final int astop, final int bstop) throws Exception{
		final Class<?> unit=Class.forName("jgi.Dedupe2$Unit");
		final Class<?> overlap=Class.forName("jgi.Dedupe2$Overlap");
		final Constructor<?> uc=unit.getDeclaredConstructor(Dedupe2.class, Read.class);
		uc.setAccessible(true);
		final Read ar=new Read(bytes(a), null, "a", 0);
		final Read br=new Read(bytes(b), null, "b", 1);
		final Object ua=uc.newInstance(tool, ar), ub=uc.newInstance(tool, br);
		final Constructor<?> oc=overlap.getDeclaredConstructor(Dedupe2.class, unit, unit,
			int.class, int.class, int.class, int.class, int.class, int.class, int.class, int.class, BandedAligner.class);
		oc.setAccessible(true);
		final BandedAligner band=BandedAligner.makeBandedAligner(7);
		final Object o=oc.newInstance(tool, ua, ub, Dedupe2.FORWARD, 0, 0, astop, bstop, astop+1, 3, 3, band);
		final Method test=method(overlap, "test", BandedAligner.class, int.class);
		final Method swap=method(overlap, "swap"), reverse=method(overlap, "reverseDirection");
		for(int strand=0; strand<2; strand++){
			for(int direction=0; direction<2; direction++){
				for(int order=0; order<2; order++){
					check((Boolean)test.invoke(o, band, 6), "Saved overlap must validate in every order and orientation");
					swap.invoke(o);
				}
				reverse.invoke(o);
			}
			if(strand==0){
				final Object changed=field(overlap, "u2").get(o);
				method(unit, "reverseComplement").invoke(changed);
				method(overlap, "flip", unit, BandedAligner.class).invoke(o, changed, band);
			}
		}
		// Corrupt the saved span, not its outside tail: the validator must reject it.
		Arrays.fill(ar.bases, (byte)'A');
		check(!(Boolean)test.invoke(o, band, 6), "Bounded validation must still detect a corrupted accepted span");
	}

	private static void boundedApi(final String a, final String b, final int astop, final int bstop){
		assert(astop<a.length() && bstop<b.length()) : "Fixture stops must leave the documented unaligned tail";
		for(int direction=0; direction<4; direction++){
			final boolean qforward=(direction==0 || direction==3), rforward=(direction<2);
			byte[] query=bytes(b), ref=bytes(a);
			if(!qforward){AminoAcid.reverseComplementBasesInPlace(query);}
			if(!rforward){AminoAcid.reverseComplementBasesInPlace(ref);}
			final int qs=qforward ? 0 : query.length-1;
			final int rs=rforward ? 0 : ref.length-1;
			final int qe=qforward ? bstop : query.length-1-bstop;
			final int re=rforward ? astop : ref.length-1-astop;
			final Probe actual=new Probe(7);
			final int distance=align(actual, direction, query, ref, qs, rs, qe, re);
			check(distance==3 && actual.lastQueryLoc==qe && actual.lastRefLoc==re,
				"The real pair must consume its saved endpoints at three edits, direction="+direction);
			final int qoffset=qforward ? 0 : qe, roffset=rforward ? 0 : re;
			final byte[] qc=Arrays.copyOfRange(query, qoffset, qforward ? qe+1 : query.length);
			final byte[] rc=Arrays.copyOfRange(ref, roffset, rforward ? re+1 : ref.length);
			final Probe cropped=new Probe(7);
			final int expected=align(cropped, direction, qc, rc, qs-qoffset, rs-roffset,
				qforward ? qc.length-1 : 0, rforward ? rc.length-1 : 0);
			check(distance==expected && actual.lastRow==cropped.lastRow && actual.lastEdits==cropped.lastEdits
				&& actual.offset()==cropped.offset() && actual.lastQueryLoc==cropped.lastQueryLoc+qoffset
				&& actual.lastRefLoc==cropped.lastRefLoc+roffset,
				"Bounded results must equal a physically cropped terminal tail, direction="+direction);
			if(qforward){Arrays.fill(query, qe+1, query.length, (byte)'N');}
			else{Arrays.fill(query, 0, qe, (byte)'N');}
			check(align(actual, direction, query, ref, qs, rs, qe, re)==expected,
				"Bases beyond the saved stop must not change revalidation");
		}
	}

	private static int align(final BandedAligner b, final int direction, final byte[] q, final byte[] r,
		final int qs, final int rs, final int qe, final int re){
		assert(direction>=0 && direction<4) : "Fixture must select a defined traversal";
		return b.alignFixed(q, r, qs, rs, qe, re, 6, true, direction==1 || direction==3);
	}

	/** Literal terminal-gap oracle: four substitutions plus three extra bases cost seven. */
	private static void fixedEndControls(){
		final byte[] query=new byte[100], ref=new byte[103];
		Arrays.fill(query, (byte)'A');
		Arrays.fill(ref, (byte)'A');
		Arrays.fill(ref, 100, 103, (byte)'G');
		for(int i=10; i<90; i+=20){query[i]='C';}
		final BandedAligner band=BandedAligner.makeBandedAligner(7);
		check(band.alignFixed(query, ref, 0, 0, 99, 102, 6, true, false)>6,
			"A four-edit prefix must not validate the seven-edit complete saved spans");
		check(band.alignFixed(query, ref, 0, 0, 99, 102, 7, true, false)==7,
			"Four substitutions and three terminal gaps must cost seven");
		Arrays.fill(query, (byte)'A');
		check(band.alignFixed(query, ref, 0, 0, 99, 102, 3, true, false)==3,
			"Valid terminal gaps must survive an extension tie favoring a shorter endpoint");
		check(band.lastQueryLoc==99 && band.lastRefLoc==102,
			"Successful fixed validation must consume both saved stops");
		// Exhaustive short binary strings use an independent unrestricted matrix.
		// A 15-column band contains every cell for these lengths.
		final BandedAligner wide=BandedAligner.makeBandedAligner(15);
		for(int n=1; n<=4; n++){
			for(int m=1; m<=4; m++){
				for(int a=0; a<(1<<n); a++){
					for(int b=0; b<(1<<m); b++){
						final byte[] q=binary(n, a), r=binary(m, b);
						final int expected=globalDistance(q, r);
						for(int direction=0; direction<4; direction++){
							final boolean qf=direction==0 || direction==3, rf=direction<2;
							final byte[] qt=q.clone(), rt=r.clone();
							if(!qf){AminoAcid.reverseComplementBasesInPlace(qt);}
							if(!rf){AminoAcid.reverseComplementBasesInPlace(rt);}
							final int actual=align(wide, direction, qt, rt, qf ? 0 : n-1,
								rf ? 0 : m-1, qf ? n-1 : 0, rf ? m-1 : 0);
							check(actual==expected, "Fixed spans must match the independent full-matrix oracle");
						}
					}
				}
			}
		}
	}

	private static byte[] binary(final int length, final int value){
		assert(length>0 && length<=4) : "Finite oracle enumerates one- to four-base strings";
		final byte[] bases=new byte[length];
		for(int i=0; i<length; i++){bases[i]=(byte)((value&(1<<i))==0 ? 'A' : 'C');}
		return bases;
	}

	/** Test-only unrestricted matrix; deliberately separate from the rolling-band implementation. */
	private static int globalDistance(final byte[] q, final byte[] r){
		assert(q.length>0 && r.length>0) : "Oracle inputs are nonempty finite strings";
		final int[][] matrix=new int[q.length+1][r.length+1];
		for(int i=0; i<=q.length; i++){matrix[i][0]=i;}
		for(int j=0; j<=r.length; j++){matrix[0][j]=j;}
		for(int i=1; i<=q.length; i++){
			for(int j=1; j<=r.length; j++){
				matrix[i][j]=Math.min(matrix[i-1][j-1]+(q[i-1]==r[j-1] ? 0 : 1),
					Math.min(matrix[i-1][j]+1, matrix[i][j-1]+1));
			}
		}
		return matrix[q.length][r.length];
	}

	@SuppressWarnings("unchecked")
	private static void workerFailure(final Dedupe2 tool) throws Exception{
		final Read read=new Read(bytes(FIRST_A), null, "intentional_worker_failure", 0);
		final Class<?> unit=Class.forName("jgi.Dedupe2$Unit");
		final Constructor<?> constructor=unit.getDeclaredConstructor(Dedupe2.class, Read.class);
		constructor.setAccessible(true);
		final Object value=constructor.newInstance(tool, read);
		read.obj=value;
		// HashThread.findOverlaps asserts that an allocated overlap list is nonempty.
		// Supply that invalid state deliberately; no production injection flag is needed.
		field(unit, "overlapList").set(value, new ArrayList<Object>());
		final Map<Long, ArrayList<Object>> map=(Map<Long, ArrayList<Object>>)field(Dedupe2.class, "codeMap").get(tool);
		final ArrayList<Object> list=new ArrayList<Object>();
		list.add(value);
		map.put((Long)field(unit, "code1").get(value), list);
		field(Dedupe2.class, "addedToMain").setLong(tool, 1);
		System.err.println("INTENTIONAL_EMPTY_OVERLAP_LIST_IN_HASH_THREAD");
		method(Dedupe2.class, "findOverlaps", Timer.class).invoke(tool, new Timer());
		// This is the old bad path: worker died but main still writes success output.
		method(Dedupe2.class, "writeOutput", String.class, Timer.class).invoke(tool, null, new Timer());
		System.err.println("WORKER_FAILURE_WAS_SILENTLY_ACCEPTED");
	}

	private static Field field(final Class<?> owner, final String name) throws Exception{
		final Field f=owner.getDeclaredField(name);
		f.setAccessible(true);
		return f;
	}
	private static Method method(final Class<?> owner, final String name, final Class<?>... types) throws Exception{
		final Method m=owner.getDeclaredMethod(name, types);
		m.setAccessible(true);
		return m;
	}
	private static byte[] bytes(final String sequence){return sequence.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean valid, final String message){
		if(!valid){throw new AssertionError(message);}
	}

	private static final class Probe extends BandedAlignerConcrete{
		Probe(final int width){super(width);}
		int offset(){return lastOffset;}
	}
	private static final String FIRST_A="GCCTATTGCTGGGGGTCTCTTAATCCTTTATCTAATCAGGGTTATAGAGTGAGAGTTACCCTTTGAACCTGACGCGGATAATCCCGCCGTAGGGAAGCCACTGATAAGCCTC";
	private static final String FIRST_B="GCCTATTGCTGGGGGTCTCTTAATCCTTTATCTAATCAGGGTTATAGAGTGAGAGTTACCCTTTGAACCTGACGCGGATAATCCCGCCGTAGGGAAGCCGATAAGCCTCACA";
	private static final String SECOND_A="TGGATGAAACAGGGGTGCCGTACGGATGGGCCGTGCGGCTGAGAGAGTCCCTTCGCACCCGATCCGGTTCGTACCGGCGTGGGAAGTTTCTACAGAACACCG";
	private static final String SECOND_B="TGGATGAAACAGGGGTGCCGTACGGATGGGCCGTGCGGCTGAGAGAGTCCCTTCGCACCCGATCCGGTTCGTACCGGCGTGGGAAGTTTCAGAACACCGCGT";
}
