package cardinality;

import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import stream.Read;

/** Deterministic sketch, merge and histogram regressions for khistmulti.
 * @author Brian Bushnell, Fischl
 */
public class KmerHistMultiTest{

	/** Runs scalar tests or writes deterministic paired FASTQ fixtures. */
	public static void main(final String[] args){
		if(args.length==1){writeFixtures(args[0]); return;}
		assert(args.length==0) : "Use no arguments for tests, or one existing fixture directory.";
		testMerge();
		testSaturation();
		testEmpty();
		testPairsAndStrands();
		testParser();
		System.out.println("KHISTMULTI_UNIT_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Tests            ----------------*/
	/*--------------------------------------------------------------*/

	/** Exercises small, exact packed, multiword and longer-than-read K values. */
	private static void testMerge(){
		final int[] ks={31, 32, 33, 64, 95, 96, 124, 127, 300, 601};
		final Random random=new Random(631);
		for(final int buckets : new int[] {64, 2048}){
			final DynamicDemiLog[] serial=KmerHistMulti.createSketches(ks, buckets, 7, 0);
			final DynamicDemiLog[][] pieces=new DynamicDemiLog[4][];
			for(int i=0; i<pieces.length; i++){pieces[i]=KmerHistMulti.createSketches(ks, buckets, 7, 0);}
			for(int i=0; i<1200; i++){
				final Read r=new Read(sequence(random, 400), null, i);
				if(i%7==0){r.bases[200]='N';}
				final int copies=1+i%5;
				for(int c=0; c<copies; c++){
					for(int j=0; j<ks.length; j++){
						serial[j].hash(r);
						pieces[(i+c)%4][j].hash(r);
					}
				}
			}
			for(int i=1; i<pieces.length; i++){KmerHistMulti.mergeSketches(pieces[0], pieces[i]);}
			for(int i=0; i<ks.length; i++){
				check(Arrays.equals(serial[i].maxArray(), pieces[0][i].maxArray()), "Maximum-score merge differs at K="+ks[i]);
				check(Arrays.equals(serial[i].counts16(), pieces[0][i].counts16()), "Depth-count merge differs at K="+ks[i]);
				check(serial[i].cardinality()==pieces[0][i].cardinality(), "Cardinality merge differs at K="+ks[i]);
				check(Arrays.equals(KmerHistMulti.histogram(serial[i], 100),
					KmerHistMulti.histogram(pieces[0][i], 100)), "Histogram merge differs at K="+ks[i]);
			}
		}
	}

	/** Saturated and pooled-tail observations must retain exactly one sampled winner. */
	private static void testSaturation(){
		final DynamicDemiLog a=DynamicDemiLog.create(64, 32, 3, 0, true);
		final DynamicDemiLog b=DynamicDemiLog.create(64, 32, 3, 0, true);
		for(int i=0; i<40000; i++){a.add(15437); b.add(15437);}
		a.add(b);
		check(KmerHistMulti.histogram(a, 65535)[65535]==1, "Merged counter must saturate, not wrap.");
		check(KmerHistMulti.histogram(a, 100)[100]==1, "Tail aggregation lost a saturated bucket.");
		for(int i=0; i<100; i++){a.add(15437);}
		check(KmerHistMulti.histogram(a, 65535)[65535]==1, "Incrementing a saturated counter wrapped.");
	}

	/** Empty, short and all-N sequences contribute no canonical kmers. */
	private static void testEmpty(){
		for(final int k : new int[] {31, 32, 64, 124, 300}){
			final DynamicDemiLog d=DynamicDemiLog.create(64, k, 0, 0, true);
			check(d.cardinality()==0, "Empty cardinality must be zero.");
			final byte[] ns=new byte[400];
			Arrays.fill(ns, (byte)'N');
			d.hash(new Read(ns, null, 1));
			d.hash(new Read(new byte[] {'A', 'C', 'G'}, null, 2));
			check(d.filledBuckets()==0, "Short/all-N sequences created kmers at K="+k);
			check(Arrays.equals(new long[11], KmerHistMulti.histogram(d, 10)), "Empty histogram is nonzero.");
		}
	}

	/** One hash of a paired root equals hashing each mate once; reverse complements agree. */
	private static void testPairsAndStrands(){
		final Random random=new Random(55);
		for(final int k : new int[] {31, 32, 95, 124, 300}){
			final DynamicDemiLog pair=DynamicDemiLog.create(2048, k, 4, 0, true);
			final DynamicDemiLog singles=DynamicDemiLog.create(2048, k, 4, 0, true);
			final DynamicDemiLog reverse=DynamicDemiLog.create(2048, k, 4, 0, true);
			final Read a=new Read(sequence(random, 420), null, 1);
			final Read b=new Read(sequence(random, 420), null, 2);
			singles.hash(a); singles.hash(b);
			a.mate=b; b.mate=a;
			pair.hash(a);
			AminoAcid.reverseComplementBasesInPlace(a.bases);
			AminoAcid.reverseComplementBasesInPlace(b.bases);
			reverse.hash(a);
			check(Arrays.equals(pair.counts16(), singles.counts16()), "Mate counted twice or omitted at K="+k);
			check(Arrays.equals(pair.maxArray(), reverse.maxArray()), "Reverse-complement hash differs at K="+k);
			check(Arrays.equals(pair.counts16(), reverse.counts16()), "Reverse-complement depths differ at K="+k);
		}
	}

	/** The K parser rejects ambiguity instead of silently dropping/rounding columns. */
	private static void testParser(){
		check(Arrays.equals(KmerHistMulti.parseK("124,33,64"), new int[] {124, 33, 64}), "Requested column order changed.");
		for(final String bad : new String[] {"", "0", "-2", "32,32", "32,", "x"}){
			boolean rejected=false;
			try{KmerHistMulti.parseK(bad);}
			catch(IllegalArgumentException e){rejected=true;}
			check(rejected, "Invalid K list accepted: "+bad);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes equivalent interleaved and split inputs with shared sequence at known multiplicities. */
	private static void writeFixtures(final String directory){
		final ByteStreamWriter paired=new ByteStreamWriter(directory+"/paired.fq", false, false, true);
		final ByteStreamWriter r1=new ByteStreamWriter(directory+"/r1.fq", false, false, true);
		final ByteStreamWriter r2=new ByteStreamWriter(directory+"/r2.fq", false, false, true);
		paired.start(); r1.start(); r2.start();
		final Random random=new Random(92);
		final byte[] quals=new byte[150];
		Arrays.fill(quals, (byte)'I');
		for(int i=0; i<200; i++){
			final byte[] a=sequence(random, 150), b=sequence(random, 150);
			for(int c=0; c<1+i%8; c++){
				for(int mate=1; mate<=2; mate++){
					final String id="@pair"+i+"_"+c+"/"+mate;
					final byte[] bases=mate==1 ? a : b;
					final ByteStreamWriter single=mate==1 ? r1 : r2;
					paired.print(id).nl().print(bases).print("\n+\n").print(quals).nl();
					single.print(id).nl().print(bases).print("\n+\n").print(quals).nl();
				}
			}
		}
		check(!(paired.poisonAndWait() | r1.poisonAndWait() | r2.poisonAndWait()), "Fixture output failed.");
	}

	/** Generates repeatable fully defined DNA without per-base objects. */
	private static byte[] sequence(final Random random, final int length){
		assert(length>0) : "Test sequences need at least one base.";
		final byte[] bases=new byte[length];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}

	/** Counts a checked invariant and fails even when Java assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks=0;
}
