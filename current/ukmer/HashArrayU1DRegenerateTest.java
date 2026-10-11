package ukmer;

import java.util.Random;

/** Regression checks for explicit-table deletion accounting and surviving keys.
 * @author Brian Bushnell, Fischl */
public class HashArrayU1DRegenerateTest {

	/** Exercises pre-zeroed entries, collisions, count thresholds and reinsertion. */
	public static void main(String[] args){
		boolean enabled=false;
		assert(enabled=true);
		if(!enabled){throw new IllegalStateException("Run regeneration regression with -ea.");}
		final boolean packed=Kmer.PACKED, masked=Kmer.MASK_CORE;
		try{
			Kmer.PACKED=true;
			Kmer.MASK_CORE=false;
			testZeroedMain();
			testNegativeKey();
			for(boolean layout : new boolean[] {false, true}){
				Kmer.PACKED=layout;
				for(boolean mask : new boolean[] {false, true}){
					Kmer.MASK_CORE=mask;
					for(int k : new int[] {32, 62, 95, 128}){
						for(int limit : new int[] {0, 2}){
							testMixed(k, limit);
						}
					}
				}
			}
		}finally{
			Kmer.PACKED=packed;
			Kmer.MASK_CORE=masked;
		}
		System.out.println("REGENERATE_ACCOUNTING_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Tests            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reproduces a main-array deletion that formerly escaped both counters. */
	private static void testZeroedMain(){
		final Kmer key=randomKmer(62, new Random(1));
		final HashArrayU1D table=newTable(key, 101);
		table.increment(key);
		check(table.size()==1 && table.victims().size()==0, "Single key must occupy the main array.");
		final int cell=table.findKmer(key);
		check(cell>=0, "Fixture key must have a real array cell before deletion.");
		table.values()[cell]=0;
		check(table.regenerate(0)==1, "Regeneration missed a pre-zeroed main-array kmer.");
		check(totalSize(table)==0, "A deleted main-array entry must leave no stored size.");
		check(table.regenerate(0)==0, "Repeated empty regeneration must not recount deletions.");
		table.increment(key);
		check(totalSize(table)==1 && table.getValue(key)==1, "Deleted cells must support fresh insertion.");
	}

	/** Guards against reviving key-bit occupancy tests for a valid all-one word. */
	private static void testNegativeKey(){
		final Kmer key=new Kmer(64);
		for(int i=0; i<64; i++){key.addRightNumeric(i<32 ? 3 : 0);}
		check(key.key()[0]==-1L, "Fixture must contain a legitimate all-one packed word.");
		final HashArrayU1D table=newTable(key, 101);
		table.set(key, 3);
		check(table.regenerate(0)==0 && table.getValue(key)==3, "All-one packed words must survive regeneration.");
		table.values()[table.findKmer(key)]=0;
		check(table.regenerate(0)==1 && totalSize(table)==0, "All-one keys must also be removable by count.");
	}

	/** Forces main-array and victim deletions, then verifies a complete oracle. */
	private static void testMixed(final int k, final int limit){
		final Random random=new Random(1000+k+limit);
		final Kmer[] keys=new Kmer[200];
		final int[] expected=new int[keys.length];
		keys[0]=randomKmer(k, random);
		final HashArrayU1D table=newTable(keys[0], 31);
		for(int i=0; i<keys.length; i++){
			if(i>0){keys[i]=randomKmer(k, random);}
			check(table.getValue(keys[i])<1, "Deterministic fixture unexpectedly repeated a canonical kmer.");
			expected[i]=1+i%5;
			table.set(keys[i], expected[i]);
		}
		check(totalSize(table)==keys.length && table.victims().size()>0, "Fixture must exercise real collision-cache entries.");
		table.initializeOwnership();
		for(int i=0; i<keys.length; i++){
			if(i%4==0){table.setOwner(keys[i], 2); expected[i]=0;}
		}
		int zeroedMain=0, zeroedVictims=0;
		for(int cell=0; cell<table.arrayLength(); cell++){
			if(table.owners.get(cell)==2){table.values()[cell]=0; zeroedMain++;}
		}
		for(KmerNodeU node : table.victims().toList()){
			if(node.owner()==2){node.set(0); zeroedVictims++;}
		}
		check(zeroedMain>0 && zeroedVictims>0, "Fixture must zero entries in both storage paths.");
		table.clearOwnership();
		long removed=0;
		for(int i=0; i<expected.length; i++){
			if(expected[i]<=limit){expected[i]=0; removed++;}
		}
		check(table.regenerate(limit)==removed, "Removal total must include pre-zeroed and thresholded entries.");
		verify(table, keys, expected);
		check(table.regenerate(0)==0, "A clean table must report no removals on a second regeneration.");
		verify(table, keys, expected);
		for(int i=0; i<keys.length; i++){
			if(expected[i]==0){table.increment(keys[i]); expected[i]=1;}
		}
		verify(table, keys, expected);
		check(table.regenerate(Integer.MAX_VALUE)==keys.length, "Removing everything must count both storage paths once.");
		check(totalSize(table)==0 && !table.walk().next(), "An emptied table must have neither entries nor stale size.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks every expected count, the stored total and an independent walk. */
	private static void verify(final HashArrayU1D table, final Kmer[] keys, final int[] expected){
		check(keys.length==expected.length, "Each fixture key needs one expected count.");
		long live=0;
		for(int i=0; i<keys.length; i++){
			final int actual=table.getValue(keys[i]);
			check(expected[i]==0 ? actual<=0 : actual==expected[i], "Surviving count or deleted-key lookup changed at index "+i+".");
			if(expected[i]>0){live++;}
		}
		check(totalSize(table)==live, "Stored size must equal the actual surviving key count.");
		long walked=0;
		final HashArrayU1D.WalkerU1D walker=table.walk();
		while(walker.next()){
			check(walker.value()>0, "Regenerated walks must never expose deleted entries.");
			walked++;
		}
		check(walked==live, "A full walk must agree with the independent key oracle.");
	}

	/** Uses the actual normalized Kmer geometry in either packing mode. */
	private static HashArrayU1D newTable(final Kmer key, final int prime){
		assert(prime>0 && key.len()>=key.kbig) : "Test tables require a complete key and positive initial size.";
		return new HashArrayU1D(new int[] {prime}, key.k, key.kbig);
	}

	/** Generates a full kmer without depending on the requested length's rounding. */
	private static Kmer randomKmer(final int k, final Random random){
		assert(k>0 && random!=null) : "Fixture generation needs a positive length and deterministic RNG.";
		final Kmer key=new Kmer(k);
		for(int i=0; i<key.kbig; i++){key.addRightNumeric(random.nextInt(4));}
		return key;
	}

	/** Includes both physical stores, as loader capacity accounting does. */
	private static long totalSize(final HashArrayU1D table){return table.size()+table.victims().size();}

	/** Counts checked facts and supplies a diagnostic assertion to negative controls. */
	private static void check(final boolean condition, final String message){
		assert(condition) : message;
		checks++;
	}

	private static long checks=0;
}
