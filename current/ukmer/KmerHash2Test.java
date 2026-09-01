package ukmer;

import java.util.HashSet;
import java.util.Random;

import dna.AminoAcid;
import shared.Tools;

/** Regression and distribution checks for Tools.splitMix64 and Kmer.xor2. */
public class KmerHash2Test {

	public static void main(String[] args){
		testSplitMix();
		testKmers();
		testMaskedTerminalFingerprint();
		testResizableTable();
		testFixedTable();
		testOwnerTombstones(true);
		testOwnerTombstones(false);
		System.out.println("KmerHash2Test PASS: "+checks+" checks");
	}

	private static void testSplitMix(){
		check(Tools.splitMix64(0)==0, "SplitMix finalizer must map zero to zero");
		check(Tools.splitMix64(1)==0x5692161D100B05E5L, "Unexpected SplitMix64(1) vector");
		check(Tools.splitMix64(-1L)==0xB4D055FCF2CBBD7BL, "Unexpected SplitMix64(-1) vector");
	}

	private static void testKmers(){
		final Random randy=new Random(1);
		final int[] ks={31, 32, 63, 64, 65, 95, 128, 300};
		for(int k : ks){
			final int samples=(k<100 ? 20000 : 5000);
			final HashSet<Long> seen=new HashSet<Long>(samples*2);
			long ones=0;
			for(int sample=0; sample<samples; sample++){
				final byte[] bases=randomBases(k, randy);
				final Kmer a=makeKmer(bases);
				final Kmer b=makeKmer(reverseComplement(bases));
				final long hash=a.xor2();
				check(hash==b.xor2(), "Reverse complements disagreed at k="+k);
				check(hash==Kmer.xor2(a.key()), "Static and cached xor2 disagreed at k="+k);
				check(seen.add(hash), "xor2 collision at k="+k+", sample="+sample);
				ones+=Long.bitCount(hash);

				final int pos=sample%k;
				final long old=AminoAcid.baseToNumber[bases[pos]];
				final long alternate=(old+1)&3;
				a.substituteBase(pos, alternate);
				bases[pos]=AminoAcid.numberToBase[(int)alternate];
				check(a.xor2()==makeKmer(bases).xor2(), "Lazy xor2 was stale after substitution at k="+k);
				a.substituteBase(pos, old);
				check(a.xor2()==hash, "xor2 was not restored at k="+k);
			}
			final double fraction=ones/(samples*63.0);
			check(fraction>0.49 && fraction<0.51, "xor2 bit balance outside 49-51% at k="+k+": "+fraction);
		}
	}

	private static void testResizableTable(){
		final HashArrayH1D table=new HashArrayH1D(new int[] {101, 211, 431});
		final Random randy=new Random(2);
		final Kmer[] kmers=new Kmer[180];
		for(int i=0; i<kmers.length; i++){
			kmers[i]=makeKmer(randomBases(128, randy));
			check(table.incrementAndReturnNumCreated(kmers[i])==1, "Resizable table rejected a new kmer");
			check(table.incrementAndReturnNumCreated(kmers[i])==0, "Resizable table recreated an existing kmer");
			check(table.getValue(kmers[i])==2, "Resizable table count mismatch");
		}
		check(table.arrayLength()==211, "Resizable table did not follow its schedule");
		check(table.size()==kmers.length, "Resizable table size mismatch");
		for(Kmer kmer : kmers){check(table.getValue(kmer)==2, "Resize lost a kmer");}

		table.initializeOwnership();
		check(table.getOwner(kmers[0])==-1, "New ownership was not clear");
		check(table.setOwner(kmers[0], 3)==3, "Ownership assignment failed");
		check(table.getOwner(kmers[0])==3, "Ownership lookup failed");
		check(table.clearOwner(kmers[0], 3), "Ownership clear failed");
		table.clearOwnership();

		check(table.regenerate(1)==0, "Regeneration removed count-2 kmers");
		check(table.regenerate(2)==kmers.length, "Regeneration failed to remove count-2 kmers");
		check(table.size()==0, "Resizable table was nonempty after regeneration");
	}

	private static void testMaskedTerminalFingerprint(){
		final boolean old=Kmer.MASK_CORE;
		Kmer.MASK_CORE=true;
		try{
			final byte[] a=new byte[95], b=new byte[95];
			java.util.Arrays.fill(a, (byte)'C');
			java.util.Arrays.fill(b, (byte)'C');
			b[94]='A';
			final Kmer ka=makeKmer(a), kb=makeKmer(b);
			check(ka.xor()==kb.xor(), "Terminal variants did not share the masked placement hash");
			check(ka.xor2()!=kb.xor2(), "Full-key fingerprint failed to distinguish terminal variants");
		}finally{
			Kmer.MASK_CORE=old;
		}
	}

	private static void testFixedTable(){
		final HashArrayH1D table=new HashArrayH1D(new int[] {101}, false);
		final Random randy=new Random(3);
		final Kmer[] retained=new Kmer[40];
		for(int i=0; i<80; i++){
			Kmer kmer=makeKmer(randomBases(300, randy));
			table.increment(kmer);
			if((i&1)==0){table.increment(kmer); retained[i/2]=kmer;}
		}
		check(table.fingerprintOnly(), "Fixed table unexpectedly stored placement hashes");
		check(table.regenerate(1)==40, "Fixed-table tombstone regeneration removed the wrong count");
		for(Kmer kmer : retained){check(table.getValue(kmer)==2, "Fixed-table regeneration lost a retained kmer");}
		for(int i=0; i<40; i++){
			Kmer kmer=makeKmer(randomBases(300, randy));
			check(table.incrementAndReturnNumCreated(kmer)==1, "Fixed table failed to reuse a tombstone");
		}
		check(table.size()==80, "Fixed table size mismatch after tombstone reuse");
		boolean failedLoudly=false;
		try{
			while(true){table.increment(makeKmer(randomBases(300, randy)));}
		}catch(IllegalStateException e){
			failedLoudly=e.getMessage().contains("cannot resize");
		}
		check(failedLoudly, "Fixed table did not fail loudly when capacity was exhausted");
	}

	private static void testOwnerTombstones(final boolean storePlacementHash){
		final HashArrayH1D table=new HashArrayH1D(new int[] {211}, storePlacementHash);
		final Random randy=new Random(storePlacementHash ? 4 : 5);
		final Kmer[] kmers=new Kmer[120];
		for(int i=0; i<kmers.length; i++){
			kmers[i]=makeKmer(randomBases(161, randy));
			check(table.incrementAndReturnNumCreated(kmers[i])==1, "Tombstone fixture rejected a new kmer");
		}
		table.initializeOwnership();
		int expected=0;
		for(int i=0; i<kmers.length; i++){
			if(i%3==0){table.setOwner(kmers[i], 2); expected++;}
			else{table.setOwner(kmers[i], 3);}
		}
		table.setOwner(kmers[0], 3);
		expected--;
		check(table.removeByOwner(2)==expected, "Ownership sweep removed the wrong number of kmers");
		table.clearOwnership();
		check(table.size()==kmers.length-expected, "Ownership sweep left the wrong live size");
		for(int i=0; i<kmers.length; i++){
			check(table.getValue(kmers[i])==(i%3==0 && i>0 ? -1 : 1), "Tombstone broke a probe-chain lookup");
		}
		for(int i=0; i<20; i++){
			check(table.incrementAndReturnNumCreated(makeKmer(randomBases(161, randy)))==1,
					"Insertion failed to reuse an ownership tombstone");
		}
	}

	private static byte[] randomBases(int length, Random randy){
		byte[] bases=new byte[length];
		for(int i=0; i<length; i++){bases[i]=AminoAcid.numberToBase[randy.nextInt(4)];}
		return bases;
	}

	private static byte[] reverseComplement(byte[] bases){
		byte[] rc=new byte[bases.length];
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			rc[i]=AminoAcid.numberToBase[AminoAcid.baseToComplementNumber[bases[j]]];
		}
		return rc;
	}

	private static Kmer makeKmer(byte[] bases){
		Kmer kmer=new Kmer(bases.length);
		for(byte b : bases){kmer.addRight(b);}
		return kmer;
	}

	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new RuntimeException(message);}
	}

	private static long checks=0;
}
