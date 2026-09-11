package assemble;

import java.util.Arrays;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import ukmer.HashArrayU1D;
import ukmer.Kmer;

/** Controlled count dictionaries isolate verification geometry, not biological accuracy.
 * Each supported word can be supplied as an independent K-base count-library read.
 * @author Fischl */
public final class LocalEditFootprintTest {

	public static void main(final String[] args){
		if(args.length!=1 || !(args[0].equals("observe") || args[0].equals("verify"))){
			throw new IllegalArgumentException("Usage: observe|verify");
		}
		final boolean packed=Kmer.PACKED,core=Kmer.MASK_CORE,changeQuality=Read.CHANGE_QUALITY;
		int failures=0;
		try{
			Kmer.PACKED=true;Kmer.MASK_CORE=false;Read.CHANGE_QUALITY=false;
			System.out.println("k\toperation\treverse\tfull_support\tapplied\texpected\taffected_windows\tlost_supported_windows\tcontract_pass");
			for(final int k:new int[]{31,62,63}){
				for(int operation=0;operation<3;operation++){
					for(final boolean full:new boolean[]{false,true}){
						for(final boolean reverse:new boolean[]{false,true}){if(!test(k,operation,full,reverse)){failures++;}}
					}
				}
			}
		}finally{Kmer.PACKED=packed;Kmer.MASK_CORE=core;Read.CHANGE_QUALITY=changeQuality;}
		System.out.println("FOOTPRINT_CONTRACT_FAILURES="+failures);
		if(args[0].equals("verify")){check(failures==0,"Full edit-footprint contract failed in "+failures+" fixtures; unsupported outside-trough words must not be introduced.");}
	}

	private static boolean test(final int k,final int operation,final boolean full,final boolean reverse){
		final int a=k+7,p=a+k/2;
		final byte[] original=new byte[4*k+30];final Random random=new Random(7000+k);
		for(int i=0;i<original.length;i++){original[i]=ALPHABET[random.nextInt(4)];}
		original[0]='A';original[original.length-1]='A';
		original[p-1]='A';original[p]='C';original[p+1]='G';
		// A...A is lexicographically before its reverse complement T...T, so
		// the test knows the corrector canonical orientation without copying its code.
		final int delta=operation==1 ? 1 : operation==2 ? -1 : 0;
		final byte[] changed=new byte[original.length+delta];
		System.arraycopy(original,0,changed,0,p);
		if(operation==0){changed[p]='T';System.arraycopy(original,p+1,changed,p+1,original.length-p-1);}
		else if(operation==1){changed[p]='T';System.arraycopy(original,p,changed,p+1,original.length-p);}
		else{System.arraycopy(original,p+1,changed,p,original.length-p-1);}
		final Counts counts=new Counts(k);
		for(int start=0;start<=original.length-k;start++){counts.addUnique(original,start,start==a ? 1 : 40);}
		final int affectedLast=operation==2 ? p-1 : p;
		final int first=full ? p-k+1 : a-1,last=full ? affectedLast : a+1+delta;
		for(int start=first;start<=last;start++){counts.addUnique(changed,start,40);}
		check(counts.depth(original,a)==1 && counts.depth(original,a-1)==40 && counts.depth(original,a+1)==40,
			"Explicit dictionary must have one low original window between supported neighbors.");
		int unsupported=0,lost=0;
		for(int start=p-k+1;start<=affectedLast;start++){
			final int before=counts.depth(original,start),after=counts.depth(changed,start);
			if(after<3){unsupported++;if(before>=3){lost++;}}
		}
		final int affected=affectedLast-(p-k+1)+1;
		check(unsupported==(full ? 0 : affected-(3+delta)) && lost==unsupported,
			"Only old trough/flank candidate words are supplied in the negative fixture; all other altered words were supported before.");
		final byte[] bases=reverse ? AminoAcid.reverseComplementBases(original) : original.clone();
		final byte[] qualities=new byte[bases.length];Arrays.fill(qualities,(byte)40);
		final Read read=new Read(bases,qualities,"footprint",0,false);
		final int applied=new LocalEditCorrector(k,counts).correctOne(read,false);
		check(applied==0 || applied==1,"One-call corrector contract must return zero or one applied edit.");
		if(applied==1){
			check(Arrays.equals(read.bases,reverse ? AminoAcid.reverseComplementBases(changed) : changed),
				"The sole supported candidate must be the explicitly constructed edit, not another event.");
			final byte[] expectedQ=new byte[changed.length];Arrays.fill(expectedQ,(byte)40);
			if(operation!=2){expectedQ[reverse ? changed.length-1-p : p]=0;}
			check(Arrays.equals(read.quality,expectedQ),"Applied edit must splice retained qualities and set only its inserted/substituted base to Q0.");
		}else{check(read.bases==bases && read.quality==qualities,"Abstention must preserve both original array identities.");}
		final int expected=full ? 1 : 0;final boolean pass=applied==expected;
		System.out.println(k+"\t"+operation+"\t"+reverse+"\t"+full+"\t"+applied+"\t"+expected+"\t"+affected+"\t"+(applied==1 ? lost : 0)+"\t"+pass);
		return pass;
	}

	private static final class Counts implements HomopolymerIndelProposal.CountLookup{
		Counts(final int k_){k=k_;key=new Kmer(k);table=new HashArrayU1D(new int[]{2003},key.k,k);}
		void setWord(final byte[] bases,final int start){
			assert(start>=0 && start+k<=bases.length) : "Fixed test word must contain exactly K existing bases.";
			key.clearFast();for(int i=start;i<start+k;i++){key.addRight(bases[i]);}
		}
		void addUnique(final byte[] bases,final int start,final int copies){
			setWord(bases,start);check(table.getValue(key)<1,"Fixture unexpectedly repeats a canonical word; original/alternate count assignments would be ambiguous.");
			for(int i=0;i<copies;i++){table.increment(key);}
		}
		int depth(final byte[] bases,final int start){setWord(bases,start);return Math.max(0,table.getValue(key));}
		@Override public int count(final Kmer word){return table.getValue(word);}
		final int k;final Kmer key;final HashArrayU1D table;
	}
	private static void check(final boolean ok,final String message){if(!ok){throw new AssertionError(message);}}
	private static final byte[] ALPHABET={'A','C','G','T'};
}
