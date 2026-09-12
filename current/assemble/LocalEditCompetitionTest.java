package assemble;

import java.util.Arrays;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import ukmer.HashArrayU1D;
import ukmer.Kmer;

/** Exact count dictionaries isolate cross-family ambiguity, not biological accuracy.
 * The legacy S-first policy and opt-in competitor check see identical counts.
 * @author Fischl */
public final class LocalEditCompetitionTest {
	public static void main(final String[] args){
		final boolean packed=Kmer.PACKED,core=Kmer.MASK_CORE,quality=Read.CHANGE_QUALITY;
		int cases=0;
		try{
			Kmer.PACKED=true;Kmer.MASK_CORE=false;Read.CHANGE_QUALITY=false;
			for(final int k:new int[]{31,62,63}){for(final int windows:new int[]{1,3}){
				for(final int delta:new int[]{-1,1}){for(int support=0;support<3;support++){
					for(final boolean reverse:new boolean[]{false,true}){test(k,windows,delta,support,reverse);cases++;}
				}}
			}}
		}finally{Kmer.PACKED=packed;Kmer.MASK_CORE=core;Read.CHANGE_QUALITY=quality;}
		System.out.println("LOCAL_EDIT_COMPETITION_OK cases="+cases);
	}
	private static void test(final int k,final int windows,final int delta,final int support,final boolean reverse){
		final int a=k+7,s=a+3,p=a+k/2;
		final byte[] original=new byte[4*k+30];final Random random=new Random(9100+k);
		for(int i=0;i<original.length;i++){original[i]=ALPHABET[random.nextInt(4)];}
		original[0]=original[original.length-1]='A';original[s]='C';
		original[p-1]='A';original[p]='C';original[p+1]='G';
		final byte[] substitution=original.clone();substitution[s]='T';
		final byte[] indel=new byte[original.length+delta];
		System.arraycopy(original,0,indel,0,p);
		if(delta>0){indel[p]='T';System.arraycopy(original,p,indel,p+1,original.length-p);}
		else{System.arraycopy(original,p+1,indel,p,original.length-p-1);}
		final Counts counts=new Counts(k);
		for(int start=0;start<=original.length-k;start++){counts.add(original,start,start==a ? 1 : 40);}
		//support0=S only, support1=indel only, support2=both distinct outcomes.
		if(support!=1){for(int start=s-k+1;start<=s;start++){counts.add(substitution,start,40);}}
		if(support!=0){for(int start=p-k+1;start<=p-(delta<0 ? 1 : 0);start++){counts.add(indel,start,40);}}
		check(counts.depth(original,a)==1,"Fixture must retain its one low original word despite adding alternate outcomes.");
		for(final boolean guard:new boolean[]{false,true}){
			final byte[] bases=reverse ? AminoAcid.reverseComplementBases(original) : original.clone();
			final byte[] qualities=new byte[bases.length];Arrays.fill(qualities,(byte)40);
			final Read read=new Read(bases,qualities,"competition",0,false);
			final LocalEditCorrector caller=new LocalEditCorrector(k,counts,windows,guard);
			final int applied=caller.correctOne(read,false);
			final boolean ambiguous=guard && support==2;
			check(applied==(ambiguous ? 0 : 1),"Competing verified edit families must abstain only when requested; k="+k+", support="+support+", guard="+guard);
			if(ambiguous){
				check(caller.ambiguousLoci==1 && read.bases==bases && read.quality==qualities,
					"Ambiguous locus must retain both original arrays, not install an intermediate substitution.");
			}else{
				final byte[] expected=support==1 ? indel : substitution;
				check(Arrays.equals(read.bases,reverse ? AminoAcid.reverseComplementBases(expected) : expected),
					"Unique supported family must still repair; legacy policy must retain S-first behavior.");
			}
		}
	}
	private static final class Counts implements HomopolymerIndelProposal.CountLookup {
		Counts(final int k_){k=k_;key=new Kmer(k);table=new HashArrayU1D(new int[]{2003},key.k,k);}
		void word(final byte[] bases,final int start){
			assert(start>=0 && start+k<=bases.length) : "Dictionary additions must use complete explicitly constructed candidate words.";
			key.clearFast();for(int i=start;i<start+k;i++){key.addRight(bases[i]);}
		}
		void add(final byte[] bases,final int start,final int copies){
			word(bases,start);
			final int existing=table.getValue(key);
			if(existing>0){check(existing==copies,"Alternate support must not overwrite the deliberately low original word.");return;}
			for(int i=0;i<copies;i++){table.increment(key);}
		}
		int depth(final byte[] bases,final int start){word(bases,start);return table.getValue(key);}
		@Override public int count(final Kmer word){return table.getValue(word);}
		final int k;final Kmer key;final HashArrayU1D table;
	}
	private static void check(final boolean ok,final String message){if(!ok){throw new AssertionError(message);}}
	private static final byte[] ALPHABET={'A','C','G','T'};
}
