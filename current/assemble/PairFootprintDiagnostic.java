package assemble;

import java.util.Arrays;
import java.util.Random;
import dna.AminoAcid;
import stream.Read;
import structures.IntList;
import ukmer.HashArrayU1D;
import ukmer.Kmer;

/** Isolated controlled-count diagnostic; not product code or biological evidence.
 * @author Fischl */
public final class PairFootprintDiagnostic {

	public static void main(String[] args){
		Kmer.PACKED=true;Kmer.MASK_CORE=false;Read.CHANGE_QUALITY=false;
		System.out.println("k\tfull_support\tpair_found\tfirst_position\tsecond_position\tbelow_threshold_final_words\tordinary_applied");
		for(int k:new int[]{31,62}){for(boolean full:new boolean[]{false,true}){test(k,full);}}
		geometry();
	}
	private static void test(final int k,final boolean full){
		final int a=k+7,p=a+k-1,q=a+1;
		final byte[] original=new byte[4*k+30];final Random random=new Random(17000+k);
		for(int i=0;i<original.length;i++){original[i]=BASES[random.nextInt(4)];}
		original[0]='A';original[original.length-1]='A';original[p]=original[q]='C';
		final byte[] first=original.clone();first[p]='T';
		final byte[] witness=first.clone();witness[q]='G';
		final Counts counts=new Counts(k);
		for(int s=0;s<=original.length-k;s++){counts.assign(original,s,s==a ? 1 : 40);}
		// The second verifier also checks its supporting right flank at a+2.
		for(int s=a;s<=p;s++){counts.assign(first,s,s==a+1 ? 0 : s<=a+2 || full ? 40 : 3);}
		for(int s=q-k+1;s<=q;s++){counts.assign(witness,s,40);}
		final IntList profile=new IntList();
		for(int s=0;s<=original.length-k;s++){profile.add(counts.depth(original,s));}
		require(profile.get(a)==1 && profile.get(a-1)==40 && profile.get(a+1)==40,"Original dictionary must have one width1 trough, not an invented profile.");
		final LocalEditPairLookahead pair=new LocalEditPairLookahead(k,counts);
		final boolean found=pair.propose(original,profile);
		require(found==full,"Full pair-footprint verification must reject weak first-only words and retain the fully supported witness.");
		if(found){
			require(pair.firstOperation==LocalSingleBaseEdit.Operation.SUBSTITUTION && pair.firstPosition==p && pair.firstBase=='T',"Actual first candidate must match the constructed substitution.");
			require(pair.secondOperation==LocalSingleBaseEdit.Operation.SUBSTITUTION && pair.secondPositionAfterFirst==q && Arrays.equals(pair.witness,witness),"Real second corrector must choose the proposed q and exact two-substitution endpoint.");
		}
		int below=0;
		for(int s=a;s<=p;s++){if(counts.depth(witness,s)<5){below++;}}
		require(below==(full ? 0 : k-3),"Final first-only words must be exactly the count3 words outside the second footprint and flank; two substitutions have no coordinate shift.");
		final Read read=new Read(original.clone(),null,"pair-footprint",0,false);
		final int ordinary=new LocalEditCorrector(k,counts).correctOne(read,false);
		require(ordinary==0 && Arrays.equals(read.bases,original),"Pair fixture must not be satisfiable by an ordinary single edit.");
		final LocalEditCorrector integrated=new LocalEditCorrector(k,counts);
		final int applied=integrated.correctOne(read,true);
		require(applied==(full ? 1 : 0) && integrated.usedPairLookahead==full && Arrays.equals(read.bases,full ? first : original),"Integrated caller must reject weak pairs and preserve supported first-edit application.");
		final byte[] rc=AminoAcid.reverseComplementBases(original);
		final byte[] qualities=new byte[rc.length];Arrays.fill(qualities,(byte)40);
		final Read reverse=new Read(rc,qualities,"pair-footprint-rc",0,false);
		final int rcApplied=new LocalEditCorrector(k,counts).correctOne(reverse,true);
		require(rcApplied==applied && Arrays.equals(reverse.bases,AminoAcid.reverseComplementBases(full ? first : original)),"Pair footprint decisions and output must commute with reverse complement.");
		if(!full){require(reverse.bases==rc && reverse.quality==qualities,"Pair abstention must preserve original sequence and quality array identity.");}
		else{for(int i=0;i<reverse.length();i++){require(reverse.quality[i]==(i==reverse.length()-1-p ? 0 : 40),"Accepted reverse pair must change only the first substituted base quality to Q0.");}}
		System.out.println(k+"\t"+full+"\t"+found+"\t"+pair.firstPosition+"\t"+pair.secondPositionAfterFirst+"\t"+below+"\t"+ordinary);
	}
	private static void geometry(){
		int cases=0,changed=0;
		final LocalSingleBaseEdit.Operation[] ops={LocalSingleBaseEdit.Operation.SUBSTITUTION,LocalSingleBaseEdit.Operation.INSERTION,LocalSingleBaseEdit.Operation.DELETION};
		for(int k:new int[]{5,31,62}){
			final int n=2*k+6;final int[] original=new int[n];for(int i=0;i<n;i++){original[i]=i;}
			for(LocalSingleBaseEdit.Operation op1:ops){for(int p:new int[]{0,1,k-1,k,n-2,n-1,n}){
				if(p==n && op1!=LocalSingleBaseEdit.Operation.INSERTION){continue;}
				final int[] first=edited(original,op1,p,-1);
				for(LocalSingleBaseEdit.Operation op2:ops){for(int q:new int[]{0,1,p-k,p-1,p,p+1,p+k,first.length-1,first.length}){
					if(q<0 || q>first.length || (q==first.length && op2!=LocalSingleBaseEdit.Operation.INSERTION)){continue;}
					final int[] last=edited(first,op2,q,-2);
					final long bounds=LocalEditPairLookahead.verificationBounds(k,last.length,op1,p,op2,q,0,k);
					final int from=(int)(bounds>>>32),to=(int)bounds;
					for(int s=0;s+k<=last.length;s++){
						boolean untouched=last[s]>=0;
						for(int j=0;j<k;j++){untouched&=last[s+j]>=0 && last[s+j]==last[s]+j;}
						if(!untouched){changed++;require(s>=from && s+k<=to,"Every materialized altered window must be inside final bounds: K="+k+" first="+op1+"@"+p+" second="+op2+"@"+q+" start="+s+" bounds="+from+".."+to);}
					}
					cases++;
				}}
			}}
		}
		require(cases>1000 && changed>1000,"Geometry panel must exercise both orders of all S/I/D edits, boundary equality, clipping and separated footprints.");
		System.out.println("PAIR_GEOMETRY_PASS cases="+cases+" changed_windows="+changed);
	}
	private static int[] edited(int[] in,LocalSingleBaseEdit.Operation op,int p,int label){
		assert(p>=0 && p<=in.length && (p<in.length || op==LocalSingleBaseEdit.Operation.INSERTION)) : "Origin-label edit must address a valid base or insertion boundary.";
		final int delta=op==LocalSingleBaseEdit.Operation.INSERTION ? 1 : op==LocalSingleBaseEdit.Operation.DELETION ? -1 : 0;
		final int[] out=new int[in.length+delta];System.arraycopy(in,0,out,0,p);
		if(delta<0){System.arraycopy(in,p+1,out,p,in.length-p-1);}
		else{out[p]=label;System.arraycopy(in,p+(delta==0 ? 1 : 0),out,p+1,in.length-p-(delta==0 ? 1 : 0));}
		return out;
	}
	private static final class Counts implements HomopolymerIndelProposal.CountLookup {
		Counts(int k_){k=k_;key=new Kmer(k);table=new HashArrayU1D(new int[]{2003},key.k,k);}
		void word(byte[] b,int start){
			assert(start>=0 && start+k<=b.length) : "Diagnostic count assignments must address complete existing K-base words.";
			key.clearFast();for(int i=start;i<start+k;i++){key.addRight(b[i]);}
		}
		void assign(byte[] b,int s,int copies){
			word(b,s);require(table.getValue(key)<1,"Unexpected canonical key collision invalidates controlled independent word counts.");
			for(int i=0;i<copies;i++){table.increment(key);}
		}
		int depth(byte[] b,int s){word(b,s);return Math.max(0,table.getValue(key));}
		@Override public int count(Kmer word){return table.getValue(word);}
		final int k;final Kmer key;final HashArrayU1D table;
	}
	private static void require(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
	private static final byte[] BASES={'A','C','G','T'};
}
