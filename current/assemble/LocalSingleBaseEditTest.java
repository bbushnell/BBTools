package assemble;

import java.util.Arrays;
import dna.AminoAcid;
import stream.Read;

/** Independent sequence/quality and reverse-complement tests for LocalSingleBaseEdit. */
public final class LocalSingleBaseEditTest {
	public static void main(final String[] args){
		final boolean oldChangeQuality=Read.CHANGE_QUALITY;
		try{
			Read.CHANGE_QUALITY=false;
			basicOperations();
			nullQuality();
			rcSymmetry();
			rejectedRequests();
			System.out.println("LOCAL_SINGLE_BASE_EDIT_TEST_OK");
		}finally{
			Read.CHANGE_QUALITY=oldChangeQuality;
		}
	}

	private static void basicOperations(){
		final byte[] source="ACGTACG".getBytes(), quality=quality(source.length);
		final int mid=source.length/2;
		for(final int p:new int[]{0,mid,source.length-1}){
			verify(source,quality,LocalSingleBaseEdit.Operation.SUBSTITUTION,p,(byte)'T');
			verify(source,quality,LocalSingleBaseEdit.Operation.DELETION,p,(byte)0);
		}
		for(final int p:new int[]{0,mid,source.length-1,source.length}){
			verify(source,quality,LocalSingleBaseEdit.Operation.INSERTION,p,(byte)'G');
		}
		verify(new byte[0],new byte[0],LocalSingleBaseEdit.Operation.INSERTION,0,(byte)'A');
		verify(new byte[]{'G'},new byte[]{93},LocalSingleBaseEdit.Operation.DELETION,0,(byte)0);
	}

	private static void nullQuality(){
		final byte[] source="GATTACA".getBytes();
		for(final LocalSingleBaseEdit.Operation op:LocalSingleBaseEdit.Operation.values()){
			final int p=op==LocalSingleBaseEdit.Operation.INSERTION ? source.length : source.length/2;
			verify(source,null,op,p,(byte)'C');
		}
	}

	private static void rcSymmetry(){
		final byte[] source="ACGTTGCAC".getBytes(), quality=quality(source.length);
		final int n=source.length;
		for(final LocalSingleBaseEdit.Operation op:LocalSingleBaseEdit.Operation.values()){
			for(final int p:op==LocalSingleBaseEdit.Operation.INSERTION ? new int[]{0,n/2,n} : new int[]{0,n/2,n-1}){
				final byte base=op==LocalSingleBaseEdit.Operation.DELETION ? 0 : (byte)'C';
				final Read forward=new Read(source.clone(),quality.clone(),"forward",0,false);
				final byte[] rcBases=AminoAcid.reverseComplementBases(source), rcQuality=reverse(quality);
				final Read rcRead=new Read(rcBases,rcQuality,"reverse",0,false);
				final LocalSingleBaseEdit f=new LocalSingleBaseEdit(), r=new LocalSingleBaseEdit();
				if(op==LocalSingleBaseEdit.Operation.DELETION){
					f.apply(forward,op,p);r.apply(rcRead,op,n-1-p);
				}else if(op==LocalSingleBaseEdit.Operation.INSERTION){
					f.apply(forward,op,p,base);r.apply(rcRead,op,n-p,complement(base));
				}else{
					f.apply(forward,op,p,base);r.apply(rcRead,op,n-1-p,complement(base));
				}
				check(Arrays.equals(forward.bases,AminoAcid.reverseComplementBases(rcRead.bases)),"RC sequence mismatch for "+op+" at "+p);
				check(Arrays.equals(forward.quality,reverse(rcRead.quality)),"RC quality mismatch for "+op+" at "+p);
			}
		}
	}

	private static void rejectedRequests(){
		final byte[] source="ACGTACG".getBytes(), quality=quality(source.length);
		final LocalSingleBaseEdit.Operation[] ops={LocalSingleBaseEdit.Operation.SUBSTITUTION,LocalSingleBaseEdit.Operation.DELETION,LocalSingleBaseEdit.Operation.INSERTION};
		for(final LocalSingleBaseEdit.Operation op:ops){
			final int[] positions=op==LocalSingleBaseEdit.Operation.INSERTION ? new int[]{-1,source.length+1} : new int[]{-1,source.length};
			for(final int p:positions){
				final Read read=new Read(source.clone(),quality.clone(),"invalid",0,false);final byte[] oldBases=read.bases,oldQuality=read.quality;
				boolean threw=false;try{apply(read,op,p,(byte)'A');}catch(IllegalArgumentException expected){threw=true;}
				check(threw && read.bases==oldBases && read.quality==oldQuality && Arrays.equals(read.bases,source) && Arrays.equals(read.quality,quality),"Invalid position mutated read.");
			}
		}
		final Read baseRead=new Read(source.clone(),quality.clone(),"invalid-base",0,false);
		for(final LocalSingleBaseEdit.Operation op:new LocalSingleBaseEdit.Operation[]{LocalSingleBaseEdit.Operation.SUBSTITUTION,LocalSingleBaseEdit.Operation.INSERTION}){
			final byte[] oldBases=baseRead.bases,oldQuality=baseRead.quality;boolean threw=false;
			try{apply(baseRead,op,2,(byte)'N');}catch(IllegalArgumentException expected){threw=true;}
			check(threw && baseRead.bases==oldBases && baseRead.quality==oldQuality,"Invalid base mutated read.");
		}
		boolean threw=false;try{new LocalSingleBaseEdit().apply(baseRead,LocalSingleBaseEdit.Operation.DELETION,2,(byte)'A');}catch(IllegalArgumentException expected){threw=true;}
		check(threw,"Deletion-with-base overload must reject the base argument.");
		threw=false;try{new LocalSingleBaseEdit().apply(baseRead,LocalSingleBaseEdit.Operation.SUBSTITUTION,2);}catch(IllegalArgumentException expected){threw=true;}
		check(threw,"No-base overload must reject non-deletion operations.");
		threw=false;try{new LocalSingleBaseEdit().apply(baseRead,null,2,(byte)'A');}catch(IllegalArgumentException expected){threw=true;}
		check(threw,"Null operation must fail loudly.");
		final Read metadata=new Read(source.clone(),quality.clone(),"metadata",0,false);metadata.setMapped(true);
		final byte[] oldBases=metadata.bases,oldQuality=metadata.quality;threw=false;
		try{new LocalSingleBaseEdit().apply(metadata,LocalSingleBaseEdit.Operation.INSERTION,2,(byte)'A');}catch(IllegalArgumentException expected){threw=true;}
		check(threw && metadata.bases==oldBases && metadata.quality==oldQuality && Arrays.equals(metadata.bases,source),"Metadata rejection mutated read.");
		final Read nullBases=new Read(source.clone(),quality.clone(),"null-bases",0,false);
		final byte[] nullQualityRef=nullBases.quality;nullBases.bases=null;threw=false;
		try{new LocalSingleBaseEdit().apply(nullBases,LocalSingleBaseEdit.Operation.DELETION,0);}catch(IllegalArgumentException expected){threw=true;}
		check(threw && nullBases.bases==null && nullBases.quality==nullQualityRef,"Null bases rejection mutated read.");
		threw=false;try{new LocalSingleBaseEdit().apply(null,LocalSingleBaseEdit.Operation.DELETION,0);}catch(IllegalArgumentException expected){threw=true;}
		check(threw,"Null read must fail loudly.");
		final Read mismatch=new Read(source.clone(),quality.clone(),"quality-mismatch",0,false);
		final byte[] mismatchBases=mismatch.bases,malformedQuality=new byte[source.length-1];mismatch.quality=malformedQuality;threw=false;
		try{new LocalSingleBaseEdit().apply(mismatch,LocalSingleBaseEdit.Operation.DELETION,2);}catch(IllegalArgumentException expected){threw=true;}
		check(threw && mismatch.bases==mismatchBases && mismatch.quality==malformedQuality,"Mismatched qualities must fail atomically.");
	}

	private static void verify(final byte[] source,final byte[] quality,final LocalSingleBaseEdit.Operation op,
			final int position,final byte base){
		final Read read=new Read(source.clone(),quality==null ? null : quality.clone(),"basic",0,false);
		final byte[] oldBases=read.bases,oldQuality=read.quality;
		if(op==LocalSingleBaseEdit.Operation.DELETION){new LocalSingleBaseEdit().apply(read,op,position);}
		else{new LocalSingleBaseEdit().apply(read,op,position,base);}
		final byte[][] expected=expected(source,quality,op,position,base);
		check(Arrays.equals(read.bases,expected[0]),"Independent base output mismatch for "+op+" at "+position);
		check(Arrays.equals(read.quality,expected[1]),"Independent quality output mismatch for "+op+" at "+position);
		check(Arrays.equals(oldBases,source) && (oldQuality==null ? quality==null : Arrays.equals(oldQuality,quality)),"Input arrays were mutated.");
		check((read.quality==null)==(quality==null),"Null-quality policy changed.");
	}

	private static byte[][] expected(final byte[] source,final byte[] quality,final LocalSingleBaseEdit.Operation op,
			final int position,final byte base){
		final int delta=op==LocalSingleBaseEdit.Operation.INSERTION ? 1 : op==LocalSingleBaseEdit.Operation.DELETION ? -1 : 0;
		final byte[] b=new byte[source.length+delta],q=quality==null ? null : new byte[b.length];
		if(op==LocalSingleBaseEdit.Operation.SUBSTITUTION){
			for(int i=0;i<source.length;i++){b[i]=i==position ? base : source[i];if(q!=null){q[i]=i==position ? 0 : quality[i];}}
		}else if(op==LocalSingleBaseEdit.Operation.INSERTION){
			for(int i=0;i<b.length;i++){if(i==position){b[i]=base;if(q!=null){q[i]=0;}}else{final int s=i<position ? i : i-1;b[i]=source[s];if(q!=null){q[i]=quality[s];}}}
		}else{
			for(int i=0;i<b.length;i++){final int s=i<position ? i : i+1;b[i]=source[s];if(q!=null){q[i]=quality[s];}}
		}
		return new byte[][]{b,q};
	}

	private static void apply(final Read read,final LocalSingleBaseEdit.Operation op,final int position,final byte base){
		if(op==LocalSingleBaseEdit.Operation.DELETION){new LocalSingleBaseEdit().apply(read,op,position);}
		else{new LocalSingleBaseEdit().apply(read,op,position,base);}
	}
	private static byte[] quality(final int length){final byte[] q=new byte[length];final int[] raw={0,1,60,93};for(int i=0;i<length;i++){q[i]=(byte)raw[i%raw.length];}return q;}
	private static byte[] reverse(final byte[] input){if(input==null){return null;}final byte[] out=new byte[input.length];for(int i=0;i<input.length;i++){out[i]=input[input.length-1-i];}return out;}
	private static byte complement(final byte base){return AminoAcid.numberToBase[AminoAcid.numberToComplement[AminoAcid.baseToNumber[base]]];}
	private static void check(final boolean ok,final String message){if(!ok){throw new AssertionError(message);}}
}
