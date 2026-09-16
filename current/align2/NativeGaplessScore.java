package align2;

/**
 * Allocation-free score of the operations produced by native genMatchNoIndels.
 * Unweighted operation-score contract; not the historical fast score or a claim
 * that every DP boundary/unknown-base convention is identical.
 * Used by a private experimental native hook; no mapper qualification.
 * @author Collei
 */
public final class NativeGaplessScore {
	private NativeGaplessScore(){}
	public static int score(MSA msa,byte[] query,byte[] reference,int start){
		if(query==null || reference==null || query.length==0){return 0;}
		byte mode=0,previous=0;int length=0,previousLength=0,total=0;
		for(int i=0;i<query.length;i++){
			long position=(long)start+i;byte r=position<0 || position>=reference.length ? (byte)'N' : reference[(int)position];
			byte q=query[i],op=q=='N' || r=='N' ? (byte)'N' : q==r ? (byte)'m' : (byte)'S';
			if(op==mode){length++;}
			else{
				if(length>0){total+=run(msa,mode,length,previous,previousLength);}
				previous=mode;previousLength=length;mode=op;length=1;
			}
		}
		return total+run(msa,mode,length,previous,previousLength);
	}
	private static int run(MSA msa,byte mode,int length,byte previous,int previousLength){
		if(mode=='m'){return msa.calcMatchScore(length);}
		if(mode=='N'){return msa.calcNocallScore(length);}
		assert(mode=='S');int score=msa.calcSubScore(length);
		if(previous=='N'){score+=msa.POINTS_SUB2()-msa.POINTS_SUB();}
		else if(previous=='m' && previousLength==1){score+=msa.POINTS_SUBR()-msa.POINTS_SUB();}
		return score;
	}
}
