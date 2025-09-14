package tracker;

import dna.AminoAcid;
import fileIO.ReadWrite;
import shared.Tools;
import stream.Read;
import structures.LongList;

/** Tracks counts of base types per position to make a consensus sequence.
 * Designed for use with BBMerge since adapters are inferred by insert size;
 * could be used with alignment too. */
public class AdapterTracker {
	
	public AdapterTracker() {
		for(int i=0; i<counts.length; i++){
			for(int j=0; j<counts[i].length; j++){
				counts[i][j]=new LongList(151);
			}
		}
	}
	
	public void storeAdapterSequence(Read r, int insert){
		reads++;
		if(r.length()<=insert) {return;}
		shortInserts++;
		
		if(looksLikePhix(r, insert)){
			phixLike++;
			if(ignorePhixAdapters){return;}
		}
		
		LongList[] lists=counts[r.pairnum()];
		byte[] bases=r.bases;
		
		for(int i=insert, j=0; i<bases.length; i++, j++){
			byte b=bases[i];
			int num=AminoAcid.baseToNumber[b];
			if(num>=0){
				lists[num].increment(j);
			}
		}
	}
	
	private boolean looksLikePhix(Read r, int insert){
		return looksLikePhix(r.bases, insert) || looksLikePhix(r.mate.bases, insert);
	}
	
	private boolean looksLikePhix(byte[] bases, int insert){
		int len=bases.length-insert;
		if(len<phixPrefix.length){return false;}
		for(int i=insert, j=0; i<bases.length && j<phixPrefix.length; i++, j++){
			byte b=bases[i];
			if(b!='N' && b!=phixPrefix[j]){
				return false;
			}
		}
//		outstream.println(new String(bases).substring(insert));
//		outstream.println(new String(phixPrefix));
		return true;
	}
	
	public boolean makeSequence() {
		seq1=seq2=null;
		seq1=toAdapterSequence(counts[0], trimPolyAorG);
		seq2=toAdapterSequence(counts[1], trimPolyAorG);
		return hasSequence();
	}
	
	public boolean hasSequence() {
		return (seq1!=null && seq1.length()>1) || (seq2!=null && seq2.length()>1);
	}
	
	public long writeAdapterConsensus(String fname){
		StringBuilder sb=new StringBuilder();
		{
			sb.append(">Read1_adapter\n");
			String adapter=toAdapterSequence(counts[0], trimPolyAorG);
			sb.append(adapter).append('\n');
		}
		if(counts.length>1){
			sb.append(">Read2_adapter\n");
			String adapter=toAdapterSequence(counts[1], trimPolyAorG);
			sb.append(adapter).append('\n');
		}
		long count=counts[0][0].get(0)+counts[0][1].get(0)+
				counts[0][2].get(0)+counts[0][3].get(0);
//		outstream.println("Adapters counted: "+count);
		ReadWrite.writeString(sb, fname);
		return count;
	}
	
	private static String toAdapterSequence(LongList[] lists, boolean trimPolyAorG){
		StringBuilder adapter=new StringBuilder();
		long max=0;
		int lastBase=-1;
		for(int i=0; true; i++){
			long a=lists[0].get(i);
			long c=lists[1].get(i);
			long g=lists[2].get(i);
			long t=lists[3].get(i);
			long sum=(a+c+g+t);
			max=Tools.max(max, sum);
			if(sum==0 || (sum<10 && sum<=max/1000) || (max>100 && sum<8)){break;}
			long thresh=(max>100 ? 4+(sum*2)/3 : (sum*2)/3);
			if(a>thresh){
				adapter.append('A');
				lastBase=i;
			}else if(c>thresh){
				adapter.append('C');
				lastBase=i;
			}else if(g>thresh){
				adapter.append('G');
				lastBase=i;
			}else if(t>thresh){
				adapter.append('T');
				lastBase=i;
			}else{
				adapter.append('N');
			}
		}
		if(lastBase<0){return "N";}

		String trimmed=trimPoly2(adapter.toString(), 'N');
		if(trimPolyAorG){
			for(int len=-1; len!=trimmed.length(); ) {
				len=trimmed.length();
				trimmed=trimPoly2(trimmed, 'G');
				trimmed=trimPoly2(trimmed, 'A');
			}
		}
		if(trimJunk){
			trimmed=trimJunk(trimmed, 6);
		}
		
//		if(lastBase>=0){
//			char A=(trimPolyAorG ? 'A' : 'N');
//			while(lastBase>=0 && (adapter.charAt(lastBase)=='N' || adapter.charAt(lastBase)==A)){lastBase--;}
//		}
		
		if(trimmed.length()<1){return "N";}
		return trimmed;
	}
	
	private static String trimPoly(String adapter, char trim){
		int lastBase=-1;
		for(int i=0; i<adapter.length(); i++){
			char c=adapter.charAt(i);
			if(AminoAcid.isFullyDefined(c)){
				lastBase=i;
			}
		}
		
		int aCount=0;
		int nCount=0;
		int count=0;
		while(lastBase>=0){
			char c=adapter.charAt(lastBase);
			if(c=='N'){nCount++;}
			else if(c==trim){aCount++;}
			else{break;}
			count++;
			lastBase--;
		}
		
		if(lastBase<0){return "N";}
		if(count==nCount || (aCount>3)){
			return adapter.substring(0, lastBase+1);
		}
		return adapter;
	}
	
	private static String trimPoly2(String adapter, char poly){
		int last=adapter.length()-1;
		int trim=0;
		while(last>=0) {
			char c=adapter.charAt(last);
//			System.err.println("c="+Character.toString(c)+", poly="+Character.toString(poly));
			if(c==poly || c=='N') {
				trim++;
				last--;
			}else{
				break;
			}
		}
		
		if(trim>3 || (trim>0 && poly=='N')) {
			adapter=adapter.substring(0, last+1);
		}
//		assert(poly=='N') : Character.toString(poly)+"\n"+adapter+"\n"+trim+"\n"+last;
		return adapter==null || adapter.length()<1 ? "N" : adapter;
	}
	
	private static String trimJunk(String s, int minScore) {
		int score=0, last=s.length()-1;
		for(; last>=0 && score<minScore; last--) {
			char c=s.charAt(last);
			if(c=='N') {
				score--;
			}else {
				score+=2;
			}
		}
		last++;
		while(last<s.length() && s.charAt(last)!='N' || (last<s.length()-1 && s.charAt(last+1)!='N')) {last++;}
		return (last<1 ? "N" : last>s.length() ? s : s.substring(0, last));
	}
	
	public void merge(AdapterTracker b){
		for(int x=0; x<counts.length; x++){
			for(int y=0; y<counts[x].length; y++){
				counts[x][y].incrementBy(b.counts[x][y]);
			}
		}
		reads+=b.reads;
		shortInserts+=b.shortInserts;
		phixLike+=b.phixLike;
	}
	
	final LongList[][] counts=new LongList[2][4];
	public String seq1=null;
	public String seq2=null;
	public long reads=0;
	public long shortInserts=0;
	public long phixLike=0;
	
	private static final byte[] phixPrefix="AGATCGGAAGAGCG".getBytes();
	public static boolean ignorePhixAdapters=false;
	public static boolean trimPolyAorG=true;
	public static boolean trimJunk=true;
	
}
