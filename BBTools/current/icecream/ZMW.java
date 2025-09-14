package icecream;

import java.util.ArrayList;

import shared.Tools;
import stream.Read;
import stream.SamLine;
import structures.IntList;

/**
 * Container for the list of reads from a single
 * PacBio ZMW.
 * @author Brian Bushnell
 * @date June 5, 2020
 */
public class ZMW extends ArrayList<Read> {
	
	/**
	 * For serialization.
	 */
	private static final long serialVersionUID = -2580124131008824113L;

	public ZMW(){super();}
	
	public ZMW(int initialSize){super(initialSize);}

	public long countBases(){
		long x=0;
		for(Read r : this){
			x+=r.length();
		}
		return x;
	}
	
	public int medianLength(boolean includeDiscarded){
		if(size()<3){return -1;}
		IntList lengths=new IntList(size()-2);
		
		for(int i=1; i<size()-1; i++){
			Read r=get(i);
			if(includeDiscarded || !r.discarded()){
				lengths.add(get(i).length());
			}
		}
		lengths.sort();
		int median=lengths.get(lengths.size/2);
		return median;
	}
	
	public int longestLength(boolean includeDiscarded){
		int max=0;
		for(Read r : this){
			if(includeDiscarded || !r.discarded()){
				max=Tools.max(max, r.length());
			}
		}
		return max;
	}
	
	public Read medianRead(boolean includeDiscarded){
		int len=medianLength(includeDiscarded);
		if(len<0){return longestRead(includeDiscarded);}
		for(int i=1; i<size()-1; i++){
			Read r=get(i);
			if((includeDiscarded || !r.discarded()) && r.length()==len){
				return r;
			}
		}
		return null;
	}
	
	public Read longestRead(boolean includeDiscarded){
		Read max=null;
		for(Read r : this){
			if((includeDiscarded || !r.discarded()) && (max==null || r.length()>max.length())){max=r;}
		}
		return max;
	}
	
	public int zid(){
		if(zid==-1){parseZID();}
		return zid;
	}
	
	private int parseZID(){
		return (size()<1 ? -1 : PBHeader.parseZMW(get(0).id));
	}
	
	public static void fixReadHeader(Read r, int leftTrim, int rightTrim){
		leftTrim=Tools.max(0, leftTrim);
		rightTrim=Tools.max(0, rightTrim);
		if(leftTrim<1 && rightTrim<1){return;}
		final int idx=r.id.lastIndexOf('/');
		if(idx>0 && idx<r.id.length()-3){
			String prefix=r.id.substring(0, idx+1);
			String suffix=r.id.substring(idx+1);
			if(suffix.indexOf('_')>0){
				String coords=suffix, comment="";
				int tab=suffix.indexOf('\t');
				if(tab<0){tab=suffix.indexOf(' ');}
				if(tab>0){
					coords=coords.substring(0, tab);
					comment=coords.substring(tab);
				}
				String[] split=Tools.underscorePattern.split(coords);
				int left=Integer.parseInt(split[0]);
				int right=Integer.parseInt(split[1]);
				left+=leftTrim;
				right-=rightTrim;
				if(left>right){left=right;}
				
				if(right-left!=r.length()){right=left+r.length();}
//				System.err.println(r.length()+", "+(right-left));
				
				r.id=prefix+left+"_"+right+comment;
				final SamLine sl=r.samline;
				if(sl!=null){
					sl.qname=r.id;
					if(sl.optional!=null){
						for(int i=0; i<sl.optional.size(); i++){
							String s=sl.optional.get(i);
							if(s.startsWith("qe:i:")){
								s="qe:i:"+right;
								sl.optional.set(i, s);
							}else if(s.startsWith("qs:i:")){
								s="qs:i:"+left;
								sl.optional.set(i, s);
							}
						}
					}
				}
			}
		}
	}
	
	public void setDiscarded(boolean b){
		for(Read r : this){
			r.setDiscarded(b);
		}
	}

	public int[] lengths() {
		final int size=size();
		int[] array=new int[size];
		for(int i=0; i<size; i++){
			Read r=get(i);
			array[i]=r==null ? -1 : r.length();
		}
		return array;
	}
	
	public float estimatePasses(){
		final int size=size();
		if(size<1){return 0;}
		else if(size==1){return 0.25f;}
		else if(size==2){return 0.5f;}
		
		int median=medianLength(true);
		int first=first().length();
		int last=last().length();

		return size-2+estimatePasses(first, median)+estimatePasses(last, median);
	}
	
	private float estimatePasses(int len, int median){
		float ratio=len/(float)median;
		//TODO: I want this to be more asymptotic
		return Tools.min(0.99f, ratio/(1+0.05f*ratio));
	}

	public boolean discarded() {
		for(Read r : this){
			if(!r.discarded()){return false;}
		}
		return true;
	}
	
	/** 
	 * Identifier assigned by streamer, not by PacBio.
	 * First identifier is 0, then 1, etc.
	 */
	public long id;
	
	/** 
	 * ZMW ID assigned by PacBio.
	 */
	private int zid=-1;

	public Read first(){return get(0);}
	public Read last(){return get(size()-1);}
	
}
