package driver;

import fileIO.TextFile;

/**
 * Summarizes match/sub/ins/del/N rates for consecutive BBMap runs
 * @author Brian Bushnell
 * @date Jan 8, 2014
 *
 */
public class SummarizeMSDIN {
	
	public static void main(String[] args){
		String fname=args[0];
		boolean M=false;
		boolean E=false;
		boolean S=true;
		boolean D=false;
		boolean I=false;
		boolean N=false;
		boolean B=false;
		boolean MS=true;

		long mcount=0;
		long ecount=0;
		long scount=0;
		long dcount=0;
		long icount=0;
		long ncount=0;
		long bcount=0;
		
		TextFile tf=new TextFile(fname);
		StringBuilder sb=new StringBuilder();
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			String[] split=s.split("\t");
			if(s.startsWith("Total time:")){
				if(B){
					if(sb.length()>0){sb.append('\t');}
					sb.append(bcount);
				}
				if(MS){
					if(sb.length()>0){sb.append('\t');}
					sb.append((mcount+scount));
				}
				if(M){
					if(sb.length()>0){sb.append('\t');}
					sb.append(mcount);
				}
				if(E){
					if(sb.length()>0){sb.append('\t');}
					sb.append(ecount);
				}
				if(S){
					if(sb.length()>0){sb.append('\t');}
					sb.append(scount);
				}
				if(D){
					if(sb.length()>0){sb.append('\t');}
					sb.append(dcount);
				}
				if(I){
					if(sb.length()>0){sb.append('\t');}
					sb.append(icount);
				}
				if(N){
					if(sb.length()>0){sb.append('\t');}
					sb.append(ncount);
				}
				System.out.println(sb);
				sb.setLength(0);
				mcount=ecount=scount=dcount=icount=ncount=bcount=0;
			}else if(s.startsWith("Match Rate:")){
				String x=split[split.length-1];
				try{mcount=(Long.parseLong(x));}catch(Exception e){}
			}else if(E && s.startsWith("Error Rate:")){
				String x=split[split.length-1];
//				if(E){
//					if(sb.length()>0){sb.append('\t');}
//					sb.append(x);
//				}
				try{ecount=(Long.parseLong(x));}catch(Exception e){}
			}else if(s.startsWith("Sub Rate:")){
				String x=split[split.length-1];
//				if(S){
//					if(sb.length()>0){sb.append('\t');}
//					sb.append(x);
//				}
				try{scount=(Long.parseLong(x));}catch(Exception e){}
			}else if(s.startsWith("Del Rate:")){
				String x=split[split.length-1];
//				if(D){
//					if(sb.length()>0){sb.append('\t');}
//					sb.append(x);
//				}
				try{dcount=(Long.parseLong(x));}catch(Exception e){}
			}else if(s.startsWith("Ins Rate:")){
				String x=split[split.length-1];
//				if(I){
//					if(sb.length()>0){sb.append('\t');}
//					sb.append(x);
//				}
				try{icount=(Long.parseLong(x));}catch(Exception e){}
			}else if(s.startsWith("N Rate:")){
				String x=split[split.length-1];
//				if(N){
//					if(sb.length()>0){sb.append('\t');}
//					sb.append(x);
//				}
				try{ncount=(Long.parseLong(x));}catch(Exception e){}
			}else if(s.startsWith("Reads Used:")){
				String x=split[split.length-1].replace("(", "").replace(" bases)", "");
//				if(B){
//					if(sb.length()>0){sb.append('\t');}
//					sb.append(x);
//				}
				try{bcount=(Long.parseLong(x));}catch(Exception e){}
			}
		}
		
	}
	
}
