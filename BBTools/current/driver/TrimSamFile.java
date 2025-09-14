package driver;

import java.util.HashSet;

import fileIO.TextFile;
import shared.Tools;
import stream.SamLine;

public class TrimSamFile {
	
	public static void main(String[] args){
		String fname=args[0];
		String scaf=args[1];
		int from=Integer.parseInt(args[2]);
		int to=Integer.parseInt(args[3]);
		TextFile tf=new TextFile(fname, false);
		HashSet<String> set=findBadLines(tf, scaf, from, to);
		tf.reset();
		printExcludingSet(tf, set);
	}
	
	
	public static HashSet<String> findBadLines(TextFile tf, String scafS, int from, int to){
		byte[] scaf=scafS.getBytes();
		HashSet<String> set=new HashSet<String>(16000);
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)!='@'){//header
				SamLine sl=new SamLine(s);
				
				if(sl.pos>=from && sl.pos<=to && Tools.equals(sl.rname(), scaf)){
					set.add(sl.qname);
				}else if(sl.pnext>=from && sl.pnext<=to && Tools.equals(sl.rnext(), scaf)){
					set.add(sl.qname);
				}else if(Tools.equals(sl.rname(), scaf) && Tools.equals(sl.rnext(), scaf) && (sl.pos<from != sl.pnext<from)){
					set.add(sl.qname);
				}else if(!sl.mapped() || !sl.nextMapped() || !sl.pairedOnSameChrom()){
					set.add(sl.qname);
				}
			}
		}
		return set;
	}
	
	
	public static void printExcludingSet(TextFile tf, HashSet<String> set){
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)=='@'){//header
				System.out.println(s);
			}else{
				SamLine sl=new SamLine(s);
				
				if(!set.contains(sl.qname)){
					System.out.println(s);
				}
			}
		}
	}
	
	
}

