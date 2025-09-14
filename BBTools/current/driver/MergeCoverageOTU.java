package driver;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import fileIO.TextFile;
import fileIO.TextStreamWriter;
import jgi.CovStatsLine;

/**
 * @author Brian Bushnell
 * @date Oct 13, 2014
 *
 */
public class MergeCoverageOTU {
	
	public static void main(String[] args){
		String a=args[0];
		String b=args[1];
		String in=null, out=null;
		if(a.toLowerCase().startsWith("in=")){
			in=a.split("=")[1];
			out=b.split("=")[1];
		}else if(a.toLowerCase().startsWith("out=")){
			in=b.split("=")[1];
			out=a.split("=")[1];
		}else{
			in=a;
			out=b;
		}
		
		TextFile tf=new TextFile(in);
		LinkedHashMap<String, CovStatsLine> map=new LinkedHashMap<String, CovStatsLine>();
		int count=0;
		ArrayList<String> headers=new ArrayList<String>();
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(count==0){
				assert(s.startsWith("#")) : "Expected a header line starting with #";
				CovStatsLine.initializeHeader(s);
			}else{
				int space=s.indexOf(' ');
				String otu=s.substring(space+1, s.indexOf('\t'));
				CovStatsLine csl=new CovStatsLine(s);
				CovStatsLine old=map.get(otu);
				if(old==null){
					map.put(otu, csl);
				}else{
					old.add(csl);
				}
			}
			count++;
		}
		tf.close();
		
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		for(String s : headers){tsw.println(s);}
		for(String s : map.keySet()){
			CovStatsLine csl=map.get(s);
			csl.id=s;
			tsw.println(csl.toString());
		}
		tsw.poisonAndWait();
	}
	
}
