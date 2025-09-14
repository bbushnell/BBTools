package driver;

import fileIO.TextFile;
import fileIO.TextStreamWriter;

/** One-off program for converting grch38 sam files to hg19 */
public class FixChr {
	
	public static void main(String[] args){
		
		String in=args[0];
		String out=args[1];
		
		TextFile tf=new TextFile(in);
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, true);
		tsw.start();
		
		String s=null;
		while((s=tf.nextLine())!=null){
			if(!s.startsWith("#")){s="chr"+s;}
			else if(s.startsWith("##contig=<ID=")){
				s="##contig=<ID=chr"+s.substring("##contig=<ID=".length());
			}
			tsw.println(s);
		}
		tf.close();
		tsw.poisonAndWait();
	}
	
}
