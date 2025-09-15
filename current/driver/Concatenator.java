package driver;

import fileIO.TextFile;
import fileIO.TextStreamWriter;

public class Concatenator {
	
	
	public static void main(String args[]){
		
		assert(args.length==2 && !args[1].contains(","));
		TextStreamWriter tsw=new TextStreamWriter(args[1], false, false, true);
		tsw.start();
		for(String s : args[0].split(",")){
			writeFile(s, tsw);
		}
		tsw.poison();
	}
	
	public static void writeFile(String fname, TextStreamWriter tsw){
		TextFile tf=new TextFile(fname, false);
		if(tsw==null){
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				System.out.println(s);
			}
		}else{
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				tsw.println(s);
			}
		}
		tf.close();
	}
	
	
	public static StringBuilder merge(String[] fnames){
		StringBuilder sb=new StringBuilder();
		
		for(int i=0; i<fnames.length; i++){
			String fname=fnames[i];
			if(fname!=null){
				TextFile tf=new TextFile(fname, false);
				String[] lines=tf.toStringLines();
				tf.close();
				for(int j=0; j<lines.length; j++){
					String s=lines[j];
					lines[j]=null;
//					if(i<2 || !s.startsWith("#")){
//						sb.append(s);
//						sb.append('\n');
//					}
					sb.append(s);
					sb.append('\n');
				}
			}
		}
		return sb;
	}
	
	
}
