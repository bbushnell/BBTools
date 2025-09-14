package driver;

import fileIO.ReadWrite;
import fileIO.TextFile;

public class TransposeTextFile {
	
	public static void main(String[] args){
		
		int skipLines=args.length>1 ? Integer.parseInt(args[1]) : 0;
		
		int minChrom=1;
		int maxChrom=22;
		
		for(int i=minChrom; i<=maxChrom; i++){
			if(args[0].contains("#")){
				process(args[0].replace("#", ""+i), skipLines);
			}else{
				process(args[0], skipLines);
				break;
			}
		}
		
	}
	
	public static void process(String fname, int skipLines){
		TextFile tf=new TextFile(fname, false);
		String[] lines=tf.toStringLines();
		tf.close();
		String[][] lines2=TextFile.doublesplitWhitespace(lines, true);
		
		StringBuilder sb=new StringBuilder(4096);
		
		int columns=lines2[skipLines].length;

		for(int column=0; column<columns; column++){
			String tab="";
			for(int row=skipLines; row<lines.length; row++){
				sb.append(tab);
				sb.append(lines2[row][column]);
				tab="\t";
			}
			sb.append("\n");
		}
		
		ReadWrite.writeString(sb, fname+".transposed");
		
	}
	
	
	
	
}
