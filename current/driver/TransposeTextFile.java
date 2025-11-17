package driver;

import fileIO.ReadWrite;
import fileIO.TextFile;

/**
 * Transposes tab-delimited text files by converting rows to columns.
 * Processes single files or multiple chromosome-numbered files using # placeholder.
 * Optionally skips header lines and outputs transposed data to .transposed files.
 * @author Brian Bushnell
 */
public class TransposeTextFile {
	
	/**
	 * Program entry point for file transposition utility.
	 * Processes files with chromosome numbering (1-22) if # placeholder is present,
	 * otherwise processes single file. Second argument specifies lines to skip.
	 * @param args Command-line arguments: [filename] [optional_skip_lines]
	 */
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
	
	/**
	 * Transposes a tab-delimited text file by converting rows to columns.
	 * Reads entire file into memory, skips specified header lines, then writes
	 * each column as a tab-separated row to .transposed output file.
	 *
	 * @param fname Input filename to transpose
	 * @param skipLines Number of header lines to skip during transposition
	 */
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
