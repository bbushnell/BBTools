package driver;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import fileIO.TextFile;

/**
 * Utility for merging two tab-delimited text files based on matching keys.
 * Combines rows from two files by matching values in specified columns,
 * preserving headers and handling missing entries gracefully.
 * @author Brian Bushnell
 */
public class MergeTextFiles2 {
	
	/**
	 * Program entry point for command-line file merging.
	 * Merges two files using first argument as file1, second as file2,
	 * matching on column 0 from file1 and column 1 from file2.
	 * @param args Command-line arguments: file1 file2
	 */
	public static void main(String[] args){
		CharSequence sb=mergeWithHeader(args[0], args[1], 0, 1);
		System.out.println(sb);
	}
	
	/**
	 * Merges two tab-delimited files based on matching keys in specified columns.
	 * Loads both files, creates lookup tables by key column, then combines all
	 * rows with matching keys. Missing entries are handled by creating placeholder
	 * rows with the matching key.
	 *
	 * @param fname1 Path to first input file
	 * @param fname2 Path to second input file
	 * @param col1 Column index in first file to use as key
	 * @param col2 Column index in second file to use as key
	 * @return StringBuilder containing merged tab-delimited output
	 */
	public static StringBuilder mergeWithHeader(String fname1, String fname2, int col1, int col2){

		TextFile tf1=new TextFile(fname1, false);
		String[][] lines1=TextFile.doublesplitTab(tf1.toStringLines(), false);
		tf1.close();
		tf1=null;
		
		TextFile tf2=new TextFile(fname2, false);
		String[][] lines2=TextFile.doublesplitTab(tf2.toStringLines(), false);
		tf2.close();
		tf2=null;

		int maxWidth1=findMaxWidth(lines1);
		int maxWidth2=findMaxWidth(lines2);
		
		Hashtable<String, String[]> table1=makeTable(lines1, col1, 1);
		Hashtable<String, String[]> table2=makeTable(lines2, col2, 1);
		
		HashSet<String> keySet=new HashSet<String>();
		keySet.addAll(table1.keySet());
		keySet.addAll(table2.keySet());
		String[] keys=keySet.toArray(new String[0]);
		Arrays.sort(keys);
		
		StringBuilder sb=new StringBuilder();
		sb.append(toString(lines1[0], lines2[0], maxWidth1, maxWidth2));
		sb.append('\n');
		
		for(String key : keys){
			String[] line1=table1.get(key);
			String[] line2=table2.get(key);
			
			if(line1==null){
				line1=new String[col1+1];
				line1[col1]=line2[col2];
			}
			
			sb.append(toString(line1, line2, maxWidth1, maxWidth2));
			sb.append('\n');
		}
		
		return sb;
	}
	
	/**
	 * Concatenates two string arrays into tab-delimited format.
	 * Outputs elements from array a followed by elements from array b,
	 * padding with tabs to specified widths and handling null arrays/elements.
	 *
	 * @param a First string array (may be null)
	 * @param b Second string array (may be null)
	 * @param alen Maximum width to output from first array
	 * @param blen Maximum width to output from second array
	 * @return StringBuilder with tab-delimited concatenated values
	 */
	private static StringBuilder toString(String[] a, String[] b, int alen, int blen){
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<alen; i++){
			if(a!=null && a.length>i && a[i]!=null){
				sb.append(a[i]);
			}
			sb.append('\t');
		}
		for(int i=0; i<blen; i++){
			if(b!=null && b.length>i && b[i]!=null){
				sb.append(b[i]);
			}
			sb.append('\t');
		}
		return sb;
	}

	/**
	 * Creates lookup table mapping key column values to complete rows.
	 * Builds hashtable using specified column as key, starting from given line
	 * to skip headers. Each key maps to its complete row array.
	 *
	 * @param lines Two-dimensional array of parsed file lines
	 * @param col Column index to use as lookup key
	 * @param firstLine Starting line index (typically 1 to skip header)
	 * @return Hashtable mapping key strings to complete row arrays
	 */
	private static Hashtable<String, String[]> makeTable(String[][] lines, int col, int firstLine) {
		Hashtable<String, String[]> table=new Hashtable<String, String[]>();
		for(int i=firstLine; i<lines.length; i++){
			String[] line=lines[i];
			table.put(line[col], line);
		}
		return table;
	}
	
	/**
	 * Finds the maximum row width in a two-dimensional string matrix.
	 * Iterates through all rows to determine the longest row length,
	 * handling null rows safely.
	 *
	 * @param matrix Two-dimensional string array to analyze
	 * @return Maximum number of columns found in any row
	 */
	private static int findMaxWidth(String[][] matrix){
		int max=0;
		for(String[] line : matrix){
			if(line!=null && max<line.length){max=line.length;}
		}
		return max;
	}
	
}
