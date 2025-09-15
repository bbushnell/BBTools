package driver;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import fileIO.TextFile;

public class MergeTextFiles {
	
	public static void main(String[] args){
		CharSequence sb=mergeWithHeader(args[0], args[1], 0, 1);
		System.out.println(sb);
	}
	
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

	private static Hashtable<String, String[]> makeTable(String[][] lines, int col, int firstLine) {
		Hashtable<String, String[]> table=new Hashtable<String, String[]>();
		for(int i=firstLine; i<lines.length; i++){
			String[] line=lines[i];
			table.put(line[col], line);
		}
		return table;
	}
	
	private static int findMaxWidth(String[][] matrix){
		int max=0;
		for(String[] line : matrix){
			if(line!=null && max<line.length){max=line.length;}
		}
		return max;
	}
	
}
