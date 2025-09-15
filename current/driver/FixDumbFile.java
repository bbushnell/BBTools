package driver;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Set;

import fileIO.ReadWrite;
import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date May 15, 2014
 *
 */
public class FixDumbFile {
	
	public static void main(String[] args){
		
		String in=args[0];
		String out=args[1];
		
		TextFile tf=new TextFile(in);
		
		LinkedHashMap<String, ArrayList<String[]>> map=new LinkedHashMap<String, ArrayList<String[]>>();
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(!s.startsWith("library_name")){
				String[] line=s.split("\t");
				String key=line[0];
				ArrayList<String[]> list=map.get(key);
				if(list==null){
					list=new ArrayList<String[]>();
					map.put(key, list);
				}else if(s.contains("\tmode\t")){
					
				}
				list.add(line);
			}
		}
		
		tf.close();
		
		StringBuilder sb=new StringBuilder();
		sb.append("library_name\trun_date");
		Set<String> keys=map.keySet();
		{
			String key0=keys.iterator().next();
			ArrayList<String[]> list0=map.get(key0);
			for(String[] term : list0){
				sb.append('\t').append(term[2]);
			}
			sb.append('\n');
		}
		
		for(String key : keys){
			ArrayList<String[]> list=map.get(key);
			String[] term0=list.get(0);
			sb.append(term0[0]);
			sb.append('\t').append(term0[1]);
			for(String[] term : list){
				sb.append('\t').append(term[3]);
			}
			sb.append('\n');
		}
		
		ReadWrite.writeString(sb, out);
		
	}
	
}
