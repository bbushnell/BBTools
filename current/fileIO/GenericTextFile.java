package fileIO;

import java.util.ArrayList;

public class GenericTextFile extends TextFile {

	public GenericTextFile(String name) {
		super(name, false);
	}
	
	

	
	public String[] toLines(){
		
		String s=null;
		ArrayList<String> list=new ArrayList<String>(4096);
		
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		
		return list.toArray(new String[list.size()]);
		
	}
	
	@Override
	public String nextLine(){
		String line=readLine();
		while(line!=null && false){
			line=readLine();
		}
		return line;
	}
	

}
