package fileIO;


public class ArrayFile extends TextFile{
	
	public static void main(String[] args){
		
		try {
			//Name of mat file
			String name=args[0];
			
			ArrayFile mat=new ArrayFile(name);
			
			String s=null;
			
			for(s=mat.readLine(); s!=null; s=mat.readLine()){
				System.out.println(s);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}
	
	
	public ArrayFile(String name){super(name, false);}
	
	@Override
	public String nextLine(){
		String line=readLine();
		char c=line.charAt(0);
		
		while(line!=null && c!='{' && c!='/'){
			line=readLine();
			c=line.charAt(0);
		}
		return line;
	}
	
	public float[] nextArray(){
		String line;
		String[] split;
		
		line=nextLine();
		if(line==null || line.startsWith("//end")){return null;}
		
		assert(line.startsWith("//name: ")) : line;
		String name=line.replace("//name: ","").trim();
		
		line=nextLine();
		assert(line.startsWith("//size: ")) : line;
		line=line.replace("//size: ","");
		int length=Integer.parseInt(line);
		
		
		float[] grid=new float[length];
		
		line=nextLine();
		assert(line.startsWith("{"));
		if(line.endsWith(",")){line=line.substring(0, line.length()-1);}
		assert(line.endsWith("}"));
		line=line.replace("{", "").replace("}", "").replace(" ", "");
		split=line.split(",");
		assert(split.length==length);
		for(int i=0; i<split.length; i++){
			grid[i]=Float.parseFloat(split[i]);
		}
		
		return grid;
	}
	
}
