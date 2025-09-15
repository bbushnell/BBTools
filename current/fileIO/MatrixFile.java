package fileIO;
import dna.Matrix;



public class MatrixFile extends TextFile{
	
	public static void main(String[] args){
		
		try {
			//Name of mat file
			String name=args[0];
			
			MatrixFile mat=new MatrixFile(name);
			
			String s=null;
			
			for(s=mat.readLine(); s!=null; s=mat.readLine()){
				System.out.println(s);
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}
	
	
	public MatrixFile(String name){super(name, false);}
	
	@Override
	public String nextLine(){
		String line=readLine();
		
		while(line!=null && line.charAt(0)!='{' && line.charAt(0)!='/'){
			line=readLine();
		}
		return line;
	}
	
	public Matrix nextMatrix(){
		String line;
		String[] split;
		
		line=nextLine();
		if(line==null || line.startsWith("//end")){return null;}
		
		assert(line.startsWith("//name: ")) : line;
		String name=line.replace("//name: ","").trim();
		
		line=nextLine();
		assert(line.startsWith("//size: ")) : line;
		line=line.replace("//size: ","");
		split=line.split("x");
		int length=Integer.parseInt(split[0]);
		int width=Integer.parseInt(split[1]);
		
		line=nextLine();
		assert(line.startsWith("//prefix: ")) : line;
		line=line.replace("//prefix: ","");
		int prefix=Integer.parseInt(line);
		
		line=nextLine();
		assert(line.startsWith("//count: ")) : line;
		line=line.replace("//count: ","");
		int count=Integer.parseInt(line);
		
		
		float[][] grid=new float[length][width];
		for(int i=0; i<length; i++){
			line=nextLine();
			
			while(line.startsWith("//")){line=nextLine();}
			
			assert(line.startsWith("{"));
			if(line.endsWith(",")){line=line.substring(0, line.length()-1);}
			assert(line.endsWith("}"));
			line=line.replace("{", "").replace("}", "").replace(" ", "");
			split=line.split(",");
			assert(split.length==width);
			for(int j=0; j<split.length; j++){
				grid[i][j]=Float.parseFloat(split[j]);
			}
		}
		
		return new Matrix(grid, prefix, name);
	}
	
}
