package fileIO;

import java.io.IOException;
import java.io.InputStream;

public class OpenFile {
	
	public static void main(String[] args){
		InputStream is=ReadWrite.getRawInputStream(args[0], false);
		byte[] line=new byte[100];
		try {
			int r=is.read(line, 0, 100);
			if(r>0){
				System.err.println("'"+new String(line, 0, r)+"'");
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
