package align2;

import java.io.File;

import fileIO.ReadWrite;
import shared.Parse;
import shared.Tools;

public class PrintTime {
	
	public static void main(String[] args){
		long millis=System.currentTimeMillis();
		
		if(args==null || args.length<1){
			System.err.println("Time:\t"+millis);
		}
		
		if(args!=null && args.length>0){
			File f=new File(args[0]);
			if(f.exists()){
				String s=ReadWrite.readString(args[0]);
//				TextFile tf=new TextFile(args[0], false, false);
//				String s=tf.nextLine();
//				tf.close();
				long old=Long.parseLong(s);
				long elapsed=millis-old;
				if(args.length<2 || Parse.parseBoolean(args[1])){
					System.out.println("Elapsed:\t"+Tools.format("%.2f", elapsed/1000d));
					if(true){
						System.err.println("Elapsed:\t"+Tools.format("%.2f", elapsed/1000d));
					}
				}
			}
			f=null;
			ReadWrite.writeString(millis+"", args[0]);
		}
	}
	
}
