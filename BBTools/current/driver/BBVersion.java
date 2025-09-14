package driver;

import shared.Shared;

public class BBVersion {
	
	public static void main(String[] args){
		System.out.println(Shared.BBTOOLS_VERSION_STRING);
		if(args.length>0){System.out.println(Shared.BBMAP_VERSION_NAME);}
	}
	
}
