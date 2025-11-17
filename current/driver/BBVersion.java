package driver;

import shared.Shared;

/**
 * Utility class for displaying BBTools version information.
 * Simple command-line tool that prints version strings from Shared constants.
 * @author Brian Bushnell
 */
public class BBVersion {
	
	/**
	 * Program entry point that prints BBTools version information.
	 * Always prints BBTOOLS_VERSION_STRING, optionally prints BBMAP_VERSION_NAME
	 * if any command-line arguments are provided.
	 * @param args Command-line arguments; presence triggers additional version output
	 */
	public static void main(String[] args){
		System.out.println(Shared.BBTOOLS_VERSION_STRING);
		if(args.length>0){System.out.println(Shared.BBMAP_VERSION_NAME);}
	}
	
}
