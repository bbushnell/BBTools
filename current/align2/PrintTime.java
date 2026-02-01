package align2;

import java.io.File;

import fileIO.ReadWrite;
import parse.Parse;
import shared.Tools;

/**
 * Simple timing utility for measuring elapsed time between program executions.
 * Stores timestamps in files and calculates time differences for performance monitoring.
 * Used within the BBTools alignment framework for benchmarking and execution tracking.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class PrintTime {
	
	/**
	 * Program entry point for timing operations.
	 * Prints current timestamp if no arguments provided.
	 * If file argument given, calculates elapsed time since last stored timestamp.
	 * Stores current timestamp in specified file for future elapsed time calculations.
	 *
	 * @param args Command-line arguments: [0] = timestamp file path, [1] = optional boolean to control output
	 */
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
