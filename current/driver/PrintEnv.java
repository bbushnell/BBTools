package driver;

import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Map;

import shared.Shared;

/**
 * Utility class for displaying system environment information.
 * Prints all environment variables in sorted order along with current timestamp
 * and local hostname.
 *
 * @author Brian Bushnell
 * @date Apr 4, 2013
 */
public class PrintEnv {
	
	public static void main(String[] args){
		
		Date d=new Date();
		System.out.println("Time: "+d.getTime()+" = "+d+"\n");
		
		Map<String, String> env=System.getenv();
		ArrayList<String> keys=new ArrayList<String>(env.keySet());
		Shared.sort(keys);
		for(String s : keys){
			System.out.println(s+"\t"+env.get(s));
		}
		try {
			java.net.InetAddress localMachine = java.net.InetAddress.getLocalHost();
			System.out.println("Hostname of local machine: " + localMachine.getHostName());
		} catch (UnknownHostException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
