package align2;

import java.util.HashSet;

import fileIO.TextFile;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Mar 14, 2013
 *
 */
public class Blacklist {
	
	public static boolean inWhitelist(Read r){
		return r==null ? false : (inWhitelist2(r) || inWhitelist2(r.mate));
	}
	
	private static boolean inWhitelist2(Read r){
		if(r==null || !r.mapped() || whitelist==null || whitelist.isEmpty()){return false;}
		byte[] name=r.getScaffoldName(false);
		return (name!=null && whitelist.contains(new String(name)));
	}
	
	public static boolean inBlacklist(Read r){
		if(r==null){return false;}
		boolean a=inBlacklist2(r);
		boolean b=inBlacklist2(r.mate);
		if(!a && !b){return false;}
		if(a){
			return b || r.mate==null || !r.mate.mapped();
		}
		return b && !r.mapped();
	}
	
	private static boolean inBlacklist2(Read r){
		if(r==null || !r.mapped() || blacklist==null || blacklist.isEmpty()){return false;}
		byte[] name=r.getScaffoldName(false);
		return (name!=null && blacklist.contains(new String(name)));
	}
	
	public static void addToBlacklist(String fname){
		addToSet(fname, true);
	}
	
	public static void addToWhitelist(String fname){
		addToSet(fname, false);
	}
	
	public static synchronized int addToSet(String fname, boolean black){
		final HashSet<String> set;
		int added=0, overwritten=0;
		if(black){
			if(blacklist==null){blacklist=new HashSet<String>(4001);}
			set=blacklist;
		}else{
			if(whitelist==null){whitelist=new HashSet<String>(4001);}
			set=whitelist;
		}
		TextFile tf=new TextFile(fname, false);
		String line=tf.nextLine();
		if(line==null){return 0;}
		final boolean fasta=(line.charAt(0)=='>');
		System.err.println("Detected "+(black ? "black" : "white")+"list file "+fname+" as "+(fasta ? "" : "non-")+"fasta-formatted.");
		while(line!=null){
			String key=null;
			if(fasta){
				if(line.charAt(0)=='>'){key=new String(line.substring(1));}
			}else{
				key=line;
			}
			if(key!=null){
				boolean b=set.add(key);
				added++;
				if(!b){
					if(overwritten==0){
						System.err.println("Duplicate "+(black ? "black" : "white")+"list key "+key);
						System.err.println("Subsequent duplicates from this file will not be mentioned.");
					}
					overwritten++;
				}
			}
			line=tf.nextLine();
		}
		if(overwritten>0){
			System.err.println("Added "+overwritten+" duplicate keys.");
		}
		return added-overwritten;
	}

	public static boolean hasBlacklist(){return blacklist!=null && !blacklist.isEmpty();}
	public static boolean hasWhitelist(){return whitelist!=null && !whitelist.isEmpty();}

	public static void clearBlacklist(){blacklist=null;}
	public static void clearWhitelist(){whitelist=null;}
	
	private static HashSet<String> blacklist=null;
	private static HashSet<String> whitelist=null;
	
}
