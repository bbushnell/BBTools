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
	
	/**
	 * Determines if a read or its mate maps to a whitelisted scaffold.
	 * Returns true if either the read or its mate is mapped to a scaffold
	 * present in the whitelist.
	 *
	 * @param r The read to check (may be null)
	 * @return true if read or mate maps to whitelisted scaffold, false otherwise
	 */
	public static boolean inWhitelist(Read r){
		return r==null ? false : (inWhitelist2(r) || inWhitelist2(r.mate));
	}
	
	/**
	 * Helper method to check if a single read maps to a whitelisted scaffold.
	 * Verifies the read is mapped and its scaffold name exists in the whitelist.
	 * @param r The read to check (may be null)
	 * @return true if read maps to whitelisted scaffold, false otherwise
	 */
	private static boolean inWhitelist2(Read r){
		if(r==null || !r.mapped() || whitelist==null || whitelist.isEmpty()){return false;}
		byte[] name=r.getScaffoldName(false);
		return (name!=null && whitelist.contains(new String(name)));
	}
	
	/**
	 * Determines if a read pair should be filtered based on blacklist criteria.
	 * Complex logic filters read pairs when one or both reads map to blacklisted
	 * scaffolds, accounting for mate mapping status.
	 *
	 * @param r The read to check (may be null)
	 * @return true if read pair should be filtered out, false otherwise
	 */
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
	
	/**
	 * Helper method to check if a single read maps to a blacklisted scaffold.
	 * Verifies the read is mapped and its scaffold name exists in the blacklist.
	 * @param r The read to check (may be null)
	 * @return true if read maps to blacklisted scaffold, false otherwise
	 */
	private static boolean inBlacklist2(Read r){
		if(r==null || !r.mapped() || blacklist==null || blacklist.isEmpty()){return false;}
		byte[] name=r.getScaffoldName(false);
		return (name!=null && blacklist.contains(new String(name)));
	}
	
	/**
	 * Loads scaffold names from a file into the blacklist.
	 * Convenience method that calls addToSet with black=true.
	 * @param fname Path to file containing scaffold names to blacklist
	 */
	public static void addToBlacklist(String fname){
		addToSet(fname, true);
	}
	
	/**
	 * Loads scaffold names from a file into the whitelist.
	 * Convenience method that calls addToSet with black=false.
	 * @param fname Path to file containing scaffold names to whitelist
	 */
	public static void addToWhitelist(String fname){
		addToSet(fname, false);
	}
	
	/**
	 * Reads scaffold names from a file and adds them to blacklist or whitelist.
	 * Supports both FASTA format (>scaffold_name) and plain text (one name per line).
	 * Thread-safe operation that auto-detects file format and reports duplicates.
	 *
	 * @param fname Path to input file containing scaffold names
	 * @param black true to add to blacklist, false for whitelist
	 * @return Number of unique scaffold names added (excludes duplicates)
	 */
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

	/** Returns true if blacklist exists and contains entries */
	public static boolean hasBlacklist(){return blacklist!=null && !blacklist.isEmpty();}
	/** Returns true if whitelist exists and contains entries */
	public static boolean hasWhitelist(){return whitelist!=null && !whitelist.isEmpty();}

	/** Clears the blacklist by setting it to null */
	public static void clearBlacklist(){blacklist=null;}
	/** Clears the whitelist by setting it to null */
	public static void clearWhitelist(){whitelist=null;}
	
	/** Set of scaffold names to exclude from processing */
	private static HashSet<String> blacklist=null;
	/** Set of scaffold names to exclusively include in processing */
	private static HashSet<String> whitelist=null;
	
}
