package align2;

import java.util.HashSet;

import fileIO.TextFile;
import stream.Read;

/**
 * Manages scaffold name filtering using blacklist and whitelist mechanisms.
 * Provides read filtering based on scaffold names to include/exclude specific
 * reference sequences during alignment processing. Supports both FASTA and
 * plain text formatted filter files.
 * Names match exactly, including spaces and descriptions; no tokenization or
 * normalization is performed. State is process-global. Load/clear it before
 * mapping workers start; synchronized loading does not make concurrent lookup
 * or clearing safe. Loading appends to existing sets, rather than replacing them.
 *
 * @author Brian Bushnell
 * @date Mar 14, 2013
 */
public class Blacklist{

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
	 * At least one mapped end must be blacklisted, and every mapped end must be
	 * blacklisted. A mapped unlisted mate protects the pair; a missing or unmapped
	 * mate does not. This differs deliberately from the whitelist's either-end rule.
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
	 * The first nonblank line selects FASTA versus literal-line input for the whole
	 * file. FASTA sequence lines are skipped and whole header text after '>' is kept.
	 * Blank lines are skipped by TextFile. Loading is serialized against other
	 * loads only; callers must exclude concurrent lookups and clear operations.
	 * A failure may leave already-added names in the global set; no rollback occurs.
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
		final TextFile tf=new TextFile(fname, false);
		try{
			String line=tf.nextLine();
			if(line==null){return 0;}
			// TextFile.readLine(true) skips empty/whitespace-only lines before charAt(0).
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
					final boolean b=set.add(key);
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
		}finally{
			// EOF does not close TextFile; this also covers empty inputs and read failures.
			tf.close();
		}
	}

	/** Returns true if a blacklist exists and contains entries. */
	public static boolean hasBlacklist(){return blacklist!=null && !blacklist.isEmpty();}
	/** Returns true if a whitelist exists and contains entries. */
	public static boolean hasWhitelist(){return whitelist!=null && !whitelist.isEmpty();}

	/** Clears global blacklist state; caller must exclude concurrent loading/lookup. */
	public static void clearBlacklist(){blacklist=null;}
	/** Clears global whitelist state; caller must exclude concurrent loading/lookup. */
	public static void clearWhitelist(){whitelist=null;}

	private static HashSet<String> blacklist=null;
	private static HashSet<String> whitelist=null;

}
