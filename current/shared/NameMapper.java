package shared;

import java.util.HashMap;
import java.util.HashSet;

/**
 * Maps exact names and, when unambiguous, names truncated at first whitespace.
 * Exact matches always take precedence.
 * Values are retained by reference. Lookups may record the first expanded query
 * for an alias, so even read-only use requires thread confinement or external
 * synchronization. Names are case-sensitive; whitespace follows
 * {@link Tools#trimToWhitespace(String)}.
 *
 * @param <T> mapped value type
 */
public class NameMapper<T>{

	/** Creates an empty mapper with an initial exact-name capacity of 16. */
	public NameMapper(){this(16);}

	/** Creates an empty mapper; capacity does not limit the number of names.
	 * @param initialSize Nonnegative initial capacity for exact-name collections */
	public NameMapper(int initialSize){
		table=new HashMap<String, T>(initialSize);
		exactNames=new HashSet<String>(initialSize);
		ambiguousAliases=new HashSet<String>();
		fallbackQueries=new HashMap<String, String>();
	}

	/**
	 * Adds or replaces an exact name and registers its first-token alias.
	 * Replacing an exact name updates its unique alias but does not clear ambiguity
	 * or previously observed fallback queries.
	 * @param name Non-null exact name, including any description
	 * @param value Non-null value retained without copying
	 * @return the previous exact value, or null
	 */
	public T put(String name, T value){
		if(name==null || value==null){throw new NullPointerException();}
		final boolean replacing=exactNames.contains(name);
		final T old=(replacing ? table.get(name) : null);
		final String alias=Tools.trimToWhitespace(name);

		if(replacing){
			table.put(name, value);
			if(!alias.equals(name) && !ambiguousAliases.contains(alias) && table.get(alias)==old){
				table.put(alias, value);
			}
			return old;
		}

		if(table.containsKey(name)){ambiguousAliases.add(name);}
		exactNames.add(name);
		table.put(name, value);
		if(alias.equals(name) || ambiguousAliases.contains(alias)){return null;}

		T prior=table.get(alias);
		if(prior==null){
			table.put(alias, value);
		}else{
			ambiguousAliases.add(alias);
			if(!exactNames.contains(alias)){table.remove(alias);}
		}
		return null;
	}

	/** Returns only a stored exact-name value; null or absent names return null. */
	public T getExact(String name){return exactNames.contains(name) ? table.get(name) : null;}

	/**
	 * Returns an exact match, then an unambiguous first-token match.
	 * If multiple stored names or multiple distinct lookup names share the same
	 * fallback token, throws rather than selecting one.
	 * A direct exact-name or stored-alias hit returns before query-history checks.
	 * @param name Query name; null returns null
	 * @return Retained value, or null when no name or alias matches
	 * @throws IllegalArgumentException If the required fallback is ambiguous
	 */
	public T get(String name){
		if(name==null){return null;}
		T value=table.get(name);
		if(value!=null){return value;}
		if(name.isEmpty()){return null;}

		String alias=Tools.trimToWhitespace(name);
		if(ambiguousAliases.contains(alias)){throw ambiguous(alias);}
		if(alias.equals(name)){return null;}
		value=table.get(alias);
		if(value==null){return null;}

		String prior=fallbackQueries.get(alias);
		if(prior==null){
			fallbackQueries.put(alias, name);
		}else if(!prior.equals(name)){
			throw ambiguous(alias);
		}
		return value;
	}

	/** Returns the number of distinct exact names, excluding generated aliases. */
	public int size(){return exactNames.size();}

	/** Builds the common failure message for stored-name and query-name collisions. */
	private static IllegalArgumentException ambiguous(String alias){
		return new IllegalArgumentException("Ambiguous sequence name '"+alias+
				"' matches multiple records. Use unique first tokens or consistent full names.");
	}

	/** Exact names and currently usable first-token aliases share this lookup table. */
	private final HashMap<String, T> table;
	/** Distinguishes exact names from generated table entries. */
	private final HashSet<String> exactNames;
	/** Tokens with multiple stored interpretations; exact table hits still win. */
	private final HashSet<String> ambiguousAliases;
	/** First expanded query observed for each successfully used fallback token. */
	private final HashMap<String, String> fallbackQueries;
}
