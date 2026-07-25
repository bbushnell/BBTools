package shared;

import java.util.HashMap;
import java.util.HashSet;

/**
 * Maps exact names and, when unambiguous, names truncated at first whitespace.
 * Exact matches always take precedence.
 *
 * @param <T> mapped value type
 */
public class NameMapper<T> {

	public NameMapper(){
		this(16);
	}

	public NameMapper(int initialSize){
		table=new HashMap<String, T>(initialSize);
		exactNames=new HashSet<String>(initialSize);
		ambiguousAliases=new HashSet<String>();
		fallbackQueries=new HashMap<String, String>();
	}

	/**
	 * Adds or replaces an exact name and registers its first-token alias.
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

	public T getExact(String name){
		return exactNames.contains(name) ? table.get(name) : null;
	}

	/**
	 * Returns an exact match, then an unambiguous first-token match.
	 * If multiple stored names or multiple distinct lookup names share the same
	 * fallback token, throws rather than selecting one.
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

	public int size(){
		return exactNames.size();
	}

	private static IllegalArgumentException ambiguous(String alias){
		return new IllegalArgumentException("Ambiguous sequence name '"+alias+
				"' matches multiple records. Use unique first tokens or consistent full names.");
	}

	private final HashMap<String, T> table;
	private final HashSet<String> exactNames;
	private final HashSet<String> ambiguousAliases;
	private final HashMap<String, String> fallbackQueries;
}
