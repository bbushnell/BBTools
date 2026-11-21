package processor;

import map.ObjectIntMap;
import shared.Tools;

/**
 * Allocates collision-free IDs in one serial mate namespace.
 * The caller supplies IDs after description trimming, normalization and decoration.
 * Identifier scope compares the prefix before the first whitespace; header scope
 * compares the complete ID. Unchanged IDs reuse the caller's String. A collision
 * appends _2, _3, etc. before a terminal mate annotation. Identifier scope keeps
 * descriptions outside the renamed key; header scope renames within the full header.
 *
 * One primitive counter is retained per distinct input key, never per generated ID.
 * Counter k means the original key and every suffix from _2 through _k are occupied,
 * either emitted by this history or skipped because an earlier input occupied them.
 * Thus a previously generated input key can be recognized from its canonical parent.
 * Serial worker clones may share this instance; concurrent calls are unsupported.
 * @author Sucrose
 */
public final class UniqueReadNameAllocator{

	/** @param fullHeaderScope_ True for the complete header, false for its identifier prefix. */
	public UniqueReadNameAllocator(final boolean fullHeaderScope_){fullHeaderScope=fullHeaderScope_;}

	/** Parses the independent uniquenamescope selector without enabling uniqueness or trimming.
	 * @param value Identifier or header, ignoring case
	 * @return False for identifier scope, true for complete-header scope
	 * @throws IllegalArgumentException For null or any unsupported value
	 */
	public static boolean parseScope(final String value){
		if("identifier".equalsIgnoreCase(value)){return false;}
		if("header".equalsIgnoreCase(value)){return true;}
		throw new IllegalArgumentException("uniquenamescope must be identifier or header: "+value);
	}

	/**
	 * Returns an unused ID, preserving any description outside the selected key.
	 * A terminal /1, /2, spaced /1 or /2, or spaced 1: or 2: stays after the suffix.
	 * Natural suffixes with leading zeroes, values below two or values beyond the
	 * generated int range are ordinary valid names. Counter exhaustion fails loudly.
	 * @param id Nonnull ID after caller-owned transformations
	 * @return Original String when unused, otherwise a renamed String
	 * @throws IllegalStateException When no further int ordinal is available for this input key
	 */
	public String allocate(final String id){
		assert(id!=null) : "The caller must supply its numeric fallback before allocating a null read ID.";
		final int whitespace=fullHeaderScope ? -1 : Tools.indexOfWhitespace(id);
		final int keyEnd=whitespace<0 ? id.length() : whitespace;
		probe.set(id, keyEnd, keyEnd, keyEnd);
		final int marker=markerStart(probe);
		final int previous=histories.get(probe);
		assert(previous==-1 || previous>=1) : "Name histories store occupied ordinals >=1; -1 alone means absent: "+previous;
		if(previous<0 && !generatedBefore(id, marker, keyEnd)){
			histories.put(new Key(id.substring(0, keyEnd)), 1);
			return id;
		}

		int ordinal=Math.max(previous, 1);
		String renamed;
		while(true){
			if(ordinal==Integer.MAX_VALUE){
				throw new IllegalStateException("Unique-name suffix counter exhausted for: "+id);
			}
			ordinal++;
			builder.setLength(0);
			builder.append(id, 0, marker).append('_').append(ordinal).append(id, marker, id.length());
			renamed=builder.toString();
			final int renamedEnd=keyEnd+(renamed.length()-id.length());
			probe.set(renamed, renamedEnd, renamedEnd, renamedEnd);
			//This candidate has exactly one canonical parent: the current input key.
			//Its ordinal exceeds that parent's counter, so only an earlier input can reserve it.
			if(histories.get(probe)<0){break;}
		}
		probe.set(id, keyEnd, keyEnd, keyEnd);
		//ObjectIntMap retains its existing key on update; never insert the mutable probe.
		histories.put(previous<0 ? new Key(id.substring(0, keyEnd)) : probe, ordinal);
		return renamed;
	}

	/** Tests the sole possible generated parent, without retaining a generated-name registry. */
	private boolean generatedBefore(final String id, final int marker, final int keyEnd){
		final int underscore=id.lastIndexOf('_', marker-1), start=underscore+1;
		if(underscore<0 || start>=marker || id.charAt(start)<'1' || id.charAt(start)>'9'){return false;}
		int ordinal=0;
		for(int i=start; i<marker; i++){
			final int digit=id.charAt(i)-'0';
			//Check before multiplication: arbitrarily long natural decimal suffixes remain valid.
			if(digit<0 || digit>9 || ordinal>(Integer.MAX_VALUE-digit)/10){return false;}
			ordinal=ordinal*10+digit;
		}
		if(ordinal<2){return false;}
		probe.set(id, underscore, marker, keyEnd);
		//Removing _k is insufficient: that parent must generate at the same boundary.
		//For example, x/1 generates x_2/1, never x/1_2; whitespace before mate tags matters too.
		if(markerStart(probe)!=underscore){return false;}
		return histories.get(probe)>=ordinal;
	}

	/** Finds a terminal mate annotation inside the selected key, including its preceding whitespace. */
	private static int markerStart(final Key key){
		assert(key.text!=null) : "Mate-marker scans require a configured identifier/header lookup.";
		final int end=key.length;
		int marker=end;
		if(end>=2 && key.charAt(end-2)=='/' && (key.charAt(end-1)=='1' || key.charAt(end-1)=='2')){
			marker=end-2;
		}else if(end>=3 && key.charAt(end-1)==':' && (key.charAt(end-2)=='1' || key.charAt(end-2)=='2') &&
			Character.isWhitespace(key.charAt(end-3))){
			marker=end-3;
		}
		if(marker<end){
			while(marker>0 && Character.isWhitespace(key.charAt(marker-1))){marker--;}
		}
		return marker;
	}

	/**
	 * Stored keys contain a stable compact String. The separate mutable lookup key
	 * projects a prefix and optional mate-annotation tail, avoiding a substring per
	 * read or parent lookup. It is never inserted for an absent input history.
	 */
	private static final class Key{
		private Key(){}
		private Key(final String text_){set(text_, text_.length(), text_.length(), text_.length());}

		/** Configures prefix+tail lookup; removing an ordinal leaves the mate annotation intact. */
		private void set(final String text_, final int prefixEnd_, final int tailStart_, final int tailEnd_){
			assert(text_!=null && prefixEnd_>=0 && prefixEnd_<=tailStart_ && tailStart_<=tailEnd_ && tailEnd_<=text_.length()) :
				"Name lookup must project ordered ranges of the original ID: "+prefixEnd_+", "+tailStart_+", "+tailEnd_;
			text=text_; prefixEnd=prefixEnd_; tailStart=tailStart_; tailEnd=tailEnd_;
			length=prefixEnd+(tailEnd-tailStart);
			if(prefixEnd==text.length() && tailStart==tailEnd){hash=text.hashCode();}
			else{
				hash=0;
				for(int i=0; i<prefixEnd; i++){hash=31*hash+text.charAt(i);}
				for(int i=tailStart; i<tailEnd; i++){hash=31*hash+text.charAt(i);}
			}
		}

		@Override
		public int hashCode(){return hash;}

		/** Reads the projected key, excluding a removed suffix between prefix and mate annotation. */
		private char charAt(final int index){
			assert(index>=0 && index<length) : "Mate-marker scans must stay within the projected name: "+index+" / "+length;
			return text.charAt(index<prefixEnd ? index : tailStart+index-prefixEnd);
		}

		@Override
		public boolean equals(final Object other){
			if(this==other){return true;}
			if(!(other instanceof Key)){return false;}
			final Key key=(Key)other;
			if(length!=key.length){return false;}
			if(prefixEnd==length){
				return text.regionMatches(0, key.text, 0, key.prefixEnd) &&
					text.regionMatches(key.prefixEnd, key.text, key.tailStart, key.tailEnd-key.tailStart);
			}
			for(int i=0; i<length; i++){
				final char a=text.charAt(i<prefixEnd ? i : tailStart+i-prefixEnd);
				final char b=key.text.charAt(i<key.prefixEnd ? i : key.tailStart+i-key.prefixEnd);
				if(a!=b){return false;}
			}
			return true;
		}

		private String text;
		private int prefixEnd, tailStart, tailEnd, length, hash;
	}

	/** Scope is immutable for the entire serial history. */
	private final boolean fullHeaderScope;
	/** Only distinct input keys are retained; values are the last occupied generated ordinal. */
	private final ObjectIntMap<Key> histories=new ObjectIntMap<Key>(Key.class);
	/** Mutable lookup scratch is borrowed for updates only, never retained as a map key. */
	private final Key probe=new Key();
	/** Reused output scratch; only the resulting ID String survives each allocation. */
	private final StringBuilder builder=new StringBuilder(64);
}
