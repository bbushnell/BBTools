package stream;

import java.util.HashMap;

import fileIO.TextFile;
import parse.LineParserS1;

/**
 * Renames reference scaffolds in a SAM/BAM stream from a 2-column old&rarr;new TSV.
 * Rewrites the RNAME and RNEXT fields of each record (via {@link #renameRecord(SamLine)})
 * and the SN: field of {@code @SQ} header lines (via {@link #renameHeaderLine(String)}).
 * Names absent from the map pass through unchanged, so partial renaming is allowed.
 *
 * <p>Coordinates are NOT altered &mdash; this is purely a naming transform, valid only
 * when the old and new names refer to the SAME assembly (e.g. UCSC "chr1" vs Ensembl "1").
 * It does NOT perform liftover between different reference builds.
 *
 * <p>Renaming reads the map without changing it, but {@link #map()} exposes the mutable
 * backing map. Callers sharing this object must publish it safely, avoid concurrent
 * map changes and coordinate mutation of each SAM record. Keep the RNAME storage
 * mode fixed and consistent with the records being renamed.
 *
 * @author UMP45
 */
public class ScaffoldRenamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads a nonempty old&rarr;new name map from a TSV. Empty lines and lines beginning
	 * with '#' are skipped; the first two tab-delimited fields are taken as old, new.
	 * Later columns are ignored, names are not trimmed and later duplicate keys replace
	 * earlier values. Supply valid scaffold names for the same reference assembly;
	 * this loader does not validate the assembly relationship or uniqueness of new names.
	 * @param tsvPath Path to the 2-column old&lt;tab&gt;new TSV
	 */
	public ScaffoldRenamer(String tsvPath){
		map=new HashMap<String, String>();
		final LineParserS1 lp=new LineParserS1('\t');
		TextFile tf=new TextFile(tsvPath);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			lp.set(line);
			int terms=lp.terms();
			//String.split discards trailing empty fields; preserve its column-count check.
			while(terms>0 && lp.length(terms-1)==0){terms--;}
			assert(terms>=2) : "Expected 2-column 'old<tab>new' TSV, got: "+line;
			map.put(lp.parseString(0), lp.parseString(1));
		}
		tf.close();
		assert(!map.isEmpty()) : "No rename pairs loaded from "+tsvPath;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Rewrites the first SN: field of an {@code @SQ} header line via the rename map.
	 * Uses an @SQ prefix check and stops at the first SN: field even if it has no mapping.
	 * A mapping hit rebuilds the line without trailing empty fields, preserving the legacy
	 * split/join result; otherwise the original string is returned.
	 * @param line A single SAM header line; null is returned unchanged
	 * @return Rebuilt line on a nonnull mapping hit, or the original reference otherwise
	 */
	public String renameHeaderLine(String line){
		if(line==null || !line.startsWith("@SQ")){return line;}
		int tab=line.indexOf('\t');
		while(tab>=0){
			final int start=tab+1, next=line.indexOf('\t', start);
			final int stop=(next<0 ? line.length() : next);
			if(stop-start>=3 && line.startsWith("SN:", start)){
				final String renamed=map.get(line.substring(start+3, stop));
				if(renamed==null){return line;}
				int end=line.length();
				while(end>0 && line.charAt(end-1)=='\t'){end--;}
				assert(end>=stop) : "SN: makes the field nonempty; trailing-tab trimming must retain it before appending the suffix: stop="+stop+", end="+end;
				return new StringBuilder(line.length()).append(line, 0, start+3).append(renamed)
						.append(line, stop, end).toString();
			}
			tab=next;
		}
		return line;
	}

	/**
	 * Rewrites the RNAME and RNEXT fields of a record via the rename map, in place.
	 * Null RNAME/RNEXT fields and RNEXT="=" (mate on the same scaffold) are left alone.
	 * The unmapped flag itself does not prevent renaming a supplied reference name.
	 * Reads RNAME via {@link SamLine#rnameS()} (mode-agnostic) and writes via the setter
	 * matching {@link SamLine#RNAME_AS_BYTES}; keep that mode consistent with the record.
	 * Byte conversions in this method use the platform default charset; SamLine setters
	 * apply their normal reference-name canonicalization. Coordinates and CIGAR are untouched.
	 * @param sl Nonnull record owned by the caller for the duration of this mutation
	 * @return true if a nonnull mapping was assigned to RNAME or RNEXT, even when its text is unchanged
	 */
	public boolean renameRecord(SamLine sl){
		boolean changed=false;
		final String rn=sl.rnameS();
		if(rn!=null){
			final String renamed=map.get(rn);
			if(renamed!=null){
				if(SamLine.RNAME_AS_BYTES){sl.setRname(renamed.getBytes());}
				else{sl.setRnameS(renamed);}
				changed=true;
			}
		}
		final byte[] rx=sl.rnext();
		if(rx!=null && !(rx.length==1 && rx[0]=='=')){
			final String renamed=map.get(new String(rx));
			if(renamed!=null){sl.setRnext(renamed.getBytes()); changed=true;}
		}
		return changed;
	}

	/**
	 * Exposes the live mutable backing map, without copying or a read-only wrapper.
	 * Changes affect subsequent renaming; callers must coordinate access themselves.
	 * @return The same old&rarr;new map retained by this renamer
	 */
	public HashMap<String, String> map(){return map;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Old&rarr;new scaffold names; the reference is final but map() exposes mutable contents. */
	private final HashMap<String, String> map;

}
