package stream;

import java.util.HashMap;

import fileIO.TextFile;

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
 * <p>Stateless per call once constructed (the map is read-only), so a single instance is
 * safe to share across the wrapper's worker threads.
 *
 * @author UMP45
 */
public class ScaffoldRenamer {

	/**
	 * Loads the old&rarr;new name map from a 2-column TSV. Blank lines and lines beginning
	 * with '#' are skipped; the first two tab-delimited fields are taken as old, new.
	 * @param tsvPath Path to the 2-column old&lt;tab&gt;new TSV
	 */
	public ScaffoldRenamer(String tsvPath){
		map=new HashMap<String, String>();
		TextFile tf=new TextFile(tsvPath);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.length()==0 || line.charAt(0)=='#'){continue;}
			String[] s=line.split("\t");
			assert(s.length>=2) : "Expected 2-column 'old<tab>new' TSV, got: "+line;
			map.put(s[0], s[1]);
		}
		tf.close();
		assert(!map.isEmpty()) : "No rename pairs loaded from "+tsvPath;
	}

	/**
	 * Rewrites the SN: field of an {@code @SQ} header line via the rename map.
	 * Non-@SQ lines, and @SQ lines whose SN is not in the map, are returned unchanged.
	 * @param line A single SAM header line
	 * @return The line with SN: renamed, or the original line if nothing changed
	 */
	public String renameHeaderLine(String line){
		if(line==null || !line.startsWith("@SQ")){return line;}
		String[] parts=line.split("\t");
		for(int i=1; i<parts.length; i++){
			if(parts[i].startsWith("SN:")){
				String old=parts[i].substring(3);
				String renamed=map.get(old);
				if(renamed!=null){
					parts[i]="SN:"+renamed;
					return String.join("\t", parts);
				}
				break;
			}
		}
		return line;
	}

	/**
	 * Rewrites the RNAME and RNEXT fields of a record via the rename map, in place.
	 * RNEXT="=" (mate on same scaffold as RNAME) and unmapped (null) fields are left alone.
	 * Reads RNAME via {@link SamLine#rnameS()} (mode-agnostic) and writes via the setter
	 * matching {@link SamLine#RNAME_AS_BYTES}.
	 * @param sl The record to rename in place
	 * @return true if RNAME or RNEXT was changed
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

	/** @return The old&rarr;new rename map (read-only use). */
	public HashMap<String, String> map(){return map;}

	/** Old&rarr;new scaffold name map. */
	private final HashMap<String, String> map;

}
