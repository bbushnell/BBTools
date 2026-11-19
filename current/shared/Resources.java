package shared;

import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;

/**
 * Handles loading of BBTools resource files with helpful download
 * instructions when resources are missing.  Wraps Data.findPath()
 * and adds a table mapping known resource filenames to their
 * download locations.
 *
 * @author Noire, Brian Bushnell
 * @date May 19, 2026
 */
public class Resources {

	/** Resolves a known optional resource's download URL without opening files or exiting. */
	public static String downloadURL(String filename){
		if(filename==null || filename.isEmpty()){throw new IllegalArgumentException("Resource filename is required");}
		final String path=(filename.startsWith("?") ? filename.substring(1) : filename).replace('\\', '/');
		if(bbgapArchive(path)!=null){return GITHUB_RELEASES;}
		final String known=RESOURCE_URLS.get(new java.io.File(path).getName());
		if(known!=null){return known;}
		if(path.startsWith("prokcc/") || path.contains("/prokcc/")){return PROKCC_ARCHIVE;}
		return null;
	}

	/** ProkCC's config, models and tables are installed together from one manual download. */
	public static String prokccDownloadInstructions(){
		return prokccDownloadInstructions(PROKCC_ARCHIVE);
	}

	/** Gives extraction instructions for the archive matching a current or legacy resource. */
	public static String prokccDownloadInstructions(String archive){
		if(!isProkccArchive(archive)){throw new IllegalArgumentException("Unknown ProkCC archive: "+archive);}
		return "Download "+archive.substring(archive.lastIndexOf('/')+1)+" from:\n  "+archive+
			"\nUntar it in the BBTools resources/ directory to create resources/prokcc/, then try again.\n";
	}

	/** Gives extraction instructions for a BBGAP compact stop-tail archive. */
	public static String bbgapDownloadInstructions(String archive){
		if(!archive.equals(BBGAP_SHORT_ARCHIVE) && !archive.equals(BBGAP_TAILS_ARCHIVE)){
			throw new IllegalArgumentException("Unknown BBGAP archive: "+archive);
		}
		final String view=archive.equals(BBGAP_SHORT_ARCHIVE) ? "short" : "tails";
		return "Download "+archive+" from:\n  "+GITHUB_RELEASES+
			"\nFrom the BBTools resources/ directory, run tar -xf "+archive+
			" to create bbgap/"+view+"/, then try again.\n";
	}

	/** Identifies the current archive and the explicitly retained previous release. */
	private static boolean isProkccArchive(String url){
		return PROKCC_ARCHIVE.equals(url) || PROKCC_LEGACY_ARCHIVE.equals(url);
	}

	/** Locates a resource file.  If fname contains commas, tries each
	 * candidate in order and returns the first one found.
	 * @param fname Resource name or comma-separated resource candidates
	 * @return Resource path */
	public static String find(String fname){
		return find(fname, true);
	}

	/** Locates a resource file, optionally exiting when no candidate is present.
	 * @param fname Resource name or comma-separated resource candidates
	 * @param exit True to exit when missing; false still prints single-name instructions,
	 * but returns null silently when no comma-separated candidate is found
	 * @return Resource path, or null if missing and exit is false */
	public static String find(String fname, boolean exit){
		if(fname.indexOf(',')>=0){
			String[] parts=fname.split(",");
			ArrayList<String> list=new ArrayList<>(parts.length);
			for(String s : parts){s=s.trim(); if(!s.isEmpty()){list.add(s);}}
			return find(list, exit);
		}
		return findSingle(fname, exit);
	}

	/** Tries each candidate in order, returns the first one found.
	 * Only complains if all are missing.
	 * @param fnames Resource candidates
	 * @return Resource path */
	public static String find(ArrayList<String> fnames){
		return find(fnames, true);
	}

	/** Tries each candidate in order, optionally exiting when all are missing.
	 * @param fnames Resource candidates
	 * @param exit True to print instructions and exit when missing
	 * @return Resource path, or null if missing and exit is false */
	public static String find(ArrayList<String> fnames, boolean exit){
		for(String f : fnames){
			String path=Data.findPath(f, false);
			if(path!=null){return path;}
		}
		if(fnames.isEmpty()){return null;}
		if(!exit){return null;}
		StringBuilder sb=new StringBuilder();
		final ArrayList<String> archives=new ArrayList<String>();
		sb.append("\nERROR: No resource found from candidates:\n");
		for(String f : fnames){
			String bare=f.startsWith("?") ? f.substring(1) : f;
			String url=downloadURL(bare);
			if(isProkccArchive(url) && !archives.contains(url)){archives.add(url);}
			sb.append("  ").append(bare);
			if(url!=null){sb.append("  (").append(url).append(')');}
			sb.append('\n');
		}
		for(String archive:archives){sb.append(prokccDownloadInstructions(archive));}
		sb.append("Place other network files in BBTools/networks/ and other resources in BBTools/resources/.\n");
		System.err.print(sb);
		System.exit(1);
		return null;
	}

	public static void throwMissing(ArrayList<String> fnames){
		if(fnames.isEmpty()){return;}
		final StringBuilder sb=new StringBuilder();
		for(String f : fnames){appendMissingMessage(sb, f);}
		throw new MissingResourceException(sb.toString());
	}

	/** Locates a single resource file.  If missing, prints download
	 * instructions and optionally exits. */
	private static String findSingle(String fname, boolean exit){
		String path=Data.findPath(fname, false);
		if(path!=null){return path;}

		StringBuilder sb=new StringBuilder();
		appendMissingMessage(sb, fname);
		System.err.print(sb);
		if(exit){
			System.exit(1);
		}
		return null;
	}

	private static void appendMissingMessage(StringBuilder sb, String fname){
		String bare=fname;
		if(bare.startsWith("?")){bare=bare.substring(1);}

		String url=downloadURL(bare);
		final String bbgapArchive=bbgapArchive(bare);
		sb.append("\nERROR: Required resource not found: ").append(bare).append('\n');
		if(isProkccArchive(url)){
			sb.append(prokccDownloadInstructions(url));
		}else if(bbgapArchive!=null){
			sb.append(bbgapDownloadInstructions(bbgapArchive));
		}else{
			if(url!=null){
				sb.append("Download it from:\n  ").append(url).append('\n');
			}else{
				sb.append("You may need to download it from:\n  ").append(SOURCEFORGE_URL).append('\n');
			}
			final String destination=Data.isNetworkFile(bare) ? "BBTools/networks/" : resourceDestination(bare);
			sb.append("Place it in ").append(destination).append(" and try again.\n");
		}
	}

	private static String resourceDestination(String bare){
		final String path=bare.replace('\\', '/');
		final int slash=path.lastIndexOf('/');
		return "BBTools/resources/"+(slash<0 ? "" : path.substring(0, slash+1));
	}

	private static String bbgapArchive(String path){
		if(path==null){return null;}
		path=path.replace('\\', '/');
		if(path.startsWith("?")){path=path.substring(1);}
		if(path.equals("bbgap/short/families.tsv") || path.equals("bbgap/short/families.tsv.gz") ||
				path.endsWith("/bbgap/short/families.tsv") || path.endsWith("/bbgap/short/families.tsv.gz")){
			return BBGAP_SHORT_ARCHIVE;
		}
		if(path.equals("bbgap/tails/families.tsv") || path.equals("bbgap/tails/families.tsv.gz") ||
				path.endsWith("/bbgap/tails/families.tsv") || path.endsWith("/bbgap/tails/families.tsv.gz")){
			return BBGAP_TAILS_ARCHIVE;
		}
		return null;
	}

	public static final class MissingResourceException extends RuntimeException {
		public MissingResourceException(String message){super(message);}
	}

	private static final String GITHUB_RELEASES="https://github.com/bbushnell/BBTools/releases/";
	private static final String GITHUB_V3982="https://github.com/bbushnell/BBTools/releases/tag/v39.82/";
	private static final String GITHUB_V3984="https://github.com/bbushnell/BBTools/releases/tag/v39.84/";
	/** Direct download prefix for the v40.00 release (size-filtered DDL + congruent spectra). */
	private static final String GITHUB_V4000="https://github.com/bbushnell/BBTools/releases/download/v40.00/";
	private static final String GITHUB_V4001="https://github.com/bbushnell/BBTools/releases/download/v40.01/";
	private static final String NERSC_URL="https://portal.nersc.gov/cfs/bbtools/";
	private static final String SOURCEFORGE_URL="https://sourceforge.net/projects/bbmap/files/Resources/";
	/** One archive contains the matching ProkCC release config, model files and tables. */
	private static final String PROKCC_ARCHIVE=SOURCEFORGE_URL+"prokcc_v1.5.tar";
	private static final String PROKCC_LEGACY_ARCHIVE=SOURCEFORGE_URL+"prokcc_v1.2.1.tar";
	private static final String BBGAP_SHORT_ARCHIVE="bbgap_short_v2.tar";
	private static final String BBGAP_TAILS_ARCHIVE="bbgap_tails_v2.tar";
	/** The size-filtered 32k DDL sketch DB is ~9.8 GB -- too large for GitHub's 2GB cap; hosted on Zenodo
	 * as a direct-download file link pinned to the v40.00 record (matches the v40.00 GitHub links above). */
	private static final String ZENODO_DDL32K="https://zenodo.org/records/21630308/files/refseqSketchDDL_k25e5b32768.tsv.gz";

	private static final HashMap<String, String> RESOURCE_URLS=new HashMap<>();
	static{
		RESOURCE_URLS.put("magqc_composite_v1.5.graft", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_graft_private_v1.5_18bit.bbnet", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_graft_heads_v1.5_18bit.bbnet", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_composite_v1.4_shrunk_18bit.bbnet.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_subnets_v1.4_shrunk_18bit.bbnets.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("composite_v1.2.1_shrunk_0.8pct_18bit.bbnet.gz", PROKCC_LEGACY_ARCHIVE);
		RESOURCE_URLS.put("magqc_subnets_v1.2.1_shrunk_57.4pct_18bit.bbnets.gz", PROKCC_LEGACY_ARCHIVE);
		RESOURCE_URLS.put("composite_d252_polished_weight18.bbnet.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("v1.bbnets.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_subnets_v1.bbnets.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_subnets_v1.full_fallback.bbnets.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_hbm_v1.rare01.a48.mqhb.bgz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_hbm_v1.rare01.hbmt.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_sidecar_v1.tsv.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("all_euk_18S_best_taxsorted.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("all_prok_16S_best_taxsorted.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b4096.tsv.gz", GITHUB_V4000+"refseqSketchDDL_k25e5b4096.tsv.gz");
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b4096_merged.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b32768.tsv.gz", ZENODO_DDL32K);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b2048.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b2048_merged.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("ribokmers.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("protein_blacklist_reduced16_k5_raw_shadow1024.faa.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("protein_blacklist_d120_k6_nonoverlap_shadow2048.faa.gz", GITHUB_RELEASES);

		RESOURCE_URLS.put("refseqA48_with_ribo.spectra.gz", GITHUB_V4001+"refseqA48_with_ribo.spectra.gz");

		RESOURCE_URLS.put("RQCFilterData.tar", NERSC_URL);
		RESOURCE_URLS.put("riboKmers20fused.fa.gz", NERSC_URL);

		RESOURCE_URLS.put("ssuSketchDDL.tsv.gz", SOURCEFORGE_URL);
	}

}
