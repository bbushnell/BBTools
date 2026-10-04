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
		if(path.startsWith("prokcc/") || path.contains("/prokcc/")){return PROKCC_ARCHIVE;}
		return RESOURCE_URLS.get(new java.io.File(path).getName());
	}

	/** ProkCC's config, models and tables are installed together from one manual download. */
	public static String prokccDownloadInstructions(){
		return "Download prokcc_v1.2.1.tar from:\n  "+PROKCC_ARCHIVE+
			"\nExtract its contents into BBTools/resources/ to create resources/prokcc/, then try again.\n";
	}

	/** Locates a resource file.  If fname contains commas, tries each
	 * candidate in order and returns the first one found. */
	public static String find(String fname){
		return find(fname, true);
	}

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
	 * Only complains if ALL are missing. */
	public static String find(ArrayList<String> fnames){
		return find(fnames, true);
	}

	public static String find(ArrayList<String> fnames, boolean exit){
		for(String f : fnames){
			String path=Data.findPath(f, false);
			if(path!=null){return path;}
		}
		if(fnames.isEmpty()){return null;}
		if(!exit){return null;}
		StringBuilder sb=new StringBuilder();
		boolean prokcc=false;
		sb.append("\nERROR: No resource found from candidates:\n");
		for(String f : fnames){
			String bare=f.startsWith("?") ? f.substring(1) : f;
			String url=downloadURL(bare);
			prokcc|=PROKCC_ARCHIVE.equals(url);
			sb.append("  ").append(bare);
			if(url!=null){sb.append("  (").append(url).append(')');}
			sb.append('\n');
		}
		if(prokcc){sb.append(prokccDownloadInstructions());}
		sb.append("Place other network files in BBTools/networks/ and other resources in BBTools/resources/.\n");
		System.err.print(sb);
		System.exit(1);
		return null;
	}

	/** Locates a single resource file.  If missing, prints download
	 * instructions and optionally exits. */
	private static String findSingle(String fname, boolean exit){
		String path=Data.findPath(fname, false);
		if(path!=null){return path;}

		String bare=fname;
		if(bare.startsWith("?")){bare=bare.substring(1);}

		String url=downloadURL(bare);
		StringBuilder sb=new StringBuilder();
		sb.append("\nERROR: Required resource not found: ").append(bare).append('\n');
		if(PROKCC_ARCHIVE.equals(url)){
			sb.append(prokccDownloadInstructions());
		}else{
			if(url!=null){
				sb.append("Download it from:\n  ").append(url).append('\n');
			}else{
				sb.append("You may need to download it from:\n  ").append(SOURCEFORGE_URL).append('\n');
			}
			final String destination=Data.isNetworkFile(bare) ? "BBTools/networks/" : "BBTools/resources/";
			sb.append("Place it in ").append(destination).append(" and try again.\n");
		}
		System.err.print(sb);
		if(exit){
			System.exit(1);
		}
		return null;
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
	private static final String PROKCC_ARCHIVE=SOURCEFORGE_URL+"prokcc_v1.2.1.tar";
	/** The size-filtered 32k DDL sketch DB is ~9.8 GB -- too large for GitHub's 2GB cap; hosted on Zenodo
	 * as a direct-download file link pinned to the v40.00 record (matches the v40.00 GitHub links above). */
	private static final String ZENODO_DDL32K="https://zenodo.org/records/21630308/files/refseqSketchDDL_k25e5b32768.tsv.gz";

	private static final HashMap<String, String> RESOURCE_URLS=new HashMap<>();
	static{
		RESOURCE_URLS.put("composite_v1.2.1_shrunk_0.8pct_18bit.bbnet.gz", PROKCC_ARCHIVE);
		RESOURCE_URLS.put("magqc_subnets_v1.2.1_shrunk_57.4pct_18bit.bbnets.gz", PROKCC_ARCHIVE);
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

		RESOURCE_URLS.put("refseqA48_with_ribo.spectra.gz", GITHUB_V4001+"refseqA48_with_ribo.spectra.gz");

		RESOURCE_URLS.put("RQCFilterData.tar", NERSC_URL);
		RESOURCE_URLS.put("riboKmers20fused.fa.gz", NERSC_URL);

		RESOURCE_URLS.put("ssuSketchDDL.tsv.gz", SOURCEFORGE_URL);
	}

}
