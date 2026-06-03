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
		sb.append("\nERROR: No resource found from candidates:\n");
		for(String f : fnames){
			String bare=f.startsWith("?") ? f.substring(1) : f;
			String url=RESOURCE_URLS.get(bare);
			sb.append("  ").append(bare);
			if(url!=null){sb.append("  (").append(url).append(')');}
			sb.append('\n');
		}
		sb.append("Download one and place it in BBTools/resources/\n");
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

		String url=RESOURCE_URLS.get(bare);
		StringBuilder sb=new StringBuilder();
		sb.append("\nERROR: Required resource not found: ").append(bare).append('\n');
		if(url!=null){
			sb.append("Download it from:\n  ").append(url).append('\n');
		}else{
			sb.append("You may need to download it from:\n  ").append(SOURCEFORGE_URL).append('\n');
		}
		sb.append("Place it in BBTools/resources/ and try again.\n");
		System.err.print(sb);
		if(exit){
			System.exit(1);
		}
		return null;
	}

	private static final String GITHUB_RELEASES="https://github.com/bbushnell/BBTools/releases/";
	private static final String GITHUB_V3982="https://github.com/bbushnell/BBTools/releases/tag/v39.82/";
	private static final String GITHUB_V3984="https://github.com/bbushnell/BBTools/releases/tag/v39.84/";
	private static final String NERSC_URL="https://portal.nersc.gov/cfs/bbtools/";
	private static final String SOURCEFORGE_URL="https://sourceforge.net/projects/bbmap/files/Resources/";

	private static final HashMap<String, String> RESOURCE_URLS=new HashMap<>();
	static{
		RESOURCE_URLS.put("all_euk_18S_best_taxsorted.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("all_prok_16S_best_taxsorted.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b4096.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b4096_merged.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b2048.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("refseqSketchDDL_k25e5b2048_merged.tsv.gz", GITHUB_RELEASES);
		RESOURCE_URLS.put("ribokmers.fa.gz", GITHUB_V3982);

		RESOURCE_URLS.put("refseqA48_with_ribo.spectra.gz", GITHUB_V3984);

		RESOURCE_URLS.put("RQCFilterData.tar", NERSC_URL);
		RESOURCE_URLS.put("riboKmers20fused.fa.gz", NERSC_URL);

		RESOURCE_URLS.put("ssuSketchDDL.tsv.gz", SOURCEFORGE_URL);
	}

}
