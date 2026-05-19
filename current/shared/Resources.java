package shared;

import java.io.File;
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

	public static String find(String fname){
		return find(fname, true);
	}

	/**
	 * Locates a resource file via Data.findPath().
	 * If missing, prints download instructions and returns null.
	 */
	public static String find(String fname, boolean exit){
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

	private static final String GITHUB_V3982="https://github.com/bbushnell/BBTools/releases/tag/v39.82/";
	private static final String GITHUB_V3984="https://github.com/bbushnell/BBTools/releases/tag/v39.84/";
	private static final String NERSC_URL="https://portal.nersc.gov/cfs/bbtools/";
	private static final String SOURCEFORGE_URL="https://sourceforge.net/projects/bbmap/files/Resources/";

	private static final HashMap<String, String> RESOURCE_URLS=new HashMap<>();
	static{
		RESOURCE_URLS.put("all_euk_18S_best_taxsorted.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("all_prok_16S_best_taxsorted.fa.gz", GITHUB_V3982);
		RESOURCE_URLS.put("refseqSketchDDL.tsv.gz", GITHUB_V3982);
		RESOURCE_URLS.put("refseqSketchDDL_merged.tsv.gz", GITHUB_V3982);
		RESOURCE_URLS.put("ribokmers.fa.gz", GITHUB_V3982);

		RESOURCE_URLS.put("refseqA48_with_ribo.spectra.gz", GITHUB_V3984);

		RESOURCE_URLS.put("RQCFilterData.tar", NERSC_URL);
		RESOURCE_URLS.put("riboKmers20fused.fa.gz", NERSC_URL);

		RESOURCE_URLS.put("ssuSketchDDL.tsv.gz", SOURCEFORGE_URL);
	}

}
