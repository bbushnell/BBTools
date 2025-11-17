package prok;

import java.util.ArrayList;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.TextStreamWriter;
import server.ServerTools;
import shared.Parse;
import shared.Tools;
import template.ThreadWaiter;

/** Crawls ncbi's ftp site to download genomes and annotations */
public class FetchProks {
	
	/**
	 * Program entry point for bacterial genome downloading.
	 * Parses command-line arguments for base address, output file, species limits,
	 * and quality selection.
	 * Creates multiple processing threads to handle genus-based parallel crawling
	 * of NCBI FTP directories.
	 *
	 * @param args Command-line arguments: [baseAddress] [output] [maxSpeciesPerGenus] [findBest]
	 */
	public static void main(String[] args){
		//ftp://ftp.ncbi.nih.gov:21/genomes/refseq/bacteria/
		
		String baseAddress=args[0];
		String out=args.length>1 ? args[1] : "stdout";
		if(args.length>2){
			maxSpeciesPerGenus=Integer.parseInt(args[2]);
			System.err.println("Set maxSpeciesPerGenus="+maxSpeciesPerGenus);
		}
		if(args.length>3){
			findBest=Parse.parseBoolean(args[3]);
			System.err.println("Set findBest="+findBest);
		}
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false, FileFormat.TEXT);
		tsw.start();

//		iterateOuter(baseAddress, tsw);
		ArrayList<String> contents=ServerTools.listDirectory(baseAddress, retries);
		
		int threads=7;
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(contents, tsw, i, threads));
		}
		for(ProcessThread pt : alpt){pt.start();}
		boolean success=ThreadWaiter.waitForThreadsToFinish(alpt);
		
		for(ProcessThread pt : alpt){
			totalSpecies+=pt.totalSpeciesT;
			totalGenus+=pt.totalGenusT;
			totalGenomes+=pt.totalGenomesT;
		}
		System.err.println("Total Genomes: "+totalGenomes);
		System.err.println("Total Species: "+totalSpecies);
		System.err.println("Total Genuses: "+totalGenus);
		
		tsw.poisonAndWait();
		assert(success);
	}
	
	/**
	 * Worker thread for processing bacterial species directories from NCBI FTP site.
	 * Each thread handles complete genera to avoid synchronization issues and maintains
	 * separate tracking of processed species and genomes.
	 * @author Brian Bushnell
	 */
	static class ProcessThread extends Thread {
		
		/**
		 * Constructs a processing thread for bacterial genome crawling.
		 *
		 * @param speciesList_ List of species directories to process
		 * @param tsw_ Output writer for generating download commands
		 * @param tid_ Thread identifier for work distribution
		 * @param threads_ Total number of processing threads
		 */
		ProcessThread(ArrayList<String> speciesList_, TextStreamWriter tsw_, int tid_, int threads_){
			speciesList=speciesList_;
			tsw=tsw_;
			tid=tid_;
			threads=threads_;
		}
		
		@Override
		public void run(){
			for(String s : speciesList){
//				if((s.hashCode()&Integer.MAX_VALUE)%threads==tid) {
//					processSpecies(s);
//				}
				
				//This way one thread handles an entire genus
				if(s!=null){
					String genus=getGenus(s);
					if(genus!=null){
						if((genus.hashCode()&Integer.MAX_VALUE)%threads==tid) {
							processSpecies(s);
						}
					}else{
						if((s.hashCode()&Integer.MAX_VALUE)%threads==tid) {
							processSpecies(s);
						}
					}
				}
			}
		}
		
		/**
		 * Processes a single bacterial species directory for genome downloading.
		 * Checks genus limits and delegates to examineSpecies for actual file discovery.
		 * Updates thread-local counters for species, genera, and genomes processed.
		 * @param species Path to the species directory on NCBI FTP
		 */
		void processSpecies(String species){
			String genus=getGenus(species);
			if(genus!=null){
				final int count=seen(genus, seen);
				
				if(maxSpeciesPerGenus<1 || count<maxSpeciesPerGenus){
					int found=examineSpecies(species, tsw);
					if(found>=1){
						totalSpeciesT++;
						totalGenomesT+=found;
						if(count==0){totalGenusT++;}
						put(genus, found, seen);
					}
				}else{
					if(verbose){System.err.println("same genus: "+species+"\n"+genus);}
				}
			}else{
				if(verbose){System.err.println("bad species: "+species+"\n"+genus);}
			}
		}
		
		/** List of species directories assigned to this thread for processing */
		final ArrayList<String> speciesList;
		/** Thread identifier for work distribution hashing */
		final int tid;
		/** Total number of processing threads for modulo distribution */
		final int threads;
		//This is OK now that threads work on a per-genus basis
		/** Thread-local map tracking species count per genus to enforce limits */
		HashMap<String, Integer> seen=new HashMap<String, Integer>();
		/** Output writer for generating download commands */
		final TextStreamWriter tsw;
		
		/** Thread-local counter of successfully processed species */
		int totalSpeciesT=0;
		/** Thread-local counter of distinct genera processed */
		int totalGenusT=0;
		/** Thread-local counter of genome files discovered for download */
		int totalGenomesT=0;
	}
	
	/**
	 * Extracts genus name from NCBI species directory path.
	 * Handles Candidatus organism naming conventions by removing the prefix.
	 * @param path Full path to species directory
	 * @return Genus name extracted from path, or null if parsing fails
	 */
	static String getGenus(String path){
		//Candidatus_Hamiltonella
		String name=path.substring(path.lastIndexOf('/')+1);
		if(name.startsWith("Candidatus_")){name=name.substring("Candidatus_".length());}
		int under=name.indexOf('_');
		if(under>0){
			return name.substring(0, under);
		}else{
			return null;
		}
	}
	
	/**
	 * Extracts species name from NCBI directory path.
	 * Handles Candidatus organism naming by removing the prefix.
	 * @param path Full path to species directory
	 * @return Species name extracted from directory path
	 */
	static String getSpecies(String path){
		//Candidatus_Hamiltonella
		String name=path.substring(path.lastIndexOf('/')+1);
		if(name.startsWith("Candidatus_")){name=name.substring("Candidatus_".length());}
		return name;
	}
	
	/**
	 * Examines a bacterial species directory for genome assemblies.
	 * Searches in priority order: reference genomes, latest assemblies, then all assemblies.
	 * Stops after finding the first suitable assembly to avoid duplicates.
	 *
	 * @param baseAddress Base FTP address for the species
	 * @param tsw Output writer for download commands
	 * @return Number of genomes found and processed
	 */
	static int examineSpecies(String baseAddress, TextStreamWriter tsw){
		if(verbose){System.err.println("examineSpecies: "+baseAddress);}
		String speciesName=getSpecies(baseAddress);
		ArrayList<String> contents=ServerTools.listDirectory(baseAddress, retries);
//		System.err.println("B: "+contents);
		int found=0;
		for(String s : contents){
//			System.err.println(s);
			if(s.contains("reference")){
//				System.err.println("Looking at '"+s+"'");
				found+=examineAssemblies(s, tsw, speciesName);
			}
		}
		if(found>0){return found;}
		for(String s : contents){
//			System.err.println(s);
			 if(s.contains("latest_assembly_versions")){
//				System.err.println("Looking at '"+s+"'");
				 found+=examineAssemblies(s, tsw, speciesName);
			}
		}
		if(found>0){return found;}
		for(String s : contents){
//			System.err.println(s);
			if(s.contains("all_assembly_versions")){
//				System.err.println("Looking at '"+s+"'");
				found+=examineAssemblies(s, tsw, speciesName);
			}
		}
		return found;
	}
	
	/**
	 * Examines assembly directory for downloadable genome files.
	 * Optionally finds best assembly by quality metrics, otherwise processes first suitable assembly.
	 *
	 * @param baseAddress FTP address of assembly directory
	 * @param tsw Output writer for generating download commands
	 * @param speciesName Species name for file naming
	 * @return Number of assemblies processed (0 or 1)
	 */
	static int examineAssemblies(String baseAddress, TextStreamWriter tsw, String speciesName){
		if(verbose){System.err.println("examineAssemblies: "+baseAddress);}
		Stats stats=null;
		if(findBest){
			stats=findBestAssembly(baseAddress);
			if(stats!=null){
				stats.name=speciesName;
				int x=examineAssembly(stats, tsw, speciesName);
				if(x>0){return x;}
			}
		}
		
		ArrayList<String> contents=ServerTools.listDirectory(baseAddress, retries);
//		System.err.println("C: "+contents);
		
		int found=0;
		for(String s : contents){
			stats=calcStats(s);
			if(stats!=null){
				stats.name=speciesName;
				found+=examineAssembly(stats, tsw, speciesName);
				if(found>0){break;}
			}
		}
		return found;
	}
	
	/** Tries to find the assembly with the longest contig */
	static Stats findBestAssembly(String baseAddress){
		if(verbose){System.err.println("findBestAssembly: "+baseAddress);}
		ArrayList<String> contents=ServerTools.listDirectory(baseAddress, retries);
//		System.err.println("C: "+contents);
		Stats best=null;
		for(String s : contents){
//			System.err.println(s);
			Stats stats=calcStats(s);
			if(stats!=null){
				if(best==null || stats.compareTo(best)>0){
					best=stats;
				}
			}
		}
		return best;
	}
	
	/**
	 * Calculates assembly statistics by parsing NCBI assembly report.
	 * Extracts contig count, genome size, longest contig, and taxonomic ID
	 * from the assembly_report.txt file.
	 *
	 * @param baseAddress FTP path to assembly directory
	 * @return Stats object with assembly metrics, or null if report unavailable
	 */
	static Stats calcStats(String baseAddress){
		if(verbose){System.err.println("calcStats: "+baseAddress);}
		ArrayList<String> contents=ServerTools.listDirectory(baseAddress, retries);
		String report=null;
		for(String s : contents){
			if(s.endsWith("_assembly_report.txt")){
				report=s;
				break;
			}
		}
		if(report==null){
			if(verbose){System.err.println("Could not find report for "+baseAddress);}
			return null;
		}
		if(verbose){System.err.println("Report: "+report);}
		ArrayList<String> data=null;
		for(int i=0; i<=retries && data==null; i++){
			try {
				data = ServerTools.readFTPFile(report);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				try {
					Thread.sleep(Tools.mid(10000, i*1000, 1000));
				} catch (InterruptedException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
		}
		if(data==null){return null;}
		int contigs=0;
		long size=0;
		long max=0;
		int taxid=-1;
		for(String s : data){
			if(s!=null && s.length()>0){
				if(s.charAt(0)=='#'){
					if(s.startsWith("# Taxid:")){
						String[] split=Tools.whitespacePlus.split(s);
						try {
							taxid=Integer.parseInt(split[split.length-1]);
						} catch (NumberFormatException e) {
							e.printStackTrace();
						}
						assert(taxid>-1) : "Bad TaxID: '"+s+"'";
					}
				}else{
					String[] split=s.split("\t");
					contigs++;
					long len;
					try {
						len=Long.parseLong(split[8]);
					} catch (NumberFormatException e) {
						len=1;
					}
					size+=len;
					max=Tools.max(max, len);
				}
			}
		}
		return new Stats(baseAddress, max, size, contigs, taxid);
	}
	
	/**
	 * Examines individual assembly for downloadable genomic and annotation files.
	 * Generates wget commands for both .fna.gz (genomic sequence) and .gff.gz (annotation) files.
	 * Supports optional file renaming and sequence ID standardization.
	 *
	 * @param stats Assembly statistics and path information
	 * @param tsw Output writer for download commands
	 * @param speciesName Species name for output file naming
	 * @return 1 if both files found and processed, 0 otherwise
	 */
	static int examineAssembly(Stats stats, TextStreamWriter tsw, String speciesName){
		if(verbose){System.err.println("examineAssembly: "+stats.path);}
		ArrayList<String> contents=ServerTools.listDirectory(stats.path, retries);
//		System.err.println("D: "+contents);
		String gff=null;
		String fna=null;
		for(String s : contents){
//			System.err.println(s);
			if(!s.contains("_from_genomic")){
				if(s.endsWith("genomic.fna.gz")){fna=s;}
				else if(s.endsWith("genomic.gff.gz")){gff=s;}
			}
		}
		if(fna!=null && gff!=null){
			System.err.println("Printing: "+fna);
			String prefix=(tidInFilename ? "tid_"+stats.taxID+"_" : "");
			
			synchronized(tsw){
				if(renameSequences){
					tsw.println("wget -q -O - "+fna+" | "
							+ "gi2taxid.sh in=stdin.fa.gz deleteinvalid zl=9 server -Xmx1g out="+prefix+speciesName+".fna.gz");
					tsw.println("wget -q -O - "+gff+" | "
							+ "gi2taxid.sh in=stdin.gff.gz deleteinvalid zl=9 server -Xmx1g out="+prefix+speciesName+".gff.gz");
				}else if(renameFiles){
					tsw.println("wget -q -O - "+fna+" > "+prefix+speciesName+".fna.gz");
					tsw.println("wget -q -O - "+gff+" > "+prefix+speciesName+".gff.gz");
				}else{
					tsw.println("wget -q "+fna);
					tsw.println("wget -q "+gff);
				}
				tsw.println();
			}
			return 1;
		}
		return 0;
	}
	
	/**
	 * Constructs FTP sub-address by combining base address with extension path.
	 * @param baseAddress Base FTP directory address
	 * @param extension Path extension to append
	 * @return Complete FTP address for subdirectory
	 */
	static String makeSubAddress(String baseAddress, String extension){
		if(!baseAddress.endsWith("/")){baseAddress=baseAddress+"/";}
		String subAddress=baseAddress+extension.substring(extension.indexOf('/')+1);
		return subAddress;
	}
	
	/**
	 * Gets count of previously processed items for the given key.
	 * @param s Key to look up in the tracking map
	 * @param map HashMap containing processing counts
	 * @return Count of previous occurrences, or 0 if not seen
	 */
	static int seen(String s, HashMap<String, Integer> map){
//		synchronized(map){
			Integer x=map.get(s);
			return x==null ? 0 : x.intValue();
//		}
	}
	/**
	 * Increments the count for a key in the tracking map.
	 * Adds the found count to any existing value for the key.
	 *
	 * @param s Key to update in the map
	 * @param found Number to add to existing count
	 * @param map HashMap containing processing counts
	 */
	static void put(String s, int found, HashMap<String, Integer> map){
//		synchronized(map){
			int present=seen(s, map);
			map.put(s, present+found);
//		}
	}
	
	/**
	 * Assembly statistics container for comparing genome quality.
	 * Implements comparison logic prioritizing assemblies with valid taxonomy,
	 * larger genome size, longer contigs, and fewer fragments.
	 * @author Brian Bushnell
	 */
	static class Stats implements Comparable<Stats>{
		
		/**
		 * Constructs assembly statistics from parsed assembly report data.
		 *
		 * @param path_ FTP path to the assembly directory
		 * @param maxContig_ Length of longest contig in the assembly
		 * @param size_ Total genome size in bases
		 * @param contigs_ Number of contigs in the assembly
		 * @param taxID_ NCBI taxonomic identifier
		 */
		public Stats(String path_, long maxContig_, long size_, int contigs_, int taxID_){
			path=path_;
			maxContig=maxContig_;
			size=size_;
			contigs=contigs_;
			taxID=taxID_;
		}

		@Override
		public int compareTo(Stats b) {//true if b is better
			if(b==null){return 1;}
			if(taxID>0 && b.taxID<1){return 1;}
			if(b.taxID>0 && taxID<1){return -1;}
			
			if(size>2*b.size){return 1;}
			if(size<2*b.size){return -1;}

			if(maxContig>b.maxContig){return 1;}
			if(maxContig<b.maxContig){return -1;}
			
			return b.contigs-contigs;
		}
		
		/** FTP path to the assembly directory */
		String path;
		/** Species name assigned for output file naming */
		String name;
		/** Length of the longest contig in this assembly */
		long maxContig;
		/** Total genome size in bases */
		long size;
		/** Total number of contigs in the assembly */
		int contigs;
		/** NCBI taxonomic identifier for this organism */
		int taxID;
	}
	
	/** Enable verbose output during processing */
	static boolean verbose=true;
//	static boolean allowSameGenus=false;
	/** Maximum number of species to download per genus for balanced sampling */
	static int maxSpeciesPerGenus=1;
	/** Whether to rename downloaded files using species names */
	static boolean renameFiles=true;
	/** Whether to standardize sequence headers using gi2taxid tool */
	static boolean renameSequences=true;
	/** Number of retry attempts for failed FTP operations */
	static int retries=40;
	/**
	 * Whether to search for best assembly by quality metrics rather than first found
	 */
	static boolean findBest=false;
	
	/** Whether to include taxonomic ID prefix in output filenames */
	static boolean tidInFilename=true;
	
//	private static HashMap<String, Integer> seen=new HashMap<String, Integer>();
	
	/** Global counter of successfully processed species across all threads */
	static int totalSpecies=0;
	/** Global counter of distinct genera processed across all threads */
	static int totalGenus=0;
	/** Global counter of genome files discovered for download across all threads */
	static int totalGenomes=0;

	/** Cached Integer object for the value 1 */
	private static final Integer one=1;
	
}
