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
	
	static class ProcessThread extends Thread {
		
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
		
		final ArrayList<String> speciesList;
		final int tid;
		final int threads;
		//This is OK now that threads work on a per-genus basis
		HashMap<String, Integer> seen=new HashMap<String, Integer>();
		final TextStreamWriter tsw;
		
		int totalSpeciesT=0;
		int totalGenusT=0;
		int totalGenomesT=0;
	}
	
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
	
	static String getSpecies(String path){
		//Candidatus_Hamiltonella
		String name=path.substring(path.lastIndexOf('/')+1);
		if(name.startsWith("Candidatus_")){name=name.substring("Candidatus_".length());}
		return name;
	}
	
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
	
	static String makeSubAddress(String baseAddress, String extension){
		if(!baseAddress.endsWith("/")){baseAddress=baseAddress+"/";}
		String subAddress=baseAddress+extension.substring(extension.indexOf('/')+1);
		return subAddress;
	}
	
	static int seen(String s, HashMap<String, Integer> map){
//		synchronized(map){
			Integer x=map.get(s);
			return x==null ? 0 : x.intValue();
//		}
	}
	static void put(String s, int found, HashMap<String, Integer> map){
//		synchronized(map){
			int present=seen(s, map);
			map.put(s, present+found);
//		}
	}
	
	static class Stats implements Comparable<Stats>{
		
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
		
		String path;
		String name;
		long maxContig;
		long size;
		int contigs;
		int taxID;
	}
	
	static boolean verbose=true;
//	static boolean allowSameGenus=false;
	static int maxSpeciesPerGenus=1;
	static boolean renameFiles=true;
	static boolean renameSequences=true;
	static int retries=40;
	static boolean findBest=false;
	
	static boolean tidInFilename=true;
	
//	private static HashMap<String, Integer> seen=new HashMap<String, Integer>();
	
	static int totalSpecies=0;
	static int totalGenus=0;
	static int totalGenomes=0;

	private static final Integer one=1;
	
}
