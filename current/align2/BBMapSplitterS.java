package align2;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;

import static align2.BBSplitter.AMBIGUOUS2_MODE;
import static align2.BBSplitter.AMBIGUOUS2_FIRST;
import static align2.BBSplitter.AMBIGUOUS2_UNSET;
import static align2.BBSplitter.AMBIGUOUS2_SPLIT;
import static align2.BBSplitter.AMBIGUOUS2_ALL;
import static align2.BBSplitter.AMBIGUOUS2_RANDOM;
import static align2.BBSplitter.AMBIGUOUS2_TOSS;
import static align2.BBSplitter.TRACK_SET_STATS;
import static align2.BBSplitter.TRACK_SCAF_STATS;
import static align2.BBSplitter.setCountTable;
import static align2.BBSplitter.scafCountTable;
import static align2.BBSplitter.overwrite;
import static align2.BBSplitter.append;
import static align2.BBSplitter.getSets;
import static align2.BBSplitter.getScaffolds;
import static align2.BBSplitter.toListNames;
import align2.BBSplitter.SetCount;

/**
 * Writer-backed splitter output for BBMapS. Routing, scaffold/set counting and
 * BAM-script generation are copied from BBSplitter; configuration and SetCount
 * tables remain the same objects used by AbstractMapper parsing and reporting.
 * Only output tables and their lifecycle belong to this class. Reference
 * preparation and classic BBSplitter remain unchanged.
 * @author Brian Bushnell, Collei
 */
public final class BBMapSplitterS {
	public static synchronized HashMap<String, Writer> makeOutputStreams(String[] args, boolean OUTPUT_READS, boolean OUTPUT_ORDERED_READS,
			int buff, boolean paired, boolean overwrite_, boolean append_, boolean ambiguous){
//		assert(false) : Arrays.toString(args);
		HashMap<String, Writer> table=new HashMap<String, Writer>();
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0];
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
//			assert(b!=null) : "Bad parameter: "+arg+"\n"+Arrays.toString(args);
			
			if(arg.indexOf('=')>0 && a.toLowerCase().startsWith("out_")){
				assert(b!=null) : "Bad parameter: "+arg+"\n"+Arrays.toString(args);
				String name=a.substring(4).replace('\\', '/');
				
				final String fname1, fname2;
				
				if(ambiguous){
					if(b.indexOf('/')>=0){
						int x=b.lastIndexOf('/');
						b=b.substring(0, x+1)+"AMBIGUOUS_"+b.substring(x+1);
					}else{
						b="AMBIGUOUS_"+b;
					}
				}
				
				if(!FileFormat.hasSamOrBamExtension(b) && ReadWrite.stripExtension(b).contains("#")){
					fname1=b.replace('#', '1');
					fname2=b.replace('#', '2');
				}else{
					fname1=b;
					fname2=null;
				}
//				assert(false) : fname1;
//				assert(!ambiguous) : fname1+", "+fname2+", "+b+", "+ambiguous;

				FileFormat ff1=FileFormat.testOutput(fname1, FileFormat.SAM, null, true, overwrite_, append_, OUTPUT_ORDERED_READS);
				FileFormat ff2=paired ? FileFormat.testOutput(fname2, FileFormat.SAM, null, true, overwrite_, append_, OUTPUT_ORDERED_READS) : null;
				Writer ros=WriterFactory.getStream(ff1, ff2, null, null, buff, null, false, Shared.threads());
				ros.start();
//				Data.sysout.println("Started output stream:\t"+t);
				table.put(name, ros);
				AbstractMapThread.OUTPUT_SAM|=ff1.samOrBam();
			}
		}
		return table.isEmpty() ? null : table;
	}
	/**
	 * @param readlist List of reads to print
	 * @param listID ID of read list, from ReadInputStream
	 * @param splitTable A temporary structure to hold sets of reads that go to the different output streams
	 * @param clearzone Min distance between best and next-best site to be considered unambiguous
	 */
	public static void printReads(ArrayList<Read> readlist, long listID, HashMap<String, ArrayList<Read>> splitTable, int clearzone){
		if(clearzone>=0 || TRACK_SET_STATS || TRACK_SCAF_STATS){
			printReadsAndProcessAmbiguous(readlist, listID, splitTable, null, clearzone);
			return;
		}
		assert((streamTable!=null && streamTable.size()>0) || (setCountTable!=null && setCountTable.size()>0) || (scafCountTable!=null && scafCountTable.size()>0));
		boolean clear=true;
		if(splitTable==null){
			splitTable=new HashMap<String, ArrayList<Read>>();
			clear=false;
		}
		
		if(!readlist.isEmpty()){
			HashSet<String> set=new HashSet<String>(8);
			for(Read r : readlist){
				if(r!=null){
					set=toListNames(r, set);
					for(String s : set){
						ArrayList<Read> alr=splitTable.get(s);
						if(alr==null){
							alr=new ArrayList<Read>();
							splitTable.put(s, alr);
						}
						//123***
						alr.add(r);
					}
					set.clear();
				}
			}
		}
		
		for(String s : streamTable.keySet()){
			ArrayList<Read> alr=splitTable.get(s);
			if(alr==null){alr=blank;}
			Writer tros=streamTable.get(s);
			tros.add(alr, listID);
		}
		if(clear){splitTable.clear();}
	}
	

	/**
	 * @param readlist List of reads to print
	 * @param listID ID of read list, from ReadInputStream
	 * @param splitTable A temporary structure to hold sets of reads that go to the different output streams
	 * @param clearzone Min distance between best and next-best site to be considered unambiguous
	 */
	public static void printReadsAndProcessAmbiguous(ArrayList<Read> readlist, long listID, HashMap<String, ArrayList<Read>> splitTable,
			HashMap<String, ArrayList<Read>> splitTableA, int clearzone){
		assert(clearzone>=0 || TRACK_SET_STATS || TRACK_SCAF_STATS);
		assert((streamTable!=null && streamTable.size()>0) || (setCountTable!=null && setCountTable.size()>0) || (scafCountTable!=null && scafCountTable.size()>0));
		boolean clear=streamTable!=null, clearA=streamTableAmbiguous!=null;
		if(splitTable==null && streamTable!=null){
			splitTable=new HashMap<String, ArrayList<Read>>();
			clear=false;
		}
		if(splitTableA==null && streamTableAmbiguous!=null){
			splitTableA=new HashMap<String, ArrayList<Read>>();
			clearA=false;
		}
		
		final HashSet<String> hss0, hss1, hss2, hss3, hsspr, hssam;
		final HashSet<String>[] hssa;
		if(TRACK_SET_STATS || streamTable!=null){
			hss0=new HashSet<String>(16);
			hss1=new HashSet<String>(16);
			hss2=new HashSet<String>(16);
			hss3=new HashSet<String>(16);
			hsspr=new HashSet<String>(16);
			hssam=new HashSet<String>(16);
			hssa=(HashSet<String>[])new HashSet[] {hss0, hss1, hss2, hss3};
		}else if(TRACK_SCAF_STATS){
			hss0=new HashSet<String>(16);
			hss1=null; hss2=null; hss3=null; hsspr=null; hssam=null; hssa=null;
		}else{
			hss0=null; hss1=null; hss2=null; hss3=null; hsspr=null; hssam=null; hssa=null;
		}
		
		for(final Read r1 : readlist){
//			System.out.println("\nProcessing read "+r1.numericID);
			final Read r2=r1==null ? null : r1.mate;
			
			if(r1!=null){addToScafCounts(r1, clearzone, hss0);} //Scafstats for read 1
			if(r2!=null){addToScafCounts(r2, clearzone, hss0);} //Scafstats for read 2
			
			if(r1!=null){

				final HashSet<String>[] sets=(TRACK_SET_STATS || streamTable!=null) ? getSets(r1, clearzone, hssa) : null;
				boolean ambiguous=false;
				if(sets!=null){
					final HashSet<String> p1=(sets[0].isEmpty() ? null : sets[0]), s1=(sets[1].isEmpty() ? null : sets[1]),
							p2=(sets[2].isEmpty() ? null : sets[2]), s2=(sets[3].isEmpty() ? null : sets[3]);
					assert(sets==hssa);
//					assert(p1!=null);
//					assert(s1!=null);
//					assert(p2!=null);
//					assert(s2!=null);

					if(p1!=null && p2!=null && !p1.equals(p2)){ambiguous=true;}
					else if(p1!=null && s1!=null && !p1.containsAll(s1)){ambiguous=true;}
					else if(p2!=null && s2!=null && !p2.containsAll(s2)){ambiguous=true;}

//					System.out.println("\nambiguous="+ambiguous);
//					System.out.println(p1);
//					System.out.println(s1);
					
					HashSet<String> primarySet=hsspr, ambigSet=hssam;
					primarySet.clear();
					ambigSet.clear();
					if(AMBIGUOUS2_MODE==AMBIGUOUS2_FIRST || AMBIGUOUS2_MODE==AMBIGUOUS2_UNSET){//pick one
						if(r2==null || r1.mapScore>=r2.mapScore){
							if(p1!=null){primarySet.addAll(p1);}
						}else{
							if(p2!=null){primarySet.addAll(p2);}
						}
					}else{//merge
						if(p1!=null){primarySet.addAll(p1);}
						if(p2!=null){primarySet.addAll(p2);}
					}
					

					if(ambiguous){
						if(AMBIGUOUS2_MODE==AMBIGUOUS2_SPLIT){
							if(primarySet!=null && s1!=null){primarySet.addAll(s1);}
							if(primarySet!=null && s2!=null){primarySet.addAll(s2);}
							ambigSet=primarySet;
							primarySet=null;
						}else if(AMBIGUOUS2_MODE==AMBIGUOUS2_ALL){
							if(primarySet!=null && s1!=null){primarySet.addAll(s1);}
							if(primarySet!=null && s2!=null){primarySet.addAll(s2);}
							ambigSet=null;
						}else if(AMBIGUOUS2_MODE==AMBIGUOUS2_RANDOM){
							throw new RuntimeException("AMBIGUOUS2_RANDOM: Not yet implemented.");
						}else if(AMBIGUOUS2_MODE==AMBIGUOUS2_TOSS){
							primarySet=null;
						}
					}
					
					if(primarySet!=null && splitTable!=null){
						for(String s : primarySet){
							ArrayList<Read> alr=splitTable.get(s);
							if(alr==null){
								alr=new ArrayList<Read>();
								splitTable.put(s, alr);
							}
							alr.add(r1);
						}
					}

					if(ambigSet!=null && splitTableA!=null){
						for(String s : ambigSet){
							ArrayList<Read> alr=splitTableA.get(s);
							if(alr==null){
								alr=new ArrayList<Read>();
								splitTableA.put(s, alr);
							}
							alr.add(r1);
						}
					}
					
					if(setCountTable!=null){
						
						primarySet=hsspr;
						primarySet.clear();
						if(p1!=null){primarySet.addAll(p1);}
						if(p2!=null){primarySet.addAll(p2);}
						if(ambiguous){
							if(s1!=null){primarySet.addAll(s1);}
							if(s2!=null){primarySet.addAll(s2);}
						}
						//	System.out.println(primarySet);
						final int incrR=r1.pairCount();
						final int incrB=r1.pairLength();

						int num=0;
						for(String s : primarySet){
							SetCount sc=setCountTable.get(s);
							assert(sc!=null) : s;
							if(ambiguous){
								synchronized(sc){
									//										System.out.println("Incrementing set "+sc);
									sc.ambiguousReads+=incrR;
									sc.ambiguousBases+=incrB;
									if(num==0){
										sc.assignedReads+=incrR;
										sc.assignedBases+=incrB;
									}
								}
							}else{
								synchronized(sc){
									//										System.out.println("Incrementing set "+sc);
									sc.mappedReads+=incrR;
									sc.mappedBases+=incrB;
									if(num==0){
										sc.assignedReads+=incrR;
										sc.assignedBases+=incrB;
									}
								}
							}
							num++;
						}
					}
					for(HashSet<String> set : sets){set.clear();}
				}
			}
		}
		if(streamTable!=null){
			for(String s : streamTable.keySet()){
//				System.err.println("Searching for "+s+" in "+splitTable.keySet());
//				System.err.println(splitTable.containsKey(s));
				ArrayList<Read> alr=splitTable.get(s);
//				System.err.println("Adding alr "+alr+"\n");
				if(alr==null){alr=blank;}
				Writer tros=streamTable.get(s);
				tros.add(alr, listID);
			}
		}
		if(streamTableAmbiguous!=null){
			for(String s : streamTableAmbiguous.keySet()){
				ArrayList<Read> alr=splitTableA.get(s);
				if(alr==null){alr=blank;}
				Writer tros=streamTableAmbiguous.get(s);
				tros.add(alr, listID);
			}
		}
		if(clear){splitTable.clear();}
		if(clearA){splitTableA.clear();}
	}
	
	/**
	 * Updates scaffold-level read count statistics.
	 * Counts mapped, ambiguous, and assigned reads per scaffold.
	 *
	 * @param r Read to count
	 * @param clearzone Score difference threshold for ambiguity determination
	 * @param hss0 Reusable set for scaffold name collection
	 */
	private static void addToScafCounts(Read r, int clearzone, HashSet<String> hss0){
		if(r==null || !r.mapped()){return;}
		assert((scafCountTable!=null)==TRACK_SCAF_STATS) : TRACK_SCAF_STATS;
		if(scafCountTable!=null){
			HashSet<String> set=getScaffolds(r, clearzone, hss0, false);
			if(set!=null && !set.isEmpty()){
				int incrRM=0;
				int incrRA=0;
				int incrBM=0;
				int incrBA=0;

				int incrRS=1+(r.mate!=null && !r.mateMapped() ? 1 : 0);
				int incrBS=r.length()+(r.mate!=null && !r.mateMapped() ? r.mateLength() : 0);
				
				
				if(r.ambiguous()){
					incrRA+=1;
					incrBA+=r.length();
					if(r.mate!=null && !r.mateMapped()){
						incrRA++;
						incrBA+=r.mateLength();
					}
				}else{
					incrRM+=1;
					incrBM+=r.length();
				}
				int num=0;
				for(String s : set){
					SetCount sc=scafCountTable.get(s);
					assert(sc!=null) : "Can't find "+s+"\nin\n"+scafCountTable.keySet()+"\n";

//					System.out.println(sc);
//					System.out.println("+ "+incrRM+", "+incrRA+", "+incrBM+", "+incrBA);
					synchronized(sc){
						//							System.out.println("Incrementing scaf "+sc);
						sc.mappedReads+=incrRM;
						sc.mappedBases+=incrBM;
						sc.ambiguousReads+=incrRA;
						sc.ambiguousBases+=incrBA;
						if(num==0){
							sc.assignedReads+=incrRS;
							sc.assignedBases+=incrBS;
						}
					}
//					System.out.println(sc);
//					System.out.println();
//					assert(false) : "\n"+incrRM+", "+incrRA+", "+incrBM+", "+incrBA+"\n"+set;
					num++;
				}
				set.clear();
			}
		}
	}
	/**
	 * Creates a shell script for converting SAM files to sorted, indexed BAM files.
	 * Generates samtools commands for each SAM/BAM file in the output streams.
	 *
	 * @param outname Output script file name
	 * @param list Additional SAM/BAM files to include
	 * @param sams Variable arguments of SAM/BAM file names
	 */
	public static void makeBamScript(String outname, ArrayList<String> list, String...sams){
		LinkedHashSet<String> set=new LinkedHashSet<String>();
		if(sams!=null){
			for(String s : sams){
				if(s!=null && (s.endsWith(".sam") || s.endsWith(".sam.gz") || s.endsWith(".bam"))){
					set.add(s);
				}
			}
		}
		if(list!=null){
			for(String s : list){
				if(s!=null && (s.endsWith(".sam") || s.endsWith(".sam.gz") || s.endsWith(".bam"))){
					set.add(s);
				}
			}
		}
		if(streamTable!=null){
			for(Writer ros : streamTable.values()){
				String s=ros.fname();
				if(s.endsWith(".sam") || s.endsWith(".sam.gz") || s.endsWith(".bam")){
					set.add(s);
				}
			}
		}
		TextStreamWriter tsw=new TextStreamWriter(outname, overwrite, append, false);
		tsw.start();
		
		String memstring=null;
		if(set.size()>0){
			tsw.println("#!/bin/bash");
			
			long mem=Runtime.getRuntime().maxMemory()/3400000;
			mem=Tools.min(100000, mem);
			if(mem<2048){memstring=mem+"M";}
			else{memstring=(mem/1024)+"G";}

			tsw.println("echo \"Note: This script is designed to run with the amount of memory detected by BBMap.\"");
			tsw.println("echo \"      If Samtools crashes, please ensure you are running on the same platform as BBMap,\"");
			tsw.println("echo \"      or reduce Samtools' memory setting (the -m flag).\"");
		}
		
		for(String sam : set){
			String bam;
			if(sam.endsWith(".sam.gz")){bam=sam.substring(0, sam.length()-6)+"bam";}
			else if(sam.endsWith(".sam")){bam=sam.substring(0, sam.length()-3)+"bam";}
			else{bam=sam;} //Hopefully, they must have outputted a bam file using samtools.
			String bam2=bam.substring(0, bam.length()-4)+"_sorted";
			String bam3=bam2+".bam";
			
			if(Data.SAMTOOLS() && !Data.SAMTOOLS_VERSION_1x){
				//do nothing
			}else{
				bam2="-o "+bam2+".bam";
			}
			
			boolean pipe=true;
			if(pipe && sam!=bam){
//				if(Data.SAMTOOLS() && !Data.SAMTOOLS_VERSION_1x){
					tsw.println("echo \"Note: Please ignore any warnings about 'EOF marker is absent'; " +
							"this is a bug in samtools that occurs when using piped input.\"");
//				}
				tsw.println("samtools view -bShu "+sam+" | samtools sort -m "+memstring+" -@ 3 - "+bam2);
			}else{
				if(sam!=bam){tsw.println("samtools view -bSh1 -o "+bam+" "+sam);}
				tsw.println("samtools sort -m "+memstring+" -@ 3 "+bam+" "+bam2);
			}
			
			tsw.println("samtools index "+bam3);
		}
		tsw.poisonAndWait();
		
		try {
			File f=new File(outname);
			f.setExecutable(true, false);
		} catch (Exception e) {
//			e.printStackTrace();
		}
	}

	private static final ArrayList<Read> blank=new ArrayList<Read>(0);
	public static HashMap<String, Writer> streamTable=null;
	public static HashMap<String, Writer> streamTableAmbiguous=null;

	/** Clear only owned output references; AbstractMapper resets shared config/statistics. */
	public static void clearStatics(){
		streamTable=null;
		streamTableAmbiguous=null;
	}
}
