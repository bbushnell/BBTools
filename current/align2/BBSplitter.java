package align2;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SiteScore;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Mar 19, 2013
 *
 */
public class BBSplitter {
	
	public static void main(String[] args){
		if(Shared.COMMAND_LINE==null){
			Shared.COMMAND_LINE=(args==null ? null : args.clone());
			Shared.BBMAP_CLASS="BBSplitter";
		}
		Timer t=new Timer();
		String[] margs=processArgs(args);
		ReadWrite.waitForWritingToFinish();
		t.stop();
		Data.sysout.println("Ref merge time:     \t"+t);
		Data.scaffoldPrefixes=true;
		if(MAP_MODE==MAP_NORMAL){
			BBMap.main(margs);
		}else if(MAP_MODE==MAP_ACC){
			BBMapAcc.main(margs);
		}else if(MAP_MODE==MAP_PACBIO){
			BBMapPacBio.main(margs);
		}else if(MAP_MODE==MAP_PACBIOSKIMMER){
			BBMapPacBioSkimmer.main(margs);
		}else{
			throw new RuntimeException();
		}
//		Data.sysout.println("\nTotal time:     \t"+t);
		
		clearStatics();
	}
	
	public static String[] processArgs(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
//		if(ReadWrite.ZIPLEVEL<2){ReadWrite.ZIPLEVEL=2;} //Should be fine for a realistic number of threads, except in perfect mode with lots of sites and a small index.
		String[] oldargs=args;
		args=remakeArgs(args);
		if(args!=oldargs){
			Data.sysout.println("Converted arguments to "+Arrays.toString(args));
		}
		
		AbstractMapper.DEFAULT_OUTPUT_FORMAT=FileFormat.FASTQ;
		ReadWrite.ZIPLEVEL=2;
		
		Timer t=new Timer();
		
		
		int ziplevel=-1;
		int build=1;
		
		LinkedHashSet<String> nameSet=new LinkedHashSet<String>();
		HashMap<String, LinkedHashSet<String>> table=new HashMap<String, LinkedHashSet<String>>();
		
		ArrayList<String> unparsed=new ArrayList<String>();
		
		String basename=null;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a0=split[0];
			String a=a0.toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("blacklist") || a.equals("ref_blacklist")){a="ref_blacklist";}
			if(a.equals("whitelist") || a.equals("ref_whitelist")){a="ref_whitelist";}
			if(a.equals("ref") || a.equals("reference")){a="ref_ref";}

			if(b!=null && (a.startsWith("ref_"))){
				String setName=a0.substring(4);
				if(setName.indexOf(',')>=0){setName=setName.replace(',', '_');}
				if(setName.indexOf('$')>=0){setName=setName.replace('$', '_');}
				nameSet.add(setName);
				if(!table.containsKey(setName)){table.put(setName, new LinkedHashSet<String>());}
				LinkedHashSet<String> set=table.get(setName);
				
				File f;
				if((f=new File(b)).exists()){
					try {
						String s=f.getCanonicalPath();
						set.add(s);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}else{
					for(String x : b.split(",")){
						f=new File(x);
						if(f.exists()){
							try {
								set.add(f.getCanonicalPath());
							} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}else{
							assert(x.startsWith("stdin")) : "Can't find file "+x;
							set.add(x);
						}
					}
				}
//				assert(false) : a+", "+b+", "+arg+", "+setName;
			}else{
				if(a.startsWith("-xmx") || a.startsWith("-xms")){
					//jvm argument; do nothing
				}else if(a.equals("path") || a.equals("root")){
					Data.setPath(b);
				}else if(a.equals("ziplevel") || a.equals("zl")){
					ReadWrite.ZIPLEVEL=Integer.parseInt(b);
					unparsed.add(args[i]);
				}else if(a.equals("build")){
					build=Integer.parseInt(b);
					unparsed.add(args[i]);
				}else if(a.equals("basename") || a.equals("pattern")){
					basename=b;
					assert(b==null || (b.indexOf('%')>=0 && (b.indexOf('%')<b.lastIndexOf('.')))) :
						"basename must contain a '%' symbol prior to file extension.";
				}else if(a.equals("append") || a.equals("app")){
					append=ReadStats.append=Parse.parseBoolean(b);
//					sysout.println("Set append to "+append);
					unparsed.add(args[i]);
				}else if(a.equals("overwrite") || a.equals("ow")){
					overwrite=Parse.parseBoolean(b);
//					Data.sysout.println("Set overwrite to "+overwrite);
					unparsed.add(args[i]);
				}else if(a.equals("verbose")){
					verbose=Parse.parseBoolean(b);
					unparsed.add(args[i]);
				}else if(a.equals("rebuild")){
					forceRebuild=Parse.parseBoolean(b);
					unparsed.add(args[i]);
				}else if(a.equals("fastawrap") || a.equals("wrap")){
					Shared.FASTA_WRAP=Integer.parseInt(b);
				}else{
					unparsed.add(args[i]);
				}
			}
		}
		
		String refname=mergeReferences(nameSet, table, build);
		ArrayList<String> outnames=gatherLists(nameSet, basename);
//		unparsed.add("scaffoldprefixes=true");
		unparsed.add("ref="+refname);
		
		String[] margs=new String[unparsed.size()+(outnames==null ? 0 : outnames.size())];
		int idx=0;
		for(int i=0; i<unparsed.size(); i++){
			margs[idx]=unparsed.get(i);
			idx++;
		}
		if(outnames!=null){
			for(int i=0; i<outnames.size(); i++){
				margs[idx]=outnames.get(i);
				idx++;
			}
		}
		
		return margs;
	}
	
	
	public static String[] remakeArgs(String[] args){
		
		LinkedHashSet<String> set=new LinkedHashSet<String>();
		HashMap<String,LinkedHashSet<String>> map=new HashMap<String,LinkedHashSet<String>>();
		int removed=0;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("mapmode") && b!=null){
				args[i]=null;
				removed++;
				if(b.equalsIgnoreCase("normal")){MAP_MODE=MAP_NORMAL;}
				else if(b.equalsIgnoreCase("accurate") || b.equalsIgnoreCase("acc")){MAP_MODE=MAP_ACC;}
				else if(b.equalsIgnoreCase("pacbio") || b.equalsIgnoreCase("pb") || b.equalsIgnoreCase("bp")){MAP_MODE=MAP_PACBIO;}
				else if(b.equalsIgnoreCase("pacbioskimmer") || b.equalsIgnoreCase("pbs") || b.equalsIgnoreCase("bps")){MAP_MODE=MAP_PACBIOSKIMMER;}
				else{throw new RuntimeException("Unknown mode: "+b);}
			}else if(a.equals("ref") && b!=null){
				args[i]=null;
				removed++;
				processRef(b, set, map);
			}
		}
		if(set.isEmpty() && removed==0){return args;}
		if(MAP_MODE==MAP_ACC){removed--;}
		String[] args2=new String[args.length+set.size()-removed];
		
		int i=0, j=0;
		if(MAP_MODE==MAP_ACC){
			args2[j]="minratio=0.4"; //Increase sensitivity in accurate mode
			j++;
		}
		for(; i<args.length; i++){
			if(args[i]!=null){
				args2[j]=args[i];
				j++;
			}
		}
		for(String key : set){
			LinkedHashSet<String> list=map.get(key);
			StringBuilder sb=new StringBuilder(200);
			sb.append("ref_"+key+"=");
			String comma="";
			for(String s : list){
				sb.append(comma);
				sb.append(s);
				comma=",";
			}
			args2[j]=sb.toString();
			j++;
		}
		return args2;
	}
	
	private static void processRef(String b, LinkedHashSet<String> set, HashMap<String,LinkedHashSet<String>> map){
		
		ArrayList<String> files=(ArrayList<String>)Tools.getFileOrFiles(b, null, true, false, false, false);
		for(String file : files){
			String name=file.replace('\\', '/');
			int x=name.lastIndexOf('/');
			if(x>=0){name=name.substring(x+1);}
			while(name.endsWith(".zip") || name.endsWith(".bz2") || name.endsWith(".gz") || name.endsWith(".gzip")){
				name=name.substring(0, name.lastIndexOf('.'));
			}
			if(name.lastIndexOf('.')>=0){
				name=name.substring(0, name.lastIndexOf('.'));
			}
			set.add(name);
			LinkedHashSet<String> list=map.get(name);
			if(list==null){
				list=new LinkedHashSet<String>();
				map.put(name, list);
			}
			list.add(file);
		}
	}
	
	public static ArrayList<String> gatherLists(LinkedHashSet<String> nameSet, String basename){
		if(basename==null){return null;}
		ArrayList<String> args=new ArrayList<String>();
		for(String name : nameSet){
			if(basename!=null){
				args.add("out_"+name+"="+(basename.replaceFirst("%", name)));
			}
		}
		return args;
	}
	
	
	public static String mergeReferences(LinkedHashSet<String> nameSet, HashMap<String, LinkedHashSet<String>> nameToFileTable, int build){
		LinkedHashSet<String> fnames=new LinkedHashSet<String>();
//		nameSet.remove("blacklist");
//		nameSet.remove("whitelist");
		addNames(fnames, nameToFileTable, "whitelist");
		for(String s : nameSet){
			if(!s.equalsIgnoreCase("blacklist") && !s.equalsIgnoreCase("whitelist")){
				addNames(fnames, nameToFileTable, s);
			}
		}
		addNames(fnames, nameToFileTable, "blacklist");
		
		final HashMap<String, LinkedHashSet<String>> fileToNameTable=new HashMap<String, LinkedHashSet<String>>();
		for(String name : nameSet){
			LinkedHashSet<String> files=nameToFileTable.get(name);
			if(files!=null){
				for(String f : files){
					LinkedHashSet<String> names=fileToNameTable.get(f);
					if(names==null){
						names=new LinkedHashSet<String>();
						fileToNameTable.put(f, names);
					}
					names.add(name);
				}
			}
		}
		
		final String root=Data.ROOT_GENOME+build;
		{
			File f=new File(root);
			if(!f.exists()){f.mkdirs();}
		}
		
		{
			final String reflist=root+"/reflist.txt";
			final String namelist=root+"/namelist.txt";
			final boolean reflistExists=new File(reflist).exists();
			boolean writeReflist=false;
			String[] oldrefs=null;
			String[] oldnames=null;
			if(reflistExists){
				TextFile tf=new TextFile(reflist, false);
				oldrefs=tf.toStringLines();
				tf.close();
				
				tf=new TextFile(namelist, false);
				oldnames=tf.toStringLines();
				tf.close();
			}
			if(fnames.size()>0){
				writeReflist=true;
				ArrayList<String> fl=new ArrayList<String>(fnames.size());
				fl.addAll(fnames);
				ArrayList<String> nl=new ArrayList<String>(nameSet.size());
				nl.addAll(nameSet);
				//TODO: Compare old to new
			}else{
				assert(oldrefs!=null) : "No reference specified, and none exists.  Please regenerate the index.";
				for(String s : oldrefs){fnames.add(s);}

				assert(oldnames!=null) : "No reference specified, and none exists.  Please regenerate the index.";
				for(String s : oldnames){nameSet.add(s);}
				
				writeReflist=false;
			}
			if(writeReflist){
				{
//					assert(false) : fnames;
//					assert(fnames.size()>0);
					TextStreamWriter tsw=new TextStreamWriter(reflist, overwrite, append, false);
					tsw.start();
					for(String s : fnames){tsw.println(s);}
					tsw.poisonAndWait();
					assert(new File(reflist).exists()) : reflist+".exists? "+new File(reflist).exists();
				}
				{
//					assert(nameSet.size()>0);
					TextStreamWriter tsw=new TextStreamWriter(namelist, overwrite, append, false);
					tsw.start();
					for(String s : nameSet){tsw.println(s);}
					tsw.poisonAndWait();
				}
			}
		}
		
		if(fnames.size()<1){
			assert(false) : "No references specified." +
					"\nTODO:  This is really annoying; I need to include reference names in some auxillary file.";
			return null;
		}else if(fnames.size()==1){
//			Data.sysout.println("Only one reference file; skipping merge.");
//			String refname=fnames.iterator().next();
//			return refname;
		}
		
		long key=0;
		for(String s : nameSet){
			key=Long.rotateLeft(key, 21);
			key=key^s.hashCode();
//			System.err.println("Hashed nameSet "+nameSet+" -> "+key);
		}
		key=(key&Long.MAX_VALUE);
		String refname0="merged_ref_"+key+".fa.gz";
		String refname=root+"/"+refname0;
		
		if(!forceRebuild){
			File f=new File(refname);
			if(f.exists()){
				//			Data.sysout.println("Merged reference file /ref/genome/"+build+"/"+refname0+" already exists; skipping merge.");
				Data.sysout.println("Merged reference file "+refname+" already exists; skipping merge.");
				return refname;
			}
//			else{
//				f=new File(root);
//				if(!f.exists()){f.mkdirs();}
//			}
		}
		//			Data.sysout.println("Creating merged reference file /ref/genome/"+build+"/"+refname0);
		Data.sysout.println("Creating merged reference file "+refname);
		
		TextStreamWriter tsw=new TextStreamWriter(refname, overwrite || forceRebuild, false, true);
		tsw.start();
		for(String fname : fnames){
			TextFile tf=new TextFile(fname, false);
			LinkedHashSet<String> listnames=fileToNameTable.get(fname);
//			assert(false) : "\n\n"+fname+"\n\n"+listnames+"\n\n"+fileToNameTable+"\n\n"+nameSet+"\n\n"+nameToFileTable+"\n\n";
			String prefix=null;
			{
				StringBuilder sb=new StringBuilder(100);
				sb.append('>');
				if(listnames!=null){
					String sep="";
					for(String s : listnames){
						sb.append(sep);
						sb.append(s);
						sep=",";
					}
				}
				sb.append('$');
				prefix=sb.toString();
			}
//			assert(false) : prefix;
//			System.err.println(prefix);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				if(prefix!=null && line.charAt(0)=='>'){
					
					tsw.print(prefix);
					tsw.println(line.substring(1));
				}else{
					tsw.println(line);
				}
			}
			tf.close();
		}
		tsw.poisonAndWait();
		
		return refname;
	}
	
	/** Returns the set of scaffold name prefixes or suffixes.
	 * 
	 * @param getPrefixes True to return prefixes (set names), false to return suffixes (scaffold names)
	 * @return A set of affixes
	 */
	public static HashSet<String> getScaffoldAffixes(boolean getPrefixes){
		final byte[][][] b3=Data.scaffoldNames;
		
		int size=(int)Tools.min((10+Data.numScaffolds*3)/2, Integer.MAX_VALUE);
		HashSet<String> set=new HashSet<String>(size);
		
		assert(b3!=null);
		for(byte[][] b2 : b3){
			if(b2!=null){
				for(byte[] bname : b2){
					if(bname!=null){
						int idx=Tools.indexOf(bname, (byte)'$');
						String prefix=null, suffix=null;
						if(idx>=0){
							if(getPrefixes){prefix=new String(bname, 0, idx);}
							else{suffix=new String(bname, idx+1, bname.length-idx-1);}
						}else{
							if(!getPrefixes){suffix=new String(bname);}
						}

						if(getPrefixes){
							if(prefix!=null){
								if(prefix.indexOf(',')>=0){
									for(String s : prefix.split(",")){
										set.add(s);
									}
								}else{
									set.add(prefix);
								}
							}
						}else{
							if(suffix!=null){
								set.add(suffix);
							}
						}
					}
				}
			}
		}
		return set;
	}
	
	public static synchronized HashMap<String, ConcurrentReadOutputStream> makeOutputStreams(String[] args, boolean OUTPUT_READS, boolean OUTPUT_ORDERED_READS,
			int buff, boolean paired, boolean overwrite_, boolean append_, boolean ambiguous){
//		assert(false) : Arrays.toString(args);
		HashMap<String, ConcurrentReadOutputStream> table=new HashMap<String, ConcurrentReadOutputStream>();
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
				ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, null, null, buff, null, false);
				ros.start();
//				Data.sysout.println("Started output stream:\t"+t);
				table.put(name, ros);
				AbstractMapThread.OUTPUT_SAM|=ff1.samOrBam();
			}
		}
		return table.isEmpty() ? null : table;
	}
	
	
	public static synchronized LinkedHashMap<String, SetCount> makeSetCountTable(){
		assert(setCountTable==null);
		HashSet<String> names=getScaffoldAffixes(true);
		setCountTable=new LinkedHashMap<String, SetCount>();
		for(String s : names){setCountTable.put(s, new SetCount(s));}
		return setCountTable;
	}
	
	
	public static synchronized LinkedHashMap<String, SetCount> makeScafCountTable(){
		assert(scafCountTable==null);
		HashSet<String> names=getScaffoldAffixes(false);
		scafCountTable=new LinkedHashMap<String, SetCount>();
		for(String s : names){scafCountTable.put(s, new SetCount(s));}
//		System.out.println("Made table "+scafCountTable);
		return scafCountTable;
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
			ConcurrentReadOutputStream tros=streamTable.get(s);
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
				ConcurrentReadOutputStream tros=streamTable.get(s);
				tros.add(alr, listID);
			}
		}
		if(streamTableAmbiguous!=null){
			for(String s : streamTableAmbiguous.keySet()){
				ArrayList<Read> alr=splitTableA.get(s);
				if(alr==null){alr=blank;}
				ConcurrentReadOutputStream tros=streamTableAmbiguous.get(s);
				tros.add(alr, listID);
			}
		}
		if(clear){splitTable.clear();}
		if(clearA){splitTableA.clear();}
	}
	
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
	
	//*********************************
	
	public static HashSet<String>[] getSets(Read r1, int clearzone, HashSet<String>[] sets){
		Read r2=r1.mate;
		if(!r1.mapped() && (r2==null || !r2.mapped())){return null;}
		
		if(sets==null){
			assert(false);
			sets=new HashSet[4];
		}else{
			for(HashSet<String> set : sets){
				assert(set==null || set.isEmpty());
			}
		}
		
		HashSet<String> primary1=sets[0], other1=sets[1], primary2=sets[2], other2=sets[3];
		if(r1.mapped()){
//			System.out.println(r1.list.size());
			SiteScore s0=r1.topSite();
			primary1=toListNames(s0, primary1);
			for(int i=1; i<r1.sites.size(); i++){
				SiteScore ss=r1.sites.get(i);
				if(ss.score+clearzone<s0.score){break;}
				other1=toListNames(ss, other1);
			}
//			System.out.println(primary1);
//			System.out.println(other1);
		}
		if(r2!=null && r2.mapped()){
			SiteScore s0=r2.topSite();
			primary2=toListNames(s0, primary2);
			for(int i=1; i<r2.sites.size(); i++){
				SiteScore ss=r2.sites.get(i);
				if(ss.score+clearzone<s0.score){break;}
				other2=toListNames(ss, other2);
			}
		}
		sets[0]=primary1;
		sets[1]=other1;
		sets[2]=primary2;
		sets[3]=other2;
		return sets;
	}
	
	
	public static HashSet<String> getScaffolds(Read r1, int clearzone, HashSet<String> set, boolean includeMate){
		Read r2=(includeMate ? r1.mate : null);
		if(!r1.mapped() && (r2==null || !r2.mapped())){return null;}
		assert(set==null || set.isEmpty());
		
		if(!r1.ambiguous() && (r2==null || !r2.ambiguous())){
			byte[] scafb1=r1.getScaffoldName(false);
			byte[] scafb2=(r2==null ? null : r2.getScaffoldName(false));
			if(scafb1==null){scafb1=scafb2;}
			if(scafb1==null){
				assert(false) : r1;
				return null;
			}
			if(scafb2==null || scafb1==scafb2){
				final String s;
				if(Data.scaffoldPrefixes){
					int idx=Tools.indexOf(scafb1, (byte)'$');
					assert(idx>=0) : idx+", "+new String(scafb1);
					s=(idx>=0 ? new String(scafb1, idx+1, scafb1.length-idx-1) : new String(scafb1));
				}else{
					s=new String(scafb1);
				}
				if(set==null){set=new HashSet<String>(1);}
//				assert(!s.contains("$")) : s+", "+Data.scaffoldPrefixes+", "+Tools.indexOf(scafb1, (byte)'$');
				set.add(s);
				return set;
			}
		}
		
		if(set==null){set=new HashSet<String>(4);}
		if(r1.mapped()){
			SiteScore s0=r1.topSite();
			for(SiteScore ss : r1.sites){
				if(ss.score+clearzone<s0.score){break;}
				byte[] b=ss.getScaffoldName(false);
				assert(b!=null);
				final String s;
				if(Data.scaffoldPrefixes){
					int idx=Tools.indexOf(b, (byte)'$');
					assert(idx>=0) : idx+", "+new String(b);
					s=(idx>=0 ? new String(b, idx+1, b.length-idx-1) : new String(b));
				}else{
					s=new String(b);
				}
				set.add(s);
			}
		}
		if(r2!=null && r2.mapped()){
			SiteScore s0=r2.topSite();
			for(SiteScore ss : r2.sites){
				if(ss.score+clearzone<s0.score){break;}
				byte[] b=ss.getScaffoldName(false);
				assert(b!=null);
				final String s;
				if(Data.scaffoldPrefixes){
					int idx=Tools.indexOf(b, (byte)'$');
					assert(idx>=0) : idx+", "+new String(b);
					s=(idx>=0 ? new String(b, idx+1, b.length-idx-1) : new String(b));
				}else{
					s=new String(b);
				}
				set.add(s);
			}
		}
		assert(set.size()>0);
		return set;
	}
	
	
	/**
	 * @param r
	 * @return A set of names of reference lists containing this read or its mate.
	 */
	public static HashSet<String> toListNames(Read r, HashSet<String> set) {
		if(r==null){return set;}
		byte[] scaf1=r.getScaffoldName(false);
		byte[] scaf2=(r.mate==null ? null : r.mate.getScaffoldName(false));
		if(scaf1==null && scaf2==null){return set;}
		
		if(set==null){set=new HashSet<String>(8);}
		else{assert(set.isEmpty());}
		
		int x=scaf1==null ? -1 : Tools.indexOf(scaf1, (byte)'$');
		if(x>=0){
			String s=new String(scaf1, 0, x);
			if(s.indexOf(',')<0){
				set.add(s);
			}else{
				for(String s2 : s.split(",")){set.add(s2);}
			}
		}
		
		x=(scaf2==null || scaf2==scaf1) ? -1 : Tools.indexOf(scaf2, (byte)'$');
		if(x>=0){
			String s=new String(scaf2, 0, x);
			if(s.indexOf(',')<0){
				set.add(s);
			}else{
				for(String s2 : s.split(",")){set.add(s2);}
			}
		}
		
		return set;
	}
	
	
	/**
	 * @param r
	 * @return A mapping of reference names to read clones.
	 */
	public static HashMap<String, Read> toNameMap(Read r, HashMap<String, Read> map) {

		if(true){throw new RuntimeException("TODO");}
		
		if(r==null){return map;}
		byte[] scaf1=r.getScaffoldName(false);
		byte[] scaf2=(r.mate==null ? null : r.mate.getScaffoldName(false));
		if(scaf1==null && scaf2==null){return map;}
		
		if(map==null){map=new HashMap<String, Read>(8);}
		else{assert(map.isEmpty());}
		
		int x=scaf1==null ? -1 : Tools.indexOf(scaf1, (byte)'$');
		if(x>=0){
			String s=new String(scaf1, 0, x);
			if(s.indexOf(',')<0){
				if(!map.containsKey(s)){
					map.put(s, cloneRead(r, s));
				}
			}else{
				for(String s2 : s.split(",")){
					if(!map.containsKey(s2)){
						map.put(s2, cloneRead(r, s2));
					}
				}
			}
		}
		
		x=(scaf2==null || scaf2==scaf1) ? -1 : Tools.indexOf(scaf2, (byte)'$');
		if(x>=0){
			String s=new String(scaf2, 0, x);
			if(s.indexOf(',')<0){
				if(!map.containsKey(s)){
					map.put(s, cloneRead(r, s));
				}
			}else{
				for(String s2 : s.split(",")){
					if(!map.containsKey(s2)){
						map.put(s2, cloneRead(r, s2));
					}
				}
			}
		}
		
		return map;
	}
	
	private static Read cloneRead(Read r, String ref){
		throw new RuntimeException("TODO");
	}

	
	/**
	 * @param r
	 * @return A set of names of reference lists containing this site.
	 */
	public static HashSet<String> toListNames(SiteScore r, HashSet<String> set) {
		if(r==null){return set;}
		byte[] scaf1=r.getScaffoldName(false);
		if(scaf1==null){return set;}
		if(set==null){set=new HashSet<String>(8);}
		int x=scaf1==null ? -1 : Tools.indexOf(scaf1, (byte)'$');
		if(x>=0){
			String s=new String(scaf1, 0, x);
			if(s.indexOf(',')<0){
				set.add(s);
			}else{
				for(String s2 : s.split(",")){set.add(s2);}
			}
		}
		return set;
	}

	private static void addNames(LinkedHashSet<String> fnames, HashMap<String, LinkedHashSet<String>> table, String setName){
		LinkedHashSet<String> set=table.get(setName);
		if(set==null){return;}
		for(String s : set){fnames.add(s);}
	}
	
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
			for(ConcurrentReadOutputStream ros : streamTable.values()){
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
	
	public static class SetCount implements Comparable<SetCount>{
		
		public SetCount(String s){
			name=s;
		}

		@Override
		public boolean equals(Object other){return equals((SetCount)other);}
		public boolean equals(SetCount other){return compareTo(other)==0;}
		
		@Override
		public int hashCode() {
			assert(false) : "This class should not be hashed.";
			return super.hashCode();
		}
		
		@Override
		public int compareTo(SetCount o) {
			if(mappedReads!=o.mappedReads){return mappedReads>o.mappedReads ? 1 : -1;}
			if(ambiguousReads!=o.ambiguousReads){return ambiguousReads>o.ambiguousReads ? 1 : -1;}
			return name.compareTo(o.name);
		}
		
		@Override
		public String toString(){
			return name+", "+mappedReads+", "+ambiguousReads+", "+mappedBases+", "+ambiguousBases;
		}
		
		public final String name;
		public long mappedReads;
		public long ambiguousReads;
		public long assignedReads;
		public long mappedBases;
		public long ambiguousBases;	
		public long assignedBases;
		
	}
	
	public static void printCounts(String fname, LinkedHashMap<String, SetCount> map, boolean header, long totalReads, boolean nzo, boolean sort){
		final ArrayList<SetCount> list=new ArrayList<SetCount>(map.size());
		for(String name : map.keySet()){
			list.add(map.get(name));
		}
		final TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, append, false);
		tsw.start();
		if(sort){
			Shared.sort(list);
			Collections.reverse(list);
		}
		
		if(header){
			tsw.print("#name\t%unambiguousReads\tunambiguousMB\t%ambiguousReads\tambiguousMB\tunambiguousReads\tambiguousReads\tassignedReads\tassignedBases\n");
		}
		final StringBuilder sb=new StringBuilder(1024);
		final double divR=100.0/(totalReads);
		final double divB=1.0/1000000;
		for(SetCount sc : list){
			if(!nzo || sc.mappedReads>0 || sc.ambiguousReads>0 || sc.assignedReads>0){
				sb.append(sc.name).append('\t');
				sb.append(Tools.format("%.5f\t", sc.mappedReads*divR));
				sb.append(Tools.format("%.6f\t", sc.mappedBases*divB));
				sb.append(Tools.format("%.5f\t", sc.ambiguousReads*divR));
				sb.append(Tools.format("%.6f\t", sc.ambiguousBases*divB));
				sb.append(sc.mappedReads).append('\t');
				sb.append(sc.ambiguousReads).append('\t');
				sb.append(sc.assignedReads).append('\t');
				sb.append(sc.assignedBases);
				sb.append('\n');
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		tsw.poisonAndWait();
	}
	
	static final void clearStatics(){


		setCountTable=null;
		 scafCountTable=null;
		streamTable=null;
		streamTableAmbiguous=null;
		
		AMBIGUOUS2_MODE=AMBIGUOUS2_UNSET;
		TRACK_SET_STATS=false;
		TRACK_SCAF_STATS=false;
		SCAF_STATS_FILE=null;
		SET_STATS_FILE=null;
		overwrite=true;
		append=false;
		verbose=false;
		forceRebuild=false;
		
		MAP_MODE=MAP_NORMAL;
	}

	public static LinkedHashMap<String, SetCount> setCountTable=null;
	public static LinkedHashMap<String, SetCount> scafCountTable=null;
	
	/**
	 * Holds named output streams.
	 */
	public static HashMap<String, ConcurrentReadOutputStream> streamTable=null;
	
	/**
	 * Holds named output streams for ambiguous (across different references) reads.
	 */
	public static HashMap<String, ConcurrentReadOutputStream> streamTableAmbiguous=null;
	public static final int AMBIGUOUS2_UNSET=0;
	public static final int AMBIGUOUS2_FIRST=1;
	public static final int AMBIGUOUS2_SPLIT=2;
	public static final int AMBIGUOUS2_TOSS=3;
	public static final int AMBIGUOUS2_RANDOM=4;
	public static final int AMBIGUOUS2_ALL=5;
	public static int AMBIGUOUS2_MODE=AMBIGUOUS2_UNSET;
	public static boolean TRACK_SET_STATS=false;
	public static boolean TRACK_SCAF_STATS=false;
	public static String SCAF_STATS_FILE=null;
	public static String SET_STATS_FILE=null;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean verbose=false;
	public static boolean forceRebuild=false;
	private static final ArrayList<Read> blank=new ArrayList<Read>(0);
	
	public static final int MAP_NORMAL=1;
	public static final int MAP_ACC=2;
	public static final int MAP_PACBIO=3;
	public static final int MAP_PACBIOSKIMMER=4;
	public static int MAP_MODE=MAP_NORMAL;
	
}
