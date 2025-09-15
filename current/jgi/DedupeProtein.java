package jgi;

import java.io.File;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.concurrent.atomic.AtomicIntegerArray;

import align2.BandedAligner;
import dna.AminoAcid;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import sort.ReadComparator;
import sort.ReadComparatorID;
import sort.ReadComparatorName;
import sort.ReadLengthComparator;
import sort.ReadQualityComparator;
import stream.ConcurrentCollectionReadInputStream;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.FastqReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongM;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Nov 30, 2017
 *
 */
public final class DedupeProtein {
	
	public static void main(String[] args){
		int rbl0=Shared.bufferLen();;
		Shared.setBufferLen(16);
		
		DedupeProtein dd=new DedupeProtein(args);
		dd.process();
		
		Shared.setBufferLen(rbl0);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public DedupeProtein(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), true);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_PIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Read.TO_UPPER_CASE=true;
		
		boolean setMcsfs=false;
		int bandwidth_=-1;
		int k_=31;
		int subset_=0, subsetCount_=1;
		boolean ascending=false;
		Parser parser=new Parser();
		
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseQTrim(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.indexOf(',')>=0 && !new File(b).exists()){
					in1=b.split(",");
				}else{
					in1=new String[] {b};
				}
			}else if(a.equals("in2")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.indexOf(',')>=0 && !new File(b).exists()){
					in2=b.split(",");
				}else{
					in2=new String[] {b};
				}
			}else if(a.equals("out") || a.equals("out")){
				out=b;
			}else if(a.equals("out2")){
				throw new RuntimeException("Dedupe does not allow 'out2'; for paired reads, output is interleaved.");
			}else if(a.equals("clusterfilepattern") || a.equals("pattern")){
				clusterFilePattern=b;
				assert(clusterFilePattern==null || clusterFilePattern.contains("%")) : "pattern must contain the % symbol.";
			}else if(a.equals("outbest")){
				outbest=b;
			}else if(a.equals("outd") || a.equals("outduplicate")){
				outdupe=b;
			}else if(a.equals("csf") || a.equals("clusterstatsfile")){
				outcsf=b;
			}else if(a.equals("dot") || a.equals("graph") || a.equals("outdot") || a.equals("outgraph")){
				outgraph=b;
			}else if(a.equals("mcsfs") || a.equals("minclustersizeforstats")){
				minClusterSizeForStats=Integer.parseInt(b);
			}else if(a.equals("mcs") || a.equals("minclustersize")){
				minClusterSize=Integer.parseInt(b);
				if(!setMcsfs){
					minClusterSizeForStats=minClusterSize;
				}
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("sort")){
				if(b==null){sort=true;}
				else if(b.equalsIgnoreCase("a") || b.equalsIgnoreCase("ascending")){
					sort=true;
					ascending=true;
				}else if(b.equalsIgnoreCase("d") || b.equalsIgnoreCase("descending")){
					sort=true;
					ascending=false;
				}else if(b.equalsIgnoreCase("length")){
					sort=true;
					comparator=ReadLengthComparator.comparator;
				}else if(b.equalsIgnoreCase("quality")){
					sort=true;
					comparator=ReadQualityComparator.comparator;
				}else if(b.equalsIgnoreCase("name")){
					sort=true;
					comparator=ReadComparatorName.comparator;
				}else if(b.equalsIgnoreCase("id")){
					sort=true;
					comparator=ReadComparatorID.comparator;
				}else{
					sort=Parse.parseBoolean(b);
				}
			}else if(a.equals("ascending")){
				ascending=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				boolean x=Parse.parseBoolean(b);
				if(x){
					ascending=true;
					sort=true;
					comparator=ReadComparatorID.comparator;
				}else{
					//do nothing
				}
			}else if(a.equals("ac") || a.equals("absorbcontainment") || a.equals("absorbcontainments") || a.equals("tc") || a.equals("testcontainment") || a.equals("containment")){
				absorbContainment=Parse.parseBoolean(b);
			}else if(a.equals("am") || a.equals("absorbmatch") || a.equals("absorbmatches") || a.equals("tm") || a.equals("testmatch")){
				absorbMatch=Parse.parseBoolean(b);
			}else if(a.equals("ao") || a.equals("absorboverlap") || a.equals("absorboverlaps") || a.equals("to") || a.equals("testoverlap")){
				absorbOverlap=Parse.parseBoolean(b);
			}else if(a.equals("fo") || a.equals("findoverlap") || a.equals("findoverlaps")){
				findOverlaps=Parse.parseBoolean(b);
			}else if(a.equals("c") || a.equals("cluster") || a.equals("clusters")){
				makeClusters=Parse.parseBoolean(b);
			}else if(a.equals("fmj") || a.equals("fixmultijoin") || a.equals("fixmultijoins")){
				fixMultiJoins=Parse.parseBoolean(b);
			}else if(a.equals("foc") || a.equals("fixoffsetcontradiction") || a.equals("fixoffsetcontradictions")){
				fixOffsetContradictions=Parse.parseBoolean(b);
			}else if(a.equals("pto") || a.equals("preventtransitiveoverlap") || a.equals("preventtransitiveoverlaps")){
				preventTransitiveOverlaps=Parse.parseBoolean(b);
			}else if(a.equals("pbr") || a.equals("pickbestrepresentative")){
				pickBestRepresentativePerCluster=Parse.parseBoolean(b);
			}else if(a.equals("mst") || a.equals("maxspanningtree")){
				maxSpanningTree=Parse.parseBoolean(b);
			}else if(a.equals("pc") || a.equals("processcluster") || a.equals("processclusters")){
				processClusters=Parse.parseBoolean(b);
			}else if(a.equals("rnc") || a.equals("renamecluster") || a.equals("renameclusters")){
				renameClusters=Parse.parseBoolean(b);
				if(renameClusters){storeName=false;}
			}else if(a.equals("rc") || a.equals("removecycles") || a.equals("removecycle")){
				removeCycles=Parse.parseBoolean(b);
			}else if(a.equals("uo") || a.equals("uniqueonly")){
				UNIQUE_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("rmn") || a.equals("requirematchingnames")){
				REQUIRE_MATCHING_NAMES=Parse.parseBoolean(b);
			}else if(a.equals("ngn") || a.equals("numbergraphnodes")){
				NUMBER_GRAPH_NODES=Parse.parseBoolean(b);
			}else if(a.equals("addpairnum")){
				ADD_PAIRNUM_TO_NAME=Parse.parseBoolean(b);
			}else if(a.equals("hashns")){
				HASH_NS=Parse.parseBoolean(b);
			}else if(a.equals("pn") || a.equals("prefixname")){
//				PREFIX_NAME=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
				assert(k_>0 && k_<11) : "k must be between 1 and 10; default is 10, and lower values are slower.";
			}else if(a.equals("minscaf") || a.equals("ms")){
				MINSCAF=FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("mlp") || a.equals("minlengthpercent")){
				minLengthPercent=Float.parseFloat(b);
			}else if(a.equals("mop") || a.equals("minoverlappercent")){
				minOverlapPercentCluster=minOverlapPercentMerge=Float.parseFloat(b);
			}else if(a.equals("mopc") || a.equals("minoverlappercentcluster")){
				minOverlapPercentCluster=Float.parseFloat(b);
			}else if(a.equals("mopm") || a.equals("minoverlappercentmerge")){
				minOverlapPercentMerge=Float.parseFloat(b);
			}else if(a.equals("mo") || a.equals("minoverlap")){
				minOverlapCluster=minOverlapMerge=Integer.parseInt(b);
			}else if(a.equals("moc") || a.equals("minoverlapcluster")){
				minOverlapCluster=Integer.parseInt(b);
			}else if(a.equals("mom") || a.equals("minoverlapmerge")){
				minOverlapMerge=Integer.parseInt(b);
			}else if(a.equals("rt") || a.equals("rigoroustransitive")){
				rigorousTransitive=Parse.parseBoolean(b);
			}else if(a.equals("e") || a.equals("maxedits") || a.equals("edits") || a.equals("edist")){
				maxEdits=Integer.parseInt(b);
			}else if(a.equals("s") || a.equals("maxsubs") || a.equals("maxsubstitutions") || a.equals("hdist")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("bw") || a.equals("bandwidth")){
				bandwidth_=Integer.parseInt(b);
			}else if(a.equals("mid") || a.equals("minidentity")){
				minIdentity=Float.parseFloat(b);
				minIdentityMult=(minIdentity==100f ? 0 : (100f-minIdentity)/100f);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
//				BandedAligner.verbose=verbose;
			}else if(a.equals("contigbreak") || (arg.contains("=") && (a.equals("n") || a.equals("-n")))){
				maxNs=Integer.parseInt(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("sn") || a.equals("storename") || a.equals("storenames") || a.equals("keepnames")){
				storeName=Parse.parseBoolean(b);
			}else if(a.equals("ssx") || a.equals("storesuffix") || a.equals("storesuffixes")){
				storeSuffix=Parse.parseBoolean(b);
			}else if(a.equals("numaffixmaps") || a.equals("nam")){
				numAffixMaps=Integer.parseInt(b);
			}else if(a.equals("amino") || a.equals("protein")){
				//do nothing
			}else if(a.equals("mac") || a.equals("maxaffixcopies")){
				maxAffixCopies=Integer.parseInt(b);
			}else if(a.equals("me") || a.equals("maxedges")){
				maxEdges=Integer.parseInt(b);
				maxEdges2=maxEdges*2;
				if(maxEdges2<1){maxEdges2=Integer.MAX_VALUE-1;}
			}else if(a.equals("ignoreaffix1") || a.equals("ia1")){
				ignoreAffix1=Parse.parseBoolean(b);
			}else if(a.equals("parsedepth") || a.equals("pd")){
				parseDepth=Parse.parseBoolean(b);
			}else if(a.equals("printlengthinedges") || a.equals("ple")){
				printLengthInEdges=Parse.parseBoolean(b);
			}else if(a.equals("depthmult") || a.equals("depthratio") || a.equals("dr")){
				depthRatio=Float.parseFloat(b);
				if(depthRatio<=0){
					parseDepth=false;
				}else{
					parseDepth=true;
					assert(depthRatio>0);
					if(depthRatio<1){depthRatio=1/depthRatio;}
				}
			}else if(a.equals("storequality") || a.equals("sq")){
				storeQuality=Parse.parseBoolean(b);
			}else if(a.equals("exact") || a.equals("ex")){
				exact=Parse.parseBoolean(b);
			}else if(a.equals("uniquenames") || a.equals("un")){
				uniqueNames=Parse.parseBoolean(b);
			}else if(a.equals("ftl") || a.equals("forcetrimleft")){
				forceTrimLeft=Integer.parseInt(b);
			}else if(a.equals("ftr") || a.equals("forcetrimright")){
				forceTrimRight=Integer.parseInt(b);
			}else if(a.equals("subset") || a.equals("sst")){
				subset_=Integer.parseInt(b);
			}else if(a.equals("subsets") || a.equals("subsetcount") || a.equals("sstc")){
				subsetCount_=Integer.parseInt(b);
			}else if(i==0 && in1==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				String c=args[i];
				if(c.indexOf(',')>=0 && !new File(c).exists()){
					in1=c.split(",");
				}else{
					in1=new String[] {c};
				}
			}else if(i==1 && out==null && arg.indexOf('=')<0 && arg.lastIndexOf('.')>0){
				out=args[i];
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		comparator.setAscending(ascending);
		
		{//Process parser fields
			Parser.processQuality();
			qTrimLeft=parser.qtrimLeft;
			qTrimRight=parser.qtrimRight;
			trimQ=parser.trimq;
			trimE=parser.trimE();
		}
		
		if(verbose){
			ReadWrite.verbose=ConcurrentGenericReadInputStream.verbose=ConcurrentReadOutputStream.verbose=ByteFile1.verbose=ByteFile2.verbose=FastqReadInputStream.verbose=true;
		}
		
		k=k_;
		k2=k-1;
		subset=subset_;
		subsetCount=subsetCount_;
		subsetMode=subsetCount>1;
		assert(subset>=0 && subset<subsetCount) : "subset="+subset+", subsetCount="+subsetCount;
		
		BandedAligner.penalizeOffCenter=true;

		if(maxSpanningTree){removeCycles=fixMultiJoins=false;}
		if(absorbOverlap){processClusters=true;}
		if(processClusters || renameClusters || maxSpanningTree){makeClusters=true;}
		if(makeClusters){findOverlaps=true;}
		if(renameClusters){uniqueNames=/*storeName=*/false;}
		
		if(bandwidth_>-1){
			bandwidth=Tools.min(bandwidth_, 2*maxEdits+1);
			customBandwidth=(bandwidth<2*maxEdits+1);
		}else{
			bandwidth=2*maxEdits+1;
			customBandwidth=false;
		}
		maxSubs=Tools.max(maxSubs, maxEdits);
		if(maxSubs>0 || minIdentity<100 || findOverlaps){storeSuffix=true;}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		for(int i=0; i<in1.length; i++){
			if(in1[i].equalsIgnoreCase("stdin") && !new File(in1[i]).exists()){in1[i]="stdin.fa";}
		}
		
//		assert(false) : Arrays.toString(in);
		
//		if(!setOut && clusterFilePattern==null){out="stdout.fa";}
//		else
//			if("stdout".equalsIgnoreCase(out) || "standarddout".equalsIgnoreCase(out)){
//			out="stdout.fa";
//			outstream=System.err;
//		}
		if(!Tools.canWrite(out, overwrite)){throw new RuntimeException("Output file "+out+" already exists, and overwrite="+overwrite);}
		
		for(int i=0; i<in1.length; i++){
			assert(!in1[i].equalsIgnoreCase(out));
		}
//		assert(false) : "\nabsorbContainment="+absorbContainment+", findOverlaps="+findOverlaps+", absorbOverlap="+absorbOverlap+"\n"+
//			"processClusters="+processClusters+", renameClusters="+renameClusters+", makeClusters="+makeClusters+", uniqueNames="+uniqueNames+"\n"+
//			"storeName="+storeName+", DISPLAY_PROGRESS="+DISPLAY_PROGRESS+", removeCycles="+removeCycles;
		if(absorbContainment || findOverlaps){
//			assert(false);
			affixMaps=new HashMap[numAffixMaps];
			for(int i=0; i<numAffixMaps; i++){
				affixMaps[i]=new HashMap<LongM, ArrayList<Unit>>(4000000);
			}
			if(affixMaps.length>0){affixMap1=affixMaps[0];}
			if(affixMaps.length>1){affixMap2=affixMaps[1];}
		}
//		assert(false) : absorbContainment+", "+(affixMap==null);
		
		if(outdupe==null){
			dupeWriter=null;
		}else{
			FileFormat ff=FileFormat.testOutput(outdupe, FileFormat.FASTA, null, true, overwrite, append, false);
			dupeWriter=new ByteStreamWriter(ff);
		}
	}
	
	public void process(){
		
		Timer t=new Timer();
		
		boolean dq0=FASTQ.DETECT_QUALITY;
		boolean ti0=FASTQ.TEST_INTERLEAVED;
		
		process2();
		
		FASTQ.DETECT_QUALITY=dq0;
		FASTQ.TEST_INTERLEAVED=ti0;
		
		t.stop();
		
		if(showSpeed){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		if(errorState){
			throw new RuntimeException("Dedupe terminated in an error state; the output may be corrupt.");
		}
	}
	
	public void process2(){
		if(dupeWriter!=null){dupeWriter.start();}
//		assert(false) : out;
		Timer t=new Timer();
		
		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			Shared.printMemory();
			outstream.println();
		}
		
		processMatches(t);
		
		forceTrimLeft=forceTrimRight=-1;
		qTrimLeft=qTrimRight=false;
		
		
		if(absorbContainment){
			processContainments(t);
		}
		
		if(dupeWriter!=null){dupeWriter.poisonAndWait();}
		
		if(findOverlaps){
			findOverlaps(t);
			
			killAffixMaps();
			
			if(processClusters || renameClusters || maxSpanningTree){codeMap=null;}
			
			if(maxSpanningTree){
				processMst(t);
			}
			
			if(processClusters){
				processClusters(t, absorbOverlap);
			}
//			if(renameClusters){
//				renameClusters(t);
//			}
//			assert(false) : (codeMap==null)+", "+(affixMap1==null)+", "+processedClusters;
		}
		
		outstream.println("Input:                  \t"+readsProcessed+" reads \t\t"+basesProcessed+" bases.");

		if(absorbMatch){
			outstream.println("Duplicates:             \t"+matches+" reads ("+Tools.format("%.2f",matches*100.0/readsProcessed)+"%) \t"+
					baseMatches+" bases ("+Tools.format("%.2f",baseMatches*100.0/basesProcessed)+"%)     \t"+collisions+" collisions.");
		}
		if(absorbContainment){
			outstream.println("Containments:           \t"+containments+" reads ("+Tools.format("%.2f",containments*100.0/readsProcessed)+"%) \t"+
					baseContainments+" bases ("+Tools.format("%.2f",baseContainments*100.0/basesProcessed)+"%)    \t"+containmentCollisions+" collisions.");
		}
		if(findOverlaps){
			outstream.println("Overlaps:               \t"+overlaps+" reads ("+Tools.format("%.2f",overlaps*100.0/readsProcessed)+"%) \t"+
					baseOverlaps+" bases ("+Tools.format("%.2f",baseOverlaps*100.0/basesProcessed)+"%)    \t"+overlapCollisions+" collisions.");
		}
//		outstream.println("Result:                 \t"+(addedToMain-containments)+" reads \t\t"+(basesProcessed-baseMatches-baseContainments)+" bases.");
		
		long outReads=(addedToMain-containments);
		if(UNIQUE_ONLY){outReads=readsProcessed-matches-containments;}
		long outBases=(basesProcessed-baseMatches-baseContainments);
		outstream.println("Result:                 \t"+outReads+" reads ("+Tools.format("%.2f",outReads*100.0/readsProcessed)+"%) \t"+
				outBases+" bases ("+Tools.format("%.2f",outBases*100.0/basesProcessed)+"%)");
		
		outstream.println("");
		
		if(out!=null || clusterFilePattern!=null || outbest!=null || outgraph!=null || outcsf!=null){
			writeOutput(outcsf, t);
		}
		
	}
	
	private void killAffixMaps(){
		if(affixMaps==null){return;}
		for(int i=0; i<numAffixMaps; i++){
			if(affixMaps[i]!=null){affixMaps[i].clear();}
			affixMaps[i]=null;
		}
		affixMap1=null;
		affixMap2=null;
		affixMaps=null;
	}
	
	private ConcurrentReadInputStream[] makeCrisArray(ArrayList<Read> list){
		final ConcurrentReadInputStream[] array;
		
		if(list!=null){
			array=new ConcurrentReadInputStream[] {new ConcurrentCollectionReadInputStream(list, null, -1)};
			array[0].start(); //This deadlocks if ConcurrentReadInputStream extends Thread rather than spawning a new thread.
		}else{
			array=new ConcurrentReadInputStream[in1.length];
			multipleInputFiles=array.length>1;
			for(int i=0; i<in1.length; i++){
				if(verbose){System.err.println("Creating cris for "+in1[i]);}

				final ConcurrentReadInputStream cris;
				{
					FileFormat ff1=FileFormat.testInput(in1[i], FileFormat.FASTA, null, !multipleInputFiles || ReadWrite.USE_UNPIGZ, true);
					FileFormat ff2=(in2==null || in2.length<=i ? null : FileFormat.testInput(in2[i], FileFormat.FASTA, null, !multipleInputFiles || ReadWrite.USE_UNPIGZ, true));
					cris=ConcurrentReadInputStream.getReadInputStream(maxReads, ff1.samOrBam(), ff1, ff2);
					cris.start();
					if(cris.paired()){
						THREADS=1;//Temp fix for losing reads when multithreaded and paired
						if(absorbContainment){
							System.err.println("Set absorbContainment to false because it is not currently supported for paired reads.");
							absorbContainment=false;
						}
					}
				}
				array[i]=cris;
			}
		}
		return array;
	}
	
	private void processMatches(Timer t){
		crisa=makeCrisArray(null);
		
		ArrayList<HashThread> alht=new ArrayList<HashThread>(THREADS);
		for(int i=0; i<THREADS; i++){alht.add(new HashThread(true, (absorbContainment|findOverlaps), absorbMatch, false, false));}
		for(HashThread ht : alht){ht.start();}
		for(HashThread ht : alht){
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			matches+=ht.matchesT;
			collisions+=ht.collisionsT;
			containments+=ht.containmentsT;
			containmentCollisions+=ht.containmentCollisionsT;
			baseContainments+=ht.baseContainmentsT;
			baseMatches+=ht.baseMatchesT;
			addedToMain+=ht.addedToMainT;
			readsProcessed+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
		}
		alht.clear();
		
		if(verbose){System.err.println("Attempting to close input streams (1).");}
		for(ConcurrentReadInputStream cris : crisa){
			errorState|=ReadWrite.closeStream(cris);
		}
		crisa=null;
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Found "+matches+" duplicates.");
			outstream.println("Finished exact matches.    Time: "+t);
			Shared.printMemory();
			if(verbose){outstream.println(affixMap1);}
			outstream.println();
			t.start();
		}
	}
	
	private void processContainments(Timer t){
		ArrayList<Read> list=new ArrayList<Read>((int)addedToMain);
		for(ArrayList<Unit> alu : codeMap.values()){
			for(Unit u : alu){
				assert(u.r.mate==null) : "Containments are not currently supported with paired reads.";
				if(u.valid() && u.r.pairnum()==0){list.add(u.r);}
			}
		}

		//	if(minLengthPercent>0){
		//		if(verbose){System.err.println("Sorting.");}
		//		Shared.sort(list, comparator);
		//		Collections.reverse(list);
		//		assert(list.isEmpty() || list.get(0).length()<=list.get(list.size()-1).length()) :
		//			list.get(0).length()+", "+list.get(list.size()-1).length();
		//	}
		
		crisa=makeCrisArray(subsetMode ? null : list);

		ArrayList<HashThread> alht=new ArrayList<HashThread>(THREADS);
		for(int i=0; i<THREADS; i++){alht.add(new HashThread(false, false, false, true, false));}

		for(HashThread ht : alht){ht.start();}
		for(HashThread ht : alht){
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			assert(ht.matchesT==0);
			assert(ht.collisionsT==0);
			assert(ht.baseMatchesT==0);
			assert(ht.addedToMainT==0);
//			assert(ht.readsProcessedT==0);
//			assert(ht.basesProcessedT==0);
			//		matches+=ht.matchesT;
			//		collisions+=ht.collisionsT;
			containments+=ht.containmentsT;
			containmentCollisions+=ht.containmentCollisionsT;
			baseContainments+=ht.baseContainmentsT;
			//		baseMatches+=ht.baseMatchesT;
			//		addedToMain+=ht.addedToMainT;
			//		readsProcessed+=ht.readsProcessedT;
			//		basesProcessed+=ht.basesProcessedT;
		}
		alht.clear();
		if(verbose){System.err.println("Attempting to close input streams (2).");}
		for(ConcurrentReadInputStream cris : crisa){
			errorState|=ReadWrite.closeStream(cris);
		}

		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Found "+containments+" contained sequences.");
			outstream.println("Finished containment.      Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		crisa=null;
		if(!findOverlaps){
			killAffixMaps();
		}

		long x=removeInvalid(list);
		list.clear();

		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Removed "+x+" invalid entries.");
			outstream.println("Finished invalid removal.  Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
	}
	
	private void findOverlaps(Timer t){
		
		ArrayList<Read> list=new ArrayList<Read>((int)addedToMain);
		for(ArrayList<Unit> alu : codeMap.values()){
			for(Unit u : alu){
				if(u.valid() && u.r.pairnum()==0){
					u.unitID=list.size();
					list.add(u.r);
					if(u.r.mate!=null){
						Unit u2=(Unit)u.r.mate.obj;
						u2.unitID=u.unitID;
					}
				}else{
					u.unitID=Integer.MAX_VALUE;
				}
			}
		}
		
		if(preventTransitiveOverlaps){
			clusterNumbers=new AtomicIntegerArray(list.size());
			for(int i=0; i<clusterNumbers.length(); i++){
				clusterNumbers.set(i, i);
			}
		}
		
		crisa=makeCrisArray(subsetMode ? null : list);
		
		ArrayList<HashThread> alht=new ArrayList<HashThread>(THREADS);
		for(int i=0; i<THREADS; i++){alht.add(new HashThread(false, false, false, false, true));}
		
		for(HashThread ht : alht){ht.start();}
		for(HashThread ht : alht){
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			assert(ht.matchesT==0);
			assert(ht.collisionsT==0);
			assert(ht.baseMatchesT==0);
			assert(ht.addedToMainT==0);
			
			overlaps+=ht.overlapsT;
			baseOverlaps+=ht.baseOverlapsT;
			overlapCollisions+=ht.overlapCollisionsT;
		}
		alht.clear();
		if(verbose){System.err.println("Attempting to close input streams (3).");}
		for(ConcurrentReadInputStream cris : crisa){
			errorState|=ReadWrite.closeStream(cris);
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Found "+overlaps+" overlaps.");
			outstream.println("Finished finding overlaps. Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		
		crisa=null;
		
		if(makeClusters){
			int intransitive=0, redundant=0;
			assert((intransitive=countIntransitive(t, list, rigorousTransitive))==0);
			assert((redundant=countRedundant(t, list))==0);
			long overlaps=countOverlaps(t, list);
			assert(intransitive==0);
			assert(redundant==0);
//			makeTransitive(t, list, rigorousTransitive);
			if(clusterQueue==null){
				clusterQueue=new ArrayDeque<ArrayList<Unit>>(list.size()/4+1);
				processedClusters=new ArrayList<ArrayList<Unit>>();
			}else{
				assert(clusterQueue.isEmpty());
			}
			makeClusters(t, list);
		}
		
		list.clear();
	}
	
	private long makeTransitive(Timer t, ArrayList<Read> list, boolean rigorous){
		assert(false) : "No longer needed.";
		long added=0;
		for(Read r : list){
			assert(r!=null);
			Unit u=(Unit) r.obj;
			assert(u!=null);
			assert(u.valid());
//			outstream.println("Considering "+r.id+"; valid="+u.valid()+", overlaps="+(u.overlapList==null ? "null" : u.overlapList.size()));
			if(u.valid()){

				if(u.overlapList!=null){
					for(Overlap o : u.overlapList){
						Unit u2=(o.u1==u ? o.u2 : o.u1);
						assert(u2!=u);
						if(u2.overlapList==null){
							u2.overlapList=new ArrayList<Overlap>(2);
							u2.overlapList.add(o);
						}else{
							boolean found=false;
							if(rigorous){
								found=u2.overlapList.contains(o);
							}else{
								for(Overlap o2 : u2.overlapList){
									if(o2.u1==u || o2.u2==u){found=true; break;}
								}
							}
							if(!found){
								added++;
								u2.overlapList.add(o);
							}
						}
					}
				}
			}
		}
		
		for(Read r : list){
			Unit u=(Unit) r.obj;
			if(u.valid()){
				assert(u.isTransitive());
			}
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Added overlaps: "+added);
			outstream.println("Made overlaps transitive.  Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		return added;
	}
	
	private int countIntransitive(Timer t, ArrayList<Read> list, boolean rigorous){
		if(!countTransitive){return 0;}
		int transitive=0, intransitive=0;
		for(Read r : list){
			assert(r!=null);
			Unit u=(Unit) r.obj;
			assert(u!=null);
			assert(u.valid());
//			outstream.println("Considering "+r.id+"; valid="+u.valid()+", overlaps="+(u.overlapList==null ? "null" : u.overlapList.size()));
			if(u.valid()){
				if(rigorous ? u.isPerfectlyTransitive() : u.isTransitive()){
					transitive++;
				}else{
					intransitive++;
				}
			}
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Intransitive:   "+intransitive+", \ttransitive: "+transitive);
			outstream.println("Checked transitivity.      Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		
		return intransitive;
	}
	
	private int countRedundant(Timer t, ArrayList<Read> list){
		if(!countRedundant){return 0;}
		int redundant=0, nonredundant=0;
		for(Read r : list){
			assert(r!=null);
			Unit u=(Unit) r.obj;
			assert(u!=null);
			assert(u.valid());
//			outstream.println("Considering "+r.id+"; valid="+u.valid()+", overlaps="+(u.overlapList==null ? "null" : u.overlapList.size()));
			if(u.valid()){
				if(u.isNonRedundant()){
					nonredundant++;
				}else{
					redundant++;
				}
			}
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Redundant:      "+redundant+",  \tnonredundant: "+nonredundant);
			outstream.println("Checked redundancy.        Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		return redundant;
	}
	
	private long countOverlaps(Timer t, ArrayList<Read> list){
		
		long overlaps=0, length=0;
		for(Read r : list){
			assert(r!=null);
			Unit u=(Unit) r.obj;
			assert(u!=null);
			assert(u.valid());
//			outstream.println("Considering "+r.id+"; valid="+u.valid()+", overlaps="+(u.overlapList==null ? "null" : u.overlapList.size()));
			if(u.valid() && u.overlapList!=null){
				for(Overlap o : u.overlapList){
					overlaps++;
					length+=o.overlapLen;
				}
			}
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Overlaps:       "+overlaps+",  \tlength: "+length);
			outstream.println("Counted overlaps.          Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		return overlaps;
	}
	
	private long fillClusterSizeMatrix(ArrayList<ArrayList<Unit>> clusters, long[][] clusterSize){
		int max=0;
		for(ArrayList<Unit> cluster : clusters){
			final int cs=Tools.min(clusterSize.length-1, cluster.size());
			{
				long reads=0, bases=0;
				for(Unit u2 : cluster){
					reads++;
					bases+=u2.length();
				}
				clusterSize[0][cs]++;
				clusterSize[1][cs]+=reads;
				clusterSize[2][cs]+=bases;
			}
			max=Tools.max(max, cluster.size());
		}
		return max;
	}
	
	private void makeClusters(Timer t, ArrayList<Read> list){
		
		final int clusterlen=70000;
		long[][] clusterSize=new long[3][clusterlen];
		int max=0;
		for(Read r : list){
			Unit u=(Unit) r.obj;

			if(!u.clustered()){
				ArrayList<Unit> cluster=u.makeCluster();
				if(cluster.size()>2){cluster.trimToSize();}
				if(cluster.size()==1 || (!processClusters && !maxSpanningTree)){processedClusters.add(cluster);}
				else{clusterQueue.add(cluster);}
				final int cs=Tools.min(clusterlen-1, cluster.size());
				{
					long reads=0, bases=0;
					for(Unit u2 : cluster){
						reads++;
						bases+=u2.length();
					}
					clusterSize[0][cs]++;
					clusterSize[1][cs]+=reads;
					clusterSize[2][cs]+=bases;
				}
				max=Tools.max(max, cluster.size());
			}
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println(toClusterSizeString(clusterSize));
			outstream.println("\nLargest:          "+max);
			outstream.println("Finished making clusters.  Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
		
		
		long x=removeInvalid(list);
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Removed "+x+" invalid entries.");
			outstream.println("Finished invalid removal.  Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
	}
	
	private String toClusterSizeString(long[][] clusterSizeMatrix){
		
		long[] clusterSize=clusterSizeMatrix[0];
		long[] clusterReads=clusterSizeMatrix[1];
		long[] clusterBases=clusterSizeMatrix[2];
		
		long totalClusters=shared.Vector.sum(clusterSize);
		
		long bigClusters=0;
		for(int i=minClusterSize; i<clusterSize.length; i++){
			bigClusters+=clusterSize[i];
		}
		
		final int spaces=18;
		final int spaces2=spaces*2, spaces3=spaces*3;

		final StringBuilder sb=new StringBuilder(100), sb2=new StringBuilder(1000);
		sb2.append("Clusters:");
		while(sb2.length()<spaces){sb2.append(' ');}
		sb2.append(totalClusters+(minClusterSize<2 ? "" : " ("+bigClusters+" of at least size "+minClusterSize+")")+"\n");
		
		sb.append("Size Range");
		while(sb.length()<spaces){sb.append(' ');}
		sb.append("Clusters");
		while(sb.length()<spaces2){sb.append(' ');}
		sb.append("Reads");
		while(sb.length()<spaces3){sb.append(' ');}
		sb.append("Bases");
		
		sb2.append('\n');
		sb2.append(sb);
		sb.setLength(0);
		
		for(int i=0; i<clusterSize.length-1; i=Tools.max(i+1, i*2)){
			int a=i+1, b=i*2;
			if(i<2){
				sb.append(a);
				while(sb.length()<spaces){sb.append(' ');}
				sb.append(clusterSize[a]);
				while(sb.length()<spaces2){sb.append(' ');}
				sb.append(clusterReads[a]);
				while(sb.length()<spaces3){sb.append(' ');}
				sb.append(clusterBases[a]);
			}else if(b>=clusterSize.length){
				long x=shared.Vector.sum(clusterSize, a, clusterSize.length-1);
				long y=shared.Vector.sum(clusterReads, a, clusterSize.length-1);
				long z=shared.Vector.sum(clusterBases, a, clusterSize.length-1);
				if(x>0){
					sb.append(a+"+");
					while(sb.length()<spaces){sb.append(' ');}
					sb.append(x);
					while(sb.length()<spaces2){sb.append(' ');}
					sb.append(y);
					while(sb.length()<spaces3){sb.append(' ');}
					sb.append(z);
				}
			}else{
				long x=shared.Vector.sum(clusterSize, a, b);
				long y=shared.Vector.sum(clusterReads, a, b);
				long z=shared.Vector.sum(clusterBases, a, b);
				if(x>0){
					sb.append(a+"-"+b);
					while(sb.length()<spaces){sb.append(' ');}
					sb.append(x);
					while(sb.length()<spaces2){sb.append(' ');}
					sb.append(y);
					while(sb.length()<spaces3){sb.append(' ');}
					sb.append(z);
				}
			}
			if(sb.length()>0){
				sb2.append('\n');
				sb2.append(sb);
				sb.setLength(0);
			}
		}
		return sb2.toString();
	}
	
	private void renameClusters(Timer t){
		assert(false) : "This is now unused; renaming is done at output time.";
		int cnum=0;
		final StringBuilder sb=new StringBuilder(64);
		for(ArrayList<Unit> alu : processedClusters){
			for(int i=0; i<alu.size(); i++){
				Unit u=alu.get(i);
				Read r=u.r;
				sb.append("Cluster ");
				sb.append(cnum);
				sb.append(",contig ");
				sb.append(i);
				if(u.offsetValid()){
					sb.append(",pos ");
					sb.append(u.offset());
				}
				r.id=sb.toString();
				sb.setLength(0);
			}
		}
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Finished cluster renaming. Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
	}
	
	private void processMst(Timer t){
		
		if(DISPLAY_PROGRESS){outstream.println("Converting to Maximum Spanning Tree.");}
		
		ArrayList<MstThread> alct=new ArrayList<MstThread>(THREADS);
		for(int i=0; i<THREADS; i++){
			alct.add(new MstThread());
		}
		
		long overlapsRemoved=0;
		long overlapBasesRemoved=0;
		long overlapsRetained=0;
		long overlapBasesRetained=0;
		
		for(MstThread ct : alct){ct.start();}
		for(MstThread ct : alct){
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			overlapsRemoved+=ct.overlapsRemovedT;
			overlapBasesRemoved+=ct.overlapBasesRemovedT;
			overlapsRetained+=ct.overlapsRetainedT;
			overlapBasesRetained+=ct.overlapBasesRetainedT;
		}
		assert(clusterQueue.isEmpty());
		if(processClusters){
			for(MstThread ct : alct){
				clusterQueue.addAll(ct.processedT);
				ct.processedT.clear();
				ct.processedT=null;
			}
		}else{
			for(MstThread ct : alct){
				processedClusters.addAll(ct.processedT);
				ct.processedT.clear();
				ct.processedT=null;
			}
			clusterQueue=null;
		}
		alct.clear();
		
		assert(affixMaps==null);
		killAffixMaps();
		
		
		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Removed "+(overlapsRemoved)+" edges ("+overlapBasesRemoved+" bases).");
			outstream.println("Retained "+(overlapsRetained)+" edges ("+overlapBasesRetained+" bases).");
			
//			outstream.println("\nAfter conversion to Maximum Spanning Tree:");
//			final int[] clusterSize=new int[8200];
//			int max=0;
//			for(ArrayList<Unit> cluster : processedClusters){
//				clusterSize[Tools.min(clusterSize.length-1, cluster.size())]++;
//				max=Tools.max(max, cluster.size());
//			}
//			outstream.println(toClusterSizeString(clusterSize));
//			outstream.println("Largest:          "+max);

			outstream.println("Finished MST conversion.   Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
	}
	
	private void processClusters(Timer t, boolean mergeClusters){
		
		ArrayList<ClusterThread> alct=new ArrayList<ClusterThread>(THREADS);
		for(int i=0; i<THREADS; i++){
			alct.add(new ClusterThread(fixMultiJoins, removeCycles, fixOffsetContradictions,
					mergeClusters, mergeClusters, mergeClusters));
		}
		
		long leafMerges=0;
		long innerMerges=0;
		long leafBaseMerges=0;
		long innerBaseMerges=0;
		
		long multiJoinFailures=0;
		long multiJoinsFound=0;
		long multiJoinBasesFound=0;
		long unitsFlipped=0;
		long overlapsFlipped=0;
		long offsetContradictoryOverlaps=0;
		long offsetContradictoryClusters=0;
		long cycleOverlaps=0;
		long cycleClusters=0;
		
		for(ClusterThread ct : alct){ct.start();}
		for(ClusterThread ct : alct){
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			leafMerges+=ct.leafMergesT;
			innerMerges+=ct.innerMergesT;
			leafBaseMerges+=ct.leafBaseMergesT;
			innerBaseMerges+=ct.innerBaseMergesT;
			
			multiJoinFailures+=ct.multiJoinFailuresT;
			multiJoinsFound+=ct.multiJoinsFoundT;
			multiJoinBasesFound+=ct.multiJoinBasesFoundT;
			unitsFlipped+=ct.unitsFlippedT;
			overlapsFlipped+=ct.overlapsFlippedT;
			offsetContradictoryOverlaps+=ct.offsetContradictoryOverlapsT;
			offsetContradictoryClusters+=ct.offsetContradictoryClustersT;
			cycleOverlaps+=ct.cycleOverlapsT;
			cycleClusters+=ct.cycleClustersT;
		}
		alct.clear();
		
		assert(affixMaps==null || affixMap1==null);
		killAffixMaps();
		
		assert(clusterQueue.isEmpty());
		clusterQueue=null;
		
		if(DISPLAY_PROGRESS){
			t.stop();
			if(fixMultiJoins){
				outstream.println("Found "+(multiJoinsFound)+" multijoins ("+multiJoinBasesFound+" bases).");
				outstream.println("Experienced "+(multiJoinFailures)+" multijoin removal failures.");
			}
			if(fixOffsetContradictions){
				outstream.println("Found "+(offsetContradictoryClusters)+" clusters ("+offsetContradictoryOverlaps+" overlaps) with contradictory offset cycles.");
			}
			outstream.println("Found "+(cycleClusters)+" clusters ("+cycleOverlaps+" overlaps) with remaining cycles.");
			if(absorbOverlap){
				outstream.println("Merged "+(leafMerges)+" leaves ("+leafBaseMerges+" bases).");
				outstream.println("Merged "+(innerMerges)+" nonleaves ("+innerBaseMerges+" bases).");
			}
			
			outstream.println("\nAfter processing clusters:");
			final long[][] clusterSize=new long[3][70000];
			final long max=fillClusterSizeMatrix(processedClusters, clusterSize);
			outstream.println(toClusterSizeString(clusterSize));
			outstream.println("\nLargest:          "+max);

			outstream.println("Finished processing.       Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
	}
	
	private long removeInvalid(ArrayList<Read> list){
		final LongM keym=new LongM();
		long removedC=0, removedP=0, removedS=0, invalid=0;
		
		for(int j=0, lim=list.size(); j<lim; j++){
			final Read r=list.get(j);
			final Unit u=(Unit)r.obj;
			
			if(!u.valid()){
				
				invalid++;
				
				if(codeMap!=null && !codeMap.isEmpty()){
					Long key=u.code1;
					ArrayList<Unit> alu=codeMap.get(key);
					if(alu!=null){
						int valid=0;
						for(int i=alu.size()-1; i>=0; i--){
							Unit u2=alu.get(i);
							if(u2==null || !u2.valid()){
								alu.remove(i);
								removedC++;
							}
							else{valid++;}
						}
						if(valid==0){codeMap.remove(key);}
					}
				}
				
				if(affixMap1!=null && !affixMap1.isEmpty()){
					{
						keym.set(u.prefix1);
						ArrayList<Unit> alu=affixMap1.get(keym);
						if(alu!=null){
							int valid=0;
							for(int i=alu.size()-1; i>=0; i--){
								Unit u2=alu.get(i);
								if(u2==null || !u2.valid()){
									alu.remove(i);
									removedP++;
								}
								else{valid++;}
							}
							if(valid==0){affixMap1.remove(keym);}
						}
					}
					if(storeSuffix){
						keym.set(u.suffix1);
						ArrayList<Unit> alu=affixMap1.get(keym);
						if(alu!=null){
							int valid=0;
							for(int i=alu.size()-1; i>=0; i--){
								Unit u2=alu.get(i);
								if(u2==null || !u2.valid()){
									alu.remove(i);
									removedS++;
								}
								else{valid++;}
							}
							if(valid==0){affixMap1.remove(keym);}
						}
					}
				}
				if(affixMap2!=null && !affixMap2.isEmpty()){
					if(u.prefix2!=-1){
						keym.set(u.prefix2);
						ArrayList<Unit> alu=affixMap2.get(keym);
						if(alu!=null){
							int valid=0;
							for(int i=alu.size()-1; i>=0; i--){
								Unit u2=alu.get(i);
								if(u2==null || !u2.valid()){
									alu.remove(i);
									removedP++;
								}
								else{valid++;}
							}
							if(valid==0){affixMap2.remove(keym);}
						}
					}
					if(storeSuffix && u.suffix2!=-1){
						keym.set(u.suffix2);
						ArrayList<Unit> alu=affixMap2.get(keym);
						if(alu!=null){
							int valid=0;
							for(int i=alu.size()-1; i>=0; i--){
								Unit u2=alu.get(i);
								if(u2==null || !u2.valid()){
									alu.remove(i);
									removedS++;
								}
								else{valid++;}
							}
							if(valid==0){affixMap2.remove(keym);}
						}
					}
					
				}
				
				list.set(j, null);
			}
		}
		
		if(invalid>0){
			Tools.condenseStrict(list);
		}
		if(verbose){
			outstream.println("Removed invalids: "+removedC+", "+removedP+", "+removedS);
		}
		return invalid;
	}
	
	
	private static ArrayList<Read> addToArray(HashMap<Long, ArrayList<Unit>> codeMap, boolean sort, boolean clear, long outNum){
		assert(outNum<=Integer.MAX_VALUE);
		if(verbose){System.err.println("Making list.");}
		ArrayList<Read> list=new ArrayList<Read>((int)outNum);
		if(verbose){System.err.println("Adding.");}
		for(ArrayList<Unit> alu : codeMap.values()){
			for(Unit u : alu){
				if(u.valid() && u.r.pairnum()==0){list.add(u.r);}
			}
			if(clear){alu.clear();}
		}
		if(clear){codeMap.clear();}
		
		if(sort){
			if(verbose){System.err.println("Sorting.");}
			Shared.sort(list, comparator);
//			if(ascending){
//				Collections.reverse(list);
//				assert(list.isEmpty() || list.get(0).length()<=list.get(list.size()-1).length()) :
//					list.get(0).length()+", "+list.get(list.size()-1).length();
//			}else{
//				assert(list.isEmpty() || list.get(0).length()>=list.get(list.size()-1).length()) :
//					list.get(0).length()+", "+list.get(list.size()-1).length();
//			}
		}
		assert(list.size()==outNum || list.size()*2L==outNum || UNIQUE_ONLY) : list.size()+", "+outNum;
		return list;
	}
	
	private void writeOutput(String clusterStatsFile, Timer t){
//		verbose=true;
//		assert(false) : (processedClusters==null)+", "+(processedClusters.isEmpty())+", "+outgraph+", "+out+", "+clusterFilePattern;
		if(processedClusters==null || processedClusters.isEmpty()){
			
			if(out!=null || clusterFilePattern!=null){

				ArrayList<Read> list=addToArray(codeMap, sort, true, addedToMain-containments);
				codeMap=null;

				if(sort){
					if(DISPLAY_PROGRESS){
						t.stop();
						outstream.println("Sorted output.             Time: "+t);
						Shared.printMemory();
						outstream.println();
						t.start();
					}
				}

				writeOutput(list);
			}
		}else{
			if(outgraph!=null){
				writeGraph(outgraph, processedClusters);
			}
			if(out!=null || clusterFilePattern!=null || clusterStatsFile!=null || outbest!=null){
				writeOutputClusters(clusterStatsFile, processedClusters);
			}
		}

		if(DISPLAY_PROGRESS){
			t.stop();
			outstream.println("Printed output.            Time: "+t);
			Shared.printMemory();
			outstream.println();
			t.start();
		}
	}
	

	
	private void writeOutput(ArrayList<Read> list){
		
		final ByteStreamWriter tsw=(out==null ? null : new ByteStreamWriter(out, overwrite, append, true));
		
		if(verbose){System.err.println("Writing from array.");}
		tsw.start();
		
		HashSet<String> names=((uniqueNames && storeName) ?
				new HashSet<String>(Tools.min(Integer.MAX_VALUE, Tools.max((int)addedToMain, (int)(addedToMain*1.35)))) : null);
		long rid=0;
		for(int x=0; x<list.size(); x++){
			Read r=list.get(x);
			list.set(x, null);
			
			if(r.mate!=null && r.pairnum()!=0){r=r.mate;}
			
			assert(r.mate==null || r.mate.discarded()==r.discarded());
			
			if(!r.discarded()){
				rid++;
				
				for(int i=0; r!=null && i<2; i++){
					if(multipleInputFiles){r.numericID=rid;}
					if(names!=null){
						String name=(r.id==null ? ""+r.numericID : r.id);
						if(names.contains(name)){
							for(long j=0; j<Integer.MAX_VALUE; j++){
								String name2=name+"_dd"+j;
								if(r.mate!=null){name2+=(" /"+(i+1));}
								if(!names.contains(name2)){
									r.id=name2;
									names.add(name2);
									break;
								}
							}
						}else{
							names.add(name);
						}
					}
					tsw.println(r);
					r.setDiscarded(true);
					r=r.mate;
				}
			}
		}
		if(verbose){System.err.println("Shutting down tsw "+tsw.fname);}
		tsw.poisonAndWait();
	}
	

	private void writeOutputClusters(String clusterStatsFile, ArrayList<ArrayList<Unit>> clist){

//		Shared.sort(clist, CLUSTER_LENGTH_COMPARATOR); //Causes a crash: java.lang.System.arraycopy(Native Method)
		clist.sort(CLUSTER_LENGTH_COMPARATOR);
		
		if(verbose){System.err.println("Writing clusters.");}
		
		final ByteStreamWriter tswAll=(out==null ? null : new ByteStreamWriter(out, overwrite, append, true));
		if(tswAll!=null){tswAll.start();}
		ByteStreamWriter tswCluster=null;
		ByteStreamWriter tswBest=null;
		
		if(outbest!=null){
			tswBest=new ByteStreamWriter(outbest, overwrite, append, true);
			tswBest.start();
		}
		
		TextStreamWriter csf=null;
		if(clusterStatsFile!=null){
			csf=new TextStreamWriter(clusterStatsFile, overwrite, false, false);
			csf.start();
			csf.print("#Name\tsize\t"+nmerLength+"-mer frequencies\n");
		}
		
		HashSet<String> names=((uniqueNames && storeName) ?
				new HashSet<String>(Tools.min(Integer.MAX_VALUE, Tools.max((int)addedToMain, (int)(addedToMain*1.35)))) : null);
		long rid=0;
		final long[] nmerCounts=new long[maxNmer+1];
		
		final StringBuilder sb=new StringBuilder(64);
		
		for(int cnum=0; cnum<clist.size(); cnum++){
			final ArrayList<Unit> alu=clist.get(cnum);
//			clist.set(cnum, null); //This breaks subsequent output processing
			
			if(alu.size()<minClusterSize){
				if(verbose){System.err.println("Ignoring small cluster "+cnum+", size "+alu.size());}

				if(csf!=null && alu.size()>=minClusterSizeForStats){
					float[] profile=makeNmerProfile(alu, nmerCounts);
					sb.append("Cluster_");
					sb.append(cnum);
					sb.append('\t');
					sb.append(alu.size());
					sb.append('\t');
					for(float f : profile){
						sb.append(Tools.format("%.5f ", f));
					}
					sb.setCharAt(sb.length()-1, '\n');
					csf.print(sb.toString());
					sb.setLength(0);
				}
			}else{
				if(verbose){System.err.println("Writing cluster "+cnum+", size "+alu.size());}
				
				if(clusterFilePattern!=null){
					if(tswCluster!=null){
						if(verbose){System.err.println("Shutting down tswCluster "+tswCluster.fname);}
						tswCluster.poisonAndWait();
						tswCluster=null;
					}
					tswCluster=new ByteStreamWriter(clusterFilePattern.replaceFirst("%", ""+cnum), overwrite, append, true);
					if(verbose){System.err.println("Starting tswCluster "+tswCluster.fname);}
					tswCluster.start();
				}

				if(csf!=null && alu.size()>=minClusterSizeForStats){
					float[] profile=makeNmerProfile(alu, nmerCounts);
					sb.append("Cluster_");
					sb.append(cnum);
					sb.append('\t');
					sb.append(alu.size());
					sb.append('\t');
					for(float f : profile){
						sb.append(Tools.format("%.5f ", f));
					}
					sb.setCharAt(sb.length()-1, '\n');
					csf.print(sb.toString());
					sb.setLength(0);
				}
				
				if(pickBestRepresentativePerCluster){
					pickBestRepresenative(alu, true);
				}
				
				if(outbest!=null){
					Unit u=pickBestRepresenative((ArrayList<Unit>)alu.clone(), false);
					tswBest.println(u.r);
					if(u.r.mate!=null){tswBest.println(u.r.mate);}
				}

				for(int contig=0; contig<alu.size(); contig++){
					final Unit u0=alu.get(contig);
					alu.set(contig, null);
					Read r=u0.r;
					if(r.mate!=null && r.pairnum()!=0){r=r.mate;}

					if(!r.discarded()){
						rid++;

						for(int i=0; r!=null && i<2; i++){
							assert(r.pairnum()==i) : i+", "+r.pairnum()+", "+(r.mate==null ? 9 : r.mate.pairnum());
							Unit u=(Unit)r.obj;
							if(verbose){System.err.println("Writing read "+r.id);}
							r.numericID=rid;
							if(renameClusters){
								sb.append("Cluster_");
								sb.append(cnum);
								sb.append(",contig_");
								sb.append(contig);
								if(u.offsetValid()){
									sb.append(",pos_");
									sb.append(u.offset());
								}
								if(r.mate!=null){sb.append(" /"+(i+1));}
								r.id=(r.id==null ? sb.toString() : r.id+"\t"+sb);
								sb.setLength(0);
							}else if(names!=null){
								String name=(r.id==null ? ""+r.numericID : r.id);
								if(names.contains(name)){
									for(long j=0; j<Integer.MAX_VALUE; j++){
										String name2=name+"_dd"+j;
										if(!names.contains(name2)){
											r.id=name2;
											names.add(name2);
											break;
										}
									}
								}else{
									names.add(name);
								}
							}
							if(tswAll!=null){tswAll.println(r);}
							if(tswCluster!=null){tswCluster.println(r);}
							r.setDiscarded(true);
							r=r.mate;
						}
					}
				}
			}
		}
		if(csf!=null){
			if(verbose){System.err.println("Shutting down csf "+csf.fname);}
			csf.poisonAndWait();
		}
		if(tswBest!=null){
			if(verbose){System.err.println("Shutting down tswBest "+tswBest.fname);}
			tswBest.poisonAndWait();
		}
		if(tswAll!=null){
			if(verbose){System.err.println("Shutting down tswAll "+tswAll.fname);}
			tswAll.poisonAndWait();
		}
		if(tswCluster!=null){
			if(verbose){System.err.println("Shutting down tswCluster "+tswCluster.fname);}
			tswCluster.poisonAndWait();
		}
	}
	

	private void writeGraph(String graphFile, ArrayList<ArrayList<Unit>> clist){
//		Shared.sort(clist, CLUSTER_LENGTH_COMPARATOR);
		clist.sort(CLUSTER_LENGTH_COMPARATOR);
		
		if(verbose){System.err.println("Writing overlap graph.");}
		
		final TextStreamWriter tsw=(graphFile==null ? null : new TextStreamWriter(graphFile, overwrite, append, true));
		if(tsw!=null){
			tsw.start();
			tsw.print("digraph G {\n");
		}
		
		for(int cnum=0; cnum<clist.size(); cnum++){
			final ArrayList<Unit> alu=clist.get(cnum);
//			clist.set(cnum, null); //This breaks subsequent output processing
//			Shared.sort(alu); //TODO: Remove

			if(alu.size()<minClusterSize){
				if(verbose){System.err.println("Ignoring small cluster "+cnum+", size "+alu.size());}
			}else{
				if(verbose){System.err.println("Processing cluster "+cnum+", size "+alu.size());}

				for(int contig=0; contig<alu.size(); contig++){
					final Unit u0=alu.get(contig);
//					alu.set(contig, null); //This breaks subsequent output processing
					Read r=u0.r;
					if(r.mate!=null && r.pairnum()!=0){r=r.mate;}
					
					{
						for(int i=0; r!=null && i<2; i++){
							assert(r.pairnum()==i) : i+", "+r.pairnum()+", "+(r.mate==null ? 9 : r.mate.pairnum());
							Unit u=(Unit)r.obj;
							if(verbose){System.err.println("Writing read "+r.id);}
							
							if(tsw!=null){
								tsw.print("\t"+toGraphName(r)+"\n");
								if(r.mate!=null && r.pairnum()==0){
									Read r2=r.mate;
									tsw.print(toGraphName(r)+" -> "+toGraphName(r2)+" [label=mate]");
								}
								if(u.overlapList!=null){
									for(Overlap o : u.overlapList){
										if(u==o.u1){
											Read r2=o.u2.r;
											tsw.print("\t"+toGraphName(r)+" -> "+toGraphName(r2)+" [label=\""+o.toLabel()+"\"]\n");
										}
									}
								}
							}
							r=r.mate;
						}
					}
				}
			}
		}
		
		if(tsw!=null){
			tsw.print("}\n");
			if(verbose){System.err.println("Shutting down tswAll "+tsw.fname);}
			tsw.poisonAndWait();
		}
	}
	
	private static String toGraphName(Read r){
		if(NUMBER_GRAPH_NODES || r.id==null){
			return r.numericID+((ADD_PAIRNUM_TO_NAME || r.mate!=null) ? "."+(r.pairnum()+1) : "");
		}else{
			return r.id.replace(' ','_').replace('\t','_');
		}
	}
	
	private Unit pickBestRepresenative(ArrayList<Unit> alu, boolean clearList){
		if(alu==null || alu.isEmpty()){return null;}
		float[] quality=new float[alu.size()];
		int[] lengths=new int[alu.size()];
		for(int i=0; i<alu.size(); i++){
			Unit u=alu.get(i);
			int len=u.r.length();
			quality[i]=u.r.expectedErrors(true, 0)/len;
			lengths[i]=len;
		}
		Arrays.sort(quality);
		Arrays.sort(lengths);
		int medianLength=lengths[lengths.length/2];
		float bestQuality=quality[0];
		
		float currentBestQuality=9999999;
		Unit best=null;
		for(int i=0; i<alu.size(); i++){
			Unit u=alu.get(i);
			int len=u.r.length();
			float deviation=Tools.absdif(len, medianLength)*1f/(medianLength+1);
			if(deviation<0.05){
				float qual=u.r.expectedErrors(true, 0)/len;
				qual=(qual+.001f)*(1+10*deviation);
				if(qual<currentBestQuality || best==null){
					currentBestQuality=qual;
					best=u;
				}
			}
		}
		if(clearList){
			alu.clear();
			alu.add(best);
		}
		return best;
	}
	
	public static long hash(byte[] bases){
		long code=bases.length;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int mode=(int)(code&31);
			assert(hashcodes[b]!=null) : "Invalid sequence character: '"+(char)b+"'";
			code=code^hashcodes[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	
	public static long hashReversed(byte[] bases){
		long code=bases.length;
		for(int i=bases.length-1; i>=0; i--){
			byte b=bases[i];
			assert(hashcodesR[b]!=null) : "Invalid sequence character: '"+(char)b+"'";
			int mode=(int)(code&31);
			code=code^hashcodesR[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	private static synchronized long[][] makeCodes(int symbols, int modes, long seed){
		Random randy=Shared.threadLocalRandom(seed);
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong();
			}
		}
		return r;
	}
	
	/** Handles IUPAC codes and invalid symbols */
	private static synchronized long[][] makeCodes2(int modes, long seed){
		long[][] r=makeCodes(128, modes, seed);
		
		for(int i=0; i<26; i++){
			char c=(char)('A'+i);
			r[Tools.toLowerCase(c)]=r[c];
		}
		return r;
	}
	
	private void addDupe(Read r){
		if(dupeWriter==null){return;}
		if(r.mate==null || r.pairnum()==0){
			synchronized(dupeWriter){
				dupeWriter.println(r);
				if(r.mate!=null){
					dupeWriter.println(r.mate);
				}
			}
		}
	}
	
	
	private final class MstThread extends Thread{

		public MstThread(){}
		
		@Override
		public void run(){
			
			ArrayList<Unit> cluster=null;
			while((cluster=nextCluster())!=null){
				makeMst(cluster);
				processedT.add(cluster);
			}
				
		}
		
		public void makeMst(ArrayList<Unit> cluster){
			assert(heap.isEmpty());
			unvisit(cluster);
			for(Unit u : cluster){
				u.flags&=~Unit.VISIT_MASK;
				Shared.sort(u.overlapList);
			}
			{
				Unit u=cluster.get(0);
				u.setVisited(true);
				heap.addAll(u.overlapList);
			}
//			assert(false) : cluster.size();
			while(!heap.isEmpty()){
				Overlap o=heap.poll();
				assert(!o.mst());
				if(!o.invalid()){
//					assert(o.u1.overlapList.contains(o)); //slow
//					assert(o.u2.overlapList.contains(o)); //slow
					assert(o.u1.visited() || o.u2.visited());
					final Unit u=(!o.u1.visited() ? o.u1 : !o.u2.visited()? o.u2 : null);
					if(u!=null){
						o.setMst(true);
						u.setVisited(true);
						overlapsRetainedT++;
						overlapBasesRetainedT+=o.overlapLen;
						for(Overlap o2 : u.overlapList){
							if(o2.mst()){
								//do nothing
							}else if(!o2.u1.visited() || !o2.u2.visited()){
								if(heap.size()>=Integer.MAX_VALUE){
									removeInvalid(heap);
								}
								heap.add(o2);
							}else if(!o2.invalid()){
								o2.setInvalid(true);
								overlapsRemovedT++;
								overlapBasesRemovedT+=o2.overlapLen;
							}
						}
					}
				}
			}
			for(Unit u : cluster){
				ArrayList<Overlap> alo=u.overlapList;
				int removed=0;
				for(int i=0; i<alo.size(); i++){
					Overlap o=alo.get(i);
					if(o.invalid()){
						assert(!o.mst());
						alo.set(i, null);
						removed++;
					}else{
						assert(o.mst());
					}
				}
				if(removed>0){
					Tools.condenseStrict(alo);
					alo.trimToSize();
				}
			}
		}
		
		private void removeInvalid(PriorityQueue<Overlap> heap){
			ArrayList<Overlap> valid=new ArrayList<Overlap>(heap.size());
			for(Overlap o : heap){
				if(!o.invalid()){
					assert(!o.u1.visited() || !o.u2.visited());
					valid.add(o);
				}
			}
			heap.clear();
			heap.addAll(valid);
		}


		public long overlapsRemovedT=0;
		public long overlapBasesRemovedT=0;
		public long overlapsRetainedT=0;
		public long overlapBasesRetainedT=0;
		
		private final PriorityQueue<Overlap> heap=new PriorityQueue<Overlap>((1<<16)-1);
		private ArrayList<ArrayList<Unit>> processedT=new ArrayList<ArrayList<Unit>>();
	}
	
	
	/**
	 * Processes clustered sets of reads.
	 * @author Brian Bushnell
	 * @date Aug 9, 2013
	 *
	 */
	private final class ClusterThread extends Thread{
		
		public ClusterThread(boolean fixMultiJoins_, boolean removeCycles_,
				boolean fixOffsetContradictions_, boolean mergeClusters_, boolean mergeLeaves_, boolean mergeInner_){
			fixMultiJoinsT=fixMultiJoins_;
			fixOffsetContradictionsT=fixOffsetContradictions_;
			mergeClustersT=mergeClusters_;
			mergeLeavesT=mergeLeaves_;
			mergeInnerT=mergeInner_;
			
//			assert(false) : fixMultiJoinsT+", "+canonicizeT+", "+fixCanonContradictionsT+", "+mergeLeavesT+", "+mergeInnerT;
			bandy=(maxEdits>0 ? BandedAligner.makeBandedAligner(bandwidth) : null);
//			assert(false) : fixMultiJoinsT+", "+canonicizeT+", "+fixCanonContradictionsT+", "+fixOffsetContradictionsT+", "+mergeClustersT+", "+removeCycles_;
		}
		
		@Override
		public void run(){
			
			final ArrayList<Unit> temp=new ArrayList<Unit>(1000);
			
			ArrayList<Unit> cluster=null;
			while((cluster=nextCluster())!=null){
				
				if(EA){
					for(Unit u : cluster){assert(u.r.mate==null) : "Cluster processing/merging is not supported for paired reads, only cluster generation.";}
				}
				
//				for(Unit u : cluster){assert(!u.visited());}
				unvisit(cluster);
				
				reorderClusterBreadthFirst(cluster);
				int multiJoinCount=findMultiJoinsInCluster(cluster, fixMultiJoinsT);
				
				if(EA){
					for(Unit u : cluster){assert(!u.visited());}
				}
				
				boolean ok=true;
				if(multiJoinCount!=0){
					assert(multiJoinCount>0);
					multiJoinsFoundT+=multiJoinCount;
					if(!fixMultiJoinsT){
						multiJoinFailuresT++;
						ok=false;
					}
				}
				
				int cycleOverlaps=0;
				if(ok){
					cycleOverlaps=findCycles(cluster, removeCycles);
					for(Unit u : cluster){assert(!u.visited());}
					if(cycleOverlaps>0){
						cycleOverlapsT+=cycleOverlaps;
						cycleClustersT++;
					}
				}
				
				int offsetContradictions=0;
				if(ok && fixOffsetContradictionsT){
					if(EA){
						for(Unit u : cluster){
							assert(!u.visited());
							assert(!u.offsetContradiction());
							assert(!u.offsetValid());
							if(u.overlapList!=null){
								for(Overlap o : u.overlapList){
									assert(!o.invalid());
									assert(!o.offsetContradiction());
									assert(o.type==FORWARD) : o;
								}
							}
						}
					}
					offsetContradictions=generateOffsetsBreadthFirst(cluster, temp);
//					System.err.println("Made offsets for cluster of size "+cluster.size()+"; contradictions = "+offsetContradictions+"; set = "+temp.size());
					temp.clear();
					for(Unit u : cluster){assert(!u.visited());}
					if(offsetContradictions>0){
						offsetContradictoryOverlapsT+=offsetContradictions;
						offsetContradictoryClustersT++;
						if(fixOffsetContradictionsT){
							if(verbose){System.err.println("Pruning cluster to remove offset contradictions.");}
							fullyPruneCluster(cluster, temp);
							if(verbose){System.err.println("Resulting size: "+cluster.size());}
							if(EA){
								for(Unit u : cluster){
									assert(!u.visited());
									assert(!u.offsetContradiction());
									assert(u.offsetValid());
									if(u.overlapList!=null){
										for(Overlap o : u.overlapList){
											assert(!o.invalid());
											assert(!o.offsetContradiction());
											assert(o.type==FORWARD) : o;
										}
									}
								}
							}
						}else{
							ok=false;
						}
					}
					if(ok){Shared.sort(cluster, UNIT_OFFSET_COMPARATOR);}
				}
				
				if(ok && absorbOverlap){
					mergeCluster(cluster);
				}
				
				processedClustersT.add(cluster);
				if(processedClustersT.size()>=threadMaxReadsToBuffer){
					synchronized(processedClusters){
						processedClusters.addAll(processedClustersT);
						processedClustersT.clear();
					}
				}
			}
			synchronized(processedClusters){
				processedClusters.addAll(processedClustersT);
				processedClustersT.clear();
			}
		}
		
		private void fullyPruneCluster(ArrayList<Unit> cluster, ArrayList<Unit> temp){
			assert(cluster.size()>1) : cluster.size();
			ArrayList<Unit> pruned=pruneCluster(cluster, true, true, temp);
			assert(temp.isEmpty());
			assert(pruned==null || pruned.size()>0);
			while(pruned!=null){
				ArrayList<Unit> subcluster=pruned;
				for(Unit u : subcluster){
					u.clearVolatileFlags();
					if(u.overlapList!=null){
						for(Overlap o : u.overlapList){
							o.clearVolatileFlags();
						}
					}
				}
				assert(subcluster.size()>0);
				pruned=pruneCluster(subcluster, false, false, temp);
				assert(temp.isEmpty());
				assert(pruned==null || pruned.size()>0);
				assert(subcluster.size()>0);
				if(subcluster.size()==1){
					processedClustersT.add(subcluster);
				}else{
					assert(subcluster.size()>1);
					synchronized(clusterQueue){
						clusterQueue.add(subcluster);
					}
				}
			}
		}
		
		/**
		 * @param cluster
		 */
		private void mergeCluster(ArrayList<Unit> cluster) {
			if(cluster.size()==1){return;}
			if(mergeLeavesT){
				mergeLeaves(cluster);
			}
			if(mergeInnerT){
				mergeInner(cluster);
			}
		}
		
		/**
		 * Finds places in the cluster where two Units are joined by multiple different Overlaps.
		 * Returns number of multijoins found.
		 * @param cluster
		 */
		private int findMultiJoinsInCluster(ArrayList<Unit> cluster, boolean resolveProblems) {
			if(cluster.size()<2){return 0;}
			int totalMultiJoins=0;
			for(Unit ua : cluster){
				ArrayList<Overlap> list=ua.overlapList;
				assert(list!=null);
				if(list.size()>1){
					Shared.sort(list);
					
					int multiJoins=0;
					for(int i=0; i<list.size(); i++){
						Overlap o=list.get(i);
						Unit ub=(o.u1==ua ? o.u2 : o.u1);
						assert(ua!=ub);
						assert(ua==o.u1 || ua==o.u2);
						if(ub.visited()){
							multiJoins++;
							multiJoinBasesFoundT+=o.overlapLen;
							if(!o.multiJoin()){o.setMultiJoin(true);}
							if(resolveProblems){list.set(i, null);}
						}else{
							ub.setVisited(true);
						}
					}
					
					if(multiJoins>0){
						totalMultiJoins+=multiJoins;
						if(resolveProblems){Tools.condenseStrict(list);}
					}
					
					for(int i=0; i<list.size(); i++){
						Overlap o=list.get(i);
						Unit ub=(o.u1==ua ? o.u2 : o.u1);
						assert(ua!=ub);
						assert(ua==o.u1 || ua==o.u2);
						assert(ub.visited());
						ub.setVisited(false);
					}
				}
				
			}
			
			return totalMultiJoins;
		}
		
		private ArrayList<Unit> pruneCluster(ArrayList<Unit> cluster, boolean pruneContradictoryNodes, boolean pruneContradictoryOverlaps, ArrayList<Unit> visited){
			if(verbose){System.err.println("pruneCluster(size="+cluster.size()+", "+pruneContradictoryNodes+", "+pruneContradictoryOverlaps+")");}
			
			//pruneContradictoryOverlaps is less strict than pruneContradictoryNodes
			assert(pruneContradictoryOverlaps || !pruneContradictoryNodes);
			
			for(Unit ua : cluster){
				assert(!ua.visited());
				assert(ua.isPerfectlyTransitive()) : ua;
				if(ua.visited()){ua.setVisited(false);}
			}
			
			int prunedOverlaps=0;
			int visits=1;
			
			{
				final Unit root=cluster.get(0);
				assert(!root.contradiction());
				root.setVisited(true);
				visited.add(root);
			}
			
			for(int i=0; i<visited.size(); i++){
				Unit ua=visited.get(i);
				
				if(ua.visited() && (!ua.contradiction() || !pruneContradictoryNodes)){
					ArrayList<Overlap> list=ua.overlapList;
					if(list!=null){
						int removed=0;
						for(int j=0; j<list.size(); j++){
							Overlap o=list.get(j);
							Unit ub=(o.u1==ua ? o.u2 : o.u1);
							assert(o.u1==ua || o.u2==ua);
							assert(ua!=ub);
							assert(ub.valid());

							assert(!o.offsetContradiction() || (ua.offsetContradiction() || ub.offsetContradiction())) :
								"\n"+o.offsetContradiction()+", "+ua.offsetContradiction()+", "+ub.offsetContradiction();

//							assert(o.contradiction()==(ua.contradiction() && ub.contradiction())) :
//								"\n"+o.canonContradiction()+", "+o.offsetContradiction()+
//								"\n"+ua.canonContradiction()+", "+ua.offsetContradiction()+
//								"\n"+ub.canonContradiction()+", "+ub.offsetContradiction();

							final boolean remove=(pruneContradictoryNodes && ub.contradiction() || (pruneContradictoryOverlaps && o.contradiction()));
							if(!remove && !ub.visited()){
								ub.setVisited(true);
								visited.add(ub);
								visits++;
							}

							if(remove){
								if(!o.invalid()){o.setInvalid(true);}
								list.set(j, null);
								removed++;
								prunedOverlaps++;
							}else{
								assert(!o.invalid());
							}
						}
						if(removed>0){Tools.condenseStrict(list);}
					}
				}
			}
			
			if(verbose){System.err.println("cluster.size()="+cluster.size()+", visits="+visits+", visited.size()="+visited.size());}
			
			final int numUnvisited=cluster.size()-visits;
			ArrayList<Unit> pruned=(numUnvisited==0 ? null : new ArrayList<Unit>(numUnvisited));
			assert(visits==visited.size());
			assert(visits>=1 && visits<=cluster.size());
			
			if(visits<cluster.size()){
				pruned=new ArrayList<Unit>(cluster.size()-visits);
				for(Unit ua : cluster){
					if(!ua.visited()){
						pruned.add(ua);
						ArrayList<Overlap> list=ua.overlapList;
						if(list!=null){
							int removed=0;
							for(int j=0; j<list.size(); j++){
								Overlap o=list.get(j);
								Unit ub=(o.u1==ua ? o.u2 : o.u1);
								assert(o.u1==ua || o.u2==ua);
								assert(ua!=ub);
								assert(ub.valid());

								if(ub.visited() || o.invalid()){
									assert(ub.visited()==o.invalid()) : "\n"+o+"\n"+ub;
									list.set(j, null);
									removed++;
								}
							}
							if(removed>0){Tools.condenseStrict(list);}
						}
					}
				}
				assert(pruned.size()==numUnvisited);
			}else{
				assert(prunedOverlaps==0) : "If this fails then I may need to mark overlaps to remove.";
			}
			for(Unit u : cluster){
				assert(u.isPerfectlyTransitive()) : u;
				if(EA){
					if(u.overlapList!=null){
						for(Overlap o : u.overlapList){assert(!o.invalid());}
					}
				}
				if(u.visited()){u.setVisited(false);}
			}
			cluster.clear();
			cluster.addAll(visited);
			cluster.trimToSize();
			
//			for(Unit u : cluster){
////				assert(u.canonicized());
//				for(Overlap o : u.overlapList){
//					assert(pruned==null || !pruned.contains(o.u1));
//					assert(pruned==null || !pruned.contains(o.u2));
//					assert(cluster.contains(o.u1));
//					assert(cluster.contains(o.u2));
//				}
//			}
//			if(pruned!=null){
//				for(Unit u : pruned){
//					for(Overlap o : u.overlapList){
//						assert(pruned.contains(o.u1));
//						assert(pruned.contains(o.u2));
//						assert(!cluster.contains(o.u1));
//						assert(!cluster.contains(o.u2));
//					}
//				}
//			}
			
			visited.clear();
			return pruned;
		}
		
		/**
		 * Cluster should already be ordered breadth-first
		 * This may fail because removing cycles could change breadth-first traversal, but if it fails, an assertion will be thrown
		 * @param cluster
		 */
		private int findCycles(ArrayList<Unit> cluster, boolean remove){
			
			{
				final Unit root=cluster.get(0);
				assert(root.length()>=cluster.get(cluster.size()-1).length());
				root.setVisited(true);
			}
			int cycles=0;
			
			for(Unit ua : cluster){
				assert(ua.visited());
				ArrayList<Overlap> list=ua.overlapList;
				if(list!=null){
					int removed=0;
					for(int i=0; i<list.size(); i++){
						Overlap o=list.get(i);
						Unit ub=(o.u1==ua ? o.u2 : o.u1);
						assert(o.u1==ua || o.u2==ua);
						assert(ua!=ub);
						assert(ub.valid());

						if(!o.visited()){
							o.setVisited(true);
							if(ub.visited()){
								if(!o.cyclic()){
									o.setCyclic(true);
									cycles++;
								}
							}else{
								ub.setVisited(true);
							}
						}
						if(remove && o.cyclic()){
							list.set(i, null);
							removed++;
						}
					}
					if(removed>0){Tools.condenseStrict(list);}
				}
			}
			
			for(Unit u : cluster){
				if(u.visited()){u.setVisited(false);}
				if(u.overlapList!=null){
					for(Overlap o : u.overlapList){
						if(o.visited()){o.setVisited(false);}
					}
				}
			}
			
			return cycles;
		}
		
		/**
		 * Cluster should already be ordered breadth-first
		 * @param cluster
		 */
		private int generateOffsetsBreadthFirst(ArrayList<Unit> cluster, ArrayList<Unit> temp){


			assert(temp!=null);
			assert(temp.isEmpty());
			{
				final Unit root=cluster.get(0);
				assert(root.length()>=cluster.get(cluster.size()-1).length());
				root.setOffset(0);
				temp.add(root);
			}
			
			int contradictions=0;
			for(int i=0; i<temp.size(); i++){
				Unit u=temp.get(i);
				assert(!u.visited()) : i;
				assert(u.offsetValid() || contradictions>0) : i+", "+temp.size()+", "+contradictions+"\n"+toString(temp);
				if(u.offsetValid() && !u.offsetContradiction()){
					contradictions+=setOffsetsNeighbors(u, temp);
					assert(contradictions==0 || (i>0 && temp.size()>2));
				}
			}
			
			int min=0;
			for(Unit u : temp){
				if(u.visited()){u.setVisited(false);}
				if(u.overlapList!=null){
					for(Overlap o : u.overlapList){
						if(o.visited()){o.setVisited(false);}
					}
				}
				if(u.offsetValid() && !u.offsetContradiction()){
					min=Tools.min(min, u.offset());
				}
			}
			
			if(verbose){
				System.err.println("min offset = "+min);
			}
			
			for(Unit u : temp){
				if(u.offsetValid()){
					if(verbose){System.err.println("Set "+u.name()+" offset from "+u.offset+" to "+(u.offset-min));}
					u.offset=u.offset-min;
				}
			}
			
			
			return contradictions;
		}
		
		/**
		 * @param root
		 */
		private int setOffsetsNeighbors(final Unit root, final ArrayList<Unit> temp) {
			if(verbose){System.err.println("\nsetOffsetsNeighbors("+root.name()+")\nroot.code1="+root.code1+"\n");}
			assert(root.valid());
			assert(!root.visited());
			assert(root.offsetValid());
			assert(!root.offsetContradiction());
			root.setVisited(true);
			if(root.overlapList==null){return 0;}
			final int contradictions=countOffsetContradictions(root, false);
			if(verbose){System.err.println("\ncontradictions="+contradictions);}
			for(Overlap o : root.overlapList){
				Unit u=(o.u1==root ? o.u2 : o.u1);
				assert(o.u1==root || o.u2==root);
				assert(root!=u);
				assert(u.valid());
				
				if(verbose){System.err.println("\nProcessing Overlap "+o);}
				if(!o.visited() && !o.offsetContradiction()){
					o.setVisited(true);
					if(!u.offsetContradiction()){
						if(verbose){System.err.println("Calling setOffset:  "+o);}
						if(!u.offsetValid()){temp.add(u);}
						boolean b=setOffset(root, u, o);
						if(verbose){System.err.println("Finished setOffset: "+o);}
						
//						if(x>0){
//							if(verbose){System.err.println("\n*********************************************");}
//							if(verbose){System.err.println("Problem detected with contig "+u.name());}
//							if(verbose){System.err.println("*********************************************\n");}
//							verbose=true;
//							int y2=countOffsetContradictions(root, false);
//							assert(contradictions==y2);
//						}
						
						assert(b) : "\n"+contradictions+", "+o.offsetContradiction()+", "+root.offsetContradiction()+", "+u.offsetContradiction()+"\n"
							+root.offsetValid()+", "+u.offsetValid()+", "+OVERLAP_TYPE_NAMES[o.type]+"\n"+b
							+fixMultiJoins; //This assertion can fail if a multijoin is present
						assert(u.offsetValid());
					}
				}
			}
			return contradictions;
		}
		
		private int countOffsetContradictions(Unit root, boolean includeKnown){
			if(verbose){System.err.println("\ncountContradictions("+root.name()+", "+includeKnown+")\nroot.code1="+root.code1+"\n");}
			assert(root.valid());
			assert(root.visited());
			assert(root.offsetValid());
//			assert(!root.offsetContradiction());
			if(root.overlapList==null){return 0;}
			int contradictions=0;
			for(Overlap o : root.overlapList){
				Unit u=(o.u1==root ? o.u2 : o.u1);
				assert(o.u1==root || o.u2==root);
				assert(root!=u);
				assert(u.valid());

				if(verbose){System.err.println("\nOverlap "+o+"\nu="+u.name()+", offsetValid="+u.offsetValid());}
				
				boolean contradictory=(u.offsetValid() && u.offset()!=calcOffset(root, u, o));
				if(verbose){System.err.println("contradictory=            \t"+contradictory);}
				if(contradictory){
					if(includeKnown || !u.offsetContradiction()){
						contradictions++;
						if(!root.offsetContradiction()){root.setOffsetContradiction(true);}
					}
					if(!o.offsetContradiction()){o.setOffsetContradiction(true);}
					if(!u.offsetContradiction()){u.setOffsetContradiction(true);}
				}
				assert(contradictory==o.offsetContradiction()) : contradictory+", "+o.offsetContradiction();
				if(verbose){
					System.err.println("root.offsetContradiction()=\t"+root.offsetContradiction());
					System.err.println("u.offsetContradiction()=   \t"+u.offsetContradiction());
					System.err.println("o.offsetContradiction()=   \t"+o.offsetContradiction());
					System.err.println("contradictions=           \t"+contradictions);
				}
			}
			if(verbose){System.err.println("Final contradictions="+contradictions+"\n");}
			return contradictions;
		}
		
		private String toString(ArrayList<Unit> cluster){
			for(int i=0; i<cluster.size(); i++){
				Unit u=cluster.get(i);
				u.r.id=""+i;
			}
			StringBuilder sb=new StringBuilder(1000);
			for(Unit u : cluster){
				sb.append(">"+u.name()+"\n");
				sb.append(new String(u.bases()));
				sb.append("\n");
			}
			sb.append("\n*****\n");
			for(Unit u : cluster){
				sb.append("\n"+u.name()+":");
				if(u.overlapList!=null){
					for(Overlap o : u.overlapList){
						Unit ub=(o.u1==u ? o.u2 : o.u1);
						sb.append(" "+ub.name());
					}
				}
			}
			sb.append("\n");
			return sb.toString();
		}
		
		private String toShortString(ArrayList<Unit> cluster){
			for(int i=0; i<cluster.size(); i++){
				Unit u=cluster.get(i);
				u.r.id=""+i;
			}
			StringBuilder sb=new StringBuilder(1000);
			for(Unit u : cluster){
				sb.append("\n"+u.name()+":");
				if(u.overlapList!=null){
					for(Overlap o : u.overlapList){
						Unit ub=(o.u1==u ? o.u2 : o.u1);
						sb.append(" "+ub.name());
					}
				}
			}
			sb.append("\n");
			return sb.toString();
		}
		
		
		/**
		 * @param root
		 * @param u2
		 * @param o
		 * @return Number of contradictions
		 */
		private int canonicize(final Unit root, final Unit u2, final Overlap o){
			if(o.type==FORWARD){return 0;}
			if(o.type==REVERSE){o.reverseDirection();}
			assert(o.type==FORWARD);
			assert(o.test(bandy, o.edits+maxEdits));
			return 0;
		}
		
		
		/**
		 * @param root
		 * @param u2
		 * @param o
		 * @return true if no contradictions
		 */
		private boolean setOffset(final Unit root, final Unit u2, final Overlap o){
			assert(root.offsetValid());
			assert(!root.offsetContradiction());
			int offset=calcOffset(root, u2, o);
			
			if(u2.offsetValid()){return u2.offset()==offset;}
			u2.setOffset(offset);
			
			if(verbose){
				System.err.println("\nroot = "+(root.name()==null ? root.r.numericID+"" : root.name())+", u2 = "+(u2.name()==null ? u2.r.numericID+"" : u2.name())
						+"\no = "+o
						+"\nroot.offset = "+root.offset()
						+"\nu2.offset = "+u2.offset());
			}
			
			return true;
		}
		
		
		private int calcOffset(final Unit root, final Unit ub, final Overlap o){
			assert(root.offsetValid());
			if(o.type==FORWARD){
				if(root==o.u1){
					int dif=o.start1-o.start2;
					if(verbose){System.err.println("root==o.u1=="+root.name()+", start1="+o.start1+"; u2==o.u2=="+ub.name()+", start2="+o.start2+", dif="+dif);}
					return root.offset+dif;
				}else{
					int dif=o.start2-o.start1;
					if(verbose){System.err.println("root==o.u2=="+root.name()+", start2="+o.start2+"; u2==o.u1=="+ub.name()+", start1="+o.start1+", dif="+dif);}
					return root.offset+dif;
				}
			}else{
				assert(false) : o;
				throw new RuntimeException("TODO");
			}
		}
		
		
		/**
		 * @param cluster
		 */
		private void mergeLeaves(ArrayList<Unit> cluster) {
			assert(false) : "TODO";
			for(Unit u : cluster){
				
			}
		}
		
		/**
		 * @param cluster
		 */
		private void mergeInner(ArrayList<Unit> cluster) {
			assert(false) : "TODO";
			for(Unit u : cluster){
				
			}
		}

		private ArrayList<ArrayList<Unit>> processedClustersT=new ArrayList<ArrayList<Unit>>(threadMaxReadsToBuffer);

		long leafMergesT=0;
		long innerMergesT=0;
		long leafBaseMergesT=0;
		long innerBaseMergesT=0;
		
		long multiJoinFailuresT=0;
		long multiJoinsFoundT=0;
		long multiJoinBasesFoundT=0;
		long unitsFlippedT=0;
		long overlapsFlippedT=0;
		long offsetContradictoryOverlapsT=0;
		long offsetContradictoryClustersT=0;
		long cycleOverlapsT=0;
		long cycleClustersT=0;
		
		private final boolean fixMultiJoinsT;
		private final boolean fixOffsetContradictionsT;
		private final boolean mergeClustersT;
		private final boolean mergeLeavesT;
		private final boolean mergeInnerT;
		private final BandedAligner bandy;
	}

	
	/**
	 * @param cluster
	 */
	private void unvisit(ArrayList<Unit> cluster) {
		for(Unit u : cluster){
			if(u.visited()){u.setVisited(false);}
		}
	}
	
	/**
	 * @param cluster
	 */
	private void reorderClusterBreadthFirst(ArrayList<Unit> cluster) {
		if(verbose){System.err.println("reorderClusterBreadthFirst");}
		
		final int size=cluster.size();
		Shared.sort(cluster); //Now it is in descending length
		final Unit root=cluster.get(0);
		assert(root.length()>=cluster.get(size-1).length()) : root.length()+", "+cluster.get(size-1).length()+", "+root.compareTo(cluster.get(size-1));
		
		ArrayList<Unit> breadthFirst=new ArrayList<Unit>(cluster.size());
		root.setVisited(true);
//		System.err.println("root="+root.name());
		breadthFirst.add(root);
		for(int i=0; i<breadthFirst.size(); i++){
			Unit u=breadthFirst.get(i);
			Shared.sort(u.overlapList); //Sorted in descending overlap length
			if(u.overlapList!=null){
				for(Overlap o : u.overlapList){
					if(!o.u1.visited()){
						//					System.err.println("Visiting "+o.u1.name());
						o.u1.setVisited(true);
						breadthFirst.add(o.u1);
					}
					if(!o.u2.visited()){
						//					System.err.println("Visiting "+o.u2.name());
						o.u2.setVisited(true);
						breadthFirst.add(o.u2);
					}
					//				System.err.println("***");
					//				System.err.println(toShortString(breadthFirst));
				}
			}
		}
		for(Unit u : cluster){
			assert(u.visited());
			if(u.visited()){u.setVisited(false);}
			if(EA){
				if(u.overlapList!=null){
					for(Overlap o : u.overlapList){assert(!o.visited());}
				}
			}
		}
//		System.err.println("***");
//		System.err.println("Final:");
//		System.err.println(toShortString(breadthFirst));
		assert(cluster.size()==breadthFirst.size());
		cluster.clear();
		cluster.addAll(breadthFirst);
	}
	

	
	/** Returns next cluster larger than 1 element.
	 * Singleton clusters are added directly to 'processed'. */
	private ArrayList<Unit> nextCluster(){
		synchronized(clusterQueue){
			ArrayList<Unit> cluster=clusterQueue.poll();
			assert(cluster==null || cluster.size()>1);
//			while(cluster!=null && cluster.size()<2){
////				unmark(cluster);
//				processedClustersT.add(cluster);
//				cluster=clusterQueue.poll();
//			}
			return cluster;
		}
	}
	
	
	/**
	 * Creates Unit objects or uses ones already attached to reads.
	 * Places them in local storage and percolates them to shared storage (codeMap), removing exact duplicates.
	 * Also hashes tips and places these in shared affixMap.
	 * Looks for containments in the affix map.
	 * @author Brian Bushnell
	 * @date Jul 24, 2013
	 *
	 */
	private final class HashThread extends Thread{
		
		public HashThread(boolean addToCodeMap_, boolean addToAffixMap_, boolean findMatches_, boolean findContainments_, boolean findOverlaps_){
			addToCodeMapT=addToCodeMap_;
			addToAffixMapT=addToAffixMap_;
			findContainmentsT=findContainments_;
			findOverlapsT=findOverlaps_;
			findMatchesT=findMatches_;
			tid=getTid();
			crisq=new ArrayDeque<ConcurrentReadInputStream>(crisa.length);
			for(int i=0; i<crisa.length; i++){
//				if(verbose){System.err.println("Adding to crisq.");}
				crisq.add(crisa[(i+tid)%crisa.length]);
			}
			bandy=(maxEdits>0 && (findOverlapsT || findContainmentsT) ? BandedAligner.makeBandedAligner(bandwidth) : null);
			
//			assert(addToCodeMapT) : "addToCodeMapT="+addToCodeMapT+", addToAffixMapT="+addToAffixMapT+", findContainmentsT="+findContainmentsT+
//			", findOverlapsT="+findOverlapsT+", findMatchesT="+findMatchesT+", convertToUpperCaseT="+convertToUpperCaseT+", numAffixMaps="+numAffixMaps;
		}
		
		@Override
		public void run(){
			
			ConcurrentReadInputStream cris=crisq.poll();
			
//			String lastName=null;
			
			while(cris!=null){
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					
					for(Read r : reads){
//						lastName=r.id;
						processReadOuter(r);
					}
					
					if(codeMapT!=null && (codeMapT.size()>threadMaxReadsToBuffer || basesStoredT>threadMaxBasesToBuffer)){
						assert(addToCodeMapT);
						long added=mergeMaps();
						addedToMainT+=added;
					}
					
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				cris.returnList(ln);
				if(codeMapT!=null && !codeMapT.isEmpty()){
					long added=mergeMaps();
					addedToMainT+=added;
				}
				cris=crisq.poll();
			}
			
//			System.err.println("Thread finished; last name = "+lastName);
			
			codeMapT=null;
			localConflictList=null;
			sharedConflictList=null;
		}
		
		/** Return true if this read was a member of this subset. */
		private boolean processReadOuter(Read r1){
			if(r1.length()<MINSCAF){return false;}
			Read r2=r1.mate;
			
			assert(r1.pairnum()==0);
			assert(r2==null || r2.pairnum()==1);
			
			if(!addToCodeMapT && r1.obj==null){
				if(r1.bases!=null && r1.length()>=MINSCAF){
					final Unit u=(r1.obj!=null ? (Unit)r1.obj : new Unit(r1));
					assert(u.r==r1 && (r1.obj==u || r1.obj==null));
					final long code=u.code1;
					r1.obj=u;
					assert(u.r==r1 && r1.obj==u);
					if(r2!=null && r2.obj==null){r2.obj=new Unit(r2);}
					
					//Check for subset membership
					final boolean inSet=u.inSet();
					if(inSet){
						final Long codeL=code;
						ArrayList<Unit> list=codeMap.get(codeL);
						boolean found=false;
						for(Unit u0 : list){
							//Replace with existing read
							if(u0.equals(u) && u0.r.numericID==r1.numericID){
								r1=u0.r;
								r2=r1.mate;
								found=true;
								break;
							}
						}
						assert(list!=null);
						if(!found){
							return false;
						}
					}
				}
			}
			boolean b=processRead(r1);
			if(r2!=null){processRead(r2);}
			return b;
		}
		
		/** Return true if this read was a member of this subset. */
		private boolean processRead(Read r){
			if(r.length()<MINSCAF){return false;}

			final boolean inSet;
			if(!storeName){r.id=null;}
			if(!storeQuality){r.quality=null;}

			if(forceTrimLeft>0 || forceTrimRight>0){//Added at request of RQC team
				if(r!=null && r.length()>0){
					TrimRead.trimToPosition(r, forceTrimLeft>0 ? forceTrimLeft : 0, forceTrimRight>0 ? forceTrimRight : r.length(), 1);
				}
			}
			if(qTrimLeft || qTrimRight){
				TrimRead.trimFast(r, qTrimLeft, qTrimRight, trimQ, trimE, 0);
			}
			if(r.length()<MINSCAF){return false;}

			readsProcessedT++;
			basesProcessedT+=r.length();

			final Unit u=(r.obj!=null ? (Unit)r.obj : new Unit(r));
			assert(u.r==r && (r.obj==u || r.obj==null));
			final long code=u.code1;

			//Check for subset membership
			inSet=u.inSet();

			r.obj=u;
			assert(u.r==r && r.obj==u);
			if(r.mate!=null && r.mate.obj==null){r.mate.obj=new Unit(r.mate);}

			if(verbose){System.err.println("Generated "+code+" for sequence "+u.name()+"\t"+new String(r.bases, 0, Tools.min(40, r.length())));}

			if(addToCodeMapT && inSet){
				final Long codeL=code;
				ArrayList<Unit> list=codeMapT.get(codeL);
				if(list==null){
					if(verbose){System.err.println("Unique.");}
					list=new ArrayList<Unit>(1);
					list.add(u);
					basesStoredT+=r.length();
					codeMapT.put(codeL, list);
				}else{
					if(verbose){System.err.println("Exists.");}
					boolean match=false;
					if(findMatchesT){
						for(Unit u2 : list){
							if(pairedEqualsRC(u, u2)){
//								if(u.r.mate!=null){
//									verbose=true;
//
//									Unit um=(Unit)u.r.mate.obj;
//									Unit u2m=(Unit)u2.r.mate.obj;
//
//									if(verbose){
//										System.err.println("********");
//										System.err.println(u.r.toFastq());
//										System.err.println(u.r.mate.toFastq());
//										System.err.println("********");
//										System.err.println(u2.r.toFastq());
//										System.err.println(u2.r.mate.toFastq());
//										System.err.println("********");
//										System.err.println(u);
//										System.err.println(u2);
//										System.err.println(um);
//										System.err.println(u2m);
//										System.err.println("********");
//										System.err.println(u.equals(u2));
//										System.err.println(u.compareTo(u2));
//										System.err.println("********");
//										System.err.println(um.equals(u2m));
//										System.err.println(um.compareTo(u2m));
//										System.err.println("********");
//									}
//
//									verbose=false;
//								}
								assert(u.r.mate==null || pairedEqualsRC((Unit)u.r.mate.obj, (Unit)u2.r.mate.obj)) :
									u.r.toFastq()+"\n"+u2.r.toFastq()+"\n"+u.r.mate.toFastq()+"\n"+u2.r.mate.toFastq()+
									"\n"+u+"\n"+u2+"\n"+u.r.mate.obj+"\n"+u2.r.mate.obj;
//								if(verbose){System.err.println("Matches "+new String(r2.bases, 0, Tools.min(40, r2.length())));}
								match=true;
								u2.absorbMatch(u);
								if(UNIQUE_ONLY){
									synchronized(u2){
										if(u2.valid()){
											matchesT++;
											baseMatchesT+=u2.length();
											u2.setValid(false);
											addDupe(u2.r);
										}
									}
								}
								break;
							}
						}
					}
					if(match){
						addDupe(r);
						matchesT++;
						baseMatchesT+=r.length();
						//							if(verbose){System.err.println("matchesT="+matchesT+", baseMatchesT="+baseMatchesT);}
					}else{
						collisionsT++;
						if(verbose){System.err.println("False collision; count = "+collisionsT);}
						list.add(u);
						basesStoredT+=r.length();
					}
				}
			}

			if(findContainmentsT){
				int x=findContainments(u);
			}

			if(findOverlapsT){
				int x=findOverlaps(u);
			}

			return inSet;
		}
		
		private int findContainments(final Unit u){
			if(minLengthPercent<=0 && maxSubs<=0 && minIdentity>=100 && !u.valid()){return 0;}
			final byte[] bases=u.bases();
			final int shift=bitsPerBase*k;
			final int shift2=shift-bitsPerBase;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			long kmer=0;
			int hits=0;
			int currentContainments=0;
			int len=0;
			
			if(bases==null || bases.length<k){return -1;}
			final LongM key=new LongM();
			
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=acidToNumber[b];
				kmer=((kmer<<bitsPerBase)|x)&mask;
				if(HASH_NS || isFullyDefined(b)){len++;}
				else{len=0;}
//				if(verbose){System.err.println("Scanning i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					key.set(kmer); //Canonical
					for(int am=0; am<affixMaps.length; am++){
						ArrayList<Unit> list=affixMaps[am].get(key);
						if(list!=null){
							for(Unit u2 : list){
								if(u!=u2 && !u.equals(u2)){
									if(u2.valid()){
										hits++;
										if(verbose){
											System.err.println("\nFound potential containment at am="+am+", i="+i+", key="+key.value()+
													", pre1="+u2.prefix1+", pre2="+u2.prefix2+
													", suf1="+u2.suffix1+", suf2="+u2.suffix2+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i, k)));
										}
										if(u.contains(u2, i, key, bandy, am)){
											synchronized(u2){
												if(u2.valid()){
													currentContainments++;
													baseContainmentsT+=u2.length();
													u2.setValid(false);
													addDupe(u2.r);
												}
											}
											if(UNIQUE_ONLY){
												synchronized(u){
													if(u.valid()){
														currentContainments++;
														baseContainmentsT+=u.length();
														u.setValid(false);
														addDupe(u.r);
													}
												}
											}

											if(verbose){System.err.println("Added containment "+u2);}
										}
									}
								}
							}
						}
					}
				}
			}
//			assert(false) : hits+", "+currentContainments+", "+baseContainments+"\n"+containmentMapT+"\n";
			
			containmentCollisionsT+=(hits-currentContainments);
//			outstream.println("hits="+hits+", currentContainments="+currentContainments);
			containmentsT+=currentContainments;
			return hits;
		}
		
		private int findOverlaps(final Unit u){
//			if(minLengthPercent<=0 && maxSubs<=0 && minIdentity>=100 && !u.valid()){return 0;}
//			if(u.overlapList!=null){u.overlapList.clear();}
			final byte[] bases=u.bases();
			final int shift=bitsPerBase*k;
			final int shift2=shift-bitsPerBase;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			long kmer=0;
			int hits=0;
			int currentOverlaps=0;
			int len=0;
			
			if(bases==null || bases.length<k){return -1;}
			final LongM key=new LongM();
			
			boolean quit=false;
			
			for(int i=0; i<bases.length && !quit; i++){
				byte b=bases[i];
				long x=acidToNumber[b];
				kmer=((kmer<<bitsPerBase)|x)&mask;
				if(HASH_NS || isFullyDefined(b)){len++;}
				else{len=0;}
//				if(verbose){System.err.println("Scanning i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					key.set(kmer); //Canonical key
					for(int am=0; am<affixMaps.length; am++){
						ArrayList<Unit> list=affixMaps[am].get(key);
						if(list!=null){//found a key collision
							for(Unit u2 : list){
								if(quit){break;}//too many edges
								int u1cluster=-1, u2cluster=-2;
								if(preventTransitiveOverlaps && u!=u2){
									u1cluster=u.determineCluster();
									u2cluster=u2.determineCluster();
								}
								if(u1cluster!=u2cluster && u!=u2 && !u.equals(u2) && u2.r!=u.r.mate){//TODO:  Not sure why identical things are banned...  possibly relates to avoiding inter-pair edges?
									if(u2.valid()){
										hits++;
										
//										boolean flag=(u.code1==-3676200394282040623L && u2.code1==-7034423913727372751L) ||
//												(u2.code1==-3676200394282040623L && u.code1==-7034423913727372751L);
										final boolean flag=false;
										if(verbose || flag){
											System.err.println("\nFound potential overlap at am="+am+", i="+i+", key="+key.value()+
													", pre1="+u2.prefix1+", pre2="+u2.prefix2+
													", suf1="+u2.suffix1+", suf2="+u2.suffix2+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i, k)));
										}
										
										final Overlap o;
										if(maxEdges>1000000000 || u.overlapList==null || u2.overlapList==null ||
												(u.overlapList.size()<maxEdges && u2.overlapList.size()<maxEdges2)){
											 o=u.makeOverlap(u2, i, key, bandy, am);

										}else{
											o=null;
											if(u.overlapList.size()>maxEdges){quit=true;}
										}
										if(o!=null){
											
											if(preventTransitiveOverlaps){
												 mergeClusterIds(u1cluster, u2cluster);
											}
											
											assert(o.test(bandy, o.edits+maxEdits)) : o;
											if(verbose || flag){System.err.println("Created overlap "+o);}

											long comp=u.length()-u2.length();
											if(comp==0){comp=u.code1-u2.code1;}
											if(comp==0){comp=u.code2-u2.code2;}
											if(comp==0){comp=u.prefix1-u2.prefix1;}
											if(comp==0){comp=u.suffix1-u2.suffix1;}
											if(comp==0){comp=(u.r.numericID-u2.r.numericID);}
											assert(comp!=0) : u+", "+u2;
											Unit ua=(comp<0 ? u : u2);
											Unit ub=(comp<0 ? u2 : u);
											assert(ua!=ub);
											if(verbose || flag){
												System.err.println("ua="+ua.code1);
												System.err.println("ub="+ub.code1);
												System.err.println("u ="+u.code1);
												System.err.println("u2="+u2.code1);
												System.err.println("u.r ="+u.r.numericID);
												System.err.println("u2.r="+u2.r.numericID);
												System.err.println("ua contains o? "+ua.alreadyHas(o));
												System.err.println("ub contains o? "+ub.alreadyHas(o));
												System.err.println("ua.list="+ua.overlapList);
												System.err.println("ub.list="+ub.overlapList);
											}
											
//											assert(ua.alreadyHas(o)==ub.alreadyHas(o));
											
											final boolean uaContainedOverlap;
											
											synchronized(ua){
												if(ua.overlapList==null){ua.overlapList=new ArrayList<Overlap>(2);}
												if(!ua.overlapList.contains(o)){
													if(EA){
														synchronized(ub){
															assert(ub.overlapList==null || !ub.overlapList.contains(o)) :
																ua.alreadyHas(o)+", "+ub.alreadyHas(o)+"\n"+o+"\n"+ub.overlapList.get(ub.overlapList.indexOf(o))+
																"\nua.list="+ua.overlapList+"\nub.list="+ub.overlapList+"\nu.code1="+u.code1+"\nu2.code1="+u2.code1;
														}
													}
													currentOverlaps++;
													baseOverlapsT+=o.overlapLen;
													ua.overlapList.add(o);
													if(verbose || flag){System.err.println("Added overlap "+o);}
													uaContainedOverlap=false;
												}else{
													if(verbose || flag){System.err.println("Already contained overlap "+o);}
													hits--;
													uaContainedOverlap=true;
												}
											}
											
											if(!uaContainedOverlap){
												synchronized(ub){
													if(ub.overlapList==null){ub.overlapList=new ArrayList<Overlap>(2);}
													assert(!ub.overlapList.contains(o));
													ub.overlapList.add(o);
													if(verbose || flag){System.err.println("Added overlap "+o);}
												}
											}else{
												if(verbose || flag){System.err.println("Already contained overlap "+o);}
											}
											

//											assert(ua.alreadyHas(o));
//											assert(ub.alreadyHas(o));
//											assert(ua.overlapList.contains(o));
//											assert(ub.overlapList.contains(o));
											if(verbose || flag){
												System.err.println("ua contains o? "+ua.alreadyHas(o));
												System.err.println("ub contains o? "+ub.alreadyHas(o));
												System.err.println("ua.list="+ua.overlapList);
												System.err.println("ub.list="+ub.overlapList);
											}
										}
									}
								}
							}
						}
					}
				}
			}
			if(EA){
				synchronized(u){
					if(u.overlapList!=null && u.overlapList.isEmpty()){
						assert(false) : "Why would this happen?";
						u.overlapList=null;
					}
				}
			}
//			assert(false) : hits+", "+currentOverlaps+", "+baseOverlaps+"\n"+overlapMapT+"\n";
			
//			assert(hits==currentOverlaps) : hits+", "+currentOverlaps;
			
			overlapCollisionsT+=(hits-currentOverlaps);
//			outstream.println("hits="+hits+", currentOverlaps="+currentOverlaps);
			overlapsT+=currentOverlaps;
			return hits;
		}
		
		/** Insert reads processed by a thread into the shared code and affix maps.
		 * If operating in subset mode, only store reads with code equal to subset mod subsetCount. */
		private long mergeMaps(){
			if(verbose){System.err.println("Merging maps.");}
			long novelReads=0, novelKeys=0;
			long collisionReads=0;
			long mergedReads=0;

			assert(localConflictList.isEmpty());
			assert(sharedConflictList.isEmpty());
			
			synchronized(codeMap){
				for(Long key : codeMapT.keySet()){
					if(codeMap.containsKey(key)){
						localConflictList.add(codeMapT.get(key));
						sharedConflictList.add(codeMap.get(key));
					}else{
						ArrayList<Unit> list=codeMapT.get(key);
						codeMap.put(key, list);
						addedList.addAll(list);
						novelReads+=list.size();
						novelKeys++;
					}
				}
			}
			
			if(verbose){System.err.println("Novel reads = "+novelReads+", conflicts = "+localConflictList.size());}
			
			for(int i=0; i<localConflictList.size(); i++){
				ArrayList<Unit> listT=localConflictList.get(i);
				ArrayList<Unit> list=sharedConflictList.get(i);
				synchronized(list){
					for(Unit u : listT){
						if(verbose){System.err.println("Processing novel unit "+u.name());}
						boolean match=false;
						if(findMatchesT){
							for(Unit u2 : list){
								if(pairedEqualsRC(u, u2)){
									//								if(verbose){System.err.println("Matches "+new String(r2.bases, 0, Tools.min(40, r2.length())));}
									u2.absorbMatch(u);
									if(UNIQUE_ONLY){
										synchronized(u2){
											if(u2.valid()){
												mergedReads++;
												baseMatchesT+=u2.length();
												u2.setValid(false);
												addDupe(u2.r);
											}
										}
									}
									match=true;
									break;
								}
							}
						}
						if(match){
							addDupe(u.r);
							mergedReads++;
							baseMatchesT+=u.length();
							if(verbose){System.err.println("matchesT="+matchesT+", baseMatchesT="+baseMatchesT);}
						}else{
							collisionReads++;
							if(verbose){System.err.println("False collision; count = "+collisionReads);}
							list.add(u);
							addedList.add(u);
						}
					}
				}
			}
			matchesT+=mergedReads;
			collisionsT+=collisionReads;
			if(verbose){System.err.println("Done Merging.");}
			if(verbose){System.err.println("mapT.size="+codeMapT.size()+", basesStoredT="+basesStoredT);}
			
			codeMapT.clear();
			localConflictList.clear();
			sharedConflictList.clear();
			
			if(!addedList.isEmpty()){
				if(addToAffixMapT){
					final LongM p=new LongM(-1, true);
					assert(affixMap1!=null || affixMap2!=null);
					if(affixMap1!=null && !ignoreAffix1){//Allows you to not use am1
						synchronized(affixMap1){
							for(Unit u : addedList){
								if(verbose){System.err.println("Processing affixes for "+u.name());}
								if(u.prefix1!=-1 || u.prefix1!=u.suffix1){
									if(verbose){System.err.println("Using prefix "+u.prefix1);}
									p.set(u.prefix1);
									ArrayList<Unit> alu=affixMap1.get(p);
									if(alu==null){
										if(verbose){System.err.println("Made new alu for "+p);}
										alu=new ArrayList<Unit>(2);
										affixMap1.put(p.iCopy(), alu);
									}
									if(alu.size()<maxAffixCopies){
										if(verbose){System.err.println("Added "+u.name());}
										alu.add(u);
									}
									if(verbose){System.err.println(affixMap1.get(p));}
								}
								if(storeSuffix && u.prefix1!=u.suffix1){
									if(verbose){System.err.println("Using suffix "+u.suffix1);}
									p.set(u.suffix1);
									ArrayList<Unit> alu=affixMap1.get(p);
									if(alu==null){
										if(verbose){System.err.println("Made new alu for "+p);}
										alu=new ArrayList<Unit>(2);
										affixMap1.put(p.iCopy(), alu);
									}
									if(alu.size()<maxAffixCopies){
										if(verbose){System.err.println("Added "+u.name());}
										alu.add(u);
									}
									if(verbose){System.err.println(affixMap1.get(p));}
								}
							}
						}
					}
					if(affixMap2!=null){
						synchronized(affixMap2){
							for(Unit u : addedList){
								if(u.prefix2!=-1 || u.prefix2!=u.suffix2){
									p.set(u.prefix2);
									ArrayList<Unit> alu=affixMap2.get(p);
									if(alu==null){
										alu=new ArrayList<Unit>(2);
										affixMap2.put(p.iCopy(), alu);
									}
									if(alu.size()<maxAffixCopies){alu.add(u);}
								}
								if(storeSuffix && u.prefix2!=u.suffix2){
									p.set(u.suffix2);
									ArrayList<Unit> alu=affixMap2.get(p);
									if(alu==null){
										alu=new ArrayList<Unit>(2);
										affixMap2.put(p.iCopy(), alu);
									}
									if(alu.size()<maxAffixCopies){alu.add(u);}
								}
							}
						}
					}
				}
			}
			
			addedList.clear();
			basesStoredT=0;
			return collisionReads+novelReads;
		}
		
		private int getTid(){
			synchronized(HashThread.class){
				int x=tcount;
				tcount++;
				return x;
			}
		}
		
		private LinkedHashMap<Long, ArrayList<Unit>> codeMapT=new LinkedHashMap<Long, ArrayList<Unit>>(threadMaxReadsToBuffer*8);
		private ArrayList<Unit> addedList=new ArrayList<Unit>(threadMaxReadsToBuffer);
		private ArrayList<ArrayList<Unit>> localConflictList=new ArrayList<ArrayList<Unit>>(threadMaxReadsToBuffer);
		private ArrayList<ArrayList<Unit>> sharedConflictList=new ArrayList<ArrayList<Unit>>(threadMaxReadsToBuffer);
		
		long matchesT=0;
		long baseMatchesT=0;
		long baseContainmentsT=0;
		long collisionsT=0;
		long containmentsT=0;
		long containmentCollisionsT=0;
		long basesStoredT=0;
		long addedToMainT=0;
		long readsProcessedT=0;
		long basesProcessedT=0;
		long overlapsT=0;
		long baseOverlapsT=0;
		long overlapCollisionsT=0;
		
		private final boolean addToCodeMapT;
		private final boolean addToAffixMapT;
		private final boolean findContainmentsT;
		private final boolean findOverlapsT;
		private final boolean findMatchesT;
//		private final boolean convertToUpperCaseT;
		private final int tid;
		private final ArrayDeque<ConcurrentReadInputStream> crisq;
		private final BandedAligner bandy;
	}
	
	public static boolean equalsRC(byte[] a, byte[] b){
		if(a==b){return true;}
		if(a==null || b==null){return false;}
		if(a.length!=b.length){return false;}
		
		for(int i=0; i<a.length; i++){
			final byte aa=a[i], bb=b[i];
			if(aa!=bb){return false;}
		}
		return true;
	}
	
	public static boolean pairedEqualsRC(Unit ua, Unit ub){
		if(verbose){System.err.println("pairedEqualsRC("+ua.name()+", "+ub.name()+")");}
		if(verbose){System.err.println("ea");}
		boolean b=equalsRC(ua, ub);
		if(verbose){System.err.println("eb");}
		if(!b){return false;}
		if(verbose){System.err.println("ec");}
		
		if(ua.r!=null && ub.r!=null){
			if(verbose){System.err.println("ed");}
			assert((ua.r.mate==null)==(ub.r.mate==null));
			if(verbose){System.err.println("ee");}
			if(ua.r.mate!=null && ub.r.mate!=null){
				if(verbose){System.err.println("ef");}
				return ua.r.pairnum()==ub.r.pairnum() && Tools.compare(ua.r.mate.bases, ub.r.mate.bases)==0;
			}
			if(verbose){System.err.println("eg");}
		}
		if(verbose){System.err.println("eh");}
		return true;
	}
	
	private static boolean equalsRC(Unit ua, Unit ub){
		if(verbose){System.err.println("equalsRC("+ua.name()+", "+ub.name()+")");}
		return ua.code1==ub.code1 && ua.code2==ub.code2 && (ua.prefix1==ub.prefix1 && ua.suffix1==ub.suffix1) && compareRC(ua, ub)==0;
	}
	
	public static int comparePairedRC(Unit ua, Unit ub){
		if(verbose){System.err.println("comparePairedRC("+ua.name()+", "+ub.name()+")");}
		int x=compareRC(ua, ub);
		if(x!=0){return x;}
		
		if(ua.r!=null && ub.r!=null && ua.r.mate!=null && ub.r.mate!=null){
			if(ua.r.pairnum()!=ub.r.pairnum()){return ua.r.pairnum()-ub.r.pairnum();}
			return compareRC((Unit)ua.r.mate.obj, (Unit)ub.r.mate.obj);
		}
		return 0;
	}
	
	//TODO
	//This is really for sorting by length.
	private static int compareRC(Unit ua, Unit ub){
		if(verbose){System.err.println("compareRC("+ua.name()+", "+ub.name()+")");}
		if(ua==ub){return 0;}
		if(verbose){System.err.println("a");}
		if(verbose){System.err.println("a1");}
		if(ua.length()!=ub.length()){return ub.length()-ua.length();}
		if(verbose){System.err.println("a2");}
		
		if(REQUIRE_MATCHING_NAMES){
			if(ua.name()!=null && ub.name()!=null){
				int x=ua.name().compareTo(ub.name());
				if(x!=0){return x;}
			}
		}
		if(verbose){System.err.println("a3");}
		
		if(ua.r==null || ub.r==null){
			if(verbose){System.err.println("b");}
			if(verbose){System.err.println("b1");}
			if(verbose){System.err.println("c");}
			if(ua.prefix1!=ub.prefix1){return ua.prefix1>ub.prefix1 ? 1 : -1;}
			if(ua.suffix1!=ub.suffix1){return ua.suffix1>ub.suffix1 ? 1 : -1;}
			if(verbose){System.err.println("e");}
			if(ua.code1!=ub.code1){return ua.code1>ub.code1 ? 1 : -1;}
			if(ua.code2!=ub.code2){return ua.code2>ub.code2 ? 1 : -1;}
			
			return ua.pairnum()-ub.pairnum();
		}
		if(verbose){System.err.println("f");}
		final byte[] a=ua.r.bases, b=ub.r.bases;
		if(a==b){return 0;}
		if(a==null || b==null){return a==null ? -1 : 1;}
		if(verbose){System.err.println("g");}
		
		if(verbose){System.err.println("h");}
		for(int i=0; i<a.length; i++){
			final byte aa=a[i], bb=b[i];
			if(aa!=bb){return aa-bb;}
		}
		
		if(verbose){System.err.println("j");}
		return ua.pairnum()-ub.pairnum();
	}
	
	private long hashTip(byte[] bases, boolean prefix, int k, int skipInitialBases){
		if(bases==null || bases.length<k){return -1;}

		final int shift=bitsPerBase*k;
		final int shift2=shift-bitsPerBase;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		
		final int start=(prefix ? 0+skipInitialBases : bases.length-k-skipInitialBases);
		final int stop=start+k;

//		if(verbose){
//			System.err.println("\n"+new String(bases));
//			System.err.println("prefix="+prefix+", start="+start+", stop="+stop);
////			System.err.print(new String(bases));
//		}
		for(int i=start; i<stop; i++){
			byte b=bases[i];
//			if(verbose){System.err.print((char)b);}
			long x=acidToNumber[b];
			kmer=((kmer<<bitsPerBase)|x)&mask;
			len++;
		}
		if(verbose){System.err.println(new String(bases, start, k)+" = "+kmer);}
		assert(len==k) : len+","+k;
		return kmer;
	}
	
	private static final int calcMaxEdits(int maxEdits, float minIdentityMult, int len){
		return minIdentityMult==0 ? maxEdits : Tools.max(maxEdits, (int)Math.round(len*minIdentityMult));
	}
	
	
	private class Overlap implements Comparable<Overlap>{
		
		public Overlap(Unit u1_, Unit u2_, int type_, int start1_, int start2_, int stop1_, int stop2_, int len_, int mismatches_, int edits_, BandedAligner bandy){
			assert(u1_!=u2_);
			if(verbose){System.err.println("\nCreating an overlap.");}
			u1=u1_;
			u2=u2_;
			type=type_;
			start1=start1_;
			start2=start2_;
			stop1=stop1_;
			stop2=stop2_;
			overlapLen=len_;
			mismatches=mismatches_;
			edits=edits_;

			assert(Tools.absdif(Tools.absdif(start1, stop1), Tools.absdif(start2, stop2))<=maxEdits) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)
				+"\n>1a\n"+new String(u1.r.bases, Tools.min(start1, stop1), Tools.max(start1, stop1)-Tools.min(start1, stop1)+1)
				+"\n>2a\n"+new String(u2.r.bases, Tools.min(start2, stop2), Tools.max(start2, stop2)-Tools.min(start2, stop2)+1);
			
			assert(start1>=0 && start1<=u1.length()) : "type "+type+": "+start1+", "+stop1+", "+u1.length()+", "+start2+", "+stop2+", "+u2.length();
			assert(stop1>=0 && stop1<=u1.length()) :   "type "+type+": "+start1+", "+stop1+", "+u1.length()+", "+start2+", "+stop2+", "+u2.length();
			assert(start2>=0 && start2<=u2.length()) : "type "+type+": "+start1+", "+stop1+", "+u1.length()+", "+start2+", "+stop2+", "+u2.length();
			assert(stop2>=0 && stop2<=u2.length()) :   "type "+type+": "+start1+", "+stop1+", "+u1.length()+", "+start2+", "+stop2+", "+u2.length();
			
			assert(type==FORWARD || type==REVERSE);
			
			if(verbose){System.err.println(this);}
			
			assert(Tools.absdif(Tools.absdif(start1, stop1), Tools.absdif(start1, stop1))<=maxEdits);
			
			assert(test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)
				+"\n>1a\n"+new String(u1.r.bases, Tools.min(start1, stop1), Tools.max(start1, stop1)-Tools.min(start1, stop1)+1)
				+"\n>2a\n"+new String(u2.r.bases, Tools.min(start2, stop2), Tools.max(start2, stop2)-Tools.min(start2, stop2)+1);
			if(verbose){System.err.println("Passed test 1.");}
			
//			bandy.verbose=true;
//			test(bandy);
//			assert(false);
			
			assert(u1!=u2);
			u1.firstInOverlap(u2);
			u2.firstInOverlap(u1);
			assert(u1.length()!=u2.length() || u1.code1!=u2.code1 || u1.code2!=u2.code2 || (u1.r!=null && u1.r.mate!=null)) : "Collision? \n"+this+"\n"+u1+"\n"+u2;
			assert(u1.firstInOverlap(u2)!=u2.firstInOverlap(u1)) :
				"\nu1.firstInOverlap(u2)="+u1.firstInOverlap(u2)+"\nu2.firstInOverlap(u1)="+u2.firstInOverlap(u1)+"\nu1="+u1+"\nu2="+u2;
			
			if(!u1.firstInOverlap(u2)){
				if(verbose){System.err.println("\nSwapping.");}
				swap();
				if(verbose){System.err.println(this);}
				
				if(EA && !customBandwidth && !test(bandy, edits+maxEdits)){
					System.err.println("\n"+this);
					swap();
					System.err.println("\n"+this);
					assert(test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)+"\n";
					System.err.println("Passed test 2a, "+bandy.lastEdits+" edits.\n");
					swap();
					System.err.println("\n"+this);
					assert(test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)
						+"\n>1a\n"+new String(u1.r.bases, Tools.min(start1, stop1), Tools.max(start1, stop1)-Tools.min(start1, stop1)+1)
						+"\n>2a\n"+new String(u2.r.bases, Tools.min(start2, stop2), Tools.max(start2, stop2)-Tools.min(start2, stop2)+1);
					System.err.println("Passed test 2b, "+bandy.lastEdits+" edits.\n");
				}
				
				assert(customBandwidth || test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)
					+"\n>1a\n"+new String(u1.r.bases, Tools.min(start1, stop1), Tools.max(start1, stop1)-Tools.min(start1, stop1)+1)
					+"\n>2a\n"+new String(u2.r.bases, Tools.min(start2, stop2), Tools.max(start2, stop2)-Tools.min(start2, stop2)+1);
				if(verbose){System.err.println("Passed test 2.");}
			}
			
			if(type==REVERSE){
				if(verbose){System.err.println("\nReversing.");}
				reverseDirection();
				if(verbose){System.err.println(this);}
				
				if(EA && !Shared.anomaly && !customBandwidth && bandy!=null && !test(bandy, edits+maxEdits)){
					Shared.anomaly=true;
					BandedAligner.verbose=true;
					System.err.println("\n********** Failed test 3, "+bandy.lastEdits+" edits. **********\n");
					reverseDirection();
					System.err.println(this);
					assert(test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)+"\n";
					System.err.println("Passed test 3a, "+bandy.lastEdits+" edits.\n");
					reverseDirection();
					System.err.println(this);
					assert(test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)
						+"\n>1a\n"+new String(u1.r.bases, Tools.min(start1, stop1), Tools.max(start1, stop1)-Tools.min(start1, stop1)+1)
						+"\n>2a\n"+new String(u2.r.bases, Tools.min(start2, stop2), Tools.max(start2, stop2)-Tools.min(start2, stop2)+1);
					System.err.println("Passed test 3b, "+bandy.lastEdits+" edits.\n");
					BandedAligner.verbose=false;
					assert(false);
				}
				
				assert(customBandwidth || test(bandy, edits+maxEdits)) : "\n"+this+"\n>1\n"+new String(u1.r.bases)+"\n>2\n"+new String(u2.r.bases)+"\n";
				if(verbose){System.err.println("Passed test 3.");}
			}
			//Now all overlaps should be FORWARD or FORWARDRC and u1 should be at least as big as u2
			assert(type==FORWARD);
			assert(u1.length()>=u2.length());
			assert(u1.firstInOverlap(u2));
			assert(!u2.firstInOverlap(u1));
			if(verbose){System.err.println("Finished overlap initialization.");}
		}

		public boolean test(BandedAligner bandy, int editLimit){
			final int last1=u1.length()-1, last2=u2.length()-1;
			if(verbose){System.err.println("Testing "+OVERLAP_TYPE_NAMES[type]+", "+start1+", "+start2);}
			if(type==FORWARD){
				assert(start1==0 || start2==0) : "start1="+start1+", stop1="+stop1+", last1="+last1+", start2="+start2+", stop2="+stop2+", last2="+last2;
				if(start2==0){
					if(verbose){System.err.println("A");}
					return u1.overlapsForward(u2, start1, start2, bandy, false, editLimit);}
				else{
					if(verbose){System.err.println("B");}
					return u2.overlapsForward(u1, start2, start1, bandy, false, editLimit);}
			}
			if(type==REVERSE){
				assert(start1==last1 || start2==last2) : "start1="+start1+", stop1="+stop1+", last1="+last1+", start2="+start2+", stop2="+stop2+", last2="+last2;
				if(start2==last2){return u1.overlapsReverse(u2, start1, start2, bandy, false, editLimit);}
				else{return u2.overlapsReverse(u1, start2, start1, bandy, false, editLimit);}
			}
			throw new RuntimeException();
		}
		
		@Override
		public boolean equals(Object o){
			return equals((Overlap)o);
		}
		
		public boolean equals(Overlap o){
			if(this==o){return true;}
			assert(o!=null) : "*A*\n"+this+"\n"+o+"\n"+u1+"\n"+u2;
			assert(u1!=null && u2!=null) : "*B*\n"+this+"\n"+o+"\n"+u1+"\n"+u2;
			assert(u1!=o.u2 || u2!=o.u1) : "*C*\n"+this+"\n"+o+"\n"+u1.firstInOverlap(u2)+"\n"+o.u1.firstInOverlap(o.u2)+"\n"+u1+"\n"+u2;
			return (u1==o.u1 && u2==o.u2 && type==o.type && start1==o.start1 && start2==o.start2 && stop1==o.stop1 && stop2==o.stop2)
					;//|| (u1==o.u2 && u2==o.u1 && type==reverseType(o.type) && start1==o.start2 && start2==o.start1);
		}
		
//		public int compareTo(Overlap o){
//			int a=compareTo2(o);
//			int b=o.compareTo2(this);
//			assert(a==-b) : "\n"+this+"\n"+o+"\na="+a+", b="+b+", equals="+this.equals(o)
//				+"\nu1.compareTo(o.u1)="+u1.compareTo(o.u1)+"\no.u1.compareTo(u1)="+o.u1.compareTo(u1)
//				+"\nu2.compareTo(o.u2)="+u2.compareTo(o.u2)+"\no.u2.compareTo(u2)="+o.u2.compareTo(u2);
//			return a;
//		}
		
		@Override
		public int compareTo(Overlap o){
			int score1=overlapLen-50*(mismatches+edits);
			int score2=o.overlapLen-50*(o.mismatches+o.edits);
			if(score1!=score2){return score2-score1;}
			if(overlapLen!=o.overlapLen){return o.overlapLen-overlapLen;}
			int x=u1.compareTo(o.u1);
			if(x!=0){return -x;}
			x=u2.compareTo(o.u2);
			if(x!=0){return -x;}
			if(type!=o.type){return type-o.type;}
			if((u1!=o.u1 || u2!=o.u2) && absorbMatch && !subsetMode){
				boolean oldv=verbose;
				verbose=true;
				System.err.println(this);
				System.err.println(o);
				System.err.println("********");
				System.err.println(u1);
				System.err.println(u2);
				System.err.println(o.u1);
				System.err.println(o.u2);
				System.err.println("********");
				System.err.println(u1.equals(o.u1));
				System.err.println("********");
				System.err.println(u2.equals(o.u2));
				System.err.println("********");
				System.err.println(u1.compareTo(o.u1));
				System.err.println("********");
				System.err.println(u2.compareTo(o.u2));
				System.err.println("********");
				verbose=oldv;
			}
			assert(!absorbMatch || (u1==o.u1 && u2==o.u2) || subsetMode) : "\n"+u1.r+"\n"+u2.r+"\n"+o.u1.r+"\n"+o.u2.r
				+"\n\n"+u1.r.mate+"\n"+u2.r.mate+"\n"+o.u1.r.mate+"\n"+o.u2.r.mate;
//			assert(false) : "\n"+this+"\n"+o+"\n>"+u1.name()+"\n"+new String(u1.bases())+"\n>"+u2.name()+"\n"+new String(u2.bases())+"\n";
			if(start1!=o.start1){return start1-o.start1;}
			if(stop1!=o.stop1){return stop1-o.stop1;}
			if(start2!=o.start2){return start2-o.start2;}
			if(stop2!=o.stop2){return stop2-o.stop2;}
			if(this.equals(o)){
				return 0;
			}else{
				//TODO: ensure this assumption is valid.
				assert(!absorbContainment || !absorbMatch || subsetMode) : "\n"+this+"\n"+o+"\n>"+u1.name()+"\n"+new String(u1.bases())+"\n>"+u2.name()+"\n"+new String(u2.bases())+"\n";
				
				if(u1.unitID!=o.u1.unitID){return u1.unitID-o.u1.unitID;}
				if(u2.unitID!=o.u2.unitID){return u2.unitID-o.u2.unitID;}
			}
			return 0;
		}
		
		@Override
		public int hashCode(){
			return u1.hashCode()^u2.hashCode()^overlapLen;
		}
		
		public void swap(){
			Unit tempu=u1;
			u1=u2;
			u2=tempu;
			int temp=start1;
			start1=start2;
			start2=temp;
			temp=stop1;
			stop1=stop2;
			stop2=temp;
//			if(type==FORWARDRC){type=REVERSERC;}
//			else if(type==REVERSERC){type=FORWARDRC;}
			//What about forward and reverse?
		}
		
		public void reverseDirection(){
			type=reverseType(type);
			int temp=start1;
			start1=stop1;
			stop1=temp;
			temp=start2;
			start2=stop2;
			stop2=temp;
		}
		
		@Override
		public String toString(){
			StringBuilder sb=new StringBuilder(80);
			sb.append("type=");
			sb.append(OVERLAP_TYPE_NAMES[type]);
			sb.append(", len=");
			sb.append(overlapLen);
			sb.append(", subs=");
			sb.append(mismatches);
			sb.append(", edits=");
			sb.append(edits);
			
			sb.append(" (");
			sb.append(u1.name()==null ? u1.r.numericID+"" : u1.name());
			if(printLengthInEdges){
				sb.append(", length=");
				sb.append(u1.length());
			}
			sb.append(", start1=");
			sb.append(start1);
			sb.append(", stop1=");
			sb.append(stop1);
			
			sb.append(") (");
			sb.append(u2.name()==null ? u2.r.numericID+"" : u2.name());
			if(printLengthInEdges){
				sb.append(", length=");
				sb.append(u2.length());
			}
			sb.append(", start2=");
			sb.append(start2);
			sb.append(", stop2=");
			sb.append(stop2);
			sb.append(")");
			return sb.toString();
		}
		
		public String toLabel(){
			StringBuilder sb=new StringBuilder(80);
			sb.append(OVERLAP_TYPE_ABBREVIATIONS[type]);
			sb.append(',');
			sb.append(overlapLen);
			sb.append(',');
			sb.append(mismatches);
			sb.append(',');
			sb.append(edits);
			
			if(printLengthInEdges){
				sb.append(',');
				sb.append(u1.length());
			}
			sb.append(',');
			sb.append(start1);
			sb.append(',');
			sb.append(stop1);

			if(printLengthInEdges){
				sb.append(',');
				sb.append(u2.length());
			}
			sb.append(',');
			sb.append(start2);
			sb.append(',');
			sb.append(stop2);
			
			return sb.toString();
		}
		
		private void setOffsetContradiction(boolean b){
			assert(b!=offsetContradiction()) : b+", "+offsetContradiction();
			if(b){flags|=OFFSET_CONTRADICTION_MASK;}
			else{flags&=~OFFSET_CONTRADICTION_MASK;}
			assert(b==offsetContradiction()) : b+", "+offsetContradiction();
		}
		
		private void setMultiJoin(boolean b){
			assert(b!=multiJoin()) : b+", "+multiJoin();
			if(b){flags|=MULTIJOIN_MASK;}
			else{flags&=~MULTIJOIN_MASK;}
			assert(b==multiJoin()) : b+", "+multiJoin();
		}
		
		private void setVisited(boolean b){
			assert(b!=visited()) : b+", "+visited();
			if(b){flags|=VISITED_MASK;}
			else{flags&=~VISITED_MASK;}
			assert(b==visited()) : b+", "+visited();
		}
		
		private void setCyclic(boolean b){
			assert(b!=cyclic()) : b+", "+cyclic();
			if(b){flags|=CYCLIC_MASK;}
			else{flags&=~CYCLIC_MASK;}
			assert(b==cyclic()) : b+", "+cyclic();
		}
		
		private void setInvalid(boolean b){
			assert(b!=invalid()) : b+", "+invalid();
			assert(b!=mst()) : b+", "+mst()+", "+invalid();
			if(b){flags|=INVALID_MASK;}
			else{flags&=~INVALID_MASK;}
			assert(b==invalid()) : b+", "+invalid();
		}
		
		private void setMst(boolean b){
			assert(b!=mst()) : b+", "+mst();
			assert(b!=invalid()) : b+", "+mst()+", "+invalid();
			if(b){flags|=MST_MASK;}
			else{flags&=~MST_MASK;}
			assert(b==mst()) : b+", "+mst();
		}
		
		public void clearVolatileFlags(){
			flags=0;
//			flags=flags&~(MULTIJOIN_MASK|VISITED_MASK|CANON_CONTRADICTION_MASK|CYCLIC_MASK|OFFSET_CONTRADICTION_MASK|INVALID_MASK);
//			assert(!canonContradiction());
//			assert(!offsetContradiction());
//			assert(!multiJoin());
//			assert(!visited());
//			assert(!cyclic());
//			assert(!invalid());
		}
		
		public boolean offsetContradiction(){return (OFFSET_CONTRADICTION_MASK&flags)==OFFSET_CONTRADICTION_MASK;}
		public boolean multiJoin(){return (MULTIJOIN_MASK&flags)==MULTIJOIN_MASK;}
		public boolean visited(){return (VISITED_MASK&flags)==VISITED_MASK;}
		public boolean cyclic(){return (CYCLIC_MASK&flags)==CYCLIC_MASK;}
		public boolean invalid(){return (INVALID_MASK&flags)==INVALID_MASK;}
		public boolean mst(){return (MST_MASK&flags)==MST_MASK;}
		public boolean contradiction(){return offsetContradiction();}

		private static final long VISITED_MASK=(1L<<0);
		private static final long MULTIJOIN_MASK=(1L<<1);
		private static final long CYCLIC_MASK=(1L<<2);
		private static final long CANON_CONTRADICTION_MASK=(1L<<3);
		private static final long OFFSET_CONTRADICTION_MASK=(1L<<4);
		private static final long INVALID_MASK=(1L<<5);
		private static final long MST_MASK=(1L<<6);
		
		Unit u1;
		Unit u2;
		int type;
		int start1;
		int start2;
		int stop1;
		int stop2;
		
		long flags=0;
		
		final int overlapLen;
		final int mismatches;
		final int edits;
	}

	/**
	 * @return
	 */
	private int determineCluster2(final int uid) {
		assert(clusterNumbers!=null);
		boolean stable=false;
		int cluster=uid;
		while(!stable){
			cluster=clusterNumbers.get(uid);
			if(cluster==0 || cluster==uid){return cluster;}
			assert(cluster<=uid);
			final int next=determineCluster2(cluster);
			if(next>=cluster){return cluster;}
			stable=clusterNumbers.compareAndSet(uid, cluster, next);
		}
		return cluster;
	}
	

	private int mergeClusterIds(int cluster1, int cluster2) {
		assert(clusterNumbers!=null);
		
//		System.err.println("Merging clusters "+cluster1+" and "+cluster2);
		
		while(cluster1!=cluster2){
			int min=Tools.min(cluster1, cluster2);
			if(cluster1!=min){
				assert(cluster1>min);
				boolean b=clusterNumbers.compareAndSet(cluster1, cluster1, min);
				if(!b){
					cluster1=determineCluster2(cluster1);
					min=Tools.min(cluster1, cluster2);
				}
			}
			if(cluster2!=min){
				assert(cluster2>min);
				boolean b=clusterNumbers.compareAndSet(cluster2, cluster2, min);
				if(!b){
					cluster2=determineCluster2(cluster2);
					min=Tools.min(cluster1, cluster2);
				}
			}
		}
//		System.err.println("Returning "+cluster1);
		return cluster1;
	}
	
	private class Unit implements Comparable<Unit>, Serializable {
		
		/**
		 * 
		 */
		private static final long serialVersionUID = 5232407003873807738L;

		public Unit(Read r_){
//			this(r_, canonical_, canonical_ ? hash(r_.bases) : hashReversed(r_.bases));
			this(r_, hash(r_.bases), hashReversed(r_.bases));
		}
		
		public Unit(Read r_, long codeF_, long codeR_){
			r=r_;
			code1=Tools.min(codeF_, codeR_);
			code2=Tools.max(codeF_, codeR_);
			long f=r.length();
			prefix1=hashTip(r.bases, true, k, 0);
			suffix1=hashTip(r.bases, false, k, 0);
			if(r.length()>bitsPerBase*k){
				prefix2=hashTip(r.bases, true, k, k);
				suffix2=hashTip(r.bases, false, k, k);
			}
//			System.err.println(r.id+", "+prefix1+", "+suffix1+", "+code1+", "+code2);
			if(r.pairnum()==1){f|=PAIRNUM_MASK;}
			flags=f;
			assert(length()==r.length());
			assert(pairnum()==r.pairnum());
			if(parseDepth){
				int[] quad=KmerNormalize.parseDepth(r.id, null);
				if(quad!=null){depth=quad[r.pairnum()];}
			}
		}
		
		int determineCluster() {
			return determineCluster2(unitID);
		}
		
		public void absorbMatch(Unit u){
			
			assert(code1==u.code1 && code2==u.code2 && length()==u.length());
			if(r==null || u.r==null){return;}
			u.r.setDiscarded(true);
			final byte[] bases1=r.bases, bases2=u.r.bases;
			final byte[] quals1=r.quality, quals2=u.r.quality;
			
			assert((r.mate==null) == (u.r.mate==null));
			
			if(r.mate!=null && !u.r.mate.discarded()){
				((Unit)r.mate.obj).absorbMatch((Unit)u.r.mate.obj);
			}
			if(quals1==null || quals2==null){return;}
			
			for(int i=0; i<bases1.length; i++){
				byte b1=bases1[i], b2=bases2[i];
				if(!isFullyDefined(b1) && isFullyDefined(b2)){bases1[i]=b2;}
				else{assert(b1==b2);}
				if(quals1!=null && quals2!=null){
					quals1[i]=Tools.max(quals1[i], quals2[i]);
				}
			}
		}
		
		public boolean alreadyHas(Overlap o){
			if(overlapList==null){return false;}
			for(int i=0; i<overlapList.size(); i++){
				Overlap o2=overlapList.get(i);
				if(o.equals(o2)){
					assert(overlapList.contains(o));
					assert(o2.equals(o));
					return true;
				}
			}
			assert(!overlapList.contains(o));
			return false;
		}
		
		public ArrayList<Unit> makeCluster() {
			assert(!visited());
			assert(!clustered());
			assert(valid());
//			assert(set.isEmpty());
			ArrayList<Unit> cluster=new ArrayList<Unit>(overlapList==null ? 1 : overlapList.size()+1);
			cluster.add(this);
			setClustered(true);
			
			int added=1;
			for(int i=0; i<cluster.size(); i++){
				Unit u=cluster.get(i);
				added+=u.visit(cluster);
			}
			
			assert(added==cluster.size());
			return cluster;
		}
		
		public int visit(ArrayList<Unit> cluster) {
			assert(!visited());
			assert(clustered());
			assert(valid());
//			assert(cluster.contains(this));
			setVisited(true);
			int added=0;
			
			if(r!=null && r.mate!=null){
				Unit u2=(Unit)r.mate.obj;
				assert(u2!=this);
				assert(u2.valid());
				if(!u2.clustered()){
					u2.setClustered(true);
					cluster.add(u2);
					added++;
				}
			}
			
			if(overlapList!=null){
				for(Overlap o : overlapList){
					Unit u2=(o.u1==this ? o.u2 : o.u1);
					assert(o.u1==this || o.u2==this);
					assert(u2!=this);
					assert(u2.valid());
					if(!u2.clustered()){
						u2.setClustered(true);
						cluster.add(u2);
						added++;
					}
				}
			}
			return added;
		}
		
		public boolean isTransitive(){
			assert(valid());
			if(overlapList==null || overlapList.size()==0){return true;}
			for(Overlap o : overlapList){
				assert(o.u1==this || o.u2==this);
				Unit u2=(o.u1==this ? o.u2 : o.u1);
				assert(u2!=this);
				if(u2.overlapList==null){
					return false;
				}else{
					boolean found=false;
					for(Overlap o2 : u2.overlapList){
						if(o2.u1==this || o2.u2==this){
							found=true; break;
						}
					}
					if(!found){return false;}
				}
			}
			return true;
		}
		
		public boolean isPerfectlyTransitive(){
			assert(valid());
			if(overlapList==null || overlapList.size()==0){return true;}
			for(Overlap o : overlapList){
				assert(o.u1==this || o.u2==this);
				Unit u2=(o.u1==this ? o.u2 : o.u1);
				assert(u2!=this);
				if(u2.overlapList==null){
					return false;
				}else{
					boolean found=false;
					for(Overlap o2 : u2.overlapList){
						if(o2==o){
							found=true; break;
						}
					}
					if(!found){return false;}
				}
			}
			return true;
		}
		
		public boolean isNonRedundant(){
			assert(valid());
			if(overlapList==null || overlapList.size()==0){return true;}
			for(int i=0; i<overlapList.size(); i++){
				Overlap a=overlapList.get(i);
				for(int j=0; j<overlapList.size(); j++){
					Overlap b=overlapList.get(j);
					if((i==j)!=(a.equals(b))){
						return false;
					}
				}
			}
			return true;
		}
		
		public boolean contains(Unit u2, int loc, LongM key, BandedAligner bandy, int tableNum) {
			if(verbose){System.err.println("contains: Considering key "+key+", unit "+u2);}
			if(u2.length()>this.length()){return false;} //Smaller things cannot contain larger things
			if(u2.length()==this.length()){//For equal size, only one can contain the other the other
				boolean pass=false;
				if(!pass && u2.unitID>this.unitID){return false;}else if(u2.unitID<this.unitID){pass=true;}
				if(!pass && u2.r.numericID>this.r.numericID){return false;}else if(u2.r.numericID<this.r.numericID){pass=true;}
				if(!pass && u2.code1>this.code1){return false;}else if(u2.code1<this.code1){pass=true;}
				if(!pass && u2.r.bases.hashCode()>this.r.bases.hashCode()){return false;}else if(u2.r.bases.hashCode()<this.r.bases.hashCode()){pass=true;}
				if(!pass){return false;}
			}
			if(minLengthPercent>0 && (u2.length()*100f/length())<minLengthPercent){return false;}
			if(verbose){System.err.println("passed tests for "+u2.length()+", "+this.length()+", "+u2.unitID+", "+this.unitID);}
			
			assert(u2.code1!=code1 || u2.code2!=code2 || u2.length()!=length() || (r!=null && r.mate!=null) || //REQUIRE_MATCHING_NAMES ||
					(u2.prefix1!=prefix1 && u2.suffix1!=suffix1)) :
						"Collision? \n"+this+"\n"+u2+"\n"+r+"\n"+u2.r;
			
			if(tableNum==0){
				if(key.value()==u2.prefix1){
					if(verbose){System.err.println("Containment A1");}
					if(containsForward(u2, loc-k2, bandy, tableNum==0)){return true;}
				}
				if(key.value()==u2.suffix1){
					if(verbose){System.err.println("Containment B1");}
					if(containsReverse(u2, loc, bandy, tableNum==0)){return true;}
				}
			}else{
				if(key.value()==u2.prefix2){
					if(verbose){System.err.println("Containment A2");}
					if(containsForward(u2, loc-k2-k, bandy, tableNum==0)){return true;}
				}
				if(key.value()==u2.suffix2){
					if(verbose){System.err.println("Containment B2");}
					if(containsReverse(u2, loc+k, bandy, tableNum==0)){return true;}
				}
			}
			return false;
		}
		
		private boolean containsForward(Unit u2, int start, BandedAligner bandy, boolean earlyExit) {
			if(start+u2.length()>length() || start<0 || start>=length()){return false;}
//			if(true){return false;}
			if(u2.r!=null){
				final byte[] a=bases(), b=u2.bases();
				int mismatches=0, maxMismatches=calcMaxEdits(maxSubs, minIdentityMult, b.length);
				
				for(int i=start, j=0; j<b.length; i++, j++){
					byte aa=a[i];
					byte bb=b[j];
					if(aa!=bb){
						if(exact || (isFullyDefined(aa) && isFullyDefined(bb))){
							if(earlyExit && j<k2){return false;}
							if((mismatches=mismatches+1)>maxMismatches){
								if(bandy==null || maxEdits<1){return false;}
								int edits=bandy.alignForward(b, a, 0, start, maxEdits, exact);
								assert(b.length<k || a.length<k || bandy.lastRow>=k2 || edits>maxEdits) : b.length+", "+k+", "+bandy.lastRow+", "+edits;
								return edits<=maxEdits && bandy.score()>4*edits;
							}
						}
					}
				}
				return true;
			}else{
				assert(false) : "TODO: Verify by hashing and checking both tips";
				return false;
			}
		}
		
		private boolean containsReverse(Unit u2, int start, BandedAligner bandy, boolean earlyExit) {
			if(start+1<u2.length() || start<0 || start>=length()){return false;}
//			if(true){return false;}
			if(u2.r!=null){
				final byte[] a=bases(), b=u2.bases();
				int mismatches=0, maxMismatches=calcMaxEdits(maxSubs, minIdentityMult, b.length);

				for(int i=start, j=b.length-1, iprefix=start-k2; j>=0; i--, j--){
					byte aa=a[i];
					byte bb=b[j];
					if(aa!=bb){
						if(exact || (isFullyDefined(aa) && isFullyDefined(bb))){
							if(earlyExit && i>iprefix){return false;}
							if((mismatches=mismatches+1)>maxMismatches){
								if(bandy==null || maxEdits<1){return false;}
								int edits=bandy.alignReverse(b, a, b.length-1, start, maxEdits, exact);
								assert(b.length<k || a.length<k || bandy.lastRow>=k2 || edits>maxEdits) : b.length+", "+k+", "+bandy.lastRow+", "+edits;
								return edits<=maxEdits && bandy.score()>4*edits;
							}
						}
					}
				}
				return true;
			}else{
				assert(false) : "TODO: Verify by hashing and checking both tips";
				return false;
			}
		}
		
		
		public boolean depthCongruent(int aa, int bb){
			if(aa<5 && bb<5){return true;}
			final int a=Tools.max(1, Tools.min(aa, bb));
			final int b=Tools.max(aa, bb);
			return a*depthRatio>=b;
		}
		
		public boolean overlaps(Unit u2, int loc, LongM key, BandedAligner bandy, int tableNum, int editLimit) {
//			return makeOverlap(u2, loc, key, bandy, earlyExit)!=null;
			
//			assert(false) : "TODO";
			if(verbose){System.err.println("overlaps: Considering key "+key+", unit "+u2);}
			if(parseDepth && !depthCongruent(depth, u2.depth)){return false;}
			if(minLengthPercent>0){
				final int len1=length(), len2=u2.length();
				if(Tools.min(len1, len2)*100f/Tools.max(len1, len2)<minLengthPercent){return false;}
			}
			assert(u2.code1!=code1 || u2.code2!=code2 || u2.length()!=length() ||
					(u2.prefix1!=prefix1 && u2.suffix1!=suffix1)) :
						"Collision? \n"+this+"\n"+u2+"\n"+r+"\n"+u2.r;
			
			
			if(tableNum==0){
				if(key.value()==u2.prefix1){
					if(verbose){System.err.println("Testing overlaps A1");}
					if(overlapsForward(u2, loc-k2, 0, bandy, tableNum==0, editLimit)){
						if(verbose){System.err.println("Found Overlap A1F");}
						return true;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
				if(key.value()==u2.suffix1){
					if(verbose){System.err.println("Testing overlaps B1");}
					if(verbose){System.err.println("Testing overlaps B1R");}
					if(overlapsReverse(u2, loc, u2.length()-1, bandy, tableNum==0, editLimit)){
						if(verbose){System.err.println("Found Overlap B1R");}
						return true;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
			}else{
				if(key.value()==u2.prefix2){
					if(verbose){System.err.println("Testing overlaps A2");}
					if(overlapsForward(u2, loc-k2-k, 0, bandy, tableNum==0, editLimit)){
						if(verbose){System.err.println("Found Overlap A2F");}
						return true;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
				if(key.value()==u2.suffix2){
					if(verbose){System.err.println("Testing overlaps B2");}
					if(overlapsReverse(u2, loc+k, u2.length()-1, bandy, tableNum==0, editLimit)){
						if(verbose){System.err.println("Found Overlap B2R");}
						return true;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
			}
			
			return false;
		}
		
		/**
		 * @param u2
		 * @param loc
		 * @param key
		 * @return
		 */
		protected Overlap makeOverlap(Unit u2, int loc, LongM key, BandedAligner bandy, int tableNum) {
			if(verbose){System.err.println("makeOverlap: Considering key "+key+", unit "+u2);}
			if(parseDepth && !depthCongruent(depth, u2.depth)){return null;}
			if(minLengthPercent>0){
				final int len1=length(), len2=u2.length();
				if(Tools.min(len1, len2)*100f/Tools.max(len1, len2)<minLengthPercent){return null;}
			}
			assert(u2.code1!=code1 || u2.code2!=code2 || u2.length()!=length() || (r!=null && r.mate!=null) ||
					(u2.prefix1!=prefix1 && u2.suffix1!=suffix1)) :
						"Collision? \n"+this+"\n"+u2+"\n"+r+"\n"+u2.r;
			
			
			Overlap o=null;
			if(tableNum==0){
				if(key.value()==u2.prefix1){
					if(verbose){System.err.println("\nTesting makeOverlap A1F");}
					if((o=makeOverlapForward(u2, loc-k2, bandy, tableNum==0))!=null){
						if(verbose){System.err.println("Made Overlap A1F");}
						return o;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
				if(key.value()==u2.suffix1){
					if(verbose){System.err.println("\nTesting makeOverlap B1R");}
					if((o=makeOverlapReverse(u2, loc, bandy, tableNum==0))!=null){
						if(verbose){System.err.println("Made Overlap B1R");}
						return o;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
			}else{
				if(key.value()==u2.prefix2){
					if(verbose){System.err.println("\nTesting makeOverlap A2F");}
					if((o=makeOverlapForward(u2, loc-k2-k, bandy, tableNum==0))!=null){
						if(verbose){System.err.println("Made Overlap A2F");}
						return o;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
				if(key.value()==u2.suffix2){
					if(verbose){System.err.println("\nTesting makeOverlap B2R");}
					if((o=makeOverlapReverse(u2, loc+k, bandy, tableNum==0))!=null){
						if(verbose){System.err.println("Made Overlap B2R");}
						return o;
					}
					if(verbose){System.err.println("No Overlap.");}
				}
			}
			return o;
		}

		private boolean overlapsForward(Unit u2, int start1, int start2, BandedAligner bandy, boolean earlyExit, int maxEdits) {
			if(verbose){System.err.println("overlapsForward(u1="+this.name()+", u2="+u2.name()+", start1="+start1+", start2="+start2+", earlyExit="+earlyExit+")");}

			final int len1=length(), len2=u2.length();
			if(start1<0){
				start2-=start1;
				start1=0;
				if(verbose){System.err.println("Modified: start1="+start1+", start2="+start2);}
			}
			int overlapLength=Tools.min(len1-start1, len2-start2);
			int overlapLength2=Tools.max(len1-start1, len2-start2);
			int stop1=start1+overlapLength-1, stop2=start2+overlapLength-1;
			if(verbose){System.err.println("Calculated stop1="+stop1+", stop2="+stop2+", overlapLength="+overlapLength);}

			if(!allowAllContainedOverlaps || overlapLength>Tools.min(len1, len2)){
				if(verbose){
					System.err.println("Side block. allowAllContainedOverlaps="+allowAllContainedOverlaps+", minOverlapCluster="+minOverlapCluster);
					System.err.println("start1="+start1+", stop1="+stop1+", len1="+len1+", start2="+start2+", stop2="+stop2+", len2="+len2+
							", overlapLen="+overlapLength+", maxEdits="+maxEdits);
				}
				if(overlapLength2<minOverlapCluster){return false;}
				if(minOverlapPercentCluster>0f && (overlapLength2*100f/Tools.min(len1, len2))<minOverlapPercentCluster){return false;}
			}

			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxSubs, minIdentityMult, overlapLength);

			if(verbose){
				System.err.println("start1="+start1+", stop1="+stop1+", len1="+len1+", start2="+start2+", stop2="+stop2+", len2="+len2+
						", overlapLen="+overlapLength+", maxMismatches="+maxMismatches+", maxEdits="+maxEdits);
			}
			
			for(int i=start1, j=start2; j<=stop2; i++, j++){
				byte aa=a[i];
				byte bb=b[j];
				if(aa!=bb){
					if(exact || (isFullyDefined(aa) && isFullyDefined(bb))){
						if(earlyExit && j<k2){return false;}
						mismatches++;
						if(verbose){System.err.println("i="+i+", j="+j+", "+(char)aa+"!="+(char)bb+", mismatches="+mismatches+"/"+maxMismatches);}
						if(mismatches>maxMismatches){
							if(bandy==null || maxEdits<1){return false;}
							if(verbose){System.err.println("Mismatches exceeded maximum, attempting banded alignment.");}
							int edits=bandy.alignForward(b, a, 0, start1, maxEdits, exact);
							assert(b.length<k || a.length<k || bandy.lastRow>=k2 || edits>maxEdits) : b.length+", "+k+", "+bandy.lastRow+", "+edits;
							stop2=bandy.lastQueryLoc;
							stop1=bandy.lastRefLoc;
							return edits<=maxEdits && bandy.score()>2*edits; //Set at 2*edits instead of 4*edits to prevent assertion errors when reversing alignment
						}
					}
				}
			}
			return true;
		}

		private boolean overlapsReverse(Unit u2, int start1, int start2, BandedAligner bandy, boolean earlyExit, int maxEdits) {
			if(verbose){System.err.println("overlapsReverse(u1="+this.name()+", u2="+u2.name()+", start1="+start1+", start2="+start2+", earlyExit="+earlyExit+")");}

			final int len1=length(), len2=u2.length();
			if(start1>=len1){
				start2-=(start1-len1+1);
				start1=len1-1;
				if(verbose){System.err.println("Modified: start1="+start1+", start2="+start2);}
			}
			final int overlapLength=Tools.min(start1+1, start2+1);
			final int overlapLength2=Tools.max(start1+1, start2+1);
			int stop1=start1-overlapLength+1, stop2=start2-overlapLength+1;
			if(verbose){System.err.println("Calculated stop1="+stop1+", stop2="+stop2+", overlapLength="+overlapLength);}

			if(!allowAllContainedOverlaps || overlapLength>Tools.min(len1, len2)){
				if(overlapLength2<minOverlapCluster){return false;}
				if(minOverlapPercentCluster>0f && (overlapLength2*100f/Tools.min(len1, len2))<minOverlapPercentCluster){return false;}
			}

			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxSubs, minIdentityMult, b.length);
			
			if(verbose){
				System.err.println("start1="+start1+", stop1="+stop1+", len1="+len1+", start2="+start2+", stop2="+stop2+", len2="+len2+
						", overlapLen="+overlapLength+", maxMismatches="+maxMismatches+", maxEdits="+maxEdits);
			}
			
			for(int i=start1, j=start2, iprefix=start1-k2; i>=stop1; i--, j--){
				byte aa=a[i];
				byte bb=b[j];
				if(aa!=bb){
					if(exact || (isFullyDefined(aa) && isFullyDefined(bb))){
						if(earlyExit && i>iprefix){return false;}
						mismatches++;
						if(verbose){System.err.println("i="+i+", j="+j+", "+(char)aa+"!="+(char)bb+", mismatches="+mismatches+"/"+maxMismatches);}
						if(mismatches>maxMismatches){
							if(bandy==null || maxEdits<1){return false;}
							if(verbose){System.err.println("Mismatches exceeded maximum, attempting banded alignment.");}
							int edits=bandy.alignReverse(b, a, b.length-1, start1, maxEdits, exact);
							assert(b.length<k || a.length<k || bandy.lastRow>=k2 || edits>maxEdits) : b.length+", "+k+", "+bandy.lastRow+", "+edits;
							stop2=bandy.lastQueryLoc;
							stop1=bandy.lastRefLoc;
							return edits<=maxEdits && bandy.score()>2*edits; //Set at 2*edits instead of 4*edits to prevent assertion errors when reversing alignment
						}
					}
				}
			}
			return true;
		}

		private Overlap makeOverlapForward(Unit u2, int start1, BandedAligner bandy, boolean earlyExit) {
			if(verbose){System.err.println("makeOverlapForward(u1="+this.name()+", u2="+u2.name()+", start="+start1+", earlyExit="+earlyExit+")");}
			final int len1=length(), len2=u2.length();
			int start2=0;
			if(start1<0){
				start2-=start1;
				start1=0;
			}
			final int overlapLength=Tools.min(len1-start1, len2-start2);
			final int overlapLength2=Tools.max(len1-start1, len2-start2);
			int stop1=start1+overlapLength-1, stop2=start2+overlapLength-1;
			if(verbose){System.err.println("Calculated stop1="+stop1+", stop2="+stop2+", overlapLength="+overlapLength);}

			if(!allowAllContainedOverlaps || overlapLength>Tools.min(len1, len2)){
				if(overlapLength<minOverlapCluster){return null;}
				if(minOverlapPercentCluster>0f && (overlapLength*100f/Tools.min(len1, len2))<minOverlapPercentCluster){return null;}
			}

			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxSubs, minIdentityMult, overlapLength);

			if(verbose){
				System.err.println("start1="+start1+", stop1="+stop1+", len1="+len1+", start2="+start2+", stop2="+stop2+", len2="+len2+
						", overlapLen="+overlapLength+", maxMismatches="+maxMismatches+", maxEdits="+maxEdits);
			}

			for(int i=start1, j=start2; j<=stop2; i++, j++){
				byte aa=a[i];
				byte bb=b[j];
				if(aa!=bb){
					if(exact || (isFullyDefined(aa) && isFullyDefined(bb))){
						if(earlyExit && j<k2){return null;}
						mismatches++;
						if(verbose){System.err.println("i="+i+", j="+j+", "+(char)aa+"!="+(char)bb+", mismatches="+mismatches+"/"+maxMismatches);}
						if(mismatches>1 && bandy!=null){
							if(maxEdits<1){return null;}
							if(verbose){System.err.println("mismatches exceeded 1, attempting banded alignment.");}
							int edits=bandy.alignForward(b, a, start2, start1, maxEdits, exact);
							if(edits>maxEdits || bandy.score()<=4*edits){
								if(verbose){System.err.println((edits>maxEdits ? "Too many edits" : "Alignment score too low")+"; returning null.");}
								return null;
							}
							assert(b.length<k || a.length<k || bandy.lastRow>=k2 || edits>maxEdits) : b.length+", "+k+", "+bandy.lastRow+", "+edits;
							stop2=bandy.lastQueryLoc;
							stop1=bandy.lastRefLoc;
//							if(bandy.lastOffset>0){//Ref longer than query
//								for(int k=0; k<bandy.lastOffset; k++){
//									if(stop1+1<=len1){stop1++;}
//									else{stop2--;}//I don't think this can happen
//								}
//							}else if(bandy.lastOffset<0){//Query longer than ref
//								for(int k=0; k>bandy.lastOffset; k--){
//									if(stop2+1<=len2){stop2++;}
//									else{stop1--;}
//								}
//							}
							return new Overlap(this, u2, FORWARD, start1, start2, stop1, stop2, overlapLength, edits, edits, bandy);
						}else if(mismatches>maxMismatches){return null;}
					}
				}
			}
			return new Overlap(this, u2, FORWARD, start1, start2, stop1, stop2, overlapLength, mismatches, 0, bandy);
		}

		private Overlap makeOverlapReverse(Unit u2, int start1, BandedAligner bandy, boolean earlyExit) {
			if(verbose){System.err.println("makeOverlapReverse(u1="+this.name()+", u2="+u2.name()+", start="+start1+", earlyExit="+earlyExit+")");}

			final int len1=length(), len2=u2.length();
			int start2=len2-1;
			if(start1>=len1){
				start2-=(start1-len1+1);
				start1=len1-1;
			}
			final int overlapLength=Tools.min(start1+1, start2+1);
			final int overlapLength2=Tools.max(start1+1, start2+1);
			int stop1=start1-overlapLength+1, stop2=start2-overlapLength+1;
			if(verbose){System.err.println("Calculated stop1="+stop1+", stop2="+stop2+", overlapLength="+overlapLength);}

			if(!allowAllContainedOverlaps || overlapLength>Tools.min(len1, len2)){
				if(overlapLength<minOverlapCluster){return null;}
				if(minOverlapPercentCluster>0f && (overlapLength*100f/Tools.min(len1, len2))<minOverlapPercentCluster){return null;}
			}

			final byte[] a=bases(), b=u2.bases();
			assert(a!=null && b!=null) : "Null bases for "+code1+" or "+u2.code1;
			int mismatches=0, maxMismatches=calcMaxEdits(maxSubs, minIdentityMult, b.length);

			if(verbose){
				System.err.println("start1="+start1+", stop1="+stop1+", len1="+len1+", start2="+start2+", stop2="+stop2+", len2="+len2+
						", overlapLen="+overlapLength+", maxMismatches="+maxMismatches+", maxEdits="+maxEdits);
			}

			for(int i=start1, j=start2, iprefix=start1-k2; i>=stop1; i--, j--){
				byte aa=a[i];
				byte bb=b[j];
				if(aa!=bb){
					if(exact || (isFullyDefined(aa) && isFullyDefined(bb))){
						if(earlyExit && i>iprefix){return null;}
						mismatches++;
						if(verbose){System.err.println("i="+i+", j="+j+", "+(char)aa+"!="+(char)bb+", mismatches="+mismatches+"/"+maxMismatches);}
						if(mismatches>1 && bandy!=null){
							if(maxEdits<1){return null;}
							if(verbose){System.err.println("mismatches exceeded 1, attempting banded alignment.");}
							int edits=bandy.alignReverse(b, a, start2, start1, maxEdits, exact);
							if(edits>maxEdits || bandy.score()<=4*edits){
								if(verbose){System.err.println((edits>maxEdits ? "Too many edits" : "Alignment score too low")+"; returning null.");}
								return null;
							}
							assert(b.length<k || a.length<k || bandy.lastRow>=k2 || edits>maxEdits) : b.length+", "+k+", "+bandy.lastRow+", "+edits;
							stop2=bandy.lastQueryLoc;
							stop1=bandy.lastRefLoc;
//							if(bandy.lastOffset>0){//Ref longer than query
//								for(int k=0; k<bandy.lastOffset; k++){
//									if(stop1>0){stop1--;}
//									else{stop2++;}//I don't think this can happen
//								}
//							}else if(bandy.lastOffset<0){//Query longer than ref
//								for(int k=0; k>bandy.lastOffset; k--){
//									if(stop2>0){stop2--;}
//									else{stop1++;}
//								}
//							}
							return new Overlap(this, u2, REVERSE, start1, start2, stop1, stop2, overlapLength, edits, edits, bandy);
						}else if(mismatches>maxMismatches){return null;}
					}
				}
			}
			return new Overlap(this, u2, REVERSE, start1, start2, stop1, stop2, overlapLength, mismatches, 0, bandy);
		}
		
		@Override
		public int compareTo(Unit b) {
			int x=comparePairedRC(this, b);
			return x;
		}
		
		@Override
		public boolean equals(Object b){return equals((Unit)b);}
		public boolean equals(Unit b){
			boolean x=pairedEqualsRC(this, b);
			return x;
		}
		
		@Override
		public int hashCode(){
			return (int)((code1^(code1>>>32))&0xFFFFFFFFL);
		}
		
		private synchronized void setValid(boolean b){
			assert(b!=valid());
//			if(!b){System.err.println("Setting invalid "+name());}
			if(b){flags&=~INVALID_MASK;}
			else{flags|=INVALID_MASK;}
			assert(b==valid());
		}
		
		private synchronized void setClustered(boolean b){
			assert(b!=clustered());
			if(b){flags|=CLUSTER_MASK;}
			else{flags&=~CLUSTER_MASK;}
			assert(b==clustered());
		}
		
		private void setVisited(boolean b){
			assert(b!=visited());
			if(b){flags|=VISIT_MASK;}
			else{flags&=~VISIT_MASK;}
			assert(b==visited());
		}
		
		private synchronized void setOffset(int x){
			offset=x;
			setOffsetValid(true);
		}
		
		private synchronized void setOffsetValid(boolean b){
			assert(!offsetValid());
			if(b){flags|=OFFSET_VALID_MASK;}
			else{flags&=~OFFSET_VALID_MASK;}
			assert(b==offsetValid());
		}
		
		private synchronized void setOffsetContradiction(boolean b){
//			assert(b!=offsetContradiction());
			assert(offsetValid());
			if(b){flags|=OFFSET_CONTRADICTION_MASK;}
			else{flags&=~OFFSET_CONTRADICTION_MASK;}
			assert(b==offsetContradiction());
		}
		
		/** Return true if 'this' should be the first Unit in the overlap object */
		public boolean firstInOverlap(Unit u2){
			assert(this!=u2) : "\n"+this.r+"\n"+u2.r;
			if(u2.length()!=length()){return u2.length()<length();}
			if(u2.code1!=code1){return u2.code1<code1;}
			if(u2.code2!=code2){return u2.code2<code2;}
			int x=compareTo(u2);
			assert(x!=0 || (r!=null && r.mate!=null));
			if(x!=0){return x>=0;}
			return r.numericID>=u2.r.numericID;
		}
		
		public final boolean inSet(){
			if(subsetCount<2){return true;}
			if(r.pairnum()>0){return ((Unit)r.mate.obj).inSet();}
			return ((code1&Long.MAX_VALUE)%subsetCount)==subset;
		}
		
		public byte[] bases(){return r==null ? null : r.bases;}

		public String name(){return r!=null ? r.id : null /*code+""*/;}
		@Override
		public String toString(){return "("+name()+","+code1+","+code2+","+length()+","+prefix1+","+suffix1+","+"c"+",d="+depth+")";}
		
		
		public final Read r;
		public final long code1;
		public final long code2;
		public long prefix1=-1;
		public long suffix1=-1;
		public long prefix2=-1;
		public long suffix2=-1;
		/** Distance of leftmost side of this read relative to leftmost side of root.
		 * Assumes everything is in 'forward' orientation. */
		public int offset=-999999999;
		public int depth=1;
//		private boolean valid=true;

		public int unitID;

		public ArrayList<Overlap> overlapList;
		
		private long flags;
		/** True if this contig should be output, false if not */
		public final boolean valid(){return (INVALID_MASK&flags)==0;}
		/** Length of this contig */
		public final int length(){return (int)(LEN_MASK&flags);}
		/** Position of this contig relative to root */
		public final int offset(){
			assert(offsetValid());
			return offset;
		}
		public int pairnum(){return (PAIRNUM_MASK&flags)==PAIRNUM_MASK ? 1 : 0;}
		
		public void clearVolatileFlags(){
			flags=flags&~(VISIT_MASK|OFFSET_VALID_MASK|OFFSET_CONTRADICTION_MASK);
			assert(!visited());
			assert(!offsetValid());
			assert(!offsetContradiction());
		}
		
		public boolean visited(){return (VISIT_MASK&flags)==VISIT_MASK;}
		public boolean clustered(){return (CLUSTER_MASK&flags)==CLUSTER_MASK;}
		public boolean offsetValid(){return (OFFSET_VALID_MASK&flags)==OFFSET_VALID_MASK;}
		public boolean offsetContradiction(){return (OFFSET_CONTRADICTION_MASK&flags)==OFFSET_CONTRADICTION_MASK;}
		public boolean contradiction(){return offsetContradiction();}

		private static final long LEN_MASK=0x7FFFFFFFL;
//		private static final long CANON_MASK=(1L<<33);
		private static final long INVALID_MASK=(1L<<34);
		private static final long VISIT_MASK=(1L<<35);
		private static final long CLUSTER_MASK=(1L<<36);
//		private static final long CANONICIZED_MASK=(1L<<37);
//		private static final long CANON_CONTRADICTION_MASK=(1L<<38);
		private static final long OFFSET_VALID_MASK=(1L<<39);
		private static final long OFFSET_CONTRADICTION_MASK=(1L<<40);
		private static final long PAIRNUM_MASK=(1L<<41);
	}
	
	private static final class UnitOffsetComparator implements Comparator<Unit> {
		
		UnitOffsetComparator(){}
		
		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(Unit a, Unit b) {
			if(a.offsetValid() && b.offsetValid()){
				int x=a.offset()-b.offset();
				if(x!=0){return x;}
			}else{
				if(a.offsetValid()){return -1;}
				if(b.offsetValid()){return 1;}
			}
			return a.compareTo(b);
		}
		
	}
	
	private static final class ClusterLengthComparator implements Comparator<ArrayList<Unit>> {
		
		ClusterLengthComparator(){}
		
		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		@Override
		public int compare(ArrayList<Unit> a, ArrayList<Unit> b) {
			if(a.size()!=b.size()){return b.size()-a.size();}
			if(a.isEmpty() && b.isEmpty()){return 0;}
			return a.get(0).compareTo(b.get(0));
		}
		
	}
	
	/** Makes a nmer (e.g., tetramer) profile of a cluster */
	private final float[] makeNmerProfile(ArrayList<Unit> alu, long[] array_){
		final int nbits=bitsPerBase*nmerLength;
		final long[] array=(array_==null ? new long[maxNmer+1] : array_);
		final int mask=~((-1)<<(nbits));
		
		long keysCounted=0;
		
		for(Unit u : alu){
			byte[] bases=u.r.bases;
			int len=0;
			int kmer=0;
			for(byte b : bases){
				int x=AminoAcid.baseToNumber[b];
				if(x<0){
					len=0;
					kmer=0;
				}else{
					kmer=((kmer<<bitsPerBase)|x)&mask;
					len++;
					if(len>=nmerLength){
						keysCounted++;
						array[kmer]++;
					}
				}
			}
		}
		
		if(keysCounted==0){keysCounted=1;}
		final float mult=1f/keysCounted;
		
		float[] r=new float[array.length];
		for(int i=0; i<array.length; i++){
			r[i]=array[i]*mult;
			array[i]=0;
		}
		return r;
	}
	
	private static boolean isFullyDefined(byte b){
		return AminoAcid.acidToNumber[b]>=0;
	}
	
	private ConcurrentReadInputStream crisa[];
	
	private final ByteStreamWriter dupeWriter;

	private String[] in1=null;
	private String[] in2=null;
	private String out=null;
	private String clusterFilePattern=null;
	private String outbest=null;
	private String outdupe=null;
	private String outcsf=null;
	private String outgraph=null;
	private int maxNs=-1;
	private long maxReads=-1;
	public boolean errorState=false;
	boolean sort=false;
//	boolean ascending=true;
	boolean absorbContainment=true;
	boolean absorbMatch=true;
	boolean findOverlaps=false;
	boolean makeClusters=false;
	boolean processClusters=false;
	boolean renameClusters=false;
	boolean absorbOverlap=false;
	boolean storeSuffix=false;
	boolean storeName=true;
	boolean storeQuality=true;
	boolean exact=true;
	boolean uniqueNames=true;
	boolean maxSpanningTree=false;
	
	boolean removeCycles=true;
	boolean fixMultiJoins=true;
	boolean fixOffsetContradictions=false;
	boolean countTransitive=false;
	boolean countRedundant=false;
	
	private boolean multipleInputFiles=false;
	private boolean rigorousTransitive=false;
	private int numAffixMaps=1;
	private int maxAffixCopies=2000000000;
	private int maxEdges=2000000000;
	private int maxEdges2=2000000000;
	private boolean allowAllContainedOverlaps=false;
//	private boolean toUpperCase=false;
	
	/** Trim left bases of the read to this position (exclusive, 0-based) */
	private int forceTrimLeft=-1;
	/** Trim right bases of the read after this position (exclusive, 0-based) */
	private int forceTrimRight=-1;
	
	private boolean qTrimLeft=false;
	private boolean qTrimRight=false;
	private float trimQ=6;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;

	int maxEdits=0;
	int maxSubs=0;
	int bandwidth=9;
	final boolean customBandwidth;
	float minIdentity=100;
	float minIdentityMult=0;
	float minLengthPercent=0;
	int minOverlapCluster=100;
	int minOverlapMerge=1;
	float minOverlapPercentCluster=0;
	float minOverlapPercentMerge=0;

	private int minClusterSize=1;
	private int minClusterSizeForStats=1;
	private boolean pickBestRepresentativePerCluster=false;

	long readsProcessed=0;
	long basesProcessed=0;
	long collisions=0;
	long containments=0;
	long baseContainments=0;
	long containmentCollisions=0;
	long matches=0;
	long baseMatches=0;
	long overlaps=0;
	long baseOverlaps=0;
	long overlapCollisions=0;
	long addedToMain=0;

	private final int subset;
	private final int subsetCount;
	private final boolean subsetMode;
	
	private final int k;
	private final int k2;
	private final boolean EA=Shared.EA();
	
//	static final boolean amino=true;
	static final int bitsPerBase=5; //amino
//	static final boolean ignoreReverseComplement=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ReadComparator comparator=ReadLengthComparator.comparator;
	
	private static int tcount=0;
	
	private LinkedHashMap<Long, ArrayList<Unit>> codeMap=new LinkedHashMap<Long, ArrayList<Unit>>(4000000);
	private HashMap<LongM, ArrayList<Unit>> affixMap1=null;
	private HashMap<LongM, ArrayList<Unit>> affixMap2=null;
	private HashMap<LongM, ArrayList<Unit>>[] affixMaps=null;
	private ArrayDeque<ArrayList<Unit>> clusterQueue=null;
	private ArrayList<ArrayList<Unit>> processedClusters=null;
	private AtomicIntegerArray clusterNumbers=null;

	private static final UnitOffsetComparator UNIT_OFFSET_COMPARATOR=new UnitOffsetComparator();
	private static final ClusterLengthComparator CLUSTER_LENGTH_COMPARATOR=new ClusterLengthComparator();
	private static final long[][] hashcodes=makeCodes2(32, 1);
	private static final long[][] hashcodesR=makeCodes2(32, 2);
	
	public static final byte[] acidToNumber=new byte[128];
	
	public static final int nmerLength=4;
	public static final int maxNmer=((int) Math.pow(20, nmerLength))-1;
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	public static boolean showSpeed=true;
	public static boolean verbose=false;
	public static boolean preventTransitiveOverlaps=false;
	public static boolean ignoreAffix1=false;
	public static boolean parseDepth=false;
	public static boolean printLengthInEdges=false;
	public static float depthRatio=2;
	public static int MINSCAF=0;
	public static int THREADS=Shared.threads();
	public static int threadMaxReadsToBuffer=4000;
	public static int threadMaxBasesToBuffer=32000000;
	public static boolean DISPLAY_PROGRESS=true;
	public static boolean UNIQUE_ONLY=false;
	public static boolean REQUIRE_MATCHING_NAMES=false;
	public static boolean NUMBER_GRAPH_NODES=true;
	public static boolean ADD_PAIRNUM_TO_NAME=true;
	public static boolean HASH_NS=false;
	
	private static int reverseType(int type){return (type+2)%4;}
	public static final int FORWARD=0;
	public static final int REVERSE=2;
	public static final String[] OVERLAP_TYPE_NAMES=new String[] {"FORWARD", "FORWARDRC", "REVERSE", "REVERSERC"};
	public static final String[] OVERLAP_TYPE_ABBREVIATIONS=new String[] {"F", "FRC", "R", "RRC"};
	
	static{//All others are 0
		for(int i=0; i<acidToNumber.length; i++){
			byte b=AminoAcid.acidToNumber[i];
			acidToNumber[i]=(b<0 ? 21 : b);
		}
	}
	
}
