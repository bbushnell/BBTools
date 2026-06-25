package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Grades a VCF against a truth set and produces concordance metrics and histograms.
 * Replaces a multi-step pipeline (filtervcf + comparevcf) with a single fast tool.
 *
 * Computes TP, FP, FN counts with optional BED region restriction, VarFilter
 * pre-filtering, and neural network scoring. Writes cumulative QUAL and NN score
 * histograms suitable for precision/recall curve plotting.
 *
 * Counting follows the "marking" contract: the truth set is normalized and split
 * into simple vars, then each truth var is marked at most once when a call matches
 * it (recording the maximum matching-call score). TP = distinct marked truth vars,
 * FN = unmarked truth vars, FP = call vars that matched no truth. This avoids
 * double-counting when a compound site splits into redundant identical components.
 * Variant matching uses VCFLine equals/hashCode — the same identity contract as
 * CompareVCF.
 *
 * @author UMP45
 * @date June 8, 2026
 */
public class GradeVCF {

	public static void main(String[] args){
		Timer t=new Timer();
		GradeVCF x=new GradeVCF(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public GradeVCF(String[] args){

		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		varFilter.clear();
		boolean setVarFilter=false;

		String netFile=null;
		boolean autoCutoff=true;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("truth") || a.equals("giab")){
				truthFile=b;
			}else if(a.equals("bed") || a.equals("bedfile")){
				bedFile=b;
			}else if(a.equals("bedweight") || a.equals("bw")){
				bedWeight=Double.parseDouble(b);
			}else if(a.equals("invertbed") || a.equals("bedinvert") || a.equals("invbed")){
				invertBed=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("minscore") || a.equals("minqual") || a.equals("minphred")){
				minScore=Double.parseDouble(b);
			}else if(a.equals("net") || a.equals("netfile")){
				netFile=b;
				useNet=(b!=null);
			}else if(a.equals("netcutoff") || a.equals("cutoff")){
				if("auto".equalsIgnoreCase(b)){
					autoCutoff=true;
				}else{
					autoCutoff=false;
					netCutoff=Float.parseFloat(b);
				}
			}else if(a.equals("usenet") || a.equals("useann") || a.equals("usenn") || a.equals("nn")){
				useNet=Parse.parseBoolean(b);
			}else if(a.equals("netmode")){
				useNet=(b!=null);
				if(b!=null){FeatureVectorMaker.setMode(b);}
			}else if(a.equals("includescore")){
				VectorUMP45.includeScore=Parse.parseBoolean(b);
			}else if(a.equals("platform")){
				if(b.equalsIgnoreCase("illumina")){VectorUMP45.platform=VectorUMP45.PLATFORM_ILLUMINA;}
				else if(b.equalsIgnoreCase("pacbio")){VectorUMP45.platform=VectorUMP45.PLATFORM_PACBIO;}
				else if(b.equalsIgnoreCase("nanopore") || b.equalsIgnoreCase("ont")){VectorUMP45.platform=VectorUMP45.PLATFORM_NANOPORE;}
				else if(b.equalsIgnoreCase("roche") || b.equalsIgnoreCase("sbx")){VectorUMP45.platform=VectorUMP45.PLATFORM_ROCHE;}
				else{VectorUMP45.platform=Integer.parseInt(b);}
			}else if(a.equals("normalize") || a.equals("leftalign") || a.equals("norm")){
				normalize=Parse.parseBoolean(b);
			}else if(a.equals("splitalleles")){
				splitAlleles=Parse.parseBoolean(b);
			}else if(a.equals("splitsubs") || a.equals("splitsnps")){
				splitSubs=Parse.parseBoolean(b);
			}else if(a.equals("split") || a.equals("sass")){
				splitAlleles=splitSubs=normalize=Parse.parseBoolean(b);
			}else if(a.equals("histnn") || a.equals("nnhist")){
				histNNFile=b;
			}else if(a.equals("histqual") || a.equals("qualhist")){
				histQualFile=b;
			}else if(a.equals("crossover") || a.equals("xover")){
				crossover=Double.parseDouble(b);
			}else if(a.equals("simdff") || a.equals("simdfeedforward")){
				simdFFEval=Parse.parseBoolean(b);
				shared.Shared.SIMD_FEED_FORWARD=simdFFEval; //Toggle the batched SIMD feed-forward used for net scoring
			}else if(a.equals("maxvcfreadthreads") || a.equals("readthreads") || a.equals("vcfthreads")){
				MAX_VCF_READ_THREADS=Math.max(1, Integer.parseInt(b)); //Tunable worker cap (was hardcoded 8); resolved as min(this, t)
			}else if(a.equals("lp") || a.equals("lineparser")){
				useLineParser=Parse.parseBoolean(b); //SIMD LineParser1 path in VCFLine (perf experiment)
			}else if(a.equals("reusevec")){
				reuseVec=Parse.parseBoolean(b); //Reuse a per-thread feature vector instead of allocating per variant
			}else if(a.equals("clearfilters")){
				if(Parse.parseBoolean(b)){
					varFilter.clear();
				}
			}else if(varFilter.parse(a, b, arg)){
				setVarFilter=true;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(!setVarFilter){varFilter=null;}

		// If NN scoring is requested but no explicit net= was given, choose the default network for this platform+ploidy.
		// (ploidy here is the arg/default value; ##ploidy from the VCF header is applied later and does not affect the
		//  current single-net choice, which covers ploidy 1-2.)
		if(useNet && netFile==null){netFile=NNChooser.choose(VectorUMP45.platform, ploidy);}

		if(netFile!=null && useNet){
			net0=CellNetParser.load(netFile);
			assert(net0!=null) : "Failed to load neural network: "+netFile;
			if(autoCutoff){netCutoff=net0.cutoff;}
			if(verbose){outstream.println("Loaded neural network: "+netFile+" (cutoff="+netCutoff+")");}
		}else{
			net0=null;
		}

		assert(FastaReadInputStream.settingsOK());

		if(in1==null){throw new RuntimeException("Error - in= is required.");}
		if(truthFile==null){throw new RuntimeException("Error - truth= is required.");}
		// Quick fix for the baseless-map corruption: a net recomputes reference-derived features
		// (contigEndDist/homopolymer) from scaffold bases, so ref= is mandatory in NN mode. See the
		// PROPER FIX note in gradeVariant — once CED/HMP are read from INFO this can be removed.
		if(net0!=null && ref==null){
			throw new RuntimeException("Error - ref= is required when net= is used "
					+"(the network recomputes contigEndDist/homopolymer features from reference bases).");
		}

		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && !Tools.testOutputFiles(overwrite, append, false, out1)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}

		if(!Tools.testInputFiles(false, true, in1, truthFile, ref, bedFile)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		ffout1=(out1!=null ? FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, true) : null);

		if(normalize && ref==null){
			throw new RuntimeException("Error - normalize requires ref=");
		}
		if(ref!=null){ScafMap.loadReference(ref, null, null, true);}

		if(bedFile!=null){
			bedMask=new BedMask(bedFile);
			outstream.println("Loaded "+bedMask.intervalsLoaded()+" BED intervals across "+bedMask.scaffolds()+" scaffolds.");
		}
		assert(bedWeight>=0) : "bedweight must be >=0 (got "+bedWeight+").";
		assert(!invertBed || bedFile!=null) : "invertbed requires a bed= file.";

		threads=Tools.min(MAX_VCF_READ_THREADS>0 ? MAX_VCF_READ_THREADS : (net0==null ? 8 : 16), Shared.threads());
	}

	private byte[] refBasesFor(String scaf){
		Scaffold sc=ScafMap.defaultScafMap().getScaffold(scaf);
		return sc==null ? null : sc.bases;
	}

	/** True if v is in the high-confidence region: inside the bed, or OUTSIDE it when invertBed is set.
	 *  With no bed loaded, everything is high-confidence (invert is ignored). */
	private boolean inHighConf(VCFLine v){
		if(bedMask==null){return true;}
		boolean inBed=bedMask.contains(v.scaf, v.pos);
		return inBed!=invertBed;
	}

	/** Concordance weight for v: 1.0 in the high-confidence region, bedWeight outside it. bedWeight=0
	 *  (default) reproduces the old hard bed restriction (off-region truth/calls contribute nothing). */
	private double weightFor(VCFLine v){
		return inHighConf(v) ? 1.0 : bedWeight;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Truth Set Loading      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads a truth VCF into a HashMap of VCFLine objects (key==value), using the same
	 * identity contract (equals/hashCode) as CompareVCF. A HashMap (not a HashSet) is
	 * used so that get(callVline) returns the canonical truth object, which can then be
	 * marked exactly once. Optionally splits multi-allelic variants, left-aligns indels,
	 * and restricts to BED regions.
	 */
	private HashMap<VCFLine,VCFLine> loadTruthSet(String fname){
		HashMap<VCFLine,VCFLine> map=new HashMap<VCFLine,VCFLine>();
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TXT, null, true, false);
		VCFFile vfile=new VCFFile(ff);
		for(Entry<VCFLine, VCFLine> e : vfile.map.entrySet()){
			VCFLine v=e.getValue();
			ArrayList<VCFLine> list=null;
			if(splitAlleles || splitSubs){list=v.split(splitAlleles, false, splitSubs);}
			if(list==null || list.isEmpty()){
				if(normalize){v.leftAlign(refBasesFor(v.scaf));}
				if(weightFor(v)>0){map.put(v, v);}
			}else{
				for(VCFLine line : list){
					if(normalize){line.leftAlign(refBasesFor(line.scaf));}
					if(weightFor(line)>0){map.put(line, line);}
				}
			}
		}
		return map;
	}

	/*--------------------------------------------------------------*/
	/*----------------         VCF Header           ----------------*/
	/*--------------------------------------------------------------*/

	private void processVcfHeader(ByteFile bf, ByteStreamWriter bsw){
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		while(line!=null && (line.length==0 || line[0]=='#')){
			if(line.length>0){
				linesProcessed++;
				bytesProcessed+=line.length;
				bb.append(line).append('\n');
				header.add(line);
				if(Tools.startsWith(line, "##contig=<ID=")){
					scafMap.addFromVcf(line);
				}else{
					String[] sp=new String(line).split("=");
					if(sp.length==2){
						String ka=sp[0], kv=sp[1];
						if(ka.equalsIgnoreCase("##ploidy")){
							ploidy=Integer.parseInt(kv);
						}else if(ka.equalsIgnoreCase("##properPairRate")){
							properPairRate=(float)Double.parseDouble(kv);
						}else if(ka.equalsIgnoreCase("##totalQualityAvg")){
							totalQualityAvg=(float)Double.parseDouble(kv);
						}else if(ka.equalsIgnoreCase("##mapqAvg")){
							totalMapqAvg=(float)Double.parseDouble(kv);
						}else if(ka.equalsIgnoreCase("##readLengthAvg")){
							readLengthAvg=(float)Double.parseDouble(kv);
						}
					}
				}
			}
			line=bf.nextLine();
		}
		if(line!=null && line.length>0){bf.pushBack(line);}
		if(bb.length()>0 && bsw!=null){
			bsw.add(bb, jobIDOffset);
			jobIDOffset++;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       MT Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<ProcessThread> spawnThreads(ByteFile bf, ByteStreamWriter bsw){
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(bf, bsw, jobIDOffset));
		}
		for(ProcessThread pt : alpt){pt.start();}
		return alpt;
	}

	private void waitForFinish(ArrayList<ProcessThread> alpt){
		boolean allSuccess=true;
		for(ProcessThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try{pt.join();}
				catch(InterruptedException e){e.printStackTrace();}
			}
			linesProcessed+=pt.linesProcessedT;
			bytesProcessed+=pt.bytesProcessedT;
			variantLinesProcessed+=pt.variantLinesProcessedT;
			// Merge each thread's marked-truth map into the global map, keeping the
			// maximum matching-call score per truth (each truth counted once).
			for(Entry<VCFLine, float[]> e : pt.matchedT.entrySet()){
				float[] inc=e.getValue();
				float[] cur=matched.get(e.getKey());
				if(cur==null){
					matched.put(e.getKey(), new float[]{inc[0], inc[1]});
				}else{
					if(inc[0]>cur[0]){cur[0]=inc[0];}
					if(inc[1]>cur[1]){cur[1]=inc[1];}
				}
			}
			for(int i=0; i<histNNFalse.length; i++){histNNFalse[i]+=pt.histNNFalseT[i];}
			for(int i=0; i<histQualFalse.length; i++){histQualFalse[i]+=pt.histQualFalseT[i];}
			allSuccess&=pt.success;
		}
		if(!allSuccess){errorState=true;}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Worker Thread         ----------------*/
	/*--------------------------------------------------------------*/

	private class ProcessThread extends Thread {

		ProcessThread(ByteFile bf_, ByteStreamWriter bsw_, long jobIDOffset_){
			bf=bf_;
			bsw=bsw_;
			offset=jobIDOffset_;
		}

		@Override
		public void run(){
			net=(net0==null ? null : net0.copy(false));
			if(net!=null){net.simdFF=simdFFEval;} //Per-eval SIMD feed-forward toggle (also gated by Shared.SIMD && SIMD_FEED_FORWARD)
			ListNum<byte[]> ln=nextListSync();
			while(ln!=null){
				ByteBuilder bb=(bsw!=null ? new ByteBuilder(4096) : null);
				for(byte[] line : ln){
					linesProcessedT++;
					processLine(line, bb);
				}
				if(bsw!=null && bb!=null){bsw.add(bb, ln.id+offset);}
				ln=nextListSync();
			}
			success=true;
			synchronized(this){notify();}
		}

		private ListNum<byte[]> nextListSync(){
			synchronized(bf){return bf.nextList();}
		}

		void processLine(byte[] line, ByteBuilder bb){
			bytesProcessedT+=line.length;
			if(line.length<1){return;}
			if(line[0]=='#'){return;}

			variantLinesProcessedT++;
			VCFLine vline=(useLineParser ? new VCFLine(line, lp) : new VCFLine(line));

			ArrayList<VCFLine> variants=null;
			if(splitAlleles || splitSubs){variants=vline.split(splitAlleles, false, splitSubs);}
			if(variants==null){
				if(normalize){vline.leftAlign(refBasesFor(vline.scaf));}
				gradeVariant(vline, line, bb);
			}else{
				for(VCFLine v2 : variants){
					if(normalize){v2.leftAlign(refBasesFor(v2.scaf));}
					gradeVariant(v2, line, bb);
				}
			}
		}

		void gradeVariant(VCFLine vline, byte[] origLine, ByteBuilder bb){
			final double weight=weightFor(vline);
			if(weight<=0){return;}

			// Optional hard pre-filter: a variant failing the VarFilter is not considered
			// at all (neither marks truth nor counts as FP). With clearfilters varFilter is null.
			if(varFilter!=null){
				boolean ok;
				try{
					// Same reference-feature concern as the scoring block below: in NN mode the
					// vector recompute needs the bases-loaded map, so use defaultScafMap (ref is
					// required when net!=null). Without a net, the composite path uses the instance
					// map (no ref required).
					ScafMap vmap=(net!=null ? ScafMap.defaultScafMap() : scafMap);
					Var v=VcfToVar.fromVCF(origLine, vmap, true, true);
					ok=varFilter.passesFilter(v, properPairRate, totalQualityAvg,
							totalMapqAvg, readLengthAvg, ploidy, vmap, net, false);
				}catch(Throwable e){
					ok=(vline.qual>=varFilter.minScore);
				}
				if(!ok){return;}
			}

			final float qualScore=(float)vline.qual;

			// NN score for the histogram axis (computed for every variant, not gated).
			float nnScore=-1;
			if(net!=null){
				try{
					// makeVector recomputes two reference-derived features from scaffold BASES:
					// contigEndDist (vec[30], capped at nScan) and homopolymerCount (vec[28]). So
					// we must pass the ref-loaded defaultScafMap, NOT the ##contig-derived instance
					// scafMap (which carries lengths but no bases). The danger is not that a baseless
					// map is unusable — it's that it silently returns a DIFFERENT value rather than
					// failing: contigEndDist skips its cap and returns the raw distance (Var.java:1927),
					// homopolymerCount returns 0, and the net then scores a corrupted vector. That is
					// why ref= is mandatory in NN mode (enforced in the constructor).
					//
					// PROPER FIX (would remove the ref requirement entirely): CED= and HMP= are ALREADY
					// written to the VCF INFO by Var.toVCF. fromVCF should parse them into cached Var
					// fields, and makeVector should prefer the cached value when present, recomputing
					// from the map only on the live calling path (where no INFO yet exists). Then
					// post-hoc scoring would need no reference at all and could never silently corrupt.
					ScafMap vmap=ScafMap.defaultScafMap();
					Var v=VcfToVar.fromVCF(origLine, vmap, true, true);
					float[] vec=(reuseVec ? VectorUMP45.makeVector(v, properPairRate, totalQualityAvg,
							totalMapqAvg, readLengthAvg, ploidy, vmap, vecBuf)
						: VectorUMP45.makeVector(v, properPairRate, totalQualityAvg,
							totalMapqAvg, readLengthAvg, ploidy, vmap));
					net.applyInput(vec);
					nnScore=net.feedForward();
				}catch(Throwable e){
					nnScore=-1;
				}
			}

			final VCFLine truth=truthMap.get(vline);
			if(truth!=null){
				// Mark the truth var, keeping the max matching-call score (NN and QUAL).
				float[] cur=matchedT.get(truth);
				if(cur==null){
					matchedT.put(truth, new float[]{nnScore, qualScore});
				}else{
					if(nnScore>cur[0]){cur[0]=nnScore;}
					if(qualScore>cur[1]){cur[1]=qualScore;}
				}
			}else{
				// FP: a call that matches no truth. Counted directly, per call, at its score.
				if(net!=null && nnScore>=0){histNNFalseT[Tools.mid(0, (int)(nnScore*NN_SCALE), NN_BINS-1)]+=weight;}
				histQualFalseT[Tools.mid(0, (int)(qualScore*QUAL_SCALE), QUAL_BINS-1)]+=weight;
			}

			// Optional passing-call output VCF, gated by the configured operating point.
			final double opScore=(net!=null ? nnScore : qualScore);
			final double opThresh=(net!=null ? netCutoff : minScore);
			if(opScore>=opThresh && bsw!=null && bb!=null){
				bb.append(origLine).append('\n');
			}
		}

		final ByteFile bf;
		final ByteStreamWriter bsw;
		CellNet net;
		final LineParser1 lp=new LineParser1((byte)'\t'); //Reusable SIMD field-splitter (used when useLineParser)
		final float[] vecBuf=new float[VectorUMP45.DIMS]; //Reusable feature vector (used when reuseVec)
		final long offset;

		long linesProcessedT=0;
		long bytesProcessedT=0;
		long variantLinesProcessedT=0;
		// Distinct truth vars this thread matched, mapped to {maxNN, maxQual} matching scores.
		HashMap<VCFLine, float[]> matchedT=new HashMap<VCFLine, float[]>();
		double[] histNNFalseT=new double[NN_BINS];
		double[] histQualFalseT=new double[QUAL_BINS];
		boolean success=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Core Processing        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){

		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		ByteFile bf=ByteFile.makeByteFile(ffin1);
		processVcfHeader(bf, bsw);

		if(scafMap.size()==0){
			throw new RuntimeException("ScafMap is empty — VCF has no ##contig lines and no ref= given.");
		}

		if(ScafMap.defaultScafMap()==null){
			ScafMap.setDefaultScafMap(scafMap, in1);
		}

		outstream.println("Loading truth set: "+truthFile);
		truthMap=loadTruthSet(truthFile);
		final long truthSize=truthMap.size();
		outstream.println("Truth variants: "+truthSize);

		ArrayList<ProcessThread> alpt=spawnThreads(bf, bsw);
		waitForFinish(alpt);

		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		// Build the truth (true-positive) histograms by binning each distinct marked truth var ONCE,
		// at the maximum score of any call that matched it, weighted by its bed weight.
		double weightedMarked=0;
		for(Entry<VCFLine, float[]> e : matched.entrySet()){
			final double w=weightFor(e.getKey());
			float[] sc=e.getValue();
			float maxNN=sc[0], maxQual=sc[1];
			if(net0!=null && maxNN>=0){histNNTruth[Tools.mid(0, (int)(maxNN*NN_SCALE), NN_BINS-1)]+=w;}
			histQualTruth[Tools.mid(0, (int)(maxQual*QUAL_SCALE), QUAL_BINS-1)]+=w;
			weightedMarked+=w;
		}

		// Weighted counts (marking contract): TP = weighted marked truth, FN = weighted unmarked truth.
		// bedWeight=0 (default) makes every truth weight 1.0 (off-region truth was never added to the map),
		// so these reduce exactly to the old integer counts.
		double weightedTruthSize=0;
		for(VCFLine tline : truthMap.keySet()){weightedTruthSize+=weightFor(tline);}
		final double callerMissed=weightedTruthSize-weightedMarked;

		// Operating-point summary at the configured threshold (NN cutoff, or minScore on QUAL).
		final boolean nnMode=(net0!=null);
		final double opThresh=(nnMode ? netCutoff : minScore);
		double tpOp=0;
		for(Entry<VCFLine, float[]> e : matched.entrySet()){
			float[] sc=e.getValue();
			double s=(nnMode ? sc[0] : sc[1]);
			if(s>=opThresh){tpOp+=weightFor(e.getKey());}
		}
		double fpOp=0;
		{
			double[] falseHist=(nnMode ? histNNFalse : histQualFalse);
			double binScale=(nnMode ? NN_SCALE : QUAL_SCALE);
			int bins=(nnMode ? NN_BINS : QUAL_BINS);
			for(int bin=0; bin<bins; bin++){
				if(bin/binScale>=opThresh){fpOp+=falseHist[bin];}
			}
		}
		final double fnOp=weightedTruthSize-tpOp;

		double precision=(tpOp+fpOp>0 ? tpOp/(tpOp+fpOp) : 0);
		double recall=(weightedTruthSize>0 ? tpOp/weightedTruthSize : 0);
		double f1=(precision+recall>0 ? 2*precision*recall/(precision+recall) : 0);
		double fpr=(tpOp+fpOp>0 ? fpOp/(tpOp+fpOp) : 0);
		double fnr=(weightedTruthSize>0 ? fnOp/weightedTruthSize : 0);

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Variant Lines In:   \t"+variantLinesProcessed);
		outstream.println("Truth Variants:     \t"+truthSize);
		outstream.printf("Caller Missed (FN): \t%.2f%n", callerMissed);
		outstream.println("Operating point:    \t"+(nnMode ? "NN>=" : "QUAL>=")+opThresh);
		outstream.println();
		outstream.printf("TP:        \t%.2f%n", tpOp);
		outstream.printf("FP:        \t%.2f%n", fpOp);
		outstream.printf("FN:        \t%.2f%n", fnOp);
		outstream.println();
		outstream.printf("Precision: \t%.6f%n", precision);
		outstream.printf("Recall:    \t%.6f%n", recall);
		outstream.printf("F1:        \t%.6f%n", f1);
		outstream.printf("FPR:       \t%.6f%n", fpr);
		outstream.printf("FNR:       \t%.6f%n", fnr);

		// FN=crossover*FP operating-point report (e.g. crossover=4 -> the FN=4*FP cutoff). Lets each eval job
		// read the optimal cutoff + FP+FN straight from GradeVCF instead of post-processing the histogram.
		if(crossover>0){
			outstream.println();
			double[] qx=findCrossover(histQualTruth, histQualFalse, QUAL_BINS, QUAL_SCALE, crossover, weightedTruthSize);
			outstream.printf("CROSSOVER_QUAL\tFN=%.1fxFP\tcutoff=%.4f\tFP=%.2f\tFN=%.2f\tFP+FN=%.2f%n",
					crossover, qx[0], qx[1], qx[2], qx[1]+qx[2]);
			if(nnMode){
				double[] nx=findCrossover(histNNTruth, histNNFalse, NN_BINS, NN_SCALE, crossover, weightedTruthSize);
				outstream.printf("CROSSOVER_NN\tFN=%.1fxFP\tcutoff=%.4f\tFP=%.2f\tFN=%.2f\tFP+FN=%.2f%n",
						crossover, nx[0], nx[1], nx[2], nx[1]+nx[2]);
			}
		}

		if(histNNFile!=null){writeHist(histNNFile, histNNTruth, histNNFalse, NN_BINS, NN_SCALE, "NN_score", callerMissed);}
		if(histQualFile!=null){writeHist(histQualFile, histQualTruth, histQualFalse, QUAL_BINS, QUAL_SCALE, "QUAL", callerMissed);}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Histogram Writers       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Writes a two-section histogram: raw per-bin true/false counts, then the cumulative
	 * TP/FP/FN curve. trueCounts holds each distinct truth var counted ONCE at the bin of
	 * its maximum matching-call score; falseCounts holds unmatched calls per bin.
	 *
	 * Cumulative FP sums high->low (FP = false calls scoring >= threshold). Cumulative FN
	 * sums low->high starting at callerMissed, which is identical to (truthSize - cumTP):
	 * every truth var is either matched at score >= threshold (TP) or not (FN).
	 */
	private void writeHist(String fname, double[] trueCounts, double[] falseCounts,
			int bins, double binScale, String scoreLabel, double callerMissed){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();

		// Raw per-bin counts
		bb.append("#").append(scoreLabel).append("\ttrueCounts\tfalseCounts\n");
		for(int bin=0; bin<bins; bin++){
			double score=bin/binScale;
			bb.append(score, 2).tab().append(trueCounts[bin], 2).tab().append(falseCounts[bin], 2).nl();
		}
		bb.nl();

		// Cumulative: FP sums high->low, FN sums low->high (starting at callerMissed)
		bb.append("#threshold\tTP\tFP\tFN\n");
		double[] cumFNArr=new double[bins];
		double cumFN=callerMissed;
		for(int bin=0; bin<bins; bin++){
			cumFNArr[bin]=cumFN;
			cumFN+=trueCounts[bin];
		}
		double cumTP=0, cumFP=0;
		for(int bin=bins-1; bin>=0; bin--){
			cumTP+=trueCounts[bin];
			cumFP+=falseCounts[bin];
			double threshold=bin/binScale;
			bb.append(threshold, 2).tab().append(cumTP, 2).tab().append(cumFP, 2).tab().append(cumFNArr[bin], 2).nl();
		}

		bsw.print(bb);
		errorState|=bsw.poisonAndWait();
	}

	/**
	 * Finds the operating point on one score axis where FN ~= crossover*FP, scanning the cumulative
	 * histogram high->low. Returns {cutoff, FP, FN} at the bin minimizing |FN - crossover*FP|.
	 * (crossover=4 -> the FN=4*FP point: four false-negatives tolerated per false-positive.)
	 */
	private static double[] findCrossover(double[] trueCounts, double[] falseCounts, int bins,
			double binScale, double crossover, double weightedTruthSize){
		double cumTP=0, cumFP=0;
		double bestCutoff=0, bestFP=0, bestFN=weightedTruthSize, bestDiff=Double.MAX_VALUE;
		for(int bin=bins-1; bin>=0; bin--){
			cumTP+=trueCounts[bin];
			cumFP+=falseCounts[bin];
			final double fp=cumFP;
			final double fn=weightedTruthSize-cumTP;
			final double diff=Math.abs(fn-crossover*fp);
			if(diff<bestDiff){bestDiff=diff; bestCutoff=bin/binScale; bestFP=fp; bestFN=fn;}
		}
		return new double[]{bestCutoff, bestFP, bestFN};
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	// Histogram resolution. NN axis: bin=(int)(nn*NN_SCALE), step 1/NN_SCALE=0.002, range 0..2.0 -> 1001 bins
	// (the net output exceeds 1.0 -- graded labels reach 1.45 -- so the FN=k*FP cutoff can lie above 1.0; the
	// extra headroom keeps high-scoring FPs separable from even-higher TPs).
	// QUAL axis: bin=(int)(qual*QUAL_SCALE), step 1/QUAL_SCALE=0.1, range 0..60 -> 601 bins.
	private static final double NN_SCALE=500.0;
	private static final double QUAL_SCALE=10.0;
	private static final int NN_BINS=1001;
	private static final int QUAL_BINS=601;

	private long linesProcessed=0;
	private long bytesProcessed=0;
	private long variantLinesProcessed=0;

	private HashMap<VCFLine,VCFLine> truthMap=null;
	// Distinct truth vars matched by any call, mapped to {maxNN, maxQual} matching scores.
	private final HashMap<VCFLine, float[]> matched=new HashMap<VCFLine, float[]>();

	// trueCounts: distinct truth binned once at max matching score (built after merge).
	private double[] histNNTruth=new double[NN_BINS];
	private double[] histQualTruth=new double[QUAL_BINS];
	// falseCounts: unmatched calls binned per call (summed across threads). Weighted by bedWeight.
	private double[] histNNFalse=new double[NN_BINS];
	private double[] histQualFalse=new double[QUAL_BINS];

	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	private long jobIDOffset=0;

	private CellNet net0=null;
	private boolean useNet=false;
	private float netCutoff=0.5f;
	/** Eval-time SIMD feed-forward toggle (simdff=t/f): stamped onto each worker's scoring net so we can
	 * measure a net's SIMD_FF-on vs SIMD_FF-off result. Also flips global Shared.SIMD_FEED_FORWARD. */
	boolean simdFFEval=false;
	/** lp=t/f: use the SIMD-accelerated LineParser1 path in VCFLine instead of the scalar tab-scan. */
	boolean useLineParser=false;
	/** reusevec=t/f: reuse a per-thread feature vector (no per-variant float[] allocation). Default on. */
	boolean reuseVec=true;

	double minScore=0;
	// FN=crossover*FP operating point. 0 (default) = off. When >0, GradeVCF reports the QUAL and NN
	// cutoffs at that point and the FP+FN there, so each eval reads the answer instead of post-processing.
	double crossover=0;
	public int ploidy=1;
	public float properPairRate=0;
	public float totalQualityAvg=30;
	public float totalMapqAvg=30;
	public float readLengthAvg=150;

	VarFilter varFilter=new VarFilter();
	BedMask bedMask=null;
	// Concordance weight for calls/truth OUTSIDE the high-conf region (in-region weight is always 1.0).
	// 0 (default) = old hard bed restriction; 0.5 = off-region counts half; 1.0 = no restriction.
	double bedWeight=0.0;
	// When true, the high-conf region is the COMPLEMENT of the bed (grade ONLY low-confidence regions).
	boolean invertBed=false;
	boolean normalize=false;
	boolean splitAlleles=false;
	boolean splitSubs=false;

	/** Override cap on concurrent VCF-processing worker threads (runtime flag maxvcfreadthreads=, aliases
	 *  readthreads=, vcfthreads=). 0 = auto: 8 with no network (read-bound; the shared path caps ~9 cores),
	 *  16 with a network (NN scoring parallelizes). Worker count = min(resolved cap, t). */
	public static int MAX_VCF_READ_THREADS=0;

	final int threads;

	private String in1=null;
	private String out1=null;
	private String truthFile=null;
	private String ref=null;
	private String bedFile=null;
	private String histNNFile=null;
	private String histQualFile=null;

	private final FileFormat ffin1;
	private final FileFormat ffout1;

	public final ScafMap scafMap=new ScafMap();

	static final ListNum<byte[]> POISON_BYTES=new ListNum<byte[]>(null, Long.MAX_VALUE, true, false);
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
}
