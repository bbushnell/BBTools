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

		threads=Tools.min(8, Shared.threads());
	}

	private byte[] refBasesFor(String scaf){
		Scaffold sc=ScafMap.defaultScafMap().getScaffold(scaf);
		return sc==null ? null : sc.bases;
	}

	private boolean inRegion(VCFLine v){
		return bedMask==null || bedMask.contains(v.scaf, v.pos);
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
				if(inRegion(v)){map.put(v, v);}
			}else{
				for(VCFLine line : list){
					if(normalize){line.leftAlign(refBasesFor(line.scaf));}
					if(inRegion(line)){map.put(line, line);}
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
			Tools.add(histNNFalse, pt.histNNFalseT);
			Tools.add(histQualFalse, pt.histQualFalseT);
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
			VCFLine vline=new VCFLine(line);

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
			if(!inRegion(vline)){return;}

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
					float[] vec=VectorUMP45.makeVector(v, properPairRate, totalQualityAvg,
							totalMapqAvg, readLengthAvg, ploidy, vmap);
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
				if(net!=null && nnScore>=0){histNNFalseT[Tools.mid(0, (int)(nnScore*100), NN_BINS-1)]++;}
				histQualFalseT[Tools.mid(0, (int)(qualScore*4), QUAL_BINS-1)]++;
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
		final long offset;

		long linesProcessedT=0;
		long bytesProcessedT=0;
		long variantLinesProcessedT=0;
		// Distinct truth vars this thread matched, mapped to {maxNN, maxQual} matching scores.
		HashMap<VCFLine, float[]> matchedT=new HashMap<VCFLine, float[]>();
		long[] histNNFalseT=new long[NN_BINS];
		long[] histQualFalseT=new long[QUAL_BINS];
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

		// Build the truth (true-positive) histograms by binning each distinct marked
		// truth var ONCE, at the maximum score of any call that matched it.
		for(float[] sc : matched.values()){
			float maxNN=sc[0], maxQual=sc[1];
			if(net0!=null && maxNN>=0){histNNTruth[Tools.mid(0, (int)(maxNN*100), NN_BINS-1)]++;}
			histQualTruth[Tools.mid(0, (int)(maxQual*4), QUAL_BINS-1)]++;
		}

		// Direct counts (marking contract): TP = marked truth, FN = unmarked truth.
		final long markedTruth=matched.size();
		final long callerMissed=truthSize-markedTruth;

		// Operating-point summary at the configured threshold (NN cutoff, or minScore on QUAL).
		final boolean nnMode=(net0!=null);
		final double opThresh=(nnMode ? netCutoff : minScore);
		long tpOp=0;
		for(float[] sc : matched.values()){
			double s=(nnMode ? sc[0] : sc[1]);
			if(s>=opThresh){tpOp++;}
		}
		long fpOp=0;
		{
			long[] falseHist=(nnMode ? histNNFalse : histQualFalse);
			double binScale=(nnMode ? 100.0 : 4.0);
			int bins=(nnMode ? NN_BINS : QUAL_BINS);
			for(int bin=0; bin<bins; bin++){
				if(bin/binScale>=opThresh){fpOp+=falseHist[bin];}
			}
		}
		final long fnOp=truthSize-tpOp;

		double precision=(tpOp+fpOp>0 ? (double)tpOp/(tpOp+fpOp) : 0);
		double recall=(truthSize>0 ? (double)tpOp/truthSize : 0);
		double f1=(precision+recall>0 ? 2*precision*recall/(precision+recall) : 0);
		double fpr=(tpOp+fpOp>0 ? (double)fpOp/(tpOp+fpOp) : 0);
		double fnr=(truthSize>0 ? (double)fnOp/truthSize : 0);

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Variant Lines In:   \t"+variantLinesProcessed);
		outstream.println("Truth Variants:     \t"+truthSize);
		outstream.println("Caller Missed (FN): \t"+callerMissed);
		outstream.println("Operating point:    \t"+(nnMode ? "NN>=" : "QUAL>=")+opThresh);
		outstream.println();
		outstream.println("TP:        \t"+tpOp);
		outstream.println("FP:        \t"+fpOp);
		outstream.println("FN:        \t"+fnOp);
		outstream.println();
		outstream.printf("Precision: \t%.6f%n", precision);
		outstream.printf("Recall:    \t%.6f%n", recall);
		outstream.printf("F1:        \t%.6f%n", f1);
		outstream.printf("FPR:       \t%.6f%n", fpr);
		outstream.printf("FNR:       \t%.6f%n", fnr);

		if(histNNFile!=null){writeHist(histNNFile, histNNTruth, histNNFalse, NN_BINS, 100.0, "NN_score", callerMissed);}
		if(histQualFile!=null){writeHist(histQualFile, histQualTruth, histQualFalse, QUAL_BINS, 4.0, "QUAL", callerMissed);}

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
	private void writeHist(String fname, long[] trueCounts, long[] falseCounts,
			int bins, double binScale, String scoreLabel, long callerMissed){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();

		// Raw per-bin counts
		bb.append("#").append(scoreLabel).append("\ttrueCounts\tfalseCounts\n");
		for(int bin=0; bin<bins; bin++){
			double score=bin/binScale;
			bb.append(score, 2).tab().append(trueCounts[bin]).tab().append(falseCounts[bin]).nl();
		}
		bb.nl();

		// Cumulative: FP sums high->low, FN sums low->high (starting at callerMissed)
		bb.append("#threshold\tTP\tFP\tFN\n");
		long[] cumFNArr=new long[bins];
		long cumFN=callerMissed;
		for(int bin=0; bin<bins; bin++){
			cumFNArr[bin]=cumFN;
			cumFN+=trueCounts[bin];
		}
		long cumTP=0, cumFP=0;
		for(int bin=bins-1; bin>=0; bin--){
			cumTP+=trueCounts[bin];
			cumFP+=falseCounts[bin];
			double threshold=bin/binScale;
			bb.append(threshold, 2).tab().append(cumTP).tab().append(cumFP).tab().append(cumFNArr[bin]).nl();
		}

		bsw.print(bb);
		errorState|=bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	private static final int NN_BINS=101;
	private static final int QUAL_BINS=241;

	private long linesProcessed=0;
	private long bytesProcessed=0;
	private long variantLinesProcessed=0;

	private HashMap<VCFLine,VCFLine> truthMap=null;
	// Distinct truth vars matched by any call, mapped to {maxNN, maxQual} matching scores.
	private final HashMap<VCFLine, float[]> matched=new HashMap<VCFLine, float[]>();

	// trueCounts: distinct truth binned once at max matching score (built after merge).
	private long[] histNNTruth=new long[NN_BINS];
	private long[] histQualTruth=new long[QUAL_BINS];
	// falseCounts: unmatched calls binned per call (summed across threads).
	private long[] histNNFalse=new long[NN_BINS];
	private long[] histQualFalse=new long[QUAL_BINS];

	public ArrayList<byte[]> header=new ArrayList<byte[]>();
	private long jobIDOffset=0;

	private CellNet net0=null;
	private boolean useNet=false;
	private float netCutoff=0.5f;

	double minScore=0;
	public int ploidy=1;
	public float properPairRate=0;
	public float totalQualityAvg=30;
	public float totalMapqAvg=30;
	public float readLengthAvg=150;

	VarFilter varFilter=new VarFilter();
	BedMask bedMask=null;
	boolean normalize=false;
	boolean splitAlleles=false;
	boolean splitSubs=false;

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
