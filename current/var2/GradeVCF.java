package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
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
 * Computes TP, FP, FN, TN counts with optional BED region restriction, VarFilter
 * pre-filtering, and neural network scoring. Writes cumulative QUAL and NN score
 * histograms suitable for precision/recall curve plotting.
 *
 * Variant matching uses VCFLine equals/hashCode — the same identity contract as
 * CompareVCF — so results agree exactly with comparevcf intersection counts.
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
	 * Loads a truth VCF into a HashSet of VCFLine objects, using the same identity
	 * contract (equals/hashCode) as CompareVCF. Optionally splits multi-allelic
	 * variants, left-aligns indels, and restricts to BED regions.
	 */
	private HashSet<VCFLine> loadTruthSet(String fname){
		HashSet<VCFLine> set=new HashSet<VCFLine>();
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TXT, null, true, false);
		VCFFile vfile=new VCFFile(ff);
		for(Entry<VCFLine, VCFLine> e : vfile.map.entrySet()){
			VCFLine v=e.getValue();
			ArrayList<VCFLine> list=null;
			if(splitAlleles || splitSubs){list=v.split(splitAlleles, false, splitSubs);}
			if(list==null || list.isEmpty()){
				if(normalize){v.leftAlign(refBasesFor(v.scaf));}
				if(inRegion(v)){set.add(v);}
			}else{
				for(VCFLine line : list){
					if(normalize){line.leftAlign(refBasesFor(line.scaf));}
					if(inRegion(line)){set.add(line);}
				}
			}
		}
		return set;
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
			tp+=pt.tpT;
			fp+=pt.fpT;
			fn+=pt.fnT;
			tn+=pt.tnT;
			synchronized(seenTruth){
				seenTruth.addAll(pt.seenTruthT);
			}
			Tools.add(histNN, pt.histNNT);
			Tools.add(histNNTruth, pt.histNNTruthT);
			Tools.add(histQual, pt.histQualT);
			Tools.add(histQualTruth, pt.histQualTruthT);
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

			boolean isTruth=truthSet.contains(vline);
			if(isTruth){seenTruthT.add(vline);}

			float qualScore=(float)vline.qual;

			boolean passFilter=true;
			if(varFilter!=null){
				try{
					Var v=VcfToVar.fromVCF(origLine, scafMap, true, true);
					passFilter=varFilter.passesFilter(v, properPairRate, totalQualityAvg,
							totalMapqAvg, readLengthAvg, ploidy, scafMap, net, false);
				}catch(Throwable e){
					passFilter=(qualScore>=varFilter.minScore);
				}
			}

			if(passFilter && qualScore<minScore){passFilter=false;}

			float nnScore=-1;
			if(passFilter && net!=null){
				try{
					Var v=VcfToVar.fromVCF(origLine, scafMap, true, true);
					float[] vec=VectorUMP45.makeVector(v, properPairRate, totalQualityAvg,
							totalMapqAvg, readLengthAvg, ploidy, scafMap);
					net.applyInput(vec);
					nnScore=net.feedForward();
					if(nnScore<netCutoff){passFilter=false;}
				}catch(Throwable e){
					passFilter=false;
				}
			}

			if(passFilter){
				if(isTruth){tpT++;}
				else{fpT++;}
			}else{
				if(isTruth){fnT++;}
				else{tnT++;}
			}

			if(nnScore>=0){
				int nnBin=Tools.mid(0, (int)(nnScore*100), NN_BINS-1);
				histNNT[nnBin]++;
				if(isTruth){histNNTruthT[nnBin]++;}
			}
			{
				int qBin=Tools.mid(0, (int)(qualScore*4), QUAL_BINS-1);
				histQualT[qBin]++;
				if(isTruth){histQualTruthT[qBin]++;}
			}

			if(passFilter && bsw!=null && bb!=null){
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
		long tpT=0, fpT=0, fnT=0, tnT=0;
		HashSet<VCFLine> seenTruthT=new HashSet<VCFLine>();
		long[] histNNT=new long[NN_BINS];
		long[] histNNTruthT=new long[NN_BINS];
		long[] histQualT=new long[QUAL_BINS];
		long[] histQualTruthT=new long[QUAL_BINS];
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
		truthSet=loadTruthSet(truthFile);
		outstream.println("Truth variants: "+truthSet.size());

		ArrayList<ProcessThread> alpt=spawnThreads(bf, bsw);
		waitForFinish(alpt);

		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		long uniqueTP=seenTruth.size();
		long callerMissed=truthSet.size()-uniqueTP;
		fn+=callerMissed;
		tp=uniqueTP;

		double precision=(tp+fp>0 ? (double)tp/(tp+fp) : 0);
		double recall=(tp+fn>0 ? (double)tp/(tp+fn) : 0);
		double f1=(precision+recall>0 ? 2*precision*recall/(precision+recall) : 0);
		double fpr=(tp+fp>0 ? (double)fp/(tp+fp) : 0);
		double fnr=(tp+fn>0 ? (double)fn/(tp+fn) : 0);

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println();
		outstream.println("Variant Lines In:   \t"+variantLinesProcessed);
		outstream.println("Truth Variants:     \t"+truthSet.size());
		outstream.println("Caller Missed (FN): \t"+callerMissed);
		outstream.println();
		outstream.println("TP:        \t"+tp);
		outstream.println("FP:        \t"+fp);
		outstream.println("FN:        \t"+fn);
		outstream.println("TN:        \t"+tn);
		outstream.println();
		outstream.printf("Precision: \t%.6f%n", precision);
		outstream.printf("Recall:    \t%.6f%n", recall);
		outstream.printf("F1:        \t%.6f%n", f1);
		outstream.printf("FPR:       \t%.6f%n", fpr);
		outstream.printf("FNR:       \t%.6f%n", fnr);

		if(histNNFile!=null){writeHistNN(histNNFile);}
		if(histQualFile!=null){writeHistQual(histQualFile);}

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Histogram Writers       ----------------*/
	/*--------------------------------------------------------------*/

	private void writeHistNN(String fname){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();
		bb.append("#threshold\tTP\tFP\tFN\n");
		long cumTP=0, cumFP=0, cumFN=fn;
		for(int bin=NN_BINS-1; bin>=0; bin--){
			long binTotal=histNN[bin];
			long binTruth=histNNTruth[bin];
			long binNonTruth=binTotal-binTruth;
			cumTP+=binTruth;
			cumFP+=binNonTruth;
			cumFN=Tools.max(0, cumFN-binTruth);
			double threshold=bin*0.01;
			bb.append(threshold, 2).tab().append(cumTP).tab().append(cumFP).tab().append(cumFN).nl();
		}
		bsw.print(bb);
		errorState|=bsw.poisonAndWait();
	}

	private void writeHistQual(String fname){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();
		bb.append("#threshold\tTP\tFP\tFN\n");
		long cumTP=0, cumFP=0, cumFN=fn;
		for(int bin=QUAL_BINS-1; bin>=0; bin--){
			long binTotal=histQual[bin];
			long binTruth=histQualTruth[bin];
			long binNonTruth=binTotal-binTruth;
			cumTP+=binTruth;
			cumFP+=binNonTruth;
			cumFN=Tools.max(0, cumFN-binTruth);
			double threshold=bin*0.25;
			bb.append(threshold, 2).tab().append(cumTP).tab().append(cumFP).tab().append(cumFN).nl();
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

	private long tp=0, fp=0, fn=0, tn=0;

	private HashSet<VCFLine> truthSet=null;
	private final HashSet<VCFLine> seenTruth=new HashSet<VCFLine>();

	private long[] histNN=new long[NN_BINS];
	private long[] histNNTruth=new long[NN_BINS];
	private long[] histQual=new long[QUAL_BINS];
	private long[] histQualTruth=new long[QUAL_BINS];

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
