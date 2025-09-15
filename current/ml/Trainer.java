package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Trains a neural network.
 * 
 * @author Brian Bushnell
 * @date February 6, 2023
 *
 */
public class Trainer implements Accumulator<WorkerThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Trainer x=new Trainer(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Trainer(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser.printExecuting=false;
			Parser.printSetThreads=false;
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
			
			commandLine="#CL "+Shared.fullCommandline(false, true);
			commands=new ArrayList<String>(1);
			commands.add(commandLine);
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;

			if(netIn==null) {netIn=parser.in1;}
			if(dataIn==null) {dataIn=parser.in2;}
//			if(netOutFinal==null) {netOutFinal=parser.out1;}
		}
		if(netOutBest==null && netOutFinal!=null) {
			netOutBest=netOutFinal;//Keeps intermediate results in case of early termination.
		}
		
		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffDataOut=FileFormat.testOutput(dataOut, FileFormat.TXT, null, true, overwrite, false, false);
		ffNetOutFinal=FileFormat.testOutput(netOutFinal, FileFormat.TXT, null, true, true, false, false);
		ffNetOutBest=FileFormat.testOutput(netOutBest, FileFormat.TXT, null, true, true, false, false);
		if(ffNetOutBest!=null && ffNetOutBest.exists() && !overwrite) {
			throw new RuntimeException(ffNetOutBest.name()+" exists and overwrite = "+overwrite);
		}
//		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffNetIn=FileFormat.testInput(netIn, FileFormat.TXT, null, true, true);
		ffDataIn=FileFormat.testInput(dataIn, FileFormat.TXT, null, true, true);
		ffValidateIn=FileFormat.testInput(validateIn, FileFormat.TXT, null, true, true);
		
//		masterQueue=new ArrayBlockingQueue<JobResults>(threads);
//		workerQueue=new ArrayBlockingQueue<JobData>(threads);

		alphaIncrease=alphaZero*(alphaMult-1.0)/(peakAlphaEpoch);
		alphaEpochs=(Tools.min(maxEpochs, 800000)-peakAlphaEpoch);
		alphaDropoff=1.0/Math.pow(alphaMult2*alphaMult, 1.0/alphaEpochs);

		randyNetSeed=(netSeed0<0 ? new Random() : new Random(netSeed0));
//		randyAnnealSeed=(annealSeed0<0 ? new Random() : new Random(netSeed0));
		
		//Can be 0.9, for example
		maxAnnealEpoch=(maxAnnealEpoch<Integer.MAX_VALUE ? maxAnnealEpoch : (int)(maxAnnealEpochMult*maxEpochs));
		fpeStart=Tools.min(maxEpochs, fpeStart);

		if(targetFPR>=0) {CellNet.compareCode=CellNet.compareFNR;}
		else if(targetFNR>=0) {CellNet.compareCode=CellNet.compareFPR;}
		else if(crossoverFpMult>0) {CellNet.compareCode=CellNet.compareFPR;}
		else{CellNet.compareCode=CellNet.compareWER;}
		
		printFPR=targetFPR<0 || forcePrintFPR;
		printFNR=targetFNR<0 || forcePrintFNR;
//		assert(false) : printFPR+", "+targetFPR+", "+forcePrintFPR;
		int[] dims=(dims0==null ? minDims : dims0);
		if(dims!=null) {
			numLayers=dims.length;
			numInputs=dims[0];
			numOutputs=dims[numLayers-1];
		}
		if(dumpEpoch<=0 && dumpRate>0){
			dumpEpoch=(int)(dumpEpochMult*maxEpochs);
		}

		if(maxLines<0){maxLines=Shared.MAX_ARRAY_LEN;}
		if(maxLinesV<0){maxLinesV=maxLines;}
		
		if(startTriage<0 && startTriageMult>=0){
			startTriage=(int)(startTriageMult*maxEpochs);
		}
		
		if(minWeightEpoch<0 && minWeightEpochMult>=0){
			minWeightEpoch=(int)(minWeightEpochMult*maxEpochs);
		}
		
		Function.normalizeTypeRates();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
//		parser.netOut="stdout";
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			//TODO: Add train vs test modes
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("train") || a.equals("training")){
				training=Parse.parseBoolean(b);
			}else if(a.equals("evaluate") || a.equals("test")){
				training=!Parse.parseBoolean(b);
			}else if(a.equals("numnetworks") || a.equals("numnets") || a.equals("nets") || a.equals("networks")){
				networksPerCycle=Parse.parseIntKMG(b);
			}else if(a.equals("cycles")){
				cycles=Parse.parseIntKMG(b);
			}else if(a.equals("scan") || a.equals("seedstoscan")){
				seedsToScan=Parse.parseIntKMG(b);
			}else if(a.equals("scanepochs") || a.equals("scanbatches")){
				scanEpochs=Parse.parseIntKMG(b);
			}else if(a.equals("scansamples")){
				scanSamples=Parse.parseIntKMG(b);
			}else if(a.equals("scanmult")){
				seedsToScanMult=Float.parseFloat(b);
//			}else if(a.equals("specialfma")){
//				CellNet.SPECIAL_FMA=Parse.parseBoolean(b);
			}else if(a.equals("balance") || a.equals("balanced")){
				if(b==null || Tools.startsWithLetter(b)) {
					balance=Parse.parseBoolean(b) ? 1.0f : 0;
				}else{
					balance=Float.parseFloat(b);
				}
			}else if(a.equals("normalizebias")){
				CellNet.NORMALIZE_BIAS=Parse.parseBoolean(b);
			}else if(a.equals("normalizeweights")){
				CellNet.NORMALIZE_WEIGHTS=Parse.parseBoolean(b);
			}else if(a.equals("normalize")){
				CellNet.NORMALIZE_BIAS=CellNet.NORMALIZE_WEIGHTS=Parse.parseBoolean(b);
			}else if(a.equals("concise")){
				CellNet.CONCISE=Parse.parseBoolean(b);
			}else if(a.equals("dense")){
				CellNet.DENSE=Parse.parseBoolean(b);
			}else if(a.equals("sparse")){
				CellNet.DENSE=!Parse.parseBoolean(b);
			}else if(a.equals("outsparse")){
				CellNet.OUT_SPARSE=Parse.parseBoolean(b);
				if(CellNet.OUT_SPARSE) {CellNet.OUT_DENSE=false;}
			}else if(a.equals("outdense")){
				CellNet.OUT_DENSE=Parse.parseBoolean(b);
				if(CellNet.OUT_DENSE) {CellNet.OUT_SPARSE=false;}
			}else if(a.equals("outhex")){
				CellNet.OUT_HEX=Parse.parseBoolean(b);
			}else if(a.equals("print") || a.equals("printstatus")){
				printStatus=Parse.parseBoolean(b);
			}else if(a.equals("machineout")){
				machineOut=Parse.parseBoolean(b);
			}else if(a.equals("printall")){
				printMode=(Parse.parseBoolean(b) ? printAll : printBest);
			}else if(a.equals("printbest")){
				printMode=(Parse.parseBoolean(b) ? printBest : printBest);
			}else if(a.equals("printfirst")){
				printMode=(Parse.parseBoolean(b) ? printFirst : printBest);
			}else if(a.equals("printaverage") || a.equals("printavg")){
				printMode=(Parse.parseBoolean(b) ? printAverage : printBest);
			}else if(a.equals("printinterval") || a.equals("printerval")){
				printInterval=Parse.parseIntKMG(b);
			}else if(a.equals("ignorebadlines")){
				DataLoader.IGNORE_BAD_LINES=Parse.parseBoolean(b);
			}else if(a.equals("lines") || a.equals("maxlines") || a.equals("samples") || a.equals("maxsamples")){
				maxLines=Parse.parseIntKMG(b);
				if(maxLines<0){maxLines=Integer.MAX_VALUE;}
				
			}else if(a.equals("vlines") || a.equals("vsamples")){
				maxLinesV=Parse.parseIntKMG(b);
				if(maxLinesV<0){maxLinesV=Integer.MAX_VALUE;}
				
			}else if(a.equals("shuffleraw")){
				shuffleRaw=Parse.parseBoolean(b);
			
			}else if(a.equals("copynet")){
				copyNetInWorkerThread=Parse.parseBoolean(b);
			}else if(a.equals("setnet")){
				setNetInWorkerThread=Parse.parseBoolean(b);
				
			}else if(a.equals("dump") || a.equals("dumprate")){
				if(b!=null && b.length()>0 && Character.isLetter(b.charAt(0))){
					boolean x=Parse.parseBoolean(b);
					dumpRate=x ? 0.5f : 0f;
				}else {
					dumpRate=Float.parseFloat(b);
				}
			}else if(a.equals("dumpepoch") || a.equals("dumpbatch")){
				if(b.indexOf('.')<0){
					dumpEpoch=Parse.parseIntKMG(b);
				}else{
					dumpEpochMult=Float.parseFloat(b);
				}
			}else if(a.equals("dumpepochmult") || a.equals("dumpbatchmult")){
				dumpEpochMult=Float.parseFloat(b);
			}else if(a.equals("partialdump")){
				if(b==null || Character.isLetter(b.charAt(0))) {
					partialDumpFraction=(Parse.parseBoolean(b) ? partialDumpFraction : 1.0f);
				}else {
					partialDumpFraction=Float.parseFloat(b);
				}
			}else if(a.equals("halfdump")){
				if(b==null || Character.isLetter(b.charAt(0))) {
					partialDumpFraction=(Parse.parseBoolean(b) ? 0.5f : 1.0f);
				}else {
					partialDumpFraction=Float.parseFloat(b);
				}
				
			}else if(a.equals("exclusive")){
				exclusive=Parse.parseBoolean(b);
			}else if(a.equals("inclusive")){
				exclusive=!Parse.parseBoolean(b);
			}else if(a.equals("validatefraction") || a.equals("vfraction")) {
				validateFraction=Float.parseFloat(b);
				
//			}else if(a.equals("shrinksubsets")){
//				shrinkSubsets=Parse.parseBoolean(b);
				
			}else if(a.equals("samplesperepoch") || a.equals("spe") || a.equals("samplesperbatch") || a.equals("spb") || a.equals("batchsize")){
//				maxSamplesPerEpoch=Parse.parseIntKMG(b);
//				if(maxSamplesPerEpoch<0){maxSamplesPerEpoch=Integer.MAX_VALUE;}
				assert(false) : "lpe disabled";
			}else if(a.equals("fractionperepoch") || a.equals("fpe") || a.equals("fractionperbatch") || a.equals("fpb")){
				fractionPerEpoch0=Float.parseFloat(b);
			}else if(a.equals("fractionperepoch2") || a.equals("fpe2") || a.equals("fractionperbatch2") || a.equals("fpb2")){
				fractionPerEpoch2=Float.parseFloat(b);
			}else if(a.equals("fpestart") || a.equals("fractionalstart")){
				fpeStart=Parse.parseIntKMG(b);
			}else if(a.equals("fptriage") || a.equals("ptriage")){
				positiveTriage=Float.parseFloat(b);
			}else if(a.equals("fntriage") || a.equals("ntriage")){
				negativeTriage=Float.parseFloat(b);
			}else if(a.equals("triage")){
				positiveTriage=negativeTriage=Float.parseFloat(b);
			}else if(a.equals("starttriage") || a.equals("triagestart")){
				startTriage=Parse.parseIntKMG(b);
			}else if(a.equals("starttriagemult") || a.equals("triagestartmult")){
				startTriageMult=Float.parseFloat(b);
			}
				
			else if(a.equals("minweightepoch") || a.equals("minweightbatch") || a.equals("weightstart")){
				minWeightEpoch=Parse.parseIntKMG(b);
			}else if(a.equals("minweightepochmult") || a.equals("minweightbatchmult") || a.equals("weightstartmult")){
				minWeightEpochMult=Float.parseFloat(b);
			}
				
			else if(a.equalsIgnoreCase("weights") || a.equalsIgnoreCase("weighted")){
				DataLoader.weighted=Parse.parseBoolean(b);
			}
				
			else if(a.equals("lowweightannealcutoff") || a.equals("lwac")){
				Cell.setLowWeightAnnealCutoff(Float.parseFloat(b));
			}else if(a.equals("startanneal") || a.equals("annealstart") || a.equals("minannealepoch")){
				minAnnealEpoch=Parse.parseIntKMG(b);
			}else if(a.equals("stopanneal") || a.equals("annealstop") || a.equals("maxannealepoch")) {
				if(b.indexOf('.')>=0){
					maxAnnealEpochMult=Float.parseFloat(b);
					maxAnnealEpoch=Integer.MAX_VALUE;
				}else{
					maxAnnealEpoch=Parse.parseIntKMG(b);
				}
			}else if(a.equals("netseed")){
				netSeed0=Long.parseLong(b);
				setNetSeed=true;
			}else if(a.equals("annealseed")){
				assert(false) : "annealseed is no longer supported.";
//				annealSeed0=Long.parseLong(b);
//				setAnnealSeed=true;
			}else if(a.equals("shuffleseed")){
				SampleSet.shuffleSeed=Long.parseLong(b);
			}else if(a.equals("seed")){
				netSeed0=Long.parseLong(b);
				setNetSeed=true;
//				netSeed0=annealSeed0=Long.parseLong(b);
//				setNetSeed=setAnnealSeed=true;
			}else if(a.equals("density") || a.equals("maxdensity")){
				maxDensity=Float.parseFloat(b);
				assert(maxDensity>0 && maxDensity<=1) : "Density must be between 0 and 1: "+maxDensity;
//				maxEdges=Long.MAX_VALUE;
			}else if(a.equals("density1")){
				maxDensity1=Float.parseFloat(b);
			}else if(a.equals("densityarray")){
				densityArray=Parse.parseFloatArray(b, ",", "x");
			}else if(a.equals("edges") || a.equals("maxedges")){
				maxEdges=Parse.parseKMG(b);
//				maxDensity=0;
			}else if(a.equals("edgeblocksize") || a.equals("blocksize") || a.equals("quant") || a.equals("ebs")){
				edgeBlockSize=Parse.parseIntKMG(b);
			}
			
			else if(a.equals("dims") || a.equals("dimensions")){
				dims0=Parse.parseIntArray(b, ",", "x");
			}else if(a.equals("mindims") || a.equals("mindimensions")){
				minDims=Parse.parseIntArray(b, ",", "x");
			}else if(a.equals("maxdims") || a.equals("maxdimensions")){
				maxDims=Parse.parseIntArray(b, ",", "x");
			}else if(a.equals("in") || a.equals("data") || a.equals("din")){
				dataIn=b;
			}else if(a.equals("validateset") || a.equals("validate") || a.equals("vin") || a.equals("vset")){
				validateIn=b;
			}else if(a.equals("net") || a.equals("network") || a.equals("netin") || a.equals("networkin")){
				netIn=b;
			}else if(a.equals("out") || a.equals("netout") || a.equals("networkout") || a.equals("outfinal") || a.equals("outf")){
				netOutFinal=b;
			}else if(a.equals("outb") || a.equals("outbest") || a.equals("netoutbest")){
				netOutBest=b;
			}else if(a.equals("epochs") || a.equals("maxepochs") || a.equals("batches")){
				maxEpochs=Parse.parseIntKMG(b);
			}else if(a.equals("subsets") || a.equals("sets")){
				subsets=Parse.parseIntKMG(b);
				if(subsets>0) {setsize=-1;}
			}else if(a.equals("setsize")){
				setsize=Parse.parseIntKMG(b);
				if(setsize>0) {subsets=-1;}
			}else if(a.equals("crossover") || a.equals("xvr")){
				if(b==null || Tools.startsWithLetter(b)) {
					boolean x=Parse.parseBoolean(b);
					if(x) {
						if(crossoverFpMult<=0) {crossoverFpMult=1;}
					}else {
						if(crossoverFpMult>0) {crossoverFpMult=-1;}
					}
				}else {
					crossoverFpMult=Float.parseFloat(b);
				}
				if(crossoverFpMult>0) {targetFPR=targetFNR=-1;}
			}else if(a.equals("fpmult")){
				crossoverFpMult=Float.parseFloat(b);
				if(crossoverFpMult>0) {targetFPR=targetFNR=-1;}
			}else if(a.equals("error") || a.equals("targeterror")){
				targetError=Float.parseFloat(b);
			}else if(a.equals("targetfpr") || a.equals("fpr")){
				targetFPR=Float.parseFloat(b);
				if(targetFPR>0) {crossoverFpMult=targetFNR=-1;}
			}else if(a.equals("targetfnr") || a.equals("fnr")){
				targetFNR=Float.parseFloat(b);
				if(targetFNR>0) {targetFPR=crossoverFpMult=-1;}
			}else if(a.equals("printtpr")){boolean printAlpha=true, printAnneal=true;
				printTPR=Parse.parseBoolean(b);
			}else if(a.equals("printtnr")){
				printTNR=Parse.parseBoolean(b);
			}else if(a.equals("printalpha") || a.equals("printalp")){
				printAlpha=Parse.parseBoolean(b);
			}else if(a.equals("printanneal") || a.equals("printann")){
				printAnneal=Parse.parseBoolean(b);
			}else if(a.equals("alpha")){
				alphaZero=Float.parseFloat(b);
			}else if(a.equals("alphamult")){
				alphaMult=Float.parseFloat(b);
			}else if(a.equals("alphamult2")){
				alphaMult2=Float.parseFloat(b);
			}else if(a.equals("anneal") || a.equals("annealstrength")){
				if(b==null || b.length()<1 || Tools.startsWithLetter(b)) {
					annealStrength0=Parse.parseBoolean(b) ? 0.04f : 0;
				}else{
					annealStrength0=Float.parseFloat(b);
				}
			}else if(a.equals("annealprob") || a.equals("annealrate")) {
				annealProb=Float.parseFloat(b);
			}else if(a.equals("annealmult2")) {
				annealMult2=Float.parseFloat(b);
//			}else if(a.equals("annealdropoff")){
//				annealDropoff0=Float.parseFloat(b);
				//assert(false) : annealDropoff;
			}else if(a.equals("biasalphamult")){
				Cell.biasAlphaMult=Float.parseFloat(b);
			}else if(a.equals("biasannealmult")){
				Cell.biasAnnealMult=Float.parseFloat(b);
			}else if(a.equals("annealbias")){
				Cell.annealBias=Parse.parseBoolean(b);
			}else if(a.equals("sortall")){
				sortAll=Parse.parseBoolean(b);
			}else if(a.equals("sort")){
				sort=Parse.parseBoolean(b);
			}else if(a.equals("sortmt") || a.equals("mtsort")){
				allowMTSort=Parse.parseBoolean(b);
			}else if(a.equals("sortinthread")){
				sortInThread=Parse.parseBoolean(b);
			}else if(a.equals("pivotsort")){
				PIVOT_SORT=Parse.parseBoolean(b);
			}else if(a.equals("launchinthread")){
				launchInThread=Parse.parseBoolean(b);
			}else if(a.equals("setlock") || a.equals("usesetlock")){
				useSetLock=Parse.parseBoolean(b);
			}else if(a.equals("shuffle")){
				SampleSet.shuffle=Parse.parseBoolean(b);
			}else if(a.equals("shuffle2") || a.equals("shufflesubset")){
				shuffleSubset=Parse.parseBoolean(b);
			}else if(a.equals("shufflemod")){
				SHUFFLEMOD=Integer.parseInt(b);
			}else if(a.equals("ordered")){
//				assert(false) : "Ordered is currently forced.";
				orderedJobs=Parse.parseBoolean(b);
			}else if(a.equals("peakalphaepoch") || a.equals("peakalpha") || a.equals("pae") || a.equals("peakalphabatch") || a.equals("peakalpha") || a.equals("pab")){
				peakAlphaEpoch=Parse.parseIntKMG(b);
			}else if(a.equals("alphadropoff")){
//				alphaDropoff=Float.parseFloat(b);
			}else if(a.equals("alphamult")){
				alphaMult=Float.parseFloat(b);
			}
			
			else if(a.equals("booleancutoff") || a.equals("cutoff") || a.equals("thresh")){
				cutoffForEvaluation=Cell.cutoffForTraining=Float.parseFloat(b);
				setCutoffForEvaluation=Cell.setCutoffForTraining=true;
			}else if(a.equals("usemidpoint")){
				Cell.useMidpoint=Parse.parseBoolean(b);
			}else if(a.equals("booleancutoffbackprop") || a.equals("cutoffbackprop") || a.equals("cutofferror")
					|| a.equals("threshbackprop") || a.equals("trainingcutoff") || a.equals("cutofftraining")){
				Cell.cutoffForTraining=Float.parseFloat(b);
				Cell.setCutoffForTraining=true;
			}else if(a.equals("cutoffeval") || a.equals("evalcutoff")){
				cutoffForEvaluation=Float.parseFloat(b);
				setCutoffForEvaluation=true;
			}
			
//			else if(a.equals("booleancutoffgoal") || a.equals("cutoffgoal") || a.equals("threshgoal")){
//				booleanCutoffGoal=Float.parseFloat(b);
//			}
			
			else if(a.equals("positiveerrormult") || a.equals("pem")){
				Cell.positiveErrorMult=Float.parseFloat(b);
			}else if(a.equals("falsepositiveerrormult") || a.equals("fpem")){
				Cell.falsePositiveErrorMult=Float.parseFloat(b);
			}else if(a.equals("excesspositiveerrormult") || a.equals("epem")){
				Cell.excessPositiveErrorMult=Float.parseFloat(b);
			}else if(a.equals("negativeerrormult") || a.equals("nem")){
				Cell.negativeErrorMult=Float.parseFloat(b);
			}else if(a.equals("falsenegativeerrormult") || a.equals("fnem")){
				Cell.falseNegativeErrorMult=Float.parseFloat(b);
			}else if(a.equals("excessnegativeerrormult") || a.equals("enem")){
				Cell.excessNegativeErrorMult=Float.parseFloat(b);
			}else if(a.equals("falsepositeverrorincr") || a.equals("fpei")){
				Cell.fpErrorIncr=Float.parseFloat(b);
			}else if(a.equals("falsenegativeerrorincr") || a.equals("fnei")){
				Cell.fnErrorIncr=Float.parseFloat(b);
			}else if(a.equals("spread")){
				Cell.spread=Float.parseFloat(b);
			}else if(a.equals("excesspivotmult") || a.equals("pivotmult") || a.equals("epm")){
				Sample.excessPivotMult=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("edgeAmplitudeIncreaseMult") || a.equals("eaim")){
				Cell.edgeAmplitudeIncreaseMult=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("edgeAmplitudeIncreaseThresh") || a.equals("eait")){
				Cell.edgeAmplitudeIncreaseThresh=Float.parseFloat(b);
			}
			
			else if(a.equalsIgnoreCase("probFlat") || a.equals("pflat")){
				CellNet.PROB_FLAT=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("probExp") || a.equals("pexp")){
				CellNet.PROB_EXP=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("explamda") || a.equals("lamda") || a.equals("lambda")){
				CellNet.EXP_LAMDA=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("randweightcap") || a.equalsIgnoreCase("weightcap") || a.equals("expcap")){
				CellNet.RAND_WEIGHT_CAP=Float.parseFloat(b);
			}
			
			else if(a.equals("convertto01") || a.equals("converttoboolean")){
				Matrix.convertTo01=Parse.parseBoolean(b);
			}else if(a.equals("outputrangemin") || a.equals("rangemin") || a.equals("minoutput")){
				Matrix.targetOutputRangeMin=Float.parseFloat(b);
				Matrix.setTargetOutputRangeMin=true;
			}else if(a.equals("outputrangemax") || a.equals("rangemax") || a.equals("maxoutput")){
				Matrix.targetOutputRangeMax=Float.parseFloat(b);
				Matrix.setTargetOutputRangeMax=true;
			}
			
			
			else if(a.equals("profile")){
				Profiler.PROFILING=Parse.parseBoolean(b);
			}
			
			else if(a.equals("final") || a.equals("finallayer") || a.equals("finaltype")){
				Cell.finalLayerType=Function.toType(b);
			}else if(a.equals("type") || a.equals("defaulttype")){
				Cell.defaultActivationType=Function.toType(b);
			}else if(Function.toType(a, false)>=0){
				int type=Function.toType(a, true);
				if(b==null){
					Cell.defaultActivationType=type;
				}else{
					Function.TYPE_RATES[type]=Float.parseFloat(b);
				}
			}
			

			else if(a.equals("maxtype")){
				Cell.MAX_TYPE=Function.toType(b);
				outstream.println(Cell.MAX_TYPE);
			}
			
			//This is all deprecated
//			else if(a.equals("sigrate")){
//				Cell.TYPE_RATES[Cell.SIGMOID]=Float.parseFloat(b);
//			}else if(a.equals("tanhrate")){
//				Cell.TYPE_RATES[Cell.TANH]=Float.parseFloat(b);
//			}else if(a.equals("rslograte")){
//				Cell.TYPE_RATES[Cell.RSLOG]=Float.parseFloat(b);
//			}else if(a.equals("swishrate")){
//				Cell.TYPE_RATES[Cell.SWISH]=Float.parseFloat(b);
//			}
//			
//			else if(a.equals("sigfinal")){
//				Cell.finalLayerType=Cell.SIGMOID;
//			}else if(a.equals("tanhfinal")){
//				Cell.finalLayerType=Cell.TANH;
//			}else if(a.equals("rslogfinal")){
//				Cell.finalLayerType=Cell.RSLOG;
//			}else if(a.equals("swishfinal")){
//				Cell.finalLayerType=Cell.SWISH;
//			}
			
//			else if(a.equals("sigmoid") || a.equals("sig")){
//				Cell.defaultActivationType=Cell.SIGMOID;
//			}else if(a.equals("tanh")){
//				Cell.defaultActivationType=Cell.TANH;
//			}else if(a.equals("rslog")){
//				Cell.defaultActivationType=Cell.RSLOG;
//			}else if(a.equals("swish")){
//				Cell.defaultActivationType=Cell.SWISH;
//			}
			
			else if(a.equals("randomtype") || a.equals("randomtypes") || a.equals("randomtyperate")){
				if(b==null || Character.isLetter(b.charAt(0))) {
					Cell.randomTypeRate=(Parse.parseBoolean(b) ? 0.5f : 0f);
				}else{
					Cell.randomTypeRate=Float.parseFloat(b);
				}
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
			
		}
		
		if(!training){
			subsets=1;
			alphaZero=0;
			annealStrength0=0;
			maxEpochs=1;
			netOutFinal=netOutBest=null;
			quiet=true;
//			seedsToScan=0;
		}
		if(Parser.silent) {quiet=true;}
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		netIn=Tools.fixExtension(netIn);
		dataIn=Tools.fixExtension(dataIn);
		validateIn=Tools.fixExtension(validateIn);
		if(netIn==null && dims0==null && (minDims==null || maxDims==null)){//TODO: This fired when dims was specified
			throw new RuntimeException("Error - an input net file or dims is required.");
		}
		if(dataIn==null && (validateIn==null || training)){throw new RuntimeException("Error - a data file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, true, netOutFinal, netOutBest, dataOut)){
			outstream.println((netOutFinal==null)+", "+netOutFinal);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+
					netOutFinal+", "+netOutFinal+", "+dataOut+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, false, dataIn, netIn)){
			//throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(netOutFinal==null && netOutBest==null && maxEpochs>1) {
			if(training) {outstream.println("Warning - no output file specified.");}
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, dataIn, netIn, netOutFinal, dataOut)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, dataIn, netIn, netOutBest, dataOut)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
//		assert(false) : "TODO";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		 //TODO: ignore the top X% most ideosyncratic cases when training (highest error magnitude) after every full pass
		//Simply swap them to the very end, and update their epoch to the future so they will stay there for a while
		
		//Reset counters
		linesProcessed=linesOut=0;
		bytesProcessed=bytesOut=0;
		
		int dataLines=training ? maxLines : Tools.min(maxLines, maxLinesV);
		
		if(validateIn!=null) {validateFraction=0;}
		SampleSet[] ssa=(dataIn==null ? null : 
			loadData(dataIn, training, dataLines, shuffleRaw, validateFraction, maxLinesV));
		
		data=(ssa==null ? null : ssa[0]);
		validateSet=(validateIn==null ? ssa[1] : 
			loadData(validateIn, false, maxLinesV, shuffleRaw && maxLinesV<Shared.MAX_ARRAY_LEN, 0, maxLinesV)[0]);
		
//		data=(dataIn==null ? null : 
//			loadData(dataIn, training, dataLines, shuffleRaw));
//		validateSet=(validateIn==null ? data.copy(maxLinesV, 0) : 
//			loadData(validateIn, false, maxLinesV, shuffleRaw && maxLinesV<Shared.MAX_ARRAY_LEN));

		if(!setCutoffForEvaluation && data!=null && Cell.useMidpoint){cutoffForEvaluation=data.outputMidpoint();}
		if(!Cell.setCutoffForTraining && data!=null && Cell.useMidpoint){Cell.cutoffForTraining=data.outputMidpoint();}
		
//		assert(false) : data+", "+validateSet;
		
		//Process the reads in separate threads
		spawnThreads();
		
		writeNetwork(finalNet, ffNetOutFinal);
		writeData(ffDataOut);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		
		//Report timing and results
		t.stop();
		
		if(!quiet) {
			long bytesProcessed=linesProcessed*4*(numInputs+numOutputs);
//		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		outstream.println(Tools.timeLinesProcessed(t, linesProcessed, 8));
//		outstream.println(Tools.thingsProcessed(t.elapsed, linesProcessed*(net0.list.size()-1), 8, "Cells"));
		
		//TODO: This is not correct for variable-width networks
		outstream.println(Tools.thingsProcessed(t.elapsed, linesProcessed*(edgesIn/networksIn), 8, "Edges"));
		
		outstream.println();
		outstream.println("Lines In:          \t"+(DataLoader.lastValidLines+DataLoader.lastInvalidLines));
//		outstream.println("Valid Lines In:    \t"+DataLoader.lastValidLines);
		if(DataLoader.lastInvalidLines>0) {
			outstream.println("Invalid Lines In:  \t"+DataLoader.lastInvalidLines);
		}
		
		if(linesOut>0) {
			outstream.println("Lines Out:         \t"+linesOut);
			outstream.println("Bytes Out:         \t"+bytesOut);
		}
		}
		
		if(machineOut && bestNet!=null && bestNet.dims!=null){
			String sh=toMachineHeader(bestNet.dims);
			String ss=toMachineOut(t);
			outstream.println(sh);
			outstream.println(ss);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	private String toMachineHeader(int[] dims) {
		final boolean printAvg=(true || networksTested>1);
		final boolean printAvgDims=(true || (maxDims!=null && printAvg));
		ByteBuilder bb=new ByteBuilder();
		bb.append("machineheader:");
		
		//99.9697571      0.010117        0.020011        0.000140        0.078500        0.381700        0.539660        0.847933
		//0.009925        0.019858        0.000140        0.100300        0.359900        0.539660        0.863700        84

		bb.append("\ttime");
		
		bb.append("\tbERR").append("\tbWER").append("\tbFPR").append("\tbFNR");
		bb.append("\tbTPR").append("\tbTNR").append("\tbCTF");
		
		if(printAvg) {
			bb.append("\taERR").append("\taWER").append("\taFPR").append("\taFNR");
			bb.append("\taTPR").append("\taTNR").append("\taCTF");
		}
		
		for(int i=1; i<dims.length-1; i++) {
			bb.tab().append("bD").append(i);
		}
		
		if(printAvgDims) {
			for(int i=1; i<dims.length-1; i++) {
				bb.tab().append("aD").append(i);
			}
		}
		
		bb.tab().append("CL");
		return bb.toString();
	}
	
	private String toMachineOut(Timer t) {
		final boolean printAvg=(true || networksTested>1);
		final boolean printAvgDims=(true || (maxDims!=null && printAvg));
		StringBuilder sb=new StringBuilder();
		sb.append("machineout:");
		sb.append('\t').append(t.timeInSeconds());

		String bdims=this.genDimsString(bestNet.dims, true);
		String adims="";
		if(printAvgDims) {
			adims=(finalNets!=null && finalNets.size()>1 ? this.genAvgDimsString(finalNets, true) : bdims);
		}
		String bstats=genPrintStringMachine(bestNet);
		
		CellNet avg=finalNet;
		if(finalNets!=null && finalNets.size()>1){avg=makeAvg(finalNets);}
		String astats=genPrintStringMachine(avg);

		sb.append('\t').append(bstats);
		if(printAvg) {sb.append('\t').append(astats);}
		sb.append('\t').append(bdims);
		if(printAvgDims) {sb.append('\t').append(adims);}
		sb.append('\t').append(commandLine);
		
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(){

		//Do anything necessary prior to processing
		annealStrength=annealStrength0;
		int annealEpochs=Tools.min(maxEpochs, maxAnnealEpoch)-minAnnealEpoch;
		annealDropoff=1.0/Math.pow(annealMult2, 1.0/annealEpochs);//annealDropoff0;
		annealDropoff0=annealDropoff;

		//Determine how many threads may be used
		final int samples=(data!=null ? data.samples.length : validateSet.samples.length);
		threads=Tools.mid(1, Shared.threads(), Tools.max(samples, validateSet.samples.length));
		jobsPerEpoch=Tools.mid(1, threads, (int)(4*(threads+(networksPerCycle-1))/(1.6*networksPerCycle)));
		int seeds=calcSeedsToScan();
		int nets=networksPerCycle*cycles;
		String ss=(seeds>nets ? seeds+" seeds, " : "");
		if(!quiet) {
			outstream.println("Using "+plural("worker thread", threads)+", "
					+plural("trainer thread", networksPerCycle)+", "+ss+plural("job", jobsPerEpoch)+"/trainer/batch"
					+", "+pluralES("batch", maxEpochs)+", "+plural("cycle", cycles)
					+" to train "+plural("network", networksPerCycle*cycles)+" total.");
		}
		
		networkQueue=new ArrayBlockingQueue<CellNet>(networksPerCycle);
		workerQueue=new ArrayBlockingQueue<JobData>(jobsPerEpoch*networksPerCycle);

		//Fill a list with WorkerThreads
		alwt=new ArrayList<WorkerThread>(threads);
		for(int i=0; i<threads; i++){
			alwt.add(new WorkerThread(i, workerQueue, cutoffForEvaluation));
		}
		
		//Start the threads and wait for them to finish
		ThreadWaiter.startThreads(alwt);

		
		finalNet=trainNetworks();
		
		//			outstream.println("Trainer finished, shutting down workers.");
		for(int i=0; i<threads; i++) {
			workerQueue.add(JobData.POISON);
		}
		
		boolean success=ThreadWaiter.waitForThreadsToFinish(alwt, this);
		mprof.printTimes();
		wprof.printTimes();
		errorState&=!success;
//		outstream.println("Master finished.");
		
		//Do anything necessary after processing

	}
	
	private int calcSeedsToScan() {
		if(seedsToScan>=0) {return seedsToScan;}
		if(netIn!=null) {return 0;}
		final int totalNets=networksPerCycle*cycles;
		long x=Tools.max(seedsToScan, (int)(Tools.min(Integer.MAX_VALUE, totalNets*seedsToScanMult)));
		return (int)Tools.min(x, Shared.MAX_ARRAY_LEN);
	}
	
	private CellNet trainNetworks() {
		bestNet=null;
		finalNet=null;
		
		final int totalNets=networksPerCycle*cycles;
		final int seedsToScan=calcSeedsToScan();
		if(seedsToScan>totalNets) {
			outstream.println("Scanning "+seedsToScan+" seeds to find the best "+totalNets+".");
			Timer t=new Timer();
			ArrayList<Seed> seedList0=scanForSeeds(seedsToScan, totalNets);
			seedList=new LongList(seedList0.size());
			for(Seed s : seedList0) {seedList.add(s.netSeed);}
			t.stopAndPrint();
		}
//		assert(false) : seedsToScan+", "+scanMult+", "+totalNets;
		
		for(int i=0; i<cycles; i++) {
			if(cycles>1) {outstream.println("Cycle "+i);}
			clearCycle();
			CellNet s=runCycle();
			if(finalNet==null || s.compareTo(finalNet)>0) {
				finalNet=s;
			}
		}
		if(!quiet) {outstream.println();}
		if(networksTested>1) {
			outstream.println("Best:\t\t"+genPrintString(bestNet, -1, false, 
					printFPR, printFNR, printTPR, printTNR, printCutoff, false, false));
			outstream.println("Final:\t\t"+genPrintString(finalNet, -1, false, 
					printFPR, printFNR, printTPR, printTNR, printCutoff, false, false));
			if(printMode==printAverage) {
			outstream.println("Average:  \t"+genPrintString(makeAvg(finalNets), -1, false, 
					printFPR, printFNR, printTPR, printTNR, printCutoff, false, false));
			}
			
			if(dims0==null) {
				outstream.println("\nBest Dims:\t"+Arrays.toString(bestNet.dims));
				outstream.println("Final Dims:\t"+Arrays.toString(finalNet.dims));
				if(printMode==printAverage) {
					outstream.println("Average Dims:  \t"+genAvgDimsString(finalNets, false));
				}
			}
			
			outstream.println("\nBest Seed:\t"+bestNet.seed);
		}else{
			if(bestNet!=null) {
				outstream.println("Best:\t\t"+genPrintString(bestNet, -1, false, 
						printFPR, printFNR, printTPR, printTNR, printCutoff, false, false));
			}
			outstream.println("Final:\t\t"+genPrintString(finalNet, -1, false, 
					printFPR, printFNR, printTPR, printTNR, printCutoff, false, false));
			if(!quiet) {outstream.println();}
		}
		return finalNet;
	}
	
	private CellNet runCycle() {
		final CellNet[] nets=fetchNetworks(netIn, networksPerCycle);
		assert(networksPerCycle==nets.length);
		ArrayList<TrainerThread> altt=new ArrayList<TrainerThread>(networksPerCycle);
		for(int i=0; i<nets.length; i++) {
			networksIn++;
			edgesIn+=nets[i].countEdges();
			altt.add(new TrainerThread(this, nets[i]));
		}
		for(TrainerThread tt : altt) {
			tt.start();
		}
		CellNet best=null;
		for(int i=0; i<networksPerCycle; i++) {
			CellNet net=fetchFromQueue(networkQueue);
			if(printMode==printAverage) {
				finalNets.add(net);
			}
			best=compare(net, best);
		}
		for(TrainerThread tt : altt) {
			accumulate(tt);
		}
		return best;
	}
	
	<X> X fetchFromQueue(ArrayBlockingQueue<X> queue) {
		X x=null;
		while(x==null) {
			try {
				x=queue.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		return x;
	}

	synchronized CellNet compareWithBest(CellNet net) {
		CellNet better=compare(net, bestNet);
//		assert(status!=null); //for testing
		if(printStatus && !quiet) {
			printStatus(net);
		}
		if(better!=bestNet) {
			bestNet=better;//TODO: Can be a setFrom instead
			writeNetwork(bestNet, ffNetOutBest);
		}
		return bestNet;
	}
	
	synchronized CellNet compare(CellNet a, CellNet b) {
		assert(a!=null);
		int x=a.compareTo(b);
		return x>0 ? a : b;
	}
	
	synchronized Status compare(Status a, Status b) {
		assert(a!=null);
		if(b==null) {return a;}
		int x=a.compareTo(b);
		return x>0 ? a : b;
	}
	
	private final void accumulate(TrainerThread tt){
		synchronized(tt) {
			mprof.accumulate(tt.mprof);
			networksTested++;
//			edgesTested+=tt.net0.countEdges();
		}
	}
	
	private final void accumulate(ScannerThread tt){
		synchronized(tt) {
			mprof.accumulate(tt.mprof);
		}
	}
	
	@Override
	public final void accumulate(WorkerThread wt){
		synchronized(wt) {
			wprof.accumulate(wt.tprof);
			linesProcessed+=wt.linesProcessedT;
			linesOut+=wt.linesOutT;
			errorState|=(!wt.success);
			errorState|=(wt.errorStateT);
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	@Deprecated
//	private LongList scanForSeeds_old(final int seedsToScan0, final int seedsToReturn0) {
//		int seedsToScan=seedsToScan0, seedsToReturn=seedsToReturn0;
//		final int stsPerThread=(seedsToScan+networksPerCycle-1)/networksPerCycle;
//		final int strPerThread=(seedsToReturn+networksPerCycle-1)/networksPerCycle;
//		
//		seedQueue=new ArrayBlockingQueue<ArrayList<Seed>>(networksPerCycle);
//		ArrayList<ScannerThread> alst=new ArrayList<ScannerThread>(networksPerCycle);
//		for(int i=0; i<networksPerCycle; i++) {
//			final long seed=(i>0 || netSeed0<0 ? randyNetSeed.nextLong()&Long.MAX_VALUE : netSeed0);
//			final int stsThisThread=Tools.min(seedsToScan, stsPerThread);
//			final int strThisThread=Tools.min(seedsToReturn, strPerThread);
//			seedsToScan-=stsThisThread;
//			seedsToReturn-=strThisThread;
//			final ScannerThread st=new ScannerThread(this, dims0, minDims, maxDims, stsThisThread, 
//					strThisThread, scanEpochs, scanSamples, seed, seedQueue);
//			alst.add(st);
//		}
//		for(ScannerThread st : alst) {
//			st.start();
//		}
//		
//		LongList seeds=new LongList(seedsToReturn0);
//		for(int i=0; i<networksPerCycle; i++) {
//			LongList seedsFromQueue=fetchFromQueue(seedQueue);
//			seeds.append(seedsFromQueue);
//		}
//		assert(seeds.size()==seedsToReturn0) : seeds.size();
//		seeds.sort();
//		seeds.shrinkToUnique();
//		
//		for(ScannerThread st : alst) {
//			accumulate(st);
//		}
//		return seeds;
//	}
	
	private ArrayList<Seed> scanForSeeds(final int seedsToScan0, final int seedsToReturn0) {
		int seedsToScan=seedsToScan0, seedsToReturn=seedsToReturn0;
		final int stsPerThread=(seedsToScan+networksPerCycle-1)/networksPerCycle;
		final int strPerThread=(seedsToReturn+networksPerCycle-1)/networksPerCycle;
		
		seedQueue=new ArrayBlockingQueue<ArrayList<Seed>>(networksPerCycle);
		ArrayList<ScannerThread> alst=new ArrayList<ScannerThread>(networksPerCycle);
		for(int i=0; i<networksPerCycle; i++) {
			final long seed=(i>0 || netSeed0<0 ? randyNetSeed.nextLong()&Long.MAX_VALUE : netSeed0);
			final int stsThisThread=Tools.min(seedsToScan, stsPerThread);
			final int strThisThread=Tools.min(seedsToReturn, strPerThread);
			seedsToScan-=stsThisThread;
			seedsToReturn-=strThisThread;
			final ScannerThread st=new ScannerThread(this, dims0, minDims, maxDims, stsThisThread, 
					strThisThread, scanEpochs, scanSamples, seed, seedQueue);
			alst.add(st);
		}
		for(ScannerThread st : alst) {
			st.start();
		}
		
		ArrayList<Seed> seeds=new ArrayList<Seed>(seedsToReturn0);
		for(int i=0; i<networksPerCycle; i++) {
			ArrayList<Seed> seedsFromQueue=fetchFromQueue(seedQueue);
			seeds.addAll(seedsFromQueue);
		}
		assert(seeds.size()==seedsToReturn0) : seeds.size();
		Collections.sort(seeds);
		
		for(ScannerThread st : alst) {
			accumulate(st);
		}
		return seeds;
	}
	
	private CellNet[] fetchNetworks(String path, int count) {
		return path==null ? createNetworks(count) : loadNetworks(path);
	}
	
	private CellNet[] createNetworks(int count){
		CellNet[] nets=new CellNet[count];
		for(int i=0; i<count; i++) {
			nets[i]=randomNetwork();
		}
		return nets;
	}
	
	private CellNet randomNetwork() {
		final long seed=(seedList!=null && netsMade<seedList.size() ? seedList.get(netsMade) : 
			(netsMade>0 || netSeed0<0 ? randyNetSeed.nextLong()&Long.MAX_VALUE : netSeed0));
//		final long aseed=(netsMade>0 || annealSeed0<0 ? 
//				randyAnnealSeed.nextLong()&Long.MAX_VALUE : annealSeed0);
		netsMade++;
		return randomNetwork(seed);
	}
	
	public final CellNet randomNetwork(long seed) {
		final int[] dims=(dims0!=null ? dims0 : makeDims(minDims, maxDims, seed));
		return randomNetwork(dims, seed, maxDensity, maxDensity1, maxEdges, edgeBlockSize, commands);
	}
	
	public static final CellNet randomNetwork(int[] dims, long seed, float maxDensity,
			float maxDensity1, long maxEdges, int edgeBlockSize, ArrayList<String> commands) {
		final long fullyConnectedEdges=fullyConnectedEdges(dims);
		float density=maxEdges<1 ? maxDensity : 
			Tools.min(maxDensity, (float)(maxEdges/(double)fullyConnectedEdges));
		float density1=(maxDensity1<=0 ? density : maxDensity1);
		assert(density>0f && density<=1f) : density+", "+maxEdges+", "+fullyConnectedEdges;
		final CellNet net=new CellNet(dims, seed, density, density1, edgeBlockSize, commands);
//		net.annealSeed=aseed;
		net.randomize();
		return net;
	}
	
	public static long fullyConnectedEdges(int[] dims) {
		long sum=0;
		for(int i=1; i<dims.length; i++) {
			long edges=dims[i-1]*dims[i];
			sum+=edges;
		}
		return sum;
	}
	
	public static final int[] makeDims(int[] minDims, int[] maxDims, long seed) {
		Random randy=Shared.threadLocalRandom(seed);
		assert(minDims.length==maxDims.length);
		assert(minDims[0]==maxDims[0]);
		assert(minDims[minDims.length-1]==maxDims[maxDims.length-1]);
		int[] dims=minDims.clone();
		for(int i=1; i<dims.length-1; i++) {
			assert(minDims[i]<=maxDims[i]);
			dims[i]=minDims[i]+randy.nextInt(maxDims[i]-minDims[i]+1);
		}
		return dims;
	}
	
	private CellNet[] loadNetworks(final String path){
		assert(path!=null);
		String[] paths=path.split(",");
		networksPerCycle=(paths==null ? networksPerCycle : paths.length);
		CellNet[] nets=new CellNet[networksPerCycle];
		for(int i=0; i<paths.length; i++) {
			nets[i]=loadNetwork(paths[i]);
		}
		return nets;
	}
	
	private CellNet loadNetwork(String path){
		final CellNet net=CellNetParser.load(path, false);
		if(setCutoffForEvaluation) {net.setCutoff(cutoffForEvaluation);}
//		dims=net.dims.clone();
		assert(numInputs<0 || net.numInputs()==numInputs);
		assert(numOutputs<0 || net.numOutputs()==numOutputs);
		numLayers=net.numLayers();
		numInputs=net.numInputs();
		numOutputs=net.numOutputs();
		net.commands.add(commandLine);
//		if(net.annealSeed<0) {net.annealSeed=annealSeed0;}
		return net;
	}
	
	private synchronized void writeNetwork(CellNet net, FileFormat ff) {
		if(ff==null) {return;}
		if(net.lastStats==null) {
			net.lastStats=genPrintStringFull(net, -1);
		}
		ByteStreamWriter bsw=makeBSW(ff);
		ByteBuilder bb=net.toBytes();
		linesOut+=net.lastLinesWritten;
		bytesOut+=bb.length();
		bsw.println(bb);
		bsw.poisonAndWait();
	}
	
	private SampleSet[] loadData(String path, boolean makeSubsets, int maxLines0, 
			boolean shuffleRaw, float splitFraction, int maxLines1){
		if(!quiet) {outstream.println("Loading "+path);}
		SampleSet[] ssa=DataLoader.load(path, maxLines0, shuffleRaw, splitFraction, maxLines1, exclusive, balance);
		for(SampleSet ss : ssa) {
			ss.makeSamples();
		}
		if(makeSubsets) {
			if(setsize>0) {
				assert(setsize>=100) : "Setsize should be at least 100";
				subsets=Tools.max(1, ssa[0].samples.length/setsize);
				outstream.println("Data was organized into "+subsets+(subsets==1 ? " set." : " sets."));
			}
			subsets=Tools.mid(1, subsets, ssa[0].samples.length);
			ssa[0].makeSubsets(subsets);
		}
		return ssa;
	}
	
	private void writeData(FileFormat ff){
		if(ff==null) {return;}
		ByteStreamWriter bsw=makeBSW(ff);
		ByteBuilder bb=new ByteBuilder();
		bb.append("dims").tab().append(data.matrix.numInputs).tab().append(data.matrix.numOutputs).nl();
		if(data.matrix.columns!=null) {
			for(String s : data.matrix.columns) {
				bb.append(s).tab();
			}
			if(bb.endsWith('\t')) {
				bb.trimLast(1).nl();
			}
		}
		for(Sample s : data.samples) {
			s.toBytes(bb);
			bsw.print(bb);
			linesOut++;
			bytesOut+=bb.length();
			bb.clear();
		}
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	
	void clearCycle() {
		lastPrintEpoch=-1;
		lastPrintTime=System.currentTimeMillis();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Printing           ----------------*/
	/*--------------------------------------------------------------*/

	private static String plural(String s, int count){return count+" "+(count==1 ? s : s+"s");}
	private static String pluralES(String s, int count){return count+" "+(count==1 ? s : s+"es");}
	
	synchronized void printStatus(CellNet net) {
		if(networksPerCycle<=1 || printMode==printAll || (printMode==printFirst && net.epoch>lastPrintEpoch)) {
			printStatusInner(net);
		}else if(printMode==printBest){
			final CellNet old=networkMap.get(net);
			CellNet best=old;
			if(old==null) {
				networkMap.put(net, net);
				best=net;
			}else{
				int x=net.compareTo(old);
				if(x>0) {
					net.count+=old.count;
					networkMap.put(net, net);
					best=net;
				}else{
					old.count+=net.count;
				}
			}
			if(best.count>=networksPerCycle){
				printStatusInner(best);
				networkMap.remove(best);
			}
		}else if(printMode==printAverage){
			ArrayList<CellNet> list=networkMap2.get(net);
			if(list==null) {
				list=new ArrayList<CellNet>(networksPerCycle);
				networkMap2.put(net, list);
			}
			list.add(net);
			if(list.size()>=networksPerCycle){
				CellNet avg=makeAvg(list);
				printStatusInner(avg);
//				Object o=networkMap.remove(list.get(0));//For some reason this doesn't work
//				assert(o==list) : o+", "+list.get(0).epoch;
				list.clear();
			}
		}
	}
	
	private ArrayList<CellNet> concentrateBest(ArrayList<CellNet> list){
		Collections.sort(list);
		Collections.reverse(list);
		final int count=(int)Math.ceil(Math.pow(list.size(), 0.6)); //Best fraction of nets
		final int start=(count>4 ? 1 : 0); //If there are enough samples, ignore the best as an outlier
		final ArrayList<CellNet> concentrate=new ArrayList<CellNet>(count);
		int x=0; //Ensure the correct number were sampled
		for(int i=start; i<count+start; i++) {
			x++;
			CellNet net=list.get(i);
			concentrate.add(net);
		}
		assert(x==count);
		return concentrate;
	}
	
	private CellNet makeAvg(ArrayList<CellNet> list0) {
		ArrayList<CellNet> list=concentrateBest(list0);
		final CellNet zero=list.get(0);
		CellNet base=new CellNet(new int[] {1,1}, 1, 1.0f, 0f, 1, null);
		base.epoch=zero.epoch;
		base.alpha=zero.alpha;
		base.annealStrength=zero.annealStrength;
		base.errorRate=0;
		base.weightedErrorRate=0;
		base.fpRate=0;
		base.fnRate=0;
		base.tpRate=0;
		base.tnRate=0;
		float cutoff=0;
		for(CellNet net : list) {
			base.errorRate+=net.errorRate;
			base.weightedErrorRate+=net.weightedErrorRate;
			base.fpRate+=net.fpRate;
			base.fnRate+=net.fnRate;
			base.tpRate+=net.tpRate;
			base.tnRate+=net.tnRate;
			cutoff+=net.cutoff;
		}
		final float mult=1.0f/list.size();
		base.errorRate*=mult;
		base.weightedErrorRate*=mult;
		base.fpRate*=mult;
		base.fnRate*=mult;
		base.tpRate*=mult;
		base.tnRate*=mult;
		base.setCutoff(cutoff*mult);
		return base;
	}
	
	private float[] genAvgDims(ArrayList<CellNet> list0) {
		ArrayList<CellNet> list=concentrateBest(list0);
		final CellNet zero=list.get(0);
		float[] dims=new float[zero.dims.length];
		for(CellNet net : list) {
			for(int i=0; i<dims.length; i++) {
				dims[i]+=net.dims[i];
			}
		}
		for(int i=0; i<dims.length; i++) {
			dims[i]/=list.size();
		}
		return dims;
	}
	
	private String genAvgDimsString(ArrayList<CellNet> list0, boolean machine) {
		float[] dims=genAvgDims(list0);
		int[] dims2=new int[dims.length];
		for(int i=0; i<dims.length; i++) {
			dims2[i]=(int)Math.round(dims[i]);
		}
		return genDimsString(dims2, machine);
	}
	
	private String genDimsString(int[] dims, boolean machine) {
		ByteBuilder bb=new ByteBuilder();
		if(machine) {
			for(int i=1; i<dims.length-1; i++) {
//				bb.append(f, 1).comma().space();
				bb.append(dims[i]).tab();
			}
			bb.setLength(bb.length-1);
		}else {
			bb.append('[');
			for(int f : dims) {
				bb.append(f).comma().space();
			}
			bb.setLength(bb.length-2);
			bb.append(']');
		}
		return bb.toString();
	}
	
	private void printStatusInner(CellNet net) {
		final long t=System.currentTimeMillis();
		final long elapsed=t-lastPrintTime;
		lastPrintTime=t;
		
		net.lastStats=genPrintStringFull(net, elapsed);
		
		String s=genPrintStringDefault(net, elapsed);
		outstream.println(s);
		lastPrintEpoch=Tools.max(net.epoch, lastPrintEpoch);
	}
	
	synchronized void printStatusOld(CellNet net, String status, int epoch) {
		assert(status!=null);
		if(epoch>lastPrintEpoch || printMode==printAll || printMode==printBest){
			outstream.println("seed="+net.seed+"\t"+status);
			lastPrintEpoch=Tools.max(epoch, lastPrintEpoch);
		}
	}

	String genPrintStringFull(CellNet net, long millis) {
		String s=genPrintString(net, millis, true, true, true, true, true, true, true, true);
		return s;
	}
	String genPrintStringDefault(CellNet stat, long millis) {
		String s=genPrintString(stat, millis, true, printFPR, printFNR, printTPR, printTNR, printCutoff, printAlpha, printAnneal);
		return s;
	}
	
	private String genPrintString(CellNet net, long millis, boolean epo,
			boolean fpr, boolean fnr, boolean tpr, boolean tnr, boolean ctf, boolean alp, boolean ann) {
		String s="ERR= %.6f\tWER= %.6f";
		s=String.format(s, net.errorRate, net.weightedErrorRate);
		if(epo) {
//			s=String.format("Epoch %d:\t", net.epoch)+s;
			int e=net.epoch;
			s="Batch "+(e<10000 ? e : ((e/1000)+"k"))+":\t"+s;
		}
		if(fpr) {
			s+=String.format("\tFPR= %.6f", net.fpRate);
		}
		if(fnr) {
			s+=String.format("\tFNR= %.6f", net.fnRate);
		}
		if(tpr) {
			s+=String.format("\tTPR= %.6f", net.tpRate);
		}
		if(tnr) {
			s+=String.format("\tTNR= %.6f", net.tnRate);
		}
		if(ctf) {
			s+=String.format("\tCTF= %.6f", net.cutoff);
		}
		if(alp) {
			s+=String.format("\tALP= %.5f", net.alpha);
		}
		if(ann) {
			s+=String.format("\tANN= %.6f", (net.epoch>=minAnnealEpoch ? net.annealStrength : 0));
		}
		if(millis>=0) {
			if(millis>=100000) {
				s+="\t"+(millis/1000)+"s";
			}else if(millis>10000){
				s+=String.format("\t%.1fs", (millis/1000f));
			}else{
				s+=String.format("\t%.2fs", (millis/1000f));
			}
		}
		return s;
	}
	
	private String genPrintStringMachine(CellNet net) {
		ByteBuilder bb=new ByteBuilder();
		
		bb.append(net.errorRate, 6);
		bb.tab().append(net.weightedErrorRate, 6);
		bb.tab().append(net.fpRate, 6);
		bb.tab().append(net.fnRate, 6);
		bb.tab().append(net.tpRate, 6);
		bb.tab().append(net.tnRate, 6);
		bb.tab().append(net.cutoff, 6);
		return bb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Streams           ----------------*/
	/*--------------------------------------------------------------*/
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final String commandLine;
	final ArrayList<String> commands;
	
	/** Primary input file path */
	private String netIn=null;

	/** Secondary input file path */
	private String dataIn=null;

	/** Secondary input file path */
	private String validateIn=null;

	/** Output data path */
	private String dataOut=null;

	/** Primary output file path */
	private String netOutFinal=null;

	/** Output file path for best intermediate network */
	private String netOutBest=null;

	private int[] dims0;
	private int[] minDims;
	private int[] maxDims;
	private int numInputs=-1, numOutputs=-1, numLayers=-1;

	SampleSet data;
	SampleSet validateSet;
	
	float validateFraction=0.1f;
	boolean exclusive=true;

	private int netsMade=0;
	private float maxDensity=1.0f;
	//For layer 1 input only (between layers 0 and 1).
	private float maxDensity1=0f;
	private int edgeBlockSize=1;
	private long maxEdges=-1;
	private float[] densityArray=null;
	private long netSeed0=-1;
//	private long annealSeed0=7;
	private boolean setNetSeed=false;
//	private boolean setAnnealSeed=false;
	
	int threads=1;
	
//	private AtomicInteger nextLine=new AtomicInteger(0);
	
//	private boolean working=true;
	
	boolean orderedJobs=true; //Without ordered, very very slight nondeterminism.
	//TODO: Ordered could be sped up a little by adding a hashmap of results as they become available
	ArrayBlockingQueue<CellNet> networkQueue;
	ArrayBlockingQueue<JobData> workerQueue;
//	ArrayBlockingQueue<LongList> seedQueue;
	ArrayBlockingQueue<ArrayList<Seed>> seedQueue;
	private final Profiler mprof=new Profiler("M", 13);
	private final Profiler wprof=new Profiler("W", 7);
	
	ArrayList<WorkerThread> alwt;
	
	boolean training=true;
	
	/*--------------------------------------------------------------*/
	
	int maxEpochs=400000;
	float targetError=-1;
	float targetFPR=-1;
	float targetFNR=-1;
	float crossoverFpMult=1;
	
	boolean sortAll=false;
	boolean sort=true;
	boolean sortInThread=false;
	boolean shuffleSubset=true;
	boolean launchInThread=false;
	boolean allowMTSort=true;
	
	static boolean PIVOT_SORT=false;
	static int SHUFFLEMOD=4;//subcycle to shuffle on. 4 seems better than 6.
	
	//setlock=t setnet=t copynet=f was 2x slower,
	//also memory use never grew over time.  These may be related since
	//a lack of memory pressure may have led to objects staying in nonlocal memory
	//...However, on 2nd and 3rd attempts, this slowness could not be replicated.
	//setlock=t setnet=f copynet=t was still fastest, but ttf uses the least memory and is only 4% slower
	boolean useSetLock=true;//Nuanced but mostly good
	static boolean setNetInWorkerThread=true;//Higher concurrency per net but more CPU usage
	static boolean copyNetInWorkerThread=false;//Lower concurrency per net but less CPU usage

//#simd=f 2385.891  16x    648m4.438s
//#tft    1439.468  12x    291m4.056s
//#ttf    1082.808  19.4x  350m20.804s //Uses least memory
//#fft    1365.212  12.6x  287m11.788s
//#ftf    1245.100  17x    354m10.348s
	
	double alphaZero=0.08f;
	double alphaMult=5.0f;
	double alphaMult2=2.8f;
	int peakAlphaEpoch=640;
	
	float annealStrength0=.003f;
	float annealProb=.225f;
	float annealMult2=800;//2000;
	double annealDropoff0=0.999975f;
	
	float maxAnnealEpochMult=0.8f;
	int minAnnealEpoch=400;
	int maxAnnealEpoch=Integer.MAX_VALUE;
	
//	static float booleanCutoffGoal=0.5f;
	static float cutoffForEvaluation=0.5f;
	static boolean setCutoffForEvaluation=false;
	
	int subsets=-1;
	int setsize=60000;

	float fractionPerEpoch0=0.08f;
	float fractionPerEpoch2=0.08f;
	int fpeStart=192;
	
	//TODO: Dump is nondeterministic, at least with halfdump enabled.  Probably need a stable sort or better compareTo.
	//Dump is mainly to get rid of low-error samples to prevent them from getting sorted repeatedly.
	float dumpRate=0.5f;//0.5 or 0.6 is good with partialDump; otherwise, 0.3 at most partialDump false.  Could also just dump everything below an error cutoff.
	private float dumpEpochMult=0.25f;
	int dumpEpoch=-1;//TODO: Could be an error threshold, or once the error starts increasing, or a fraction
	float partialDumpFraction=0.9f;
	final boolean shrinkSubsets=true; //False gives very poor results, presumably due to reduced variety
	
	//For BBMerge ntriage=0.002, ptriage=0.0001 is good (for 20k samples per subset).
	float positiveTriage=0.0001f;
	float negativeTriage=0.0005f;
	int startTriage=-1;
	float startTriageMult=0.2f;
	
	int minWeightEpoch=-1;
	float minWeightEpochMult=0.03f;
	
	/*--------------------------------------------------------------*/

	boolean forcePrintFPR=false, forcePrintFNR=false;
	final boolean printFPR, printFNR;
	boolean printTPR=false, printTNR=false, printCutoff=true;
	boolean printAlpha=false, printAnneal=false;
	
	double annealStrength;
	double annealDropoff;
	

	final double alphaIncrease;//=alphaZero*(alphaMult-1.0)/(peakAlphaEpoch);
	final int alphaEpochs;//=(Tools.min(maxEpochs, 800000)-peakAlphaEpoch);
	final double alphaDropoff;//=1.0/Math.pow(alphaMult2, 1.0/alphaEpochs);
	
	/*--------------------------------------------------------------*/
	
	int jobsPerEpoch;
	int networksPerCycle=1;
	int cycles=1;
	long networksTested=0;
	
	private LongList seedList;
	private int seedsToScan=-1;
	private float seedsToScanMult=0;//10 worked best for BBMerge; needs more testing;
	private int scanEpochs=1024;
	private int scanSamples=80000;
	
	/*--------------------------------------------------------------*/
	
	//Cycle Fields
	int lastPrintEpoch=-1;
	
	private CellNet finalNet=null;
	private CellNet bestNet=null;
	
	/*--------------------------------------------------------------*/
	
	int printInterval=10000;//TODO: Lower is slower but better; decouple print and evaluate
	boolean quiet=false;
	boolean printStatus=true;
	boolean machineOut=false;
	int printMode=printAverage;
	long lastPrintTime=-1;
	private static final int printFirst=0, printAll=1, printBest=2, printAverage=3;

	HashMap<CellNet, CellNet> networkMap=new HashMap<CellNet, CellNet>();
	HashMap<CellNet, ArrayList<CellNet>> networkMap2=new HashMap<CellNet, ArrayList<CellNet>>();
	ArrayList<CellNet> finalNets=new ArrayList<CellNet>();
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;

	private long edgesIn=0;
	private long networksIn=0;
	
	private int maxLines=Shared.MAX_ARRAY_LEN;
	int maxLinesV=-1;
	private boolean shuffleRaw=false;
	private float balance=0.4f;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Network Input File */
	private final FileFormat ffNetIn;

	/** Data Input File */
	private final FileFormat ffDataIn;
	/** Data Input File */
	private final FileFormat ffValidateIn;
	/** Data Output File */
	private final FileFormat ffDataOut;
	/** Network Output File */
	private final FileFormat ffNetOutFinal;
	/** Intermediate Output File */
	private final FileFormat ffNetOutBest;
//	/** Optional Output File for Junk */
//	private final FileFormat ffoutInvalid;
	
//	final Random randyAnnealSeed;
	final Random randyNetSeed;
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
